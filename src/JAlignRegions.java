import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/*
 * usage
 *
 * align_regions -D '.PARA' -d '.End of Sentence'  file1 file2
 *
 * outputs two files: file1.al & file2.al
 *
 * regions are delimited by the -D and -d args
 *
 * the program is allowed to delete -d delimiters as necessary in order
 * align the files, but it cannot change -D delimiters.
 */
public class JAlignRegions {

    private static final int BIG_DISTANCE = 2500;

    // foreign characters per English character
    public static final double DEFAULT_FOREIGN_CHARS_PER_ENG_CHAR = 1d;
    // variance per English character
    public static final double DEFAULT_VAR_PER_ENG_CHAR = 6.8d;

    private final double foreignCharsPerEngChar;
    private final double varPerEngChar;
    private final Charset encoding;
    /**
     * Affects how the length of a region is determined.
     * <ul>
     * <li>true: Counts text length by number of bytes in the system default
     * encoding. This is how the C version does it.
     * <li>false: Counts text length by number of codepoints. This seems more
     * faithful to the intention behind the algorithm.
     * </ul>
     */
    private final boolean countBytes;

    public JAlignRegions() {
        this(DEFAULT_FOREIGN_CHARS_PER_ENG_CHAR, DEFAULT_VAR_PER_ENG_CHAR);
    }

    public JAlignRegions(double foreignCharsPerEngChar, double varPerEngChar) {
        this(foreignCharsPerEngChar, varPerEngChar, Charset.defaultCharset(), false);
    }

    public JAlignRegions(double foreignCharsPerEngChar, double varPerEngChar, Charset encoding, boolean countBytes) {
        this.foreignCharsPerEngChar = foreignCharsPerEngChar;
        this.varPerEngChar = varPerEngChar;
        this.encoding = encoding;
        this.countBytes = countBytes;
    }

    public static class Alignment {
        public final int x1;
        public final int y1;
        public final int x2;
        public final int y2;
        public final int d;

        public Alignment(int x1, int y1, int x2, int y2, int d) {
            this.x1 = x1;
            this.y1 = y1;
            this.x2 = x2;
            this.y2 = y2;
            this.d = d;
        }
    }

    /*
     * seq_align by Mike Riley
     * Sequence alignment routine.
     * This version allows for contraction/expansions.
     * 
     * x and y are sequences of objects, represented as non-zero ints, to be
     * aligned.
     * 
     * dist_funct(x1, y1, x2, y2) is a distance function of 4 args:
     * 
     *   dist_funct(x1, y1, 0, 0) gives cost of substitution of x1 by y1.
     *   dist_funct(x1, 0, 0, 0) gives cost of deletion of x1.
     *   dist_funct(0, y1, 0, 0) gives cost of insertion of y1.
     *   dist_funct(x1, y1, x2, 0) gives cost of contraction of (x1,x2) to y1.
     *   dist_funct(x1, y1, 0, y2) gives cost of expansion of x1 to (y1,y2).
     *   dist_funct(x1, y1, x2, y2) gives cost to match (x1,x2) to (y1,y2).
     * 
     * align is the alignment, with (align[i].x1, align[i].x2) aligned with
     * (align[i].y1, align[i].y2). Zero in align[].x1 and align[].y1 correspond
     * to insertion and deletion, respectively. Non-zero in align[].x2 and
     * align[].y2 correspond to contraction and expansion, respectively.
     * align[].d gives the distance for that pairing.
     * 
     * The function returns the length of the alignment.
     */
    private static List<Alignment> seqAlign(int[] x, int[] y, IDistanceFunction distFunc) {
        int[][] distances = new int[x.length + 1][y.length + 1];
        int[][] pathX = new int[x.length + 1][y.length + 1];
        int[][] pathY = new int[x.length + 1][y.length + 1];
        List<Alignment> ralign = new ArrayList<Alignment>();

        for (int j = 0; j <= y.length; j++) {
            for (int i = 0; i <= x.length; i++) {
                int d1 = i > 0 && j > 0 ? // substitution
                        distances[i - 1][j - 1] + distFunc.calculate(x[i - 1], y[j - 1], 0, 0)
                        : Integer.MAX_VALUE;
                int d2 = i > 0 ? // deletion
                        distances[i - 1][j] + distFunc.calculate(x[i - 1], 0, 0, 0)
                        : Integer.MAX_VALUE;
                int d3 = j > 0 ? // insertion
                        distances[i][j - 1] + distFunc.calculate(0, y[j - 1], 0, 0)
                        : Integer.MAX_VALUE;
                int d4 = i > 1 && j > 0 ? // contraction
                        distances[i - 2][j - 1] + distFunc.calculate(x[i - 2], y[j - 1], x[i - 1], 0)
                        : Integer.MAX_VALUE;
                int d5 = i > 0 && j > 1 ? // expansion
                        distances[i - 1][j - 2] + distFunc.calculate(x[i - 1], y[j - 2], 0, y[j - 1])
                        : Integer.MAX_VALUE;
                int d6 = i > 1 && j > 1 ? // melding
                        distances[i - 2][j - 2] + distFunc.calculate(x[i - 2], y[j - 2], x[i - 1], y[j - 1])
                        : Integer.MAX_VALUE;

                int dmin = d1;
                if (d2 < dmin) {
                    dmin = d2;
                }
                if (d3 < dmin) {
                    dmin = d3;
                }
                if (d4 < dmin) {
                    dmin = d4;
                }
                if (d5 < dmin) {
                    dmin = d5;
                }
                if (d6 < dmin) {
                    dmin = d6;
                }

                if (dmin == Integer.MAX_VALUE) {
                    distances[i][j] = 0;
                } else if (dmin == d1) {
                    distances[i][j] = d1;
                    pathX[i][j] = i - 1;
                    pathY[i][j] = j - 1;
                } else if (dmin == d2) {
                    distances[i][j] = d2;
                    pathX[i][j] = i - 1;
                    pathY[i][j] = j;
                } else if (dmin == d3) {
                    distances[i][j] = d3;
                    pathX[i][j] = i;
                    pathY[i][j] = j - 1;
                } else if (dmin == d4) {
                    distances[i][j] = d4;
                    pathX[i][j] = i - 2;
                    pathY[i][j] = j - 1;
                } else if (dmin == d5) {
                    distances[i][j] = d5;
                    pathX[i][j] = i - 1;
                    pathY[i][j] = j - 2;
                } else { // dmin == d6
                    distances[i][j] = d6;
                    pathX[i][j] = i - 2;
                    pathY[i][j] = j - 2;
                }
            }
        }
        
        for (int oi, oj, i = x.length, j = y.length; i > 0 || j > 0; i = oi, j = oj) {
            oi = pathX[i][j];
            oj = pathY[i][j];
            int di = i - oi;
            int dj = j - oj;
            
            if (di == 1 && dj == 1) { // substitution
                ralign.add(new Alignment(x[i - 1],
                        y[j - 1],
                        0,
                        0,
                        distances[i][j] - distances[i - 1][j - 1]));
            } else if (di == 1 && dj == 0) { // deletion
                ralign.add(new Alignment(x[i - 1],
                        0,
                        0,
                        0,
                        distances[i][j] - distances[i - 1][j]));
            } else if (di == 0 && dj == 1) { // insertion
                ralign.add(new Alignment(0,
                        y[j - 1],
                        0,
                        0,
                        distances[i][j] - distances[i][j - 1]));
            } else if (dj == 1) { // contraction
                ralign.add(new Alignment(x[i - 2],
                        y[j - 1],
                        x[i - 1],
                        0,
                        distances[i][j] - distances[i - 2][j - 1]));
            } else if (di == 1) { // expansion
                ralign.add(new Alignment(x[i - 1],
                        y[j - 2],
                        0,
                        y[j - 1],
                        distances[i][j] - distances[i - 1][j - 2]));
            } else { // di == 2 && dj == 2 : melding
                ralign.add(new Alignment(x[i - 2],
                        y[j - 2],
                        x[i - 1],
                        y[j - 1],
                        distances[i][j] - distances[i - 2][j - 2]));
            }
        }

        Collections.reverse(ralign);
        return ralign;
    }

    private interface IDistanceFunction {
        int calculate(int x1, int y1, int x2, int y2);
    }
    
    /*
     * Returns the area under a normal distribution from -inf to z standard
     * deviations
     */
    private static double pnorm(double z) {
        double t = 1 / (1 + 0.2316419 * z);
        double pd = 1 - 0.3989423 * Math.exp(-z * z / 2)
                * ((((1.330274429 * t - 1.821255978) * t + 1.781477937) * t - 0.356563782) * t + 0.319381530) * t;
        /* see Gradsteyn & Rhyzik, 26.2.17 p932 */
        return pd;
    }

    /*
     * Return -100 * log probability that an English sentence of length len1 is
     * a translation of a foreign sentence of length len2. The probability is
     * based on two parameters, the mean and variance of number of foreign
     * characters per English character.
     */
    private int match(int len1, int len2) {

        if (len1 == 0 && len2 == 0) {
            return 0;
        }
        double mean = (len1 + len2 / foreignCharsPerEngChar) / 2;
        double z = (foreignCharsPerEngChar * len1 - len2) / Math.sqrt(varPerEngChar * mean);

        // Need to deal with both sides of the normal distribution
        if (z < 0) {
            z = -z;
        }
        double pd = 2 * (1 - pnorm(z));

        if (pd > 0) {
            return (int) (-100 * Math.log(pd));
        } else {
            return BIG_DISTANCE;
        }
    }

    private final IDistanceFunction TWO_SIDE_DISTANCE = new IDistanceFunction() {
        @Override
        public int calculate(int x1, int y1, int x2, int y2) {
            int penalty21 = 230; // -100 * log([prob of 2-1 match] / [prob of 1-1 match])
            int penalty22 = 440; // -100 * log([prob of 2-2 match] / [prob of 1-1 match])
            int penalty01 = 450; // -100 * log([prob of 0-1 match] / [prob of 1-1 match])
            
            if (x2 == 0 && y2 == 0) {
                if (x1 == 0) { // insertion
                    return match(x1, y1) + penalty01;
                } else if (y1 == 0) { // deletion
                    return match(x1, y1) + penalty01;
                } else { // substitution
                    return match(x1, y1);
                }
            } else if (x2 == 0) { // expansion
                return match(x1, y1 + y2) + penalty21;
            } else if (y2 == 0) { // contraction
                return match(x1 + x2, y1) + penalty21;
            } else { // melding
                return match(x1 + x2, y1 + y2) + penalty22;
            }
        }
    };

    /*
     * return an array of strings, one string for each line of the file
     */
    private List<String> readlines(String filename) throws IOException {
        List<String> lines = new ArrayList<String>();
        FileInputStream fis = null;
        InputStreamReader isr = null;
        BufferedReader br = null;
        try {
            fis = new FileInputStream(filename);
            isr = new InputStreamReader(fis, encoding);
            br = new BufferedReader(isr);
            String line;
            while ((line = br.readLine()) != null) {
                lines.add(line);
            }
            fis.close();
            isr.close();
            br.close();
        } finally {
            if (fis != null) {
                fis.close();
            }
            if (isr != null) {
                isr.close();
            }
            if (br != null) {
                br.close();
            }
        }
        return lines;
    }

    private static void printRegion(PrintWriter out, List<String> lines, int score) {
        for (String line : lines) {
            out.println(line);
        }
    }

    private int lengthOfARegion(List<?> lines) {
        int result = lines.size();
        
        for (Object o : lines) {
            result += lengthOfAString(o.toString());
        }
        return result;
    }

    private int lengthOfARegion(String line) {
        // This count will differ from the tokenized count by virtue of
        // including spaces instead of the number of tokens.
        return lengthOfAString(line);
    }

    private int lengthOfAString(String s) {
        return countBytes ? encoding.encode(s).limit() : s.codePointCount(0, s.length());
    }

    private int[] regionLengths(List<?> regions) {
        int[] result = new int[regions.size()];

        for (int i = 0; i < regions.size(); i++) {
            Object region = regions.get(i);
            int len;
            if (region instanceof String) {
                len = lengthOfARegion((String) region);
            } else if (region instanceof List<?>) {
                len = lengthOfARegion((List<?>) region);
            } else {
                throw new IllegalArgumentException("Unsupported region type: " + region.getClass().getName());
            }
            result[i] = len;
        }
        return result;
    }

    private static List<List<String>> findSubRegions(List<String> lines, String delimiter) {
        List<List<String>> result = new ArrayList<List<String>>();

        List<String> temp = new ArrayList<String>();
        for (String line : lines) {
            if (delimiter.equals(line)) {
                result.add(temp);
                temp = new ArrayList<String>();
            } else {
                temp.add(line);
            }
        }

        return result;
    }

    public static class Bead<T> {
        public final Alignment alignment;
        public final List<T> sourceLines;
        public final List<T> targetLines;
        private final List<String> sourceDebug;
        private final List<String> targetDebug;

        public Bead(Alignment alignment, List<T> sourceLines, List<T> targetLines,
                List<String> sourceDebug, List<String> targetDebug) {
            this.alignment = alignment;
            this.sourceLines = Collections.unmodifiableList(sourceLines);
            this.targetLines = Collections.unmodifiableList(targetLines);
            this.sourceDebug = Collections.unmodifiableList(sourceDebug);
            this.targetDebug = Collections.unmodifiableList(targetDebug);
        }
    }

    public <T> List<Bead<T>> align(List<T> softRegions1, List<T> softRegions2) {
        int[] len1 = regionLengths(softRegions1);
        int[] len2 = regionLengths(softRegions2);

        List<Alignment> align = seqAlign(len1, len2, TWO_SIDE_DISTANCE);
        List<Bead<T>> result = new ArrayList<Bead<T>>();

        int prevx = 0, prevy = 0, ix = 0, iy = 0;
        for (int i = 0; i < align.size(); i++) {
            Alignment a = align.get(i);
            if (a.x2 > 0) {
                ix++;
            } else if (a.x1 == 0) {
                ix--;
            }
            if (a.y2 > 0) {
                iy++;
            } else if (a.y1 == 0) {
                iy--;
            }
            if (a.x1 == 0 && a.y1 == 0 && a.x2 == 0 && a.y2 == 0) {
                ix++;
                iy++;
            }
            ix++;
            iy++;

            List<T> sourceLines = new ArrayList<T>();
            List<T> targetLines = new ArrayList<T>();
            List<String> sourceDebug = new ArrayList<String>();
            List<String> targetDebug = new ArrayList<String>();

            for (; prevx < ix; prevx++) {
                sourceDebug.add("ix=" + ix + " prevx=" + prevx);
                sourceLines.add(softRegions1.get(prevx));
            }

            for (; prevy < iy; prevy++) {
                targetDebug.add("iy=" + iy + " prevy=" + prevy);
                targetLines.add(softRegions2.get(prevy));
            }

            result.add(new Bead<T>(a, sourceLines, targetLines, sourceDebug, targetDebug));
        }

        return result;
    }

    public static void main(String[] args) throws Exception {
        boolean verbose = false;
        boolean debug = false;
        String softDelimiter = null;
        String hardDelimiter = null;
        String filename1 = null;
        String filename2 = null;
        for (int i = 0; i < args.length; i++) {
            String arg = args[i];
            if ("-v".equals(arg)) {
                verbose = true;
            } else if ("-V".equals(arg)) {
                verbose = true;
                debug = true;
            } else if ("-d".equals(arg)) {
                softDelimiter = args[++i];
            } else if ("-D".equals(arg)) {
                hardDelimiter = args[++i];
            } else if (filename1 == null) {
                filename1 = arg;
            } else if (filename2 == null) {
                filename2 = arg;
            } else {
                System.err.print("usage: java JAlignRegions [d (soft delimiter)] [D (hard delimiter)]");
                System.exit(2);
            }
        }

        if (softDelimiter == null || hardDelimiter == null || filename1 == null || filename2 == null) {
            err("wrong number of arguments");
        }
        
        String outFilename1 = filename1 + ".al";
        PrintWriter out1 = null;
        try {
            out1 = new PrintWriter(outFilename1);
        } catch (Exception ex) {
            System.err.println("can't open " + outFilename1);
            System.exit(2);
        }
        
        String outFilename2 = filename2 + ".al";
        PrintWriter out2 = null;
        try {
            out2 = new PrintWriter(outFilename2);
        } catch (Exception ex) {
            System.err.println("can't open " + outFilename2);
            System.exit(2);
        }

        JAlignRegions aligner = new JAlignRegions();

        List<String> lines1 = aligner.readlines(filename1);
        List<String> lines2 = aligner.readlines(filename2);

        List<List<String>> hardRegions1 = findSubRegions(lines1, hardDelimiter);
        List<List<String>> hardRegions2 = findSubRegions(lines2, hardDelimiter);

        if (hardRegions1.size() != hardRegions2.size()) {
            System.err.println("align_regions: input files do not contain the same number " + "of hard regions ("
                    + hardDelimiter + ").");
            System.err.println(filename1 + " has " + hardRegions1.size() + " and " + filename2 + " has "
                    + hardRegions2.size() + ".");
            System.exit(2);
        }

        for (int j = 0; j < hardRegions1.size(); j++) {
            List<List<String>> softRegions1 = findSubRegions(hardRegions1.get(j), softDelimiter);
            List<List<String>> softRegions2 = findSubRegions(hardRegions2.get(j), softDelimiter);
            
            if (debug) {
                out1.println("number of soft regions=" + softRegions1.size());
                out2.println("number of soft regions=" + softRegions2.size());
            }

            List<Bead<List<String>>> beads = aligner.align(softRegions1, softRegions2);
            for (int i = 0; i < beads.size(); i++) {
                Bead<List<String>> bead = beads.get(i);

                if (debug) {
                    String out = "n=" + beads.size() + " i=" + i + " x1=" + bead.alignment.x1 + " y1="
                            + bead.alignment.y1 + " x2=" + bead.alignment.x2 + " y2=" + bead.alignment.y2;
                    out1.println(out);
                    out2.println(out);
                }
                if (verbose) {
                    String out = ".Score " + bead.alignment.d;
                    out1.println(out);
                    out2.println(out);
                }

                for (int k = 0; k < bead.sourceLines.size(); k++) {
                    if (debug) {
                        out1.print(bead.sourceDebug.get(k) + " ");
                    }
                    printRegion(out1, bead.sourceLines.get(k), bead.alignment.d);
                }
                out1.println(softDelimiter);

                for (int k = 0; k < bead.targetLines.size(); k++) {
                    if (debug) {
                        out2.print(bead.targetDebug.get(k) + " ");
                    }
                    printRegion(out2, bead.targetLines.get(k), bead.alignment.d);
                }
                out2.println(softDelimiter);
            }
            out1.println(hardDelimiter);
            out2.println(hardDelimiter);
        }

        if (out1.checkError()) {
            System.err.println("Problem writing output file: " + outFilename1);
            System.exit(2);
        }
        if (out2.checkError()) {
            System.err.println("Problem writing output file: " + outFilename2);
            System.exit(2);
        }

        out1.close();
        out2.close();
    }

    private static void err(String msg) {
        System.err.println("**ERROR**: " + msg);
        System.exit(2);
    }
}
