import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.List;

public class JAlignRegions {

    private static final int BIG_DISTANCE = 2500;

    private static class Alignment {
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

    private static List<Alignment> seqAlign(int[] x, int[] y, IDistanceFunction distFunc) {
        int[][] distances = new int[x.length + 1][y.length + 1];
        int[][] pathX = new int[x.length + 1][y.length + 1];
        int[][] pathY = new int[x.length + 1][y.length + 1];
        List<Alignment> ralign = new ArrayList<Alignment>();

        for (int j = 0; j <= y.length; j++) {
            for (int i = 0; i < x.length; i++) {
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

        return ralign;
    }

    private interface IDistanceFunction {
        int calculate(int x1, int y1, int x2, int y2);
    }
    
    private static double pnorm(double z) {
        double t = 1 / (1 + 0.2316419 * z);
        double pd = 1 - 0.3989423 * Math.exp(-z * z / 2)
                * ((((1.330274429 * t - 1.821255978) * t + 1.781477937) * t - 0.356563782) * t + 0.319381530) * t;
        return pd;
    }

    private static int match(int len1, int len2) {
        // foreign characters per English character
        double foreignCharsPerEngChar = 1d;
        // variance per English character
        double varPerEngChar = 6.8d;

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

    private static final IDistanceFunction TWO_SIDE_DISTANCE = new IDistanceFunction() {
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

    private static List<String> readlines(String filename) throws IOException {
        List<String> lines = new ArrayList<String>();
        FileInputStream fis = null;
        InputStreamReader isr = null;
        BufferedReader br = null;
        try {
            fis = new FileInputStream(filename);
            isr = new InputStreamReader(fis, "UTF-8");
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

    private static int lengthOfARegion(List<String> lines) {
        int result = 0;
        
        for (String line : lines) {
            result += line.codePointCount(0, line.length());
        }
        return result;
    }

    private static int[] regionLengths(List<List<String>> regions) {
        int[] result = new int[regions.size()];

        for (int i = 0; i < regions.size(); i++) {
            result[i] = lengthOfARegion(regions.get(i));
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

    public static void main(String[] args) throws Exception {
        boolean verbose = false;
        boolean debug = false;
        String softDelimiter = null;
        String hardDelimiter = null;
        String filename1 = null;
        String filename2 = null;
        for (int i = 1; i < args.length; i++) {
            String arg = args[i];
            if ("-v".equals(arg)) {
                verbose = true;
            } else if ("-V".equals(arg)) {
                verbose = true;
                debug = true;
            } else if ("-d".equals(arg)) {
                softDelimiter = args[i++];
            } else if ("-D".equals(arg)) {
                hardDelimiter = args[i++];
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

        List<String> lines1 = readlines(filename1);
        List<String> lines2 = readlines(filename2);

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

            int[] len1 = regionLengths(softRegions1);
            int[] len2 = regionLengths(softRegions2);

            List<Alignment> align = seqAlign(len1, len2, TWO_SIDE_DISTANCE);

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

                if (debug) {
                    String out = "n=" + align.size() + " i=" + i + " x1=" + a.x1 + " y1=" + a.y1 + " x2=" + a.x2
                            + " y2=" + a.y2;
                    out1.println(out);
                    out2.println(out);
                }
                if (verbose) {
                    String out = ".Score " + a.d;
                    out1.println(out);
                    out2.println(out);
                }

                for (; prevx < ix; prevx++) {
                    if (debug) {
                        out1.println("ix=" + ix + " prevx=" + prevx);
                    }
                    printRegion(out1, softRegions1.get(prevx), a.d);
                }
                out1.println(softDelimiter);

                for (; prevy < iy; prevy++) {
                    if (debug) {
                        out1.println("iy=" + iy + " prevy=" + prevy);
                    }
                    printRegion(out2, softRegions2.get(prevy), a.d);
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
