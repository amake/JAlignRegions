import sys
import os
from codecs import open

def main(*files):
    for in_file in files:
        with open(in_file, encoding='utf-8') as inpt:
            lines = inpt.read().splitlines()
        out_file = os.path.basename(in_file) + '.out'
        print 'Outputting to:', out_file
        with open(out_file, 'w', encoding='utf-8') as output:
            for line in lines:
                if line.startswith('#'):
                    output.write('.PARA\n')
                else:
                    for word in line.split():
                        output.write(word + '\n')
                    output.write('.End of Sentence\n')

if __name__ == '__main__':
    main(sys.argv[1:])