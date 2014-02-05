from os.path import basename
from ibidas import *;

class LatexFile(object):
    def __init__(self, filename):
        self.f = open(filename, 'w')
        self.write_header()

    def write_header(self):
        self.f.write(r"""
\documentclass[a4paper,10pt]{article}
\usepackage[utf8x]{inputenc}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage[top=3cm, bottom=3cm, left=2cm, right=2cm]{geometry}
""")

    
    def write_title(self, title, author):
        self.f.write("\\title{%s}\n" % title)
        self.f.write("\\author{%s}\n" % author)
        self.f.write("\\begin{document}\n")
        self.f.write("\\maketitle\n\n\n")
        self.f.write("\\tableofcontents\n\n");
                                             
    sectionlevels = ['section','subsection','subsubsection','paragraph','subparagraph']
    def start_section(self, sectionname, level=0):
        s = self.sectionlevels[level]        
        self.f.write("\n\n\n\\%s{%s}\n" % (s, sectionname))

    def texcape(self, txt):
        return txt.replace("_","\_")

    def add_text(self, text):
        self.f.write("%s\n" % self.texcape(text));

    def write_rep(self, R, caption):
        D = zip(*R());
        cols = 'l'*len(R.Names);
        self.f.write("\\begin{center}\n");
        self.f.write("\\begin{longtable}{%s}\n" % cols);
        self.f.write("\\caption{%s} \\\\ \\hline\n\n" % (self.texcape(caption)));
        self.f.write(" & ".join(["\\textbf{%s}" % self.texcape(n) for n in R.Names]) + '\\\\ \\hline\n');
        self.f.write("\\endfirsthead\n\n");
        self.f.write("\\multicolumn{%d}{c}{{\\bfseries \\tablename \\thetable{} -- continued from previous page}} \\\\ \\hline\n" % len(R.Names));
        self.f.write(" & ".join(["\\textbf{%s}" % self.texcape(n) for n in R.Names]) + '\\\\ \\hline\n');
        self.f.write("\\endhead\n\n");
        self.f.write("\\multicolumn{%d}{c}{{\\bfseries \\tablename \\thetable{} -- continued on next page}} \\\\ \\hline\n" % len(R.Names));
        self.f.write("\\endfoot\n\n");
        self.f.write("\\hline \\hline\n");
        self.f.write("\\endlastfoot\n\n");
        for row in D:
            self.f.write(" & ".join([ self.texcape(str(elem)) for elem in row]) + ' \\\\ \n');
        #efor
        self.f.write("\\end{longtable}\n");
        self.f.write("\\end{center}\n\n");

    def include_figure(self, filename, refname, caption = None, width=1.0):
        self.f.write("\\begin{figure}[htp]\n")
        self.f.write("\\noindent\\makebox[\\textwidth]{%\n")
        filename = basename(filename)
        self.f.write('\\includegraphics[width=%g\\textwidth]{%s}}\n' % (width, filename))
        
        if not caption is None:
            self.f.write("\\caption{%s}\n" % self.texcape(caption))
        
        self.f.write("\\label{%s}\n" % refname)
        self.f.write("\\end{figure}\n")

    def clear_page(self):
        self.f.write("\\clearpage\n\n")

    def end_document(self):
        self.f.write("\\end{document}\n")
        self.f.close()

