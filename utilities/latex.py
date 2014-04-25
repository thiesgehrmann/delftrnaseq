import os.path
from ibidas import *
import numpy as np

class LatexFile(object):
    def __init__(self, filename) :
        self.f = open(filename, 'w')
        dirsplit = os.path.split(filename)
        self.basedir = dirsplit[0]
        self.write_header()

    def write_header(self) :
        self.f.write(r"""
\documentclass[a4paper,10pt]{article}
\usepackage[utf8x]{inputenc}
\usepackage{graphicx}
\usepackage{longtable}
\usepackage{booktabs}
\usepackage{amsmath}
\usepackage{verbatim}
\usepackage{hyperref}
\usepackage{morefloats}
\usepackage[top=3cm, bottom=3cm, left=2cm, right=2cm]{geometry}
\newcommand{\ra}[1]{\renewcommand{\arraystretch}{#1}}
""")

    
    def write_title(self, title, author) :
        self.f.write("\\title{%s}\n" % title)
        self.f.write("\\author{%s}\n" % author)
        self.f.write("\\begin{document}\n")
        self.f.write("\\maketitle\n\n\n")
        self.f.write("\\tableofcontents\n\n");
                                             
    sectionlevels = ['section','subsection','subsubsection','paragraph','subparagraph']
    def start_section(self, sectionname, level=0) :
        s = self.sectionlevels[level]        
        self.f.write("\n\n\n\\%s{%s}\n" % (s, sectionname))

    def texcape(self, txt) :
        return txt.replace("_","\_")

    def add_text(self, text, texcape = True) :
        if texcape :
            self.f.write("%s\n" % self.texcape(text));
        else :
            self.f.write("%s\n" % text);

    def write_table(self, columns, data, caption, alignment = None) :
        n_columns = len(columns)
        self.f.write('\\begin{table}[h]\n')
        self.f.write('\\caption{%s}\n' % self.texcape(caption))
        self.f.write('\\centering\n')
        self.f.write('\\ra{1.2}\n')
        if alignment is None :
            self.f.write('\\begin{tabular}{%s}\n' % ("r" * n_columns))
        else :
            self.f.write('\\begin{tabular}{%s}\n' % alignment)
        self.f.write('\\toprule\n')
        for i in xrange(n_columns) :
            col = columns[i]
            if isinstance(col, int) :
                str_rep = '$%s$' % format(col, ',d').replace(',', '\,')
            elif isinstance(col, float) :
                str_rep = '$%.2g$' % col
            else :
                str_rep = self.texcape(str(col))
            self.f.write('%s' % str_rep)
            if i != n_columns - 1 :
                self.f.write(' & ')
        self.f.write('\\tabularnewline\n')
        self.f.write('\\midrule\n')
        data = np.atleast_2d(data)
        for i in xrange(data.shape[0]) :
            for j in xrange(n_columns) :
                obj = data[i][j]
                if isinstance(obj, int) :
                    str_rep = '$%s$' % format(obj, ',d').replace(',', '\,')
                elif isinstance(obj, float) :
                    str_rep = '$%.2g$' % obj
                else :
                    str_rep = self.texcape(str(obj))
                self.f.write('%s' % str_rep)
                if j < n_columns - 1 :
                    self.f.write(' & ')
            self.f.write('\\tabularnewline\n')
        self.f.write('\\bottomrule\n')
        self.f.write('\\end{tabular}\n')
        self.f.write('\\end{table}\n')

    def write_rep(self, R, caption) :
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
        for row in D :
            self.f.write(" & ".join([ self.texcape(str(elem)) for elem in row]) + ' \\\\ \n');
        #efor
        self.f.write("\\end{longtable}\n");
        self.f.write("\\end{center}\n\n");

    def include_figure(self, filename, refname, caption = None, width = 1.0, additional_options = None) :
        self.f.write("\\begin{figure}[htp]\n")
        self.f.write("\\noindent\\makebox[\\textwidth]{%\n")
        if os.path.isabs(filename) :
            filename = os.path.relpath(filename, self.basedir)
        if additional_options is None :
            self.f.write('\\includegraphics[width=%g\\textwidth]{%s}}\n' % (width, filename))
        else :
            self.f.write('\\includegraphics[%s,width=%g\\textwidth]{%s}}\n' % (additional_options, width, filename))
        
        if not caption is None :
            self.f.write("\\caption{%s}\n" % self.texcape(caption))
        
        self.f.write("\\label{%s}\n" % refname)
        self.f.write("\\end{figure}\n")

    def clear_page(self) :
        self.f.write("\\clearpage\n\n")

    def end_document(self) :
        self.f.write("\\end{document}\n")
        self.f.close()

