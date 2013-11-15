from os.path import basename
class LatexFile(object):
    def __init__(self, filename):
        self.f = open(filename, 'w')
        self.write_header()

    def write_header(self):
        self.f.write(r"""
\documentclass[a4paper,10pt]{article}
\usepackage[utf8x]{inputenc}
\usepackage{graphicx}""")

    
    def write_title(self, title, author):
        self.f.write("\\title{%s}\n" % title)
        self.f.write("\\author{%s}\n" % author)
        self.f.write("\\begin{document}\n")
        self.f.write("\\maketitle")
                                             
    sectionlevels = ['section','subsection','subsubsection','paragraph','subparagraph']
    def start_section(self, sectionname, level=0):
        s = self.sectionlevels[level]        
        self.f.write("\n\\%s{%s}\n" % (s, sectionname))

    def texcape(self, txt):
        return txt.replace("_","\_")

    def include_figure(self, filename, refname, caption = None, width=None):
        self.f.write("\\begin{figure}[htp]\n")
        self.f.write("\\noindent\\makebox[\\textwidth]{%\n")
        filename = basename(filename)
        if width is None:
            self.f.write('\\includegraphics[width=\\textwidth]{%s}}\n' % filename)
        else:
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


