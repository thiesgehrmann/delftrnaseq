delftrnaseq
===========

RNA-Seq pipeline used at the Delft Bioinformatics Lab for paired and single end reads.

![An outline of what the pipeline can (or will) be able to produce](/readme/pipeline.png)

Dependencies
=============

  * IBIDAS: https://github.com/mhulsman/ibidas
  * Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
  * STAR: https://code.google.com/p/rna-star/
  * SAMTOOLS: http://samtools.sourceforge.net/
  * CUFFLINKS: http://cufflinks.cbcb.umd.edu/

Further functionality is provided by installing:
  * TRINTITY: http://trinityrnaseq.sourceforge.net/
  * BLAST+: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
  * HTSeq: http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html
  * FastQC: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  * Matplotlib: http://matplotlib.org/
  * LATEX: http://www.latex-project.org/

Usage
=========

Given a configuration file (such as 'pipeline_conf_example.py') (most functionality is exemplified there)

```shell
> ./pipeline_conf_example.py # Produces a makefile 'Makefile', and a configuration file 'conf'
```
You can get a status of each step by issuing:
```shell
> make status
```
And you can run any one of these steps with:
```shell¸
> make TRIMMOMATIC
```

If you start a step and you have not completed the dependencies of that step, it will run those first. E.g.:

Installation
===============

Clone the delftrnaseq repository to somewhere:

```shell
> git clone https://github.com/thiesgehrmann/delftrnaseq.git
```

In each configuration file you create, the first line should reflect the location of 'pipeline.py' in your current installation.
