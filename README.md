delftrnaseq
===========

RNA-Seq pipeline used at the Delft Bioinformatics Lab.

![An outline of what the pipeline can (or will) be able to produce](/delftrnaseq.png)

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
  * Matplotlib: http://matplotlib.org/
  * LATEX: http://www.latex-project.org/

Usage
=========

Given a configuration file (such as 'pipeline_conf_example.py'):

```shell
> ./pipeline_conf_example.py # Produces a makefile 'Makefile', and a configuration file 'conf'
```
You can get a status of each step by issuing:
```shell
> make status
|[TRIMMOMATIC]:         INCOMPLETE
[STAR_GG]:             INCOMPLETE
[STAR_AL]:             INCOMPLETE
[POST_STAR_AL_BAM]:    INCOMPLETE
[POST_STAR_AL_SORT]:   INCOMPLETE
[PRE_CUFFLINKS_MERGE]: INCOMPLETE
[PRE_CUFFLINKS_SORT]:  INCOMPLETE
[GENOME_ANNOT_FORMAT]: INCOMPLETE
[CUFFLINKS]:           INCOMPLETE
[CUFFDIFF]:            INCOMPLETE
[CUFFLINKS_INDIV]:     INCOMPLETE
[TRINITY]:             INCOMPLETE
[TRINITY_ORF]:         INCOMPLETE
[UNMAPPED_BLAST]:      INCOMPLETE
[UNMAPPED]:            INCOMPLETE
[CUFFDIFF_COMBINE]:    INCOMPLETE
[QUALITYREPORT]:       INCOMPLETE
[POSTANALYSIS]:        INCOMPLETE
```
And you can run any one of these steps with:
```shel
> make TRIMMOMATIC
> make status
[TRIMMOMATIC]:         COMPLETE
[STAR_GG]:             INCOMPLETE
[STAR_AL]:             INCOMPLETE
[POST_STAR_AL_BAM]:    INCOMPLETE
[POST_STAR_AL_SORT]:   INCOMPLETE
[PRE_CUFFLINKS_MERGE]: INCOMPLETE
[PRE_CUFFLINKS_SORT]:  INCOMPLETE
[GENOME_ANNOT_FORMAT]: INCOMPLETE
[CUFFLINKS]:           INCOMPLETE
[CUFFDIFF]:            INCOMPLETE
[CUFFLINKS_INDIV]:     INCOMPLETE
[TRINITY]:             INCOMPLETE
[TRINITY_ORF]:         INCOMPLETE
[UNMAPPED_BLAST]:      INCOMPLETE
[UNMAPPED]:            INCOMPLETE
[CUFFDIFF_COMBINE]:    INCOMPLETE
[QUALITYREPORT]:       INCOMPLETE
[POSTANALYSIS]:        INCOMPLETE
```

If you start a step and you have not completed the dependencies of that step, it will run those first:
```shell
> make CUFFDIFF
> make status
[TRIMMOMATIC]:         COMPLETE
[STAR_GG]:             COMPLETE
[STAR_AL]:             COMPLETE
[POST_STAR_AL_BAM]:    COMPLETE
[POST_STAR_AL_SORT]:   COMPLETE
[PRE_CUFFLINKS_MERGE]: COMPLETE
[PRE_CUFFLINKS_SORT]:  COMPLETE
[GENOME_ANNOT_FORMAT]: COMPLETE
[CUFFLINKS]:           COMPLETE
[CUFFDIFF]:            COMPLETE
[CUFFLINKS_INDIV]:     INCOMPLETE
[TRINITY]:             INCOMPLETE
[TRINITY_ORF]:         INCOMPLETE
[UNMAPPED_BLAST]:      INCOMPLETE
[UNMAPPED]:            INCOMPLETE
[CUFFDIFF_COMBINE]:    INCOMPLETE
[QUALITYREPORT]:       INCOMPLETE
[POSTANALYSIS]:        INCOMPLETE
```
