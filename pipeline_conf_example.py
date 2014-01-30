#!/home/nfs/thiesgehrmann/groups/w/phd/tasks/rnaseq_pipeline/delftrnaseq/pipeline.py

class config(PIPELINECONF):
  #inst_loc

  __max_threads__ = 12;

  email = "thiesgehrmann@gmail.com"
  jobname = "SCHCO_RNASEQ"

  location = "./conf";
  workdir  = None;
  outdir   = "/data/tmp/thiesgehrmann/rnaseq_schco3";
  #makefile="";
  blast_db = "/data/tmp/thiesgehrmann/blastdb/nt.nal";

  # find "`pwd`" | grep R1.fastq | sort > r1.list; find "`pwd`" | grep R2.fastq | sort > r2.list; paste r1.list r2.list | awk -F'\t' '{printf "(\047%s\047, \047%s\047),\n", $1, $2;}'; rm r1.list r2.list
  samples = [
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-01_ATCACG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-01_ATCACG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-02_CGATGT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-02_CGATGT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-03_TTAGGC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-03_TTAGGC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-04_TGACCA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-04_TGACCA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-05_ACAGTG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-05_ACAGTG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-06_GCCAAT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-06_GCCAAT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-07_CAGATC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-07_CAGATC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-08_ACTTGA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-08_ACTTGA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-09_GATCAG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-09_GATCAG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-10_TAGCTT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-10_TAGCTT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-11_GGCTAC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-11_GGCTAC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-12_CTTGTA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-12_CTTGTA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-13_AGTCAA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-13_AGTCAA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-14_AGTTCC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-14_AGTTCC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-15_ATGTCA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-15_ATGTCA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-16_CCGTCC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-16_CCGTCC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-17_GTCCGC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-17_GTCCGC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-18_GTGAAA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-18_GTGAAA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-19_GTGGCC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-19_GTGGCC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-20_GTTTCG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V86ACXX_102137-20_GTTTCG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-21_CGTACG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-21_CGTACG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-22_GAGTGG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-22_GAGTGG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-23_ACTGAT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-23_ACTGAT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-24_ATTCCT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-24_ATTCCT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-25_ATCACG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-25_ATCACG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-26_CGATGT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-26_CGATGT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-27_TTAGGC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-27_TTAGGC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-28_TGACCA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-28_TGACCA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-29_ACAGTG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-29_ACAGTG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-30_GCCAAT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-30_GCCAAT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-32_CAGATC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-32_CAGATC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-33_ACTTGA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-33_ACTTGA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-34_GATCAG_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-34_GATCAG_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-35_TAGCTT_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-35_TAGCTT_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-36_GGCTAC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-36_GGCTAC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-37_CTTGTA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-37_CTTGTA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-38_AGTCAA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-38_AGTCAA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-39_AGTTCC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-39_AGTTCC_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-40_ATGTCA_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-40_ATGTCA_L007_R2.fastq'),
('/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-41_CCGTCC_L007_R1.fastq', '/data/tmp/thiesgehrmann/data/schco_rnaseq_data/RawData/C1V8EACXX_102137-41_CCGTCC_L007_R2.fastq')
]
  sample_names  = [ "4.8-1", "4.8-2", "4.8-3", "4.8-4", "wc1-1", "wc1-2", "wc1-3", "wc1-4", "wc2-1", "wc2-2", "wc2-3", "wc2-4", "hom1-1", "hom1-2", "hom1-3", "hom1-4", "hom2-1", "hom2-2", "hom2-3", "hom2-4", "fst3-1", "fst3-2", "fst3-3", "fst3-4", "fst4-1", "fst4-2", "fst4-3", "fst4-4", "bri1-1", "bri1-2", "bri1-4", "gat1-1", "gat1-2", "gat1-3", "gat1-4", "c2h2-1", "c2h2-2", "c2h2-3", "c2h2-4", "bri1-3" ];
  sample_labels = [ 0,       0,       1,       1,       2,       2,       3,       3,       4,       4,       5,       5,       6,        6,        7,        7,        8,        8,        9,        9,        10,       10,       11,       11,       12,       12,       13,       13,       14,       14,       15,       16,       16,       17,       17,       18,       18,       19,       19,       15       ];
  label_names   = [ "4.8-1", "4.8-2", "wc1-1", "wc1-2", "wc2-1", "wc2-2", "hom1-1", "hom1-2", "hom2-1", "hom2-2","fst3-1", "fst3-2", "fst4-1", "fst4-2", "bri1-1", "bri1-2", "gat1-1", "gat1-2", "c2h2-1", "c2h2-2" ];

#  samples       = samples[20:];
#  sample_names  = sample_names[20:];
#  sample_labels = sample_labels[20:];
#  label_names   = label_names[20:];
  

  genome       = '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/genome.fasta';
  genome_annot = '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/cleaned_gene_model_proteinid.gff';

  #############################################################################
  # Trimmomatic options                                                       #
  #############################################################################
  #trimmomatic_opts="";
  trimmomatic_trim="ILLUMINACLIP:/home/nfs/thiesgehrmann/groups/w/phd/data/illumina_adapters.fasta:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

  #############################################################################
  # STAR GENOME GENERATION OPTIONS                                            #
  #############################################################################
  #star_gg_opts = "";

  build_splice_db = False;

  #############################################################################
  # STAR ALIGNMENT OPTIONS                                                    #
  #############################################################################
  #star_al_opts = "";

  #############################################################################
  # CUFFDIFF OPTIONS                                                          #
  #############################################################################

    # Versus wildtype without phenotype
  cuffdiff_cmp = [ (0, i) for i in xrange(2, 20, 2) ] + \
                 [ (1, i) for i in xrange(3, 20, 2) ];

    # Versus wildtype with phenotype
  cuffdiff_cmp.extend([ (0,2), (0,3), (0,4), (0,5),  (0,8),  (0,9),  (0,12), (0,13), (0,14), (1,15), \
                        (0,1), (0,6), (1,7), (0,10), (1,11), (0,16), (1,17), (0,18), (1,19) ]);

    # T1 vs T2
  cuffdiff_cmp.extend([ (i-1, i) for i in xrange(1, 20, 2) ]);
  #cuffdiff_cmp = None;

#eclass
