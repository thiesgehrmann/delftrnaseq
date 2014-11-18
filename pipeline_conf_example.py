#!/home/nfs/thiesgehrmann/groups/w/phd/tasks/rnaseq_pipeline/delftrnaseq/pipeline.py

#### Running on Markov
####
####

class config(PIPELINECONF):

    # Specify the number of threads
  __max_threads__ = 12;

    # Emails to receive updates at
  email = "thiesgehrmann@gmail.com"
    # Job title
  jobname = "SCHCO_RNASEQ"

    # Where to store the configuration file for the pipeline
  location = "./conf";
    # Working directory (None is default)
  workdir  = None;
    # Where to output the results
  outdir   = "/data/tmp/thiesgehrmann/rnaseq_schco3";
    # Makefile name (Default = Makefile)
  #makefile="";
    # local blast DB
  blast_db = "/data/tmp/thiesgehrmann/blastdb/nt.nal";
  #server=markov

    # Paired end reads (False for Single End reads)
  PE = True;
    # Sample files structure [ (left_end1, right_end1) ... (left_endn, right_endn) ]
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

    # A name for each sample
  sample_names  = [ "4.8-1", "4.8-2", "4.8-3", "4.8-4", "wc1-1", "wc1-2", "wc1-3", "wc1-4", "wc2-1", "wc2-2", "wc2-3", "wc2-4", "hom1-1", "hom1-2", "hom1-3", "hom1-4", "hom2-1", "hom2-2", "hom2-3", "hom2-4", "fst3-1", "fst3-2", "fst3-3", "fst3-4", "fst4-1", "fst4-2", "fst4-3", "fst4-4", "bri1-1", "bri1-2", "bri1-4", "gat1-1", "gat1-2", "gat1-3", "gat1-4", "c2h2-1", "c2h2-2", "c2h2-3", "c2h2-4", "bri1-3" ];
    # Replicate groups
  sample_labels = [ 0,       0,       1,       1,       2,       2,       3,       3,       4,       4,       5,       5,       6,        6,        7,        7,        8,        8,        9,        9,        10,       10,       11,       11,       12,       12,       13,       13,       14,       14,       15,       16,       16,       17,       17,       18,       18,       19,       19,       15       ];
    # Replicate group names
  label_names   = [ "4.8-1", "4.8-2", "wc1-1", "wc1-2", "wc2-1", "wc2-2", "hom1-1", "hom1-2", "hom2-1", "hom2-2","fst3-1", "fst3-2", "fst4-1", "fst4-2", "bri1-1", "bri1-2", "gat1-1", "gat1-2", "c2h2-1", "c2h2-2" ];

    # Genome sequence file
  genome       = '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/genome.fasta';
    # Gene model file
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

    # Build a splice DB for star?
  build_splice_db = True;

  #############################################################################
  # STAR ALIGNMENT OPTIONS                                                    #
  #############################################################################
  #star_al_opts    = "--runThreadN %d --alignIntronMax 1500 --alignIntronMin 10 " % __max_threads__;
  #star_preal_opts = "--runThreadN %d --alignIntronMax 1500 --alignIntronMin 10 " % __max_threads__;

  #############################################################################
  # CUFFDIFF OPTIONS                                                          #
  #############################################################################


    # Replicate groups to compare with cuffdiff
  cuffdiff_cmp = [(1, 3), (10, 11), (0, 14), (1, 17), (1, 15), (8, 9), (0, 16), (18, 19), (0, 10), (0, 3), (1, 11), (16, 17), (2, 3), (6, 7), (12, 13), (1, 5), (0, 4), (4, 5), (1, 13), (0, 18), (0, 12), (1, 19), (14, 15), (1, 9), (0, 8), (0, 1), (0, 13), (0, 6), (1, 7), (0, 9), (0, 5), (0, 2)];

    # cuffdiff test to perform. =  cds|gene|isoform|tss
  #cuffdiff_test_type         = 'isoform';

  #cuffdiff_opts = "";
   
  #############################################################################
  # CUFFDIFF COMBINE OPTIONS                                                  #
  #############################################################################

  #############################################################################
  # ANALYSIS OPTIONS                                                          #
  #############################################################################
  #perform_quality_report     = True;
  #perform_analysis           = True;
  #check_contamination        = False;


    # Enrichment analysis annotation files
  annotation_files = [ '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/Schco3_GeneCatalog_proteins_20130812_GO_CUT.tab',  '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/IPR_annot.tsv' ];
    # Enrichmeny analysis annotation names
  annotation_names = [ 'GO' , 'IPR' ];

    # Filter the postanalysis with these filters
  analysis_filter       = [None, '/home/nfs/thiesgehrmann/groups/w/phd/data/schco3/pred_transfacs.tsv' ];
    # Name the filters
  analysis_filter_names = [ "NoFilter", 'TransFacs' ];

    # Show venn diagrams of differencially expressed gene overlaps
  analysis_venn                    = [ (16,20,24), (19,8,6), (6,27,8), (19,23,15) ];
    # Show different venn diagrams for up and down regulated genes
  analysis_venn_updown_split       = True;
  #analysis_enrichment_updown_split = True;

  #############################################################################
  # DENSE GENOME ISOFORM ANALYSIS                                             #
  #############################################################################

  isoform_dense_analysis = True;
  isoform_dense_build_splice_db = True;
#eclass
