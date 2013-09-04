#!./pipeline.py

class config(PIPELINECONF):
  #inst_loc

  email   = "m.hulsman@tudelft.nl";
  jobname = "FENGFENG";

  location = "./conf";
  workdir = None;
  outdir  = "/data/tmp/marchulsman/fengfeng_run/outdir"
  #makefile="";
  blast_db = '/data/tmp/marchulsman/blast_db/nt.nal'    
  


# echo '['; find | grep fastq | cut -c1-35 | sort | uniq | while read x; do x=`basename "$x"`; echo "('`pwd`/${x}1.fastq', '`pwd`/${x}2.fastq'),"; done; echo ']'
  samples = [
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-01_ATCACG_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-01_ATCACG_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-02_CGATGT_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-02_CGATGT_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-03_TTAGGC_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-03_TTAGGC_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-04_TGACCA_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-04_TGACCA_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-05_ACAGTG_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-05_ACAGTG_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-06_GCCAAT_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-06_GCCAAT_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-07_CAGATC_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-07_CAGATC_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-08_ACTTGA_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-08_ACTTGA_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-09_GATCAG_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-09_GATCAG_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-10_TAGCTT_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-10_TAGCTT_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-11_GGCTAC_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-11_GGCTAC_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-12_CTTGTA_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-12_CTTGTA_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-13_AGTCAA_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-13_AGTCAA_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-14_AGTTCC_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-14_AGTTCC_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-15_ATGTCA_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-15_ATGTCA_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-16_CCGTCC_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-16_CCGTCC_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-17_GTCCGC_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-17_GTCCGC_L003_R2.fastq'),
('/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-18_GTGAAA_L003_R1.fastq', '/data/tmp/marchulsman/data_feng/C1V35ACXX_102179-18_GTGAAA_L003_R2.fastq'),
]
  sample_names = [ "F1_a","F3_a", "F5_a", "F1_b", "F3_b", "F5_b", "F1_c","F3_c","F5_c","N1_a","N3_a", "N5_a","N1_b","N3_b","N5_b", "N1_c","N3_c","N5_c"];
  
  sample_labels = [ 0, 1, 2, 0, 1, 2, 0, 1, 2, 3, 4, 5, 3, 4, 5, 3, 4, 5];
  label_names   = [ "F1", "F3", "F5", "N1", "N3", "N5" ];

  genome = '/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/assembly/n402_atcc.unpadded.fasta' 
  genome_annot = '/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/annotations/n402_annotation.gff'

  #############################################################################
  # Trimmomatic options                                                       #
  #############################################################################
  #trimmomatic_opts="";
  trimmomatic_trim="ILLUMINACLIP:/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/aniger_rnaseq_dec2012/trim/illumina.fasta:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

  #############################################################################
  # STAR GENOME GENERATION OPTIONS                                            #
  #############################################################################
  #star_gg_opts = "";

  #############################################################################
  # STAR ALIGNMENT OPTIONS                                                    #
  #############################################################################
  #star_al_opts = "";
  
  cuffdiff_cmp = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (2, 5)]; 
  
#eclass
