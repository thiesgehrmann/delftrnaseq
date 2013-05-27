#!./pipeline.py

class config(PIPELINECONF):
  #inst_loc

  email   = "thiesgehrmann@gmail.com";
  jobname = "test";

  workdir = None;
  outdir  = "rnapipeline_test_data/outdir"
  #makefile="";


  # echo '['; find | grep fastq | cut -c1-35 | uniq | while read x; do x=`basename "$x"`; echo "('`pwd`/${x}1.fastq', '`pwd`/${x}2.fastq'),"; done; echo ']'
  samples=[("my_data/102020-01_R1.short.fastq", "my_data/102020-01_R2.short.fastq")];
  sample_names = [ "sample1"];

  genome = [ '/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/assembly/n402_atcc.unpadded.fasta' ] 
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

  
#eclass
