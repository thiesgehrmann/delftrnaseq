#!/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/software/delftrnaseq/pipeline.py

class config(PIPELINECONF):
  #inst_loc

  email   = "m.hulsman@tudelft.nl";
  jobname = "RJTIME";

  location = "./conf";
  workdir = None;
  outdir  = "/data/tmp/marchulsman/run_rj"
  #makefile="";
  blast_db = '/data/tmp/marchulsman/blast_db/nt.nal'    
  


# echo '['; find | grep fastq | cut -c1-35 | sort | uniq | while read x; do x=`basename "$x"`; echo "('`pwd`/${x}1.fastq', '`pwd`/${x}2.fastq'),"; done; echo ']'
[
('/data/tmp/marchulsman/data_rj/102020-01_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-01_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-01_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-01_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-02_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-02_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-02_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-02_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-03_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-03_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-03_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-03_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-04_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-04_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-04_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-04_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-05_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-05_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-05_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-05_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-06_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-06_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-06_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-06_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-07_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-07_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-07_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-07_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-09_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-09_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-09_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-09_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-10_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-10_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-10_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-10_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-12_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-12_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-12_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-12_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-13_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-13_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-13_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-13_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-14_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-14_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-14_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-14_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-17_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-17_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-17_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-17_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-19_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-19_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-19_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-19_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-22_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-22_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-22_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-22_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-23_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-23_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-23_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-23_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-24_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-24_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-24_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-24_R2.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-28_R1.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-28_R1.fastq2.fastq'),
('/data/tmp/marchulsman/data_rj/102020-28_R2.fastq1.fastq', '/data/tmp/marchulsman/data_rj/102020-28_R2.fastq2.fastq'),
]

  sample_names =  ["n402_d3_r3","dhex_d7_r3","dhex_d3_r1","dhex_d7_r1","n402_d5_r1","n402_d5_r3","n402_d7_r3","n402_d9_r3","n402_d3_r2","n402_d7_r1","dhex_d3_r2","dhex_d7_r2","n402_d5_r2","n402_d3_r1","n402_d7_r2","n402_d9_r2","n402_d9_r1","dhex_d3_r3"]
  sample_labels = [0,            5,           4,           5,          1,            1,           2,           3,           0,           2,          4,            5,           1,           0,           2,           3,           3,           4]
  label_names   = [ "n402_d3", "n402_d5", "n402_d7", "n402_d9", "dhex_d3", "dhex_d7" ];

  genome = '/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/assembly/n402_atcc.unpadded.fasta' 
  genome_annot = '/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/annotations/results/n402_annotations.gff'

  annots_file = '/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/n402_sequence/annotation_data/annots_v2.dat'

  #############################################################################
  # Trimmomatic options                                                       #
  #############################################################################
  #trimmomatic_opts="";
  trimmomatic_trim="ILLUMINACLIP:/tudelft.net/staff-groups/ewi/insy/DBL/marchulsman/projects/aniger_rnaseq_dec2012/trim/illumina.fasta:2:40:15 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

  #############################################################################
  # STAR GENOME GENERATION OPTIONS                                            #
  #############################################################################
  #adapt sjdbOverhang to read length - 1
  def star_gg_opts(self):
    return "--runThreadN %d --sjdbOverhang 99" % self.__max_threads__;

  #############################################################################
  # STAR ALIGNMENT OPTIONS                                                    #
  #############################################################################
  def star_preal_opts(self):
    return "--runThreadN %d  --alignIntronMax 5000 --alignMatesGapMax 5000" % self.__max_threads__;
  
  #splice_sites_minmax_overhang = 36

  def star_al_opts(self):
    return "--runThreadN %d  --alignIntronMax 5000 --alignMatesGapMax 5000" % self.__max_threads__;
  

  cuffdiff_cmp = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 4), (2, 5), (3, 4), (3, 5), (4, 5)]; 
  
#eclass
