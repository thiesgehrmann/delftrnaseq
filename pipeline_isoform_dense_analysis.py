#!/usr/bin/env python
from pipeline_common import *;
from utilities import splicing_statistics as sstats;
from utilities import enrichment as enrich;
from utilities import multi_test_corr as mtc;
from utilities import latex;
from ibidas import *;

C = init_conf()

prg="STAR"

cmds = [];

#old_gff, old_fasta = cor(C.__isoform_dense_genome_split_output__);
orig_gff  = cor(C.__genome_annot_format_output__);
new_gff   = cor(C.__isoform_dense_genome_unsplit_output__);
attr_name = cor(C.__isoform_dense_attr_name__);
outdir    = cor(C.__isoform_dense_genome_analysis_outdir__);



cmd = "mkdir -p %s" % outdir;

retval = run_cmd(cmd);
if retval != 0:
  sys.exit(retval);
#fi


GO = Read(orig_gff);
NG = Read(new_gff);

IDX = sstats.index_gff3_id(GO);
  
  # For each exon in each isoform, classify what kind of variant it is
D = sstats.diagnose_annotation(IDX, NG, attr_name, display=False);
  # Gather statistics
S = sstats.gather_diagnosis_statistics(D);
  # Print the statistics
T, F = sstats.draw_statistics(S, outdir, outfmt='svg');

Export(T, '%s/output_counts.tsv' % outdir);
Export(T, '%s/output_counts.dat' % outdir);
Export(D, '%s/diagnosis.tsv' % outdir);
Save(D,   '%s/diagnosis.dat' % outdir);
#Export(sstats.gff_add_diagnosis(D, NG), '%s/diagnosed_annotation.gff' % outdir);

cmd = "gffread -EF %s/diagnosed_annotation.gff -o %s/diagnoised_annotation.cleaned.gff" % (outdir, outdir);

#retval = run_cmd(cmd);
#if retval != 0:
#  sys.exit(retval);
#fi


diagnosis = D

# Draw figures for interesting genes
if C.isoform_dense_genome_analysis_genes_of_interest != None:
  genes_of_interest = C.isoform_dense_genome_analysis_genes_of_interest;

  for gene, geneid in genes_of_interest:
    out_name = '%s/isoforms_of_%s_%s.pdf' % (outdir, gene, str(geneid));
    sstats.draw_gene_isoforms(diagnosis, geneid, out_name, outfmt='pdf');
  #efor
#fi

# Enrichment

annotation_files, annotation_names = ( cor(C.isoform_dense_genome_annotation_files), cor(C.isoform_dense_genome_annotation_names) );

enrich_table_gene    = sstats.gather_per_gene_statistics(diagnosis,    func=lambda d: any(d));
enrich_table_isoform = sstats.gather_per_isoform_statistics(diagnosis, func=lambda d: any(d));

counts = T

def organize_isoform_anal(D, func):

  M = D.GroupBy(_.alt_gene).To(_.orig_gene, Do=_.Get(_.orig_gene[0])).Without(_.exon_start, _.exon_end, _.orig_exon_coverage);
  M = M.To(_.retention, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt5, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt3,  Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.skipped, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.new_exon, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.identical, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt_gene, Do= _ / 'id');

  return M.Copy();
#edef

def organize_gene_anal(D, func):

  M = D.GroupBy(_.orig_gene).Without(_.alt_gene, _.exon_start, _.exon_end, _.orig_exon_coverage);
  M = M.To(_.retention, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt5, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.alt3,  Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.skipped, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.new_exon, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.identical, Do=_.Array().Each(lambda x: func(x)));
  M = M.To(_.orig_gene, Do= _ / 'id');

  return M.Copy()
#edef
for filter, filter_name in zip(cor(C.isoform_dense_genome_analysis_filter) , cor(C.isoform_dense_genome_analysis_filter_names) ):

  print filter, filter_name

  if filter == None:
    table_filt = D;
  else:
    filter = Read(filter);
    table_filt = D[_.orig_gene.In(filter.Get(0))];
  #fi

  enrich_table_gene_filter    = organize_gene_anal(table_filt,    func=lambda d: any(d))
  enrich_table_isoform_filter = organize_isoform_anal(table_filt, func=lambda d: any(d));

  l = latex.LatexFile('%s/isoform_analysis_filter=%s.tex' % ( outdir, filter_name))
  l.write_title("Isoform analysis for dense genomes for ``%s''" % C.title, C.author);
  l.add_text("Filter for this analysis: %s.\n" % str(filter_name));

  S = sstats.gather_diagnosis_statistics(table_filt);
    # Print the statistics
  T, F = sstats.draw_statistics(S, outdir, suffix='filter=%s' % filter_name, outfmt='pdf');

  l.start_section('Basic Statistics');
  l.write_rep(T, 'Counts of various splicing events');
  for i, (filename, title) in enumerate(F):
    l.include_figure(filename, 'statistics%s' % str(i),  title);
  #efor

  if C.isoform_dense_genome_analysis_genes_of_interest != None:
    l.start_section('Isoforms found for specific genes');
    for gene, geneid in C.isoform_dense_genome_analysis_genes_of_interest:
      l.include_figure('isoforms_of_%s_%s.pdf' % ( gene, str(geneid)), 'isoform_%s' % str(gene), "Isoforms discovered for %s" % gene)
    #efor
  #fi

  for annot, annot_name in zip(cor(C.isoform_dense_genome_annotation_files), cor(C.isoform_dense_genome_annotation_names) ):
 
    l.start_section('Enrichment of GO terms in sets of genes');
 
    annot_table = Read(annot);
  
    for table, table_name in zip([enrich_table_gene_filter, enrich_table_isoform_filter], ['genes', 'isoforms']):

      for field in [ 'retention', 'alt5', 'alt3', 'skipped', 'new_exon' ]:
        enrich_table = table.Get(_.id, field) | Match(0,0, jointype='left', merge_same='equi')| annot_table.GroupBy(0).Get(0,1);

        print field, annot, table_name
        enrichment = enrich.fast_enrich(enrich_table, annot_table.Without(0).Unique(0), 0.05);
        l.write_rep(enrichment[_.qvalue < cor(C.__isoform_dense_genome_analysis_enrichment_alpha__)].Sort(_.qvalue), '%s annotations significantly enriched in %s which have isoforms that have a %s' % (annot_name, table_name, field));
      #efor
    #efor
  #efor

  l.end_document()

#efor

sys.exit(run_seq_cmds(cmds));

