
library(tibble)

# https://github.com/oncokb/oncokb-annotator/blob/a80ef0ce937c287778c36d45bf1cc8397539910c/AnnotatorCore.py#L118
consequence_map = tribble(
  ~mutationType, ~consequence_final_coding,
  '3\'Flank',  'any',
  '5\'Flank',  'any',
  'Targeted_Region',  'inframe_deletion', 'inframe_insertion',
  'COMPLEX_INDEL',  'inframe_deletion', 'inframe_insertion',
  'ESSENTIAL_SPLICE_SITE',  'feature_truncation',
  'Exon skipping',  'inframe_deletion',
  'Frameshift deletion',  'frameshift_variant',
  'Frameshift insertion',  'frameshift_variant',
  'FRAMESHIFT_CODING',  'frameshift_variant',
  'Frame_Shift_Del',  'frameshift_variant',
  'Frame_Shift_Ins',  'frameshift_variant',
  'Fusion',  'fusion',
  'Indel',  'frameshift_variant', 'inframe_deletion', 'inframe_insertion',
  'In_Frame_Del',  'inframe_deletion',
  'In_Frame_Ins',  'inframe_insertion',
  'Missense',  'missense_variant',
  'Missense_Mutation',  'missense_variant',
  'Nonsense_Mutation',  'stop_gained',
  'Nonstop_Mutation',  'stop_lost',
  'Splice_Site',  'splice_region_variant',
  'Splice_Site_Del',  'splice_region_variant',
  'Splice_Site_SNP',  'splice_region_variant',
  'splicing',  'splice_region_variant',
  'Translation_Start_Site',  'start_lost',
  'vIII deletion',  'any',

  # Karissa Added,
  'Splice_Region',  'splice_region_variant'
)


usethis::use_data(consequence_map, overwrite = TRUE)
