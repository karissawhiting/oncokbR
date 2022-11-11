library(tibble)

mutations_hgvsg <- tibble::tribble(
  ~NCBI_Build, ~Hugo_Symbol, ~Variant_Classification,          ~Tumor_Sample_Barcode,     ~HGVSp_Short,         ~HGVSp,             ~HGVSg, ~Chromosome, ~Start_Position, ~End_Position, ~Reference_Allele, ~Tumor_Seq_Allele1, ~Tumor_Seq_Allele2,
     "GRCh38",       "CUL1",     "Missense_Mutation", "TCGA-A6-2672-01A-01W-0833-10",        "p.Y466S",    "Tyr466Ser",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",       "AKT3",     "Nonsense_Mutation",              "TCGA-05-4417-01",        "p.E182*",      "Glu182*",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",     "PIK3CA",     "Missense_Mutation",              "TCGA-02-0033-01",        "p.E542K",    "Glu542Lys", "3:g.179218294G>A",          3L,      179218294L,    179218294L,               "G",                "A",                "A",
     "GRCh38",      "FGFR3",     "Missense_Mutation",              "TCGA-05-4417-01",        "p.V271M",    "Val271Met",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",       "EGFR",     "Missense_Mutation",              "TCGA-06-0155-01",        "p.H304Y",    "His304Tyr",  "7:g.55155850C>T",          7L,       55155850L,     55155850L,               "C",                "T",                "T",
     "GRCh38",       "PTEN",     "Missense_Mutation",              "TCGA-06-0155-01",        "p.C136R",    "Cys136Arg", "10:g.87933165T>C",         10L,       87933165L,     87933165L,               "T",                "C",                "C",
     "GRCh38",      "FGFR2",     "Missense_Mutation",              "TCGA-02-0033-01",        "p.Q212K",    "Gln121Lys",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",        "ATM",     "Missense_Mutation",              "TCGA-05-4417-01",       "p.L2890R",   "Leu2890Arg",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",       "KRAS",     "Missense_Mutation",              "TCGA-05-4417-01",         "p.G12C",     "Gly12Cys", "12:g.25245351C>A",         12L,       25245351L,     25245351L,               "C",                "A",                "A",
     "GRCh38",        "RB1",     "Nonsense_Mutation",              "TCGA-02-0033-01",        "p.Q702*",      "Gln702*",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",       "TP53",     "Missense_Mutation",              "TCGA-02-0033-01",        "p.R248Q",    "Arg248Gln",  "17:g.7674220C>T",         17L,        7674220L,      7674220L,               "C",                "T",                "T",
     "GRCh38",        "NF1",           "Splice_Site",              "TCGA-02-0033-01", "p.X1445_splice", "X1445_splice", "17:g.31259031G>A",         17L,       31259031L,     31259031L,               "G",                "A",                "A",
     "GRCh38",      "STK11",     "Missense_Mutation",              "TCGA-05-4417-01",        "p.H168R",    "His168Arg",                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA,
     "GRCh38",      "MYD88",     "Missense_Mutation",              "TCGA-05-4417-01",          "M219T",             NA,                 NA,          NA,              NA,            NA,                NA,                 NA,                 NA
  )



usethis::use_data(mutations_hgvsg, overwrite = TRUE)
