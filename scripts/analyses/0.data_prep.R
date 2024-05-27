# Data prep

# This script is used for data cleaning and wrangling, resulting in an .RData file for loading when the project is used. Do not run more than once unless there are changes to the underlying data structures.

genes <- c("16S", "ITS") %>% 
  setNames(., .)

COUNTS <- lapply(genes, \(s) {
  counts(read.biom(str_glue("data/{s}.feature_table.biom")))
})

tmp.taxonomy <- list.files("data", "taxonomy", full.names = T)
taxa_levels <- c("Domain", 
                 "Phylum", 
                 "Class", 
                 "Order", 
                 "Family", 
                 "Genus", 
                 "Species") %>% 
  setNames(., .)
TAXONOMY <- lapply(tmp.taxonomy, \(s) {
  tb <- read_tsv(s)
  separate(tb, "Taxon", taxa_levels, ";")
})
names(TAXONOMY) <- str_extract(tmp.taxonomy, "[0-9A-Z]+.[^_]+")

PHYLOGENY <- lapply(genes, \(s) read.tree(str_glue("data/{s}.rooted_tree.nwk")))

METADATA <- read_csv("data/sample_data.csv", col_types = "cfff")

# Reorder columns according to metadata
COUNTS <- lapply(COUNTS, \(m) {
  cn <- METADATA$sample[METADATA$sample %in% colnames(m)]
  m[, cn]
})

RCOUNTS <- lapply(COUNTS, \(m) {
  min_count <- min(colSums(m))
  t(
    vegan::rrarefy(t(m), min_count)
  )
})

fig_formats <- c("png", "tiff", "svg")
results_dir <- "./results"

rm(tmp.taxonomy)
save.image()