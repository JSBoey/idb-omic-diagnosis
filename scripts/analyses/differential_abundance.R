# Differential abundance

# Comparison between status

comp_A <- map(COUNTS, \(m) {
  covariates <- dplyr::filter(METADATA, sample %in% colnames(m))
  g <- as.character(covariates$status)
  
  x <- aldex.clr(m, g, mc.samples = 1e3, verbose = TRUE, gamma = 0.25)
  tt <- aldex.ttest(x, hist.plot = FALSE, paired.test = FALSE, verbose = TRUE)
  eff <- aldex.effect(x, CI = TRUE, paired.test = FALSE, verbose = TRUE)
  
  data.frame(tt, eff)
  
})

par(mfrow = c(1, 2))
walk2(comp_A, names(comp_A), \(x, nm) {
  aldex.plot(x, type = "volcano", test = "welch", 
             main = str_glue("Volcano plot for {nm}"))
})

walk2(comp_A, names(comp_A), \(df, gene) {
  write_csv(rownames_to_column(df, "ASV"), 
            str_glue("{results_dir}/aldex2_status.ASV.{gene}.csv"))
})

# At ASV level, there are no differentially abundant ASVs between health status, but there are a few taxa that are differentially abundant between individuals.

# Comparison between status at genus level

comp_B <-  map2(COUNTS, TAXONOMY[c("16S.Greengenes2", "ITS.UNITE")], \(m, tx) {
  covariates <- dplyr::filter(METADATA, sample %in% colnames(m))
  g <- as.character(covariates$status)
  
  m <- left_join(
    as_tibble(m, rownames = "ASV"),
    dplyr::select(tx, -Confidence), 
    by = join_by("ASV" == "Feature ID")
  ) %>% 
    group_by(Genus) %>% 
    summarise(across(where(is.numeric), sum)) %>% 
    dplyr::filter(!(is.na(Genus) | grepl("^g__$", Genus))) %>% 
    column_to_rownames("Genus") %>% 
    as.matrix()
  
  x <- aldex.clr(m, g, mc.samples = 1e3, verbose = TRUE, gamma = 0.25)
  tt <- aldex.ttest(x, hist.plot = FALSE, paired.test = FALSE, verbose = TRUE)
  eff <- aldex.effect(x, CI = TRUE, paired.test = FALSE, verbose = TRUE)

  data.frame(tt, eff)
  
})

par(mfrow = c(1, 2))
walk2(comp_B, names(comp_B), \(x, nm) {
  aldex.plot(x, type = "volcano", test = "welch", 
             main = str_glue("Volcano plot for {nm}"))
})

walk2(comp_B, names(comp_A), \(df, gene) {
  write_csv(rownames_to_column(df, "ASV"), 
            str_glue("{results_dir}/aldex2_status.genus.{gene}.csv"))
})




