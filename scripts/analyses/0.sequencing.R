# Sequencing

col_names <- c("file", "depth")
raw_reads <- read_csv("data/raw_reads.csv", col_names = col_names)
trimmed_reads <- list.files(path = "data/", 
                            pattern = "_trimmed_reads.csv", 
                            full.names = TRUE) %>% 
  map(\(i) read_csv(i, col_names = col_names)) %>% 
  list_rbind()

raw_reads <- mutate(raw_reads, 
                    sample = str_extract(file, "^[^_]+"),
                    gene = str_extract(file, "(16S|ITS)"),
                    pair = str_replace(file, ".*L1.([12]).fq.gz", "\\1")) %>% 
  dplyr::select(-c(pair, file)) %>% 
  distinct()

trimmed_reads <- mutate(trimmed_reads,
                        sample = str_extract(file, "^[^_]+"),
                        gene = str_extract(file, "(16S|ITS)"),
                        pair = str_replace(file, ".*\\.([12]).fq.gz", "\\1")) %>% 
  dplyr::select(-c(pair, file)) %>% 
  distinct()

all_reads <- left_join(raw_reads, trimmed_reads, 
                       by = c("sample", "gene"), 
                       suffix = c("_raw", "_trimmed")) %>% 
  pivot_longer(where(is.numeric), names_to = "process", values_to = "depth")

fig <- ggplot(all_reads, aes(x = sample, 
                      y = depth/1e5,
                      fill = process)) +
  geom_col(position = "dodge") +
  labs(y = "Number of reads (x 10,000)",
       x = "Sample",
       fill = "") +
  scale_fill_discrete(labels = c("Raw", "Trimmed")) +
  facet_grid(process ~ gene) +
  coord_flip()

walk(fig_formats, \(fmt) {
  filename <- str_glue("{results_dir}/figure.sequence_depth.{fmt}")
  ggsave(filename, 
         plot = fig, 
         width = 3, 
         height = 3.5,
         scale = 2,
         dpi = 600)  
})
