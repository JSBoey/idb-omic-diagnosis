# Beta diversity

dist_matrix <- list(
  "UniFrac" = map2(RCOUNTS, PHYLOGENY, \(m, tr) {
    rbiom::unifrac(m, weighted = TRUE, tree = tr)
  }),
  "Bray_Curtis" = map(RCOUNTS, \(m) {
    vegdist(t(m))
  }),
  "Robust_Aitchison" = map(COUNTS, \(m) {
    vegdist(t(m), method = "robust.aitchison")
  })
)

ordination <- map2(dist_matrix, names(dist_matrix), \(d, nm) {
  if (nm == "Robust_Aitchison") {
    rclr <- map(COUNTS, \(m) decostand(t(m), method = "rclr"))
    o <- map(rclr, \(m) rda(m ~ 1))
    return(o)
  }
  
  map(d, \(D) metaMDS(D, try = 1e3-1, trymax = 1e4-1, autotransform = FALSE))
}) %>% 
  list_flatten()

permanova <- map_depth(dist_matrix, 2, \(d) {
  data <- dplyr::filter(METADATA, sample %in% labels(d))
  
  if (!all(data$sample == labels(d))) 
    stop("Order samples properly!")
  
  variable <- c("status", "individual", "collection_time") %>% 
    setNames(., .)
  
  map(variable, \(v) {
    f <- paste0("d ~ ", v)
    broom::tidy(
      adonis2(as.formula(f), data, permutations = 1e3-1)
    )
  }) %>% 
    list_rbind(names_to = "variable")
}) %>% 
  map(\(l) list_rbind(l, names_to = "gene")) %>% 
  list_rbind(names_to = "distance")

write_csv(permanova, str_glue("{results_dir}/table.permanova.csv"))

dispersion <- map_depth(dist_matrix, 2, \(d) {
  data <- dplyr::filter(METADATA, sample %in% labels(d))
  
  if (!all(data$sample == labels(d))) 
    stop("Order samples properly!")
  
  variable <- c("status", "individual", "collection_time") %>% 
    setNames(., .)
  
  map(variable, \(v) {
    betadisper(d, data[[v]], bias.adjust = TRUE, sqrt.dist = TRUE)
  })
}) %>% 
  map(list_flatten) %>% 
  list_flatten()

dispersion_test <- map(dispersion, \(disp) {
  as_tibble(permutest(disp)$tab, rownames = "measure") 
}) %>% 
  list_rbind(names_to = "variable") %>% 
  mutate(
    "distance_metric" = str_remove(variable, "_(ITS|16S)_.*$"),
    "gene" = str_extract(variable, "(ITS|16S)"),
    "variable" = str_remove(variable, "^.*_(ITS|16S)_")
  )

write_csv(dispersion_test, str_glue("{results_dir}/table.dispersion.csv"))

# Plot

plot_data_ordination <- map(ordination, \(ord) {
  
  xy <- scores(ord, display = "sites", tidy = TRUE) %>% 
    rownames_to_column("sample") %>% 
    left_join(METADATA) %>% 
    dplyr::select(-c("score", "label")) %>% 
    na.omit()
  
  names(xy)[c(2, 3)] <- c("Axis1", "Axis2")
  xy
  
}) %>% 
  list_rbind(names_to = "variable") %>% 
  mutate(
    distance_metric = str_remove(variable, "_(ITS|16S)$"),
    distance_metric = factor(distance_metric, levels = c("Bray_Curtis",
                                                         "UniFrac",
                                                         "Robust_Aitchison")),
    gene = str_extract(variable, "(ITS|16S)")
  )

plot_ordination <- ggplot(plot_data_ordination, aes(x = Axis1,
                                                    y = Axis2,
                                                    colour = status)) +
  geom_point() +
  labs(x = "Axis 1", 
       y = "Axis 2",
       colour = "") +
  scale_colour_discrete(labels = c("control" = "Control",
                                   "UC" = "UC")) +
  facet_wrap(distance_metric ~ gene, 
             ncol = 2, 
             nrow = 3,
             scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

walk(fig_formats, \(fmt) {
  filename <- str_glue("{results_dir}/figure.ordination.{fmt}")
  ggsave(filename, 
         plot = plot_ordination, 
         width = 2.2, 
         height = 3,
         scale = 2.5,
         dpi = 600)  
})


plot_data_dispersion <- map(dispersion, \(disp, nm) {
  
  xy <- disp %$% tibble(
    "sample" = names(distances),
    "distance" = distances,
    "group" = group
  )
  
}) %>% 
  list_rbind(names_to = "variable") %>% 
  mutate(
    distance_metric = str_remove(variable, "_(ITS|16S)_.*$"),
    distance_metric = factor(distance_metric, levels = c("Bray_Curtis",
                                                         "UniFrac",
                                                         "Robust_Aitchison")),
    gene = str_extract(variable, "(ITS|16S)"),
    variable = str_remove(variable, "^.*_(ITS|16S)_"),
    variable = factor(variable, levels = c("status",
                                           "collection_time",
                                           "individual"))
  )

plot_dispersion <- ggplot(plot_data_dispersion, aes(x = group, 
                                                    y = distance)) +
  geom_boxplot() +
  labs(x = "Grouping",
       y = "Ordination distance") +
  facet_grid(distance_metric ~ gene + variable, 
             space = "free_x",
             scales = "free") +
  theme_bw() +
  theme(panel.grid = element_blank())

walk(fig_formats, \(fmt) {
  filename <- str_glue("{results_dir}/figure.dispersion.{fmt}")
  ggsave(filename, 
         plot = plot_dispersion, 
         width = 4, 
         height = 2.5,
         scale = 2.5,
         dpi = 600)  
})
