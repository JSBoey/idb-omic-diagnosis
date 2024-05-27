# Alpha diversity

# TODO: Chao diversity

alpha_diversity <- map(RCOUNTS, \(mt) {
  tibble(
    "sample" = colnames(mt),
    "Richness" = specnumber(mt, MARGIN = 2),
    "Simpson" = diversity(mt, index = "invsimpson", MARGIN = 2),
    "Evenness" = Simpson/Richness
  ) %>% 
    left_join(METADATA)
})

alpha_diversity_tests <- map(alpha_diversity, \(tb) {
  index <- c("Richness", "Simpson", "Evenness") %>% 
    setNames(., .)
  variable <- c("status", "individual", "collection_time")
  
  as_tibble(
    expand.grid("Metric" = index, 
                "Variable" = variable)
  ) %>% 
    mutate(
      test = map2(Metric, Variable, \(i, j) {
        g <- as.factor(tb[[j]])
        x <- tb[[i]]
        
        kruskal.test(x ~ g)
      }),
      statistic = map_dbl(test, "statistic"),
      p = map_dbl(test, "p.value")
    )
})

alpha_diversity_plot <- map2(alpha_diversity, names(alpha_diversity), 
                             \(tb, nm) {
  tb <- pivot_longer(tb, c("Richness", "Simpson", "Evenness"), 
                       names_to = "index",
                       values_to = "value") %>% 
    mutate(index = factor(index, levels = c("Richness", "Simpson", "Evenness")))
  ggplot(tb, aes(x = status, y = value)) +
    geom_boxplot(alpha = 0.25, 
                 outliers = FALSE) +
    geom_point(aes(colour = individual, 
                   shape = collection_time,
                   group = individual),
               size = 2,
               position = position_dodge(width = 0.75, 
                                         preserve = "total")) +
    facet_wrap(~ index, scales = "free_y") +
    labs(title = nm) +
    scale_x_discrete(name = "status",
                     labels = c("control" = "Control",
                                "UC" = "UC")) +
    scale_colour_brewer(palette = "Paired",
                        name = "Individual") +
    scale_shape_discrete(name = "Collection time") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank())
})

fig <- wrap_plots(alpha_diversity_plot, nrow = 2, ncol = 1, guides = "collect")

walk(fig_formats, \(fmt) {
  filename <- str_glue("{results_dir}/figure.alpha_diversity.{fmt}")
  ggsave(filename, 
         plot = fig, 
         width = 3, 
         height = 2,
         scale = 2.75,
         dpi = 600)  
})





