# Taxonomy

# 1. Are 16S taxonomies between databases congruent?

match_taxa <- function(x, y) {
  taxa <- c(x, y)
  
  if (anyNA(taxa)) return(-1)
  
  if (x == y) return(1)
  
  p <- taxa[which.min(sapply(taxa, str_length))]
  s <- taxa[taxa != p]
  as.numeric(str_detect(s, p))
  
}

taxonomy_match <- list(
  "Greengenes2" = as.list(TAXONOMY[["16S.Greengenes2"]][, -c(1, 9)]),
  "SILVA" = as.list(TAXONOMY[["16S.SILVA"]][, -c(1, 9)])
) %$% 
  map2(Greengenes2, SILVA, Vectorize(match_taxa)) %>% 
  as.data.frame(row.names = TAXONOMY[["16S.SILVA"]][[1]])

par(mfrow = c(2, 2))
barplot(colSums(taxonomy_match) / nrow(taxonomy_match), 
        main = "Agreement between Greengenes2 and SILVA 138.1")
plot(density(TAXONOMY[["16S.Greengenes2"]]$Confidence), 
     col = "forestgreen",
     xlab = "Confidence score",
     main = "Prediction confidence")
lines(density(TAXONOMY[["16S.SILVA"]]$Confidence), 
      col = "red")
boxplot(
  data.frame("Greengenes2" = TAXONOMY[["16S.Greengenes2"]]$Confidence,
             "SILVA" = TAXONOMY[["16S.SILVA"]]$Confidence),
  ylab = "Confidence score"
)
wilcox.test(TAXONOMY[["16S.Greengenes2"]]$Confidence, 
            TAXONOMY[["16S.SILVA"]]$Confidence,
            alternative = "two.sided")

## Conclusion: Going with Greengenes2 for better prediction confidence and it is a newer database with links to genomic information on GTDB

# 2. Relative abundance of taxa

taxa_summarised_counts <- map2(
  COUNTS,
  TAXONOMY[c("16S.Greengenes2", "ITS.UNITE")],
  \(x, y) {
    # Join tables
    x <- as_tibble(x, rownames = "ASV")
    y <- dplyr::select(y, -Confidence)
    xy <- left_join(x, y, by = join_by("ASV" == "Feature ID"))
    
    # Summarise
    lapply(taxa_levels, \(s) {
      S <- str_to_lower(str_sub(s, 1, 1))
      group_by(xy, .data[[s]]) %>% 
        mutate(
          {{ s }} := ifelse(
            is.na(.data[[s]]) | .data[[s]] == str_glue("{S}__"), 
            str_glue("Unclassified_{s}"), 
            .data[[s]]
          )
        ) %>%
        summarise(across(where(is.numeric), sum))
    })
  }
)

taxa_summarised_abundance <- map_depth(
  taxa_summarised_counts, 2, \(tb) {
    mutate(tb, 
           across(where(is.numeric), proportions))
  }
)

relative_abundance_plot <- map_depth(
  taxa_summarised_counts, 2, \(tb) {
    s <- names(tb)[1]
    sy <- sym(s)
    Ar <- mutate(tb, 
                 across(where(is.numeric), proportions))
    Ar_long <- pivot_longer(Ar, 
                            where(is.numeric), 
                            names_to = "sample", 
                            values_to = "relative_abundance") %>% 
      mutate({{ s }} := str_remove(.data[[s]], "[a-z]__")) %>% 
      left_join(METADATA)
    
    Ar_mean <- group_by(Ar_long, status, .data[[s]]) %>% 
      summarise(avg_abd = mean(relative_abundance)) %>% 
      arrange(avg_abd)
    
    top_taxa <- if (nrow(Ar_mean) > 10) {
      slice_max(Ar_mean, order_by = avg_abd, n = 10) %>% 
        arrange(avg_abd)
    } else {
      Ar_mean
    }
    
    fig_tb <- Ar_long[Ar_long[[s]] %in% unique(top_taxa[[s]]), ]
    fig_tb[[s]] <- factor(fig_tb[[s]], levels = unique(top_taxa[[s]]))
    fig_tb$sample <- factor(fig_tb$sample, levels = METADATA$sample)
    
    n_taxa <- length(unique(top_taxa[[s]]))
    main_colour <- RColorBrewer::brewer.pal(12, "Paired")
    secondary_colour <- RColorBrewer::brewer.pal(12, "Set3")
    taxa_fill <- if (n_taxa > 12) {
      c(main_colour, sample(secondary_colour, n_taxa - 12))
    } else {
      main_colour
    }

    fig <- ggplot(fig_tb, aes(x = sample,
                              y = relative_abundance,
                              fill = !!sy)) +
      geom_col() +
      scale_y_continuous(expand = c(0, 0), 
                         limits = c(0, 1)) +
      scale_fill_manual(values = taxa_fill, 
                        labels = \(s) str_replace_all(s, "_", " ")) +
      labs(x = "Sample",
           y = "Relative abundance") +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

    list(
      "average_abundance" = Ar_mean,
      "top_taxa" = top_taxa,
      "figure" = fig
    )
  }
)

plot_list <- list(
  "16S_Phylum" = relative_abundance_plot[["16S"]]$Phylum$figure,
  "16S_Genus" = relative_abundance_plot[["16S"]]$Genus$figure,
  "ITS_Phylum" = relative_abundance_plot[["ITS"]]$Phylum$figure,
  "ITS_Genus" = relative_abundance_plot[["ITS"]]$Genus$figure
) %>% 
  map2(., names(.), \(gg, nm) gg + ggtitle(str_replace(nm, "_", " ")))

fig <- wrap_plots(plot_list, ncol = 2, nrow = 2, byrow = TRUE)

walk(fig_formats, \(fmt) {
  filename <- str_glue("{results_dir}/figure.taxonomy.{fmt}")
  ggsave(filename, 
         plot = fig, 
         width = 6, 
         height = 5.5,
         scale = 2,
         dpi = 600)  
})

relative_abundance_plot[["16S"]]$Phylum$figure
relative_abundance_plot[["16S"]]$Genus$average_abundance
relative_abundance_plot[["16S"]]$Genus$top_taxa

map(relative_abundance_plot[["16S"]], ~ length(unique(.x$top_taxa[[2]])))



