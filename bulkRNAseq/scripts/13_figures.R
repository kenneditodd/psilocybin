library(ggplot2)
library(dplyr)
library(patchwork)
library(scales)

setwd(".")
files <- list.files("../refs/metascape_output/")

for (file in files) {
  # read table
  df <- read.delim2(paste0("../refs/metascape_output/", file))
  
  # reformat table
  colnames(df) <- tolower(colnames(df))
  df$gene_count <- sapply(strsplit(df$genes, "\\|"), length)
  df$gene_count <- as.numeric(df$gene_count)
  df$log10.p. <- as.numeric(as.character(df$log10.p.)) * -1
  
  # Set limits
  min_p <- min(df$log10.p., na.rm = TRUE)
  max_p <- max(df$log10.p., na.rm = TRUE)
  max_x <- ceiling(max(df$gene_count, na.rm = TRUE) / 5) * 5
  
  # Evenly spaced color breaks
  red_colors <- c("#FFE5E5", "#FF7F7F", "#B20000")    # light → medium → dark red
  blue_colors <- c("#E5F0FF", "#7FBFFF", "#004C99")   # light → medium → dark blue
  color_values <- rescale(c(min_p, (min_p + max_p)/2, max_p), to = c(0, 1))
  
  # Prepare UP plot
  df_up <- df %>%
    filter(direction == "UP") %>%
    arrange(gene_count) %>%
    mutate(pathway = factor(pathway, levels = unique(pathway)))
  p_up <- ggplot(df_up, aes(y = pathway, x = gene_count, fill = log10.p.)) +
    geom_col(width = 0.7) +
    scale_fill_gradientn(colors = red_colors, values = color_values, limits = c(min_p, max_p)) +
    scale_x_continuous(
      limits = c(0, max_x),
      breaks = seq(0, max_x, by = 5)
    ) +
    labs(y = NULL, x = NULL, fill = "Up-regulated Genes -log10(p)", title = gsub(".tsv", "", file)) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
  
  # Prepare DOWN plot
  df_down <- df %>%
    filter(direction == "DOWN") %>%
    arrange(gene_count) %>%
    mutate(pathway = factor(pathway, levels = unique(pathway)))
  p_down <- ggplot(df_down, aes(y = pathway, x = gene_count, fill = log10.p.)) +
    geom_col(width = 0.7) +
    scale_fill_gradientn(colors = blue_colors, values = color_values, limits = c(min_p, max_p)) +
    scale_x_continuous(
      limits = c(0, max_x),
      breaks = seq(0, max_x, by = 5)
    ) +
    labs(y = NULL, x = "Gene Count", fill = "Down-regulated Genes -log10(p)") +
    theme_minimal()

  
  # Stack with patchwork
  final_plot <- (p_up / p_down) +
    plot_layout(
      heights = c(nrow(df_up), nrow(df_down)),  # match # of bars
      guides = "collect"
    ) &
    theme(plot.margin = margin(2, 5, 0, 5),
          legend.position = "right")  # or "bottom"
  
  # Compute exact number of bars
  n_up <- nrow(df_up)
  n_down <- nrow(df_down)
  total_bars <- n_up + n_down
  
  # Fixed height per bar
  height_per_bar <- 0.4
  legend_buffer <- 1.5  # inches for legend height (works well for vertical legends)
  
  # Total height = bars + space for legend
  ggsave_height <- total_bars * height_per_bar + legend_buffer
  
  # Save
  path <- paste0("../results/figures/", gsub(".tsv", ".pdf", file))
  ggsave(path, final_plot, width = 15, height = ggsave_height, units = "in")
}
