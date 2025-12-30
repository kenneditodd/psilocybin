# server.R
library(shiny)
library(dplyr)
library(ggplot2)
library(ggtranscript)

cb18 <- c(
  "#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442", "#999999",
  "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#661100", "#AA4499",
  "#882255", "#000000"
)

quant_long <- readRDS("quant_long.rds")
ann <- readRDS("annotation.rds")

server <- function(input, output, session) {
  
  quant_long <- readRDS(normalizePath("quant_long.rds", mustWork = TRUE))
  ann <- readRDS(normalizePath("annotation.rds", mustWork = TRUE))
  
  get_fill_palette <- function(n) cb18[seq_len(min(n, length(cb18)))]
  
  observe({
    gene_choices <- quant_long %>%
      dplyr::pull(gene_name) %>%
      unique() %>%
      sort()
    
    updateSelectizeInput(
      session,
      inputId = "gene_name",
      choices = gene_choices,
      selected = if ("Dbp" %in% gene_choices) "Dbp" else gene_choices[1],
      server = TRUE
    )
  })
  
  plot_df <- reactive({
    req(input$gene_name)
    
    quant_long %>%
      filter(gene_name == input$gene_name) %>%
      mutate(
        transcript_label = if_else(
          is.na(transcript_name) | transcript_name == "",
          target_id,
          transcript_name
        )
      )
  })
  
  plot_ensembl_protein_coding <- function(goi.gtf) {
    goi <- goi.gtf$gene_name[1]
    
    goi.gtf <- subset(goi.gtf, transcript_type == "protein_coding")
    goi.exons <- subset(goi.gtf, type == "exon")
    goi.cds <- subset(goi.gtf, type == "CDS")
    
    goi.exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_name
      )) +
      geom_range(
        fill = "gray",
        height = 0.25
      ) +
      geom_range(
        data = goi.cds,
        fill = "purple"
      ) +
      geom_intron(
        data = to_intron(goi.exons, "transcript_name"),
        aes(strand = strand)
      ) +
      ggtitle(label = paste(goi, "Ensembl Annotation")) +
      scale_x_continuous(name = paste("Chromosome", goi.gtf$seqnames[1], "Position")) +
      scale_y_discrete(name = "Transcript Name") +
      theme_classic() +
      theme(text = element_text(size = 20))
  }
  
  plot_ensembl_all <- function(goi.gtf) {
    goi <- goi.gtf$gene_name[1]
    goi.exons <- subset(goi.gtf, type == "exon")
    
    goi.exons %>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = transcript_name
      )) +
      geom_range(aes(fill = transcript_type)) +
      geom_intron(
        data = to_intron(goi.exons, "transcript_name"),
        aes(strand = strand)
      ) +
      ggtitle(label = paste(goi, "Ensembl Annotation")) +
      scale_x_continuous(name = paste("Chromosome", goi.gtf$seqnames[1], "Position")) +
      scale_y_discrete(name = "Transcript Name") +
      theme_classic() +
      theme(text = element_text(size = 20))
  }
  
  plot_ensembl_table <- function(goi.gtf) {
    goi.gtf <- subset(goi.gtf, type == "transcript")
    goi.gtf
  }
  
  output$isoform_tpm_bar <- renderPlot({
    df <- plot_df()
    req(nrow(df) > 0)
    
    tpm_col <- switch(
      input$tpm_level,
      "group"  = "tpm_group",
      "group2" = "tpm_group2"
    )
    
    grp_col <- switch(
      input$tpm_level,
      "group"  = "group",
      "group2" = "group2"
    )
    
    y_lab <- switch(
      input$tpm_level,
      "group"  = "TPM (group avg)",
      "group2" = "TPM (group2 avg)"
    )
    
    gene_id_val <- df$gene_id[!is.na(df$gene_id)][1]
    plot_title <- paste0(input$gene_name, ":", gene_id_val)
    
    group_levels <- c(
      "S.8h", "S.24h", "S.7d",
      "L.8h", "L.24h", "L.7d",
      "H.8h", "H.24h", "H.7d"
    )
    
    group2_levels <- as.vector(
      rbind(
        paste0(group_levels, ".F"),
        paste0(group_levels, ".M")
      )
    )
    
    if (grp_col == "group") {
      df$group <- factor(df$group, levels = group_levels)
    } else {
      df$group2 <- factor(df$group2, levels = group2_levels)
    }
    
    ggplot(
      df,
      aes(
        x = transcript_label,
        y = .data[[tpm_col]],
        fill = .data[[grp_col]]
      )
    ) +
      geom_col(position = position_dodge(width = 0.85), width = 0.8) +
      labs(
        title = plot_title,
        x = "Isoform (transcript_name)",
        y = y_lab,
        fill = input$tpm_level
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title  = element_text(face = "bold"),
        panel.grid.major.x = element_blank()
      ) +
      scale_fill_manual(values = cb18)
  }, res = 120)
  
  output$ensembl_transcripts <- renderPlot({
    req(input$gene_name)
    
    goi <- subset(ann, gene_name == input$gene_name)
    req(nrow(goi) > 0)
    
    if (input$tx_mode == "protein coding only") {
      plot_ensembl_protein_coding(goi)
    } else {
      plot_ensembl_all(goi)
    }
  }, res = 120)
  
  output$ensembl_table <- renderDataTable({
    req(input$gene_name)
    
    goi <- subset(ann, gene_name == input$gene_name)
    req(nrow(goi) > 0)
    
    plot_ensembl_table(goi)
  })
  
}
