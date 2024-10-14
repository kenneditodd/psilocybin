# Load libraries
library(dplyr)
library(ggplot2)
library(ggrepel)
library(plotly)
library(shiny)
library(stringr)

# Load input options
comparisons <- readr::read_tsv("input_options.tsv")
group1_options <- comparisons$group1
geneOptions <- readRDS("gene_options.rds")  # Assume this is a character vector of gene names
cpm <- read.table("counts/CPM_postfiltering.tsv")
theme_set(theme_bw())

# Function to assign colors based on FDRq and logFC thresholds
assign_colors <- function(data, fdrq, lfc) {
  data %>%
    mutate(color_adjpval = case_when(
      adj.P.Val < fdrq & logFC > lfc ~ "up-regulated",      # Upregulated (logFC > lfc)
      adj.P.Val < fdrq & logFC < -lfc ~ "down-regulated",   # Downregulated (logFC < -lfc)
      TRUE ~ "not significant"                              # Non-significant
    ))
}

# Function to generate color palette
get_color_palette <- function(unique_colors) {
  color_map <- c("up-regulated" = "red", "down-regulated" = "blue", "not significant" = "gray")
  return(color_map[unique_colors])  # Adjust palette to the number of unique colors
}

createVolcano <- function(fdrq, lfc, group1, group2, labels = NULL, dynamic = FALSE) {
  
  # Read the differential expression data
  path <- paste0("DEGs/", group1, "_vs_", group2, "_FDRq_1.00_LFC_0.00.tsv")
  data <- readr::read_tsv(file = path)
  
  # Check if the column 'gene_name' exists
  if (!"gene_name" %in% colnames(data)) {
    stop("The 'gene_name' column does not exist in the dataset.")
  }
  
  # Assign colors based on FDRq and logFC
  data <- assign_colors(data, fdrq, lfc)
  
  # Get unique color categories from the data
  unique_colors <- unique(data$color_adjpval)
  
  # Set color palette, make sure it matches the number of unique categories
  color_palette <- get_color_palette(unique_colors)
  
  # Filter top 10 genes by adj.P.Val
  top_genes <- data %>%
    filter(adj.P.Val < fdrq & (logFC > lfc | logFC < -lfc)) %>%
    arrange(adj.P.Val) %>%
    head(10)
  
  # Top 10 down-regulated genes (lowest negative logFC)
  top_downregulated <- data %>%
    filter(adj.P.Val < fdrq & logFC < -lfc) %>%
    arrange(logFC) %>%
    head(10)
  
  # Top 10 up-regulated genes (highest positive logFC)
  top_upregulated <- data %>%
    filter(adj.P.Val < fdrq & logFC > lfc) %>%
    arrange(desc(logFC)) %>%
    head(10)
  
  # Combine all top genes and remove duplicates
  all_top_genes <- unique(rbind(top_genes, top_downregulated, top_upregulated))
  
  # Exclude the user-selected gene from the top gene labels to avoid duplicate labeling
  if (!is.null(labels)) {
    all_top_genes <- all_top_genes %>% filter(gene_name != labels)
  }
  
  # Generate plot
  p <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val), color = color_adjpval, 
                        text = paste("Gene:", gene_name, 
                                     "<br>Adj.P.Val:", round(adj.P.Val, 4),
                                     "<br>LogFC:", round(logFC, 2)))) + # Add `text` aesthetic for tooltips
    geom_point() +
    scale_color_manual(values = color_palette, breaks = c("up-regulated", "down-regulated", "not significant")) +
    labs(x = "log2(Fold Change)", 
         y = "-log10(FDRq)", 
         title = paste0(group1, " vs ", group2),
         color = "User Defined Thresholds") +
    theme_minimal() +
    theme(
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      plot.title = element_text(size = 20, color = "black", face="bold"),
      legend.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 15, color = "black")
    ) +
    
    # Add vertical lines at log2(FC) cutoff and -log2(FC) cutoff
    geom_vline(xintercept = c(-lfc, lfc), linetype = "dashed", color = "darkgray") +
    
    # Add horizontal line at -log10(FDRq cutoff)
    geom_hline(yintercept = -log10(fdrq), linetype = "dashed", color = "darkgray")
  
  # For static plots, annotate the top 10 genes if any meet the thresholds
  if (!dynamic && nrow(all_top_genes) > 0) {
    p <- p + geom_text_repel(data = all_top_genes, 
                             aes(label = gene_name), 
                             max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
                             fontface="italic",
                             size = 4, 
                             color = "black")
  }
  
  # If a specific gene is provided, label it in static mode only
  if (!is.null(labels) && labels %in% data$gene_name && !dynamic) {
    print(paste("Highlighting gene:", labels))  # Debugging
    p <- p + geom_text_repel(
      data = data %>% filter(gene_name == labels),  # Filter data for the specific gene
      aes(label = gene_name), 
      size = 5, 
      fontface = "italic",
      max.overlaps = getOption("ggrepel.max.overlaps", default = 30),
      color = "black"
    )
  } else if (!is.null(labels) && !dynamic) {
    print(paste("Gene", labels, "not found in dataset"))
  }
  
  return(p)
}


# function to plot boxplot
plotBoxplot <- function(counts, gene) {
  
  # create df
  names <- colnames(counts)
  group <- str_match(names,"([SLH7248hd\\.]+)\\.[MF]\\.[0-9]+")[,2]
  value <- as.numeric(subset(counts, rownames(counts) == gene))
  df <- data.frame(sample = names,
                   group = factor(group, levels = c("S.8h","L.8h","H.8h",
                                                    "S.24h","L.24h","H.24h",
                                                    "S.7d","L.7d","H.7d")),
                   value = as.vector(value))
  df <- df[order(df$group),]
  rownames(df) <- 1:nrow(df)
  
  # Visualize the distribution of genes detected per sample via boxplot
  group_colors <- c("pink","lightyellow","lightblue",
                    "red","gold","cornflowerblue",
                    "darkred", "gold3","navy")
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0)) +
    ggtitle(gene) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    ylab("CPM") + xlab("Group") +
    theme(
      axis.text = element_text(size = 15, color = "black"),
      axis.title = element_text(size = 15, color = "black"),
      plot.title = element_text(size = 20, color = "black", face = "bold"),
      legend.title = element_text(size = 15, color = "black"),
      legend.text = element_text(size = 15, color = "black")
    )
  
} # end boxplot function

# Shiny app structure
ui <- fluidPage(
  titlePanel("Psilocybin Project 1 Bulk RNA-sequencing"),
  sidebarLayout(
    sidebarPanel(
      selectInput("group1", "Select Group 1", choices = group1_options, selected = "H.8h"),  # Default to H.8h
      uiOutput("group2_ui"),  # This will be dynamically generated based on group1
      numericInput("fdrq_cutoff", "FDRq Cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
      numericInput("lfc_cutoff", "log2(FC) Cutoff", value = 0.2, min = 0, max = 5, step = 0.1),
      # Use selectizeInput for gene selection with 5 suggestions as user types for the Volcano Plot
      selectizeInput("gene_label", 
                     "Gene to Highlight (Volcano Plot)", 
                     choices = geneOptions,
                     options = list(maxOptions = 5), 
                     selected = "Fmo2"),  # Default to Fmo2 for the Volcano Plot
      radioButtons("plot_type", "Select Plot Type", choices = c("Static", "Dynamic")),
      actionButton("generate_plot", "Generate Plot")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Volcano Plot",
                 p("Using data from all_samples/both_sexes/DEGs_with_hbb"),
                 
                 # Conditional Panels for Static and Dynamic plot
                 conditionalPanel(
                   condition = "input.plot_type == 'Static'",
                   plotOutput("volcano_plot_static")
                 ),
                 conditionalPanel(
                   condition = "input.plot_type == 'Dynamic'",
                   plotlyOutput("volcano_plot_dynamic")
                 )
        ),
        
        # CPM Tab (removed redundant gene selection input from this tab)
        tabPanel("CPM", 
                 p("Using CPM post-filtering data"),
                 plotOutput("cpm_data")  # Output for the CPM boxplot
        )
      )
    )
  )
)

# Server logic
server <- function(input, output, session) {
  # Dynamically generate group2 options based on group1 selection
  output$group2_ui <- renderUI({
    group2_options <- comparisons %>% filter(group1 == input$group1) %>% pull(group2)
    selectInput("group2", "Select Group 2", choices = group2_options, selected = "S.8h")  # Default to S.8h
  })
  
  # Plot rendering
  observeEvent(input$generate_plot, {
    print(paste("Selected Gene:", input$gene_label))  # Debugging to check the input
    
    if (input$plot_type == "Static") {
      output$volcano_plot_static <- renderPlot({
        createVolcano(input$fdrq_cutoff, input$lfc_cutoff, input$group1, input$group2, input$gene_label, dynamic = FALSE)
      })
      output$volcano_plot_dynamic <- NULL
    } else {
      output$volcano_plot_dynamic <- renderPlotly({
        p <- createVolcano(input$fdrq_cutoff, input$lfc_cutoff, input$group1, input$group2, dynamic = TRUE)
        ggplotly(p, tooltip = "text")  # Ensure the tooltip shows the custom text aesthetic
      })
      output$volcano_plot_static <- NULL
    }
  })
  
  # CPM data display - use plotBoxplot function to plot boxplot for the selected gene
  output$cpm_data <- renderPlot({
    req(input$gene_label)  # Make sure a gene is selected
    
    # Call the plotBoxplot function to generate the boxplot
    plotBoxplot(cpm, input$gene_label)
  })
  
}

shinyApp(ui = ui, server = server)
