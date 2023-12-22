# load input options
comparisons <- readr::read_tsv("input_options.tsv")
cpm <- read.delim2(file = "counts/CPM_prefiltering.tsv")
group1_options <- comparisons$group1
group_order <- c("sal.8h","sal.24h","sal.7d",
                 "psilo.low.8h","psilo.low.24h","psilo.low.7d",
                 "psilo.high.8h","psilo.high.24h","psilo.high.7d")
group_colors <- c("gray90","gray60","gray30","lightblue","cornflowerblue","blue",
                          "sienna1","red2","red4")

# function to plot boxplot
plotBoxplot <- function(counts, gene) {
  
  # create df
  names <- colnames(counts)
  group <- str_match(names,"(.+).[MF].[0-9]+")[,2]
  value <- as.numeric(subset(counts, rownames(counts) == gene))
  df <- data.frame(sample = names,
                   group = factor(group),
                   value = as.vector(value))
  rownames(df) <- 1:nrow(df)
  
  # Visualize the distribution of genes detected per sample via boxplot
  df$group <- factor(df$group, levels = group_order)
  ggplot(df, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter() + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0)) +
    ggtitle(gene) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    ylab("CPM") + xlab("Group")
  
} # end boxplot function


# function to plot bar graph
plotBar <- function(counts, gene) {
  
  # create df
  names <- colnames(counts)
  group <- str_match(names,"(.+).[MF].[0-9]+")[,2]
  value <- as.numeric(subset(counts, rownames(counts) == gene))
  df <- data.frame(sample = names,
                   group = factor(group),
                   value = as.vector(value))
  rownames(df) <- 1:nrow(df)
  df$group <- factor(df$group, levels = group_order)
  
  # Visualize the distribution of genes detected per sample via boxplot
  b <- ggplot(df, aes(x = factor(names, levels = unique(names)), y = value, fill = group)) +
    geom_bar(stat = "identity") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust = 0)) +
    ggtitle(gene) +
    scale_color_manual(values = group_colors) +
    scale_fill_manual(values = group_colors) +
    ylab("CPM") + xlab("Sample")
  b
} # end boxplot function


# server
server <- function(input, output) {
  
  observe({
    group2_options <- comparisons |>
      filter(group1 == input$group1) |>
      pull(group2) |>
      sort()
    
    updateSelectInput(
      inputId = "group2", 
      choices = group2_options
    )
  })
  
  output$volcano <- renderPlotly({
    # set variables
    req(input$fdrq,input$lfc,input$group1,input$group2)
    fdrq <- input$fdrq
    lfc <- input$lfc
    
    # read file
    path <- paste0("DEGs/",input$group1,"_vs_",input$group2,"_FDRq_1.00.tsv")
    data <- readr::read_tsv(file = path)
    
    # assign colors
    color_values <- vector()
    max <- nrow(data)
    for(row in 1:max){
      if (data$adj.P.Val[row] < fdrq){
        if (data$logFC [row] > lfc){
          color_values <- c(color_values, 1) # 1 when logFC > input.lfc and FDRq < input.fdrq
        }
        else if (data$logFC[row] < -lfc){
          color_values <- c(color_values, 2) # 2 when logFC < input.lfc and FDRq < input.fdrq
        }
        else {
          color_values <- c(color_values, 3) # 3 when input.lfc NOT met and FDRq < input.fdrq
        }
      }
      else{
        color_values <- c(color_values, 3) # 3 when input.fdrq NOT met
      }
    }
    data$color_adjpval <- factor(color_values)
    
    # plot title
    plot_title <- paste0(input$group1, " vs. ", input$group2, ", FDRq < ", fdrq, ", Log2(FC) = ", lfc)
    
    # subset genes to label
    up <- data[data$color_adjpval == 1,]
    up15 <- up[1:15,]
    down <- data[data$color_adjpval == 2,]
    down15 <- down[1:15,]
    
    # set manual colors
    if (!1 %in% unique(data$color_adjpval)) {
      my_colors <- c("blue","gray")
    } else if (!2 %in% unique(data$color_adjpval)) {
      my_colors <- c("red","gray")
    } else if (!1 %in% unique(data$color_adjpval) && !2 %in% unique(data$color_adjpval)) {
      my_colors <- c("gray")
    } else {
      my_colors <- c("red","blue","gray")
    }
    
    # set significance threshold
    hadjpval <- (-log10(max(data$P.Value[data$adj.P.Val < fdrq], na.rm=TRUE)))
    
    # plot
    p <-
      ggplot(data = data, 
             aes(label = gene_name,
                 label2 = adj.P.Val)) +
      geom_point(aes(x = logFC, y = -log10(P.Value), color = color_adjpval),
                 alpha = 0.8, size = 2) +  # create scatterplot, alpha makes points transparent
      theme_bw() +  # set color theme
      theme(legend.position = "none") +  # no legend
      scale_color_manual(values = my_colors) +  # set factor colors
      labs(
        title = plot_title, # main title
        x = "Log2(FC)", # x-axis title
        y = "-Log10(p-value)" # y-axis title
      ) +
      theme(axis.title.x = element_text(size = 15),
            axis.text.x = element_text(size = 15)) +
      theme(axis.title.y = element_text(size = 15),
            axis.text.y = element_text(size = 15)) +
      theme(plot.title = element_text(size = 15)) +
      geom_hline(yintercept = hadjpval,  #  horizontal line
                 colour = "#000000",
                 linetype = "dashed") +
      geom_vline(xintercept = lfc,  #  vertical line
                 colour = "#000000",
                 linetype = "dashed") +
      geom_vline(xintercept = -lfc,  #  vertical line
                 colour = "#000000",
                 linetype = "dashed")
    plotly::ggplotly(p, tooltip = c("gene_name","logFC","adj.P.Val"))
  })
  
  # box plot
  output$boxplot <- renderPlot({
    req(input$goi)
    plotBoxplot(counts = cpm, gene = input$goi)
  })
  
  # bar plot
  output$bar <- renderPlot({
    req(input$goi)
    plotBar(counts = cpm, gene = input$goi)
  })
  
  output$manhattan <- renderPlotly({
    req(input$gprofilerInput)
    # read list
    data <- as.data.frame(read.table(paste0("gprofilerInput/",input$gprofilerInput), sep = "\t"))
    data <- data$V1
    # query
    gost.up <- gost(query = data,
                    organism = "hsapiens",
                    ordered_query = TRUE)
    # plot
    gostplot(gostres = gost.up, capped = FALSE)
  })
  
}
