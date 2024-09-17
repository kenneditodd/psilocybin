# load libraries
library(dplyr)
library(ggrepel)
library(gprofiler2)
library(plotly)
library(shiny)
library(stringr)

# load input options
comparisons <- readr::read_tsv("input_options.tsv")
genes <- readRDS("genes.rds")
group1_options <- comparisons$group1
gprofilerInputFiles <- readRDS("gprofilerInputFiles.rds")

# user interface
ui <- fluidPage(
  
  # application title
  titlePanel("Psilocybin Project 1 Bulk RNA-Sequencing"),
  
  tabsetPanel(
    
    ############################### VOLCANO TAB ################################
    tabPanel(
      title = "volcano",
      # Sidebar with two select input options 
      sidebarLayout(
        sidebarPanel(
          selectInput(inputId = "group1",
                      label = "Select group 1",
                      choices = group1_options,
                      selected = "psilo.low.7d"),
          selectInput(inputId = "group2",
                      label = "Select group 2",
                      choices = c(),
                      selected = NULL),
          numericInput(inputId = "fdrq",
                       label = "Select FDRq cutoff",
                       value = 0.05,
                       min = 0,
                       max = 1,
                       step = 0.01),
          numericInput(inputId = "lfc",
                       label = "Select log2(FC) cutoff",
                       value = 0,
                       step = 0.1)
        ),
        
        # plot
        mainPanel(
          h5("Samples: all samples"),
          h5("Filtering: strict, all samples in one treatment (sal or psi) must have CPM >= 1"),
          h5("Model: group + sex + Hbb-bs"),
          plotlyOutput("volcano")
        )
        
      ) # end sidebarLayout
    ), # end tabPanel - volcano
      
    ################################# CPM TAB ##################################
    tabPanel(
      title = "CPM",
      # Sidebar with a slider input for number of bins 
      sidebarLayout(
        sidebarPanel(
          # Populate gene options based on selected tissue
          selectizeInput(inputId = "goi",
                         label = "Select a gene",
                         choices = genes,
                         selected = "Hbb-bs",
                         options = list(maxOptions = 5))
        ), # end sidebarPanel
        
        # Show a plot of the generated distribution
        mainPanel(
          h3("The postfiltering CPM table is used for these graphs."),
          plotOutput("boxplot"),
          br(),
          plotOutput("bar")
        )
      ) # end SidebarLayout
    ), # end tabPanel - CPM
    
    ############################### GPROFILER TAB ###############################
    tabPanel(
      title = "gprofiler2",
      sidebarLayout(
        sidebarPanel(
          selectInput(
            inputId = "gprofilerInput",
            label = "Select gene list for input",
            choice = gprofilerInputFiles)
        ), # end sidebarPanel
        mainPanel(
          h3("Gprofiler2"),
          p("Gprofiler2 performs a functional enrichment analysis on an input gene list. 
            It maps genes to known functional information sources and detects statistically significantly enriched terms. 
            In addition to Gene Ontology, pathways from KEGG Reactome and WikiPathways are included; miRNA targets from miRTarBase and regulatory motif matches from TRANSFAC; tissue specificity from Human Protein Atlas; protein complexes from CORUM and human disease phenotypes from Human Phenotype Ontology.
            GO hierarchy: MF = molecular function, BP = biological process, CC = cellular component"),
          br(),
          p("May take a minute to load. Some gene lists contain too few genes to 
            return significant results. One unique feature of Gprofiler2 over 
            metascape is it takes into account the order of importance in the 
            gene list. So, it takes into account that the genes at the top of the 
            list are more significant than the genes in the bottom of the list."),
          br(),
          p("The plot below is called a Manhattan plot. The size of the points 
            refer to how many genes are in the pathway which is also listed in 
            parenthesis if you hover over a point. The x-axis has various 
            databases. The y-axis shows the adjusted p-value from the functional 
            enrichment analysis query."),
          br(),
          p("Using FDRq < 0.1 and |LFC| > 0.2"),
          plotlyOutput("manhattan")
        )
      ) # sidebarPanel
    ) # end gprofiler2 tabPanel
  ) # end of tabsetPanel
  

) # end of fluidPage
