library(shiny)

ui <- fluidPage(
  titlePanel("Isoform TPM Viewer"),
  
  br(),
  
  fluidRow(
    style = "background-color: #F8F8F8",
    
    column(
      width = 3,
      selectizeInput(
        inputId = "gene_name",
        label = "Select gene:",
        choices = NULL,
        selected = "Dbp",
        options = list(placeholder = "Type to search a gene...", maxOptions = 5000)
      )
    ),
    
    column(
      width = 3,
      radioButtons(
        inputId  = "tpm_level",
        label    = "TPM var to aggregate by (plot 1 option)",
        choices  = c("group", "group2"),
        selected = "group",
        inline   = TRUE
      )
    ),
    
    column(
      width = 6,
      radioButtons(
        inputId  = "tx_mode",
        label    = "Transcript types (plot 2 option)",
        choices  = c("protein coding only", "all transcript types"),
        selected = "protein coding only",
        inline   = TRUE
      )
    )
  ),
  
  br(),
  
  fluidRow(
    column(width = 12, plotOutput("isoform_tpm_bar", height = "420px"))
  ),
  
  br(),
  
  fluidRow(
    column(width = 12, plotOutput("ensembl_transcripts", height = "75vh"))
  ),
  
  br(),
  fluidRow(
    column(width = 12, dataTableOutput("ensembl_table"))
  )
  
)
