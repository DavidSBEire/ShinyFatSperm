library(shiny)
library(shinydashboard)
library(ggVennDiagram)
library(ggiraph)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(eulerr)

# Load and prepare datasets
hfda_acute_df <- read.csv("HFD_acute.csv", stringsAsFactors = FALSE, check.names = FALSE) %>% 
  select(Gene, `Gene Name`, miRNA, Human)

hfda_chronic_df <- read.csv("HFD_chronic.csv", stringsAsFactors = FALSE, check.names = FALSE) %>% 
  select(Gene, `Gene Name`, miRNA, Human)

bodycomp_df  <- read.csv("Bodycomposition.csv", stringsAsFactors = FALSE, check.names = FALSE) %>% 
  select(Gene, `Gene Name`, Phenotype)
glucose_df   <- read.csv("Glucose.csv",         stringsAsFactors = FALSE, check.names = FALSE) %>% 
  select(Gene, `Gene Name`, Phenotype)
insulin_df   <- read.csv("Insulin.csv",         stringsAsFactors = FALSE, check.names = FALSE) %>% 
  select(Gene, `Gene Name`, Phenotype)



# Build a Gene → GeneName lookup (first occurrence)
bodycomp_names <- bodycomp_df  %>% select(Gene, GeneName = `Gene Name`) %>% distinct()
glucose_names  <- glucose_df %>% select(Gene, GeneName = `Gene Name`) %>% distinct()
insulin_names  <- insulin_df    %>% select(Gene, GeneName = `Gene Name`) %>% distinct()
hfda_acute_names     <- hfda_acute_df    %>% select(Gene, GeneName = `Gene Name`) %>% distinct()
hfda_chronic_names     <- hfda_chronic_df    %>% select(Gene, GeneName = `Gene Name`) %>% distinct()

gene_info_acute <- bind_rows(bodycomp_names, glucose_names, insulin_names, hfda_acute_names) %>%
  distinct(Gene, .keep_all = TRUE)

gene_info_chronic <- bind_rows(bodycomp_names, glucose_names, insulin_names, hfda_acute_names) %>%
  distinct(Gene, .keep_all = TRUE)

# IPA input tables
acute_IPA_df   <- read.csv("acute_IPA.csv",   stringsAsFactors = FALSE, check.names = FALSE)
chronic_IPA_df <- read.csv("chronic_IPA.csv", stringsAsFactors = FALSE, check.names = FALSE)

# Combined heatmap files
acute_pathways_df   <- read.csv("acute_pathways.csv",   stringsAsFactors = FALSE, check.names = FALSE)
acute_functions_df  <- read.csv("acute_functions.csv",  stringsAsFactors = FALSE, check.names = FALSE)
chronic_pathways_df <- read.csv("chronic_pathways.csv", stringsAsFactors = FALSE, check.names = FALSE)
chronic_functions_df<- read.csv("chronic_functions.csv",stringsAsFactors = FALSE, check.names = FALSE)

##############################################
# UI
##############################################
ui <- navbarPage(
  title    = "ShinyFatSperm",
  id       = "tabs",
  position = "fixed-top",
  inverse  = TRUE,
  header   = tags$head(tags$style(HTML("
    .navbar { background-color: #000000 !important; }
    .navbar .navbar-brand,
    .navbar .navbar-nav > li > a { color: white !important; }
    body { padding-top: 70px; }
    .custom-button-style { width:100%; background:#000000; color:white; border-radius:10px; padding:20px; text-align:center; }
    .custom-icon { font-size:3em; color:white; display:block; margin:auto; }
    .icon-heading, .icon-description { color:white; }
  "))),
  
  ########################################
  # Homepage
  ########################################
  tabPanel(
    title = "Homepage",
    value = "home",
    fluidRow(
      column(12, align="center", style="margin-top:100px;",
             img(src = "Logo.png", 
                 style = "max-width: 300px; height: 240px; margin-bottom: 10px;border-radius:30%;"),
             h2("Welcome to ShinyFatSperm"),
             p("Explore the potential links between paternal sperm epigenomes and obesity phenotypes.")
      )
    ),
    fluidRow(style="margin:30px 0;",
             column(3,
                    actionButton("to_Acute", class="custom-button-style",
                                 div(icon("bolt", class="custom-icon"),
                                     strong("Acute Overlap", class="icon-heading"),
                                     p("HFD acute vs phenotypes", class="icon-description")
                                 )
                    )
             ),
             column(3,
                    actionButton("to_Chronic", class="custom-button-style",
                                 div(icon("clock", class="custom-icon"),
                                     strong("Chronic Overlap", class="icon-heading"),
                                     p("HFD chronic vs phenotypes", class="icon-description")
                                 )
                    )
             ),
             column(3,
                    actionButton("to_IPA", class="custom-button-style",
                                 div(icon("microscope", class="custom-icon"),
                                     strong("In-Silico Data", class="icon-heading"),
                                     p("Pathway & function enrichments", class="icon-description")
                                 )
                    )
             ),
             column(3,
                    actionButton("to_Team", class="custom-button-style",
                                 div(icon("info-circle", class="custom-icon"),
                                     strong("Meet the Team", class="icon-heading"),
                                     p("Get in contact", class="icon-description")
                                 )
                    )
             )
    ),
    # Footer / references
    fluidRow(
      column(12, align = "center", style = "margin-top: 5px; color: #666;",
             p(
               "If you use this application or its data in your work, please cite the following:",
               br(), 
               a(
                 HTML("Laurent, K., et al. (2026) Epigenetic Legacy: The Role of Sperm miRNAs in the Paternal Inheritance of Diabetes and Obesity Development, <i>Diabetes/Metabolism Research and Reviews</i>"),
                 href = "https://onlinelibrary.wiley.com/doi/10.1002/dmrr.70157", target = "_blank"
               ),
               br(), br(),  # adds spacing
               p(
                 HTML('Full code is available via <a href="https://github.com/DavidSBEire/ShinyFatSperm" target="_blank">GitHub</a>')
               )
             )
      )
    )
  ),
  
  ########################################
  # Acute Overlap
  ########################################
  tabPanel(
    title = "Acute Overlap",
    value = "acute",
    sidebarLayout(
      sidebarPanel(width=3,
                   style = "background: #3A3A3A;
                            color: #FFFFFF;
                            border: 1px solid #000000;
                            border-radius: 12px;
                            padding: 15px;",
                   strong("Instructions:"),
                   tags$ul(
                     style = "color: #FFFFFF; margin-top: 10px;",
                     tags$li("Click a Venn intersect to explore the overlaps"),
                     tags$li("Download data as CSV or Excel file"),
                     tags$li(
                       tags$em("In-silico"),
                       " Workflow shown below"
                     )
                   ),
                   # Human-mapped filter checkbox
                   tags$hr(style = "border-color: #666;"),
                   checkboxInput(
                     inputId = "acute_human_only",
                     label   = "Show human-mapped genes only",
                     value   = FALSE
                   ),
                   img(src="Workflow_Acute.png", width="100%", style="margin-top:20px;")
      ),
      mainPanel(width=9,
                girafeOutput("vennAcute", width="700px", height="700px"),
                hr(),
                DTOutput("geneTableAcute"),
      )
    )
  ),
  
  ########################################
  # Chronic Overlap
  ########################################
  tabPanel(
    title = "Chronic Overlap",
    value = "chronic",
    sidebarLayout(
      sidebarPanel(width=3,
                   style = "background: #3A3A3A;
                            color: #FFFFFF;
                            border: 1px solid #000000;
                            border-radius: 12px;
                            padding: 15px;",
                   strong("Instructions:"),
                   tags$ul(
                     style = "color: #FFFFFF; margin-top: 10px;",
                     tags$li("Click a Venn intersect to explore the overlaps"),
                     tags$li("Download data as CSV or Excel file"),
                     tags$li(
                       tags$em("In-silico"),
                       " Workflow shown below"
                     )
                   ),
                   # Human-mapped filter checkbox
                   tags$hr(style = "border-color: #666;"),
                   checkboxInput(
                     inputId = "chronic_human_only",
                     label   = "Show human-mapped genes only",
                     value   = FALSE
                   ),
                   
                   img(src="Workflow_Chronic.png", width="100%", style="margin-top:20px;")
      ),
      mainPanel(width=9,
                girafeOutput("vennChronic", width="700px", height="700px"),
                hr(),
                DTOutput("geneTableChronic"),
      )
    )
  ),
  
  
  ########################################
  # In-Silico placeholder
  ########################################
  tabPanel(
    title = "In-Silico Data",
    value = "In-Silico Data",
    tags$div(
      style = "
      padding: 10px;
      background: #3A3A3A;
      color: #FFFFFF;
      border: 1px solid #FFFFFF;
      border-radius: 12px;
      margin-bottom: 15px;
    ",
      strong("Instructions:"),
      tags$ul(
        tags$li("Select your phenotype of interest — your choice live updates both Acute & Chronic analyses"),
        tags$li("Download the full annotated table, including all mapped pathways & functions"),
        tags$li("Search by gene symbol or filter by molecule type — this will update your download"),
        tags$li("Full download takes a moment, please be patient")
      )
    ),
    
    #Global phenotype selector
    fluidRow(
      column(4, offset = 4,
             selectInput(
               inputId = "inSilico_pheno",
               label   = "Choose phenotype group:",
               choices = c("Bodycomposition", "Glucose", "Insulin"),
               selected = "Bodycomposition",
               width = "100%"
             )
      )
    ),
    
    #Two columns: Acute vs Chronic
    fluidRow(
      
      #Acute HFD column
      column(6,
             wellPanel(
               style = "background-color: #F6D8ED;
                 border: 1px solid #9F5590;
                 border-radius: 8px;
                 padding: 15px;",
             h4("Acute Overlap"),
             fluidRow(
               column(4,
                      textInput(
                        inputId = "acute_search",
                        label   = "Search by Gene symbol:",
                        value   = ""
                      ),
                      selectInput(
                        inputId = "acute_type",
                        label   = "Filter by Type:",
                        choices = sort(unique(acute_IPA_df$Type)),
                        multiple = TRUE,
                        width = "100%"
                      )
               ),
               column(8, align="centre",
                      img(src = "Acute.png", width = "30%", style = "margin-bottom:10px;")
               )
             ),
             DTOutput("acute_table"),
             downloadButton("downloadAcute", "Download Acute CSV"),
             fluidRow(
               column(6, plotlyOutput("acute_pathways_heatmap", height = "300px")),
               column(6, plotlyOutput("acute_functions_heatmap", height = "300px"))
             )
             )
      ),
      
      #Chronic HFD column
      column(6,
             wellPanel(
               style = "background-color: #D3EEEA;
                 border: 1px solid #33887E;
                 border-radius: 8px;
                 padding: 15px;",
             h4("Chronic Overlap"),
             fluidRow(
               column(4,
                      textInput(
                        inputId = "chronic_search",
                        label   = "Search by Gene symbol:",
                        value   = ""
                      ),
                      selectInput(
                        inputId = "chronic_type",
                        label   = "Filter by Type:",
                        choices = sort(unique(chronic_IPA_df$Type)),
                        multiple = TRUE,
                        width = "100%"
                      )
               ),
               column(8, align="center",
                      img(src = "Chronic.png", width = "40%", style = "margin-bottom:10px;")
               )
             ),
             DTOutput("chronic_table"),
             downloadButton("downloadChronic", "Download Chronic CSV"),
             fluidRow(
               column(6, plotlyOutput("chronic_pathways_heatmap", height = "300px")),
               column(6, plotlyOutput("chronic_functions_heatmap", height = "300px"))
             )
      )
      )
    )
  ),
  
  ########################################
  # Team
  ########################################
  tabPanel("Team",
           # Section header
           fluidRow(
             column(12, align = "center",
                    h2("About the Team"),
                    p("Our team consists of experts in reproductive biology, epigenetics, and bioinformatics. Click below to get in touch or learn more about us.")
             )
           ),
           
           # Three team members side-by-side
           fluidRow(
             # Dr. David Skerrett-Byrne
             column(4, align = "center",
                    div(class = "team-member", style = "margin-bottom: 30px;",
                        img(src = "David.jpg", height = "300px", style = "border-radius:50%;"),
                        h4("Dr. David Skerrett-Byrne"),
                        p(em("NMHRC Research Fellow"), style = "margin:5px 0;"),
                        p("Institute of Experimental Genetics", style = "margin:2px 0;"),
                        p("Helmholtz Zentrum München",        style = "margin:2px 0;"),
                        p("Hunter Medical Research Institute",style = "margin:2px 0;"),
                        p("University of Newcastle (Australia)",style = "margin:2px 0;"),
                        div(
                          style = "margin-top:10px;",
                          a(icon("linkedin"), href = "https://www.linkedin.com/in/david-skerrett-byrne-17596737/",
                            target = "_blank", style = "margin-right:10px;"),
                          a(icon("envelope"), href = "mailto:David.Skerrett-Byrne@helmholtz-munich.de",
                            style = "margin-right:10px;"),
                          a(icon("globe"), href = "https://www.newcastle.edu.au/profile/david-skerrett-byrne",
                            target = "_blank")
                        )
                    )
             ),
             
             # Katharina Laurent
             column(4, align = "center",
                    div(class = "team-member", style = "margin-bottom: 30px;",
                        img(src = "Katharina.jpg", height = "300px", style = "border-radius:50%;"),
                        h4("Katharina Laurent"),
                        p(em("Doctoral Researcher"), style = "margin:5px 0;"),
                        p("Gene Regulation and Epigenetics Group", style = "margin:2px 0;"),
                        p("Institute of Experimental Genetics", style = "margin:2px 0;"),
                        p("Helmholtz Zentrum München",         style = "margin:2px 0;"),
                        div(
                          style = "margin-top:10px;",
                          a(icon("linkedin"), href = "https://www.linkedin.com/in/katharina-laurent-a24380206/",
                            target = "_blank", style = "margin-right:10px;"),
                          a(icon("envelope"), href = "mailto:katharina.laurent@helmholtz-munich.de",
                            style = "margin-right:10px;"),
                          a(icon("globe"), href = "https://www.helmholtz-munich.de/en/ieg/research-groups/gene-regulation-and-epigenetics",
                            target = "_blank")
                        )
                    )
             ),
             
             # Prof. Johannes Beckers
             column(4, align = "center",
                    div(class = "team-member", style = "margin-bottom: 30px;",
                        img(src = "Johannes.jpg", height = "300px", style = "border-radius:50%;"),
                        h4("Prof. Johannes Beckers"),
                        p(em("Principal Investigator, Gene Regulation and Epigenetics Group"), style = "margin:5px 0;"),
                        p(em("Deputy Director of the Institute of Experimental Genetics"), style = "margin:5px 0;"),
                        p("Helmholtz Zentrum München",        style = "margin:2px 0;"),
                        div(
                          style = "margin-top:10px;",
                          a(icon("linkedin"), href = "https://www.linkedin.com/in/prof-dr-johannes-beckers-executive-mba-b2aa542b/",
                            target = "_blank", style = "margin-right:10px;"),
                          a(icon("envelope"), href = "mailto:johannes.beckers@helmholtz-munich.de",
                            style = "margin-right:10px;"),
                          a(icon("globe"), href = "https://www.helmholtz-munich.de/en/pi-3-21",
                            target = "_blank")
                        )
                    )
             )
             
           ),
           
           # Footer / references
           fluidRow(
             column(12, align = "center", style = "margin-top: 5px; color: #666;",
                    h2("Referencing"),
                    p("If you use this application or its data in your work, please cite the following:"),
                    br(), 
                    a(
                      HTML("Laurent, K., et al. (2026) Epigenetic Legacy: The Role of Sperm miRNAs in the Paternal Inheritance of Diabetes and Obesity Development, <i>Diabetes/Metabolism Research and Reviews</i>"),
                      href = "https://onlinelibrary.wiley.com/doi/10.1002/dmrr.70157", target = "_blank"
                    ),
                    br(), br(),  # adds spacing
                    p(
                      HTML('Full code is available via <a href="https://github.com/DavidSBEire/ShinyFatSperm" target="_blank">GitHub</a>')
                    )
             )
           )
  )
)

##############################################
# Server
##############################################

server <- function(input, output, session) {

  
  ##########################
  # Navigation (Homepage) 
  ##########################
  observeEvent(input$to_Acute,   updateTabsetPanel(session, "tabs", "acute"))
  observeEvent(input$to_Chronic, updateTabsetPanel(session, "tabs", "chronic"))
  observeEvent(input$to_IPA,     updateTabsetPanel(session, "tabs", "In-Silico Data"))
  observeEvent(input$to_Team,    updateTabsetPanel(session, "tabs", "Team"))
  
  ##################################
  # Overlap Tab: Venn & Gene Table #
  ##################################
  
  # Extract unique Gene vectors
  genesHFDacute_r <- reactive({
    df <- hfda_acute_df
    if (isTRUE(input$acute_human_only)) {
      df <- df %>% filter(Human == "✓")
    }
    unique(df$Gene)
  })
  
  genesHFDchronic_r <- reactive({
    df <- hfda_chronic_df
    if (isTRUE(input$chronic_human_only)) {
      df <- df %>% filter(Human == "✓")
    }
    unique(df$Gene)
  })
  genesBodycomp   <- unique(bodycomp_df$Gene)
  genesGlucose    <- unique(glucose_df$Gene)
  genesInsulin    <- unique(insulin_df$Gene)
  
  makeVenn <- function(g1, g2, g3, g4) {
    xlist <- list(HFD = g1, Bodycomposition = g2, Glucose = g3, Insulin = g4)
    vo    <- Venn(xlist)
    vd    <- process_data(vo)
    edges <- venn_regionedge(vd)   %>% filter(name != "")
    labs  <- venn_regionlabel(vd)  %>% filter(name != "")
    sets  <- venn_setlabel(vd)
    p <- ggplot() +
      geom_polygon_interactive(
        data = edges,
        aes(x = X, y = Y, group = id, fill = name, data_id = name),
        color = "black", size = 0.6
      ) +
      scale_fill_manual(
        values = c(
          Bodycomposition = "#ABABAB",
          Glucose         = "#ABABAB",
          Insulin         = "#ABABAB",
          HFD             = "#6E8C9A"
        ), guide = FALSE
      ) +
      geom_text(data = sets, aes(X, Y, label = name), color = "black", size = 6, fontface = "italic") +
      geom_text(data = labs, aes(X, Y, label = count), size = 5) +
      coord_fixed() +
      theme_void(base_size = 14) +
      theme(plot.margin = margin(20, 20, 20, 20))
    girafe(
      ggobj     = p,
      width_svg  = 7,
      height_svg = 7,
      options    = list(
        opts_hover(css = "fill: rgba(52,151,143,0.4); cursor:pointer;"),
        opts_selection(
          type       = "single",
          only_shiny = FALSE,
          css        = "fill: rgba(250,190,0,0.4); stroke:#CD9C02;"
        )
      )
    )
  }
  
  membershipAcute <- reactive({
    allG <- unique(c(genesHFDacute_r(), genesBodycomp, genesGlucose, genesInsulin))
    df <- data.frame(
      Gene            = allG,
      HFD             = allG %in% genesHFDacute_r(),
      Bodycomposition = allG %in% genesBodycomp,
      Glucose         = allG %in% genesGlucose,
      Insulin         = allG %in% genesInsulin,
      stringsAsFactors = FALSE
    )
    left_join(df, gene_info_acute, by = "Gene")
  })
  
  output$vennAcute <- renderGirafe({
    req(membershipAcute())
    makeVenn(genesHFDacute_r(), genesBodycomp, genesGlucose, genesInsulin)
  })
  
  filteredAcute <- reactive({
    sel_str <- input$vennAcute_selected
    req(!is.null(sel_str) && nzchar(sel_str))
    sel   <- intersect(strsplit(sel_str, "/", TRUE)[[1]],
                       c("HFD","Bodycomposition","Glucose","Insulin"))
    other <- setdiff(c("HFD","Bodycomposition","Glucose","Insulin"), sel)
    
    # base membership + gene names
    df0 <- membershipAcute() %>% as.data.frame()
    keep <- Reduce(`&`, lapply(sel, function(x) df0[[x]]), init = TRUE) &
      !Reduce(`|`, lapply(other, function(x) df0[[x]]), init = FALSE)
    out  <- df0[keep, c("Gene","GeneName")]
    
    # add miRNA if HFD selected
    if ("HFD" %in% sel) {
      mi <- hfda_acute_df %>%
        filter(Gene %in% out$Gene) %>%
        group_by(Gene) %>%
        summarize(miRNA = paste(unique(miRNA), collapse = "; "), .groups="drop")
      out <- left_join(out, mi, by="Gene")
    }
    
    # now collapse all selected Phenotypes into one column
    pheno_df <- bind_rows(
      bodycomp_df   %>% mutate(Source = "Bodycomposition"),
      glucose_df    %>% mutate(Source = "Glucose"),
      insulin_df    %>% mutate(Source = "Insulin")
    )
    phenos <- pheno_df %>%
      filter(Gene %in% out$Gene, Source %in% sel) %>%
      group_by(Gene) %>%
      summarize(
        Phenotypes = paste(unique(Phenotype), collapse = "; "),
        .groups="drop"
      )
    
    left_join(out, phenos, by="Gene")
  })
  
  membershipChronic <- reactive({
    allG <- unique(c(genesHFDchronic_r(), genesBodycomp, genesGlucose, genesInsulin))
    df <- data.frame(
      Gene            = allG,
      HFD             = allG %in% genesHFDchronic_r(),
      Bodycomposition = allG %in% genesBodycomp,
      Glucose         = allG %in% genesGlucose,
      Insulin         = allG %in% genesInsulin,
      stringsAsFactors = FALSE
    )
    left_join(df, gene_info_chronic, by = "Gene")
  })
  
  output$vennChronic <- renderGirafe({
    req(membershipChronic())
    makeVenn(genesHFDchronic_r(), genesBodycomp, genesGlucose, genesInsulin)
  })
  
  
  filteredChronic <- reactive({
    sel_str <- input$vennChronic_selected
    req(!is.null(sel_str) && nzchar(sel_str))
    sel   <- intersect(strsplit(sel_str, "/", TRUE)[[1]],
                       c("HFD","Bodycomposition","Glucose","Insulin"))
    other <- setdiff(c("HFD","Bodycomposition","Glucose","Insulin"), sel)
    
    df0 <- membershipChronic() %>% as.data.frame()
    keep <- Reduce(`&`, lapply(sel, function(x) df0[[x]]), init = TRUE) &
      !Reduce(`|`, lapply(other, function(x) df0[[x]]), init = FALSE)
    out  <- df0[keep, c("Gene","GeneName")]
    
    if ("HFD" %in% sel) {
      mi <- hfda_chronic_df %>%
        filter(Gene %in% out$Gene) %>%
        group_by(Gene) %>%
        summarize(miRNA = paste(unique(miRNA), collapse = "; "), .groups="drop")
      out <- left_join(out, mi, by="Gene")
    }
    
    pheno_df <- bind_rows(
      bodycomp_df   %>% mutate(Source = "Bodycomposition"),
      glucose_df    %>% mutate(Source = "Glucose"),
      insulin_df    %>% mutate(Source = "Insulin")
    )
    phenos <- pheno_df %>%
      filter(Gene %in% out$Gene, Source %in% sel) %>%
      group_by(Gene) %>%
      summarize(
        Phenotypes = paste(unique(Phenotype), collapse = "; "),
        .groups="drop"
      )
    
    left_join(out, phenos, by="Gene")
  })
  
  #Download tables
  output$geneTableAcute <- renderDT({
    datatable(
      filteredAcute(),
      rownames   = FALSE,
      extensions = "Buttons",
      options    = list(
        dom     = "Bfrtip",
        buttons = list(
          list(extend   = "copy",  title = NULL),
          list(extend   = "csv",     filename = "HFD_Acute_Overlap"),
          list(extend   = "excel",   filename = "HFD_Acute_Overlap")
        ),
        order      = list(list(0, "asc")),
        pageLength = 5
      )
    )
  }, server = FALSE)
  
  output$geneTableChronic <- renderDT({
    datatable(
      filteredChronic(),
      rownames   = FALSE,
      extensions = "Buttons",
      options    = list(
        dom     = "Bfrtip",
        buttons = list(
          list(extend   = "copy",  title = NULL),
          list(extend   = "csv",     filename = "HFD_Chronic_Overlap"),
          list(extend   = "excel",   filename = "HFD_Chronic_Overlap")
        ),
        order      = list(list(0, "asc")),
        pageLength = 5
      )
    )
  }, server = FALSE)
  
  ##################################
  # In-Silico Tab: Heatmaps & Tables
  ##################################
  
  truncate_label <- function(x, max_chars = 30) {
    vapply(x, function(lbl) {
      if (nchar(lbl) > max_chars) paste0(substr(lbl, 1, max_chars), "…") else lbl
    }, character(1), USE.NAMES = FALSE)
  }
  
  truncate_text <- function(x, n = 30) {
    vapply(x, function(s) {
      if (nchar(s, type="chars") > n) {
        paste0(substr(s, 1, n), "…")
      } else {
        s
      }
    }, character(1), USE.NAMES = FALSE)
  }
  
  membershipISAcute <- reactive({
    allG <- unique(c(genesHFDacute_r(), genesBodycomp, genesGlucose, genesInsulin))
    data.frame(
      Gene            = allG,
      Bodycomposition = allG %in% genesBodycomp,
      Glucose         = allG %in% genesGlucose,
      Insulin         = allG %in% genesInsulin,
      stringsAsFactors = FALSE
    )
  })
  
  membershipISChronic <- reactive({
    allG <- unique(c(genesHFDchronic_r(), genesBodycomp, genesGlucose, genesInsulin))
    data.frame(
      Gene            = allG,
      Bodycomposition = allG %in% genesBodycomp,
      Glucose         = allG %in% genesGlucose,
      Insulin         = allG %in% genesInsulin,
      stringsAsFactors = FALSE
    )
  })
  
  subset_heatmap <- function(hm_df, gene_vec) {
    want <- toupper(gene_vec)
    keep <- sapply(hm_df$Genes, function(gstr) {
      any(toupper(strsplit(gstr, ";")[[1]]) %in% want)
    })
    hm_df[keep, , drop = FALSE]
  }
  
  base_theme <- theme_minimal() +
    theme(
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background  = element_rect(fill = "transparent", color = NA),
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title   = element_blank(),
      plot.title   = element_text(hjust = 0.5)
    )
  
  filteredAcuteIS <- reactive({
    ph <- req(input$inSilico_pheno)
    genes_keep <- membershipISAcute()[[ph]]
    df <- acute_IPA_df %>% filter(Gene %in% membershipISAcute()$Gene[genes_keep])
    if (nzchar(input$acute_search)) df <- df %>% filter(grepl(input$acute_search, Gene, ignore.case = TRUE))
    if (length(input$acute_type) > 0)    df <- df %>% filter(Type %in% input$acute_type)
    df
  })
  
  filteredChronicIS <- reactive({
    ph <- req(input$inSilico_pheno)
    genes_keep <- membershipISChronic()[[ph]]
    df <- chronic_IPA_df %>% filter(Gene %in% membershipISChronic()$Gene[genes_keep])
    if (nzchar(input$chronic_search)) df <- df %>% filter(grepl(input$chronic_search, Gene, ignore.case = TRUE))
    if (length(input$chronic_type) > 0)    df <- df %>% filter(Type %in% input$chronic_type)
    df
  })
  
  annotate_gene <- function(gene, hm_df, name_col) {
    gene_u <- toupper(gene)
    hits <- vapply(hm_df$Genes, function(gstr) {
      genes_in_row <- strsplit(gstr, ";")[[1]]
      genes_in_row <- toupper(trimws(genes_in_row))
      if (gene_u %in% genes_in_row) TRUE else FALSE
    }, logical(1))
    
    found <- hm_df[[name_col]][hits]
    if (length(found) == 0) return("")
    paste(unique(found), collapse = "; ")
  }
  
  output$acute_table <- renderDT({
    df <- filteredAcuteIS() %>%
      select(-Type) %>% 
      mutate(
        miRNA     = truncate_text(miRNA,     30),
        Phenotypes = truncate_text(Phenotypes, 30)
      )
    datatable(
      df,
      rownames   = FALSE,
      extensions = "Buttons",
      options    = list(dom = "ft", order = list(list(0, "asc")), pageLength = 3)
    )
  }, server = FALSE)
  
  output$downloadAcute <- downloadHandler(
      filename = function() paste0("acute_IPA_", input$inSilico_pheno, "_", Sys.Date(), ".csv"),
      content  = function(file) {
        filteredAcuteIS() %>%
          rowwise() %>%
          mutate(
            Pathways  = annotate_gene(Gene, acute_pathways_df,  "Pathway"),
            Functions = annotate_gene(Gene, acute_functions_df, "Function")
          ) %>%
          ungroup() %>%
          write.csv(file, row.names=FALSE)
      }
    )
  
  output$chronic_table <- renderDT({
    df <- filteredChronicIS() %>%
      select(-Type) %>% 
      mutate(
        miRNA     = truncate_text(miRNA,     30),
        Phenotypes = truncate_text(Phenotypes, 30)
      )
    datatable(
      df,
      rownames   = FALSE,
      extensions = "Buttons",
      options    = list(dom = "ft", order = list(list(0, "asc")), pageLength = 3)
    )
  }, server = FALSE)
  
  output$downloadChronic <- downloadHandler(
    filename = function() paste0("chronic_IPA_", input$inSilico_pheno, "_", Sys.Date(), ".csv"),
    content  = function(file) {
      filteredChronicIS() %>%
        rowwise() %>%
        mutate(
          Pathways  = annotate_gene(Gene, chronic_pathways_df,  "Pathway"),
          Functions = annotate_gene(Gene, chronic_functions_df, "Function")
        ) %>%
        ungroup() %>%
        write.csv(file, row.names=FALSE)
    }
  )
  
  output$acute_pathways_heatmap <- renderPlotly({
    df <- acute_pathways_df %>% filter(Phenotype == input$inSilico_pheno) %>%
      subset_heatmap(filteredAcuteIS()$Gene) %>%
      slice_max(`p-value (-log10)`, n = 10, with_ties = FALSE) %>%
      mutate(DisplayPathway = truncate_label(Pathway))
    max_val <- max(df$`p-value (-log10)`, na.rm = TRUE)
    gg <- ggplot(df, aes(x = 1, y = reorder(DisplayPathway, `p-value (-log10)`), fill = `p-value (-log10)`)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "#9CECFB", high = "#0052D4", limits = c(min(df$`p-value (-log10)`), max_val)) +
      labs(title = "Acute Pathways") + base_theme
    ggplotly(gg, tooltip = c("y","fill","Genes"))
  })
  
  output$acute_functions_heatmap <- renderPlotly({
    df <- acute_functions_df %>% filter(Phenotype == input$inSilico_pheno) %>%
      subset_heatmap(filteredAcuteIS()$Gene) %>%
      slice_max(`p-value (-log10)`, n = 10, with_ties = FALSE) %>%
      mutate(DisplayFunction = truncate_label(Function))
    max_val <- max(df$`p-value (-log10)`, na.rm = TRUE)
    gg <- ggplot(df, aes(x = 1, y = reorder(DisplayFunction, `p-value (-log10)`), fill = `p-value (-log10)`)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "#ffc0cb", high = "#800080", limits = c(min(df$`p-value (-log10)`), max_val)) +
      labs(title = "Acute Functions") + base_theme
    ggplotly(gg, tooltip = c("y","fill","Genes"))
  })
  
  output$chronic_pathways_heatmap <- renderPlotly({
    df <- chronic_pathways_df %>% filter(Phenotype == input$inSilico_pheno) %>%
      subset_heatmap(filteredChronicIS()$Gene) %>%
      slice_max(`p-value (-log10)`, n = 10, with_ties = FALSE) %>%
      mutate(DisplayPathway = truncate_label(Pathway))
    max_val <- max(df$`p-value (-log10)`, na.rm = TRUE)
    gg <- ggplot(df, aes(x = 1, y = reorder(DisplayPathway, `p-value (-log10)`), fill = `p-value (-log10)`)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "#9CECFB", high = "#0052D4", limits = c(min(df$`p-value (-log10)`), max_val)) +
      labs(title = "Chronic Pathways") + base_theme
    ggplotly(gg, tooltip = c("y","fill","Genes"))
  })
  
  output$chronic_functions_heatmap <- renderPlotly({
    df <- chronic_functions_df %>% filter(Phenotype == input$inSilico_pheno) %>%
      subset_heatmap(filteredChronicIS()$Gene) %>%
      slice_max(`p-value (-log10)`, n = 10, with_ties = FALSE) %>%
      mutate(DisplayFunction = truncate_label(Function))
    max_val <- max(df$`p-value (-log10)`, na.rm = TRUE)
    gg <- ggplot(df, aes(x = 1, y = reorder(DisplayFunction, `p-value (-log10)`), fill = `p-value (-log10)`)) +
      geom_tile(color = "white") +
      scale_fill_gradient(low = "#FFC0CB", high = "#800080", limits = c(min(df$`p-value (-log10)`), max_val)) +
      labs(title = "Chronic Functions") + base_theme
    ggplotly(gg, tooltip = c("y","fill","Genes"))
  })
  
}

shinyApp(ui, server)
