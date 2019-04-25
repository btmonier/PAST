library(PAST)
library(dplyr)
library(ggplot2)
library(shiny)
library(shinydashboard)
library(gridExtra)

# user interface
ui <- dashboardPage(
  title = "PAST",
  
  # set the title
  dashboardHeader(title = textOutput("title")),
  
  # create the sidebar
  dashboardSidebar(
    
    # create the menu
    sidebarMenu(
      
      # set the analysis title
      textInput("title", "Analysis Title", value = "New Analysis"),
      
      #create a selector for analysis type
      selectInput(
        "type",
        "Analysis Type:",
        choices = c("new", "saved"),
        selected = "new"
      ),
      
      # set the menu based on analysis type (see server)
      menuItemOutput("menuitem"),
      
      # create a plot submenu 
      menuItem(
        "Plot",
        
        # create a selector for filtering type
        selectInput(
          "filter_type",
          "Filter Parameter:",
          choices = c("p-value", "q-value"),
          selected = "p-value"
        ),
        
        # create a slider to specify filtering level
        sliderInput(
          "significance_cutoff",
          "Pathway Signficance Filter",
          min = 0,
          max = 0.3,
          value = 0.05,
          step = 0.01
        )
      )
    )
  ),
  
  # create the body
  shinydashboard::dashboardBody(fluidRow(
    
    # create the table view
    box(
      title = textOutput("box_title_table"),
      width = 4,
      status = "primary",
      solidHeader = TRUE,
      DT::dataTableOutput("pathways"),
      style = "height:82vh; overflow-y: scroll;"
    ),
    
    # create the graph view
    box(
      title = textOutput("box_title_plot"),
      status = "primary",
      color = "red",
      solidHeader = TRUE,
      width = 8,
      plotOutput("plots", height = "auto"),
      style = "height:82vh; overflow-y: scroll;"
    ),
    
    # add the download button
    downloadButton("download_data", 
                   "Download Results", 
                   style = "margin-left: 15px")
  ))
)

# server functions
server <- function(input, output) {
  options(shiny.maxRequestSize = 600 * 1024 ^ 2)
  
  # set the title
  output$title <- renderText({
    input$title
  })
  
  # build the menu based on new/saved
  output$menuitem <- renderMenu({
    
    # new analysis menu
    if (input$type == "new")
      return(
        menuItem(
          # title
          "Parameters",
          tabName =  "live",
          
          # file inputs
          fileInput("association_file", "Association File"),
          menuItem("Association Column Names",
                   textInput("association_trait", 
                             "Trait Column Name", 
                             value = "Trait"),
                   textInput("association_marker", 
                             "Marker Column Name", 
                             value = "Marker"),
                   textInput("association_locus", 
                             "Locus Column Name", 
                             value = "Locus"),
                   textInput("association_site", 
                             "Site Column Name", 
                             value = "Site"),
                   textInput("association_p", 
                             "P-value Column Name", 
                             value = "p"),
                   textInput("association_marker_R2", 
                             "Marker R^2 Column Name", 
                             value = "marker_R2")),
          fileInput("effects_file", "Effects File"),
          menuItem("Effects Column Names",
                   textInput("effects_trait", 
                             "Trait Column Name", 
                             value = "Trait"),
                   textInput("effects_marker", 
                             "Marker Column Name", 
                             value = "Marker"),
                   textInput("effects_locus", 
                             "Locus Column Name", 
                             value = "Locus"),
                   textInput("effects_site", 
                             "Site Column Name", 
                             value = "Site"),
                   textInput("effects_effect", 
                             "Effect Column Name", 
                             value = "Effect")),
          fileInput("LD_file", "Linkage Disequilibrium File"),
          menuItem("LD Column Names",
                   textInput("LD_locus_1", 
                             "Locus1 Column Name", 
                             value = "Locus1"),
                   textInput("LD_position_1", 
                             "Position1 Column Name", 
                             value = "Position1"),
                   textInput("LD_site_1", 
                             "Locus Column Name", 
                             value = "Site1"),
                   textInput("LD_position_2", 
                             "Site Column Name", 
                             value = "Position2"),
                   textInput("LD_site_2", 
                             "Site Column Name", 
                             value = "Site2"),
                   textInput("LD_distance", 
                             "Site Column Name", 
                             value = "Dist_bp"),
                   textInput("LD_R2", 
                             "Effect Column Name", 
                             value = "R.2")),
          fileInput("genes_file", "Genes File"),
          fileInput("pathway_file", "Pathways File"),
          
          # number of cores
          numericInput("num_cores", "Number of Cores", value = 2),
          
          # analysis type
          selectInput("mode", "Mode:",
                      choices = c("increasing", "decreasing")),
          
          # advanced options
          menuItem(
            "Advanced",
            
            # window size to search for genes
            numericInput("window", "Window Size", value = 1000),
            
            # r^2 for LD
            sliderInput(
              "r_squared_cutoff",
              "R^2 cutoff:",
              min = 0,
              max = 1,
              value = 0.8,
              step = 0.05
            ),
            
            # minimum number of genes in a pathway
            numericInput("gene_cutoff", "Gene Cutoff", value = 5),
            
            # number of times to sample effects for generating null distribution
            numericInput("sample", "Effects", value = 1000)
          )
        )
      )
    # saved analysis menu
    else
      (return (menuItem("Saved", fileInput("load_file", "Data"))))
  })
  
  # download handler
  output$download_data <- downloadHandler(
    
    # make filename = analysis_title.zip
    filename = function() {
      paste(input$title, ".zip", sep = "")
    },
    # set up content
    content = function(filename) {
      # get filter type and cutoff
      filter_type <- input$filter_type
      significance_cutoff <- input$significance_cutoff
      
      # set up empty vector of files and move to tempdir()
      fs <- c()
      setwd(tempdir())
      
      # get pathway significance data based on analysis type
      # R Shiny will either calculate this if something has changed or use 
      # what's already on the screen in most cases
      if (input$type == "new") {
        pathways <- find_pathways()
      } else {
        pathways <- load_pathways_from_file()
      }
      
      # write full pathways file
      write.table(pathways, paste0(input$title, ".pathways.tsv"), sep = "\t")
      
      # filter pathways
      if (filter_type == "p-value") {
        filtered_pathways <- dplyr::filter(pathways, 
                                           pvalue < significance_cutoff)
      } else {
        filtered_pathways <- dplyr::filter(pathways, 
                                           qvalue < significance_cutoff)
      }
      
      # write filtered pathways file
      write.table(
        filtered_pathways,
        paste0(input$title, ".pathways.filtered.tsv"),
        sep = "\t"
      )
      
      # store pathways files
      fs <-
        c(
          paste0(input$title, ".pathways.tsv"),
          paste0(input$title, ".pathways.filtered.tsv")
        )
      
      # break until we have pathways information
      if (is.null(pathways))
        return(NULL)
      
      # set up rugplots by arranged by pathway_number and filtering
      rugplots_data <- pathways %>% dplyr::arrange(.data$pathway_number)
      if (filter_type == "p-value") {
        rugplots_data <- dplyr::filter(rugplots_data, 
                                       pvalue < significance_cutoff)
      } else {
        rugplots_data <- dplyr::filter(rugplots_data, 
                                       qvalue < significance_cutoff)
      }
      
      # split based on pathway number
      rugplots_split <- split(rugplots_data, rugplots_data$pathway_number)
      
      # for each pathway, draw rugplot
      for (rank in names(rugplots_split)) {
        
        # get data
        temp_data <- rugplots_split[[rank]]
        
        # title is "PWY-ID - Pathway Description"
        title <- paste0(unique(as.character(temp_data$pathway_id)),
                        " - ",
                        unique(as.character(temp_data$pathway_name)))
        
        # intercept should be at rank of maximum enrichment score
        intercept <-
          temp_data %>% dplyr::arrange(desc(.data$running_enrichment_score)) %>%
          dplyr::select(rank)
        intercept <- intercept[, 1][1]
        
        # set up the rugplot
        rugplot <-
          ggplot(temp_data, aes(x = rank, y = running_enrichment_score)) +
          geom_line(stat = "identity") +
          geom_rug(sides = "t", position = "jitter") +
          geom_vline(xintercept = intercept,
                     color = "black",
                     linetype = "longdash") +
          ggtitle(title) +
          labs(x = "Gene Rank", y = "Running Enrichment Score") +
          scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000)) +
          theme(axis.text = element_text (color = "black"),
                panel.background = element_rect (color = "black", fill = "pink")
          )
        
        # set up the output path for the rugplot
        path <- paste0(input$title,
                       ".",
                       unique(as.character(temp_data$pathway_id)),
                       ".png")
        
        # add the rugplot to the files to be zipped and save it in tempdir()
        fs <- c(fs, path)
        ggsave(path, rugplot)
      }
      
      # zip the file using the name and files specified earlier
      zip(zipfile = filename, files = fs)
    },
    contentType = "application/zip"
  )
  
  # set the table box title based on whether we are doing a new analysis
  # or loading saved data
  output$box_title_table <- renderText({
    if (input$type == "new") {
      return(paste0(input$mode, " pathways"))
    }
    else {
      return(paste0(input$load_file$name, " pathways"))
    }
  })
  
  # set the plot box title based on whether we are doing a new analysis
  # or loading saved data
  output$box_title_plot <- renderText({
    if (input$type == "new") {
      return(paste0(input$mode, " plots"))
    }
    else {
      return(paste0(input$load_file$name, " plots"))
    }
  })
  
  # get GWAS data
  get_gwas_data <- reactive({
    association_file <- input$association_file
    effects_file <- input$effects_file
    association_columns = c(input$association_trait,
                            input$association_marker,
                            input$association_locus,
                            input$association_site,
                            input$association_p,
                            input$association_marker_R2)
    
    effects_columns = c(input$effects_trait,
                        input$effects_marker,
                        input$effects_locus,
                        input$effects_site,
                        input$effects_effect)
    
    if (is.null(association_file) | is.null(effects_file))
      return(NULL)
    load_GWAS_data(association_file$datapath, 
                   effects_file$datapath,
                   association_columns,
                   effects_columns)
  })
  
  # load the LD data
  get_LD <- reactive({
    LD_file <- input$LD_file
    if (is.null(LD_file))
      return(NULL)
    LD_columns = c(input$LD_locus_1,
                   input$LD_position_1,
                   input$LD_site_1,
                   input$LD_position_2,
                   input$LD_site_2,
                   input$LD_distance,
                   input$LD_R2)
    load_LD(LD_file$datapath)
  })
  
  # assign SNPs to genes
  assign_genes <- reactive({
    all_data <- get_gwas_data()
    LD <- get_LD()
    genes_file <- input$genes_file
    if (is.null(genes_file) | is.null(LD))
      return(NULL)
    window <- input$window
    r_squared_cutoff <- input$r_squared_cutoff
    num_cores <- input$num_cores
    assign_SNPs_to_genes(all_data,
                         LD,
                         genes_file$datapath,
                         window,
                         r_squared_cutoff,
                         num_cores)
  })
  
  # find pathway significance
  find_pathways <- reactive({
    genes <- assign_genes()
    pathway_file <- input$pathway_file
    if (is.null(pathway_file) | is.null(genes))
      return(NULL)
    gene_cutoff <- input$gene_cutoff
    mode <- input$mode
    sample <- input$sample
    num_cores <- input$num_cores
    find_pathway_significance(genes,
                              pathway_file$datapath,
                              gene_cutoff,
                              mode,
                              sample,
                              num_cores)
  })
  
  # load pathway sifnificance data from file
  load_pathways_from_file <- reactive({
    pathway_file <- input$load_file
    if (is.null(pathway_file))
      return(NULL)
    read.table(pathway_file$datapath, header = TRUE, sep = "\t")
  })
  
  # render pathway table
  output$pathways <- DT::renderDataTable({
    
    # determine source of data based on analysis type
    if (input$type == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    
    # get filter and cut-off
    filter_type <- input$filter_type
    significance_cutoff <- input$significance_cutoff
    if (is.null(pathways))
      return(NULL)
    
    # get the pathway, p-value, and q-value
    pathways <- pathways %>%
      dplyr::select(.data$pathway_id, .data$pvalue, .data$qvalue) %>%
      unique()
    
    # filter pathways
    if (filter_type == "p-value") {
      pathways <- dplyr::filter(pathways, pvalue < significance_cutoff)
    } else {
      pathways <- dplyr::filter(pathways, qvalue < significance_cutoff)
    }
    
    # render data table without paging and no rownames
    DT::datatable(pathways,
                  rownames = FALSE,
                  options = list(paging = FALSE))
  })
  
  # render plots
  output$plots <- renderPlot({
    
    # determine source of data based on analysis type
    if (input$type == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    
    # get filter and cut-off
    filter_type <- input$filter_type
    significance_cutoff <- input$significance_cutoff
    if (is.null(pathways))
      return(NULL)
    
    # filter pathways
    rugplots_data <- pathways %>% dplyr::arrange(.data$pathway_number)
    if (filter_type == "p-value") {
      rugplots_data <- dplyr::filter(rugplots_data, 
                                     pvalue < significance_cutoff)
    } else {
      rugplots_data <- dplyr::filter(rugplots_data,
                                     qvalue < significance_cutoff)
    }
    
    # split rugplots data into a list of pathways
    rugplots_split <- split(rugplots_data, rugplots_data$pathway_number)
    
    # set up a list to store plot instructions
    plots_list <- list()
    
    # for each pathway, draw rugplot
    for (rank in names(rugplots_split)) {
      
      # get data
      temp_data <- rugplots_split[[rank]]
      
      # title is "PWY-ID - Pathway Description"
      title <- paste0(unique(as.character(temp_data$pathway_id)),
                      " - ",
                      unique(as.character(temp_data$pathway_name)))
      
      # intercept should be at rank of maximum enrichment score
      intercept <- temp_data %>%
        dplyr::arrange(desc(.data$running_enrichment_score)) %>%
        dplyr::select(.data$rank)
      intercept <- intercept[, 1][1]
      
      # set up the rugplot
      rugplot <-
        ggplot(temp_data, aes(x = rank, y = running_enrichment_score)) +
        geom_line(stat = "identity") +
        geom_rug(sides = "t", position = "jitter") +
        geom_vline(xintercept = intercept,
                   color = "black",
                   linetype = "longdash") +
        ggtitle(title) +
        labs(x = "Gene Rank", y = "Running Enrichment Score") +
        scale_x_continuous(breaks = c(0, 5000, 10000, 15000, 20000, 25000)) +
        theme(axis.text = element_text (color = "black"),
              panel.background = element_rect (color = "black", fill = "pink")
        )
      
      # store the rugplot in the list
      plots_list[[rank]] <- rugplot
    }
    
    # draw a blank if no rugplots were plotted, otherwise, draw rugplots in grid
    if (nrow(rugplots_data) == 0) {
      ggplot() + geom_blank()
    } else {
      columns <- floor(sqrt(length(plots_list)))
      do.call("grid.arrange", c(plots_list, ncol = columns))
    }
  } # complicated way to set the height based on the number of columns
  , height = reactive({
    
    # get number of pathways being plotted
    if (input$type == "new") {
      pathways <- find_pathways()
    } else {
      pathways <- load_pathways_from_file()
    }
    
    # filter
    filter_type <- input$filter_type
    significance_cutoff <- input$significance_cutoff
    
    # return 100 if there are no pathways to keep the function from breaking
    if (is.null(pathways))
      return(100)

    # get rugplots data for column calculation
    rugplots_data <- pathways %>% dplyr::arrange(.data$pathway_number)
    if (filter_type == "p-value") {
      rugplots_data <- dplyr::filter(rugplots_data, 
                                     pvalue < significance_cutoff)
    } else {
      rugplots_data <- dplyr::filter(rugplots_data,
                                     qvalue < significance_cutoff)
    }
    
    # return 100 if there are no pathways to keep the function from breaking  
    if (nrow(rugplots_data) == 0)
      return(100)
    
    # split rugplots data into a list of pathways
    rugplots_split <- split(rugplots_data, rugplots_data$pathway_number)
    
    # get number of columns
    columns <- floor(sqrt(length(names(rugplots_split))))

    # return height
    return(length(names(rugplots_split)) / columns * 200)
  }))
}

# run the app
shinyApp(ui = ui, server = server)