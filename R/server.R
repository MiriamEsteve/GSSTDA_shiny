#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(grid)
# BiocManager::install("InteractiveComplexHeatmap", force=T)
library(InteractiveComplexHeatmap)
library(shinyjs)
library(visNetwork)

source("func.R")


options(shiny.maxRequestSize=30*1024^2)
# Define server logic required to draw a histogram
function(input, output, session) {
  ##   ##   ##   ##   ##   ##   ## Data Descriptive Panel   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##

  ## input$file is a data frame and contains the details around the name, size and temp location of the files uploaded
  # this reactive output display the content of the input$file dataframe
  output$filedf <- renderTable({
    if(is.null(input$file)){return ()}
    # the file1 input data frame object that contains the file attributes
    subset(input$file, select = c("name", "size", "type"))
    })


  ## Below code to display the structure of the input file object
  output$fileob <- renderPrint({
    if(is.null(input$file)){return ()}
    str(input$file)
  })

  ## Side bar select input widget coming through renderUI()
  # Following code displays the select input widget with the list of file loaded by the user
  output$selectfile <- renderUI({
    if(is.null(input$file)) {return()}
    fluidRow( list(hr(),
         helpText("Select the files for which you need to see data and summary stats"),
         selectInput("Select", "Select", choices=input$file$name)
    ))
  })

  ## Summary Stats code ##
  # this reactive output contains the summary of the dataset and display the summary in table format
  output$summ <- renderPrint({
    if(is.null(input$file)){return()}
    if(input$rowNames == TRUE) row_names <- 1 else row_names <- NULL
    summary(read.table(file=input$file$datapath[input$file$name==input$Select],
                       sep=input$sep,
                       header = input$header,
                       stringsAsFactors = input$stringAsFactors, row.names=row_names))})

  ## Dataset code ##
  # This reactive output contains the dataset and display the dataset in table format
  output$table <- renderTable({
    if(is.null(input$file)){return()}
    if(input$rowNames == TRUE) row_names <- 1 else row_names <- NULL

    if(input$Select == "full_data.csv"){
      data <- head(read.table(file=input$file$datapath[input$file$name==input$Select], sep=input$sep, header = input$header, stringsAsFactors = input$stringAsFactors, row.names=row_names))
      data <- select(data, colnames(data)[1:10])
    }
    else{
      data <- read.table(file=input$file$datapath[input$file$name==input$Select], sep=input$sep, header = input$header, stringsAsFactors = input$stringAsFactors, row.names=row_names)
      nrows <- trunc(nrow(data)/2)
      data <- head(data, nrows)
    }
    data
  })

  output$dataset_info <- renderPrint({
    if(input$rowNames == TRUE) row_names <- 1 else row_names <- NULL

    if(input$Select == "full_data.csv")
      cat("It shows the 10 first columns and rows of the data set")
    else{
      data <- read.table(file=input$file$datapath[input$file$name==input$Select], sep=input$sep, header = input$header, stringsAsFactors = input$stringAsFactors, row.names=row_names)
      nrows <- trunc(nrow(data)/2)
      cat("It shows the", nrows , " first rows of the data set")
    }

  })

  ## MainPanel tabset renderUI code ##
  # the following renderUI is used to dynamically generate the tabsets when the file is loaded.
  # Until the file is loaded, app will not show the tabset.
  output$tb <- renderUI({
    if(is.null(input$file)) {return()}
    else
      tabsetPanel(
        tabPanel("Dataset", textOutput("dataset_info"), tableOutput("table")),
        tabPanel("Summary Stats", verbatimTextOutput("summ")))
  })

  ## ERRORS MESSAGES
  myText <- "ERROR. Load data files."
  myText.color="red"

  ##   ##   ##   ##   ##   ##  DGSA panel ##  ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##

  ## Side bar select input widget coming through renderUI()
  output$select_control_tag <- renderUI({
    # add validate here
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
    {
      data <- read_files(input)
      control_tag_opt <- unique(data$case_tag)
      selectInput("select_control_tag", paste("What is the tag of the healthy patient (value in the case_tag)? (",
                                              control_tag_opt[[1]], " or ", control_tag_opt[[2]], "): " , sep=""), choices=control_tag_opt,
                  selected=control_tag_opt[[1]])
    }
})

  ## DGSA panel ##
  output$dgsa_text <- renderPrint({
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
      cat("Click button to execute DGSA process")
  })

  dgsa_button <- observeEvent(input$dgsa_button, {
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
    {
      data <- read_files(input)
      output$dgsa_text <- renderPrint({cat("\nThe pre-process DGSA is started")})

      DGSA_object <- DGSA_app(data$full_data, data$survival_time, data$survival_event, data$case_tag,
                               input$select_yes_no, input$select_control_tag)
      DGSA_info <- results_DGSA_app(DGSA_object[["matrix_disease_component"]])
      output$dgsa_text <- renderPrint({cat("\nThe pre-process DGSA is finished\n")})
      showNotification("\nThe pre-process DGSA is finished\n")


      # Plot heatmaps
      col_fun = circlize::colorRamp2(c(-4, 0,4),
                                     c("red", "black", "green"))
      ha = HeatmapAnnotation(Group = case_tag,
                             annotation_legend_param = list(
                               Group = list(title = "Group", at = c(unique(case_tag)))
                             ))
      ht1 = Heatmap(DGSA_info,
                    cluster_columns = T,
                    col = col_fun,
                    cluster_rows = F,
                    heatmap_legend_param = list(
                      title = "Expression",
                      at = c(-4, 0, 4)),
                    top_annotation = ha)

      ht1 = draw(ht1)

      makeInteractiveComplexHeatmap(input, output, session, ht1, "heatmap_1")
      dev.off()

    }

})

  ##   ##   ##   ##   ##   ##   geneSelection panel ##  ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##

  ## Side bar select input widget coming through renderUI()
  output$select_control_tag2 <- renderUI({
    # add validate here
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
    {
      data <- read_files(input)
      control_tag_opt <- unique(data$case_tag)
      selectInput("select_control_tag2", paste("What is the tag of the healthy patient (value in the case_tag)? (",
                                              control_tag_opt[[1]], " or ", control_tag_opt[[2]], "): " , sep=""), choices=control_tag_opt,
                  selected=control_tag_opt[[1]])
    }
  })

  output$geneSelection_text <- renderPrint({
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
      cat("Click button to execute gene selection process")
  })

  ## geneSelection panel ##
  output$cox_matrix <- renderTable({
    as.table(geneSelection_object[["cox_all_matrix"]])
  })
  output$gene_select <- renderTable({
    as.vector(geneSelection_object[["genes_selected"]])
  })
  output$genes_disease <- renderTable({
    as.table(geneSelection_object[["genes_disease_component"]])
  })
  output$filter <- renderTable({
    as.table(geneSelection_object[["filter_values"]])
  })

  output$gene <- renderUI({
    # add validate here
    if(input$select_dgsa == "no"){
      shinyjs::hide("gen_select_type2")
      shinyjs::hide("percent_gen_select2")
    }
    else
    {
      fluidRow(
        selectInput("gen_select_type2", "Select gene type", c("Abs", "Top_Bot")),
        numericInput("percent_gen_select2", "Introduce percent of gene to be selected", value = 10, min=10, max=100)
      )
    }
  })

  geneSelection_button <- observeEvent(input$geneSelection_button, {
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
    {
      data <- read_files(input)

      if(input$select_dgsa == "yes"){
        DGSA_object <- DGSA_app(data$full_data, data$survival_time, data$survival_event, data$case_tag,
                                input$select_yes_no2, input$select_control_tag2)

        # Create a Progress object
        progress <- shiny::Progress$new()

        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        progress$set(message = "\nCalculating the matrix of Zcox from DGSA data.", value = 0)

        geneSelection_object <<- geneSelection_app(progress, DGSA_object, input$gen_select_type2, input$percent_gen_select2, input$select_yes_no2, input$select_control_tag2)

        output$geneSelection_text <- renderPrint({cat("\nThe gene selection process from DGSA data is finished\n")})
        showNotification("\nThe gene selection process from DGSA data is finished\n")

      }
      else{
        # Create a Progress object
        progress <- shiny::Progress$new()
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        progress$set(message = "Calculating the matrix of Zcox from files uploaded.", value = 0)


        #Select gene from data object
        geneSelection_object <<- geneSelection_app(progress, data, input$gen_select_type2, input$percent_gen_select2, input$select_yes_no2, input$select_control_tag2)

        output$geneSelection_text <- renderPrint({cat("\nThe gene selection process from files uploaded is finished\n")})
        showNotification("\nThe gene selection process from files uploaded is finished\n")

      }

      # Output panel
      output$geneSelection_tb <- renderUI({
        tabsetPanel(
          tabPanel("cox_all_matrix", tableOutput("cox_matrix")),
          tabPanel("genes_selected", tableOutput("gene_select"))
          #tabPanel("genes_disease_component", tableOutput("genes_disease")),
          #tabPanel("filter_values", tableOutput("filter"))
        )
      })

    }

  })

  ##   ##   ##   ##   ##   ##   Mapper panel ##  ##   ##   ##   ##   ##   ##   ##   ##   ##   ##   ##
  ## Side bar select input widget coming through renderUI()
  output$select_control_tag3 <- renderUI({
    # add validate here
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
    {
      data <- read_files(input)
      control_tag_opt <- unique(data$case_tag)
      selectInput("select_control_tag3", paste("What is the tag of the healthy patient (value in the case_tag)? (",
                                               control_tag_opt[[1]], " or ", control_tag_opt[[2]], "): " , sep=""), choices=control_tag_opt,
                  selected=control_tag_opt[[1]])
    }
  })

  output$mapper_text <- renderPrint({
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
      cat("Click button to execute mapper process")
  })

  output$gene_mapper <- renderUI({
    # add validate here
    if(input$select_dgsa_gene == "no"){
      shinyjs::hide("gen_select_type3")
      shinyjs::hide("percent_gen_select3")
    }
    else
    {
      fluidRow(
         selectInput("gen_select_type3", "Select gene type", c("Abs", "Top_Bot")),
         numericInput("percent_gen_select3", "Introduce percent of gene to be selected", value = 10, min=10, max=100)
      )
    }
  })

 # Output mapper

  output$plot_mapper <- renderVisNetwork({
    plot_mapper(mapper_object)
  })

  #output$interval_data <- renderPrint({
  #  print(mapper_object[["interval_data"]])
  #})
  #output$sample_in_level <- renderPrint({
  #  print(mapper_object[["sample_in_level"]])
  #})
  #output$interval_data <- renderPrint({
  #  print(mapper_object[["interval_data"]])
  #})
  #output$clustering_all_levels <- renderPrint({
  #  print(mapper_object[["clustering_all_levels"]])
  #})

  mapper_button <- observeEvent(input$mapper_button, {
    if(is.null(input$file)){
      validate(
        need(is.null(data), myText) # display custom message in need
      )
      showNotification(myText)
    }
    else
    {
      data <- read_files(input)

      if(input$select_dgsa_gene == "yes"){

        DGSA_object <- DGSA_app(data$full_data, data$survival_time, data$survival_event, data$case_tag,
                                input$select_yes_no3, input$select_control_tag3)

        # Create a Progress object
        progress <- shiny::Progress$new()
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())
        progress$set(message = "\nCalculating the matrix of Zcox from DGSA data.", value = 0)

        geneSelection_object <- geneSelection_app(progress, DGSA_object, input$gen_select_type3, input$percent_gen_select3, input$select_yes_no3, input$select_control_tag3)

        showNotification("\nThe gene selection process from DGSA data is finished\n")


        mapper_object <<- mapper_app(geneSelection_object[["genes_disease_component"]],
                                    geneSelection_object[["filter_values"]],
                                    input$num_intervals,
                                    input$percent_overlap,
                                    input$distance_type,
                                    input$clustering_type,
                                    input$num_bins_when_clustering,
                                    input$linkage_type,
                                    input$optimal_clustering_mode,
                                    input$select_yes_no3,
                                    "checked")

        output$mapper_text <- renderPrint({cat("\nThe mapper process from DGSA data and gene selection is finished\n")})
        showNotification("\nThe mapper process from DGSA data and gene selection is finished\n")
      }
      else{
        mapper_object <<- mapper_app(geneSelection_object[["genes_disease_component"]],
                                     geneSelection_object[["filter_values"]],
                                     input$num_intervals,
                                     input$percent_overlap,
                                     input$distance_type,
                                     input$clustering_type,
                                     input$num_bins_when_clustering,
                                     input$linkage_type,
                                     input$optimal_clustering_mode,
                                     input$select_yes_no3)

        showNotification("\nThe mapper process from files uploaded is finished\n")
      }

      # Output panel
      output$mapper_tb <- renderUI({
        tabsetPanel(
          tabPanel("plot_mapper", visNetworkOutput("plot_mapper"))
          #tabPanel("interval_data", textOutput("interval_data")),
          #tabPanel("sample_in_level",textOutput("sample_in_level")),
          #tabPanel("clustering_all_levels", textOutput("clustering_all_levels"))
        )
      })

    }

  })


}
