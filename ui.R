library(shiny)
library(shinythemes)
library(InteractiveComplexHeatmap)
library(GSSTDA)

# Define UI ----
ui <- fluidPage(theme = shinytheme("journal"),
                titlePanel("GSSTDA: Gene Structure Survival using Topological Data Analysis"),
                sidebarPanel(
                  fileInput("file","Upload the 4 CSV file:
                            full_data.csv, survival_time.csv, survival_event.csv and case_tag.csv",
                            multiple = TRUE,  accept = ".csv"), # fileinput() function is used to get the file upload control option
                  helpText("Default max. file size is 5MB"),
                  helpText("Select the read.table parameters below"),
                  checkboxInput(inputId = 'header', label = 'Header', value = TRUE),
                  checkboxInput(inputId = "stringAsFactors", "stringAsFactors", FALSE),
                  checkboxInput(inputId = "rowNames", "row.names", TRUE),
                  radioButtons(inputId = 'sep', label = 'Separator', choices = c(Comma=',',Semicolon=';',Tab='\t', Space=''), selected = ','),
                  uiOutput("selectfile"),
                  fluidRow("File Objects information", tableOutput("filedf"))

                ),
                mainPanel(
                  tabsetPanel(
                    tabPanel("Data",
                             fluidRow(h3("Data descriptive")),
                             uiOutput("tb")

                    ),
                    tabPanel("DGSA",
                             fluidRow(
                               h3("First step of the process: DGSA"),
                               p("This analysis, developed by Nicolau et al. (2007) is independent of the rest of the process and can be used with the data for further analysis other than mapper.
                 It allows the calculation of the 'disease component' which consists of, through linear models, eliminating the part of the data
                 that is considered normal or healthy and keeping only the component that is due to the disease.")
                             ),
                             fluidRow(sidebarPanel(
                               fluidRow(selectInput("select_yes_no", "Are the columns of the data set the patient and the rows the genes?", c("yes", "no"), selected = "yes")),
                               fluidRow(uiOutput("select_control_tag"))
                             )),
                             fluidRow(actionButton("dgsa_button", "Update Data Model")),
                             fluidRow(textOutput("dgsa_text")),
                             fluidRow(InteractiveComplexHeatmapOutput("heatmap_1"))

                    ),
                    tabPanel("Gene selection",
                             fluidRow(
                               h3("Second step of the process: Select the genes within the DGSA object created in the previous step and calcute the values of the filtering functions."),
                               p("After performing a survival analysis of each gene, this function selects the genes to be used in the mapper according to both their variability within the database and their relationship with survival. Subsequently, with the genes selected, the values of the filtering functions are calculated for each patient. The filter function allows to summarise each vector of each individual in a single data. This function takes into account the survival associated with each gene.")
                             ),
                             fluidRow(sidebarPanel(
                              fluidRow(selectInput("select_yes_no2", "Are the columns of the data set the patient and the rows the genes?", c("yes", "no"), selected = "yes")),
                              fluidRow(uiOutput("select_control_tag2")),
                              fluidRow(selectInput("select_dgsa", "Do you want select the genes using DGSA process?", c("yes", "no"), selected = "yes")),
                              fluidRow(selectInput("gen_select_type", "Select gene type", c("Abs", "Top_Bot"))),
                              fluidRow(numericInput("percent_gen_select", "Introduce percent of gene to be selected", value = 10, min=10, max=100))
                             )),
                             fluidRow(actionButton("geneSelection_button", "Update Data Model")),
                             fluidRow(textOutput("geneSelection_text")),
                             uiOutput("geneSelection_tb")
                    ),
                    tabPanel("Mapper",
                             fluidRow(
                               h3("Third step of the process: Create the mapper object with disease component matrix with only the selected genes and the filter function obtained in the gene selection step."),
                               p("Mapper condenses the information of high-dimensional datasets into a combinatory graph that is referred to as the skeleton of the dataset. To do so, it divides the dataset into different levels according to its value of the filtering function. These levels overlap each other. Within each level, an independent clustering is performed using the input matrix and the indicated distance type. Subsequently, clusters from different levels that share patients with each other are joined by a vertex.
                                  This function is independent from the rest and could be used without having done DGSA and gene selection")

                              ),
                             sidebarPanel(
                               fluidRow(selectInput("select_yes_no3", "Are the columns of the data set the patient and the rows the genes?", c("yes", "no"), selected = "yes")),
                               fluidRow(uiOutput("select_control_tag3")),
                               fluidRow(selectInput("select_dgsa_gene", "Do you want select the genes using DGSA and gene selection process?", c("yes", "no"), selected = "yes")),
                               fluidRow(uiOutput("gene_mapper")),
                             ),

                             fluidRow(sidebarPanel(
                               fluidRow(numericInput("num_intervals", "Introduce num intervals to be selected", value = 5, min=10, max=100)),
                               fluidRow(numericInput("percent_overlap", "Introduce percent of overlap to be selected", value = 40, min=10, max=100)),
                               fluidRow(selectInput("distance_type", "Select distance type", c("cor","euclidean"))),
                               fluidRow(selectInput("clustering_type", "Select clustering type", c("hierarchical","PAM"))),
                               fluidRow(selectInput("optimal_clustering_mode", "Select optimal cluster number method", c("standard","silhouette"))),
                               fluidRow(selectInput("linkage_type", "Select linkage type", c("single","average","complete")))
                             )),
                             fluidRow(actionButton("mapper_button", "Update Data Model")),
                             fluidRow(textOutput("mapper_text")),
                             uiOutput("mapper_tb")
                    )
                  )
                )
)
