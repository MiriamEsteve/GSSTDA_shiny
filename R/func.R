library(GSSTDA)
library(dplyr)
library(shiny)
library(shinyalert)

source("One_D_Mapper_app.R")

# Function to read files uploaded
read_files <- function(input){

  if(input$rowNames == TRUE) row_names <- 1 else row_names <- NULL

  data <- list("full_data" = NULL,
            "survival_event" = NULL,
            "survival_time" = NULL,
            "case_tag" = NULL
  )
  for (i in 1:length(input$file$name)){

    if(input$file$name[i] == "survival_event.csv"){
      data$survival_event <- as.vector(unlist(read.table(file=input$file$datapath[input$file$name=="survival_event.csv"],
                                           sep=input$sep, header = input$header, stringsAsFactors = input$stringAsFactors,
                                        row.names=row_names), use.names = FALSE))
    }
    if(input$file$name[i] == "survival_time.csv"){
      data$survival_time <-  as.vector(unlist(read.table(file=input$file$datapath[input$file$name=="survival_time.csv"], sep=input$sep,
                                          header = input$header, stringsAsFactors = input$stringAsFactors,
                                          row.names=row_names), use.names = FALSE))
    }
    if(input$file$name[i] == "case_tag.csv"){
      data$case_tag <-  as.vector(unlist(read.table(file=input$file$datapath[input$file$name=="case_tag.csv"], sep=input$sep,
                                     header = input$header, stringsAsFactors = input$stringAsFactors,
                                     row.names=row_names), use.names = FALSE))
    }
    if(input$file$name[i] == "full_data.csv"){
      data$full_data <- as.matrix(read.table(file=input$file$datapath[input$file$name=="full_data.csv"], sep=input$sep,
                                      header = input$header, stringsAsFactors = input$stringAsFactors, row.names=row_names))
    }

  }
  return(data)
}

## Function to check full_data with col_row_check
check_full_data_app <- function(full_data, col_row_check,  na.rm = TRUE){
  # If this function has been executed don't do nothing
  if(na.rm == "checked"){
    return(full_data)
  }

  #Read the data set
  yes_no <- col_row_check

  if(yes_no == "no" | yes_no == "n"){
    #Transpose the data set. Columns = patient and rows = genes
    full_data <- t(full_data)
  }
  #Convert full_data to matrix type
  full_data <- as.matrix(full_data)

  #Omit NAN's values
  if (na.rm == TRUE){
    nrow_ini = nrow(full_data)
    # Remove rows (genes) with NA's values
    full_data <- full_data[rowSums(is.na(full_data))==0,]

    # Show a simple modal
    msg <- paste(nrow(full_data) - nrow_ini, " missing values and NaN's are omitted in the genes (rows)")
    showNotification(msg, duration = NULL)
  }
  return(full_data)
}

## Function to check vectors with col_row_check
check_vectors_app <- function(full_data, survival_time, survival_event, case_tag, control_tag, na.rm = TRUE){
  ncol_full_data <- ncol(full_data)

  # Check if the arguments are vectors; a valid type of data; and the vectors are the same dimension as a full_data
  if(!is.vector(survival_time) | !is.numeric(survival_time) | length(survival_time) != ncol_full_data){
    msg <- "survival_time must be a valid values vector and its length must be the same as the number of patients (columns) of the full_data."
    shinyalert::shinyalert(title = msg, type = "error")
  }

  # Omit NAN's values in checking
  if(!is.vector(survival_event) | !(length(unique(stats::na.omit(survival_event))) == 2 & is.numeric(stats::na.omit(survival_event))) | length(survival_event) != ncol_full_data){
    msg <- "survival_event must be a valid values vector. Only two type of event (0 or 1). Also, its length must be the same as the number of patients (columns) of the full_data."
    shinyalert::shinyalert(title = msg, type = "error")
  }

  #If exits NAN's values remove it and check if it contain only two cases and it has the same dimension (columns) as full_data
  if(!is.vector(case_tag)){
    msg <- "case_tag must be a valid values vector."
    shinyalert::shinyalert(title = msg, type = "error")
  }
  if( any(is.na(case_tag))){
    without_nan_patient <- which(!is.na(case_tag))
    case_tag <- case_tag[without_nan_patient]
    full_data <- full_data[,without_nan_patient]
    survival_event <- survival_event[without_nan_patient]
    survival_time <- survival_time[without_nan_patient]
    ncol_full_data <- ncol(full_data)
    msg <- "NAN's values in patient was removed in case_tag, full_data, survival_time and survival_event"
    showNotification(msg, duration = NULL)

  }

  if(length(unique(case_tag)) != 2){
    msg <- "case_tag must has only two type of tags."
    shinyalert::shinyalert(title = msg, type = "error")

  }

  if(length(case_tag) != ncol_full_data){
    msg <- "The length of case_tag must be the same as the number of patients (columns) of the full_data."
    shinyalert::shinyalert(title = msg, type = "error")
}


  return(list(full_data, survival_event, survival_time, case_tag))
}

## Function to check gene selection
check_gene_selection_app <- function(num_genes, gen_select_type, percent_gen_select){
  #Convert text to lowercase
  gen_select_type <- tolower(gen_select_type)
  #Check gen_select_type
  gen <- c("top_bot","abs")
  if(!gen_select_type %in% gen){
    msg <- paste("Invalid gene selection type selected. Choose one of the folowing: ", paste(gen, collapse = ", "))
    shinyalert::shinyalert(title = msg, type = "error")
  }

  #Number of genes to be selected in gene_selection_surv function
  num_gen_select <- trunc((percent_gen_select/100) * num_genes)

  return(num_gen_select)
}

# Function to check filter values
check_filter_values_app <- function(full_data, filter_values, na.rm = TRUE){
  # Check if filter_values is a vector
  if(!is.vector(filter_values)){
    msg <- "filter_values must be a valid values vector"
    shinyalert::shinyalert(title = msg, type = "error")
  }

  #Check if the names of the filter_values are the same as the cols of full_data.
  if(!setequal(names(filter_values), colnames(full_data))){
    msg <- "The name of the filter_values must be the same as the patient name of the full_data (or genes_disease_component)."
    shinyalert::shinyalert(title = msg, type = "error")
  }

  #Omit NAN's values
  if (na.rm == TRUE){
    # Remove rows (subjects) and their filter values with NA's values
    filter_values <- filter_values[colnames(full_data)]
    # Remove filter values and respective rows with NA's values
    filter_values <- stats::na.omit(filter_values)
    full_data <- full_data[,names(filter_values)]
  }
  return(list(full_data, filter_values))
}

# Function to check mapper arguments
check_arg_mapper_app <- function(full_data, filter_values, distance_type, clustering_type, option, linkage_type, na.rm = TRUE){
  #Check distance_type
  distances <- c("cor","euclidean")
  if(!distance_type %in% distances){
    msg <- paste("Invalid distance selected. Choose one of the folowing: ", paste(distances, collapse = ", "))
    shinyalert::shinyalert(title = msg, type = "error")
  }

  #Check clustering_type
  clust_types <- c("hierarchical","PAM")
  if(!clustering_type %in% clust_types){
    msg <- paste("Invalid clustering method selected. Choose one of the folowing: ", paste(clust_types,collapse = ", "))
    shinyalert::shinyalert(title = msg, type = "error")
  }

  optimal_clustering_mode <- "silhouette"

  if(clustering_type == "hierarchical"){

    if(option != "standard"){
      optimal_clustering_mode <- "silhouette"
    }
  }

  msg <- paste("The optimal clustering mode is '", optimal_clustering_mode, "' by default")
  showNotification(msg, duration = NULL)

  #Check linkage_type
  link_types <- c("single","average","complete")
  if(!linkage_type %in% link_types){
    msg <- paste("Invalid linkage method selected. Choose one of the folowing: ", paste(link_types,collapse = ", "))
    shinyalert::shinyalert(title = msg, type = "error")
  }

  # Check if filter_values == [] the filter_values is not calculated yet. So, we checked only the others args
  if(length(filter_values) != 0 & na.rm != "checked"){
    full_data_and_filter_values <- check_filter_values_app(full_data, filter_values)
    full_data <- full_data_and_filter_values[[1]]
    filter_values <- full_data_and_filter_values[[2]]
  }

  return(list(full_data, filter_values, optimal_clustering_mode))
}

## Function Disease-Specific Genomic Analysis (DGSA) with col_row_check, control_tag
DGSA_app <- function(full_data,  survival_time, survival_event, case_tag, col_row_check, control_tag, na.rm = TRUE){

  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data_app(full_data, col_row_check, na.rm)

  #Select the control_tag
  return_check <- check_vectors_app(full_data, survival_time, survival_event, case_tag, control_tag, na.rm)
  full_data <- return_check[[1]]
  survival_event <- return_check[[2]]
  survival_time <- return_check[[3]]
  case_tag <- return_check[[4]]


  ################### BLOCK I: Pre-process. DGSA (using "NT" control_tag) ##############################
  #   Select the normal tissue data gene expression matrix.
  normal_tiss <- full_data[,which(case_tag == control_tag)]

  #   Obtain the gene expression matrix containing the flattened version of the vectors.
  matrix_flatten_normal_tiss <- flatten_normal_tiss(normal_tiss)
  #   Obtain the normal space
  normal_space <- denoise_rectangular_matrix(matrix_flatten_normal_tiss)
  #   Obtain the disease component of the normal_space
  matrix_disease_component <- generate_disease_component(full_data, normal_space)

  ############################################  Create the object #########################################
  DGSA_object <- list("full_data" = full_data,
                      "control_tag" = control_tag,
                      "case_tag" = case_tag,
                      "survival_event" = survival_event,
                      "survival_time" = survival_time,
                      "normal_space" = normal_space,
                      "matrix_disease_component" = matrix_disease_component)

  class(DGSA_object) <- "DGSA_object"

  return(DGSA_object)
}

## Function that calculates the 100 genes with the highest variability in the matrix
#' disease component between samples and use them to draw the heat map.
results_DGSA_app <- function(matrix_disease_component){
  genes_sd <- apply(matrix_disease_component,1,stats::sd)
  selected_genes_sd <- names(genes_sd[order(genes_sd,decreasing = T)])[1:100]
  selected_matrix_disease_component_sd <- matrix_disease_component[selected_genes_sd,]

  # DT::datatable(selected_matrix_disease_component_sd)
  return(selected_matrix_disease_component_sd)
}


# Function  Gene selection and calculation of filter function values.
geneSelection_app <- function(progress, data_object, gen_select_type,
                          percent_gen_select, col_row_check, control_tag, na.rm = TRUE){
  UseMethod("geneSelection_app")
}

# Private function to select Gene without DGSA process
geneSelection_app.DGSA_object <- function(updateProgress, data_object, gen_select_type, percent_gen_select, col_row_check, control_tag, na.rm = TRUE){

  progress$inc(0.2, detail = "")

  matrix_disease_component <- data_object[["matrix_disease_component"]]
  #Check and obtain gene selection (we use in the gene_select_surv)
  num_gen_select <- check_gene_selection_app(nrow(matrix_disease_component),
                                         gen_select_type, percent_gen_select)

  control_tag <- data_object[["control_tag"]]
  survival_event <- data_object[["survival_event"]]
  survival_time <- data_object[["survival_time"]]
  case_tag <- data_object[["case_tag"]]

  progress$inc(0.2, detail = "\nCheking gene selection")


  control_tag_cases <- which(case_tag == control_tag)
  geneSelection_object <- gene_selection(matrix_disease_component, survival_time, survival_event,
                                         control_tag_cases, gen_select_type, num_gen_select)
  progress$inc(0.6, detail = "\nSelecting genes")


  return(geneSelection_object)
}

# Private function to select Gene without DGSA process
geneSelection_app.default <- function(progress, data_object, gen_select_type, percent_gen_select, col_row_check, control_tag, na.rm = TRUE){
  # If we were passed a progress update function, call it
  progress$inc(0.2, detail = "")


  full_data <- data_object[["full_data"]]
  survival_event <- data_object[["survival_event"]]
  survival_time <- data_object[["survival_time"]]
  case_tag <- data_object[["case_tag"]]

  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data_app(full_data, col_row_check, na.rm)

  #Select the control_tag
  return_check <- check_vectors_app(full_data, survival_time, survival_event, case_tag, control_tag, na.rm)
  full_data <- return_check[[1]]
  survival_event <- return_check[[2]]
  survival_time <- return_check[[3]]
  case_tag <- return_check[[4]]

  progress$inc(0.2, detail = "\nChecking vectors")

  #Select the control_tag. This do it inside of the DGSA function
  #Check and obtain gene selection (we use in the gene_select_surv)
  num_gen_select <- check_gene_selection_app(nrow(full_data), gen_select_type, percent_gen_select)

  progress$inc(0.2, detail = "\nChecking gene selection")

  control_tag_cases <- which(case_tag == control_tag)

  geneSelection_object <- gene_selection(full_data, survival_time, survival_event,
                                         control_tag_cases, gen_select_type, num_gen_select)
  progress$inc(0.4, detail = "\nSelecting genes")


  return(geneSelection_object)
}

# Private function to mapper
one_D_Mapper_app <- function(mapper_object_ini){

  full_data <- mapper_object_ini[["full_data"]]
  filter_values <- mapper_object_ini[["filter_values"]]

  #Getting intervals.
  interval_data <- get_intervals_One_D(filter_values, mapper_object_ini[["num_intervals"]], mapper_object_ini[["percent_overlap"]])

  #Getting samples on each interval.
  samp_in_lev <- samples_in_levels(interval_data, filter_values)

  #Clustering all levels.
  test_clust_all_levels <- clust_all_levels(full_data, samp_in_lev, mapper_object_ini[["distance_type"]], mapper_object_ini[["clustering_type"]],
                                            mapper_object_ini[["linkage_type"]], mapper_object_ini[["optimal_clustering_mode"]],  mapper_object_ini[["num_bins_when_clustering"]])
  #Transforming levels into nodes.
  node_samples <- levels_to_nodes(test_clust_all_levels)

  #Computing adjacency matrix.
  adj_matrix_out <- compute_node_adjacency(node_samples)

  node_sizes = unlist(lapply(node_samples,length))
  # average of the filter function of each node
  node_average_filt = lapply(node_samples,function(x,y) mean(y[x]),filter_values)

  # additional parameters
  n_nodes <- length(node_sizes)
  av_node_size <- mean(node_sizes)
  sd_node_size <- stats::sd(node_sizes)

  adj_mat <- adj_matrix_out
  upper_tri <- upper.tri(adj_mat)
  lower_tri <- lower.tri(adj_mat)

  n_connections <- sum(adj_mat[upper_tri] == 1)
  prop_connections <- n_connections/length(adj_mat[upper_tri])

  adj_mat[lower_tri] <- t(adj_mat)[lower_tri]
  diag(adj_mat) <- 0
  n_ramifications <- colSums(adj_mat)-2
  n_ramifications[n_ramifications %in% c(-1,-2)] <- 0
  n_ramifications <- sum(n_ramifications)

  #Generating the object of the output data
  mapper_object <- list("interval_data" = interval_data,
                        "sample_in_level" = samp_in_lev,
                        "clustering_all_levels" = test_clust_all_levels,
                        "node_samples" = node_samples,
                        "node_sizes" = node_sizes,
                        "node_average_filt" = node_average_filt,
                        "adj_matrix" = adj_matrix_out,
                        "n_sizes" = n_nodes,
                        "average_nodes"= av_node_size,
                        "standard_desviation_nodes" = sd_node_size,
                        "number_connections" = n_connections,
                        "proportion_connections" = prop_connections,
                        "number_ramifications" = n_ramifications)

  class(mapper_object) <- "mapper_object"
  return(mapper_object)
}

# Function mapper
mapper_app <- function(full_data, filter_values, num_intervals, percent_overlap,
                   distance_type, clustering_type,
                   num_bins_when_clustering, linkage_type ,
                   optimal_clustering_mode, col_row_check, na.rm=TRUE){
  # Don't call by GSSTDA function
  if (na.rm != "checked"){
    # Check the full_data introduces
    full_data <- check_full_data_app(full_data, col_row_check,na.rm)
  }
  # Check mapper arguments
  check_return <- check_arg_mapper_app(full_data, filter_values, distance_type, clustering_type, optimal_clustering_mode,
                                   linkage_type)

  full_data <- check_return[[1]]
  filter_values <- check_return[[2]]
  optimal_clustering_mode <- check_return[[3]]


  mapper_object_ini <- list("full_data" = full_data,
                            "filter_values" = filter_values,
                            "num_intervals" = num_intervals,
                            "percent_overlap" = percent_overlap/100,
                            "distance_type" = distance_type,
                            "optimal_clustering_mode" = optimal_clustering_mode,
                            "num_bins_when_clustering" = num_bins_when_clustering,
                            "clustering_type" = clustering_type,
                            "linkage_type" = linkage_type,
                            "optimal_clustering_mode" = optimal_clustering_mode)

  class(mapper_object_ini) <- "mapper_initialization"

  mapper_object <- one_D_Mapper_app(mapper_object_ini)

  return(mapper_object)
}


