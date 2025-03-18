library(shiny)
library(Seurat)
library(dplyr)
library(DT)
library(ggplot2)
library(plotly)
library(gridExtra)
library(shinyjs)

#Set working directories for local and for deployment on scilifelab serve
#setwd('/srv/shiny-server')
#setwd('/Users/stefanebmeyer/nbis/support_projects/O_Andersson_2024/scripts/scilifelab_serve_shiny2/shiny_december24/')


ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  tags$head(
    tags$style(HTML("
        .scrolling-panel { 
            overflow-y: auto;  
        }
    "))
  ),
  
  titlePanel("Data Exploration - Enteroendocrine Cells in Zebrafish"),
  
  sidebarLayout(
    sidebarPanel(
      width=3,
      selectInput(
        inputId = "seurat_obj",
        label = "Select Dataset:",
        choices = c("2023" = "complete_2023",
                    "2023 ablated" = "2023_ablated",
                    "2023 non-ablated" = "2023_nonablated",
                    "2024" = "complete_2024",
                    "2024 ablated" = "2024_ablated",
                    "2024 non-ablated" = "2024_nonablated",
                    "2024 fed" = "2024_fed",
                    "2024 starved" = "2024_starved", 
                    "progenitors complete" = "complete_progenitors", 
                    "progenitors DMSO" = "progenitors_DMSO", 
                    "progenitors SMER28" = "progenitors_SMER28")
      ),
      radioButtons(
        inputId = "resolution",
        label = "Louvain cluster Resolution:",
        choices = c("0.1" = 'leiden_cluster.01',
                    "0.3" = 'leiden_cluster.03',
                    "0.6" = 'leiden_cluster.06',
                    "0.9" = 'leiden_cluster.09',
                    "1.2" = 'leiden_cluster.12',
                    "custom" = 'custom_resolution')
      ),
      radioButtons(
        inputId = 'clust_annot',
        label = 'Cluster/Dotplot annotation:',
        choices = c('Cluster numbers' = 'leiden_cluster.01', 
                    'Cell type' = 'cell_type')
      ),
      
      textInput(
        inputId = "genes",
        label = "Enter Gene Names (comma separated):",
        value = "tph1b, pyyb, neurod1, insl1a, ghrl"
      ),
      actionButton("update", "Update Plot")
    ),
    
    mainPanel(
      uiOutput("tabs_ui")
    )
  )
)

server <- function(input, output, session) {
  # Load Seurat objects into reactive values for use in other reactive expressions
  seurat_objects <- reactiveValues()
  
  observe({
    
    #all_objs <- readRDS('/Users/stefanebmeyer/nbis/support_projects/O_Andersson_2024/objects_november24/all_objects_int.rds')
    #Below path is for scilifelab serve
    all_objs <- readRDS('/home/data/all_objects_int.rds')
    
    #Read in distinct seurat objects
    seurat_objects$complete_2023 <- all_objs$enteroendocrine_2023$seurat_objects$complete_2023
    seurat_objects$'2023_ablated' <- all_objs$enteroendocrine_2023$seurat_objects$ablated
    seurat_objects$'2023_nonablated' <- all_objs$enteroendocrine_2023$seurat_objects$'non-ablated'
    
    seurat_objects$complete_2024 <- all_objs$enteroendocrine_2024$seurat_objects$complete_2024
    seurat_objects$'2024_ablated' <- all_objs$enteroendocrine_2024$seurat_objects$ablated
    seurat_objects$'2024_nonablated' <- all_objs$enteroendocrine_2024$seurat_objects$'non-ablated'
    seurat_objects$'2024_fed' <- all_objs$enteroendocrine_2024$seurat_objects$fed
    seurat_objects$'2024_starved' <- all_objs$enteroendocrine_2024$seurat_objects$starved
    
    seurat_objects$complete_progenitors <- all_objs$progenitors_2024$seurat_objects$complete_progenitors_2024
    seurat_objects$progenitors_DMSO <- all_objs$progenitors_2024$seurat_objects$DMSO
    seurat_objects$progenitors_SMER28 <- all_objs$progenitors_2024$seurat_objects$SMER28
    
    #Read in DEGs
    seurat_objects$DEGs_2023 <- all_objs$enteroendocrine_2023$DEGs
    seurat_objects$DEGs_2024 <- all_objs$enteroendocrine_2024$DEGs
    seurat_objects$DEGs_progenitors <- all_objs$progenitors_2024$DEGs
    
    #DEGs between conditions for MEECs 2024
    seurat_objects$DEGs_2024_beta_cells <- all_objs$enteroendocrine_2024$DEGs$complete_2024$DEGs_beta_06
    seurat_objects$DEGs_2024_feeding <- all_objs$enteroendocrine_2024$DEGs$complete_2024$DEGs_feed_06
    
    #DEGs between conditions for progenitors 2024
    seurat_objects$DEGs_treat <- all_objs$progenitors_2024$DEGs$progenitors_complete_2024$DEGs_treat_09$all
    seurat_objects$DEGs_treat_day1 <- all_objs$progenitors_2024$DEGs$progenitors_complete_2024$DEGs_treat_09$day1
    seurat_objects$DEGs_treat_day3 <- all_objs$progenitors_2024$DEGs$progenitors_complete_2024$DEGs_treat_09$day3
    seurat_objects$DEGs_treat_day5 <- all_objs$progenitors_2024$DEGs$progenitors_complete_2024$DEGs_treat_09$day5
    
    #Add value to observe when a cluster is clicked in clusterplot
    seurat_objects$selected_cluster <- reactiveVal(0)
  })  
  
  #Reset resolution and cluster annotation upon changing the dataset
  observeEvent(input$seurat_obj, {
    updateRadioButtons(session, "resolution", selected = 'leiden_cluster.01')
    # Update the 'clust_annot' radio button to 'Cluster numbers'
    updateRadioButtons(session, "clust_annot", selected = 'leiden_cluster.01')
  })
  
  #Set resolution to corresponding value and disable selection when cell type is selected as annotation
  observeEvent(input$clust_annot, {
    if (input$clust_annot == "cell_type") {
      if (input$seurat_obj %in% c('complete_progenitors','progenitors_DMSO','progenitors_SMER28')){
        updateRadioButtons(session, "resolution", selected = 'leiden_cluster.09')
      }
      else {
        updateRadioButtons(session, "resolution", selected = 'custom_resolution')
      }
      shinyjs::disable("resolution")
    } else {
      updateRadioButtons(session, "resolution", selected = 'leiden_cluster.01')
      shinyjs::enable("resolution")
    }
  })
  
  #Enable annotation by cell type only if dataset has annotations
  observeEvent(input$seurat_obj, {
    if (input$seurat_obj %in% c('complete_progenitors','complete_2024')){
      shinyjs::enable('clust_annot')
    }
    else {
      shinyjs::disable('clust_annot')
    }
  })
  
  #Trigger update button upon changing the seurat object
  observeEvent(input$seurat_obj, {
    seurat_objects$selected_cluster(0)
    shinyjs::click('update')
  })
  
  #Trigger update button upon changing the resolution
  observeEvent(input$resolution, {
    seurat_objects$selected_cluster(0)
    shinyjs::click('update')
  })
  
  #Observe clicks on the clustered UMAP
  observeEvent(event_data('plotly_click', source = 'cluster_umap'), {
    click_data <- event_data('plotly_click', source = 'cluster_umap')
    print(paste0('Click_data: ', click_data))
    if (!is.null(click_data) && !is.na(click_data$curve)) {
      seurat_objects$selected_cluster(click_data$curve)
    } else {
      seurat_objects$selected_cluster(0)  # Default to cluster 0
    }
  }, ignoreInit = TRUE)
  
  # Create a reactive expression for a subset of genes
  gene_data <- eventReactive(input$update, {
    req(input$genes) # Require gene input to proceed
    genes <- unlist(strsplit(input$genes, ",[ ]*"))
    
    #Check that seurat object exists
    req(input$seurat_obj)
    
    seurat_obj <- if (input$seurat_obj == "complete_2023") seurat_objects$complete_2023 
    else if (input$seurat_obj == "2023_ablated") seurat_objects$'2023_ablated'
    else if (input$seurat_obj == "2023_nonablated") seurat_objects$'2023_nonablated'
    else if (input$seurat_obj == "complete_2024") seurat_objects$complete_2024
    else if (input$seurat_obj == "2024_ablated") seurat_objects$'2024_ablated'
    else if (input$seurat_obj == "2024_nonablated") seurat_objects$'2024_nonablated'
    else if (input$seurat_obj == "2024_fed") seurat_objects$'2024_fed'
    else if (input$seurat_obj == "2024_starved") seurat_objects$'2024_starved'
    else if (input$seurat_obj == "complete_progenitors") seurat_objects$complete_progenitors
    else if (input$seurat_obj == "progenitors_DMSO") seurat_objects$progenitors_DMSO
    else if (input$seurat_obj == "progenitors_SMER28") seurat_objects$progenitors_SMER28
    
    resolution <- input$resolution
    #print warning if the selected resolution does not exist in the respective object
    print(resolution)
    if (! resolution %in% colnames(seurat_objects[[input$seurat_obj]]@meta.data)){
      print('Checking resolution')
      showNotification(
        paste0('Resolution ', input$resolution, ' not found in ', input$seurat_obj, '!'),
        type='error',
        duration=10)
      
      #Set resolution to 01
      #Reset resolution and cluster annotation upon changing the dataset
      observeEvent(input$seurat_obj, {
        updateRadioButtons(session, "resolution", selected = 'leiden_cluster.01')
      })
    }
    
    #Set object identities to resolution
    Idents(seurat_obj) <- seurat_obj@meta.data[,colnames(seurat_obj[[resolution]])[1]]
    
    degs <- if (grepl('2023', input$seurat_obj)) seurat_objects$DEGs_2023
    else if (grepl('progenitor', input$seurat_obj)) seurat_objects$DEGs_progenitors
    else if (grepl('2024', input$seurat_obj)) seurat_objects$DEGs_2024
    
    
    degs_beta <- seurat_objects$DEGs_2024_beta_cells
    degs_feeding <- seurat_objects$DEGs_2024_feeding
    degs_treat <- seurat_objects$DEGs_treat
    degs_treat_day1 <- seurat_objects$DEGs_treat_day1
    degs_treat_day3 <- seurat_objects$DEGs_treat_day3
    degs_treat_day5 <- seurat_objects$DEGs_treat_day5
    
    return(list(genes = genes, seurat_obj = seurat_obj, resolution = resolution, degs=degs,
                degs_beta = degs_beta, degs_feeding = degs_feeding, degs_treat = degs_treat, degs_treat_day1 = degs_treat_day1,
                degs_treat_day3 = degs_treat_day3, degs_treat_day5 = degs_treat_day5))
  }, ignoreNULL = FALSE)
  
  #Function to be called reactively to generate QC plots per cluster
  generate_cluster_plots <- function(seurat_obj, cluster, resolution) {
    print(resolution)
    if (input$clust_annot=='cell_type'){
      cell_type_name <- seurat_obj@meta.data$cell_type[seurat_obj[[resolution]]==cluster][1]
      print(cell_type_name)
      seurat_obj <- subset(seurat_obj, subset = cell_type == cell_type_name)
    }
    else {
      seurat_obj <- subset(seurat_obj, idents = cluster)
    }
    
    # Violin plot for mitochondrial genes
    vln_mt <- VlnPlot(seurat_obj, features = c("percent_mt", "percent_ribo", "nFeature_RNA", "nCount_RNA"), ncol = 4, log = TRUE)
    # Feature plot for S phase
    feature_s <- FeaturePlot(seurat_obj, features = c("nFeature_RNA", "percent_mt"), ncol=2, cols=c('blue', 'red'))
    
    #Also generate frequency tables
    if (input$seurat_obj %in% c('complete_progenitors','progenitors_DMSO','progenitors_SMER28')){
      df1 <- as.data.frame(table(seurat_obj@meta.data$treatment))
      names(df1) <- c('Treatment', 'n cells')
      
      df2 <- as.data.frame(table(seurat_obj@meta.data$day))
      names(df2) <- c('Day', 'n cells')
      
    }
    else {
      df1 <- as.data.frame(table(seurat_obj@meta.data$beta_cells))
      names(df1) <- c('Beta-cell state', 'n cells')
      print(df1)
      
      df2 <- as.data.frame(table(seurat_obj@meta.data$feeding))
      names(df2) <- c('Feeding state', 'n cells')
    }
    
    return(list(vln_mt = vln_mt, feature_s = feature_s, df1 = df1, df2 = df2))
  }
  
  #Generate function to create QC plots based on cluster input
  cluster_plots <- reactive({
    
    genes_data <- gene_data()
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    req(genes_data)
    
    #Set resolution manually in case clust_annot == 'cell_type'
    if (input$clust_annot == 'cell_type'){
      if (input$seurat_obj %in% c('complete_2024')){
        resolution = 'custom_resolution'
      }
      else if (input$seurat_obj %in% c('complete_progenitors')){
        resolution = 'leiden_cluster.09'
      }
    }
    #Set seurat object identities to resolution
    print(paste0('Setting identities to ', colnames(seurat_obj[[resolution]])[1]))
    Idents(seurat_obj) <- seurat_obj@meta.data[,colnames(seurat_obj[[resolution]])[1]]
    
    print(paste0('Selected cluster: ',seurat_objects$selected_cluster()))
    selected_cluster_val <- seurat_objects$selected_cluster()
    print(paste0('selected_cluster_val: ',selected_cluster_val))
    
    clst_plts <- generate_cluster_plots(seurat_obj, selected_cluster_val, resolution)
    
    return(clst_plts)
  })
  
  #Create list of UI output tabs
  output$tabs_ui <- renderUI({
    tabs <- list(
      tabPanel('QC',
               h4("Data Quality Control Plots"),
               br(),
               h5('Percent of mitochondrial and ribiosomal genes per sample, and count of total reads per cell and sample
                      (nCount_RNA) and number of expressed genes per cell and sample (nFeature_RNA).'),
               plotOutput(outputId = 'qc_mito', width = "90%", height = "700px"),
               br(),
               h5('FeatureScatter plots'),
               plotOutput(outputId = 'qc_featurescatter', width = "100%", height = "400px"),
               br(),
               h5('Expressed genes on (integrated) UMAP'),
               plotOutput(outputId = 'qc_nfeature_rna', width = "90%", height = "600px"),
               br(),
               h5('Cell Cycle Score on (integrated) UMAP'),
               plotOutput(outputId = 'qc_cellcycle', width = "90%", height = "1100px")
      ),
      tabPanel('Overview',
               h4("Data/Metadata Overview"),
               p("UMAP plots with cells colored by different metadata."),
               uiOutput('dynamicOverview')
      ),
      tabPanel('Clustering',
               h4("Leiden Clusters and Dot Plot"),
               p("Shows the different leiden clusters at the selected resolution. The dot plot shows gene expression relative to other clusters."),
               plotlyOutput(outputId = "cluster_umap", width = "100%", height = "600px"),
               p(),
               p(),
               #dynamically define height for dotplot
               plotOutput(outputId = "genePlot", width = "100%", height = paste0(length(input$genes)*50+450, 'px')),
               h5('QC plots per cluster'), 
               p('Click on a cluster in the UMAP to display the respective QC plots and stats'),
               plotOutput(outputId = "clst_vln", width='100%', height='400px'),
               plotOutput(outputId = "clst_expr", width='100%', height='600px'),
               fluidRow(
                 column(6, DT::dataTableOutput(outputId = 'cdtn1')),
                 column(6, DT::dataTableOutput(outputId = 'cdtn2'))
               ),
               p()
      ),
      tabPanel('Expression Plots',
               h4('Feature/Expression Plots'),
               p('Each cell is colored by relative expression of the respective gene.'),
               div(class='fixed-size-plot', uiOutput("dynamicFeaturePlot"))
      ),
      tabPanel('DEG Table',
               h4('Differentially Expressed Gene Table'),
               p(HTML('Differentially expressed genes between one cluster and all others. Columns: <br>
                    <b>avg_log2FC</b>: Average log2 fold Change in expression <br>
                    <b>pct.1</b>: Percentage of cells in respective cluster expressing the gene<br> 
                    <b>pct.2</b>: percentage of all other cells expressing the gene<br>
                    <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                    <b>cluster</b>: cluster number<br>
                    <b>gene</b>: gene name')),
               DT::dataTableOutput(outputId = 'deg_table'))
    )
    
    if (grepl('ete_2024', input$seurat_obj) & input$resolution=='custom_resolution') {
      tabs <- append(tabs, list(
        tabPanel('DEGs 2024 feeding',
                 h4('Differentially Expressed Genes for selected clusters (leiden clustering resolution 0.6) between conditions "fed" and "starved" for variable "feeding"'),
                 p(HTML('Differentially expressed genes between "starved" and "fed" for respective cluster. Columns:<br>
                      <b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "starved"<br>
                      <b>pct.1</b>: Percentage of cells in "starved" group expressing the gene<br>
                      <b>pct.2</b>: percentage of cells in "fed" group expressing the gene<br>
                      <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                      <b>cluster</b>: cluster number')),
                 DT::dataTableOutput(outputId = 'deg_table_2024_feeding')
        ),
        tabPanel('DEGs 2024 beta-cells',
                 h4('Differentially expressed genes for selected clusters (leiden clustering resolution 0.6) between conditions "ablated" and "non-ablated" for variable "beta_cells"'),
                 p(HTML('Differentially expressed genes between "ablated" and "non-ablated" for respective cluster. Columns:<br>
                      <b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "ablated"<br>
                      <b>pct.1</b>: Percentage of cells in "ablated" group expressing the gene<br>
                      <b>pct.2</b>: percentage of cells in "non-ablated" group expressing the gene<br>
                      <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                      <b>cluster</b>: cluster number')),
                 DT::dataTableOutput(outputId = 'deg_table_2024_beta')
        )
      ))
    }
    
    else if (grepl('ete_progenitors', input$seurat_obj) & input$resolution=='leiden_cluster.09') {
      tabs <- append(tabs, list(
        tabPanel('DEGs progenitors treatment all timepoints',
                 h4('Differentially Expressed Genes for selected clusters (leiden clustering resolution 0.9) between conditions "DMSO" and "SMER28" for variable "treatment"'),
                 p(HTML('Differentially expressed genes between "DMSO" and "SMER28" for respective cluster. Columns:<br>
                      <b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "SMER28"<br>
                      <b>pct.1</b>: Percentage of cells in "SMER28" group expressing the gene<br>
                      <b>pct.2</b>: percentage of cells in "DMSO" group expressing the gene<br>
                      <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                      <b>cluster</b>: cluster number')),
                 DT::dataTableOutput(outputId = 'deg_table_2024_treat')
        ),
        tabPanel('DEGs progenitors treatment day1',
                 h4('Differentially Expressed Genes for selected clusters (leiden clustering resolution 0.9) between conditions "DMSO" and "SMER28" for variable "treatment"'),
                 p(HTML('Differentially expressed genes between "SMER28" and "DMSO" for respective cluster. Columns:<br>
                      <b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "SMER28"<br>
                      <b>pct.1</b>: Percentage of cells in "SMER28" group expressing the gene<br>
                      <b>pct.2</b>: percentage of cells in "DMSO" group expressing the gene<br>
                      <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                      <b>cluster</b>: cluster number')),
                 DT::dataTableOutput(outputId = 'deg_table_2024_treat1')
        ),
        tabPanel('DEGs progenitors treatment day3',
                 h4('Differentially Expressed Genes for selected clusters (leiden clustering resolution 0.9) between conditions "DMSO" and "SMER28" for variable "treatment"'),
                 p(HTML('Differentially expressed genes between "SMER28" and "DMSO" for respective cluster. Columns:<br>
                      <b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "SMER28"<br>
                      <b>pct.1</b>: Percentage of cells in "SMER28" group expressing the gene<br>
                      <b>pct.2</b>: percentage of cells in "DMSO" group expressing the gene<br>
                      <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                      <b>cluster</b>: cluster number')),
                 DT::dataTableOutput(outputId = 'deg_table_2024_treat3')
        ),
        tabPanel('DEGs progenitors treatment day5',
                 h4('Differentially Expressed Genes for selected clusters (leiden clustering resolution 0.9) between conditions "DMSO" and "SMER28" for variable "treatment"'),
                 p(HTML('Differentially expressed genes between "SMER28" and "DMSO" for respective cluster. Columns:<br>
                      <b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "SMER28"<br>
                      <b>pct.1</b>: Percentage of cells in "SMER28" group expressing the gene<br>
                      <b>pct.2</b>: percentage of cells in "DMSO" group expressing the gene<br>
                      <b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                      <b>cluster</b>: cluster number')),
                 DT::dataTableOutput(outputId = 'deg_table_2024_treat5')
        )
      ))
    }
    
    do.call(tabsetPanel, tabs)
    
  })
  
  output$dynamicOverview <- renderUI({
    
    if (grepl('progenitor', input$seurat_obj)){
      tagList(
        plotOutput(outputId = 'treatment_umap', width = "80%", height = "500px"),
        plotOutput(outputId = 'day_umap', width = "80%", height = "500px"),
        plotOutput(outputId = 'sample_umap', width = "80%", height = "500px")
      )
    }
    else {
      tagList(
        plotOutput(outputId = 'beta_umap', width = "80%", height = "500px"),
        plotOutput(outputId = 'feed_umap', width = "80%", height = "500px"),
        plotOutput(outputId = 'sample_umap', width = "80%", height = "500px")
      )
    }
  })
  
  output$clst_vln <- renderPlot({
    req(cluster_plots())
    cluster_plots()$vln_mt
  })
  
  output$clst_expr <- renderPlot({
    req(cluster_plots())
    cluster_plots()$feature_s
  })
  
  output$cdtn1 <- DT::renderDataTable({
    req(cluster_plots())
    DT::datatable(cluster_plots()$df1, options = list(pageLength = 3, searching = FALSE, 
                                                      lengthChange = FALSE, 
                                                      info = FALSE, 
                                                      paging = FALSE))
  })
  
  output$cdtn2 <- DT::renderDataTable({
    req(cluster_plots())
    DT::datatable(cluster_plots()$df2, options = list(pageLength = 3, searching = FALSE, 
                                                      lengthChange = FALSE, 
                                                      info = FALSE, 
                                                      paging = FALSE))
  })
  
  # Render the plot for the selected gene expression
  output$genePlot <- renderPlot({
    # Ensure that the gene data is available
    genes_data <- gene_data()
    
    genes <- genes_data$genes
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    if (input$clust_annot=='cell_type' & 'cell_type' %in% colnames(seurat_obj@meta.data)){
      grouping='cell_type'
      label_size=10
      legend_size=10
    }
    else {
      grouping=resolution
      label_size=14
      legend_size=14
    }
    
    if (grouping=='cell_type'){
      
      #Sort seurat object in reverse in order not to cut off cell names in app
      values <- seurat_obj@meta.data$cell_type
      sorted_values <- sort(unique(values), decreasing = TRUE)
      sorted_obj <- seurat_obj
      sorted_obj[['cell_type']] <- factor(values, levels = sorted_values)
      
      DotPlot(sorted_obj, features = genes, group.by = grouping) + coord_flip() +
        theme(axis.text.x = element_text(angle=45, hjust = 1),
              plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt"))
    }
    else {
      DotPlot(seurat_obj, features = genes, group.by = grouping) + coord_flip()
    }
    
  })
  
  # Render the UMAP plot for the selected resolution
  output$cluster_umap <- renderPlotly({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    if (input$clust_annot=='cell_type' & 'cell_type' %in% colnames(seurat_obj@meta.data)){
      #Set identities and grouping for the UMAP
      if (input$seurat_obj %in% c('complete_2024')){
        Idents(seurat_obj) <- 'custom_resolution'
        
        # Order cell types based on the numeric cluster ID
        ordered_cell_types <- seurat_obj@meta.data$cell_type[order(as.numeric(seurat_obj$custom_resolution))]
        # Remove duplicates while preserving order
        ordered_cell_types <- unique(ordered_cell_types)
        # Update the factor levels
        seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = ordered_cell_types)
        
        grouping='cell_type'
        label_size=10
        legend_size=10
      }
      else if (input$seurat_obj %in% c('complete_progenitors')){
        Idents(seurat_obj) <- 'leiden_cluster.09'
        
        # Order cell types based on the numeric cluster ID
        ordered_cell_types <- seurat_obj@meta.data$cell_type[order(as.numeric(seurat_obj$leiden_cluster.09))]
        # Remove duplicates while preserving order
        ordered_cell_types <- unique(ordered_cell_types)
        # Update the factor levels
        seurat_obj$cell_type <- factor(seurat_obj$cell_type, levels = ordered_cell_types)
        
        grouping='cell_type'
        label_size=10
        legend_size=10
      }
    }
    else {
      
      Idents(seurat_obj) <- seurat_obj@meta.data[,colnames(seurat_obj[[resolution]])[1]]
      grouping=resolution
      label_size=14
      legend_size=14
    }
    
    #create cluster dimplot based on input data
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved', 'complete_progenitors',
                                'progenitors_DMSO', 'progenitors_SMER28')){
      clusterplot <- DimPlot(seurat_obj, reduction='umap_integrated', group.by = grouping, label = TRUE, label.size= label_size, repel = TRUE)}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      clusterplot <- DimPlot(seurat_obj, reduction='umap', group.by = grouping, label = TRUE, label.size= label_size, repel = TRUE)}
    
    #Modify plotly object to make the plot interactive
    if (!is.null(clusterplot)){
      plotly_obj <- ggplotly(source='cluster_umap') %>% event_register('plotly_click')
      plotly_obj <- config(plotly_obj, displaylogo = FALSE, modeBarButtonsToRemove=c('pan2d', 'hoverCompareCartesian', 'drawrect', 'select2d')) %>%
        layout(dragmode = 'select', legend=list(font=list(size=legend_size)), title=list(text='Louvain clusters', size=7))
      
      #For each cluster in the data. go through the list and modify the text
      for (i in 1:(length(plotly_obj$x$data)-1)){
        htext <- plotly_obj$x$data[[i]]$text
        elements <- unlist(strsplit(htext, '<br />'))
        if (input$clust_annot=='cell_type' & 'cell_type' %in% colnames(seurat_obj@meta.data)){
          cluster_info <- grep('^cell_type', elements, value = TRUE)
          new_info <- gsub('cell_type: ', '', cluster_info) 
        }
        else {
          if(! input$resolution == 'custom_resolution'){
            cluster_info <- grep('^leiden', elements, value = TRUE)
            new_info <- gsub('leiden_cluster\\.([0-9]+): ', 'Cluster: ', cluster_info)
          }
          else {
            cluster_info <- grep('^custom', elements, value = TRUE)
            new_info <- gsub('custom_resolution: ', 'Cluster: ', cluster_info)
          }
        }
        plotly_obj$x$data[[i]]$text <- new_info
      } 
      
      plotly_obj
    }
    
  })
  
  #Conditional umap plots colored by conditions
  #Beta-cell state
  output$beta_umap <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved')){
      DimPlot(seurat_obj, reduction='umap_integrated', group.by = 'beta_cells')}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      DimPlot(seurat_obj, reduction='umap', group.by = 'beta_cells')}
    else {
      return(NULL)
    }
  })
  
  #sample_id
  output$sample_umap <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved', 'complete_progenitors',
                                'progenitors_DMSO', 'progenitors_SMER28')){
      DimPlot(seurat_obj, reduction='umap_integrated', group.by = 'sample_id')}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      DimPlot(seurat_obj, reduction='umap', group.by = 'sample_id')}
    else {
      return(NULL)
    }
  })
  
  
  #Feeding state
  output$feed_umap <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved')){
      DimPlot(seurat_obj, reduction='umap_integrated', group.by = 'feeding')}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      DimPlot(seurat_obj, reduction='umap', group.by = 'feeding')}
    else {
      return(NULL)
    }
  })
  
  #Featureplot showing expression of selected genes on umap
  output$featureplot <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved', 'complete_progenitors',
                                'progenitors_DMSO', 'progenitors_SMER28')){
      FeaturePlot(seurat_obj, reduction='umap_integrated', dims = 1:2, features=genes_data$genes, ncol = 1)}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      FeaturePlot(seurat_obj, reduction='umap', dims = 1:2, features=genes_data$genes, ncol = 1)}
  })
  
  #treatment for progenitor cells
  output$treatment_umap <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    DimPlot(seurat_obj, reduction='umap_integrated', dims = 1:2, group.by = 'treatment', ncol = 1)
    
  })
  
  #day for progenitor cells
  output$day_umap <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    DimPlot(seurat_obj, reduction='umap_integrated', dims = 1:2, group.by = 'day', ncol = 1)
    
  })
  
  #Render feature plot to show expression of different genes on UMAP
  output$dynamicFeaturePlot <- renderUI({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    
    # Calculating the dynamic height based on the number of genes
    plotHeight <- 500 * length(genes_data$genes)
    
    # Generating the 'plotOutput' with dynamic height
    plotOutput("featureplot", height = paste0(plotHeight, 'px'), width = '85%')
  })
  
  #Create a table output for displaying the DEGs
  output$deg_table <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs
    
    #Make special adjustment for progenitor subsets, as these are named differently
    if (dataset=="progenitors_DMSO"){
      dataset = 'prog_DMSO'
    }
    else if (dataset == 'progenitors_SMER28'){
      dataset = 'prog_SMER28'
    }
    else if (dataset == 'complete_progenitors'){
      dataset = 'progenitors_complete_2024'
    }
    
    dataframe <- genes_data$degs[[dataset]][[resolution]]
    
    
    #Render dataframe with DT
    DT::datatable(dataframe, options = list(pageLength = 60))
  })
  
  #Create a table output for displaying the DEGs from 2024 data for cluster resolution 0.6, based on feeding
  output$deg_table_2024_feeding <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs_feeding
    
    #Render dataframe with DT
    DT::datatable(deg_list, options = list(pageLength = 60))
  })
  
  #Create a table output for displaying the DEGs from 2024 data for cluster resolution 0.6, based on beta_cells
  output$deg_table_2024_beta <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs_beta
    
    #Render dataframe with DT
    DT::datatable(deg_list, options = list(pageLength = 60))
  })
  
  output$deg_table_2024_treat <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs_treat
    
    #Render dataframe with DT
    DT::datatable(deg_list, options = list(pageLength = 60))
  })
  
  #Create a table output for displaying the DEGs from 2024 data for cluster resolution 0.6, based on feeding
  output$deg_table_2024_treat1 <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs_treat_day1
    
    #Render dataframe with DT
    DT::datatable(deg_list, options = list(pageLength = 60))
  })
  
  #Create a table output for displaying the DEGs from 2024 data for cluster resolution 0.6, based on feeding
  output$deg_table_2024_treat3 <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs_treat_day3
    
    #Render dataframe with DT
    DT::datatable(deg_list, options = list(pageLength = 60))
  })
  
  #Create a table output for displaying the DEGs from 2024 data for cluster resolution 0.6, based on feeding
  output$deg_table_2024_treat5 <- renderDataTable({
    # Ensuring gene_data is available
    genes_data <- gene_data()
    req(genes_data)
    resolution <- genes_data$resolution
    dataset <- input$seurat_obj
    deg_list <- genes_data$degs_treat_day5
    
    #Render dataframe with DT
    DT::datatable(deg_list, options = list(pageLength = 60))
  })
  
  output$qc_mito <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    
    #create Violin plot based on input data
    VlnPlot(seurat_obj, features=c('percent_mt', 'percent_ribo', 'nFeature_RNA', 'nCount_RNA'), group.by = 'sample_id', ncol=2)
  })
  
  output$qc_featurescatter <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    p1 <-FeatureScatter(seurat_obj, feature1 = 'nCount_RNA', feature2 = 'nFeature_RNA', group.by='sample_id') +
      theme(legend.position = "none", 
            axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    p2 <- FeatureScatter(seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent_mt', group.by='sample_id') +
      theme(legend.position = "none", 
            axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    p3 <- FeatureScatter(seurat_obj, feature1 = 'nCount_RNA', feature2 = 'percent_ribo', group.by='sample_id') +
      theme(legend.position = "none", 
            axis.text.x = element_text(angle=45, hjust=1, vjust=1))
    
    grid.arrange(grobs = list(p1, p2, p3), ncol=3, nrow=1)
  })
  
  output$qc_nfeature_rna <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved', 'complete_progenitors',
                                'progenitors_DMSO', 'progenitors_SMER28')){
      FeaturePlot(seurat_obj, reduction='umap_integrated', dims = 1:2, features='nFeature_RNA', ncol = 1)}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      FeaturePlot(seurat_obj, reduction='umap', dims = 1:2, features='nFeature_RNA', ncol = 1)}
  })
  
  output$qc_cellcycle <- renderPlot({
    genes_data <- gene_data()
    req(genes_data) # Ensure that the gene data is available
    seurat_obj <- genes_data$seurat_obj
    resolution <- genes_data$resolution
    
    if (is.null(genes_data)) {
      return()
    }
    if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                '2024_nonablated', '2023_starved', 'complete_progenitors',
                                'progenitors_DMSO', 'progenitors_SMER28')){
      FeaturePlot(seurat_obj, reduction='umap_integrated', dims = 1:2, features=c('S.Score', 'G2M.Score'), ncol = 1)}
    else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
      FeaturePlot(seurat_obj, reduction='umap', dims = 1:2, features=c('S.Score', 'G2M.Score'), ncol = 1)}
  })
  
}

shinyApp(ui = ui, server = server)
