library(shiny)
library(Seurat)
library(dplyr)
library(DT)
library(ggplot2)
library(plotly)

#Set working directories for local and for deployment on scilifelab serve
#setwd('/srv/shiny-server')
#setwd('/Users/stefanebmeyer/nbis/support_projects/O_Andersson_2024/scripts/scilifelab_serve_shiny2/shiny_december24/')


ui <- fluidPage(
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
                label = "Leiden cluster Resolution:",
                choices = c("0.1" = 'leiden_cluster.01',
                            "0.3" = 'leiden_cluster.03',
                            "0.6" = 'leiden_cluster.06',
                            "0.9" = 'leiden_cluster.09',
                            "1.2" = 'leiden_cluster.12')
            ),
            textInput(
                inputId = "genes",
                label = "Enter Gene Names (comma separated):",
                value = "tph1b, pyyb, neurod1, insl1a, ghrl"
            ),
            actionButton("update", "Update Plot")
        ),
        
        mainPanel(
            tabsetPanel(
                tabPanel('Overview',
                    h4("Data/Metadata Overview"),
                    p("UMAP plots with cells colored by different metadata."),
                    uiOutput('dynamicOverview')
                ),
                tabPanel('Clustering',
                    h4("Leiden Clusters and Dot Plot"),
                    p("Shows the different leiden clusters at the selected resolution. The dot plot shows gene expression relative to other clusters."),
                    plotlyOutput(outputId = "cluster_umap", width = "90%", height = "600px"),
                    plotOutput(outputId = "genePlot", width = "100%", height = "500px")
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
                    DT::dataTableOutput(outputId = 'deg_table')
                )
                #tabPanel('DEGs 2024 feeding',
                #h4('Differentially Expressed Genes for selected clusters (leiden clustering resolution 0.6) between conditions "fed" and "starved" for variable "feeding"'),
                #p(HTML('Differentially expressed genes between "starved" and "fed" for respective cluster. Columns:<br>
                #<b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "starved"<br>
                #<b>pct.1</b>: Percentage of cells in "starved" group expressing the gene<br>
                #<b>pct.2</b>: percentage of cells in "fed" group expressing the gene<br>
                #<b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                #<b>cluster</b>: cluster number')),
                #    DT::dataTableOutput(outputId = 'deg_table_2024_feeding')
                #),
                #tabPanel('DEGs 2024 beta-cells',
                #h4('Differentially expressed genes for selected clusters (leiden clustering resolution 0.6) between conditions "ablated" and "non-ablated" for variable "beta_cells"'),
                #p(HTML('Differentially expressed genes between "ablated" and "non-ablated" for respective cluster. Columns:<br>
                #<b>avg_log2FC</b>: Average log2 fold Change in expression, positive values indicate higher expression in "ablated"<br>
                #<b>pct.1</b>: Percentage of cells in "ablated" group expressing the gene<br>
                #<b>pct.2</b>: percentage of cells in "non-ablated" group expressing the gene<br>
                #<b>p_val_adjust</b>: Adjusted p-value after correction for multiple testing<br>
                #<b>cluster</b>: cluster number')),
                #    DT::dataTableOutput(outputId = 'deg_table_2024_beta')
                #)
            )
        )
    )
)

server <- function(input, output) {
    # Load Seurat objects into reactive values for use in other reactive expressions
    seurat_objects <- reactiveValues()
    
    observe({
      
        #all_objs <- readRDS('/Users/stefanebmeyer/nbis/support_projects/O_Andersson_2024/scripts/scilifelab_serve_shiny2/shiny_december24/data/all_objects_int.rds')
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
        #seurat_objects$DEGs_2024_beta_cells <- readRDS('/home/data/2024_leiden06_beta_DEGs.rds')
        #seurat_objects$DEGs_2024_feeding <- readRDS('/home/data/2024_leiden06_feed_DEGs.rds')
    })  
    
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

        degs <- if (grepl('2023', input$seurat_obj)) seurat_objects$DEGs_2023
                else if (grepl('progenitor', input$seurat_obj)) seurat_objects$DEGs_progenitors
                else if (grepl('2024', input$seurat_obj)) seurat_objects$DEGs_2024
                

        degs_beta <- seurat_objects$DEGs_2024_beta_cells
        degs_feeding <- seurat_objects$DEGs_2024_feeding

        return(list(genes = genes, seurat_obj = seurat_obj, resolution = resolution, degs=degs,
        degs_beta = degs_beta, degs_feeding = degs_feeding))
    }, ignoreNULL = FALSE)

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

        DotPlot(seurat_obj, features = genes, group.by = resolution) + coord_flip()
            
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

        #create cluster dimplot based on input data
        if (input$seurat_obj %in% c('complete_2023', 'complete_2024', '2024_fed',
                                    '2024_nonablated', '2023_starved', 'complete_progenitors',
                                    'progenitors_DMSO', 'progenitors_SMER28')){
            clusterplot <- DimPlot(seurat_obj, reduction='umap_integrated', group.by = resolution, label = TRUE, label.size= 7, repel = TRUE)}
        else if (input$seurat_obj %in% c('2024_ablated', '2024_starved', '2023_ablated', '2023_nonablated')){
            clusterplot <- DimPlot(seurat_obj, reduction='umap', group.by = resolution, label = TRUE, label.size= 7, repel = TRUE)}
        

        #Modify plotly object to make the plot interactive
        if (!is.null(clusterplot)){
            plotly_obj <- ggplotly(clusterplot)
            plotly_obj <- config(plotly_obj, displaylogo = FALSE, modeBarButtonsToRemove=c('pan2d', 'hoverCompareCartesian', 'drawrect', 'select2d')) %>%
                layout(dragmode = 'select', legend='toggleothers', title=list(text='Louvain clusters', size=7))

            #For each cluster in the data. go through the list and modify the text
            for (i in 1:(length(plotly_obj$x$data)-1)){
                htext <- plotly_obj$x$data[[i]]$text
                elements <- unlist(strsplit(htext, '<br />'))
                cluster_info <- grep('^leiden', elements, value = TRUE)
                new_info <- gsub('leiden_cluster\\.([0-9]+): ', 'Cluster: ', cluster_info)
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

    }

shinyApp(ui = ui, server = server)
