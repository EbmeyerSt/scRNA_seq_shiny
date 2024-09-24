library(shiny)
library(Seurat)
library(dplyr)
library(DT)
library(ggplot2)

#Set working directories for local and for deployment on scilifelab serve
#setwd('/srv/shiny-server')
#setwd('/Users/stefanebmeyer/nbis/support_projects/O_Andersson_2024/scripts/scilifelab_serve_shiny/shiny_app')


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
                            "2024 starved" = "2024_starved")
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
                    p("UMAP plots with cells colored by different metadata. Notice that for all 2024 datasets,
                    integration based on the sequencing lane was done, and the integrated UMAP is shown."),
                    plotOutput(outputId = 'beta_umap', width = "80%", height = "500px"),
                    plotOutput(outputId = 'feed_umap', width = "80%", height = "500px"),
                    plotOutput(outputId = 'sample_umap', width = "80%", height = "500px"),
                    plotOutput(outputId = 'lane_umap', width = "80%", height = "500px")
                ),
                tabPanel('Clustering',
                    h4("Leiden Clusters and Dot Plot"),
                    p("Shows the different leiden clusters at the selected resolution. The dot plot shows gene exprrssion relative to other clusters."),
                    plotOutput(outputId = "cluster_umap", width = "90%", height = "600px"),
                    plotOutput(outputId = "genePlot", width = "100%", height = "500px")
                ),
                tabPanel('Expression Plots',
                    h4('Feature/Expression Plots'),
                    p('Each cell is colored by relative expression of the respective gene.'),
                    div(class='fixed-size-plot', uiOutput("dynamicFeaturePlot"))
                ),
                tabPanel('DEG Table',
                    h4('Differentially Expressed Gene Table'),
                    p('Differentially expressed genes between one cluster and all others. Columns: **avg_log2FC**: Average log2 fold Change in expression,
                    **pct.1**: Percentage of cells in respective cluster expressing the gene, **pct.2***: percentage of all other cells expressing the gene,
                    **p_val_adjust**: Adjusted p-value after correction for multiple testing, **cluster**: cluster number, **gene**: gene name.'),
                    DT::dataTableOutput(outputId = 'deg_table')
                )
            )
        )
    )
)

server <- function(input, output) {
    # Load Seurat objects into reactive values for use in other reactive expressions
    seurat_objects <- reactiveValues()
    
    observe({

        #Read in distinct seurat objects
        seurat_objects$complete_2023 <- readRDS("/home/data/2023_complete_filtered_norm_scaled.rds")
        seurat_objects$'2023_ablated' <- readRDS("/home/data/2023_ablated_filtered_norm_scaled.rds")
        seurat_objects$'2023_nonablated' <- readRDS("/home/data/2023_nonablated_filtered_norm_scaled.rds")
        
        seurat_objects$complete_2024 <- readRDS("/home/data/2024_complete_filtered_norm_scaled.rds")
        seurat_objects$'2024_ablated' <- readRDS("/home/data/2024_ablated_filtered_norm_scaled.rds")
        seurat_objects$'2024_nonablated' <- readRDS("/home/data/2024_nonablated_filtered_norm_scaled.rds")
        seurat_objects$'2024_fed' <- readRDS("/home/data/2024_fed_filtered_norm_scaled.rds")
        seurat_objects$'2024_starved' <- readRDS("/home/data/2024_starved_filtered_norm_scaled.rds")
    
        #Read in DEGs
        seurat_objects$DEGs_2023 <- readRDS('/home/data/DEGs_2023.rds')
        seurat_objects$DEGs_2024 <- readRDS('/home/data/DEGs_2024.rds')
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

        resolution <- input$resolution

        degs <- if (grepl('2023', input$seurat_obj)) seurat_objects$DEGs_2023
                else if (grepl('2024', input$seurat_obj)) seurat_objects$DEGs_2024

        return(list(genes = genes, seurat_obj = seurat_obj, resolution = resolution, degs=degs))
    }, ignoreNULL = FALSE)


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
    output$cluster_umap <- renderPlot({
        genes_data <- gene_data()
        req(genes_data) # Ensure that the gene data is available
        seurat_obj <- genes_data$seurat_obj
        resolution <- genes_data$resolution

        if (is.null(genes_data)) {
            return()
        }

        if (input$seurat_obj %in% c('complete_2023', '2023_ablated', '2023_nonablated')){
        DimPlot(seurat_obj, reduction='umap', group.by = resolution, label = TRUE, label.size = 7, repel = TRUE)}

        else if (input$seurat_obj %in% c('complete_2024', '2024_ablated', '2024_nonablated', '2024_fed', '2024_starved')){
        DimPlot(seurat_obj, reduction = "umap_harmony", group.by = resolution, label = TRUE, label.size= 7, repel = TRUE)}

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

        if (input$seurat_obj %in% c('complete_2023', '2023_ablated', '2023_nonablated')){
            DimPlot(seurat_obj, reduction='umap', group.by = 'beta_cells')}
        else if (input$seurat_obj %in% c('complete_2024', '2024_ablated', '2024_nonablated', '2024_fed', '2024_starved')){
            DimPlot(seurat_obj, reduction='umap_harmony', group.by = 'beta_cells')}
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

        if (input$seurat_obj %in% c('complete_2023', '2023_ablated', '2023_nonablated')){
            DimPlot(seurat_obj, reduction='umap', group.by = 'sample_id')}
        else if (input$seurat_obj %in% c('complete_2024', '2024_ablated', '2024_nonablated', '2024_fed', '2024_starved')){
            DimPlot(seurat_obj, reduction='umap_harmony', group.by = 'sample_id')}
        else {
            return(NULL)
        }
    })

    #lane
    output$lane_umap <- renderPlot({
        genes_data <- gene_data()
        req(genes_data) # Ensure that the gene data is available
        seurat_obj <- genes_data$seurat_obj
        resolution <- genes_data$resolution

        if (is.null(genes_data)) {
            return()
        }

        if (input$seurat_obj %in% c('complete_2023', '2023_ablated', '2023_nonablated')){
            DimPlot(seurat_obj, reduction='umap', group.by = 'lane')}
        else if (input$seurat_obj %in% c('complete_2024', '2024_ablated', '2024_nonablated', '2024_fed', '2024_starved')){
            DimPlot(seurat_obj, reduction='umap_harmony', group.by = 'lane')}
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

        if (input$seurat_obj %in% c('complete_2023', '2023_ablated', '2023_nonablated')){
            DimPlot(seurat_obj, reduction='umap', group.by = 'feeding')}
        else if (input$seurat_obj %in% c('complete_2024', '2024_ablated', '2024_nonablated', '2024_fed', '2024_starved')){
            DimPlot(seurat_obj, reduction='umap_harmony', group.by = 'feeding')}
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
        if (input$seurat_obj %in% c('complete_2023', '2023_ablated', '2023_nonablated')){
            FeaturePlot(seurat_obj, reduction='umap', dims = 1:2, features=genes_data$genes, ncol = 1)}
        else if (input$seurat_obj %in% c('complete_2024', '2024_ablated', '2024_nonablated', '2024_fed', '2024_starved')){
            FeaturePlot(seurat_obj, reduction='umap_harmony', dims = 1:2, features=genes_data$genes, ncol = 1)}
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

        dataframe <- deg_list[[dataset]][[resolution]]

        #Render dataframe with DT
        DT::datatable(dataframe, options = list(pageLength = 60))
    })

    }

shinyApp(ui = ui, server = server)