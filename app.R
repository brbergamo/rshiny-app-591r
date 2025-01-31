library(shiny)
library(bslib)
library(tidyverse)
library(ggplot2)
library(colourpicker)
library(DESeq2)
library(DT)
library(gplots)
library(fgsea)

# Set the maximum file upload size to 100 MB to handle large datasets
options(shiny.maxRequestSize = 100 * 1024^2)

# Define the User Interface (UI) of the app
ui <- navbarPage(
  title = "Rshiny Report",  # Title of the app
  theme = bs_theme(version = 5),  # Use Bootstrap 5 for consistent styling
  
  # Tab for uploading and viewing sample information
  nav_panel(
    "Samples",
    sidebarLayout(
      sidebarPanel(
        # File input to upload sample metadata
        fileInput("samples", "Upload sample info:", 
                  buttonLabel = "Upload...", accept = c(".csv", ".tsv")),
        # Button to process the uploaded sample metadata
        actionButton("submit_samples", "Submit")
      ),
      mainPanel(
        tabsetPanel(
          # Display summary of uploaded sample metadata
          tabPanel("Summary", tableOutput("summaryTable")),
          # Display full sample metadata table
          tabPanel("Table", DTOutput("dataTable")),
          # Plot based on sample metadata (e.g., histogram of conditions)
          tabPanel("Plots", plotOutput("histogramPlot"))
        )
      )
    )
  ),
  
  # Tab for uploading and analyzing count data
  nav_panel(
    "Counts",
    sidebarLayout(
      sidebarPanel(
        # File input to upload count data
        fileInput("counts", "Upload Counts:", buttonLabel = "Upload...", accept = c(".csv", ".tsv")),
        # Button to process the uploaded count data
        actionButton("submit_counts", "Submit"),
        # Slider to filter genes based on variance percentile
        sliderInput("variancePercentile", "Select Variance Percentile:",
                    min = 0, max = 100, value = 95, step = 1),
        # Slider to filter genes based on the number of non-zero samples
        sliderInput("nonZeroSamples", "Minimum Non-Zero Samples:",
                    min = 1, max = 100, value = 5, step = 1)
      ),
      mainPanel(
        h3("Counts data"),
        tabsetPanel(
          # Display summary of filters applied to count data
          tabPanel("Filter Summary", tableOutput("filterSummary")),
          # Diagnostic plots (e.g., median vs variance)
          tabPanel("Diagnostic Plots", 
                   plotOutput("medianVsVariance"),
                   plotOutput("medianVsZeros")),
          # Heatmap of filtered count data
          tabPanel("Heatmap", plotOutput("heatmapPlot")),
          # PCA analysis and visualization
          tabPanel("PCA",
                   sidebarLayout(
                     sidebarPanel(
                       h4("PCA Options"),
                       # Dropdowns to select principal components for PCA scatter plot
                       selectInput("pc_x", "X-axis PC:", choices = c("PC1", "PC2"), selected = "PC1"), 
                       selectInput("pc_y", "Y-axis PC:", choices = c("PC1", "PC2"), selected = "PC2"), 
                       actionButton("run_pca", "Run PCA")  # Button to trigger PCA analysis
                     ),
                     mainPanel(
                       plotOutput("pcaScatterPlot")  # Display PCA scatter plot
                     )
                   ))
        )
      )
    )
  ),
  
  # Tab for differential expression analysis (DE)
  nav_panel(
    "DE",
    sidebarLayout(
      sidebarPanel(
        h4("Differential Expression Analysis"),
        # File input to upload differential expression results
        fileInput("deresults", "Upload Differential Expression Results:",
                  buttonLabel = "Upload...", accept = c(".csv", ".tsv")),
        # Numeric input to set adjusted p-value cutoff for filtering DE results
        numericInput("padj_cutoff", "P-adjust Cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
        # Button to process differential expression results
        actionButton("submit_deresults", "Submit")
      ),
      mainPanel(
        h3("DE analysis content goes here"),
        tabsetPanel(
          # Display DE results as a table
          tabPanel("Table", DTOutput("DETable")),
          # Display volcano plot of DE results
          tabPanel("Volcano Plot", plotOutput("dePlot"))
        )
      )
    )
  ),
  
  # Tab for Gene Set Enrichment Analysis (GSEA)
  nav_panel(
    "GSEA",
    sidebarLayout(
      sidebarPanel(
        h4("Gene Set Enrichment Analysis"),
        # File input to upload differential expression results for GSEA
        fileInput("deresults", "Upload Differential Expression Results:",
                  buttonLabel = "Upload...", accept = c(".csv", ".tsv")),
        # Button to execute GSEA analysis
        actionButton("run_fgsea", "Run GSEA")
      ),
      mainPanel(
        tabsetPanel(
          # Bar plot of top pathways from GSEA
          tabPanel("Top Pathways Plot", 
                   sidebarLayout(
                     sidebarPanel(
                       # Slider to select the number of pathways to display in the plot
                       sliderInput("topPathways", "Number of Top Pathways to Display:",
                                   min = 1, max = 50, value = 10, step = 1)
                     ),
                     mainPanel(plotOutput("fgseaPlot"))
                   )),
          # Table of filtered FGSEA results
          tabPanel("Filter Table", 
                   sidebarLayout(
                     sidebarPanel(
                       h4("Filter FGSEA Results"),
                       # Slider to filter FGSEA results by adjusted p-value
                       sliderInput("padjFilter", "P-adjust Cutoff:", min = 0, max = 1, value = 0.05, step = 0.01),
                       # Radio buttons to filter FGSEA results by NES direction
                       radioButtons("nesFilter", "NES Direction:",
                                    choices = list("All" = "all", "Positive" = "positive", "Negative" = "negative"),
                                    selected = "all"),
                       # Button to download filtered FGSEA results as a CSV file
                       downloadButton("downloadTable", "Download Filtered Table")
                     ),
                     mainPanel(
                       h3("Filtered FGSEA Results Table"),
                       DTOutput("filteredFgseaTable")  # Display filtered FGSEA results in a table
                     ))),
          # Scatter plot of NES vs. -log10(adjusted p-value)
          tabPanel("Scatter Plot",
                   sidebarLayout(
                     sidebarPanel(
                       h4("Scatter Plot Settings"),
                       # Slider to filter gene sets for the scatter plot by adjusted p-value
                       sliderInput("scatterPadjFilter", "P-adjust Cutoff:", min = 0, max = 1, value = 0.05, step = 0.01)
                     ),
                     mainPanel(
                       h3("Scatter Plot of NES vs. -log10(padj)"),
                       plotOutput("nesScatterPlot")  # Display scatter plot of NES vs. -log10(padj)
                     )
                   ))
        )
      )
    )
  )
)

# Define server logic for the app
server <- function(input, output, session) {
  
  samplesData <- reactiveVal(NULL)  # Reactive value to store sample metadata
  pca_result <- reactiveVal(NULL)  # Reactive value to store PCA results
  pca_variance <- reactiveVal(NULL)  # Reactive value to store PCA variance explained
  countsData <- reactiveVal(NULL)  # Reactive value to store count data
  fgseaResults <- reactiveVal(NULL)  # Reactive value to store GSEA results
  
  # Event to process uploaded sample metadata
  observeEvent(input$submit_samples, {
    req(input$samples)  # Ensure a file is uploaded
    df <- read.delim(input$samples$datapath, header = TRUE)  # Read the uploaded file
    samplesData(df)  # Store the data in a reactive value
  })
  
  # Render summary of sample metadata
  output$summaryTable <- renderTable({
    req(samplesData())  # Ensure data is available
    df <- samplesData()
    tibble::tibble(
      ColumnName = names(df),  # Column names
      Type = sapply(df, class),  # Data types of columns
      Values = sapply(df, function(column) paste(unique(column), collapse = ", "))  # Unique values
    )
  })
  
  # Render full sample metadata table
  output$dataTable <- renderDT({
    req(samplesData())
    datatable(samplesData(), 
              options = list(pageLength = 10, autoWidth = TRUE), 
              rownames = FALSE)
  })
  
  # Histogram plot for sample metadata
  output$histogramPlot <- renderPlot({
    req(samplesData())
    df <- samplesData()
    ggplot(df, aes(x = Condition, fill = Condition)) + 
      geom_histogram(stat = "count") + 
      theme_minimal() + 
      labs(title = "Condition Frequency", x = "Condition", y = "Count")
  })
  
  
  # Event to process uploaded count data
  observeEvent(input$submit_counts, {
    req(input$counts)  # Ensure a file is uploaded
    df_counts <- read.delim(input$counts$datapath, header = TRUE)  # Read the uploaded count data file
    countsData(df_counts)  # Store the count data in a reactive value for further processing
  })
  
  # Render a summary table of filtering results
  output$filterSummary <- renderTable({
    req(countsData())  # Ensure count data is available
    
    df <- countsData()  # Retrieve the reactive count data
    total_genes <- nrow(df)  # Total number of genes (rows in the data)
    total_samples <- ncol(df) - 1  # Total number of samples (columns excluding gene IDs)
    
    # Filtering process based on user inputs
    filtered_data <- df %>%
      mutate(Variance = apply(select(., -1), 1, var)) %>%  # Calculate variance for each gene
      filter(Variance >= quantile(Variance, input$variancePercentile / 100, na.rm = TRUE)) %>%  # Filter by variance percentile
      mutate(NonZeroCount = apply(select(., where(is.numeric)), 1, function(row) sum(row > 0))) %>%  # Count non-zero samples for each gene
      filter(NonZeroCount >= input$nonZeroSamples)  # Filter by minimum non-zero samples
    
    passing_genes <- nrow(filtered_data)  # Number of genes passing filters
    not_passing_genes <- total_genes - passing_genes  # Number of genes not passing filters
    passing_percentage <- round((passing_genes / total_genes) * 100, 2)  # Percentage of passing genes
    not_passing_percentage <- round((not_passing_genes / total_genes) * 100, 2)  # Percentage of non-passing genes
    
    # Return a summary as a tibble
    tibble::tibble(
      Metric = c("Number of Samples", 
                 "Total Number of Genes", 
                 "Number of Genes Passing Filter", 
                 "% of Genes Passing Filter", 
                 "Number of Genes Not Passing Filter", 
                 "% of Genes Not Passing Filter"),
      Value = c(total_samples, 
                total_genes, 
                passing_genes, 
                paste0(passing_percentage, "%"), 
                not_passing_genes, 
                paste0(not_passing_percentage, "%"))
    )
  })
  
  # Plot: Median Count vs Variance
  output$medianVsVariance <- renderPlot({
    req(countsData())  # Ensure count data is available
    
    df <- countsData()  # Retrieve the reactive count data
    
    # Calculate metrics for each gene
    df_metrics <- df %>%
      mutate(
        Variance = apply(select(., -1), 1, var),  # Variance for each gene
        MedianCount = apply(select(., -1), 1, median),  # Median count for each gene
        NonZeroCount = apply(select(., -1), 1, function(row) sum(row > 0))  # Number of non-zero samples per gene
      ) %>%
      mutate(
        PassFilter = Variance >= quantile(Variance, input$variancePercentile / 100, na.rm = TRUE) &  # Check if variance passes filter
          NonZeroCount >= input$nonZeroSamples  # Check if non-zero count passes filter
      )
    
    # Plot Median Count vs Variance with filter highlights
    ggplot(df_metrics, aes(x = MedianCount, y = Variance, color = PassFilter)) +
      geom_point(alpha = 0.6) +  # Add points with transparency
      scale_color_manual(values = c("TRUE" = "darkblue", "FALSE" = "lightblue")) +  # Custom color scale
      scale_x_log10() + scale_y_log10() +  # Logarithmic scales for better visualization
      labs(
        title = "Median Count vs Variance",
        x = "Median Count (log scale)",
        y = "Variance (log scale)",
        color = "Pass Filter"  # Legend title
      ) +
      theme_minimal()  # Minimal theme for clean aesthetics
  })
  
  # Plot: Median Count vs number of zeros
  output$medianVsZeros <- renderPlot({
    req(countsData())
    
    df <- countsData()
    
    
    df_metrics <- df %>%
      mutate(
        Variance = apply(select(., -1), 1, var),
        MedianCount = apply(select(., -1), 1, median),
        NonZeroCount = apply(select(., -1), 1, function(row) sum(row > 0))
      ) %>%
      mutate(
        PassFilter = Variance >= quantile(Variance, input$variancePercentile / 100, na.rm = TRUE) &
          NonZeroCount >= input$nonZeroSamples
      )
    
    
    ggplot(df_metrics, aes(x = MedianCount, y = ncol(df) - 1 - NonZeroCount, color = PassFilter)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("TRUE" = "darkgreen", "FALSE" = "lightgreen")) +
      scale_x_log10() +
      labs(
        title = "Median Count vs Number of Zeros",
        x = "Median Count (log scale)",
        y = "Number of Zeros",
        color = "Pass Filter"
      ) +
      theme_minimal()
  })
  
  
  
  # Heatmap of filtered counts
  output$heatmapPlot <- renderPlot({
    req(countsData())  # Ensure count data is available
    
    # Apply filtering to count data based on user criteria
    df <- countsData() %>%
      mutate(
        Variance = apply(select(., -1), 1, var),  # Calculate variance for each gene
        NonZeroCount = apply(select(., -1), 1, function(row) sum(row > 0))  # Count non-zero samples for each gene
      ) %>%
      filter(
        Variance >= quantile(Variance, input$variancePercentile / 100, na.rm = TRUE),  # Filter by variance percentile
        NonZeroCount >= input$nonZeroSamples  # Filter by minimum non-zero samples
      ) %>%
      select(-Variance, -NonZeroCount)  # Remove auxiliary columns used for filtering
    
    # Convert the filtered data to a matrix format for heatmap generation
    counts_matrix <- as.matrix(df[, -1])  # Exclude gene IDs
    rownames(counts_matrix) <- df[[1]]  # Set gene IDs as row names
    
    # Apply log transformation to handle skewed counts
    log_transformed <- log1p(counts_matrix)  # log(x + 1) to handle zeros
    
    # Generate heatmap
    heatmap.2(
      x = log_transformed,
      scale = "row",  # Scale data by rows
      dendrogram = "both",  # Cluster rows and columns
      trace = "none",  # Remove trace lines
      col = colorRampPalette(c("blue", "white", "red"))(50),  # Custom color palette
      key = TRUE,  # Show color key
      key.title = "Expression",  # Title for the color key
      key.xlab = "Log-scaled Counts",  # Label for the color key
      main = "Clustered Heatmap of Filtered Counts",  # Title of the heatmap
      cexRow = 0.5,  # Font size for row labels
      cexCol = 0.7  # Font size for column labels
    )
  })
  
  
  
  #pca
  
  observeEvent(input$run_pca, {
    # Check that the required data is available
    req(countsData(), samplesData())
    
    # Process the counts data:
    counts <- countsData() %>%
      # Compute the variance for each gene (row) after excluding the first column (gene IDs)
      mutate(Variance = apply(select(., -1), 1, var)) %>%
      # Filter by the variance percentile defined by the user input
      filter(Variance >= quantile(Variance, input$variancePercentile / 100, na.rm = TRUE)) %>%
      # Count the number of samples with a non-zero expression value for each gene
      mutate(NonZeroCount = apply(select(., where(is.numeric)), 1, function(row) sum(row > 0))) %>%
      # Filter by the minimum number of non-zero samples as defined by the user input
      filter(NonZeroCount >= input$nonZeroSamples) %>%
      # Remove auxiliary columns before proceeding
      select(-Variance, -NonZeroCount)
    
    # Retrieve sample information
    sample_info <- samplesData()
    
    # Prepare data for PCA
    # Remove the gene ID column from counts and transpose so that samples are rows and genes are columns
    counts_matrix <- counts[, -1]
    rownames(counts_matrix) <- counts[[1]] # Set gene IDs as row names
    counts_matrix <- t(as.matrix(counts_matrix))
    
    # Run PCA on the processed counts data
    pca <- prcomp(counts_matrix, scale = TRUE)
    
    # Calculate the percentage of variance explained by each principal component
    pca_variance(round(100 * pca$sdev^2 / sum(pca$sdev^2), 2))
    
    # Convert PCA results into a data frame
    pca_scores <- as.data.frame(pca$x)
    pca_scores$Sample <- rownames(pca_scores)  # Add sample names to the results
    
    # Merge the PCA results with the sample metadata
    merged_data <- inner_join(pca_scores, sample_info, by = "Sample")
    
    # Store the merged PCA results for downstream use
    pca_result(merged_data)
    
    # Dynamically update the dropdown options for selecting PCs to plot
    pc_choices <- paste0("PC", 1:ncol(pca$x))
    updateSelectInput(session, "pc_x", choices = pc_choices, selected = "PC1")
    updateSelectInput(session, "pc_y", choices = pc_choices, selected = "PC2")
  })
  
  
  # Render PCA scatter plot
  output$pcaScatterPlot <- renderPlot({
    req(pca_result(), input$pc_x, input$pc_y)
    
    # Use merged data
    merged_data <- pca_result()
    variance <- pca_variance()
    

    
    ggplot(merged_data, aes(x = !!sym(input$pc_x), y = !!sym(input$pc_y), color = Condition)) +
      geom_point(alpha = 0.6, size = 3) +
      labs(
        title = "PCA Scatter Plot (Samples Colored by Group)",
        x = paste0(input$pc_x, " (", variance[as.integer(gsub("PC", "", input$pc_x))], "% variance)"),
        y = paste0(input$pc_y, " (", variance[as.integer(gsub("PC", "", input$pc_y))], "% variance)"),
        color = "Condition"
      ) +
      theme_minimal()
  })
  
  
  # DE
  
  output$DETable <-renderDT({
    req(input$deresults)
    
    df <- read_delim(input$deresults$datapath)
    datatable(df, 
              options = list(pageLength = 10, autoWidth = TRUE), 
              rownames = FALSE)
  })
  
  
  
  output$dePlot <- renderPlot({
    req(input$deresults,input$padj_cutoff)
    
    df <- read_delim(input$deresults$datapath)
    
    df <- df %>%
      mutate(Significant = padj < input$padj_cutoff)
    
    ggplot(df, aes(x = log2FoldChange, y = -log10(padj),  color = Significant)) +
      geom_point(alpha = 0.6, size = 3) +
      theme_minimal() + scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"), name = paste0("P-adjust <", input$padj_cutoff))
    
  })
 
  
  # Run GSEA
  
  observeEvent(input$run_fgsea, {
    req(input$deresults)
    
    # Read differential expression results
    res <- read_delim(input$deresults$datapath)
    
    # Prepare ranked genes
    ranked_genes <- res %>%
      filter(!is.na(log2FoldChange) & !is.na(symbol)) %>%
      arrange(desc(log2FoldChange)) %>%
      pull(log2FoldChange)
    names(ranked_genes) <- res %>%
      filter(!is.na(log2FoldChange) & !is.na(symbol)) %>%
      arrange(desc(log2FoldChange)) %>%
      pull(symbol)
    ranked_genes <- ranked_genes[!duplicated(names(ranked_genes))]
    
    # Load pathways
    pathways <- gmtPathways("h.all.v2024.1.Hs.symbols.gmt")
    
    # Run GSEA
    fgsea_results <- fgsea(pathways = pathways, stats = ranked_genes)
    fgsea_table <- as_tibble(fgsea_results) %>%
      arrange(padj)
    
    # Store results in fgseaResults
    fgseaResults(fgsea_table)
  })
  
    
  # gsea plot  
  output$fgseaPlot <- renderPlot({
    req(fgseaResults())
    
    # Select top pathways based on the slider input
    top_pathways <- fgseaResults() %>%
      arrange(padj) %>%
      head(input$topPathways)
    
    # Create the barplot
    ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = NES)) +
      geom_bar(stat = "identity", width = 0.8) +
      coord_flip() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
      labs(
        title = "Top Enriched Pathways",
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)",
        fill = "NES"
      ) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12)
      )
  })
  

  
  # filter fgsea results
  filteredFgseaResults <- reactive({
    req(fgseaResults())
    
    # Apply p-adjust filter
    filtered_table <- fgseaResults() %>%
      filter(padj <= input$padjFilter)
    
    # Apply NES direction filter
    if (input$nesFilter == "positive") {
      filtered_table <- filtered_table %>% filter(NES > 0)
    } else if (input$nesFilter == "negative") {
      filtered_table <- filtered_table %>% filter(NES < 0)
    }
    
    filtered_table
  })


  # Render the filtered table
  output$filteredFgseaTable <- renderDT({
    req(filteredFgseaResults())
    
    filteredFgseaResults() %>%
      mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", "))) %>%
      datatable(
        options = list(pageLength = 10, autoWidth = TRUE, order = list(list(1, "asc"))), 
        rownames = FALSE
      )
  })
  
  # Download button functionality
  output$downloadTable <- downloadHandler(
    filename = function() {
      paste("filtered_fgsea_results_", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      # Convert list columns to strings
      processed_data <- filteredFgseaResults() %>%
        mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse = ", ")))
      
      # Write the processed data to the file
      write.csv(processed_data, file, row.names = FALSE)
    }
  )
  
  
  
  
  
  # Reactive table for scatter plot filtering
  scatterFilteredResults <- reactive({
    req(fgseaResults())
    fgseaResults() %>%
      mutate(log_padj = -log10(padj))  # Calculate -log10(padj) for y-axis
  })
  
  # Scatter plot
  output$nesScatterPlot <- renderPlot({
    req(scatterFilteredResults())
    
    filtered_results <- scatterFilteredResults()
    
    # Highlight gene sets below the adjusted p-value threshold
    filtered_results <- filtered_results %>%
      mutate(
        Highlight = ifelse(padj <= input$scatterPadjFilter, "Significant", "Not Significant")
      )
    
    ggplot(filtered_results, aes(x = NES, y = log_padj, color = Highlight)) +
      geom_point(size = 3, alpha = 0.6) +
      scale_color_manual(values = c("Significant" = "blue", "Not Significant" = "grey")) +
      labs(
        title = "Scatter Plot of NES vs. -log10(adjusted p-value)",
        x = "Normalized Enrichment Score (NES)",
        y = "-log10(adjusted p-value)",
        color = "Gene Set Significance"
      ) +
      theme_minimal() +
      theme(
        legend.position = "top",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 12)
      )
  })
  
  
  
}

# Run the application
shinyApp(ui = ui, server = server)

