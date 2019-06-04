library(shiny)
library(data.table)
library(DT)
library(RColorBrewer)
library(umap)
library(reticulate)

# Had to do this to set the bioconductor repo, and also installed MASS with
# devtools::install_version("MASS", "7.3-51.1")
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.8/bioc"))

#reticulate::virtualenv_create(envname = "python_environment", python = "python3")
#reticulate::use_virtualenv("python_environment", required = TRUE)
#reticulate::virtualenv_install("python_environment", packages = c('pandas','numpy'))
#reticulate::use_virtualenv("python_environment", required = TRUE)
#py_install('gsea_incontext')
#reticulate::py_config()
#print(reticulate::py_config())
source_python('run_gsea_incontext.py')

shinyServer(function(input, output) {
  
  shinyjs::useShinyjs()
  
  # Initialize cache
  VALUES <- reactiveValues(
    APP_DATA = list(
      ALL_ANNOT = ALL_ANNOT,
      ALL_CELLS = ALL_CELLS,
      ALL_DRUGS = ALL_DRUGS,
      ALL_TISSUES = ALL_TISSUES
    ),
    EXPT_DATA = list(
      ORIG442 = list(
        UMAP_2D = NULL,
        UMAP_3D = NULL
      ),
      CMAP01 = list(
        UMAP_2D = NULL,
        UMAP_3D = NULL
      ),
      NCI60 = list(
        UMAP_2D = NULL,
        UMAP_3D = NULL
      ),
      INTERSECTED = list(
        UMAP_2D = NULL,
        UMAP_3D = NULL,
        NES_UMAPS = list(`All collections` = all_umap)
      )
    ),
    USER_INPUT = list(
      SELECTED_EXPTS = NULL,
      SAVED_BG_EXPTS = NULL,
      SELECTED_ROW = NULL,
      GMT = NULL
    ),
    RESULTS = list(
      TOP_PATHWAYS = NULL,
      FGSEA_RESULTS = NULL,
      FGSEA_RANKS = NULL,
      FGSEA_PATHWAYS = NULL,
      ORDERED_NEIGHBORS = NULL,
      INCONTEXT_CSV = F
    ),
    DATA_FOR_DOWNLOAD = list(
      SAMPLE_ANNOTATION_TABLE = NULL,
      GSEA_TABLE = NULL,
      INCONTEXT_TABLE = NULL,
      JOINED_TABLE = NULL
    )
  )
  
  # --------- TAB 1: BACKGROUND DATA EXPLORER ---------------------------------
  
  # ------------ TOP ROW ------------------------------------------------------
  
  # ----------------- LEFT SIDEBAR --------------------------------------------
  
  # Info modal
  observeEvent(input$info_filters, {
    showModal(modalDialog(
      title = 'Background data explorer',
      HTML('In order to run the GSEA-InContext algorithm ('),
      a('Powers et al. Bioinformatics, 2018', href = 'https://sites.google.com/site/jamesccostello4/'),
      '), you first need to select a "background set" of experiments that you will be comparing your own experiment to. This tab allows you to explore the 5,718 experiments we curated for this purpose.',
      strong('Tailor your background set of experiments in order to address the biological question of interest.'),
      HTML('<ul><li>First, using the left sidebar menu, select the dataset/platform and use the additional fields below to subset of the dataset by drug, tissue, or cell line.'),
      HTML('<ul><li>As filters are added, the UMAP plot will automatically highlight the subset of selected experiments. The barplot on the far right will also change, showing the percentage of the selected experiments having each of the Hallmarks gene sets significantly positively or negatively enriched.</li></ul>'),
      HTML('</li><li>Alternatively, instead of using the filters in the sidebar, a subset of experiments can also be selected by clicking and dragging on the UMAP plot.'),
      HTML('</li><li>When you are satisfied with the subset of experiments selected, click "Save" (lower left) to save that subset of background experiments for running GSEA-InContext on the "Run GSEA-InContext" tab.</li></ul>'),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # User-defined input filters
  filter_p <- reactive(input$filter_platform)
  filter_d <- reactive(input$filter_drug)
  filter_c <- reactive(input$filter_cell)
  filter_t <- reactive(input$filter_tissue)
  fix_barplot_order <- reactive(input$fix_barplot)
  active_gmt <- reactive(input$filter_gmt)
  
  # Clear filters when the Clear button is clicked
  observeEvent(input$button_clear, {
    VALUES$USER_INPUT$SELECTED_EXPTS <- NULL
    shinyjs::reset('tab1_side_panel')
    #js$resetSelected() # defined in ui.R
  })
  
  # ----------------- RIGHT / MAIN PANEL --------------------------------------
  
  # Get pre-computed UMAP data based on the platform selected
  # Just using the 2D ones for now, but could use 3D too in the future
  umap_data <- reactive({
    if (filter_p() == 'GEO 442 (HGU133 Plus 2.0)'){
      if (is.null(VALUES$EXPT_DATA$ORIG442$UMAP_2D)){
        VALUES$EXPT_DATA$ORIG442$UMAP_2D <- readRDS('data/Rank_UMAPs/Orig442_ranks_UMAP.rds')
      }
      dat = VALUES$EXPT_DATA$ORIG442$UMAP_2D$layout
    } else if (filter_p() == 'CMap Build 01 (HGU133A)'){
      if (is.null(VALUES$EXPT_DATA$CMAP01$UMAP_2D)){
        VALUES$EXPT_DATA$CMAP01$UMAP_2D <- readRDS('data/Rank_UMAPs/CMap01_ranks_UMAP.rds')
      }
      dat = VALUES$EXPT_DATA$CMAP01$UMAP_2D$layout
    } else if (filter_p() == 'NCI-60 (HGU133A 2.0)'){
      if (is.null(VALUES$EXPT_DATA$NCI60$UMAP_2D)){
        VALUES$EXPT_DATA$NCI60$UMAP_2D <- readRDS('data/Rank_UMAPs/NCI-60_ranks_UMAP.rds')
      }
      dat = VALUES$EXPT_DATA$NCI60$UMAP_2D$layout
    } else {
      if (is.null(VALUES$EXPT_DATA$INTERSECTED$UMAP_2D)){
        VALUES$EXPT_DATA$INTERSECTED$UMAP_2D <- readRDS('data/Rank_UMAPs/intersected_ranks_UMAP.rds')
      }
      dat = VALUES$EXPT_DATA$INTERSECTED$UMAP_2D$layout
    }
    annot = VALUES$APP_DATA$ALL_ANNOT  # regardless of data set, use master table of annotations
    row.names(annot) = annot$rnk_list
    dat = as.data.frame(dat)
    row.names(dat) = gsub('_gene', '', row.names(dat))
    names(dat) = c('V1', 'V2')
    return(list(dat = dat, annot = annot))
  })
  
  # Save manually selected points to reactive values (this allows them to be reset)
  observe({
    VALUES$USER_INPUT$SELECTED_EXPTS <- event_data('plotly_selected', source = 'A')
    # Output for event_data() is 4 columns: curveNumber, pointNumber, x, y
  })
  
  # If points were manually selected, highlight those, otherwise
  # select expts to highlight on the UMAP plot based on user filters
  get_selected <- reactive(VALUES$USER_INPUT$SELECTED_EXPTS)
  
  expts_to_highlight <- reactive({
    
    if (!is.null(get_selected())){
      # If points were manually selected, round the x and y coordinates
      # to join them with the umap data 
      # (Using pointNumber would require knowing the tissues (curveNumbers)
      # and would be a pain as far as I can tell)
      s = as.data.frame(get_selected())
      s$V1 = as.numeric(round(s$x, 5))
      s$V2 = as.numeric(round(s$y, 5))
      dat = data_for_umap_plot()
      dat$V1 = as.numeric(round(dat$V1, 5))
      dat$V2 = as.numeric(round(dat$V2, 5))
      s = join(s, dat)
      return(s$rnk_list)
      
    } else {
      # Filter by platform
      if (filter_p() == 'All'){
        expts_to_plot = VALUES$APP_DATA$ALL_ANNOT
      } else{
        platform = dataset_platforms[filter_p()]
        expts_to_plot = VALUES$APP_DATA$ALL_ANNOT[VALUES$APP_DATA$ALL_ANNOT$platform == platform,]
      }
      # Filter by drug
      if (length(filter_d()) > 0){
        if (any(expts_to_plot$drug %in% filter_d())){
          expts_to_plot = expts_to_plot[expts_to_plot$drug %in% filter_d(),]
        } else{
          return(NULL)
        }
      }
      # Filter by cell line
      if (length(filter_c()) > 0){
        if (any(expts_to_plot$cell_line %in% filter_c())){
          expts_to_plot = expts_to_plot[expts_to_plot$cell_line %in% filter_c(),]
        } else{
          return(NULL)
        }
      }
      # Filter by tissue type
      if (length(filter_t()) > 0){
        if (any(expts_to_plot$tissue %in% filter_t())){
          expts_to_plot = expts_to_plot[expts_to_plot$tissue %in% filter_t(),]
        } else{
          return(NULL)
        }
      }
      return(expts_to_plot$rnk_list)
    }
  })
  
  # Get the relevant UMAP data for data set and join with sample annotations
  data_for_umap_plot <- reactive({
    data_for_plot = umap_data()
    umap_values = data.frame(data_for_plot$dat)
    sample_annots = data_for_plot$annot
    umap_values$rnk_list = sample_annots[row.names(umap_values), 'rnk_list']
    umap_values = join(umap_values, sample_annots)
    return(umap_values)
  })
  
  # Create title for the UMAP plot
  output$umap_plot_title <- renderText({
    return(paste0('UMAP plot highlighting ', length(expts_to_highlight()), ' experiments'))
  })
  
  # Title barplot based on gmt selected
  output$title_barplot <- renderText({
    return(paste0('GSEA results for selected experiments - ', active_gmt()))
  })
  
  # Plot UMAP and highlight selected expts (if any)
  output$umap_plot <- renderPlotly({
    umap_values = data_for_umap_plot()
    if (!is.null(expts_to_highlight())){
      umap_values$highlight = ifelse(umap_values$rnk_list %in% expts_to_highlight(),
                                     'yes', 'no')
    } else{
      umap_values$highlight = 'no'
    }
    
    if (sum(umap_values$highlight == 'no') == 0){
      # If no expts were highlighted, plot them all with color background and black outline
      p <- plot_ly(umap_values,
                   source = 'A',
                   type = 'scatter', mode = 'markers',
                   x = ~V1,
                   y = ~V2,
                   color = ~tissue,
                   colors = tissue_cols,
                   text = ~paste('ID: ', rnk_list, 
                                 '<br>Tissue: ', tissue,
                                 '<br>Drug: ', drug,
                                 '<br>Dose: ', expt_dose),
                   marker = list(
                     size = 5,
                     line = list(
                       color = 'black',
                       width = 1
                     )
                   )
      ) %>% layout(dragmode = 'select') 
    } else{
      # If some expts were highlighted, plot the non-highlighted ones with no outline
      p <- plot_ly(umap_values,
                   source = 'A',
                   type = 'scatter', mode = 'markers',
                   x = ~V1,
                   y = ~V2,
                   color = ~tissue,
                   colors = tissue_cols,
                   text = ~paste('ID: ', rnk_list, 
                                 '<br>Tissue: ', tissue,
                                 '<br>Drug: ', drug,
                                 '<br>Dose: ', expt_dose),
                   marker = list(
                     opacity = 0.5,
                     size = 5,
                     line = list(
                       color = 'transparent',
                       width = 1
                     )
                   )
      ) %>%
        # Then plot the highlighted expts and give them a black outline
        add_trace(
          x = umap_values[umap_values$highlight == 'yes', 'V1'],
          y = umap_values[umap_values$highlight == 'yes', 'V2'],
          text = paste('ID: ', umap_values[umap_values$highlight == 'yes', 'rnk_list'],
                       '<br>Tissue: ', umap_values[umap_values$highlight == 'yes', 'tissue'],
                       '<br>Drug: ', umap_values[umap_values$highlight == 'yes', 'drug'],
                       '<br>Dose: ', umap_values[umap_values$highlight == 'yes', 'expt_dose']),
          name = 'User-filtered',
          marker = list(
            color = tissue_cols[umap_values[umap_values$highlight == 'yes', 'tissue']],
            size = 5,
            line = list(
              color = 'black',
              width = 1
            )
          ),
          inherit = F,
          type = 'scatter', mode = 'markers'
        ) %>% layout(dragmode = 'select')
    }
  })
  
  # Get relevant data for barplot
  data_for_barplot <- reactive({
    tmp = BARPLOT_DATA
    g = active_gmt()
    if (g != 'All collections'){
      tmp = tmp[tmp$Collection %in% names(which(GMTS == g)),
                names(tmp) %in% c('Collection', 'Term', expts_to_highlight())]
    }
    VALUES$DATA_FOR_DOWNLOAD$GSEA_TABLE <- tmp
    return(tmp)
  })
  
  # Display barplot based on all / selected experiments
  output$barplot <- renderPlotly({
    
    # Get the current subset of data to plot in barplot
    all_barplot_data = data_for_barplot()
    all_barplot_data$Term = all_barplot_data$Collection = NULL
    
    # Calculate % of selected expts with each pathway pos/neg enriched
    n = ncol(all_barplot_data)
    data_for_plot = data.frame(gene_sets = row.names(all_barplot_data),
                               pos = apply(all_barplot_data, 1, function(x){
                                  sum(x > 0, na.rm = T)/n
                                }), 
                               neg = apply(all_barplot_data, 1, function(x){
                                  -1*sum(x < 0, na.rm = T)/n
                                }))
    data_for_plot = data_for_plot[order(data_for_plot$pos - abs(data_for_plot$neg)),]
    data_for_plot$gene_sets = gsub('HALLMARK_', '', data_for_plot$gene_sets)
    data_for_plot$pos = round(data_for_plot$pos, 3)*100
    data_for_plot$neg = round(data_for_plot$neg, 3)*100
    row.names(data_for_plot) = data_for_plot$gene_sets
    
    # If there are more than 50 gene sets, just show the top and bottom 25 so the plot is sane
    if (nrow(data_for_plot) > 50){
      data_for_plot = data_for_plot[c(1:25,
                                      (nrow(data_for_plot)-26):nrow(data_for_plot)),]
    }
    
    # Save top 3 positive and top 3 negative pathways
    tops = c(data_for_plot$pos[length(data_for_plot$pos):(length(data_for_plot$pos)-2)], 
             data_for_plot$neg[1:3])
    names(tops) = data_for_plot$gene_sets[c(length(data_for_plot$pos):(length(data_for_plot$pos)-2), 
                                            1:3)]
    VALUES$RESULTS$TOP_PATHWAYS <- tops
    
    # Optionally, retain original pathway ordering instead of changing based on subset
    if (fix_barplot_order()){
      data_for_plot = unique(na.omit(data_for_plot[FIXED_BARPLOT_PATHWAYS,]))
    }
  
    # Plot positive bars
    plot_ly(data_for_plot, source = 'B',
            x = ~gene_sets, y = ~pos, type = 'bar', 
            name = 'Positively enriched (q < .05)', 
            marker = list(color = cols[1])) %>%
      # Add negative bars below
      add_bars(y = ~neg, name = 'Negatively enriched (q < .05)', 
               marker = list(color = cols[10])) %>%
      layout(xaxis = list(title = 'Gene set (hover for name)',
                          ticktext = list(''),
                          tickvals = list(1:50),
                          tickmode = 'array',
                          automargin = TRUE), 
             yaxis = list(title = paste0('% experiments (n =', n, ' selected)'),
                          range = c(-100,100)), 
             barmode = 'overlay')
  })
  
  # Format the text that will be displayed under the plots
  output$message_top_pathways <- renderText({
    if (is.null(expts_to_highlight())){
      return('Oops! No experiments match the criteria selected. Click the CLEAR ALL FILTERS button on the left to try again!')
    } else{
      tops = VALUES$RESULTS$TOP_PATHWAYS
      return(paste0('Top positively enriched:\n\t', names(tops)[1], 
                    ' (', round(as.numeric(tops[1]),3), '% of experiments), ',
                    names(tops)[2], 
                    ' (', round(as.numeric(tops[2]),3), '% of experiments), ',
                    names(tops)[3], 
                    ' (', round(as.numeric(tops[3]),3), '% of experiments)',
                    '\nTop negatively enriched:\n\t', names(tops)[4],
                    ' (', -1*round(as.numeric(tops[4]),3), '% of experiments), ',
                    names(tops)[5],
                    ' (', -1*round(as.numeric(tops[5]),3), '% of experiments), ',
                    names(tops)[6],
                    ' (', -1*round(as.numeric(tops[6]),3), '% of experiments)'))
    }
  })
  
  # ------------ BOTTOM ROW ---------------------------------------------------
  
  # ----------------- LEFT SIDEBAR --------------------------------------------
  
  # Display number of experiments selected
  output$message_n_rnks_selected <- renderText({
    paste0('Experiments currently selected: ', length(expts_to_highlight()))
  })
  
  # Save highlighted expts for background set
  saved_rnks_selected <- eventReactive(input$button_save, {
    # Save list of rnk IDs
    VALUES$USER_INPUT$SAVED_BG_EXPTS <- expts_to_highlight()
    # Save subsetted annotation table for easy downloading later
    tmp = VALUES$APP_DATA$ALL_ANNOT[VALUES$USER_INPUT$SAVED_BG_EXPTS,]
    if (any(is.na(tmp$rnk_list))){
      print('NA values in the rnk lists saved!')
    }
    VALUES$DATA_FOR_DOWNLOAD$SAMPLE_ANNOTATION_TABLE <- tmp
    return(paste0('Saved ', length(expts_to_highlight()),
                  ' experiments to use as background set'))
  })
  
  # Show message when save button clicked
  output$message_saved_rnks_selected <- renderText({
    saved_rnks_selected()
  })
  
  # Download table of sample annotations of selected background experiments
  output$download_rnks_selected <- downloadHandler(
    filename = function() {
      platforms = ifelse(filter_p() == 'All', 'allPlatforms', paste0(filter_p(), sep = '_'))
      drugs = ifelse(length(filter_d() > 0), paste0(filter_d(), sep = '_'), 'allDrugs')
      tissues = ifelse(length(filter_t() > 0), paste0(filter_t(), sep = '_'), 'allTissues')
      return(paste0(paste0(c('Background_set', platforms, drugs, tissues), collapse = '_'), '.csv'))
    },
    content = function(file) {
      write.csv(VALUES$DATA_FOR_DOWNLOAD$SAMPLE_ANNOTATION_TABLE, 
                file, row.names = F)
    }
  )
  
  # Download table of NES values of selected background experiments
  output$download_nes_selected <- downloadHandler(
    filename = function() {
      return('GSEA_NES_for_selected_expts.csv')
    },
    content = function(file) {
      write.csv(VALUES$DATA_FOR_DOWNLOAD$GSEA_TABLE, file,
                row.names = F)
    }
  )
  
  # ----------------- RIGHT / MAIN PANEL --------------------------------------
  
  # Display experiment annotations for selected expts in table
  output$summary_table <- renderDataTable({
    datatable(ALL_ANNOT[expts_to_highlight(),], 
              selection = 'none',
              rownames = F,
              colnames = c('Comparison ID', 'Cell Line', 'Drug', 'Tissue',
                           'Platform', 'Timepoint', 'Treatment Dose'))
  })
  
  # Temporarily replace with incontext results to test
  #output$summary_table <- renderDataTable({
  #  dat = read.csv('out/gsea_incontext.incontext.gene_sets.report.csv', stringsAsFactors = F)
  #  datatable(dat)
  #})
  
  # --------- TAB 2: RUN GSEA-INCONTEXT ---------------------------------------
  
  # ----------------- LEFT SIDEBAR --------------------------------------------
  
  # Info modal
  observeEvent(input$info_k, {
    showModal(modalDialog(
      title = 'Upload & Compare',
      HTML('In order to run the GSEA-InContext algorithm ('),
      a('Powers et al. Bioinformatics, 2018', href = 'https://sites.google.com/site/jamesccostello4/'),
      '), you first need to select a "background set" of experiments that you will be comparing your own experiment to. This tab allows you to upload your own experiment and explore its k nearest neighbors, of the 5,718 experiments we curated.',
      strong('Your background set of experiments will be determined by UMAP similarity.'),
      HTML('<ul><li>First, using the left sidebar menu, select the dataset/platform to begin with.</li>'),
      HTML('<li>Upload your ranked list (two-column format: GeneSymbol, logFC)</li>'),
      HTML('<li>Adjust the k parameter to increase or decrease the number of neighbors selected.</li>'),
      HTML('<ul><li>As k is adjusted, the UMAP plot will automatically highlight the subset of selected experiments. The barplot on the far right will also change, showing the percentage of the selected experiments having each of the Hallmarks gene sets significantly positively or negatively enriched.</li></ul>'),
      HTML('</li><li>When you are satisfied with the subset of experiments selected, click "Save" (lower left) to save that subset of background experiments for running GSEA-InContext on the "Run GSEA-InContext" tab.</li></ul>'),
      easyClose = TRUE,
      footer = NULL
    ))
  })
  
  # Validate user input ranked file
  input_rnk <- reactive(input$input_file$datapath)
  input_filename <- reactive(input$input_file$name)
  validated_data <- reactive(makeValidRnk(input_rnk()))
  kk <- reactive(input$k)
  
  # Select the gmt file
  observeEvent(input$filter_gmt2, {
    VALUES$USER_INPUT$GMT <- gmtPathways(as.character(GMT_PATHS[input$filter_gmt2]))
  })
  
  get_gmt_path <- reactive(return(as.character(GMT_PATHS[input$filter_gmt2])))
  
  # Calculate fgsea when the user inputs a file
  observeEvent(input$input_file, {
    pathways = VALUES$USER_INPUT$GMT
    ranks = validated_data()
    ranks = ranks$df
    names(ranks) = c('V1', 'V2')
    ranks = setNames(ranks$V2, ranks$V1)
    VALUES$RESULTS$FGSEA_RANKS <- ranks
    VALUES$RESULTS$FGSEA_PATHWAYS <- pathways
    res = fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
    VALUES$RESULTS$FGSEA_RESULTS <- res
  })
  
  # If they already uploaded a file, use that one
  # Otherwise they can upload one on this tab
  # TODO can clean up this logic and just check if values$FINAL_RNK is null
  
  # Display the number of background experiments the user has set
  output$message_bg_to_use <- renderText({
    return(paste0('Using ', length(VALUES$USER_INPUT$SAVED_BG_EXPTS),
                  ' experiments as background set'))
  })
  
  get_preranked <- reactive(return(VALUES$RESULTS$FGSEA_RESULTS))
  get_bg_expts <- reactive(return(VALUES$USER_INPUT$SAVED_BG_EXPTS))
  
  # When run button is clicked, get preranked results, compute InContext, and join them
  observeEvent(input$button_run_incontext,{
    
    if (is.null(input_rnk())){
      return(NULL)
    } else{
      
      # Get background experiment ranks
      bg_expts = get_bg_expts()
      all_rnks = readRDS('data/intersected_compiled_reranked_rnks.rds')
      all_rnks$Gene = NULL
      
      # Only keep genes in both the user's data and the bg expts
      user_ranks = validated_data()
      user_ranks = user_ranks$df
      genes_in_common = sort(intersect(user_ranks$GeneSymbol, row.names(all_rnks)))
      user_ranks = user_ranks[user_ranks$GeneSymbol %in% genes_in_common,]
      all_rnks = all_rnks[row.names(all_rnks) %in% genes_in_common, names(all_rnks) %in% bg_expts]
      print(dim(all_rnks))
      
      # Convert the rank indicies into a matrix of ordered gene names per experiment (column)
      # This is very slow right now - can we use do.call or something?
      for (i in 1:ncol(all_rnks)){
        print(i)
        if (i == 1){
          gene_mat = data.frame(row.names(all_rnks)[order(all_rnks[,i])])
        } else{
          gene_mat = cbind(gene_mat, data.frame(row.names(all_rnks)[order(all_rnks[,i])]))
        }
      }
      gene_mat = t(gene_mat)
      
      # These test paths work: (first test gene sets, then rnk, then b)
      #rnk = "gsea_incontext_test_data/GSE5145_DEG_Expt1_Control_vs_Group1_gene.rnk"
      #gene_sets = "data/gene_sets/hallmarks.gmt"
      #background_csv = "gsea_incontext_test_data/MCF7_22_background_lists_permuted_x100.csv"
      
      # Write out files whose paths will be passed to run_gsea_incontext()
      rnk = 'out/incontext_user_rnk.rnk'
      write.table(user_ranks, rnk, sep = '\t',
                  row.names = F, col.names = F, quote = F)
      gene_sets = get_gmt_path()
      print(gene_sets)
      
      # TEMP: if the background set is > 100 expts, just save that csv
      #if (length(bg_expts) > 100){
        #VALUES$RESULTS$INCONTEXT_CSV <- T
      background_csv = 'out/incontext_bg_input.csv'
      write.table(gene_mat, background_csv, sep = ',', row.names = F, col.names = F, quote = F)
      #}
      
      # ELSE: Create a csv that is based on the beta binomial
      
      # Run GSEA-InContext
      incontext_results = run_gsea_incontext(rnk, gene_sets, background_csv, outdir = 'out',
                                             permutation_num = nrow(gene_mat))
      dat = read.csv('out/gsea_incontext.incontext.gene_sets.report.csv',
                                                           stringsAsFactors = F)
      dat = dat[,c(1:3,5)]
      names(dat) = c('Term', 'ES', 'NES', 'InContext_qval')
      
      # Get preranked results, which were computed when the user uploaded their file
      preranked_table = get_preranked()
      preranked_table = preranked_table[,c(1,3)]
      names(preranked_table) = c('Term', 'Preranked_qval')
      
      # Save data frame for later plotting and stuff
      VALUES$DATA_FOR_DOWNLOAD$INCONTEXT_TABLE <- join(dat, preranked_table, type = 'full')
    }
  })
  
  # ----------------- RIGHT / MAIN PANEL --------------------------------------
  
  # Show instructions before a plot has been generated
  output$message_run_incontext <- renderText({
    if(is.null(VALUES$DATA_FOR_DOWNLOAD$INCONTEXT_TABLE)){
      return('Upload a .rnk file then click "RUN GSEA-INCONTEXT"')
    } else{
      return(NULL)
    }
  })
  
  # Save selected row to reactive values (this allows them to be reset)
  observe({
    VALUES$USER_INPUT$SELECTED_ROW <- input$gsea_incontext_table_rows_selected
  })
  
  # If user has clicked on a row in the table, show the GSEA enrichment plot 
  # otherwise show the plot for the top pathway by q-value
  output$fgsea_plot <- renderPlot({
    if (is.null(VALUES$DATA_FOR_DOWNLOAD$INCONTEXT_TABLE)){
      return(NULL)
    } else{
      fgsea_results = VALUES$RESULTS$FGSEA_RESULTS
      ranks = VALUES$RESULTS$FGSEA_RANKS
      pathways = VALUES$RESULTS$FGSEA_PATHWAYS
      
      if (is.null(VALUES$USER_INPUT$SELECTED_ROW)){
        # If no row is selected, show the top most significant
        top_gene_set = as.character(VALUES$DATA_FOR_DOWNLOAD$INCONTEXT_TABLE[1, 1])
        return(plotEnrichment(pathways[[top_gene_set]],
                              ranks) + labs(title = top_gene_set))
      } else{
        # Show the gene set that was clicked on
        selected_gene_set = as.character(VALUES$DATA_FOR_DOWNLOAD$INCONTEXT_TABLE[VALUES$USER_INPUT$SELECTED_ROW, 1])
        return(plotEnrichment(pathways[[selected_gene_set]],
                              ranks) + labs(title = selected_gene_set))
      }
    }
  })
  
  get_current_results <- reactive(return(VALUES$DATA_FOR_DOWNLOAD$INCONTEXT_TABLE))
  
  # Plot GSEA q-vals vs GSEA-InContext q-vals
  output$qq_plot <- renderPlotly({
    if (is.null(get_current_results())){
      return(NULL)
    } else{
      dat = get_current_results()
      dat = dat[order(dat$InContext_qval, decreasing = F),]
      dat = na.omit(dat)  # the NAs are from the full join I'm going on GSEAPreranked + GSEA-InContext
      dat[dat$Preranked_qval == 0, 'Preranked_qval'] = 1e-5
      dat[dat$InContext_qval == 0, 'InContext_qval'] = 1e-5
      dat$Preranked_qval_log10 = -log10(dat$Preranked_qval)
      dat$InContext_qval_log10 = -log10(dat$InContext_qval)
      
      # Save for download
      VALUES$DATA_FOR_DOWNLOAD$JOINED_TABLE <- dat
      
      # Color by quadrant
      dat$Quadrant = ifelse(dat$Preranked_qval < .05 & dat$InContext_qval < .05, names(quadrant_cols)[1],
                            ifelse(dat$Preranked_qval < .05 & dat$InContext_qval >= .05, names(quadrant_cols)[2],
                                   ifelse(dat$Preranked_qval >= .05 & dat$InContext_qval < .05, names(quadrant_cols)[3],
                                          names(quadrant_cols)[4])))
      p <- plot_ly(dat,
                   source = 'Q',
                   type = 'scatter', 
                   mode = 'markers',
                   x = ~Preranked_qval_log10,
                   y = ~InContext_qval_log10,
                   text = ~paste('Pathway: ', Term),
                   color = ~Quadrant,
                   colors = quadrant_cols,
                   marker = list(
                     #symbol = TODO
                     size = 5,
                     line = list(
                       color = 'black',
                       width = 1
                     )
                   ))
      p <- add_segments(p, x = 0, xend = 5, y = -log10(.05), yend = -log10(.05),
                        mode = 'lines', inherit = F, hoverinfo = 'none', showlegend = F,
                        line = list(color = 'grey')) %>%
        add_segments(x = -log10(.05), xend = -log10(.05), y = 0, yend = 5,
                     mode = 'lines', inherit = F, hoverinfo = 'none', showlegend = F,
                     line = list(color = 'grey')) %>%
        layout(xaxis = list(title = 'GSEA Preranked q-value',
                            showline = F,
                            showgrid = F), 
               yaxis = list(title = 'GSEA-InContext q-value',
                            showline = F,
                            showgrid = F))
    }
  })
  
  # Display the results table under the plots
  output$gsea_incontext_table <- DT::renderDataTable({
    if (is.null(get_current_results())){
      return(NULL)
    } else{
      dat = get_current_results()
      dat = dat[order(dat$InContext_qval, decreasing = F),]
      return(DT::datatable(dat, selection = 'single'))
    }
  })
  
  output$download_gseaincontext <- downloadHandler(
    filename = function() {
      return('GSEA-InContext_and_GSEAPreranked_results.csv')
    },
    content = function(file) {
      write.csv(VALUES$DATA_FOR_DOWNLOAD$JOINED_TABLE, file,
                row.names = F)
    }
  )
  
  # --------- TAB 3: UPLOAD AND COMPARE ---------------------------------------
  
  # ------------ TOP ROW ------------------------------------------------------
  
  # ----------------- LEFT SIDEBAR --------------------------------------------
  
  # Info modal
  
  # User-defined input filters
  fix_barplot_order3 <- reactive(input$fix_barplot3)
  active_gmt3 <- reactive(input$filter_gmt3)
  
  # Validate user input ranked file
  input_rnk3 <- reactive(input$input_file3$datapath)
  input_filename3 <- reactive(input$input_file3$name)
  validated_data3 <- reactive(makeValidRnk(input_rnk3()))
  kk <- reactive(input$k)
  
  # Read in active .gmt file
  observeEvent(input$filter_gmt3, {
    VALUES$USER_INPUT$GMT <- gmtPathways(as.character(GMT_PATHS[input$filter_gmt3]))
  })
  
  # Calculate fgsea when the user inputs a file
  observeEvent(input$input_file3, {
    pathways = VALUES$USER_INPUT$GMT
    ranks = validated_data3()
    ranks = ranks$df
    names(ranks) = c('V1', 'V2')
    ranks = setNames(ranks$V2, ranks$V1)
    VALUES$RESULTS$FGSEA_RANKS <- ranks
    VALUES$RESULTS$FGSEA_PATHWAYS <- pathways
    res = fgsea(pathways, ranks, minSize=15, maxSize=500, nperm=1000)
    VALUES$RESULTS$FGSEA_RESULTS <- res
  })
  
  # ----------------- RIGHT / MAIN PANEL --------------------------------------
  
  # Title the plot based on value of k
  output$title_k <- renderText({
    if (is.null(input_rnk())){
      return(paste0(c('UMAP plot (upload a file to see neighbors)')))
    } else {
      return(paste0(c('UMAP plot (k = ', kk(), ' neighbors)')))
    }
  })
  
  # Title barplot based on gmt selected
  output$title_barplot3 <- renderText({
    return(paste0('GSEA results for selected experiments - ', active_gmt3()))
  })
  
  # Show UMAP plot based on the selected gene set collection
  umap_data3 <- reactive({
    g = active_gmt3()
    umaps_list = VALUES$EXPT_DATA$INTERSECTED$NES_UMAPS
    if (!hasName(umaps_list, g)){
      u = readRDS(as.character(UMAP_PATHS[g]))
      VALUES$EXPT_DATA$INTERSECTED$NES_UMAPS[[length(VALUES$EXPT_DATA$INTERSECTED$NES_UMAPS)+1]] <- u
      names(VALUES$EXPT_DATA$INTERSECTED$NES_UMAPS)[length(VALUES$EXPT_DATA$INTERSECTED$NES_UMAPS)] <- g
    }
    return(VALUES$EXPT_DATA$INTERSECTED$NES_UMAPS[[g]])
  })
  
  # Plot pre-computed data, or new data with user's expt mapped onto it
  get_gsea_results <- reactive(VALUES$RESULTS$FGSEA_RESULTS)
  
  data_for_plotting3 <- reactive({
    u = umap_data3()
    xy_for_plot = as.data.frame(u$layout)
    gene_sets = colnames(u$data)
    names(xy_for_plot) = c('V1', 'V2')
    xy_for_plot$rnk_list = row.names(xy_for_plot)
    xy_for_plot = join(xy_for_plot, VALUES$APP_DATA$ALL_ANNOT)
    
    if (!is.null(input_rnk3())){
      # When a valid file has been uploaded, project onto the UMAP
      user_data = as.data.frame(get_gsea_results())
      row.names(user_data) = user_data$pathway
      user_data = as.data.frame(user_data[gene_sets, 'NES'])
      row.names(user_data) = gene_sets
      user_data = na.omit(user_data)
      predicted_pt = predict(u, t(user_data))
      
      # Save matrix of all xy values, including user's uploaded expt
      # (Add tissue label to user's expt so it shows up in plot)
      xy_for_plot = rbind.fill(xy_for_plot, 
                               data.frame(rnk_list = 'UserInput',
                                          tissue = 'UserInput',
                                          V1 = as.numeric(predicted_pt)[1],
                                          V2 = as.numeric(predicted_pt)[2]))
      
      # Compute knn and save ordered neighbor list
      ks = get.knnx(xy_for_plot[,c('V1', 'V2')], query = predicted_pt, 
                    k = nrow(xy_for_plot))
      
      # Order by distance to uploaded expt
      neighbors_in_order = xy_for_plot$rnk_list[ks$nn.index]
      VALUES$RESULTS$ORDERED_NEIGHBORS <- neighbors_in_order
      row.names(xy_for_plot) = xy_for_plot$rnk_list
      xy_for_plot = xy_for_plot[neighbors_in_order,]
    }
    row.names(xy_for_plot) = xy_for_plot$rnk_list
    return(xy_for_plot)
  })
  
  get_neighbors <- reactive({
    i = kk()
    neighbors = VALUES$RESULTS$ORDERED_NEIGHBORS
    return(neighbors[1:i])
  })
  
  # Identify expts to highlight on the UMAP based on user input k
  expts_to_highlight3 <- reactive({
    if (!is.null(input_rnk3())){
      return(get_neighbors())
    } else{
      return(VALUES$APP_DATA$ALL_ANNOT$rnk_list)
    }
  })
  
  # Plot UMAP and highlight selected expts (if any)
  output$umap_plot3 <- renderPlotly({
    
    # Get the UMAP values and determine which points to highlight
    to_plot = data_for_plotting3()
    expts = expts_to_highlight3()
    to_plot$highlight = ifelse(to_plot$rnk_list %in% expts,
                               'yes', 'no')
    
    if (sum(to_plot$highlight == 'no') == 0){
      p <- plot_ly(to_plot,
                   source = 'C',
                   type = 'scatter', mode = 'markers',
                   x = ~V1,
                   y = ~V2,
                   color = ~tissue,
                   colors = tissue_cols,
                   text = ~paste('ID: ', rnk_list, 
                                 '<br>Tissue: ', tissue,
                                 '<br>Drug: ', drug),
                   marker = list(
                     size = 5,
                     line = list(
                       color = 'black',
                       width = 1
                     )
                   )
      ) 
    } else{
      p <- plot_ly(to_plot,
                   source = 'C',
                   type = 'scatter', mode = 'markers',
                   x = ~V1,
                   y = ~V2,
                   color = ~tissue,
                   colors = tissue_cols,
                   text = ~paste('ID: ', rnk_list, 
                                 '<br>Tissue: ', tissue,
                                 '<br>Drug: ', drug),
                   marker = list(
                     opacity = 0.5,
                     size = 5,
                     line = list(
                       color = 'transparent',
                       width = 1
                     )
                   )
      ) %>%
        add_trace(
          x = to_plot[to_plot$highlight == 'yes', 'V1'],
          y = to_plot[to_plot$highlight == 'yes', 'V2'],
          text = paste('ID: ', to_plot[to_plot$highlight == 'yes', 'rnk_list'],
                       '<br>Tissue: ', to_plot[to_plot$highlight == 'yes', 'tissue'],
                       '<br>Drug: ', to_plot[to_plot$highlight == 'yes', 'drug']),
          name = 'User-filtered',
          marker = list(
            color = tissue_cols[to_plot[to_plot$highlight == 'yes', 'tissue']],
            size = 5,
            line = list(
              color = 'black',
              width = 1
            )
          ),
          inherit = F,
          type = 'scatter', mode = 'markers'
        )
    }
  })
  
  # Get relevant data for barplot
  data_for_barplot3 <- reactive({
    tmp = BARPLOT_DATA
    g = active_gmt3()
    if (g != 'All collections'){
      tmp = tmp[tmp$Collection %in% names(which(GMTS == g)),
                names(tmp) %in% c('Collection', 'Term', expts_to_highlight3())]
    }
    return(tmp)
  })
  
  # Display barplot based on all / selected experiments
  output$barplot3 <- renderPlotly({
    
    # Get the current subset of data to plot in barplot
    all_barplot_data = data_for_barplot3()
    all_barplot_data$Term = all_barplot_data$Collection = NULL
    
    # Calculate % of selected expts with each pathway pos/neg enriched
    n = ncol(all_barplot_data)
    data_for_plot = data.frame(gene_sets = row.names(all_barplot_data),
                               pos = apply(all_barplot_data, 1, function(x){
                                 sum(x > 0, na.rm = T)/n
                               }), 
                               neg = apply(all_barplot_data, 1, function(x){
                                 -1*sum(x < 0, na.rm = T)/n
                               }))
    data_for_plot = data_for_plot[order(data_for_plot$pos - abs(data_for_plot$neg)),]
    data_for_plot$gene_sets = gsub('HALLMARK_', '', data_for_plot$gene_sets)
    data_for_plot$pos = round(data_for_plot$pos, 3)*100
    data_for_plot$neg = round(data_for_plot$neg, 3)*100
    row.names(data_for_plot) = data_for_plot$gene_sets
    
    # If there are more than 50 gene sets, just show the top and bottom 25 so the plot is sane
    if (nrow(data_for_plot) > 50){
      data_for_plot = data_for_plot[c(1:25,
                                      (nrow(data_for_plot)-26):nrow(data_for_plot)),]
    }
    
    # Save top 3 positive and top 3 negative pathways
    tops = c(data_for_plot$pos[length(data_for_plot$pos):(length(data_for_plot$pos)-2)], 
             data_for_plot$neg[1:3])
    names(tops) = data_for_plot$gene_sets[c(length(data_for_plot$pos):(length(data_for_plot$pos)-2), 
                                            1:3)]
    VALUES$RESULTS$TOP_PATHWAYS <- tops
    
    # Optionally, retain original pathway ordering instead of changing based on subset
    if (fix_barplot_order3()){
      data_for_plot = unique(na.omit(data_for_plot[FIXED_BARPLOT_PATHWAYS,]))
    }
    
    # Plot positive bars
    plot_ly(data_for_plot, source = 'B',
            x = ~gene_sets, y = ~pos, type = 'bar', 
            name = 'Positively enriched (q < .05)', 
            marker = list(color = cols[1])) %>%
      # Add negative bars below
      add_bars(y = ~neg, name = 'Negatively enriched (q < .05)', 
               marker = list(color = cols[10])) %>%
      layout(xaxis = list(title = 'Gene set (hover for name)',
                          ticktext = list(''),
                          tickvals = list(1:50),
                          tickmode = 'array',
                          automargin = TRUE), 
             yaxis = list(title = paste0('% experiments (n =', n, ' selected)'),
                          range = c(-100,100)), 
             barmode = 'overlay')
  })
  
  # Format the text that will be displayed under the plots
  output$message_top_pathways3 <- renderText({
    tops = VALUES$RESULTS$TOP_PATHWAYS
    return(paste0('Top positively enriched:\n\t', names(tops)[1], 
                    ' (', round(as.numeric(tops[1]),3), '% of experiments), ',
                    names(tops)[2], 
                    ' (', round(as.numeric(tops[2]),3), '% of experiments), ',
                    names(tops)[3], 
                    ' (', round(as.numeric(tops[3]),3), '% of experiments)',
                    '\nTop negatively enriched:\n\t', names(tops)[4],
                    ' (', -1*round(as.numeric(tops[4]),3), '% of experiments), ',
                    names(tops)[5],
                    ' (', -1*round(as.numeric(tops[5]),3), '% of experiments), ',
                    names(tops)[6],
                    ' (', -1*round(as.numeric(tops[6]),3), '% of experiments)'))
  })
  
  # Title for summary table
  output$title_summary_table3 <- renderText({
    if (!is.null(input_rnk3())){
      return(paste0('Experiment annotations for ', kk(), 
        ' neighbors, in order of increasing distance to uploaded experiment'))
    } else{
      return('Experiment annotations for all experiments')
    }
  })
  
  # Display experiment annotations for selected expts in table
  output$summary_table3 <- renderDataTable({
    rnks = expts_to_highlight3()
    dat = VALUES$APP_DATA$ALL_ANNOT[rnks[rnks != 'UserInput'],]
    VALUES$DATA_FOR_DOWNLOAD$SAMPLE_ANNOTATION_TABLE <- dat
    datatable(dat, 
              selection = 'none',
              rownames = F,
              colnames = c('Comparison ID', 'Cell Line', 'Drug', 'Tissue',
                           'Platform', 'Timepoint', 'Treatment Dose'))
  })
  
  # Display number of experiments selected
  
  # Save selected neighbors
  
  # Saved neighbor message
  
  # Download table of background experiments
  
  # ------------ BOTTOM ROW ---------------------------------------------------
  
  # ----------------- LEFT SIDEBAR
  
  # TODO THESE TWO NEED TO BE UPDATED TO PULL TAB 3 DATA
  # Download table of sample annotations of selected background experiments
  output$download_annotations_selected3 <- downloadHandler(
    filename = function() {
      return('neighbors_sample_annotations.csv')
    },
    content = function(file) {
      write.csv(VALUES$DATA_FOR_DOWNLOAD$SAMPLE_ANNOTATION_TABLE, 
                file, row.names = F)
    }
  )
  
  # Download table of NES values of selected background experiments
  output$download_nes_selected3 <- downloadHandler(
    filename = function() {
      return('GSEA_NES_for_selected_expts.csv')
    },
    content = function(file) {
      write.csv(VALUES$DATA_FOR_DOWNLOAD$GSEA_TABLE, file,
                row.names = F)
    }
  )
})