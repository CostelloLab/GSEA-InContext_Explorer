library(shiny)
library(shinythemes)
library(shinycssloaders)
library(plotly)
library(data.table)
library(DT)
source('dropbox.R')
source('helpers.R')
# Order or the sourcing is 1) ui imports, 2) helpers, 3) rest of ui, 4) server

ui <- fluidPage(
  
  theme = shinytheme('paper'),
  title = 'GSEA-InContext Explorer',
  br(),

  # -------------- TOP ROW -----------------------------------------------
  fluidRow(
    column(3,
           h3('GSEA-InContext Explorer'),
           br(),
           
           # -------------- LEFT SIDEBAR ----------------------------
           wellPanel(
             
             # ------------------- TAB 1 -------------
             conditionalPanel(
               'input.tab === "User Guide"',
               h5('Contact'),
               HTML('Rani Powers and James Costello<br>University of Colorado Anschutz Medical Campus<br>james.costello@cuanschutz.edu<br>'),
               a('Costello Lab website', href='https://sites.google.com/site/jamesccostello4/')
             ),
             
             # ------------------- TAB 2 -------------
             conditionalPanel(
               'input.tab === "Define background set"',
               
               # for resetting this panel with server.button_clear
               shinyjs::useShinyjs(),
               id = 'tab1_side_panel',
               #extendShinyjs(text = "shinyjs.resetSelected = function() { Shiny.onInputChange('.clientValue-plotly_selected-A', 'null'); }",
               #              functions = c("resetSelected")),
               
               h5('Filter experiments'),
               selectInput('filter_platform', label = 'Data set',
                           selected = 'All',
                           choices = c('All', 
                                       'Powers 442 (HGU133 Plus 2.0)',
                                       'CMap Build 01 (HGU133A)',
                                       'NCI-60 (HGU133A 2.0)')),
               selectizeInput('filter_drug', label = 'Drugs (leave blank to include all)',
                              multiple = T,
                              choices = ALL_DRUGS),
               selectizeInput('filter_cell', label = 'Cell Lines (leave blank to include all)',
                              multiple = T,
                              choices = ALL_CELLS),
               selectizeInput('filter_tissue', label = 'Tissues (leave blank to include all)',
                              multiple = T,
                              choices = ALL_TISSUES),
               hr(),
               h5('Set barplot parameters'),
               selectInput('filter_gmt', label = 'Gene set collection',
                           selected = 'MSigDB: Hallmarks',
                           choices = c(as.character(GMTS), 'All collections')),
               checkboxInput('fix_barplot', 'Fix barplot ordering based on full data set'),
               helpText('Note: click and drag on the UMAP plot to ignore the above filters and manually select points instead'),
               actionButton('button_clear', 'Clear all filters', icon = icon('refresh')),
               div(actionLink('info_filters', 'Need help?', icon('info-circle')),
                   style = 'margin-top: 8px;')
             ),
             
             # ------------------- TAB 3 -------------
             conditionalPanel(
               'input.tab === "Run GSEA-InContext"',
               h5('Begin analysis'),
               radioButtons('select_new', 'Start new analysis or look-up existing?',
                            choices = c('New', 'Existing'), selected = 'New'),
               uiOutput('analysis_ui'),
               div(actionLink('info_run', 'Need help?', icon('info-circle')),
                   style = 'margin-top: 8px;')
             ),
             
             # ------------------- TAB 4 -------------
             conditionalPanel(
               'input.tab === "Explore pathways"',
               h5('Upload ranked list'),
               helpText('Columns: GeneSymbol, logFC'),
               fileInput('input_file3', 'Upload'),
               h5('Other parameters'),
               selectInput('filter_gmt3', label = 'Gene set collection',
                           selected = 'All collections',
                           choices = c('All collections',
                                       as.character(GMTS))),
               numericInput('k', 'Select k neighbors', value = 30, min = 1, max = 5717),
               checkboxInput('fix_barplot3', 'Fix barplot ordering based on full data set'),
               actionLink('info_k', 'Need help?', icon('info-circle'))
             )
          ),
             style="overflow-x: scroll; overflow-y: scroll"),
    
    # -------------- RIGHT / MAIN PANEL ----------------------------
    column(9,
           tabsetPanel(
             id = 'tab',
             
             # ------------------- TAB 1 -------------
             tabPanel('User Guide',
                      h5('Welcome!'),
                      strong('About the GSEA-InContext Explorer tool'),
                      HTML('<br>This tool allow users to explore a compendium of <i>in vitro</i> gene expression experiments and run the GSEA-InContext algorithm.'),
                      h5('Instructions'),
                      strong('1. Define background set'),
                      br(),
                      'Use the filters on the "Define background set tab" to select a relevant subset of experiments. As filter criteria are added, the barplot on the right side will automatically update to show the positively and negatively enriched gene sets in the selected experiments. Once you are satisfied with the selected set of background experiments, click the "Save Background Set" button to save the selected experiments for use in the GSEA-InContext algorithm.',
                      br(),
                      strong('2. Run GSEA-InContext'),
                      br(),
                      'Upload your differential expression results as a .rnk format file. Confirm that the selected background experiments are correct, and click the "Run GSEA-InContext" button to begin analysis.',
                      br(),
                      strong('3. Explore Pathway similarity'),
                      br(),
                      'On the "Explore pathway similarity" tab, you can upload an experiment in a .rnk format file and visualize how it compares to other experiments in the compendium.'
             ),
             
             # ------------------- TAB 2 -------------
             tabPanel('Define background set',
                      br(),
                      fluidRow(
                        column(6,
                               h5(textOutput('umap_plot_title')),
                               'Each point represents one differential experiment. UMAP is computed on the ranked lists for each experiment.',
                               withSpinner(plotlyOutput('umap_plot'))),
                        column(6,
                               h5(textOutput('title_barplot')),
                               'Barplot of 50 most positively and negatively enriched gene sets in the collection',
                               withSpinner(plotlyOutput('barplot'))
                          )
                      ),
                      br(),
                      verbatimTextOutput('message_top_pathways')
             ),
             
             # ------------------- TAB 3 -------------
             tabPanel('Run GSEA-InContext',
                      br(),
                      fluidRow(
                        column(5,
                               h5('Enrichment plot'),
                               h6(span(textOutput('message_run_incontext'), style='color:2196F3')),
                               withSpinner(plotOutput('fgsea_plot'))),
                        column(7,
                               h5('GSEA Preranked vs GSEA-InContext q-values'),
                               withSpinner(plotlyOutput('qq_plot')))
                      )
             ),
             
             # ------------------- TAB 4 -------------
             tabPanel('Explore pathways',
                      br(),
                      fluidRow(
                        column(6,
                               h5(textOutput('title_k')),
                               'Each point represents one differential experiment. UMAP is computed on the NES values for each experiment, for the selected gene set collection.',
                               withSpinner(plotlyOutput('umap_plot3'))),
                        column(6,
                               h5(textOutput('title_barplot3')),
                               'Barplot of 50 most positively and negatively enriched gene sets in the collection',
                               withSpinner(plotlyOutput('barplot3')))
                      ),
                      br(),
                      verbatimTextOutput('message_top_pathways3')
                      #verbatimTextOutput('message_umap_genes')
             )
         )
      )
  ),
  
  # -------------- BOTTOM ROW --------------------------------------------
  fluidRow(
    column(3,
           
           # -------------- LEFT SIDEBAR ----------------------------
           wellPanel(
             
             # ------------------- TAB 1 -------------
             conditionalPanel(
               'input.tab === "User Guide"',
               h5('References'),
               'Powers, et al. GSEA-InContext: identifying novel and common patterns in expression experiments. Bioinformatics, Volume 34, Issue 13, 1 July 2018, Pages i555–i564'
             ),
             
             # ------------------- TAB 2 -------------
             conditionalPanel(
               'input.tab === "Define background set"',
               h5(span(textOutput('message_n_rnks_selected'), style='color:#2196F3')),
               helpText('Add additional filter above, if desired, then click below to save these experiments as your background set'),
               actionButton('button_save', 'Save background set', icon = icon('check-square')),
               span(textOutput('message_saved_rnks_selected'), style='color:#2196F3'),
               helpText('Click to download table shown to the right'),
               downloadButton('download_rnks_selected', 'Download annotation file'),
               helpText('Click to download GSEA NES values for selected experiments'),
               downloadButton('download_nes_selected', 'Download GSEA results')
             ),
             
             # ------------------- TAB 3 -------------
             conditionalPanel(
               'input.tab === "Run GSEA-InContext"',
               downloadButton('download_gseaincontext', 'Download GSEA-InContext results'),
               helpText('Download the table shown to the right')
             ),
             
             # ------------------- TAB 4 -------------
             conditionalPanel(
               'input.tab === "Explore pathways"',
             #  h5(span(textOutput('message_neighbors_selected'), style="color:#2196F3")),
               helpText('Adjust k above, if desired, then click below to save the annotations for the highlighted neighboring experiments.'),
               #actionButton('button_save2', 'Save background set', icon = icon('check-square')),
               #span(textOutput('message_saved_neighbors_selected'), style='color:#2196F3'),
               br(),
               downloadButton('download_annotations_selected3', 'Download annotation file')
             )
        )
    ),
    
    # -------------- RIGHT / MAIN PANEL ----------------------------
    column(9,
           
           # ------------------- TAB 1 -------------
           conditionalPanel(
             'input.tab === "Define background set"',
             h5('Annotations for selected experiments'),
             withSpinner(DT::dataTableOutput('summary_table'))
           ),
           
           # ------------------- TAB 2 -------------
           conditionalPanel(
             'input.tab === "Run GSEA-InContext"',
             h5('All GSEA Preranked and GSEA-InContext results'),
             withSpinner(dataTableOutput('gsea_incontext_table'))
           ),
           
           # ------------------- TAB 3 -------------
           conditionalPanel(
             'input.tab === "Explore pathways"',
             h5(textOutput('title_summary_table3')),
             withSpinner(dataTableOutput('summary_table3'))
           )
      )
    )
)