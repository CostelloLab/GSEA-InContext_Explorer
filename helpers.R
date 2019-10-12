library(RColorBrewer)
library(plyr)
library(FNN)
library(fgsea)
library(umap)

cols = brewer.pal(11, 'Spectral')

# Master table of all experiment annotations
ALL_ANNOT <- get_app_data_dropbox('All_expt_annotations.rds')
ALL_DRUGS <- sort(unique(ALL_ANNOT$drug))
ALL_CELLS <- sort(unique(ALL_ANNOT$cell_line))
ALL_TISSUES <- sort(unique(ALL_ANNOT$tissue))

# Read barplot data and original barplot ordering
BARPLOT_DATA <- get_app_data_dropbox('All_nes_q05.rds')
FIXED_BARPLOT_PATHWAYS <- get_app_data_dropbox('fixed_barplot_pathways.rds')

# Possible gene set collections
GMTS <- c('MSigDB C1: positional gene sets', 
          'MSigDB C2: curated pathways (BioCarta)',
          'MSigDB C2: curated pathways (perturbations)', 
          'MSigDB C2: curated pathways (KEGG)',
          'MSigDB C2: curated pathways (Reactome)', 
          'MSigDB C3: motifs (microRNA targets)',
          'MSigDB C3: motifs (TF targets)', 
          'MSigDB C4: computational (cancer gene neighborhoods)',
          'MSigDB C4: computational (cancer modules)',
          'MSigDB C5: GO gene sets (biological process)',
          'MSigDB C5: GO gene sets (cellular component)',
          'MSigDB C5: GO gene sets (molecular function)',
          'MSigDB C6: oncogenic signatures',
          'DSigDB D1: approved drugs',
          'MSigDB: Hallmarks')
names(GMTS) = sort(unique(BARPLOT_DATA$Collection))

GMT_PATHS <- sort(list_app_data_dropbox('gene_sets', full.names = T))
names(GMT_PATHS) = c('All collections', as.character(GMTS))

all_umap <- get_app_data_dropbox('NES_UMAPs/allCollections_UMAP.rds')
UMAP_PATHS <- sort(list_app_data_dropbox('NES_UMAPs', full.names = T))
names(UMAP_PATHS) = c('All collections', as.character(GMTS))

# Colors for tissue UMAP
tissue_cols = c(brewer.pal(11, 'Spectral'),
                brewer.pal(10, 'Paired'),
                brewer.pal(5, 'Set2'),
                'yellow')
names(tissue_cols) = c(ALL_TISSUES, 'UserInput')

# Colors for GSEA Preranked vs InContext plot
quadrant_cols = c('limegreen', 'goldenrod', 'indianred1', 'grey')
names(quadrant_cols) = c('Both sig', 'Preranked only', 'InContext only', 'Not sig')

# Data set and platform mapping
dataset_platforms = c('HGU133 Plus 2.0',
                      'HGU133A',
                      'HGU133A 2.0')
names(dataset_platforms) = c('GEO 442 (HGU133 Plus 2.0)',
                             'CMap Build 01 (HGU133A)',
                             'NCI-60 (HGU133A 2.0)')

# NES and QVals
#gsea_nes = read.csv('data/Original_GSEAPreranked_hallmarks_NES.csv', 
#                    stringsAsFactors = F)
#gsea_qval = read.csv('data/Original_GSEAPreranked_hallmarks_qvals_ALL.csv',
#                     stringsAsFactors = F)
#incontext_nes = read.csv('data/Modified_GSEAPreranked_hallmarks_NES.csv',
#                         stringsAsFactors = F)
#incontext_qval = read.csv('data/Modified_GSEAPreranked_hallmarks_qvals_ALL.csv',
#                          stringsAsFactors = F)

makeValidRnk <- function(filepath){
  if (is.null(filepath)){
    results = list(err = '<- Please select a .rnk file to upload',
                   df = NULL)
  } else {
    dat = read.table(filepath, stringsAsFactors = F, sep = '\t',
                     col.names = c('GeneSymbol', 'logFC'))
    if (is.character(dat$GeneSymbol) & is.numeric(dat$logFC)){
      results = list(err = '',
                     df = na.omit(dat[order(dat$logFC, decreasing = T),]))
    } else{
      results = list(err = 'Input rnk does not match the format requirements',
                     df = NULL)
    }
    return(results)
  }
}

prepUploadedData <- function(rnk, bg){
  # rnk = user's validated 2-column ranked list
  # bg = is a list with 'rank_indices' (genes x expts) and 'annot' (expts x metadata)
  
  # Subset the user's data and bg data to genes in common
  user_genes = row.names(rnk) = rnk$GeneSymbol
  bg_genes = row.names(bg$rank_indices)
  genes_in_common = sort(intersect(user_genes, bg_genes))
  subset_rnk = rnk[genes_in_common,]
  subset_rnk = subset_rnk[order(subset_rnk$logFC, decreasing = T),]
  subset_bg = bg$rank_indices[genes_in_common,]
  
  # Rank user's genes
  user_ranks = sapply(genes_in_common, function(x){
    which(row.names(subset_rnk) == x)
  })
  
  # Rerank bg genes - this is currently very slow...
  bg_ranks = apply(subset_bg, 2, function(x){
    order(order(x))
  })
  
  # Remove duplicates
  # (Putting user input first here guarantees it stays, 
  # even if it's the same as a rnk list in the bg set)
  mat_for_umap = cbind(UserInput = user_ranks, bg_ranks)
  mat_for_umap2 = mat_for_umap[,!duplicated(t(mat_for_umap))]
  return(mat_for_umap2)
}
