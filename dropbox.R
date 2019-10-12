# Functions for reading and writing to Dropbox
library(rdrop2)
#token <- drop_auth()
#saveRDS(token, "tokenfile.RDS")
drop_auth(rdstoken = "tokenfile.RDS")

dropbox_root_read <- 'Shiny_GSEA-InContext_Explorer/app_data/'
dropbox_root_write <- 'Shiny_GSEA-InContext_Explorer/user_data/'

if (!dir.exists('tmp')){
  dir.create('tmp')
}

list_app_data_dropbox <- function(path, full.names = F){
  tab = as.data.frame(drop_dir(paste0(dropbox_root_read, path)))
  if (full.names){
    tab$name = paste0(path, '/', tab$name)
  }
  return(tab$name)
}

get_app_data_dropbox <- function(f){
  ff = strsplit(f, '/')[[1]]
  ff = ff[length(ff)]
  ext = strsplit(ff, '')[[1]]
  ext = paste0(ext[(length(ext)-3):(length(ext))], collapse = '')
  print(ff)
  if (!file.exists(paste0('tmp/', ff))){
    success = drop_download(path = paste0(dropbox_root_read, f), 
                            local_path = paste0('tmp/', ff),
                            progress = F)
  }
  if (ext == '.rds'){
    dat = readRDS(paste0('tmp/', ff))
  } else if (ext == '.gmt'){
    dat = gmtPathways(paste0('tmp/', ff))
  }
  return(dat)
}

read_UMAP_dropbox <- function(dataset_name){
  if (dataset_name == 'Powers 442 (HGU133 Plus 2.0)'){
    f = 'Orig442_ranks_UMAP.rds'
  } else if (dataset_name == 'CMap Build 01 (HGU133A)'){
    f = 'CMap01_ranks_UMAP.rds'
  } else if (dataset_name == 'NCI-60 (HGU133A 2.0)'){
    f = 'NCI-60_ranks_UMAP.rds'
  } else {
    f = 'intersected_ranks_UMAP.rds'
  }
  if (!file.exists(paste0('tmp/', f))){
    drop_download(path = paste0(dropbox_root_read, 'rank_UMAPs/', f), 
                  local_path = paste0('tmp/', f),
                  progress = F)
  }
  dat = readRDS(paste0('tmp/', f))
  return(dat)
}

initialize_analysis_dropbox <- function(bg_expts){
  id = paste0(sample(c(LETTERS, 0:9), 6, replace = T), collapse = '')
  # TODO check that id doesn't already exist in Dropbox
  # TODO create folder in Dropbox
  dropbox_path = paste0(dropbox_root, 'user_data/', id, '/')
  # Save csv of bg_expts to the folder
  write.csv(bg_expts, 
            paste0(dropbox_path, '/bg_expts.csv'), 
            col.names = 'bg_expts', row.names = F)
  return(id)
}

save_analysis_dropbox <- function(id, user_rnk, results_table){
  dropbox_path = paste0(dropbox_root, 'user_data/', id, '/')
  # Save user_rnk
  # Save results_table
  
}

