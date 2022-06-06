#A script to be sourced at the beginning of each working scripts to load common libraries, working paths, and helpful general lookup tables/datasets 

#working paths
sys_path <- ifelse(Sys.info()["sysname"]=="Windows", "G:/Shared drives/",ifelse(Sys.info()["sysname"]=="Linux", "/home/jason/Documents/Gdrive_sync/emlab_shared/", "~/Google Drive/Shared drives/"))
# Path to our emLab's data folder
emlab_data_path <- paste0(sys_path,"emlab/data")
# Path to this project's folder
project_path <- paste0(sys_path,"emlab/projects/current-projects/blue-prosperity-coalition")
#Path to project's data folder
project_data_path <- file.path(project_path,"broader-research/mpa-bycatch/data")
# path to project figures
project_figure_path <- file.path(project_path, "broader-research", "mpa-bycatch", "figures")


# Edit working paths if it's a non-emLabber
# setup should be: 
# |- Google Drive
#   |- My Drive
#     |- mpa-bycatch
#       |- data
#       |- figures
if(!dir.exists(sys_path)) { 
  sys_path <- ifelse(Sys.info()["sysname"]=="Windows", "G:/My Drive/", # windows
                     "~/Google Drive/My Drive/") # mac
  emlab_data_path <- NULL
  project_path <- paste0(sys_path, "mpa-bycatch")
  project_data_path <- file.path(project_path, "data")
  project_figure_path <- file.path(project_path, "figures")
  }
