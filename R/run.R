#' It runs the model (as time series prediction) in MATLAB
#'
#' @description
#' `run` opens the MATLAB, run the model and close it when is done
#'
#' @details
#' This is a function run the model and MATLAB using the input variables set before
#' **important** it needs MATLAB installed in the machine and the functions
#' @seealso [setup(), set_static_inputs(), set_ts_inputs()] be previously ran
#' with the same site_name and run_name
#'
#' @param patch the patch to the STEMMUS_SCOPE model directory, default "D:/model/STEMMUS_SCOPE/",
#' @param site_name,run_name the name of the location and the name of run. The last can be used to name runs with different model parameters or settings,
#' @param cores set 1,2,4,6,8 cores for ech run **experimental**
#' @return run the model in MATLAB that will return the resuls in the folder ~/output/site_name_run_name_time/
#'
#' @examples
#' \dontrun{
#'
#' n = 10
#' for (i in 1:10) {
#'
#' run_Matlab(patch = "D:/model/STEMMUS_SCOPE/",
#'            site_name = site_names[i],
#'            run_name = "DWD_36")
#'            svMisc::progress(i, n, progress.bar = TRUE, init = T) # time delay
#'            Sys.sleep(1800) # time delay in seconds 7200/60
#'            if (i == n) message("Done!")
#' }
#'
#' }
#'
#' @export

run <- function(patch = "D:/model/STEMMUS_SCOPE/",
                site_name = NA,
                run_name = NA,
                cores = NULL

){

# change the config file ----
config_file <- utils::read.table(, file = paste0(patch, "input/runs/", site_name,  "_", run_name, "/", site_name,"_", run_name, "_", "config.txt"))

config_file$V1[[12]] <- paste0("OutputPath=",patch,"output/", site_name, "_", run_name, "_", format(Sys.time(), "%Y%b%d_%H%M"), "/")

utils::write.table(config_file, file = paste0(patch, "input/runs/", site_name, "_", run_name, "/", site_name,"_", run_name, "_", "config.txt"),
            quote = FALSE, col.names = FALSE, row.names = FALSE)

# create MATLAB path to find the new input directory ----
CFG = paste0(patch, "input/runs/", site_name, "_", run_name, "/", site_name,"_", run_name, "_", "config.txt")

# write the one line text file into the run directory ----
# (will be overwrite every time the run_setup is run)

utils::write.table(CFG, file = paste0(patch, "input/runs/", "path.txt"),
            sep = " ", col.names = F, row.names = F, quote = F, eol = "", append=F)

if(missing(cores)){

# run Matlab simulation and close it after finishing ----

system('matlab -useStartupFolderPref -r "STEMMUS_SCOPE; exit"')

} else {

  if(cores == 1) {
    system('matlab -useStartupFolderPref -r "parpool(1); STEMMUS_SCOPE; exit"')
  }

  if(cores == 2) {
    system('matlab -useStartupFolderPref -r "parpool(2); STEMMUS_SCOPE; exit"')
  }

  if(cores == 4) {
    system('matlab -useStartupFolderPref -r "parpool(4); STEMMUS_SCOPE; exit"')
  }

  if(cores == 6) {
    system('matlab -useStartupFolderPref -r "parpool(6); STEMMUS_SCOPE; exit"')
  }

  if(cores == 8) {
    system('matlab -useStartupFolderPref -r "parpool(8); STEMMUS_SCOPE; exit"')
  }

}

return(print("A Matlab window will open and close automatically after run"))
}

