
###################### 
## Create folder "./req_lib/"
## that stores package installations.
###################### 
wd <- getwd()
subfolder_new <- paste0("req_lib/")

if (!dir.exists(subfolder_new)) { ## IF directory doesnt exist: create it.
       dir.create(subfolder_new)
}



###################### 
## Installation through CRAN: example
###################### 
if (!require(mvtnorm, lib = req_lib_dir)) { ## IF installed: load/IF NOT installed:
    .libPaths(req_lib_dir)                  ## Set library path to local directory "./req_lib/".
    install.packages("mvtnorm",             ## Install from CRAN to "./req_lib/"
                     lib = req_lib_dir,
                     repos = "https://archive.linux.duke.edu/cran/")
    library(mvtnorm, lib = req_lib_dir)     ## Load package from "./req_lib/"
}

###################### 
## Installation through source code: example
###################### 
if (!require(equalCovs, lib = req_lib_dir)) { ## IF installed: load/IF NOT installed:
    .libPaths(req_lib_dir)                    ## Set library path to local directory "./req_lib/".
    install.packages(                         ## Install from source file tar.gz to "./req_lib/"
        'equalCovs_1.0.tar.gz',               ## Note: make sure tar.gz file is saved in wd.
        repos=NULL, type='source',
        lib = req_lib_dir)
    library(equalCovs, lib = req_lib_dir)     ## Load package from "./req_lib/"
}

###################### 
## Note 1: If this code is ran multiple times,
##        the IFs are such that the packages
##        install only once, and in later runs
##        they only load the package without 
##        reinstalling.
###################### 
## Note 2: Only the packages that are sourced
##        with source files tar.gz have controlled
##        versions. The code for the packages loaded
##        from CRAN do have versioning. 
