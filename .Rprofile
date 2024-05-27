dynamicRequire <- function(libs, quietly = TRUE) {
  
  libs_to_install <- character(0)
  
  for (lib in libs) {
    if(!require(lib, character.only = TRUE, quietly = quietly)) {
      libs_to_install <- c(libs_to_install, lib)
    }
  }
  
  if (length(libs_to_install) > 0) {
    utils::install.packages(libs_to_install)
    sapply(libs_to_install, \(lib) {
      suppressPackageStartupMessages(library(
        lib, character.only = TRUE, quietly = quietly
      ))
    })
  }
  
  sapply(libs_to_install, \(lib) {
    suppressPackageStartupMessages(library(
      lib, character.only = TRUE, quietly = quietly
    ))
  })
  
  invisible(libs)
}

libs <- c(
  # Data wrangling
  "tidyverse", "magrittr",
  # Ecology
  "vegan", "rbiom",
  # Differential abundance
  "ALDEx2",
  # Visualisation
  "scales", "patchwork"
)

dynamicRequire(libs)

# Additionally requires BiocManager::install("rhdf5")

if (!require(BiocManager)) {
  install.packages("BiocManager")
}

if (!require(rhdf5)) {
  BiocManager::install('rhdf5')
}

library(rhdf5)

rm(libs)
