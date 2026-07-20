


# ==============================================================================
# R Package Installation and Loading Function (with GitHub Search)
# ==============================================================================

#' Install and Load R Packages from Multiple Sources
#'
#' This function iterates through a vector of package names. If a package is
#' not installed, it attempts to install it from CRAN, then Bioconductor. If
#' still not found, it will search GitHub for a repository with that name and
#' prompt the user to interactively select from the top results for installation.
#'
#' @param packages A character vector or list of package names.
#'   To install directly from GitHub without searching, use the 'username/reponame' format.
#'
#' @return Loads the requested packages into the current session. Prints status
#'   messages and continues processing even if one package fails.
#'
#' @examples
#' # packages <- c("ggplot2", "limma", "tidyverse/dplyr", "hrbrthemes")
#' # install_and_load(packages)
install_and_load <- function(packages) {
  
  # Loop through each package in the provided vector/list
  for (pkg in packages) {
    message(paste("\n--- Processing package:", pkg, "---"))
    
    actual_name <- basename(pkg)
    
    # 1. Check if the package is already installed and load it
    if (requireNamespace(actual_name, quietly = TRUE)) {
      message(paste("✅ Package '", actual_name, "' is already installed. Loading now.", sep = ""))
      library(actual_name, character.only = TRUE)
      next # Skip to the next package in the loop
    }
    
    message(paste("-> Package '", actual_name, "' not found. Attempting installation...", sep = ""))
    is_remote <- grepl("/", pkg)
    
    if (is_remote) {
      # Direct installation from GitHub/GitLab
      tryCatch({
        if (!requireNamespace("remotes", quietly = TRUE)) {
          message("--> Installing 'remotes' to handle remote installation...")
          install.packages("remotes", quiet = TRUE, repos = "https://cran.rstudio.com/")
        }
        message(paste("--> Installing '", pkg, "' from remote source...", sep = ""))
        remotes::install_github(pkg, quiet = TRUE)
      }, error = function(e) {
        warning(paste("--> Remote installation failed for '", pkg, "'.", sep = ""))
      })
      
    } else {
      # Installation from repositories (CRAN, Bioconductor, then search GitHub)
      # 2a. Try CRAN
      tryCatch({
        message(paste("--> Attempt 1/3: Installing '", pkg, "' from CRAN...", sep = ""))
        install.packages(pkg, quiet = TRUE, repos = "https://cran.rstudio.com/")
      }, error = function(e) {})
      
      # 2b. Try Bioconductor if not found on CRAN
      if (!requireNamespace(pkg, quietly = TRUE)) {
        tryCatch({
          message(paste("--> Attempt 2/3: Installing '", pkg, "' from Bioconductor...", sep = ""))
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", quiet = TRUE, repos = "https://cran.rstudio.com/")
          }
          # Suppress updates to avoid lengthy builds during script run
          BiocManager::install(pkg, ask = FALSE, update = FALSE)
        }, error = function(e) {})
      }
      
      # 2c. Search GitHub if not found on CRAN or Bioconductor
      if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("--> Attempt 3/3: Searching for '", pkg, "' on GitHub...", sep = ""))
        tryCatch({
          # Ensure 'gh' package is installed for API access
          if (!requireNamespace("gh", quietly = TRUE)) {
            message("--> Installing 'gh' package to search GitHub API...")
            install.packages("gh", quiet = TRUE, repos = "https://cran.rstudio.com/")
          }
          
          # Search GitHub repositories, sorted by stars
          search_query <- paste0(pkg, " in:name language:R")
          results <- gh::gh("/search/repositories", q = search_query, per_page = 5, sort = "stars", order = "desc")
          
          if (results$total_count > 0) {
            repos <- sapply(results$items, function(item) item$full_name)
            choice_msg <- paste("Found repositories for '", pkg, "'. Please choose one to install:", sep="")
            
            # Use menu() for interactive selection in the console
            choice <- utils::menu(c(repos, "None of the above"), title = choice_msg)
            
            if (choice > 0 && choice <= length(repos)) {
              chosen_repo <- repos[choice]
              message(paste("--> User selected:", chosen_repo))
              
              if (!requireNamespace("remotes", quietly = TRUE)) {
                install.packages("remotes", quiet = TRUE, repos = "https://cran.rstudio.com/")
              }
              remotes::install_github(chosen_repo, quiet = TRUE)
            } else {
              message("--> No package selected. Skipping GitHub installation.")
            }
          } else {
            message("--> No R packages found on GitHub matching that name.")
          }
        }, error = function(e) {
          warning(paste("--> GitHub search failed. Error:", e$message))
        })
      }
    }
    
    # --- Final Loading Attempt ---
    if (requireNamespace(actual_name, quietly = TRUE)) {
      message(paste("✅ Successfully installed and loaded '", actual_name, "'.", sep = ""))
      library(actual_name, character.only = TRUE)
    } else {
      warning(paste("❌ Failed to install package '", actual_name, "' from any source.", sep = ""))
    }
  } # End of for loop
}

# ==============================================================================
# Example Usage
# ==============================================================================
# To run these examples, source this file and then call the function in your console.
#
# source("install_and_load.R")
#
# # Example: Install a vector of packages, including one to be found on GitHub
# message("\n--- Testing a vector of packages ---")
#
packages_to_install <- c(
  "ggplot2",                                 # From CRAN
  "limma",                                   # From Bioconductor
  "tidyverse/dplyr",                         # Direct from GitHub
  "hrbrthemes"                               # Will trigger GitHub search
)



