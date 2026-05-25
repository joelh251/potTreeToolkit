#' Supertree Utilities Class
#'
#' @description
#' Collection of utilities for building a supertree
#' @importFrom R6 R6Class
#' @importFrom treeio read.beast
#' @import ggplot2
#' @import stringr
#' @importFrom ape as.phylo
#' @import dplyr
#'
#' @export
supertreeUtils <- R6Class("supertreeUtils",
  public = list(
    #' @field constituent_trees Named list of trees making up supertree as phylo objects
    constituent_trees = NULL,
    #' @field constituent_tree_files Named list of trees making up supertree as beast objects
    constituent_tree_files = NULL,
    #' @field mean_posteriors Dataframe of mean posterior probabilities of constituent trees
    mean_posteriors = NULL,

    #' @description
    #' Initialises supertreeUtils Class
    #' @param tree_dir Location of consituent trees of supertree
    #' @return
    #' constituent_trees
    #'
    #' constituent_tree_files
    initialize = function(tree_dir) {
      self$constituent_trees <- list()
      self$constituent_tree_files <- list()

      filenames <- list.files(tree_dir, pattern = "\\.tree$", full.names = TRUE)
      for (filename in filenames) {
        file <- read.beast(filename)
        tree <- as.phylo(file)
        index <- as.character(as.numeric(str_extract(basename(filename), "\\d+")))
        self$constituent_trees[[index]] <- tree
        self$constituent_tree_files[[index]] <- file
      }
      self$constituent_trees <- self$constituent_trees[order(as.numeric(names(self$constituent_trees)))]
      self$constituent_tree_files <- self$constituent_tree_files[order(as.numeric(names(self$constituent_tree_files)))]
    },

    #' @description
    #' Calculates mean posterior probability of nodes in each tree
    #' @return
    #' mean_posteriors
    mean_posterior_support = function() {
      self$mean_posteriors <- data.frame(index = as.numeric(names(self$constituent_tree_files)),
                                         mean_posterior = NA)

      for (i in self$mean_posteriors$index) {
        tree_data <- as_tibble(self$constituent_tree_files[[as.character(i)]])
        posteriors <- tree_data$posterior
        mean_posterior <- mean(posteriors, na.rm = TRUE)
        self$mean_posteriors$mean_posterior[self$mean_posteriors$index == i] <- mean_posterior
      }
    }
  )
)
