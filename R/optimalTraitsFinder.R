#' Optimal Traits Finder Class
#'
#' @description
#' Finds the optimal no. traits for a phylogeny
#' @importFrom R6 R6Class
#' @importFrom treeio read.beast
#' @importFrom phangorn RF.dist
#' @import ggplot2
#' @import stringr
#' @importFrom ape as.phylo
#'
#' @export
optimalTraitsFinder <- R6Class("optimalTraitsFinder",
  public = list(
    #' @field trees Named list of trees as phylo objects
    trees = NULL,
    #' @field tree_files Named list of trees as beast objects
    tree_files = NULL,
    #' @field RF_distances Dataframe of pairwise RF distances
    RF_distances = NULL,
    #' @field mean_posteriors Dataframe of mean posterior probabilities for each tree
    mean_posteriors = NULL,

    #' @description
    #' Initialises optimalTraitsFinder class
    #' @param tree_dir Location of tree files
    #' @return
    #' trees
    #'
    #' tree_files
    initialize = function(tree_dir) {
      self$trees <- list()
      self$tree_files <- list()
      filenames = list.files(tree_dir, pattern="*.tree", full.names = TRUE)
      for (filename in filenames) {
        file <- read.beast(filename)
        tree <- as.phylo(file)
        PCs <- as.character(as.numeric(str_extract(basename(filename), "(?<=_|^)\\d+(?=_|\\.|$)")) * 5)
        self$trees[[PCs]] <- tree
        self$tree_files[[PCs]] <- file
      }
      self$trees <- self$trees[order(as.numeric(names(self$trees)))]
    },

    #' @description
    #' Calculates pairwise normalised RF distances for input trees
    #' @return
    #' RF_distances
    calculate_RF_distances = function() {
      self$RF_distances <- expand.grid(X=names(self$trees),
                                       Y=names(self$trees))

      self$RF_distances$RF <- sapply(1:nrow(self$RF_distances), function(i) {
        RF.dist(self$trees[[self$RF_distances$X[i]]],
                self$trees[[self$RF_distances$Y[i]]],
                normalize  = TRUE,
                rooted     = TRUE,
                check.labels = TRUE)
      })
    },

    #' @description
    #' Plots a heatmap of RF distances
    #' @return
    #' Heatmap of RF distances
    generate_RF_heatmap = function() {
      plot <- ggplot(self$RF_distances, aes(x=X, y=Y, fill = RF)) +
        geom_tile() +
        theme_bw() +
        labs(
          x = "PCs",
          y = "PCs"
        )

      print(plot)
    },

    #' @description
    #' Calculates mean posterior probability of nodes in each tree
    #' @return
    #' mean_posteriors
    mean_posterior_support = function() {
      self$mean_posteriors <- data.frame(PCs = as.numeric(names(self$tree_files)),
                                         mean_posterior = NA)

      for (PC in self$mean_posteriors$PCs) {
        tree_data <- as_tibble(self$tree_files[[as.character(PC)]])
        posteriors <- tree_data$posterior
        mean_posterior <- mean(posteriors, na.rm = TRUE)
        self$mean_posteriors$mean_posterior[self$mean_posteriors$PCs == PC] <- mean_posterior
      }
    }
  )
)
