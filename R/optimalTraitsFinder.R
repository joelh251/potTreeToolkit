library(phangorn)
library(ggtree)
library(R6)
library(treeio)
library(phytools)
library(stringr)
library(ggplot2)
#' Optimal Traits Finder Class
#'
#' @description
#' Finds the optimal no. traits for a phylogeny
#' @importFrom R6 R6Class
#' @importFrom treeio read.beast
#' @importFrom phangorn RF.dist
#' @import ggplot2
#' @importFrom ape as.phylo
#'
#' @export
optimalTraitsFinder <- R6Class("optimalTraitsFinder",
  public = list(
    trees = NULL,
    RF_distances = NULL,
    initialize = function(tree_dir) {
      self$trees <- list()
      filenames = list.files(tree_dir, pattern="*.tree", full.names = TRUE)
      for (filename in filenames) {
        file <- read.beast(filename)
        tree <- as.phylo(file)
        PCs <- as.character(as.numeric(str_extract(basename(filename), "(?<=_|^)\\d+(?=_|\\.|$)")) * 5)
        self$trees[[PCs]] <- tree
      }
      self$trees <- self$trees[order(as.numeric(names(self$trees)))]
    },
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
    generate_RF_heatmap = function() {
      plot <- ggplot(self$RF_distances, aes(x=X, y=Y, fill = RF)) +
        geom_tile() +
        theme_bw()

      print(plot)
    }
  )
)

thing <- optimalTraitsFinder$new("C:/users/joelh/Documents/GitHub/pots-project/optimal_traits/trees")
thing$calculate_RF_distances()
thing$generate_RF_heatmap()
thing$RF_distances
thing$trees
length(thing$trees)
