#' Tree Processor Class
#'
#' @description
#' Read and plots trees
#' @importFrom R6 R6Class
#' @importFrom treeio read.beast
#' @import ape
#' @import ggtree
#' @import ggimage
#' @importFrom ggtreeExtra geom_fruit
#' @import dplyr
#'
#' @export
treeProcessor <- R6Class("treeProcessor",
  public = list(
    #' @field data Beast formatted tree
    data = NULL,

    #' @description
    #' Initialises the treeProcessor class
    #' @param filename Path to tree file
    initialize = function(filename) {
      self$data <- read.beast(filename)
    },

    #' @description
    #' Plots a tree with images as tip labels
    #' @param image_dir Directory to find pot images
    #' @param plotname Filename for output
    #' @return
    #' Tree with images as tip labels
    tree_with_images = function(image_dir, plotname = "tree.pdf") {
      tree <- as.phylo(self$data)

      images <- data.frame(label = tree$tip.label)
      images$path <- paste0(image_dir, '/', images$label, '.png')

      plot <- ggtree(tree, layout="rectangular") +
        geom_fruit(
          data     = images,
          geom     = geom_image,
          mapping  = aes(y = label, image = path),
          offset   = 0,
          size     = 0.01
        )

      ggsave(
        filename = plotname,
        bg = "white",
        width = 10,
        height = 10
      )
    },

    #' @description
    #' Plots a tree with images as tip labels
    #'
    #' Collapses poorly supported nodes into multichotomies
    #' @param image_dir Directory to find pot images
    #' @param posterior_threshold Threshold for node support (default = 0.5)
    #' @param plotname Filename for output
    #' @return
    #' Tree with images as tip labels
    collapsed_posterior_tree = function(image_dir,
                                        posterior_threshold=0.5,
                                        plotname="collapsed_tree.pdf") {
      tree <- as_tibble(self$data)
      tree2 <- as.phylo(tree)
      tips <- length(tree2$tip.label)

      images <- data.frame(label = tree2$tip.label)
      images$path <- paste0(image_dir, '/', images$label, '.png')

      tree2$node.label <- tree %>%
        filter(node > tips) %>%
        arrange(node) %>%
        pull(posterior)

      tree3 <- tree2
      poor_nodes <- which(as.numeric(tree2$node.label) < posterior_threshold) + tips
      poor_edges <- which(tree3$edge[, 2] %in% poor_nodes)
      tree3$edge.length[poor_edges] <- 0
      tree3 <- di2multi(tree3, tol = 1e-8)

      plot <- ggtree(tree3, layout="rectangular") +
        geom_fruit(
          data     = images,
          geom     = geom_image,
          mapping  = aes(y = label, image = path),
          offset   = 0,
          size     = 0.01
        )

      ggsave(
        filename = plotname,
        plot = plot,
        bg = "white",
        width = 10,
        height = 10
      )
    }
  )
)
