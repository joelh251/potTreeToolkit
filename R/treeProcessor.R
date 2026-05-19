#' Tree Processor Class
#'
#' @description
#' Read and plots trees
#' @importFrom R6 R6Class
#' @importFrom treeio read.beast
#' @importFrom ape as.phylo
#' @import ggtree
#' @import ggimage
#'
#' @export
treeProcessor <- R6Class("treeProcessor",
  public = list(
    data = NULL,
    initialize = function(filename) {
      self$data <- read.beast(filename)
    },
    tree_with_images = function(image_dir) {
      tree <- as.phylo(self$data)

      images <- data.frame(label = tree$tip.label)
      images$path <- paste0(image_dir, '/', images$label, '.png')

      plot <- ggtree(tree, layout="circular") +
        geom_fruit(
          data     = images,
          geom     = geom_image,
          mapping  = aes(y = label, image = path),
          offset   = 0,
          size     = 0.01
        )
      print(plot)
    }
  )
)
