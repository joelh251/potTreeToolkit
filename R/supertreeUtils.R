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
#' @import phangorn
#' @import ggimage
#' @importFrom ggtreeExtra geom_fruit
#' @import ggtree
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
    #' @field constituent_tree_posteriors Named list of posterior probabilities of tres making up supertree
    constituent_tree_posteriors = NULL,
    #' @field supertree Supertree phylo object
    supertree = NULL,
    #' @field ancestral_supertree Supertree with location ancestral state reconstruction
    ancestral_supertree = NULL,

    #' @description
    #' Initialises supertreeUtils Class
    #' @param tree_dir Location of consituent trees of supertree
    #' @param supertree_path Path to supertree .tre file
    #' @return
    #' constituent_trees
    #'
    #' constituent_tree_files
    initialize = function(tree_dir, supertree_path) {
      self$constituent_trees <- list()
      self$constituent_tree_files <- list()

      supertree <- read.tree(supertree_path)
      supertree$tip.label <- gsub("'", "", supertree$tip.label)
      self$supertree <- as.phylo(supertree)

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
    },

    #' @description
    #' Creates a named list of posterior probabilities of each tree
    #' @return
    #' constituent_tree_posteriors
    get_posteriors = function() {
      posteriors <- list()
      for (name in names(self$constituent_tree_files)) {
        tree_data <- as_tibble(self$constituent_tree_files[[name]])
        posteriors[[name]] <- tree_data$posterior

      self$constituent_tree_posteriors <- posteriors
      }
    },

    #' @description
    #' Performs ancestral state reconstruction of location
    #' @param generate_plot Logical, determines whether tree is plotted
    #' @param save_path Save location of plot (must be .pdf)
    #' @param image_dir Location of images for tree tip labels
    #' @return
    #' ancestral_supertree
    ancestral_location_reconstruction = function(generate_plot=FALSE, save_path="supertree_anc.pdf", image_dir="") {

      # 0 = italiote, 1 = greek
      location_data <- setNames(
        ifelse(grepl("italiote", self$supertree$tip.label, ignore.case = TRUE), "Italiote", "Greek"),
        self$supertree$tip.label
      )


      trait_phydat <- phyDat(
        as.matrix(location_data),
        type = "USER",
        levels = levels(as.factor(location_data))
      )


      anc.pars <- anc_pars(self$supertree, trait_phydat, type = "ACCTRAN")

      self$ancestral_supertree <- self$supertree

      unlisted_states <- unlist(anc.pars$state)
      self$ancestral_supertree$node.label <- ifelse(unlisted_states == 2, "Italiote", "Greek")

      if (generate_plot) {

        images <- data.frame(label = self$ancestral_supertree$tip.label)
        images$path <- paste0(image_dir, '/', images$label, '.png')

        plot <- ggtree(self$ancestral_supertree, layout = "circular") +
          geom_nodepoint(aes(color = label),
                         size    = 1.5,
                         na.rm = TRUE
          ) +
          ggtreeExtra::geom_fruit(
            data     = images,
            geom     = geom_image,
            mapping  = aes(y = label, image = path),
            offset   = 0.02,
            size     = 0.008
          ) +
          theme(legend.position = "none")

        ggsave(
          plot = plot,
          filename = save_path,
          bg = "white",
          width = 10,
          height = 10
        )
      }


    }
  )
)
