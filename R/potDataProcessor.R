#' Pot Data Processor Class
#'
#' @description
#' Process pot data and generate basic plots
#' @importFrom R6 R6Class
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom scales hue_pal
#' @importFrom readr write_tsv
#' @importFrom ape as.phylo
#' @importFrom ape write.nexus.data
#'
#' @export
potDataProcessor <- R6Class("potDataProcessor",
  public = list(
    #' @field data Raw pot data
    data = NULL,
    #' @field character_data Formatted PCA data for each pot, only PCs with variance > 0
    character_data = NULL,
    #' @field taxa_data Min and max ages of each pot
    taxa_data = NULL,
    #' @field PCA_data All PCs
    PCA_data = NULL,
    #' @field scaled_taxa_data Ages scaled relative to youngest taxon
    scaled_taxa_data = NULL,
    #' @field raw_pca_output Raw output of prcomp function
    raw_pca_output = NULL,
    #' @field pair_sample_character_data Subset of character data
    pair_sample_character_data = NULL,
    #' @field pair_sample_taxa_data Subset of taxa data
    pair_sample_taxa_data = NULL,
    #' @field dendrogram Ward tree of pots as a phylo object
    dendrogram = NULL,

    #' @description
    #' Initialises potDataProcessor class
    #' @param data Dataframe of pot data
    initialize = function(data) {
      self$data <- data
    },

    #' @description
    #' Transforms input data into a RevBayes-suitable format.
    #'
    #' Optionally saves data as a .nex and .tsv file.
    #' @param save_dir Directory to save outputs to
    #' @param character_data_filename Filename without extension
    #' @param taxa_data_filename Filename without extension
    #' @return
    #' character_data (optionally saved as a .nex file)
    #'
    #' taxa_data (optionally saved as a .tsv file)
    #'
    #' PCA_data
    #'
    #' raw_pca_output
    generate_revbayes_data = function(save_dir="",
                                      character_data_filename="character_data",
                                      taxa_data_filename="taxa_data") {
      # Do PCA
      shapes <- self$data[, c("type", "shape", "place", "x", "y", "point_order")]
      wide <- shapes%>%
        pivot_wider(id_cols = c(type, shape, place), names_from = point_order, values_from = c(x, y))
      pca <- prcomp(wide[, -c(1, 2, 3)], scale. = TRUE)
      self$raw_pca_output <- pca

      # Format PCA nicer
      scores <- as.data.frame(pca$x)
      scores$id <- wide$type
      scores$shape <- wide$shape
      scores$place <- wide$place
      scores$shape_place <- paste(scores$shape, scores$place)
      self$PCA_data <- scores

      # Generate taxon data
      names <- unique(self$data$type)
      max_age <- sapply(names, function(z) self$data$date_earliest[self$data$type == z][1])
      max_age <- sapply(max_age, private$convertDate)
      min_age <- sapply(names, function(z) self$data$date_latest[self$data$type == z][1])
      min_age <- sapply(min_age, private$convertDate)
      taxa_data <- data.frame(taxon = names, min_age = min_age, max_age = max_age)
      taxa_data <- na.omit(taxa_data)

      # FBD can't handle age range of 0 - must increase interval
      no_age_range <- taxa_data$taxon[(taxa_data$max_age-taxa_data$min_age) == 0]
      if (length(no_age_range) > 0) {
        warning("The following taxa have had their age range expanded +- 25 years due to having an age range of 0:\n",
                paste(no_age_range, collapse = "\n"),
                call. = FALSE)
        taxa_data <- taxa_data %>%
          mutate(min_age = ifelse(taxon %in% no_age_range, min_age - 25, min_age)) %>%
          mutate(max_age = ifelse(taxon %in% no_age_range, max_age + 25, max_age))
      }

      # Only save characters which explain a proportion of variance > 0
      prop_var <- summary(pca)$importance[2,]
      useful_characters <- length(prop_var[prop_var > 0])
      character_data <- scores[,1:useful_characters]
      rownames(character_data) <- as.character(scores$id)

      # Filter for only shared taxa between datasets
      shared_taxa <- intersect(rownames(character_data), taxa_data$taxon)
      character_data <- as.matrix(character_data)
      character_data    <- character_data[rownames(character_data) %in% shared_taxa, ]
      self$taxa_data    <- taxa_data[taxa_data$taxon %in% shared_taxa, ]
      self$character_data <- character_data

      # Save data
      if (save_dir != "") {
        character_data_path <- paste0(save_dir, "/", character_data_filename, ".nex")
        taxa_data_path <- paste0(save_dir, "/", taxa_data_filename, ".tsv")

        write.nexus.data(self$character_data, file = character_data_path, format = "continuous")
        write_tsv(self$taxa_data, taxa_data_path)
      }
    },

    #' @description
    #' Creates a plot of PCA results
    #' @param cols A vector of 2 column names to plot
    #' @param add_elipses Logical, dictates whether elipses are drawn around clusters
    #' @return
    #' 2 axis PCA plot
    PC_plot = function(cols=c("PC1", "PC2"), add_elipses=FALSE){
      # Generate data if not already generated
      if (is.null(self$PCA_data)) {
        self$generate_revbayes_data()
        warning("Input data not processed. Running self$generate_revbayes_data()")
      }

      # Make the plot
      plot <- ggplot(self$PCA_data, aes(x=.data[[cols[1]]], y=.data[[cols[2]]], colour = shape)) +
        geom_point()  +
        theme_bw() +
        coord_fixed(ratio = 1) +
        theme(
          legend.position = "bottom"
        )

      # Optionally add elipses
      if (add_elipses) {
        plot <- plot + stat_ellipse()
      }

      print(plot)
    },

    #' @description
    #' Generates a Ward tree of pots
    #' @param n Number of PCs to use in distance calculations
    #' @return
    #' Dendrogram plot
    #'
    #' dendrogram - phylo object of dendrogram
    plot_dendrogram = function(n=10){
      # Generate data if not already generated
      if (is.null(self$PCA_data)) {
        self$generate_revbayes_data()
        warning("Input data not processed. Running self$generate_revbayes_data()")
      }

      # Make there are enough PCs in the data
      if (n > ncol(self$PCA_data)-4) {
        stop("Supplied number of PCs greater than number of PCs in data")
      }

      # Construct ward tree - does not scale PCs
      cluster <- as.matrix(self$PCA_data[, 1:n])
      rownames(cluster) <- self$PCA_data$id
      d <- dist(cluster, method = "euclidean")
      hc <- hclust(d, method = "ward.D2")
      tree <- as.phylo(hc)

      self$dendrogram <- tree

      # Colours tip labels based on a priori shape assignment
      metadata <- data.frame(
        label = tree$tip.label,
        type = self$PCA_data$shape[match(tree$tip.label, self$PCA_data$id)])

      p <- ggtree(self$dendrogram) %<+% metadata +
        geom_tiplab(aes(colour = type), size = 2)

      print(p)
    },

    #TODO code
    #' @description
    #' UNFINISHED
    #'
    #' Generates an image of a single pot profile
    #' @param data REMOVE
    #' @param name Name of pot
    #' @param save_dir Location to save image (don't include filename)
    #' @return
    #' Image of a single pot
    #'
    #' Optional - saves image to .png file
    generate_pot_image = function(data, name, save_dir = "") {
      #TODO adjust code
      groups      <- sort(unique(data$shape))   # ggplot sorts levels alphabetically
      n           <- length(groups)
      default_cols <- setNames(hue_pal()(n), groups)

      df <- subset(data, type == name)
      plot <- ggplot(df, aes(x=x, y=y, group=type, fill=shape))+
        geom_polygon()+
        xlim(-1.5,1.5)+
        ylim(-1.5, 1.5)+
        coord_fixed(1) +
        theme(aspect_ratio=1)+
        scale_fill_manual(values = default_cols) +
        theme_void() +
        theme(
          panel.background  = element_rect(fill = "transparent", colour = NA),
          plot.background   = element_rect(fill = "transparent", colour = NA),
          legend.position = "none",
          plot.margin = margin(0,0,0,0))

      if (save_dir != "") {
        #save it
        #filename is pot name + .png (maybe allow specification of savetype)
      }
      else {
        print(plot)
      }
    },

    #' @description
    #' Generates images of all pots in data and saves to .png files
    #' @param save_dir Location to save pot images (don't include filenames)
    #' @return
    #' A directory of pot images
    generate_pot_images = function(save_dir) {
      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

      # Create colour scale of pot shapes
      groups      <- sort(unique(self$data$shape))
      n           <- length(groups)
      pot_palette <- setNames(hue_pal()(n), groups)

      # Generate pot images
      self$data %>%
        group_by(type) %>%
        arrange(desc(point_order), .by_group = TRUE) %>%
        group_walk(private$pot_image, save_dir = save_dir, palette = pot_palette)
    },

    #' @description
    #' UNFINISHED
    #'
    #' Creates a figure with all pot profiles
    #' @return
    #' Figure of all pot profiles
    #'
    #' Optional - saves to a .pdf file
    summarise_pot_profiles = function() {
      #TODO code it you dummy
    },

    #' @description
    #' Scales all pot ages relative to the youngest taxon
    #' @param save_dir Location to save output
    #' @param filename Name of output file (don't include extension)
    #' @return
    #' scaled_taxa_data
    #'
    #' Optional - saves output to a .tsv file
    scale_taxa_data = function(save_dir="", filename="scaled_taxa_data") {
      # Generate data if not already generated
      if (is.null(self$taxa_data)) {
        self$generate_revbayes_data()
        warning("Input data not processed. Running self$generate_revbayes_data()")
      }

      self$scaled_taxa_data <- self$taxa_data[self$taxa_data$taxon %in% rownames(self$character_data),]

      minimum_age <- min(self$scaled_taxa_data$min_age)

      self$scaled_taxa_data$min_age <- self$scaled_taxa_data$min_age - minimum_age
      self$scaled_taxa_data$max_age <- self$scaled_taxa_data$max_age - minimum_age

      if (save_dir != "") {
        write_tsv(self$scaled_taxa_data, paste0(save_dir, "/", filename, ".tsv"))
      }
    },

    #' @description
    #' Uses HCA to group pots into pairs by greates shape similarity
    #' @param n_PCs Number of PCs to feed into HCA (default = 15)
    #' @return
    #' Pair column in data with pair identity of all pots
    make_pot_pairs = function(n_PCs=15) {
      # Generate data if not already generated
      if (is.null(self$PCA_data)) {
        self$generate_revbayes_data()
        warning("Input data not processed. Running self$generate_revbayes_data()")
      }

      self$data$pair <- NA

      self_ref <- self
      private$cluster_pairs(.data = self$data, n_PCs = n_PCs, self_ref = self_ref)

      message("Generated ", length(unique(self$data$pair)), " pairs")

    },

    #' @description
    #' Randomly samples (without replacement) a given number of pot pairs into
    #' a subset of the total data.
    #'
    #' Automatically numbers output files if generating multiple subsets.
    #' @param sets Number of subsets to form
    #' @param n Number of pairs to sample
    #' @param save_dir Location to save output files
    #' @param taxa_data_filename Name for taxa data file (don't include extension)
    #' @param character_data_filename Name for character data file (don't include extension)
    #' @param n_PCs Number of PCs to use in pot pairing algorithm
    #' @return
    #' pair_sample_character_data
    #'
    #' pair_sample_taxa_data
    #'
    #' Optional - saves outputs to a directory
    pot_pair_sample = function(sets = 1,
                               n=25,
                               save_dir = "",
                               taxa_data_filename="pair_sample_taxa_data",
                               character_data_filename="pair_sample_character_data",
                               n_PCs=15) {
      if (!"pair" %in% colnames(self$data)) {
        self$make_pot_pairs(n_PCs = n_PCs)
      }

      if (sets == 1) {
        pairs <- sample(unique(self$data$pair), n, replace=FALSE)
        pair_names <- unique(self$data$type[self$data$pair %in% pairs])

        self$pair_sample_character_data <- self$character_data[rownames(self$character_data) %in% pair_names,]
        self$pair_sample_taxa_data <- self$taxa_data[self$taxa_data$taxon %in% pair_names,]

        if (save_dir != "") {
          char_filepath <- paste0(save_dir, "/", character_data_filename, ".nex")
          taxa_filepath <- paste0(save_dir, "/", taxa_data_filename, ".tsv")

          write.nexus.data(self$pair_sample_character_data, file = char_filepath, format = "continuous")
          write_tsv(self$pair_sample_taxa_data, taxa_filepath)
        }
      }

      if (sets > 1) {
        for (i in 1:sets) {
          pairs <- sample(unique(self$data$pair), n, replace=FALSE)
          pair_names <- unique(self$data$type[self$data$pair %in% pairs])

          pair_sample_character_data <- self$character_data[rownames(self$character_data) %in% pair_names,]
          pair_sample_taxa_data <- self$taxa_data[self$taxa_data$taxon %in% pair_names,]

          if (save_dir != "") {
            char_filepath <- paste0(save_dir, "/", character_data_filename, "_", i, ".nex")
            taxa_filepath <- paste0(save_dir, "/", taxa_data_filename, "_", i, ".tsv")

            write.nexus.data(pair_sample_character_data, file = char_filepath, format = "continuous")
            write_tsv(pair_sample_taxa_data, taxa_filepath)
          }
        }
      }
    }
  ),
  private = list(
    #' description
    #' Converts dates from BCE/AD to years before present
    #' param date Date
    convertDate = function(date) {
      if (is.na(date))
      {
        return(NA)
      }
      date <- as.numeric(date)
      if (date <= 0) {
        newDate <- 1950+abs(date)
      }
      else {
        newDate <- 1950-date
      }
      return(newDate)
    },

    #' description
    #' Generates a single pot image.
    #'
    #' Designed to work with dplyr
    #' param .x Data
    #' param .y Grouping variables
    #' param save_dir Location to save images
    #' param palette Colour palette
    #' return
    #' Pot image
    pot_image = function(.x, .y, save_dir, palette) {
      .x <- .x %>% arrange(point_order)

      plot <- ggplot(.x, aes(x=x, y=y, fill=shape))+
        geom_polygon()+
        coord_fixed(1, xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
        scale_fill_manual(values = palette) +
        theme_void() +
        theme(
          panel.background  = element_rect(fill = "transparent", colour = NA),
          plot.background   = element_rect(fill = "transparent", colour = NA),
          legend.position = "none",
          plot.margin = margin(0,0,0,0),
          aspect.ratio=1)

        name <- .y$type
        file = paste0(save_dir, "/", name, ".png")
        ggsave(
          filename = file,
          plot     = plot,
          bg       = "transparent",
          width    = 5,
          height   = 5,
          dpi      = 300
        )
      },

    #' description
    #' Divides taxa into time bins
    #' param bin Size of time bin in years
    form_time_bins = function(bin = 100) {
      self$data <- self$data %>%
        mutate(time_bin = floor(date_earliest / bin) * bin)
    },

    #' description
    #' HCA pair clustering algorithm
    #' param .data Data
    #' param n_PCs Number of PCs to cluster by (default = 15)
    #' param self_ref Allows function to write to class
    #' return
    #' Pair identities
    cluster_pairs = function(.data, n_PCs = 15, self_ref) {

      all_types <- unique(.data$type)

      # If only one type exists, mark everything as unpaired and exit
      if (length(all_types) < 2) {
        self_ref$data$pair <- "unpaired"
        return()
      }

      pca_subset <- self_ref$PCA_data %>%
        filter(id %in% all_types) %>%
        select(id, 1:n_PCs)

      cluster_matrix <- as.matrix(pca_subset[, -1])
      rownames(cluster_matrix) <- pca_subset$id

      d  <- dist(cluster_matrix, method = "euclidean")
      hc <- hclust(d, method = "ward.D2")

      ordered_taxa <- pca_subset$id[hc$order]

      n_taxa   <- length(ordered_taxa)
      pair_idx <- 1L

      i <- 1L
      while (i <= n_taxa) {
        if (i + 1L <= n_taxa) {
          pair_taxa  <- ordered_taxa[i:(i + 1L)]
          pair_label <- as.character(pair_idx)
          pair_idx   <- pair_idx + 1L
          i          <- i + 2L
        } else {
          pair_taxa  <- ordered_taxa[i]
          pair_label <- "unpaired"
          i          <- i + 1L
        }

        self_ref$data$pair[self_ref$data$type %in% pair_taxa] <- pair_label
      }
    }
  )
)
