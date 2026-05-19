#' Pot Data Processor Class
#'
#' @description
#' Process pot data and generate basic plots
#' @importFrom R6 R6Class
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom scales hue_pal
#' @importFrom readr read_tsv
#' @importFrom ape as.phylo
#' @importFrom ape write.nexus.data
#'
#' @export
potDataProcessor <- R6Class("potDataProcessor",
  public = list(
    data = NULL,
    character_data = NULL,
    taxa_data = NULL,
    PCA_data = NULL,
    scaled_taxa_data = NULL,
    raw_pca_output = NULL,
    initialize = function(data) {
      self$data <- data
    },
    generate_revbayes_data = function(save_dir="", character_data_filename="character_data", taxa_data_filename="taxa_data") {
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

      no_age_range <- taxa_data$taxon[(taxa_data$max_age-taxa_data$min_age) == 0]
      if (length(no_age_range) > 0) {
        warning("The following taxa have been excluded due to having an age range of 0:\n",
                paste(no_age_range, collapse = "\n"),
                call. = FALSE)
        taxa_data <- taxa_data[!(taxa_data$taxon %in% no_age_range),]
      }

      # Only save characters which explain a proportion of variance > 0
      prop_var <- summary(pca)$importance[2,]
      useful_characters <- length(prop_var[prop_var > 0])

      character_data <- scores[,1:useful_characters]
      rownames(character_data) <- as.character(scores$id)
      # filter for only shared taxa between datasets
      shared_taxa <- intersect(rownames(character_data), taxa_data$taxon)
      character_data <- as.matrix(character_data)
      character_data    <- character_data[rownames(character_data) %in% shared_taxa, ]
      self$taxa_data    <- taxa_data[taxa_data$taxon %in% shared_taxa, ]
      self$character_data <- character_data

      if (save_dir != "") {
        character_data_path <- paste0(save_dir, "/", character_data_filename, ".nex")
        taxa_data_path <- paste0(save_dir, "/", taxa_data_filename, ".tsv")

        write.nexus.data(self$character_data, file = character_data_path, format = "continuous")
        write_tsv(self$taxa_data, taxa_data_path)
      }

    },
    PC_plot = function(cols=c("PC1", "PC2"), add_elipses=FALSE){
      if (is.null(self$PCA_data)) {
        self$generate_revbayes_data()
        warning("Input data not processed. Running self$generate_revbayes_data()")
      }

      plot <- ggplot(self$PCA_data, aes(x=.data[[cols[1]]], y=.data[[cols[2]]], colour = shape)) +
        geom_point()  +
        theme_bw() +
        coord_fixed(ratio = 1) +
        theme(
          legend.position = "bottom"
        )
      if (add_elipses) {
        plot <- plot + stat_ellipse()
      }

      print(plot)
    },
    dendrogram = function(n=10){
      if (is.null(self$PCA_data)) {
        self$generate_revbayes_data()
        warning("Input data not processed. Running self$generate_revbayes_data()")
      }

      if (n > ncol(self$PCA_data)-4) {
        stop("Supplied number of PCs greater than number of PCs in data")
      }

      cluster <- as.matrix(self$PCA_data[, 1:n])
      rownames(cluster) <- self$PCA_data$id

      cluster <- scale(cluster)
      d <- dist(cluster, method = "euclidean")
      hc <- hclust(d, method = "ward.D2")
      tree <- as.phylo(hc)

      metadata <- data.frame(
        label = tree$tip.label,
        type = self$PCA_data$shape)

      p <- ggtree(tree) %<+% metadata +
        geom_tiplab(aes(colour = type), size = 2)

      print(p)

    },
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
        theme_void() +                              # removes axes, gridlines, etc.
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
    generate_pot_images = function(save_dir) {
      groups      <- sort(unique(self$data$shape))
      n           <- length(groups)
      pot_palette <- setNames(hue_pal()(n), groups)
      #TODO adjust code

      if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

      self$data %>%
        group_by(type) %>%
        arrange(desc(point_order), .by_group = TRUE) %>%
        group_walk(private$pot_image, save_dir = save_dir, palette = pot_palette)

    },
    summarise_pot_profiles = function() {
      #TODO code it you dummy
    },
    scale_taxa_data = function(save_dir="", filename="scaled_taxa_data") {
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
    }
  ),
  private = list(
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
      }
  )
)
