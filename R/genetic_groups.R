#' @title Identify Genetic Groups
#' @description Identifies genetic groups from a pairwise relatedness matrix using a graph-based approach.
#' @param relate A data.frame with columns: `ind1`, `ind2`, and the estimator column. Based on the output format of the `coancestry` function of the `related` package (Pew et al., 2015).
#' @param estimator The name of the column in `relate` containing relatedness estimates.
#' @param threshold The minimum relatedness value to consider for grouping (default: 0.4).
#' @param samples An `sf` object with spatial data for individuals.
#' @return An `sf` object with individuals and their assigned genetic groups.
#' @importFrom igraph graph_from_adjacency_matrix plot
#' @importFrom dplyr rows_patch filter mutate
#' @export
#' @examples
#' relate <- read.csv("path/to/relatedness_matrix.csv", sep = ";", dec = ",")
#' samples <- st_read("path/to/spatial_data.gpkg")
#' obs <- genetic_groups(relate, "wang", 0.4, samples)
genetic_groups <- function(relate, estimator, threshold = 0.4, samples) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required but not installed.")
  }

  # list of individuals
  indiv_list <- as.data.frame(unique(c(relate$ind1, relate$ind2)))

  # Filter and format the relatedness matrix
  relate <- relate[, c("ind1", "ind2", estimator)] %>%
    filter(.data[[estimator]] > threshold)

  # Add diagonals
  lvl1 <- levels(as.factor(relate$ind1))
  lvl2 <- levels(as.factor(relate$ind2))
  ind <- unique(c(lvl1, lvl2))
  tab <- data.frame(ind1 = ind, ind2 = ind, wang = NA)
  relate <- rbind(relate, tab)

  # Create the similarity matrix
  relate <- xtabs(relate[, 3] ~ relate[, 2] + relate[, 1],
                  addNA = TRUE, drop.unused.levels = FALSE, sparse = TRUE)

  # Create the graph and detect groups
  g <- graph_from_adjacency_matrix(relate, weighted = TRUE, diag = FALSE, mode = "lower")
  plot(g, main = "Graph of Genetic Relationships")
  dec <- as.data.frame(components(g)$membership)
  colnames(dec) <- c("group")
  dec$Individual <- rownames(dec)

  # Handle individuals without a group
  colnames(indiv_list) <- c("Individual")
  indiv_without_group <- subset(indiv_list, !(indiv_list$Individual %in% dec$Individual))
  if (nrow(indiv_without_group) > 0) {
    indiv_without_group <- indiv_without_group %>%
      mutate(group = max(dec$group) + row_number())
    groups <- merge(x = indiv_list, y = dec, by = "Individual", all.x = TRUE)
    groups <- rows_patch(groups, indiv_without_group, by = "Individual")
  } else {
    groups <- merge(x = indiv_list, y = dec, by = "Individual", all.x = TRUE)
  }

  # Merge with spatial data
  obs <- merge(x = samples, y = groups, by = "Individual")
  return(obs)
}
