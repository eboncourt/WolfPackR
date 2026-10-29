#' @title Identify Genetic Groups
#' @description Identifies genetic groups from pairwise relatedness: pairs above `threshold` become the edges
#'   of a graph (weighted by relatedness), then groups are found either as connected components or with a
#'   community detection algorithm.
#' @param relate A data.frame with columns `ind1`, `ind2` and the estimator column, e.g. the output of
#'   [relatedness()] or of the `coancestry` function of the `related` package (Pew et al., 2015).
#'   One row per pair is enough, in any order.
#' @param estimator The name of the column in `relate` containing relatedness estimates.
#' @param threshold The minimum relatedness value to consider for grouping (default: 0.4).
#'   Pick it according to the estimator: for QG or LR, 0.25 is the value expected for second-degree relatives
#'   (0.5 for first-degree ones), so pairs above it are likely close kin; for a LOD, 3 is a classic choice.
#'   Must be >= 0 when the estimator can be negative (QG, LR): a relatedness <= 0 can't be a kin link.
#' @param method How groups are built from the graph:
#'   * `"components"` (default, behaviour of earlier versions): connected components, i.e. two individuals
#'     end up together as soon as a chain of links joins them. Simple, but a single spurious link merges two families;
#'   * `"leiden"`, `"louvain"`, `"walktrap"`: community detection (igraph), which looks for densely linked
#'     groups and is much less sensitive to isolated links. Leiden and Louvain maximise modularity.
#' @param seed Random seed for the stochastic methods (Leiden, Louvain), so results are reproducible.
#'   Your own random number stream is left untouched. `NULL` to skip.
#' @param plot Draw the graph of genetic relationships (default: FALSE).
#' @return A data.frame with one row per individual: `Individual` and its `group` (integer).
#'   Individuals without any link above the threshold get a group of their own.
#' @importFrom igraph graph_from_adjacency_matrix plot.igraph components cluster_leiden cluster_louvain cluster_walktrap membership V
#' @importFrom dplyr rows_patch filter mutate row_number %>% select across .data
#' @importFrom stats xtabs
#' @export
#' @examples
#' relate <- data.frame(ind1 = c("W1", "W1", "W1", "W2", "W2", "W3"), ind2 = c("W1", "W2", "W3", "W2", "W3", "W3"), indicator = c(1.0, 0.5, 0.2, 1.0, 0.3, 1.0))
#' result <- genetic_groups(relate, estimator = "indicator", threshold = 0.4)
#'
#' # With relatedness computed from genotypes and community detection
#' data("genotypes_Carpathians")
#' rel <- relatedness(genotypes_Carpathians)
#' groups <- genetic_groups(rel, estimator = "QG", threshold = 0.25, method = "leiden")
genetic_groups <- function(relate, estimator, threshold = 0.4,
                           method = c("components", "leiden", "louvain", "walktrap"), seed = 1, plot = FALSE) {
  method <- match.arg(method)

  # list of individuals
  indiv_list <- data.frame(Individual = unique(c(as.character(relate$ind1), as.character(relate$ind2))))

  # Filter and format the relatedness matrix
  relate <- relate %>%
    mutate(across(c(ind1, ind2), as.character)) %>%  # Convert to character to avoid problems with factors in dplyr
    select(ind1, ind2, {{estimator}}) %>%
    filter(.data[[estimator]] > threshold)

  # Links need a positive weight: a relatedness <= 0 can't mean kinship, and igraph versions don't agree
  # on what to do with such edges (igraph 2 silently drops them when symmetrising, igraph 1.x keeps them)
  if (any(relate[[estimator]] <= 0)) stop("Kin links need positive weights: use a threshold >= 0.")

  # Nobody above the threshold: everyone is on their own
  if (nrow(relate) == 0) {
    indiv_list$group <- seq_len(nrow(indiv_list))
    return(indiv_list[order(indiv_list$Individual), ])
  }

  # Add diagonals so that every retained individual shows up as a node
  ind <- sort(unique(c(relate$ind1, relate$ind2)))
  tab <- data.frame(ind1 = ind, ind2 = ind, wang = 1)
  names(tab)[3] <- estimator
  relate <- rbind(relate, tab)

  # Build the square similarity matrix. Rows and columns share the same levels,
  # so a pair ends up in the same spot whatever the order of ind1 / ind2
  relate <- xtabs(relate[, 3] ~ factor(relate[, 2], levels = ind) + factor(relate[, 1], levels = ind),
                  addNA = TRUE, drop.unused.levels = FALSE, sparse = TRUE)

  # Create the graph
  # mode = "max": a pair counts as linked whichever triangle it sits in (coancestry gives one row per pair,
  # in no particular order, so "lower" silently dropped about half of them)
  g <- graph_from_adjacency_matrix(relate, weighted = TRUE, diag = FALSE, mode = "max")
  if (plot) plot.igraph(g, main = "Graph of Genetic Relationships")

  # Detect groups: connected components, or communities (weighted by relatedness)
  if (!is.null(seed) && method %in% c("leiden", "louvain")) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (is.null(old_seed)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old_seed, envir = globalenv()))
    set.seed(seed)
  }
  memb <- switch(method,
    components = components(g)$membership,
    leiden     = membership(cluster_leiden(g, objective_function = "modularity", n_iterations = -1)),
    louvain    = membership(cluster_louvain(g)),
    walktrap   = membership(cluster_walktrap(g))
  )
  dec <- data.frame(group = as.integer(memb), Individual = V(g)$name)

  # Handle individuals without a group
  indiv_without_group <- subset(indiv_list, !(indiv_list$Individual %in% dec$Individual))
  if (nrow(indiv_without_group) > 0) {
    indiv_without_group <- indiv_without_group %>%
      mutate(group = max(dec$group) + dplyr::row_number())
    groups <- merge(x = indiv_list, y = dec, by = "Individual", all.x = TRUE)
    groups <- rows_patch(groups, indiv_without_group, by = "Individual")
  } else {
    groups <- merge(x = indiv_list, y = dec, by = "Individual", all.x = TRUE)
  }

  return(groups)
}
