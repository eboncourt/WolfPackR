#' Identify Spatial Subgroups
#'
#' This function identifies spatial subgroups within genetic groups using Minimum Convex Polygons (MCP)
#' and clustering based on spatial overlap.
#'
#' @param obs An sf object containing the observations, with columns for individual IDs and coordinates.
#' @param genetic_group The name of the column in `obs` that contains the genetic group information.
#' @param mcp.percent The percentile for the MCP calculation (e.g., 95 for 95% MCP).
#' @param buffer_radius The buffer radius to use if the MCP cannot be calculated.
#' @param max_iterations The maximum number of iterations for subgroup identification.
#'
#' @return A data.frame with updated subgroup assignments.
#'
#' @examples
#' # Example usage:
#' obs_sf <- st_as_sf(obs, coords = c("Longitude", "Latitude"), crs = 4326)
#' subgroups_result <- subgroups(obs = obs_sf, genetic_group = "group", mcp.percent = 95, buffer_radius = 1000)
#'
#' @export
subgroups <- function(obs, genetic_group, mcp.percent = 100, buffer_radius = 1, max_iterations = 20) {
  
  # Initialize results
  results <- data.frame(Individual = obs$Individual, Genetic_Group = obs[[genetic_group]], Subgroup = obs[[genetic_group]], stringsAsFactors = FALSE)

  # Iterate for each genetic group
  groups <- unique(obs[[genetic_group]])

  for (gengroup in groups) {
    group_obs <- obs[obs[[genetic_group]] == gengroup, ]
    group_individuals <- unique(group_obs$Individual)

    # Calculate initial MCP for the genetic group
    mcp_group <- calculate_mcp(group_individuals, obs, mcp.percent, buffer_radius)

    # If the MCP cannot be calculated, all individuals are marked as "Lone Individual"
    if (is.null(mcp_group)) {
      results$Subgroup[results$Individual %in% group_individuals] <- "Lone Individual"
      next
    }

    # Identify spatial subgroups
    subgroups <- list(group_individuals)
    subgroup_mcps <- list(mcp_group)

    stable <- FALSE
    iteration <- 0

    while (!stable && iteration < max_iterations) {
      stable <- TRUE
      new_subgroups <- list()
      new_subgroup_mcps <- list()

      # For each existing subgroup
      for (i in seq_along(subgroups)) {
        current_subgroup <- subgroups[[i]]
        current_mcp <- if (i <= length(subgroup_mcps)) subgroup_mcps[[i]] else NULL

        # If the subgroup has fewer than 3 individuals, keep it as is
        if (length(current_subgroup) < 3) {
          new_subgroups <- c(new_subgroups, list(current_subgroup))
          new_subgroup_mcps <- c(new_subgroup_mcps, list(current_mcp))
          next
        }

        # Calculate individual MCPs
        individual_mcps <- lapply(current_subgroup, function(ind) {
          ind_obs <- obs[obs$Individual == ind, ]
          calculate_mcp(c(ind), obs, mcp.percent, buffer_radius)
        })
        names(individual_mcps) <- current_subgroup

        # Identify consistent subgroups
        overlap_matrix <- matrix(0, nrow = length(current_subgroup), ncol = length(current_subgroup))
        for (j in seq_along(current_subgroup)) {
          for (k in seq_along(current_subgroup)) {
            if (j != k) {
              overlap_value <- calculate_overlap(individual_mcps[[j]], individual_mcps[[k]])
              if (is.numeric(overlap_value) && length(overlap_value) == 1) {
                overlap_matrix[j, k] <- overlap_value
              } else {
                overlap_matrix[j, k] <- 0
              }
            }
          }
        }

        # Cluster individuals into spatial subgroups
        diag(overlap_matrix) <- 1
        graph <- graph_from_adjacency_matrix(overlap_matrix, weighted = TRUE, mode = "undirected", diag = FALSE)
        clusters <- clusters(graph)

        # Create new subgroups
        for (cluster_id in unique(clusters$membership)) {
          cluster_members <- current_subgroup[clusters$membership == cluster_id]
          new_subgroups <- c(new_subgroups, list(cluster_members))
          new_mcp <- calculate_mcp(cluster_members, obs, mcp.percent, buffer_radius)
          new_subgroup_mcps <- c(new_subgroup_mcps, list(new_mcp))
        }
      }

      # Check stability
      if (length(new_subgroups) != length(subgroups)) {
        stable <- FALSE
      } else {
        same_groups <- TRUE
        for (i in seq_along(new_subgroups)) {
          if (i > length(subgroups) || !setequal(new_subgroups[[i]], subgroups[[i]])) {
            same_groups <- FALSE
            break
          }
        }
        if (!same_groups) {
          stable <- FALSE
        }
      }

      subgroups <- new_subgroups
      subgroup_mcps <- new_subgroup_mcps
      iteration <- iteration + 1
    }

    # Assign spatial subgroups to results
    for (i in seq_along(subgroups)) {
      if (length(subgroups[[i]]) >= 3) {
        results$Subgroup[results$Individual %in% subgroups[[i]]] <- paste0(gengroup, "_Subgroup_", i)
      } else {
        results$Subgroup[results$Individual %in% subgroups[[i]]] <- "Lone Individual"
      }
    }
  }

  return(results)
}
