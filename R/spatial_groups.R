#' @title Identify spatial groups or subgroups
#'
#' This function identifies spatial groups or subgroups within groups using Minimum Convex Polygons (MCP)
#' and clustering based on spatial overlap.
#'
#' @param obs An sf object containing the observations, with columns for individual IDs and coordinates.
#' @param group The name of the column in `obs` that contains the group information. If not provided, all individuals are treated as a single group.
#' @param percentile The percentile for the MCP calculation (e.g., 95 for 95% MCP).
#' @param buffer_radius The buffer radius to use if the MCP cannot be calculated.
#' @param max_iterations The maximum number of iterations for subgroup identification.
#'
#' @return A data.frame with updated subgroup assignments.
#'
#' @examples
#' # Example usage:
#' obs_sf <- st_as_sf(obs, coords = c("Longitude", "Latitude"), crs = 4326)
#' subgroups_result <- spatial_groups(obs = obs_sf, group = group, percentile = 95, buffer_radius = 1000)
#' # Without group
#' subgroups_result <- spatial_groups(obs = obs_sf, percentile = 95, buffer_radius = 1000)
#'
#' @export
spatial_groups <- function(obs, group = NULL, percentile = 100, buffer_radius = 1, max_iterations = 20) {
  # Initialize results
  results <- data.frame(
    Individual = obs$Individual,
    Group = NA, 
    Subgroup = "All_Individuals",
    stringsAsFactors = FALSE
  )

  # If group is provided, use it; otherwise, treat all individuals as a single group
  if (!is.null(group)) {
    group_name <- deparse(substitute(group))
    results$Group <- obs[[group_name]]
    groups <- unique(obs[[group_name]])
  } else {
    groups <- "All_Individuals"
    results$Group <- "All_Individuals"
  }

  for (gengroup in groups) {
    if (!is.null(group)) {
      group_obs <- obs %>% dplyr::filter(.data[[group_name]] == gengroup)
    } else {
      group_obs <- obs
    }

    group_individuals <- unique(group_obs$Individual)

    # Calculate initial MCP for the group
    mcp_group <- mcp_sf(group_obs, percentile, buffer_radius)

    # If the MCP cannot be calculated, all individuals are marked as "Lone Individual"
    if (is.null(mcp_group)) {
      results$Subgroup[results$Individual %in% group_individuals] <- "Lone Individual"
      next
    }

    # Identify spatial subgroups
    subgroups_list <- list(group_individuals)
    subgroup_mcps <- list(mcp_group)

    stable <- FALSE
    iteration <- 0

    while (!stable && iteration < max_iterations) {
      stable <- TRUE
      new_subgroups <- list()
      new_subgroup_mcps <- list()

      # For each existing subgroup
      for (i in seq_along(subgroups_list)) {
        current_subgroup <- subgroups_list[[i]]
        current_mcp <- if (i <= length(subgroup_mcps)) subgroup_mcps[[i]] else NULL

        # If the subgroup has fewer than 3 individuals, keep it as is
        if (length(current_subgroup) < 3) {
          new_subgroups <- c(new_subgroups, list(current_subgroup))
          new_subgroup_mcps <- c(new_subgroup_mcps, list(current_mcp))
          next
        }

        # Calculate individual MCPs
        individual_mcps <- lapply(current_subgroup, function(ind) {
          ind_obs <- obs %>% dplyr::filter(Individual == ind)
          mcp_sf(ind_obs, percentile, buffer_radius)
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
        graph <- igraph::graph_from_adjacency_matrix(overlap_matrix, weighted = TRUE, mode = "undirected", diag = FALSE)
        clusters <- igraph::components(graph)

        # Create new subgroups
        for (cluster_id in unique(clusters$membership)) {
          cluster_members <- current_subgroup[clusters$membership == cluster_id]
          new_subgroups <- c(new_subgroups, list(cluster_members))
          new_mcp <- mcp_sf(obs %>% dplyr::filter(Individual %in% cluster_members), percentile, buffer_radius)
          new_subgroup_mcps <- c(new_subgroup_mcps, list(new_mcp))
        }
      }

      # Check stability
      if (length(new_subgroups) != length(subgroups_list)) {
        stable <- FALSE
      } else {
        same_groups <- TRUE
        for (i in seq_along(new_subgroups)) {
          if (i > length(subgroups_list) || !setequal(new_subgroups[[i]], subgroups_list[[i]])) {
            same_groups <- FALSE
            break
          }
        }
        if (!same_groups) {
          stable <- FALSE
        }
      }

      subgroups_list <- new_subgroups
      subgroup_mcps <- new_subgroup_mcps
      iteration <- iteration + 1
    }

    # Assign spatial subgroups to results
    for (i in seq_along(subgroups_list)) {
      if (length(subgroups_list[[i]]) >= 2) {
        if (!is.null(group)) {
          results$Subgroup[results$Individual %in% subgroups_list[[i]]] <- paste0(gengroup, "_Subgroup_", i)
        } else {
          results$Subgroup[results$Individual %in% subgroups_list[[i]]] <- paste0("Subgroup_", i)
        }
      } else {
        results$Subgroup[results$Individual %in% subgroups_list[[i]]] <- "Lone Individual"
      }
    }
  }

  return(results)
}
