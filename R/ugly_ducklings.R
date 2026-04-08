#' @title Identify and Assign Ugly Ducklings and Lone Individuals
#' @description Finds individuals that are genetically linked but spatially isolated, and assigns them to the most suitable pack.
#' @param obs An `sf` object with individuals and their groups (genetic or spatial groups). The first column must contain the individuals' ID.
#' @param group Name of the column containing the identifier of groups of individuals, which may be genetic or spatial groups.
#' @param min_overlap The minimum overlap percentage to consider an individual as integrated (default: 0.5).
#' @param buffer The buffer distance around the MCP of the group to consider for inclusion (default: 0).
#' @param mcp.percent The percentile of points (distance from centroid) excluded before calculating the MCP (default: 100).
#' @return A data.frame containing all individuals, their genetic group, and their assigned pack.
#' @importFrom sf st_buffer st_within st_geometry st_cast st_combine st_intersection
#' @importFrom dplyr filter
#' @export
#' @examples
#' data("relate")
#' data("samples")
#' gengroups <- genetic_groups(relate, threshold = 0.4, estimator = "wang")
#' obs <- merge(samples, gengroups, by = "Individual")
#' ud_df <- ugly_ducklings(obs, group, min_overlap = 0.5, buffer = 1000, mcp.percent = 95)

ugly_ducklings <- function(obs, group, min_overlap = 0.5, buffer = 0, mcp.percent = 100) {
  group_name <- deparse(substitute(group))

  # Initialize results data.frame
  results <- data.frame(
    Individual = obs[[1]],
    Genetic_Group = obs[[group_name]],
    Pack = as.character(obs[[group_name]]),
    stringsAsFactors = FALSE
  )

  # Step 1: Assign status to each individual (In Group or Lone Individual)
  obs$status <- "In Group"
  groups <- unique(obs[[group_name]])
  colnames(obs)[1] <- "Individual"

  for (gengroup in groups) {
    group_obs <- obs %>% filter(.data[[group_name]] == gengroup)
    group_individuals <- unique(group_obs$Individual)

    for (ind in group_individuals) {
      ind_obs <- group_obs %>% filter(Individual == ind)
      rest_obs <- group_obs %>% filter(Individual != ind)

      if (nrow(ind_obs) == 1) {
        if (nrow(rest_obs) >= 3) {
          mcp_rest <- mcp_sf(rest_obs, percentile = mcp.percent)
          if (!is.null(mcp_rest)) {
            buffered_mcp <- st_buffer(mcp_rest, dist = buffer)
            point_within_mcp <- st_within(st_geometry(ind_obs), buffered_mcp, sparse = FALSE)
            if (!all(point_within_mcp[, 1])) {
              obs$status[obs$Individual == ind] <- "Lone Individual"
            }
          } else {
            obs$status[obs$Individual == ind] <- "Lone Individual"
          }
        } else {
          obs$status[obs$Individual == ind] <- "Lone Individual"
        }
      } else if (nrow(ind_obs) == 2) {
        if (nrow(rest_obs) >= 3) {
          mcp_rest <- mcp_sf(rest_obs, percentile = mcp.percent)
          if (!is.null(mcp_rest)) {
            buffered_mcp <- st_buffer(mcp_rest, dist = buffer)
            line_segment <- st_cast(st_combine(st_geometry(ind_obs)), "LINESTRING")
            intersection <- tryCatch(
              { st_intersection(line_segment, buffered_mcp) },
              error = function(e) NULL
            )
            if (is.null(intersection) || length(intersection) == 0) {
              overlap_ratio <- 0
            } else {
              intersection_length <- sum(as.numeric(st_length(intersection)))
              total_length <- as.numeric(st_length(line_segment))
              overlap_ratio <- intersection_length / total_length
            }
            if (overlap_ratio <= min_overlap) {
              obs$status[obs$Individual == ind] <- "Lone Individual"
            }
          } else {
            obs$status[obs$Individual == ind] <- "Lone Individual"
          }
        } else {
          obs$status[obs$Individual == ind] <- "Lone Individual"
        }
      }
    }
  }

  # Step 2: Identify Ugly Ducklings among individuals with >= 3 points
  for (gengroup in groups) {
    group_obs <- obs %>% filter(.data[[group_name]] == gengroup & status != "Lone Individual")

    # Check that group_obs$Individual is not empty
    if (nrow(group_obs) > 0) {
      effective_ind <- as.data.frame(table(group_obs$Individual))

      # Check that effective_ind is not empty
      if (nrow(effective_ind) > 0) {
        effective_ind <- effective_ind[effective_ind$Freq >= 3, ]
        n_effective_ind <- nrow(effective_ind)

        if (n_effective_ind > 0) {
          for (k in 1:n_effective_ind) {
            ind_id <- as.character(effective_ind$Var1[k])
            pack_obs <- obs %>% filter(.data[[group_name]] == gengroup & Individual != ind_id & status != "Lone Individual")
            if (nrow(pack_obs) >= 3) {
              mcp_pack <- mcp_sf(pack_obs, percentile = mcp.percent)
              if (!is.null(mcp_pack)) {
                ind_obs <- obs %>% filter(Individual == ind_id)
                mcp_ind <- mcp_sf(ind_obs, percentile = mcp.percent)
                if (!is.null(mcp_ind)) {
                  a <- as.numeric(st_area(st_intersection(mcp_pack, mcp_ind))) / as.numeric(st_area(mcp_ind))
                  if (length(a) == 0 || is.nan(as.numeric(a)) || as.numeric(a) <= min_overlap) {
                    results$Pack[results$Individual == ind_id] <- "Ugly Duckling"
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  # Step 3: Assign Ugly Ducklings and Lone Individuals to the most suitable pack
  for (i in 1:nrow(results)) {
    ind_id <- results$Individual[i]
    if (results$Pack[i] %in% c("Ugly Duckling", "Lone Individual")) {
      ind_obs <- obs %>% filter(Individual == ind_id)
      mcp_ind <- mcp_sf(ind_obs, percentile = mcp.percent)
      if (!is.null(mcp_ind)) {
        max_overlap <- 0
        for (gengroup in groups) {
          pack_obs <- obs %>% filter(.data[[group_name]] == gengroup & status != "Lone Individual" & Individual != ind_id)
          if (nrow(pack_obs) >= 3) {
            mcp_pack <- mcp_sf(pack_obs, percentile = mcp.percent)
            if (!is.null(mcp_pack)) {
              buffered_mcp_pack <- st_buffer(mcp_pack, dist = buffer)
              a <- as.numeric(st_area(st_intersection(buffered_mcp_pack, mcp_ind))) / as.numeric(st_area(mcp_ind))
              if (length(a) > 0 && !is.nan(as.numeric(a)) && as.numeric(a) > max_overlap) {
                max_overlap <- a
                best_pack <- as.character(gengroup)
              }
            }
          }
        }
        if (exists("best_pack")) {
          results$Pack[results$Individual == ind_id] <- best_pack
        } else {
          results$Pack[results$Individual == ind_id] <- as.character(results$Genetic_Group[i])
        }
      }
    }
  }

  # Mark Lone Individuals in results
  results$Pack[results$Individual %in% obs$Individual[obs$status == "Lone Individual"]] <- "Lone Individual"

  # Adding the new_group column
  results$new_group <- results$Pack

  # Assignment of a unique identifier for each ‘Lone Individual’
  lone_individuals <- which(results$new_group == "Lone Individual")
  if (length(lone_individuals) > 0) {
    unique_ids <- paste0("Lone_", seq_len(length(lone_individuals)))
    results$new_group[lone_individuals] <- unique_ids
  }

  return(unique(results))
}
