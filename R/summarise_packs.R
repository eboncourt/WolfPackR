#' @title Summarise Packs and Individuals
#' @description Summarises individuals and their packs, including putative dominant individuals.
#' @param obs An `sf` object with individuals and their packs.
#' @param pack The name of the column in `obs` that contains the pack information.
#' @param sex_column The name of the column in obs that contains sex information (default: "Sex").
#' @param male_pattern A regex pattern to identify males in Individual IDs (default: "M").
#' @param female_pattern A regex pattern to identify females in Individual IDs (default: "F").
#' @return A list containing detailed summaries of packs and individuals.
#' @importFrom dplyr filter
#' @importFrom sf st_convex_hull st_union st_area
#' @export
#' @examples
#' summary <- summarise_packs(obs, pack = Pack)
summarise_packs <- function(obs, pack, sex_column = "Sex", male_pattern = "M", female_pattern = "F") {
  pack_name <- deparse(substitute(pack))

  # Vérifier que la colonne pack existe
  if (!pack_name %in% names(obs)) {
    stop(paste("La colonne", pack_name, "n'existe pas dans les données."))
  }

  # Count the number of points per individual directly from obs
  point_counts <- as.data.frame(table(obs$Individual))
  colnames(point_counts) <- c("Individual", "PointCount")

  # Merge point_counts with obs to get the number of points per individual
  obs_with_points <- merge(obs, point_counts, by = "Individual", all.x = TRUE)

  # Create a summary list for each pack
  packs <- unique(obs_with_points[[pack_name]])
  pack_summaries <- list()

  for (pack in packs) {
    pack_members <- obs_with_points[obs_with_points[[pack_name]] == pack, ]
    total_individuals <- length(unique(pack_members$Individual))
    individual_names <- unique(pack_members$Individual)
    point_counts_list <- point_counts[point_counts$Individual %in% individual_names, ]

    # Total points in the pack
    total_points <- sum(point_counts_list$PointCount, na.rm = TRUE)

    # Determine putative dominant individuals (male and female with the most points)
    putative_dominant_male <- NA
    putative_dominant_female <- NA

    # Calculate territory area based on convex hull
    territory_area <- NA

    if (pack != "Lone Individual") {
      # Calculate territory area based on convex hull
      pack_obs <- obs_with_points %>% dplyr::filter(.data[[pack_name]] == pack)
      mcp_pack <- st_convex_hull(st_union(pack_obs))
      territory_area <- st_area(mcp_pack)

      # Check if sex information is available in obs
      if (sex_column %in% colnames(obs)) {
        pack_obs <- obs[obs$Individual %in% individual_names, ]
        males <- pack_obs[pack_obs[[sex_column]] == "M", ]
        females <- pack_obs[pack_obs[[sex_column]] == "F", ]

        if (nrow(males) > 0) {
          male_point_counts <- point_counts[point_counts$Individual %in% males$Individual, ]
          max_points_male <- max(male_point_counts$PointCount)
          putative_dominant_male <- as.character(male_point_counts$Individual[male_point_counts$PointCount == max_points_male])
          if (length(putative_dominant_male) > 1) putative_dominant_male <- paste(putative_dominant_male, collapse = ", ")
        }

        if (nrow(females) > 0) {
          female_point_counts <- point_counts[point_counts$Individual %in% females$Individual, ]
          max_points_female <- max(female_point_counts$PointCount)
          putative_dominant_female <- as.character(female_point_counts$Individual[female_point_counts$PointCount == max_points_female])
          if (length(putative_dominant_female) > 1) putative_dominant_female <- paste(putative_dominant_female, collapse = ", ")
        }
      } else {
        # If no sex column, try to infer from Individual ID using grepl
        males <- individual_names[grepl(male_pattern, individual_names, ignore.case = FALSE)]
        females <- individual_names[grepl(female_pattern, individual_names, ignore.case = FALSE)]

        if (length(males) > 0) {
          male_point_counts <- point_counts[point_counts$Individual %in% males, ]
          max_points_male <- max(male_point_counts$PointCount)
          putative_dominant_male <- as.character(male_point_counts$Individual[male_point_counts$PointCount == max_points_male])
          if (length(putative_dominant_male) > 1) putative_dominant_male <- paste(putative_dominant_male, collapse = ", ")
        }

        if (length(females) > 0) {
          female_point_counts <- point_counts[point_counts$Individual %in% females, ]
          max_points_female <- max(female_point_counts$PointCount)
          putative_dominant_female <- as.character(female_point_counts$Individual[female_point_counts$PointCount == max_points_female])
          if (length(putative_dominant_female) > 1) putative_dominant_female <- paste(putative_dominant_female, collapse = ", ")
        }
      }
    }

    # Create a summary for the pack
    pack_summary <- list(
      Pack = pack,
      TotalIndividuals = total_individuals,
      TotalPoints = total_points,
      Individuals = individual_names,
      PointCounts = point_counts_list,
      PutativeDominantMale = putative_dominant_male,
      PutativeDominantFemale = putative_dominant_female,
      TerritoryArea = territory_area
    )

    pack_summaries[[as.character(pack)]] <- pack_summary
  }

  return(pack_summaries)
}
