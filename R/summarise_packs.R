#' @title Summarise Packs and Individuals
#' @description Summarises individuals and their packs, including putative dominant individuals.
#' @param obs An `sf` object with individuals and their genetic groups.
#' @param ud_df A data.frame containing individuals, their genetic group, and their assigned pack.
#' @param sex_column The name of the column in obs that contains sex information (default: "Sex").
#' @param male_pattern A regex pattern to identify males in Individual IDs (default: "M").
#' @param female_pattern A regex pattern to identify females in Individual IDs (default: "F").
#' @return A list containing detailed summaries of packs and individuals.
#' @examples
#' summary <- summarise_packs(obs, ud_df, sex_column = "Sex")
summarise_packs <- function(obs, ud_df, sex_column = "Sex", male_pattern = "M", female_pattern = "F") {
  # Count the number of points per individual directly from obs
  point_counts <- as.data.frame(table(obs$Individual))
  colnames(point_counts) <- c("Individual", "PointCount")
  
  # Merge ud_df with point_counts to get the number of points per individual
  ud_df_with_points <- merge(ud_df, point_counts, by = "Individual", all.x = TRUE)
  
  # Create a summary list for each pack
  packs <- unique(ud_df$Pack)
  pack_summaries <- list()
  
  for (pack in packs) {
    pack_members <- ud_df_with_points[ud_df_with_points$Pack == pack, ]
    total_individuals <- length(unique(pack_members$Individual))
    individual_names <- unique(pack_members$Individual)
    point_counts_list <- point_counts[point_counts$Individual %in% individual_names, ]
    genetic_groups <- unique(pack_members$Genetic_Group)
    
    # Total points in the pack
    total_points <- sum(point_counts_list$PointCount, na.rm = TRUE)
    
    # Determine putative dominant individuals (male and female with the most points)
    putative_dominant_male <- NA
    putative_dominant_female <- NA
    
    # Calculate territory area based on convex hull
    territory_area <- NA
    
    if (pack != "Lone Individual") {
      
      # Calculate territory area based on convex hull
      obs_with_ud <- merge(obs, ud_df, by = "Individual", all.x = TRUE)
      pack_obs <- obs_with_ud %>% filter(Pack == pack)
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
        
        if (nrow(femelles) > 0) {
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
      GeneticGroups = genetic_groups,
      PutativeDominantMale = putative_dominant_male,
      PutativeDominantFemale = putative_dominant_female,
      TerritoryArea = territory_area
    )
    
    pack_summaries[[pack]] <- pack_summary
  }
  
  return(pack_summaries)
}