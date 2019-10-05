# Simp function extract "effect" and "size columns from data frame

simp <- function(x) {
  return(x[,which(colnames(x) == "effect"):which(colnames(x) == "size")])}

#' Compute Average Treatment Effect
#'
#' \code{ATE} computes average treatment effect of the matched subsamples by
#' a weighted average of the estimated treatment effects in each matched group.
#' The weight is the number of matched units.
#'
#' @param FLAME_object object returned by applying the FLAME algorithm
#'   (\code{\link{FLAME_bit}}, \code{\link{FLAME_PostgreSQL}}, or
#'   \code{\link{FLAME_SQLite}})
#' @examples
#' data(toy_data)
#' result <- FLAME::FLAME_bit(data = toy_data, holdout = toy_data)
#' FLAME::ATE(result)
#' @return average treatment effect (ATE) of the matched subsamples
#' @export

rep_each_by <- function(x, n_reps) {
  out <- NULL
  for (i in seq_along(x)) {
    out <- c(out, rep(x[i], n_reps[i]))
  }
  return(out)
}

ATE <- function(FLAME_object) {
#  mg_sizes <- # Get sizes of each matched group
#    lapply(FLAME_object[['matched_group']], function(x) x$size) %>%
#    unlist() 
  
#  ATE <-
#    FLAME_object[['matched_data']] %>% 
#    filter(matched != 0) %>% # Exclude unmatched units
#    mutate(mg_id = unlist(mapply(rep, 1:length(mg_sizes), mg_sizes))) %>% # Give units IDs corresponding to their MG
#    group_by(mg_id) %>% # For each of these matched groups
#    summarize(CATE = mean(outcome[treated == 1]) - mean(outcome[treated == 0]),
#              mg_size = n()) %>% # Get the CATE and MG size
#    summarise(ATE = sum(CATE * mg_size, na.rm = TRUE) / sum(mg_size)) # And take the weighted average
#  return(ATE$ATE)
#}
  # for (i in 1:length(mg_sizes)) {
  #   curr_matched <- filter(matched_data, matched == i)
  #   
  # }
  # 
  # CATEs <- 
  #   FLAME_object %>% 
  #   pluck('matched_data') %>% 
  #   group_by(matched) %>% 
  #   summarise(CATE = mean(outcome[treated == 1]) - mean(outcome[treated == 0])) %>% 
  #   select(CATE)
  # 
  # # Get summary data frame with effects and size from all matched units
  CATE_df <- do.call(rbind,lapply(FLAME_object[[2]],simp))
   
  effect <- CATE_df[,which(colnames(CATE_df) == "effect")]
  size <- CATE_df[,which(colnames(CATE_df) == "size")]
  
  return(sum(effect * size)/sum(size))
}



