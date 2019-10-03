# gen_connected_sgs <- function(g, v_not_considered, sgs_so_far, neighbors, threshold) {
#   if (length(sgs_so_far) >= threshold) {
#     return(list(sgs_so_far))
#   }
#   if (length(sgs_so_far) == 0) {
#     candidates <- v_not_considered
#   }
#   else {
#     candidates <- intersect(v_not_considered, neighbors)
#   }
#   if (length(candidates) == 0) {
#     return(list(sgs_so_far))
#   } 
#   else {
#     v <- candidates[1]
#     v_not_considered <- setdiff(v_not_considered, v)
#     return(c(gen_connected_sgs(g, v_not_considered,
#                                sgs_so_far,
#                                neighbors, threshold),
#              gen_connected_sgs(g, v_not_considered, 
#                                union(sgs_so_far, v),
#                                union(neighbors, neighbors(g, v)), threshold)))
#   }
# }
# 
# threshold_list_subgraphs = function(V, k) {
#   if (k == 'max') {
#     k <- min(length(V), k)
#   }
#   else {###### I know this is stupid. 
#     k <- min(length(V), k)
#   }
#   sgs = c()
#   for (j in 1:k) {
#     if (j == 1 && length(V) == 1) {
#       sgs = c(sgs, list(V)) # Can probably simplify
#     } else if (j <= length(V)) {
#       sgs = c(sgs, combn(V, j, simplify = FALSE))
#     }
#   }
#   sgs
# }
# 
# threshold_all_neighborhood_subgraphs = function(G, k = 'max') {
#   # For each vertex i, generates sgs[[i]] which is all possible combinations of
#   # ≤ k vertices out of i's neighbors
#   # k <- sort(ego_size(G), decreasing = TRUE)[2]
#   sgs = list()
#   for (i in V(G)) {
#     # if (i %% 10 == 0) {print(i)}
#     neighbs <- neighbors(G, i)
#     if (length(neighbs) == 0)
#       sgs[[i]] = numeric(0)
#     else
#       sgs[[i]] = threshold_list_subgraphs(neighbs, k)
#   }
#   sgs
# }
# 
# g <- erdos.renyi.game(80, 0.15)
# neighbs <- neighborhood(g)
# n <- length(V(g))
# 
# start1 <- Sys.time()
# all_neighb_subgraphs <- vector(mode = 'list', length = n)
#   for (i in 1:n) {
#     tmp <- gen_connected_sgs(g, neighbors(g, i), sgs_so_far = c(), c(), threshold = 5)
#     all_neighb_subgraphs[[i]] <- tmp[2:length(tmp)]
#   }
#   # all_feats <- gen_all_features(g, all_neighb_subgraphs)
#   # dta <- gen_data(all_feats)
# print(Sys.time() - start1)
# 
# start2 <- Sys.time()
# out <- threshold_all_neighborhood_subgraphs(g, 5)
# print(Sys.time() - start2)
#   # dta2 <- gen_data(gen_all_features(g, threshold_all_neighborhood_subgraphs(g, 5)))
# beep()

my_combn <- function(x, m) {
  if (length(x) == 1) {
    return(list(x))
  }
  return(combn(x, m, simplify = FALSE))
}

require(Rcpp)
require(RcppArmadillo)
require(igraph)
require(magrittr)
sourceCpp('subgraph_enumerate.cpp')
G <- erdos.renyi.game(20, 0.07)
# A <- matrix(c(0, 0, 0, 0, 0, 
#               0, 0, 0, 1, 1, 
#               0, 0, 0, 1, 0, 
#               0, 1, 1, 0, 1, 
#               0, 1, 0, 1, 0), 
#             nrow = 5)
A <- get.adjacency(G, type = 'both', sparse = FALSE)
n <- dim(A)[1]
# Z <- c(1, 1, 1, 1, 1)
Z <- rep(1, n)
out <- get_node_subgraph_counts(A, Z)
max_len <- max(vapply(out, length, numeric(1)))
out %<>% 
  sapply(out, function(x) c(x, rep(0, max_len - length(x)))) %>%
  as.data.frame()

system.time({
  G <- erdos.renyi.game(25, 0.1)
  A <- get.adjacency(G, type = 'both', sparse = FALSE)
  for (i in 1:500) {
    n <- dim(A)[1]
    Z <- rep(1, n)
    dta <- get_node_subgraph_counts(A, Z)
  }
})
beep()












G <- erdos.renyi.game(29, 0.1)
A <- get.adjacency(G, type = 'both', sparse = FALSE)
n <- dim(A)[1]
Z <- rep(1, n)
dta <- get_node_subgraph_counts(A, Z)















