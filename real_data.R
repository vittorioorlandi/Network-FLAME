require(igraph)
source('network_flame_sims.R')
source('FLAME_bit.R')

# Gives graph
G <- graph_from_adjacency_matrix(as.matrix(A), mode = 'undirected')
G <- induced_subgraph(G, vertices_with_treatment_info)

# Enumerates all possible subgraphs and puts into dataframe
all_subgraphs = threshold_all_neighborhood_subgraphs(G, 'max')
all_features = gen_all_features(G, all_subgraphs)
dta = gen_data(all_features)

# Add covariate information
## Check that order of covariates X is same as order of subgraph counts
dta <- cbind(dta, X)

# Convert everything to factor 
dta = data.frame(sapply(dta, factor), stringsAsFactors = T)
dta = dta[, order(sapply(dta, function(x) length(levels(x))), decreasing = T)]
dta <- data.frame(sapply(dta, factor), stringsAsFactors = T)

# Adds outcome and treatment
dta$outcome = Y # Y should be numeric
dta$treated = factor(Z) # Z should be binary vector

flame_out <- FLAME_bit(dta, dta, A = A, network_lik_weight = 0, iterate_FLAME = TRUE)
ATE_out <- ATE(flame_out)