require(Rcpp)
require(RcppArmadillo)
require(igraph)
require(magrittr)
source('network_flame_sims.R')
source('FLAME_bit.R')
sourceCpp('subgraph_enumerate.cpp')
my_combn <- function(x, m) {
  if (length(x) == 1) {
    return(list(x))
  }
  return(combn(as.integer(x), m, simplify = FALSE))
}

# all_dat <- load('application.RData')

## Load in data, name the adjacency matrix A, outcome Y, treatment Z, categorical covariates X
# Navigate to where the RProject is
setwd('/Users/vittorioorlandi/Desktop/Network FLAME/Network-FLAME/')
A <- 
  read.csv('./Data/adj_allVillageRelationships_vilno_1.csv',
           header = FALSE) %>%
  as.matrix()

demographics <- read.csv('./Data/village_1.csv')

units_with_treatment_info <- demographics$adjmatrix_key
A <- A[units_with_treatment_info, units_with_treatment_info]
Y <- demographics$Y
Z <- demographics$Z
X <- demographics[, which(!colnames(demographics) %in% c('adjmatrix_key', 'Y', 'Z'))]

untreated <- which(Z == 0)

# To drop control -- control edges 
# A[untreated, untreated] <- 0 

# To drop any edges involving a control individual 
A[untreated, ] <- 0
A[, untreated] <- 0

n <- dim(A)[1]

# Brute force symmetry test because isSymmetric.matrix(A) outputs FALSE
# for (i in 1:n) {
#   for (j in 1:n) {
#     if (A[i, j] != A[j, i]) {
#       print(c(i, j))
#     }
#   }
# }

# Gives graph
G <- graph_from_adjacency_matrix(A, mode = 'undirected')
# G <- induced_subgraph(G, units_with_treatment_info)

# Enumerates all possible subgraphs and puts into dataframe
# all_subgraphs = threshold_all_neighborhood_subgraphs(G, 5)
threshold <- 0
all_subgraphs <- get_neighb_subgraphs(A, threshold)
all_features = gen_all_features(G, all_subgraphs)$feats
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
