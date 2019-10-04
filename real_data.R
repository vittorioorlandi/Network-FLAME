require(Rcpp)
require(RcppArmadillo)
require(igraph)
require(magrittr)
source('network_flame_sims.R')
source('FLAME_bit.R')
source('ATE.R')

my_combn <- function(x, m) {
  if (length(x) == 1) {
    return(list(x))
  }
  return(combn(as.integer(x), m, simplify = FALSE))
}
sourceCpp('subgraph_enumerate.cpp')

# all_dat <- load('application.RData')

## Load in data, name the adjacency matrix A, outcome Y, treatment Z, categorical covariates X
# Navigate to where the RProject is
#setwd('/Users/vittorioorlandi/Desktop/Network FLAME/Network-FLAME/')
setwd("/Users/musaidawan/Dropbox/Duke/Projects/Data for Network paper/Network-FLAME/")
village_codes <- setdiff(c(1:77), c(13,22))
outcome <- matrix(NA, nrow = 75, ncol = 3)

for (qs in c(1,2,3)) {
  for (val in village_codes) {
    
    qs <- 1
    val <- 1
    print(val)
    A <-
      read.csv(paste('./Data/Adjency/adj_andRelationships_vilno_',val,'.csv', sep = ""),
               header = FALSE) %>%
      as.matrix()

    demographics <- read.csv(paste('./Data/characteristics_',qs,'/village_',val,'.csv',sep =""))

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

    # Enumerates all possible subgraphs and puts into dataframe
    ## Old method; avoid using if you can
    #all_subgraphs <- threshold_all_neighborhood_subgraphs(G, 3)

    all_subgraphs <- get_neighb_subgraphs(A)
    all_features = gen_all_features(G, all_subgraphs)
    dta = gen_data(all_features)


    # Add covariate information
    ## Check that order of covariates X is same as order of subgraph counts
    dta <- cbind(dta, X)

    # Convert everything to factor
    dta = data.frame(sapply(dta, factor), stringsAsFactors = T)

    # Adds outcome and treatment
    dta$outcome = Y # Y should be numeric
    dta$treated = factor(Z) # Z should be binary vector

    # drop cols with no variation
    drop_these <- which(lapply(dta, function(x) length(unique(x))) == 1)
    tmp <- dta[, -drop_these]

    # cols with missing values: // should be no missing values, data pre-processing
    # lapply(tmp, function(x) sum(is.na(x)) == 1)

    # impute missing val by median, lower
    #tmp$ration_color_2[is.na(tmp$ration_color_2)] <- as.factor(median(as.numeric(tmp$ration_color_2),na.rm = TRUE))

    # FLAME
    flame_out <- FLAME_bit(tmp, tmp, A = A, network_lik_weight = 0, iterate_FLAME = TRUE)
    ATT_out <- flame_out$ATT
    if (val == 1) {
      covs_list_out <- flame_out$covs_list
    }
    outcome[val][qs] <- ATT_out
  }
}
