require(Matrix)
require(igraph)
require(dplyr)
require(magrittr)
require(stringr)

standardize <- function(x, type = '0-1', limits = c(0, 1)) {
  if (type == 'center') {
    return((x - mean(x, na.rm = T)) / (ifelse(sd(x) == 0, 1, sd(x))))
  } else if (type == '0-1') {
    if (length(unique(x)) == 1) {return(x)}
    return((limits[2] - limits[1]) / (max(x) - min(x)) * (x - max(x)) + limits[2])
  } else {
    stop("I don't know how to standardize this!")
  }
}

count_degree = function(G){
  dm = max(degree(G))
  X = matrix(0, length(V(G)), dm)
  for(i in V(G)){
    isg = induced_subgraph(G, neighbors(G, i))
    dgs = table(degree(isg))
    X[i, as.integer(names(dgs)) + 1] = dgs
  }
  X[, which(colSums(X) > 0)]
}

count_feature <- function(A, G, feature) {
  # Takes in:
  #   The (treated) adjacency matrix A
  #   The (treated) graph G
  #   A string, feature, specifying which feature should be counted
  # Valid features are:
  #   'betweenness', 'closeness', 'degree', 'triangle',
  #   'kstar(n)' (which finds n-stars)
  # Returns:
  #   A vector of length n_units with the value of the feature for each unit
  n <- vcount(G)
  if (feature == 'betweenness') {
    return(betweenness(G, directed = FALSE, normalized = TRUE))
  }
  if (feature == 'closeness') {
    return(closeness(G, normalized = TRUE))
  }
  if (feature == 'degree') {
    return(rowSums(A))
  }
  else {
    feature_counts <- vector(mode = 'numeric', length = n)
  }
  for (i in 1:n) {
    isg <- induced_subgraph(G, neighbors(G, i))
    if (feature == 'triangle') {
      feature_counts[i] <- length(E(isg))
    } else if (str_detect(feature, 'kstar')) {
        feature_counts[i] <- choose(length(V(isg)),
                                  as.numeric(str_extract(feature, '[:digit:]+')))
    } else if (str_detect(feature, 'degree-neigh')) {
        degree_thresh <- str_sub(feature, 1, (str_locate(feature, '-')[1] - 1))
        feature_counts[i] <- sum(degree(isg) >= degree_thresh)
    } else {
      stop("I don't know what feature this is!")
    }
  }
  return(feature_counts)
}

threshold_list_subgraphs = function(V, k) {
  if (k == 'max') {
    k <- min(length(V), k)
  }
  else {
    k <- min(length(V), k)
  }
  sgs = c()
  for (j in 1:k) {
    if (j == 1 && length(V) == 1) {
      sgs = c(sgs, list(V)) # Can probably simplify
    } else if (j <= length(V)) {
      sgs = c(sgs, combn(V, j, simplify = FALSE))
    }
  }
  # for (v in combn(V, k, simplify = FALSE)){
  #   sgs = c(sgs, list.flatten(list_subgraphs(v)))
  # }
  sgs
}

threshold_all_neighborhood_subgraphs = function(G, k = 'max') {
  # For each vertex i, generates sgs[[i]] which is all possible combinations of
  # ≤ k vertices out of i's neighbors
  # k <- sort(ego_size(G), decreasing = TRUE)[2]
  sgs = list()
  for (i in V(G)) {
    if (i%%10 == 0) {
      print(i)
    }
    neighbs <- neighbors(G, i)
    if (length(neighbs) == 0)
      sgs[[i]] = numeric(0)
    else
      sgs[[i]] = threshold_list_subgraphs(neighbs, k)
  }
}

find_sg_type = function(G, sg, sg_types) {
  # Given a subgraph sg and sg_types, returns
  # i if sg is isomorphic to the i'th subgraph in sg_types
  # and NA otherwise
  if (length(sg_types) == 0)
    return(NA)
  sG = induced_subgraph(G, sg)
  coloring <- vertex_attr(sG, 'color')
  # if (!is.null(coloring) & sum(coloring) == 0) { # There's at least one untreated individual
  #   return(-1)
  # }
  # speeds things up, if possible, when not coloring
  method <- ifelse(is.null(coloring), 'auto', 'vf2')
  for (i in 1:length(sg_types)) {
    if (is_isomorphic_to(sG, sg_types[[i]], method = method))
      return(i)
  }
  return(NA)
}

gen_features = function(G, sgs, sg_types = list()) {
  # sgs is a list, each entry of which is a subset of vertices in the
  # neighborhood of a single vertex in G
  feats = list()

  for (g in sgs) {
    type = find_sg_type(G, g, sg_types)

    # If we haven't yet seen this type of subgraph...
    if (is.na(type)) {
      # Introduce a new type and classify it as such
      type = length(sg_types) + 1
      sg_types[[type]] = induced_subgraph(G,g)
      feats[[type]] = 0
    }
    if (type == -1) { # wholly untreated subgraph
      next
    }
    if (type > length(feats) || is.null(feats[[type]]))
      feats[[type]] = 0

    # One more observation of this subgraph type
    feats[[type]] = feats[[type]] + 1
  }
  return(list(feats, sg_types))
}

gen_all_features = function(G, sgs) {
  sg_types = list()
  feats = list()
  for (i in 1:length(sgs)) { # For each unit
    res = gen_features(G, sgs[[i]], sg_types)
    sg_types = res[[2]]
    feats[[i]] = res[[1]]
  }
  return(list(feats = feats,
              types = sg_types))
}

gen_data = function(feats) {
  dta = NULL
  n_feats = max(sapply(feats, length))
  for (i in feats){
    lst = i
    if(length(lst) < n_feats)
      lst[[n_feats]] = 0
    vs = unlist(lapply(lst, function(x) ifelse(is.null(x), 0, x)))
    dta = rbind(dta, vs)
  }
  data.frame(dta, row.names = 1:length(feats))
}

new_sg_counter <- function(g, threshold) {
  n <- length(V(g))
  all_neighb_subgraphs <- vector(mode = 'list', length = n)
  for (i in 1:n) {
    tmp <- gen_connected_sgs(g, neighbors(g, i), sgs_so_far = c(), c(), threshold)
    all_neighb_subgraphs[[i]] <- tmp[2:length(tmp)]
  }
  return(all_neighb_subgraphs)
}

ATE_est_error <- function(Y, f, estimator_type, G, A, ATE, Z,
                          feature_counts, network_lik_weight, treat_prob, iterate_FLAME, threshold) {
  # Takes in:
  #  The outcome vector Y
  #  The interference vector f
  #  The estimator type, one of:
  #    'true', 'first_eigenvector', 'all_eigenvectors', 'degree_dist', 'naive', 'FLAME'
  #  The graph G
  #  The (true) average treatment effect ATE
  # Returns absolute error between the estimated ATE and the true value, a scalar

  n <- length(Y)
  total_treatment_effect <- 0
  if (estimator_type == 'true') {
    error <- abs(true_est(Y, Z, f) - ATE)
  } else if (estimator_type == 'stratified') {
    error <- abs(strat_naive(A, Z, Y) - ATE)
  } else if (estimator_type == 'SANIA') {
    error <- abs(ind_SANIA(G, Z, Y, mean(Z)) - ATE)
  } else if (estimator_type == 'first_eigenvector') {
    error <- abs(first_eig(A, Z, Y) - ATE)
  } else if (estimator_type == 'all_eigenvectors') {
    error <- abs(all_eigs(A, Z, Y) - ATE)
  } else if (estimator_type == 'FLAME') {

    # all_subgraphs <- threshold_all_neighborhood_subgraphs(G, threshold)
    all_subgraphs <- get_neighb_subgraphs(G, threshold)
    # all_subgraphs <- new_sg_counter(G, threshold)
    features_and_graphs <- gen_all_features(G, all_subgraphs)
    all_features <- features_and_graphs$feats
    sg_types <- features_and_graphs$types
    dta <- gen_data(all_features)
    dta <- data.frame(sapply(dta, factor), stringsAsFactors = T)
    # dta = dta[, order(sapply(dta, function(x) length(levels(x))), decreasing = T)]
    # dta <- data.frame(sapply(dta, factor), stringsAsFactors = T)
    dta$outcome <- Y
    dta$treated <- factor(Z)
    # browser()
    flame_out <- FLAME_bit(dta, dta, A = A,
                           network_lik_weight = network_lik_weight, iterate_FLAME = iterate_FLAME)
    error <- abs(ATE(flame_out) - ATE)
  } else if (estimator_type == 'FLAM') { #ignore me

    matching_features <- colnames(feature_counts)

    # Make the features into a dataframe suitable for FLAME
    # including making covariates into factors
    feature_counts %<>%
      apply(2, factor) %>%
      data.frame(stringsAsFactors = T) %>%
      `colnames<-`(matching_features)

    # Could be unnecessary
    feature_counts <-
      feature_counts[, order(sapply(feature_counts, function(x) {
        length(levels(x))}),
        decreasing = T)]

    # Add numeric outcome for FLAME
    feature_counts$outcome <- Y

    # Add factor treatment for FLAME
    feature_counts$treated <- factor(Z)
    # browser()
    # save('feature_counts', file = 'data_causing_problems.RData')
    flame_out <- FLAME_bit(feature_counts, feature_counts) # pass in adj mat
    error <- abs(ATE(flame_out) - ATE)
  } else if (estimator_type == 'degree_dist') {
    X <- count_degree(G)
    dta <- data.frame(apply(X, 2,  as.character))
    dta = dta[, order(sapply(dta, function(x) length(levels(x))), decreasing = T)]
    dta$outcome = Y
    dta$treated = factor(Z)
    flame_res = FLAME_bit(dta, dta, AZ)
    error <- abs(ATE(flame_res) - ATE)
  } else if (estimator_type == 'naive') {
    error <- abs(mean(Y[Z == 1]) - mean(Y[Z == 0]) - ATE)
  } else {
    stop("I don't know what this estimator is!")
  }
  return(error)
}

simulate_network_matching <- function(sim_type = 'ER',
                                      n_sims = 100,
                                      n_units = 100,
                                      treat_prob = 0.5,
                                      interference_type = 'drop_mutual_untreated_edges',
                                      n_treated = floor(n_units / 2),
                                      erdos_renyi_p = 0.07,
                                      standardization_type = '0-1',
                                      interference_features,
                                      interference_parameters = NULL,
                                      matching_features = c(),
                                      estimators,
                                      network_lik_weight = 1,
                                      coloring = TRUE,
                                      n_blocks = 5,
                                      iterate_FLAME = FALSE,
                                      multiplicative = FALSE,
                                      threshold) {

  # interference_features are the features that should be used to construct
  # the interference function, their prevalences weighted by the corresponding
  # entries in interference_parameters
  # matching_features are the features to be matched on via FLAME
  # valid features include:
  #   'degree', 'kstar(n)', which denotes an n-star, 'triangle'

  if (network_lik_weight < 0) {
    stop('network_lik_weight should be positive, as it weighs the (negative) log-likelihood')
  }
  # browser()
  stopifnot(length(interference_features) == length(interference_parameters))
  if ('degree_dist' %in% matching_features & length(matching_features) > 1) {
    stop('If matching on degree sequence, cannot match on other features as well')
  }

  # To store errors of the estimators
  abs_error <- matrix(, nrow = n_sims, ncol = length(estimators))
  colnames(abs_error) <- estimators

  # Parameters used to generate outcome
  alpha <- rnorm(n_units, 0, 1)
  beta <- rnorm(n_units, 5, 1)
  ATE <- mean(beta)

  all_features <- base::union(interference_features, matching_features)
  n_features <- length(all_features)

  if (sim_type == 'SBM') {
    interference_parameters %<>%
      lapply(function(x) runif(n_blocks, x[1], x[2])) %>%
      unlist() %>%
      rep(each = n_units / n_blocks) %>%
      matrix(nrow = n_units)
  }
  pmat = Matrix(runif(n_blocks ^ 2, 0, 0.1),
                nrow = n_blocks,
                ncol = n_blocks)
  pmat = forceSymmetric(pmat)
  bsizes = rep(n_units / n_blocks, n_blocks)

  for (sim in 1:n_sims) {

    # Generate an undirected Erdos-Renyi or SBM graph
    if (sim_type == 'ER') {
      G <- erdos.renyi.game(n_units, erdos_renyi_p, directed = FALSE)
    } else if (sim_type == 'SBM') {
      G <- sbm.game(n_units, pmat, bsizes, directed = F, loops = F)
    } else {
      stop("I don't recognize this type of simulation!")
    }

    if (sim %% 25 == 0) {
      message(sim)
    }

    # Randomly assign treatment to n_treated units
    # Z <- rep(0, n_units)
    # Z[sample(1:n_units, n_treated)] = 1

    # Assign treatment randomly with probability treat_prob
    Z <- sample(c(0, 1), size = n_units, replace = TRUE, prob = c(1 - treat_prob, treat_prob))
    n_treated <- sum(Z)

    if (coloring) {
      G <- set_vertex_attr(G, 'color', value = Z)
    }

    # Get the graph adjacency matrix
    A <- get.adjacency(G, type= 'both', sparse = FALSE)

    feature_counts <- matrix(0, nrow = n_units, ncol = n_features)
    colnames(feature_counts) <- all_features

    if (interference_type == 'drop_mutual_untreated_edges') {
      AZ <- A
      untreated <- which(Z == 0)
      to_drop <- combn(untreated, 2)
      for (i in 1:ncol(to_drop)) {
        AZ[to_drop[1, i], to_drop[2, i]] <- 0
      }
      AZ %<>% forceSymmetric()
    } else if (interference_type == 'drop_untreated_edges') {
      untreated <- which(Z == 0)
      AZ <- A
      AZ[untreated, ] <- 0
      AZ[, untreated] <- 0
    } else {
      print('Assigning interference to the whole, untouched graph!')
      AZ <- A
    }
    GZ <- graph_from_adjacency_matrix(as.matrix(AZ), mode = 'undirected')
    for (i in 1:n_features) {
      feature_counts[, i] <- count_feature(AZ, GZ, all_features[i])
    }
    std_feature_counts <- apply(feature_counts, 2, standardize, standardization_type)
# check MARGIN!!
    if (sim_type == 'ER') {
      f <- sweep(std_feature_counts[, interference_features],
                 MARGIN = 2,
                 STATS = interference_parameters,
                 FUN = '*')
    } else {
      f <- std_feature_counts[, interference_features] * interference_parameters
    }

    # f <- rowSums(f) + 7 * ATE
    # ff <- rowSums(f)
    f <- rowSums(f)
    # browser()
    if (multiplicative) {
      tmp <- std_feature_counts[, interference_features]
      cts <- tmp[, 1] * tmp[, 2]
      f <- prod(interference_parameters) * cts
    }
    # f <- standardize(Z * (rowSums(f) + 7 * ATE), limits = c(0, mean(beta)))
    Y = alpha + beta * Z + f

    for (j in 1:length(estimators)) {
      abs_error[sim, j] <-
        ATE_est_error(Y, f, estimators[j], G, A, ATE, Z,
                      feature_counts, network_lik_weight, treat_prob, iterate_FLAME, threshold)
    }
  }

  mean_abs_error <- colMeans(abs_error)
  sd_abs_error <- apply(abs_error, 2, sd)

  # print(sprintf('For %d units, %s coloring....', n_units, c('without', 'with')[1 + coloring]))
  print(sprintf('With p = %.2f,', erdos_renyi_p))
  print(interference_parameters)
  print(interference_features)
  for (i in seq_along(estimators)) {
    'The mean absolute error of the %s estimator was %.2f, with a standard deviation of %.2f' %>%
      sprintf(estimators[i], mean_abs_error[i], sd_abs_error[i]) %>%
      print()
  }
  return(abs_error)
}
#

# interference_features <- list(c('kstar(2)', 'degree'),
#                               c('degree', 'kstar(4)', '3-degree-neighb',
#                                 'betweenness', 'closeness', 'triangle', 'kstar(2)'),
#                               c('degree', 'kstar(4)', '3-degree-neighb',
#                                 'betweenness', 'closeness', 'triangle', 'kstar(2)'))
#
# interference_params <- list(c(10, 10),
#                             c(1, 1, 1, 1, 1, 1, 1),
#                             c(1, 1, 1, 10, 5, 1, 1))
# interference_params <- interference_params[[1]]
# interference_features <- interference_features[[1]]
#
# # ## TEST 4
# # interference_params <- list(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1))
# #
# # ## TEST 5
# # interference_params <- list(c(5, 10), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(5, 10), c(5, 10))
# #
# # ## TEST 6
# # interference_params <- list(c(0, 1), c(0, 1), c(0, 1), c(5, 10), c(5, 10), c(0, 1), c(0, 1))
# er_p <- c(.05, .07, .1)
# interference_features <- list(c('kstar(2)', 'degree'),
#                               c('kstar(2)', 'degree'),
#                               c('degree', 'kstar(4)', '3-degree-neighb',
#                                 'betweenness', 'closeness', 'triangle', 'kstar(2)'),
#                               c('degree', 'kstar(4)', '3-degree-neighb',
#                                 'betweenness', 'closeness', 'triangle', 'kstar(2)'),
#                               c('degree', 'kstar(4)', '3-degree-neighb',
#                                 'betweenness', 'closeness', 'triangle', 'kstar(2)'))
#
# interference_params <- list(c(1, 1),
#                             c(1, -1),
#                             c(-1, 1, 1, 1, 1, -1, 1),
#                             c(1, 1, 1, 10, -5, 1, 1),
#                             c(-5, 1, 1, 10, -5, 1, 10))
#
# for (i in seq_along(er_p)) {
#   for (j in seq_along(interference_features)) {
#     simulate_network_matching(sim_type = 'ER',
#                               n_sims = 50,
#                               n_units = 50,
#                               erdos_renyi_p = er_p[i],
#                               standardization_type = 'center',
#                               interference_type = 'drop_mutual_untreated_edges',
#                               estimators = c('true',
#                                              'first_eigenvector',
#                                              'all_eigenvectors',
#                                              'FLAME',
#                                              'naive',
#                                              'stratified',
#                                              'SANIA'),
#                               interference_features = interference_features[[j]],
#                               interference_parameters = interference_params[[j]],
#                               coloring = FALSE,
#                               network_lik_weight = 0,)
#   }
# }
#
# require(beepr)
# beep()
