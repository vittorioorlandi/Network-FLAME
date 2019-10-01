## The file from which to run simulations
## The function to call is simulate_network_matching, which takes in arguments:
#  sim_type = 'ER': either 'ER' for Erdos-Renyi graph simulation or 'SBM' for stochastic block model generation
#  n_sims = 50: the number of simulations to run
#  n_units = 50: the number of units in the network
#   Careful making this too large as it will result in a massive number of subgraphs / computation time
#  n_blocks = 5: the number of communities when performing SBM generation
#  n_treated = floor(n_units / 2): the number of treated units
#   for now, use treat_prob (below) instead 
#  treat_prob = 0.5: the probability that each unit is treated
#   Use this instead of n_treated because I've derived the SANIA estimator in this setting 
#  erdos_renyi_p = 0.07; the p parameter in ER graph simulation.
#   Careful making this too large as it will result in a massive number of subgraphs / computation time
#  standardization_type = 'center': how to standardize feature counts when producing interference
#   Either 'center' which standardizes to mean 0, sd 1, or '0-1' which standardizes to be in [0, 1]
#  estimators: the estimators to use. Any of:
#   'true': nearest neighbor on true interference
#   'naive': difference in means
#   'all_eigenvectors': weird eigenvector Mahalanobis-like nearest neighbor thing
#   'first_eigenvectors': the above but only with the largest eigenvector
#   'FLAME': our approach
#   'stratified': the Sussman Airoldi 2017 stratified naive estimator 
#   'SANIA': the Sussman Airoldi 2017 SANIA minimum variance linear unbiased estimator
#  coloring: whether to take treatment into account when testing for isomorphism
#   for now, let's just go ahead and leave this as FALSE, though results should be similar
#  network_lik_weight: how much to weigh the maximized network log likelihood when 
#   computing match quality. This should be positive.
#   Would recommend making this 0 until I find a better way to determine it; 
#   otherwise, will probably be useless at best and counterproductive at worst
#   interference_type = 'drop_mutual_untreated_edges': determines what features
#   to count when computing / assigning interference. One of:
#   'drop_untreated_edges': drops any edge involving a control unit
#   'drop_mutual_untreated_edges': drops any edges between 2 control units
#   I recommend the second 
#  interference_features: the features to count when assigning interference. Any of:
#   'triangle' to count triangles
#   'kstar(n)' to count n-stars (sorry for ugly syntax; previously for ergm package compatibility)
#   'degree' to count degree
#   'betweenness' to compute vertex betweenness 
#   'closeness' to compute closeness centrality
#   'k-degree-neighb' to count number of neighbors with degree >= k (sorry for ugly syntax)
#  interference_parameters: the values by which to weight the above features.
#   For sim_type == 'ER':
#     A vector of numerics equal in length to length(interference_features)
#     e.g. interference_features = c('degree', 'kstar(2)') and interference_parameters = c(1, 2)
#       implies interference = 1 * degree count + 2 * 2-star count
#   For sim_type == 'SBM':
#     A list equal in length to length(interference_features)
#     Each entry of the list is a vector of length 2 with lower and upper bounds, respectively,
#       from which the weight for the respective interference feature is to be sampled uniformly
#     E.g. interference_features = c('degree', 'kstar(2)') and
#       interference_parameters = list(c(0, 1), c(2, 5)) implies
#         interference = gamma1 * degree count + gamma2 * 2-star count, where
#           gamma1 ~ U(0, 1) and gamma2 ~ U(2, 5)
#  iterate_flame = FALSE: 
#   a boolean for whether to perform FLAME (TRUE) or simply do 1-round of exact matching (FALSE)
#   the latter is much more efficient (obviously) and recommended for dense / large networks


interference_features <- list(c('kstar(2)', 'degree'),
                              c('degree', 'kstar(4)', '3-degree-neighb', 
                                'betweenness', 'closeness', 'triangle', 'kstar(2)'),
                              c('degree', 'kstar(4)', '3-degree-neighb', 
                                'betweenness', 'closeness', 'triangle', 'kstar(2)'))

interference_params <- list(c(10, 10), 
                            c(1, 1, 1, 1, 1, 1, 1), 
                            c(1, 1, 1, 10, 5, 1, 1))
interference_params <- interference_params[[1]]
interference_features <- interference_features[[1]]

# ## TEST 4
# interference_params <- list(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1))
# 
# ## TEST 5
# interference_params <- list(c(5, 10), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(5, 10), c(5, 10))
# 
# ## TEST 6
# interference_params <- list(c(0, 1), c(0, 1), c(0, 1), c(5, 10), c(5, 10), c(0, 1), c(0, 1))

er_p <- c(.05, .07, .1)
interference_features <- list(c('kstar(2)', 'degree'),
                              c('kstar(2)', 'degree'),
                              c('degree', 'kstar(4)', '3-degree-neighb', 
                                'betweenness', 'closeness', 'triangle', 'kstar(2)'),
                              c('degree', 'kstar(4)', '3-degree-neighb', 
                                'betweenness', 'closeness', 'triangle', 'kstar(2)'),
                              c('degree', 'kstar(4)', '3-degree-neighb', 
                                'betweenness', 'closeness', 'triangle', 'kstar(2)'))

interference_params <- list(c(1, 1), 
                            c(1, -1),
                            c(-1, 1, 1, 1, 1, -1, 1), 
                            c(1, 1, 1, 10, -5, 1, 1), 
                            c(-5, 1, 1, 10, -5, 1, 10))
interference_features <- c('triangle', 'degree')
interference_params <- c(10, 0)
simulate_network_matching(sim_type = 'ER', 
                            n_sims = 50,
                            n_units = 50,
                            erdos_renyi_p = 0.07,
                            standardization_type = 'center',
                            interference_type = 'drop_mutual_untreated_edges',
                            estimators = c('true',
                                           'first_eigenvector',
                                           'all_eigenvectors',
                                           'FLAME',
                                           'naive', 
                                           'stratified', 
                                           'SANIA'),
                            interference_features = interference_features,
                            interference_parameters = interference_params,
                            coloring = FALSE, 
                            network_lik_weight = 0, 
                          iterate_FLAME = FALSE)

require(beepr)
beep()

