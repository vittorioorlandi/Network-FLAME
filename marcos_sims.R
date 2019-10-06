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

# ## TEST 4
# interference_params <- list(c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(0, 1))
# 
# ## TEST 5
# interference_params <- list(c(5, 10), c(0, 1), c(0, 1), c(0, 1), c(0, 1), c(5, 10), c(5, 10))
# 
# ## TEST 6
# interference_params <- list(c(0, 1), c(0, 1), c(0, 1), c(5, 10), c(5, 10), c(0, 1), c(0, 1))

set.seed(42069)
setwd('~/Dropbox/Duke/projects/learning_interference/Network-FLAME/')
source('network_flame_sims.R')
#require(plyr)
require(ggplot2)
require(reshape2)

##################################################################################################
# Multiplicative Interference
##################################################################################################
sim_name = 'ER(50, 0.1)_degree_betweenness_triangle_multiplicative'
interference_features <- list(c('degree', 'triangle'), 
                              c('degree', 'betweenness'), 
                              c('triangle', 'betweenness'), 
                              c('triangle', 'kstar(4)'))
interference_params <- list(c(5, 1), c(5, 1), c(5, 1), c(5, 1))
out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                n_units = 50,
                n_treated = 25,
                erdos_renyi_p = 0.05,
                standardization_type = 'center',
                interference_type = 'drop_mutual_untreated_edges',
                estimators = c('true',
                               'first_eigenvector',
                               'all_eigenvectors',
                               'naive', 
                               'FLAME',
                               'stratified', 
                               'SANIA'),
                coloring = TRUE, 
                network_lik_weight = 0.5, 
                iterate_FLAME = FALSE,
                multiplicative = TRUE, 
                threshold = 2)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(interference_features = interference_features[[i]], 
                                    interference_parameters=interference_params[[i]]), 
                                    settings))
}
save(out_all, file=paste(sim_name, '.RData', sep=""))
require(beepr)
beep()


#### Plotting code
ggdata = NULL
for (i in 1:length(interference_params)){
  ggdata = rbind(ggdata, data.frame(as.data.frame(out_all[[i]]$results), setting=i))
}

ggdata = melt(ggdata, id.vars = 'setting')
levels(ggdata$variable) = c('True', 'First\n Eigenvector', 'All\n Eigenvectors',
                            'Naive',  'FLAME-Networks',  'Stratified', 'SANIA')

#this is for coloring
leq = rep(NA, length(interference_params) * 7)
i = 1
for(m in 1:length(levels(ggdata$variable))){
  for(k in 1:length(interference_params)){
    means = tapply(ggdata[ggdata$setting==k, 'value'], 
                   ggdata[ggdata$setting==k, 'variable'], mean)
    leq[i] = (means[m] <= means['FLAME-Networks'])
    i = 1 + i
  }
}
repcol = rep(leq, each=50)
ggplot(ggdata, aes(x=variable, y=value)) + 
  geom_violin(aes(fill=repcol), draw_quantiles = 0.5 ) + 
  geom_hline(data = data.frame(y = tapply(ggdata[ggdata$variable=='FLAME-Networks', 'value'], 
                                          ggdata[ggdata$variable=='FLAME-Networks', 'setting'], mean), 
                               setting=1:length(interference_params)), 
             aes(yintercept=y), linetype=2) + 
  scale_fill_brewer(palette = 'Set1') + 
  xlab('') + ylab('Mean Absolute Error') + 
  ggtitle('Simulation Results: Multiplicative Interference') + 
  facet_grid(setting ~., scales = 'free_y') +
  theme_bw() + theme(legend.position='None', plot.title = element_text(hjust = 0.5), 
                     text = element_text(color='black', size=16)) 

ggsave(paste(sim_name, '.png', sep=''), width=10, height=8, units = 'in', device = 'png', dpi=300)

##################################################################################################
# Additive Interference
##################################################################################################
sim_name = 'ER(50, 0.1)_like_presentation'
interference_features <- list(c('degree', 'triangle'),
                              c('degree', 'triangle'),
                              c('degree', 'triangle', 'kstar(2)', 'kstar(4)',
                                '3-degree-neighb', 'betweenness', 'closeness'),
                              c('triangle', 'kstar(2)', 'kstar(4)',
                                '3-degree-neighb', 'betweenness', 'closeness'))
interference_params <- list(c(0, 10),
                            c(10, 10),
                            c(-5, 1, 10, 1, 1, 1, -1),
                            c(1, 10, 1, 1, 1, -1))

out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                n_units = 50,
                n_treated = 25,
                erdos_renyi_p = 0.1,
                standardization_type = 'center',
                interference_type = 'drop_mutual_untreated_edges',
                estimators = c('true',
                               'first_eigenvector',
                               'all_eigenvectors',
                               'naive', 
                               'FLAME',
                               'stratified', 
                               'SANIA'),
                coloring = TRUE, 
                network_lik_weight = 0.5, 
                iterate_FLAME = FALSE,
                multiplicative = FALSE, 
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(
                                    interference_features = interference_features[[i]],
                                    interference_parameters = interference_params[[i]]), 
                                    settings))
}
save(out_all, file=paste(sim_name, '.RData', sep=""))
require(beepr)
beep()


#### Plotting code
ggdata = NULL
for (i in 1:length(interference_params)){
  ggdata = rbind(ggdata, data.frame(as.data.frame(out_all[[i]]$results), setting=i))
}

ggdata = melt(ggdata, id.vars = 'setting')
levels(ggdata$variable) = c('True', 'First\n Eigenvector', 'All\n Eigenvectors',
                            'Naive',  'FLAME-Networks',  'Stratified', 'SANIA')
#this is for coloring
leq = rep(NA, length(interference_params) * 7)
i = 1
for(m in 1:7){
  for(k in 1:length(interference_params)){
    means = tapply(ggdata[ggdata$setting==k, 'value'], 
                   ggdata[ggdata$setting==k, 'variable'], mean)
    leq[i] = (means[m] <= means['FLAME-Networks'])
    i = 1 + i
  }
}
repcol = rep(leq, each=50)
ggplot(ggdata, aes(x=variable, y=value)) + 
  geom_violin(aes(fill=repcol), draw_quantiles = 0.5 ) + 
  geom_hline(data = data.frame(y = tapply(ggdata[ggdata$variable=='FLAME-Networks', 'value'], 
                                          ggdata[ggdata$variable=='FLAME-Networks', 'setting'], mean), 
                               setting=1:length(interference_params)), 
             aes(yintercept=y), linetype=2) + 
  scale_fill_brewer(palette = 'Set1') + 
  xlab('') + ylab('Mean Absolute Error') + 
  ggtitle('Simulation Results: Additive Interference') + 
  facet_grid(setting ~., scales = 'free_y') +
  theme_bw() + theme(legend.position='None', plot.title = element_text(hjust = 0.5), 
                     text = element_text(color='black', size=16)) 

ggsave(paste(sim_name, '.png', sep=''), width=10, height=8, units = 'in', device = 'png', dpi=300)

##################################################################################################
# Degree vs Triangle
##################################################################################################
sim_name = 'ER(100, 0,05)_degree_triangle'
interference_features <- c('degree', 'triangle')
interference_params <- list(c(5, 0), c(4, 1), c(3, 2), c(2, 3), c(1, 4), c(0, 5))
out_all <- vector(mode = 'list', length = length(interference_params))
settings = list(sim_type = 'ER', 
                n_sims = 50,
                n_units = 75,
                n_treated = 25,
                erdos_renyi_p = 0.03,
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
                coloring = TRUE, 
                network_lik_weight = 0.5, 
                iterate_FLAME = FALSE,
                multiplicative = FALSE, 
                threshold = 5)
for (i in 1:length(interference_params)) {
  print(paste('Setting', i, 'of', length(interference_params)))
  out_all[[i]] <- list(settings)
  out_all[[i]]$results <- do.call(simulate_network_matching, 
                                  c(list(interference_parameters=interference_params[[i]]), 
                                    settings))
}
save(out_all, file=paste(sim_name, '.RData', sep=""))
require(beepr)
beep()

ggdata = NULL
for (i in 1:length(interference_params)){
  ggdata = rbind(ggdata, data.frame(as.data.frame(out_all[[i]]$results), setting=i))
}

ggdata = data.frame(mean = c(tapply(ggdata$FLAME, ggdata$setting, mean), 
                             tapply(ggdata$stratified, ggdata$setting, mean)),
                    lo = c(tapply(ggdata$FLAME, ggdata$setting, 
                                  function(x) sort(x)[length(x) * 0.25]), 
                           tapply(ggdata$stratified, ggdata$setting, 
                                  function(x) sort(x)[length(x) * 0.25])),
                    hi = c(tapply(ggdata$FLAME, ggdata$setting, 
                                  function(x) sort(x)[length(x) * 0.75]), 
                           tapply(ggdata$stratified, ggdata$setting, 
                                  function(x) sort(x)[length(x) * 0.75])),
                    setting = rep(0:5, 2), method=rep(c('FLAME', 'Stratified'), each=6))

ggplot(ggdata, aes(x=setting, y=mean, color=method, fill=method, ymin=lo, ymax=hi)) + 
  geom_ribbon(alpha=0.2, color=NA) + geom_line(size=1) +  geom_point() + 
  xlab(TeX('$\\gamma$')) + ylab('Mean Absolute Error') + 
  scale_color_brewer(palette='Set1') + scale_fill_brewer(palette='Set1') + 
  ggtitle(TeX('Simulation Results: $f$ = ($(5-\\gamma)d_i$ + $\\gamma\\Delta_i$', output='expression')) + 
  theme_bw() + theme(legend.title = element_blank(), plot.title = element_text(hjust = 0.5), 
                     text = element_text(color='black', size=16), 
                     legend.position = c(0.9,0.89), 
                     legend.background = element_rect(colour = 'black')) 
ggsave(paste(sim_name, '.png', sep=''), width=10, height=5, units = 'in', device = 'png', dpi=300)

