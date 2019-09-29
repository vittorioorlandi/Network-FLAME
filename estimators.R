n_fun <- function(z, d, Z, treated_degree) {
  sum(Z == z & treated_degree == d)
}

C_fun <- function(d, Z, treated_degree) {
  n_0d <- n_fun(0, d, Z, treated_degree)
  n_1d <- n_fun(1, d, Z, treated_degree)
  if (n_0d * n_1d == 0) {
    return(0)
  }
  return(1 / (1 / n_0d + 1 / n_1d))
}
  
strat_degree <- function(A, Z) {
  d <- colSums(A)
  treated_degree <- colSums(A * Z)
  n <- dim(A)[1]
  all_C <- vapply(1:n, function(i) {C_fun(d[i], Z, treated_degree)}, 
                  FUN.VALUE = numeric(1))
  C_weight <- all_C / sum(all_C)
  weights <- sapply(1:n, function(i) {
    a <- C_weight[i]
    b <- (2 * Z[i] - 1) / n_fun(Z[i], treated_degree[i], Z, treated_degree)
    a * b
  })
  weights
}

g <- erdos.renyi.game(n = 20, p = 0.15)
A <- get.adjacency(g, type= 'both', sparse = FALSE)
Z <- sample(c(0, 1), 20, replace = TRUE)
print(strat_degree(A, Z))