#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]
using namespace Rcpp;

// #include <iostream>
// #include <vector>
// 
// using namespace std;
// 
// vector<int> people;
// vector<int> combination;
// 
// void go(int offset, int k) {
//   if (k == 0) {
//     return combination;
//   }
//   for (int i = offset; i <= people.size() - k; ++i) {
//     combination.push_back(people[i]);
//     go(i+1, k-1);
//     combination.pop_back();
//   }
// }
// 
// int main() {
//   int n = 5, k = 3;
//   
//   for (int i = 0; i < n; ++i) { people.push_back(i+1); }
//   go(0, k);
//   
//   return 0;
// }



// [[Rcpp::export]]
List nth_neighbor(int j, IntegerMatrix A, IntegerVector Z) {
  // Load relevant packages
	Environment igraph("package:igraph");
  Environment utils("package:utils");
  // Load relevant functions from those packages
  Function combn = utils["combn"];
	Function graph_from_adjacency = igraph["graph_from_adjacency_matrix"];
	Function get_neighbors = igraph["neighbors"];
	Function is_iso = igraph["is_isomorphic_to"];
	
	// Get number of treated individuals
	int n_treated = sum(Z);
	int n = A.nrow();
	
	// Get the graph (will need this for determining isomorphisms)
	SEXP g = graph_from_adjacency(Named("adjmatrix", A), Named("mode", "undirected"));
	// how to set vertex attributes? 
	List all_subgraph_types = List::create();
	int sg_num = -1; // default value for not being isormophic to a previously seen graph
	int n_seen_subgraphs; // total # non-isomorphic subgraphs seen 
	int counts_so_far; // defined later 
	List all_nodes_subgraph_counts = List::create();
	for (int i = 0; i < n; ++i) { // for each unit
	  IntegerVector subgraph_counts (); 
	  // Get neighbors 
	  IntegerVector neighbors = get_neighbors(Named("graph", g), Named("v", i)); // can in theory just go from the adj mat
	  // for each subgraph size we're considering (< n_treated bc otherwise won't be entirely treated graph)
  	for (int j = 1; j < n_treated; ++j) { // Start at 1 because we don't care about 0 node subgraphs
  	  // would be faster (and harder) computing combinations ourselves and then doing what follows one at a time
  	  List all_combs = combn(Named("x", neighbors), Named("m", j), Named("simplify", false));	
  	  int n_combs = all_combs.length()
  	  for (int k = 0; k < n_combs; ++k) {
  	    // For now assuming coloring and the fact that we only want fully treated graphs
  	    IntegerVector current_subgraph = all_combs[[k]]; // nodes in current subgraph
  	    if (sum(Z[current_subgraph]) != current_subgraph.length()) { // If not every node treated
  	      continue;
  	    }
  	    n_seen_subgraphs = all_subgraph_types.length(); 
  	    if (n_seen_subgraphs == 0) {
  	      sg_num = 0; // First graph will go at 0th position of all_subgraph_types
  	      subgraph_counts.insert(sg_num, 1) // Start the count at 1 for first subgraph ever seen
  	    }
  	    else {
  	      neighborhood_sg = induced_subgraph(Named("graph", g), Named("vids", current_subgraph))
  	      for (l = 0; l < n_seen_subgraphs; ++l) {
  	        if (is_iso(Named("graph1", all_subgraph_types[[l]]), 
                       Named("graph2", neighborhood_sg), 
                       Named("method", "vf2"))) {
              sg_num = l;
  	          break;
  	        }
  	      }
  	      if (sg_num == -1) { // haven't seen this graph before
            all_subgraph_types.push_back(neighborhood_sg); // Add it to the list
  	        subgraph_counts.push_back(1); // Start my count at 1
  	      }
          else {
            counts_so_far = subgraph_counts[sg_num]; // of this subgraph
            all_subgraph_types.insert(sg_num, counts_so_far + 1); 
          }
  	    }
  	  }
  	}
  	all_nodes_subgraph_counts.push_back(subgraph_counts);
	}
	return all_nodes_subgraph_counts;
}

/*** R
set.seed(20130810)
library(igraph)

A <- matrix(c(0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0), nrow = 5)
nth_neighbor(1, A)
*/