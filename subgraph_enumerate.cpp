// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// // [[Rcpp::export]]
// IntegerVector tester2() {
//   IntegerVector v;
//   v.push_back(3);
//   v.push_back(5);
//   v.insert(1, 10);
//   cout << v; 
//   return v;
// }

// [[Rcpp::export]]
List get_node_subgraph_counts(IntegerMatrix A, IntegerVector Z) {
  // Load relevant packages
	Environment igraph("package:igraph");
  Environment utils("package:utils");
  // Load relevant functions from those packages
  Function combn = utils["combn"];
	Function graph_from_adjacency = igraph["graph_from_adjacency_matrix"];
	Function get_neighbors = igraph["neighbors"];
	Function is_iso = igraph["is_isomorphic_to"];
  Function induced_subgraph = igraph["induced_subgraph"];
	// Get number of treated individuals
	// int n_treated = sum(Z);
	// int n_treated = std::accumulate(Z.begin(), Z.end(), 0);
	int n = A.nrow();

	// Get the graph (will need this for determining isomorphisms)
	SEXP g = graph_from_adjacency(Named("adjmatrix", A), Named("mode", "undirected"));
	// how to set vertex attributes?
	List all_subgraph_types;
	int n_seen_subgraphs; // total # non-isomorphic subgraphs seen (length of above)
	int sg_num = -1; // default value for not being isormophic to a previously seen graph

	// What we'll be returning: a list with a vector for each unit
  // corresponding to its (neighborhood) subgraph counts
	List all_nodes_subgraph_counts;

	for (int i = 1; i <= 1; ++i) { // for each unit
	  IntegerVector subgraph_counts; // Where we'll store the unit's counts
	  // Get neighbors
	  IntegerVector neighbors = get_neighbors(Named("graph", g), Named("v", i)); // can in theory just go from the adj mat
	  // for each subgraph size we're considering (< n_treated bc otherwise won't be entirely treated graph)
	  // in the case where color matters, should have j <= min(n_treated, neighbors.length())
  	for (int j = 1; j <= neighbors.length(); ++j) { // Start at 1 because we don't care about 0 node subgraphs
  	  // would be faster (and harder) computing combinations ourselves and then doing what follows one at a time
  	  List all_combs = combn(Named("x", neighbors), Named("m", j), Named("simplify", false));
  	  int n_combs = all_combs.length(); // number of subgraphs with j many nodes
  	  for (int k = 0; k < n_combs; ++k) { // for each of these j-node subgraphs
  	    sg_num = -1; // resets 'not seen' flag
  	    // For now assuming coloring and the fact that we only want fully treated graphs
  	    // Ideally would initialize this beforehand but I don't know its length...
  	    IntegerVector current_subgraph = all_combs[k]; // nodes in current subgraph
  	    // not sure if reinitializing each time is super slow; check that.
  	    SEXP neighborhood_sg = induced_subgraph(Named("graph", g), Named("vids", current_subgraph));
  	    // check that subgraph indices don't start at 0
  	    // ok temporarily don't do this and assume it's all treated
  	    // if (sum(Z[current_subgraph - 1]) != current_subgraph.length()) { // If not every node treated. -1 bc 0 based indexing.
  	    //   continue; // move on to the next j-node subgraph
  	    // }
  	    n_seen_subgraphs = all_subgraph_types.length();
  	    if (n_seen_subgraphs == 0) { // Haven't seen any graphs at all
  	      subgraph_counts.push_back(1); // Start the count at 1 for first subgraph ever seen.
  	      all_subgraph_types.push_back(neighborhood_sg); // Add this first one to the list
  	    }
  	    else { // We've seen at least 1 graph before
  	      for (int l = 0; l < n_seen_subgraphs; ++l) { // Loop through previously seen subgraphs
  	        if (is_iso(Named("graph1", all_subgraph_types[l]),
                       Named("graph2", neighborhood_sg))) { // take coloring into account!!
              sg_num = l;
  	          if (j == 2) {
  	            cout << sg_num;
  	            all_nodes_subgraph_counts.push_back(neighborhood_sg);
  	            all_nodes_subgraph_counts.push_back(all_subgraph_types[l]);
  	            return all_nodes_subgraph_counts;
  	          }
  	          break;
  	        }
  	      }
  	      if (sg_num == -1) { // haven't seen this graph before
  	        subgraph_counts.push_back(1); // Start my count at 1
            all_subgraph_types.push_back(neighborhood_sg); // Add it to the list
  	      }
          else {
            subgraph_counts[sg_num] += 1; 
          }
  	    }
  	  }
  	}
  	all_nodes_subgraph_counts.push_back(subgraph_counts);
	}
	return all_nodes_subgraph_counts;
}
