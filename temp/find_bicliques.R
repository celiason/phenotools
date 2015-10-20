# find bicliques

# temporary function

samp <- sample(1:1000, 10)

M <- tmp$data[1:10, samp]

M
library(igraph)

g <- graph_from_incidence_matrix(M)

plot(g, layout=layout.bipartite)


# function to compute maximum biclique

# Algorithm: MBEA (Starred lines apply to iMBEA)
# procedure biclique_find(G, L, R, P, Q);
# G: a bipartite graph G = (U ∪ V , E);
# L: set of vertices ∈ U that are common neighbors of vertices in R, initially L = U;
# R: set of vertices ∈ V belonging to the current biclique, initially empty;
# P: set of vertices ∈ V that can be added to R, initially
# P = V, sorted by non-decreasing order of common neighborhood size;
# Q: set of vertices used to determine maximality, initially empty;
* i ← 0; // Position of selected candidate
in P

biclique_find <- function(igraph)
# v = vertices with type
# u = vertices matching to?
# edges
# vertices set 1
L <- V(igraph)
R <- NULL
# vertices set 2
P <- 
Q <- NULL



while P ̸= ∅ do
* x ← P[i + +]; // Select next candidate
from P in order
R′ ←R∪{x};
L′ ←{u∈L|(u,x)∈E(G)};
       // Observation 1: extend biclique
* L′ ← L \ L′; C ← {x};
P′ ← ∅; Q′ ← ∅;// Create new sets
// Observation 2: check maximality is_maximal ← TRUE;
forall the v in Q do
N[v]← {u ∈ L′ | (u,v) ∈ E(G)};
// Observation 4: end of branch if |N[v] | = |L′| then
is_maximal ← FALSE; break;
elseif|N[v]|>0thenQ′ ←Q′∪{v};
if is_maximal = TRUE then forall the v in P, v ̸= x do
N[v]←{u∈L′ |(u,v)∈E(G)};//Getthe neighbors of v in L′
if |N[v] | = |L′| then R′ ← R′ ∪ {v};
                // Observation 3: expand
to maximal
*S←{u∈L′ |(u,v)∈E(G)}; * if |S| = 0 then C ← C ∪ {v};
                // Observation 5: further
                pruning
* else if |N[v]| > 0 then
* P′ ← P′ ∪ {v} // Insert v into P′
                in non-decreasing order of
                common neighborhood size
PRINT(L′, R′); // Report maximal biclique
if P′ ̸= ∅ then biclique_find(G, L′, R′, P′, Q′);
// Move C from candidate set to former candidate set
* Q ← Q ∪ C; P ← P \ C;


