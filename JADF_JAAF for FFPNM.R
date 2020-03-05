
###################################################################### 
# The following functions require R libraries "network".
# One way to acquire these libraries is to install the "statnet"         
# suite of packages (www.statnetproject.org). 
###################################################################### 
# Function to compute the feature "JADF" discussed in the paper.
# G is a network object
# node1 denotes the ID of a drug
# node2 denotes the ID of an ADE

# This function returns a vector denotes the value of feature "JADF" 
#####################################################################

compute.jaccard.aa.drug.feature.sum = function(G, node1, node2) {
    N1 = get.neighborhood(G, node1)
    N2 = get.neighborhood(G, node2)
    N2 = setdiff(N2, node1)           
    n.neighbors = length(N2)
    D = degree(G, gmode="graph")
    D.node1 = D[node1]
    result.ja.drug = numeric(3)
    result.ja.drug[1:3] <- NA
    names(result.ja.drug) = c('ja_drug_min','ja_drug_max','ja_drug_sum')
    if (n.neighbors == 0) {
       result.ja.drug['ja_drug_min'] = 0
       result.ja.drug['ja_drug_max'] = 0
       result.ja.drug['ja_drug_mean'] = 0
    } else {
        ja.vector = numeric(n.neighbors)
        for (i in 1:n.neighbors) {
        neighbor.i = N2[i]       
        N1.i = get.neighborhood(G, neighbor.i)
        intersection.i = intersect(N1, N1.i)
        union.i = union(N1, N1.i)
        D.i=D[neighbor.i]
        ja.vector[i] = (length(intersection.i)/length(union.i))/log(D.i)+1/log(D.i)
        }
        result.ja.drug['ja_drug_min'] = min(ja.vector)
        result.ja.drug['ja_drug_max'] = max(ja.vector)
        result.ja.drug['ja_drug_sum'] = sum(ja.vector)
      }
  return(result.ja.drug['ja_drug_sum'])
}

###################################################################### 
# Function to compute the feature "JAAF" discussed in the paper.
# G is a network object
# node1 denotes the ID of a drug
# node2 denotes the ID of an ADE
#
# This function returns a vector denotes the value of feature "JAAF" 
#####################################################################

compute.jaccard.aa.ADE.feature.sum = function(G, node1, node2){
  N1 = get.neighborhood(G, node1)
  N2 = get.neighborhood(G, node2)
  N1 = setdiff(N1, node2)             
  n.neighbors = length(N1)
  D = degree(G, gmode="graph")
  D.node2 = D[node2] 
  result.ja.ade = numeric(3)
  result.ja.ade[1:3] <- NA
  names(result.ja.ade) =
    c('ja_ade_min','ja_ade_max','ja_ade_sum')
  if (n.neighbors == 0) {
    result.ja.ade['ja_ade_min'] = 0
    result.ja.ade['ja_ade_max'] = 0
    result.ja.ade['ja_ade_mean'] = 0
  } else {
    ja.vector = numeric(n.neighbors)
    for (i in 1:n.neighbors) {
      neighbor.i = N1[i]    
      N2.i = get.neighborhood(G, neighbor.i)  
      intersection.i = intersect(N2, N2.i)    
      union.i = union(N2, N2.i)
      D.i=D[neighbor.i] 
      ja.vector[i] = (length(intersection.i)/length(union.i))/log(D.i)+1/log(D.i) 
    }
    result.ja.ade['ja_ade_min'] = min(ja.vector)
    result.ja.ade['ja_ade_max'] = max(ja.vector)
    result.ja.ade['ja_ade_sum'] = sum(ja.vector)
  }
  return(result.ja.ade['ja_ade_sum'])
}
