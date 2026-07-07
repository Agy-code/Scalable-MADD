The following libraries are needed here: 
Rfast: For distance matrix calculation
Caret: Used in creating fold IDs for cross-validation

Here, the main function is NN_MADD: a Nearest Neighbor (NN) classifier that uses MADD distances to classify observations.
It may use either the entire training set or a subset of it as a reference set. This function uses the
functions given in: NN with MADD and MADD_sc, Sampling from a DPP, Class-wise DPP index selection, and Cross-validation as 
helper function.
  
