The following libraries are needed here: 

Rfast: For distance matrix calculation

Caret: Used in creating fold IDs for cross-validation

The main function here is the NN_MADD classifier. It utilizes the other functions as helper functions. The inputs are:

trainX: The features of the training data
trainY: The class labels of the training data
testX: The features of the test data
testY: The class labels of the test data (if available)
use_reference: This is used to indicate if MADD based Nearest neighbor classifier uses the full training data (NN-MADD) or
               its scalable version (NN-MADDsc)
reference.ind: A vector of indices considered as a reference set,
reference_k: 
reference_k_choices: NULL,
do_cv =
Dtrain = NULL,
memory.eff = "NO",
n_tol_dpp = 2000,
use_RFF = FALSE,
D_RFF = 500,
dpp_method = "kDPP"
cv_sample_method = "Deterministic"
return_MADD = FALSE
                                
