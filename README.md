The following libraries are needed: 

Rfast: For distance matrix calculation

Caret: Used in creating fold IDs for cross-validation

The main function here is the NN_MADD classifier. It utilizes the other functions as helper functions. The inputs are:

trainX: The features of the training data

trainY: The class labels of the training data

testX: The features of the test data

testY: The class labels of the test data (if available)

use_reference: This indicates whether the MADD-based Nearest neighbor classifier uses the full training data (NN-MADD) or  its scalable version (NN-MADDsc). TRUE implies it will use references and implement the scalable version of MADD. 

reference.ind: A vector of indices considered as a reference set provided by the user.

reference_k: No reference.ind is provided, but the user specifies the number of references to be used per class.

do_cv: Whether cross-validation (CV) is to be used or not. TRUE implies CV is to be done.

reference_k_choices: If do_cv= TRUE, then the range of reference sizes based on which CV is performed.

Dtrain: Distance matrix of the training data points.

memory.eff: Indicates whether the user wants to use the memory-efficient approach, which avoids storing the entire distance matrix for the training data points. This is comparatively slower. The default is "NO".

n_tol_dpp: If the sample size of any class crosses n_tol_dpp, it uses the RFF version of k-DPP. The default choice is 2000.

use_RFF: If the user wants to forcefully use the RFF approximation.

D_RFF: Value of D used for RFF, default is 500.

dpp_method: "kDPP"

cv_sample_method: "Deterministic"

return_MADD: If the user wants to print the MADD matrix, then return_MADD = TRUE. The default is FALSE.

                                
