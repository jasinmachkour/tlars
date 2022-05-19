// export_tlars_cpp.cpp

#include "tlars_cpp.h"

RCPP_MODULE(tlars_cpp){
    Rcpp::class_<tlars_cpp>("tlars_cpp")

    // Constructors
    .constructor<arma::mat, arma::vec , bool , bool , bool , int, std::string>("Creates a new object of the class tlars_cpp.")
    .constructor<Rcpp::List>("Re-creates an object of the class tlars_cpp based on a list of class variables that is obtained via get_all().")

    // Methods
    .method("execute_lars_step", &tlars_cpp::execute_lars_step, "Void. Executes lars steps until a stopping-condition is satisfied.")

    // Output Getters
    .method("get_beta", &tlars_cpp::get_beta,"Returns the estimate of the beta vector.")
    .method("get_beta_path", &tlars_cpp::get_beta_path,"Returns a a matrix with the estimates of the beta vectors at all steps.")
    .method("get_num_active",&tlars_cpp::get_num_active, "Returns the number of active predictors.")
    .method("get_num_active_dummies",&tlars_cpp::get_num_active_dummies, "Returns the number of dummy variables that have been included.")
    .method("get_num_dummies",&tlars_cpp::get_num_dummies, "Returns the number of dummy predictors.")
    .method("get_actions", &tlars_cpp::get_actions, "Returns the indices of added/removed variables along the solution path.")
    .method("get_df", &tlars_cpp::get_df, "Returns the degrees of freedom at each step which is given by number of active variables (+1 if intercept is true).")
    .method("get_R2", &tlars_cpp::get_R2, "Returns the R^2 statistic at each step.")
    .method("get_RSS", &tlars_cpp::get_RSS, "Returns the residual sum of squares at each step.")
    .method("get_Cp", &tlars_cpp::get_Cp, "Returns the Cp-statistic at each step.")
    .method("get_lambda",&tlars_cpp::get_lambda, "Returns the lambda-values (penalty parameters) at each step along the solution path.")
    .method("get_entry", &tlars_cpp::get_entry, "Returns the first entry/selection steps of the predictors along the solution path.")
    .method("get_norm_X", &tlars_cpp::get_norm_X, "Returns the L2-norm of the predictors.")
    .method("get_meanX", &tlars_cpp::get_mean_X, "Returns the sample means of the predictors.")
    .method("get_mean_y", &tlars_cpp::get_mean_y, "Returns the sample mean of the response y.")
    .method("get_all", &tlars_cpp::get_all, "Returns all class variables: This list can be used as an input to the constructor to re-create an object of class tlars_cpp.")

    // Input variables
    .field("X", &tlars_cpp::X, "Real valued predictor matrix.")
    .field("y", &tlars_cpp::y, "Response vector.")
    .field("verbose", &tlars_cpp::verbose, "Logical. If TRUE progress in computations is shown.")
    .field("intercept", &tlars_cpp::intercept, "Logical. If TRUE an intercept is included.")
    .field("standardize", &tlars_cpp::standardize, "Logical. If TRUE the predictors are standardized and the response is centered.")
    .field("num_dummies", &tlars_cpp::num_dummies, "Number of dummies that are appended to the predictor matrix.")
    .field("type", &tlars_cpp::type, "Type of used algorithm (currently possible choices: 'lar' or 'lasso').")
    ;
}
