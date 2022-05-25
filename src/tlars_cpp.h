// tlars_cpp.h

#ifndef tlars_cpp_H
#define tlars_cpp_H

#include <vector>
#include <RcppArmadillo.h>

/**
 * Class that implements the Terminating-LARS (T-LARS) algorithm.
 *
 */

class tlars_cpp
{
public:

    // Constructors
    tlars_cpp(arma::mat X, arma::vec y, bool verbose, bool intercept, bool standardize, int num_dummies, std::string type);
    tlars_cpp(Rcpp::List lars_state);

    // Methods
    void execute_lars_step(int T_stop, bool early_stop);

    // Output Getters
    std::vector<double> get_beta();
    std::list<std::vector<double>> get_beta_path();
    int get_num_active();
    int get_num_active_dummies();
    int get_num_dummies();
    std::list<int> get_actions();
    std::list<int> get_df();
    std::list<double> get_R2();
    std::list<double> get_RSS();
    arma::vec get_Cp();
    arma::vec get_lambda();
    std::vector<int> get_entry();
    double get_mean_y();
    arma::vec get_norm_X();
    arma::vec get_mean_X();
    Rcpp::List get_all();

    // State variables
    arma::mat X;
    arma::vec y;
    bool verbose;
    bool intercept;
    bool standardize;
    int num_dummies;
    std::string type;



private:

    // Methods
    void initialize_values();
    void initialize_values(Rcpp::List lars_state);
    void update_decomp(arma::mat new_X, arma::mat old_X);
    void remove_var_from_decomp(int removal_index);
    arma::vec solve_upper_triangular(arma::mat upperT_X, arma::vec vec_b);
    arma::vec solve_lower_triangular(arma::mat lowerT_X, arma::vec vec_b);
    arma::mat cholesky_decomp(arma::mat square_matrix);
    arma::vec double_list_to_vector(std::list<double> double_list);
    arma::vec int_list_to_vector(std::list<int> int_list);
    void update_df();

    // State variables
    int n;
    int p;
    int effective_n;
    std::list<int> active_pred;
    int count_active_pred;
    std::list<int> new_pred;
    int count_new_pred;
    std::list<int> inactive_pred;
    int count_inactive_pred;
    Rcpp::LogicalVector ignored_pred;
    int count_ignored_pred;
    arma::vec mean_x;
    arma::vec norm_x;
    double mean_y;
    arma::vec corr_predictors;
    Rcpp::LogicalVector pos_corr_predictors;
    double ssy;
    arma::vec residuals;
    int max_steps;
    std::list<std::vector<double>> beta_state;
    std::list<double> RSS;
    double RSS_next;
    std::list<double> R2;
    double R2_next;
    arma::vec lambda;
    std::vector<int> first_in;
    arma::mat active_data_decomp;
    int active_data_rank;
    arma::mat A;
    arma::vec w;
    arma::vec Gi1;
    arma::vec a;
    arma::vec u;
    double gamhat;
    double max_gam1;
    double max_gam2;
    std::list<double> gamrat;
    std::list<double> gamhat_list;
    bool drop;
    std::list<int> drop_ind;
    arma::vec sign_vec;
    std::string step_type;
    int i;
    int j;
    int counter;
    int count_dummies;
    std::list<int>::iterator it;
    std::list<int>::iterator inner_it;
    std::list<double>::iterator double_it;
    int k;
    bool early_stop;
    arma::vec gamhat1;
    arma::vec gamhat2;
    arma::mat mod_X_matrix;
    std::vector<double> next_beta;
    arma::mat old_active_data_decomp;
    arma::vec active_beta;
    arma::vec gam_lasso;
    double machine_prec;
    std::list<int> actions;
    std::list<int> df;
};

#endif /* tlars_cpp_h */
