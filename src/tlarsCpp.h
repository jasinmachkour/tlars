//  tlarsCpp.h

#ifndef tlarsCpp_H
#define tlarsCpp_H

#include <vector>
#include <RcppArmadillo.h>

using namespace Rcpp;

class tlarsCpp
{
public:


    // Constructors
    tlarsCpp(arma::mat X, arma::vec y , bool verbose, bool intercept, bool normalize, int L_val, std::string type);
    tlarsCpp(Rcpp::List lars_state);

    // Methods
    void executeLarsStep(int T_stop, bool earlyStop);
    void initializeValues(arma::mat X, arma::vec y , bool verbose, bool intercept, bool normalize, std::string type);
    void initializeValues(Rcpp::List lars_state);
    void update_decomp(arma::mat newX, arma::mat oldX);
    void remove_var_from_decomp(int removal_index);
    arma::vec solveUpperTriangular(arma::mat upperT_X, arma::vec vec_b);
    arma::vec solveLowerTriangular(arma::mat lowerT_X, arma::vec vec_b);
    arma::mat choleskyDecomp(arma::mat square_matrix);
    arma::vec doubleListToVector(std::list<double> doubleList);
    arma::vec intListToVector(std::list<int> intList);
    void updateDf();


    // Output Getters
    Rcpp::List getAll();
    arma::mat getActive_data_decomp();
    int getActive_data_rank();
    std::string getType();
    arma::vec getLambda();
    arma::vec getR2();
    arma::vec getRSS();
    arma::vec getCp();
    std::list<int> getActions();
    std::list<int> getDf();
    std::vector<int> getEntry();
    arma::vec getGamrat();
    arma::vec getGamHat();
    std::list<std::vector<double>> getBeta();
    std::vector<double> getLastBeta();
    double getMu();
    arma::vec getNormX();
    arma::vec getMeanX();

// private:

    //'state variables':
        int n;
        int p;
        int effective_n;
        std::list<int> active_pred;
        int count_active_pred;
        std::list<int> new_pred;
        int count_new_pred;
        std::list<int> inactive_pred;
        int count_inactive_pred;
        LogicalVector ignored_pred;
        int count_ignored_pred;
        arma::vec norm_x;
        arma::vec mean_x;
        double mean_y;
        arma::vec corr_predictors;
        LogicalVector pos_corr_predictors;
        double ssy;
        arma::vec residuals;
        int max_steps;
        std::list<std::vector<double>> beta_state;
        std::list<double> RSS;
        double RSS_next;
        std::list<double> R2;
        double R2_next;
        arma::vec lambda;
        arma::mat X;
        arma::vec y;
        std::vector<int> integerInputs;
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
        bool verbose;
        int L_val;
        bool normalize;
        bool intercept;
        std::string type;
        std::string step_type;
        int i;
        int j;
        int counter;
        int count_knockoffs;
        std::list<int>::iterator it;
        std::list<int>::iterator inner_it;
        std::list<double>::iterator double_it;
        int k;
        bool e_stop;
        arma::vec gamhat1;
        arma::vec gamhat2;
        arma::mat modX_matrix;
        std::vector<double> next_beta;
        arma::mat old_active_data_decomp;
        arma::vec active_beta;
        arma::vec gam_lasso;
        double machinePrec;
        std::list<int> actions;
        std::list<int> df;


};

#endif /* tlarsCpp_h */
