//  tlars_cpp.cpp

#include "tlars_cpp.h"
#include <limits>
#include <iterator>
#include <iostream>



// Constructors

/** Constructor for a new tlars_cpp-object
 *
 * Creates a new object of the class tlars_cpp.
 *
 * @param X Real valued Predictor matrix.
 * @param y Response vector.
 * @param verbose Logical. If TRUE progress in computations is shown.
 * @param intercept Logical. If TRUE an intercept is included.
 * @param standardize Logical. If TRUE the predictors are standardized and the response is centered.
 * @param num_dummies Number of dummies that are appended to the predictor matrix.
 * @param type Type of used algorithm (currently possible choices: 'lar' or 'lasso').
 */
tlars_cpp::tlars_cpp(arma::mat X, arma::vec y, bool verbose, bool intercept, bool standardize, int num_dummies, std::string type)
{
    this->X = X;
    this->y = y;
    this->verbose = verbose;
    this->intercept = intercept;
    this->standardize = standardize;
    this->num_dummies = num_dummies;
    this->type = type;
    initialize_values();
}

/** Constructor for a tlars_cpp-object with previous LARS-state as an input
 *
 * Re-creates an object of the class tlars_cpp based on a list of class variables that is obtained via get_all().
 *
 * This constructor is necessary for "warm re-starts" of the solution path at the last stopping point.
 * It is also necessary when at risk of losing the object pointer when, e.g., the
 * R sessions on multiple parallel workers are closed which leads to the invalidity of the pointers.
 * The object can then be re-created when the list with all class variables was extracted from the
 * previous object using get_all().
 *
 * @param lars_state Input list that was extracted from a previous tlars_cpp object using get_all().
 */
tlars_cpp::tlars_cpp(Rcpp::List lars_state)
{
    initialize_values(lars_state);
}



// Getters

/** Returns the estimate of the beta vector
 *
 * @return beta
 */
std::vector<double> tlars_cpp::get_beta()
{
    std::vector<double> last_beta = beta_state.back();
    for(int i = 0; i<p; i++)
    {
        last_beta[i] = last_beta.at(i)/norm_x(i);
    }
    return last_beta;
}

/** Returns a a matrix with the estimates of the beta vectors at all steps
 *
 * @return beta_path
 */
std::list<std::vector<double>> tlars_cpp::get_beta_path()
{
    std::list<std::vector<double>> beta;
    std::list<std::vector<double>>::iterator it;
    std::vector<double> curr_beta(p);
    for (it = beta_state.begin(); it != beta_state.end(); ++it)
    {
        curr_beta = *it;
        for(int i = 0; i<p; i++)
        {
            curr_beta[i] = curr_beta.at(i)/norm_x(i);

        }
        beta.push_back(curr_beta);
    }
    return beta;
}

/** Returns the number of active predictors
 *
 * @return num_active
 */
int tlars_cpp::get_num_active()
{
    return count_active_pred - count_dummies;
}

/** Returns the number of dummy predictors that have been included
 *
 * @return num_active_dummies
 */
int tlars_cpp::get_num_active_dummies()
{
    return count_dummies;
}

/** Returns the number of dummy predictors
 *
 * @return num_dummies
 */
 int tlars_cpp::get_num_dummies()
{
    return num_dummies;
}

/** Returns the indices of added/removed variables along the solution path
 *
 * @return actions
 */
std::list<int> tlars_cpp::get_actions()
{
    return actions;
}

/** Returns the degrees of freedom at each step which is
 * given by number of active variables (+1 if intercept is true)
 *
 * @return df
 */
std::list<int> tlars_cpp::get_df()
{
    update_df();
    return df;
}

/** Returns the R^2 statistic at each step
 *
 * @return R2
 */
std::list<double> tlars_cpp::get_R2()
{
    return R2;
}

/** Returns the residual sum of squares at each step
 *
 * @return RSS
 */
std::list<double> tlars_cpp::get_RSS()
{
    return RSS;
}

/** Returns the Cp-statistic at each step
 *
 * @return Cp
 */
arma::vec tlars_cpp::get_Cp()
{
    update_df();
    double rss_big = RSS.back();
    double df_big = n-df.back();
    double sigma2;
    if(rss_big>machine_prec && df_big>machine_prec)
    {
        sigma2 = rss_big/df_big;
    }
    else
    {
        sigma2 = NAN;
    }
    return double_list_to_vector(RSS) / sigma2 - n + 2*int_list_to_vector(df);
}

/** Returns the lambda-values (penalty parameters) at each
 * step along the solution path
 *
 * @return lambda
 */
arma::vec tlars_cpp::get_lambda()
{
    return lambda;
}

/** Returns the first entry/selection steps of the predictors
 * along the solution path
 *
 * @return first_in
 */
std::vector<int> tlars_cpp::get_entry()
{
    return first_in;
}

/** Returns the L2-norm of the predictors
 *
 * @return norm_x
 */
arma::vec tlars_cpp::get_norm_X()
{
    return norm_x;
}

/** Returns the sample means of the predictors
 *
 * @return mean_x
 */
arma::vec tlars_cpp::get_mean_X()
{
    return mean_x;
}

/** Returns the sample mean of the response y
 *
 * @return mean_y
 */
double tlars_cpp::get_mean_y()
{
    return mean_y;
}

/** Returns all class variables: This list can be used as an input to the constructor to re-create an object of class tlars_cpp
 *
 * @return lars_state
 */
Rcpp::List tlars_cpp::get_all()
{
    Rcpp::List l1 =  Rcpp::List::create(
                         Rcpp::Named("n") = n,
                         Rcpp::Named("p") = p,
                         Rcpp::Named("effective_n") = effective_n,
                         Rcpp::Named("active_pred") = active_pred,
                         Rcpp::Named("count_active_pred") = count_active_pred,
                         Rcpp::Named("new_pred") = new_pred,
                         Rcpp::Named("count_new_pred") = count_new_pred,
                         Rcpp::Named("inactive_pred") = inactive_pred,
                         Rcpp::Named("count_inactive_pred") = count_inactive_pred,
                         Rcpp::Named("ignored_pred") = ignored_pred,
                         Rcpp::Named("count_ignored_pred") = count_ignored_pred,
                         Rcpp::Named("norm_x") = norm_x,
                         Rcpp::Named("mean_x") = mean_x,
                         Rcpp::Named("mean_y") = mean_y,
                         Rcpp::Named("corr_predictors") = corr_predictors,
                         Rcpp::Named("pos_corr_predictors") = pos_corr_predictors,
                         Rcpp::Named("ssy") = ssy,
                         Rcpp::Named("residuals") = residuals,
                         Rcpp::Named("max_steps") = max_steps,
                         Rcpp::Named("beta_state") = beta_state
                     );

    Rcpp::List l2 =  Rcpp::List::create(
                         Rcpp::Named("RSS") = RSS,
                         Rcpp::Named("RSS_next") = RSS_next,
                         Rcpp::Named("R2") = R2,
                         Rcpp::Named("R2_next") = R2_next,
                         Rcpp::Named("lambda") = lambda,
                         Rcpp::Named("X") = X,
                         Rcpp::Named("y") = y,
                         Rcpp::Named("first_in") = first_in,
                         Rcpp::Named("active_data_decomp") = active_data_decomp,
                         Rcpp::Named("active_data_rank") = active_data_rank,
                         Rcpp::Named("A") = A,
                         Rcpp::Named("w") = w,
                         Rcpp::Named("Gi1") = Gi1,
                         Rcpp::Named("a") = a,
                         Rcpp::Named("u") = u,
                         Rcpp::Named("gamhat") = gamhat,
                         Rcpp::Named("max_gam1") = max_gam1,
                         Rcpp::Named("max_gam2") = max_gam2,
                         Rcpp::Named("gamrat") = gamrat,
                         Rcpp::Named("gamhat_list") = gamhat_list
                     );

    Rcpp::List l3 =  Rcpp::List::create(
                         Rcpp::Named("drop") = drop,
                         Rcpp::Named("drop_ind") = drop_ind,
                         Rcpp::Named("sign_vec") = sign_vec,
                         Rcpp::Named("verbose") = verbose,
                         Rcpp::Named("num_dummies") = num_dummies,
                         Rcpp::Named("standardize") = standardize,
                         Rcpp::Named("intercept") = intercept,
                         Rcpp::Named("type") = type,
                         Rcpp::Named("step_type") = step_type,
                         Rcpp::Named("count_dummies") = count_dummies,
                         Rcpp::Named("k") = k,
                         Rcpp::Named("early_stop") = early_stop,
                         Rcpp::Named("gamhat1") = gamhat1,
                         Rcpp::Named("gamhat2") = gamhat2,
                         Rcpp::Named("mod_X_matrix") = mod_X_matrix,
                         Rcpp::Named("next_beta") = next_beta,
                         Rcpp::Named("old_active_data_decomp") = old_active_data_decomp,
                         Rcpp::Named("active_beta") = active_beta,
                         Rcpp::Named("gam_lasso") = gam_lasso,
                         Rcpp::Named("machine_prec") = machine_prec
                     );

    Rcpp::List l4 =  Rcpp::List::create(
                         Rcpp::Named("actions") = actions,
                         Rcpp::Named("df") = df
                     );

    return Rcpp::List::create(
               Rcpp::Named("l1") = l1,
               Rcpp::Named("l2") = l2,
               Rcpp::Named("l3") = l3,
               Rcpp::Named("l4") = l4
           );
}



// Methods

/** Executes all necessary pre-processing steps
 *
 */
void tlars_cpp::initialize_values()
{

    // initialize dimensions p and sample size n
    n = X.n_rows;
    p = X.n_cols;

    // set machine precision
    machine_prec = std::numeric_limits<float>::denorm_min();


    //effective n by 1 reduced if intercept is true
    effective_n = n;
    if (intercept==true)
    {
        effective_n  = n-1;
    }

    // initialize dummy counter
    count_dummies = 0;

    // initialize the list that lists all predictors (all are inactive at the start)
    count_active_pred = 0;
    count_new_pred = 0;
    count_inactive_pred = p;
    for (i=0; i<p; i++)
    {
        inactive_pred.push_back(i);
    }

    // if intercept is true, remove the mean in the data X and in the output y
    mean_x = arma::zeros<arma::vec>(p);
    mean_y = 0;

    for (i=0; i<p; i++)
    {
        double dim_mean = 0;
        for (j=0; j<n; j++)
        {
            dim_mean = dim_mean + X(j,i);
        }
        mean_x(i) = dim_mean/n;
        if(intercept)
        {
            X.col(i) = X.col(i) - mean_x(i);
        }
    }
    for(j=0; j<n; j++)
    {
        mean_y = mean_y + y(j);
    }
    mean_y = mean_y/n;
    if(intercept)
    {
        y = y - mean_y;
    }


    // If standardize is true:
    // 1. If the variance of the signal is below the threshold epsilon, the predictor is ignored.
    // 2. The signal is standardized.
    Rcpp::LogicalVector ignored_pred(p);
    count_ignored_pred = 0;
    norm_x = arma::ones<arma::vec>(p);
    if (standardize == true)
    {
        for (it = inactive_pred.begin(); it != inactive_pred.end(); ++it)
        {
            double squared_sum = 0;
            for(j=0; j<n; j++)
            {
                squared_sum = squared_sum + pow(X(j,*it),2);
            }
            norm_x(*it) = sqrt(squared_sum);
            if (norm_x(*it)/sqrt(n)< machine_prec)
            {
                norm_x(*it) = machine_prec*sqrt(n);
                ignored_pred(*it) = true;
                count_ignored_pred++;
            }
            else
            {
                X.col(*it) = X.col(*it)/norm_x(*it);
            }
        }
        if (count_ignored_pred>0)
        {
            for(i=0; i<p; i++)
            {
                if(ignored_pred(i) == true)
                {
                    inactive_pred.remove(i);
                    count_inactive_pred--;
                }
            }
            if(verbose)
            {
                Rcpp::Rcout << count_ignored_pred << " predictor(s) dropped because of low variance \n";
            }
        }
    }

    // Initialize vector with correlations of the predictor data with y
    corr_predictors = (y.t() * X).t();
    Rcpp::LogicalVector pos_corr_predictors(p);

    // Initialize summed squared response and summed squared residuals
    ssy = dot(y, y);


    // Initialize residuals
    residuals = y;

    // Initialize maximum number of steps
    if (p<effective_n)
    {
        max_steps = 8*p;
    }
    else
    {
        max_steps = 8*effective_n;
    }

    //Initialize the first zero-entry of the beta-vector list
    std::vector<double> zero_vector(p,0);
    beta_state.push_back(zero_vector);

    //Initialize some statistical measures
    RSS.push_back(ssy);
    R2.push_back(0);

    //Initialize lambda vector
    lambda = arma::zeros<arma::vec>(max_steps);

    //Initialize vector that documents parameters entering the model
    std::vector<int> first_in(p);

    //Initialize R-matrix and the rank of the model
    active_data_decomp.set_size(1,1);
    active_data_rank = 0;

    //Initialize some loop-parameters
    k = 0;
    early_stop = false;
    drop = false;

    // Initialize remaining parameters
    A.set_size(1,1);

    // Initialize type for algorithm steps

    step_type = type;

    this->first_in = first_in;
    this->ignored_pred = ignored_pred;

}

/** Initializes values while recreating an object with a previous
 *  obtained list using get_all()
 *
 */
void tlars_cpp::initialize_values(Rcpp::List lars_state)
{

    // Extract inner lists from outer list
    Rcpp::List l1 = lars_state["l1"];
    Rcpp::List l2 = lars_state["l2"];
    Rcpp::List l3 = lars_state["l3"];
    Rcpp::List l4 = lars_state["l4"];

    // initialize all variables
    n = l1["n"];
    p = l1["p"];
    effective_n = l1["effective_n"];
    active_pred = Rcpp::as<std::list<int>>(l1["active_pred"]);
    count_active_pred = l1["count_active_pred"];
    new_pred = Rcpp::as<std::list<int>>(l1["new_pred"]);
    count_new_pred = l1["count_new_pred"];
    inactive_pred = Rcpp::as<std::list<int>>(l1["inactive_pred"]);
    count_inactive_pred = l1["count_inactive_pred"];
    ignored_pred = l1["ignored_pred"];
    count_ignored_pred = l1["count_ignored_pred"];
    norm_x = Rcpp::as<arma::vec>(l1["norm_x"]);
    mean_x = Rcpp::as<arma::vec>(l1["mean_x"]);
    mean_y = l1["mean_y"];
    corr_predictors = Rcpp::as<arma::vec>(l1["corr_predictors"]);
    pos_corr_predictors = l1["pos_corr_predictors"];
    ssy = l1["ssy"];
    residuals = Rcpp::as<arma::vec>(l1["residuals"]);
    max_steps = l1["max_steps"];
    beta_state = Rcpp::as<std::list<std::vector<double>>>(l1["beta_state"]);

    RSS = Rcpp::as<std::list<double>>(l2["RSS"]);
    RSS_next = l2["RSS_next"];
    R2 = Rcpp::as<std::list<double>>(l2["R2"]);
    R2_next = l2["R2_next"];
    lambda = Rcpp::as<arma::vec>(l2["lambda"]);
    X = Rcpp::as<arma::mat>(l2["X"]);
    y = Rcpp::as<arma::vec>(l2["y"]);
    first_in = Rcpp::as<std::vector<int>>(l2["first_in"]);
    active_data_decomp = Rcpp::as<arma::mat>(l2["active_data_decomp"]);
    active_data_rank = l2["active_data_rank"];
    A = Rcpp::as<arma::mat>(l2["A"]);
    w = Rcpp::as<arma::vec>(l2["w"]);
    Gi1 = Rcpp::as<arma::vec>(l2["Gi1"]);
    a = Rcpp::as<arma::vec>(l2["a"]);
    u = Rcpp::as<arma::vec>(l2["u"]);
    gamhat = l2["gamhat"];
    max_gam1 = l2["max_gam1"];
    max_gam2 = l2["max_gam2"];
    gamrat = Rcpp::as<std::list<double>>(l2["gamrat"]);
    gamhat_list = Rcpp::as<std::list<double>>(l2["gamhat_list"]);

    drop = l3["drop"];
    drop_ind = Rcpp::as<std::list<int>>(l3["drop_ind"]);
    sign_vec = Rcpp::as<arma::vec>(l3["sign_vec"]);
    verbose = l3["verbose"];
    num_dummies = l3["num_dummies"];
    standardize = l3["standardize"];
    intercept = l3["intercept"];
    type = Rcpp::as<std::string>(l3["type"]);
    step_type = Rcpp::as<std::string>(l3["step_type"]);
    count_dummies = l3["count_dummies"];
    k = l3["k"];
    early_stop = l3["early_stop"];
    gamhat1 = Rcpp::as<arma::vec>(l3["gamhat1"]);
    gamhat2 = Rcpp::as<arma::vec>(l3["gamhat2"]);
    mod_X_matrix = Rcpp::as<arma::mat>(l3["mod_X_matrix"]);
    next_beta = Rcpp::as<std::vector<double>>(l3["next_beta"]);
    old_active_data_decomp = Rcpp::as<arma::mat>(l3["old_active_data_decomp"]);
    active_beta = Rcpp::as<arma::vec>(l3["active_beta"]);
    gam_lasso = Rcpp::as<arma::vec>(l3["gam_lasso"]);
    machine_prec = l3["machine_prec"];

    actions = Rcpp::as<std::list<int>>(l4["actions"]);
    df = Rcpp::as<std::list<int>>(l4["df"]);

}

/** Executes T-LARS steps until a stopping-condition is satisfied
 *
 * @param T_stop Number of included dummies after which the random experiments (i.e., forward selection processes) are stopped.
 * @param early_stop Logical. If TRUE, then the forward selection process is stopped after T_stop dummies have been included. Otherwise
 * the entire solution path is computed.
 */
void tlars_cpp::execute_lars_step(int T_stop, bool early_stop)
{

    // Determine the index of the first dummy
    int dummy_ind = p - num_dummies;

    //Begin LARS-algorithm
    while (k < max_steps&&
            count_inactive_pred > 0 &&
            count_active_pred < effective_n &&
            (count_dummies < T_stop || early_stop == false))
    {
        //Obtain correlations of all inactive predictors
        arma::vec corr_inactive(count_inactive_pred);
        counter = 0;
        for (it = inactive_pred.begin(); it != inactive_pred.end(); ++it)
        {
            corr_inactive(counter) = corr_predictors(*it);
            counter++;
        }
        //Obtain maximum correlation over all inactive predictors
        double corr_max_inactive = max(abs(corr_inactive));



        //Break if maximum correlation is too low
        if (corr_max_inactive<100*machine_prec)
        {
            if(verbose == true)
            {
                Rcpp::Rcout << "Stopped because of too low correlations \n";
            }
            break;
        }
        //Set ideal lambda for this step equal to the maximum correlation
        lambda(k) = corr_max_inactive;

        //
        if(drop==false)
        {

            // Define predictors that will enter the set of active predictors
            new_pred.clear();
            for (it = inactive_pred.begin(); it != inactive_pred.end(); ++it)
            {
                if(corr_predictors(*it)>= corr_max_inactive - machine_prec ||
                        corr_predictors(*it)<= -corr_max_inactive + machine_prec)
                {
                    new_pred.push_back(*it);
                }
            }

            // For every new predictor do:
            for (it = new_pred.begin(); it!= new_pred.end(); it++)
            {
                arma::mat oldX(n,count_active_pred);
                counter= 0;
                // Create oldX which is the predictor matrix X of only the active predictors
                for (inner_it = active_pred.begin(); inner_it!= active_pred.end(); inner_it++)
                {
                    oldX.col(counter) = X.col(*inner_it);
                    counter++;
                }
                // Check for rank including a new predictor
                old_active_data_decomp = active_data_decomp;
                update_decomp(X.col(*it), oldX);
                // If the new predictor is linear dependent on the previous ones, ignore new predictor.
                if(active_data_rank == count_active_pred)
                {
                    active_data_decomp = old_active_data_decomp;
                    ignored_pred(*it) = true;
                    count_ignored_pred++;

                    if(verbose)
                    {
                        Rcpp::Rcout << "Predictor dropped because of linear dependency \n";
                    }

                    // Else add the new predictor to the active set
                }
                else
                {
                    // Did the new predictor enter the set for the first time?
                    if (first_in[*it]==0)
                    {
                        first_in[*it] = k;
                    }
                    active_pred.push_back(*it);
                    actions.push_back(*it+1);

                    // Add 1 to the dummy counter if the corresponding predictor was a dummy.
                    if (*it>=dummy_ind)
                    {
                        count_dummies++;
                    }
                    count_active_pred++;
                }
                counter = 0;
                for (inner_it = inactive_pred.begin(); inner_it!= inactive_pred.end(); inner_it++)
                {
                    if (*inner_it== *it)
                    {
                        corr_inactive.shed_rows(counter,counter);
                    }
                    counter++;
                }
                inactive_pred.remove(*it);
                count_inactive_pred--;
            }
        }
        // Calculate sign-vector
        counter = 0;
        sign_vec.resize(count_active_pred);
        for (it = active_pred.begin(); it!= active_pred.end(); it++)
        {
            if (corr_predictors(*it) >= 0)
                sign_vec(counter) = 1;
            else
                sign_vec(counter) = -1;
            counter++;
        }
        // Calculate Lars Step
        Gi1 = solve_upper_triangular(active_data_decomp,solve_lower_triangular(active_data_decomp.t(),sign_vec));

        A = Gi1.t() * sign_vec;
        A = sqrt(1/A);
        w = (A*Gi1.t()).t();
        mod_X_matrix.resize(n,count_active_pred);
        counter=0;
        for (it = active_pred.begin(); it!= active_pred.end(); it++)
        {
            mod_X_matrix.col(counter) = X.col(*it);
            counter++;
        }
        u = mod_X_matrix*w;
        if(count_active_pred >= effective_n || count_active_pred >= p - count_ignored_pred)
            gamhat = corr_max_inactive/A(0,0);
        else
        {
            mod_X_matrix.resize(n,count_inactive_pred);
            counter = 0;
            for (it = inactive_pred.begin(); it!= inactive_pred.end(); it++)
            {
                mod_X_matrix.col(counter) = X.col(*it);
                counter++;
            }
            a = (u.t()*mod_X_matrix).t();
            gamhat1 = (corr_max_inactive - corr_inactive)/(A(0,0) - a);
            gamhat2 = (corr_max_inactive + corr_inactive)/(A(0,0) + a);
            max_gam1 = gamhat1.max();
            max_gam2 = gamhat2.max();
            if (max_gam1 >= max_gam2 && max_gam1>=0)
            {
                gamhat = max_gam1;
            }
            else if(max_gam2 >=0)
            {
                gamhat= max_gam2;
            }
            else
            {
                if(verbose)
                {
                    Rcpp::Rcout << "Warning: No positve Gamma. Exiting tlars \n";
                }
                break;
            }
            for(i=0; i<count_inactive_pred; i++)
            {
                if(gamhat1(i)<gamhat && gamhat1(i)>=machine_prec)
                {
                    gamhat=gamhat1(i);
                }
                if(gamhat2(i)<gamhat && gamhat2(i)>=machine_prec)
                {
                    gamhat=gamhat2(i);
                }
            }
        }
        next_beta = beta_state.back();
        // check if variables need to be removed
        if(step_type=="lasso")
        {
            drop = false;
            active_beta.set_size(active_pred.size());
            counter = 0;
            for (it = active_pred.begin(); it!= active_pred.end(); it++)
            {
                active_beta(counter) = next_beta[*it];
                counter ++;
            }
            gam_lasso = -active_beta/w;
            for(i=0; i<count_active_pred; i++)
            {
                if(gam_lasso(i)<gamhat && gam_lasso(i)>=machine_prec)
                {
                    drop = true;
                    gamhat = gam_lasso(i);
                }
            }
        }
        residuals = residuals - gamhat*u;
        corr_predictors = (residuals.t() * X).t();
        gamrat.push_back(gamhat*A(0,0)/corr_max_inactive);
        gamhat_list.push_back(gamhat);
        counter=0;
        for(it = active_pred.begin(); it!= active_pred.end(); it++)
        {
            next_beta[*it] = next_beta[*it] + gamhat*w(counter);
            counter++;
        }
        if (drop == true)
        {
            counter = 0;
            it = active_pred.begin();
            while (it!= active_pred.end())
            {
                if(gamhat == gam_lasso(counter))
                {
                    remove_var_from_decomp(counter);
                    count_active_pred--;
                    inactive_pred.push_back(*it);
                    actions.push_back(-*it-1);
                    count_inactive_pred++;
                    next_beta[*it] = 0;
                    counter++;
                    active_pred.erase(it++);
                    //++it;
                    if (*it>=dummy_ind)
                    {
                        count_dummies--;
                    }
                }
                else
                {
                    counter ++;
                    ++it;
                }

            }
        }
        beta_state.push_back(next_beta);
        k++;


        // Calculate some outputs
        RSS_next = 0;
        for(j=0; j<n; j++)
        {
            RSS_next = RSS_next + pow(residuals(j),2);
        }
        RSS.push_back(RSS_next);
        R2_next = 1 - RSS_next/ssy;
        R2.push_back(R2_next);
    }
}

/** Add a predictor to the Cholesky-decomposition of the active data
 *
 * @param new_X Predictor to be added.
 * @param old_X Matrix containing all active predictors.
 */
void tlars_cpp::update_decomp(arma::mat new_X, arma::mat old_X)
{
    double xtx=0;
    double norm_xnew;
    int dim = active_data_decomp.n_cols;
    int j;
    for(j=0; j<n; j++)
    {
        xtx = xtx + pow(new_X(j),2);
    }
    norm_xnew= sqrt(xtx);
    if(active_data_rank == 0)
    {
        active_data_rank = 1;
        active_data_decomp(0,0) = norm_xnew;
    }
    else
    {
        arma::vec Xtx = (new_X.t() * old_X).t();
        arma::vec r;
        r = solve_lower_triangular(active_data_decomp.t(), Xtx);
        double rpp = pow(norm_xnew, 2);
        for(j=0; j<active_data_rank; j++)
        {
            rpp = rpp - pow(r(j),2);
        }
        // Check for machine singularity
        if(rpp<=machine_prec) rpp = machine_prec;
        else
        {
            rpp = sqrt(rpp);
            active_data_rank++;
            active_data_decomp.resize(dim+1,dim+1);
            for (j=0; j<active_data_rank-1; j++)
            {
                active_data_decomp(j,dim) = r(j);
            }
            active_data_decomp(dim,dim) = rpp;
        }
    }

}

/** Remove a variable from the Cholesky decomposition
 *
 * @param removal_index Index of the variable to be removed.
 */
void tlars_cpp::remove_var_from_decomp(int removal_index)
{
    if(removal_index == active_data_rank-1)
    {
        active_data_decomp.shed_cols(removal_index, removal_index);
        active_data_decomp.shed_rows(count_active_pred-1,count_active_pred-1);
        active_data_rank--;
    }
    else
    {
        arma::mat modified_submatrix;
        arma::mat modified_submatrix2;
        modified_submatrix = active_data_decomp.submat(removal_index+1, removal_index+1, count_active_pred-1, count_active_pred-1);
        modified_submatrix = modified_submatrix.t()*modified_submatrix;
        modified_submatrix2= active_data_decomp.submat(removal_index, removal_index+1, removal_index, count_active_pred-1);
        modified_submatrix2 = modified_submatrix2.t()*modified_submatrix2;
        modified_submatrix = cholesky_decomp(modified_submatrix + modified_submatrix2);
        active_data_decomp.shed_cols(removal_index, removal_index);
        active_data_decomp.shed_rows(count_active_pred-1,count_active_pred-1);
        active_data_decomp.submat(removal_index, removal_index, count_active_pred-2, count_active_pred-2) = modified_submatrix;
        active_data_rank--;
    }
}

/** Solves A*y = x for y where A is an upper triangular matrix
 *
 * @param upper_t Upper triangular matrix A.
 * @param x Vector on the right side of the equation.
 *
 * @return Solution for vector y.
 */
arma::vec tlars_cpp::solve_upper_triangular(arma::mat upper_t, arma::vec x)
{
    int n_solve = x.n_elem;
    int max_index_solve = n_solve-1;
    for (i=max_index_solve; i>=0; i--)
    {
        for(j=max_index_solve; j>i; j--)
        {
            x(i) = x(i) - upper_t(i,j) * x(j);
        }
        x(i) = x(i)/upper_t(i,i);
    }
    return x;
}

/** Solves A*y = x for y where A is a lower triangular matrix
 *
 * @param lower_t Lower triangular matrix A.
 * @param x Vector on the right side of the equation.
 *
 * @return Solution for vector y.
 */
arma::vec tlars_cpp::solve_lower_triangular(arma::mat lower_t, arma::vec x)
{
    int n_solve = x.n_elem;
    for (i=0; i<n_solve; i++)
    {
        for(j=0; j<i; j++)
        {
            x(i) = x(i) - lower_t(i,j) * x(j);
        }
        x(i) = x(i)/lower_t(i,i);
    }
    return x;
}

/** Computes the Cholesky-decomposition of a square matrix
 *
 * @param square_matrix A square matrix.
 *
 * @return Cholesky-decomposition of square_matrix.
 */
arma::mat tlars_cpp::cholesky_decomp(arma::mat square_matrix)
{
    int mat_size = square_matrix.n_rows;
    arma::mat lower_triangular(mat_size,mat_size);
    double sum;
    int h;
    for(i=0; i<mat_size; i++)
    {

        for (j=0; j<=i; j++)
        {
            sum = 0;
            for (h=0; h<j; h++)
            {
                sum = sum + lower_triangular(i,h)*lower_triangular(j,h);
            }
            if(j!=i)
            {
                lower_triangular(i,j) = (square_matrix(i,j)-sum)/lower_triangular(j,j);
            }
            else
            {
                lower_triangular(j,j) = std::sqrt(square_matrix(j,j)-sum);
            }
        }
    }
    return lower_triangular.t();
}

/** Updates the degrees of freedom
 *
 */
void tlars_cpp::update_df()
{
    df.clear();
    counter=0;
    if(intercept)
    {
        counter++;
    }
    df.push_back(counter);
    for (it = actions.begin(); it!= actions.end(); it++)
    {
        if (*it > 0)
        {
            counter++;
        }
        else
        {
            counter--;
        }
        df.push_back(counter);
    }
}

/** Converts a std::list<double> into an arma::vec
 *
 * @param double_list A std::list<double>.
 *
 * @return Transformed list of type arma::vec.
 */
arma::vec tlars_cpp::double_list_to_vector(std::list<double> double_list)
{
    arma::vec output(double_list.size());
    counter= 0;
    for (double_it = double_list.begin(); double_it!= double_list.end(); double_it++)
    {
        output(counter) = *double_it;
        counter++;
    }
    return output;
}

/** Converts a std::list<int> into an arma::vec
 *
 * @param int_list A std::list<int>.
 *
 * @return Transformed list of type arma::vec.
 */
arma::vec tlars_cpp::int_list_to_vector(std::list<int> int_list)
{
    arma::vec output(int_list.size());
    counter= 0;
    for (it = int_list.begin(); it!= int_list.end(); it++)
    {
        output(counter) = *it;
        counter++;
    }
    return output;
}
