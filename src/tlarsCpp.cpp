//  tlarsCpp.cpp

#include <limits>
#include "tlarsCpp.h"
#include <vector>
#include <iterator>
#include <RcppArmadillo.h>
#include <iostream>

using namespace Rcpp;


// Constructors
tlarsCpp::tlarsCpp(arma::mat X, arma::vec y, bool verbose, bool intercept, bool normalize, int L_val, std::string type)
{
    this->X = X;
    this->y = y;
    this->verbose = verbose;
    this->intercept = intercept;
    this->normalize = normalize;
    this->L_val = L_val;
    this->type = type;
    this->initializeValues(X, y, verbose, intercept, normalize, type);
}

tlarsCpp::tlarsCpp(Rcpp::List lars_state)
{
    this->initializeValues(lars_state);
}

// Getters
Rcpp::List tlarsCpp::getAll()
{
    Rcpp::List l1 =  Rcpp::List::create(
                         Rcpp::Named("n") = n,
                         _("p") = p,
                         _("effective_n") = effective_n,
                         _("active_pred") = active_pred,
                         _("count_active_pred") = count_active_pred,
                         _("new_pred") = new_pred,
                         _("count_new_pred") = count_new_pred,
                         _("inactive_pred") = inactive_pred,
                         _("count_inactive_pred") = count_inactive_pred,
                         _("ignored_pred") = ignored_pred,
                         _("count_ignored_pred") = count_ignored_pred,
                         _("norm_x") = norm_x,
                         _("mean_x") = mean_x,
                         _("mean_y") = mean_y,
                         _("corr_predictors") = corr_predictors,
                         _("pos_corr_predictors") = pos_corr_predictors,
                         _("ssy") = ssy,
                         _("residuals") = residuals,
                         _("max_steps") = max_steps,
                         _("beta_state") = beta_state
                     );

    Rcpp::List l2 =  Rcpp::List::create(
                         Rcpp::Named("RSS") = RSS,
                         _("RSS_next") = RSS_next,
                         _("R2") = R2,
                         _("R2_next") = R2_next,
                         _("lambda") = lambda,
                         _("X") = X,
                         _("y") = y,
                         _("integerInputs") = integerInputs,
                         _("first_in") = first_in,
                         _("active_data_decomp") = active_data_decomp,
                         _("active_data_rank") = active_data_rank,
                         _("A") = A,
                         _("w") = w,
                         _("Gi1") = Gi1,
                         _("a") = a,
                         _("u") = u,
                         _("gamhat") = gamhat,
                         _("max_gam1") = max_gam1,
                         _("max_gam2") = max_gam2,
                         _("gamrat") = gamrat
                     );

    Rcpp::List l3 =  Rcpp::List::create(
                         Rcpp::Named("gamhat_list") = gamhat_list,
                         _("drop") = drop,
                         _("drop_ind") = drop_ind,
                         _("sign_vec") = sign_vec,
                         _("verbose") = verbose,
                         _("L_val") = L_val,
                         _("normalize") = normalize,
                         _("intercept") = intercept,
                         _("type") = type,
                         _("step_type") = step_type,
                         _("count_knockoffs") = count_knockoffs,
                         _("k") = k,
                         _("e_stop") = e_stop,
                         _("gamhat1") = gamhat1,
                         _("gamhat2") = gamhat2,
                         _("modX_matrix") = modX_matrix,
                         _("next_beta") = next_beta,
                         _("old_active_data_decomp") = old_active_data_decomp,
                         _("active_beta") = active_beta,
                         _("gam_lasso") = gam_lasso
                     );

    Rcpp::List l4 =  Rcpp::List::create(
                         Rcpp::Named("machinePrec") = machinePrec,
                         _("actions") = actions,
                         _("df") = df
                     );

    return Rcpp::List::create(
               Rcpp::Named("l1") = l1,
               _("l2") = l2,
               _("l3") = l3,
               _("l4") = l4
           );
}

arma::mat tlarsCpp::getActive_data_decomp()
{
    return active_data_decomp;
}

int tlarsCpp::getActive_data_rank()
{
    return active_data_rank;
}

std::string tlarsCpp::getType()
{
    return type;
}

arma::vec tlarsCpp::getLambda()
{
    return lambda;
}

arma::vec tlarsCpp::getR2()
{
    return doubleListToVector(R2);
}

arma::vec tlarsCpp::getRSS()
{
    return doubleListToVector(RSS);
}

arma::vec tlarsCpp::getCp()
{
    updateDf();
    double rss_big = RSS.back();
    double df_big = n-df.back();
    double sigma2;
    if(rss_big>machinePrec && df_big>machinePrec)
    {
        sigma2 = rss_big/df_big;
    }
    else
    {
        sigma2 = NAN;
    }
    return doubleListToVector(RSS) / sigma2 - n + 2*intListToVector(df);
}


std::list<int> tlarsCpp::getActions()
{
    return actions;
}

std::list<int> tlarsCpp::getDf()
{
    updateDf();
    return df;
}

std::vector<int> tlarsCpp::getEntry()
{
    return first_in;
}

arma::vec tlarsCpp::getGamrat()
{
    return doubleListToVector(gamrat);
}

arma::vec tlarsCpp::getGamHat()
{
    return doubleListToVector(gamhat_list);
}

std::list<std::vector<double>> tlarsCpp::getBeta()
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

std::vector<double> tlarsCpp::getLastBeta()
{
    std::vector<double> last_beta = beta_state.back();
    for(int i = 0; i<p; i++)
    {
        last_beta[i] = last_beta.at(i)/norm_x(i);
    }
    return last_beta;
}

double tlarsCpp::getMu()
{
    return mean_y;
}

arma::vec tlarsCpp::getNormX()
{
    return norm_x;
}

arma::vec tlarsCpp::getMeanX()
{
    return mean_x;
}


// Methods
void tlarsCpp::initializeValues(arma::mat X, arma::vec y, bool verbose, bool intercept, bool normalize, std::string type)
{

    // initialize dimensions p and sample size n
    n = X.n_rows;
    p = X.n_cols;

    // set machine precision
    machinePrec = std::numeric_limits<float>::denorm_min();


    //effective n by 1 reduced if intercept is true
    effective_n = n;
    if (intercept==true)
    {
        effective_n  = n-1;
    }

    // initialize knockoff counter
    count_knockoffs = 0;

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
    if (intercept==true)
    {
        for (i=0; i<p; i++)
        {
            double dim_mean = 0;
            for (j=0; j<n; j++)
            {
                dim_mean = dim_mean + X(j,i);
            }
            mean_x(i) = dim_mean/n;
            X.col(i) = X.col(i) - mean_x(i);
        }
        for(j=0; j<n; j++)
        {
            mean_y = mean_y + y(j);
        }
        mean_y = mean_y/n;
        y = y - mean_y;
    }

    // If normalize is true:
    // 1. If the variance of the signal is below the threshold epsilon, the predictor is ignored.
    // 2. The signal is normalized.
    LogicalVector ignored_pred(p);
    count_ignored_pred = 0;
    norm_x = arma::ones<arma::vec>(p);
    if (normalize == true)
    {
        for (it = inactive_pred.begin(); it != inactive_pred.end(); ++it)
        {
            double squared_sum = 0;
            for(j=0; j<n; j++)
            {
                squared_sum = squared_sum + pow(X(j,*it),2);
            }
            norm_x(*it) = sqrt(squared_sum);
            if (norm_x(*it)/sqrt(n)< machinePrec)
            {
                norm_x(*it) = machinePrec*sqrt(n);
                ignored_pred.push_back(*it);
                count_ignored_pred++;
                inactive_pred.remove(*it);
                count_inactive_pred--;
            }
            else
            {
                X.col(*it) = X.col(*it)/norm_x(*it);
            }
        }
        if (count_ignored_pred>0 && verbose==true)
        {
            // std::cout << count_ignored_pred << " predictors with too low power dropped while initializing";
        }
    }

    // Initialize vector with correlations of the predictor data with y
    corr_predictors = (y.t() * X).t();
    LogicalVector pos_corr_predictors(p);

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
    e_stop = false;
    drop = false;

    // Initialize remaining parameters
    A.set_size(1,1);

    // Initialize type for algorithm steps

    this->step_type = type;

    this->X = X;
    this->first_in = first_in;
    this->ignored_pred = ignored_pred;

}

void tlarsCpp::initializeValues(Rcpp::List lars_state)
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
    integerInputs = Rcpp::as<std::vector<int>>(l2["integerInputs"]);
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

    gamhat_list = Rcpp::as<std::list<double>>(l3["gamhat_list"]);
    drop = l3["drop"];
    drop_ind = Rcpp::as<std::list<int>>(l3["drop_ind"]);
    sign_vec = Rcpp::as<arma::vec>(l3["sign_vec"]);
    verbose = l3["verbose"];
    L_val = l3["L_val"];
    normalize = l3["normalize"];
    intercept = l3["intercept"];
    type = Rcpp::as<std::string>(l3["type"]);
    step_type = Rcpp::as<std::string>(l3["step_type"]);
    count_knockoffs = l3["count_knockoffs"];
    k = l3["k"];
    e_stop = l3["e_stop"];
    gamhat1 = Rcpp::as<arma::vec>(l3["gamhat1"]);
    gamhat2 = Rcpp::as<arma::vec>(l3["gamhat2"]);
    modX_matrix = Rcpp::as<arma::mat>(l3["modX_matrix"]);
    next_beta = Rcpp::as<std::vector<double>>(l3["next_beta"]);
    old_active_data_decomp = Rcpp::as<arma::mat>(l3["old_active_data_decomp"]);
    active_beta = Rcpp::as<arma::vec>(l3["active_beta"]);
    gam_lasso = Rcpp::as<arma::vec>(l3["gam_lasso"]);

    machinePrec = l4["machinePrec"];
    actions = Rcpp::as<std::list<int>>(l4["actions"]);
    df = Rcpp::as<std::list<int>>(l4["df"]);

}

void tlarsCpp::executeLarsStep(int T_stop, bool earlyStop)
{

    // Determine the index of the first knockoff // and initialize knockoff counter
    int knockoff_ind = p - L_val;
    // int count_knockoffs = 0;
    //Begin LARS-algorithm
    while (k < max_steps&&
            count_inactive_pred > 0 &&
            count_active_pred < effective_n &&
            (count_knockoffs < T_stop || earlyStop == false))
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
        if (corr_max_inactive<100*machinePrec)
        {
            if(verbose == true)
            {
                // std::cout << "Stopped because of too low correlations";
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
                if(corr_predictors(*it)>= corr_max_inactive - machinePrec ||
                        corr_predictors(*it)<= -corr_max_inactive + machinePrec)
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

                    // Add 1 to the knockoff counter if the corresponding predictor was a knockoff.
                    if (*it>=knockoff_ind)
                    {
                        count_knockoffs++;
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
        Gi1 = solveUpperTriangular(active_data_decomp,solveLowerTriangular(active_data_decomp.t(),sign_vec));
        A = Gi1.t() * sign_vec;
        A = sqrt(1/A);
        w = (A*Gi1.t()).t();
        modX_matrix.resize(n,count_active_pred);
        counter=0;
        for (it = active_pred.begin(); it!= active_pred.end(); it++)
        {
            modX_matrix.col(counter) = X.col(*it);
            counter++;
        }
        u = modX_matrix*w;
        if(count_active_pred >= effective_n || count_active_pred >= p - count_ignored_pred)
            gamhat = corr_max_inactive/A(0,0);
        else
        {
            modX_matrix.resize(n,count_inactive_pred);
            counter = 0;
            for (it = inactive_pred.begin(); it!= inactive_pred.end(); it++)
            {
                modX_matrix.col(counter) = X.col(*it);
                counter++;
            }
            a = (u.t()*modX_matrix).t();
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
                // std::cout << "Error: No positve Gamma";
                break;
            }
            for(i=0; i<count_inactive_pred; i++)
            {
                if(gamhat1(i)<gamhat && gamhat1(i)>=machinePrec)
                {
                    gamhat=gamhat1(i);
                }
                if(gamhat2(i)<gamhat && gamhat2(i)>=machinePrec)
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
                if(gam_lasso(i)<gamhat && gam_lasso(i)>=machinePrec)
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
                    if (*it>=knockoff_ind)
                    {
                        count_knockoffs--;
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


void tlarsCpp::update_decomp(arma::mat new_X, arma::mat old_X)
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
        r = solveLowerTriangular(active_data_decomp.t(), Xtx);
        double rpp = pow(norm_xnew, 2);
        for(j=0; j<active_data_rank; j++)
        {
            rpp = rpp - pow(r(j),2);
        }
        // Check for machine singularity
        if(rpp<=machinePrec) rpp = machinePrec;
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

// Remove a variable at removal_index in the cholesky decomposition
void tlarsCpp::remove_var_from_decomp(int removal_index)
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
        modified_submatrix = choleskyDecomp(modified_submatrix + modified_submatrix2);
        active_data_decomp.shed_cols(removal_index, removal_index);
        active_data_decomp.shed_rows(count_active_pred-1,count_active_pred-1);
        active_data_decomp.submat(removal_index, removal_index, count_active_pred-2, count_active_pred-2) = modified_submatrix;
        active_data_rank--;
    }
}

arma::vec tlarsCpp::solveUpperTriangular(arma::mat upperT_R, arma::vec x)
{
    int n_solve = x.n_elem;
    int max_index_solve = n_solve-1;
    for (i=max_index_solve; i>=0; i--)
    {
        for(j=max_index_solve; j>i; j--)
        {
            x(i) = x(i) - upperT_R(i,j) * x(j);
        }
        x(i) = x(i)/upperT_R(i,i);
    }
    return x;
}

arma::vec tlarsCpp::solveLowerTriangular(arma::mat lowerT_R, arma::vec x)
{
    int n_solve = x.n_elem;
    for (i=0; i<n_solve; i++)
    {
        for(j=0; j<i; j++)
        {
            x(i) = x(i) - lowerT_R(i,j) * x(j);
        }
        x(i) = x(i)/lowerT_R(i,i);
    }
    return x;
}

arma::mat tlarsCpp::choleskyDecomp(arma::mat square_matrix)
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

arma::vec tlarsCpp::doubleListToVector(std::list<double> doubleList)
{
    arma::vec output(doubleList.size());
    counter= 0;
    for (double_it = doubleList.begin(); double_it!= doubleList.end(); double_it++)
    {
        output(counter) = *double_it;
        counter++;
    }
    return output;
}

arma::vec tlarsCpp::intListToVector(std::list<int> intList)
{
    arma::vec output(intList.size());
    counter= 0;
    for (it = intList.begin(); it!= intList.end(); it++)
    {
        output(counter) = *it;
        counter++;
    }
    return output;
}


void tlarsCpp::updateDf()
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

