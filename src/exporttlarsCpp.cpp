// exporttlarsCpp.cpp

#include <RcppArmadillo.h>
#include "tlarsCpp.h"
using namespace Rcpp;
using namespace RcppArmadillo;


RCPP_MODULE(tlarsCpp){
class_<tlarsCpp>("tlarsCpp")
    .constructor<arma::mat, arma::vec , bool , bool , bool , int, std::string>()
    .constructor<Rcpp::List>()
    .method("getAll", &tlarsCpp::getAll)
    .method("getBeta", &tlarsCpp::getBeta)
    .method("getLastBeta", &tlarsCpp::getLastBeta)
    //.field("R_rank", &tlarsCpp::R_rank)
    .method("getActive_data_decomp", &tlarsCpp::getActive_data_decomp)
    .method("getActive_data_rank", &tlarsCpp::getActive_data_rank)
    .method("getR2", &tlarsCpp::getR2)
    .method("getRSS", &tlarsCpp::getRSS)
    .method("getActions", &tlarsCpp::getActions)
    .method("getDf", &tlarsCpp::getDf)
    .method("getCp", &tlarsCpp::getCp)
    .method("getEntry", &tlarsCpp::getEntry)
    .method("getGamrat", &tlarsCpp::getGamrat)
    .method("getGamHat", &tlarsCpp::getGamHat)
    .method("getMu", &tlarsCpp::getMu)
    .method("getNormX", &tlarsCpp::getNormX)
    .method("getMeanX", &tlarsCpp::getMeanX)
    .method("executeLarsStep", &tlarsCpp::executeLarsStep);
}



