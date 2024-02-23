## Changes in tlars version 1.0.1 2024-02-24

* Package update:

  - Replace the-annals-of-statistics.csl with ieee.csl to fix a minor vignette re-building error on some linux flavors.

## Changes in tlars version 1.0.0 2024-02-20

* Package update (stable version):

  - Update supporting scientific papers.
  - Polish previous version.
  - Drop C++11 specification.

## Changes in tlars version 0.0.1 2022-07-11

* First release:

  - Polish pre-release version.

## Changes in tlars version 0.0.0.9000 2022-05-18

* Initial pre-release version contains C++ class: 

	- tlars_cpp.

* It includes two major functions for initializing a T-LARS model and executing the T-LARS algorithm:

	- tlars_model(); and
	- tlars().

* It includes two functions for visualizing the T-LARS results:

	- plot.Rcpp_tlars_cpp(); and
	- print.Rcpp_tlars_cpp(). 
	
* Vignette illustrates how the tlars package is used.
