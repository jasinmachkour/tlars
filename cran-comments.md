## Package update
This is an update of the tlars R package:

  * Version is increased to 1.0.0 (stable version).
  * Supporting scientific papers were updated.
  * C++11 specification was dropped.
  * DESCRIPTION, NEWS, README were updated.
  * Package was polished (vignettes, etc.).
  

## Resubmission 2
This is a resubmission that addresses the comments on the previous submission. The following changes were made:

 * A bug in resetting the user's options() in the README file and the vignette was fixed.
 * Two unit tests testing that tlars_model() and tlars() also work in low-dimensional settings were added.
 
## Resubmission 1
This is a resubmission that addresses the comments on the initial submission. The following changes were made:

 * The title is reduced to 65 characters.
 * Package names are written in single quotes.
 * References in the description of the DESCRIPTION file are written in the formats authors (year) <doi:...> and authors (year) <arXiv:...>.
 * \value has been added to all .Rd files.
 * User's options() are reset in the README file, examples, and the vignette.

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.

## Downstream dependencies
There are currently no downstream dependencies for this package.
