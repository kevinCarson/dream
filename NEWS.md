# dream 2.1.5 (2026-07-17)

## Minor Changes
* The output from the `gof_rem()` function is now an `dream_gof` S3 object. The following
methods are included for this object class: `plot()`, `summary()`, and `print()`. 
* Added a new relational event sequence `data.frame` object to help users get started
with one-mode relational event analysis. The new object is named `simulated.res`.
* Happy dreaming!



# dream 2.1.4 (2026-06-24)

## Minor Changes
* Set a default value for the `plot()` S3 method for `dream_rem` objects, so that the
user does not have to initially specify a value. 
* Updated the `predict()` S3 method for `dream_rem` objects to now allow
updated `dream_sequence` objects to based the predicted event rates on. 
* Updated the `print()` and `summary()` S3 methods for `dream_sequence` objects to 
relay additional information about the processed relational event sequence to the user. 
* Fixed a small bug in the `simulate_rem_seq()` function for the three-cycles statistic when 
parallel processing is requested. 


# dream 2.1.3 (2026-06-10)

## Minor Changes
* Fixed a bug in the `gof_rem()` function that resulted in the proportion of receivers
correctly predicted to be incorrect, and added the `rseed` argument to allow
the user to set the random seed.
* Add the `cumulative` option to the `create_res()` function as an additional 
option. Case-control sampling is now indicated by the argument `case_control` 
in the `create_res()` function. 


# dream 2.1.2 (2026-06-02)

## Minor Changes
* Added three new functions for the examination of relational event model fits. First, `gof_rem()`
computes the proportion of correctly predicted dyadic events for a `dream_rem` relational event 
model fit. Secondly, we added the `plot()` and `residuals()` S3 methods for the
`dream_rem` object class. Please see the function documentation pages for more information
on these new functions. Happy dreaming!



# dream 2.1.1 (2026-05-27)

## Major Changes
* The package has undergone a new API for the analysis of relational event sequences 
aimed at reducing redundant user-inputs and increasing useability. Importantly, the 
package now contains two S3 object classes: `dream_sequence` and `dream_rem`, 
where `dream_sequence` is the object class created after using `create_res()` in the 
dream package to generate the post-processing relational event 
sequence (i.e., (sampled) realized events and sampled control/null events) or 
after using `dream_sequence()`, a new function to create the S3 object based 
upon a user-generated processing sequence. Moreover, in prior versions of 
the package, users had to input 6 arguments into the functions to compute 
various sufficient network statistics (e.g., four-cycles, repetition), these 
6 arguments are now pulled from the `dream_sequence` data object. Additionally, 
the functions that compute the sufficient network statistics add the vector of 
computed statistics to the user-inputted `dream_sequence` object as an additional
element in the `statistics` list. The `estimate_rem()` function now relies 
upon the user-created `dream_sequence` objects and extracts the relational event 
model covariates from the data object. Finally, `dream_rem` S3 objects now contain
a new set of functions to extract relevant model information: `vcov()`, `logLik()`, 
`coef()`, and `predict()`. 
* `create_riskset_dynamic()` and `create_riskset_constant()` have been replaced  
by the new function `create_res()`, which internally is the same as the aforementioned functions
and generate risksets based upon the `riskset` argument. We encourage users
to read the help page for the `create_res()` to see the updated set of full 
riskset possibilities. 
* As mentioned above, the `dreamstats_` series of functions have been updated
to reduce the number of arguments supplied by the user. Previously, the following
arguments needed to be specified: `time`, `sender`, `receiver`, `eventID`, 
`observed`, and `sampled`. However, the new updated functions internally create
these based upon the `data` argument which is an S3 `dream_sequence` object. 
* The newly created S3 `dream_sequence` object, which is created using the 
`create_res()` function, also has a function named `dream_sequence()` that
can create a `dream_sequence` object based upon the user-inputted values. In 
addition, the `as.data.frame()` has been supplied to extract the post-processing
relational event sequence and the computed statistics. Users may find this function
useful if they would like to estimate the relational event model by a different 
modeling function such as the `clogit()` function in the `survival` package or
a user-created function. 
* The new package version also includes a new suite of functions that expands the 
set of exogenous statistics for the estimation of relational event models. These are
`dreamstats_actor()`, `dreamstats_actorfe()`, `dreamstats_dyadic()`, 
`dreamstats_dyadfe()`, and `dreamstats_event()`. These functions, respectively, 
allow users to (1) add actor-level time-varying and time-invariant exogenous
covariates, (2) add actor-level (sender and receiver) fixed effects, (3) add 
dyadic-level time-varying and time-invariant exogenous covariates, (4) add
dyadic-level fixed effects, and (3) add event-level exogenous covariates. 
* Lastly, the `estimate_rem()` function now allows for the estimation of 
interval timing relational event models. In terms of estimation routines, the
c++ functions to compute the model log-likelihoods, hessians, and gradients 
are now computed using `RcppArmadillo` with the `arma` set of functions. 
* All functions that were deprecated in `dream` 1.1.1 were removed, in addition
to the past function arguments for the `dreamstats_` functions. 



# dream 1.1.1 (2026-03-24)

## Major Changes

* All `remstats_`-related functions were renamed to `dreamstats_` to provide a more
consistently named API for the `dream` package. As a result, all `remstats_` original functions are now 
deprecated starting on version 1.1.1, but remain available in this update. 
* `create_riskset()` was renamed to `create_riskset_constant()`. 
* All functions that were deprecated in `dream` 1.0.0 were removed.


# dream 1.0.1 (2026-01-30)

## Major Changes

* Implemented a new function to create time-dynamic risk sets for one- and two-mode
relational event sequences. The new function is named `create_riskset_dynamic()`. 

## Minor Changes

* Updated the `remstats_`-related `C++` files to delete past network edges (i.e., past 
event weights) that are smaller than the past event cut off value (i.e., `dyadic_weight`) when this feature is 
selected by the user. That is, when the `dyadic_weight` argument is non-zero.

# dream 1.0.0 (2026-01-19)

## Major Changes

* Added a `NEWS.md` file to track changes to the package.
* Two key changes occured in the new version of `dream` package. First, the majority 
of functions now depend upon c++ routines via the `Rccp` package to improve 
computational runtime. For instance, the new function to compute the repetition
network statistic, `remstats_repetition()`, is around 90 times faster than the previous 
function, `computeRepetition()`, in **dream** version 0.0.1. Secondly, the `dream` package
has an updated API to allow better categorization for the functions. Please run `?dream` for 
more details on the new API. Based upon this transition, the original functions are now 
deprecated starting on version 1.0.0, but remain available in this update. 
* The following functions were deprecated for `remstats_triads()`: `computeISP()`, 
`computeITP()`, `computeOSP()`, `computeTriads`, and `computeOTP()`. 
* The following functions were deprecated for `remstats_degree()`: `computeSenderOutdegree()`, 
`computeReceiverOutdegree()`, `computeSenderIndegree()`, and `computeReceiverIndegree()`. 
* `estimateREM()` was deprecated for `estimate_rem_logit()`.
* `simulateREseq()` was deprecated for `simulate_rem_seq()`.
* `computeFourCycles()` was deprecated for `remstats_fourcycles()`.
* `computeRemDyadCut()` was deprecated for `remstats_dyadcut()`.
* `computePersistence()` was deprecated for `remstats_persistence()`.
* `computePersistence()` was deprecated for `remstats_persistence()`.
* `processOMEventSeq()` and `processTMEventSeq()` were deprecated for `createriskset()`.
* `computePrefAttach()` was deprecated for `remstats_prefattachment()`.
* `computeRepetition()` was deprecated for `remstats_repetition()`.
* `computeRepetition()` was deprecated for `remstats_repetition()`.
* `remExpWeights()` was deprecated and will not be replaced in the current version as the 
updated functions do not require it. 
* `remExpWeights()` was deprecated and will not be replaced in the current version as the 
updated functions do not require it. 
* `computeTMDegree()` was deprecated for `netstats_tm_degreecent()`.
* `computeTMEgoDis()` was deprecated for `netstats_tm_egodistance()`.
* `computeBCConstraint()` was deprecated for `netstats_tm_constraint()`.
* `computeBCES()` was deprecated for `netstats_tm_effective()`.
* `computeBCRedund()` was deprecated for `netstats_tm_redundancy()`.
* `computeBurtsConstraint()` was deprecated for `netstats_om_constraint()`.
* `computeBurtsES()` was deprecated for `netstats_om_effective()`.
* `computeHomFourCycles()` was deprecated for `netstats_tm_homfourcycles()`.
* `computeLealBrokerage()` was deprecated for `netstats_om_pib()`.
* `computeNPaths()` was deprecated for `netstats_om_nwalks()`.
* `computeTMDens()` was deprecated for `netstats_tm_density()`.


