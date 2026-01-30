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


