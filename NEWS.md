# dream 1.0.0 (2026-01-09)

## Major Changes

* Added a `NEWS.md` file to track changes to the package.
* Superseded the following functions in the prior version: computeISP(), 
computeITP(), computeOSP(), and computeOTP(). The new function that compute 
triadic formation statistics is now triadicstats(). Please see the help package for the
new function. 
* Superseded the following functions in the prior version: computeSenderIndegree(), 
computeSenderOutdegree(), computeReceiverOutdegree(), and computeReceiverIndegree(). The new function that compute 
triadic formation statistics is now degreestats(). Please see the help package for the
new function. 
* estimateREM() is now named remlogit().
* simulateREseq() is now named remsimulate().
* computeFourCycles() is now named fourcycles().
* computeRemDyadCut() is now named remdyadcut().
* computePersistence() is now named persistence().
* Superseded the following functions in the prior version: processOMEventSeq(), 
and processTMEventSeq(). The new function that creates risk sets is now createriskset(). Please see the help package for the
new function. 
* computePrefAttach() is now named prefattachment().
* computeRepetition() is now named repetition().
* computeTMDegree() is now named tmdegreecent().
* computeTMEgoDis() is now named tmegodistance().
* computeTriads() was removed as it returned the same values as computeOTP(). 
* remExpWeights() was removed as it is now longer required. 
* computeBCConstraint() is now named bctmconstraint().
* computeBCES() is now named bctmeffective().
* computeBCRedund() is now named bctmredundancy().
* computeBurtsConstraint() is now named burtsconstraint().
* computeBurtsES() is now named burtseffective().
* computeHomFourCycles() is now named homfourcycles().
* computeLealBrokerage() is now named lealpib().
* computeNPaths() is now named nwalks().
* computeTMDens() is now named tmdensity().
* computeTMEgoDis() is now named tmegodistance().

* A key change to the new version of **dream** is that most functions now use c++ functions
via the *Rccp* package to improve computational speed. For instance, the new version of 
repetition() that sources to c++ is around 90 times faster than the original coding of 
computeRepetition() in **dream** version 0.0.1. 


