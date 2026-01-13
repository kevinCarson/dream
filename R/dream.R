#' dream: A Package for Dynamic Relational Event Analysis and Modeling
#'
#' The dream package provides users with helpful functions for relational
#' and event analysis. In particular, dream provides users with helper functions
#' for large relational event analysis, such as recently proposed sampling
#' procedures for creating relational risk sets. Alongside the set of functions
#' for relational event analysis, this package includes functions for the
#' structural analysis of one- and two-mode networks, such as network constraint
#' and effective size measures. This package was developed with support from the
#' National Science Foundation’s (NSF) Human Networks and Data Science Program (HNDS)
#' under award number 2241536 (PI: Diego F. Leal). Any opinions, findings,
#' and conclusions, or recommendations expressed in this material are those
#' of the authors and do not necessarily reflect the views of the NSF.
#'
#'
#' @section dream functions:
#'
#' The dream package 'API' is structured into six categories, where the prefix identifies what category the specific function corresponds to (see below):
#'
#'  - `remstats_`
#'  - `netstats_om_`
#'  - `netstats_tm_`
#'  - `estimate_`
#'  - `simulate_`
#'  - `create_`
#'
#'The `remstats_` functions compute relational/network statistics for relational event sequences.
#'For instance, `remstats_fourcycles` computes the four-cycles network statistic for a two-mode
#'relational event sequence. The `create_` function creates a risk-set for one- and two-mode
#'relational event sequences based on a set of sampling procedures. The `netstats_om_` series of functions compute
#'static network statics for one-mode networks
#'(i.e., `netstats_om_pib` computes [Leal (2025)](https://journals.sagepub.com/doi/10.1177/00491241251322517) measure for
#'potential for intercultural brokerage). The `netstats_om_` set of functions compute static network
#'statics for two-mode networks (i.e., `netstats_om_effective`
#'computes [Burchard and Cornwell (2018)](https://www.sciencedirect.com/science/article/abs/pii/S0378873317302241) measure for two-mode
#'ego effective size). The `estimate_` functions estimate relational event models for relational event sequences. Currently,
#'the only function in this set is `estimate_rem_logit`, which estimates the ordinal timing relational event model and,
#'under certain conditions, can estimate a Cox-proportional hazard model for exact timing relational event
#'models (see [Bianchi et al. (2024)](https://www.annualreviews.org/content/journals/10.1146/annurev-statistics-040722-060248)
#'and [Butts (2008)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2008.00203.x) for more information on
#'these models). Finally, the `simulate_` functions simulate one-mode relational event sequences based upon
#'results of a relational event model.
#'
#'
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @name dream
#' @useDynLib dream, .registration=TRUE
"_PACKAGE"
NULL
#> NULL
