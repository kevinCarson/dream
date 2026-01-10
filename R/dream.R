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
#' The functions in dream can be grouped into four useful categories:
#'
#' - Create Dynamic Risk Sets for (Large) Relational Event Models with \code{\link{createriskset}}.
#' - Compute Network Statistics for (Large) Relational Event Models
#'   - Functions: \code{\link{degreestats}}, \code{\link{triadicstats}}, \code{\link{fourcycles}},
#'   \code{\link{persistence}}, \code{\link{prefattachment}}, \code{\link{recency}},
#'   \code{\link{reciprocity}}, \code{\link{remdyadcut}}, and \code{\link{repetition}}.
#' - Estimate and Simulate (Large) Relational Event Models
#'   - Functions: \code{\link{remlogit}} and \code{\link{remsimulate}}.
#' - Compute One- and Two-Mode Network Structural Measures
#'   - Functions: \code{\link{bctmconstraint}}, \code{\link{bctmeffective}}, \code{\link{bctmredundancy}},
#'    \code{\link{burtsconstraint}}, \code{\link{burtseffective}}, \code{\link{homfourcycles}},
#'    \code{\link{lealpib}}, \code{\link{nwalks}}, \code{\link{tmdegreecent}},
#'    \code{\link{tmdensity}}, and \code{\link{tmegodistance}}.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @name dream
#' @useDynLib dream, .registration=TRUE
"_PACKAGE"
NULL
#> NULL
