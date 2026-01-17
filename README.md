
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dream: An R Package for Dynamic Relational Event Analysis and Modeling

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/dream?color=blue)](https://cran.r-project.org/package=dream)
[![](http://cranlogs.r-pkg.org/badges/grand-total/dream?color=red)](https://cran.r-project.org/package=dream)
[![CRAN
checks](https://badges.cranchecks.info/summary/dream.svg)](https://cran.r-project.org/web/checks/check_results_dream.html)

[![R-CMD-check](https://github.com/kevinCarson/dream/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kevinCarson/dream/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `dream` package provides users with helpful functions for (large)
relational event modeling/analysis. For an introduction to relational
events analysis/modeling, see [Butts
(2008)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2008.00203.x),
[Butts et
al. (2023)](https://www.cambridge.org/core/journals/network-science/article/relational-event-models-in-network-science/D3C9C3D18FC9C3FCC8561CCC4FEE747A),
[Duxbury
(2023)](https://methods.sagepub.com/book/mono/longitudinal-network-models/toc#_),
and [Bianchi et
al. (2024)](https://www.annualreviews.org/content/journals/10.1146/annurev-statistics-040722-060248).
In particular, `dream` provides users with helper functions for large
relational event analysis, such as recently proposed sampling procedures
for creating relational risk sets (i.e., sampling from the observed
event sequence [(Lerner and Lomi
2020)](https://www.cambridge.org/core/journals/network-science/article/reliability-of-relational-event-model-estimates-under-sampling-how-to-fit-a-relational-event-model-to-360-million-dyadic-events/B4286F370CD3A1A4ED30DF5120F04897),
case-control sampling [(Vu et
al. 2015)](https://www.sciencedirect.com/science/article/abs/pii/S0378873315000477)).
Alongside the set of functions for relational event analysis, this
package includes functions for the structural analysis of one- and
two-mode networks, such as network constraint and effective size
measures.

This package was developed with support from the National Science
Foundation’s (NSF) Human Networks and Data Science Program (HNDS) under
award number 2241536 (PI: Diego F. Leal). Any opinions, findings, and
conclusions, or recommendations expressed in this material are those of
the authors and do not necessarily reflect the views of the NSF.

## Authors

**Kevin A. Carson**  
- Author & Maintainer  
- PhD Candidate at the [University of Arizona School of
Sociology](https://sociology.arizona.edu/)  
- Email: <kacarson@arizona.edu>  
- Website: <https://kevincarson.github.io/>

**Diego F. Leal**  
- Author & Maintainer  
- Associate Professor at the [University of Arizona School of
Sociology](https://sociology.arizona.edu/)  
- Email: <dflc@arizona.edu>  
- Website: <https://www.diegoleal.info/index.html>

## Installation

You can install the stable verison of `dream` from
[CRAN](https://cran.r-project.org/package=dream) via:

``` r
install.packages("dream")
```

You can install the development version of `dream` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kevinCarson/dream")
```

## The `dream` Package API

The dream package ‘API’ is structured into six categories, where the
prefix identifies what category the specific function corresponds to
(see below):

- `remstats_`
- `netstats_om_`
- `netstats_tm_`
- `estimate_`
- `simulate_`
- `create_`

The `remstats_` functions compute relational/network statistics for
relational event sequences. For instance, `remstats_fourcycles` computes
the four-cycles network statistic for a two-mode relational event
sequence. The `create_` function creates a risk-set for one- and
two-mode relational event sequences based on a set of sampling
procedures. The `netstats_om_` series of functions compute static
network statics for one-mode networks (i.e., `netstats_om_pib` computes
[Leal
(2025)](https://journals.sagepub.com/doi/10.1177/00491241251322517)
measure for potential for intercultural brokerage). The `netstats_tm_`
set of functions compute static network statistics for two-mode networks
(i.e., `netstats_tm_effective` computes [Burchard and Cornwell
(2018)](https://www.sciencedirect.com/science/article/abs/pii/S0378873317302241)
measure for two-mode ego effective size). The `estimate_` function
estimates a relational event models for relational event sequences.
Currently, the only function in this set is `estimate_rem_logit`, which
estimates the ordinal timing relational event model and, under certain
conditions, can estimate a Cox-proportional hazard model for exact
timing relational event models (see [Bianchi et
al. (2024)](https://www.annualreviews.org/content/journals/10.1146/annurev-statistics-040722-060248)
and [Butts
(2008)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2008.00203.x)
for more information on these models). Finally, the `simulate_`
functions simulate one-mode relational event sequences based upon
results of a relational event model.

## Estimating an (Ordinal) Relational Event Model in `dream`

### Sampling from the Observed Events and Case-Control Sampling

This is a basic example which shows how to sample from the observed
events and employ the case-control sampling technique for large
relational event models (see [Butts
2008](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2008.00203.x))
following [Lerner and Lomi
(2020)](https://www.cambridge.org/core/journals/network-science/article/reliability-of-relational-event-model-estimates-under-sampling-how-to-fit-a-relational-event-model-to-360-million-dyadic-events/B4286F370CD3A1A4ED30DF5120F04897)
and [Vu et
al. (2015)](https://www.sciencedirect.com/science/article/abs/pii/S0378873315000477).
Then based upon the post-processing event sequence, the example computes
a set of standard network statistics for two-mode relational event
models. Lastly, the examples estimates an ordinal timing relational
event model. The event sequence included in this example is based a
subset (i.e., the first 100,000 events) of the 2018 Wikipedia
article-edit event sequence used in [Lerner and Lomi
(2020)](https://www.cambridge.org/core/journals/network-science/article/reliability-of-relational-event-model-estimates-under-sampling-how-to-fit-a-relational-event-model-to-360-million-dyadic-events/B4286F370CD3A1A4ED30DF5120F04897).
Across five replications, the average execution time for the example
below on a standard MacBook Air was 46.392 seconds.

``` r
library(dream)
data("WikiEvent2018.first100k")
WikiEvent2018.first100k$time <- as.numeric(WikiEvent2018.first100k$time)
### Creating the EventSet By Employing Case-Control Sampling With M = 10 and
### Sampling from the Observed Event Sequence with P = 0.01
EventSet <- create_riskset(
  type = "two-mode",
  time = WikiEvent2018.first100k$time, # The Time Variable
  eventID = WikiEvent2018.first100k$eventID, # The Event Sequence Variable
  sender = WikiEvent2018.first100k$user, # The Sender Variable
  receiver = WikiEvent2018.first100k$article, # The Receiver Variable
  p_samplingobserved = 0.10, # The Probability of Selection
  n_controls = 10, # The Number of Controls to Sample from the Full Risk Set
  seed = 9999) # The Seed for Replication

post.processing.riskset <- EventSet[EventSet$sampled == 1,] #only those sampled events! 
```

``` r
nrow(post.processing.riskset) #the total number of post-processing events
#> [1] 110000
table(post.processing.riskset$observed) # 0 = null events; 1 = observed events
#> 
#>      0      1 
#> 100000  10000
```

Based on the above results, the post-processing event sequence contains
110,000 post-processing events, that is, 10,000 observed events and 10
control events per observed events (i.e., 100,000 null events).

### A Miniature Replication of Lerner and Lomi (2020)

``` r
# computing the inertia statistic with the exponential weights and a halflife
# value of 30 days
post.processing.riskset$repetition <- remstats_repetition(
   time = EventSet$time,
   sender = EventSet$sender,
   receiver = EventSet$receiver,
   sampled = EventSet$sampled,
   observed = EventSet$observed,
   halflife = 2.592e+09, 
   dyadic_weight = 0,
   exp_weight_form = FALSE)

# computing the sender outdegree statistic with the exponential weights and a halflife
# value of 30 days
post.processing.riskset$sender.outdegree <- remstats_degree(
   formation = "sender-outdegree",
   time = EventSet$time,
   sender = EventSet$sender,
   receiver = EventSet$receiver,
   sampled = EventSet$sampled,
   observed = EventSet$observed,
   halflife = 2.592e+09, 
   dyadic_weight = 0,
   exp_weight_form = FALSE)

# computing the receiver indegree statistic with the exponential weights and a halflife
# value of 30 days
post.processing.riskset$receiver.indegree <- remstats_degree(
   formation = "receiver-indegree",
   time = EventSet$time,
   sender = EventSet$sender,
   receiver = EventSet$receiver,
   sampled = EventSet$sampled,
   observed = EventSet$observed,
   halflife = 2.592e+09, 
   dyadic_weight = 0,
   exp_weight_form = FALSE)

# computing the four-cycles statistic with the exponential weights and a halflife
# value of 30 days
post.processing.riskset$fourcycles <- remstats_fourcycles(
   time = EventSet$time,
   sender = EventSet$sender,
   receiver = EventSet$receiver,
   sampled = EventSet$sampled,
   observed = EventSet$observed,
   halflife = 2.592e+09, 
   dyadic_weight = 0,
   exp_weight_form = FALSE)

# Estimating the ordinal relational event model! 
lerner.lomi.rem <- estimate_rem_logit(observed ~ 
                            repetition +
                            sender.outdegree + 
                            receiver.indegree + 
                            receiver.indegree:sender.outdegree +
                            fourcycles,
                            event.cluster = post.processing.riskset$eventID,
                            newton.rhapson=FALSE,
                            data = post.processing.riskset)
#> Extracting user-provided data.
#> Prepping data for numerical optimization.
#> Starting optimzation for parameters.
```

``` r
summary(lerner.lomi.rem)
#> Ordinal Timing Relational Event Model
#> 
#> Call:
#> estimate_rem_logit(formula = observed ~ repetition + sender.outdegree + 
#>     receiver.indegree + receiver.indegree:sender.outdegree + 
#>     fourcycles, event.cluster = post.processing.riskset$eventID, 
#>     data = post.processing.riskset, newton.rhapson = FALSE)
#> 
#>  n events: 10000 null events: 1e+05 
#> 
#> Coefficients:
#>                                    Estimate Std. Error z value  Pr(>|z|)    
#> repetition                           7.1365     0.3970  17.976 < 2.2e-16 ***
#> sender.outdegree                     0.0219     0.0005  44.825 < 2.2e-16 ***
#> receiver.indegree                    0.1580     0.0093  16.995 < 2.2e-16 ***
#> fourcycles                           1.3161     0.0465  28.334 < 2.2e-16 ***
#> sender.outdegree:receiver.indegree  -0.0009     0.0001 -10.218 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Null Likelihood: -23978.95 Model Likelihood: -7558.12 
#> 
#> Likelihood Ratio Test: 32841.67  with df: 5 p-value: 0 
#> 
#> AIC 15126.24 BIC 15162.29
```

## Questions, Comments, or Suggestions!

If you have any questions, comments, or suggestions please feel free to
open an issue!
