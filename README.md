
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dream: An R Package for Dynamic Relational Event Analysis and Modeling

<!-- badges: start -->

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

You can install the developmental version of `dream` from
[GitHub](https://github.com/) with:

``` r
remotes::install_github("kevinCarson/dream")
```

## The `dream` Package API

The dream package ‘API’ is structured into six categories, where the
prefix identifies what category the specific function corresponds to
(see below):

- `dreamstats_`
- `netstats_om_`
- `netstats_tm_`
- `estimate_rem`
- `simulate_rem_seq`
- `create_res`’

The `dreamstats_` functions compute relational/network statistics for
relational event sequences. For instance, `dreamstats_fourcycles`
computes the four-cycles network statistic for a two-mode relational
event sequence. The `create_res` function creates a risk-set for one-
and two-mode relational event sequences based on a set of sampling
procedures. The `netstats_om_` series of functions compute static
network statics for one-mode networks (i.e., `netstats_om_pib` computes
[Leal
(2025)](https://journals.sagepub.com/doi/10.1177/00491241251322517)
measure for potential for intercultural brokerage). The `netstats_om_`
set of functions compute static network statics for two-mode networks
(i.e., `netstats_om_effective` computes [Burchard and Cornwell
(2018)](https://www.sciencedirect.com/science/article/abs/pii/S0378873317302241)
measure for two-mode ego effective size). The `estimate_rem` functions
estimate relational event models for relational event sequences. The
function estimates the interval and ordinal timing relational event
model and, under certain conditions, can estimate a Cox-proportional
hazard model for exact timing relational event models (see [Bianchi et
al. (2024)](https://www.annualreviews.org/content/journals/10.1146/annurev-statistics-040722-060248)
and [Butts
(2008)](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-9531.2008.00203.x)
for more information on these models). Finally, the `simulate_rem_seq`
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

``` r
library(dream)
data("WikiEvent2018.first100k", package = "dream")
WikiEvent2018.first100k$time <- as.numeric(WikiEvent2018.first100k$time)
### Creating the EventSet By Employing Case-Control Sampling With M = 10 and
### Sampling from the Observed Event Sequence with P = 0.01
processed <- create_res(
  type = "two-mode",
  ordinal = TRUE,
  riskset = "cumulative",
  time = WikiEvent2018.first100k$time, # The Time Variable
  sender = WikiEvent2018.first100k$user, # The Sender Variable
  receiver = WikiEvent2018.first100k$article, # The Receiver Variable
  p_samplingobserved = 0.20, # The Probability of Selection
  case_control=TRUE,
  n_controls = 10, # The Number of Controls to Sample from the Full Risk Set
  seed = 9999) # The Seed for Replication
```

``` r
processed #printing the summary information
#> Processed Relational Event Sequence: 
#>  -> Relational event sequence type: two-mode 
#>  -> Relational event sequence timing: ordinal 
#>  -> Number of senders: 4041 
#>  -> Number of receivers: 24936 
#>  -> Number of realized events: 1e+05 
#>  -> Sampling from the realized event sequence?: yes 
#>  -> The probability of sampling from the realized event sequence: 0.2 
#>  -> Number of sampled realized events: 20000  
#>  -> The risk set definition: cumulative 
#>  -> Case-control sampling of control events?: yes 
#>  -> Case-control sampling m: 10 
#>  -> Number of non-realized (i.e., control) events: 2e+05  
#>  -> The total number of processed realized and non-realized events: 3e+05
```

Based on the above results, the post-processing event sequence contains
110,000 post-processing events, that is, 10,000 observed events and 10
control events per observed events (i.e., 100,000 null events).

### A Miniature Replication of Lerner and Lomi (2020)

``` r
# computing the inertia statistic with the exponential weights and a halflife
# value of 30 days
processed <- dreamstats_repetition(data = processed, 
                                   halflife = 2.592e+09, 
                                   dyadic_weight = 0.01)

# computing the sender outdegree statistic with the exponential weights and a halflife
# value of 30 days
processed <- dreamstats_degree(formation = "sender-outdegree",
                               data = processed,
                               halflife = 2.592e+09, 
                               dyadic_weight = 0.01)

# computing the receiver indegree statistic with the exponential weights and a halflife
# value of 30 days
processed <- dreamstats_degree(formation ="receiver-indegree",
                               data = processed,
                               halflife = 2.592e+09, 
                               dyadic_weight = 0.01)

# computing the four-cycles statistic with the exponential weights and a halflife
# value of 30 days
processed <- dreamstats_fourcycles(data = processed, 
                                   halflife = 2.592e+09, 
                                   dyadic_weight = 0.01)

#transforming the variables following Lerner and Lomi (2020) 
processed$statistics$four.cycles <- log1p(processed$statistics$four.cycles)
processed$statistics$sender.outdegree <- log1p(processed$statistics$sender.outdegree)
processed$statistics$receiver.indegree <- log1p(processed$statistics$receiver.indegree)
processed$statistics$repetition <- log1p(processed$statistics$repetition)

#some realized events have the same time point, so to make them ordered
#in the sequence of the realized relational event sequence per Lerner and Lomi (2020)
processed$processed_sequence$time <- processed$processed_sequence$eventID

# Estimating the ordinal relational event model! 
lerner.lomi.rem <- estimate_rem(formula= ~ repetition + sender.outdegree + 
                            receiver.indegree + four.cycles +
                            receiver.indegree:sender.outdegree,
                            data = processed)
```

``` r
summary(lerner.lomi.rem)
#> Ordinal Timing Relational Event Model
#> 
#> Call:
#> estimate_rem(formula = ~repetition + sender.outdegree + receiver.indegree + 
#>     four.cycles + receiver.indegree:sender.outdegree, data = processed)
#> 
#> Relational Event Sequence Information:
#> The number of realized events = 20000 
#> The number of control events = 2e+05 
#> 
#> Coefficients:
#>                                     Estimate Std. Error z value Pr(>|z|)    
#> repetition                          8.843061   0.379085 23.3274  < 2e-16 ***
#> sender.outdegree                    1.391764   0.015132 91.9769  < 2e-16 ***
#> receiver.indegree                   1.206558   0.044584 27.0624  < 2e-16 ***
#> four.cycles                         0.328161   0.079055  4.1511  3.3e-05 ***
#> sender.outdegree:receiver.indegree -0.082913   0.020910 -3.9651  7.3e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> REM Fit Information:
#> Null Model Likelihood = -47957.91 
#> Full Model Likelihood = -7578.563 
#> Likelihood Ratio Test: 80758.69 (df=5; p-value: 0)
#> AIC = 15167.13 
#> BIC = 15206.64 
#> Number of Newton Iterations = 8
```

## Questions, Comments, or Suggestions!

If you have any questions, comments, or suggestions please feel free to
open an issue!
