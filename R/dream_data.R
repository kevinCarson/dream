# the master file for dream R package functions related to data files included


#' Wikipedia Edit Event Sequence 2018
#'
#'The first 100,000 events of the (two-mode) Wikipedia edit event sequence, where
#'an event is described as a Wikipedia user editing a Wikipedia article. The user
#'column represents the unique event senders, the article column represents the unique
#'event receivers (targets), and the time variable is in milliseconds.
#'
#' @format ## `WikiEvent2018.first100k`
#' The first 100,000 events of the Wikipedia edit event sequence, where an event
#' is described as a Wikipedia user editing a Wikipedia article. The user column
#' represents the unique event senders, the article column represents the unique
#' event receivers (targets), and the time variable is in milliseconds.
#' \describe{
#'   \item{user}{the column that represents the unique event senders.}
#'   \item{article}{the article column represents the unique event receivers.}
#'   \item{time}{the event time variable in milliseconds.}
#'   \item{eventID}{the numerical id for each event in the event sequence}
#' }
#'
#' @usage data(WikiEvent2018.first100k)
#' @source <https://zenodo.org/records/1626323>
#' @source  Lerner, Jurgen and Alessandro Lomi. 2020. "Reliability of relational event model estimates
#' under sampling: how to fit a relational event model to 360 million dyadic events."
#' *Network Science* 8(1):97-135. (DOI: https://doi.org/10.1017/nws.2019.57)

"WikiEvent2018.first100k"




#' Davis Southern Women's Dataset
#'
#'
#' @format ## `southern.women`
#' Two-Mode affliation matrix from Davis et al.(1941) Southern Women study. 18 women x 14 events.
#' Dataset is taken from the networkdata R package (Almquist 2014)
#' @usage data(southern.women)
#' @source  Almquist, Zach. 2014. _networkdata: Lin Freeman's Network Data Collection_. R package version 0.01, <https://github.com/Z-co/networkdata>.
#' @source  Brieger, Ronald. 1974. "Duality of Persons and Groups." *Social Forces* 53(2): 181-190.
#' @source  Davis, Allison, Burleigh B. Gardner, and Mary R. Gardner. 1941. *Deep South: A Social Anthropological Study of Caste and Class*. University of Chicago Press.
#
"southern.women"


