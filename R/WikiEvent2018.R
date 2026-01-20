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
