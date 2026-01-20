#######
# The Creation of Interpersonal Cultural Brokerage Statistics Based on
# Reference: Leal, Diego F. Forthcoming. "Cultural Brokerage and Diffusion Across Bright Symbolic
# Boundaries. A New Measure to Locate Brokers Spanning Cultural Holes"
# Reference: Gould, Roger V., and Roberto M. Fernandez. 1989. "Structures of Mediation: A Formal Approach to Brokerage in Transaction Networks"
#            Sociological Methodology, Vol. 19, pp. 89-126
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 09-21-24

#' @title Compute Potential for Intercultural Brokerage (PIB) Based on Leal (2025)
#' @name netstats_om_pib
#' @param net The one-mode adjacency matrix.
#' @param g.mem The vector of membership values that the brokerage scores will be based on.
#' @param symmetric TRUE/FALSE. TRUE indicates that network matrix will be treated as symmetric. FALSE indicates that the network matrix will be treated as asymmetric. Set to TRUE by default.
#' @param triad.type The string value (or vector) that indicates what specific triadic (star) structures the potential for cultural brokerage will be computed for. Possible values are "ANY", "OTS", "ITS", "MTS" (see the details section). The function defaults to “ANY”.
#' @param count TRUE/FALSE. TRUE indicates that the number of culturally brokered open triangles will be returned. FALSE indicates that the proportion of culturally brokered open triangles to all open triangles will be returned (see the details section). Set to TRUE by default.
#' @param isolate If count = FALSE, the numerical value that will be given to isolates. This value is set to NA by default, as 0/0 is undefined. The user can specify this value!
#' @import Rcpp
#' @return The vector of interpersonal cultural brokerage values for the one-mode network.
#' @export
#'
#'
#'
#'
#'
#' @description
#' `r lifecycle::badge("stable")`
#'
#' Following Leal (2025), this function calculates node’s Potential for
#' Intercultural Brokerage (PIB) in a one-mode network, that is, brokerage
#' based on nodes’ distinct group memberships. For example, users can examine
#' PIB based on actors’ gender. The option count determines what is returned
#' by the function. If count is TRUE, then the count of ‘culturally’ dissimilar
#' pairs brokered by ego is included (i.e., ego’s total count of brokered open
#' triangles where the alters at the two endpoints of said open triangles
#' are ‘culturally’ dissimilar from one another). If count is FALSE, the
#' proportion of ego’s brokered open triangles where the endpoints
#' are ‘culturally’ dissimilar out of all of ego’s brokered open
#' triangles (regardless of the cultural identity of the alters) is returned. The
#' formula for computing interpersonal brokerage is presented in the details section.
#'
#' @details
#' Following Leal (2025), the formula for interpersonal brokerage is:
#'
#' \deqn{ \text{PIB}_i = \sum_{j < k} \frac{S_{jik}}{S_{jk}} m_{jk}, \quad S_{jik} \neq 0 \text{ and } i \neq j \neq k }
#'
#' where:
#' \itemize{
#'   \item \eqn{S_{jik} = 1} if there is an (un)directed two-path connecting actors *j* and *k* through actor *i*; 0 otherwise.
#'   \item \eqn{m_{jk} = 1} if actors *j* and *k* are on different sides of a symbolic boundary; 0 otherwise.
#'   \item Following Gould (1989), \eqn{S_{jik}} represents the total number of two-paths between actors *j* and *k*.
#' }
#' If the network is non-symmetric (i.e., the user specified symmetric = FALSE), then
#' the function can compute the cultural brokerage scores for different star structures. The possible
#' values are: "ANY", which computes the scores for all structures, where a tie exists
#' between *i* and *j*, *j* and *k*, and one does not exist between *i* and *k*. "OTS" computes the
#' values for outgoing two-stars (i<-j->k or the 021D triad according to the M.A.N. notation; see Wasserman and Faust 1994), where j is the broker. "ITS" computes the
#' values for incoming two-stars (i->j<-k or the 021U triad according to the M.A.N. notation; see Wasserman and Faust 1994), where j is the broker. "MTS" computes PIB for
#' mixed triadic structures (i<-j<-k or i->j->k or the 021C triad according to the M.A.N. notation; see Wasserman and Faust 1994). If not specified, the function defaults to the
#' "ANY" category. This function can also compute all of the formations at once.
#'
#' @author Kevin A. Carson <kacarson@arizona.edu>, Diego F. Leal <dflc@arizona.edu>
#' @references
#' Gould, Roger. 1989. "Power and Social Structure in Community Elites." *Social Forces* 68(2): 531-552.
#'
#' Leal, Diego F. 2025. "Locating Cultural Holes Brokers in Diffusion Dynamics Across Bright Symbolic Boundaries." *Sociological Methods & Research* \doi{10.1177/00491241251322517}
#'
#' Wasserman, Stanley and Katherine Faust. 1994. *Social Network Analysis: Methods and Applications*. Cambridge: Cambridge University Press.
#'
#' @examples
#'
#' # For this example, we recreate Figure 3 in Leal (2025)
#' LealNet <- matrix( c(
#'  0,1,0,0,0,0,0,
#'  1,0,1,1,0,0,0,
#'  0,1,0,0,1,1,0,
#'  0,1,0,0,1,0,0,
#'  0,0,1,1,0,0,0,
#'  0,0,1,0,0,0,1,
#'  0,0,0,0,0,1,0),
#'  nrow = 7, ncol = 7, byrow = TRUE)
#'
#' colnames(LealNet) <- rownames(LealNet) <- c("A", "B", "C","D",
#'                                            "E", "F", "G")
#' categorical_variable <- c(0,0,1,0,0,0,0)
#' #These values are exactly the same as reported by Leal (2025)
#' netstats_om_pib(LealNet,
#'    symmetric = TRUE,
#'    g.mem = categorical_variable)
#'
#'
#'
netstats_om_pib <- function(net,
                 g.mem,
                 symmetric = TRUE,
                 triad.type = NULL, #can be of c("ANY", "OTS", "ITS", "MTS")
                 count = TRUE,
                 isolate = NA){

  if(symmetric == FALSE){
    message("Note: the network is specified as symmetric by default. However, you
           specified the network as asymmetric. Therefore, if traid.type was
           not specified, the function defaults to 'ANY'.")}
  if(symmetric == FALSE & is.null(triad.type)){ #if the triad type was not specifed,
    #simply compute all of the possible starst
    triad.type <- "ANY"
  }
  if(sum(!(triad.type %in% c("ANY", "OTS", "ITS", "MTS"))) != 0){
    base::stop("an unknown type was specifed in the triad.type argument. Please
               update this and retry!")
  }
  if(symmetric){
    triad.type <- "ANY" # this does not matter, as there are no structural differences in symmetric networks.
  }

  ntypes <- length(triad.type)
  results <- matrix(0,
                    nrow = nrow(net),
                    ncol = length(ntypes))
  colnames(results) <- triad.type
  for(i in 1:ntypes){
    results[,i] <-  pibcpp(net = net,
                           gmem = g.mem,
                           symmetric = symmetric,
                           traidtype = triad.type[i],
                           count=count,
                           isolate=isolate)
  }
  if(is.null(rownames(net))){
    rownames(results) <- paste0("Node_", 1:nrow(net))
  }else{
    rownames(results) <- rownames(net)

  }
  return(results)

}
