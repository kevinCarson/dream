#######
# The Creation of Interpersonal Cultural Brokerage Statistics Based on
# Reference: Leal, Diego F. Forthcoming. "Cultural Brokerage and Diffusion Across Bright Symbolic
# Boundaries. A New Measure to Locate Brokers Spanning Cultural Holes"
# Reference: Gould, Roger V., and Roberto M. Fernandez. 1989. "Structures of Mediation: A Formal Approach to Brokerage in Transaction Networks"
#            Sociological Methodology, Vol. 19, pp. 89-126
## Code written by Kevin Carson (kacarson@arizona.edu) and Deigo Leal (https://www.diegoleal.info/)
## Last Updated: 09-21-24

#' @title Compute Potential for Intercultural Brokerage (PIB) Based on Leal (2025)
#' @name computeLealBrokerage
#' @param net The one-mode adjacency matrix.
#' @param g.mem The vector of membership values that the brokerage scores will be based on.
#' @param symmetric TRUE/FALSE. TRUE indicates that network matrix will be treated as symmetric. FALSE indicates that the network matrix will be treated as asymmetric. Set to TRUE by default.
#' @param triad.type The string value (or vector) that indicates what specific triadic (star) structures the potential for cultural brokerage will be computed for. Possible values are "ANY", "OTS", "ITS", "MTS" (see the details section). The function defaults to “ANY”.
#' @param count TRUE/FALSE. TRUE indicates that the number of culturally brokered open triangles will be returned. FALSE indicates that the proportion of culturally brokered open triangles to all open triangles will be returned (see the details section). Set to TRUE by default.
#' @param isolate If count = FALSE, the numerical value that will be given to isolates. This value is set to NA by default, as 0/0 is undefined. The user can specify this value!
#' @return The vector of interpersonal cultural brokerage values for the one-mode network.
#' @export
#'
#'
#'
#'
#'
#' @description
#' `r lifecycle::badge("deprecated")`
#'
#' `computeLealBrokerage()` has been deprecated starting on version 1.0.0 of the `dream` package. Please use the `netstats_om_pib()` function and see the `NEWS.md` file for more details.
#'
#'
#' Following Leal (2025), this function calculates node’s Potential
#' for Intercultural Brokerage (PIB) in a one-mode network, that is, brokerage
#' based on nodes’ distinct group memberships. For example, users can examine PIB
#' based on actors’ gender. The option count determines what is returned by
#' the function. If count is TRUE, then the count of ‘culturally’ dissimilar pairs
#' brokered by ego is included (i.e., ego’s total count of brokered open triangles
#' where the alters at the two endpoints of said open triangles are ‘culturally’
#' dissimilar from one another). If count is FALSE, the proportion of ego’s brokered
#' open triangles where the endpoints are ‘culturally’ dissimilar out of all of ego’s
#' brokered open triangles (regardless of the cultural identity of the alters) is
#' returned. The formula for computing interpersonal brokerage is presented in
#' the details section.
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
#' values for incoming two-stars (i->j<-k or the 021U triad according to the M.A.N. notation; see Wasserman and Faust 1994 ), where j is the broker. "MTS" computes PIB for
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
#' computeLealBrokerage(LealNet,
#'    symmetric = TRUE,
#'    g.mem = categorical_variable)
#'
#'
#'
computeLealBrokerage <- function(net,
                 g.mem,
                 symmetric = TRUE,
                 triad.type = NULL, #can be of c("ANY", "OTS", "ITS", "MTS")
                 count = TRUE,
                 isolate = NA){

  lifecycle::deprecate_warn("1.0.0", "computeLealBrokerage()", "netstats_om_pib()")

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

  if(count == TRUE){

    ## all potential triads (i , j  , k)  ; j -> ego , (i,k) alters, and treating the ties as directed (i.e., i -> k or k -> i) is a tie
    if(symmetric == FALSE){

      starsT <- length(triad.type)
      results <- matrix(0,
                        nrow = nrow(net),
                        ncol = length(triad.type))
      colnames(results) <- triad.type
      for(type in 1:starsT){

        if(triad.type[type] == "OTS"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          structural.scores <- rep(0, nrow(net))
          #for all alters i
          #for all alters i
          for(i in 1:nrow(net)){

            #for all alters k
            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[j,k] == 1   ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[j,i] == 1  ) ){

                    # now let's check for redundancy within the network
                    red.actor <- 0


                    for(z in 1:nrow(net)){ #for z in all other alters!

                      if( ( (z != k) & (z != i)  & (z != j) )  & #z is not i, j, or k (members of the original triad)

                          #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                          ( net[z,k] == 1 ) &

                          #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                          ( net[z,i] == 1 )){


                        red.actor <- red.actor + 1 #if this alter also is a cultural broker, increase the score by 1
                      }
                    }
                    if(red.actor == 0){ #if there are no other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + 1 #give the jth actor a full 1 broker score
                    }
                    if(red.actor != 0){#if there are other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + (1 /( red.actor + 1) ) #give the jth actor a score of  1 / the number of actors being a broker

                    }
                  }

                }}} }
          results[,which(colnames(results) == "OTS")] <- structural.scores
        }
        if(triad.type[type] == "ANY"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          structural.scores <- rep(0, nrow(net))
          #for all alters i
          for(i in 1:nrow(net)){

            #for all alters k
            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[j,k] == 1 | net[k,j] == 1    ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[j,i] == 1 | net[i,j] == 1    )){

                    # now let's check for redundancy within the network
                    red.actor <- 0


                    for(z in 1:nrow(net)){ #for z in all other alters!

                      if( ( (z != k) & (z != i)  & (z != j) )  & #z is not i, j, or k (members of the original triad)

                          #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                          ( net[z,k] == 1 | net[k,z] == 1    ) &

                          #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                          ( net[z,i] == 1 | net[i,z] == 1    )){


                        red.actor <- red.actor + 1 #if this alter also is a cultural broker, increase the score by 1
                      }



                    }

                    if(red.actor == 0){ #if there are no other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + 1 #give the jth actor a full 1 broker score

                    }

                    if(red.actor != 0){#if there are other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + (1 /( red.actor + 1) ) #give the jth actor a score of  1 / the number of actors being a broker

                    }
                  }

                }}} }
          results[,which(colnames(results) == "ANY")] <- structural.scores
        }
        if(triad.type[type] == "ITS"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          structural.scores <- rep(0, nrow(net))
          #for all alters i
          #for all alters i
          for(i in 1:nrow(net)){

            #for all alters k
            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[k,j] == 1   ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[i,j] == 1  ) ){

                    # now let's check for redundancy within the network
                    red.actor <- 0


                    for(z in 1:nrow(net)){ #for z in all other alters!

                      if( ( (z != k) & (z != i)  & (z != j) )  & #z is not i, j, or k (members of the original triad)

                          #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                          ( net[k,z] == 1 ) &

                          #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                          ( net[i,z] == 1 )){


                        red.actor <- red.actor + 1 #if this alter also is a cultural broker, increase the score by 1
                      }
                    }
                    if(red.actor == 0){ #if there are no other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + 1 #give the jth actor a full 1 broker score
                    }
                    if(red.actor != 0){#if there are other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + (1 /( red.actor + 1) ) #give the jth actor a score of  1 / the number of actors being a broker

                    }
                  }

                }}} }
          results[,which(colnames(results) == "ITS")] <- structural.scores
        }
        if(triad.type[type] == "MTS"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          structural.scores <- rep(0, nrow(net))
          #for all alters i
          for(i in 1:nrow(net)){

            #for all alters k
            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( (( (j != k) & (j != i)  )  &

                       #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                       ( net[k,j] == 1 & net[j,i] == 1    )) |
                      (( (j != k) & (j != i)  )  &

                       #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                       (net[j,k] == 1 & net[i,j] == 1   )) ){

                    # now let's check for redundancy within the network
                    red.actor <- 0


                    for(z in 1:nrow(net)){ #for z in all other alters!

                      if( (( (z != k) & (z != i)  & (z != j) )  & #z is not i, j, or k (members of the original triad)

                           #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                           ( net[k,z] == 1 & net[z,i] == 1    )) |
                          (( (z != k) & (z != i)  & (z != j) )  &
                           #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                           (net[z,k] == 1 & net[i,z] == 1   ))){


                        red.actor <- red.actor + 1 #if this alter also is a cultural broker, increase the score by 1
                      }
                    }
                    if(red.actor == 0){ #if there are no other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + 1 #give the jth actor a full 1 broker score
                    }
                    if(red.actor != 0){#if there are other actors brokering the tie
                      structural.scores[j] <- structural.scores[j] + (1 /( red.actor + 1) ) #give the jth actor a score of  1 / the number of actors being a broker

                    }
                  }

                }}} }
          results[,which(colnames(results) == "MTS")] <- structural.scores
        }

      }

    }
    ## all potential triads (i , j  , k)  ; j -> ego , (i,k) alters, and treating the ties as undirected (i.e., i -> j and j -> i; SYMMETRIC) is a tie
    if(symmetric == TRUE ){
      #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
      structural.scores <- rep(0, nrow(net))
      #for all alters i
      for(i in 1:nrow(net)){

        #for all alters k
        for(k in 1:nrow(net)){

          #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
          if( ( (i != k) & (i > k) )  &

              #if j -> k, or k -> j,        i and k must not be  friends
              ( net[i,k] == 0 |  net[k,i] == 0    ) &

              #if k and i do not share the same group membership, ie the endpoints are
              #different boundary points (for instance k is male and i is female)

              #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
              (g.mem[k] != g.mem[i])){

            for(j in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (j != k) & (j != i)  )  &

                  #if j -> k, or k -> j,        j and k must be symmetric friends
                  ( net[j,k] == 1 & net[k,j] == 1    ) &

                  #if j -> i, or i -> j,        j and i must be symmetric friends
                  ( net[j,i] == 1 & net[i,j] == 1    )){

                # now let's check for redundancy within the network
                red.actor <- 0


                for(z in 1:nrow(net)){ #for z in all other alters!

                  if( ( (z != k) & (z != i)  & (z != j) )  & #z is not i, j, or k (members of the original triad)

                      #if j -> k, or k -> j,        z and k must be symmetric friends
                      ( net[z,k] == 1 & net[k,z] == 1    ) &

                      #if j -> i, or i -> j,        z and i must be symmetric friends
                      ( net[z,i] == 1 & net[i,z] == 1    )){


                    red.actor <- red.actor + 1 #if this alter also is a cultural broker, increase the score by 1
                  }



                }

                if(red.actor == 0){ #if there are no other actors brokering the tie
                  structural.scores[j] <- structural.scores[j] + 1 #give the jth actor a full 1 broker score

                }

                if(red.actor != 0){#if there are other actors brokering the tie
                  structural.scores[j] <- structural.scores[j] + (1 /( red.actor + 1) ) #give the jth actor a score of  1 / the number of actors being a broker

                }
              }

            }}} }
    }

    if(symmetric == TRUE ){
      if(is.null(rownames(net))){
        names(structural.scores)  <- NULL
      }else{
        names(structural.scores) <- rownames(net)
      }
      return(structural.scores)
    }

    if(symmetric == FALSE ){
      if(is.null(rownames(net))){
        rownames(results) <- 1:nrow(net)
      }else{
        rownames(results) <- rownames(net)
      }
      return(results)
    }


  }


  if(count == FALSE){

    if(symmetric == FALSE){
      starsT <- length(triad.type)
      results <- matrix(0,
                        nrow = nrow(net),
                        ncol = length(triad.type))
      colnames(results) <- triad.type
      for(type in 1:starsT){

        if(triad.type[type] == "OTS"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          cul.scores <- rep(0, nrow(net))
          struct.scores <- rep(0, nrow(net))

          #for all alters i
          #for all alters i
          for(i in 1:nrow(net)){

            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[j,k] == 1   ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[j,i] == 1  ) ){

                    # now let's check for redundancy within the network
                    cul.scores[j] <- cul.scores[j] + 1

                  }

                }}

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) ){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[j,k] == 1   ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[j,i] == 1  ) ){

                    # now let's check for redundancy within the network
                    struct.scores[j] <- struct.scores[j] + 1

                  }

                }}






            }
          }

          structural.scores <- cul.scores / struct.scores
          structural.scores[is.nan(structural.scores)] <-   isolate
          results[,which(colnames(results) == "OTS")] <- structural.scores
        }




        if(triad.type[type] == "ANY"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          cul.scores <- rep(0, nrow(net))
          struct.scores <- rep(0, nrow(net))

          for(i in 1:nrow(net)){

            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[j,k] == 1 | net[k,j] == 1    ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[j,i] == 1 | net[i,j] == 1    ) ){

                    # now let's check for redundancy within the network
                    cul.scores[j] <- cul.scores[j] + 1

                  }

                }}

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) ){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[j,k] == 1 | net[k,j] == 1    ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[j,i] == 1 | net[i,j] == 1    ) ){

                    # now let's check for redundancy within the network
                    struct.scores[j] <- struct.scores[j] + 1

                  }

                }}






            }
          }
          structural.scores <- cul.scores / struct.scores
          structural.scores[is.nan(structural.scores)] <-   isolate
          results[,which(colnames(results) == "ANY")] <- structural.scores
        }






        if(triad.type[type] == "ITS"){
          cul.scores <- rep(0, nrow(net))
          struct.scores <- rep(0, nrow(net))

          for(i in 1:nrow(net)){

            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[k,j] == 1   ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[i,j] == 1  ) ){

                    # now let's check for redundancy within the network
                    cul.scores[j] <- cul.scores[j] + 1

                  }

                }}

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) ){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( ( (j != k) & (j != i)  )  &

                      #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                      ( net[k,j] == 1   ) &

                      #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                      ( net[i,j] == 1  ) ){

                    # now let's check for redundancy within the network
                    struct.scores[j] <- struct.scores[j] + 1

                  }

                }}






            }
          }
          structural.scores <- cul.scores / struct.scores
          structural.scores[is.nan(structural.scores)] <-   isolate
          results[,which(colnames(results) == "ITS")] <- structural.scores
        }
        if(triad.type[type] == "MTS"){
          #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
          cul.scores <- rep(0, nrow(net))
          struct.scores <- rep(0, nrow(net))

          for(i in 1:nrow(net)){

            for(k in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) &

                  #if k and i do not share the same group membership, ie the endpoints are
                  #different boundary points (for instance k is male and i is female)

                  #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
                  (g.mem[k] != g.mem[i])){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( (( (j != k) & (j != i)  )  &

                       #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                       ( net[k,j] == 1 & net[j,i] == 1    )) |
                      (( (j != k) & (j != i)  )  &
                       #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                       (net[j,k] == 1 & net[i,j] == 1   )) ){

                    # now let's check for redundancy within the network
                    cul.scores[j] <- cul.scores[j] + 1

                  }

                }}

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (i != k) & (i > k) )  &

                  #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                  ( net[i,k] == 0 &  net[k,i] == 0    ) ){

                for(j in 1:nrow(net)){

                  #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
                  if( (( (j != k) & (j != i)  )  &

                       #if j -> k, or k -> j,        j and k must be at least asymmetric friends
                       ( net[k,j] == 1 & net[j,i] == 1    )) |
                      (( (j != k) & (j != i)  )  &
                       #if j -> i, or i -> j,        j and i must be at least asymmetric friends
                       (net[j,k] == 1 & net[i,j] == 1   ) )){

                    # now let's check for redundancy within the network
                    struct.scores[j] <- struct.scores[j] + 1

                  }

                }}






            }
          }
          structural.scores <- cul.scores / struct.scores
          structural.scores[is.nan(structural.scores)] <-   isolate
          results[,which(colnames(results) == "MTS")] <- structural.scores
        }

      }
    }

    ## all potential triads (i , j  , k)  ; j -> ego , (i,k) alters, and treating the ties as undirected (i.e., i -> j and j -> i; SYMMETRIC) is a tie
    if(symmetric == TRUE){
      #creating a vector to hold the different brokerage scores of each ego (STRUCTURAL REDUNDANCY)
      cul.scores <- rep(0, nrow(net))
      struct.scores <- rep(0, nrow(net))
      #for all alters i
      for(i in 1:nrow(net)){

        #for all alters k
        for(k in 1:nrow(net)){

          #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
          if( ( (i != k) & (i > k) )  &

              #if j -> k, or k -> j,        i and k must not be  friends
              ( net[i,k] == 0 |  net[k,i] == 0    ) &

              #if k and i do not share the same group membership, ie the endpoints are
              #different boundary points (for instance k is male and i is female)

              #DO NOTE: this function always for brokering for the liaison, representative, and gatekeeper broker roles
              (g.mem[k] != g.mem[i])){

            for(j in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (j != k) & (j != i)  )  &

                  #if j -> k, or k -> j,        j and k must be symmetric friends
                  ( net[j,k] == 1 & net[k,j] == 1    ) &

                  #if j -> i, or i -> j,        j and i must be symmetric friends
                  ( net[j,i] == 1 & net[i,j] == 1    )){

                # now let's check for redundancy within the network
                cul.scores[j] <- cul.scores[j] + 1


              }

            }
          }



          #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
          if( ( (i != k) & (i > k) )  &

              #if j -> k, or k -> j,        i and k must not be  friends
              ( net[i,k] == 0 |  net[k,i] == 0    ) ){

            for(j in 1:nrow(net)){

              #if i, j, and k are different alters, and going through the upper triangle of the friendship matrix
              if( ( (j != k) & (j != i)  )  &

                  #if j -> k, or k -> j,        j and k must be symmetric friends
                  ( net[j,k] == 1 & net[k,j] == 1    ) &

                  #if j -> i, or i -> j,        j and i must be symmetric friends
                  ( net[j,i] == 1 & net[i,j] == 1    )){

                # now let's check for redundancy within the network
                struct.scores[j] <- struct.scores[j] + 1


              }

            }
          }
        }
      }
      structural.scores <- cul.scores / struct.scores
      structural.scores[is.nan(structural.scores)] <-   isolate
    }

    if(symmetric == TRUE ){
      ifelse(is.null(rownames(net)),
             names(structural.scores)  <- NULL,
             names(structural.scores) <- rownames(net))
      return(structural.scores)
    }

    if(symmetric == FALSE ){
      ifelse(is.null(rownames(net)),
             rownames(results) <- 1:nrow(net),
             rownames(results) <- rownames(net))
      return(results)
    }

  }
}
