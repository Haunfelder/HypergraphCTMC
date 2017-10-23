This project contains all code to replicate simulations in
`Hypergraph Modeling of the Evolution of Team Relationships with Applications to Academic Funding Collaboration`
by Ryan Haunfelder, Bailey Fosdick, and Haonan Wang.

The hypergraph data object is stored as a data frame with three columns.

    lambda_true = c(1/250, 1/200, 1/200, 1/300, 1/300, 1/150 , 1/150)
    names(lambda_true) = c("lambda00","lambda10","lambda20","lambda30","lambda01","lambda11","lambda21")

    trans_probs_true = c(0.45, 0.45, 0.2, 0.45, 0.25, 0.3, 0.35, 0.35, 0.2, 0.3)
    names(trans_probs_true) = c("p01|00","p10|00","p11|10","p20|10","p11|01",
                              "p21|11","p31|21","p31|30",
                              "p21|20","p30|20")
    pars = c(lambda_true,trans_probs_true)

    triad_sim(pars, 2000)

    ##   dyads triads TimeBetween  cumtime
    ## 1     0      0    489.7654 489.7654
    ## 2     0      1   1510.2346      Inf
