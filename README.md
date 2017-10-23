This project contains all code to replicate simulations in
`Hypergraph Modeling of the Evolution of Team Relationships with Applications to Academic Funding Collaboration`
by Ryan Haunfelder, Bailey Fosdick, and Haonan Wang.

Given holding parameters for the eight subgraphs, transition
probabilities between them, and a viewing length, the `triad_sim`
function generates a single observation from the CTMC model. The dyads
and triads column give each state of the simulated chain. Each chain
starts at the initial state
*h*<sub>00</sub>
, so the first row will always have a dyad and triad value of 0.

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

The simulations are

<!-- html table generated in R 3.4.1 by xtable 1.8-2 package -->
<!-- Sun Oct 22 19:06:45 2017 -->
<table border="1">
<tr>
<th>
</th>
<th>
\_{00}
</th>
<th>
\_{10}
</th>
<th>
\_{20}
</th>
<th>
\_{30}
</th>
<th>
\_{01}
</th>
<th>
\_{11}
</th>
<th>
\_{21}
</th>
</tr>
<tr>
<td align="right">
Simulation 1
</td>
<td align="right">
250.00
</td>
<td align="right">
200.00
</td>
<td align="right">
200.00
</td>
<td align="right">
300.00
</td>
<td align="right">
300.00
</td>
<td align="right">
150.00
</td>
<td align="right">
150.00
</td>
</tr>
<tr>
<td align="right">
Simulation 2
</td>
<td align="right">
500.00
</td>
<td align="right">
100.00
</td>
<td align="right">
100.00
</td>
<td align="right">
150.00
</td>
<td align="right">
150.00
</td>
<td align="right">
200.00
</td>
<td align="right">
200.00
</td>
</tr>
<tr>
<td align="right">
Simulation 3
</td>
<td align="right">
200.00
</td>
<td align="right">
100.00
</td>
<td align="right">
100.00
</td>
<td align="right">
150.00
</td>
<td align="right">
400.00
</td>
<td align="right">
450.00
</td>
<td align="right">
450.00
</td>
</tr>
<tr>
<td align="right">
Simulation 4
</td>
<td align="right">
200.00
</td>
<td align="right">
400.00
</td>
<td align="right">
400.00
</td>
<td align="right">
450.00
</td>
<td align="right">
100.00
</td>
<td align="right">
150.00
</td>
<td align="right">
150.00
</td>
</tr>
</table>
<!-- html table generated in R 3.4.1 by xtable 1.8-2 package -->
<!-- Sun Oct 22 19:06:46 2017 -->
<table border="1">
<tr>
<th>
</th>
<th>
p\_{01|00}
</th>
<th>
p\_{10|00}
</th>
<th>
p\_{11|10}
</th>
<th>
p\_{20|10}
</th>
<th>
p\_{11|01}
</th>
<th>
p\_{21|11}
</th>
<th>
p\_{31|21}
</th>
<th>
p\_{31|30}
</th>
<th>
p\_{21|20}
</th>
<th>
p\_{30|20}
</th>
</tr>
<tr>
<td align="right">
Simulation 1
</td>
<td align="right">
0.45
</td>
<td align="right">
0.45
</td>
<td align="right">
0.20
</td>
<td align="right">
0.45
</td>
<td align="right">
0.25
</td>
<td align="right">
0.30
</td>
<td align="right">
0.35
</td>
<td align="right">
0.35
</td>
<td align="right">
0.20
</td>
<td align="right">
0.30
</td>
</tr>
<tr>
<td align="right">
Simulation 2
</td>
<td align="right">
0.10
</td>
<td align="right">
0.30
</td>
<td align="right">
0.35
</td>
<td align="right">
0.60
</td>
<td align="right">
0.90
</td>
<td align="right">
0.90
</td>
<td align="right">
0.80
</td>
<td align="right">
0.80
</td>
<td align="right">
0.25
</td>
<td align="right">
0.65
</td>
</tr>
<tr>
<td align="right">
Simulation 3
</td>
<td align="right">
0.10
</td>
<td align="right">
0.40
</td>
<td align="right">
0.05
</td>
<td align="right">
0.70
</td>
<td align="right">
0.50
</td>
<td align="right">
0.50
</td>
<td align="right">
0.50
</td>
<td align="right">
0.80
</td>
<td align="right">
0.05
</td>
<td align="right">
0.75
</td>
</tr>
<tr>
<td align="right">
Simulation 4
</td>
<td align="right">
0.60
</td>
<td align="right">
0.20
</td>
<td align="right">
0.60
</td>
<td align="right">
0.10
</td>
<td align="right">
0.70
</td>
<td align="right">
0.70
</td>
<td align="right">
0.70
</td>
<td align="right">
0.60
</td>
<td align="right">
0.50
</td>
<td align="right">
0.40
</td>
</tr>
</table>
