This project contains all code to replicate simulations in
`Hypergraph Modeling of the Evolution of Team Relationships with Applications to Academic Funding Collaboration`
by Ryan Haunfelder, Bailey Fosdick, and Haonan Wang.

Given holding parameters for the eight subgraphs, transition
probabilities between them, and a viewing length, the `triad_sim`
function generates a single observation from the CTMC model. The dyads
and triads column give each state of the simulated chain. Each chain
starts at the initial state *h*<sub>00</sub>, so the first row will
always have a dyad and triad value of 0.

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

The four simulation parameters are summarized in the table below.

<table>
<colgroup>
<col width="13%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
<col width="12%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"> </th>
<th align="center"><span class="math inline"><em>λ</em><sub>00</sub></span></th>
<th align="center"><span class="math inline"><em>λ</em><sub>10</sub></span></th>
<th align="center"><span class="math inline"><em>λ</em><sub>20</sub></span></th>
<th align="center"><span class="math inline"><em>λ</em><sub>30</sub></span></th>
<th align="center"><span class="math inline"><em>λ</em><sub>01</sub></span></th>
<th align="center"><span class="math inline"><em>λ</em><sub>11</sub></span></th>
<th align="center"><span class="math inline"><em>λ</em><sub>21</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong>Simulation 1</strong></td>
<td align="center">250</td>
<td align="center">200</td>
<td align="center">200</td>
<td align="center">300</td>
<td align="center">300</td>
<td align="center">150</td>
<td align="center">150</td>
</tr>
<tr class="even">
<td align="center"><strong>Simulation 2</strong></td>
<td align="center">500</td>
<td align="center">100</td>
<td align="center">100</td>
<td align="center">150</td>
<td align="center">150</td>
<td align="center">200</td>
<td align="center">200</td>
</tr>
<tr class="odd">
<td align="center"><strong>Simulation 3</strong></td>
<td align="center">200</td>
<td align="center">100</td>
<td align="center">100</td>
<td align="center">150</td>
<td align="center">400</td>
<td align="center">450</td>
<td align="center">450</td>
</tr>
<tr class="even">
<td align="center"><strong>Simulation 4</strong></td>
<td align="center">200</td>
<td align="center">400</td>
<td align="center">400</td>
<td align="center">450</td>
<td align="center">100</td>
<td align="center">150</td>
<td align="center">150</td>
</tr>
</tbody>
</table>

<table style="width:100%;">
<colgroup>
<col width="11%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
<col width="8%" />
</colgroup>
<thead>
<tr class="header">
<th align="center"> </th>
<th align="center"><span class="math inline"><em>p</em><sub>01|00</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>10|00</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>11|10</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>20|10</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>11|01</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>21|11</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>31|21</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>31|30</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>21|20</sub></span></th>
<th align="center"><span class="math inline"><em>p</em><sub>30|20</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="center"><strong>Simulation 1</strong></td>
<td align="center">0.45</td>
<td align="center">0.45</td>
<td align="center">0.2</td>
<td align="center">0.45</td>
<td align="center">0.25</td>
<td align="center">0.3</td>
<td align="center">0.35</td>
<td align="center">0.35</td>
<td align="center">0.2</td>
<td align="center">0.3</td>
</tr>
<tr class="even">
<td align="center"><strong>Simulation 2</strong></td>
<td align="center">0.1</td>
<td align="center">0.3</td>
<td align="center">0.35</td>
<td align="center">0.6</td>
<td align="center">0.9</td>
<td align="center">0.9</td>
<td align="center">0.8</td>
<td align="center">0.8</td>
<td align="center">0.25</td>
<td align="center">0.65</td>
</tr>
<tr class="odd">
<td align="center"><strong>Simulation 3</strong></td>
<td align="center">0.1</td>
<td align="center">0.4</td>
<td align="center">0.05</td>
<td align="center">0.7</td>
<td align="center">0.5</td>
<td align="center">0.5</td>
<td align="center">0.5</td>
<td align="center">0.8</td>
<td align="center">0.05</td>
<td align="center">0.75</td>
</tr>
<tr class="even">
<td align="center"><strong>Simulation 4</strong></td>
<td align="center">0.6</td>
<td align="center">0.2</td>
<td align="center">0.6</td>
<td align="center">0.1</td>
<td align="center">0.7</td>
<td align="center">0.7</td>
<td align="center">0.7</td>
<td align="center">0.6</td>
<td align="center">0.5</td>
<td align="center">0.4</td>
</tr>
</tbody>
</table>
