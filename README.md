# remix

## Notes

* `calc_weights` will return type 1 weights AND type 2 weights, if requested 
using `type` and `keep.type1`. This is because it's quite computationally 
demanding, but the steps in the computation of the type 1 weights can be reused 
in the computation of the type 2 weights and we often want both for comparison.
* PSIS was originally applied in `calc_weights`, but again we wanted to compute 
all weights: type 1, type 2, type 1 with PSIS and type 2 with PSIS, and there 
is lots of shared computation in these. So now the PSIS transformation of the 
weights is applied using function `pareto_smooth`, called AFTER `calc_weights`.
