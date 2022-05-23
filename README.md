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
* The fractionated prior concept: this originates in the published literature 
for the consensus and DPE algorithms. It is a theoretical requirement but 
perhaps not a practical one. If you don't use the fractionated prior, your 
samples will not be from the posterior distribution implied by your original 
choice of prior distribution. However, it is not always possible to perform the 
fractionation of the prior distribution, and if the prior is not very 
informative it possibly does not matter very much.
* `R` package `BayesLogit` has been reduced to containing only the `rpg` 
function. I will need to add my own Gibbs sampler.

## Spark notes

* I have JRE version 11 installed. When first launching a Spark cluster with 
`sparklyr` I got: "Java 11 is only supported for Spark 3.0.0+" (I had Spark 
2.4.3 installed, by default). So installed Spark 3.1.2 using 
`sparklyr::spark_install` (tried 3.2.x but got a bad URL error).
* Spark Web UI accessed via http://localhost:4040
* Currently in vignettes I am using `sdf_copy_to`, but `sdf_repartition` might 
be better (requires Spark v2 or greater).
* Launching a Spark cluster with exactly `M` executors: setting:

```
config$spark.dynamicAllocation.enabled <- "false"
config$spark.executor.instances <- M
```

is not sufficient. The configs required in local mode are different from those 
on the cluster.
