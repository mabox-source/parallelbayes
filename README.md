# parallelbayes

This package implements the Laplace enriched multiple importance estimator of 
[[1]](#1). This is an importance weighting scheme for Monte Carlo estimation of 
posterior expectations using samples drawn from multiple proposal 
distributions. The idea is that this be used to reconstruct posterior inference 
when the data are partitioned and a sampling algorithm, such as MCMC, was used 
to target the local posteriors conditioned on parts of data.

Also implemented are the consensus Monte Carlo algorithm of [[2]](#2) and the 
NDPE and SDPE algorithms of [[3]](#3),


## Installation

The latest version can be installed from GitHub by cloning the repository and 
using the `devtools` package, or by using the `remotes` package:

```
install.packages("remotes")
remotes::install_github("mabox-source/parallelbayes")
```


## Notes

* `lemie.weights` will return type 1 weights AND type 2 weights, if requested 
using `type` and `keep.type1`. This is because it's quite computationally 
demanding, but the steps in the computation of the type 1 weights can be reused 
in the computation of the type 2 weights and we often want both for comparison.
* PSIS was originally applied in `lemie.weights`, but again we wanted to 
compute all weights: type 1, type 2, type 1 with PSIS and type 2 with PSIS, and 
there is lots of shared computation in these. So now the PSIS transformation of 
the weights is applied using function `pareto_smooth`, called AFTER 
`lemie.weights`.
* The fractionated prior concept: this originates in the published literature 
for the consensus and DPE algorithms. It is a theoretical requirement but 
perhaps not a practical one. If you don't use the fractionated prior, your 
samples will not be from the posterior distribution implied by your original 
choice of prior distribution. However, it is not always possible to perform the 
fractionation of the prior distribution, and if the prior is not very 
informative it possibly does not matter very much.
* `R` package `BayesLogit` has been reduced to containing only the `rpg` 
function. I have added my own Gibbs sampler, which makes use of `rpg`.
* In the consensus and DPE algorithms, samples are pooled so that you end up 
with only $\frac{H^{total}}{M}$ samples. The LEMIE (and naive pooling) methods 
do not: you end up with $H^{total}$ samples. This is reflected in greater ESS 
(see the vignettes).
* I made type 1 weights an argument of the `lemie.weights` function, so that we 
can compute the type 2 weights after computing the type 1 weights only. This 
allows us to get timings for both algorithms without wasting much time (much 
computation is shared between the type 1 and type 2 weights). It will also 
speed up the type 2 weight computation if `loglik` is also supplied as an 
argument. This is why both `loglik` and `w.type_1` can be returned by 
`lemie.weights`.

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


## References

<a id="1">[1]</a> Box, M. (2022). Importance Sampling Methods for Bayesian Inference with Partitioned Data. *arXiv preprint arXiv:2210.06620*.<br/>
<a id="2">[2]</a> Scott, S. L., Blocker, A. W., Bonassi, F. V., Chipman, H. A., George, E. I., and McCulloch, R. E. (2016). Bayes and big data: The consensus monte carlo algorithm. *International Journal of Management Science and Engineering Management*, 11(2):78–88.<br/>
<a id="3">[3]</a> Neiswanger, W., Wang, C., and Xing, E. (2013). Asymptotically exact, embarrassingly parallel mcmc. *arXiv preprint arXiv:1311.4780*.<br/>
