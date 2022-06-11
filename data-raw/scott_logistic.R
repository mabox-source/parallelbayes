
# Create the small data set of binary outcomes and predictor variables used for 
# logistic regression in Scott et al 2016.

n <- 2755 + 2753 + 1186 + 717 + 1173 + 305 + 301 + 706 + 32 + 17 + 24 + 10 + 2 + 13 + 2 + 4

scott_logistic <- data.frame(
  y = c(rep(0, 2755 - 266), rep(1, 266), rep(0, 2753 - 116), rep(1, 116), rep(0, 1186 - 34), rep(1, 34), rep(0, 717 - 190), rep(1, 190), rep(0, 1173 - 61), rep(1, 61), rep(0, 305 - 37), rep(1, 37), rep(0, 301 - 68), rep(1, 68), rep(0, 706 - 119), rep(1, 119), rep(0, 32 - 18), rep(1, 18), rep(0, 17 - 13), rep(1, 13), rep(0, 24 - 18), rep(1, 18), rep(0, 10 - 8), rep(1, 8), rep(0, 2 - 2), rep(1, 2), rep(0, 13 - 7), rep(1, 7), rep(0, 2 - 2), rep(1, 2), rep(0, 4 - 3), rep(1, 3)),
  x1 = rep(1, n),
  x2 = c(rep(0, 2755), rep(0, 2753), rep(0, 1186), rep(1, 717), rep(0, 1173), rep(1, 305), rep(1, 301), rep(1, 706), rep(0, 32), rep(0, 17), rep(0, 24), rep(1, 10), rep(1, 2), rep(0, 13), rep(1, 2), rep(1, 4)),
  x3 = c(rep(0, 2755), rep(0, 2753), rep(1, 1186), rep(0, 717), rep(1, 1173), rep(1, 305), rep(1, 301), rep(0, 706), rep(0, 32), rep(1, 17), rep(0, 24), rep(0, 10), rep(1, 2), rep(1, 13), rep(1, 2), rep(0, 4)),
  x4 = c(rep(1, 2755), rep(0, 2753), rep(0, 1186), rep(1, 717), rep(1, 1173), rep(0, 305), rep(1, 301), rep(0, 706), rep(0, 32), rep(1, 17), rep(1, 24), rep(1, 10), rep(0, 2), rep(0, 13), rep(1, 2), rep(0, 4)),
  x5 = c(rep(0, 2755), rep(0, 2753), rep(0, 1186), rep(0, 717), rep(0, 1173), rep(0, 305), rep(0, 301), rep(0, 706), rep(1, 32), rep(1, 17), rep(1, 24), rep(1, 10), rep(1, 2), rep(1, 13), rep(1, 2), rep(1, 4))
)

save(scott_logistic, file = "../data/scott_logistic.RData")
