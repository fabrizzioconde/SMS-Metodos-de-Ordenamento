# empirical CAR scores for diabetes data
# load care library
library("care")
data(efron2004)
xnames = colnames(efron2004$x)
n = dim(efron2004$x)[1]
car1 = carscore(efron2004$x, efron2004$y, lambda=0)
car1
car2 = carscore(efron2004$x, efron2004$y, lambda=0, diagonal = TRUE) # Correlação Margianal
car2
car3 = carscore(efron2004$x, efron2004$y, diagonal = FALSE) # CARSCORE
car3
# compare orderings
# variables ordered by squared CAR scores
xnames[order(car1^2, decreasing=TRUE)]
xnames[order(car2^2, decreasing=TRUE)]
# "bmi" "s5" "bp" "s3" "s4" "s6" "sex" "age" "s2" "s1"
# compare with ordering by t-scores / partial correlations
pcor = pcor.shrink(cbind(efron2004$y,efron2004$x), lambda=0, verbose=FALSE)[-1,1]
xnames[order(pcor^2, decreasing=TRUE)]
# "bmi" "bp" "s5" "sex" "s1" "s2" "s4" "s6" "s3" "age"
# compare with ordering by marginal correlations
#mcor = cor(efron2004$y,efron2004$x)
mcor = carscore(efron2004$x, efron2004$y, diagonal=TRUE, lambda=0)
xnames[order(mcor^2, decreasing=TRUE)]
