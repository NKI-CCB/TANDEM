# TANDEM 1.0.2
Fixed a bug in relative.contributions() that decreased the relative contributions for correlated data types (by not taking into account the correlated part). Thanks to Wouter Touw for pointing this out!

# TANDEM 1.0.1
Fixed a bug that occurred when using lambda.upstream="lambda.min" or lambda.downstream="lambda.min" in tandem(). All regression coefficients were correct, but the intercept was estimated using lambda.1se. This also affected results from predict() or nested.cv() using TANDEM fits where lambda.upstream="lambda.min" or lambda.downstream="lambda.min", negatively affecting the predictive performance.

In addition, I've added an option to set lambda="lambda.min" for the glmnet model in relative.contributions().


# TANDEM 1.0.0
First release version!
