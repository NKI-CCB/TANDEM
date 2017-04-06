# TANDEM 1.0.1
Fixed a bug that occurred when using lambda.upstream="lambda.min" or lambda.downstream="lambda.min" in tandem(). All regression coefficients were correct, but the intercept was estimated using lambda.1se. This also affected results from predict() or nested.cv() using TANDEM fits where lambda.upstream="lambda.min" or lambda.downstream="lambda.min", negatively affecting the predictive performance.

In addition, I've added an option to set lambda="lambda.min" for the glmnet model in relative.contributions().


# TANDEM 1.0.0
First release version!
