#' Fits a TANDEM model by performing a two-stage regression
#'
#' Fits a TANDEM model by performing a two-stage regression. In the first stage, all upstream features (x[,upstream]) are regressed
#' on the output y. In the second stage, the downstream features (x[,!upstream]) are regressed on the residuals of the first stage.
#' In both stages Elastic Net regression (as implemented in cv.glmnet() from the glmnet package) is used to perform the regression.
#'
#' @param x A feature matrix, where the rows correspond to samples and the columns to features.
#' @param y A vector containing the response.
#' @param upstream A boolean vector that indicates for each feature whether it's upstream (TRUE) or downstream (FALSE).
#' @param family The family parameter that's passed to cv.glmnet(). Currently, only family='gaussian' is supported.
#' @param nfolds Number of cross-validation folds (default is 10) used to determine the optimal lambda in cv.glmnet().
#' @param foldid An optional vector indicating in which cross-validation fold each sample should be. Overrides nfolds when used.
#' @param lambda_upstream For the first stage (using the upstream features), should glmnet use lambda.min or lambda.1se? Default is lambda.1se.
#' @param lambda_downstream For the second stage (using the downstream features), should glmnet use lambda.min or lambda.1se? Default is lambda.1se.
#' @param ... Other parameters that are passed to cv.glmnet().
#' @return A tandem-object.
#' @examples
#' # unpack example data
#' x = example_data$x
#' y = example_data$y
#' upstream = example_data$upstream
#'
#' # fit a tandem model, determine the coefficients and create a prediction
#' fit = tandem(x, y, upstream, alpha=0.5)
#' beta = coef(fit)
#' y_hat = predict(fit, newx=x)
#' @export
tandem = function(x, y, upstream, family="gaussian", nfolds=10, foldid=NULL, lambda_upstream="lambda.1se", lambda_downstream="lambda.1se", ...) {
  if(class(x)!="matrix" | !class(x[1,1]) %in% c("numeric","integer"))
    stop("x needs to be a numeric matrix")
  if(!class(y) %in% c("numeric", "integer"))
    stop("y needs to be a numeric vector")
  if(class(upstream)!="logical")
    stop("upstream needs to a logical index vector, integer index vectors are currently not supported")
  if(nrow(x)!=length(y))
    stop("Number of samples in x and y don't match")
  if(ncol(x)!=length(upstream))
    stop("Number of features in x and upstream don't match")
  if(family!="gaussian")
    stop("Currently only the glmnet family='guassian' is supported")
  if(any(is.na(x)))
    stop("NAs in x are not allowed")
  if(any(is.na(y)))
    stop("NAs in y are not allowed")
  if(nrow(x)<nfolds)
    stop("nfolds should be smaller than the number of samples")
  if(!lambda_upstream %in% c("lambda.1se", "lambda.min"))
    stop("lambda_upstream should be either lambda.1se or lambda.min")
  if(!lambda_downstream %in% c("lambda.1se", "lambda.min"))
    stop("lambda_downstream should be either lambda.1se or lambda.min")

  if(is.null(foldid)) {
    n = length(y)
    foldid = ceiling(sample(1:n)/n * nfolds)
  }

  fit1 = glmnet::cv.glmnet(x[,upstream], y, foldid=foldid, ...)
  residuals = y - stats::predict(fit1, newx=x[,upstream], s=lambda_upstream)
  fit2 = glmnet::cv.glmnet(x[,!upstream], residuals, foldid=foldid, ...)

  beta0 =stats::coef(fit1, s=lambda_upstream)[1] +stats::coef(fit2, s=lambda_downstream)[1]
  beta = matrix(NA, ncol(x), 1)
  beta[upstream,] = as.matrix(stats::coef(fit1, s=lambda_upstream)[-1,,drop=F])
  beta[!upstream,] = as.matrix(stats::coef(fit2, s=lambda_downstream)[-1,,drop=F])
  rownames(beta) = colnames(x)
  beta = Matrix::Matrix(beta)
  fit = list(beta0=beta0, beta=beta)
  class(fit) = "tandem"

  return(fit)
}

#' Creates a prediction using a tandem-object
#'
#' Creates a prediction using a tandem-object.
#'
#' @param object A tandem-object, as returned by tandem()
#' @param newx A feature matrix, where the rows correspond to samples and the columns to features.
#' @param ... Not used. Other arguments for predict().
#' @return The predicted response vector.
#' @examples
#' # unpack example data
#' x = example_data$x
#' y = example_data$y
#' upstream = example_data$upstream
#'
#' # fit a tandem model, determine the coefficients and create a prediction
#' fit = tandem(x, y, upstream, alpha=0.5)
#' beta = coef(fit)
#' y_hat = predict(fit, newx=x)
#' @export
predict.tandem = function(object, newx, ...) {
  x = newx

  if(class(object)!="tandem")
    stop("object needs to be a tandem-object")
  if(class(x)!="matrix" | !class(x[1,1]) %in% c("numeric","integer"))
    stop("x needs to be a numeric matrix")
  if(nrow(object$beta)!=ncol(x))
    stop("Number of features in the fitted tandem object does not match number of features in x")
  if(any(is.na(x)))
    stop("NAs in x are not allowed")

  beta0 = object$beta0
  beta = object$beta
  y = x%*%beta + beta0

  return(y)
}

#' Estimating predictive performance via nested cross-validation
#'
#' Performs a nested cross-validation to assess the predictive performance. The inner loop is used to determine the optimal lambda
#' (as in cv.glmnet) and the outer loop is used to asses the predictive performance in an unbiased way.
#'
#' @param x A feature matrix, where the rows correspond to samples and the columns to features.
#' @param y A vector containing the response.
#' @param upstream A logical index vector that indicates for each feature whether it's upstream (TRUE) or downstream (FALSE).
#' @param method Indicates whether the nested cross-validation is performed on TANDEM or on the classic approach (glmnet). Should be either "tandem" or "glmnet".
#' @param family The family parameter that's passed to cv.glmnet(). Currently, only family='gaussian' is supported.
#' @param nfolds Number of cross-validation folds (default is 10) used in the outer cross-validation loop.
#' @param nfolds_inner Number of cross-validation folds (default is 10) used to determine the optimal lambda in the inner cross-validation loop.
#' @param foldid An optional vector indicating in which cross-validation fold each sample should be in the outer cross-validation loop. Overrides nfolds when used.
#' @param lambda_upstream Only used when method='tandem'. For the first stage (using the upstream features), should glmnet use lambda.min or lambda.1se? Default is lambda.1se.
#' @param lambda_downstream Only used when method='tandem'. For the second stage (using the downstream features), should glmnet use lambda.min or lambda.1se? Default is lambda.1se.
#' @param lambda_glmnet Only used when method='glmnet'. Should glmnet use lambda.min or lambda.1se? Default is lambda.1se.
#' @param ... Other parameters that are passed to cv.glmnet().
#' @return The predicted response vector y_hat and the mean-squared error (MSE).
#' @examples
#' # unpack example data
#' x = example_data$x
#' y = example_data$y
#' upstream = example_data$upstream
#'
#' # assess the prediction error in a nested cv-loop
#' # fix the seed to have the same foldids between the two methods
#' set.seed(1)
#' cv_tandem = nested.cv(x, y, upstream, method="tandem", alpha=0.5)
#' set.seed(1)
#' cv_glmnet = nested.cv(x, y, upstream, method="glmnet", alpha=0.5)
#' barplot(c(cv_tandem$mse, cv_glmnet$mse), ylab="MSE", names=c("TANDEM", "Classic approach"))
#' @export
nested.cv = function(x, y, upstream, method="tandem", family="gaussian", nfolds=10, nfolds_inner=10, foldid=NULL,
                     lambda_upstream="lambda.1se", lambda_downstream="lambda.1se", lambda_glmnet="lambda.1se",...) {
  if(class(x)!="matrix" | !class(x[1,1]) %in% c("numeric","integer"))
    stop("x needs to be a numeric matrix")
  if(!class(y) %in% c("numeric", "integer"))
    stop("y needs to be a numeric vector")
  if(class(upstream)!="logical")
    stop("upstream needs to a logical index vector, integer index vectors are currently not supported")
  if(nrow(x)!=length(y))
    stop("Number of samples in x and y don't match")
  if(ncol(x)!=length(upstream))
    stop("Number of features in x and upstream don't match")
  if(family!="gaussian")
    stop("Currently only the glmnet family='guassian' is supported")
  if(any(is.na(x)))
    stop("NAs in x are not allowed")
  if(any(is.na(y)))
    stop("NAs in y are not allowed")
  if(!method %in% c("tandem", "glmnet"))
    stop("method should equal 'tandem' or 'glmnet'")
  if(!lambda_upstream %in% c("lambda.1se", "lambda.min"))
    stop("lambda_upstream should be either lambda.1se or lambda.min")
  if(!lambda_downstream %in% c("lambda.1se", "lambda.min"))
    stop("lambda_downstream should be either lambda.1se or lambda.min")
  if(!lambda_glmnet %in% c("lambda.1se", "lambda.min"))
    stop("lambda_glmnet should be either lambda.1se or lambda.min")

  if(is.null(foldid)) {
    n = length(y)
    foldid = ceiling(sample(1:n)/n * nfolds)
  }

  y_hat = rep(NA, length(y))
  for(i in 1:nfolds) {
    ind = foldid==i
    x_train = x[!ind,]
    y_train = y[!ind]
    x_test = x[ind,]

    n_i = length(y_train)
    foldid_i = ceiling(sample(1:n_i)/n_i * nfolds_inner)

    if(method=="tandem") {
      fit = tandem(x_train, y_train, upstream, family=family, foldid=foldid_i, lambda_upstream=lambda_upstream, lambda_downstream=lambda_downstream, ...)
      y_hat[ind] = as.vector(predict.tandem(fit, newx=x_test))
    } else {
      fit = glmnet::cv.glmnet(x_train, y_train, family=family, foldid=foldid_i, ...)
      y_hat[ind] = as.vector(stats::predict(fit, newx=x_test, s=lambda_glmnet))
    }
  }
  mse = mean((y-y_hat)^2)

  return(list(mse=mse, y_hat=y_hat))
}

#' Returns the regression coefficients from a TANDEM fit
#'
#' Returns the regression coefficients from a TANDEM fit.
#'
#' @param object A tandem-object, as returned by tandem()
#' @param ... Not used. Other arguments for coef().
#' @return The regression coefficients.
#' @examples
#' # unpack example data
#' x = example_data$x
#' y = example_data$y
#' upstream = example_data$upstream
#'
#' # fit a tandem model, determine the coefficients and create a prediction
#' fit = tandem(x, y, upstream, alpha=0.5)
#' beta = coef(fit)
#' y_hat = predict(fit, newx=x)
#' @export
coef.tandem = function(object, ...) {
  if(class(object)!="tandem")
    stop("object needs to be a tandem-object")

  if(length(rownames(object$beta))==0) {
    rownames(object$beta) = 1:nrow(object$beta)
  }
  beta = rbind(object$beta0,object$beta)
  rownames(beta)[1] = "(Intercept)"
  return(beta)
}

#' Determine the relative contribution per data type
#'
#' For each data type, determine its relative contribution to the overall prediction.
#'
#' @param fit Either a tandem-object or a cv.glmnet-object
#' @param x The feature matrix used to train the fit, where the rows correspond to samples and the columns to features.
#' @param data_types A vector of the same length as the number of features, that indicates for each feature to which data
#' type it belongs. This vector doesn't need to correspond to the 'upstream' vector used in tandem(). For example, the upstream
#' features be spread across various data types (such as mutation, CNA, methylation and cancer type) and the downstream features
#' could be gene expression.
#' @param lambda_glmnet Only used when fit is a cv.glmnet object. Should glmnet use lambda.min or lambda.1se? Default is lambda.1se.
#' Note that for TANDEM objects, the lambda_upstream and lambda_downstream parameters should be specified during the tandem() call, as
#' they are used while fitting the model.
#' @return A vector that indicates the relative contribution per data type. These numbers sum up to one.
#' @examples
#' ## simple example
#' # unpack example data
#' x = example_data$x
#' y = example_data$y
#' upstream = example_data$upstream
#' data_types = example_data$data_types
#'
#' # fit TANDEM model
#' fit = tandem(x, y, upstream, alpha=0.5)
#'
#' # assess the relative contribution of upstream and downstream features
#' contr = relative.contributions(fit, x, data_types)
#' barplot(contr, ylab="Relative contribution", ylim=0:1)
#'
#' ## comparing TANDEM and classic model (glmnet)
#' # unpack example data
#' x = example_data$x
#' y = example_data$y
#' upstream = example_data$upstream
#' data_types = example_data$data_types
#'
#' # fix the cv folds, to facilitate a comparison between models
#' set.seed(1)
#' n = nrow(x)
#' nfolds = 10
#' foldid = ceiling(sample(1:n)/n * nfolds)
#'
#' # fit both a TANDEM and a classic model (glmnet)
#' fit = tandem(x, y, upstream, alpha=0.5)
#' library(glmnet)
#' fit2 = cv.glmnet(x, y, alpha=0.5, foldid=foldid)
#'
#' # assess the relative contribution of upstream and downstream features
#' # using both methods
#' contr_tandem = relative.contributions(fit, x, data_types)
#' contr_glmnet = relative.contributions(fit2, x, data_types)
#' par(mfrow=c(1,2))
#' barplot(contr_tandem, ylab="Relative contribution", main="TANDEM", ylim=0:1)
#' barplot(contr_glmnet, ylab="Relative contribution", main="Classic approach", ylim=0:1)
#' par(mfrow=c(1,1))
#' @export
relative.contributions = function(fit, x, data_types, lambda_glmnet="lambda.1se") {
  if(class(fit)!="tandem" & class(fit)!="cv.glmnet")
    stop("fit needs to be etiher a tandem-object or a cv.glmnet-object")
  if(ncol(x)!=length(data_types))
    stop("Data types should be a vector of length equal to the number of columns in x")

  x = scale(x, center=T, scale=F)
  p = ncol(x)
  betas = list()
  for(i in unique(data_types)) {
    beta = rep(0,p)
    ind = data_types==i
    if(class(fit)=="tandem") {
      beta[ind] = coef.tandem(fit)[-1][ind]
    } else if(class(fit)=="cv.glmnet") {
      beta[ind] =stats::coef(fit, s=lambda_glmnet)[-1][ind]
    }
    i = as.character(i)
    betas[[i]] = beta
  }

  y_hat = list()
  for(i in unique(data_types)) {
    i = as.character(i)
    y_hat[[i]] = x%*%betas[[i]]
  }

  a = rep(NA, length(unique(data_types)))
  names(a) = unique(data_types)
  for(i in unique(data_types)) {
    i = as.character(i)
    a[i] = sum(y_hat[[i]]^2)
    for(j in unique(data_types)) {
      j = as.character(j)
      if(i == j)
        next
      a[i] = a[i] + t(y_hat[[i]])%*%y_hat[[j]]
    }
    a[i] = abs(a[i])
  }
  a = a / sum(a)

  return(a)
}
