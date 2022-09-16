#' Internal bootstraping validation logistic model
#'
#' @description Validate logistic regression using bootstrap.
#' @param fit Object from glm fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param gamma Value of "gamma" parameter selected for relaxed model
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @importFrom stats predict glm formula
#' @export
vboot.glm <- function(fit, x = NULL, y = NULL, s = NULL, gamma = NULL, nfolds = NULL, B = 200, cv_replicates = NULL, n_cores = max(1, parallel::detectCores() - 1), save_ps=F){
    x <- fit$data[,!names(fit$data) %in% names(fit$model[1])]
    y <- fit$model[1][,1]

    orig_p <- predict(fit, type = "response")
    orig_measures <- score_binary(y = y, p = orig_p)
    # Making index to bootstrap
    bootstrap <- function(x, y, save_ps=save_ps){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]
        databoot <- data.frame(yboot, xboot)
        names(databoot)[1] <- names(fit$model)[1]
        # Fit the model using bootstrap dataset
        boot_fit <- glm(formula(fit), data = databoot, family = fit$family)
        Cb_boot_p <- predict.glm(boot_fit, newdata = databoot, type = "response")
        Cb_boot <- score_binary(y = yboot, p = Cb_boot_p)

        # fit bootstrap model to the original dataset
        Cb_orig_p <-  predict.glm(boot_fit, newdata = fit$data, type = "response")
        Cb_orig <- score_binary(y = y, p = Cb_orig_p)

        if (save_ps){
            ps = data.frame("yboot"=yboot, "Cb_boot_p"=Cb_boot_p, "Cb_orig_p"=Cb_orig_p)
        } else{
            ps = NULL
        }

        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, ps = ps))
    }

    if(n_cores > 1){
        CBOOT <- future.apply::future_replicate(B, tryCatch(bootstrap(x, y, save_ps = save_ps), error = function(e) NA))
    }
    else{
        CBOOT <- pbapply::pbreplicate(B, bootstrap(x, y, save_ps = save_ps))
    }

    # Optimism
    #Optimisms <- stats::na.omit(do.call(rbind, CBOOT[1,])) - stats::na.omit(do.call(rbind, CBOOT[2,]))
    Optimisms <- do.call(rbind, CBOOT[1,]) - do.call(rbind, CBOOT[2,])
    # eff_B <- length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    # O <- sapply(Optimisms,function(x) eff_B^-1 * sum(x))
    eff_B <- apply(Optimisms, 2, function(x) sum(!is.na(x))) # length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    O <- sapply(Optimisms, mean, na.rm=T)

    # Adjusted Optimist
    if(any(eff_B < B)) warning(paste("Due to sample size effective bootstraps samples are <", B, "see vboot$eff_B"))
    measures_adj <- orig_measures + O*-1
    output <- list(measures = data.frame(measures = names(orig_measures),
                                         Original_value = as.numeric(orig_measures),
                                         Mean_Optimism = O,
                                         Adjusted_value = as.numeric(measures_adj),
                                         row.names = NULL),
                   Effective_Bootstraps = eff_B,
                   Optimisms = Optimisms,
                   yorig = if(save_ps) y else NULL,
                   porig = if(save_ps) orig_p else NULL,
                   ps = if(save_ps) CBOOT[3,] else NULL
    )
    class(output) <- "bootVal"
    return(output)
}


#' Internal bootstraping validation linear model
#'
#' @description Validate linear regression using bootstrap.
#' @param fit Object from lm fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param gamma Value of "gamma" parameter selected for relaxed model
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @importFrom stats var predict lm formula
#' @export
vboot.lm <- function(fit, x = NULL, y = NULL, s = NULL, gamma = NULL, nfolds = NULL, B = 200, cv_replicates = NULL, n_cores = max(1, parallel::detectCores() - 1), save_ps=F){
    x <- fit$data[,!names(fit$data) %in% names(fit$model[1])]
    y <- fit$model[1][,1]

    orig_p <- predict(fit, type = "response")
    orig_measures <- score_linear(y = y, p = orig_p)
    # Making index to bootstrap
    bootstrap <- function(x, y, save_ps=save_ps){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]
        databoot <- data.frame(yboot, xboot)
        names(databoot)[1] <- names(fit$model)[1]
        # Fit the model using bootstrap dataset
        boot_fit <- lm(formula(fit), data = databoot)
        Cb_boot_p <- predict.lm(boot_fit, newdata = databoot, type = "response")
        Cb_boot <- score_linear(y = yboot, p = Cb_boot_p)

        # fit bootstrap model to the original dataset
        Cb_orig_p <-  predict.lm(boot_fit, newdata = fit$data, type = "response")
        Cb_orig <- score_linear(y = y, p = Cb_orig_p)

        if (save_ps){
            ps = data.frame("yboot"=yboot, "Cb_boot_p"=Cb_boot_p, "Cb_orig_p"=Cb_orig_p)
        } else{
            ps = NULL
        }

        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, ps = ps))
    }

    if(n_cores > 1){
        CBOOT <- future.apply::future_replicate(B, tryCatch(bootstrap(x, y, save_ps = save_ps), error = function(e) NA))
    }
    else{
        CBOOT <- pbapply::pbreplicate(B, bootstrap(x, y, save_ps = save_ps))
    }

    # Optimism
    #Optimisms <- stats::na.omit(do.call(rbind, CBOOT[1,])) - stats::na.omit(do.call(rbind, CBOOT[2,]))
    Optimisms <- do.call(rbind, CBOOT[1,]) - do.call(rbind, CBOOT[2,])
    # eff_B <- length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    # O <- sapply(Optimisms,function(x) eff_B^-1 * sum(x))
    eff_B <- apply(Optimisms, 2, function(x) sum(!is.na(x))) # length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    O <- sapply(Optimisms, mean, na.rm=T)

    # Adjusted Optimist
    if(any(eff_B < B)) warning(paste("Due to sample size effective bootstraps samples are <", B, "see vboot$eff_B"))
    measures_adj <- orig_measures + O*-1
    output <- list(measures = data.frame(measures = names(orig_measures),
                                         Original_value = as.numeric(orig_measures),
                                         Mean_Optimism = O,
                                         Adjusted_value = as.numeric(measures_adj),
                                         row.names = NULL),
                   Effective_Bootstraps = eff_B,
                   Optimisms = Optimisms,
                   yorig = if(save_ps) y else NULL,
                   porig = if(save_ps) orig_p else NULL,
                   ps = if(save_ps) CBOOT[3,] else NULL
    )
    class(output) <- "bootVal"
    return(output)
}

#' Score predicted probabilities
#'
#' @description Score predicted probabilities against observed outcomes
#' @param y vector containing binary outcome
#' @param p vector containing predicted probabilities
#' @return A `data.frame` containing binomial deviance, mean-squared error
#'    (brier score), mean absolute error, AUC, intercept and slope from a
#'    logistic calibration model (glm(y ~ qlogis(p), family="binomial"))
#' @importFrom pROC roc
#' @export
score_binary <- function(y, p){
    if (any(is.infinite(qlogis(p)))){
        y2 = y[!is.infinite(qlogis(p))]
        p2 = p[!is.infinite(qlogis(p))]
        cal <- stats::glm(y2 ~ qlogis(p2), family = 'binomial')
        warning("Preidicted probabilities of 0/1 removed for logistic caibration")
    } else{
        cal <- stats::glm(y ~ qlogis(p), family = 'binomial')
    }

    bd = mean(-2*(y*log(p) + (1 - y)*log(1-p)))

    scores <- data.frame("Binomial Deviance" = if(is.nan(bd)) NA else bd,
                         "Mean-Squared Error"  = mean((y - p)^2),
                         "Mean absolute Error"= mean(abs(y - p)),
                         "AUC" = suppressMessages(pROC::auc(response=y, predictor=p)[1]),
                         "Calib Intercept" = cal$coefficients[[1]],
                         "Calib Slope" = cal$coefficients[[2]]
    )
    return(scores)
}

#' Score predictions from linear model
#'
#' @description Score predictions from linear model
#' @param y vector containing outcome
#' @param p vector containing predictions
#' @return A `data.frame` containing R2, MSE, MAE, intercept and slope from a
#'    linear calibration model (lm(y ~ p))
#' @export
score_linear <- function(y, p){
    cal = lm(y ~ p)
    scores <- data.frame("R2" = 1 - (var(p - y) / var(y)),
                         "Mean-Squared Error"  = mean((y - p)^2),
                         "Mean absolute Error"= mean(abs(y - p)),
                         "Calib Intercept" = cal$coefficients[[1]],
                         "Calib Slope" = cal$coefficients[[2]]
    )
    return(scores)
}

#' Calculate AIC with a glmnetfit object
#'
#' @description Calculate AIC for glmnet
#' @param fit object fit with glmnet
#' @return
#' @export
AIC.glmnet <- function(fit){
    LL = fit$nulldev - deviance(fit)
    k = fit$df
    n = fit$nobs
    # with small sample size correction
    AIC = -LL+2*k+2*k*(k+1)/(n-k-1)

    return(AIC)
}


#' Internal bootstraping validation logistic glmnet model
#'
#' @description Validate glmnet logistic regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. Should be a factor with two levels
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param gamma Value of "gamma" parameter selected for relaxed model
#' @param nfolds Number of folds for cross validation as in cv.glmnet
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step in 'cv.glmnet'
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @param AIC_select use AIC (instead of cross vaidation) to select lambda? This should have been done to select the original lambdq (default = FALSE)
#' @param save_ps save the predicted probabilities from the bootstrapped models? This is used for non-linear calibration curves (default = FALSE)
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet
#' @importFrom pROC roc
#' @importFrom pbapply pbreplicate
#' @importFrom stats median coef na.omit
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @importFrom future.apply future_replicate
#' @export
vboot.glmnet <- function(fit, x, y, s, gamma = NULL, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1), AIC_select=FALSE, save_ps = FALSE){

    if (is.null(fit$call$alpha)){
        # running glmnet w/o specifying alpha leads to NULL
        fit$call$alpha <- 1
    }

    family = fit$call$family
    if (is.call(family)){
        family = as.character(family)[1]
        warning("Note: only family='binomial' or 'gaussian' are supported. Other calls (e.g., binomial(link='probit')) are not currently supported")
    }

    if (!family %in% c("gaussian", "binomial")) stop("glmnet family not supported")

    if (family == "binomial"){
        orig_p = plogis(glmnet::predict.glmnet(object = fit, newx = x, s = s, type = "link")[,1]) # type = "response" doesn't work as expected glmnet_4.1-4
        orig_measures <- score_binary(y = y, p = orig_p)
    } else if (family == "gaussian"){
        orig_p = glmnet::predict.glmnet(object = fit, newx = x, s = s, type = "response")[,1]
        orig_measures <- score_linear(y = y, p = orig_p)
    }

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, AIC_select = AIC_select, family=family, save_ps=save_ps){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        if (!AIC_select){
            # Fit the model using bootstrap dataset
            cv.glmnet_b <- suppressWarnings(pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha,
                                                                                                           family = family, nfolds = nfolds)$lambda.1se, error = function(e) NA)))
            l <- median(cv.glmnet_b, na.rm = TRUE)
            boot_fit <- suppressWarnings(tryCatch(glmnet::glmnet(xboot, yboot, alpha = alpha, family = family), error = function(e) NA))
        } else{
            boot_fit <- glmnet::glmnet(x = xboot, y = yboot, alpha = alpha, family = family)
            l <- boot_fit$lambda[which.min(AIC.glmnet(boot_fit))]
        }
        if (family == "binomial"){
            # eval boot model on boot data
            Cb_boot_p <- plogis(glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "link")[,1])
            Cb_boot <- score_binary(y = yboot, p = Cb_boot_p)

            # eval boot model on orig data
            Cb_orig_p <- plogis(glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "link")[,1])
            Cb_orig <- score_binary(y = y, p = Cb_orig_p)
        } else if (family == "gaussian"){
            # eval boot model on boot data
            Cb_boot_p <- glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "response")[,1]
            Cb_boot <- score_linear(y = yboot, p = Cb_boot_p)

            # eval boot model on orig data
            Cb_orig_p <- glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "response")[,1]
            Cb_orig <- score_linear(y = y, p = Cb_orig_p)
        }

        #selected variables
        boot_coefs <- glmnet::coef.glmnet(boot_fit, s = l)[-1, ]
        var_sel <- names(boot_coefs[boot_coefs!=0])

        if (save_ps){
            ps = data.frame("yboot"=yboot, "Cb_boot_p"=Cb_boot_p, "Cb_orig_p"=Cb_orig_p)
        } else{
            ps = NULL
        }

        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel, ps = ps))
    }
    if(n_cores > 1){
        CBOOT <- future.apply::future_replicate(B, tryCatch(bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds, AIC_select = AIC_select, family = family, save_ps = save_ps), error = function(e) NA))
    }
    else{
        CBOOT <- pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds, AIC_select = AIC_select, family = family, save_ps = save_ps))
        }

    # Optimism
    #Optimisms <- stats::na.omit(do.call(rbind, CBOOT[1,])) - stats::na.omit(do.call(rbind, CBOOT[2,]))
    Optimisms <- do.call(rbind, CBOOT[1,]) - do.call(rbind, CBOOT[2,])
    # eff_B <- length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    # O <- sapply(Optimisms,function(x) eff_B^-1 * sum(x))
    eff_B <- apply(Optimisms, 2, function(x) sum(!is.na(x))) # length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    O <- sapply(Optimisms, mean, na.rm=T)

    # Adjusted Optimist
    if(any(eff_B < B)) warning(paste("Due to sample size effective bootstraps samples are <", B, "see vboot$eff_B"))
    measures_adj <- orig_measures + O*-1
    output <- list(measures = data.frame(measures = names(orig_measures),
                                         Original_value = as.numeric(orig_measures),
                                         Mean_Optimism = O,
                                         Adjusted_value = as.numeric(measures_adj),
                                         row.names = NULL),
                   varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B,
                   Optimisms = Optimisms,
                   yorig = if(save_ps) y else NULL,
                   porig = if(save_ps) orig_p else NULL,
                   ps = if(save_ps) CBOOT[4,] else NULL
                   )
    class(output) <- "bootVal"
    return(output)
}

#' Internal bootstraping validation for relaxed logistic glmnet model
#'
#' @description Validate glmnet relaxed logistic regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. Should be a factor with two levels
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param gamma Value of "gamma" parameter selected for relaxed model
#' @param nfolds Number of folds for cross validation as in cv.glmnet
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step in 'cv.glmnet'
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet
#' @importFrom pROC roc
#' @importFrom pbapply pbreplicate
#' @importFrom stats median coef na.omit
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @importFrom future.apply future_replicate
#' @importFrom MASS kde2d
#' @export
vboot.relaxed <- function(fit, x, y, s, gamma, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    orig_asses <- glmnet::assess.glmnet(fit$relaxed, newx = x, newy = y, s = s, gamma = gamma ,family = fit$call$family)
    orig_measures <- data.frame("Binomial Deviance" = orig_asses$deviance,
                                "Misclassification Error" = orig_asses$class,
                                "Mean-Squared Error"  = orig_asses$mse,
                                "Mean absolute Error"= orig_asses$mae,
                                "AUC" = orig_asses$auc)

    # Making index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- suppressWarnings(pbapply::pbreplicate(cv_replicates, tryCatch(glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha,
                                                                                                       family = "binomial", nfolds = nfolds, relax = TRUE)$relaxed[c(4,6)], error = function(e) NA)))
        density2d <- MASS::kde2d(na.omit(unlist(cv.glmnet_b[1,])), na.omit(unlist(cv.glmnet_b[2,])), h = 1)
        position <- which(density2d$z == max(density2d$z), arr.ind = TRUE)
        l <- density2d$x[position[1]]
        g <- density2d$y[position[2]]
        boot_fit <- suppressWarnings(tryCatch(glmnet::glmnet(xboot, yboot, alpha = fit$call$alpha, family = "binomial", relax = TRUE), error = function(e) NA))
        Cb_boot_asses <- glmnet::assess.glmnet(boot_fit$relaxed, newx = xboot, newy = yboot, s = l, gamma = g, family = fit$call$family)
        Cb_boot <- data.frame("Binomial Deviance" = Cb_boot_asses$deviance,
                              "Misclassification Error" = Cb_boot_asses$class,
                              "Mean-Squared Error"  = Cb_boot_asses$mse,
                              "Mean absolute Error"= Cb_boot_asses$mae,
                              "AUC" = Cb_boot_asses$auc)

        #selected variables
        var_sel <- colnames(x)[glmnet::predict.glmnet(boot_fit$relaxed,newx = xboot, s = l, gamma = g, type = "nonzero")[[1]]]

        # fit bootstrap model to the original dataset
        Cb_orig_asses <- glmnet::assess.glmnet(boot_fit$relaxed, newx = x, newy = y, s = l, gamma = g)
        Cb_orig <- data.frame("Binomial Deviance" = Cb_orig_asses$deviance,
                              "Misclassification Error" = Cb_orig_asses$class,
                              "Mean-Squared Error"  = Cb_orig_asses$mse,
                              "Mean absolute Error"= Cb_orig_asses$mae,
                              "AUC" = Cb_orig_asses$auc)

        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel))
    }
    if(n_cores > 1){
        CBOOT <- future.apply::future_replicate(B, tryCatch(bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds), error = function(e) NA))
    }
    else{
        CBOOT <- pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds))
    }

    # Optimist
    Optimisms <- stats::na.omit(do.call(rbind, CBOOT[1,])) - stats::na.omit(do.call(rbind, CBOOT[2,]))
    eff_B <- length(stats::na.omit(do.call(rbind, CBOOT[2,]))[,1])
    O <- sapply(Optimisms,function(x) eff_B^-1 * sum(x))

    # Adjusted Optimist
    if(eff_B  != B) warning(paste("Due to sample size effective bootstraps samples are",eff_B))
    measures_adj <- orig_measures + O*-1
    output <- list(measures = data.frame(measures = c("Binomial Deviance",
                                                      "Misclassification Error",
                                                      "Mean-Squared Error",
                                                      "Mean absolute Error",
                                                      "AUC"),
                                         Original_value = as.numeric(c(orig_asses$deviance,
                                                                       orig_asses$class,
                                                                       orig_asses$mse,
                                                                       orig_asses$mae,
                                                                       orig_asses$auc)),
                                         Mean_Optimism = O,
                                         Adjusted_value = c(measures_adj$Binomial.Deviance,
                                                            measures_adj$Misclassification.Error,
                                                            measures_adj$Mean.Squared.Error,
                                                            measures_adj$Mean.absolute.Error,
                                                            measures_adj$AUC),
                                         row.names = NULL),
                   varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B,
                   Optimisms = Optimisms)
    class(output) <- "bootVal"
    return(output)
}

#' Internal bootstraping validation linear glmnet model
#'
#' @description Validate glmnet linear regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. Should be numeric
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'
#' @param gamma Value of "gamma" parameter selected for relaxed model
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @importFrom glmnet cv.glmnet glmnet predict.glmnet
#' @importFrom pROC roc
#' @importFrom pbapply pbreplicate
#' @importFrom stats median var coef
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
vboot.elnet <- function(fit, x, y, s, gamma = NULL, nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1)){
    orig_predict <- glmnet::predict.glmnet(fit, newx = x, s=s, type = "response")
    orig_r2 <- 1 - (var(orig_predict - y) / var(y))
    # Create index to bootstrap
    bootstrap <- function(x, y, alpha = fit$call$alpha, nfolds = nfolds, B = B){
        index <- sample(1:nrow(x), replace = TRUE)
        xboot <- x[index, ]
        yboot <- y[index]

        # Fit the model using bootstrap dataset
        cv.glmnet_b <- pbapply::pbreplicate(cv_replicates, glmnet::cv.glmnet(xboot, yboot, alpha = fit$call$alpha, family = "gaussian", nfolds = nfolds)$lambda.1se)
        l <- median(cv.glmnet_b)
        boot_fit <- glmnet::glmnet(xboot, yboot, alpha = fit$call$alpha, family = "gaussian")
        boot_predict <- glmnet::predict.glmnet(boot_fit, newx = xboot, s = l, type = "response")
        Cb_boot <- 1 - (var(boot_predict - yboot) / var(yboot))
        #selected variables
        var_sel <- names(coef(fit, s = s)[coef(fit, s = s)[,1] != 0,])[-1]

        # fit bootstrap model to the original dataset
        bootorig_predict <-  glmnet::predict.glmnet(boot_fit, newx = x, s = l, type = "response")
        Cb_orig <- 1 - (var(bootorig_predict - y) / var(y))
        return(list(Cb_boot = Cb_boot, Cb_orig = Cb_orig, var_sel = var_sel))
    }
    if( n_cores > 1){
        cl <- parallel::makeCluster(n_cores)
        parallel::clusterExport(cl, varlist = c("B", "x", "y", "fit", "nfolds", "bootstrap"), envir = environment())
        CBOOT <- parallel::parSapply(cl, 1:B, function(i) bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds))
        parallel::stopCluster(cl)
        closeAllConnections()
    }
    else{
        CBOOT <- suppressWarnings(pbapply::pbreplicate(B, bootstrap(x, y, alpha = fit$call$alpha, nfolds = nfolds)))
    }
    # Optimist
    Optimisms <- as.numeric(stats::na.omit(unlist(CBOOT[1, ])) - stats::na.omit(unlist(CBOOT[2, ])))
    eff_B <- length( stats::na.omit(unlist(CBOOT[2, ])))
    O <- B^-1 * sum(unlist(CBOOT[1,]) - unlist(CBOOT[2, ]))
    # Adjusted Optimist
    Oadj <- as.numeric(orig_r2 - O)
    #output
    output <- list(round(data.frame(Original_R2 = as.numeric(orig_r2), Optimism = O, Validated_R2 = Oadj),3), varImportance = CBOOT[3,],
                   Effective_Bootstraps = eff_B, Optimisms = Optimisms)
    class(output) <- "bootVal"
    return(output)

}

#' Print function
#'
#' @description internal function to print vboot object
#' @param x A vboot object
#' @param ... Further arguments of generic function
#' @return validation
#' @export
print.bootVal <- function(x, ...){
    print(x[[1]])
}


#' Plot for repeat_cv results
#'
#' @description Plots a grid of slices from the estimates of the repeat_cv function
#' @param x A vboot object
#' @param order order plot by importance
#' @param n Variables to be displayed
#' @param ... further arguments passed to plot
#' @importFrom graphics barplot
#' @return validation
#' @export
plot.bootVal <- function(x, order = TRUE, n = length(unique(unlist(x$varImportance))), ...){
    raw_table <- table(unlist(x[[2]]))[1:n]/x$Effective_Bootstraps
    if(order) raw_table <- raw_table[order(unlist(x[[2]]), decreasing = TRUE)]
    barplot(raw_table, horiz = TRUE,
            las = 1, xlim = c(0,1), main = "Importance variable", xlab = "Inclusion porportion", ...)
}


#' Plot non-linear (lowess) calibration curve
#'
#' @description Plots non-linear calibration curve using lowess smoother
#'    Also plots optimism corrected calibration curve using method described in Harrell (2015) p. 270
#'    modified from `rms::calibrate`.
#' @param x A bootVal object (must have been run with save_ps = TRUE)
#'
#' @references Harrell Jr, F. E. (2015). Regression modeling strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. Springer.
#' @return A plot of apparent and optimism corrected calibration curves
#' @return a list of original lowess smooth and optimisms
#' @export
calib.bootVal <- function(x){
    if (is.null(x$yorig) | is.null(x$ps) | is.null(x$porig))
        stop("model predictions and/or yorig missing. Run vboot with save_ps = TRUE")

    smooth <- function(p, y, porig){
        # adapted from https://github.com/harrelfe/rms/blob/master/R/calibrate.default.s
        smo = lowess(p, y, iter = 0)
        intsmo = approx(x = smo$x, y = smo$y, xout=porig, ties=function(x)x[1])$y
        return(intsmo)
    }

    orig_smooth = lowess(x$porig, x$yorig, iter = 0)

    optimisms = sapply(x$ps, function(xx){
        boot_smooth = tryCatch(smooth(p = xx$Cb_boot_p, y = xx$yboot, porig = orig_smooth$x), error=function(e) NA)
        bootorig_smooth = tryCatch(smooth(p = xx$Cb_orig_p, y = x$yorig, porig = orig_smooth$x), error=function(e) NA)

        # distance to 45^o
        optimism = (boot_smooth - orig_smooth$x) - (bootorig_smooth - orig_smooth$x)
    })

    missings = apply(optimisms, 2, function(x) mean(is.na(x)))
    if (any(missings > .25)) warning("Large proportion missing")

    O = apply(optimisms, MARGIN = 1, FUN = mean, na.rm=T)

    plot(NA, type = "l", xlim=c(0, 1),# max(orig_smooth$x)),
         ylim=c(0, 1),
         xlab="Predicted Probability", ylab="Actual Probability")
    abline(0,1, lwd=4, col="lightgrey")
    lines(orig_smooth)

    lines(orig_smooth$x, orig_smooth$y - O, lty=2)

    legend(x = 0, y = 1, legend = c("Apparent", "Optimism Corrected"), lty=c(1,2))

    h = hist(x$porig, plot = F, breaks = 100)

    points(x = h$mids, y = .1*h$density/max(h$density), type = "h")

    return(invisible(list("orig_smooth" = orig_smooth, "optimism"= O)))
}


#' Generic function for bootstrap validation
#'
#' @description Validate 'glmnet' linear, logistic or cox regression using bootstrap.
#' @param fit Object from glmnet fit
#' @param x A matrix of the predictors, each row is an observation vector.
#' @param y A vector of response variable. It should be quantitative for lineal regression, a factor with two levels for logistic regression or a two-column matrix with columns named 'time' and 'status' for cox regression.
#' @param s Value of the penalty parameter "lambda" selected from the original 'cv.glmnet'.
#' @param gamma Value of "gamma" parameter selected for relaxed model
#' @param nfolds Number of folds for cross validation as in 'cv.glmnet'.
#' @param B Number of bootsrap samples
#' @param cv_replicates Number of replicates for the cross-validation step
#' @param n_cores number of cores to use in parallel. Default detectCores()-1
#' @references Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). Regularization Paths for Generalized Linear Models via Coordinate Descent. Journal of Statistical Software, 33(1), 1-22. URL http://www.jstatsoft.org/v33/i01/.
#' @references Noah Simon, Jerome Friedman, Trevor Hastie, Rob Tibshirani (2011). Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software, 39(5), 1-13. URL http://www.jstatsoft.org/v39/i05/.
#' @references Harrell Jr, F. E. (2015). Regression modeling strategies: with applications to linear models, logistic and ordinal regression, and survival analysis. Springer.
#' @references Gordon C.S. Smith, Shaun R. Seaman, Angela M. Wood, Patrick Royston, Ian R. White (2014). Correcting for Optimistic Prediction in Small Data Sets, American Journal of Epidemiology, Volume 180, Issue 3, 1 August 2014, Pages 318-324, https://doi.org/10.1093/aje/kwu140
#' @importFrom pROC roc multiclass.roc
#' @importFrom survival survfit Surv
#' @importFrom pbapply pbreplicate
#' @importFrom stats median var predict.glm predict.lm glm lm coef formula
#' @importFrom graphics barplot
#' @importFrom parallel parSapply makeCluster detectCores clusterExport stopCluster
#' @export
#' @examples
#' # Create the data
#' set.seed(25)
#' x <- matrix(rnorm(80),ncol=4)
#' y <- x[,4]*0.8+x[,3]*0.4+rnorm(20)
#' # Fit glmnet model
#' fit_enet <- glmnet::glmnet(x, y, alpha = 0.5)
#' # Bootstrap validation
#' vboot(fit_enet, x, y, nfolds = 3, B = 2, s = 0.5, cv_replicates = 5, n_cores = 1)

vboot <- function(fit, x, y, s, gamma = FALSE,nfolds = 5, B = 200, cv_replicates = 100, n_cores = max(1, parallel::detectCores() - 1), AIC_select=FALSE, save_ps=FALSE){
    UseMethod("vboot")
}


