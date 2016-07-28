demogSmooth.wrapper = function(...)
{
  tryCatch((result = demogSmooth(...)),
           error = function(e) {
             cat("\n\nERROR\n\n")
             print(e)
             dots <- list(...)
             if(length(dots) != 0) {
               cat("Parameters:\n")
               print(dots)
             }
           })
  return(result$result + ifelse(is.null(result$yearsEffect), 0, result$yearsEffect) + ifelse(is.null(result$cohortEffect), 0, result$cohortEffect))
}


#' Estimates smoothing parameters
#'
#' @param data Demographic data presented as a matrix.
#' @param effects Controls if the cohort and period effects are taking into account.
#' @param cornerLength Sets the smallest length of a digonal to be considered for cohort effects. The first diagonal is at the lowerest left conner of data matrix.
#' @param affdDiagonals Diagonals to be used for cohort effects.
#' @param affdYears Years to be used for period effects.
#' @param parameters Optional vector with some initial values to replace values in init parameter.
#' @param lower Lowest possible values for the optimisation procedure.
#' @param upper Highest possible values for the optimisation procedure.
#' @param init Initial values for the optimisation procedure.
#' @param reltol Relative tolerance parameter to be supplied to optim function (standard optimiser for R).
#' @param trace Controls if tracing is on.
#' @param control The control data passed directly to rq.fit.sfn method (quantreg package).
#' @return A vector of optimal smoothing parameters.
#' @examples
#' # library(demography)
#' # m = log(fr.mort$rate$female[1:30, 150:160])
#' # parameters = estPar(m)
#' @references \url{https://business.monash.edu/econometrics-and-business-statistics/research/publications/ebs/two-dimensional_smoothing_of_mortality_rates..pdf}
#' @author Alex Dokumentov
#' @export

estPar = function(data,
                  effects = TRUE,
                  cornerLength = 7,
                  affdDiagonals = NULL,
                  affdYears = NULL,
                  parameters = NULL,
                  lower = head(c(0.01, 0.01, 0.01, 2.0, 0.01, 2.0, 0.01), 3 + effects*4),
                  upper = head(c(1.2,  1.8,  1.2,  12,  0.4,  12,  0.4), 3 + effects*4),
                  init =  head(c(0.1,  0.1,  0.2,  4,   0.05, 4,   0.05), 3 + effects*4),
                  reltol = 0.001,
                  trace = F,
                  control = list(nnzlmax = 1000000, nsubmax = 2000000, tmpmax = 200000))
{
  counter = 0

  f = function(x)
  {
    if(trace) {
      print("Function f(x). x:")
      print(x)
    } else {
      counter <<- (counter + 1) %% 4
      cat("\r"); cat(paste0(c("\\","|","/","-")[counter + 1], "   "))
    }
    if(length(x) != length(lower)) stop("Error: length(x) != length(lower)")
    xx = pmin(pmax(x, lower), upper)
    if(sum(abs(xx - x)) > 0) {
      if(trace) {
        print("Beyond lower or upper bounds")
        print("smoothCv result: Inf")
      }
      return(Inf)
    }
    lambdaaa <- x[1]
    lambdayy <- x[2]
    lambdaay <- x[3]
    lambdaYearsEffect <- x[4]
    thetaYearsEffect <- x[5]
    lambdaCohortEffect <- x[6]
    thetaCohortEffect <- x[7]
    # time = system.time({
    cv <- smoothCv(demogSmooth.wrapper, data = data,
                   lambda = 1, lambdaaa = lambdaaa, lambdayy = lambdayy, lambdaay = lambdaay,
                   lambdaYearsEffect = lambdaYearsEffect, thetaYearsEffect = thetaYearsEffect,
                   lambdaCohortEffect = lambdaCohortEffect, thetaCohortEffect = thetaCohortEffect,
                   cornerLength = cornerLength, effects = effects,
                   affdDiagonals = affdDiagonals, affdYears = affdYears,
                   control = control)
    # })
    if(trace) {
      # print(time)
      print("smoothCv result:")
      print(cv$MAE)
    }
    return(cv$MAE)
  }

  if(!is.null(parameters)) {
    for(i in 1L:min(length(parameters), length(init))) {
      init[i] = min(max(parameters[i], lower[i]), upper[i])
    }
  }

  result = optim(par = init, fn = f, control = list(reltol = reltol))
  return(list(par = result$par, cv = result$value))
}

estAA = function(data,
                 lower = 0.2,
                 upper = 10,
                 step = 0.2,
                 trace = trace,
                 control = list(nnzlmax = 1000000, nsubmax = 2000000, tmpmax = 200000))
{
  lambdas = seq(lower, upper, by = step)
  cv = 0
  for(i in seq_along(lambdas)) {
    if(!trace) {cat("\r "); cat(paste0(c("\\","|","/","-")[i %% 4 + 1], "   "))}
    cv[i] <- smoothCv(demogSmooth.wrapper, data = data,
                      lambda = 1, lambdaaa = lambdas[i], lambdayy = 0, lambdaay = 0,
                      effects = FALSE, control = control)$MAE
    if(trace) {print(paste("lambdaAA:", lambdas[i])); print(cv[i])}
  }
  return(lambdas[which.min(cv)])
}

estYY = function(data,
                 lower = 0.2,
                 upper = 10,
                 step = 0.2,
                 trace = trace,
                 control = list(nnzlmax = 1000000, nsubmax = 2000000, tmpmax = 200000))
{
  lambdas = seq(lower, upper, by = step)
  cv = 0
  for(i in seq_along(lambdas)) {
    if(!trace) {cat("\r "); cat(paste0(c("\\","|","/","-")[i %% 4 + 1], "   "))}
    cv[i] <- smoothCv(demogSmooth.wrapper, data = data,
                      lambda = 1, lambdaaa = 0, lambdayy = lambdas[i], lambdaay = 0,
                      effects = FALSE, control = control)$MAE
    if(trace) {print(paste("lambdaYY:", lambdas[i])); print(cv[i])}
  }
  return(lambdas[which.min(cv)])
}

#' Smoothes demographic data using automatically estimated parameters and optionally
#' taking into account period and cohort effects
#'
#' If period and cohort effects are taken into account (effects = TRUE) the method uses all
#' available years and diagonals for estimation of the period and cohort effects.
#'
#' @seealso \code{\link{twoStepDemogSmooth}} might give slightly better performance.
#' @param data Demographic data presented as a matrix.
#' @param effects Controls if the cohort and period effects are taking into account.
#' @param cornerLength Sets the smallest length of a digonal to be considered for cohort effects.
#' @param affdDiagonals Diagonals to be used for cohort effects. The first diagonal is at the lowerest left conner of data matrix.
#' @param affdYears Years to be used for period effects.
#' @param lower Lowest possible values for the optimisation procedure.
#' @param upper Highest possible values for the optimisation procedure.
#' @param init Initial values for the optimisation procedure.
#' @param reltol Relative tolerance parameter to be supplied to optim function (standard optimiser for R).
#' @param parameters Optional parameters (output of estPar function). If not provied the parameters are estimated inside the function.
#' @param trace Controls if tracing is on.
#' @param control The control data passed directly to rq.fit.sfn method (quantreg package).
#' @return A list with four components: smooth surface, period effects, cohort effects and parameters
#' used for smoothing (passed as a parameter or estimated).
#' @examples
#' # library(demography)
#' # m = log(fr.mort$rate$female[1:30, 150:160])
#' # Show(m)
#' # sm = autoDemogSmooth(m)
#' # plot(sm)
#' # plot(sm, "period")
#' # plot(sm, "cohort")
#' @references \url{https://business.monash.edu/econometrics-and-business-statistics/research/publications/ebs/two-dimensional_smoothing_of_mortality_rates..pdf}
#' @author Alex Dokumentov
#' @export

autoDemogSmooth = function(data,
                           effects = TRUE,
                           cornerLength = 7,
                           affdDiagonals = NULL,
                           affdYears = NULL,
                           lower = head(c(0.01, 0.01, 0.01, 2.0, 0.001, 2.0, 0.001), 3 + effects*4),
                           upper = head(c(1.2,  1.8,  1.2,  12,  0.4,  12,  0.4), 3 + effects*4),
                           init =  head(c(0.1,  0.1,  0.2,  4,   0.01, 4,   0.01), 3 + effects*4),
                           reltol = 0.001,
                           parameters = NULL,
                           trace = F,
                           control = list(nnzlmax = 1000000, nsubmax = 2000000, tmpmax = 200000))
{
  if(missing(parameters)) {
    parameters = estPar(data,
                        effects = effects,
                        cornerLength = cornerLength,
                        affdDiagonals = affdDiagonals,
                        affdYears = affdYears,
                        lower = lower,
                        upper = upper,
                        init =  init,
                        reltol = reltol,
                        trace = trace,
                        control = control)$par
  }
  result = demogSmooth(data,
                       lambda = 1,
                       lambdaaa = parameters[1],
                       lambdayy = parameters[2],
                       lambdaay = parameters[3],
                       lambdaYearsEffect = parameters[4],
                       thetaYearsEffect = parameters[5],
                       lambdaCohortEffect = parameters[6],
                       thetaCohortEffect = parameters[7],
                       cornerLength = cornerLength,
                       effects = effects,
                       control = control)
  result$parameters = parameters
  return(result)
}

my.t.test = function(x, alternative = "two.sided")
{
  x = as.numeric(na.omit(x))
  if(length(x) >= 2) {
    return(t.test(x, alternative = alternative)$p.value)
  }
  else {
    return(NA)
  }
}

getAffected = function(resid, p.value = 0.05)
{
  d = diags(resid)
  p.values.t.1 = vapply(2:(nrow(d)-1), function(i) my.t.test(x=d[i,]), c(p.value=0))
  affdDiagonals1 = (2:(nrow(d)-1))[p.values.t.1 <= 0.05]
  p.values.t.2 = vapply(3:(nrow(d)-2), function(i) my.t.test(x=d[i,-1]*d[i,-ncol(d)], alternative="greater"), c(p.value=0))
  affdDiagonals2 = (3:(nrow(d)-2))[p.values.t.2 <= 0.05]
  affdDiagonals = sort(union(affdDiagonals1, affdDiagonals2))

  p.values.t.3 = vapply(1:ncol(resid), function(j) my.t.test(x=resid[,j]), c(p.value=0))
  affdYears1 = (1:ncol(resid))[p.values.t.3 <= 0.05]
  p.values.t.4 = vapply(1:ncol(resid), function(j) my.t.test(x=resid[-1,j]*resid[-nrow(resid),j], alternative="greater"), c(p.value=0))
  affdYears2 = (1:ncol(resid))[p.values.t.4 <= 0.05]
  affdYears = sort(union(affdYears1, affdYears2))

  return(list(affdYears = affdYears, affdDiagonals = affdDiagonals))
}


#' Smoothes demographic data using automatically estimated parameters and
#' taking into account only significant period and cohort effects.
#'
#' It is a heuristic procedure which tries to figure out positions of
#' period and cohort effects in the data. It also uses a few steps to estimate
#' model's parameters. The procedure is supposed to outperform \code{\link{autoDemogSmooth}} slightly.
#'
#' @param data Demographic data presented as a matrix.
#' @param p.value P-value used to test effects for significance.
#' The lower the value the less diagonals and years will be used to find cohort and period effects.
#' @param cornerLength Sets the smallest length of a digonal to be considered for cohort effects.
#' @param lower Lowest possible values for the optimisation procedure.
#' @param upper Highest possible values for the optimisation procedure.
#' @param init Initial values for the optimisation procedure.
#' @param reltol Relative tolerance parameter to be supplied to optim function (standard optimiser for R).
#' @param trace Controls if tracing is on.
#' @param control The control data passed directly to rq.fit.sfn method (quantreg package).
#' @return A list with six components: smooth surface, period effects, cohort effects, parameters
#' used for smoothing, diagonals used for cohort effects and years used for period effects.
#' @examples
#' # library(demography)
#' m = log(fr.mort$rate$female[1:30, 120:139])
#' # Show(m)
#' # sm = twoStepDemogSmooth(m)
#' # plot(sm)
#' # plot(sm, "surface")
#' # plot(sm, "period")
#' # plot(sm, "cohort")
#' @references \url{https://business.monash.edu/econometrics-and-business-statistics/research/publications/ebs/two-dimensional_smoothing_of_mortality_rates..pdf}
#' @author Alex Dokumentov
#' @export

twoStepDemogSmooth = function(data,
                              p.value = 0.05,
                              cornerLength = 7,
                              lower = c(0.01, 0.01, 0.01, 1.0, 0.001, 1.0, 0.001),
                              upper = c(1.2,  1.8,  1.2,  12,  0.4,  12,  0.4),
                              init =  c(0.1,  0.1,  0.2,  4,   0.001, 4,   0.001),
                              reltol = 0.001,
                              trace = F,
                              control = list(nnzlmax=1000000, nsubmax = 2000000, tmpmax = 200000))
{
  lambdayy = estYY(data,
                   lower = lower[2],
                   upper = upper[2],
                   step = abs(upper[2]-lower[2])/20,
                   trace = trace,
                   control = control)
  resid = smoothCv(demogSmooth.wrapper,
                   data = data,
                   lambda = 1,
                   lambdaaa = 0,
                   lambdayy = lambdayy,
                   lambdaay = 0,
                   effects = FALSE,
                   control = control)$cvResiduals
  lambdaYearsEffect = estAA(resid,
                           lower = lower[4],
                           upper = upper[4],
                           step = abs(upper[4]-lower[4])/20,
                           trace = trace,
                           control = control)
  result1 = demogSmooth(resid,
                       lambda = 1,
                       lambdaaa = lambdaYearsEffect,
                       lambdayy = 0,
                       lambdaay = 0,
                       effects = F,
                       control = control)
  vals = colSums(abs(result1$result))
  dataNA = data
  colsNA = NULL
  k = min(ceiling(length(vals)*0.15) , sum(vals > 2*mean(vals)))
  if(k > 0) {
    n <- length(vals)
    colsNA = which(vals >= sort(vals, partial = n-k+1)[n-k+1])
    dataNA[,colsNA] = NA
  }
  parametersNA = estPar(dataNA,
                        effects = FALSE,
                        lower = lower[1:3],
                        upper = upper[1:3],
                        init = init[1:3],
                        reltol = reltol,
                        trace = trace,
                        control = control)$par
  resid2 = smoothCv(demogSmooth.wrapper,
                   data = data,
                   lambda = 1,
                   lambdaaa = parametersNA[1],
                   lambdayy = parametersNA[2],
                   lambdaay = parametersNA[3],
                   effects = FALSE,
                   control = control)$cvResiduals
  affd = getAffected(resid2, p.value = p.value)
  affd$affdYears = sort(union(affd$affdYears, colsNA))
  # Estimating also taking into account period and cohort effects
  parameters = estPar(data,
                      effects = TRUE,
                      affdDiagonals = affd$affdDiagonals,
                      affdYears = affd$affdYears,
                      parameters = c(parametersNA,lambdaYearsEffect,init[5],lambdaYearsEffect,init[7]),
                      cornerLength = cornerLength,
                      lower = lower,
                      upper = upper,
                      init = init,
                      reltol = reltol,
                      trace = trace,
                      control = control)$par
  result = demogSmooth(data,
                       lambda = 1,
                       lambdaaa = parameters[1],
                       lambdayy = parameters[2],
                       lambdaay = parameters[3],
                       lambdaYearsEffect = parameters[4],
                       thetaYearsEffect = parameters[5],
                       lambdaCohortEffect = parameters[6],
                       thetaCohortEffect = parameters[7],
                       cornerLength = cornerLength,
                       effects = TRUE,
                       affdDiagonals = affd$affdDiagonals,
                       affdYears = affd$affdYears,
                       control = control)
  residuals = result$original - result$result - result$yearsEffect - result$cohortEffect
  # Removing very small period effects:
  affd$affdYears = which(colSums(result$yearsEffect) > 2 * mean(abs(residuals)))
  # Reestimating
  parameters = estPar(data,
                      effects = TRUE,
                      affdDiagonals = affd$affdDiagonals,
                      affdYears = affd$affdYears,
                      parameters = parameters,
                      cornerLength = cornerLength,
                      lower = lower,
                      upper = upper,
                      init = init,
                      reltol = reltol,
                      trace = trace,
                      control = control)$par
  result = demogSmooth(data,
                       lambda = 1,
                       lambdaaa = parameters[1],
                       lambdayy = parameters[2],
                       lambdaay = parameters[3],
                       lambdaYearsEffect = parameters[4],
                       thetaYearsEffect = parameters[5],
                       lambdaCohortEffect = parameters[6],
                       thetaCohortEffect = parameters[7],
                       cornerLength = cornerLength,
                       effects = TRUE,
                       affdDiagonals = affd$affdDiagonals,
                       affdYears = affd$affdYears,
                       control = control)
  result$parameters = parameters
  result$affdDiagonals = affd$affdDiagonals
  result$affdYears = affd$affdYears
  return(result)
}

