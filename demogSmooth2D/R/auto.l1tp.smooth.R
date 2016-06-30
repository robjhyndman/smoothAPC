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
    penalty = sum(abs(xx - x))
    x = xx
    if(trace) print(x)
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
      print(cv[2])
      print("Penalty:")
      print(penalty)
    }
    return(cv[2] + penalty)
  }

  if(!is.null(parameters)) {
    for(i in 1L:min(length(parameters), length(init))) {
      init[i] = min(max(parameters[i], lower[i]), upper[i])
    }
  }

  result = optim(par = init, fn = f, control = list(reltol = reltol))
  return(result$par)
}

estYY = function(data,
                 lower = 0.2,
                 upper = 10,
                 step = 0.2,
                 control = list(nnzlmax = 1000000, nsubmax = 2000000, tmpmax = 200000))
{
  lambdas = seq(lower, upper, by = step)
  cv = 0
  for(i in seq_along(lambdas)) {
    cat("\r "); cat(paste0(c("\\","|","/","-")[i %% 4 + 1], "   "))
    cv[i] <- smoothCv(demogSmooth.wrapper, data = data,
                      lambda = 1, lambdaaa = 0.01*lambdas[i], lambdayy = lambdas[i], lambdaay = 0.01*lambdas[i],
                      effects = FALSE, control = control)[2]
  }
  return(lambdas[which.min(cv)])
}


#' Smoothes demographic data using automatically estimated parameters and optionally
#' taking into account period and cohort effects
#'
#' If period and cohort effects are taken into account (effects = TRUE) the method uses all
#' available years and diagonals for estimation of the period and cohort effects.
#'
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
#' # sm = autoDemogSmooth(m)
#' # Show(m)
#' # Show(sm$result)
#' # Show(sm$yearsEffect)
#' # Show(sm$cohortEffect)
#' @export

autoDemogSmooth = function(data,
                           effects = TRUE,
                           cornerLength = 7,
                           affdDiagonals = NULL,
                           affdYears = NULL,
                           lower = head(c(0.01, 0.01, 0.01, 2.0, 0.01, 2.0, 0.01), 3 + effects*4),
                           upper = head(c(1.2,  1.8,  1.2,  12,  0.4,  12,  0.4), 3 + effects*4),
                           init =  head(c(0.1,  0.1,  0.2,  4,   0.05, 4,   0.05), 3 + effects*4),
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
                        control = control)
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
  result[[4]] = parameters
  names(result)[4] = "parameters"
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
#' taking into account only significant period and cohort effects
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
#' # m = log(fr.mort$rate$female[1:30, 150:160])
#' # sm = twoStepDemogSmooth(m)
#' # Show(m)
#' # Show(sm$result)
#' # Show(sm$yearsEffect)
#' # Show(sm$cohortEffect)
#' @export

twoStepDemogSmooth = function(data,
                              p.value = 0.05,
                              cornerLength = 7,
                              lower = c(0.01, 0.01, 0.01, 2.0, 0.01, 2.0, 0.01),
                              upper = c(1.2,  1.8,  1.2,  12,  0.4,  12,  0.4),
                              init =  c(0.1,  0.1,  0.2,  4,   0.05, 4,   0.05),
                              reltol = 0.001,
                              trace = F,
                              control = list(nnzlmax=1000000, nsubmax = 2000000, tmpmax = 200000))
{
  lambdayy = estYY(data,
                   lower = lower[2],
                   upper = upper[2],
                   step = abs(upper[2]-lower[2])/20,
                   control = control)
  result0 = demogSmooth(data,
                        lambda = 1,
                        lambdaaa = 0,
                        lambdayy = lambdayy,
                        lambdaay = 0,
                        effects = FALSE,
                        control = control)
  resid = data - result0$result
  affd = getAffected(resid, p.value = p.value)

  parameters = estPar(data,
                      effects = TRUE,
                      affdDiagonals = affd$affdDiagonals,
                      affdYears = affd$affdYears,
                      parameters = c(0.1*lambdayy, 0.1*lambdayy, lambdayy),
                      cornerLength = cornerLength,
                      lower = lower,
                      upper = upper,
                      init =  init,
                      reltol = reltol,
                      trace = trace,
                      control = control)
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
  result[[4]] = parameters
  names(result)[4] = "parameters"
  result[[5]] = affd$affdDiagonals
  names(result)[5] = "affdDiagonals"
  result[[6]] = affd$affdYears
  names(result)[6] = "affdYears"
  return(result)
}
