require(compiler)

cvPoints.nc = function(k, mini, maxi, maxj)
{
  k = (k+2) %% 5
  result = list()
  n = 1
  for(i in mini:maxi)
    for(j in 1:maxj)
      if((k + i + 2*j) %% 5 == 0)
      {
        result[[n]] = c(i,j)
        n = n + 1
      }
  return(result)
}

cvPoints = cmpfun(cvPoints.nc)

partialSSE.nc = function(data1, data2, points)
{
  SSE = 0
  for(p in points) {
    SSE = SSE + (data1[p[1],p[2]] - data2[p[1],p[2]])^2
  }
  return(SSE)
}

partialSSE = cmpfun(partialSSE.nc)

partialSAE.nc = function(data1, data2, points)
{
  SAE = 0
  for(p in points) {
    SAE = SAE + abs(data1[p[1],p[2]] - data2[p[1],p[2]])
  }
  return(SAE)
}

partialSAE = cmpfun(partialSAE.nc)

smoothCv.nc = function(smoothFun, data, upToAgeInd = dim(data)[1], fromAgeInd = 1, folds = 0:4, ...)
{
  nAges = dim(data)[1]
  nYears = dim(data)[2]
  SSE = 0
  SAE = 0
  l = 0
  for(k in folds) {
    lmWithNAs = data
    naPoints = cvPoints(k, max(1, fromAgeInd), min(nAges, upToAgeInd), nYears)
    fltPoints = list()
    n = 1
    for(i in naPoints) {
      if(!is.na(lmWithNAs[i[1],i[2]])) {
        lmWithNAs[i[1],i[2]] = NA
        fltPoints[[n]] = i
        n = n + 1
      }
    }
    result = smoothFun(lmWithNAs, ...)
    SSE = SSE + partialSSE(data, result, fltPoints)
    SAE = SAE + partialSAE(data, result, fltPoints)
    l = l + length(fltPoints)
  }
  return(c(MSE = SSE/l, MAE = SAE/l))
}

smoothCv = cmpfun(smoothCv.nc)
