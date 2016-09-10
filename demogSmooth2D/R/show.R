#' @importFrom rgl open3d
#' @importFrom rgl lines3d
#' @importFrom rgl rgl.bringtotop
#' @importFrom rgl persp3d
#' @importFrom rgl title3d
#' @importFrom colorspace hex
#' @importFrom colorspace polarLUV
#' @importFrom grDevices rainbow
#' @importFrom graphics axis
#' @importFrom graphics filled.contour
#' @importFrom graphics title
#' @importFrom methods new
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom utils head


plot3D = function(x, y, z,
                title="",
                labs = c("Age", "Time", "Value"),
                aspect = c(1, 1, 0.6),
                grid,
                gridCol = "red",
                type = 2,
                color.pallete = function(n) rainbow(n, start=0.0, end=0.7))
{
  aspectX = x[length(x)] - x[1]
  aspectY = y[length(y)] - y[1]
  mn = min(z, na.rm=TRUE)
  mx = max(z, na.rm=TRUE)
  if(type == 1)
    c = color.pallete(dim(z)[1]*dim(z)[2])
  else if(type == 2)
  {
    if(mx == mn) {
      c = "green"
    }
    else {
      c = (z - mn)/(mx - mn)
      r = color.pallete(65536)
      c = r[round(c * 65535) + 1]
    }
  }
  open3d(windowRect=c(10,35,810,835))
  persp3d(x = x, y = y, z = z, aspect = c(aspectX * aspect[1], aspectY * aspect[2], aspectX * aspect[3]), xlab = labs[1], ylab = labs[2], zlab = labs[3], col=c)
  title3d(main=title)

  eps1 = 0.001 * abs(mx - mn)
  if(!missing(grid) && !is.na(grid[1])) {
    xGrid = seq(1, length(x), grid[1])
    for(xInd in xGrid) {
      lines3d(x[xInd], y, z[xInd,] + eps1, col = gridCol)
      lines3d(x[xInd], y, z[xInd,] - eps1, col = gridCol)
    }
  }
  if(!missing(grid) && !is.na(grid[2])) {
    yGrid = seq(1, length(y), grid[2])
    for(yInd in yGrid) {
      lines3d(x, y[yInd], z[,yInd] + eps1, col = gridCol)
      lines3d(x, y[yInd], z[,yInd] - eps1, col = gridCol)
    }
  }
  rgl.bringtotop()
}

#' Presents demographic data using 3D surface and/or a heatmap.
#'
#' @param x Result of smoothing (object of class \code{sm2D}).
#' @param component "smooth", "period", "cohort", "residuals" or "original"
#' @param labs Vector of labels for X, Y and Z axes.
#' @param types Vector of plot types to plot. Possible types are \code{"2D"} and \code{"3D"}. Default value for the parameter is \code{c("3D", "2D")}.
#' @param ... Other parameters. They are currently ignored.
#' @examples
#' \dontrun{
#'
#' library(demography)
#' m = log(fr.mort$rate$female[1:30, 150:160])
#' sm = autoDemogSmooth(m)
#'
#' plot(sm)
#' plot(sm, "surface")
#' plot(sm, "cohort")
#' plot(sm, "period")
#' plot(sm, "residuals")
#' plot(sm, "original")
#'
#' }
#' @author Alexander Dokumentov
#' @export

plot.sm2D = function(x,
                     component = c("all", "surface", "period", "cohort", "residuals", "original"),
                     labs = c("Age", "Time", NA),
                     types = c("3D", "2D"),
                     ...)
{
  if(!(component[1] %in% c("all", "surface", "period", "cohort", "residuals", "original")))
    stop("Incorrect component.")
  data = switch (which(component[1] == c("all","surface", "period", "cohort", "residuals", "original")),
    {
      combined = x$result
      if(!is.null(x$yearsEffect)) combined = combined + x$yearsEffect
      if(!is.null(x$cohortEffect)) combined = combined + x$cohortEffect
      combined
    },
    x$result,
    x$yearsEffect,
    x$cohortEffect,
    {
      residuals = x$original - x$result
      if(!is.null(x$yearsEffect)) residuals = residuals - x$yearsEffect
      if(!is.null(x$cohortEffect)) residuals = residuals - x$cohortEffect
      residuals
    },
    x$original
  )
  if(is.null(data)) stop(paste0('Component "', component[1], '" cannot be extracted.'))
  labs[3] = ifelse(is.na(labs[3]), gsub("(^[[:alpha:]])", "\\U\\1", component[1], perl=TRUE), labs[3])
  plot.matrix(x = data, labs = labs, types = types)
}

#' Presents matrix as 3D surface and/or a heatmap.
#'
#' @param x Matrix to plot.
#' @param labs Vector of lables for X, Y and Z axes.
#' @param types Vector of plot types to plot. Possible types are \code{"2D"} and \code{"3D"}. Default value for the parameter is \code{c("3D", "2D")}.
#' @param ... Other parameters. They are currently ignored.
#' @examples
#' \dontrun{
#'
#' plot(matrix(rnorm(100),10,10))
#' plot(matrix(1:100,10,10), c("Dimension 1", "Dimension 2", "Value"))
#'
#' }
#' @author Alexander Dokumentov
#' @export

plot.matrix = function(x, labs = c("X", "Y", "Z"), types = c("2D", "3D"), ...)
{
  if("3D" %in% types) plot3D(1:dim(x)[1], 1:dim(x)[2], x, labs = labs)
  if("2D" %in% types) plot2D(1:dim(x)[1], 1:dim(x)[2], x, labs = labs)
}

my.colors =
  function(n, h = c(260, 0), c = 80, l = c(20, 90), power = .7,
          fixup = TRUE, gamma = NULL, ...)
{
  if(!is.null(gamma))
    warning("'gamma' is deprecated and has no effect")
  if(n < 1)
    return(character(0))
  h <- rep(h, length.out = 2)
  c <- c[1]
  l <- rep(l, length.out = 2)
  power <- rep(power, length.out = 2)
  rval <- seq(1, -1, length = n)
  rval <- hex(polarLUV(L = l[2] - diff(l) * abs(rval)^power[2],
                       C = c * abs(rval)^power[1], H = ifelse(rval > 0, h[1],
                                                              h[2])), fixup = fixup, ...)
  return(rval)
}

plot2D = function(ages, years, z, labs=c("X", "Y", "Z"))
{
  if(max(z)>0 && min(z)<0) {
    filled.contour(ages, years, z, zlim = c(-max(abs(z)), max(abs(z))), color.palette = my.colors,
                   plot.title = title(main = labs[3], xlab = labs[1], ylab = labs[2]),
                   plot.axes = {axis(1, seq(ages[1],ages[length(ages)],10)); axis(2, years[seq(1,100,5)])})
  } else {
    filled.contour(ages, years, z, color.palette = function(n) rainbow(n, start=0.0, end=0.7),
                   plot.title = title(main = labs[3], xlab = labs[1], ylab = labs[2]),
                   plot.axes = {axis(1, seq(ages[1],ages[length(ages)],10)); axis(2, years[seq(1,100,5)])})
  }
}
