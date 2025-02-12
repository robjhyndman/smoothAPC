#' @importFrom rgl open3d
#' @importFrom rgl lines3d
#' @importFrom rgl rgl.bringtotop
#' @importFrom rgl persp3d
#' @importFrom rgl title3d
#' @importFrom rgl axis3d
#' @importFrom colorspace hex
#' @importFrom colorspace polarLUV
#' @importFrom colorspace heat_hcl
#' @importFrom grDevices rainbow
#' @importFrom graphics axis
#' @importFrom graphics filled.contour
#' @importFrom graphics title
#' @importFrom methods new
#' @importFrom stats na.omit
#' @importFrom stats optim
#' @importFrom utils head

# Base function for 2d and 3d plotting of matrices

plot_matrix <- function(x, labs = c("X", "Y", "Z"), types = c("2D", "3D"), color.palette = c("default", "special"), ...) {
  if (max(x) > 0 && min(x) < 0) {
    zlim <- c(-max(abs(x)), max(abs(x)))
    if (inherits(color.palette, "function")) {
      palette <- color.palette
    } else {
      palette <- my.colors
    }
  } else {
    zlim <- c(min(x), max(x))
    palette <- function(n) heat_hcl(n)
    if (inherits(color.palette, "character") && color.palette[1] == "special") {
      palette <- function(n) rainbow(n, start = 0.0, end = 0.7)
    } else if (inherits(color.palette, "function")) {
      palette <- color.palette
    }
  }
  ages <- rownames(x)
  if (is.null(ages)) ages <- 1:dim(x)[1]
  years <- colnames(x)
  if (is.null(years)) years <- 1:dim(x)[2]
  if ("3D" %in% types) plot3D(ages, years, x, labs = labs, color.palette = palette)
  if ("2D" %in% types) {
    filled.contour(seq_along(ages), seq_along(years), x,
      zlim = zlim, color.palette = palette,
      plot.title = title(main = labs[3], xlab = labs[1], ylab = labs[2]),
      plot.axes = {
        axis(side = 1, at = seq(1, length(ages), 5), labels = ages[seq(1, length(ages), 5)])
        axis(side = 2, at = seq(1, length(years), 5), labels = years[seq(1, length(years), 5)])
      }
    )
  }
}

# Base function for 2d and 3d plotting of smAPC objects

plot_smAPC <- function(
    x,
    component = c("all", "surface", "period", "cohort", "residuals", "original"),
    labs = c("Age", "Time", NA),
    types = c("3D", "2D"),
    color.palette = c("default", "special"),
    ...) {
  if (!(component[1] %in% c("all", "surface", "period", "cohort", "residuals", "original"))) {
    stop("Incorrect component.")
  }
  data <- switch(which(component[1] == c("all", "surface", "period", "cohort", "residuals", "original")),
    {
      combined <- x$result
      if (!is.null(x$yearsEffect)) combined <- combined + x$yearsEffect
      if (!is.null(x$cohortEffect)) combined <- combined + x$cohortEffect
      combined
    },
    x$result,
    x$yearsEffect,
    x$cohortEffect,
    {
      residuals <- x$original - x$result
      if (!is.null(x$yearsEffect)) residuals <- residuals - x$yearsEffect
      if (!is.null(x$cohortEffect)) residuals <- residuals - x$cohortEffect
      residuals
    },
    x$original
  )
  if (is.null(data)) stop(paste0('Component "', component[1], '" cannot be extracted.'))
  labs[3] <- ifelse(is.na(labs[3]), gsub("(^[[:alpha:]])", "\\U\\1", component[1], perl = TRUE), labs[3])
  plot_matrix(x = data, labs = labs, types = types, color.palette = color.palette)
}

# Default color palette for heat maps
my.colors <-
  function(n, h = c(260, 0), c = 80, l = c(20, 90), power = .7,
           fixup = TRUE, gamma = NULL, ...) {
    if (!is.null(gamma)) {
      warning("'gamma' is deprecated and has no effect")
    }
    if (n < 1) {
      return(character(0))
    }
    h <- rep(h, length.out = 2)
    c <- c[1]
    l <- rep(l, length.out = 2)
    power <- rep(power, length.out = 2)
    rval <- seq(1, -1, length = n)
    rval <- hex(polarLUV(
      L = l[2] - diff(l) * abs(rval)^power[2],
      C = c * abs(rval)^power[1], H = ifelse(rval > 0, h[1],
        h[2]
      )
    ), fixup = fixup, ...)
    return(rval)
  }

# 2d plotting functions --------------------------------------------------------

#' Methods for producing heatmap plots for matrix and smAPC objects.
#'
#' Displays a matrix or an estimated demographic model of class \code{smAPC} as a heatmap.
#'
#' @param x Either a matrix or the result of fitting a smoothed demographic model (object of class \code{smAPC}).
#' @param labs Vector of labels for X and Y axes.
#' @param color.palette Character string \code{"default"} or \code{"special"}
#' or a function accepting one argument and returning a color palette
#' (for example \code{\link[grDevices]{rainbow}}).
#' @param main Title for the plot.
#' @param ... Other parameters. They are currently ignored.
#' @examples
#' plot(matrix(rnorm(100), 10, 10), main = "Noise")
#' plot(matrix(1:100, 10, 10), c("Dimension 1", "Dimension 2"), main = "Value")
#'
#' m <- log(demography::fr.mort$rate$female[1:30, 150:160])
#' plot(m)
#' plot(m, color.palette = "special")
#' plot(m, color.palette = rainbow)
#' @author Alexander Dokumentov
#' @export

plot.matrix <- function(x, labs = c("X", "Y"), color.palette = c("default", "special"), main = "", ...) {
  plot_matrix(x = x, labs = c(labs, main), types = "2D", color.palette = color.palette, ...)
}

#' @param component Character string specifying which component to plot.
#' @rdname plot.matrix
#' @examples
#'
#' \dontrun{
#' m <- log(demography::fr.mort$rate$female[1:30, 150:160])
#' sm <- autoSmoothAPC(m)
#'
#' plot(sm)
#' plot(sm, "surface")
#' plot(sm, "cohort")
#' plot(sm, "period")
#' plot(sm, "residuals")
#' plot(sm, "original", main = "Original data")
#' }
#' @author Alexander Dokumentov
#' @export

plot.smAPC <- function(
    x,
    component = c("all", "surface", "period", "cohort", "residuals", "original"),
    labs = c("Age", "Time"),
    color.palette = c("default", "special"),
    main = "",
    ...) {
  plot_smAPC(x = x, component = component, labs = c(labs, main), types = "2D", color.palette = color.palette, ...)
}

# 3d plotting functions --------------------------------------------------------

plot3D <- function(
    x, y, z,
    title = "",
    labs = c("Age", "Time", "Value"),
    aspect = c(1, 1, 0.6),
    color.palette = function(n) rainbow(n, start = 0.0, end = 0.7)) {
  aspectX <- length(x)
  aspectY <- length(y)
  mn <- min(z, na.rm = TRUE)
  mx <- max(z, na.rm = TRUE)
  if (mx == mn) {
    c <- "green"
  } else {
    c <- (z - mn) / (mx - mn)
    r <- color.palette(65536)
    c <- r[round(c * 65535) + 1]
  }
  open3d(windowRect = c(10, 35, 810, 835))
  persp3d(x = seq_along(x), y = seq_along(y), z = z, aspect = c(aspectX * aspect[1], aspectY * aspect[2], aspectX * aspect[3]), xlab = labs[1], ylab = labs[2], zlab = labs[3], col = c, axes = FALSE)
  axis3d(edge = "x", at = seq(1, length(x), 5), labels = x[seq(1, length(x), 5)])
  axis3d(edge = "y", at = seq(1, length(y), 5), labels = y[seq(1, length(y), 5)])
  axis3d(edge = "z")
  title3d(main = title)
  rgl.bringtotop()
}

#' Methods for producing 3d surface plots
#'
#' Plots a matrix or an estimated demographic model of class \code{smAPC} as a 3D surface.
#'
#' @param x Either a matrix or the result of fitting a smoothed demographic model (object of class \code{smAPC}).
#' @param labs Vector of labels for X, Y and Z axes.
#' @param color.palette Character string \code{"default"} or \code{"special"}
#' or a function accepting one argument and returning a color palette
#' (for example \code{\link[grDevices]{rainbow}}).
#' @param ... Other parameters. They are currently ignored.
#'
#' @export

plot3d <- function(x, ...) UseMethod("plot3d")

#' @rdname plot3d
#' @examples
#' \dontrun{
#' plot3d(matrix(rnorm(100), 10, 10))
#' plot3d(matrix(1:100, 10, 10), c("Dimension 1", "Dimension 2", "Value"))
#'
#' m <- log(demography::fr.mort$rate$female[1:30, 150:160])
#' plot3d(m)
#' plot3d(m, color.palette = "special")
#' plot3d(m, color.palette = rainbow)
#' }
#' @author Alexander Dokumentov
#' @export

plot3d.matrix <- function(x, labs = c("X", "Y", "Z"), color.palette = c("default", "special"), ...) {
  plot_matrix(x = x, labs = labs, types = "3D", color.palette = color.palette, ...)
}

#' @param component Character string specifying which component to plot.
#' @rdname plot3d
#' @examples
#'
#' \dontrun{
#' m <- log(demography::fr.mort$rate$female[1:30, 150:160])
#' sm <- autoSmoothAPC(m)
#'
#' plot3d(sm)
#' plot3d(sm, "surface", color.palette = "special")
#' plot3d(sm, "cohort")
#' plot3d(sm, "period")
#' plot3d(sm, "residuals")
#' plot3d(sm, "original", color.palette = rainbow)
#' }
#' @author Alexander Dokumentov
#' @export

plot3d.smAPC <- function(
    x,
    component = c("all", "surface", "period", "cohort", "residuals", "original"),
    labs = c("Age", "Time", NA),
    color.palette = c("default", "special"),
    ...) {
  plot_smAPC(x = x, component = component, labs = labs, types = "3D", color.palette = color.palette, ...)
}
