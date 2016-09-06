library(demography)
m = log(fr.mort$rate$female[1:30, 150:160])
sm = autoDemogSmooth(m)

plot(sm)
plot(sm, "surface", types = "2D")
plot(sm, "cohort")
plot(sm, "period")
plot(sm, "residuals")
plot(sm, "original")

m = matrix(rnorm(225), 15, 15)
plot(m, types = "2D")

plot(matrix(1:100,10,10), c("Dimension 1", "Dimension 2", "Value"), type = "2D")
