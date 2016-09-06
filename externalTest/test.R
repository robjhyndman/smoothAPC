library(demography)
m = log(fr.mort$rate$female[1:30, 150:160])
sm = autoDemogSmooth(m)

plot(sm)
plot(sm, "cohort")
plot(sm, "period")
plot(sm, "residuals")
plot(sm, "original")

