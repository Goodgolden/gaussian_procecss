# Create simulated progesterone curves
#
# Only one year measured on each subject
#
# in addition to cycle day, the outcome (progesterone) is dependent
# 	on known covariates: age, time-to-menopause (ttm)
#
# Written by Nichole Carlson
#
# A general formula is at end of code


dir <- "C:/Users/mawhinns/Dropbox/students/randy_debashis/simulations/"
file.run <- paste(dir, "prog_sim.R", sep = "")


### Create age range of observations -------------------------------------------
age <- seq(45, 55, by = 1)

### Create time to menopause for each age --------------------------------------
Nsub <- NULL
obs.ttm <- NULL
menoage <- seq(44, 60, by = 1)

MeanMenoAge <- 51
StdMenoAge <- 4.5
menodist <- dnorm(menoage, MeanMenoAge, StdMenoAge)
cummenodist <- 1 - pnorm(menoage, MeanMenoAge, StdMenoAge)

#### Even age distribution
Nagegrp <- round(900 / length(age), 0)

### Create the time to menopause for each individual
for (i in 1:length(age)) {
  menoprob <- (pnorm(
    menoage[i:length(menoage)] + 1,
    MeanMenoAge, StdMenoAge
  ) -
    pnorm(
      menoage[i:length(menoage)],
      MeanMenoAge, StdMenoAge
    )) / cummenodist[i]
  ttm.1 <- c(seq(0, 5, by = 1), rep(5, 12 - i))
  menoinfo <- cbind(menoage[i:length(menoage)] + 1, menoprob, ttm.1)
  probmeno <- by(menoprob, ttm.1, sum)
  probmeno <- as.vector(probmeno)
  Nmeno <- round(probmeno * Nagegrp)
  Nsub <- c(Nsub, sum(Nmeno))
  obs.ttm <- c(obs.ttm, rep(c(0:5), Nmeno))
}

### Create age distribution at cross-sectional look of data
ageall <- rep(age, Nsub)

#### Create length(day) time series for each subject

### Make random like ttm code above
CycleLength <- 28
Day <- seq(1, CycleLength, by = 1.0) ### should modify to be random for each subject

LutealLength <- 14 ### Make random like ttm code above####

ind <- c(rep(0, CycleLength - LutealLength), rep(1, LutealLength)) ### Indictor whether in the Luteal Phase

### Set parameters in the model
Beta0 <- 0.5
CurveAgeInt <- 23.5
CurveAgeSlope <- -0.20
CurveTTMSlope <- 2.4
REStd <- 0.2

PeakLoc <- CycleLength - LutealLength / 2
Denom <- (LutealLength / 2)^2
CurveAge <- CurveAgeInt + CurveAgeSlope * ageall
CurveTTM <- CurveTTMSlope * (obs.ttm - 5)
CurveCoeff <- CurveAge + CurveTTM
CurveCoeffRE <- rnorm(length(ageall), CurveCoeff, REStd)

### Concentration profiles -----------------------------------------------------
prog <- NULL
lprogerr <- NULL
for (i in 1:length(ageall)) {
  subj.prog <- Beta0 +
    (-CurveCoeffRE[i] / Denom *
      (Day - PeakLoc)^2 +
      CurveCoeffRE[i]) * ind

  sub.lprogerr <- log(subj.prog) +
    rnorm(length(subj.prog), 0, sqrt(0.02))
  prog <- c(prog, subj.prog)
  lprogerr <- c(lprogerr, sub.lprogerr)
}

#### Get baseline age and time to menopause repeated
#### for stacked data formulation
daylength <- rep(28, length(ageall))
agestacked <- rep(ageall, daylength)
ttmstacked <- rep(obs.ttm, daylength)
daystacked <- rep(Day, length(ageall))
CurveAgeStacked <- rep(CurveAge, daylength)
CurveTTMStacked <- rep(CurveTTM, daylength)
CurveCoeffStacked <- rep(CurveCoeff, daylength)
CurveCoeffStackedRE <- rep(CurveCoeffRE, daylength)
Id <- seq(1, length(ageall), by = 1)
SubjectId <- rep(Id, daylength)

Simdata <-
  cbind(
    SubjectId,
    agestacked,
    ttmstacked,
    daystacked,
    prog, lprogerr,
    CurveCoeffStacked,
    CurveCoeffStackedRE,
    CurveAgeStacked,
    CurveTTMStacked
  ) %>%
  as.data.frame()

write.csv((Simdata),
  file = paste(dir, "progsim.csv", sep = "")
)
write.csv((Simdata[Simdata$SubjectId == 1, ]),
  file = paste(dir, "progsim_one.csv", sep = "")
)

plot(Simdata$daystacked,
  Simdata$prog,
  col = 2
)
lines(Simdata$daystacked,
  exp(Simdata$lprogerr),
  type = "o"
)

plot(
  Simdata$daystacked[Simdata$SubjectId == 1],
  Simdata$prog[Simdata$SubjectId == 1]
)

#### general formula, linear change with ttm, peak at day 21
prog <- 1.0 + 0.1 * (age - 55) +
  (-2.4 * ttm / 49 * (Day - 21)^2 + 2.4 * ttm) * ind
lprogerr <- log(prog) + rnorm(28, 0, sqrt(0.02))
