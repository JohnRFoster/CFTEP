#--------------------------------
#
# Fit Occupancy models and generate various outputs
#
# By: Ryan S. Miller
#
# Fit Basic Occupancy models using covariates for detection and occupancy
#
# Last Updated: 3 Nov 2020
#
#----------------------------------

#--Clean Workspace
#rm(list=ls(all=TRUE))
#gc()

#---- Load some useful aliases and functions
# NOTE: requires an active internet connection
# install.packages("devtools") # Install your first time through
require(devtools)
source_gist(9216051) # Rprofile.R
source_gist(9216061) # various.R

# Loading required packages
toLoad = c(
  "unmarked",
  "lubridate",
  "plyr",
  "dplyr",
  "AICcmodavg",
  "car",
  "ggplot2",
  "grid"
)
instant_pkgs(toLoad)

# Load a modification of an`unmarked` function and custom plotting function

code.path <- "C:/Documents/Project Documents/CattleFeverTick/Code/"

source(paste0(code.path, "plotOccu.R"))
source(paste0(code.path, "update_formatMult.R"))
source(paste0(code.path, "newmodSel.R"))


#----Load Packages for Session----
packages.vec <- c(
  "sp",
  "plyr",
  "raster",
  "wiqid",
  "unmarked",
  "AICcmodavg",
  "DMwR"
)

# Check if installed and install if not
if (
  FALSE %in% unique(is.element(packages.vec, rownames(installed.packages())))
) {
  install.packages(setdiff(packages.vec, rownames(installed.packages())))
}
# load
lapply(packages.vec, require, character.only = TRUE)
#----END Load libraries----

#---- Set Paths

model.dat.path <- "C:/Documents/Project Documents/CattleFeverTick/Model.Data/"
figure.path <- "C:/Documents/Project Documents/CattleFeverTick/Figures/"

#---- END Paths

#----READ DATA----

#--Read Occurance Data
#occ.raw<-read.csv(paste0(model.dat.path,"dat.occ.cov.2020-11-02.csv"))
#summary(occ.raw)
#nrow(occ.raw)

occ.dat <- csvToUMF(
  paste0(model.dat.path, "dat.occ.cov.2020-11-03.csv"),
  long = TRUE,
  type = "unmarkedFrameOccu"
)
summary(occ.dat)

#-- Scale and center predictors

#Observations
#obsCovs(occ.dat)$month <- scale(obsCovs(occ.dat)$month, center = TRUE, scale = TRUE)
obsCovs(occ.dat)$qty_inspected.std <- scale(
  obsCovs(occ.dat)$qty_inspected,
  center = TRUE,
  scale = TRUE
)
obsCovs(occ.dat)$qty_infested.std <- scale(
  obsCovs(occ.dat)$qty_infested,
  center = TRUE,
  scale = TRUE
)
obsCovs(occ.dat)$density.std <- scale(
  obsCovs(occ.dat)$density,
  center = TRUE,
  scale = TRUE
)

#Site
#siteCovs(occ.dat)$qty_added <- scale(siteCovs(occ.dat)$qty_added, center = TRUE, scale = TRUE)
#siteCovs(occ.dat)$pasture_size <- scale(siteCovs(occ.dat)$pasture_size, center = TRUE, scale = TRUE)

#---- END Read Data ----

#---- Fit Occupancy Models from MacKenzie et. al (2003) ----

#-- Naive occupancy estimate
# Numerator sums the number of sites with at least 1 detected BACS
naive.occ <- sum(apply(occ.dat@y, 1, sum) > 0) / nrow(occ.dat@y)

# Double right-hand side formula describing covariates of
# detection and occupancy in that order

#---- Intercept only model (null)
m0 <- occu(
  ~1 ~ 1,
  data = occ.dat,
  method = "BFGS",
  se = TRUE,
  control = list(trace = 1, maxit = 1e4)
)
summary(m0)

#--Back Transform Estimates
det.est0 <- backTransform(m0, "det")
psi.est0 <- backTransform(m0, "state")

confint(psi.est0)


#---- Detection Model
m.det <- occu(
  ~ species + qty_inspected.std + density.std ~ 1,
  data = occ.dat,
  method = "BFGS",
  se = TRUE,
  control = list(trace = 3, maxit = 1e4)
)
summary(m.det)


#---- Full Model - Detection and Occupancy predictors
m.full <- occu(
  ~ month + species + qty_inspected + density ~ 1,
  data = occ.dat,
  method = "BFGS",
  se = TRUE,
  control = list(trace = 3, maxit = 1e4)
)
summary(m.full)

#---- END Fit Occupancy Models ----

#---- Compare Models ----

# Model comparison
modList <- list(m0, m.det, m.full)
names(modList) <- c(
  "p(.)psi(.)",
  "p(month+species+qty_inspected+density)psi(.)",
  "p(month+species+qty_inspected+density)psi(pasture_qty_acres)"
)
fmList <- fitList(fits = modList)

# Model selection
ms <- newmodSel(fmList, IC = c("AICc"))
ms

#---- END Compare Models ----

#---- Evaluate Best Model ----

#Best model
mod <- m.det

# Goodness of fit test (MacKenzie & Bailey 2004)
gof_boot <- mb.gof.test(mod, nsim = 100, plot.hist = TRUE)
gof_boot

# values of c-hat > 1 indicate overdispersion (variance > mean),
# but � values much higher than 1 (i.e., > 4) probably indicate
# lack-of-fit. In cases of moderate overdispersion, one usually
# multiplies the variance-covariance matrix of the estimates by c-hat.
# As a result, the SE�s of the estimates are inflated (c-hat is also known
# as a variance inflation factor)

#The parametric bootstrap can be used to check the adequacy of model ???t. Here we use a ??2 statistic appropriate for binary data.

chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y > 1] <- 1
  sr <- fm@sitesRemoved
  if (length(sr) > 0) {
    y <- y[-sr, , drop = FALSE]
  }
  fv <- fitted(fm, na.rm = TRUE)
  y[is.na(fv)] <- NA
  sum((y - fv)^2 / (fv * (1 - fv)), na.rm = TRUE)
} #END Function

pb <- parboot(mod, statistic = chisq, nsim = 100, parallel = FALSE)

pb


# Function returning three fit-statistics.
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2)
  chisq <- sum((observed - expected)^2 / expected)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2)
  out <- c(SSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(out)
} #END Function

pb <- parboot(fm, fitstats, nsim = 25, report = 1)
plot(pb, main = "")

#---- END Model Evaluation ----

#---- Generate Model Outputs ----

#--Assign model to use
mod <- m.det

#---- Make Model Estimates ----

#--Expected probability that a site was occupied (population-level estimate)
# Coefficients: c(state, 0,...) 0 inidicates to use the mean value of the variable
pop.est <- backTransform(linearComb(mod, coefficients = c(1), type = 'state'))

#--Probability of detection given a site is occupied
# Coefficients: c(det, 0,...) 0 inidicates to use the mean value of the variable
n.var <- nrow(summary(mod)$det)
det.est <- backTransform(linearComb(
  mod,
  coefficients = c(1, rep(0, n.var - 1)),
  type = 'det'
))


#--Make Estimates

#--Detection estimate
estimate <- mod@estimates[2]@estimates
det.beta.mat <- cbind.data.frame(estimate, confint(mod, type = 'det'))

det.est <- cbind.data.frame(
  det.est = det.est@estimate,
  t(antilogit(confint(mod, type = 'det')[1, ]))
)

#--Get site level detection estimates
det.vec <- getP(mod)

#Make site detection means
det.site.mu <- rowMeans(det.vec, na.rm = TRUE)

#Generate bootstrap
det.stat <- bootstrap.CI(det.site.mu)

#--State estimate
estimate <- mod@estimates[1]@estimates
state.beta.mat <- cbind.data.frame(estimate, confint(mod, type = 'state'))

state.est <- cbind.data.frame(
  det.est = pop.est@estimate,
  t(antilogit(confint(mod, type = 'state')[1, ]))
)


#--Confidence intervals for projected occupancy (this takes awhile)
#mod <- nonparboot(mod, B=5)
#occ.est<-cbind.data.frame(projected=projected(mod)[2,], SE=mod@projected.mean.bsse[2,])

#---- End Make Model Estimates ----

#---- PLOTs ----

#---- Estimate Pasture Level Probability of Detection on Successive Surveys

p.star <- est.prob.detection.repeated.sampling(det.stat, n.surveys = 20)

pdf(
  paste0(figure.path, "FIG.prob.detection.many.surveys.", Sys.Date(), ".pdf"),
  width = 8,
  height = 8
)
boxplot(
  p.star$p.star ~ p.star$x,
  las = 1,
  ylab = "Probability of Detection (P star)",
  xlab = "Number of Surveys",
  outline = FALSE
)
abline(h = 0.95, lty = 2, lwd = 2)
dev.off()

#---- END

#---- Plot Detection Probability by Year

#--Generate Predictions for each site/observation
pred.det <- predict(mod, type = 'det')

#--Merge with observation data
pred.det <- cbind.data.frame(occ.dat@obsCovs, pred.det)

pdf(
  paste0(figure.path, "FIG.prob.detection.annual.", Sys.Date(), ".pdf"),
  width = 8,
  height = 8
)
boxplot(
  pred.det$Predicted ~ pred.det$year,
  las = 1,
  ylab = "Probability of Detection on Single Inspection",
  xlab = "",
  outline = FALSE
)
abline(h = 0.95, lty = 2, lwd = 2)
dev.off()

#---- END Plot Detection Probability by Year

#---- Detection Probability by Species

var.plt <- "qty_inspected.std"
var.mean <- "density.std"

#--Set Variable to hold at mean
x.mean.var <- mean(occ.dat@obsCovs[, var.mean], na.rm = TRUE)

#--Set variable of interest
x <- occ.dat@obsCovs[, var.plt]
x.var <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE) + 3, 1)

#--Unscale Values
center <- attr(x, "scaled:center")
scale <- attr(x, "scaled:scale")

unscale <- function(x) {
  t(t(x) * scale + center)
}


#--Make Plot

pdf(
  paste0(figure.path, "FIG.Det.Probability.by.species.num.inspected.pdf"),
  width = 8,
  height = 8
)

month.vec <- c(1, 3, 6, 9)
var.vec <- c("Bovine", "Whitetail", "Nilgai", "Equine", "Other Species")
lab.vec <- c("Bovine", "Whitetail", "Nilgai", "Equine", "Other")

for (i in 1:length(var.vec)) {
  newData <- data.frame(
    x.mean = x.mean.var,
    species = factor(
      var.vec[i],
      levels = c(
        "Bovine",
        "Equine",
        "None",
        "Whitetail",
        "Nilgai",
        "Other Species"
      )
    ),
    var.plt = x.var
  )
  colnames(newData)[1] <- var.mean
  colnames(newData)[3] <- var.plt

  tmp <- predict(mod, type = 'det', newdata = newData, appendData = TRUE)

  tmp$x <- unscale(tmp[, var.plt])

  if (i == 1) {
    plot(
      1,
      xlim = c(-100, max(tmp$x)),
      ylim = c(0, 1),
      type = "n",
      axes = FALSE,
      xlab = "Number of Animals Inspected",
      pch = 20,
      ylab = "Detection Probability",
      cex.lab = 1.25,
      cex.main = 1.75
    )
  }

  lines(tmp$x, tmp$Predicted, col = "black", lwd = 2)

  text(x = 0, y = tmp$Predicted[1], labels = lab.vec[i], pos = 2)
} #END Plot

axis(side = 1)
axis(side = 2, las = 2)

#x<-mean(occ.raw$qty_inspected,na.rm=TRUE)
#text(x=x, y=0,
#     label="Mean number of\nanimals inspected per visit",pos=4, col="red")

#abline(v=x,col="red")

mtext(
  side = 3,
  adj = 0,
  "a) Functional relationship between detection probability and animals inspected"
)

dev.off()

#---- END Detection Probability by Species

#---- Detection Probability by Species

var.plt <- "density.std"
var.mean <- "qty_inspected.std"

#--Set Variable to hold at mean
x.mean.var <- mean(occ.dat@obsCovs[, var.mean], na.rm = TRUE)

#--Set variable of interest
x <- occ.dat@obsCovs[, var.plt]
x.var <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE) + 3, 1)

#--Unscale Values
center <- attr(x, "scaled:center")
scale <- attr(x, "scaled:scale")

unscale <- function(x) {
  t(t(x) * scale + center)
}


#--Make Plot

pdf(
  paste0(figure.path, "FIG.Det.Probability.by.species.density.pdf"),
  width = 8,
  height = 8
)

month.vec <- c(1, 3, 6, 9)
var.vec <- c("Bovine", "Whitetail", "Nilgai", "Equine", "Other Species")
lab.vec <- c("Bovine", "Whitetail", "Nilgai", "Equine", "Other")

for (i in 1:length(var.vec)) {
  newData <- data.frame(
    x.mean = x.mean.var,
    species = factor(
      var.vec[i],
      levels = c(
        "Bovine",
        "Equine",
        "None",
        "Whitetail",
        "Nilgai",
        "Other Species"
      )
    ),
    var.plt = x.var
  )
  colnames(newData)[1] <- var.mean
  colnames(newData)[3] <- var.plt

  tmp <- predict(mod, type = 'det', newdata = newData, appendData = TRUE)

  tmp$x <- unscale(tmp[, var.plt])

  if (i == 1) {
    plot(
      1,
      xlim = c(-100, max(tmp$x)),
      ylim = c(0, 1),
      type = "n",
      axes = FALSE,
      xlab = "Number of Animals Inspected",
      pch = 20,
      ylab = "Detection Probability",
      cex.lab = 1.25,
      cex.main = 1.75
    )
  }

  lines(tmp$x, tmp$Predicted, col = "black", lwd = 2)

  text(x = 0, y = tmp$Predicted[1], labels = lab.vec[i], pos = 2)
} #END Plot

axis(side = 1)
axis(side = 2, las = 2)

#x<-mean(occ.raw$qty_inspected,na.rm=TRUE)
#text(x=x, y=0,
#     label="Mean number of\nanimals inspected per visit",pos=4, col="red")

#abline(v=x,col="red")

mtext(
  side = 3,
  adj = 0,
  "a) Functional relationship between detection probability and animal density"
)

dev.off()

#---- END Detection Probability by Species

#---- (NOT WORKING) Estimate Latent Proportion of Sites Occupied (PAO) ----
re <- ranef(mod)
EBUP <- bup(re, stat = "mode")
CI <- confint(re, level = 0.95)

rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / nrow(CI))


#--Bind with site level data
pao.dat <- cbind.data.frame(re@post, ebup, occ.dat@siteCovs)

cnty.vec <- unique(occ.dat@siteCovs$county_name)

for (i in 1:length(cnty.vec)) {
  val.loc <- which(pao.dat$county_name == cnty.vec[i])

  tmp <- post.df[val.loc, ]
  if (nrow(tmp) == 1) {
    next
  }

  tmp <- cbind.data.frame(county = cnty.vec[i], p = tmp[tmp$N == 1, "p"])

  if (i == 1) {
    latent.est <- tmp
  }
  if (i > 1) {
    latent.est <- rbind.data.frame(latent.est, tmp)
  }
}

all.est <- latent.est
all.est$county <- "All"

latent.est <- rbind.data.frame(latent.est, all.est)


cnt <- plyr::count(latent.est$county)
x <- cnt[cnt$freq > 3, "x"]

plt.dat <- latent.est[latent.est$county %in% x, ]
plt.dat$county <- as.character(plt.dat$county)


boxplot(
  p ~ county,
  data = plt.dat,
  horizontal = TRUE,
  axes = FALSE,
  xlab = "",
  ylab = ""
)
axis(side = 2, at = seq(1, length(x), 1), labels = x, las = 2)
axis(side = 1)
mtext(side = 1, "Estimate Latent Proportion of Sites Occupied", line = 2)


#-------------------
#---- FUNCTIONS ----

bootstrap.CI <- function(in.dat) {
  sampmean <- function(y, i) mean(y[i])
  bootmean <- boot(data = in.dat, statistic = sampmean, R = 10000)
  out.ci <- boot.ci(bootmean, conf = .95, type = c("norm"))

  out <- list(mu = bootmean$t, summary = out.ci)
  return(out)
} #END Function


antilogit <- function(x) {
  exp(x) / (1 + exp(x))
}


#-- FNC: Estimate Pasture Level Probability of Detection on Successive Surveys

est.prob.detection.repeated.sampling <- function(det.stat, n.surveys) {
  #--Controls
  n.iter <- length(det.stat$mu)
  n.surveys <- 20

  #--Simulate for testing
  #psi<- rbeta(n=n.iter, shape1=2, shape2=5)
  psi <- det.stat$mu

  #--Storage
  p.star <- array(NA, dim = c(n.iter, n.surveys))

  #--Set up x data (surveys)
  for (i in 1:n.surveys) {
    y <- rep(i, n.iter)
    if (i == 1) {
      x <- y
    }
    if (i > 1) {
      x <- cbind(x, y)
    }
  }

  #--Estimate Prob Detection
  for (i in 1:n.iter) {
    for (j in 1:n.surveys) {
      p.star[i, j] <- 1 - (1 - psi[i])^j
    } #END j
  } #END i

  out <- list(x = x, p.star = p.star)

  return(out)
} #END Function
