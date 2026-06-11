#-------------------
#
# By: Ryan Miller
#
# Dynamic Occupancy models
#
# Last updated: 15 Jan 2021
#
#------------------

## Clean Workspace
rm(list = ls(all = TRUE))
gc()


round(memory.limit() / 2^20, 2)
memory.size(max = TRUE)


#---- Set Directories (may need to change as needed) ----
code.path <- "C:/Documents/Project Documents/CattleFeverTick/Code/"
model.dat.path <- "C:/Documents/Project Documents/CattleFeverTick/Model.Data/"
model.results.path <- "C:/Documents/Project Documents/CattleFeverTick/Model.Results/"
figure.path <- "C:/Documents/Project Documents/CattleFeverTick/Figures/"

setwd(code.path)
#---- END Directories ----

#---- Load Packages for Session ----

#-- Load some useful aliases and functions
source("C:/Documents/R CODE/Rprofile/Rprofile.R")
source("C:/Documents/R CODE/Rprofile/Various.R")

#----Load Additional Packages for Session----
toLoad = c(
  "unmarked",
  "lubridate",
  "plyr",
  "dplyr",
  "AICcmodavg",
  "ROCR",
  "car",
  "ggplot2",
  "grid",
  "sp",
  "plyr",
  "raster",
  "wiqid",
  "unmarked",
  "reshape",
  "operators",
  "boot",
  "stringr"
)
instant_pkgs(toLoad)

#-- Load a modification of an`unmarked` function and custom plotting function
source("./newmodSel.R")
source("./plotOccu.R")
source("./update_formatMult.R") #Fix error in formatMulti function in unmarked package
source("./utils.R")
source(paste0(code.path, "Supporting.Functions.R"))

#---- END Load libraries ----

#---- READ DATA ----

#--Occurance Data via Surveillance

file.name <- "dat.occ.cov.sp.covariates.2021-02-05.csv"

occ.dat <- read.csv(paste0(model.dat.path, file.name))

length(unique(occ.dat$pasture_name))

#--Remove some of the prems with only 1 sample
cnt <- plyr::count(occ.dat[, c("pasture_name")])
cnt <- cnt[cnt$freq == 1, "x"]

#--Randomly select sites with only 1 observation
remove.vec <- sample(cnt, size = length(cnt) * .80, replace = FALSE)
occ.dat <- occ.dat[occ.dat$pasture_name %!in% remove.vec, ]

#--Drop some variables with large numbers of NA
occ.dat <- occ.dat[,
  colnames(occ.dat) %!in%
    c("lag.1.count", "lag.1.rate", "FIPS", "frm", "inv", "mean.farm.size")
]

summary(occ.dat)

#--Subset data for code development
#cnt<-plyr::count(occ.dat[,c("pasture_name")])
#cnt<-cnt[cnt$freq>=5,"x"]
#remove.vec<-sample(cnt, size=length(cnt)*.50,replace=FALSE)
#occ.dat<-occ.dat[occ.dat$pasture_name %!in% remove.vec,]
#write.csv(occ.dat, paste0(model.dat.path,"subset.", file.name))

#--Fix issue with county (temporary fix)
vec <- c("253-602a, 253-602A", "253-598", "214-542")
occ.dat[occ.dat$pasture_name %in% vec, "county_name"] <- "Zapata"


#--Reorder for formateMulti function
#Column 1 = year number
#Column 2 = site name or number
#Column 3 = julian date or chronological sample number during year
#Column 4 = observations (y)
#Column 5 to Final Column = covariates

col.names <- colnames(occ.dat)
col.names <- col.names[
  col.names %!in% c("pasture_name", "sample.year", "y", "year")
]

occ.dat <- occ.dat[, c("year", "pasture_name", "sample.year", "y", col.names)]


#---- Set Factors and Reference Categories ----

#--Add county grouping variable
tmp <- plyr::count(occ.dat$cnty.index)
tmp[order(tmp$freq), ]

#County list
cnty.vec <- c(
  "Kinney",
  "Maverick",
  "Dimmit",
  "Webb",
  "Zapata",
  "Starr",
  "Hidalgo",
  "Cameron",
  "Jim Hogg",
  "Willacy",
  "Brooks",
  "La Salle",
  "Uvalde",
  "Zavala",
  "Jim Wells"
)

#Set counties to "other"
occ.dat$cnty.index <- occ.dat$county_name
occ.dat[occ.dat$county_name %!in% cnty.vec, "cnty.index"] <- "other"

#--Set factor relative to other counties
occ.dat$cnty.index <- factor(occ.dat$cnty.index)
occ.dat <- within(occ.dat, cnty.index <- relevel(cnty.index, ref = 8))

#--Set inspection factor relative to 14 day inspection
occ.dat$inspection_type <- factor(occ.dat$inspection_type)
occ.dat <- within(occ.dat, inspection_type <- relevel(inspection_type, ref = 1))

#---- END Set Factors

#---- Add Site Level Unique ID for Linking Data ----

occ.dat$unk.id <- as.numeric(as.factor(occ.dat$pasture_name))

#-- Make site level characteristics dataframe
site.char <- unique(occ.dat[, c(
  "pasture_name",
  "county_name",
  "pasture_longitude",
  "pasture_latitude",
  "cnty.index",
  "unk.id"
)])
#-- END

#---- Make unmarked model object (takes a a long time) ----
occ.dat <- formatMult(occ.dat)

#--Save/read RDS to reduce reloading during future seasons
#saveRDS(occ.dat, file=paste0(model.dat.path,"formated.occ.dat.",Sys.Date(),".RDS"))

#saveRDS(occ.dat, file=paste0(model.dat.path,"subset.formated.occ.dat.",Sys.Date(),".RDS"))
occ.dat <- readRDS(paste0(model.dat.path, "formated.occ.dat.2021-02-10.RDS"))


#---- MAKE Site LUT ----

#-- Extract sites
site.dat <- siteCovs(occ.dat)

#-- Identify sites never infested
num.events.pos <- rowSums(occ.dat@y, na.rm = TRUE)
site.dat <- cbind.data.frame(site.dat, num.events.pos)

#-- Identify number of events per site
num.events <- apply(occ.dat@y, MARGIN = 1, function(x) length(which(!is.na(x))))
site.dat <- cbind.data.frame(site.dat, num.events)

#-- Assign sites never positive
site.dat$never.pos <- 0
site.dat[site.dat$num.events.pos == 0, "never.pos"] <- 1

#-- Add pasture id
site.dat$pasture_name <- rownames(site.dat)

rownames(site.dat) <- as.numeric(seq(1, nrow(site.dat), 1))

#-- END Make Site LUT

#---- END Make Unmarked Model ----

#--Plot and summerize data
#plot(occ.dat)
#summary(occ.dat)

#---- Scale and Center predictors ----

#Observations
col.names <- c(
  "qty_inspected",
  "qty_infested",
  "month",
  "apr.infest.prev",
  "prop.new.animals",
  "site.dens.bovine",
  "site.dens.equine",
  "site.mean.pas.dist",
  "site.sd.pas.dist",
  "num.adj.pas",
  "lag.count",
  "lag.rate",
  "tick.suit",
  "med.dist.pos",
  "sd.dist.pos"
)
for (i in 1:length(col.names)) {
  obsCovs(occ.dat)[, col.names[i]] <- as.numeric(as.character(obsCovs(
    occ.dat
  )[, col.names[i]]))
  obsCovs(occ.dat)[, paste0(col.names[i], ".std")] <- scale(
    obsCovs(occ.dat)[, col.names[i]],
    center = TRUE,
    scale = TRUE
  )
}

#Sites
col.names <- c(
  "site.pasture_qty_acres",
  "tick.suit",
  "site.mean.pas.dist",
  "site.sd.pas.dist",
  "site.med.dist.pos",
  "site.sd.dist.pos",
  "site.num.adj.pas",
  "site.lag.count",
  "site.lag.rate",
  "mx.distance"
)
for (i in 1:length(col.names)) {
  siteCovs(occ.dat)[, col.names[i]] <- as.numeric(as.character(siteCovs(
    occ.dat
  )[, col.names[i]]))
  siteCovs(occ.dat)[, paste0(col.names[i], ".std")] <- scale(
    siteCovs(occ.dat)[, col.names[i]],
    center = TRUE,
    scale = TRUE
  )
}

#---- END Data ----

#---- Fit Occupancy Models from MacKenzie et. al (2003) ----

#--Set up occupancy models

# psiformula = first-season's occupancy
# gammaformula = colonization rate  (site level)
# epsilonformula = extinction rate  (site level)
# pformula = detection probability

#--Model with covariates for detection probability and colonization

#--Covariates to include
col.vars <- c("tick.suit.std", "site.lag.rate.std", "site.med.dist.pos.std")
ext.vars <- c("tick.suit.std")
p.vars <- c(
  "as.factor(season)",
  "qty_inspected.std",
  "site.dens.bovine.std",
  "species",
  "inspection_type",
  "cnty.index"
)

#--Make formulas
col.form <- as.formula(paste("~ ", paste(col.vars, collapse = "+")))
ext.form <- as.formula(paste("~ ", paste(ext.vars, collapse = "+")))
p.form <- as.formula(paste("~ ", paste(p.vars, collapse = "+")))

#--Fit Model
m.full <- colext(
  psiformula = ~1,
  gammaformula = col.form,
  epsilonformula = ext.form,
  pformula = p.form,
  data = occ.dat,
  method = "BFGS",
  se = TRUE,
  control = list(trace = 1, maxit = 1e4)
)
summary(m.full)

#---- END Fit Models ----

#---- Conduct Model Selection ----

#-- Dredge all possible combinations of the full occupancy model
#occ.dredge <- dredge(m.full)

#-- Model comparison to explore the results
#mc <- as.data.frame(occ.dredge) %>%
#  select(starts_with("psi(p"), df, AICc, delta, weight)

#-- Shorten names for printing
#names(mc) <- names(mc) %>%
#  str_extract("(?<=psi\\(pland_[0-9]{2}_)[a-z_]+") %>%
#  coalesce(names(mc))

#-- Print model selection table
#mutate_all(mc, ~ round(., 3)) %>%
#  head(18) %>%
#  knitr::kable()

# select models with the most support for model averaging (< 2.5 delta aicc)
#occ.dredge.delta <- get.models(occ.dredge, subset = weight >= .1)

# average models based on model weights
#occ.avg <- model.avg(occ.dredge.delta, fit = TRUE)

# model averaged coefficients for occupancy and detection probability
#coef(occ.avg)

# model comparison to explore the results for detection
#md <- as.data.frame(occ.dredge) %>%
#  select(starts_with("p("), df, AICc, delta, weight)

# shorten names for printing
#names(md) <- names(md) %>%
#  str_extract("(?<=p\\(pland_[0-9]{2}_)[a-z_]+") %>%
#  coalesce(names(md))

#coef(occ_avg) %>%
#  enframe()

#---- END Model Selection ----

#---- Compare Models ----

# Model comparison
#modList <- list(m0, m.det, m.full)
#names(modList) <- c("p(.)psi(.)gamma(.)epsilon(.)",
#                    "p(qty_inspected+month+inspection_type)psi(.)gamma(.)epsilon(.)",
#                    "p(qty_inspected+month+inspection_type)psi(.)gamma(qty_added+pasture_size)epsilon(.)")
#fmList <- fitList(fits = modList)

# Model selection
#ms <- modSel(fmList)

#---- END Compare Models ----

#---- Assign Best Model and Generate SE

#Assign Best model
mod <- m.full

#--Generate bootstrap standard errors for smoothed trajectory
# Bootstrap sample information propagates through to derived quantities
# below if method="nonparboot" is used

mod <- nonparboot(mod, B = 10) # This takes a while and should be at min 1000 for final model!

#--Save model to avoid having to refit in future sessions (can take some time)

saveRDS(
  mod,
  file = paste0(model.results.path, "fit.col.ext.model.", Sys.Date(), ".RDS")
)

#mod<-readRDS(paste0(model.results.path,"fit.col.ext.model.2021-02-11.RDS"))

#---- End Assign Best Model ----

#---- Evaluate Model ----

# Cross Validation
#cross.val.stat<-crossVal(mod, method="Kfold", folds=10, holdoutPct=0.25, statistic=RMSE_MAE)

# Assess goodness-of-fit
#parboot(mod)
#plot(mod)

# Goodness of fit test (MacKenzie & Bailey 2004)
gof_boot <- unmarked::mb.gof.test(mod, nsim = 20, plot.hist = TRUE)
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

mod.chisq <- parboot(mod, statistic = chisq, nsim = 10, parallel = FALSE)


# Function returning three fit-statistics.
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2, na.rm = TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm = TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm = TRUE)
  out <- c(SSE = sse, Chisq = chisq, freemanTukey = freeTuke)
  return(out)
}

pb <- parboot(mod, fitstats, nsim = 10, report = 1)
plot(pb, main = "")


# Finite-sample inference for a derived parameter.
# Population size in sampled area

Nhat <- function(fm) {
  sum(bup(ranef(fm, K = 10)))
}

pb.N <- parboot(mod, Nhat, nsim = 10, report = 5)

# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(mod, K = 10)))


#---- ROC Curve Evaluation ----

#-- Get Fitted and Observed Data
x <- fitted(mod)
y <- getY(mod)

#-- Subset by county
cnty.vec <- unique(site.dat$county_name)
row.vec <- as.numeric(rownames(site.dat))

x.list <- vector(mode = "list")
y.list <- vector(mode = "list")

for (i in 1:length(cnty.vec)) {
  row.index <- row.vec[which(site.dat$county_name %in% cnty.vec[i])]

  y.vec <- c(y[row.index, ])
  y.vec <- y.vec[is.na(y.vec) == FALSE]

  x.vec <- c(x[row.index, ])
  x.vec <- x.vec[is.na(x.vec) == FALSE]

  if (sum(y.vec) > 0) {
    x.list[[i]] <- x.vec
    y.list[[i]] <- y.vec
  } #END Logical
} #END Loop

x.list <- x.list[lapply(x.list, length) > 0]
y.list <- y.list[lapply(y.list, length) > 0]

x.vec <- c(x)
y.vec <- c(y)
x.vec <- x.vec[is.na(x.vec) == FALSE]
y.vec <- y.vec[is.na(y.vec) == FALSE]


#tmp<-roc(response=x.vec,predictor=y.vec)

pred <- prediction(x.list, y.list)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")

plot(perf, col = "gray")

pred <- prediction(x.vec, y.vec)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")

plot(perf, col = "black", add = TRUE, lwd = 3)

df <- cbind.data.frame(
  x = as.vector(perf@x.values),
  y = as.vector(perf@y.values)
)
colnames(df) <- c("x", "y")

auc.val <- DescTools::AUC(x = df$x, y = df$y)

text(x = .5, y = .5, labels = paste0("AUC = ", round(auc.val, 3)), cex = 1.5)


#---- END ROC Evaluation ----

#---- END Evaluate Model ----

#---- Generate Model Outputs ----

#-- Make Model Estimates

ci.method = "normal" #Profile should be used for final results but takes awhile to run

#-- Expected extinction rate (population-level estimate)
# Coefficients: c(ext, 0,...) 0 inidicates to use the mean value of the variable
n.var <- nrow(summary(mod)$ext)
ext.est <- backTransform(linearComb(
  mod,
  coefficients = c(1, rep(0, n.var - 1)),
  type = 'ext'
))

ext.est <- cbind.data.frame(
  ext.est = ext.est@estimate,
  t(antilogit(confint(mod, type = 'ext', method = ci.method)[1, ]))
)


#-- Expected colinization rate (population-level estimate)
# Coefficients: c(col, 0,...) 0 inidicates to use the mean value of the variable
n.var <- nrow(summary(mod)$col)
col.est <- backTransform(linearComb(
  mod,
  coefficients = c(1, rep(0, n.var - 1)),
  type = 'col'
))

col.est <- cbind.data.frame(
  det.est = col.est@estimate,
  t(antilogit(confint(mod, type = 'col', method = ci.method)[1, ]))
)

# Coefficient Table (logit-scale)
estimate <- mod@estimates[2]@estimates
col.beta.mat <- cbind.data.frame(estimate, confint(mod, type = 'col'))

# Convert to non-logit scale
col.beta.mat <- antilogit(col.beta.mat)


#-- Probability of detection given a site is occupied
# Coefficients: c(det, 0,...) 0 inidicates to use the mean value of the variable
n.var <- nrow(summary(mod)$det)
det.est <- backTransform(linearComb(
  mod,
  coefficients = c(1, rep(0, n.var - 1)),
  type = 'det'
))

det.est <- cbind.data.frame(
  det.est = det.est@estimate,
  t(antilogit(confint(mod, type = 'det', method = ci.method)[1, ]))
)

# Coefficient Table (logit-scale)
estimate <- mod@estimates[4]@estimates
det.beta.mat <- cbind.data.frame(
  estimate,
  confint(mod, type = 'det', method = ci.method)
)

# Convert to non-logit scale
det.beta.mat <- antilogit(det.beta.mat)


#--Get site level detection estimates
det.vec <- getP(mod)

#Make site detection means and sd
det.site.mu <- rowMeans(det.vec, na.rm = TRUE)
det.site.sd <- apply(det.vec, 1, sd, na.rm = TRUE)

det.site.mu.all <- cbind.data.frame(site.dat, det.site.mu)

#Generate bootstrap
#det.stat<-bootstrap.CI(det.site.mu.all$det.site.mu)

#Make site detection for previously infested
det.site.mu.pos <- det.site.mu.all[det.site.mu.all$never.pos == 0, ]

#Generate bootstrap
#det.stat.pos<-bootstrap.CI(det.site.mu)

### NEED to rework

#--Make Site level detection estimates with county
site.dat <- cbind.data.frame(siteCovs(occ.dat), det.site.mu)

site.dat <- cbind.data.frame(site.dat, det.site.sd)

site.lut <- unique(obsCovs(occ.dat)[, c(
  "county_name",
  "cnty.index",
  "pasture_longitude",
  "pasture_latitude"
)])

site.lut <- merge(
  site.dat,
  site.lut,
  by = c("pasture_longitude", "pasture_latitude"),
  all.x = TRUE
)


#--Confidence intervals for projected occupancy (needs work)
occ.est <- cbind.data.frame(
  projected = projected(mod)[2, ],
  smoothed = smoothed(mod)[2, ],
  SE = mod@projected.mean.bsse[2, ]
)

occ.est$l95ci <- occ.est$projected +
  qt(c(0.025), sampleSize(mod) - 1) * occ.est$SE
occ.est$u95ci <- occ.est$projected +
  qt(c(0.975), sampleSize(mod) - 1) * occ.est$SE

#Hack for CI less than zero (need to fix)
occ.est[occ.est$l95ci < 0, "l95ci"] <- 0


#--Equilibrium occupancy
equil.occ <- col.beta.mat$estimate[1] / (col.beta.mat$estimate[1] + ext.est[1])


#--Write outputs

write.csv(
  equil.occ,
  paste0(write.path, "Equilibrium.occupancy.", Sys.Date(), ".csv"),
  row.names = FALSE
)

write.csv(
  occ.est,
  paste0(write.path, "occupancy.probability.population.", Sys.Date(), ".csv"),
  row.names = FALSE
)


#-- Make apparent infestation rate

# Estimates of conditional occupancy distribution at each site and time point
re <- ranef(mod)

#--Draw from the posterior predictive distribution
ppd <- posteriorSamples(re, nsims = 10000)

#--Generate means
ppd.mean <- posterior.mean(ppd)
#ppd.cnt<-posterior.cnt.success(ppd)

#--Merge with sites
ppd.mean <- cbind.data.frame(site.dat, ppd.mean)

#--Never Infested
ppd.never <- ppd.mean[ppd.mean$never.pos == 1, ]
ppd.never <- ppd.never[, (ncol(ppd.never) - 6):ncol(ppd.never)]
mean.never <- colMeans(ppd.never)

ci.lower <- apply(ppd.never, MARGIN = 2, FUN = function(x) {
  mean(x) + qt(c(0.025), length(x) - 1) * std.error(x)
})
ci.upper <- apply(ppd.never, MARGIN = 2, FUN = function(x) {
  mean(x) + qt(c(0.975), length(x) - 1) * std.error(x)
})

mean.never <- cbind.data.frame(
  mean = mean.never,
  ci.lower = ci.lower,
  ci.upper = ci.upper
)


#--Infested
ppd.pos <- ppd.mean[ppd.mean$never.pos == 0, ]
ppd.pos <- ppd.pos[, (ncol(ppd.pos) - 6):ncol(ppd.pos)]
mean.pos <- colMeans(ppd.pos)

ci.pos <- t(apply(ppd.pos, MARGIN = 2, FUN = function(x) {
  quantile(x, probs = c(.25, .75))
}))


ci.lower <- apply(ppd.pos, MARGIN = 2, FUN = function(x) {
  mean(x) + qt(c(0.025), length(x) - 1) * std.error(x)
})
ci.upper <- apply(ppd.pos, MARGIN = 2, FUN = function(x) {
  mean(x) + qt(c(0.975), length(x) - 1) * std.error(x)
})

mean.pos <- cbind.data.frame(
  mean = mean.pos,
  ci.lower = ci.lower,
  ci.upper = ci.upper
)


#--Calculate Turn Over (not working yet)
turnover <- function(fm) {
  psi.hat <- plogis(coef(fm, type = "psi"))
  if (length(psi.hat) > 1) {
    stop("this function only works if psi is scalar")
  }
  T <- unmarked::getData(fm)@numPrimary
  tau.hat <- numeric(T - 1)
  gamma.hat <- plogis(coef(fm, type = "col"))
  phi.hat <- 1 - plogis(coef(fm, type = "ext"))
  if (length(gamma.hat) != T - 1 | length(phi.hat) != T - 1) {
    stop("this function only works if gamma and phi T-1 vectors")
  }
  for (t in 2:T) {
    psi.hat[t] <- psi.hat[t - 1] *
      phi.hat[t - 1] +
      (1 - psi.hat[t - 1]) * gamma.hat[t - 1]
    tau.hat[t - 1] <- gamma.hat[t - 1] * (1 - psi.hat[t - 1]) / psi.hat[t]
  }
  return(tau.hat)
} #END

pb <- parboot(m1, statistic = turnover, nsim = 2)
turnCI <- cbind(
  pb@t0,
  t(apply(pb@t.star, 2, quantile, probs = c(0.025, 0.975)))
)
colnames(turnCI) <- c("tau", "lower", "upper")


#---- End Make Model Estimates ----

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


#---- FNC - Plot regression coeff

plot.betas <- function(plt.dat) {
  x.labs <- pretty(as.vector(as.matrix(plt.dat)))

  x.min <- min(x.labs)
  x.max <- max(x.labs)

  plt.dat$name <- rownames(plt.dat)
  plt.dat$y <- seq(1, nrow(plt.dat), 1)

  plot(
    x = plt.dat$estimate,
    y = plt.dat$y,
    xlim = c(x.min, x.max),
    pch = 19,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )

  abline(v = 0, col = "gray")

  segments(
    x0 = plt.dat$`0.025`,
    x1 = plt.dat$`0.975`,
    y0 = plt.dat$y,
    y1 = plt.dat$y,
    lwd = 1.5
  )
  points(x = plt.dat$estimate, y = plt.dat$y, pch = 19)

  axis(side = 1)
  axis(side = 2, at = plt.dat$y, labels = plt.dat$name, las = 2)
} #END Plot Betas


generate.annual.subsets <- function(occ.raw, year.val) {
  occ.raw <- occ.raw[occ.raw$year == year.val, ]

  col.names <- colnames(occ.raw)
  col.names <- col.names[
    col.names %!in% c("pasture_name", "jdate", "y", "year")
  ]

  occ.raw <- occ.raw[, c("year", "pasture_name", "jdate", "y", col.names)]

  #--Add Site Year
  occ.raw$site.year <- occ.raw$year

  #--Make unmarked model object (takes a a long time)
  occ.raw <- formatMult(occ.raw)

  #Observations
  col.names <- c(
    "qty_inspected",
    "qty_infested",
    "apr.infest.prev",
    "prop.new.animals",
    "density",
    "med.dist.pos",
    "sd.dist.pos"
  )
  for (i in 1:length(col.names)) {
    obsCovs(occ.raw)[, col.names[i]] <- as.numeric(as.character(obsCovs(
      occ.raw
    )[, col.names[i]]))
    obsCovs(occ.raw)[, paste0(col.names[i], ".std")] <- scale(
      obsCovs(occ.raw)[, col.names[i]],
      center = TRUE,
      scale = TRUE
    )
  }

  #Sites
  col.names <- c(
    "pasture_qty_acres",
    "mx.distance",
    "site.med.dist.pos",
    "site.sd.dist.pos",
    "qty_added"
  )
  for (i in 1:length(col.names)) {
    siteCovs(occ.raw)[, col.names[i]] <- as.numeric(as.character(siteCovs(
      occ.raw
    )[, col.names[i]]))
    siteCovs(occ.raw)[, paste0(col.names[i], ".std")] <- scale(
      siteCovs(occ.raw)[, col.names[i]],
      center = TRUE,
      scale = TRUE
    )
  }

  return(occ.raw)
} #END
