

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
rm(list=ls(all=TRUE))
gc()


round(memory.limit()/2^20, 2)
memory.size(max=TRUE)



#---- Set Directories (may need to change as needed) ----
code.path<-"C:/Documents/Project Documents/CattleFeverTick/Code/"
model.dat.path<-"C:/Documents/Project Documents/CattleFeverTick/Model.Data/"
model.results.path<-"C:/Documents/Project Documents/CattleFeverTick/Model.Results/"
figure.path<-"C:/Documents/Project Documents/CattleFeverTick/Figures/"

setwd(code.path)
#---- END Directories ----



#---- Load Packages for Session ----

#-- Load some useful aliases and functions
source("C:/Documents/R CODE/Rprofile/Rprofile.R")
source("C:/Documents/R CODE/Rprofile/Various.R")

#----Load Additional Packages for Session----
toLoad = c("unmarked", "lubridate", "plyr", "dplyr", "AICcmodavg",
          "ROCR", "car", "ggplot2", "grid","sp","plyr","raster","wiqid","unmarked","reshape","operators","boot","stringr")
instant_pkgs(toLoad)

#-- Load a modification of an`unmarked` function and custom plotting function
source("./newmodSel.R")
source("./plotOccu.R")
source("./update_formatMult.R") #Fix error in formatMulti function in unmarked package
source("./utils.R")
source(paste0(code.path,"Supporting.Functions.R"))

#---- END Load libraries ----



#---- READ DATA ----

#--Occurance Data via Surveillance

file.name <- "dat.occ.cov.sp.covariates.2021-02-05.csv"

occ.dat<-read.csv( paste0(model.dat.path, file.name) )

length(unique(occ.dat$pasture_name))

#--Remove some of the prems with only 1 sample
cnt<-plyr::count(occ.dat[,c("pasture_name")])
cnt<-cnt[cnt$freq==1,"x"]

#--Randomly select sites with only 1 observation
remove.vec<-sample(cnt, size=length(cnt)*.80,replace=FALSE)
occ.dat<-occ.dat[occ.dat$pasture_name %!in% remove.vec,]

#--Drop some variables with large numbers of NA
occ.dat<-occ.dat[,colnames(occ.dat) %!in% c("lag.1.count","lag.1.rate","FIPS","frm","inv","mean.farm.size")]

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

col.names<-colnames(occ.dat)
col.names<-col.names[col.names %!in% c("pasture_name","sample.year","y","year")]

occ.dat<-occ.dat[,c("year","pasture_name","sample.year","y",col.names)]


#---- Set Factors and Reference Categories ----

#--Add county grouping variable
tmp<-plyr::count(occ.dat$cnty.index)
tmp[order(tmp$freq),]

#County list
cnty.vec<-c("Kinney", "Maverick", "Dimmit", "Webb", "Zapata", "Starr", "Hidalgo", "Cameron",
            "Jim Hogg", "Willacy", "Brooks", "La Salle", "Uvalde", "Zavala", "Jim Wells")

#Set counties to "other"
occ.dat$cnty.index <- occ.dat$county_name
occ.dat[occ.dat$county_name %!in% cnty.vec,"cnty.index"] <- "other"

#--Set factor relative to other counties
occ.dat$cnty.index <- factor(occ.dat$cnty.index)
occ.dat <- within(occ.dat, cnty.index <- relevel(cnty.index, ref=8))

#--Set inspection factor relative to 14 day inspection
occ.dat$inspection_type <- factor(occ.dat$inspection_type)
occ.dat <- within(occ.dat, inspection_type <- relevel(inspection_type, ref=1))

#---- END Set Factors


#---- Add Site Level Unique ID for Linking Data ----

occ.dat$unk.id<-as.numeric(as.factor(occ.dat$pasture_name))

#-- Make site level characteristics dataframe
site.char <- unique(occ.dat[,c("pasture_name","county_name",
                               "pasture_longitude","pasture_latitude",
                               "cnty.index","unk.id")])
#-- END


#---- Make unmarked model object (takes a a long time) ----
occ.dat<-formatMult(occ.dat)

#--Save/read RDS to reduce reloading during future seasons
#saveRDS(occ.dat, file=paste0(model.dat.path,"formated.occ.dat.",Sys.Date(),".RDS"))

#saveRDS(occ.dat, file=paste0(model.dat.path,"subset.formated.occ.dat.",Sys.Date(),".RDS"))
occ.dat<-readRDS(paste0(model.dat.path,"formated.occ.dat.2021-02-10.RDS"))



#---- MAKE Site LUT ----

#-- Extract sites 
site.dat<-siteCovs(occ.dat)

#-- Identify sites never infested
num.events.pos<-rowSums(occ.dat@y, na.rm=TRUE)
site.dat<-cbind.data.frame(site.dat,num.events.pos)

#-- Identify number of events per site
num.events<-apply(occ.dat@y, MARGIN=1, function(x) length(which(!is.na(x))))
site.dat<-cbind.data.frame(site.dat,num.events)

#-- Assign sites never positive
site.dat$never.pos <- 0
site.dat[site.dat$num.events.pos==0,"never.pos"] <- 1

#-- Add pasture id
site.dat$pasture_name <- rownames(site.dat)

rownames(site.dat)<-as.numeric(seq(1,nrow(site.dat),1))

#-- END Make Site LUT




#---- END Make Unmarked Model ----


#--Plot and summerize data
#plot(occ.dat)
#summary(occ.dat)


#---- Scale and Center predictors ----

#Observations
col.names<-c("qty_inspected","qty_infested","month",
             "apr.infest.prev","prop.new.animals", "site.dens.bovine", "site.dens.equine",
             "site.mean.pas.dist","site.sd.pas.dist",
             "num.adj.pas",
             "lag.count","lag.rate",
             "tick.suit","med.dist.pos","sd.dist.pos")
for(i in 1:length(col.names)){
  obsCovs(occ.dat)[,col.names[i]]<-as.numeric(as.character(obsCovs(occ.dat)[,col.names[i]]))
  obsCovs(occ.dat)[,paste0(col.names[i],".std")] <- scale(obsCovs(occ.dat)[,col.names[i]], center = TRUE, scale = TRUE)
}

#Sites
col.names<-c("site.pasture_qty_acres","tick.suit",
             "site.mean.pas.dist","site.sd.pas.dist",
             "site.med.dist.pos","site.sd.dist.pos","site.num.adj.pas",
             "site.lag.count","site.lag.rate","mx.distance")
for(i in 1:length(col.names)){
  siteCovs(occ.dat)[,col.names[i]]<-as.numeric(as.character(siteCovs(occ.dat)[,col.names[i]]))
  siteCovs(occ.dat)[,paste0(col.names[i],".std")] <- scale(siteCovs(occ.dat)[,col.names[i]], center = TRUE, scale = TRUE)
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
col.vars<-c("tick.suit.std", "site.lag.rate.std", "site.med.dist.pos.std")
ext.vars<-c("tick.suit.std")
p.vars<-c("as.factor(season)","qty_inspected.std","site.dens.bovine.std","species","inspection_type","cnty.index")

#--Make formulas
col.form <- as.formula(paste("~ ", paste(col.vars, collapse= "+")))
ext.form <- as.formula(paste("~ ", paste(ext.vars, collapse= "+")))
p.form <- as.formula(paste("~ ", paste(p.vars, collapse= "+")))

#--Fit Model
m.full <- colext(psiformula= ~1, 
                 gammaformula = col.form, 
                 epsilonformula = ext.form, 
                 pformula = p.form, 
                 data = occ.dat,
                 method="BFGS", se=TRUE, 
                 control = list(trace=1, maxit = 1e4))
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
mod<-m.full

#--Generate bootstrap standard errors for smoothed trajectory
# Bootstrap sample information propagates through to derived quantities
# below if method="nonparboot" is used

mod <- nonparboot(mod, B = 2)  # This takes a while and should be at min 1000 for final model!

#--Save model to avoid having to refit in future sessions (can take some time)
#saveRDS(mod, file=paste0(model.results.path,"fit.col.ext.model.",Sys.Date(),".RDS"))

#mod<-readRDS(paste0(model.results.path,"fit.col.ext.model.2021-02-11.RDS"))

#---- End Assign Best Model ----




#---- Evaluate Model ----

# Cross Validation
#cross.val.stat<-crossVal(mod, method="Kfold", folds=10, holdoutPct=0.25, statistic=RMSE_MAE)

# Assess goodness-of-fit
#parboot(mod)
#plot(mod)

# Goodness of fit test (MacKenzie & Bailey 2004)
gof_boot <- unmarked::mb.gof.test(mod, nsim=20, plot.hist=TRUE)
gof_boot

# values of c-hat > 1 indicate overdispersion (variance > mean),
# but … values much higher than 1 (i.e., > 4) probably indicate
# lack-of-fit. In cases of moderate overdispersion, one usually
# multiplies the variance-covariance matrix of the estimates by c-hat.
# As a result, the SE’s of the estimates are inflated (c-hat is also known
# as a variance inflation factor)

#The parametric bootstrap can be used to check the adequacy of model ???t. Here we use a ??2 statistic appropriate for binary data.

chisq <- function(fm) {
  umf <- getData(fm)
  y <- getY(umf)
  y[y>1] <- 1
  sr <- fm@sitesRemoved
  if(length(sr)>0)
    y <- y[-sr,,drop=FALSE]
  fv <- fitted(fm, na.rm=TRUE)
  y[is.na(fv)] <- NA
  sum((y-fv)^2/(fv*(1-fv)), na.rm=TRUE) }#END Function

mod.chisq <- parboot(mod, statistic=chisq, nsim=10, parallel=FALSE)


# Function returning three fit-statistics.
fitstats <- function(fm) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2, na.rm=TRUE)
  chisq <- sum((observed - expected)^2 / expected, na.rm=TRUE)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=TRUE)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

pb <- parboot(mod, fitstats, nsim=10, report=1)
plot(pb, main="")


# Finite-sample inference for a derived parameter.
# Population size in sampled area

Nhat <- function(fm) {
  sum(bup(ranef(fm, K=10)))
}

pb.N <- parboot(mod, Nhat, nsim=10, report=5)

# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(mod, K=10)))


#---- ROC Curve Evaluation ----

#-- Get Fitted and Observed Data
x<-fitted(mod)
y<-getY(mod)

#-- Subset by county
cnty.vec<-unique(site.dat$county_name)
row.vec<-as.numeric(rownames(site.dat))

x.list <- vector(mode = "list")
y.list <- vector(mode = "list")

for(i in 1:length(cnty.vec)){
  row.index<-row.vec[which(site.dat$county_name %in% cnty.vec[i])]
  
  y.vec<-c(y[row.index,])
  y.vec<-y.vec[is.na(y.vec)==FALSE]
  
  x.vec<-c(x[row.index,])
  x.vec<-x.vec[is.na(x.vec)==FALSE]
  
  if(sum(y.vec)>0){
    x.list[[i]]<-x.vec
    y.list[[i]]<-y.vec
  }#END Logical
}#END Loop

x.list<-x.list[lapply(x.list,length)>0]
y.list<-y.list[lapply(y.list,length)>0]

x.vec<-c(x)
y.vec<-c(y)
x.vec<-x.vec[is.na(x.vec)==FALSE]
y.vec<-y.vec[is.na(y.vec)==FALSE]


#tmp<-roc(response=x.vec,predictor=y.vec)


pred<-prediction(x.list, y.list)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")

plot(perf, col="gray")

pred<-prediction(x.vec, y.vec)
perf <- performance(pred, measure = "tpr", x.measure = "fpr")

plot(perf, col="black",add=TRUE,lwd=3)

df<-cbind.data.frame(x=as.vector(perf@x.values),y=as.vector(perf@y.values))
colnames(df)<-c("x","y")

auc.val<-DescTools::AUC(x=df$x,y=df$y)

text(x=.5,y=.5, labels=paste0("AUC = ", round(auc.val,3)), cex=1.5)


#---- END ROC Evaluation ----

#---- END Evaluate Model ----





#---- Generate Model Outputs ----

#-- Make Model Estimates

ci.method="normal" #Profile should be used for final results but takes awhile to run

#-- Expected extinction rate (population-level estimate)
# Coefficients: c(ext, 0,...) 0 inidicates to use the mean value of the variable
n.var<-nrow(summary(mod)$ext)
ext.est<-backTransform(linearComb(mod, coefficients=c(1,rep(0,n.var-1)), type='ext')) 

ext.est<-cbind.data.frame(ext.est=ext.est@estimate, t(antilogit(confint(mod, type='ext',method=ci.method)[1,])))



#-- Expected colinization rate (population-level estimate)
# Coefficients: c(col, 0,...) 0 inidicates to use the mean value of the variable
n.var<-nrow(summary(mod)$col)
col.est<-backTransform(linearComb(mod, coefficients=c(1,rep(0,n.var-1)), type='col')) 

col.est<-cbind.data.frame(det.est=col.est@estimate, t(antilogit(confint(mod, type='col',method=ci.method)[1,])))

# Coefficient Table (logit-scale)
estimate<-mod@estimates[2]@estimates
col.beta.mat<-cbind.data.frame(estimate,confint(mod, type='col'))

# Convert to non-logit scale 
col.beta.mat<-antilogit(col.beta.mat)



#-- Probability of detection given a site is occupied
# Coefficients: c(det, 0,...) 0 inidicates to use the mean value of the variable
n.var<-nrow(summary(mod)$det)
det.est<-backTransform(linearComb(mod, coefficients=c(1,rep(0,n.var-1)), type='det')) 

det.est<-cbind.data.frame(det.est=det.est@estimate, t(antilogit(confint(mod, type='det',method=ci.method)[1,])))

# Coefficient Table (logit-scale)
estimate<-mod@estimates[4]@estimates
det.beta.mat<-cbind.data.frame(estimate,confint(mod, type='det',method=ci.method))

# Convert to non-logit scale 
det.beta.mat<-antilogit(det.beta.mat)



#--Get site level detection estimates
det.vec<-getP(mod)

#Make site detection means and sd
det.site.mu<-rowMeans(det.vec, na.rm=TRUE)
det.site.sd<-apply(det.vec,1, sd, na.rm = TRUE)

det.site.mu.all<-cbind.data.frame(site.dat,det.site.mu)

#Generate bootstrap
#det.stat<-bootstrap.CI(det.site.mu.all$det.site.mu)


#Make site detection for previously infested
det.site.mu.pos<-det.site.mu.all[det.site.mu.all$never.pos==0,]

#Generate bootstrap
#det.stat.pos<-bootstrap.CI(det.site.mu)




### NEED to rework

#--Make Site level detection estimates with county
site.dat<-cbind.data.frame(siteCovs(occ.dat),det.site.mu)

site.dat<-cbind.data.frame(site.dat,det.site.sd)

site.lut<-unique(obsCovs(occ.dat)[,c("county_name","cnty.index","pasture_longitude", "pasture_latitude")])

site.lut<-merge(site.dat,site.lut, by=c("pasture_longitude", "pasture_latitude"),all.x=TRUE)











#--Confidence intervals for projected occupancy (needs work)
occ.est<-cbind.data.frame(projected=projected(mod)[2,], smoothed=smoothed(mod)[2,],SE=mod@projected.mean.bsse[2,])

occ.est$l95ci<-occ.est$projected + qt( c(0.025), sampleSize(mod) - 1) * occ.est$SE
occ.est$u95ci<-occ.est$projected + qt( c(0.975), sampleSize(mod) - 1) * occ.est$SE

#Hack for CI less than zero (need to fix)
occ.est[occ.est$l95ci<0,"l95ci"]<-0


#--Equilibrium occupancy
equil.occ<-col.beta.mat$estimate[1]/(col.beta.mat$estimate[1] + ext.est[1])










#-- Make apparent infestation rate







# Estimates of conditional occupancy distribution at each site and time point
re <- ranef(mod)

#--Draw from the posterior predictive distribution
ppd <- posteriorSamples(re, nsims=10000)

#--Generate means
ppd.mean <- posterior.mean(ppd)
#ppd.cnt<-posterior.cnt.success(ppd)


#--Merge with sites
ppd.mean<-cbind.data.frame(site.dat, ppd.mean)

#--Never Infested
ppd.never <- ppd.mean[ppd.mean$never.pos==1,]
ppd.never<-ppd.never[,(ncol(ppd.never)-6):ncol(ppd.never)]
mean.never<-colMeans(ppd.never)

ci.lower<-apply(ppd.never, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
ci.upper<-apply(ppd.never, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))

mean.never<-cbind.data.frame(mean=mean.never,ci.lower=ci.lower,ci.upper=ci.upper)


#--Infested
ppd.pos <- ppd.mean[ppd.mean$never.pos==0,]
ppd.pos<-ppd.pos[,(ncol(ppd.pos)-6):ncol(ppd.pos)]
mean.pos<-colMeans(ppd.pos)

ci.pos<-t(apply(ppd.pos, MARGIN=2, FUN=function(x) quantile(x, probs=c(.25,.75))))


ci.lower<-apply(ppd.pos, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
ci.upper<-apply(ppd.pos, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))

mean.pos<-cbind.data.frame(mean=mean.pos,ci.lower=ci.lower,ci.upper=ci.upper)










# Best Unbiased Predictors
tmp<-bup(re, stat="mean")           # Posterior mean
bup(re, stat="mode")           # Posterior mode

bup.ci<-unmarked::confint(re, level=0.9, type="state") # 90% CI


bup.ci[100,,]















#--Calculate Turn Over (not working yet)
turnover <- function(fm) {
  psi.hat <- plogis(coef(fm, type="psi"))
  if(length(psi.hat) > 1)
    stop("this function only works if psi is scalar")
  T <- unmarked::getData(fm)@numPrimary
  tau.hat <- numeric(T-1)
  gamma.hat <- plogis(coef(fm, type="col"))
  phi.hat <- 1 - plogis(coef(fm, type="ext"))
  if(length(gamma.hat) != T-1 | length(phi.hat) != T-1)
    stop("this function only works if gamma and phi T-1 vectors")
  for(t in 2:T) { psi.hat[t] <- psi.hat[t-1]*phi.hat[t-1] + (1-psi.hat[t-1])*gamma.hat[t-1]
  tau.hat[t-1] <- gamma.hat[t-1]*(1-psi.hat[t-1]) / psi.hat[t] }
  return(tau.hat)
  }#END

pb <- parboot(m1, statistic=turnover, nsim=2)
turnCI <- cbind(pb@t0, t(apply(pb@t.star, 2, quantile, probs=c(0.025, 0.975))))
colnames(turnCI) <- c("tau", "lower", "upper")




#---- End Make Model Estimates ----







#---- Plot Regression Coeff

pdf(paste0(figure.path,"Regression.Betas.",Sys.Date(),".pdf"), width=11,height=8)

mat<-matrix(c(1,4,
              2,4,
              3,5),ncol=2,nrow=3,byrow=TRUE)

layout(mat, widths=c(1,1,1), heights=c(1,1,1))

par(mar=c(3,10,3,1))

#--Detection
names.vec<-c("(Intercept)","Equine","Nilgai","Other_Wildlife","Whitetail")
vars.vec<-c("(Intercept)","speciesEquine","speciesNilgai",
           "speciesOther_Wildlife","speciesWhitetail")
plot.betas(det.beta.mat, names.vec, vars.vec)

mtext(side=3, line=1, adj=0, "a) Detection: species regression coefficents", cex=1.2)

#--Inspection type
vars.vec<-rownames(det.beta.mat[grep(rownames(det.beta.mat),pattern="inspection_type"),])
names.vec<-str_remove(vars.vec, pattern="inspection_type")
vars.vec<-c("(Intercept)",vars.vec)
names.vec<-c("(Intercept)",names.vec)

plot.betas(det.beta.mat, names.vec, vars.vec)

mtext(side=3, line=1, adj=0, "b) Detection: inspection type regression coefficents", cex=1.1)

#--prem characteristics
names.vec<-c("(Intercept)","Season","Number inspected","Density bovine")
vars.vec<-c("(Intercept)","as.factor(season)1","qty_inspected.std","site.dens.bovine.std")
plot.betas(det.beta.mat, names.vec, vars.vec)

mtext(side=3, line=1, adj=0, "c) Detection: prem characteristics regression coefficents", cex=1.1)

#--County
vars.vec<-rownames(det.beta.mat[grep(rownames(det.beta.mat),pattern="cnty.index"),])
names.vec<-str_remove(vars.vec, pattern="cnty.index")
vars.vec<-c("(Intercept)",vars.vec)
names.vec<-c("(Intercept)",names.vec)

plot.betas(det.beta.mat, names.vec, vars.vec)

mtext(side=3, line=1, adj=0, "d) Detection: county regression coefficents", cex=1.1)

#--Colinization
names.vec<-c("(Intercept)","Tick Suitability","Adjacent infestion rate",
             "Distance to infested")
vars.vec<-c("(Intercept)","tick.suit.std","site.lag.rate.std","site.med.dist.pos.std")
plot.betas(col.beta.mat, names.vec, vars.vec)

mtext(side=3, line=1, adj=0, "e) Colinization: regression coefficents", cex=1.1)

dev.off()

#---- END Plot ----






#---- Generate Annual Estimates of COL, EXT, DETECTION ----

#--Prep data for prediction
occ.raw<-read.csv(paste0(model.dat.path,file.name))

#Year vec
year.vec<-unique(occ.raw$year)
year.vec<-year.vec[is.na(year.vec)==FALSE]
year.vec<-sort(as.numeric(as.character(year.vec)))
year.vec<-year.vec[year.vec %!in% c(2020,2021)]

#--Observations
col.names<-unique(c(col.vars,ext.vars,p.vars))

#Remove factors
col.names<-col.names[-grep(col.names, pattern="as.fact")]
col.names<-str_remove_all(col.names, pattern=".std")

factor.vars<-c("cnty.index","species","inspection_type","as.factor(season)")

col.scale <- col.names[col.names %!in% factor.vars]

#--Scale Variables
for(i in 1:length(col.scale)){
  if(class(occ.raw[,col.scale[i]])!="character"){
    occ.raw[,paste0(col.scale[i],".std")] <- scale(occ.raw[,col.scale[i]], center=TRUE, scale=TRUE)
  }#END Logical
}#END Loop
  
#-- END Prep data

#-- Set Factors and Reference Categories

#--Add county grouping variable
tmp<-plyr::count(occ.raw$cnty.index)
tmp[order(tmp$freq),]

#County list
cnty.vec<-c("Kinney", "Maverick", "Dimmit", "Webb", "Zapata", "Starr", "Hidalgo", "Cameron",
            "Jim Hogg", "Willacy", "Brooks", "La Salle", "Uvalde", "Zavala", "Jim Wells")

#Set counties to "other"
occ.raw$cnty.index <- occ.raw$county_name
occ.raw[occ.raw$county_name %!in% cnty.vec,"cnty.index"] <- "other"

#--Set factor relative to other counties
occ.raw$cnty.index <- factor(occ.raw$cnty.index)
occ.raw <- within(occ.raw, cnty.index <- relevel(cnty.index, ref=8))

#--Set inspection factor relative to 14 day inspection
occ.raw$inspection_type <- factor(occ.raw$inspection_type)
occ.raw <- within(occ.raw, inspection_type <- relevel(inspection_type, ref=1))

#--Set season factor
occ.raw$season <- factor(occ.raw$season)

#-- END Set Factors


#-- Predict ext, col

for(i in 1:length(year.vec)){
  #Occupancy
  tmp <- predict(mod, type="psi", newdata=occ.raw[occ.raw$year==year.vec[i],])
  tmp <- unique(tmp)
  tmp<-apply(tmp, MARGIN=2, FUN=median) #Not the right way to do this but place holder
  tmp<-data.frame(year=year.vec[i], t(tmp))
  if(i==1){E.occ<-tmp}
  if(i>1){E.occ<-rbind(E.occ,tmp)}
  
  #Extinction
  tmp <- predict(mod, type='ext', newdata=occ.raw[occ.raw$year==year.vec[i],])
    tmp <- unique(tmp)
    tmp<-apply(tmp, MARGIN=2, FUN=median) #Not the right way to do this but place holder
    tmp<-data.frame(year=year.vec[i], t(tmp))
  if(i==1){E.ext<-tmp}
  if(i>1){E.ext<-rbind(E.ext,tmp)}
  
  #Colinization
  tmp <- predict(mod, type='col', newdata=occ.raw[occ.raw$year==year.vec[i],])
    tmp <- unique(tmp)
    tmp<-apply(tmp, MARGIN=2, FUN=median) #Not the right way to do this but place holder
    tmp<-data.frame(year=year.vec[i], t(tmp))
  if(i==1){E.col<-tmp}
  if(i>1){E.col<-rbind(E.col,tmp)}

  #Detection
  tmp <- predict(mod, type='det', newdata=occ.raw[occ.raw$year==year.vec[i],])
    tmp <- unique(tmp)
    tmp<-apply(tmp, MARGIN=2, FUN=median) #Not the right way to do this but place holder
    tmp<-data.frame(year=year.vec[i], t(tmp))
  if(i==1){E.det<-tmp}
  if(i>1){E.det<-rbind.data.frame(E.det,tmp)}
}#END

#---- END Predict ----






#---- PLOT Detection, EXT, COL ----

#-- Generate column means

#--Prep data for prediction
occ.raw<-read.csv(paste0(model.dat.path,file.name))


#-- Set Factors and Reference Categories

#--Add county grouping variable
tmp<-plyr::count(occ.raw$cnty.index)
tmp[order(tmp$freq),]

#County list
cnty.vec<-c("Kinney", "Maverick", "Dimmit", "Webb", "Zapata", "Starr", "Hidalgo", "Cameron",
            "Jim Hogg", "Willacy", "Brooks", "Uvalde", "Zavala", "Jim Wells")

#Set counties to "other"
occ.raw$cnty.index <- occ.raw$county_name
occ.raw[occ.raw$county_name %!in% cnty.vec,"cnty.index"] <- "other"

#--Set factor relative to other counties
occ.raw$cnty.index <- factor(occ.raw$cnty.index)
occ.raw <- within(occ.raw, cnty.index <- relevel(cnty.index, ref=8))

#--Set inspection factor relative to 14 day inspection
occ.raw$inspection_type <- factor(occ.raw$inspection_type)
occ.raw <- within(occ.raw, inspection_type <- relevel(inspection_type, ref=1))

#--Set season factor
occ.raw$season <- factor(occ.raw$season)

#-- END Set Factors



#--Column names
col.names<-unique(c(col.vars,ext.vars,p.vars))

#Alter names
col.names<-str_remove_all(col.names, pattern=".std")
col.names[col.names=="as.factor(season)"] <- "season"


#--Subset to only numeric
tmp<-occ.raw[,col.names]
col.class<-unlist(lapply(tmp,class))

col.numeric <- names(col.class[col.class %!in% c("character","factor")])
col.char <- names(col.class[col.class %in% c("character","factor")])


#-- Scale 
occ.raw[,col.numeric] <- scale(occ.raw[,col.numeric])


#-- Aggregate Data

#Generate formula
agg.formula <- formula(paste0("cbind(", paste(col.numeric, collapse=","),") ~ ",paste(c("FIPS","county_name","year", col.char), collapse="+")))

#Aggregate data
pred.dat <- aggregate(agg.formula, data=occ.raw, FUN=mean)

#Alter columns
colnames(pred.dat)[which(colnames(pred.dat) %in% col.numeric)] <- paste0(col.numeric,".std")

pred.dat$'as.factor(season)' <- pred.dat$season


#Occupancy
e.occ <- predict(mod, type="psi", newdata=pred.dat)
colnames(e.occ) <- paste0("psi.",colnames(e.occ))
  
e.ext <- predict(mod, type="ext", newdata=pred.dat)
colnames(e.ext) <- paste0("ext.",colnames(e.ext))
  
e.col <- predict(mod, type="col", newdata=pred.dat)
colnames(e.col) <- paste0("col.",colnames(e.col))
  
e.det <- predict(mod, type="det", newdata=pred.dat)
colnames(e.det) <- paste0("det.",colnames(e.det))
  
  
#Bind all predictions
x<-list(e.occ, e.ext, e.col, e.det)
preds <- do.call("cbind.data.frame", x)

pred.cols <- colnames(preds)

preds <- cbind.data.frame(pred.dat, preds)

  
#-- Aggregate Data
  
det.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ cnty.index+year, data=preds, FUN=mean)
col.preds <- aggregate(cbind(col.Predicted,col.SE,col.lower,col.upper) ~ cnty.index+year, data=preds, FUN=mean)
ext.preds <- aggregate(cbind(ext.Predicted,ext.SE,ext.lower,ext.upper) ~ cnty.index+year, data=preds, FUN=mean)
  

x <- preds[preds$species=="Bovine" & preds$inspection_type=="Scratch",]
det.cattle.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ cnty.index+year, data=x, FUN=mean)

x <- preds[preds$species=="Nilgai" & preds$inspection_type=="Scratch",]
det.nilgai.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ cnty.index+year, data=x, FUN=mean)

x <- preds[preds$species=="Whitetail" & preds$inspection_type=="Scratch",]
det.deer.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ cnty.index+year, data=x, FUN=mean)

x <- preds[preds$species=="Equine" & preds$inspection_type=="Scratch",]
det.equine.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ cnty.index+year, data=x, FUN=mean)

plyr::count(preds[,c("species","inspection_type")])

#-- Generate Plots
  
pdf(paste0(figure.path,"County.Detection.Probability",Sys.Date(),".pdf"), width=11,height=8)
  plot.county(dat=det.cattle.preds, col.plt="det", grp.var="cnty.index")
dev.off()
  
pdf(paste0(figure.path,"County.Extinction.Rate",Sys.Date(),".pdf"), width=11,height=8)
    plot.county(dat=ext.preds, col.plt="ext", grp.var="cnty.index")
dev.off()
  
pdf(paste0(figure.path,"County.Colinization.Rate",Sys.Date(),".pdf"), width=11,height=8)
    plot.county(dat=col.preds, col.plt="col", grp.var="cnty.index")
dev.off()



#-- PLOT Colinization Versus Extinction




det.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ cnty.index+year, data=preds, FUN=mean)
col.preds <- aggregate(cbind(col.Predicted,col.SE,col.lower,col.upper) ~ cnty.index+year, data=preds, FUN=mean)
ext.preds <- aggregate(cbind(ext.Predicted,ext.SE,ext.lower,ext.upper) ~ cnty.index+year, data=preds, FUN=mean)





grp.var="cnty.index"


dat<-col.preds

vec<-unique(dat[,grp.var])

year.vec<-seq(min(dat$year),max(dat$year),1)


par(mfrow=c(4,4))


for(i in 1:length(vec)){
  
  
  col.plt="ext"
  dat<-ext.preds
  
  col.val<-"#238b45"
  
  tmp.plt<-dat[dat[,grp.var] == vec[i],]
  
  plot(x=tmp.plt$year, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, axes=FALSE, ylab="", xlab="",
       ylim=c(0,1), xlim=c(year.vec[1],year.vec[length(year.vec)]), col=col.val)
  arrows(tmp.plt$year,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0.03, col=col.val)
  
  
  col.plt="col"
  dat<-col.preds
  tmp.plt<-dat[dat[,grp.var] == vec[i],]
  
  col.val<-"#a50f15"
  
  points(x=tmp.plt$year, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0.03, col=col.val)
  
  
  
  axis(side=1)
  axis(side=2,las=2)
  
  mtext(side=3,line=1,vec[i],adj=0)
 
  
  #abline(h=equil.occ,col="red")
  
  
   
}





#-- END




#-- Plot Detection All Types Counties

pdf(paste0(figure.path,"County.Detection.Probability.All.Species.",Sys.Date(),".pdf"), width=11,height=10)

dat<-det.preds
dat2<-det.cattle.preds
dat3<-det.nilgai.preds
dat4<-det.deer.preds
dat5<-det.equine.preds

vec<-unique(dat[,grp.var])
vec<-sort(vec)
  
year.vec<-seq(min(dat$year),max(dat$year),1)

col.plt="det"

par(mfrow=c(3,4), mar=c(4,5,3,1))

for(i in 1:length(vec)){
  
  x.lim<-c(year.vec[1]-1,year.vec[length(year.vec)])
  
  tmp.plt<-dat[dat[,grp.var] == vec[i],]
  col.val="#525252"
  
  plot(x=tmp.plt$year, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, axes=FALSE, ylab="", xlab="",
       ylim=c(0,1), xlim=x.lim, col=col.val)
  arrows(tmp.plt$year,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  abline(h=0,col="light gray")
  
  tmp.plt<-dat2[dat2[,grp.var] == vec[i],]
  col.val="#cc4c02"
  
  points(x=tmp.plt$year+.2, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year+.2,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year+.2,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  tmp.plt<-dat3[dat3[,grp.var] == vec[i],]
  col.val="#1d91c0"
  
  points(x=tmp.plt$year-.2, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year-.2,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year-.2,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  tmp.plt<-dat4[dat4[,grp.var] == vec[i],]
  col.val="#253494"
  
  points(x=tmp.plt$year-.4, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year-.4,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year-.4,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  tmp.plt<-dat5[dat5[,grp.var] == vec[i],]
  col.val="black"
  
  points(x=tmp.plt$year+.4, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year+.4,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year+.4,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  axis(side=1)
  axis(side=2,las=2)
  
  mtext(side=3,line=1,vec[i],adj=0)
  
  if(i %in% c(1,5,9)){
    mtext(side=2,line=2.5,"Detection Probability")
  }
  
  if(i==1){
  text(x=c(2013,2013,2013,2013), y=c(1,.94,.86,.8,.75),
       c("Cattle","Nilgai","Whitetail","Equine","All methods/species"),
       col=c("#cc4c02","#1d91c0","#253494","black","#525252"),pos=4, font=2)
  }
  
}#END Loop

dev.off()









x <- preds[preds$species=="Bovine" & preds$inspection_type=="Scratch",]
det.cattle.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ year, data=x, FUN=mean)

x <- preds[preds$species=="Nilgai" & preds$inspection_type=="Scratch",]
det.nilgai.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ year, data=x, FUN=mean)

x <- preds[preds$species=="Whitetail" & preds$inspection_type=="Scratch",]
det.deer.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ year, data=x, FUN=mean)

x <- preds[preds$species=="Equine" & preds$inspection_type=="Scratch",]
det.equine.preds <- aggregate(cbind(det.Predicted,det.SE,det.lower,det.upper) ~ year, data=x, FUN=mean)




pdf(paste0(figure.path,"Detection.Probability.Scratch.All.Species.",Sys.Date(),".pdf"), width=11,height=10)

dat2<-det.cattle.preds
dat3<-det.nilgai.preds
dat4<-det.deer.preds
dat5<-det.equine.preds

vec<-unique(dat[,grp.var])
vec<-sort(vec)

year.vec<-seq(min(dat$year),max(dat$year),1)

col.plt="det"

  x.lim<-c(year.vec[1]-1,year.vec[length(year.vec)])
  
  tmp.plt<-dat2
  col.val="#cc4c02"
  
  plot(x=tmp.plt$year, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, axes=FALSE, ylab="", xlab="",
       ylim=c(0,1), xlim=x.lim, col=col.val)
  arrows(tmp.plt$year,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  abline(h=0,col="light gray")
  
  tmp.plt<-dat3
  col.val="#1d91c0"
  
  points(x=tmp.plt$year-.1, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year-.1,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year-.1,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  tmp.plt<-dat4
  col.val="#253494"
  
  points(x=tmp.plt$year-.2, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year-.2,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year-.2,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  tmp.plt<-dat5
  col.val="black"
  
  points(x=tmp.plt$year+.1, y=tmp.plt[,paste0(col.plt,".Predicted")], pch=19, col=col.val)
  arrows(tmp.plt$year+.1,
         tmp.plt[,paste0(col.plt,".lower")], tmp.plt$year+.1,
         tmp.plt[,paste0(col.plt,".upper")], code=3, angle=90, length=0, col=col.val)
  
  axis(side=1)
  axis(side=2,las=2)
  
  mtext(side=3,line=1,"a) Detection probability (pasture level) for scratched animals",adj=0,cex=1.1)
  
  
  mtext(side=2,line=2.5,"Detection Probability")
  
  
  text(x=c(2013,2013,2013,2013), y=c(1,.94,.86,.8),
         c("Cattle","Nilgai","Whitetail","Equine"),
         col=c("#cc4c02","#1d91c0","#253494","#525252"),pos=4, font=2)
  
dev.off()

#---- END Plot ----








#---- Plot Changes in Occupancy and Detection Through Time ----

pdf(paste0(figure.path,"Occupancy.Detection.Probability",Sys.Date(),".pdf"), width=11,height=8)

op <- par(mfrow=c(2,1), mai=c(0.6, 0.6, 0.6, 0.1))

vec<-seq(min(year.vec),(min(year.vec)+nrow(occ.est))-1,1)

plot(vec, occ.est$projected, pch=19, axes=FALSE, ylab="", xlab="",
     ylim=c(0,1), col=4)

abline(h=equil.occ, col="#67000d")
text(y=equil.occ-.03, x=vec[1], label="Equilibrium Occupancy", col="#67000d", pos=4)


arrows(vec, occ.est$l95ci, vec, occ.est$u95ci, code=3, angle=90, length=0.03, col=4)

axis(side=1)
axis(side=2,las=2)

mtext(side=2,line=3,"Probability")
mtext(side=3,line=0,"a) Occupancy probability (pasture scale)",adj=0)

#abline(h=equil.occ,col="gray")

plot(year.vec, E.det$Predicted, pch=19, axes=FALSE, ylab="", xlab="",
     ylim=c(0,1), col=4)
arrows(year.vec, E.det$lower, year.vec, E.det$upper, code=3, angle=90, length=0.03, col=4)

axis(side=1)
axis(side=2,las=2)

mtext(side=2,line=3,"Probability")
mtext(side=3,line=0,"b) Detection probability (pasture scale)",adj=0)

#abline(h=equil.occ,col="gray")

dev.off()

#---- END Plot ----








#---- Plot 

pdf(paste0(figure.path,"Occ.Probability.pdf"), width=11,height=8)
plot(occ.est$projected,ylim=c(0.04,.07), pch=16, ylab="Occupancy Probability",xlab="Year")
abline(h=equil.occ,col="red")
text(y=equil.occ, x=1,label="Equilibrium Occupancy", col="red", pos=4)
dev.off()










#---- Generate Latent Proportion of Occupied Sites ----

#--Empirical Bayes Prediction of Latent Proportion of Sites Occupied (PAO)
re <- ranef(mod)
e.bup <- colSums(bup(re, stat="mode"))
ci.vals <- confint(re, level=0.95)
rbind(PAO = c(Estimate = sum(e.bup), rowMeans(colSums(ci.vals))) / nrow(re@post)) 




# Empirical Bayes estimates of number of sites occupied in each year
re <- ranef(mod)
modes <- colSums(bup(re, stat="mode"))
plot(1:7, modes, xlab="Year", ylab="Sites occupied", ylim=c(0, 70))

## Find bootstrap standard errors for smoothed trajectory
fm <- nonparboot(mod, B = 100)  # This takes a while!
fm@smoothed.mean.bsse

## get the trajectory estimates
smoothed(mod)
projected(mod)

## try yearly transition rates
yearlySiteCovs(umf) <- data.frame(year = factor(rep(1:7, numSites(umf))))

(fm.yearly <- colext(psiformula = ~ 1,
                     gammaformula = ~ year,
                     epsilonformula = ~ year,
                     pformula = ~ JulianDate + I(JulianDate^2), umf,
                     control = list(trace=1, maxit=1e4)))


#--Parameter Estimates

backTransform(linearComb(fm2, coefficients = c(1,0,0), type = 'det'))



psi<-backTransform(mod, type="psi")
backTransform(mod, type="det")
backTransform(mod, type="ext")
backTransform(mod, type="col")

confint(backTransform(mod, type="psi"))


mod@projected.mean


#--Confidence intervals for projected occupancy
mod <- nonparboot(mod, B = 5)
occ.est<-cbind.data.frame(projected=projected(mod)[2,], SE=mod@projected.mean.bsse[2,])

#occ.est$projected + qt( c(0.05, 0.95), length(occ.est$projected) - 1) * occ.est$SE


#--Equilibrium occupancy
equil.occ<-backTransform(mod, type="col")@estimate/(backTransform(mod, type="col")@estimate + backTransform(mod, type="ext")@estimate)


#---- Plot 

pdf(paste0(figure.path,"Occ.Probability.pdf"), width=11,height=8)
plot(occ.est$projected,ylim=c(0.04,.07), pch=16, ylab="Occupancy Probability",xlab="Year")
abline(h=equil.occ,col="red")
text(y=equil.occ, x=1,label="Equilibrium Occupancy", col="red", pos=4)
dev.off()

#---- Plot Detection Probability

nd <- data.frame(("sampled"=(seq(1, 1000,1))))
colnames(nd)<-c("qty_inspected.std")
pred.p <- predict(mod, type='det', newdata=nd)
print(predictions<-cbind(pred.p, nd))

pdf(paste0(figure.path,"Detection.Probability.pdf"), width=11,height=8)

plot(1, xlim=c(0,1000), ylim=c(0,1), type="n", axes=T, xlab="Number of Animals Inspected per Visit",
     pch=20, ylab="Detection Probability", 
     cex.lab=1.25, cex.main=1.75)

abline(v=mean(det.dat.num.samp,na.rm=TRUE), col="red")

text(x=mean(det.dat.num.samp,na.rm=TRUE), y=0,
     label="Mean number of\nanimals inspected",pos=4, col="red")

lines(predictions$sampled, predictions$Predicted, col="black", lwd=2)
lines(predictions$sampled, predictions$lower, lty=2, col="black")
lines(predictions$sampled, predictions$upper, lty=2, col="black")

mtext(side=3,adj=0,"a) Functional relationship between detection probability and animals inspected")

dev.off()



#---- PLOT Functional Relationships ----








#---- Estimate Pasture Level Probability of Detection on Successive Surveys

#--Generate data
p.star<-est.prob.detection.repeated.sampling(det.stat, n.surveys=20)

#--Plot
pdf(paste0(figure.path,"FIG.prob.detection.many.surveys.",Sys.Date(),".pdf"),width=8,height=8)

par(mar=c(4,4,4,1))

boxplot(p.star$p.star ~ p.star$x, las=1, ylab="",
        xlab="", outline=FALSE, axes=FALSE,
        ylim=c(0,1))

axis(side=1)
axis(side=2, at=seq(0,1,.1), las=2)

mtext(side=1, line=2.5, "Number of Surveys")
mtext(side=2, line=2.5, "Probability of Detection (P star)")

mtext(side=3, adj=0, line=0, "a) Probability of detection with increasing number of surveys")


abline(h=0.95, lty=2, lwd=2)

dev.off()

#---- END 














#---- Estimate Pasture Level Probability of Detection on Successive Surveys

in.dat<-det.site.mu.all

in.dat2<-det.site.mu.pos

colnames(in.dat)[which(colnames(in.dat)=="det.site.mu")] <- "mu"
colnames(in.dat2)[which(colnames(in.dat2)=="det.site.mu")] <- "mu"

vec<-unique(in.dat$cnty.index)
vec<-vec[is.na(vec)==FALSE]
vec<-vec[vec!="other"]
vec<-sort(vec)

n.surveys=20

#--Plot
pdf(paste0(figure.path,"Detection.Probability.County.",Sys.Date(),".pdf"), width=11,height=8)

mat<-matrix(c(1,2,3,4,
              5,6,7,8,
              9,10,11,12),ncol=4,nrow=3,byrow=TRUE)

layout(mat, widths=c(1,1,1,1), heights=c(1,1,1))

par(mar=c(4,4,4,1))

for(i in 1:length(vec)){

x<-in.dat[in.dat$cnty.index==vec[i],]
  
x<-bootstrap.CI(x$mu)

p.star<-est.prob.detection.repeated.sampling(x, n.surveys=n.surveys)

#plt.dat<-boxplot(p.star$p.star ~ p.star$x, plot=FALSE)
plt.dat<-boxplot(p.star$p.star ~ p.star$x, range=0,plot=FALSE)

plot(y=plt.dat$stats[3,], x=plt.dat$names, ylab="",
        xlab="", axes=FALSE,
        ylim=c(0,1), xlim=c(0,n.surveys),col="white")

col.val<-"#252525"

shade.col<-transparent(orig.col = col.val, trans.val = .25, maxColorValue = 255)

Epi::matshade(x=as.numeric(plt.dat$names),
              y=as.matrix(t(plt.dat$stats[c(3,1,5),])),
              col=col.val, col.shade=shade.col, lwd=2)

if(i==1){
  text(x=c(n.surveys), y=c(.15),
       labels=c("Any Pasture"),
       pos=2, col=col.val)
}


x<-in.dat2[in.dat2$cnty.index==vec[i],]

if(nrow(x)>0){

x<-bootstrap.CI(x$mu)

p.star<-est.prob.detection.repeated.sampling(x, n.surveys=n.surveys)

#plt.dat<-boxplot(p.star$p.star ~ p.star$x, plot=FALSE)
plt.dat<-boxplot(p.star$p.star ~ p.star$x, range=0,plot=FALSE)


col.val<-"#67000d"

shade.col<-transparent(orig.col = col.val, trans.val = .25, maxColorValue = 255)

Epi::matshade(x=as.numeric(plt.dat$names),
              y=as.matrix(t(plt.dat$stats[c(3,1,5),])),
              col=col.val, col.shade=shade.col, lwd=2)

}

if(i==1){
  text(x=c(n.surveys), y=c(.05),
       labels=c("Previously Infested Pasture"),
       pos=2, col=col.val)
}


#--Add 
axis(side=1)
axis(side=2, at=seq(0,1,.1), las=2)

mtext(side=1, line=2.5, "Number of Surveys")
mtext(side=2, line=2.5, "Probability of Detection")
mtext(side=3, adj=0, line=0, vec[i])




abline(h=0.95, lty=2, lwd=2, col="red")
}#END Loop

dev.off()

#---- END 


Epi::matshade(x=as.numeric(plt.dat$names),
              y=as.matrix(t(plt.dat$stats[c(3,1,5),])), col="#525252", col.shade="#d9d9d9")











#---- Plot Detection Probability by Year

#--Generate Predictions for each site/observation
pred.det <- predict(mod, type='det')

#--Merge with observation data
pred.det<-cbind.data.frame(occ.dat@obsCovs, pred.det)

pdf(paste0(figure.path,"FIG.prob.detection.annual.",Sys.Date(),".pdf"),width=8,height=8)

  boxplot(pred.det$Predicted ~ pred.det$year, las=1, ylab="Probability of Detection on Single Inspection", xlab="", outline=FALSE)
  abline(h=0.95, lty=2, lwd=2)

dev.off()

#---- END Plot Detection Probability by Year











#---- Detection Probability by Species

var.plt<-"density.std"
var.mean<-"qty_inspected.std"

#--Set Variable to hold at mean
x.mean.var<-mean(occ.dat@obsCovs[,var.mean],na.rm=TRUE)

#--Set variable of interest
x<-occ.dat@obsCovs[,var.plt]
x.var<-seq(min(x,na.rm=TRUE),max(x,na.rm=TRUE)+3,1)

#--Unscale Values
center <- attr(x,"scaled:center")
scale <- attr(x,"scaled:scale")

unscale<-function(x){t(t(x) * scale + center)}


#--Make Plot

pdf(paste0(figure.path,"FIG.Det.Probability.by.species.density.pdf"), width=8,height=8)

month.vec<-c(1,3,6,9)
var.vec<-c("Bovine","Whitetail","Nilgai","Equine","Other Species")
lab.vec<-c("Bovine","Whitetail","Nilgai","Equine","Other")

for(i in 1:length(var.vec)){
  newData <- data.frame(x.mean=x.mean.var,
                        species = factor(var.vec[i], levels=c("Bovine","Equine","None","Whitetail","Nilgai","Other Species")),
                        var.plt=x.var)
  colnames(newData)[1]<-var.mean
  colnames(newData)[3]<-var.plt
  
  tmp<-predict(mod, type = 'det', newdata = newData, appendData=TRUE)
  
  
  tmp$x<-unscale(tmp[,var.plt])
  
  if(i==1){
    plot(1, xlim=c(-100,max(tmp$x)), ylim=c(0,1), type="n", axes=FALSE, xlab="Number of Animals Inspected",
         pch=20, ylab="Detection Probability", 
         cex.lab=1.25, cex.main=1.75)
  }
  
  lines(tmp$x, tmp$Predicted, col="black", lwd=2)
  
  text(x=0,y=tmp$Predicted[1],labels=lab.vec[i],pos=2)
  
}#END Plot

axis(side=1)
axis(side=2,las=2)

#x<-mean(occ.raw$qty_inspected,na.rm=TRUE)
#text(x=x, y=0,
#     label="Mean number of\nanimals inspected per visit",pos=4, col="red")

#abline(v=x,col="red")

mtext(side=3,adj=0,"a) Functional relationship between detection probability and animal density")

dev.off()

#---- END Detection Probability by Species










#---- PLOT Predicted Pasture Occupancy ----

# Estimates of conditional occupancy distribution at each site and time point
re <- ranef(mod)

#--Draw from the posterior predictive distribution
ppd <- posteriorSamples(re, nsims=100000)

#--Generate means
ppd.mean <- posterior.mean(ppd)
#ppd.cnt<-posterior.cnt.success(ppd)

#--Merge with sites
ppd.mean<-cbind.data.frame(site.dat, ppd.mean)

#--Never Infested
ppd.never <- ppd.mean[ppd.mean$never.pos==1,]
ppd.never<-ppd.never[,(ncol(ppd.never)-6):ncol(ppd.never)]
mean.never<-colMeans(ppd.never)

ci.lower<-apply(ppd.never, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
ci.upper<-apply(ppd.never, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))

mean.never<-cbind.data.frame(mean=mean.never,ci.lower=ci.lower,ci.upper=ci.upper)


#--Infested
ppd.pos <- ppd.mean[ppd.mean$never.pos==0,]
ppd.pos<-ppd.pos[,(ncol(ppd.pos)-6):ncol(ppd.pos)]
mean.pos<-colMeans(ppd.pos)

ci.pos<-t(apply(ppd.pos, MARGIN=2, FUN=function(x) quantile(x, probs=c(.25,.75))))


ci.lower<-apply(ppd.pos, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
ci.upper<-apply(ppd.pos, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))

mean.pos<-cbind.data.frame(mean=mean.pos,ci.lower=ci.lower,ci.upper=ci.upper)






pdf(paste0(figure.path,"Probability.Occupancy.Global.",Sys.Date(),".pdf"), width=8,height=6)

  par(mar=c(3,4,2,1.5))
  
  #-- All pastures
  ppd.subset <- ppd.mean
  ppd.subset<-ppd.subset[,(ncol(ppd.subset)-6):ncol(ppd.subset)]
  subset.mean<-colMeans(ppd.subset)
  
  ci.lower<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
  ci.upper<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))
  
  stat.dat<-cbind.data.frame(mean=subset.mean,ci.lower=ci.lower,ci.upper=ci.upper)
  
  tmp.plt<-stat.dat
  col.val="#525252"
  
  year.vec<-seq(2014,2020,1)
  
  shift.val<-0.1
  
  plot(x=year.vec-shift.val, y=tmp.plt$mean, pch=19, axes=FALSE, ylab="", xlab="",
       ylim=c(0,1), col=col.val)
  
  arrows(year.vec-shift.val,
         tmp.plt[,2], year.vec-shift.val,
         tmp.plt[,3], code=3, angle=90, length=0, col=col.val)
  
  
  #Pastures that have ever been infested
  ppd.subset <- ppd.mean[ppd.mean$never.pos==0,]
  ppd.subset<-ppd.subset[,(ncol(ppd.subset)-6):ncol(ppd.subset)]
  subset.mean<-colMeans(ppd.subset)
  
  ci.lower<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
  ci.upper<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))
  
  stat.dat<-cbind.data.frame(mean=subset.mean,ci.lower=ci.lower,ci.upper=ci.upper)
  
  
  tmp.plt<-stat.dat
  col.val="#cc4c02"
  
  points(x=year.vec+shift.val, y=tmp.plt$mean, pch=19, col=col.val)
  
  arrows(year.vec+shift.val,
         tmp.plt[,2], year.vec+shift.val,
         tmp.plt[,3], code=3, angle=90, length=0, col=col.val)
  
  axis(side=1)
  axis(side=2, las=2)
  
  mtext(side=2, line=2.5,"Probability of Pasture Infestation")
  text(x=c(2014,2014), y=c(0.98,.93), labels=c("Previously Infested Pasture","Any Pasture"),
         col=c("#cc4c02","#525252"), pos=4, font=2)
  text(x=2014.5,y=equil.occ,pos=1,cex=.8,labels="Equilibrium Occupancy",col="red")

  mtext(side=3, adj=0, line=0, "a) Probability a pasture is infested", cex=1.2)
  
  abline(h=equil.occ,col="red")
  
dev.off()

#----- END END -----












#---- Probability Occupancy County Level ----

pdf(paste0(figure.path,"Probability.Occupancy.County.",Sys.Date(),".pdf"), width=11,height=10)



ppd.mean[ppd.mean$cnty.index %in% "La Salle","cnty.index"] <- "other"

vec<-sort(unique(ppd.mean$cnty.index))

mat<-matrix(c(1,2,3,4,
              5,6,7,8,
              9,10,11,12), byrow=TRUE, ncol=4, nrow=3)

layout(mat, heights=c(1,1,1,1), widths=c(1,1,1,1))

for(i in 1:length(vec)){

par(mar=c(3,4,2,1.5))
  
#-- All pastures
ppd.subset <- ppd.mean[ppd.mean$cnty.index==vec[i],]
ppd.subset<-ppd.subset[,(ncol(ppd.subset)-6):ncol(ppd.subset)]
subset.mean<-colMeans(ppd.subset)

ci.lower<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
ci.upper<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))

stat.dat<-cbind.data.frame(mean=subset.mean,ci.lower=ci.lower,ci.upper=ci.upper)

tmp.plt<-stat.dat
col.val="#525252"

year.vec<-seq(2014,2020,1)

shift.val<-0.1

plot(x=year.vec-shift.val, y=tmp.plt$mean, pch=19, axes=FALSE, ylab="", xlab="",
     ylim=c(0,1), col=col.val)

arrows(year.vec-shift.val,
       tmp.plt[,2], year.vec-shift.val,
       tmp.plt[,3], code=3, angle=90, length=0, col=col.val)


#Pastures that have ever been infested
ppd.subset <- ppd.mean[ppd.mean$cnty.index==vec[i] & ppd.mean$never.pos==0,]
ppd.subset<-ppd.subset[,(ncol(ppd.subset)-6):ncol(ppd.subset)]
subset.mean<-colMeans(ppd.subset)

ci.lower<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.025), length(x) - 1) * std.error(x))
ci.upper<-apply(ppd.subset, MARGIN=2, FUN=function(x) mean(x) + qt( c(0.975), length(x) - 1) * std.error(x))

stat.dat<-cbind.data.frame(mean=subset.mean,ci.lower=ci.lower,ci.upper=ci.upper)


tmp.plt<-stat.dat
col.val="#cc4c02"

points(x=year.vec+shift.val, y=tmp.plt$mean, pch=19, col=col.val)

arrows(year.vec+shift.val,
       tmp.plt[,2], year.vec+shift.val,
       tmp.plt[,3], code=3, angle=90, length=0, col=col.val)

axis(side=1)
axis(side=2, las=2)

if(i %in% c(1,5,9)) {mtext(side=2, line=2.5,"Probability of Pasture Infestation", cex=.9)}

if(i == 1){
  text(x=c(2014,2014), y=c(0.98,.93), labels=c("Previously Infested Pasture","Any Pasture"),
       col=c("#cc4c02","#525252"), pos=4, font=2)
  text(x=2015.2,y=equil.occ,pos=1,cex=.8,labels="Equilibrium Occupancy",col="red")
}

mtext(side=3, adj=0, line=0, vec[i])

abline(h=equil.occ,col="red")

}#END Loop

dev.off()

#----- END END -----




EBUP <- bup(re, stat="mode")






# Best Unbiased Predictors
tmp<-bup(re, stat="mean")           # Posterior mean
bup(re, stat="mode")           # Posterior mode

bup.ci<-unmarked::confint(re, level=0.9, type="state") # 90% CI


bup.ci[100,,]




plot(re, subset=site %in% c(1:10), layout=c(5, 2))

plot(re, subset=year, layout=c(5, 2))


























#-------------------
#---- FUNCTIONS ----


bootstrap.CI <- function(in.dat){
  sampmean <- function(y,i) mean(y[i])
  bootmean <- boot(data=in.dat,statistic=sampmean,R=10000)
  out.ci<-boot.ci(bootmean, conf=.95,type=c("norm"))
  
  out<-list(mu=bootmean$t,summary=out.ci)
  return(out)
}#END Function


antilogit <- function(x) { exp(x) / (1 + exp(x) ) }




#-- FNC: Estimate Pasture Level Probability of Detection on Successive Surveys

est.prob.detection.repeated.sampling <- function(det.stat, n.surveys){
  
  #--Controls
  n.iter<-length(det.stat$mu)
  n.surveys<-20
  
  #--Simulate for testing
  #psi<- rbeta(n=n.iter, shape1=2, shape2=5)
  psi <- det.stat$mu
  
  #--Storage
  p.star<-array(NA, dim=c(n.iter, n.surveys))
  
  #--Set up x data (surveys)
  for(i in 1:n.surveys){
    y<-rep(i,n.iter)
    if(i==1){x<-y}
    if(i>1){x<-cbind(x,y)}
  }
  
  #--Estimate Prob Detection
  for(i in 1:n.iter){
    for(j in 1:n.surveys){
      p.star[i,j]<- 1 - (1 - psi[i])^j
    }#END j
  }#END i
  
  out<-list(x=x, p.star=p.star)
  
  return(out)
  
}#END Function



#---- FNC - Plot regression coeff 

plot.betas<-function(plt.dat){
  
  x.labs<-pretty(as.vector(as.matrix(plt.dat)))
  
  x.min<-min(x.labs)
  x.max<-max(x.labs)
  
  plt.dat$name<-rownames(plt.dat)
  plt.dat$y<-seq(1,nrow(plt.dat),1)
  
  plot(x=plt.dat$estimate, y=plt.dat$y,
       xlim=c(x.min,x.max),pch=19,
       axes=FALSE, xlab="",ylab="")
  
  abline(v=0,col="gray")
  
  segments(x0=plt.dat$`0.025`,x1=plt.dat$`0.975`,y0=plt.dat$y,y1=plt.dat$y,lwd=1.5)
  points(x=plt.dat$estimate, y=plt.dat$y,pch=19)
  
  axis(side=1)
  axis(side=2, at=plt.dat$y, labels=plt.dat$name, las=2)
  
}#END Plot Betas









generate.annual.subsets<-function(occ.raw, year.val){

  occ.raw<-occ.raw[occ.raw$year==year.val,]


  col.names<-colnames(occ.raw)
  col.names<-col.names[col.names %!in% c("pasture_name","jdate","y","year")]
  
  occ.raw<-occ.raw[,c("year","pasture_name","jdate","y",col.names)]
  
  #--Add Site Year
  occ.raw$site.year<-occ.raw$year
  
  #--Make unmarked model object (takes a a long time)
  occ.raw<-formatMult(occ.raw)
  
  
  #Observations
  col.names<-c("qty_inspected","qty_infested","apr.infest.prev","prop.new.animals", "density","med.dist.pos","sd.dist.pos")
  for(i in 1:length(col.names)){
    obsCovs(occ.raw)[,col.names[i]]<-as.numeric(as.character(obsCovs(occ.raw)[,col.names[i]]))
    obsCovs(occ.raw)[,paste0(col.names[i],".std")] <- scale(obsCovs(occ.raw)[,col.names[i]], center = TRUE, scale = TRUE)
  }
  
  #Sites
  col.names<-c("pasture_qty_acres","mx.distance","site.med.dist.pos","site.sd.dist.pos","qty_added")
  for(i in 1:length(col.names)){
    siteCovs(occ.raw)[,col.names[i]]<-as.numeric(as.character(siteCovs(occ.raw)[,col.names[i]]))
    siteCovs(occ.raw)[,paste0(col.names[i],".std")] <- scale(siteCovs(occ.raw)[,col.names[i]], center = TRUE, scale = TRUE)
  }

return(occ.raw)
}#END






