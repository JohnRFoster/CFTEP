#-------------
#
# Supporting Functions
#
# By Ryan Miller
#
#-------------

#--FNC - Function to convert from wide to long
wide.to.long <- function(in.dat, mast.col, pattern.str) {
  #Reduce to columns of interest
  col.names <- colnames(in.dat)
  name.vec <- col.names[grep(x = col.names, pattern = pattern.str)]
  tmp <- in.dat[, c(mast.col, name.vec)]

  #Reshape from wide to long
  tmp.long <- gather(
    in.dat[, c(mast.col, name.vec)],
    key = Species,
    value = value,
    name.vec[1]:name.vec[length(name.vec)],
    factor_key = FALSE
  )

  #Fix Species names
  tmp.long$Species <- str_remove(
    string = tmp.long$Species,
    pattern = paste0(pattern.str, ".")
  )

  #Rename column
  colnames(tmp.long)[ncol(tmp.long)] <- pattern.str

  return(tmp.long)
} #END Function


#--FNC - Rename column names
alter.column.names <- function(in.dat, pattern.str) {
  #Get Column names
  col.names <- colnames(in.dat)
  name.vec <- col.names[grep(x = col.names, pattern = pattern.str)]

  #Change Names
  name.vec <- str_remove(string = name.vec, pattern = paste0("_", pattern.str))
  name.vec <- paste0(pattern.str, ".", name.vec)

  #Assign to DF
  col.names[grep(x = col.names, pattern = paste0("_", pattern.str))] <- name.vec

  return(col.names)
} #END Function


#FNC - Rename columns
rename.columns <- function(in.dat, org.col, new.col) {
  length(new.col) == length(org.col)

  for (i in 1:length(org.col)) {
    colnames(in.dat)[which(colnames(in.dat) == org.col[i])] <- new.col[i]
  } #END Loop

  return(in.dat)
} #END Function


#FNC - Return only columns that have at least 1 data value
completeFun <- function(in.dat, col.names) {
  for (i in 1:length(col.names)) {
    vec <- complete.cases(in.dat[, col.names[i]])

    vec[vec == TRUE] <- 1
    vec[vec == FALSE] <- 0

    if (i == 1) {
      out <- as.data.frame(vec)
    }
    if (i > 1) {
      out <- cbind.data.frame(out, vec)
    }
  }
  vec <- rowSums(out)

  vec <- vec > 0

  in.dat <- in.dat[vec == TRUE, ]

  return(in.dat)
} #END Function


#-- FNC - Median distance to infested pastures in each year
median.distance.to.infested <- function(dat, type) {
  if (type == "year") {
    dat$i.loop <- dat$year
    i.vec <- unique(dat$i.loop)
  }
  if (type == "season") {
    dat$i.loop <- paste0(dat$year, dat$season)
    i.vec <- unique(dat$i.loop)
  }

  pas.vec <- unique(dat$pasture_name)

  pb <- txtProgressBar(min = 0, max = length(i.vec), style = 3)

  #Loop over year
  for (i in 1:length(i.vec)) {
    infest.dat <- dat[
      dat$y == 1 & dat$i.loop == i.vec[i],
      c("pasture_name", "pasture_longitude", "pasture_latitude")
    ]
    infest.dat <- unique(infest.dat)

    #Loop over pastures
    for (j in 1:length(pas.vec)) {
      if (
        nrow(dat[dat$pasture_name == pas.vec[j] & dat$i.loop == i.vec[i], ]) !=
          0
      ) {
        xy <- unique(dat[
          dat$pasture_name == pas.vec[j],
          c("pasture_longitude", "pasture_latitude")
        ])

        #Calc distances
        pnt.dist <- pointDistance(
          p1 = xy,
          p2 = infest.dat[
            infest.dat$pasture_name != pas.vec[j],
            c("pasture_longitude", "pasture_latitude")
          ],
          lonlat = TRUE
        )

        #Identify outliers and remove
        tmp <- boxplot(pnt.dist, plot = FALSE)$out

        if (length(tmp) != 0) {
          cut.pnt <- min(tmp[tmp > median(pnt.dist)]) - 0.01
          val.med <- median(pnt.dist[pnt.dist < cut.pnt])
          val.sd <- sd(pnt.dist[pnt.dist < cut.pnt])
        }
        if (length(tmp) == 0) {
          val.med <- median(pnt.dist)
          val.sd <- sd(pnt.dist)
        }

        #Assign median values
        dat[
          dat$pasture_name == pas.vec[j] & dat$i.loop == i.vec[i],
          "med.dist.pos"
        ] <- val.med
        dat[
          dat$pasture_name == pas.vec[j] & dat$i.loop == i.vec[i],
          "sd.dist.pos"
        ] <- val.sd
      } #END Logical
    } #END Loop

    setTxtProgressBar(pb, i)
  } #END Loop

  dat <- dat[, colnames(dat) %!in% c("i.loop")]

  return(dat)
} #-- END Distance to infested pastures


#-- Reduce to locations by year

generate.lag.rate <- function(dat) {
  year.vec <- unique(dat$year)

  pb <- txtProgressBar(min = 0, max = length(year.vec), style = 3)

  for (i in 1:length(year.vec)) {
    # Subset data by year
    tmp <- dat[dat$year == year.vec[i], ]

    # Generate unique values
    tmp <- unique(tmp[, c(
      "pasture_name",
      "pasture_latitude",
      "pasture_longitude",
      "year",
      "y"
    )])

    tmp <- aggregate(
      y ~ pasture_name + pasture_latitude + pasture_longitude + year,
      data = tmp,
      FUN = sum
    )

    # Address pastures with same lat/lon
    tmp$pasture_latitude <- (runif(n = nrow(tmp), min = 1, max = 100) /
      1000000) +
      tmp$pasture_latitude
    tmp$pasture_longitude <- (runif(n = nrow(tmp), min = 1, max = 100) /
      1000000) +
      tmp$pasture_longitude

    # Generate polygons
    df.sp <- tmp

    names(df.sp)[1] <- "Pasture_Name"

    coordinates(df.sp) <- ~ pasture_longitude + pasture_latitude

    vp <- voronoipolygons(df.sp)

    names(vp) <- "pasture_name"

    vp <- merge(vp, tmp, by = "pasture_name")

    #-- Add Spatial Relationship to Adjacents

    #--Create Neighbor List
    require(spdep)
    W.nb <- poly2nb(vp, row.names = rownames(vp@data))
    W.list <- nb2listw(W.nb, style = "W")

    #--END Make Neighbor List

    rownames(vp@data) <- vp$pasture_name

    vp$num.adj.pas <- NA
    vp$lag.count <- NA
    vp$lag.rate <- NA

    for (j in 1:length(W.list$neighbours)) {
      row.vec <- W.list$neighbours[[j]]
      tmp <- vp@data[row.vec, ]

      vp@data[j, ]$lag.rate <- sum(tmp$y) / length(tmp$y)
      vp@data[j, ]$lag.count <- sum(tmp$y)
      vp@data[j, ]$num.adj.pas <- length(tmp$y)
    } #END Loop

    if (i == 1) {
      out <- vp@data
    }
    if (i > 1) {
      out <- rbind.data.frame(out, vp@data)
    }

    setTxtProgressBar(pb, i)
  } #END Loop

  return(out)
} #END Function


#-- FNC - Make voronoi polygons

voronoipolygons <- function(layer) {
  require(deldir)
  crds = layer@coords
  z = deldir(crds[, 1], crds[, 2])
  w = tile.list(z)
  polys = vector(mode = 'list', length = length(w))
  require(sp)
  for (i in seq(along = polys)) {
    pcrds = cbind(w[[i]]$x, w[[i]]$y)
    pcrds = rbind(pcrds, pcrds[1, ])
    polys[[i]] = Polygons(list(Polygon(pcrds)), ID = as.character(i))
  }
  SP = SpatialPolygons(polys)
  voronoi = SpatialPolygonsDataFrame(
    SP,
    data = data.frame(
      Pasture_Name = layer$Pasture_Name,
      row.names = sapply(slot(SP, 'polygons'), function(x) slot(x, 'ID'))
    )
  )
  return(voronoi)
} #END Function


#-- Reduce to locations by year

generate.adj.pasture.dist <- function(dat) {
  year.vec <- unique(dat$year)
  year.vec <- year.vec[is.na(year.vec) == FALSE]

  pb <- txtProgressBar(min = 0, max = length(year.vec), style = 3)

  for (i in 1:length(year.vec)) {
    # Subset data by year
    tmp <- dat[dat$year == year.vec[i], ]

    # Generate unique values
    tmp <- unique(tmp[, c(
      "pasture_name",
      "pasture_latitude",
      "pasture_longitude",
      "year",
      "y"
    )])

    tmp <- aggregate(
      y ~ pasture_name + pasture_latitude + pasture_longitude + year,
      data = tmp,
      FUN = sum
    )

    # Address pastures with same lat/lon
    tmp$pasture_latitude <- (runif(n = nrow(tmp), min = 1, max = 100) /
      1000000) +
      tmp$pasture_latitude
    tmp$pasture_longitude <- (runif(n = nrow(tmp), min = 1, max = 100) /
      1000000) +
      tmp$pasture_longitude

    # Generate polygons
    df.sp <- tmp

    names(df.sp)[1] <- "Pasture_Name"

    coordinates(df.sp) <- ~ pasture_longitude + pasture_latitude

    vp <- voronoipolygons(df.sp)

    names(vp) <- "pasture_name"

    vp <- merge(vp, tmp, by = "pasture_name", all.x = TRUE)

    #-- Add Spatial Relationship to Adjacents

    #--Create Neighbor List
    require(spdep)
    W.nb <- poly2nb(vp, row.names = rownames(vp@data))
    W.list <- nb2listw(W.nb, style = "W")

    #--END Make Neighbor List

    #-- Add columns
    rownames(vp@data) <- vp$pasture_name

    vp$mean.pas.dist <- NA
    vp$sd.pas.dist <- NA

    #-- Generate distances
    for (j in 1:length(W.list$neighbours)) {
      row.vec <- W.list$neighbours[[j]]
      tmp <- vp@data[row.vec, ]

      xy <- unique(vp@data[j, c("pasture_longitude", "pasture_latitude")])
      to.xy <- vp@data[row.vec, c("pasture_longitude", "pasture_latitude")]

      pnt.dist <- pointDistance(p1 = xy, p2 = to.xy, lonlat = TRUE)

      vp@data[j, ]$mean.pas.dist <- mean(pnt.dist)
      vp@data[j, ]$sd.pas.dist <- sd(pnt.dist)
    } #END Loop

    if (i == 1) {
      out <- vp@data
    }
    if (i > 1) {
      out <- rbind.data.frame(out, vp@data)
    }

    setTxtProgressBar(pb, i)
  } #END Loop

  return(out)
} #END Function


#-- FNC - Generate time lagged counts and rates for neighborhood infested

generate.lag.1.rate <- function(dat) {
  # Generate Lag Data
  lag.dat <- generate.lag.rate(dat)

  # Unique Pasture Names
  pas.names <- unique(lag.dat$pasture_name)

  pb <- txtProgressBar(min = 0, max = length(pas.names), style = 3)

  for (i in 1:length(pas.names)) {
    year.vec <- lag.dat[lag.dat$pasture_name %in% pas.names[i], "year"]

    year.vec <- sort(year.vec)

    lag.dat$lag.1.count <- NA
    lag.dat$lag.1.rate <- NA

    if (length(year.vec) > 1) {
      for (j in 2:length(year.vec)) {
        cnt <- lag.dat[
          lag.dat$pasture_name %in%
            pas.names[i] &
            lag.dat$year == year.vec[j - 1],
          "lag.count"
        ]
        rate <- lag.dat[
          lag.dat$pasture_name %in%
            pas.names[i] &
            lag.dat$year == year.vec[j - 1],
          "lag.rate"
        ]

        lag.dat[
          lag.dat$pasture_name %in% pas.names[i] & lag.dat$year == year.vec[j],
          "lag.1.count"
        ] <- cnt
        lag.dat[
          lag.dat$pasture_name %in% pas.names[i] & lag.dat$year == year.vec[j],
          "lag.1.rate"
        ] <- rate
      } #END Loop over years
    } #END Logical

    tmp <- lag.dat[lag.dat$pasture_name %in% pas.names[i], ]

    if (i == 1) {
      out <- tmp
    }
    if (i > 1) {
      out <- rbind.data.frame(out, tmp)
    }

    setTxtProgressBar(pb, j)
  } #END Loop over pastures

  return(out)
} #END Function


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

plot.betas <- function(plt.dat, names.vec, vars.vec) {
  plt.dat <- plt.dat[rownames(plt.dat) %in% vars.vec, ]

  x.labs <- pretty(as.vector(as.matrix(plt.dat)))

  x.min <- min(x.labs)
  x.max <- max(x.labs)

  if (x.max < 1) {
    x.max <- 1
  }

  plt.dat$name <- names.vec
  plt.dat$y <- seq(1, nrow(plt.dat), 1)

  plt.dat <- plt.dat[order(plt.dat$estimate), ]
  plt.dat[plt.dat$name != "(Intercept)", "y"] <- seq(2, nrow(plt.dat), 1)
  plt.dat[plt.dat$name == "(Intercept)", "y"] <- 1

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


#-- FNC - Generate Neighbor Mean

neighbor.mean <- function(in.dat, col.name) {
  pb <- txtProgressBar(min = 0, max = length(year.vec), style = 3)

  # Generate unique values
  tmp <- unique(in.dat[, c(
    "Pasture_Name",
    "Pasture_Latitude",
    "Pasture_Longitude"
  )])

  # Address pastures with same lat/lon
  tmp$Pasture_Latitude <- (runif(n = nrow(tmp), min = 1, max = 100) / 1000000) +
    tmp$Pasture_Latitude
  tmp$Pasture_Longitude <- (runif(n = nrow(tmp), min = 1, max = 100) /
    1000000) +
    tmp$Pasture_Longitude

  # Generate polygons
  df.sp <- tmp

  coordinates(df.sp) <- ~ Pasture_Longitude + Pasture_Latitude

  vp <- voronoipolygons(df.sp)

  vp <- merge(vp, tmp, by = "Pasture_Name", all.x = TRUE)

  #-- Add Spatial Relationship to Adjacents

  #--Create Neighbor List
  require(spdep)
  W.nb <- poly2nb(vp, row.names = rownames(vp@data))
  W.list <- nb2listw(W.nb, style = "W")

  #--END Make Neighbor List

  #-- Add columns
  rownames(vp@data) <- vp$Pasture_Name

  pas.vec <- unique(in.dat[is.na(in.dat[, col.name]) == TRUE, "Pasture_Name"])

  #Loop over pastures
  for (i in 1:length(pas.vec)) {
    row.vec <- W.list$neighbours[[i]]

    vec <- vp@data[row.vec, "Pasture_Name"]

    vec <- vec[vec != pas.vec[i]]

    tmp <- in.dat[in.dat$Pasture_Name %in% vec, c("Pasture_Name", col.name)]
    tmp <- unique(tmp)

    val <- mean(tmp[, col.name], na.rm = TRUE)

    in.dat[in.dat$Pasture_Name == pas.vec[i], col.name] <- val
  } #END Loop

  setTxtProgressBar(pb, i)

  return(in.dat)
} #END Function


#-- FNC - Plot Function for County Level Estimates
plot.county <- function(dat, col.plt = "det", grp.var = "cnty.index") {
  vec <- unique(dat[, grp.var])

  year.vec <- seq(min(dat$year), max(dat$year), 1)

  par(mfrow = c(4, 4))

  for (i in 1:length(vec)) {
    tmp.plt <- dat[dat[, grp.var] == vec[i], ]

    plot(
      x = tmp.plt$year,
      y = tmp.plt[, paste0(col.plt, ".Predicted")],
      pch = 19,
      axes = FALSE,
      ylab = "",
      xlab = "",
      ylim = c(0, 1),
      xlim = c(year.vec[1], year.vec[length(year.vec)]),
      col = 4
    )
    arrows(
      tmp.plt$year,
      tmp.plt[, paste0(col.plt, ".lower")],
      tmp.plt$year,
      tmp.plt[, paste0(col.plt, ".upper")],
      code = 3,
      angle = 90,
      length = 0.03,
      col = 4
    )

    axis(side = 1)
    axis(side = 2, las = 2)

    mtext(side = 3, line = 1, vec[i], adj = 0)
  } #END Loop
} #END Function


#-- FNC - Function to generate posterior means for each site/year

posterior.mean <- function(ppd) {
  out <- matrix(ncol = dim(ppd@samples)[2], nrow = dim(ppd@samples)[1])

  for (i in 1:dim(out)[1]) {
    for (j in 1:dim(out)[2]) {
      out[i, j] <- sum(ppd@samples[i, j, ]) / length(ppd@samples[i, j, ])
    } #END Columns
  } #END Rows

  return(out)
} #END FNC


#-- FNC - Function to generate posterior means for each site/year

posterior.cnt.success <- function(ppd) {
  out <- matrix(ncol = dim(ppd@samples)[2], nrow = dim(ppd@samples)[1])

  for (i in 1:dim(out)[1]) {
    for (j in 1:dim(out)[2]) {
      out[i, j] <- sum(ppd@samples[i, j, ])
    } #END Columns
  } #END Rows

  return(out)
} #END FNC
