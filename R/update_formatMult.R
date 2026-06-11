library(R.utils)
library(reshape)

formatMult <- function(df.in) {
  years <- sort(unique(df.in[[1]]))
  nY <- length(years)
  df.obs <- list()
  nsamp <- numeric()
  maxsamp <- max(table(df.in[[1]], df.in[[2]])) # the maximum samples/yr
  for (t in 1:nY) {
    df.t <- df.in[df.in[[1]] == years[t], ] # subset for current year
    df.t <- df.t[, -1] # remove year column
    df.t <- unmarked:::dateToObs(df.t)
    nsamp <- max(df.t$obsNum)
    if (nsamp < maxsamp) {
      newrows <- df.t[1:(maxsamp - nsamp), ] # just a placeholder
      newrows[, "obsNum"] <- ((nsamp + 1):maxsamp)
      newrows[, 3:(ncol(df.t) - 1)] <- NA
      df.t <- rbind(df.t, newrows)
    }
    df.obs <- rbind(df.obs, cbind(year = years[t], df.t))
  }
  dfnm <- colnames(df.obs)
  nV <- length(dfnm) - 1 # last variable is obsNum

  ### Identify variables that are not factors
  # Include julian date/visit in search as it is added back in later
  fac <- sapply(df.obs[, c(3, 5:nV)], is.factor)
  nonfac <- names(df.obs[, c(3, 5:nV)])[!fac]

  # create y matrix using reshape
  expr <- substitute(
    recast(
      df.obs,
      var1 ~ year + obsNum + variable,
      id.var = c(dfnm[2], "year", "obsNum"),
      measure.var = dfnm[4]
    ),
    list(var1 = as.name(dfnm[2]))
  )
  y <- as.matrix(eval(expr)[, -1])

  # create obsdata with reshape
  # include date (3rd col) and other measured vars
  expr <- substitute(
    recast(
      df.obs,
      newvar ~ year + obsNum ~ variable,
      id.var = c(dfnm[2], "year", "obsNum"),
      measure.var = dfnm[c(3, 5:nV)]
    ),
    list(newvar = as.name(dfnm[2]))
  )
  obsvars <- eval(expr)

  rownames(y) <- dimnames(obsvars)[[1]]
  colnames(y) <- dimnames(obsvars)[[2]]
  y <- as.matrix(y)

  obsvars.list <- arrToList(obsvars)

  # Return any non-factors to the correct mode
  if (length(nonfac) >= 1) {
    modes <- apply(df.obs[, nonfac], 2, mode)
    for (i in 1:length(nonfac)) {
      mode(obsvars.list[[nonfac[i]]]) <- modes[i]
    }
  }

  obsvars.list <- lapply(obsvars.list, function(x) as.vector(t(x)))
  obsvars.df <- as.data.frame(obsvars.list)

  ## check for siteCovs
  obsNum <- ncol(y)
  M <- nrow(y)
  site.inds <- matrix(1:(M * obsNum), M, obsNum, byrow = TRUE)
  siteCovs <- sapply(obsvars.df, function(x) {
    obsmat <- matrix(x, M, obsNum, byrow = TRUE)
    l.u <- apply(obsmat, 1, function(y) {
      row.u <- unique(y)
      length(row.u[!is.na(row.u)])
    })
    ## if there are 0 or 1 unique vals per row, we have a sitecov
    if (all(l.u %in% 0:1)) {
      u <- apply(obsmat, 1, function(y) {
        row.u <- unique(y)
        ## only remove NAs if there are some non-NAs.
        if (!all(is.na(row.u))) {
          row.u <- row.u[!is.na(row.u)]
        }
        row.u
      })
      u
    }
  })
  siteCovs <- as.data.frame(siteCovs[!sapply(siteCovs, is.null)])
  if (nrow(siteCovs) == 0) {
    siteCovs <- NULL
  }

  ## only check non-sitecovs
  obsvars.df2 <- as.data.frame(obsvars.df[,
    !(names(obsvars.df) %in%
      names(siteCovs))
  ])
  names(obsvars.df2) <- names(obsvars.df)[
    !(names(obsvars.df) %in%
      names(siteCovs))
  ]

  yearlySiteCovs <- sapply(obsvars.df2, function(x) {
    obsmat <- matrix(x, M * nY, obsNum / nY, byrow = TRUE)
    l.u <- apply(obsmat, 1, function(y) {
      row.u <- unique(y)
      length(row.u[!is.na(row.u)])
    })
    ## if there are 0 or 1 unique vals per row, we have a sitecov
    if (all(l.u %in% 0:1)) {
      u <- apply(obsmat, 1, function(y) {
        row.u <- unique(y)
        ## only remove NAs if there are some non-NAs.
        if (!all(is.na(row.u))) {
          row.u <- row.u[!is.na(row.u)]
        }
        row.u
      })
      u
    }
  })
  yearlySiteCovs <- as.data.frame(yearlySiteCovs[
    !sapply(yearlySiteCovs, is.null)
  ])
  if (nrow(yearlySiteCovs) == 0) {
    yearlySiteCovs <- NULL
  }

  # Extract siteCovs and yearlySiteCovs from obsvars
  finalobsvars.df <- as.data.frame(obsvars.df[,
    !(names(obsvars.df) %in%
      c(names(siteCovs), names(yearlySiteCovs)))
  ])
  names(finalobsvars.df) <- names(obsvars.df)[
    !(names(obsvars.df) %in%
      c(names(siteCovs), names(yearlySiteCovs)))
  ]

  umf <- unmarkedMultFrame(
    y = y,
    siteCovs = siteCovs,
    obsCovs = finalobsvars.df,
    yearlySiteCovs = yearlySiteCovs,
    numPrimary = nY
  )
  return(umf)
}

reassignInPackage("formatMult", "unmarked", formatMult)
