newformatLong <- function(dfin, species = NULL, timename = "JulianDate", type) {
	if (missing(type)) {
		stop("type must be supplied")
	}
	nc <- ncol(dfin)
	dfin[[nc + 1]] <- dfin[[2]]
	names(dfin)[nc + 1] <- "Date"
	timefac <- is.factor(dfin[[2]])
	if (!is.null(species)) {
		dfin$y <- ifelse(dfin$species == species, dfin$y, 0)
		dfin$y[is.na(dfin$y)] <- 0
		dfin$species = NULL
	}
	names(dfin)[3] <- "y"
	dfin <- unmarked:::dateToObs(dfin)
	dfnm <- colnames(dfin)
	nV <- length(dfnm) - 1
	expr <- substitute(
		recast(
			dfin,
			newvar ~ obsNum + variable,
			id.var = c(dfnm[1], "obsNum"),
			measure.var = dfnm[3]
		),
		list(newvar = as.name(dfnm[1]))
	)
	y <- as.matrix(eval(expr)[, -1])
	attr(y, "class") <- "matrix"

	## At this point the data frame is recast into an array of matrices of dimension R (# sites) x
	## J (max # of sampling periods per site)
	## The problem is that recast (specifically the melt) doesn't play well with a combination of
	## factor (categorical variables) and integer/numeric variables
	## Here, I perform the operations separately for fac/char, if present, and int/num and recombine

	## Factors/characters first; note that melt converts factors to character by default
	facs <- sapply(dfin[4:nV], function(x) any(is.character(x), is.factor(x)))
	if (sum(facs) > 0) {
		expr <- substitute(
			recast(
				dfin,
				newvar ~ obsNum ~ variable,
				id.var = c(dfnm[1], "obsNum"),
				measure.var = dfnm[4:nV][facs]
			),
			list(newvar = as.name(dfnm[1]))
		)
		facvars <- eval(expr)
		if (timefac) {
			which.date <- which(dimnames(facvars)$variable == "Date")
		}
		dimnames(facvars)$variable[which.date] <- timename
	}

	## Non-factors
	if (length(facs == 0) > 0) {
		expr <- substitute(
			recast(
				dfin,
				newvar ~ obsNum ~ variable,
				id.var = c(dfnm[1], "obsNum"),
				measure.var = dfnm[4:nV][!facs]
			),
			list(newvar = as.name(dfnm[1]))
		)
		nonfacvars <- eval(expr)
		if (!timefac) {
			which.date <- which(dimnames(nonfacvars)$variable == "Date")
			dimnames(nonfacvars)$variable[which.date] <- timename
		}

		if (sum(facs) > 0) {
			obsvars.matlist <- c(
				unmarked:::arrToList(facvars),
				unmarked:::arrToList(nonfacvars)
			)
		} else {
			obsvars.matlist <- unmarked:::arrToList(nonfacvars)
		}
	} else {
		obsvars.matlist <- unmarked:::arrToList(facvars)
	}

	## obsvars.matlist is a list of length nV - 3 matrices of dimensions M (nSites) x sample periods
	## For each matrix, check if ALL rows have unique length of 1

	which.siteCovs <- which(
		lapply(obsvars.matlist, function(x) {
			vec <- apply(x, 1, function(i) length(unique(i)))
			same <- all(vec == 1)
			same
		}) ==
			1
	)
	dfin.siteCovs <- which(dfnm %in% names(which.siteCovs))
	sitevars.df <- unique(dfin[, c(1, dfin.siteCovs)])[, -1, drop = FALSE]

	## Now proceed with observation level variables
	obsvars.veclist <- lapply(obsvars.matlist[-which.siteCovs], function(x) {
		as.vector(t(x))
	})
	obsvars.df <- data.frame(obsvars.veclist)

	do.call(type, list(y = y, siteCovs = sitevars.df, obsCovs = obsvars.df))
}
