newshow <- function(object) {
    out <- as(object, "data.frame")
    rownames(out) <- out$model
    ICs <- grep("AIC", names(out), value = TRUE)
    out <- out[, c("nPars", ICs[1], "delta", ICs[2], "cumltvWt", "Rsq")]
    if (all(is.na(out$Rsq))) {
        out$Rsq <- NULL
    }
    print(format(out, digits = 2, nsmall = 2))
}

environment(newshow) <- environment(occu)

setMethod("show", signature(object = "unmarkedModSel"), newshow)

newmodSel <- function(object, IC = c("AICc", "AIC"), ...) {
    IC <- match.arg(IC)
    .local <- function(object, nullmod = NULL) {
        if (!is.character(nullmod) && !is.null(nullmod)) {
            stop(
                "nullmod must be character name of null model fit in the fitlist."
            )
        }
        fits <- object@fits
        estList <- lapply(fits, coef, altNames = TRUE)
        seList <- lapply(fits, function(x) {
            se <- tryCatch(
                sqrt(diag(vcov(x, altNames = TRUE))),
                error = function(e) simpleError("Hessian is singular.")
            )
            if (identical(class(se)[1], "simpleError")) {
                cat(se$message, fill = TRUE)
                se <- rep(NA, length(coef(x)))
            }
            return(se)
        })
        eNames <- sort(unique(unlist(sapply(estList, names))))
        seNames <- paste("SE", eNames, sep = "")
        eseNames <- character(l <- length(c(eNames, seNames)))
        eseNames[seq(1, l, by = 2)] <- eNames
        eseNames[seq(2, l, by = 2)] <- seNames
        cNames <- c("model", "formula", eseNames)
        out <- data.frame(matrix(
            NA,
            ncol = length(cNames),
            nrow = length(fits)
        ))
        colnames(out) <- cNames
        out$model <- names(fits)
        out$formula <- sapply(fits, function(x) {
            f <- as.character(x@formula)
            f <- paste(f[2], "~", f[3])
            f
        })
        for (i in 1:length(eNames)) {
            out[, eNames[i]] <- sapply(estList, function(x) x[eNames[i]])
            out[, seNames[i]] <- sapply(seList, function(x) x[eNames[i]])
        }
        out$Converge <- sapply(fits, function(x) x@opt$convergence)
        out$CondNum <- sapply(fits, function(x) unmarked:::cn(x))
        out$negLogLike <- sapply(fits, function(x) x@negLogLike)
        out$nPars <- sapply(fits, function(x) length(coef(x)))
        out$n <- sapply(fits, function(x) sampleSize(x))
        if (IC == "AICc") {
            #out$AIC <- sapply(fits, function(x) x@AIC)
            out$AICc <- sapply(fits, function(x) x@AIC) + #AIC
                (2 * out$nPars * (out$nPars - 1)) / (out$n - out$nPars - 1)
            out$delta <- out$AICc - min(out$AICc)
            out$AICcwt <- exp(-out$delta / 2)
            out$AICcwt <- format(
                round(out$AICcwt / sum(out$AICcwt), 3),
                nsmall = 3
            )
        } else {
            out$AIC <- sapply(fits, function(x) x@AIC)
            out$delta <- out$AIC - min(out$AIC)
            out$AICwt <- exp(-out$delta / 2)
            out$AICwt <- out$AICwt / sum(out$AICwt)
        }
        out$Rsq <- NA
        if (!is.null(nullmod)) {
            if (is.na(match(nullmod, names(fits)))) {
                stop(paste("No fit named", nullmod, "was found in fits."))
            }
            nullmod <- fits[[nullmod]]
            out$Rsq <- sapply(fits, nagR2, nullmod)
        }
        if (IC == "AICc") {
            out <- out[order(out$AICc), ]
            out$cumltvWt <- cumsum(out$AICcwt)
        } else {
            out <- out[order(out$AIC), ]
            out$cumltvWt <- cumsum(out$AICwt)
        }
        msout <- new(
            "unmarkedModSel",
            Full = out,
            Names = rbind(Coefs = eNames, SEs = seNames)
        )
        return(msout)
    }
    .local(object, ...)
}
