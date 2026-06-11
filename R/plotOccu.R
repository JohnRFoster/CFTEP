plotOccu <- function(
    occuFit,
    data,
    plotVars = "",
    varNames = "",
    modelVars = "",
    covar_ct = covar_meds,
    covar_var = covar_mads,
    ncols = 2,
    output = TRUE,
    fn = "./occupancy_plot.png",
    panel_width = 4
) {
    #plotVars: covariates to plot, in desired plotting order
    #varNames: covariate labels for the plots, in same order as plotVars
    #modelVars: other covariates in the model but not to be plotted

    # Check for plotVars
    if (length(plotVars) == 1 && plotVars[1] == "") {
        stop("You must enter at least one covariate to plot.")
    }

    # Check for names
    if (length(varNames) == 1 && varNames == "") {
        varNames <- plotVars
    }

    # List to hold ggplot objects
    figure <- vector(mode = "list", length = length(plotVars))
    names(figure) <- plotVars

    # String of all variables used in the model, not just those plotted
    modelVars <- c(plotVars, modelVars)

    for (var in plotVars) {
        # Specify new data for prediction and intervals
        # other covariates are held at their median (i.e., scaled value of 0)
        uniqueVals <- sort(unique(data[, var]))
        newDat <- data.frame(
            var_scaled = uniqueVals,
            combinations(length(modelVars) - 1)
        )
        names(newDat) <- c(var, modelVars[-which(modelVars == var)])

        # Predicting over the range of current variable, all others at medians
        newDat <- predict(
            occuFit,
            newdata = newDat,
            type = "state",
            appendData = TRUE
        )

        # Return covariates to original scale for interpretability
        for (i in seq(ncol(newDat))) {
            cur_var <- names(newDat)[i]
            if (cur_var %in% names(covar_ct)) {
                newDat[, cur_var] <-
                    (newDat[, cur_var] * covar_var[cur_var] + covar_ct[cur_var])
            }
        }

        rugDat <- data.frame(
            var = data[, var] * covar_var[var] + covar_ct[var],
            Predicted = runif(length(data[, var]))
        ) # nonsense to make ggplot happy

        plotDat <- data.frame(var = newDat[, var], newDat)

        # Create plot

        # Get some parameters for labels
        xmin <- min(plotDat$var)
        xmax <- max(plotDat$var)
        x_range <- xmax - xmin

        p <- ggplot(plotDat, aes(x = var, y = Predicted)) +
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.3) +
            geom_line(size = 1.25)

        # Tidy the plot
        p <- p +
            xlab(varNames[which(plotVars == var)]) +
            geom_rug(
                data = rugDat,
                size = 0.1,
                sides = "b",
                position = "jitter",
                alpha = 0.3
            ) +
            theme(
                legend.position = "none",
                plot.margin = unit(c(0.1, 0.2, 0.1, 0), "cm")
            )
        #  theme(legend.justification=c(0.5,1), legend.position=c(0.5, 1),
        #        legend.key = element_blank(), legend.direction = "horizontal",
        #        legend.key.width = unit(0.05, units = "npc"))

        # Label panels
        index <- which(plotVars == var)
        #p <- p + annotate("text", x = min(plotDat$var), y = 1,
        #                  label = LETTERS[index], hjust = 0, vjust=0.6,
        #                  size = 8)

        # Axis manipulation for multipanel plot
        p <- p +
            if (index %in% seq(1, length(plotVars), ncols)) {
                scale_y_continuous(
                    "Predicted occupancy",
                    limits = c(0, 1),
                    breaks = seq(0, 1, 0.2)
                )
            } else {
                scale_y_continuous(
                    "",
                    limits = c(0, 1),
                    breaks = seq(0, 1, 0.2)
                )
            }

        figure[[var]] <- p
    }

    ncols <- ifelse(length(plotVars) == 1, 1, ncols)
    fig <- multiplot(
        plotlist = figure,
        layout = matrix(1:length(plotVars), ncol = ncols, byrow = T)
    )

    if (output) {
        if (ncols < 2) {
            png(
                file = fn,
                width = panel_width,
                height = panel_width * ceiling(length(plotVars) / ncols),
                units = "in",
                res = 600
            )
        } else {
            png(
                file = fn,
                width = panel_width * ncols,
                height = panel_width * ceiling(length(plotVars) / ncols),
                units = "in",
                res = 600
            )
        }

        fig
        dev.off()
    } else {
        fig
    }
}
