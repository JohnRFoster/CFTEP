#--FNC - Get NASS Data
get.nass.data <- function(source.vec, commodity.vec, year.vec, state.vec) {
  #--Loop Over Survey and Census
  for (s in 1:length(source.vec)) {
    #--Loop Over States
    for (j in 1:length(state.vec)) {
      #--Loop Over Commodities
      for (i in 1:length(commodity.vec)) {
        #--Querry to use
        params <- list(
          "source_desc" = source.vec[s],
          "commodity_desc" = commodity.vec[i],
          "agg_level_desc" = "COUNTY",
          "state_alpha" = state.vec[j]
        )

        if (rnassqs::nassqs_record_count(params) > 0) {
          #--Loop Over Years
          if (source.vec[s] == "CENSUS") {
            time.vec <- c(1997, 2002, 2007, 2012, 2017)
          }
          if (source.vec[s] == "SURVEY") {
            time.vec <- year.vec
          }

          for (y in 1:length(time.vec)) {
            params <- list(
              "source_desc" = source.vec[s],
              "commodity_desc" = commodity.vec[i],
              "agg_level_desc" = "COUNTY",
              "state_alpha" = state.vec[j],
              "year" = time.vec[y]
            )

            records <- rnassqs::nassqs_record_count(params)
            assertthat::assert_that(as.integer(records$count) <= 50000)

            if (records > 0) {
              raw <- rnassqs::nassqs_GET(params)
              parsed <- rnassqs::nassqs_parse(raw, as = 'data.frame')
            } #END Logical

            if (exists("parsed") == TRUE & exists("tmp.dat") == FALSE) {
              tmp.dat <- parsed
            }
            if (exists("parsed") == TRUE & exists("tmp.dat") == TRUE) {
              tmp.dat <- rbind.data.frame(tmp.dat, parsed)
            }
            if (exists("parsed") == TRUE) {
              rm(parsed)
            }
          } #END Year (Y) LOOP

          if (exists("tmp.dat") == TRUE & exists("year.dat") == FALSE) {
            year.dat <- tmp.dat
          }
          if (exists("tmp.dat") == TRUE & exists("year.dat") == TRUE) {
            year.dat <- rbind.data.frame(year.dat, tmp.dat)
          }
          if (exists("tmp.dat") == TRUE) {
            rm(tmp.dat)
          }
        } #END Logical
      } #END Commodity (I) Loop

      if (exists("year.dat") == TRUE & exists("com.dat") == FALSE) {
        com.dat <- year.dat
      }
      if (exists("year.dat") == TRUE & exists("com.dat") == TRUE) {
        com.dat <- rbind.data.frame(com.dat, year.dat)
      }
      if (exists("year.dat") == TRUE) {
        rm(year.dat)
      }

      print(paste0(
        state.vec[j],
        "  Percent complete: ",
        round(j / length(state.vec) * 100, digits = 1)
      ))
    } #END State (J) Loop

    if (exists("com.dat") == TRUE & exists("source.dat") == FALSE) {
      source.dat <- com.dat
    }
    if (exists("com.dat") == TRUE & exists("source.dat") == TRUE) {
      source.dat <- rbind.data.frame(source.dat, com.dat)
    }
    if (exists("com.dat") == TRUE) {
      rm(com.dat)
    }
  } #END Source (S) Loop

  #--Remove commas from value (ignore character values)
  source.dat$Value.numeric <- as.numeric(gsub(",", "", source.dat$Value))

  #--Clean up Values
  source.dat$Value <- trimws(source.dat$Value, which = c("both"))

  #--Add FIPS Code
  source.dat$FIPS <- paste0(source.dat$state_fips_code, source.dat$county_code)

  rm(com.dat, year.dat, tmp.dat, parsed)

  return(source.dat)
} #END Function


#--FNC - Impute missing values

impute.missing <- function(dat) {
  require(imputeTS)

  #--Convert FIPS codes
  dat$FIPS.numeric <- as.integer(dat$FIPS)

  #----Find and fill in missing values

  #--Find withheld data
  num.start <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  missing.loc <- dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ]
  missing.loc <- missing.loc[
    order(missing.loc$FIPS.numeric, missing.loc$year),
  ]
  missing.loc <- plyr::count(missing.loc[, c(
    "FIPS.numeric",
    "source_desc",
    "domaincat_desc",
    "commodity_desc",
    "statisticcat_desc",
    "short_desc",
    "unit_desc"
  )])

  #--Loop over FIPS
  for (i in 1:nrow(missing.loc)) {
    val.vec <- dat[
      dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
        dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
        dat$short_desc == missing.loc[i, "short_desc"] &
        dat$unit_desc == missing.loc[i, "unit_desc"],
    ]

    x <- val.vec[, "Value.numeric"]

    if (length(x[is.na(x) == FALSE]) >= 2) {
      x <- na_interpolation(ts(val.vec[, "Value.numeric"]))

      dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"],
        "Value.numeric"
      ] <- x
    } #END Logical
  }

  num.left <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  print(paste0("Imputed: ", num.start - num.left))
  print(paste0("Missing: ", num.left))

  return(dat)
} #END Function


#--FNC - Impute missing values

impute.missing.inventory <- function(
  dat,
  statisticcat_desc.val,
  short_desc.val
) {
  require(imputeTS)

  #--Convert FIPS codes
  dat$FIPS.numeric <- as.integer(dat$FIPS)

  #----Find and fill in missing values

  #--Find withheld data
  num.start <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  missing.loc <- dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ]
  missing.loc <- missing.loc[
    order(missing.loc$FIPS.numeric, missing.loc$year),
  ]
  missing.loc <- plyr::count(missing.loc[, c(
    "FIPS.numeric",
    "year",
    "domaincat_desc",
    "commodity_desc",
    "statisticcat_desc",
    "short_desc",
    "unit_desc"
  )])

  #--Loop over FIPS
  for (i in 1:nrow(missing.loc)) {
    val.vec <- dat[
      dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
        dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
        dat$short_desc == missing.loc[i, "short_desc"] &
        dat$unit_desc == missing.loc[i, "unit_desc"],
    ]

    x <- val.vec[, "Value.numeric"]
    x <- x[is.na(x) == FALSE]

    if (length(x) != 0) {
      x <- mean(x, na.rm = TRUE)

      n <- length(dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ])

      y <- rpois(n = n, lambda = x)

      dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ] <- y
    } #END Logical
  } #END Loop

  num.left <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  print(paste0("Imputed: ", num.start - num.left))
  print(paste0("Missing: ", num.left))

  return(dat)
} #END Function


#--FNC - Impute missing values

spatial.impute.missing <- function(
  dat,
  adj.file,
  statisticcat_desc.val,
  short_desc.val
) {
  require(imputeTS)

  #--Convert FIPS codes
  dat$FIPS.numeric <- as.integer(dat$FIPS)

  #----Find and fill in missing values

  #--Find withheld data
  num.start <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  missing.loc <- dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ]
  missing.loc <- missing.loc[
    order(missing.loc$FIPS.numeric, missing.loc$year),
  ]
  missing.loc <- plyr::count(missing.loc[, c(
    "FIPS.numeric",
    "year",
    "domaincat_desc",
    "commodity_desc",
    "statisticcat_desc",
    "short_desc",
    "unit_desc"
  )])

  #--Loop over FIPS
  for (i in 1:nrow(missing.loc)) {
    #Get Adjacent fips codes
    adj.vec <- adj.file[
      adj.file$FIPS == missing.loc[i, "FIPS.numeric"],
      "FIPS.neighbor"
    ]

    val.vec <- dat[
      dat$FIPS.numeric %in%
        adj.vec &
        dat$year == missing.loc[i, "year"] &
        dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
        dat$short_desc == missing.loc[i, "short_desc"] &
        dat$unit_desc == missing.loc[i, "unit_desc"],
    ]

    x <- val.vec[, "Value.numeric"]
    x <- x[is.na(x) == FALSE]

    if (length(x) != 0) {
      x <- mean(x, na.rm = TRUE)

      n <- length(dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$year == missing.loc[i, "year"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ])

      y <- rpois(n = n, lambda = x)

      dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$year == missing.loc[i, "year"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ] <- y
    } #END Logical
  } #END Loop

  num.left <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  print(paste0("Imputed: ", num.start - num.left))
  print(paste0("Missing: ", num.left))

  return(dat)
} #END Function


#--FNC - Impute missing values

impute.missing.using.state.mean <- function(
  dat,
  statisticcat_desc.val,
  short_desc.val
) {
  require(imputeTS)

  #--Convert FIPS codes
  dat$FIPS.numeric <- as.integer(dat$FIPS)

  #--Make copy of original
  org.dat <- dat

  #----Find and fill in missing values

  #--Find withheld data
  num.start <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc",
      "state_name"
    )
  ])

  missing.loc <- dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc",
      "state_name"
    )
  ]
  missing.loc <- missing.loc[
    order(missing.loc$FIPS.numeric, missing.loc$year),
  ]
  missing.loc <- plyr::count(missing.loc[, c(
    "FIPS.numeric",
    "year",
    "domaincat_desc",
    "commodity_desc",
    "statisticcat_desc",
    "short_desc",
    "unit_desc",
    "state_name"
  )])

  #--Loop over FIPS
  for (i in 1:nrow(missing.loc)) {
    #--Generate State Values

    #Find FIPS with operators fewer than 3
    loc.tab <- org.dat[
      org.dat$state_name == missing.loc[i, "state_name"] &
        org.dat$year == missing.loc[i, "year"] &
        org.dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        org.dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        org.dat$unit_desc == "OPERATIONS",
    ]

    loc.tab <- loc.tab[loc.tab$Value.numeric < 3, "FIPS.numeric"]

    #Generate vector of values
    val.vec <- dat[
      dat$FIPS.numeric %in%
        loc.tab &
        dat$year == missing.loc[i, "year"] &
        dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
        dat$short_desc == missing.loc[i, "short_desc"] &
        dat$unit_desc == missing.loc[i, "unit_desc"],
    ]

    x <- val.vec[, "Value.numeric"]
    x <- x[is.na(x) == FALSE]

    if (length(x) != 0) {
      x <- mean(x, na.rm = TRUE)

      n <- length(dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$year == missing.loc[i, "year"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ])

      y <- rpois(n = n, lambda = x)

      dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$year == missing.loc[i, "year"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ] <- y
    } #END Logical
  } #END Loop

  num.left <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "source_desc",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc",
      "state_name"
    )
  ])

  print(paste0("Imputed: ", num.start - num.left))
  print(paste0("Missing: ", num.left))

  return(dat)
} #END Function


#--FNC - Impute missing values

impute.missing.using.national.mean <- function(
  dat,
  statisticcat_desc.val,
  short_desc.val
) {
  require(imputeTS)

  #--Convert FIPS codes
  dat$FIPS.numeric <- as.integer(dat$FIPS)

  #--Make copy of original
  org.dat <- dat

  #----Find and fill in missing values

  #--Find withheld data
  num.start <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc",
      "state_name"
    )
  ])

  missing.loc <- dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc",
      "state_name"
    )
  ]
  missing.loc <- missing.loc[
    order(missing.loc$FIPS.numeric, missing.loc$year),
  ]
  missing.loc <- plyr::count(missing.loc[, c(
    "FIPS.numeric",
    "year",
    "domaincat_desc",
    "commodity_desc",
    "statisticcat_desc",
    "short_desc",
    "unit_desc",
    "state_name"
  )])

  #--Loop over FIPS
  for (i in 1:nrow(missing.loc)) {
    #--Generate State Values

    #Find FIPS with operators fewer than 3
    loc.tab <- org.dat[
      org.dat$year == missing.loc[i, "year"] &
        org.dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        org.dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        org.dat$unit_desc == "OPERATIONS",
    ]

    loc.tab <- loc.tab[loc.tab$Value.numeric < 3, "FIPS.numeric"]

    #Generate vector of values
    val.vec <- dat[
      dat$FIPS.numeric %in%
        loc.tab &
        dat$year == missing.loc[i, "year"] &
        dat$commodity_desc == missing.loc[i, "commodity_desc"] &
        dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
        dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
        dat$short_desc == missing.loc[i, "short_desc"] &
        dat$unit_desc == missing.loc[i, "unit_desc"],
    ]

    x <- val.vec[, "Value.numeric"]
    x <- x[is.na(x) == FALSE]

    if (length(x) != 0) {
      x <- mean(x, na.rm = TRUE)

      n <- length(dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$year == missing.loc[i, "year"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ])

      y <- rpois(n = n, lambda = x)

      dat[
        dat$FIPS.numeric == missing.loc[i, "FIPS.numeric"] &
          dat$year == missing.loc[i, "year"] &
          dat$commodity_desc == missing.loc[i, "commodity_desc"] &
          dat$domaincat_desc == missing.loc[i, "domaincat_desc"] &
          dat$statisticcat_desc == missing.loc[i, "statisticcat_desc"] &
          dat$short_desc == missing.loc[i, "short_desc"] &
          dat$unit_desc == missing.loc[i, "unit_desc"] &
          is.na(dat$Value.numeric) == TRUE,
        "Value.numeric"
      ] <- y
    } #END Logical
  } #END Loop

  num.left <- nrow(dat[
    dat$Value %in%
      c("(D)", "(S)") &
      is.na(dat$Value.numeric) == TRUE &
      dat$statisticcat_desc == statisticcat_desc.val &
      dat$short_desc %in% short_desc.val &
      dat$sector_desc == "ANIMALS & PRODUCTS",
    c(
      "FIPS.numeric",
      "year",
      "domaincat_desc",
      "commodity_desc",
      "statisticcat_desc",
      "short_desc",
      "unit_desc"
    )
  ])

  print(paste0("Imputed: ", num.start - num.left))
  print(paste0("Missing: ", num.left))

  return(dat)
} #END Function
