compute_fossil_threshold <- function(case) {
    ## emissions factors taken from Mohr et al (2015), Table 4
    emissions_conversion <- c(
        'coal_anthracite' = 2122,
        'coal_bituminous' = 2026,
        'coal_sub_bituminous' = 1510,
        'coal_lignite' = 1126,
        'coal_semianthracite' = 2107,
        'coal_bit_and_subbit' = 1768,
        'coal_black' = 2107,
        'coal_brown' = 1318,
        'oil_conventional' = 434.2,
        'oil_tight' = 434.2,
        'oil_bitumen' = 434.2,
        'oil_extra_heavy' = 434.2,
        'oil_kerogen' = 610,
        'gas_conventional' = 54.6,
        'gas_CBM' = 54.6,
        'gas_shale' = 54.6,
        'gas_tight' = 54.6,
        'gas_hydrates' = 54.6
    )

    ## URR mass scenario estimates taken from Mohr et al (2015), Table 2
    
    low_URR <- c(
        'coal_anthracite' = 14,
        'coal_bituminous' = 454,
        'coal_sub_bituminous' = 34,
        'coal_lignite' = 89,
        'coal_semianthracite' = 1,
        'coal_bit_and_subbit' = 1,
        'coal_black' = 61,
        'coal_brown' = 9,
        'oil_conventional' = 2479,
        'oil_tight' = 340,
        'oil_bitumen' = 423,
        'oil_extra_heavy' = 250,
        'oil_kerogen' = 2,
        'gas_conventional' = 10595,
        'gas_CBM' = 830,
        'gas_shale' = 1286,
        'gas_tight' = 551,
        'gas_hydrates' = 0
    )

    base_URR <- c(
        'coal_anthracite' = 17,
        'coal_bituminous' = 616,
        'coal_sub_bituminous' = 138,
        'coal_lignite' = 306,
        'coal_semianthracite' = 2,
        'coal_bit_and_subbit' = 1,
        'coal_black' = 61,
        'coal_brown' = 20,
        'oil_conventional' = 2547,
        'oil_tight' = 369,
        'oil_bitumen' = 457,
        'oil_extra_heavy' = 302,
        'oil_kerogen' = 769,
        'gas_conventional' = 12512,
        'gas_CBM' = 1045,
        'gas_shale' = 6679,
        'gas_tight' = 1867,
        'gas_hydrates' = 4383
    )

     hi_URR <- c(
        'coal_anthracite' = 14,
        'coal_bituminous' = 814,
        'coal_sub_bituminous' = 159,
        'coal_lignite' = 588,
        'coal_semianthracite' = 3,
        'coal_bit_and_subbit' = 2,
        'coal_black' = 117,
        'coal_brown' = 22,
        'oil_conventional' = 3686,
        'oil_tight' = 644,
        'oil_bitumen' = 664,
        'oil_extra_heavy' = 606,
        'oil_kerogen' = 1937,
        'gas_conventional' = 21466,
        'gas_CBM' = 1852,
        'gas_shale' = 6723,
        'gas_tight' = 2406,
        'gas_hydrates' = 12036
     )
    ## select URR estimate based on case
    URR <- switch(case, 'low' = low_URR, 'base' = base_URR, 'high' = hi_URR)
    fuel_dat <- data.frame(URR=URR, conv_factor=emissions_conversion / 1000) # divide by 1000 to convert from Mt CO2e to Gt CO2e
    ## compute thresholds by technology
    hi_row_idx <- grepl('coal', rownames(fuel_dat))
    hi_thresh <- sum(with(fuel_dat[hi_row_idx, ], URR * conv_factor))
    lo_thresh <- sum(with(fuel_dat[!hi_row_idx, ], URR * conv_factor))
    c(hi_thresh, lo_thresh)
}
