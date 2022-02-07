# mh-execute

## Note on ithimr

This package is hosted on Github. Please install it using this command:

> remotes::install_github("ithim/ithim-r")

# GENERAL MH CONVENTIONS

save r: saveRds(rdsfile, version = 2)

city region names: <https://github.com/metahit/mh-execute/blob/master/inputs/mh_regions_lad_lookup.csv>

mh-execute: everything it calls, nothing it doesn't call

# Structure

Loop through all scenarios (`cyc_scen` and `car_scen`) in `global_scen`:

-   Set `NSAMPLES` (sample size) - currently set to 4

-   Set parameters with their distributions (for instance `MMET_CYCLING`)

-   Set default mode speeds

-   Get `DR` relationships from `ithimr` package (NOTE: Will be replaced by `drpa` package)

-   Read `DISEASE_INVENTORY` from `inputs/dose_response/disease_outcomes_lookup.csv`

-   Read `BACKGROUND_POLLUION_TABLE` from `inputs/background-air-pollution/1_apmeans.csv`

-   Read `age categories` from `inputs/scenarios/190330_sp_ind_codebook.xlsx`

-   Read scenario and city specific `emissions_distances` and `pa_distances` from `inputs/distances/`

-   Read `city_total_distances` by mode, from `inputs/distances/mode_road_city.csv`

-   Read `injury_table` from `inputs/injury/processed_injuries_9.Rds`

-   Read `city_regions_table` from `inputs/mh_regions_lad_lookup.csv` and filter for only those regions mentioned in the `injury_table`

-   Read `scenarios` files per for the `global_scen` for each `LAD` (Local Authority District)

-   Set `ithim_setup_parameters`

-   Loop through all `city_regions`

    -   Read `GBD` data

    -   Read `inh_distances`

    -   Filter `LAD` for the current `city_region`

    -   Create `per person summary` (`pp_summary`) for all scenarios by summarizing total mode specific duration per person. Also summarize, where present, `work_ltpa_marg_mmet`

    -   Summarise `DIST` for motorised modes (`car`, `motorcycle` and `bus`) for each `LAD` from `emission_distances`

    -   Set city-specific parameters for background `PM2.5` and `transport_share`, which are read from `BACKGROUND_POLLUION_TABLE`

    -   For each sample, do the following:

        -   Create `SiN` coefficients for `CAS_EXPONENT` and `STR_EXPONENT` from `parameters`
        -   Initialize `VEHICLE_INVENTORY` for the modes
        -   Using `DIST` and `pp_summary` to calculate `scenario_pm_calculations`
        -   Using `pm_conc_pp` calcualate `RR_AP_calculations` using `gen_ap_rr`
        -   Calculate `total_mmet` using `pp_summary`
        -   Calculate `RR_PA_calculations` using `ithimr::gen_pa_rr`
        -   Calculate combined `AP` and `PA` relative risk `RR_PA_AP_calculations` - by multiplying, using `combined_rr_ap_pa`
        -   

# 
