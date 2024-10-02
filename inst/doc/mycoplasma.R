## ----include = FALSE--------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 6,
  fig.path = ""
)


## ----setup------------------------------------------------
library(poems)
library(epizootic)
nsims <- 1
parallel_cores <- 1
timesteps <- 23


## ----region-----------------------------------------------
library(raster)
epizootic::finch_region
region <- Region$new(template_raster = finch_region)
raster::plot(region$region_raster, colNA = "blue",
             main = "Study Region")


## ----bsl--------------------------------------------------
data("bsl_raster")
raster::plot(bsl_raster, colNA = "blue",
             main = "Breeding Season Length (days)")


## ----carrying-capacity------------------------------------
data("habitat_suitability")
raster::plot(habitat_suitability, colNA = "blue",
             main = "Habitat Suitability")


## ----initial-abundance------------------------------------
data("initial_abundance")
region$raster_from_values(initial_abundance[1, ]) |>
  raster::plot(colNA = "blue", main = "Susceptible Juveniles")


## ----fixed parameters-------------------------------------
model_template <- DiseaseModel$new(
  simulation_function = "disease_simulator",
  region = region,
  time_steps = timesteps,
  populations = region$region_cells,
  replicates = 1,
  stages = 2, # life cycle stages
  compartments = 4, # disease compartments
  seasons = 2, # seasons in a year
  mortality_unit = list(c(1, 1, 0, 0, 1, 1, 0, 0), # is mortality seasonal
                        c(1, 1, 0, 0, 1, 1, 0, 0)), # or daily?
  fecundity_unit = rep(1, 8), # is fecundity seasonal or daily?
  fecundity_mask = rep(c(0, 1), 4), # which stages/compartments reproduce
  transmission_unit = rep(0, 4), # Is transmission rate seasonal or daily?
  transmission_mask = c(1, 1, 0, 0, 1, 1, 0, 0), # which compartments can become                                                  # infected
  recovery_unit = rep(0, 4), # is recovery rate seasonal or daily?
  recovery_mask = c(0, 0, 1, 1, 0, 0, 1, 1), # which compartments can recover
  breeding_season_length = bsl_raster,
  hs = habitat_suitability,
  abundance = initial_abundance,
  dispersal_type = "stages", # indicates that life cycle stages disperse
                             # differently
  simulation_order = list(c("transition", "season_functions", "results"),
                          c("dispersal", "season_functions", "results")),
  results_selection = c("abundance"), # what do want included in the result?
  results_breakdown = "segments", # are the results broken down by life cycle
                                  # stage, disease compartment, or both?
  season_functions = list(siri_model_summer, siri_model_winter),
  verbose = FALSE,
  random_seed = 648,
  attribute_aliases = list(
    mortality_Sj_summer = "mortality$summer$a",
    mortality_Sa_summer = "mortality$summer$b",
    mortality_I1j_summer = "mortality$summer$c",
    mortality_I1a_summer = "mortality$summer$d",
    mortality_Rj_summer = "mortality$summer$e",
    mortality_Ra_summer = "mortality$summer$f",
    mortality_I2j_summer = "mortality$summer$g",
    mortality_I2a_summer = "mortality$summer$h",
    mortality_Sj_winter = "mortality$winter$a",
    mortality_Sa_winter = "mortality$winter$b",
    mortality_I1j_winter = "mortality$winter$c",
    mortality_I1a_winter = "mortality$winter$d",
    mortality_Rj_winter = "mortality$winter$e",
    mortality_Ra_winter = "mortality$winter$f",
    mortality_I2j_winter = "mortality$winter$g",
    mortality_I2a_winter = "mortality$winter$h",
    beta_Sj_summer = "transmission$summer$a",
    beta_Sa_summer = "transmission$summer$b",
    beta_Rj_summer = "transmission$summer$c",
    beta_Ra_summer = "transmission$summer$d",
    beta_Sj_winter = "transmission$winter$a",
    beta_Sa_winter = "transmission$winter$b",
    beta_Rj_winter = "transmission$winter$c",
    beta_Ra_winter = "transmission$winter$d",
    recovery_I1j_summer = "recovery$summer$a",
    recovery_I1a_summer = "recovery$summer$b",
    recovery_I2j_summer = "recovery$summer$c",
    recovery_I2a_summer = "recovery$summer$d",
    recovery_I1j_winter = "recovery$winter$a",
    recovery_I1a_winter = "recovery$winter$b",
    recovery_I2j_winter = "recovery$winter$c",
    recovery_I2a_winter = "recovery$winter$d",
    dispersal1 = "dispersal$a",
    dispersal2 = "dispersal$b"
  )
)


## ----test consistency and completeness--------------------
model_template$is_complete()
model_template$is_consistent()


## ----lhs--------------------------------------------------
lhs_generator <- LatinHypercubeSampler$new()

# Dispersal parameters
lhs_generator$set_beta_parameter("dispersal_p_juv", alpha = 9.834837,
                                 beta = 2.019125)
lhs_generator$set_beta_parameter("dispersal_p_adult", alpha = 1.5685,
                                 beta = 2.365266)
lhs_generator$set_truncnorm_parameter("dispersal_r_juv", lower = 0, upper = 1500,
                                      mean = 725.9071, sd = sqrt(204006.6))
lhs_generator$set_normal_parameter("dispersal_r_adult", mean = 679.4172,
                                   sd = sqrt(18594.59))
lhs_generator$set_uniform_parameter("dispersal_source_n_k_cutoff", lower = 0,
                                    upper = 1)
lhs_generator$set_uniform_parameter("dispersal_source_n_k_threshold", lower = 0,
                                    upper = 1)
lhs_generator$set_uniform_parameter("dispersal_target_n_k_cutoff", lower = 0,
                                    upper = 1)
lhs_generator$set_uniform_parameter("dispersal_target_n_k_threshold", lower = 0,
                                    upper = 1)

# Population growth parameters
lhs_generator$set_uniform_parameter("abundance_threshold", lower = 0, upper = 25, decimals = 0)
lhs_generator$set_uniform_parameter("initial_release", lower = 5, upper = 50, decimals = 0)
lhs_generator$set_uniform_parameter("density_max", lower = 186000, upper = 310000, decimals = 0)
lhs_generator$set_poisson_parameter("fecundity", lambda = 8.509018)

# Transmission parameters
lhs_generator$set_uniform_parameter("beta_Sa_winter", lower = 0, upper = 0.07588)
lhs_generator$set_uniform_parameter("beta_Sa_summer", lower = 0, upper = 0.007784)
lhs_generator$set_triangular_parameter("Sj_multiplier", lower = 0, upper = 8.5,
                                       mode = 3)
lhs_generator$set_beta_parameter("beta_I2_modifier", alpha = 1.547023,
                                 beta = 0.4239236)

# Mortality parameters
lhs_generator$set_beta_parameter("mortality_Sj_winter", alpha = 3.962104,
                                 beta = 2.228683)
lhs_generator$set_beta_parameter("mortality_Sa_winter", alpha = 21.89136,
                                 beta = 19.59278)
lhs_generator$set_beta_parameter("mortality_Sj_summer", alpha = 14.51403,
                                 beta = 21.53632)
lhs_generator$set_beta_parameter("mortality_I1j_summer", alpha = 2.756404,
                                 beta = 62.47181)
lhs_generator$set_beta_parameter("mortality_I1j_winter", alpha = 2.756404,
                                 beta = 62.47181)
lhs_generator$set_beta_parameter("mortality_I1a_summer", alpha = 1.771183,
                                 beta = 27.19457)
lhs_generator$set_beta_parameter("mortality_I1a_winter", alpha = 1.678424,
                                 beta = 41.15975)
lhs_generator$set_beta_parameter("mortality_I2_modifier", alpha = 1.033367,
                                 beta = 3.505319)

# Recovery parameters
lhs_generator$set_beta_parameter("recovery_I1", alpha = 9.347533,
                                 beta = 620.1732)
lhs_generator$set_beta_parameter("recovery_I2", alpha = 1.181112,
                                 beta = 29.18489)

# How many birds are infected in DC at timestep 1?
lhs_generator$set_uniform_parameter("infected_t1", lower = 1, upper = 20, decimals = 0)

sample_data <- lhs_generator$generate_samples(number = nsims,
                        random_seed = 630)
sample_data$sample <- 1:nsims
sample_data$mortality_Sa_summer <- 0
sample_data$mortality_I2j_summer <- sample_data$mortality_I2_modifier * sample_data$mortality_I1j_summer
sample_data$mortality_I2j_winter <- sample_data$mortality_I2_modifier * sample_data$mortality_I1j_winter
sample_data$mortality_I2a_winter <- sample_data$mortality_I2_modifier * sample_data$mortality_I1a_winter
sample_data$mortality_I2a_summer <- sample_data$mortality_I2_modifier * sample_data$mortality_I1a_summer
sample_data$mortality_Rj_summer <- sample_data$mortality_Sj_summer
sample_data$mortality_Ra_summer <- sample_data$mortality_Sa_summer
sample_data$mortality_Rj_winter <- sample_data$mortality_Sj_winter
sample_data$mortality_Ra_winter <- sample_data$mortality_Sa_winter
sample_data


## ----spatial correlation----------------------------------
env_corr <- SpatialCorrelation$new(region = region,
                                   amplitude = 0.99,
                                   breadth = 850,
                                   distance_scale = 1000)
env_corr$calculate_compact_decomposition(decimals = 4)


## ----adult dispersal--------------------------------------
b_lookup <- data.frame(d_max = -Inf, b = 0:904)
for (i in 2:904) {
  b_lookup$d_max[i] <- which.max(exp(-1*(1:1501)/b_lookup$b[i]) <= 0.19)
}
b_lookup$d_max[905] <- 1501

adult_dispersal_gen <- DispersalGenerator$new(
  region = region,
  spatial_correlation = env_corr,
  distance_classes = seq(10, 1500, 10),
  distance_scale = 1000, # km
  dispersal_function_data = b_lookup,
  inputs = c("dispersal_p_adult",
             "dispersal_r_adult"),
  attribute_aliases = list(dispersal_r_adult = "dispersal_max_distance",
                           dispersal_p_adult = "dispersal_proportion",
                           dispersal_adult = "dispersal_data"),
  decimals = 3
)
# This stage is computationally intensive
distance_matrix <- adult_dispersal_gen$calculate_distance_matrix()
adult_dispersal_gen$calculate_distance_data(distance_matrix = distance_matrix)


## ----juvenile dispersal-----------------------------------
juvenile_dispersal_gen <- DispersalGenerator$new(
  region = region,
  spatial_correlation = env_corr,
  distance_classes = seq(10, 1500, 10),
  distance_scale = 1000, # km
  dispersal_function_data = b_lookup,
  decimals = 3,
  inputs = c("dispersal_p_juv",
             "dispersal_r_juv"),
  attribute_aliases = list(dispersal_r_juv = "dispersal_max_distance",
                 dispersal_p_juv = "dispersal_proportion",
                 dispersal_source_n_k_cutoff = "dispersal_source_n_k$cutoff",
                 dispersal_juv = "dispersal_data"),
  decimals = 3
)
juvenile_dispersal_gen$distance_data <- adult_dispersal_gen$distance_data


## ----capacity-gen-----------------------------------------
capacity_gen <- Generator$new(
  description = "capacity",
  region = region,
  generate_rasters = FALSE,
  time_steps = timesteps,
  generative_requirements = list(
    initial_abundance = "function",
    carrying_capacity = "function"
  ),
  inputs = c("density_max", "hs", "abundance", "infected_t1"),
  outputs = c("initial_abundance", "carrying_capacity")
)
# Here we subset the habitat suitability raster to have only the region cells,
# and we add the burn in. Also, we tell the generator to generate the
# carrying_capacity based on "density_max" and "hs".
capacity_gen$add_function_template(
  param = "carrying_capacity",
  function_def = function(params) {
    hs_matrix <- params$hs |> as.matrix() |>
      _[params$region$region_indices, 1:params$time_steps]
    hs_matrix[!is.finite(hs_matrix)] <- 0
    # round the density values
    round(params$density_max * hs_matrix)
  },
  call_params = c("density_max", "hs", "region",
                  "time_steps")
)
# Here we tell the generator what function to use to generate initial_abundance
# based on the number of finches first infected in Washington, D.C.
capacity_gen$add_function_template(
  param = "initial_abundance",
  function_def = function(params) {
    abundance <- params$abundance
    infected_adults <- round(runif(1, 0, params$infected_t1))
    infected_juv <- params$infected_t1 - infected_adults
    abundance[3, 3531] <- infected_juv
    abundance[4, 3531] <- infected_adults
    return(abundance)
  },
  call_params = c("abundance", "infected_t1")
)

test_capacity <- capacity_gen$generate(input_values = list(density_max = 186000, infected_t1 = 5, hs = habitat_suitability, abundance = initial_abundance))

raster::plot(
  region$raster_from_values(test_capacity[[1]][1,]),
  main = "Initial abundance of susceptible juveniles",
  colNA = "blue"
)


## ----simulate---------------------------------------------
handler <- SimulationHandler$new(
  sample_data = sample_data,
  model_template = model_template,
  generators = list(juvenile_dispersal_gen,
                    adult_dispersal_gen,
                    capacity_gen),
  parallel_cores = parallel_cores,
  results_dir = tempdir()
)
sim_log <- handler$run()
sim_log$summary


## ----results----------------------------------------------
results <- qs::qread(file.path(tempdir(), "sample_1_results.qs"))
str(results)


## ----time-series------------------------------------------
plot(y = results$all$abundance_segments$stage_1_compartment_2[, 1],
     x = 1994:2016,
     ylab = "Abundance", xlab = "Year",
     main = "Juveniles Infected for the First Time (Summer)")


## ----map--------------------------------------------------
region$raster_from_values(results$abundance_segments$stage_2_compartment_1[ , 23, 1]) |>
  raster::plot(colNA = "blue", main = "Susceptible Adults in Summer 2016")

