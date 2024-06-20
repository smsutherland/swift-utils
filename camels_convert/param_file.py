def make_param_file(run_name, A_SN1, A_SN2, A_AGN1, A_AGN2, Omega_m, Omega_b, **kwargs):
    h = 0.6711
    Omega_cdm = Omega_m - Omega_b
    Omega_lambda = 1.0 - Omega_m
    return f"""
# Define some meta-data about the simulation
MetaData:
  run_name:   {run_name}

# Define the system of units to use internally. 
InternalUnitSystem:
  UnitMass_in_cgs:     1.98841e43    # 10^10 M_sun in grams
  UnitLength_in_cgs:   3.08567758e24 # Mpc in centimeters
  UnitVelocity_in_cgs: 1e5           # km/s in centimeters per second
  UnitCurrent_in_cgs:  1             # Amperes
  UnitTemp_in_cgs:     1             # Kelvin

# Cosmological parameters
Cosmology:
  h:              {h}        # Reduced Hubble constant
  a_begin:        0.0078125     # Initial scale-factor of the simulation
  a_end:          1.0           # Final scale factor of the simulation
  Omega_cdm:      {Omega_cdm}     # Cold Dark Matter density parameter
  Omega_lambda:   {Omega_lambda}         # Dark-energy density parameter
  Omega_b:        {Omega_b}     # Baryon density parameter

# Parameters governing the time integration
TimeIntegration:
  dt_min:     1e-10 # The minimal time-step size of the simulation (in internal units).
  dt_max:     1e-2  # The maximal time-step size of the simulation (in internal units).
  
# Parameters governing the snapshots
Snapshots:
  basename:            snapshot # Common part of the name of output files
  output_list_on:      1
  output_list:         ./output_list.txt
  #recording_triggers_part: [1.0227e-4, 1.0227e-5]   # Recording starts 100M and 10M years before a snapshot
  #recording_triggers_bpart: [1.0227e-4, 1.0227e-5]   # Recording starts 100M and 10M years before a snapshot
  recording_triggers_part: [-1, -1]   # No particle recording
  recording_triggers_bpart: [-1, -1]   # No particle recording

# Parameters governing the conserved quantities statistics
Statistics:
  delta_time:           1.01
  scale_factor_first:   0.01

# Parameters for the self-gravity scheme
Gravity:
  eta:                         0.025     # Constant dimensionless multiplier for time integration.
  MAC:                         geometric # Use the geometric opening angle condition
  theta_cr:                    0.7       # Opening angle (Multipole acceptance criterion)
  use_tree_below_softening:    0
  mesh_side_length:            256
  comoving_DM_softening:         0.0005 # Comoving softening for DM
  max_physical_DM_softening:     0.0005 # Physical softening for DM
  comoving_baryon_softening:     0.0005 # Comoving softening for baryons
  max_physical_baryon_softening: 0.0005 # Physical softening for baryons

# Parameters for the hydrodynamics scheme
SPH:
  resolution_eta:                    1.2348   # Target smoothing length in units of the mean inter-particle separation (1.2348 == 48Ngbs with the cubic spline kernel).
  h_min_ratio:                       0.25     # Minimal smoothing length in units of softening.
  h_max:                             0.5      # Maximal smoothing length in co-moving internal units.
  CFL_condition:                     0.2      # Courant-Friedrich-Levy condition for time integration.
  minimal_temperature:               100.0    # (internal units)
  initial_temperature:               268.7    # (internal units)
  particle_splitting:                1        # Particle splitting is ON
  particle_splitting_mass_threshold: 8.0e-3   # (internal units, i.e. 8e7 Msun ~ 4x initial gas particle mass)

# Parameters of the stars neighbour search
Stars:
  resolution_eta:        1.1642   # Target smoothing length in units of the mean inter-particle separation
  h_tolerance:           7e-3
  luminosity_filename:   ./photometry

BlackHoles:
  resolution_eta:        2.0
  h_tolerance:           7e-3
  h_max:                 0.00588                # 2.94 kpc (2 kpc/h in Simba)
  h_min:                 0.00147               # 0.294 kpc (10% of 2kpc/h in Simba)

# Parameters for the Friends-Of-Friends algorithm
FOF:
  basename:                        fof_output  # Filename for the FOF outputs.
  min_group_size:                  16          # The minimum no. of particles required for a group.
  linking_length_ratio:            0.05         # Linking length in units of the main inter-particle separation.
  seed_black_holes_enabled:        1           # Enable seeding of black holes in FoF groups
  black_hole_seed_host_mass_Msun:  5.8e8      # Minimal group mass in which to seed a black hole (in solar masses).
  cold_gas_temperature_threshold:  1.0e5      # Cold gas below this temperature in K
  cold_gas_n_H_threshold_cgs:      0.13       # Cold gas above this density in H/cc
  scale_factor_first:              0.05        # Scale-factor of first FoF black hole seeding calls.
  delta_time:                      1.01     # Scale-factor ratio between consecutive FoF black hole seeding calls.

Scheduler:
  max_top_level_cells:   16
  cell_split_size:       200
  
Restarts:
  onexit:       1
  delta_hours:  0.5

# Parameters related to the initial conditions
InitialConditions:
  file_name:  ic.hdf5
  periodic:   1
  cleanup_h_factors: 0               # Remove the h-factors inherited from Gadget
  cleanup_velocity_factors: 0        # Remove the sqrt(a) factor in the velocities inherited from Gadget
  generate_gas_in_ics: 0             # Generate gas particles from the DM-only ICs
  cleanup_smoothing_lengths: 0       # Since we generate gas, make use of the (expensive) cleaning-up procedure.
  remap_ids: 1                       # Re-map the IDs to [1, N] to avoid collision problems when splitting

# Impose primoridal metallicity
EAGLEChemistry:
  init_abundance_metal:     0.
  init_abundance_Hydrogen:  0.752
  init_abundance_Helium:    0.248
  init_abundance_Carbon:    0.0
  init_abundance_Nitrogen:  0.0
  init_abundance_Oxygen:    0.0
  init_abundance_Neon:      0.0
  init_abundance_Magnesium: 0.0
  init_abundance_Silicon:   0.0
  init_abundance_Iron:      0.0

# PS2020Cooling cooling parameters
PS2020Cooling:
  dir_name:                ./UV_dust1_CR1_G1_shield1.hdf5 # Location of the cooling tables
  H_reion_z:               7.5               # Redshift of Hydrogen re-ionization (Planck 2018)
  H_reion_eV_p_H:          2.0
  He_reion_z_centre:       3.5               # Redshift of the centre of the Helium re-ionization Gaussian
  He_reion_z_sigma:        0.5               # Spread in redshift of the  Helium re-ionization Gaussian
  He_reion_eV_p_H:         2.0               # Energy inject by Helium re-ionization in electron-volt per Hydrogen atom
  delta_logTEOS_subgrid_properties: 0.3      # delta log T above the EOS below which the subgrid properties use Teq assumption
  rapid_cooling_threshold:          0.333333 # Switch to rapid cooling regime for dt / t_cool above this threshold.

# SIMBA star formation parameters
SIMBAStarFormation:
  SF_threshold:                      Subgrid      # Zdep (Schaye 2004) or Subgrid
  SF_model:                          PressureLaw  # PressureLaw (Schaye et al. 2008) or SchmidtLaw
  KS_normalisation:                  1.515e-4     # The normalization of the Kennicutt-Schmidt law in Msun / kpc^2 / yr.
  KS_exponent:                       1.4          # The exponent of the Kennicutt-Schmidt law.
  min_over_density:                  10000.0        # The over-density above which star-formation is allowed.
  KS_high_density_threshold_H_p_cm3: 1e8          # Hydrogen number density above which the Kennicut-Schmidt law changes slope in Hydrogen atoms per cm^3.
  KS_high_density_exponent:          1.4          # Slope of the Kennicut-Schmidt law above the high-density threshold.
  threshold_temperature1_K:          100.0        # When using subgrid-based SF threshold, subgrid temperature below which gas is star-forming.
  threshold_temperature2_K:          1.e5        # When using subgrid-based SF threshold, subgrid temperature below which gas is star-forming if also above the density limit.
  threshold_number_density_H_p_cm3:  0.13           # When using subgrid-based SF threshold, subgrid number density above which gas is star-forming if also below the second temperature limit.
  H2_model:                          KMT    # 'Thresh' sets fH2=1; 'KMT' computes fH2 from KG11; 'Grackle' uses Grackle-computed fH2
  clumping_factor_scaling:           1.0        # scaling with resolution of KMT clumping factor
  
# Parameters for the EAGLE "equation of state"
SIMBAEntropyFloor:
  Jeans_density_threshold_H_p_cm3: 0.1      # Physical density above which the SIMBA Jeans limiter entropy floor kicks in expressed in Hydrogen atoms per cm^3.
  Jeans_over_density_threshold:    1000.       # Overdensity above which the SIMBA Jeans limiter entropy floor can kick in.
  Jeans_temperature_norm_K:        1.e4       # Temperature of the SIMBA Jeans limiter entropy floor at the density threshold expressed in Kelvin.
  Jeans_gamma_effective:           1.0       # Slope the of the SIMBA Jeans limiter entropy floor
  Cool_density_threshold_H_p_cm3: 1.e-7       # Physical density for the SIMBA Cool limiter entropy floor kicks in expressed in Hydrogen atoms per cm^3, between this value up to Jeans_density_threshold_H_p_cm3.
  Cool_over_density_threshold:    1.        # Overdensity above which the SIMBA Cool limiter entropy floor can kick in.
  Cool_temperature_norm_K:        100.        # Temperature of the SIMBA Cool limiter entropy floor at the density threshold expressed in Kelvin. (NOTE: This is below the min T and hence this floor does nothing)
  Cool_gamma_effective:           1.3333         # Slope the of the SIMBA Cool limiter entropy floor

# SIMBA feedback model (based on EAGLE-thermal)
SIMBAFeedback:
  use_SNII_feedback:                    1               # Global switch for SNII thermal (stochastic) feedback.
  use_SNIa_feedback:                    1               # Global switch for SNIa thermal (continuous) feedback.
  use_AGB_enrichment:                   1               # Global switch for enrichement from AGB stars.
  use_SNII_enrichment:                  1               # Global switch for enrichement from SNII stars.
  use_SNIa_enrichment:                  1               # Global switch for enrichement from SNIa stars.
  filename:                             ./yieldtables/  # Path to the directory containing the EAGLE yield tables.
  IMF_min_mass_Msun:                    0.1             # Minimal stellar mass considered for the Chabrier IMF in solar masses.
  IMF_max_mass_Msun:                  100.0             # Maximal stellar mass considered for the Chabrier IMF in solar masses.
  SNII_min_mass_Msun:                   8.0             # Minimal mass considered for SNII stars in solar masses.
  SNII_max_mass_Msun:                 100.0             # Maximal mass considered for SNII stars in solar masses.
  SNII_energy_erg:                      1.0e51          # Energy of one SNII explosion in ergs.
  SNIa_DTD:                             Exponential     # Functional form of the SNIa delay time distribution.
  SNIa_DTD_delay_Gyr:                   0.04            # Stellar age after which SNIa start in Gyr (40 Myr corresponds to stars ~ 8 Msun).
  SNIa_DTD_exp_timescale_Gyr:           2.0             # Time-scale of the exponential decay of the SNIa rates in Gyr.
  SNIa_DTD_exp_norm_p_Msun:             0.002           # Normalisation of the SNIa rates in inverse solar masses.
  SNIa_energy_erg:                     1.0e51           # Energy of one SNIa explosion in ergs.
  AGB_ejecta_velocity_km_p_s:          10.0             # Velocity of the AGB ejectas in km/s.
  stellar_evolution_age_cut_Gyr:        0.1             # Stellar age in Gyr above which the enrichment is down-sampled.
  stellar_evolution_sampling_rate:       10             # Number of time-steps in-between two enrichment events for a star above the age threshold.
  SNII_yield_factor_Hydrogen:           1.0             # (Optional) Correction factor to apply to the Hydrogen yield from the SNII channel.
  SNII_yield_factor_Helium:             1.0             # (Optional) Correction factor to apply to the Helium yield from the SNII channel.
  SNII_yield_factor_Carbon:             1.0             # (Optional) Correction factor to apply to the Carbon yield from the SNII channel.
  SNII_yield_factor_Nitrogen:           1.0             # (Optional) Correction factor to apply to the Nitrogen yield from the SNII channel.
  SNII_yield_factor_Oxygen:             1.0             # (Optional) Correction factor to apply to the Oxygen yield from the SNII channel.
  SNII_yield_factor_Neon:               1.0             # (Optional) Correction factor to apply to the Neon yield from the SNII channel.
  SNII_yield_factor_Magnesium:          1.0             # (Optional) Correction factor to apply to the Magnesium yield from the SNII channel.
  SNII_yield_factor_Silicon:            1.0             # (Optional) Correction factor to apply to the Silicon yield from the SNII channel.
  SNII_yield_factor_Iron:               1.0             # (Optional) Correction factor to apply to the Iron yield from the SNII channel.
  FIRE_velocity_normalization:          {1.6*A_SN2}
  FIRE_velocity_slope:                  0.12
  FIRE_eta_normalization:               {9.0*A_SN1}
  FIRE_eta_break_Msun:                  5.2e9
  FIRE_eta_lower_slope:                 -0.317
  FIRE_eta_upper_slope:                 -0.761
  early_wind_suppression_enabled:       1
  early_stellar_mass_norm_Msun:         2.9e8
  early_wind_suppression_scale_factor:  0.25
  early_wind_suppression_slope:         2.0
  minimum_galaxy_stellar_mass_Msun:     6.4e8             # Minimum mass to consider galaxy for SF. Simba: 6.4e8
  kick_velocity_scatter:                0.0               # Made to be 0 for no scatter
  wind_decouple_time_factor:            0.02
  SN_energy_scale:                      100               # never limit wind velocity due to available supernova energy

# Simba AGN model
SIMBAAGN:
  subgrid_seed_mass_Msun:             {1.0e4/h}           # Black hole subgrid mass at creation time in solar masses.
  use_multi_phase_bondi:              0               # Compute Bondi rates per neighbour particle?
  use_subgrid_bondi:                  0               # Compute Bondi rates using the subgrid extrapolation of the gas properties around the BH?
  with_angmom_limiter:                0               # Are we applying the Rosas-Guevara et al. (2015) viscous time-scale reduction term?
  with_boost_factor:                  0               # Are we using the model from Booth & Schaye (2009)?
  fixed_T_above_EoS_dex:              0.3             # Distance above the entropy floor for which we use a fixed sound-speed
  fixed_T_near_EoS_K:                 8000            # Fixed temperature assumed to compute the sound-speed of gas on the entropy floor in the Bondy-Hoyle accretion term
  radiative_efficiency:               0.1             # Fraction of the accreted mass that gets radiated.
  use_nibbling:                       1               # Continuously transfer small amounts of mass from all gas neighbours to a black hole [1] or stochastically swallow whole gas particles [0]?
  min_gas_mass_for_nibbling_Msun:          5.0e5           # Minimum mass for a gas particle to be nibbled from [M_Sun]. Only used if use_nibbling is 1.
  max_eddington_fraction:             3.              # Maximal allowed accretion rate in units of the Eddington rate.
  eddington_fraction_for_recording:   0.1             # Record the last time BHs reached an Eddington ratio above this threshold.
  coupling_efficiency:                0.1             # Fraction of the radiated energy that couples to the gas in feedback events.
  AGN_feedback_model:                 MinimumDistance # Feedback modes: Random, Isotropic, MinimumDistance, MinimumDensity
  AGN_use_deterministic_feedback:     1               # Deterministic (reservoir) [1] or stochastic [0] AGN feedback?
  use_variable_delta_T:               1               # Switch to enable adaptive calculation of AGN dT [1], rather than using a constant value [0].
  AGN_with_locally_adaptive_delta_T:  1               # Switch to enable additional dependence of AGN dT on local gas density and temperature (only used if use_variable_delta_T is 1).
  AGN_delta_T_mass_norm:              3e8             # Normalisation temperature of AGN dT scaling with BH subgrid mass [K] (only used if use_variable_delta_T is 1).
  AGN_delta_T_mass_reference:         1e8             # BH subgrid mass at which the normalisation temperature set above applies [M_Sun] (only used if use_variable_delta_T is 1).
  AGN_delta_T_mass_exponent:          0.666667        # Power-law index of AGN dT scaling with BH subgrid mass (only used if use_variable_delta_T is 1).
  AGN_delta_T_crit_factor:            1.0             # Multiple of critical dT for numerical efficiency (Dalla Vecchia & Schaye 2012) to use as dT floor (only used if use_variable_delta_T and AGN_with_locally_adaptive_delta_T are both 1).
  AGN_delta_T_background_factor:      0.0             # Multiple of local gas temperature to use as dT floor (only used if use_variable_delta_T and AGN_with_locally_adaptive_delta_T are both 1).
  AGN_delta_T_min:                    1e7             # Minimum allowed value of AGN dT [K] (only used if use_variable_delta_T is 1).
  AGN_delta_T_max:                    1e9             # Maximum allowed value of AGN dT [K] (only used if use_variable_delta_T is 1).
  AGN_delta_T_K:                      3.16228e8       # Change in temperature to apply to the gas particle in an AGN feedback event [K] (used if use_variable_delta_T is 0 or AGN_use_nheat_with_fixed_dT is 1 AND to initialise the BHs).
  AGN_use_nheat_with_fixed_dT:        0               # Switch to use the constant AGN dT, rather than the adaptive one, for calculating the energy reservoir threshold.
  AGN_use_adaptive_energy_reservoir_threshold: 0      # Switch to calculate an adaptive AGN energy reservoir threshold.
  AGN_num_ngb_to_heat:                1.              # Target number of gas neighbours to heat in an AGN feedback event (only used if AGN_use_adaptive_energy_reservoir_threshold is 0).
  max_reposition_mass:                1e30            # Maximal BH mass considered for BH repositioning in solar masses (large number implies we always reposition).
  max_reposition_distance_ratio:      4.0             # Maximal distance a BH can be repositioned, in units of the softening length.
  with_reposition_velocity_threshold: 0               # Should we only reposition to particles that move slowly w.r.t. the black hole?
  set_reposition_speed:               0               # Should we reposition black holes with (at most) a prescribed speed towards the potential minimum?
  threshold_major_merger:             0.333           # Mass ratio threshold to consider a BH merger as 'major'
  threshold_minor_merger:             0.1             # Mass ratio threshold to consider a BH merger as 'minor'
  merger_threshold_type:              EscapeVelocity               # Type of velocity threshold for BH mergers (0: v_circ at kernel edge, 1: v_esc at actual distance, with softening, 2: v_esc at actual distance, no softening).
  merger_max_distance_ratio:          3.0             # Maximal distance over which two BHs can merge, in units of the softening length.
  minimum_timestep_Myr:               0.1             # Minimum of the accretion-limited time-step length.
  jet_heating_velocity_threshold:     2000.0
  jet_velocity:                       {7000.0*A_AGN2}
  scale_jet_temperature_with_mass:    1        # use Tjet ~ MBH^2/3.
  jet_temperature:                    1.0e8    # T of jet ejected particles; if scaled, this is T at MBH=1.e9
  eddington_fraction_lower_boundary:  0.2
  jet_mass_min_Msun:                  3.1e7
  jet_mass_spread_Msun:               0.0
  environment_temperature_cut:        1.0e5
  with_potential_correction:          1
  wind_momentum_flux:                 {20.0*A_AGN1}
  f_accretion:                        0.1       # Bondi accretion rate is multiplied by this
  torque_accretion_norm:              1.0      # Torque accretion rate is multiplied by this
  xray_heating_velocity_threshold:    9999.0    # Set high to turn off x-ray feedback
  xray_maximum_heating_factor:        1000000.0    # OPTIONAL: Default 1000.0 * u_part,gas
  xray_kinetic_fraction:              0.5       # OPTIONAL: Default 0.5
  xray_heating_n_H_threshold_cgs:     0.13      # OPTIONAL: Default 0.13 cm^-3
  xray_heating_T_threshold_cgs:       3.1e5     # OPTIONAL: Default 5.0e5
  xray_shutoff_cooling:               0         # shut off cooling for dyn time after being kicked/heated by X-ray feedback
  dt_accretion_factor:                0.05      # OPTIONAL: Default 1.0, timestep limiter 5% growth
  xray_radiation_loss:                1.0       # radiative loss factor for X-ray feedback (1=no loss, 0=no xray fb)
  xray_f_gas_limit:                   0.5       # X-ray feedback only active if cold gas frac within BH kernel is lower than this, and linearly scales with lowering f_gas.
  suppress_growth:                    2         # 0=None, 1=Hopkins+21, 2=Simba-style exponential
  sigma_crit_Msun_pc2:                3000      # critical mass surface density for Hopkins+21 suppression, should be 3000
  sigma_crit_resolution_factor:       0.1       # fudge factor to multiply sigma_crit_Msun_pc2 
  bh_characteristic_suppression_mass: 5.e6      # the exponential e-folding mass for Simba BH suppression (see black_holes.h)
  bondi_rate_limiting_bh_mass:        1.e10     # Bondi rate cannot exceed the rate for BH of this mass.
  wind_decouple_time_factor:          1.e-4
  bondi_fraction_for_jet:             1.1       # Setting this above 1 to guarantee for sure that this fraction is never exceeded.
  jet_velocity_scaling_with_mass:     0
  jet_velocity_max_multiplier:        1
  jet_velocity_spread_alpha:          1         # (optional) turn off random jet velocity spread
  jet_velocity_spread_beta:           0         # (optional) turn off random jet velocity spread
"""
