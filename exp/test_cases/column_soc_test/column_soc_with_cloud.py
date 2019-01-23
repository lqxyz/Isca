import os
import numpy as np
from isca import ColumnSocratesCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE
import sys
sys.path.append('../column_test_case/')
from scm_interp_routine import scm_interp, global_average_lat_lon

NCORES = 1
NUM_LEVELS = 31

base_dir = os.path.dirname(os.path.realpath(__file__))
cb = ColumnSocratesCodeBase.from_directory(GFDL_BASE)
cb.compile()

ds = scm_interp(filename=os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc'),
               varname='ozone_1990', nlevels=NUM_LEVELS)
global_average_lat_lon(ds, 'ozone_1990_interp')

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('column', 'ps', time_avg=True)
diag.add_field('column', 'bk')
diag.add_field('column', 'pk')
diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('column', 'sphum', time_avg=True)
diag.add_field('column', 'ucomp', time_avg=True)
diag.add_field('column', 'vcomp', time_avg=True)
diag.add_field('column', 'temp', time_avg=True)
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'dt_ug_diffusion', time_avg=True)
diag.add_field('atmosphere', 'dt_vg_diffusion', time_avg=True)

diag.add_field('socrates', 'soc_tdt_lw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_sw', time_avg=True)
diag.add_field('socrates', 'soc_tdt_rad', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_surf_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_olr', time_avg=True)
diag.add_field('socrates', 'soc_toa_sw', time_avg=True)
diag.add_field('socrates', 'soc_flux_lw', time_avg=True)
diag.add_field('socrates', 'soc_flux_sw', time_avg=True)
diag.add_field('socrates', 'soc_ozone', time_avg=True)

diag.add_field('cloud_simple', 'cf', time_avg=True)
diag.add_field('cloud_simple', 'reff_rad', time_avg=True)
diag.add_field('cloud_simple', 'frac_liq', time_avg=True)
diag.add_field('cloud_simple', 'qcl_rad', time_avg=True)
diag.add_field('cloud_simple', 'simple_rhcrit', time_avg=True)
diag.add_field('cloud_simple', 'rh_in_cf', time_avg=True)

cf_diag_formula_names = ['sundqvist', 'smith']
for cf_diag_formula_name in cf_diag_formula_names:
    print cf_diag_formula_name
    # create an Experiment object to handle the configuration of model parameters
    exp = Experiment('column_soc_with_cloud_'+cf_diag_formula_name, codebase=cb)

    exp.diag_table = diag
    #Empty the run directory ready to run
    exp.clear_rundir()

    #Define values for the 'core' namelist
    exp.namelist = Namelist({
        'main_nml':{
        'days'   : 30,
        'hours'  : 0,
        'minutes': 0,
        'seconds': 0,
        'dt_atmos':360,
        'current_date' : [1,1,1,0,0,0],
        'calendar' : 'thirty_day'
        },

        'socrates_rad_nml': {
            'stellar_constant': 1370.,
            'lw_spectral_filename': os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_lw_ga7'),
            'sw_spectral_filename': os.path.join(GFDL_BASE,'src/atmos_param/socrates/src/trunk/data/spectra/ga7/sp_sw_ga7'),
            'do_read_ozone': False,
            #'ozone_file_name': 'ozone_1990',
            #'ozone_field_name': 'ozone_1990',
            'dt_rad': 3600,
            'store_intermediate_rad': True,
            'chunk_size': 1,
            'use_pressure_interp_for_half_levels': False,
            'tidally_locked': False,
            'solday': 90,
            'do_scm_ozone': True,
            'scm_ozone': np.squeeze(ds.ozone_1990_interp_area_av.mean('time').values).tolist(),
        },

        'atmosphere_nml': {
            'idealized_moist_model': True
        },

        'column_nml': {
            'lon_max': 1, # number of columns in longitude, default begins at lon=0.0
            'lat_max': 1, # number of columns in latitude, precise
                        # latitude can be set in column_grid_nml if only 1 lat used.
            'num_levels': 31,  # number of levels
            'initial_sphum': 1e-6,
        },

        'column_grid_nml': {
            'lat_value': np.rad2deg(np.arcsin(1/np.sqrt(3))) # set latitude to that which causes insolation in frierson p2 radiation to be insolation / 4.
            #'global_average': True # don't use this option at the moment
        },

        # set initial condition, NOTE: currently there is not an option to read in initial condition from a file.
        'column_init_cond_nml': {
            'initial_temperature': 264., # initial atmospheric temperature
            'surf_geopotential': 0.0, # applied to all columns
            'surface_wind': 5. # as described above
        },

        'idealized_moist_phys_nml': {
            'do_damping': False, # no damping in column model, surface wind prescribed
            'turb':True,
            'mixed_layer_bc':True, # need surface, how is this trying to modify the wind field? ****
            'do_simple': True,
            'roughness_mom': 3.21e-05,
            'roughness_heat':3.21e-05,
            'roughness_moist':3.21e-05,
            'two_stream_gray': False,
            'do_socrates_radiation': True,
            'convection_scheme': 'SIMPLE_BETTS_MILLER',
            'do_cloud_simple': True,
        },

        'cloud_simple_nml': {
            'simple_cca':0.0,
            'rhcsfc': 0.95,
            'rhc700': 0.7,
            'rhc200': 0.3,
            'cf_diag_formula_name': cf_diag_formula_name,
        },

        'qe_moist_convection_nml': {
            'rhbm': 0.7, # rh criterion for convection
            'Tmin': 160., # min temperature for convection scheme look up tables
            'Tmax': 350.  # max temperature for convection scheme look up tables
        },

        'lscale_cond_nml': {
            'do_simple': True, # only rain
            'do_evap': False,  # no re-evaporation of falling precipitation
        },

        'surface_flux_nml': {
            'use_virtual_temp': True, # use virtual temperature for BL stability
            'do_simple': True,
            'old_dtaudv': True
        },

        'vert_turb_driver_nml': { # DONT WANT TO USE THIS, BUT NOT DOING SO IS STOPPING MIXED LAYER FROM WORKING
            'do_mellor_yamada': False,     # default: True
            'do_diffusivity': True,        # default: False
            'do_simple': True,             # default: False
            'constant_gust': 0.0,          # default: 1.0
            'use_tau': False
        },

        'mixed_layer_nml': {
            'tconst': 285.,
            'prescribe_initial_dist': False,
            'evaporation': True,
            'depth': 2.5,                          #Depth of mixed layer used
            'albedo_value': 0.30,                  #Albedo value used
        },

        'sat_vapor_pres_nml': {
            'do_simple': True,
        },

        # FMS Framework configuration
        'diag_manager_nml': {
            'mix_snapshot_average_fields': False  # time avg fields are labelled with time in middle of window
        },
        'fms_nml': {
            'domains_stack_size': 600000                        # default: 0
        },
        'fms_io_nml': {
            'threading_write': 'single',                         # default: multi
            'fileset_write': 'single',                           # default: multi
        },
    })

    # Run the simulation
    exp.run(1, use_restart=False, num_cores=NCORES, overwrite_data=False)
    for i in range(2, 25):
        exp.run(i, num_cores=NCORES, overwrite_data=False)
