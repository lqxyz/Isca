import os
import numpy as np
from isca import IscaCodeBase, DiagTable, Experiment, Namelist, GFDL_BASE

NCORES = 32 
base_dir = os.path.dirname(os.path.realpath(__file__))
cb = IscaCodeBase.from_directory(GFDL_BASE)

cb.compile()  # compile the source code to working directory $GFDL_WORK/codebase

#Tell model how to write diagnostics
diag = DiagTable()
diag.add_file('atmos_monthly', 30, 'days', time_units='days')

#Tell model which diagnostics to write
diag.add_field('atmosphere', 'precipitation', time_avg=True)
diag.add_field('atmosphere', 'convection_rain', time_avg=True)
diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
diag.add_field('atmosphere', 'rh', time_avg=True)

diag.add_field('mixed_layer', 't_surf', time_avg=True)
diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
diag.add_field('mixed_layer', 'flux_t', time_avg=True)
diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True)
diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)

diag.add_field('dynamics', 'ps', time_avg=True)
diag.add_field('dynamics', 'bk')
diag.add_field('dynamics', 'pk')
diag.add_field('dynamics', 'temp', time_avg=True)
diag.add_field('dynamics', 'sphum', time_avg=True)
diag.add_field('dynamics', 'ucomp', time_avg=True)
diag.add_field('dynamics', 'vcomp', time_avg=True)
diag.add_field('dynamics', 'omega', time_avg=True)
diag.add_field('dynamics', 'height', time_avg=True)
diag.add_field('dynamics', 'ucomp_temp', time_avg=True)
diag.add_field('dynamics', 'vcomp_temp', time_avg=True)
diag.add_field('dynamics', 'omega_temp', time_avg=True)
diag.add_field('dynamics', 'sphum_u', time_avg=True)
diag.add_field('dynamics', 'sphum_v', time_avg=True)
diag.add_field('dynamics', 'sphum_w', time_avg=True)
diag.add_field('dynamics', 'ucomp_height', time_avg=True)
diag.add_field('dynamics', 'vcomp_height', time_avg=True)
diag.add_field('dynamics', 'omega_height', time_avg=True)
diag.add_field('dynamics', 'slp', time_avg=True)

diag.add_field('rrtm_radiation', 'co2', time_avg=True)
diag.add_field('rrtm_radiation', 'flux_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'flux_lw', time_avg=True)
diag.add_field('rrtm_radiation', 'olr', time_avg=True)
diag.add_field('rrtm_radiation', 'toa_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_sw', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_lw', time_avg=True)
diag.add_field('rrtm_radiation', 'tdt_rad', time_avg=True)
diag.add_field('rrtm_radiation', 'coszen', time_avg=True)


#ALBEDO = 0.3
for ALBEDO in [0.38, 0.3, 0.27, 0.33]:
    print("Albedo is " + str(ALBEDO))
    exp = Experiment('rrtm_albedo_'+str(ALBEDO)+'_diag_mse', codebase=cb)
    exp.inputfiles = [os.path.join(GFDL_BASE,'input/rrtm_input_files/ozone_1990.nc')]
    exp.diag_table = diag
    #Empty the run directory ready to run
    exp.clear_rundir()

    #Define values for the 'core' namelist
    exp.namelist = namelist = Namelist({
        'main_nml':{
         'days'   : 30,
         'hours'  : 0,
         'minutes': 0,
         'seconds': 0,
         'dt_atmos':720,
         'current_date' : [1,1,1,0,0,0],
         'calendar' : 'thirty_day'
        },

        'idealized_moist_phys_nml': {
            'do_damping': True,
            'turb':True,
            'mixed_layer_bc':True,
            'do_virtual' :False,
            'do_simple': True,
            'roughness_mom':3.21e-05,
            'roughness_heat':3.21e-05,
            'roughness_moist':3.21e-05,                
            'two_stream_gray': False,     #Use grey radiation
            'do_rrtm_radiation': True,     #Use RRTM radiation
            'convection_scheme': 'SIMPLE_BETTS_MILLER', #Use the simple Betts Miller convection scheme from Frierson
        },

        'vert_turb_driver_nml': {
            'do_mellor_yamada': False,     # default: True
            'do_diffusivity': True,        # default: False
            'do_simple': True,             # default: False
            'constant_gust': 0.0,          # default: 1.0
            'use_tau': False
        },
        
        'diffusivity_nml': {
            'do_entrain':False,
            'do_simple': True,
        },

        'surface_flux_nml': {
            'use_virtual_temp': False,
            'do_simple': True,
            'old_dtaudv': True    
        },

        'atmosphere_nml': {
            'idealized_moist_model': True
        },

        'mixed_layer_nml': {
            'tconst' : 285.,
            'prescribe_initial_dist':True,
            'evaporation':True,   
            'depth': 10,                          #Depth of mixed layer used
            'albedo_value': ALBEDO,                  #Albedo value used             
        },

        'qe_moist_convection_nml': {
            'rhbm':0.7,
            'Tmin':160.,
            'Tmax':350.   
        },

        'betts_miller_nml': {
           'rhbm': .7   , 
           'do_simp': False, 
           'do_shallower': True, 
        },
        
        'lscale_cond_nml': {
            'do_simple':True,
            'do_evap':True
        },
        
        'sat_vapor_pres_nml': {
            'do_simple':True
        },
        
        'damping_driver_nml': {
            'do_rayleigh': True,
            'trayfric': -0.25,              # neg. value: time in *days*
            'sponge_pbottom':  5000.,           #Bottom of the model's sponge down to 50hPa (units are Pa)
            'do_conserve_energy': True,             
        },

        'rrtm_radiation_nml': {
            #'do_seasonal': True,
            'solday': 90,
            'do_read_ozone':True,
            'ozone_file':'ozone_1990',
            'dt_rad': 3600., #Set RRTM radiation timestep to 3600 seconds, meaning it runs every 5 atmospheric timesteps
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

        'spectral_dynamics_nml': {
            'damping_order': 4,             
            'water_correction_limit': 200.e2,
            'reference_sea_level_press':1.0e5,
            'num_levels':25,               #How many model pressure levels to use
            'valid_range_t':[100.,800.],
            'initial_sphum':[2.e-6],
            'vert_coord_option':'input', #Use the vertical levels from Frierson 2006
            'surf_res':0.5,
            'scale_heights' : 11.0,
            'exponent':7.0,
            'robert_coeff':0.03
        },
        'vert_coordinate_nml': {
            'bk': [0.000000, 0.0117665, 0.0196679, 0.0315244, 0.0485411, 0.0719344, 0.1027829, 0.1418581, 0.1894648, 0.2453219, 0.3085103, 0.3775033, 0.4502789, 0.5244989, 0.5977253, 0.6676441, 0.7322627, 0.7900587, 0.8400683, 0.8819111, 0.9157609, 0.9422770, 0.9625127, 0.9778177, 0.9897489, 1.0000000],
            'pk': [0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000],
           }
    })

    #Lets do a run!
    if __name__=="__main__":
        exp.run(1, use_restart=False, num_cores=NCORES)
        for i in range(2, 241):
            exp.run(i, num_cores=NCORES)
