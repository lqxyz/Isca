import numpy as np
import os
from gfdl.experiment import Experiment, DiagTable
import f90nml

#Define our base experiment to compile
base_dir=os.getcwd()
GFDL_BASE        = os.environ['GFDL_BASE']


for albedo in [0.27, 0.33, 0.38]:

    print ''
    print '--------------------------------------------------------------------------------------'
    print '                                      Albedo is:'+str(albedo)
    print '--------------------------------------------------------------------------------------'
    print ''

    baseexp = Experiment('byrne_estimate_forcing_albedo_'+str(albedo)+'_co2_360', overwrite_data=False)
    #baseexp = Experiment('byrne_estimate_forcing_albedo_0.27_co2_360', overwrite_data=False)
    baseexp.inputfiles = [os.path.join(base_dir,'input/sst_zonal_mean_byrne.nc')]

    baseexp.namelist['mixed_layer_nml']['albedo_value'] = albedo                  #Albedo value used

#Tell model how to write diagnostics
    diag = DiagTable()
    diag.add_file('atmos_monthly', 30, 'days', time_units='days')

    #Tell model which diagnostics to write
    #diag.add_field('atmosphere', 'precipitation', time_avg=True)
    #diag.add_field('atmosphere', 'convection_rain', time_avg=True)
    #diag.add_field('atmosphere', 'condensation_rain', time_avg=True)
    #diag.add_field('atmosphere', 'rh', time_avg=True)
    diag.add_field('mixed_layer', 't_surf', time_avg=True)
    #diag.add_field('mixed_layer', 'flux_lhe', time_avg=True)
    #diag.add_field('mixed_layer', 'flux_t', time_avg=True)
    #diag.add_field('mixed_layer', 'flux_q', time_avg=True) # No such outpu variable in mixed_layer
    #diag.add_field('mixed_layer', 'flux_r', time_avg=True)
    #diag.add_field('mixed_layer', 'flux_u', time_avg=True)
    #diag.add_field('mixed_layer', 'flux_v', time_avg=True)
    #diag.add_field('mixed_layer', 'ml_heat_cap', time_avg=True)
    #diag.add_field('mixed_layer', 'flux_oceanq', time_avg=True)
    #diag.add_field('dynamics', 'temp', time_avg=True)
    #diag.add_field('dynamics', 'ucomp', time_avg=True)
    #diag.add_field('dynamics', 'vcomp', time_avg=True)
    #diag.add_field('dynamics', 'ucomp_temp', time_avg=True)
    #diag.add_field('dynamics', 'vcomp_temp', time_avg=True)
    #diag.add_field('dynamics', 'omega', time_avg=True)
    #diag.add_field('dynamics', 'height', time_avg=True)
    #diag.add_field('dynamics', 'height_half', time_avg=True)
    #diag.add_field('dynamics', 'sphum', time_avg=True)
    #diag.add_field('dynamics', 'pk', time_avg=True)
    #diag.add_field('dynamics', 'bk', time_avg=True)
    #diag.add_field('dynamics', 'slp', time_avg=True)
    #diag.add_field('dynamics', 'ps', time_avg=True)
    #diag.add_field('two_stream', 'co2', time_avg=True)
    diag.add_field('two_stream', 'flux_sw', time_avg=True)
    diag.add_field('two_stream', 'flux_lw', time_avg=True)
    diag.add_field('two_stream', 'swdn_sfc', time_avg=True)
    diag.add_field('two_stream', 'swdn_toa', time_avg=True)
    diag.add_field('two_stream', 'lwup_sfc', time_avg=True)
    diag.add_field('two_stream', 'lwdn_sfc', time_avg=True)
    diag.add_field('two_stream', 'olr', time_avg=True)
    diag.add_field('two_stream', 'net_lw_surf', time_avg=True)
    #diag.add_field('two_stream', 'tdt_rad', time_avg=True)
    #diag.add_field('two_stream', 'lw_dtrans', time_avg=True)

    baseexp.use_diag_table(diag)

#Turn off the full, slow radiation scheme compilation

    baseexp.disable_rrtm()

#Compile model if not already compiled
    baseexp.compile()

#Empty the run directory ready to run
    baseexp.clear_rundir()

#Define values for the 'core' namelist
    baseexp.namelist['main_nml'] = f90nml.Namelist({
         'days'   : 30,
         'hours'  : 0,
         'minutes': 0,
         'seconds': 0,
         'dt_atmos':720,
         'current_date' : [0001,1,1,0,0,0],
         'calendar' : 'thirty_day'
    })

    baseexp.namelist['mixed_layer_nml']['depth'] = 10.                         #Depth of mixed layer used
#baseexp.namelist['mixed_layer_nml']['albedo_value'] = 0.27                  #Albedo value used

# ******************Read SST from file******************#
    baseexp.namelist['mixed_layer_nml']['do_read_sst'] = True
    baseexp.namelist['mixed_layer_nml']['do_sc_sst'] = True
    baseexp.namelist['mixed_layer_nml']['sst_file'] = 'sst_zonal_mean_byrne'

#baseexp.namelist['mixed_layer_nml']['do_qflux'] = False                    # Q-flux
    baseexp.namelist['spectral_dynamics_nml']['num_levels'] = 25               #How many model pressure levels to use
    baseexp.namelist['idealized_moist_phys_nml']['two_stream_gray'] = True     #Use the simple radiation scheme
    baseexp.namelist['idealized_moist_phys_nml']['do_rrtm_radiation'] = False  #Turn off the full radiation scheme
    baseexp.namelist['idealized_moist_phys_nml']['lwet_convection'] = True     #Use the simple convection scheme
    baseexp.namelist['idealized_moist_phys_nml']['do_bm'] = False              #Turn off the more complex convection scheme

    baseexp.namelist['spectral_dynamics_nml']['vert_coord_option'] = 'input'   #Use the vertical levels from Frierson 2007
    baseexp.namelist['damping_driver_nml']['sponge_pbottom'] = 5000.           #Bottom of the model's sponge down to 50hPa
    baseexp.namelist['damping_driver_nml']['trayfric'] = -0.25                 #Drag timescale for model's sponge

    baseexp.namelist['two_stream_gray_rad_nml']['carbon_conc'] = 360.            #Set a constant CO2 value
    baseexp.namelist['two_stream_gray_rad_nml']['rad_scheme'] = 'byrne'        #Select radiation scheme to use
    baseexp.namelist['two_stream_gray_rad_nml']['bog_b'] = 1997.9              #radiation scheme parameter setting strength of water-vapour feedback
    baseexp.namelist['two_stream_gray_rad_nml']['bog_a'] = 0.1627125           #radiation scheme parameter setting amount of absorption by well-mixed gases
    baseexp.namelist['two_stream_gray_rad_nml']['do_seasonal'] = True          #Required to use perpetual equinox
    baseexp.namelist['two_stream_gray_rad_nml']['solday'] = 90                 #Set the calendar day used for perpetual equinox
    baseexp.namelist['two_stream_gray_rad_nml']['equinox_day'] = 0.75          #A calendar parameter for setting perpetual equinox

    baseexp.namelist['qflux_nml']['qflux_amp'] = 30.                           #Prescribe some ocean heat transport

    #Lets do a run!
    baseexp.runmonth(1, use_restart=False,num_cores=8, light=False)
    for i in range(2,241):
        baseexp.runmonth(i, num_cores=8, light=False)
