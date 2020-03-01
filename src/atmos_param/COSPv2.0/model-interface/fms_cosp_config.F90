module fms_cosp_config_mod
  ! This module is modified from
  ! https://github.com/CFMIP/COSPv2.0/blob/master/driver/src/cosp2_test.f90
  use cosp_kinds,          only: wp   
  USE MOD_COSP_CONFIG,     ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE,SR_BINS,    &
                                 N_HYDRO,RTTOV_MAX_CHANNELS,numMISRHgtBins,               &
                                 cloudsat_DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,         &
                                 CFODD_NDBZE,      CFODD_NICOD,                           &
                                 CFODD_BNDRE,      CFODD_NCLASS,                          &
                                 CFODD_DBZE_MIN,   CFODD_DBZE_MAX,                        &
                                 CFODD_ICOD_MIN,   CFODD_ICOD_MAX,                        &
                                 CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH,                      &
                                 WR_NREGIME,                                              &
                                 numMODISTauBins,numMODISPresBins,                        &
                                 numMODISReffIceBins,numMODISReffLiqBins,                 &
                                 numISCCPTauBins,numISCCPPresBins,numMISRTauBins,         &
                                 ntau,modis_histTau,tau_binBounds,                        &
                                 modis_histTauEdges,tau_binEdges,                         &
                                 modis_histTauCenters,tau_binCenters,ntauV1p4,            &
                                 tau_binBoundsV1p4,tau_binEdgesV1p4, tau_binCentersV1p4,  &
                                 grLidar532_histBsct,atlid_histBsct,vgrid_zu,vgrid_zl,    & 
                                 Nlvgrid_local  => Nlvgrid,                               &
                                 vgrid_z_local  => vgrid_z,cloudsat_preclvl

  ! ================== Input namelist fields ================== !
  ! Values are from https://github.com/CFMIP/COSPv2.0/blob/master/driver/run/cosp2_input_nl.txt

  implicit none

  integer ::                       & !
      Ncolumns=20,                 & ! Number of subcolumns
      Npoints_it=16,               & !5000,& ! Number of gridpoints to be processed in one iteration
      Nlvgrid=40,                  & ! Number of vertical levels for statistical outputs (USE_VGRID=.true.)
      surface_radar=0,             & ! surface=1/spaceborne=0
      cloudsat_use_gas_abs=1,      & ! Include gaseous absorption (1=yes/0=no)
      cloudsat_do_ray=0,           & ! Calculate output Rayleigh (1=yes/0=no)
      lidar_ice_type=0,            & ! Ice particle shape in lidar calculations 
                                     ! (0=ice-spheres/1=ice-non-spherical)
      overlap=3,                   & ! Overlap type: 1=max, 2=rand, 3=max/rand
      isccp_topheight=1,           & ! ISCCP cloud top height 
                                     !
                                     ! 1 = adjust top height using both a computed
                                     ! infrared brightness temperature and the visible
                                     ! optical depth to adjust cloud top pressure. Note
                                     ! that this calculation is most appropriate to compare
                                     ! to ISCCP data during sunlit hours.
                                     !
                                     ! 2 = do not adjust top height, that is cloud top
                                     ! pressure is the actual cloud top pressure
                                     ! in the model
                                     !
                                     ! 3 = adjust top height using only the computed
                                     ! infrared brightness temperature. Note that this
                                     ! calculation is most appropriate to compare to ISCCP
                                     ! IR only algortihm (i.e. you can compare to nighttime
                                     ! ISCCP data with this option)
  
      isccp_topheight_direction=2, & ! ISCCP cloud top height direction
                                     ! direction for finding atmosphere pressure level
                                     ! with interpolated temperature equal to the radiance
                                     ! determined cloud-top temperature
                                     ! 1 = find the *lowest* altitude (highest pressure) level
                                     ! with interpolated temperature equal to the radiance
                                     ! determined cloud-top temperature
                                     ! 2 = find the *highest* altitude (lowest pressure) level
                                     ! with interpolated temperature equal to the radiance
                                     ! determined cloud-top temperature. This is the
                                     ! default value since V4.0 of the ISCCP simulator.
                                     ! ONLY APPLICABLE IF top_height EQUALS 1 or 3
 
      rttov_platform=1,            & ! RTTOV: Satellite platform
      rttov_satellite=15,          & ! RTTOV: Satellite
      rttov_instrument=5,          & ! RTTOV: Instrument
      rttov_Nchannels=3              ! RTTOV: Number of channels to be computed

  real(wp) ::                      & !
      cloudsat_radar_freq=94.0,    & ! CloudSat radar frequency (GHz)
      cloudsat_k2=-1,              & ! |K|^2, -1=use frequency dependent default
      rttov_ZenAng=50.0,           & ! RTTOV: Satellite Zenith Angle
      !co2=5.241e-04,               & ! CO2 mixing ratio
      ch4=9.139e-07,               & ! CH4 mixing ratio
      n2o=4.665e-07,               & ! n2o mixing ratio
      co=2.098e-07                   ! co mixing ratio
  logical ::                       & !
      use_vgrid=.false.,            & ! Use fixed vertical grid for outputs?
      csat_vgrid=.false.,           & ! CloudSat vertical grid? (if .true. then the CloudSat standard grid
                                     ! is used for the outputs. USE_VGRID needs also be .true.)
      use_precipitation_fluxes=.false.          ! True if precipitation fluxes are input to the algorithm 
  integer,dimension(RTTOV_MAX_CHANNELS) :: &
      rttov_Channels ! =(/1, 2, 3/)             ! RTTOV: Channel numbers
  real(wp),dimension(RTTOV_MAX_CHANNELS) :: &
      rttov_Surfem !=(/0.0, 0.0, 0.0/)          ! RTTOV: Surface emissivity
  character(len=64) :: &
      cloudsat_micro_scheme='MMF_v3_single_moment'   ! Microphysical scheme used in cloudsat radar simulator
                                                      !'MMF_v3.5_two_moment'

  ! COSP time stepping and spatial sampling
  ! **** Copied from socrates_config_mod ****
  integer :: dt_cosp = 0               ! COSP timestep - every step if dt_cosp<dt_atmos
  integer :: dt_cosp_avg = -1          ! If averaging, over what time? dt_cosp_avg=dt_cosp if dt_cosp_avg<=0
  integer :: solday = 0                ! if >0, do perpetual run corresponding to day of the year = solday \in [0,days per year]
  logical :: do_cosp_time_avg = .true. ! Average coszen for SW radiation over dt_cosp?
  real    :: equinox_day = 0.75        ! fraction of the year defining NH autumn equinox \in [0,1]
  real    :: del_sol = 1.4
  real    :: del_sw  = 0.0
  logical :: tidally_locked = .false.
  logical :: frierson_solar_rad = .false.

  logical :: do_read_ozone = .false. ! read ozone from an external file?
  character(len=256) :: ozone_file_name = 'ozone' !Name of file containing ozone field - n.b. don't need to include '.nc'
  character(len=256) :: ozone_field_name = 'ozone' !Name of ozone variable in ozone file
  logical :: input_o3_file_is_mmr = .true. ! Does the ozone input file contain values as a mass mixing ratio (set to true) or a volume mixing ratio (set to false)?
  logical :: do_read_co2 = .false. ! read ozone from an external file?
  character(len=256) :: co2_file_name = 'co2' !Name of file containing co2 field - n.b. don't need to include '.nc'
  character(len=256) :: co2_field_name = 'co2' !Name of co2 variable in co2 file  
  real :: co2_ppmv = 300. !Default CO2 concentration in PPMV ! COSP default value is about 345 ppm (5.241e-04kg/kg)
  logical :: input_co2_mmr=.false. !Socrates wants input concentrations as mmr not vmr, so need to make sure input data supplied is converted if necessary


  namelist/cosp_input_nml/ overlap, isccp_topheight, isccp_topheight_direction,     &
      npoints_it, ncolumns, use_vgrid, Nlvgrid, csat_vgrid,                         &
      cloudsat_radar_freq, surface_radar, cloudsat_use_gas_abs,cloudsat_do_ray,     &
      cloudsat_k2, cloudsat_micro_scheme, lidar_ice_type, use_precipitation_fluxes, &
      rttov_platform, rttov_satellite, rttov_Instrument, rttov_Nchannels,           &
      rttov_Channels, rttov_Surfem, rttov_ZenAng, ch4, n2o, co,                     & ! co2
      dt_cosp, dt_cosp_avg, solday, do_cosp_time_avg, tidally_locked,               &
      frierson_solar_rad, equinox_day, del_sol, del_sw,                             &
      do_read_ozone, ozone_file_name, ozone_field_name, input_o3_file_is_mmr,       &
      do_read_co2, co2_file_name, co2_field_name, input_co2_mmr, co2_ppmv
  
  ! ================== Output namelist ================== !
  logical :: &
      !- CloudSat
      Lcfaddbze94=.true., &
      Ldbze94=.true., &
     
      !- CALIPSO
      Latb532=.true., &
      LcfadLidarsr532=.true., &
      Lclcalipso=.true., &
      Lclhcalipso=.true., &
      Lcllcalipso=.true., &
      Lclmcalipso=.true., &
      Lcltcalipso=.true., &
      LparasolRefl=.true., &
     
      ! CALIPSO phase diagnostics
      Lclcalipsoliq=.true., &
      Lclcalipsoice=.true., &
      Lclcalipsoun=.true., &
      Lclcalipsotmp=.true., &
      Lclcalipsotmpliq=.true., &
      Lclcalipsotmpice=.true., &
      Lclcalipsotmpun=.true., &
      Lclhcalipsoliq=.true., &
      Lcllcalipsoliq=.true., &
      Lclmcalipsoliq=.true., &
      Lcltcalipsoliq=.true., &
      Lclhcalipsoice=.true., &
      Lcllcalipsoice=.true., &
      Lclmcalipsoice=.true., &
      Lcltcalipsoice=.true., &
      Lclhcalipsoun=.true., &
      Lcllcalipsoun=.true., &
      Lclmcalipsoun=.true., &
      Lcltcalipsoun=.true., &
      
      ! CALIPSO OPAQ diagnostics
      Lclopaquecalipso=.true., &
      Lclthincalipso=.true., &
      Lclzopaquecalipso=.true., &
      Lclcalipsoopaque=.true., &
      Lclcalipsothin=.true., &
      Lclcalipsozopaque=.true., &
      Lclcalipsoopacity=.true., &
      Lclopaquetemp=.true.,  &
      Lclthintemp=.true.,  &
      Lclzopaquetemp=.true.,  &
      Lclopaquemeanz=.true.,  &
      Lclthinmeanz=.true.,  &
      Lclthinemis=.true.,  &
      Lclopaquemeanzse=.true., &
      Lclthinmeanzse=.true.,  &
      Lclzopaquecalipsose=.true., &
      
      ! GROUND LIDAR diagnostics
      LlidarBetaMol532gr=.true., &
      LcfadLidarsr532gr=.true., &
      Latb532gr=.true., &
      LclgrLidar532=.true., &
      LclhgrLidar532=.true., &
      LcllgrLidar532=.true., &
      LclmgrLidar532=.true., &
      LcltgrLidar532=.true., &
      
      ! ATLID diagnostics
      LlidarBetaMol355=.true., &
      LcfadLidarsr355=.true., &
      Latb355=.true., &
      Lclatlid=.true., &
      Lclhatlid=.true., &
      Lcllatlid=.true., &
      Lclmatlid=.true., &
      Lcltatlid=.true., &
      
      !- ISCCP
      Lalbisccp=.true., &
      Lboxptopisccp=.true., &
      Lboxtauisccp=.true., &
      Lpctisccp=.true., &
      Lclisccp=.true., &
      Ltauisccp=.true., &
      Lcltisccp=.true., &
      Lmeantbisccp=.true., &
      Lmeantbclrisccp=.true., &
      
      !- MISR
      LclMISR=.true., &

      !- Use lidar and radar
      Lclcalipso2=.true., &
      Lcltlidarradar=.true., &
      Lcloudsat_tcc=.true., &
      Lcloudsat_tcc2=.true., &

      !- These are provided for debugging or special purposes
      Lfracout=.true., &
      LlidarBetaMol532=.true., &

      !- MODIS
      Lcltmodis=.true., &
      Lclwmodis=.true., &
      Lclimodis=.true., &
      Lclhmodis=.true., &
      Lclmmodis=.true., &
      Lcllmodis=.true., &
      Ltautmodis=.true., &
      Ltauwmodis=.true., &
      Ltauimodis=.true., &
      Ltautlogmodis=.true., &
      Ltauwlogmodis=.true., &
      Ltauilogmodis=.true., &
      Lreffclwmodis=.true., &
      Lreffclimodis=.true., &
      Lpctmodis=.true., &
      Llwpmodis=.true., &
      Liwpmodis=.true., &
      Lclmodis=.true., &
      
      !- RTTOV &
      Ltbrttov=.true., &
      
      ! -CLOUDSAT precipitation frequency/occurence diagnostics
      Lptradarflag0=.true., &
      Lptradarflag1=.true., &
      Lptradarflag2=.true., &
      Lptradarflag3=.true., &
      Lptradarflag4=.true., &
      Lptradarflag5=.true., &
      Lptradarflag6=.true., &
      Lptradarflag7=.true., &
      Lptradarflag8=.true., &
      Lptradarflag9=.true., &
      Lradarpia=.true., &
    
      !- CloudSat+MODIS joint diagnostics
      Lwr_occfreq=.true., &
      Lcfodd=.true.

  namelist/cosp_output_nml/Lcfaddbze94,Ldbze94,Latb532,LcfadLidarsr532,Lclcalipso,&
                Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,LparasolRefl,     &
                Lclcalipsoliq,Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,           &
                Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,Lclhcalipsoliq, &
                Lcllcalipsoliq,Lclmcalipsoliq,Lcltcalipsoliq,Lclhcalipsoice,      &
                Lcllcalipsoice,Lclmcalipsoice,Lcltcalipsoice,Lclhcalipsoun,       &
                Lcllcalipsoun,Lclmcalipsoun,Lcltcalipsoun,Lclopaquecalipso,       &
                Lclthincalipso,Lclzopaquecalipso,Lclcalipsoopaque,Lclcalipsothin, &
                Lclcalipsozopaque,Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,    &
                Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,Lclthinemis,           &
                Lclopaquemeanzse,Lclthinmeanzse,Lclzopaquecalipsose,              &
                LlidarBetaMol532gr,LcfadLidarsr532gr,Latb532gr,LclgrLidar532,     &
                LclhgrLidar532,LcllgrLidar532,LclmgrLidar532,LcltgrLidar532,      &
                LlidarBetaMol355,LcfadLidarsr355,Latb355,Lclatlid,                &
                Lclhatlid,Lcllatlid,Lclmatlid,Lcltatlid,Lalbisccp,Lboxptopisccp,  &
                Lboxtauisccp,Lpctisccp,Lclisccp,Ltauisccp,Lcltisccp,Lmeantbisccp, &
                Lmeantbclrisccp,LclMISR,Lclcalipso2,Lcltlidarradar,               &
                Lcloudsat_tcc, Lcloudsat_tcc2, Lfracout,                          &
                LlidarBetaMol532,Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,         &
                Lclmmodis,Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,             &
                Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,Lreffclwmodis,          &
                Lreffclimodis,Lpctmodis,Llwpmodis,Liwpmodis,Lclmodis,Ltbrttov,    &
                Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,          &
                Lptradarflag4,Lptradarflag5,Lptradarflag6,Lptradarflag7,          &
                Lptradarflag8,Lptradarflag9,Lradarpia,                            &
                Lwr_occfreq, Lcfodd      

end module fms_cosp_config_mod
