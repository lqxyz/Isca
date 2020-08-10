module fms_cosp_interface_mod
  ! This module is modified from
  ! https://github.com/CFMIP/COSPv2.0/blob/master/driver/src/cosp2_test.f90

  ! COSP calculation interface module
  ! Takes FMS time, temperature, and pressure, cloud properties
  ! Outputs FMS cloud
  ! QL, 02/2020

#ifdef INTERNAL_FILE_NML
  use mpp_mod, only: input_nml_file
#else
  use fms_mod, only: open_namelist_file, close_file
#endif

  use diag_manager_mod,     only: register_diag_field, send_data, diag_axis_init
  use time_manager_mod,     only: time_type, OPERATOR(+), OPERATOR(-), OPERATOR(/=), &
                                  length_of_day, length_of_year, get_time, set_time
  use fms_mod,              only: stdlog, FATAL, WARNING, NOTE, error_mesg
  use interpolator_mod,     only: interpolator_init, interpolator_end, &
                                  interpolate_type, interpolator, ZERO
  use astronomy_mod,        only: astronomy_init, diurnal_solar
  use constants_mod,        only: grav, pi, rdgas, gas_constant ! wtmco2, wtmozone

  ! COSP modules
  use cosp_kinds,           only: wp
  use mod_cosp,             only: cosp_init, cosp_optical_inputs, cosp_column_inputs, &
                                  cosp_outputs, cosp_cleanUp, cosp_simulator
  use mod_quickbeam_optics, only: size_distribution, hydro_class_init, quickbeam_optics, &
                                  quickbeam_optics_init, gases
  use quickbeam,            only: radar_cfg
  use cosp_phys_constants,  only: amO3, amCO2
  use mod_rng,              only: rng_state, init_rng
  use mod_scops,            only: scops
  use mod_prec_scops,       only: prec_scops
  use mod_cosp_utils,       only: cosp_precip_mxratio
  use cosp_optics,          only: cosp_simulator_optics,lidar_optics,modis_optics,  &
                                  modis_optics_partition
  use mod_cosp_stats,       only: cosp_change_vertical_grid

  ! Sef defined modules
  use fms_cosp_config_mod  ! Flags to select which simulator to use
  use fms_cosp_io_mod,      only:  map_point_to_ll  ! ,map_ll_to_point

  implicit none

  type(size_distribution)   :: sd             ! Hydrometeor description
  type(radar_cfg)           :: rcfg_cloudsat  ! Radar configuration
  type(cosp_outputs)        :: cospOUT        ! COSP simulator outputs
  type(cosp_optical_inputs) :: cospIN         ! COSP optical (or derived?) fields needed by simulators
  type(cosp_column_inputs)  :: cospstateIN    ! COSP model fields needed by simulators
  character(len=256), dimension(100) :: cosp_status

  character(len=10), parameter  :: mod_name = 'cosp'
  real :: dt_last !Time of last COSP calculation - used to tell whether it is time to recompute COSP or not
  type(interpolate_type),save   :: o3_interp, co2_interp  ! use external file for ozone and co2

  ! ================== Local variables ================== !
  logical :: &
        lsingle     = .true., & ! true if using mmf_v3_single_moment cloudsat microphysical scheme (default)
        ldouble     = .false., & ! true if using mmf_v3.5_two_moment cloudsat microphysical scheme
        lisccp      = .false., & ! local on/off switch for simulators (used by initialization)
        lmodis      = .false., &
        lmisr       = .false., &
        lcalipso    = .false., &
        lgrlidar532 = .false., &
        latlid      = .false., &
        lcloudsat   = .false., &
        lrttov      = .false., &
        lparasol    = .false.

    ! Indices to address arrays of LS and CONV hydrometeors
    integer,parameter :: &
        I_LSCLIQ = 1, & ! Large-scale (stratiform) liquid
        I_LSCICE = 2, & ! Large-scale (stratiform) ice
        I_LSRAIN = 3, & ! Large-scale (stratiform) rain
        I_LSSNOW = 4, & ! Large-scale (stratiform) snow
        I_CVCLIQ = 5, & ! Convective liquid
        I_CVCICE = 6, & ! Convective ice
        I_CVRAIN = 7, & ! Convective rain
        I_CVSNOW = 8, & ! Convective snow
        I_LSGRPL = 9    ! Large-scale (stratiform) groupel

   ! Stratiform and convective clouds in frac_out (scops output).
   integer, parameter :: &
        I_LSC = 1, & ! Large-scale clouds
        I_CVC = 2    ! Convective clouds

   ! Microphysical settings for the precipitation flux to mixing ratio conversion
   real(wp),parameter,dimension(N_HYDRO) :: &
                  ! LSL   LSI      LSR       LSS   CVL  CVI      CVR       CVS       LSG
        N_ax    = (/-1., -1.,     8.e6,     3.e6, -1., -1.,     8.e6,     3.e6,     4.e6/),&
        N_bx    = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
        alpha_x = (/-1., -1.,      0.0,      0.0, -1., -1.,      0.0,      0.0,      0.0/),&
        c_x     = (/-1., -1.,    842.0,     4.84, -1., -1.,    842.0,     4.84,     94.5/),&
        d_x     = (/-1., -1.,      0.8,     0.25, -1., -1.,      0.8,     0.25,      0.5/),&
        g_x     = (/-1., -1.,      0.5,      0.5, -1., -1.,      0.5,      0.5,      0.5/),&
        a_x     = (/-1., -1.,    524.0,    52.36, -1., -1.,    524.0,    52.36,   209.44/),&
        b_x     = (/-1., -1.,      3.0,      3.0, -1., -1.,      3.0,      3.0,      3.0/),&
        gamma_1 = (/-1., -1., 17.83725, 8.284701, -1., -1., 17.83725, 8.284701, 11.63230/),&
        gamma_2 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/),&
        gamma_3 = (/-1., -1.,      2.0,      2.0, -1., -1.,      2.0,      2.0,      2.0/),&
        gamma_4 = (/-1., -1.,      6.0,      6.0, -1., -1.,      6.0,      6.0,      6.0/)

  integer :: Npoints  ! Number of gridpoints ==> (nlat x nlon)
  integer :: Nlevels  ! Number of model vertical levels

  ! ************************* Output variables ***********************************!
  integer ::  id_cltcalipso_sat, id_cllcalipso_sat, id_clmcalipso_sat,  &
              id_clhcalipso_sat
  integer ::  id_cltcalipso, id_cllcalipso, id_clmcalipso, id_clhcalipso, &
              id_cltlidarradar, id_cltisccp, id_ctpisccp, id_tauisccp, &
              id_tbisccp, id_tbclrisccp, &
              id_betamol532, &
              id_albisccp, id_clcalipso, id_clcalipso2, &
              id_clcalipso_sat, id_clcalipso2_sat, &
              id_clcalipso_mdl, id_clcalipso2_mdl, &
              id_boxtauisccp, id_boxptopisccp, id_parasolrefl, &
              id_parasolrefl_sat, &
              id_sampling_sat, id_location_sat, id_lat_sat, id_lon_sat
  integer ::  id_tclmodis, id_lclmodis, id_iclmodis, id_ttaumodis, &
              id_ltaumodis, id_itaumodis, id_tlogtaumodis, &
              id_llogtaumodis, id_ilogtaumodis, id_lremodis, &
              id_badlremodis, id_badiremodis, &
              id_locldmodis, id_mdcldmodis, id_hicldmodis, &
              id_iremodis, id_ctpmodis, id_lwpmodis, id_iwpmodis
  integer, allocatable, dimension(:) :: id_dbze94, id_cloudsatcfad, &
                                        id_cloudsatcfad_sat, &
                                        id_atb532, id_calipsosrcfad, &
                                        id_calipsosrcfad_sat, &
                                        id_cloud_type, id_boxtauisccp_n, &
                                        id_boxptopisccp_n, &
                                        id_taumodis_n, id_ptopmodis_n, &
                                        id_badsizemodis_n, &
                                        id_sizemodis_n, id_phasemodis_n
  integer, allocatable, dimension(:) :: id_cloudsatcfad_mdl, &
                                        id_calipsosrcfad_mdl
  real  :: missing_value = R_UNDEF
  integer :: geomode

  contains

  subroutine fms_cosp_init(is, ie, js, je, num_levels, axes, Time, lonb, latb, delta_t_atmos) !, do_cloud_simple)

    implicit none

    ! arguments
    integer, intent(in), dimension(4) :: axes
    type(time_type), intent(in)       :: Time, delta_t_atmos
    integer, intent(in)               :: is, ie, js, je, num_levels
    real, intent(in) , dimension(:,:) :: lonb, latb
    ! logical, intent(in)               :: do_cloud_simple

    integer :: io, stdlog_unit
    integer :: res, time_step_seconds
    real    :: day_in_s_check

    ! Reads COSP option from 'cosp_input_nml' and 'cosp_output_nml',
    ! which are defined in 'fms_cosp_config_mod.F90'.
#ifdef INTERNAL_FILE_NML
    read (input_nml_file, nml=cosp_input_nml, iostat=io)
    read (input_nml_file, nml=cosp_output_nml, iostat=io)
#else
    if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      read (nml_unit, cosp_input_nml, iostat=io)
      call close_file(nml_unit)
    endif

    if ( file_exist('input.nml') ) then
      nml_unit = open_namelist_file()
      read (nml_unit, cosp_output_nml, iostat=io)
      call close_file(nml_unit)
    endif
#endif
    stdlog_unit = stdlog()
    write(stdlog_unit, cosp_input_nml)
    write(stdlog_unit, cosp_output_nml)

    !Initialise astronomy
    call astronomy_init

    !Initialise variables related to timestep
    call get_time(delta_t_atmos, time_step_seconds)

    if (dt_cosp .le. 0.) then
        dt_cosp = time_step_seconds !Make sure that dt_cosp is set if it is not specified in the namelist
    endif

    dt_last = -real(dt_cosp) !make sure we are calling COSP at the first time step

    if (dt_cosp .gt. time_step_seconds) then
        res = mod(dt_cosp, time_step_seconds)

        if (res.ne.0) then
            call error_mesg( 'fms_cosp_init', &
                'dt_cosp must be an integer multiple of dt_atmos', FATAL)
        endif

        day_in_s_check = length_of_day()
        res = mod(int(day_in_s_check), dt_cosp)

        if (res.ne.0) then
            call error_mesg( 'fms_cosp_init', &
                'dt_cosp does not fit into one day an integer number of times', WARNING)
        endif
    endif

    if(dt_cosp_avg .le. 0) dt_cosp_avg = dt_cosp

    if(do_read_ozone)then
      call interpolator_init(o3_interp, trim(ozone_file_name)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
    endif

    if(do_read_co2)then
      call interpolator_init(co2_interp, trim(co2_file_name)//'.nc', lonb, latb, data_out_of_bounds=(/ZERO/))
    endif

    ! if (mod((size(lonb,1)-1)*(size(latb,1)-1), chunk_size) .ne. 0) then
    !   call error_mesg( 'socrates_init', &
    !     'chunk_size must equally divide number of points per processor, which it currently does not.', FATAL)
    ! endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Initialize COSP
    !*This only needs to be done the first time that COSP is called.*
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    !!!! calculation !!!
    Nlevels = num_levels
    Npoints = (ie-is+1) * (je-js+1)

    call select_simulator()

    ! Initialize quickbeam_optics, also if two-moment radar microphysics scheme is wanted...
    if (cloudsat_micro_scheme == 'MMF_v3.5_two_moment')  then
      ldouble = .true.
      lsingle = .false.
    endif
    call quickbeam_optics_init()

    ! Initialize the distributional parameters for hydrometeors in radar simulator
    call hydro_class_init(lsingle, ldouble, sd)

    !write(*,*) 'QL: before call cosp init...'
    ! Initialize COSP simulator
    !  SUBROUTINE COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, LgrLidar532,
    !      Latlid, Lparasol, Lrttov,     &
    !     cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs, cloudsat_do_ray,   &
    !     isccp_top_height, isccp_top_height_direction, surface_radar, rcfg, lusevgrid, &
    !     luseCSATvgrid, Nvgrid, Nlevels, cloudsat_micro_scheme)
    ! from cosp_input_nml except 'rcfg_cloudsat'
    call cosp_init(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, LgrLidar532, Latlid,     &
        Lparasol, Lrttov,                                                               &
        cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs,                         &
        cloudsat_do_ray, isccp_topheight, isccp_topheight_direction, surface_radar,     &
        rcfg_cloudsat, use_vgrid, csat_vgrid, Nlvgrid, Nlevels, cloudsat_micro_scheme)
    !write(*,*) 'QL: after call cosp init...'

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Construct output derived type.
    ! *NOTE* The "construct/destroy" subroutines are local to this module and should be
    !        modified for your configuration. E.g. it may be overkill to query each field.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call construct_cosp_outputs(Lpctisccp, Lclisccp, Lboxptopisccp, Lboxtauisccp,         &
        Ltauisccp, Lcltisccp, Lmeantbisccp, Lmeantbclrisccp, Lalbisccp, LclMISR,          &
        Lcltmodis, Lclwmodis, Lclimodis, Lclhmodis, Lclmmodis, Lcllmodis, Ltautmodis,     &
        Ltauwmodis, Ltauimodis, Ltautlogmodis, Ltauwlogmodis, Ltauilogmodis,              &
        Lreffclwmodis, Lreffclimodis, Lpctmodis, Llwpmodis, Liwpmodis, Lclmodis, Latb532, &
        Latb532gr, Latb355, LlidarBetaMol532, LlidarBetaMol532gr, LlidarBetaMol355,       &
        LcfadLidarsr532, LcfadLidarsr532gr, LcfadLidarsr355, Lclcalipso2,                 &
        Lclcalipso, LclgrLidar532, Lclatlid, Lclhcalipso, Lcllcalipso, Lclmcalipso,       &
        Lcltcalipso, LclhgrLidar532, LcllgrLidar532, LclmgrLidar532, LcltgrLidar532,      &
        Lclhatlid, Lcllatlid, Lclmatlid, Lcltatlid, Lcltlidarradar,  Lcloudsat_tcc,       &
        Lcloudsat_tcc2, Lclcalipsoliq,                                                    &
        Lclcalipsoice, Lclcalipsoun, Lclcalipsotmp, Lclcalipsotmpliq, Lclcalipsotmpice,   &
        Lclcalipsotmpun, Lcltcalipsoliq, Lcltcalipsoice, Lcltcalipsoun, Lclhcalipsoliq,   &
        Lclhcalipsoice, Lclhcalipsoun, Lclmcalipsoliq, Lclmcalipsoice, Lclmcalipsoun,     &
        Lcllcalipsoliq, Lcllcalipsoice, Lcllcalipsoun, Lclopaquecalipso, Lclthincalipso,  &
        Lclzopaquecalipso, Lclcalipsoopaque, Lclcalipsothin, Lclcalipsozopaque,           &
        Lclcalipsoopacity, Lclopaquetemp, Lclthintemp, Lclzopaquetemp, Lclopaquemeanz,    &
        Lclthinmeanz, Lclthinemis, Lclopaquemeanzse, Lclthinmeanzse, Lclzopaquecalipsose, &
        LcfadDbze94, Ldbze94, Lparasolrefl,                                               &
        Ltbrttov, Lptradarflag0,Lptradarflag1,Lptradarflag2,Lptradarflag3,Lptradarflag4,  &
        Lptradarflag5,Lptradarflag6,Lptradarflag7,Lptradarflag8,Lptradarflag9,Lradarpia,  &
        Lwr_occfreq, Lcfodd,                                                              &
        Npoints, Ncolumns, Nlevels, Nlvgrid_local, rttov_Nchannels, cospOUT)

    !!!! ************ Define ids of output variables ************ !!!!
    call diag_field_init(Time, axes)

    !--------------------------------------------------------------------
    !   variable geomode indicates that the grid (i,j) => (lon,lat)
    !--------------------------------------------------------------------
    ! See am3/src/atmos_param/cosp/cosp_driver.F90 for details
    geomode = 2 ! geomode=2 for (lon,lat) mode; geomode=3 for (lat,lon) mode.

    if (js==1) then
      call error_mesg( 'fms_cosp_init', 'COSP v2.0 is initialized.', NOTE)
    endif

  end subroutine fms_cosp_init


  subroutine fms_cosp_interface(rad_lat, rad_lon, &
      fms_p_full, fms_p_half, fms_z_full, fms_z_half, fms_temp, &
      fms_spec_hum, fms_rh, fms_tot_cld_amt, fms_conv_cld_amt, &
      fms_mmr_ls_liq, fms_mmr_ls_ice, fms_mmr_cc_liq, fms_mmr_cc_ice, &
      fms_flux_ls_rain, fms_flux_ls_snow, fms_flux_ls_graupel, &
      fms_flux_cc_rain, fms_flux_cc_snow, fms_dtau_strat, fms_dtau_conv, &
      fms_dlw_emissivity_strat, fms_dlw_emissivity_conv, fms_mmr_ozone, &
      fms_mmr_co2, fms_hydrometeor_Reff, fms_t_surf, fms_landmask, &
      fms_u_wind, fms_v_wind, fms_sunlit, fms_z_surf, fms_sfc_lw_emissivity)

    implicit none

    real(wp), dimension(:,:),     intent(in) :: rad_lon, rad_lat
    real(wp), dimension(:,:,:),   intent(in) :: fms_p_full, fms_z_full,        &
        fms_temp, fms_spec_hum, fms_rh, fms_tot_cld_amt, fms_conv_cld_amt,     &
        fms_mmr_ls_liq, fms_mmr_ls_ice, fms_mmr_cc_liq, fms_mmr_cc_ice,        &
        fms_dtau_strat, fms_dtau_conv,fms_dlw_emissivity_strat,                &
        fms_dlw_emissivity_conv, fms_mmr_ozone, fms_mmr_co2, fms_flux_ls_rain, &
        fms_flux_ls_snow, fms_flux_ls_graupel, fms_flux_cc_rain, fms_flux_cc_snow
    real(wp), dimension(:,:,:),   intent(in) :: fms_p_half, fms_z_half
    real(wp), dimension(:,:,:,:), intent(in) :: fms_hydrometeor_Reff !(lat,lon,lev,nhydro)
    real(wp), dimension(:,:),     intent(in) ::  fms_t_surf, fms_z_surf, &
                                  fms_landmask, fms_u_wind, fms_v_wind, fms_sunlit
    real(wp), intent(in) :: fms_sfc_lw_emissivity

    ! --------------- Send to COSP --------------- !
    real(wp), dimension(Npoints) :: lon, lat
    real(wp), dimension(Npoints, Nlevels) :: p, zlev, T, sh, rh, tca, cca,    &
              mr_lsliq, mr_lsice, mr_ccliq, mr_ccice, dtau_s, dtau_c, dem_s,  &
              dem_c, mr_ozone, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
    real(wp), dimension(Npoints, Nlevels+1) :: ph, zlev_half
    real(wp), dimension(Npoints, Nlevels, N_HYDRO) :: Reff
    real(wp), dimension(Npoints) :: skt, landmask, u_wind, v_wind, sunlit, surfelev
    real(wp) :: emsfc_lw, mr_co2

    ! local variables
    integer :: iChunk, nChunks, start_idx, end_idx, nPtsPerIt, ij

    integer :: si, sj, sk

    si = size(fms_temp,1)
    sj = size(fms_temp,2)
    sk = size(fms_temp,3)

    ! Reshap the FMS variables to COSP version
    lat = reshape(rad_lat, (/si*sj/))
    lon = reshape(rad_lon, (/si*sj/))

    p = reshape(fms_p_full, (/si*sj, sk/))
    ph = reshape(fms_p_half, (/si*sj, sk+1/)) !! sk or sk+1, check again !!!!!

    ! *********** z_full is geopotential height??? should be height!!!********!
    zlev = reshape(fms_z_full, (/si*sj, sk/))
    zlev_half = reshape(fms_z_half, (/si*sj, sk+1/)) !! sk or sk+1, check again !!!!!

    T = reshape(fms_temp, (/si*sj, sk/))
    sh = reshape(fms_spec_hum, (/si*sj, sk/))
    rh = reshape(fms_rh, (/si*sj, sk/))

    ! Cloud properties
    tca = reshape(fms_tot_cld_amt, (/si*sj, sk/))
    cca = reshape(fms_conv_cld_amt, (/si*sj, sk/))
    mr_lsliq = reshape(fms_mmr_ls_liq, (/si*sj, sk/))
    mr_lsice = reshape(fms_mmr_ls_ice, (/si*sj, sk/))
    mr_ccliq = reshape(fms_mmr_cc_liq, (/si*sj, sk/))
    mr_ccice = reshape(fms_mmr_cc_ice, (/si*sj, sk/))

    fl_lsrain = reshape(fms_flux_ls_rain, (/si*sj, sk/))
    fl_lssnow = reshape(fms_flux_ls_snow, (/si*sj, sk/))
    fl_lsgrpl = reshape(fms_flux_ls_graupel, (/si*sj, sk/))
    fl_ccrain = reshape(fms_flux_cc_rain, (/si*sj, sk/))
    fl_ccsnow = reshape(fms_flux_cc_snow, (/si*sj, sk/))

    dtau_s = reshape(fms_dtau_strat, (/si*sj, sk/))
    dtau_c = reshape(fms_dtau_conv, (/si*sj, sk/))
    dem_s = reshape(fms_dlw_emissivity_strat, (/si*sj, sk/))
    dem_c = reshape(fms_dlw_emissivity_conv, (/si*sj, sk/))

    mr_ozone = reshape(fms_mmr_ozone, (/si*sj, sk/))
    mr_co2 = sum(fms_mmr_co2)/(si*sj*sk)
    !write(*,*) 'QL, mr_co2', mr_co2

    ! fms_hydrometeor_Reff(nlat, nlon, nlev, N_HYDRO)
    Reff = reshape(fms_hydrometeor_Reff,  (/si*sj, sk, N_HYDRO/))

    skt = reshape(fms_t_surf, (/si*sj/))
    landmask = reshape(fms_landmask, (/si*sj/))
    u_wind = reshape(fms_u_wind, (/si*sj/))
    v_wind = reshape(fms_v_wind, (/si*sj/))
    sunlit = reshape(fms_sunlit, (/si*sj/))  ! check again!
    surfelev = reshape(fms_z_surf, (/si*sj/))
    emsfc_lw = fms_sfc_lw_emissivity

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Break COSP up into pieces and loop over each COSP 'chunk'.
    ! nChunks = # Points to Process (nPoints) / # Points per COSP iteration (nPoints_it)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nChunks = nPoints/nPoints_it+1
    if (nPoints .eq. nPoints_it) nChunks = 1
    !write(*,*) 'QL, nChunks, nPoints', nChunks, nPoints

    do iChunk=1,nChunks
      !write(*,*) 'QL, iChunk', iChunk
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Determine indices for "chunking" (again, if necessary)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (nChunks .eq. 1) then
        start_idx = 1
        end_idx   = nPoints
        nPtsPerIt = nPoints
      else
        start_idx = (iChunk-1)*nPoints_it+1
        end_idx   = iChunk*nPoints_it
        if (end_idx .gt. nPoints) end_idx=nPoints
        nPtsPerIt = end_idx-start_idx+1
      endif

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Construct COSP input types
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (iChunk .eq. 1) then
        call construct_cospIN(Nptsperit, nColumns, nLevels, cospIN)
        call construct_cospstateIN(Nptsperit, nLevels, rttov_nChannels, cospstateIN)
      endif
      if (iChunk .eq. nChunks) then
        call destroy_cospIN(cospIN)
        call destroy_cospstateIN(cospstateIN)
        call construct_cospIN(Nptsperit, nColumns, nLevels, cospIN)
        call construct_cospstateIN(Nptsperit, nLevels, rttov_nChannels, cospstateIN)
      endif

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Populate input types with model fields.
      ! Here the 3D sample model fields (temperature,pressure,etc...) are ordered from the
      ! surface-2-TOA, whereas COSP expects all fields to be ordered from TOA-2-SFC. So the
      ! vertical fields are flipped prior to storing to COSP input type.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Should get input from model....
      cospIN%emsfc_lw         = emsfc_lw
      cospIN%rcfg_cloudsat    = rcfg_cloudsat
      cospstateIN%hgt_matrix  = zlev(start_idx:end_idx,1:Nlevels) ! km ! It should be m, right??
      cospstateIN%sunlit      = sunlit(start_idx:end_idx)         ! 0-1
      cospstateIN%skt         = skt(start_idx:end_idx)            ! K
      cospstateIN%surfelev    = surfelev(start_idx:end_idx)       ! m
      ! In the mask: 1 for land, 0 for sea (refer to lines 463-468 in quickbeam.F90)
      cospstateIN%land        = landmask(start_idx:end_idx)       ! 0-1 (*note* model specific)
      cospstateIN%qv          = sh(start_idx:end_idx,1:Nlevels)   ! kg/kg
      cospstateIN%at          = T(start_idx:end_idx,1:Nlevels)    ! K
      ! Pressure at interface (nlevels+1). Set uppermost interface to 0.
      cospstateIN%pfull       = p(start_idx:end_idx,1:Nlevels)    ! Pa
      cospstateIN%phalf(:,2:Nlevels+1) = ph(start_idx:end_idx,2:Nlevels+1) ! Pa
      cospstateIN%phalf(:,1)  = 0._wp
      ! Height at interface (nlevels+1). Set lowermost interface to 0.
      cospstateIN%hgt_matrix_half(:,1:Nlevels) = zlev_half(start_idx:end_idx,2:Nlevels+1) !m not km, right?
      cospstateIN%hgt_matrix_half(:,Nlevels+1) = 0._wp
      ! Update o3 and co2
      cospstateIN%o3          = mr_ozone(start_idx:end_idx,1:Nlevels) ! kg/kg
      cospstateIN%co2         = mr_co2                                ! kg/kg

      ! write(*,*) 'QL init IN test level order, p:',  cospstateIN%pfull(2,1:Nlevels)
      ! write(*,*) 'ph=', ph(2,1:Nlevels+1)
      ! write(*,*) 'QL init IN test level order, p_half:', cospstateIN%phalf(2,1:Nlevels+1)
      ! write(*,*) 'QL init IN test level order, zlev:',  cospstateIN%hgt_matrix(2,1:Nlevels)
      ! write(*,*) 'QL zlev_half(2,1:Nlevels+1)', zlev_half(2,1:Nlevels+1)
      ! write(*,*) 'QL init IN test level order, zlev_half:', cospstateIN%hgt_matrix_half(2,1:Nlevels+1)

      !write(*,*) 'QL: before call subsample_and_optics...'
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Generate subcolumns and compute optical inputs.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      call subsample_and_optics(nPtsPerIt,nLevels,nColumns,N_HYDRO,overlap,                    &
          use_precipitation_fluxes,lidar_ice_type,sd,                                          &
          tca(start_idx:end_idx,Nlevels:1:-1),cca(start_idx:end_idx,Nlevels:1:-1),             &
          fl_lsrain(start_idx:end_idx,Nlevels:1:-1),fl_lssnow(start_idx:end_idx,Nlevels:1:-1), &
          fl_lsgrpl(start_idx:end_idx,Nlevels:1:-1),fl_ccrain(start_idx:end_idx,Nlevels:1:-1), &
          fl_ccsnow(start_idx:end_idx,Nlevels:1:-1),mr_lsliq(start_idx:end_idx,Nlevels:1:-1),  &
          mr_lsice(start_idx:end_idx,Nlevels:1:-1),mr_ccliq(start_idx:end_idx,Nlevels:1:-1),   &
          mr_ccice(start_idx:end_idx,Nlevels:1:-1),Reff(start_idx:end_idx,Nlevels:1:-1,:),     &
          dtau_c(start_idx:end_idx,nLevels:1:-1),dtau_s(start_idx:end_idx,nLevels:1:-1),       &
          dem_c(start_idx:end_idx,nLevels:1:-1),dem_s(start_idx:end_idx,nLevels:1:-1),         &
          cospstateIN,cospIN)

      !write(*,*) 'QL: after call subsample_and_optics...'
      !write(*,*) 'QL: before call COSP_SIMULATOR...'
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Call COSP
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      cosp_status = COSP_SIMULATOR(cospIN, cospstateIN, cospOUT, start_idx, end_idx, .false.)
      !write(*,*) 'QL: after call COSP_SIMULATOR...'
      do ij=1,size(cosp_status,1)
        if (cosp_status(ij) .ne. '') print*,trim(cosp_status(ij))
      end do

    end do  ! iChunk
  end subroutine fms_cosp_interface


  subroutine fms_run_cosp(Time, Time_diag, is, ie, js, je, rad_lat, rad_lon, p_full_in,  &
                p_half_in, z_full_in, z_half_in, temp_in, q_in, rh_in, cf_rad, reff_rad, &
                qcl_rad, t_surf_in, landmask_in, u_wind_in, v_wind_in, z_surf_in,        &
                ls_cloud_absorptivity, cnv_cloud_absorptivity,                           &
                ls_cloud_extinction,   cnv_cloud_extinction)
    implicit none

    ! Input time
    type(time_type), intent(in)        :: Time, Time_diag
    integer, intent(in)                :: is, ie, js, je   ! might use, otherwise delete them
    real, intent(in), dimension(:,:)   :: rad_lat, rad_lon, t_surf_in, landmask_in, &
                                          u_wind_in, v_wind_in, z_surf_in
    real, intent(in), dimension(:,:,:) :: temp_in, p_full_in, q_in, rh_in, z_full_in,    &
                                          cf_rad, reff_rad, qcl_rad,                     &
                                          ls_cloud_absorptivity, cnv_cloud_absorptivity, &
                                          ls_cloud_extinction, cnv_cloud_extinction
    real, intent(in), dimension(:,:,:) :: p_half_in, z_half_in

    ! -------------- variables passed to the COSP interface  -------------- !
    real(wp), dimension(size(temp_in,1),size(temp_in,2),size(temp_in,3)) ::  &
                    p_full_cosp, z_full_cosp, temp_cosp, q_cosp,             &
                    rh_cosp, cf_rad_cosp, reff_rad_cosp, qcl_rad_cosp,       &
                    ls_cloud_absorptivity_cosp, cnv_cloud_absorptivity_cosp, &
                    ls_cloud_extinction_cosp, cnv_cloud_extinction_cosp,     &
                    ls_cloud_tau_cosp, cnv_cloud_tau_cosp
    real(wp), dimension(size(temp_in,1),size(temp_in,2),size(temp_in,3)+1) :: p_half_cosp, z_half_cosp
    real(wp), dimension(size(temp_in,1),size(temp_in,2)) :: &
        t_surf_cosp, landmask_cosp, u_wind_cosp, v_wind_cosp, z_surf_cosp
    real(wp), dimension(size(temp_in,1),size(temp_in,2),size(temp_in,3)) :: &
        tot_cld_amt, conv_cld_amt, mmr_ls_liq, mmr_ls_ice, mmr_cc_liq, mmr_cc_ice, dtau_strat, &
        dtau_conv, dlw_emissivity_strat, dlw_emissivity_conv, ozone_cosp, co2_cosp, &
        flux_ls_rain, flux_ls_snow, flux_ls_graupel, flux_cc_rain, flux_cc_snow
    real(wp), dimension(size(temp_in,1),size(temp_in,2)) :: rad_lat_cosp, rad_lon_cosp, sunlit_cosp
    real(wp), dimension(size(temp_in,1),size(temp_in,2),size(temp_in,3),N_HYDRO) :: hydrometeor_Reff
    real(wp) :: sfc_lw_emissivity

    ! -------------- other local variables -------------- !
    integer :: seconds, days, year_in_s
    real    :: r_seconds, r_days, r_total_seconds, frac_of_day, frac_of_year, gmt, &
               time_since_ae, dt_cosp_radians, day_in_s, r_solday, r_dt_cosp_avg, rrsun
    real, dimension(size(temp_in,1),size(temp_in,2)) :: coszen, sunlit_in, p2, fracsun, dp
    real, dimension(size(temp_in,1),size(temp_in,2),size(temp_in,3)) :: ozone_in, co2_in
    type(time_type) :: Time_loc
    integer :: k

    !write (*,*) 'QL, go into run_cosp now'

    ! check if we really want to recompute COSP
    ! alarm
    call get_time(Time,seconds,days)
    r_days = real(days)
    r_seconds = real(seconds)
    r_total_seconds = r_seconds + r_days*86400.

    if(r_total_seconds - dt_last .ge. dt_cosp) then
        dt_last = r_total_seconds
    else
      ! here should add output??
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Output
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      !if (js==1) write(*,*) 'QL else-clause, before output, dt_last', is, ie, js, je, dt_last, r_total_seconds
      !call output_cosp_fields(Time_diag, is, ie, js, je)
      call output_cosp_fields(Time_diag, size(temp_in,1), size(temp_in,2))
      !if (js==1) write(*,*) 'QL  else-clause, after output', dt_last, r_total_seconds

      return !not time yet
    end if

    ! ************* Sort out the time related input variables ************** !
    !make sure we run perpetual when solday > 0)
    if(solday > 0)then
      Time_loc = set_time(seconds,solday)
    else
      Time_loc = Time
    endif

    !Set tide-locked flux if tidally-locked = .true. Else use diurnal-solar
    !to calculate insolation from orbit!
    if (tidally_locked.eq..true.) then
      coszen = COS(rad_lat(:,:))*COS(rad_lon(:,:))
      WHERE (coszen < 0.0) coszen = 0.0
      rrsun = 1 ! needs to be set, set to 1 so that stellar_radiation is unchanged
    elseif (frierson_solar_rad .eq. .true.) then
      p2     = (1. - 3.*sin(rad_lat(:,:))**2)/4.
      coszen = 0.25 * (1.0 + del_sol * p2 + del_sw * sin(rad_lat(:,:)))
      rrsun  = 1 ! needs to be set, set to 1 so that stellar_radiation is unchanged
    else
      ! compute zenith angle
      call get_time(Time_loc, seconds, days)
      call get_time(length_of_year(), year_in_s)
      day_in_s = length_of_day()

      r_seconds=real(seconds)
      r_days=real(days)
      r_total_seconds=r_seconds+(r_days*86400.)
      frac_of_day = r_total_seconds / day_in_s

      if(solday > 0) then
        r_solday=real(solday)
        frac_of_year=(r_solday*day_in_s)/year_in_s
      else
        frac_of_year = r_total_seconds / year_in_s
      endif
      gmt = abs(mod(frac_of_day, 1.0)) * 2.0 * pi
      time_since_ae = modulo(frac_of_year-equinox_day, 1.0) * 2.0 * pi

      if(do_cosp_time_avg) then
        r_dt_cosp_avg = real(dt_cosp_avg)
        dt_cosp_radians = (r_dt_cosp_avg/day_in_s)*2.0*pi
        call diurnal_solar(rad_lat, rad_lon, gmt, time_since_ae, coszen, fracsun, rrsun, dt_cosp_radians)
      else
        ! Seasonal Cycle: Use astronomical parameters to calculate insolation
        call diurnal_solar(rad_lat, rad_lon, gmt, time_since_ae, coszen, fracsun, rrsun)
      end if
    endif

    ! Copied from lines 5737-5741 am3/src/atmos_param/sea_esf_rad/cloudrad_diagnostics.F90
    ! AM3 source code can be obtained from https://www.gfdl.noaa.gov/am3/#acquire
    ! Compute sunlit falg
    sunlit_in = 0.0
    where(coszen > 1.e-06) sunlit_in = 1

    ! Input ozone and co2
    ozone_in = 0.0

    !get ozone
    if(do_read_ozone)then
        call interpolator(o3_interp, Time_diag, p_half_in, ozone_in, trim(ozone_field_name))
        if (input_o3_file_is_mmr==.false.) then
            ozone_in = ozone_in * amO3 / (1.e3 * gas_constant / rdgas )
            ! COSP expects all abundances to be mass mixing ratio. So if input file is volume mixing ratio,
            ! it must be converted to mass mixing ratio using the molar masses of dry air and ozone
            ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol
            ! from kg/mol by multiplying by 1000. This conversion is necessary because wtmozone is in g/mol.
        endif
    endif

    if (input_co2_mmr==.false.) then
        co2_in = co2_ppmv * 1.e-6 * amCO2 / (1.e3 * gas_constant / rdgas )
    else
        co2_in = co2_ppmv * 1.e-6 !No need to convert if it is already a mmr
    endif

    !get co2
    if(do_read_co2) then
      call interpolator(co2_interp, Time_diag, p_half_in, co2_in, trim(co2_field_name))
      if (input_co2_mmr==.false.) then
          co2_in = co2_in * 1.e-6 * amCO2 / (1.e3 * gas_constant / rdgas )
          ! Molar mass of dry air calculated from gas_constant / rdgas, and converted into g/mol
          ! from kg/mol by multiplying by 1000. This conversion is necessary because wtmco2 is in g/mol.
      endif
    endif

    !write(*,*) 'QL: after interpolate co2'

    ! cast into COSP wp format
    rad_lat_cosp = real(rad_lat, kind(wp))
    rad_lon_cosp = real(rad_lon, kind(wp))
    p_full_cosp = real(p_full_in, kind(wp))
    p_half_cosp = real(p_half_in, kind(wp))
    z_full_cosp = real(z_full_in, kind(wp))
    z_half_cosp = real(z_half_in, kind(wp))
    temp_cosp = real(temp_in, kind(wp))
    q_cosp = real(q_in, kind(wp))
    rh_cosp = real(rh_in, kind(wp))
    cf_rad_cosp = real(cf_rad, kind(wp))
    reff_rad_cosp = real(reff_rad, kind(wp))
    qcl_rad_cosp = real(qcl_rad, kind(wp))

    ls_cloud_absorptivity_cosp = real(ls_cloud_absorptivity, kind(wp))
    cnv_cloud_absorptivity_cosp = real(cnv_cloud_absorptivity, kind(wp))
    ls_cloud_extinction_cosp = real(ls_cloud_extinction, kind(wp))
    cnv_cloud_extinction_cosp = real(cnv_cloud_extinction, kind(wp))

    t_surf_cosp = real(t_surf_in, kind(wp))
    landmask_cosp = real(landmask_in, kind(wp))
    u_wind_cosp = real(u_wind_in, kind(wp))
    v_wind_cosp = real(v_wind_in, kind(wp))
    z_surf_cosp = real(z_surf_in, kind(wp))
    sunlit_cosp = real(sunlit_in, kind(wp))
    ozone_cosp = real(ozone_in, kind(wp))
    co2_cosp = real(co2_in, kind(wp))

    ! Cloud related variables
    tot_cld_amt = cf_rad_cosp
    conv_cld_amt = 0.0
    mmr_ls_liq = qcl_rad_cosp / (1.0 - qcl_rad_cosp)
    mmr_ls_ice = 0.0
    mmr_cc_liq = 0.0
    mmr_cc_ice = 0.0
    flux_ls_rain = 0.0 !ls_rain_cosp
    flux_ls_snow = 0.0 !ls_snow_cosp
    flux_ls_graupel = 0.0
    flux_cc_rain = 0.0 !conv_rain_cosp
    flux_cc_snow = 0.0 !conv_snow_cosp

    ! change extinction to optical depth
    ! rho*k*dz = rho*k*dp/(rho*g) = k*dp/g
    ! check the sign, pos or neg?
    do k=Nlevels+1,2,-1
      dp = p_half_cosp(:,:,k) - p_half_cosp(:,:,k-1)
      ! write(*,*) 'QL, k, phalf:', sum(p_half_cosp(:,:,k) )/size(p_half_cosp(:,:,k) ), &
      !                             sum(p_half_cosp(:,:,k-1) )/size(p_half_cosp(:,:,k-1) )
      ! write(*,*) 'QL dp', k, sum(dp)/size(dp), &
      !     sum(ls_cloud_extinction_cosp(:,:,k))/size(ls_cloud_extinction_cosp(:,:,k))
      ls_cloud_tau_cosp(:,:,k) = ls_cloud_extinction_cosp(:,:,k) * dp / grav
      cnv_cloud_tau_cosp(:,:,k) = cnv_cloud_extinction_cosp(:,:,k) * dp / grav
      where ( ls_cloud_tau_cosp(:,:,k) < 0)  ls_cloud_tau_cosp(:,:,k) = 0.0
      where (cnv_cloud_tau_cosp(:,:,k) < 0) cnv_cloud_tau_cosp(:,:,k) = 0.0
      ! write(*,*) 'QL ls tau', k, sum(ls_cloud_tau_cosp(:,:,k))/size(ls_cloud_tau_cosp(:,:,k))
      ! write(*,*) 'QL cnv tau', k, sum(cnv_cloud_tau_cosp(:,:,k))/size(cnv_cloud_tau_cosp(:,:,k))
    enddo

    dtau_strat = ls_cloud_tau_cosp
    dtau_conv = cnv_cloud_tau_cosp
    dlw_emissivity_strat = ls_cloud_absorptivity_cosp
    dlw_emissivity_conv = cnv_cloud_absorptivity_cosp

    !do i = 1,N_HYDRO
    !  hydrometeor_Reff(:,:,:,i) = reff_rad_cosp
    !end do
    !hydrometeor_Reff(:,:,:,1) = reff_rad_cosp
    hydrometeor_Reff = 0.0

    sfc_lw_emissivity = 1.0

    !write(*,*) 'QL: before call cosp interface'
    call fms_cosp_interface(rad_lat_cosp, rad_lon_cosp, &
              p_full_cosp, p_half_cosp, z_full_cosp, z_half_cosp, temp_cosp, &
              q_cosp, rh_cosp, tot_cld_amt, conv_cld_amt, &
              mmr_ls_liq, mmr_ls_ice, mmr_cc_liq, mmr_cc_ice, &
              flux_ls_rain, flux_ls_snow, flux_ls_graupel, &
              flux_cc_rain, flux_cc_snow, &
              dtau_strat, dtau_conv, dlw_emissivity_strat, &
              dlw_emissivity_conv, ozone_cosp, co2_cosp, hydrometeor_Reff, &
              t_surf_cosp, landmask_cosp, u_wind_cosp, v_wind_cosp, sunlit_cosp, &
              z_surf_cosp, sfc_lw_emissivity)
    !write(*,*) 'QL: after call cosp interface'

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Output
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !if (js==1) write(*,*) 'QL before output, dt_last', is, ie, js, je, dt_last, r_total_seconds
    !call output_cosp_fields(Time_diag, is, ie, js, je)
    call output_cosp_fields(Time_diag, size(temp_in,1), size(temp_in,2))
    !if (js==1) write(*,*) 'QL after output', dt_last, r_total_seconds

  end subroutine fms_run_cosp


  subroutine fms_cosp_end()
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Free up memory
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call destroy_cosp_outputs(cospOUT)
    call destroy_cospIN(cospIN)
    call destroy_cospstateIN(cospstateIN)
    call cosp_cleanUp()

    if(do_read_ozone) call interpolator_end(o3_interp)
    if(do_read_co2)   call interpolator_end(co2_interp)
  end subroutine fms_cosp_end


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Which simulators need to be run? Look at which outputs are requested.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine select_simulator()
    if (Lpctisccp .or. Lclisccp .or. Lboxptopisccp .or.  Lboxtauisccp .or. Ltauisccp .or. &
      Lcltisccp .or. Lmeantbisccp .or. Lmeantbclrisccp .or. Lalbisccp) Lisccp = .true.

    if (LclMISR) Lmisr = .true.

    if (Lcltmodis .or. Lclwmodis .or. Lclimodis .or. Lclhmodis .or. Lclmmodis .or.        &
      Lcllmodis .or. Ltautmodis .or. Ltauwmodis .or. Ltauimodis .or. Ltautlogmodis .or. &
      Ltauwlogmodis .or. Ltauilogmodis .or. Lreffclwmodis .or. Lreffclimodis .or.       &
      Lpctmodis .or. Llwpmodis .or. Liwpmodis .or. Lclmodis) Lmodis = .true.

    if (Lclcalipso2 .or. Lclcalipso .or.  Lclhcalipso .or. Lcllcalipso .or. Lclmcalipso   &
      .or. Lcltcalipso .or. Lcltlidarradar .or. Lclcalipsoliq .or. Lclcalipsoice .or.   &
      Lclcalipsoun .or. Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsotmpice .or.  &
      Lclcalipsotmpun .or. Lcltcalipsoliq .or. Lcltcalipsoice .or. Lcltcalipsoun .or.   &
      Lclhcalipsoliq .or. Lclhcalipsoice .or. Lclhcalipsoun .or. Lclmcalipsoliq .or.    &
      Lclmcalipsoice .or. Lclmcalipsoun .or. Lcllcalipsoliq .or. Lcllcalipsoice .or.    &
      Lcllcalipsoun .or. LlidarBetaMol532 .or. LcfadLidarsr532 .or. Lcltlidarradar .or. &
      Lcltlidarradar .or. Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso   &
      .or. Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or.             &
      Lclcalipsoopacity .or. Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp .or.    &
      Lclopaquemeanz .or. Lclthinmeanz .or. Lclthinemis .or. Lclopaquemeanzse .or.      &
      Lclthinmeanzse .or. Lclzopaquecalipsose) Lcalipso = .true.

    if (LlidarBetaMol532gr .or. LcfadLidarsr532gr .or. Latb532gr .or. LclgrLidar532 .or.  &
      LclhgrLidar532 .or. LcllgrLidar532 .or. LclmgrLidar532 .or. LcltgrLidar532)   &
      LgrLidar532 = .true.

    if (LlidarBetaMol355 .or. LcfadLidarsr355 .or. Latb355 .or. Lclatlid .or.           &
      Lclhatlid .or. Lcllatlid .or. Lclmatlid .or. Lcltatlid)                           &
      Latlid = .true.

    if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar) Lcloudsat = .true.

    if (LcfadDbze94 .or. Ldbze94 .or. Lcltlidarradar .or. Lptradarflag0 .or. Lptradarflag1 &
      .or. Lptradarflag2 .or. Lptradarflag3 .or. Lptradarflag4 .or. Lptradarflag5 .or.  &
      Lptradarflag6 .or. Lptradarflag7 .or. Lptradarflag8 .or. Lptradarflag9 .or.       &
      Lradarpia) Lcloudsat = .true.

    if (Lparasolrefl) Lparasol = .true.

    if (Ltbrttov) Lrttov = .true.

  end subroutine select_simulator


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Register output variables according to cospOUT flags
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine diag_field_init(Time, axes)

    type(time_type), intent(in) :: Time
    integer, dimension(4), intent(in) :: axes

    ! ISCCP simulator outputs
    if (associated(cospOUT%isccp_totalcldarea)) then
      id_cltisccp = register_diag_field &
        (mod_name, 'cltisccp', axes(1:2), Time, &
        'ISCCP Total Cloud Fraction', &
        '%', missing_value=missing_value) ! mask_variant=.true., 
    endif

    ! CALIPSO simulator outputs
    if (associated(cospOUT%calipso_cldlayer)) then
      ! Low-level
      id_cllcalipso = register_diag_field &
        (mod_name, 'cllcalipso', axes(1:2), Time, &
        'CALIPSO Low Level Cloud Fraction', &
        '%', missing_value=missing_value)
      ! Mid-level
      id_clmcalipso = register_diag_field &
        (mod_name, 'clmcalipso', axes(1:2), Time, &
        'CALIPSO Mid Level Cloud Fraction', &
        '%', missing_value=missing_value)
      ! High-level
      id_clhcalipso = register_diag_field &
        (mod_name, 'clhcalipso', axes(1:2), Time, &
        'CALIPSO High Level Cloud Fraction', &
        '%', missing_value=missing_value)
      ! Total
      id_cltcalipso = register_diag_field &
        (mod_name, 'cltcalipso', axes(1:2), Time, &
        'CALIPSO Total Cloud Fraction', &
        '%', missing_value=missing_value)
    endif

  end subroutine diag_field_init


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Output the COSP variables
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine output_cosp_fields(Time_diag, Nlon, Nlat)
    type(time_type), intent(in) :: Time_diag
    integer, intent(in) :: Nlat, Nlon
    logical :: used
    real, dimension(Nlon,Nlat) :: y2save

    ! Cast data from type real(wp) to real.

    if (id_cltisccp > 0) then
      call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%isccp_totalcldarea), y2=y2save)
      used = send_data(id_cltisccp, y2save, Time_diag) !, mask=y2save/=missing_value)
    endif

    if (id_cllcalipso > 0) then
      call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,1)), y2=y2save)
      used = send_data(id_cllcalipso, y2save, Time_diag)
    endif
    
    if (id_clmcalipso > 0) then
      call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,2)), y2=y2save)
      used = send_data(id_clmcalipso, y2save, Time_diag)
    endif

    if (id_clhcalipso > 0) then
      call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,3)), y2=y2save)
      used = send_data(id_clhcalipso, y2save, Time_diag)
    endif
    
    if (id_cltcalipso > 0) then
      call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,4)), y2=y2save)
      used = send_data(id_cltcalipso, y2save, Time_diag)
    endif
  
  end subroutine output_cosp_fields


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE subsample_and_optics
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine subsample_and_optics(nPoints, nLevels, nColumns, nHydro, overlap,           &
    use_precipitation_fluxes, lidar_ice_type, sd, tca, cca, fl_lsrainIN, fl_lssnowIN,    &
    fl_lsgrplIN, fl_ccrainIN, fl_ccsnowIN, mr_lsliq, mr_lsice, mr_ccliq, mr_ccice,       &
    reffIN, dtau_c, dtau_s, dem_c, dem_s, cospstateIN, cospIN)
    ! Inputs
    integer,intent(in) :: nPoints, nLevels, nColumns, nHydro, overlap, lidar_ice_type
    real(wp),intent(in),dimension(nPoints,nLevels) :: tca, cca, mr_lsliq, mr_lsice, mr_ccliq, &
          mr_ccice, dtau_c, dtau_s, dem_c, dem_s, fl_lsrainIN, fl_lssnowIN, fl_lsgrplIN,      &
          fl_ccrainIN, fl_ccsnowIN
    real(wp),intent(in),dimension(nPoints,nLevels,nHydro) :: reffIN
    logical,intent(in) :: use_precipitation_fluxes
    type(size_distribution),intent(inout) :: sd

    ! Outputs
    type(cosp_optical_inputs),intent(inout) :: cospIN
    type(cosp_column_inputs),intent(inout)  :: cospstateIN

    ! Local variables
    type(rng_state),allocatable,dimension(:) :: rngs  ! Seeds for random number generator
    integer,dimension(:),allocatable :: seed
    integer,dimension(:),allocatable :: cloudsat_preclvl_index
    integer :: i,j,k
    real(wp) :: zstep
    real(wp),dimension(:,:), allocatable :: &
          ls_p_rate, cv_p_rate, frac_ls, frac_cv, prec_ls, prec_cv,g_vol
    real(wp),dimension(:,:,:),  allocatable :: &
          frac_prec, MODIS_cloudWater, MODIS_cloudIce, fracPrecipIce, fracPrecipIce_statGrid,&
          MODIS_watersize,MODIS_iceSize, MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce
    real(wp),dimension(:,:,:,:),allocatable :: &
          mr_hydro, Reff, Np
    real(wp),dimension(nPoints,nLevels) :: &
          column_frac_out, column_prec_out, fl_lsrain, fl_lssnow, fl_lsgrpl, fl_ccrain, fl_ccsnow
    real(wp),dimension(nPoints,nColumns,Nlvgrid_local) :: tempOut
    logical :: cmpGases=.true.

    if (Ncolumns .gt. 1) then
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Generate subcolumns for clouds (SCOPS) and precipitation type (PREC_SCOPS)
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! RNG used for subcolumn generation
      allocate(rngs(nPoints),seed(nPoints))
      seed(:)=0
      seed = int(cospstateIN%phalf(:,Nlevels+1))  ! In case of NPoints=1
      ! *NOTE* Chunking will change the seed
      if (NPoints .gt. 1) seed=int((cospstateIN%phalf(:,Nlevels+1)-minval(cospstateIN%phalf(:,Nlevels+1)))/      &
          (maxval(cospstateIN%phalf(:,Nlevels+1))-minval(cospstateIN%phalf(:,Nlevels+1)))*100000) + 1
      call init_rng(rngs, seed)

      ! Call scops
      call scops(NPoints,Nlevels,Ncolumns,rngs,tca,cca,overlap,cospIN%frac_out,0)
      deallocate(seed,rngs)

      ! Sum up precipitation rates
      allocate(ls_p_rate(nPoints,nLevels),cv_p_rate(nPoints,Nlevels))
      if(use_precipitation_fluxes) then
        ls_p_rate(:,1:nLevels) = fl_lsrainIN + fl_lssnowIN + fl_lsgrplIN
        cv_p_rate(:,1:nLevels) = fl_ccrainIN + fl_ccsnowIN
      else
        ls_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow) + mixing_ratio (groupel)
        cv_p_rate(:,1:nLevels) = 0 ! mixing_ratio(rain) + mixing_ratio(snow)
      endif

      ! Call PREC_SCOPS
      allocate(frac_prec(nPoints,nColumns,nLevels))
      call prec_scops(nPoints,nLevels,nColumns,ls_p_rate,cv_p_rate,cospIN%frac_out,frac_prec)
      deallocate(ls_p_rate,cv_p_rate)

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Compute fraction in each gridbox for precipitation  and cloud type.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Allocate
      allocate(frac_ls(nPoints,nLevels),prec_ls(nPoints,nLevels),                       &
              frac_cv(nPoints,nLevels),prec_cv(nPoints,nLevels))

      ! Initialize
      frac_ls(1:nPoints,1:nLevels) = 0._wp
      prec_ls(1:nPoints,1:nLevels) = 0._wp
      frac_cv(1:nPoints,1:nLevels) = 0._wp
      prec_cv(1:nPoints,1:nLevels) = 0._wp
      do j=1,nPoints
        do k=1,nLevels
            do i=1,nColumns
              if (cospIN%frac_out(j,i,k)  .eq. 1)  frac_ls(j,k) = frac_ls(j,k)+1._wp
              if (cospIN%frac_out(j,i,k)  .eq. 2)  frac_cv(j,k) = frac_cv(j,k)+1._wp
              if (frac_prec(j,i,k) .eq. 1)  prec_ls(j,k) = prec_ls(j,k)+1._wp
              if (frac_prec(j,i,k) .eq. 2)  prec_cv(j,k) = prec_cv(j,k)+1._wp
              if (frac_prec(j,i,k) .eq. 3)  prec_cv(j,k) = prec_cv(j,k)+1._wp
              if (frac_prec(j,i,k) .eq. 3)  prec_ls(j,k) = prec_ls(j,k)+1._wp
            enddo
            frac_ls(j,k)=frac_ls(j,k)/nColumns
            frac_cv(j,k)=frac_cv(j,k)/nColumns
            prec_ls(j,k)=prec_ls(j,k)/nColumns
            prec_cv(j,k)=prec_cv(j,k)/nColumns
        enddo
      enddo

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Assign gridmean mixing-ratios (mr_XXXXX), effective radius (ReffIN) and number
      ! concentration (not defined) to appropriate sub-column. Here we are using scops.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      allocate(mr_hydro(nPoints,nColumns,nLevels,nHydro),                               &
              Reff(nPoints,nColumns,nLevels,nHydro),                                   &
              Np(nPoints,nColumns,nLevels,nHydro))

      ! Initialize
      mr_hydro(:,:,:,:) = 0._wp
      Reff(:,:,:,:)     = 0._wp
      Np(:,:,:,:)       = 0._wp
      do k=1,nColumns
        ! Subcolumn cloud fraction
        column_frac_out = cospIN%frac_out(:,k,:)

        ! LS clouds
        where (column_frac_out == I_LSC)
            mr_hydro(:,k,:,I_LSCLIQ) = mr_lsliq
            mr_hydro(:,k,:,I_LSCICE) = mr_lsice
            Reff(:,k,:,I_LSCLIQ)     = ReffIN(:,:,I_LSCLIQ)
            Reff(:,k,:,I_LSCICE)     = ReffIN(:,:,I_LSCICE)
        ! CONV clouds
        elsewhere (column_frac_out == I_CVC)
            mr_hydro(:,k,:,I_CVCLIQ) = mr_ccliq
            mr_hydro(:,k,:,I_CVCICE) = mr_ccice
            Reff(:,k,:,I_CVCLIQ)     = ReffIN(:,:,I_CVCLIQ)
            Reff(:,k,:,I_CVCICE)     = ReffIN(:,:,I_CVCICE)
        end where

        ! Subcolumn precipitation
        column_prec_out = frac_prec(:,k,:)

        ! LS Precipitation
        where ((column_prec_out == 1) .or. (column_prec_out == 3) )
            Reff(:,k,:,I_LSRAIN) = ReffIN(:,:,I_LSRAIN)
            Reff(:,k,:,I_LSSNOW) = ReffIN(:,:,I_LSSNOW)
            Reff(:,k,:,I_LSGRPL) = ReffIN(:,:,I_LSGRPL)
            ! CONV precipitation
        elsewhere ((column_prec_out == 2) .or. (column_prec_out == 3))
            Reff(:,k,:,I_CVRAIN) = ReffIN(:,:,I_CVRAIN)
            Reff(:,k,:,I_CVSNOW) = ReffIN(:,:,I_CVSNOW)
        end where
      enddo

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Convert the subcolumn mixing ratio and precipitation fluxes from gridbox mean
      ! values to fraction-based values.
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Initialize
      fl_lsrain(:,:) = 0._wp
      fl_lssnow(:,:) = 0._wp
      fl_lsgrpl(:,:) = 0._wp
      fl_ccrain(:,:) = 0._wp
      fl_ccsnow(:,:) = 0._wp
      do k=1,nLevels
        do j=1,nPoints
            ! In-cloud mixing ratios.
            if (frac_ls(j,k) .ne. 0.) then
              mr_hydro(j,:,k,I_LSCLIQ) = mr_hydro(j,:,k,I_LSCLIQ)/frac_ls(j,k)
              mr_hydro(j,:,k,I_LSCICE) = mr_hydro(j,:,k,I_LSCICE)/frac_ls(j,k)
            endif
            if (frac_cv(j,k) .ne. 0.) then
              mr_hydro(j,:,k,I_CVCLIQ) = mr_hydro(j,:,k,I_CVCLIQ)/frac_cv(j,k)
              mr_hydro(j,:,k,I_CVCICE) = mr_hydro(j,:,k,I_CVCICE)/frac_cv(j,k)
            endif
            ! Precipitation
            if (use_precipitation_fluxes) then
              if (prec_ls(j,k) .ne. 0.) then
                  fl_lsrain(j,k) = fl_lsrainIN(j,k)/prec_ls(j,k)
                  fl_lssnow(j,k) = fl_lssnowIN(j,k)/prec_ls(j,k)
                  fl_lsgrpl(j,k) = fl_lsgrplIN(j,k)/prec_ls(j,k)
              endif
              if (prec_cv(j,k) .ne. 0.) then
                  fl_ccrain(j,k) = fl_ccrainIN(j,k)/prec_cv(j,k)
                  fl_ccsnow(j,k) = fl_ccsnowIN(j,k)/prec_cv(j,k)
              endif
            else
              if (prec_ls(j,k) .ne. 0.) then
                  mr_hydro(j,:,k,I_LSRAIN) = mr_hydro(j,:,k,I_LSRAIN)/prec_ls(j,k)
                  mr_hydro(j,:,k,I_LSSNOW) = mr_hydro(j,:,k,I_LSSNOW)/prec_ls(j,k)
                  mr_hydro(j,:,k,I_LSGRPL) = mr_hydro(j,:,k,I_LSGRPL)/prec_ls(j,k)
              endif
              if (prec_cv(j,k) .ne. 0.) then
                  mr_hydro(j,:,k,I_CVRAIN) = mr_hydro(j,:,k,I_CVRAIN)/prec_cv(j,k)
                  mr_hydro(j,:,k,I_CVSNOW) = mr_hydro(j,:,k,I_CVSNOW)/prec_cv(j,k)
              endif
            endif
        enddo
      enddo
      deallocate(frac_ls,prec_ls,frac_cv,prec_cv)

      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      ! Convert precipitation fluxes to mixing ratios
      !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if (use_precipitation_fluxes) then
        ! LS rain
        call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSRAIN), n_bx(I_LSRAIN),         &
              alpha_x(I_LSRAIN), c_x(I_LSRAIN),   d_x(I_LSRAIN),   g_x(I_LSRAIN),       &
              a_x(I_LSRAIN),   b_x(I_LSRAIN),   gamma_1(I_LSRAIN), gamma_2(I_LSRAIN),   &
              gamma_3(I_LSRAIN), gamma_4(I_LSRAIN), fl_lsrain,                          &
              mr_hydro(:,:,:,I_LSRAIN), Reff(:,:,:,I_LSRAIN))
        ! LS snow
        call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 1._wp,  n_ax(I_LSSNOW),  n_bx(I_LSSNOW),       &
              alpha_x(I_LSSNOW), c_x(I_LSSNOW),  d_x(I_LSSNOW),  g_x(I_LSSNOW),         &
              a_x(I_LSSNOW),   b_x(I_LSSNOW),   gamma_1(I_LSSNOW),  gamma_2(I_LSSNOW),  &
              gamma_3(I_LSSNOW), gamma_4(I_LSSNOW), fl_lssnow,                          &
              mr_hydro(:,:,:,I_LSSNOW), Reff(:,:,:,I_LSSNOW))
        ! CV rain
        call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVRAIN),  n_bx(I_CVRAIN),        &
              alpha_x(I_CVRAIN), c_x(I_CVRAIN),   d_x(I_CVRAIN),   g_x(I_CVRAIN),       &
              a_x(I_CVRAIN),   b_x(I_CVRAIN),   gamma_1(I_CVRAIN), gamma_2(I_CVRAIN),   &
              gamma_3(I_CVRAIN), gamma_4(I_CVRAIN), fl_ccrain,                          &
              mr_hydro(:,:,:,I_CVRAIN), Reff(:,:,:,I_CVRAIN))
        ! CV snow
        call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 2._wp, n_ax(I_CVSNOW),  n_bx(I_CVSNOW),        &
              alpha_x(I_CVSNOW),  c_x(I_CVSNOW),   d_x(I_CVSNOW),   g_x(I_CVSNOW),      &
              a_x(I_CVSNOW),   b_x(I_CVSNOW),   gamma_1(I_CVSNOW), gamma_2(I_CVSNOW),   &
              gamma_3(I_CVSNOW), gamma_4(I_CVSNOW), fl_ccsnow,                          &
              mr_hydro(:,:,:,I_CVSNOW), Reff(:,:,:,I_CVSNOW))
        ! LS groupel.
        call cosp_precip_mxratio(nPoints, nLevels, nColumns, cospstateIN%pfull,        &
              cospstateIN%at, frac_prec, 1._wp, n_ax(I_LSGRPL),  n_bx(I_LSGRPL),        &
              alpha_x(I_LSGRPL), c_x(I_LSGRPL),   d_x(I_LSGRPL),   g_x(I_LSGRPL),       &
              a_x(I_LSGRPL),   b_x(I_LSGRPL),   gamma_1(I_LSGRPL),  gamma_2(I_LSGRPL),  &
              gamma_3(I_LSGRPL), gamma_4(I_LSGRPL), fl_lsgrpl,                          &
              mr_hydro(:,:,:,I_LSGRPL), Reff(:,:,:,I_LSGRPL))
        deallocate(frac_prec)
      endif

    else
      cospIN%frac_out(:,:,:) = 1
      allocate(mr_hydro(nPoints,1,nLevels,nHydro),Reff(nPoints,1,nLevels,nHydro),       &
              Np(nPoints,1,nLevels,nHydro))
      mr_hydro(:,1,:,I_LSCLIQ) = mr_lsliq
      mr_hydro(:,1,:,I_LSCICE) = mr_lsice
      mr_hydro(:,1,:,I_CVCLIQ) = mr_ccliq
      mr_hydro(:,1,:,I_CVCICE) = mr_ccice
      Reff(:,1,:,:)            = ReffIN
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 11 micron emissivity
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp) then
      call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dem_c,dem_s,    &
                                cospIN%emiss_11)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 0.67 micron optical depth
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp .or. Lmisr .or. Lmodis) then
      call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,dtau_c,dtau_s,  &
                                cospIN%tau_067)
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! LIDAR Polarized optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !write(*,*) 'QL before call  any(cospIN%tautot_calipso .lt. 0)',  any(cospIN%tautot_calipso .lt. 0)
    !write(*,*) 'QL, Lcalipso=', Lcalipso
    if (Lcalipso) then
      call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 532, .false.,     &
          mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
          mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
          ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
          cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_calipso,                    &
          cospIN%betatot_calipso, cospIN%tau_mol_calipso, cospIN%tautot_calipso,         &
          cospIN%tautot_S_liq, cospIN%tautot_S_ice, cospIN%betatot_ice_calipso,          &
          cospIN%betatot_liq_calipso, cospIN%tautot_ice_calipso, cospIN%tautot_liq_calipso)
      !write(*,*) 'QL after call  any(cospIN%tautot_calipso .lt. 0)', any(cospIN%tautot_calipso .lt. 0)
    endif

    if (LgrLidar532) then
      call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 532, .true.,      &
          mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
          mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
          ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
          cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_grLidar532,                 &
          cospIN%betatot_grLidar532, cospIN%tau_mol_grLidar532, cospIN%tautot_grLidar532)
    endif

    if (Latlid) then
      call lidar_optics(nPoints, nColumns, nLevels, 4, lidar_ice_type, 355, .false.,     &
          mr_hydro(:,:,:,I_LSCLIQ),  mr_hydro(:,:,:,I_LSCICE), mr_hydro(:,:,:,I_CVCLIQ), &
          mr_hydro(:,:,:,I_CVCICE), ReffIN(:,:,I_LSCLIQ), ReffIN(:,:,I_LSCICE),          &
          ReffIN(:,:,I_CVCLIQ), ReffIN(:,:,I_CVCICE), cospstateIN%pfull,                 &
          cospstateIN%phalf, cospstateIN%at, cospIN%beta_mol_atlid, cospIN%betatot_atlid,&
          cospIN%tau_mol_atlid, cospIN%tautot_atlid)
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! CLOUDSAT RADAR OPTICS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (lcloudsat) then

      ! Compute gaseous absorption (assume identical for each subcolun)
      allocate(g_vol(nPoints,nLevels))
      g_vol(:,:)=0._wp
      do i=1,nPoints
        do j=1,nLevels
            if (rcfg_cloudsat%use_gas_abs == 1 .or. (rcfg_cloudsat%use_gas_abs == 2 .and. j .eq. 1)) then
              g_vol(i,j) = gases(cospstateIN%pfull(i,j), cospstateIN%at(i,j),cospstateIN%qv(i,j),rcfg_cloudsat%freq)
            endif
            cospIN%g_vol_cloudsat(i,:,j)=g_vol(i,j)
        end do
      end do

      ! Loop over all subcolumns
      allocate(fracPrecipIce(nPoints,nColumns,nLevels))
      fracPrecipIce(:,:,:) = 0._wp
      do k=1,nColumns
        call quickbeam_optics(sd, rcfg_cloudsat, nPoints, nLevels, R_UNDEF,   &
              mr_hydro(:,k,:,1:nHydro)*1000._wp, Reff(:,k,:,1:nHydro)*1.e6_wp,&
              Np(:,k,:,1:nHydro), cospstateIN%pfull, cospstateIN%at,          &
              cospstateIN%qv, cospIN%z_vol_cloudsat(1:nPoints,k,:),           &
              cospIN%kr_vol_cloudsat(1:nPoints,k,:))

        ! At each model level, what fraction of the precipitation is frozen?
        where(mr_hydro(:,k,:,I_LSRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_LSSNOW) .gt. 0 .or. &
              mr_hydro(:,k,:,I_CVRAIN) .gt. 0 .or. mr_hydro(:,k,:,I_CVSNOW) .gt. 0 .or. &
              mr_hydro(:,k,:,I_LSGRPL) .gt. 0)
            fracPrecipIce(:,k,:) = (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + &
                mr_hydro(:,k,:,I_LSGRPL)) / &
                (mr_hydro(:,k,:,I_LSSNOW) + mr_hydro(:,k,:,I_CVSNOW) + mr_hydro(:,k,:,I_LSGRPL) + &
                mr_hydro(:,k,:,I_LSRAIN)  + mr_hydro(:,k,:,I_CVRAIN))
        elsewhere
            fracPrecipIce(:,k,:) = 0._wp
        endwhere
      enddo

      ! Regrid frozen fraction to Cloudsat/Calipso statistical grid
      allocate(fracPrecipIce_statGrid(nPoints,nColumns,Nlvgrid_local))
      fracPrecipIce_statGrid(:,:,:) = 0._wp
      call cosp_change_vertical_grid(Npoints, Ncolumns, Nlevels, cospstateIN%hgt_matrix(:,Nlevels:1:-1), &
          cospstateIN%hgt_matrix_half(:,Nlevels:1:-1), fracPrecipIce(:,:,Nlevels:1:-1), Nlvgrid_local,  &
          vgrid_zl(Nlvgrid_local:1:-1), vgrid_zu(Nlvgrid_local:1:-1), fracPrecipIce_statGrid(:,:,Nlvgrid_local:1:-1))

      ! Find proper layer above de surface elevation to compute precip flags in Cloudsat/Calipso statistical grid
      allocate(cloudsat_preclvl_index(nPoints))
      cloudsat_preclvl_index(:) = 0._wp
      ! Compute the zstep distance between two atmopsheric layers
      zstep = vgrid_zl(1)-vgrid_zl(2)
      ! Computing altitude index for precip flags calculation (one layer above surfelev layer)
      cloudsat_preclvl_index(:) = cloudsat_preclvl - floor( cospstateIN%surfelev(:)/zstep )

      ! For near-surface diagnostics, we only need the frozen fraction at one layer.
      do i=1,nPoints
        cospIN%fracPrecipIce(i,:) = fracPrecipIce_statGrid(i,:,cloudsat_preclvl_index(i))
      enddo

    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! MODIS optics
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lmodis) then
      allocate(MODIS_cloudWater(nPoints,nColumns,nLevels),                               &
              MODIS_cloudIce(nPoints,nColumns,nLevels),                                  &
              MODIS_waterSize(nPoints,nColumns,nLevels),                                 &
              MODIS_iceSize(nPoints,nColumns,nLevels),                                   &
              MODIS_opticalThicknessLiq(nPoints,nColumns,nLevels),                       &
              MODIS_opticalThicknessIce(nPoints,nColumns,nLevels))
      ! Cloud water
      call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,               &
          mr_hydro(:,:,:,I_CVCLIQ),mr_hydro(:,:,:,I_LSCLIQ),MODIS_cloudWater)
      ! Cloud ice
      call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,               &
          mr_hydro(:,:,:,I_CVCICE),mr_hydro(:,:,:,I_LSCICE),MODIS_cloudIce)
      ! Water droplet size
      call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,               &
          Reff(:,:,:,I_CVCLIQ),Reff(:,:,:,I_LSCLIQ),MODIS_waterSize)
      ! Ice crystal size
      call cosp_simulator_optics(nPoints,nColumns,nLevels,cospIN%frac_out,               &
          Reff(:,:,:,I_CVCICE),Reff(:,:,:,I_LSCICE),MODIS_iceSize)

      ! Partition optical thickness into liquid and ice parts
      call modis_optics_partition(nPoints, nLevels, nColumns, MODIS_cloudWater,          &
          MODIS_cloudIce, MODIS_waterSize, MODIS_iceSize, cospIN%tau_067,                &
          MODIS_opticalThicknessLiq, MODIS_opticalThicknessIce)

      ! Compute assymetry parameter and single scattering albedo
      call modis_optics(nPoints, nLevels, nColumns, MODIS_opticalThicknessLiq,           &
          MODIS_waterSize*1.0e6_wp, MODIS_opticalThicknessIce,                           &
          MODIS_iceSize*1.0e6_wp, cospIN%fracLiq, cospIN%asym, cospIN%ss_alb)

      ! Deallocate memory
      deallocate(MODIS_cloudWater,MODIS_cloudIce,MODIS_WaterSize,MODIS_iceSize,          &
          MODIS_opticalThicknessLiq,MODIS_opticalThicknessIce,mr_hydro,                  &
          Np,Reff)
    endif
  end subroutine subsample_and_optics


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospIN(npoints,ncolumns,nlevels,y)
    ! Inputs
    integer,intent(in) :: &
          npoints,  & ! Number of horizontal gridpoints
          ncolumns, & ! Number of subcolumns
          nlevels     ! Number of vertical levels
    ! Outputs
    type(cosp_optical_inputs),intent(out) :: y

    ! Dimensions
    y%Npoints  = Npoints
    y%Ncolumns = Ncolumns
    y%Nlevels  = Nlevels
    y%Npart    = 4
    y%Nrefl    = PARASOL_NREFL
    allocate(y%frac_out(npoints, ncolumns, nlevels))

    if (Lmodis .or. Lmisr .or. Lisccp) then
      allocate(y%tau_067(npoints, ncolumns, nlevels), &
              y%emiss_11(npoints, ncolumns, nlevels))
    endif
    if (Lcalipso) then
      allocate(y%betatot_calipso(npoints,    ncolumns, nlevels),&
              y%betatot_ice_calipso(npoints, ncolumns, nlevels),&
              y%betatot_liq_calipso(npoints, ncolumns, nlevels),&
              y%tautot_calipso(npoints,      ncolumns, nlevels),&
              y%tautot_ice_calipso(npoints,  ncolumns, nlevels),&
              y%tautot_liq_calipso(npoints,  ncolumns, nlevels),&
              y%beta_mol_calipso(npoints,              nlevels),&
              y%tau_mol_calipso(npoints,               nlevels),&
              y%tautot_S_ice(npoints,        ncolumns),         &
              y%tautot_S_liq(npoints,        ncolumns))
    endif

    if (LgrLidar532) then
      allocate(y%beta_mol_grLidar532(npoints,          nlevels),&
              y%betatot_grLidar532(npoints,  ncolumns, nlevels),&
              y%tau_mol_grLidar532(npoints,            nlevels),&
              y%tautot_grLidar532(npoints,   ncolumns, nlevels))
    endif

    if (Latlid) then
      allocate(y%beta_mol_atlid(npoints,           nlevels),&
              y%betatot_atlid(npoints,   ncolumns, nlevels),&
              y%tau_mol_atlid(npoints,             nlevels),&
              y%tautot_atlid(npoints,    ncolumns, nlevels))
    endif

    if (Lcloudsat) then
      allocate(y%z_vol_cloudsat(npoints, ncolumns, nlevels),&
              y%kr_vol_cloudsat(npoints, ncolumns, nlevels),&
              y%g_vol_cloudsat(npoints,  ncolumns, nlevels),&
              y%fracPrecipIce(npoints,   ncolumns))
    endif
    if (Lmodis) then
      allocate(y%fracLiq(npoints,        ncolumns, nlevels),&
              y%asym(npoints,            ncolumns, nlevels),&
              y%ss_alb(npoints,          ncolumns, nlevels))
    endif

  end subroutine construct_cospIN


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE construct_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine construct_cospstateIN(npoints,nlevels,nchan,y)
    ! Inputs
    integer,intent(in) :: &
          npoints, & ! Number of horizontal gridpoints
          nlevels, & ! Number of vertical levels
          nchan      ! Number of channels
    ! Outputs
    type(cosp_column_inputs),intent(out) :: y

    allocate(y%sunlit(npoints),y%skt(npoints),y%land(npoints),y%at(npoints,nlevels),      &
              y%pfull(npoints,nlevels),y%phalf(npoints,nlevels+1),y%qv(npoints,nlevels),  &
              y%o3(npoints,nlevels),y%hgt_matrix(npoints,nlevels),y%u_sfc(npoints),       &
              y%v_sfc(npoints),y%lat(npoints),y%lon(nPoints),y%emis_sfc(nchan),           &
              y%cloudIce(nPoints,nLevels),y%cloudLiq(nPoints,nLevels),y%surfelev(npoints),&
              y%fl_snow(nPoints,nLevels),y%fl_rain(nPoints,nLevels),y%seaice(npoints),    &
              y%tca(nPoints,nLevels),y%hgt_matrix_half(npoints,nlevels+1))

  end subroutine construct_cospstateIN


  ! ######################################################################################
  ! SUBROUTINE construct_cosp_outputs
  !
  ! This subroutine allocates output fields based on input logical flag switches.
  ! ######################################################################################
  subroutine construct_cosp_outputs(Lpctisccp,Lclisccp,&
                                  Lboxptopisccp,Lboxtauisccp,Ltauisccp,Lcltisccp,      &
                                  Lmeantbisccp,Lmeantbclrisccp,Lalbisccp,LclMISR,      &
                                  Lcltmodis,Lclwmodis,Lclimodis,Lclhmodis,Lclmmodis,   &
                                  Lcllmodis,Ltautmodis,Ltauwmodis,Ltauimodis,          &
                                  Ltautlogmodis,Ltauwlogmodis,Ltauilogmodis,           &
                                  Lreffclwmodis,Lreffclimodis,Lpctmodis,Llwpmodis,     &
                                  Liwpmodis,Lclmodis,Latb532,Latb532gr,Latb355,        &
                                  LlidarBetaMol532,LlidarBetaMol532gr,LlidarBetaMol355,&
                                  LcfadLidarsr532,LcfadLidarsr532gr,LcfadLidarsr355,   &
                                  Lclcalipso2,Lclcalipso,LclgrLidar532,Lclatlid,      &
                                  Lclhcalipso,Lcllcalipso,Lclmcalipso,Lcltcalipso,     &
                                  LclhgrLidar532,LcllgrLidar532,LclmgrLidar532,     &
                                  LcltgrLidar532,Lclhatlid,Lcllatlid,Lclmatlid,       &
                                  Lcltatlid,Lcltlidarradar,Lcloudsat_tcc,            &
                                  Lcloudsat_tcc2,Lclcalipsoliq,              &
                                  Lclcalipsoice,Lclcalipsoun,Lclcalipsotmp,            &
                                  Lclcalipsotmpliq,Lclcalipsotmpice,Lclcalipsotmpun,   &
                                  Lcltcalipsoliq,Lcltcalipsoice,Lcltcalipsoun,         &
                                  Lclhcalipsoliq,Lclhcalipsoice,Lclhcalipsoun,         &
                                  Lclmcalipsoliq,Lclmcalipsoice,Lclmcalipsoun,         &
                                  Lcllcalipsoliq,Lcllcalipsoice,Lcllcalipsoun,         &
                                  Lclopaquecalipso,Lclthincalipso,Lclzopaquecalipso,   &
                                  Lclcalipsoopaque,Lclcalipsothin,Lclcalipsozopaque,   &
                                  Lclcalipsoopacity,Lclopaquetemp,Lclthintemp,         &
                                  Lclzopaquetemp,Lclopaquemeanz,Lclthinmeanz,          &
                                  Lclthinemis,Lclopaquemeanzse,Lclthinmeanzse,         &
                                  Lclzopaquecalipsose,LcfadDbze94,Ldbze94,Lparasolrefl,&
                                  Ltbrttov, Lptradarflag0,Lptradarflag1,Lptradarflag2, &
                                  Lptradarflag3,Lptradarflag4,Lptradarflag5,           &
                                  Lptradarflag6,Lptradarflag7,Lptradarflag8,           &
                                  Lptradarflag9,Lradarpia,Lwr_occfreq,Lcfodd,          &
                                  Npoints,Ncolumns,Nlevels,Nlvgrid,Nchan,x)
    ! Inputs
    logical,intent(in) :: &
        Lpctisccp,        & ! ISCCP mean cloud top pressure
        Lclisccp,         & ! ISCCP cloud area fraction
        Lboxptopisccp,    & ! ISCCP CTP in each column
        Lboxtauisccp,     & ! ISCCP optical epth in each column
        Ltauisccp,        & ! ISCCP mean optical depth
        Lcltisccp,        & ! ISCCP total cloud fraction
        Lmeantbisccp,     & ! ISCCP mean all-sky 10.5micron brightness temperature
        Lmeantbclrisccp,  & ! ISCCP mean clear-sky 10.5micron brightness temperature
        Lalbisccp,        & ! ISCCP mean cloud albedo
        LclMISR,          & ! MISR cloud fraction
        Lcltmodis,        & ! MODIS total cloud fraction
        Lclwmodis,        & ! MODIS liquid cloud fraction
        Lclimodis,        & ! MODIS ice cloud fraction
        Lclhmodis,        & ! MODIS high-level cloud fraction
        Lclmmodis,        & ! MODIS mid-level cloud fraction
        Lcllmodis,        & ! MODIS low-level cloud fraction
        Ltautmodis,       & ! MODIS total cloud optical thicknes
        Ltauwmodis,       & ! MODIS liquid optical thickness
        Ltauimodis,       & ! MODIS ice optical thickness
        Ltautlogmodis,    & ! MODIS total cloud optical thickness (log10 mean)
        Ltauwlogmodis,    & ! MODIS liquid optical thickness (log10 mean)
        Ltauilogmodis,    & ! MODIS ice optical thickness (log10 mean)
        Lreffclwmodis,    & ! MODIS liquid cloud particle size
        Lreffclimodis,    & ! MODIS ice particle size
        Lpctmodis,        & ! MODIS cloud top pressure
        Llwpmodis,        & ! MODIS cloud liquid water path
        Liwpmodis,        & ! MODIS cloud ice water path
        Lclmodis,         & ! MODIS cloud area fraction
        Latb532,          & ! CALIPSO attenuated total backscatter (532nm)
        Latb532gr,        & ! GROUND LIDAR @ 532NM attenuated total backscatter (532nm)
        Latb355,          & ! ATLID attenuated total backscatter (355nm)
        LlidarBetaMol532, & ! CALIPSO molecular backscatter (532nm)
        LlidarBetaMol532gr,&! GROUND LIDAR @ 532NM molecular backscatter (532nm)
        LlidarBetaMol355, & ! ATLID molecular backscatter (355nm)
        LcfadLidarsr532,  & ! CALIPSO scattering ratio CFAD
        LcfadLidarsr532gr,& ! GROUND LIDAR @ 532NM scattering ratio CFAD
        LcfadLidarsr355,  & ! ATLID scattering ratio CFAD
        Lclcalipso2,      & ! CALIPSO cloud fraction undetected by cloudsat
        Lclcalipso,       & ! CALIPSO cloud area fraction
        LclgrLidar532,    & ! GROUND LIDAR @ 532NM cloud area fraction
        Lclatlid,         & ! ATLID cloud area fraction
        Lclhcalipso,      & ! CALIPSO high-level cloud fraction
        Lcllcalipso,      & ! CALIPSO low-level cloud fraction
        Lclmcalipso,      & ! CALIPSO mid-level cloud fraction
        Lcltcalipso,      & ! CALIPSO total cloud fraction
        LclhgrLidar532,   & ! GROUND LIDAR @ 532NM high-level cloud fraction
        LcllgrLidar532,   & ! GROUND LIDAR @ 532NM low-level cloud fraction
        LclmgrLidar532,   & ! GROUND LIDAR @ 532NM mid-level cloud fraction
        LcltgrLidar532,   & ! GROUND LIDAR @ 532NM total cloud fraction
        Lclhatlid,        & ! ATLID high-level cloud fraction
        Lcllatlid,        & ! ATLID low-level cloud fraction
        Lclmatlid,        & ! ATLID mid-level cloud fraction
        Lcltatlid,        & ! ATLID total cloud fraction
        Lcltlidarradar,   & ! CALIPSO-CLOUDSAT total cloud fraction
        Lcloudsat_tcc,    & !
        Lcloudsat_tcc2,   & !
        Lclcalipsoliq,    & ! CALIPSO liquid cloud area fraction
        Lclcalipsoice,    & ! CALIPSO ice cloud area fraction
        Lclcalipsoun,     & ! CALIPSO undetected cloud area fraction
        Lclcalipsotmp,    & ! CALIPSO undetected cloud area fraction
        Lclcalipsotmpliq, & ! CALIPSO liquid cloud area fraction
        Lclcalipsotmpice, & ! CALIPSO ice cloud area fraction
        Lclcalipsotmpun,  & ! CALIPSO undetected cloud area fraction
        Lcltcalipsoliq,   & ! CALIPSO liquid total cloud fraction
        Lcltcalipsoice,   & ! CALIPSO ice total cloud fraction
        Lcltcalipsoun,    & ! CALIPSO undetected total cloud fraction
        Lclhcalipsoliq,   & ! CALIPSO high-level liquid cloud fraction
        Lclhcalipsoice,   & ! CALIPSO high-level ice cloud fraction
        Lclhcalipsoun,    & ! CALIPSO high-level undetected cloud fraction
        Lclmcalipsoliq,   & ! CALIPSO mid-level liquid cloud fraction
        Lclmcalipsoice,   & ! CALIPSO mid-level ice cloud fraction
        Lclmcalipsoun,    & ! CALIPSO mid-level undetected cloud fraction
        Lcllcalipsoliq,   & ! CALIPSO low-level liquid cloud fraction
        Lcllcalipsoice,   & ! CALIPSO low-level ice cloud fraction
        Lcllcalipsoun,    & ! CALIPSO low-level undetected cloud fraction
        Lclopaquecalipso, & ! CALIPSO opaque cloud cover (2D Map)
        Lclthincalipso,   & ! CALIPSO thin cloud cover (2D Map)
        Lclzopaquecalipso,& ! CALIPSO z_opaque altitude (opaque clouds only, 2D Map)
        Lclcalipsoopaque, & ! CALIPSO opaque cloud profiles 3D fraction
        Lclcalipsothin,   & ! CALIPSO thin cloud profiles 3D fraction
        Lclcalipsozopaque,& ! CALIPSO z_opaque 3D fraction
        Lclcalipsoopacity,& ! CALIPSO opacity 3D fraction
        Lclopaquetemp,    & ! CALIPSO opaque cloud temperature
        Lclthintemp,      & ! CALIPSO thin cloud temperature
        Lclzopaquetemp,   & ! CALIPSO z_opaque temperature
        Lclopaquemeanz,   & ! CALIPSO opaque cloud altitude
        Lclthinmeanz,     & ! CALIPSO thin cloud altitude
        Lclthinemis,      & ! CALIPSO thin cloud emissivity
        Lclopaquemeanzse,   & ! CALIPSO opaque cloud altitude with respect to SE
        Lclthinmeanzse,     & ! CALIPSO thin cloud altitude with respect to SE
        Lclzopaquecalipsose,& ! CALIPSO z_opaque altitude with respect to SE
        LcfadDbze94,      & ! CLOUDSAT radar reflectivity CFAD
        Ldbze94,          & ! CLOUDSAT radar reflectivity
        LparasolRefl,     & ! PARASOL reflectance
        Ltbrttov,         & ! RTTOV mean clear-sky brightness temperature
        Lptradarflag0,    & ! CLOUDSAT
        Lptradarflag1,    & ! CLOUDSAT
        Lptradarflag2,    & ! CLOUDSAT
        Lptradarflag3,    & ! CLOUDSAT
        Lptradarflag4,    & ! CLOUDSAT
        Lptradarflag5,    & ! CLOUDSAT
        Lptradarflag6,    & ! CLOUDSAT
        Lptradarflag7,    & ! CLOUDSAT
        Lptradarflag8,    & ! CLOUDSAT
        Lptradarflag9,    & ! CLOUDSAT
        Lradarpia,        & ! CLOUDSAT
        Lwr_occfreq,      & ! CloudSat+MODIS joint diagnostics
        Lcfodd              ! CloudSat+MODIS joint diagnostics

    integer,intent(in) :: &
        Npoints,         & ! Number of sampled points
        Ncolumns,        & ! Number of subgrid columns
        Nlevels,         & ! Number of model levels
        Nlvgrid,         & ! Number of levels in L3 stats computation
        Nchan              ! Number of RTTOV channels

    ! Outputs
    type(cosp_outputs),intent(out) :: &
        x           ! COSP output structure

    ! ISCCP simulator outputs
    if (Lboxtauisccp)    allocate(x%isccp_boxtau(Npoints,Ncolumns))
    if (Lboxptopisccp)   allocate(x%isccp_boxptop(Npoints,Ncolumns))
    if (Lclisccp)        allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
    if (Lcltisccp)       allocate(x%isccp_totalcldarea(Npoints))
    if (Lpctisccp)       allocate(x%isccp_meanptop(Npoints))
    if (Ltauisccp)       allocate(x%isccp_meantaucld(Npoints))
    if (Lmeantbisccp)    allocate(x%isccp_meantb(Npoints))
    if (Lmeantbclrisccp) allocate(x%isccp_meantbclr(Npoints))
    if (Lalbisccp)       allocate(x%isccp_meanalbedocld(Npoints))

    ! MISR simulator
    if (LclMISR) then
      allocate(x%misr_fq(Npoints, numMISRTauBins, numMISRHgtBins))
      ! *NOTE* These 3 fields are not output, but were part of the v1.4.0 cosp_misr, so
      !        they are still computed. Should probably have a logical to control these
      !        outputs.
      allocate(x%misr_dist_model_layertops(Npoints, numMISRHgtBins))
      allocate(x%misr_meanztop(Npoints))
      allocate(x%misr_cldarea(Npoints))
    endif

    ! MODIS simulator
    if (Lcltmodis)     allocate(x%modis_Cloud_Fraction_Total_Mean(Npoints))
    if (Lclwmodis)     allocate(x%modis_Cloud_Fraction_Water_Mean(Npoints))
    if (Lclimodis)     allocate(x%modis_Cloud_Fraction_Ice_Mean(Npoints))
    if (Lclhmodis)     allocate(x%modis_Cloud_Fraction_High_Mean(Npoints))
    if (Lclmmodis)     allocate(x%modis_Cloud_Fraction_Mid_Mean(Npoints))
    if (Lcllmodis)     allocate(x%modis_Cloud_Fraction_Low_Mean(Npoints))
    if (Ltautmodis)    allocate(x%modis_Optical_Thickness_Total_Mean(Npoints))
    if (Ltauwmodis)    allocate(x%modis_Optical_Thickness_Water_Mean(Npoints))
    if (Ltauimodis)    allocate(x%modis_Optical_Thickness_Ice_Mean(Npoints))
    if (Ltautlogmodis) allocate(x%modis_Optical_Thickness_Total_LogMean(Npoints))
    if (Ltauwlogmodis) allocate(x%modis_Optical_Thickness_Water_LogMean(Npoints))
    if (Ltauilogmodis) allocate(x%modis_Optical_Thickness_Ice_LogMean(Npoints))
    if (Lreffclwmodis) allocate(x%modis_Cloud_Particle_Size_Water_Mean(Npoints))
    if (Lreffclimodis) allocate(x%modis_Cloud_Particle_Size_Ice_Mean(Npoints))
    if (Lpctmodis)     allocate(x%modis_Cloud_Top_Pressure_Total_Mean(Npoints))
    if (Llwpmodis)     allocate(x%modis_Liquid_Water_Path_Mean(Npoints))
    if (Liwpmodis)     allocate(x%modis_Ice_Water_Path_Mean(Npoints))
    if (Lclmodis) then
      allocate(x%modis_Optical_Thickness_vs_Cloud_Top_Pressure(nPoints, numModisTauBins, numMODISPresBins))
      allocate(x%modis_Optical_thickness_vs_ReffLIQ(nPoints, numMODISTauBins, numMODISReffLiqBins))
      allocate(x%modis_Optical_Thickness_vs_ReffICE(nPoints, numMODISTauBins, numMODISReffIceBins))
    endif

    ! LIDAR simulator
    if (LlidarBetaMol532) allocate(x%calipso_beta_mol(Npoints,Nlevels))
    if (Latb532)          allocate(x%calipso_beta_tot(Npoints,Ncolumns,Nlevels))
    if (LcfadLidarsr532)  then
      allocate(x%calipso_srbval(SR_BINS+1))
      allocate(x%calipso_cfad_sr(Npoints,SR_BINS,Nlvgrid))
      allocate(x%calipso_betaperp_tot(Npoints,Ncolumns,Nlevels))
    endif
    if (Lclcalipso)       allocate(x%calipso_lidarcld(Npoints,Nlvgrid))
    if (Lclhcalipso .or. Lclmcalipso .or. Lcllcalipso .or. Lcltcalipso) then
      allocate(x%calipso_cldlayer(Npoints,LIDAR_NCAT))
    endif
    if (Lclcalipsoice .or. Lclcalipsoliq .or. Lclcalipsoun) then
      allocate(x%calipso_lidarcldphase(Npoints,Nlvgrid,6))
    endif
    if (Lclcalipsotmp .or. Lclcalipsotmpliq .or. Lclcalipsoice .or. Lclcalipsotmpun .or. Lclcalipsotmpice) then
      allocate(x%calipso_lidarcldtmp(Npoints,LIDAR_NTEMP,5))
    endif
    if (Lcllcalipsoice .or. Lclmcalipsoice .or. Lclhcalipsoice .or.                   &
      Lcltcalipsoice .or. Lcllcalipsoliq .or. Lclmcalipsoliq .or.                   &
      Lclhcalipsoliq .or. Lcltcalipsoliq .or. Lcllcalipsoun  .or.                   &
      Lclmcalipsoun  .or. Lclhcalipsoun  .or. Lcltcalipsoun) then
      allocate(x%calipso_cldlayerphase(Npoints,LIDAR_NCAT,6))
    endif
    if (Lclopaquecalipso .or. Lclthincalipso .or. Lclzopaquecalipso) then
      allocate(x%calipso_cldtype(Npoints,LIDAR_NTYPE))
    endif
    if (Lclopaquetemp .or. Lclthintemp .or. Lclzopaquetemp) then
      allocate(x%calipso_cldtypetemp(Npoints,LIDAR_NTYPE))
    endif
    if (Lclopaquemeanz .or. Lclthinmeanz) then
      allocate(x%calipso_cldtypemeanz(Npoints,2))
    endif
    if (Lclopaquemeanzse .or. Lclthinmeanzse .or. Lclzopaquecalipsose) then
      allocate(x%calipso_cldtypemeanzse(Npoints,3))
    endif
    if (Lclthinemis) then
      allocate(x%calipso_cldthinemis(Npoints))
    endif
    if (Lclcalipsoopaque .or. Lclcalipsothin .or. Lclcalipsozopaque .or. Lclcalipsoopacity) then
      allocate(x%calipso_lidarcldtype(Npoints, Nlvgrid, LIDAR_NTYPE+1))
    endif
    ! These 2 outputs are part of the calipso output type, but are not controlled by an
    ! logical switch in the output namelist, so if all other fields are on, then allocate
    if (LlidarBetaMol532 .or. Latb532      .or. LcfadLidarsr532 .or. Lclcalipso  .or.  &
      Lclcalipsoice    .or. Lclcalipsoliq  .or. Lclcalipsoun    .or. Lclcalipso2 .or.  &
      Lclhcalipso      .or. Lclmcalipso    .or. Lcllcalipso     .or. Lcltcalipso .or.  &
      Lclcalipsotmp    .or. Lclcalipsoice  .or. Lclcalipsotmpun .or.                   &
      Lclcalipsotmpliq .or. Lcllcalipsoice .or. Lclmcalipsoice  .or.                   &
      Lclhcalipsoice   .or. Lcltcalipsoice .or. Lcllcalipsoliq  .or.                   &
      Lclmcalipsoliq   .or. Lclhcalipsoliq .or. Lcltcalipsoliq  .or.                   &
      Lcllcalipsoun    .or. Lclmcalipsoun  .or. Lclhcalipsoun   .or. Lcltcalipsoun) then
      allocate(x%calipso_tau_tot(Npoints,Ncolumns,Nlevels))
      allocate(x%calipso_temp_tot(Npoints,Nlevels))
    endif

    ! GROUND LIDAR @ 532NM simulator
    if (LlidarBetaMol532gr) allocate(x%grLidar532_beta_mol(Npoints, Nlevels))
    if (Latb532gr)          allocate(x%grLidar532_beta_tot(Npoints, Ncolumns, Nlevels))
    if (LcfadLidarsr532gr) then
      allocate(x%grLidar532_srbval(SR_BINS+1))
      allocate(x%grLidar532_cfad_sr(Npoints, SR_BINS, Nlvgrid))
    endif
    if (LclgrLidar532)     allocate(x%grLidar532_lidarcld(Npoints, Nlvgrid))
    if (LclhgrLidar532 .or. LclmgrLidar532 .or. LcllgrLidar532 .or. LcltgrLidar532) then
      allocate(x%grLidar532_cldlayer(Npoints, LIDAR_NCAT))
    endif

    ! ATLID simulator
    if (LlidarBetaMol355) allocate(x%atlid_beta_mol(Npoints, Nlevels))
    if (Latb355)          allocate(x%atlid_beta_tot(Npoints, Ncolumns, Nlevels))
    if (LcfadLidarsr355) then
      allocate(x%atlid_srbval(SR_BINS+1))
      allocate(x%atlid_cfad_sr(Npoints,SR_BINS,Nlvgrid))
    endif
    if (Lclatlid)     allocate(x%atlid_lidarcld(Npoints, Nlvgrid))
    if (Lclhatlid .or. Lclmatlid .or. Lcllatlid .or. Lcltatlid) then
      allocate(x%atlid_cldlayer(Npoints, LIDAR_NCAT))
    endif

    ! PARASOL
    if (Lparasolrefl) then
      allocate(x%parasolPix_refl(Npoints, Ncolumns, PARASOL_NREFL))
      allocate(x%parasolGrid_refl(Npoints, PARASOL_NREFL))
    endif

    ! Cloudsat simulator
    if (Ldbze94)        allocate(x%cloudsat_Ze_tot(Npoints, Ncolumns,Nlevels))
    if (LcfadDbze94)    allocate(x%cloudsat_cfad_ze(Npoints, cloudsat_DBZE_BINS, Nlvgrid))
    if (Lptradarflag0 .or. Lptradarflag1 .or. Lptradarflag2 .or. Lptradarflag3 .or. &
      Lptradarflag4 .or. Lptradarflag5 .or. Lptradarflag6 .or. Lptradarflag7 .or. &
      Lptradarflag8 .or. Lptradarflag9) then
      allocate(x%cloudsat_precip_cover(Npoints, cloudsat_DBZE_BINS))
    endif
    if (Lradarpia) allocate(x%cloudsat_pia(Npoints))

    ! Combined CALIPSO/CLOUDSAT fields
    if (Lclcalipso2)    allocate(x%lidar_only_freq_cloud(Npoints,Nlvgrid))
    if (Lcltlidarradar) allocate(x%radar_lidar_tcc(Npoints))
    if (Lcloudsat_tcc) allocate(x%cloudsat_tcc(Npoints))
    if (Lcloudsat_tcc2) allocate(x%cloudsat_tcc2(Npoints))

    ! RTTOV
    if (Ltbrttov) allocate(x%rttov_tbs(Npoints,Nchan))

    ! Joint MODIS/CloudSat Statistics
    if (Lwr_occfreq)  allocate(x%wr_occfreq_ntotal(Npoints,WR_NREGIME))
    if (Lcfodd)       allocate(x%cfodd_ntotal(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS))

  end subroutine construct_cosp_outputs


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospIN(y)
    type(cosp_optical_inputs),intent(inout) :: y
    if (allocated(y%tau_067))             deallocate(y%tau_067)
    if (allocated(y%emiss_11))            deallocate(y%emiss_11)
    if (allocated(y%frac_out))            deallocate(y%frac_out)
    if (allocated(y%beta_mol_calipso))    deallocate(y%beta_mol_calipso)
    if (allocated(y%tau_mol_calipso))     deallocate(y%tau_mol_calipso)
    if (allocated(y%betatot_calipso))     deallocate(y%betatot_calipso)
    if (allocated(y%betatot_ice_calipso)) deallocate(y%betatot_ice_calipso)
    if (allocated(y%betatot_liq_calipso)) deallocate(y%betatot_liq_calipso)
    if (allocated(y%tautot_calipso))      deallocate(y%tautot_calipso)
    if (allocated(y%tautot_ice_calipso))  deallocate(y%tautot_ice_calipso)
    if (allocated(y%tautot_liq_calipso))  deallocate(y%tautot_liq_calipso)
    if (allocated(y%tautot_S_liq))        deallocate(y%tautot_S_liq)
    if (allocated(y%tautot_S_ice))        deallocate(y%tautot_S_ice)
    if (allocated(y%z_vol_cloudsat))      deallocate(y%z_vol_cloudsat)
    if (allocated(y%kr_vol_cloudsat))     deallocate(y%kr_vol_cloudsat)
    if (allocated(y%g_vol_cloudsat))      deallocate(y%g_vol_cloudsat)
    if (allocated(y%asym))                deallocate(y%asym)
    if (allocated(y%ss_alb))              deallocate(y%ss_alb)
    if (allocated(y%fracLiq))             deallocate(y%fracLiq)
    if (allocated(y%beta_mol_grLidar532)) deallocate(y%beta_mol_grLidar532)
    if (allocated(y%betatot_grLidar532))  deallocate(y%betatot_grLidar532)
    if (allocated(y%tau_mol_grLidar532))  deallocate(y%tau_mol_grLidar532)
    if (allocated(y%tautot_grLidar532))   deallocate(y%tautot_grLidar532)
    if (allocated(y%beta_mol_atlid))      deallocate(y%beta_mol_atlid)
    if (allocated(y%betatot_atlid))       deallocate(y%betatot_atlid)
    if (allocated(y%tau_mol_atlid))       deallocate(y%tau_mol_atlid)
    if (allocated(y%tautot_atlid))        deallocate(y%tautot_atlid)
    if (allocated(y%fracPrecipIce))      deallocate(y%fracPrecipIce)
  end subroutine destroy_cospIN


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cospstateIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cospstateIN(y)
    type(cosp_column_inputs),intent(inout) :: y

    if (allocated(y%sunlit))          deallocate(y%sunlit)
    if (allocated(y%skt))             deallocate(y%skt)
    if (allocated(y%land))            deallocate(y%land)
    if (allocated(y%at))              deallocate(y%at)
    if (allocated(y%pfull))           deallocate(y%pfull)
    if (allocated(y%phalf))           deallocate(y%phalf)
    if (allocated(y%qv))              deallocate(y%qv)
    if (allocated(y%o3))              deallocate(y%o3)
    if (allocated(y%hgt_matrix))      deallocate(y%hgt_matrix)
    if (allocated(y%u_sfc))           deallocate(y%u_sfc)
    if (allocated(y%v_sfc))           deallocate(y%v_sfc)
    if (allocated(y%lat))             deallocate(y%lat)
    if (allocated(y%lon))             deallocate(y%lon)
    if (allocated(y%emis_sfc))        deallocate(y%emis_sfc)
    if (allocated(y%cloudIce))        deallocate(y%cloudIce)
    if (allocated(y%cloudLiq))        deallocate(y%cloudLiq)
    if (allocated(y%seaice))          deallocate(y%seaice)
    if (allocated(y%fl_rain))         deallocate(y%fl_rain)
    if (allocated(y%fl_snow))         deallocate(y%fl_snow)
    if (allocated(y%tca))             deallocate(y%tca)
    if (allocated(y%hgt_matrix_half)) deallocate(y%hgt_matrix_half)
    if (allocated(y%surfelev))        deallocate(y%surfelev)

  end subroutine destroy_cospstateIN


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE destroy_cosp_outputs
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine destroy_cosp_outputs(y)
    type(cosp_outputs),intent(inout) :: y

    ! Deallocate and nullify
    if (associated(y%calipso_beta_mol))          then
      deallocate(y%calipso_beta_mol)
      nullify(y%calipso_beta_mol)
    endif
    if (associated(y%calipso_temp_tot))          then
      deallocate(y%calipso_temp_tot)
      nullify(y%calipso_temp_tot)
    endif
    if (associated(y%calipso_betaperp_tot))      then
      deallocate(y%calipso_betaperp_tot)
      nullify(y%calipso_betaperp_tot)
    endif
    if (associated(y%calipso_beta_tot))          then
      deallocate(y%calipso_beta_tot)
      nullify(y%calipso_beta_tot)
    endif
    if (associated(y%calipso_tau_tot))           then
      deallocate(y%calipso_tau_tot)
      nullify(y%calipso_tau_tot)
    endif
    if (associated(y%calipso_lidarcldphase))     then
      deallocate(y%calipso_lidarcldphase)
      nullify(y%calipso_lidarcldphase)
    endif
    if (associated(y%calipso_lidarcldtype))     then
      deallocate(y%calipso_lidarcldtype)
      nullify(y%calipso_lidarcldtype)
    endif
    if (associated(y%calipso_cldlayerphase))     then
      deallocate(y%calipso_cldlayerphase)
      nullify(y%calipso_cldlayerphase)
    endif
    if (associated(y%calipso_lidarcldtmp))       then
      deallocate(y%calipso_lidarcldtmp)
      nullify(y%calipso_lidarcldtmp)
    endif
    if (associated(y%calipso_cldlayer))          then
      deallocate(y%calipso_cldlayer)
      nullify(y%calipso_cldlayer)
    endif
    if (associated(y%calipso_cldtype))          then
      deallocate(y%calipso_cldtype)
      nullify(y%calipso_cldtype)
    endif
    if (associated(y%calipso_cldtypetemp))      then
      deallocate(y%calipso_cldtypetemp)
      nullify(y%calipso_cldtypetemp)
    endif
    if (associated(y%calipso_cldtypemeanz))     then
      deallocate(y%calipso_cldtypemeanz)
      nullify(y%calipso_cldtypemeanz)
    endif
    if (associated(y%calipso_cldtypemeanzse))   then
      deallocate(y%calipso_cldtypemeanzse)
      nullify(y%calipso_cldtypemeanzse)
    endif
    if (associated(y%calipso_cldthinemis))      then
      deallocate(y%calipso_cldthinemis)
      nullify(y%calipso_cldthinemis)
    endif
    if (associated(y%calipso_lidarcld))         then
      deallocate(y%calipso_lidarcld)
      nullify(y%calipso_lidarcld)
    endif
    if (associated(y%calipso_srbval))            then
      deallocate(y%calipso_srbval)
      nullify(y%calipso_srbval)
    endif
    if (associated(y%calipso_cfad_sr))          then
      deallocate(y%calipso_cfad_sr)
      nullify(y%calipso_cfad_sr)
    endif
    if (associated(y%grLidar532_beta_mol))     then
      deallocate(y%grLidar532_beta_mol)
      nullify(y%grLidar532_beta_mol)
    endif
    if (associated(y%grLidar532_beta_tot))     then
      deallocate(y%grLidar532_beta_tot)
      nullify(y%grLidar532_beta_tot)
    endif
    if (associated(y%grLidar532_cldlayer))     then
      deallocate(y%grLidar532_cldlayer)
      nullify(y%grLidar532_cldlayer)
    endif
    if (associated(y%grLidar532_lidarcld))     then
      deallocate(y%grLidar532_lidarcld)
      nullify(y%grLidar532_lidarcld)
    endif
    if (associated(y%grLidar532_cfad_sr))      then
      deallocate(y%grLidar532_cfad_sr)
      nullify(y%grLidar532_cfad_sr)
    endif
    if (associated(y%grLidar532_srbval))       then
      deallocate(y%grLidar532_srbval)
      nullify(y%grLidar532_srbval)
    endif
    if (associated(y%atlid_beta_mol))           then
      deallocate(y%atlid_beta_mol)
      nullify(y%atlid_beta_mol)
    endif
    if (associated(y%atlid_beta_tot))           then
      deallocate(y%atlid_beta_tot)
      nullify(y%atlid_beta_tot)
    endif
    if (associated(y%atlid_cldlayer))           then
      deallocate(y%atlid_cldlayer)
      nullify(y%atlid_cldlayer)
    endif
    if (associated(y%atlid_lidarcld))           then
      deallocate(y%atlid_lidarcld)
      nullify(y%atlid_lidarcld)
    endif
    if (associated(y%atlid_cfad_sr))            then
      deallocate(y%atlid_cfad_sr)
      nullify(y%atlid_cfad_sr)
    endif
    if (associated(y%atlid_srbval))             then
      deallocate(y%atlid_srbval)
      nullify(y%atlid_srbval)
    endif
    if (associated(y%parasolPix_refl))           then
      deallocate(y%parasolPix_refl)
      nullify(y%parasolPix_refl)
    endif
    if (associated(y%parasolGrid_refl))          then
      deallocate(y%parasolGrid_refl)
      nullify(y%parasolGrid_refl)
    endif
    if (associated(y%cloudsat_Ze_tot))           then
      deallocate(y%cloudsat_Ze_tot)
      nullify(y%cloudsat_Ze_tot)
    endif
    if (associated(y%cloudsat_cfad_ze))          then
      deallocate(y%cloudsat_cfad_ze)
      nullify(y%cloudsat_cfad_ze)
    endif
    if (associated(y%cloudsat_precip_cover))     then
      deallocate(y%cloudsat_precip_cover)
      nullify(y%cloudsat_precip_cover)
    endif
    if (associated(y%cloudsat_pia))              then
      deallocate(y%cloudsat_pia)
      nullify(y%cloudsat_pia)
    endif
    if (associated(y%cloudsat_tcc))              then
      deallocate(y%cloudsat_tcc)
      nullify(y%cloudsat_tcc)
    endif
    if (associated(y%cloudsat_tcc2))             then
      deallocate(y%cloudsat_tcc2)
      nullify(y%cloudsat_tcc2)
    endif
    if (associated(y%radar_lidar_tcc))           then
      deallocate(y%radar_lidar_tcc)
      nullify(y%radar_lidar_tcc)
    endif
    if (associated(y%cloudsat_tcc))              then
      deallocate(y%cloudsat_tcc)
      nullify(y%cloudsat_tcc)
    endif
    if (associated(y%cloudsat_tcc2))             then
      deallocate(y%cloudsat_tcc2)
      nullify(y%cloudsat_tcc2)
    endif
    if (associated(y%lidar_only_freq_cloud))     then
      deallocate(y%lidar_only_freq_cloud)
      nullify(y%lidar_only_freq_cloud)
    endif
    if (associated(y%isccp_totalcldarea))        then
      deallocate(y%isccp_totalcldarea)
      nullify(y%isccp_totalcldarea)
    endif
    if (associated(y%isccp_meantb))              then
      deallocate(y%isccp_meantb)
      nullify(y%isccp_meantb)
    endif
    if (associated(y%isccp_meantbclr))           then
      deallocate(y%isccp_meantbclr)
      nullify(y%isccp_meantbclr)
    endif
    if (associated(y%isccp_meanptop))            then
      deallocate(y%isccp_meanptop)
      nullify(y%isccp_meanptop)
    endif
    if (associated(y%isccp_meantaucld))          then
      deallocate(y%isccp_meantaucld)
      nullify(y%isccp_meantaucld)
    endif
    if (associated(y%isccp_meanalbedocld))       then
      deallocate(y%isccp_meanalbedocld)
      nullify(y%isccp_meanalbedocld)
    endif
    if (associated(y%isccp_boxtau))              then
      deallocate(y%isccp_boxtau)
      nullify(y%isccp_boxtau)
    endif
    if (associated(y%isccp_boxptop))             then
      deallocate(y%isccp_boxptop)
      nullify(y%isccp_boxptop)
    endif
    if (associated(y%isccp_fq))                  then
      deallocate(y%isccp_fq)
      nullify(y%isccp_fq)
    endif
    if (associated(y%misr_fq))                   then
      deallocate(y%misr_fq)
      nullify(y%misr_fq)
    endif
    if (associated(y%misr_dist_model_layertops)) then
      deallocate(y%misr_dist_model_layertops)
      nullify(y%misr_dist_model_layertops)
    endif
    if (associated(y%misr_meanztop))             then
      deallocate(y%misr_meanztop)
      nullify(y%misr_meanztop)
    endif
    if (associated(y%misr_cldarea))              then
      deallocate(y%misr_cldarea)
      nullify(y%misr_cldarea)
    endif
    if (associated(y%rttov_tbs))                 then
      deallocate(y%rttov_tbs)
      nullify(y%rttov_tbs)
    endif
    if (associated(y%modis_Cloud_Fraction_Total_Mean))                      then
      deallocate(y%modis_Cloud_Fraction_Total_Mean)
      nullify(y%modis_Cloud_Fraction_Total_Mean)
    endif
    if (associated(y%modis_Cloud_Fraction_Ice_Mean))                        then
      deallocate(y%modis_Cloud_Fraction_Ice_Mean)
      nullify(y%modis_Cloud_Fraction_Ice_Mean)
    endif
    if (associated(y%modis_Cloud_Fraction_Water_Mean))                      then
      deallocate(y%modis_Cloud_Fraction_Water_Mean)
      nullify(y%modis_Cloud_Fraction_Water_Mean)
    endif
    if (associated(y%modis_Cloud_Fraction_High_Mean))                       then
      deallocate(y%modis_Cloud_Fraction_High_Mean)
      nullify(y%modis_Cloud_Fraction_High_Mean)
    endif
    if (associated(y%modis_Cloud_Fraction_Mid_Mean))                        then
      deallocate(y%modis_Cloud_Fraction_Mid_Mean)
      nullify(y%modis_Cloud_Fraction_Mid_Mean)
    endif
    if (associated(y%modis_Cloud_Fraction_Low_Mean))                        then
      deallocate(y%modis_Cloud_Fraction_Low_Mean)
      nullify(y%modis_Cloud_Fraction_Low_Mean)
    endif
    if (associated(y%modis_Optical_Thickness_Total_Mean))                   then
      deallocate(y%modis_Optical_Thickness_Total_Mean)
      nullify(y%modis_Optical_Thickness_Total_Mean)
    endif
    if (associated(y%modis_Optical_Thickness_Water_Mean))                   then
      deallocate(y%modis_Optical_Thickness_Water_Mean)
      nullify(y%modis_Optical_Thickness_Water_Mean)
    endif
    if (associated(y%modis_Optical_Thickness_Ice_Mean))                     then
      deallocate(y%modis_Optical_Thickness_Ice_Mean)
      nullify(y%modis_Optical_Thickness_Ice_Mean)
    endif
    if (associated(y%modis_Optical_Thickness_Total_LogMean))                then
      deallocate(y%modis_Optical_Thickness_Total_LogMean)
      nullify(y%modis_Optical_Thickness_Total_LogMean)
    endif
    if (associated(y%modis_Optical_Thickness_Water_LogMean))                then
      deallocate(y%modis_Optical_Thickness_Water_LogMean)
      nullify(y%modis_Optical_Thickness_Water_LogMean)
    endif
    if (associated(y%modis_Optical_Thickness_Ice_LogMean))                  then
      deallocate(y%modis_Optical_Thickness_Ice_LogMean)
      nullify(y%modis_Optical_Thickness_Ice_LogMean)
    endif
    if (associated(y%modis_Cloud_Particle_Size_Water_Mean))                 then
      deallocate(y%modis_Cloud_Particle_Size_Water_Mean)
      nullify(y%modis_Cloud_Particle_Size_Water_Mean)
    endif
    if (associated(y%modis_Cloud_Particle_Size_Ice_Mean))                   then
      deallocate(y%modis_Cloud_Particle_Size_Ice_Mean)
      nullify(y%modis_Cloud_Particle_Size_Ice_Mean)
    endif
    if (associated(y%modis_Cloud_Top_Pressure_Total_Mean))                  then
      deallocate(y%modis_Cloud_Top_Pressure_Total_Mean)
      nullify(y%modis_Cloud_Top_Pressure_Total_Mean)
    endif
    if (associated(y%modis_Liquid_Water_Path_Mean))                         then
      deallocate(y%modis_Liquid_Water_Path_Mean)
      nullify(y%modis_Liquid_Water_Path_Mean)
    endif
    if (associated(y%modis_Ice_Water_Path_Mean))                            then
      deallocate(y%modis_Ice_Water_Path_Mean)
      nullify(y%modis_Ice_Water_Path_Mean)
    endif
    if (associated(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure))        then
      deallocate(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
      nullify(y%modis_Optical_Thickness_vs_Cloud_Top_Pressure)
    endif
    if (associated(y%modis_Optical_thickness_vs_ReffLIQ))                   then
      deallocate(y%modis_Optical_thickness_vs_ReffLIQ)
      nullify(y%modis_Optical_thickness_vs_ReffLIQ)
    endif
    if (associated(y%modis_Optical_thickness_vs_ReffICE))                   then
      deallocate(y%modis_Optical_thickness_vs_ReffICE)
      nullify(y%modis_Optical_thickness_vs_ReffICE)
    endif
    if (associated(y%cfodd_ntotal)) then
      deallocate(y%cfodd_ntotal)
      nullify(y%cfodd_ntotal)
    endif
    if (associated(y%wr_occfreq_ntotal)) then
      deallocate(y%wr_occfreq_ntotal)
      nullify(y%wr_occfreq_ntotal)
    endif

  end subroutine destroy_cosp_outputs

end module fms_cosp_interface_mod
