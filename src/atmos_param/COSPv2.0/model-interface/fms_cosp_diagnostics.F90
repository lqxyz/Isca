module fms_cosp_diags_mod
! Modified from https://github.com/NOAA-GFDL/AM4/blob/master/src/atmos_param/cosp/cosp_diagnostics.F90

use mpp_mod,                  only: input_nml_file
use fms_mod,                  only: open_namelist_file, open_file,  &
                                    close_file, error_mesg, FATAL, &
                                    file_exist, mpp_pe, mpp_root_pe,   &
                                    check_nml_error, write_version_number,&
                                    stdlog
use time_manager_mod,         only: set_date, time_type, operator (+), &
                                    operator(-), operator(<),    &
                                    operator(>), operator(<=), &
                                    operator(>=),  get_date, print_date, &
                                    get_calendar_type, NOLEAP, &
                                    assignment(=), set_time
use diag_grid_mod,            only: get_local_indexes2
use diag_manager_mod,         only: register_diag_field, send_data, diag_axis_init

use fms_cosp_config_mod
use fms_cosp_io_mod,          only:  map_point_to_ll

use cosp_kinds,               only: wp
use mod_cosp,                 only: cosp_outputs
use netcdf
USE MOD_COSP_CONFIG,          only: R_UNDEF, LIDAR_NCAT, SR_BINS, PARASOL_NREFL, cloudsat_DBZE_BINS, &
                                    numMODISReffIceBins, numMODISReffLiqBins, ntau, tau_binBounds, tau_binCenters, &
                                    tau_binEdges, npres, pres_binBounds, pres_binCenters, pres_binEdges, nhgt,     &
                                    hgt_binBounds, hgt_binCenters, hgt_binEdges, vgrid_z,      &
                                    reffICE_binCenters, reffLIQ_binCenters, cloudsat_binCenters, PARASOL_SZA,      &
                                    calipso_binCenters, grLidar532_binCenters, atlid_binCenters,                   &
                                    CFODD_NDBZE, CFODD_HISTDBZE, CFODD_HISTDBZEcenters,                            &
                                    CFODD_NICOD, CFODD_HISTICOD, CFODD_HISTICODcenters ! Nlvgrid,

IMPLICIT NONE

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
              id_boxtauisccp, id_boxptopisccp, & ! id_clisccp,
              id_parasolrefl, id_parasolrefl_sat, &
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
  integer, dimension(npres) :: id_clisccp

  ! For axes
  integer :: id_column_idx, id_lev_idx, id_levSat_idx, id_bnds, &
             id_tau7, id_tau7_bnds, id_pres7, id_pres7_bnds, &
             id_hgt16, id_hgt16_bnds

  integer, dimension(18)  :: cosp_axes
  real  :: missing_value = R_UNDEF

  ! ---------------- End of declaration of variables --------------

contains

!######################################################################

subroutine cosp_diag_field_init(Time, axes, Nlevels, Ncolumns, Nlvgrid, cospOUT, mod_name)

  type(time_type), intent(in)       :: Time
  integer, dimension(4), intent(in) :: axes
  integer, intent(in)               :: Nlevels, Ncolumns, Nlvgrid
  type(cosp_outputs), intent(in)    :: cospOUT   ! COSP simulator outputs
  character(len=10),  intent(in)    :: mod_name

  real :: level_ax(Nlevels)
  real :: column_ax(Ncolumns)
  real :: lvgrid_ax(Nlvgrid)
  real :: bnds(2) = (/1.0, 2.0/)
  real :: SR_bins_ax(SR_BINS)

  integer, dimension(3) :: columnindx = (/1,2,5/)
  integer, dimension(3) :: tauindx = (/1,2,8/)

  integer :: i
  character(len=2) :: chvers
  character(len=8) :: chvers2, chvers3

  ! ======================= Define cosp axes ======================= !

  !--------------------------------------------------------------------
  ! define the varisous axes needed for this data.
  !--------------------------------------------------------------------
   cosp_axes(1:4) = axes(1:4)

  !--------------------------------------------------------------------
  ! a subcolumn counter:
  !--------------------------------------------------------------------
  do i=1,Ncolumns
    column_ax(i) = real(i)
  end do
  cosp_axes(5) = diag_axis_init('column_idx', column_ax, &
          'subcol', 'n', 'subcolumn number', set_name=mod_name)

  !--------------------------------------------------------------------
  ! a level counter:
  !--------------------------------------------------------------------
  do i=1,Nlevels
    level_ax(i) = real(i)
  end do
  cosp_axes(6) = diag_axis_init('lev_idx', level_ax, &
        'levels', 'n', 'level number', set_name=mod_name)

  !--------------------------------------------------------------------
  ! a lvgrid counter:
  !--------------------------------------------------------------------
  do i=1,Nlvgrid
    lvgrid_ax(i) = real(i)
  end do
  cosp_axes(7) = diag_axis_init('levStat_idx', lvgrid_ax, &
        'vgrid_levels', 'n', 'vgrid level number', set_name=mod_name)

  cosp_axes(8) = diag_axis_init('tau7', real(tau_binCenters), &
        'tau_binCenters', 'n', 'isccp tau category', set_name=mod_name)
  cosp_axes(9) = diag_axis_init('bnds', bnds, &
        'bnds', 'n', 'isccp bnds', set_name=mod_name)
  cosp_axes(10)= diag_axis_init('pres7', real(pres_binCenters), &
        'pres_binCenters', 'n', 'isccp pres category', set_name=mod_name)
  cosp_axes(11) = diag_axis_init('hgt16', real(hgt_binCenters), &
        'hgt_binCenters', 'n', 'MISR cloud-top height bins', set_name=mod_name)

  do i=1,SR_BINS
    SR_bins_ax(i) = real(i)
  end do
  cosp_axes(12) = diag_axis_init('SR_bins_idx', SR_bins_ax, &
        'SR_bins', 'n', 'Number of bins (backscattering coefficient) '// &
        'in CALOPSO LIDAR simulator', set_name=mod_name)

  cosp_axes(13) = diag_axis_init('PARASOL_NREFL', real(PARASOL_SZA), &
        'PARASOL_SZA', 'n', 'PARASOL simulator: Number of angles in LUT', set_name=mod_name)
  cosp_axes(14) = diag_axis_init('cloudsat_DBZE_BINS', real(cloudsat_binCenters), &
        'cloudsat_binCenters', 'n', 'CLOUDSAT: dBZe bins in histogram (cfad)', set_name=mod_name)
  cosp_axes(15) = diag_axis_init('RELIQ_MODIS', real(reffLIQ_binCenters), &
        'reffLIQ_binCenters', 'n', 'MODIS: Effective radius bin centers for LIQ', set_name=mod_name)
  cosp_axes(16) = diag_axis_init('REICE_MODIS', real(reffICE_binCenters), &
        'reffICE_binCenters', 'n', 'MODIS: Effective radius bin centers for ICE', set_name=mod_name)
  cosp_axes(17) = diag_axis_init('CFODD_NDBZE', real(CFODD_HISTDBZEcenters), &
        'CFODD_HISTDBZEcenters', 'n', 'CFODD dBZe bins', set_name=mod_name)
  cosp_axes(18) = diag_axis_init('CFODD_NICOD', real(CFODD_HISTICODcenters), &
        'CFODD_HISTICODcenters', 'n', 'CFODD ICOD bins', set_name=mod_name)

  ! ======================= Define outputs ======================= !
  ! ! Joint-histogram axis
  ! ! Tau
  ! id_tau7 = register_diag_field &
  !     (mod_name, 'tau7', cosp_axes(8), Time, &
  !     'cloud ptical depth bin centers', &
  !     '1', mask_variant=.true., missing_value=missing_value)

  ! Tau Edges
  id_tau7_bnds = register_diag_field &
      (mod_name, 'tau7_bnds', cosp_axes((/9,8/)), Time, &
      'cloud optical depth bin edges', &
      '1', mask_variant=.true., missing_value=missing_value)

  ! ! Pressure
  ! id_pres7 = register_diag_field &
  !     (mod_name, 'pres7', cosp_axes(10), Time, &
  !     'air pressure bin centers', &
  !     '1', mask_variant=.true., missing_value=missing_value)

  ! Pressure Edges
  id_pres7_bnds = register_diag_field &
      (mod_name, 'pres7_bnds', cosp_axes(9:10), Time, &
      'air pressure bin edges', &
      '1', mask_variant=.true., missing_value=missing_value)

  ! ! Height
  ! id_hgt16 = register_diag_field &
  !     (mod_name, 'hgt16', cosp_axes(11), Time, &
  !     'altitude bin centers', &
  !     '1', mask_variant=.true., missing_value=missing_value)

  ! Height Edges
  id_hgt16_bnds = register_diag_field &
      (mod_name, 'hgt16_bnds', cosp_axes((/9,11/)), Time, &
      'altitude bin edges', &
      '1', mask_variant=.true., missing_value=missing_value)

  ! ! Levels
  ! id_lev_idx = register_diag_field &
  !     (mod_name, 'lev_idx', cosp_axes(6), Time, &
  !     'level indices', &
  !     '1', mask_variant=.true., missing_value=missing_value)

  ! ! Levels for statistical diagnostics (lidar and radar)
  ! id_levSat_idx = register_diag_field &
  !     (mod_name, 'levStat_idx', cosp_axes(7), Time, &
  !     'level indices', &
  !     '1', mask_variant=.true., missing_value=missing_value)

  ! ! Subcolumms
  ! id_column_idx = register_diag_field &
  !     (mod_name, 'cosp_scol', cosp_axes(5), Time, &
  !     'subcolumn indices', &
  !     '1', mask_variant=.true., missing_value=missing_value)

  ! ! Bnds
  ! id_bnds = register_diag_field &
  !     (mod_name, 'bnds', cosp_axes(9), Time, &
  !     'bounds', &
  !     '1', mask_variant=.true., missing_value=missing_value)


  ! ======================= ISCCP simulator outputs ======================= !

  ! if (Lclisccp)        allocate(x%isccp_fq(Npoints,numISCCPTauBins,numISCCPPresBins))
  if (associated(cospOUT%isccp_totalcldarea)) then
    id_cltisccp = register_diag_field &
      (mod_name, 'cltisccp', axes(1:2), Time, &
      'ISCCP Total Cloud Fraction', &
      '%', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_meanptop)) then
    id_ctpisccp = register_diag_field &
      (mod_name, 'ctpisccp', axes(1:2), Time, &
      'Mean Cloud Top Pressure *CPCT as Calculated by the ISCCP Simulator', &
      'Pa', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_meantaucld)) then
    id_tauisccp = register_diag_field &
      (mod_name, 'tauisccp', axes(1:2), Time, &
      'Mean Optical Depth *CPCT as Calculated by the ISCCP Simulator', &
      'dimensionless', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_meantb)) then
    id_tbisccp = register_diag_field &
      (mod_name, 'tbisccp', axes(1:2), Time, &
      'Mean All-sky 10.5 micron brightness temp -- ISCCP Simulator', &
      'deg K', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_meantbclr)) then
    id_tbclrisccp = register_diag_field &
      (mod_name, 'tbclrisccp', axes(1:2), Time, &
      'Mean Clr-sky 10.5 micron brightness temp -- ISCCP Simulator', &
      'deg K', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_meanalbedocld)) then
    id_albisccp = register_diag_field &
      (mod_name, 'albisccp', axes(1:2), Time, &
      'Mean Cloud Albedo -- ISCCP Simulator', &
      '1', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_boxtau)) then
    ! cosp_scol
    id_boxtauisccp = register_diag_field &
      (mod_name, 'boxtauisccp', cosp_axes(columnindx), Time, &
      'ISCCP Subcolumn Optical Depth', &
      '1', mask_variant=.true., missing_value=missing_value)
  endif
  if (associated(cospOUT%isccp_boxptop)) then
    id_boxptopisccp = register_diag_field &
      (mod_name, 'boxptopisccp', cosp_axes(columnindx), Time, &
      'ISCCP Subcolumn Cloud Top Pressure', &
      'Pa', mask_variant=.true., missing_value=missing_value)
  endif
  ! if (associated(cospOUT%isccp_fq)) then
  !   id_clisccp = register_diag_field &
  !     (mod_name, 'clisccp', cosp_axes((/1,2,8,10/)), Time, &
  !     'ISCCP joint-PDF of cloud top pressure and optical depth', &
  !     '%', mask_variant=.true., missing_value=missing_value)
  ! endif

  if (associated(cospOUT%isccp_fq)) then
    do i=1,npres
      write(chvers, '(i1)') i
      write(chvers2, '(i6)') INT(pres_binEdges(1,i)*1.0e-02)
      write(chvers3, '(i6)') INT(pres_binEdges(2,i)*1.0e-02)
      id_clisccp(i) = register_diag_field &
        (mod_name, 'clisccp_'// trim(chvers), cosp_axes(tauindx), Time, &
        'ISCCP Cld Frac for clouds between ' // trim(chvers2) // ' and' // trim(chvers3) // ' hPa', &
        '%', mask_variant=.true., missing_value=missing_value)
    end do

  endif

  ! ======================= CALIPSO simulator outputs ======================= !

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

end subroutine cosp_diag_field_init

!####################################################################

subroutine cosp_output_fields(Time_diag, Nlon, Nlat, Ncolumns, Nlevels, geomode, cospOUT)
  type(time_type), intent(in)    :: Time_diag
  integer, intent(in)            :: Nlat, Nlon, Ncolumns, Nlevels
  integer, intent(in)            :: geomode
  type(cosp_outputs), intent(in) :: cospOUT ! COSP simulator outputs

  logical :: used

  real, dimension(Nlon,Nlat)            :: y2
  real, dimension(Nlon,Nlat,Nlevels)    :: y3_levs
  real, dimension(Nlon,Nlat,Ncolumns)   :: y3_cols
  real, dimension(Nlon,Nlat,ntau,npres) :: y4_tau_pres  ! (ntau, npres) or (npres, ntau)??

  integer :: i

  ! real :: level_ax(Nlevels)
  ! real :: column_ax(Ncolumns)
  ! real :: lvgrid_ax(Nlvgrid)
  ! real :: bnds(2) = (/1.0, 2.0/)
  ! real :: SR_bins_ax(SR_BINS)
  !
  ! do i=1,Nlevels
  !   level_ax(i) = real(i)
  ! end do
  ! do i=1,Ncolumns
  !   column_ax(i) = real(i)
  ! end do
  ! do i=1,Nlvgrid
  !   lvgrid_ax(i) = real(i)
  ! end do
  ! do i=1,SR_BINS
  !   SR_bins_ax(i) = real(i)
  ! end do

  ! ! Cast data from type real(wp) to real.

  ! ! ======================================= Axis related ======================================= !
  ! if (id_tau7 > 0) then
  !   used = send_data(id_tau7, real(tau_binCenters), Time_diag)
  ! endif

  if (id_tau7_bnds > 0) then
    used = send_data(id_tau7_bnds, real(tau_binEdges), Time_diag, &
          mask=real(tau_binEdges)/=missing_value)
  endif

  ! if (id_pres7 > 0) then
  !   used = send_data(id_pres7, real(pres_binCenters), Time_diag)
  ! endif

  if (id_pres7_bnds > 0) then
    used = send_data(id_pres7_bnds, real(pres_binEdges), Time_diag, &
          mask=real(pres_binEdges)/=missing_value)
  endif

  ! if (id_hgt16 > 0) then
  !   used = send_data(id_hgt16, real(hgt_binCenters), Time_diag)
  ! endif

  if (id_hgt16_bnds > 0) then
    used = send_data(id_hgt16_bnds, real(hgt_binEdges), Time_diag, &
          mask=real(hgt_binEdges)/=missing_value)
  endif

  ! if (id_lev_idx > 0) then
  !   used = send_data(id_lev_idx, level_ax, Time_diag)
  ! endif

  ! if (id_levSat_idx > 0) then
  !   used = send_data(id_levSat_idx, lvgrid_ax, Time_diag)  !vgrid_z,
  ! endif

  ! if (id_column_idx > 0) then
  !   used = send_data(id_column_idx, column_ax, Time_diag)  !vgrid_z,
  ! endif

  ! if (id_bnds > 0) then
  !   used = send_data(id_bnds, bnds, Time_diag)  !vgrid_z,
  ! endif


  ! =======================================  ISCCP ======================================= !
  ! 2d arrays (i,j):
  if (id_cltisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%isccp_totalcldarea), y2=y2)
    used = send_data(id_cltisccp, y2, Time_diag, mask=y2/=missing_value)
  endif
  if (id_ctpisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%isccp_totalcldarea), y2=y2)
    used = send_data(id_ctpisccp, y2, Time_diag, mask=y2/=missing_value)
  endif
  if (id_tauisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%isccp_meantaucld), y2=y2)
    used = send_data(id_tauisccp, y2, Time_diag, mask=y2/=missing_value)
  endif
  if (id_tbisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%isccp_totalcldarea), y2=y2)
    used = send_data(id_tbisccp, y2, Time_diag, mask=y2/=missing_value)
  endif
  if (id_tbclrisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%isccp_totalcldarea), y2=y2)
    used = send_data(id_tbclrisccp, y2, Time_diag, mask=y2/=missing_value)
  endif

  ! 3d arrays (i,j,columns):
  if (id_boxtauisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x2=real(cospOUT%isccp_boxtau), y3=y3_cols)
    used = send_data(id_boxtauisccp, y3_cols, Time_diag, mask=y3_cols/=missing_value)
  endif
  if (id_boxptopisccp > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x2=real(cospOUT%isccp_boxptop), y3=y3_cols)
    used = send_data(id_boxptopisccp, y3_cols, Time_diag, mask=y3_cols/=missing_value)
  endif

  ! FMS send data does not support 4d data
  ! 4d array (i,j, isccp_tau,isccp_press):
  ! if (id_clisccp > 0) then
  !   call map_point_to_ll(Nlon, Nlat, geomode, x3=real(cospOUT%isccp_fq), y4=y4_tau_pres)
  !   used = send_data(id_clisccp, y4_tau_pres, Time_diag) ! mask=y4/=missing_value)
  ! endif


  if (id_clisccp(1) > 0) then
     call map_point_to_ll(Nlon, Nlat, geomode, x3=real(cospOUT%isccp_fq), y4=y4_tau_pres)
    do i=1,npres
      used = send_data(id_clisccp(i), y4_tau_pres(:,:,:,i), Time_diag, mask=y4_tau_pres(:,:,:,i)/=missing_value) ! mask=y4/=missing_value)
    enddo
  endif

  ! =======================================  CALIPSO ======================================= !

  if (id_cllcalipso > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,1)), y2=y2)
    used = send_data(id_cllcalipso, y2, Time_diag)
  endif

  if (id_clmcalipso > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,2)), y2=y2)
    used = send_data(id_clmcalipso, y2, Time_diag, mask=y2/=missing_value)
  endif

  if (id_clhcalipso > 0) then
    call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,3)), y2=y2)
    used = send_data(id_clhcalipso, y2, Time_diag, mask=y2/=missing_value)
  endif

  ! if (id_cltcalipso > 0) then
  !   call map_point_to_ll(Nlon, Nlat, geomode, x1=real(cospOUT%calipso_cldlayer(:,4)), y2=y2)
  !   used = send_data(id_cltcalipso, y2, Time_diag)
  ! endif

end subroutine cosp_output_fields

end module fms_cosp_diags_mod

