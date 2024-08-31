!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!-----------------------------------------------------------------
! This subroutine prepares information required for SVM to predict
! brightnbess temperature.
! This subroutine is called by subroutine noahmp36_getsnwdpred
!
! Written by Yonghwan Kwon (01 Aug. 2017)
!-----------------------------------------------------------------

subroutine noahmp36_LIS_SVM(n, k, swe_obsgrid, slwc_obsgrid, tp1_obsgrid, &
                           DTB_10H_36H, DTB_10V_36V, &
                           DTB_18H_36H, DTB_18V_36V, &
                           TB_10H, TB_10V, &
                           TB_18H, TB_18V, &
                           TB_36H, TB_36V)

! !USES:
  use ESMF
  use LIS_logMod,   only: LIS_logunit, LIS_verify
  use LIS_coreMod
  use noahmp36_lsmMod
  use LIS_mpiMod
  use mwSVM_routines
  use LIS_DAobservationsMod

  implicit none
! !ARGUMEMTS:
  integer, intent(in)                    :: n
  integer, intent(in)                    :: k 
  real,dimension(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n)) :: swe_obsgrid, slwc_obsgrid, tp1_obsgrid, tsurf_obsgrid
  real,dimension(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n)) :: DTB_10H_36H, DTB_10V_36V, DTB_18H_36H, DTB_18V_36V
  real,dimension(LIS_rc%obs_ngrid(k),LIS_rc%nensem(n)) :: TB_10H, TB_10V, TB_18H, TB_18V, TB_36H, TB_36V 

!EOP
  integer                                :: i,m,obsgrid_i,svm_i,svmlat_i
  integer                                :: svm_lon_id, svm_lat_id
  integer                                :: nlat_ease, nlon_ease
  integer                                :: num_of_inputs
  integer                                :: myid, suffix_id, ierr
  integer                                :: year, month, day, dofyr, mon
  integer, dimension(12)                 :: num_days_month
  real                                   :: predict_10H_36H, predict_10V_36V
  real                                   :: predict_18H_36H, predict_18V_36V
  real                                   :: predict_10H, predict_10V
  real                                   :: predict_18H, predict_18V
  real                                   :: predict_36H, predict_36V
  real                                   :: obs_lat, obs_lon
  real                                   :: svm_lon_interval
  real                                   :: aa, bb
  real, dimension(:),     allocatable    :: EASE_latsb     !EASE Grid south bound
  character(len=9)                       :: training_period
  character(len=2)                       :: training_target
  character(len=9)                       :: svm_domain
  character(len=4)                       :: scaling_option
  character(len=3)                       :: DOY
  character(len=4)                       :: YYYY
  character(len=2)                       :: MM, DD
  character(len=2)                       :: str_num_of_inputs, str_myid
  character(len=6)                       :: str_ind_obsgrid
  character(:), allocatable              :: case_directory, svm_directory
  character(:), allocatable              :: tmpfilename, Tb_output_filename
  logical                                :: forest_decoup_option, atm_decoup_option
  logical                                :: Tb_output_file_exist
  logical                                :: debug

  ! debug option
  debug = .false.

  !case_directory = '/discover/nobackup/projects/hma/users/ykwon/LIS/kyh_case/noahmp36_SynDA_kyh_02/'
  case_directory = '/discover/nobackup/projects/hma/users/ykwon/LIS/kyh_case/noahmp36_HMA/'
  !svm_directory  = '/discover/nobackup/ykwon/SVM/'
  svm_directory  = '/discover/nobackup/projects/hma/users/ykwon/SVM/'

  !---------------------------------------------------------------------
  ! Define some information for SVM
  ! (These need to be included in lis.config file later)
  ! Example,
  ! num_of_inputs == 1,2,3,4,5,6,7,8,9,10,11
  ! training_period == 'fortnight', 'month', 'seasonal'
  ! training_target == 'Tb','db'
  ! scaling_option == 'none','linear','standardization','unitvector'

  !nlat_ease        = 584    !EASE v2
  nlat_ease        = 586    !EASE v1
  !nlon_ease        = 1388   !EASE v2
  nlon_ease        = 1383   !EASE v1
  !nSVM_grid        = nlat_ease * nlon_ease  
  svm_lon_interval = 0.2594

  !account_name         = 'yhkwon'
  num_of_inputs        = 3       !4
  training_period      = 'fortnight'
  training_target      = 'db'   !'Tb'   !'db'
  svm_domain           = 'Indus'  !'SynDA'    !'WColorado'   !'Indus'
  scaling_option       = 'none'
  forest_decoup_option = .false.
  atm_decoup_option    = .false.
  !---------------------------------------------------------------------

  ! read EASE_Grid latitude information from text file
  allocate(EASE_latsb(nlat_ease))
  !open(22,file="EASE_Grid_v2_25km_southbound.txt",form="formatted",status="old",action="read")  !EASE v2
  open(22,file="EASE_Grid_southbound.txt",form="formatted",status="old",action="read")   !EASE v1
  do svmlat_i = 1, nlat_ease
     read(22,*) EASE_latsb(svmlat_i)
     if (debug) then
        !write(LIS_logunit,*) EASE_latsb(svmlat_i)
     endif
  enddo
  close(22)

  ! get model date and time information
  year  = LIS_rc%yr
  month = LIS_rc%mo
  day   = LIS_rc%da

  if (debug) then
     write(LIS_logunit,*) 'year=', year
     write(LIS_logunit,*) 'month=', month
     write(LIS_logunit,*) 'day=', day
  endif

  ! get day of year
  aa = real(year) / 4.
  bb = aa - floor(aa)

  if (debug) then
     write(LIS_logunit,*) 'aa=', aa
     write(LIS_logunit,*) 'bb=', bb
  endif

  if (bb == 0) then
     num_days_month = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  else
     num_days_month = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)
  endif

  if (debug) then
     write(LIS_logunit,*) 'num_days_month=', num_days_month
  endif 

  if (month == 1) then
     dofyr = day
  else
     dofyr = 0
     do mon = 1,month-1
        dofyr = dofyr + num_days_month(mon)
     enddo
     dofyr = dofyr + day
  endif
  if (debug) then
     write(LIS_logunit,*) 'dofyr=', dofyr
  endif  

#if (defined SPMD)
  call MPI_COMM_RANK( LIS_mpi_comm, myid, ierr )
#else
  ierr = 0
  myid = 0
#endif

  write(DOY, '(I3.3)') dofyr
  write(YYYY, '(I4.4)') year
  write(MM, '(I2.2)') month
  write(DD, '(I2.2)') day
  write(str_num_of_inputs, '(I2.2)') num_of_inputs

  do i=1,LIS_rc%obs_ngrid(k)
     do m=1,LIS_rc%nensem(n)

        ! find latitude and longitude of observation and matching it to SVM grid
        ! This is a temporary procedure before finilizing this routine.
        !----------------------------------
        obs_lat = LIS_obs_domain(n,k)%lat(i)
        obs_lon = LIS_obs_domain(n,k)%lon(i)

        if (debug) then
           write(LIS_logunit,*) 'obs_lat=', obs_lat
           write(LIS_logunit,*) 'obs_lon=', obs_lon
        endif

        svm_lon_id = ceiling((180 + obs_lon) / svm_lon_interval)
        if (svm_lon_id > nlon_ease) then
           svm_lon_id = nlon_ease
        endif

        do svmlat_i = 1, nlat_ease
           if (obs_lat >= EASE_latsb(svmlat_i)) then
              svm_lat_id = svmlat_i
              exit
           endif
        enddo

        svm_i = svm_lon_id + (nlon_ease * (svm_lat_id - 1))

        if (debug) then
           write(LIS_logunit,*) 'svm_lon_id=', svm_lon_id
           write(LIS_logunit,*) 'svm_lat_id=', svm_lat_id
           write(LIS_logunit,*) 'svm_i=', svm_i
        endif
        !---------------------------------- 

        suffix_id = myid + 10
        obsgrid_i = i
        write(str_myid, '(I2.2)') suffix_id
        write(str_ind_obsgrid, '(I6.6)') obsgrid_i

        if (debug) then
           write(LIS_logunit,*) 'obsgrid_i=', obsgrid_i
        endif

        tmpfilename = case_directory//'OUTPUT/SVM_inputs/LIS_svm_' &
                      //DOY//'vy'//YYYY//'_'//str_num_of_inputs//'_obsgid'//str_ind_obsgrid//'_'//str_myid     
  
        if (debug) then
           write(LIS_logunit,*) 'tmpfilename=', tmpfilename
        endif
        if (swe_obsgrid(i,m) > 0.01) then    !TB is estimated and thus SWE is updated
                                             !only when prior SWE is greater than the specified threshold SWE value (0.01 m)
                                             !SVM has been trained using SWE larger than this thereshold

           open(unit=60+myid, file=tmpfilename, action='write')

           ! Assemble right formatted model states for mwSVM to read
           ! swe_obsgrid: Snow water equivalent (m)
           ! slwc_obsgrid: Snow liquid water content (kg/m2)
           ! tp1_obsgrid: top layer soil temperature
           ! tsurf_obsgrid: surface temperature
           ! rhosn1_obsgrid: snow density (top layer)
           ! rhosn2_obsgrid: snow density (mid layer)
           ! rhosn3_obsgrid: snow density (bottom layer)
           ! tpsn3_obsgrid: bottom snow temperature
           ! tpsn1_obsgrid: top snow temperature
           ! tair_obsgrid: air temperature

           if (debug) then
              write(LIS_logunit,*) 'i=', i
              write(LIS_logunit,*) 'm=', m
              write(LIS_logunit,*) 'swe_obsgrid=', swe_obsgrid(i,m)
              write(LIS_logunit,*) 'slwc_obsgrid=', slwc_obsgrid(i,m)
              write(LIS_logunit,*) 'tp1_obsgrid=', tp1_obsgrid(i,m)
           endif

           select case (num_of_inputs)
              case (1)
                 !write(60+myid,  '(I1.1, I2.1, A1, F18.7)') 0, 1,':',swe_catch/100.
              case (2)
                 !write(60+myid, '(I1.1, 2(I2.1, A1, F18.7))') 0, 1,':',swe_catch/100.,  &
                 !                                                2,':',slwc_catch*1000.
              case (3)
                 write(60+myid, '(I1.1, 3(I2.1, A1, F18.7))') 0, 1,':',swe_obsgrid(i,m)*10., &
                                                                 2,':',slwc_obsgrid(i,m)*1.,3,':',tp1_obsgrid(i,m)*0.01
              case (4)
                 write(60+myid, '(I1.1, 3(I2.1, A1, F18.7))') 0, 1,':',swe_obsgrid(i,m)*10., &
                                                                 2,':',slwc_obsgrid(i,m)*1.,3,':',tp1_obsgrid(i,m)*0.01, &
                                                                 4,':',tsurf_obsgrid(i,m)*0.001
           end select
    
        ! Predict TB for each SVM grid cells
           
           if (debug) then
              write(LIS_logunit,*) 'SWE is greater than zero; swe=', swe_obsgrid(i,m)
           endif
 
           call output_svmTb(year,dofyr,case_directory,svm_directory,num_of_inputs,       &
                             training_period,training_target,                             &
                             svm_domain, scaling_option,                                  &
                             forest_decoup_option,                                        &
                             atm_decoup_option, myid, obsgrid_i,svm_i,                    &
                             predict_10H_36H, predict_10V_36V,                            &
                             predict_18H_36H, predict_18V_36V,                            &
                             predict_10H,                                                 &
                             predict_10V,                                                 &
                             predict_18H,                                                 &
                             predict_18V,                                                 &
                             predict_36H,                                                 &
                             predict_36V)

           DTB_10H_36H(i,m) = predict_10H_36H
           DTB_10V_36V(i,m) = predict_10V_36V
           DTB_18H_36H(i,m) = predict_18H_36H
           DTB_18V_36V(i,m) = predict_18V_36V
           TB_10H(i,m)      = predict_10H
           TB_10V(i,m)      = predict_10V
           TB_18H(i,m)      = predict_18H
           TB_18V(i,m)      = predict_18V
           TB_36H(i,m)      = predict_36H
           TB_36V(i,m)      = predict_36V
        else
           DTB_10H_36H(i,m) = -9999
           DTB_10V_36V(i,m) = -9999
           DTB_18H_36H(i,m) = -9999
           DTB_18V_36V(i,m) = -9999
           TB_10H(i,m)      = -9999
           TB_10V(i,m)      = -9999
           TB_18H(i,m)      = -9999
           TB_18V(i,m)      = -9999
           TB_36H(i,m)      = -9999
           TB_36V(i,m)      = -9999
        endif
        close(60+myid, status='delete')

        if (debug) then
           if (DTB_10V_36V(i,m) > -9999 ) then
              write(LIS_logunit,*) 'DTB_10H_36H=', DTB_10H_36H(i,m)
              write(LIS_logunit,*) 'DTB_10V_36V=', DTB_10V_36V(i,m)
              write(LIS_logunit,*) 'DTB_18H_36H=', DTB_18H_36H(i,m)
              write(LIS_logunit,*) 'DTB_18V_36V=', DTB_18V_36V(i,m)
              write(LIS_logunit,*) 'TB_10H=', TB_10H(i,m)
              write(LIS_logunit,*) 'TB_10V=', TB_10V(i,m)
              write(LIS_logunit,*) 'TB_18H=', TB_18H(i,m)
              write(LIS_logunit,*) 'TB_18V=', TB_18V(i,m)
              write(LIS_logunit,*) 'TB_36H=', TB_36H(i,m)
              write(LIS_logunit,*) 'TB_36V=', TB_36V(i,m)
           endif
        endif

     enddo   !m=1,LIS_rc%nensem(n) 
  enddo   !i=1,LIS_rc%obs_ngrid(k)

end subroutine noahmp36_LIS_SVM
