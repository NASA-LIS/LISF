!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readAGRMETdata
! \label{readAGRMETdata}
!
! !INTERFACE: 
subroutine readAGRMETdata(source)
! 
! !USES:   
  use grib_api
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use LVT_histDataMod
  use AGRMET_dataMod

  implicit none
  
!
! !INPUT PARAMETERS: 
  integer, intent(in)      :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  
!  This subroutine provides the data reader for AGRMET output data. 
!  The plugin processes the land surface model output variables of
!  surface fluxes, snow and model forcings. The routine also spatially
!  interpolates the AGRMET data to the LVT target grid 
!  and resolution. 
! 
!  NOTES: 
!   The AGRMET output is available at 3-hourly intervals. So 
!   the comparisons against model data should use at least a 
!   3 hour averaging interval. The fluxes (qle, qh, qg) and 
!   snow depth are not available before 2002/05/29
! 
!   The arguments are: 
!   \begin{description}
!    \item[source] source of the data stream (1 or 2)
!   \end{description}
!
! !REVISION HISTORY: 
!  09 Dec 2010: Sujay Kumar, Initial Specification
!  09 Nov 2017: Eric Kemp (SSAI).  Refactored code to only read GRIB fields 
!               that LVT can actually use.  Also reduced memory footprint,
!               replaced one dynamically allocatable array with an automatic
!               array, and replaced some 1-d arrays with scalar integers.
!               This greatly reduces run-time.
!  13 Nov 2017: Eric Kemp (SSAI).  Added specific humidity forcing.
!  02 Nov 2018: Eric Kemp (SSAI).  Added support for n1280e domain.
! 
!EOP

  !integer,   parameter         :: nc = 1440, nr=600
  integer :: nc,nr
  character*200                :: filename
  logical                      :: file_exists
  integer                      :: nunit,ufn,iret,ierr
  integer                      :: c,r
  integer                      :: ftn
  integer                      :: count1
  real                         :: &
       swd(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nswd(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       tair(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       ntair(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       rainf(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nrainf(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       lwd(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nlwd(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       qle(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nqle(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       qh(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nqh(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       qg(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nqg(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       swe(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nswe(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       snod(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nsnod(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       sm1(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nsm1(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       sm2(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nsm2(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       sm3(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nsm3(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       sm4(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nsm4(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       st1(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nst1(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       st2(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nst2(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       st3(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nst3(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       st4(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: &
       nst4(agrmetdata(source)%nc*agrmetdata(source)%nr)
  real                         :: &
       qair(agrmetdata(source)%nc*agrmetdata(source)%nr) ! EMK Specific humidity
  integer                      :: &
       nqair(agrmetdata(source)%nc*agrmetdata(source)%nr)

  logical*1                    :: &
       lb(agrmetdata(source)%nc*agrmetdata(source)%nr)
  integer                      :: swd_topt, lwd_topt,rainf_topt,tair_topt
  integer                      :: qle_topt, qh_topt, qg_topt
  integer                      :: swe_topt, snod_topt
  integer                      :: sm_topt,st_topt
  integer                      :: qair_topt ! EMK

  integer                      :: swd_index, lwd_index,rainf_index, tair_index
  integer                      :: qle_index, qh_index, qg_index
  integer                      :: swe_index, snod_index
  integer                      :: sm_index, st_index
  integer                      :: qair_index ! EMK

  real                         :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  integer                      :: yr1, mo1, da1,hr1,mn1,ss1
  integer                      :: yr2, mo2, da2,hr2,mn2,ss2
  type(ESMF_Time)              :: time1
  type(ESMF_Time)              :: time2
  type(ESMF_TimeInterval)      :: lis_ts
  integer                      :: nvars
  real                         :: &
       var(agrmetdata(source)%nc*agrmetdata(source)%nr) ! EMK

  integer                      :: index
  integer                      :: igrib
  integer                      :: pid  ! EMK
  integer                      :: tid  ! EMK
  integer                      :: lid  ! EMK

  integer                      :: status

  nc = agrmetdata(source)%nc
  nr = agrmetdata(source)%nr

!These are hardcoded for now. 1-instantaneous, 7-timeavged
  swd_topt = 7
  lwd_topt = 7
  rainf_topt =133 
  tair_topt = 7 
  qle_topt = 7
  qh_topt = 7
  qg_topt = 7 
  swe_topt = 1
  snod_topt = 1
  sm_topt = 7
  st_topt = 7
  qair_topt = 7 ! EMK

!gribids 
  swd_index = 145
  lwd_index = 144
  rainf_index = 61
  tair_index = 11
  qle_index = 121
  qh_index = 122
  qg_index = 155
  swe_index = 65
  snod_index = 66
  sm_index = 201
  st_index = 85
 ! tskin_index = 148
  qair_index = 51 ! EMK

  yr1 = LVT_rc%dyr(source)
  mo1 = LVT_rc%dmo(source)
  da1 = LVT_rc%dda(source)
  hr1 = LVT_rc%dhr(source)
  mn1 = 0
  ss1 = 0
  
  swd = 0 
  nswd = 0 
  
  lwd = 0 
  nlwd = 0 
  
  rainf = 0 
  nrainf = 0 
  
  tair = 0 
  ntair = 0 
  
  qle = 0 
  nqle = 0 
  
  qh = 0 
  nqh = 0 
  
  qg = 0 
  nqg = 0
  
  swe = 0 
  nswe = 0 
  
  snod = 0 
  nsnod = 0 
  
  sm1 = 0 
  nsm1 = 0 
  sm2 = 0
  nsm2 = 0
  sm3 = 0 
  nsm3 = 0 
  sm4 = 0
  nsm4 = 0 
  
  st1 = 0 
  nst1 = 0 
  st2 = 0
  nst2 = 0
  st3 = 0 
  nst3 = 0 
  st4 = 0
  nst4 = 0 

  qair = 0   ! EMK
  nqair = 0  ! EMK

  call ESMF_TimeSet(time1,yy=yr1, mm=mo1, dd=da1, &
       h=hr1,m=mn1,s=ss1,calendar=LVT_calendar, rc=status)
  call LVT_verify(status)  

  call ESMF_TimeIntervalSet(lis_ts, s = LVT_rc%ts, &
       rc=status)
  call LVT_verify(status)  
  
!  time2 = time1 -lis_ts
  time2 = time1

!read files from time1 to time2
  do while (time2.le.time1) 

     call ESMF_TimeGet(time2,yy=yr2, mm=mo2, dd=da2, &
          h=hr2,m=mn2,s=ss2,calendar=LVT_calendar, rc=status)
     call LVT_verify(status)  

     call create_agrmetdata_filename(source, &
          yr2, mo2, da2, hr2, filename)

     inquire(file=trim(filename),exist=file_exists)
     
     !EMK
     if (.not. file_exists) then
        write(LVT_logunit,*) '[INFO] Cannot find file : ',trim(filename)
     end if

     if(file_exists) then 
        call grib_open_file(ftn,trim(filename),'r',iret)
        if(iret.ne.0) then 
           write(LVT_logunit,*) '[ERR] Stopping routine: File not opened : ',iret
           call LVT_endrun
        endif
        write(LVT_logunit,*) '[INFO] Reading AGRMET ',trim(filename)
        call grib_count_in_file(ftn,nvars)
        
        ! EMK...Initialize before use!
        pid = 0
        tid = 0
        lid = 0
        var(:) = 0

        do index=1,nvars
           call grib_new_from_file(ftn,igrib,iret)
           call LVT_verify(iret,'grib_new_from_file failed in readAgrmetdata')
           
           call grib_get(igrib,"indicatorOfParameter",pid,iret)
           call LVT_verify(iret, & 
                'grib_get failed for indicatorOfParameter in readAgrmetdata')
           
           call grib_get(igrib, "timeRangeIndicator",tid,iret)
           call LVT_verify(iret,&
                'grib_get failed for timeRangeIndicator in readAgrmetdata')
           
           call grib_get(igrib, "bottomLevel",lid,iret)
           call LVT_verify(iret,&
                'grib_get failed for bottomLevel in readAgrmetdata')
           
           ! EMK...Only read the data field if LVT can use it.
           if(pid.eq.swd_index.and.tid.eq.swd_topt) then 
              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       swd(c+(r-1)*nc) = swd(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nswd(c+(r-1)*nc) = nswd(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.lwd_index.and.tid.eq.lwd_topt) then 
              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       lwd(c+(r-1)*nc) = lwd(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nlwd(c+(r-1)*nc) = nlwd(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.rainf_index.and.tid.eq.rainf_topt) then 
              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       rainf(c+(r-1)*nc) = rainf(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nrainf(c+(r-1)*nc) = nrainf(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.tair_index.and.tid.eq.tair_topt) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       tair(c+(r-1)*nc) = tair(c+(r-1)*nc)+var(c+(r-1)*nc)
                       ntair(c+(r-1)*nc) = ntair(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.qle_index.and.tid.eq.qle_topt) then 

             call grib_get(igrib,"values",var,iret)
             call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       qle(c+(r-1)*nc) = qle(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nqle(c+(r-1)*nc) = nqle(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.qh_index.and.tid.eq.qh_topt) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       qh(c+(r-1)*nc) = qh(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nqh(c+(r-1)*nc) = nqh(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.qg_index.and.tid.eq.qg_topt) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       qg(c+(r-1)*nc) = qg(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nqg(c+(r-1)*nc) = nqg(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.swe_index.and.tid.eq.swe_topt) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       swe(c+(r-1)*nc) = swe(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nswe(c+(r-1)*nc) = nswe(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo
           elseif(pid.eq.snod_index.and.tid.eq.snod_topt) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       snod(c+(r-1)*nc) = snod(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nsnod(c+(r-1)*nc) = nsnod(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.sm_index.and.tid.eq.sm_topt.and.&
                lid.eq.10) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       sm1(c+(r-1)*nc) = sm1(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nsm1(c+(r-1)*nc) = nsm1(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.sm_index.and.tid.eq.sm_topt.and.&
                lid.eq.40) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       sm2(c+(r-1)*nc) = sm2(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nsm2(c+(r-1)*nc) = nsm2(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.sm_index.and.tid.eq.sm_topt.and.&
                lid.eq.100) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       sm3(c+(r-1)*nc) = sm3(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nsm3(c+(r-1)*nc) = nsm3(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.sm_index.and.tid.eq.sm_topt.and.&
                lid.eq.200) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       sm4(c+(r-1)*nc) = sm4(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nsm4(c+(r-1)*nc) = nsm4(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.st_index.and.tid.eq.st_topt.and.&
                lid.eq.10) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       st1(c+(r-1)*nc) = st1(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nst1(c+(r-1)*nc) = nst1(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.st_index.and.tid.eq.st_topt.and.&
                lid.eq.40) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       st2(c+(r-1)*nc) = st2(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nst2(c+(r-1)*nc) = nst2(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.st_index.and.tid.eq.st_topt.and.&
                lid.eq.100) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       st3(c+(r-1)*nc) = st3(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nst3(c+(r-1)*nc) = nst3(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           elseif(pid.eq.st_index.and.tid.eq.st_topt.and.&
                lid.eq.200) then 

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       st4(c+(r-1)*nc) = st4(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nst4(c+(r-1)*nc) = nst4(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           ! EMK...Support specific humidity
           elseif(pid.eq.qair_index.and.tid.eq.qair_topt) then

              call grib_get(igrib,"values",var,iret)
              call LVT_verify(iret,'grib_get failed for values in readAgrmetdata')
              do r=1,nr
                 do c=1,nc
                    if(var(c+(r-1)*nc).ne.9999.0) then 
                       qair(c+(r-1)*nc) = qair(c+(r-1)*nc)+var(c+(r-1)*nc)
                       nqair(c+(r-1)*nc) = nqair(c+(r-1)*nc)+1
                    endif
                 enddo
              enddo

           endif

           call grib_release(igrib,iret)
           call LVT_verify(iret,'grib_release failed in readAgrmetdata')

        enddo

        call grib_close_file(ftn,iret)
        
        write(LVT_logunit,*) '[INFO] Successfully processed ',trim(filename)
     endif
     time2 = time2 + agrmetdata(source)%ts
  end do

  call interp_agrmetvar(source, nc,nr, swd, nswd, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_swdownforc,source,varfield,vlevel=1,units="W/m2")

  call interp_agrmetvar(source,nc,nr, lwd, nlwd, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_lwdownforc,source,varfield,vlevel=1,units="W/m2")

  call interp_agrmetvar(source,nc,nr, rainf, nrainf, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_totalprecip,source,varfield,vlevel=1,units="kg/m2")

  call interp_agrmetvar(source,nc,nr, tair, ntair, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_tairforc,source,varfield,vlevel=1,units="K")

  call interp_agrmetvar(source,nc,nr, qle, nqle, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_qle,source,varfield,vlevel=1,units="W/m2")

  call interp_agrmetvar(source,nc,nr, qh, nqh, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_qh,source,varfield,vlevel=1,units="W/m2")

  call interp_agrmetvar(source,nc,nr, qg, nqg, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_qg,source,varfield,vlevel=1,units="W/m2",dir='DN')

  call interp_agrmetvar(source,nc,nr, swe, nswe, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_swe,source,varfield,vlevel=1,units="kg/m2")

  call interp_agrmetvar(source,nc,nr, snod, nsnod, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth,source,varfield,vlevel=1,units="m")

  call interp_agrmetvar(source,nc,nr, sm1, nsm1, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,varfield,units="m3/m3",&
       vlevel = 1)

  call interp_agrmetvar(source,nc,nr, sm2, nsm2, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,varfield,units="m3/m3",&
       vlevel = 2)

  call interp_agrmetvar(source,nc,nr, sm3, nsm3, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,varfield,units="m3/m3",&
       vlevel = 3)

  call interp_agrmetvar(source,nc,nr, sm4, nsm4, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soilmoist,source,varfield,units="m3/m3",&
       vlevel = 4)

  call interp_agrmetvar(source,nc,nr, st1, nst1, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp,source,varfield,units="K",&
       vlevel = 1)

  call interp_agrmetvar(source,nc,nr, st2, nst2, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp,source,varfield,units="K",&
       vlevel = 2)

  call interp_agrmetvar(source,nc,nr, st3, nst3, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp,source,varfield,units="K",&
       vlevel = 3)

  call interp_agrmetvar(source,nc,nr, st4, nst4, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_soiltemp,source,varfield,units="K",&
       vlevel = 4)

  ! EMK...Support specific humidity
  call interp_agrmetvar(source,nc,nr, qair, nqair, varfield)
  call LVT_logSingleDataStreamVar(LVT_MOC_qairforc,source,varfield,vlevel=1,&
       units="kg/kg")

end subroutine readAGRMETdata

!BOP
! 
! !ROUTINE: interp_agrmetvar
!  \label{interp_agrmetvar}
!
! !INTERFACE: 
subroutine interp_agrmetvar(source,nc, nr, var_input, nvar_input, var_output)
! 
! !USES:   
  use LVT_logMod, only: LVT_logunit
  use LVT_coreMod,   only : LVT_rc
  use AGRMET_dataMod, only : agrmetdata

  implicit none
!
! !ARGUMENTS: 
! 
  integer            :: source
  integer            :: nc
  integer            :: nr
  real               :: var_input(nc*nr)
  integer            :: nvar_input(nc*nr)
  logical*1          :: lb(nc*nr)
  real               :: var_output(LVT_rc%lnc, LVT_rc%lnr)
!
! !DESCRIPTION: 
! 
!   This subroutine spatially interpolates the AGRMET variable to the
!   LVT output grid and resolution, using a bilinear interpolation
!   approach. 
! 
!   The arguments are: 
!   \begin{description}
!    \item[source] source of the data stream (1 or 2)
!    \item[nc]      number of columns in the input (AGRMET) grid
!    \item[nr]      number of rows in the input (AGRMET) grid
!    \item[var_input] input variable to be interpolated
!    \item[lb]        input bitmap (true//false)
!    \item[var_output] resulting interplated field
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  integer            :: iret
  integer            :: c,r
  logical*1          :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real               :: go(LVT_rc%lnc*LVT_rc%lnr)

  var_output = LVT_rc%udef
  lb = .false. 
  do r=1,nr
     do c=1,nc
        if(nvar_input(c+(r-1)*nc).gt.0) then 
           var_input(c+(r-1)*nc) = &
                var_input(c+(r-1)*nc)/&
                nvar_input(c+(r-1)*nc)
           lb(c+(r-1)*nc) = .true.
        else
           var_input(c+(r-1)*nc) = LVT_rc%udef
        endif
     enddo
  enddo

  call bilinear_interp(LVT_rc%gridDesc,lb, var_input,&
       lo,go,nc*nr,LVT_rc%lnc*LVT_rc%lnr,&
       agrmetdata(source)%rlat,agrmetdata(source)%rlon,&
       agrmetdata(source)%w11,agrmetdata(source)%w12,&
       agrmetdata(source)%w21,agrmetdata(source)%w22,&
       agrmetdata(source)%n11,agrmetdata(source)%n12,&
       agrmetdata(source)%n21,agrmetdata(source)%n22,LVT_rc%udef,iret)
  
  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        var_output(c,r) = go(c+(r-1)*LVT_rc%lnc)
     enddo
  enddo

end subroutine interp_agrmetvar
!BOP
! 
! !ROUTINE: create_agrmetdata_filename
! \label{create_agrmetdata_filename}
!
! !INTERFACE: 
subroutine create_agrmetdata_filename(source, yr, mo, da, hr, filename)
! !USES:   
  use AGRMET_dataMod
  use LVT_logMod, only: LVT_endrun, LVT_logunit

  implicit none

! !INPUT PARAMETERS: 
! 
  integer,    intent(in)    :: source
  integer           :: yr
  integer           :: mo
  integer           :: da
  integer           :: hr
  character(len=*)  :: filename

!
! !DESCRIPTION: 
!  This subroutine creates the AGRMET filename based on the given 
!  date (year, month, day and hour)
! 
!  The arguments are: 
!  \begin{description}
!   \item[source]      source of the data stream (1/2)
!   \item[yr]        year of data
!   \item[mo]        month of data
!   \item[da]        day of data
!   \item[hr]        hour of data
!   \item[filename]  Name of the AGRMET file
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

  character*4       :: fyr
  character*2       :: fmo
  character*2       :: fda
  character*2       :: fhr

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr

  if (trim(agrmetdata(source)%gridname) == "GLOBAL") then
     ! Old 0.25 deg deterministic run
     filename = trim(agrmetdata(source)%odir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//'/'&
       //'PS.AFWA_SC.'//&
       trim(agrmetdata(source)%security_class)//'_DI.'//&
       trim(agrmetdata(source)%distribution_class)//&
       '_DC.'//trim(agrmetdata(source)%data_category)//'_GP.LIS_GR'//&
       '.C0P25DEG_AR.'//trim(agrmetdata(source)%area_of_data)//&
       '_PA.03-HR-SUM_DD.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'_DT.'//trim(fhr)//'00_DF.GR1'
  else if (trim(agrmetdata(source)%gridname) == "n1280e") then
     ! In-house Bratseth run matching GALWEM n1280e domain
     filename = trim(agrmetdata(source)%odir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//'/'&
       //'PS.AFWA_SC.'//&
       trim(agrmetdata(source)%security_class)//'_DI.'//&
       trim(agrmetdata(source)%distribution_class)//&
       '_DC.'//trim(agrmetdata(source)%data_category)//'_GP.LIS_GR'//&
       '.C0P09DEG_AR.'//trim(agrmetdata(source)%area_of_data)//&
       '_PA.03-HR-SUM_DD.'//&
       trim(fyr)//trim(fmo)//trim(fda)//'_DT.'//trim(fhr)//'00_DF.GR1'
       
  else
     write(LVT_logunit,*) &
          '[ERR] Internal error, unknown AGRMET data gridname ', &
          trim(agrmetdata(source)%gridname)
     call LVT_endrun()
  end if
end subroutine create_agrmetdata_filename
  
