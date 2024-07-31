!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LVT_misc.h"
!BOP
! 
! !ROUTINE: readGLERLObs
! \label{readGLERLObs}
!
! !INTERFACE: 
subroutine readGLERLObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc, LVT_domain
  use LVT_timeMgrMod,  only : LVT_calendar
  use LVT_logMod
  use LVT_histDataMod
  use GLERL_dataMod,   only : GLERLObs
  use map_utils

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The GLERL output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  character*100           :: filename
  logical                 :: file_exists
  integer                 :: nid, ios
  integer                 :: stn_row, stn_col
  real                    :: col, row
  integer                 :: nc,nr
  real                    :: qle(12),watert(12)
  character*4             :: fyr
  integer                 :: ftn,c,r,t,kk,k
  type(ESMF_Time)         :: currTime
  type(ESMF_TimeInterval) :: ts
  integer                 :: cyr, cmo, cda, chr, cmn, css
  integer                 :: status
  real                    :: qle_2d(LVT_rc%lnc,LVT_rc%lnr)
  real                    :: watert_2d(LVT_rc%lnc,LVT_rc%lnr)


  if(GLERLobs(source)%yr.ne.LVT_rc%dyr(source).or.&
       LVT_rc%resetFlag(source)) then      

     LVT_rc%resetFlag(source) = .false. 

     GLERLobs(source)%yr = LVT_rc%dyr(source)

     GLERLobs(source)%qle = LVT_rc%udef
     GLERLobs(source)%watert = LVT_rc%udef

     write(fyr, '(i4.4)') LVT_rc%dyr(source)
     do k=1,GLERLobs(source)%nlocs
        filename = trim(GLERLobs(source)%odir)//'/Lake'//&
             trim(GLERLobs(source)%lake_locname(k))//'_Evap_'//&
             trim(fyr)//'.dat'
        
        inquire(file=trim(filename),exist=file_exists) 
        
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading GLERL Evap file ',trim(filename)
           ftn = LVT_getNextUnitNumber()
           
           open(ftn,file=trim(filename), form='formatted')
           read(ftn,*) qle
           call LVT_releaseUnitNumber(ftn)
           
           GLERLobs(source)%qle(k,:) = qle
           
        endif
     enddo
     do k=1,GLERLobs(source)%nlocs
        filename = trim(GLERLobs(source)%odir)//'/Lake'//&
             trim(GLERLobs(source)%lake_locname(k))//'_WaterTemps_'//&
             trim(fyr)//'.dat'
        
        inquire(file=trim(filename),exist=file_exists) 
        
        if(file_exists) then 
           write(LVT_logunit,*) '[INFO] Reading GLERL Evap file ',trim(filename)
           ftn = LVT_getNextUnitNumber()
           
           open(ftn,file=trim(filename), form='formatted')
           read(ftn,*) watert
           call LVT_releaseUnitNumber(ftn)
           
           GLERLobs(source)%watert(k,:) = watert
        endif
     enddo
  endif
  !log the data only at the change of a month. 
  if(LVT_rc%d_nmo(source).ne.GLERLobs(source)%mo) then 
     qle_2d     = LVT_rc%udef
     watert_2d  = LVT_rc%udef
     cmo = GLERLobs(source)%mo
     write(LVT_logunit,*) '[INFO] Reading GLERL data for ',cmo
     
     do k=1,GLERLobs(source)%nlocs
        call latlon_to_ij(LVT_domain%lvtproj,GLERLobs(source)%lake_lat(k),&
             GLERLobs(source)%lake_lon(k),col,row)
        stn_col = nint(col)
        stn_row = nint(row)
        
        qle_2d(stn_col, stn_row) = GLERLobs(source)%qle(k,cmo)
        watert_2d(stn_col, stn_row) = GLERLobs(source)%watert(k,cmo) + 273.15
     enddo

     GLERLobs(source)%mo = LVT_rc%d_nmo(source)

  else
     qle_2d  = LVT_rc%udef
     watert_2d  = LVT_rc%udef
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source, qle_2d,vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_AVGSURFT,source, watert_2d,vlevel=1,units="K")

end subroutine readGLERLObs
