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
! !ROUTINE: readMODISsportLAIObs
! \label{readMODISsportLAIObs}
!
! !INTERFACE: 
subroutine readMODISsportLAIObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_logMod
  use LVT_histDataMod
  use MODISsportLAIobsMod
          
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
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  7 Mar 2015: Sujay Kumar, Initial Specification
! 
!EOP

  integer                :: c,r, tindex
  integer                :: flag
  integer                :: ftn
  character*100          :: fname
  logical                :: file_exists
  integer                :: iret  
  integer                :: nc, nr
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)

  varfield = LVT_rc%udef

  nc = MODISsportLAIobs(source)%nc
  nr = MODISsportLAIobs(source)%nr

  LVT_rc%resetFlag(source) = .false. 
  MODISsportLAIobs(source)%startFlag = .false. 
  
  call create_MODISsportLAI_filename(MODISsportLAIobs(Source)%odir, &
       LVT_rc%dyr(source),LVT_rc%ddoy(source),&
       fname)
  
  inquire(file=trim(fname),exist=file_exists) 
  
  ! Check if both files exist:
  if( file_exists ) then 
     write(LVT_logunit,*) '[INFO] Reading MODIS SPORT LAI file ',trim(fname)
     
     call read_MODISsport_LAI_data(source, fname, varfield)
     
  else
     write(LVT_logunit,*)'[WARN] MODIS SPORT LAI file missing: ',trim(fname)
     varfield  = LVT_rc%udef
  endif
  
  call LVT_logSingleDataStreamVar(LVT_MOC_LAI,source,varfield,&
       vlevel=1,units="-")
  
end subroutine readMODISsportLAIObs

!BOP
! 
! !ROUTINE: read_MODISsport_LAI_data
! \label{read_MODISsport_LAI_data}
!
! !INTERFACE:
subroutine read_MODISsport_LAI_data(source, fname, laiobs_ip)
! 
! !USES:   

  use LVT_coreMod
  use LVT_logMod
  use LVT_timeMgrMod
  use MODISsportLAIobsMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
  integer                       :: source
  character (len=*)             :: fname
  real                          :: laiobs_ip(LVT_rc%lnc*LVT_rc%lnr)

! !OUTPUT PARAMETERS:
!
!
! !DESCRIPTION: 
!  This subroutine reads the GLASS LAI file and applies the data
!  quality flags to filter the data. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]            index of the nest
!  \item[fname]        name of the RTGLASS AMSR-E file
!  \item[laiobs\_ip]    soil moisture data processed to the LIS domain
! \end{description}
!
!
!EOP
  integer                 :: ftn
  integer                 :: c,r,iret
  real                    :: lai_in(MODISsportLAIobs(source)%nc*MODISsportLAIobs(source)%nr)
  logical*1               :: lai_data_b(MODISsportLAIobs(source)%nc*MODISsportLAIobs(source)%nr)
  logical*1               :: laiobs_b_ip(LVT_rc%lnc*LVT_rc%lnr)


  ftn = LVT_getNextUnitNumber()
  open(ftn,file=trim(fname),form='unformatted',access='direct',&
       convert="little_endian",&
       recl=MODISsportLAIobs(source)%nc*MODISsportLAIobs(source)%nc*4)
  read(ftn,rec=1) lai_in
  call LVT_releaseUnitNumber(ftn)

  lai_data_b = .false. 

  do r=1,MODISsportLAIobs(source)%nr
     do c=1,MODISsportLAIobs(source)%nc
        if(lai_in(c+(r-1)*MODISsportLAIobs(source)%nc).gt.0) then 
           lai_data_b(c+(r-1)*MODISsportLAIobs(source)%nc) =  .true. 
        else
           lai_data_b(c+(r-1)*MODISsportLAIobs(source)%nc) = .false. 
        endif
     enddo
  enddo

  if(LVT_isAtAfinerResolution(MODISsportLAIobs(source)%datares)) then
!--------------------------------------------------------------------------
! Interpolate to the LVT analysis domain
!-------------------------------------------------------------------------- 
     call bilinear_interp(LVT_rc%gridDesc,&
          lai_data_b, lai_in, laiobs_b_ip, laiobs_ip, &
          MODISsportLAIobs(source)%nc*MODISsportLAIobs(source)%nr, &
          LVT_rc%lnc*LVT_rc%lnr, &
          MODISsportLAIobs(source)%rlat,MODISsportLAIobs(source)%rlon,&
          MODISsportLAIobs(source)%w11,MODISsportLAIobs(source)%w12,&
          MODISsportLAIobs(source)%w21,MODISsportLAIobs(source)%w22,&
          MODISsportLAIobs(source)%n11,MODISsportLAIobs(source)%n12,&
          MODISsportLAIobs(source)%n21,MODISsportLAIobs(source)%n22,LVT_rc%udef,iret)

  else
     call upscaleByAveraging(MODISsportLAIobs(source)%nc*MODISsportLAIobs(source)%nr,&
          LVT_rc%lnc*LVT_rc%lnr, &
          LVT_rc%udef, MODISsportLAIobs(source)%n11,&
          lai_data_b,lai_in, laiobs_b_ip, laiobs_ip)
  endif

end subroutine read_MODISsport_LAI_data


!BOP
! !ROUTINE: create_MODISsportLAI_filename
! \label{create_MODISsportLAI_filename}
! 
! !INTERFACE: 
subroutine create_MODISsportLAI_filename(ndir,  yr, doy, filename)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filename
  integer           :: yr, doy
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the GLASS LAI filename
!  based on the time and date 
! 
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the GLASS LAI data directory
!  \item[yr]  current year
!  \item[mo]  current doy
!  \item[filename] Generated GLASS LAI filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=3) :: fdoy
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fdoy, fmt='(i3.3)') doy

  filename = trim(ndir)//'/'//trim(fyr)//'/MLAI_FLT_'//&
       trim(fyr)//trim(fdoy)//'.dat'

end subroutine create_MODISsportLAI_filename


