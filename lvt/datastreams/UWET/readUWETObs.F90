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
! !ROUTINE: readUWETObs
! \label{readUWETObs}
!
! !INTERFACE: 
subroutine readUWETObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit, LVT_verify, &
       LVT_getNextUnitNumber, LVT_releaseUnitNumber
  use LVT_histDataMod
  use UWET_obsMod, only : UWETObs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,  intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The UWET output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  19 Jul 2013: Sujay Kumar, Initial Specification
! 
!EOP

  character*100          :: filename
  integer                :: ftn 
  logical                :: file_exists
  integer                :: nid, ios
  integer                :: c1,r1,line
  real    :: qle(uwetobs(source)%nc,uwetobs(source)%nr)
  real    :: qle1d(uwetobs(source)%nc*uwetobs(source)%nr)
  logical*1 :: li(uwetobs(source)%nc*uwetobs(source)%nr)
  real                   :: lat1,lon1
  integer                :: c,r,t,kk
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: gridDesc(6)

  if((Uwetobs(Source)%mo.ne.LVT_rc%d_nmo(source)).or.&
       LVT_rc%resetFlag(source)) then

     LVT_rc%resetFlag(source) = .false. 
     if(uwetobs(source)%startFlag) then 
        Uwetobs(Source)%yr = LVT_rc%dyr(source)
        Uwetobs(Source)%mo = LVT_rc%dmo(source)
        uwetobs(source)%startflag = .false.
     endif

     Uwetobs(Source)%qle = LVT_rc%udef

     call create_uwet_filename(Uwetobs(Source)%odir, Uwetobs(Source)%yr,&
          Uwetobs(Source)%mo,filename)

     inquire(file=trim(filename),exist=file_exists)
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading UWET LH file ',trim(filename)
        
        gridDesc = 0 
        gridDesc(1) = uwetobs(source)%gridDesc(4)
        gridDesc(2) = uwetobs(source)%gridDesc(5)
        gridDesc(3) = uwetobs(source)%gridDesc(7)
        gridDesc(4) = uwetobs(source)%gridDesc(8)
        gridDesc(5) = 0.05
        gridDesc(6) = 0.05

        ftn = LVT_getNextUnitNumber()
        open(ftn,file=trim(filename),form='unformatted',recl=4,&
             access='direct', convert="little_endian",iostat=ios)
        do r=1,uwetobs(source)%nr
           do c=1,uwetobs(source)%nc
              lat1 = uwetobs(source)%gridDesc(4)+(r-1)*0.05
              lon1 = uwetobs(source)%gridDesc(5)+(c-1)*0.05
              r1 = nint((lat1-25.025)/0.05)+1
              c1 = nint((lon1+124.975)/0.05)+1
              line = c1+(r1-1)*1160
              read(ftn,rec=line) qle(c,r)
           enddo
        enddo

        call LVT_releaseUnitNumber(ftn)
        
        li = .false. 
!values
        do r=1,uwetobs(source)%nr
           do c=1,uwetobs(source)%nc
              if(qle(c,r).ne.-9999.0) then 
                 qle1d(c+(r-1)*uwetobs(source)%nc) = qle(c,r) *2.5E6/(30*24*60*60)
                 li(c+(r-1)*uwetobs(source)%nc) = .true. 
              else
                 qle1d(c+(r-1)*uwetobs(source)%nc) = -9999.0
                 li(c+(r-1)*uwetobs(source)%nc) = .false. 
              endif
           enddo
        enddo
        
        call upscaleByAveraging(uwetobs(source)%nc*uwetobs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef,&
             Uwetobs(Source)%n11,li,qle1d,lo,Uwetobs(Source)%qle)
           
     endif

     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(lo(c+(r-1)*LVT_rc%lnc)) then 
              varfield(c,r) = Uwetobs(Source)%qle(c+(r-1)*LVT_rc%lnc)          
           else
              varfield(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
     
     Uwetobs(Source)%yr = LVT_rc%d_nyr(source)
     Uwetobs(Source)%mo = LVT_rc%d_nmo(source)
  else
     varfield  = LVT_rc%udef
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source, varfield,vlevel=1,units="W/m2")
  
end subroutine readUWETObs

!BOP
! 
! !ROUTINE: create_uwet_filename
! \label{create_uwet_filename}
!
! !INTERFACE: 
subroutine create_uwet_filename(odir,yr,mo,filename)
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION:
! 
! This routine creates a timestamped filename for UWET_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the UWET_LH file
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
  character(len=*)             :: odir
  integer                      :: yr
  integer                      :: mo
  character(len=*)             :: filename
!
!EOP

  character*4             :: fyr
  character*2             :: fmo

  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  
  filename = trim(odir)//'/UW_ET_Monthly_'//trim(fyr)//trim(fmo)//&
       '_0.05.bin'
  
end subroutine create_uwet_filename


