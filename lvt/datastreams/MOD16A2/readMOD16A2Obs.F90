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
! !ROUTINE: readMOD16A2Obs
! \label{readMOD16A2Obs}
!
! !INTERFACE: 
subroutine readMOD16A2Obs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_logMod,      only : LVT_logunit, LVT_verify, &
       LVT_getNextUnitNumber, LVT_releaseUnitNumber
  use LVT_histDataMod
  use MOD16A2_obsMod, only : MOD16A2Obs

  implicit none
!
! !INPUT PARAMETERS: 
  integer,     intent(in)    :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The MOD16A2 output is available at monthly intervals. So 
!   the comparisons against model data should use at least a 
!   24 hour (1day) averaging interval. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  10 Dec 2010: Sujay Kumar, Initial Specification
! 
!EOP

  character*100          :: filename
  integer                :: ftn 
  logical                :: file_exists
  integer                :: nid, ios
  integer                :: c1,r1,line
  real    :: qle(mod16a2obs(source)%nc,mod16a2obs(source)%nr)
  real    :: qle1d(mod16a2obs(source)%nc*mod16a2obs(source)%nr)
  logical*1 :: li(mod16a2obs(source)%nc*mod16a2obs(source)%nr)
  real                   :: lat1,lon1
  integer                :: c,r,t,kk
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: gridDesc(6)

  if(Mod16a2obs(Source)%mo.ne.LVT_rc%d_nmo(source).or.&
       LVT_rc%resetFlag(source)) then

     LVT_rc%resetFlag(source) = .false. 
!     if(mod16a2obs(source)%startFlag) then 
!        if(LVT_rc%o_nmo.ne.1) then        
!           Mod16a2obs(Source)%yr = LVT_rc%o_nyr
!           Mod16a2obs(Source)%mo = LVT_rc%o_nmo-1
!        else
!           Mod16a2obs(Source)%yr = LVT_rc%o_nyr-1
!           Mod16a2obs(Source)%mo = 12
!        endif

!        mod16a2obs(source)%startflag = .false.
!     endif
!enable this for runmode = 0
     if(mod16a2obs(source)%startFlag) then 
        Mod16a2obs(Source)%yr = LVT_rc%dyr(source)
        Mod16a2obs(Source)%mo = LVT_rc%dmo(source)
        mod16a2obs(source)%startflag = .false.
     endif

     Mod16a2obs(Source)%qle = LVT_rc%udef

     call create_mod16a2_filename(Mod16a2obs(Source)%odir, Mod16a2obs(Source)%yr,&
          Mod16a2obs(Source)%mo,filename)

     inquire(file=trim(filename),exist=file_exists) 
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading MOD16A2 LH file ',trim(filename)
        
        gridDesc = 0 
        gridDesc(1) = mod16a2obs(source)%gridDesc(4)
        gridDesc(2) = mod16a2obs(source)%gridDesc(5)
        gridDesc(3) = mod16a2obs(source)%gridDesc(7)
        gridDesc(4) = mod16a2obs(source)%gridDesc(8)
        gridDesc(5) = 0.01
        gridDesc(6) = 0.01

        ftn = LVT_getNextUnitNumber()
        open(ftn,file=trim(filename),form='unformatted',recl=4,&
             access='direct',iostat=ios)
        do r=1,mod16a2obs(source)%nr
           do c=1,mod16a2obs(source)%nc
              lat1 = mod16a2obs(source)%gridDesc(4)+(r-1)*0.01
              lon1 = mod16a2obs(source)%gridDesc(5)+(c-1)*0.01
              r1 = nint((lat1+59.995)/0.01)+1
              c1 = nint((lon1+179.995)/0.01)+1
              line = c1+(r1-1)*36000
              read(ftn,rec=line) qle(c,r)
           enddo
        enddo

        call LVT_releaseUnitNumber(ftn)
        
        li = .false. 
!values
        do r=1,mod16a2obs(source)%nr
           do c=1,mod16a2obs(source)%nc
              if(qle(c,r).ne.-9999.0) then 
! Divide by 86400 (number of seconds in a day)
! to convert units from J/m^s/day to W/m^2.
                 qle1d(c+(r-1)*mod16a2obs(source)%nc) = qle(c,r) / 86400.0
                 li(c+(r-1)*mod16a2obs(source)%nc) = .true. 
              else
                 qle1d(c+(r-1)*mod16a2obs(source)%nc) = -9999.0
                 li(c+(r-1)*mod16a2obs(source)%nc) = .false. 
              endif
           enddo
        enddo
        
        call upscaleByAveraging(mod16a2obs(source)%nc*mod16a2obs(source)%nr,&
             LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef,&
             Mod16a2obs(Source)%n11,li,qle1d,lo,Mod16a2obs(Source)%qle)
     endif

     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           if(lo(c+(r-1)*LVT_rc%lnc)) then 
              varfield(c,r) = Mod16a2obs(Source)%qle(c+(r-1)*LVT_rc%lnc)          
           else
              varfield(c,r) = LVT_rc%udef
           endif
        enddo
     enddo
     
     Mod16a2obs(Source)%yr = LVT_rc%d_nyr(source)
     Mod16a2obs(Source)%mo = LVT_rc%d_nmo(source)
  else
     varfield  = LVT_rc%udef
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,varfield,vlevel=1,units="W/m2")
  
end subroutine readMOD16A2Obs

!BOP
! 
! !ROUTINE: create_mod16a2_filename
! \label{create_mod16a2_filename}
!
! !INTERFACE: 
subroutine create_mod16a2_filename(odir,yr,mo,filename)
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
! This routine creates a timestamped filename for MOD16A2_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the MOD16A2_LH file
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
  
  filename = trim(odir)//'/Y'//trim(fyr)//'/M'//trim(fmo)//&
       '/LE.1gd4r'
  
end subroutine create_mod16a2_filename


