!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! 
! !ROUTINE: readNatSFObs
! \label{readNatSFObs}
!
! !INTERFACE: 
  subroutine readNatSFObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,        only : LVT_rc, LVT_domain
  use LVT_histDataMod
  use LVT_timeMgrMod,     only : LVT_calendar, LVT_tick
  use LVT_logMod,         only : LVT_verify, LVT_getNextUnitNumber, &
       LVT_releaseUnitNumber, LVT_logunit
  use NatSF_obsMod, only : NatSFobs
  use map_utils

  implicit none


!
! !INPUT PARAMETERS: 
! 
  integer, intent(in)    :: source
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  11 May 2011: Sujay Kumar, Initial Specification
! 
!EOP
!--------------------------------------------------------------------------

  character*100       :: filename
  character*20        :: datestring, qstring
  integer             :: stn_col, stn_row, c,r
  real                :: col,row
  integer             :: tind
  real                :: offset
  real                :: discharge
  integer             :: i,t,kk
  integer             :: ftn
  character*4         :: fyr
  character*2         :: fmo
  integer             :: iloc
  integer             :: ios,ios1,ios2,status
  integer             :: yr, mo, da
  character*100       :: line
  logical             :: file_exists
  type(ESMF_Time)     :: obsTime,obsTime1
  real                :: q(LVT_rc%lnc,LVT_rc%lnr)
  integer             :: nq(LVT_rc%lnc,LVT_rc%lnr)

  q = 0
  nq = 0 
!every new year read the data, for each station and store it in memory

  if(NatSFobs(source)%mo.ne.LVT_rc%d_nmo(source) &
       .or.NatSFobs(source)%startFlag.or.&
       LVT_rc%resetFlag(source)) then

     LVT_rc%resetFlag(source) = .false. 
     NatSFobs(source)%q  = LVT_rc%udef

     if(NatSFobs(source)%startFlag) then 
        if(LVT_rc%d_nmo(source).ne.1) then 
           NatSFobs(source)%yr = LVT_rc%d_nyr(source)
           NatSFobs(source)%mo = LVT_rc%d_nmo(source) -1
        else
           NatSFobs(source)%yr = LVT_rc%d_nyr(source)-1
           NatSFobs(source)%mo = 12
        endif
        NatSFobs(source)%startFlag = .false.
     endif

     do i = 1, NatSFobs(source)%n_stns
        write(fyr,fmt='(i4.4)') LVT_rc%dyr(source)
        write(fmo,fmt='(i2.2)') LVT_rc%dmo(source)

        filename = trim(NatSFobs(source)%odir)//'/'//&
             trim(Natsfobs(source)%stn_name(i))//&
             '.'//trim(fyr)//trim(fmo)//'.dat'

        inquire(file=trim(filename), exist=file_exists) 
        if(file_exists) then            
           write(LVT_logunit,*) '[INFO] reading Nat streamflow file ',trim(filename)
           ftn = LVT_getNextUnitNumber()

           open(ftn,file=trim(filename),form='formatted')
           read(ftn,*,iostat=ios) yr, mo, discharge
           NatSFobs(source)%q(i) = discharge
           call LVT_releaseUnitNumber(ftn)
        end if
     enddo

     q = LVT_rc%udef
     do i=1,NatSFobs(source)%n_stns
        call latlon_to_ij(LVT_domain%lvtproj, &
             NatSFobs(source)%stnlat(i),NatSFobs(source)%stnlon(i),&
             col,row)
        stn_col = nint(col)
        stn_row = nint(row)
        if(stn_col.gt.0.and.stn_row.gt.0.and.&
             stn_col.le.LVT_rc%lnc.and.stn_row.le.LVT_rc%lnr) then 
           if(NatSFobs(source)%q(i).ne.LVT_rc%udef) then 
              q(stn_col, stn_row) = NatSFobs(source)%q(i)
           endif
        endif
     enddo

     NatSFobs(source)%yr = LVT_rc%d_nyr(source)
     NatSFobs(source)%mo = LVT_rc%d_nmo(source)

  else
     q = LVT_rc%udef

  endif
 
  call LVT_logSingleDataStreamVar(LVT_MOC_streamflow,source,&
       q,vlevel=1,units="m3/s")

end subroutine readNatSFObs

