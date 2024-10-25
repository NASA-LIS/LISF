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
! !ROUTINE: readFLUXNETmteObs
! \label{readFLUXNETmteObs}
!
! !INTERFACE: 
subroutine readFLUXNETmteObs(source)
! 
! !USES:   
  use ESMF
  use LVT_coreMod,     only : LVT_rc
  use LVT_timeMgrMod,  only : LVT_calendar
  use LVT_logMod,      only : LVT_logunit, LVT_verify
  use LVT_histDataMod
  use FLUXNETmte_obsMod, only : FLUXNETmteObs

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
!
! !INPUT PARAMETERS: 
  integer, intent(in) :: source
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
!  NOTES: 
!   The FLUXNETmte output is available at monthly intervals. So 
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
  logical                :: file_exists
  integer                :: nid, ios
  integer                :: qleid, qhid,rowId, colId
  integer                :: nc,nr
  real,  allocatable     :: qle(:,:,:)
  real,  allocatable     :: qle1d(:)
  real,  allocatable     :: qh(:,:,:)
  real,  allocatable     :: qh1d(:)
  logical*1, allocatable :: li(:)
  integer                :: c,r,t,kk
  type(ESMF_Time)        :: currTime
  type(ESMF_TimeInterval) :: ts
  integer                 :: cyr, cmo, cda, chr, cmn, css
  integer                 :: status
  logical*1              :: lo(LVT_rc%lnc*LVT_rc%lnr)
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: varfield1(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: br(LVT_rc%lnc,LVT_rc%lnr)
  real                   :: ef(LVT_rc%lnc,LVT_rc%lnr)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
  if((FLUXNETmteobs(Source)%yr.ne.LVT_rc%dyr(source)).or.&
       LVT_rc%resetFlag(source)) then      

     LVT_rc%resetFlag(source) = .false. 

     FLUXNETmteobs(Source)%yr = LVT_rc%dyr(source)

     FLUXNETmteobs(Source)%qle = LVT_rc%udef
     FLUXNETmteobs(Source)%qh = LVT_rc%udef
     
     call create_fluxnet_lh_filename(FLUXNETmteobs(Source)%odir, &
          LVT_rc%dyr(source),filename)
     inquire(file=trim(filename),exist=file_exists) 
     
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading FLUXNET MTE LH file ',trim(filename)
        
        
        ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
        call LVT_verify(ios, 'Error opening file '//trim(filename))
        
        ios = nf90_inq_varid(nid,'EnsembleLEcor_May09',qleid)
        call LVT_verify(ios, 'Error nf90_inq_varid: qle')
        ! dimensions
        ios = nf90_inq_dimid(nid,'latitude',rowId)
        call LVT_verify(ios, 'Error nf90_inq_dimid: latitude')
        
        ios = nf90_inquire_dimension(nid,rowId, len=nr)
        call LVT_verify(ios, 'Error nf90_inquire_dimension: row')
        
        ios = nf90_inq_dimid(nid,'longitude',colId)
        call LVT_verify(ios, 'Error nf90_inq_dimid: cols')
        
        ios = nf90_inquire_dimension(nid,colId, len=nc)
        call LVT_verify(ios, 'Error nf90_inquire_dimension: col')
        
        allocate(qle(nc,nr,12))
        allocate(qle1d(nc*nr))
        allocate(li(nc*nr))
        li = .true. 
        
        !values
        ios = nf90_get_var(nid,qleid, qle)
        call LVT_verify(ios, 'Error nf90_get_var: qle')
        
        ios = nf90_close(nid)
        call LVT_verify(ios, 'Error in nf90_close')
        
        do t=1,12
           qle1d = -9999.0
           do r=1,nr
              do c=1,nc
                 if(isNaN(qle(c,nr-r+1,t))) then 
                    qle1d(c+(r-1)*nc) = LVT_rc%udef
                    li(c+(r-1)*nc) = .false. 
                 elseif(qle(c,nr-r+1,t).lt.0) then 
                    qle1d(c+(r-1)*nc) = LVT_rc%udef
                    li(c+(r-1)*nc) = .false. 
                 else
                    qle1d(c+(r-1)*nc) = qle(c,nr-r+1,t)*0.01*1E6/86400
                    li(c+(r-1)*nc) = .true. 
                 endif
              enddo
           enddo
           
           call neighbor_interp(LVT_rc%gridDesc,li,qle1d,&
                lo, FLUXNETmteobs(Source)%qle(:,t), nc*nr, &
                LVT_rc%lnc*LVT_rc%lnr,&
                FLUXNETmteobs(Source)%rlat, FLUXNETmteobs(Source)%rlon, &
                FLUXNETmteobs(Source)%n11,LVT_rc%udef, ios)
        enddo
           
        deallocate(qle)
        deallocate(qle1d)
        deallocate(li)
     endif
     
     call create_fluxnet_sh_filename(FLUXNETmteobs(Source)%odir,&
          LVT_rc%dyr(source),filename)
     
     inquire(file=trim(filename),exist=file_exists) 
     
     if(file_exists) then 
        write(LVT_logunit,*) '[INFO] Reading FLUXNET MTE SH file ',trim(filename)
        
        
        ios = nf90_open(path=trim(filename),mode=NF90_NOWRITE,ncid=nid)
        call LVT_verify(ios, 'Error opening file '//trim(filename))
        
        ios = nf90_inq_varid(nid,'EnsembleHcor_Aug09',qhid)
        call LVT_verify(ios, 'Error nf90_inq_varid: qh')
        ! dimensions
        ios = nf90_inq_dimid(nid,'latitude',rowId)
        call LVT_verify(ios, 'Error nf90_inq_dimid: latitude')
        
        ios = nf90_inquire_dimension(nid,rowId, len=nr)
        call LVT_verify(ios, 'Error nf90_inquire_dimension: row')
        
        ios = nf90_inq_dimid(nid,'longitude',colId)
        call LVT_verify(ios, 'Error nf90_inq_dimid: cols')
        
        ios = nf90_inquire_dimension(nid,colId, len=nc)
        call LVT_verify(ios, 'Error nf90_inquire_dimension: col')
        
        allocate(qh(nc,nr,12))
        allocate(qh1d(nc*nr))
        allocate(li(nc*nr))
        li = .true. 
        !values
        ios = nf90_get_var(nid,qhid, qh)
        call LVT_verify(ios, 'Error nf90_get_var: qh')
        
        ios = nf90_close(nid)
        call LVT_verify(ios, 'Error in nf90_close')
        
        do t=1,12
           qh1d = -9999.0
           do r=1,nr
              do c=1,nc
                 if(isNaN(qh(c,nr-r+1,t))) then 
                    qh1d(c+(r-1)*nc) = LVT_rc%udef
                    li(c+(r-1)*nc) = .false. 
                 elseif(qh(c,nr-r+1,t).lt.0) then 
                    qh1d(c+(r-1)*nc) = LVT_rc%udef
                    li(c+(r-1)*nc) = .false. 
                 else
                    qh1d(c+(r-1)*nc) = qh(c,nr-r+1,t)*0.01*1E6/86400
                    li(c+(r-1)*nc) = .true. 
                 endif
              enddo
           enddo
           
           call neighbor_interp(LVT_rc%gridDesc,li,qh1d,&
                lo, FLUXNETmteobs(Source)%qh(:,t), nc*nr, &
                LVT_rc%lnc*LVT_rc%lnr,&
                FLUXNETmteobs(Source)%rlat, FLUXNETmteobs(Source)%rlon, &
                FLUXNETmteobs(Source)%n11,LVT_rc%udef, ios)
        enddo
           
        deallocate(qh)
        deallocate(qh1d)
        deallocate(li)
     endif
  endif
#endif

  !log the data only at the change of a month. 
  if(LVT_rc%d_nmo(source).ne.FLUXNETmteobs(Source)%mo) then 
     cmo = FLUXNETmteobs(Source)%mo
     write(LVT_logunit,*) '[INFO] Reading FLUXNET MTE data for ',cmo
     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           varfield(c,r) = FLUXNETmteobs(Source)%qle(c+(r-1)*LVT_rc%lnc, &
                cmo)
        enddo
     enddo

     do r=1, LVT_rc%lnr
        do c=1, LVT_rc%lnc
           varfield1(c,r) = FLUXNETmteobs(Source)%qh(c+(r-1)*LVT_rc%lnc, &
                cmo)
        enddo
     enddo
     
     FLUXNETmteobs(Source)%mo = LVT_rc%d_nmo(source)

  else
     varfield  = LVT_rc%udef
     varfield1 = LVT_rc%udef
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_QLE,source,varfield,&
       vlevel=1,units="W/m2")
  call LVT_logSingleDataStreamVar(LVT_MOC_QH,source,varfield1,&
       vlevel=1,units="W/m2")

  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(varfield(c,r).ne.-9999.0.and.varfield1(c,r).ne.-9999.0) then 
           if(varfield(c,r).ne.0) then 
              br(c,r) = varfield1(c,r)/varfield(c,r)
           else
              br(c,r) = LVT_rc%udef
           endif
        else
           br(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_BR,source,br,vlevel=1,units="-")


  do r=1,LVT_rc%lnr
     do c=1,LVT_rc%lnc
        if(varfield(c,r).ne.-9999.0.and.varfield1(c,r).ne.-9999.0) then 
           if((varfield(c,r)+varfield1(c,r)).ne.0) then 
              ef(c,r) = varfield(c,r)/(varfield(c,r)+varfield1(c,r))
           else
              ef(c,r) = LVT_rc%udef
           endif
        else
           ef(c,r) = LVT_rc%udef
        endif
     enddo
  enddo
  call LVT_logSingleDataStreamVar(LVT_MOC_EF,source,ef,vlevel=1,units="-")
  
end subroutine readFLUXNETmteObs

!BOP
! 
! !ROUTINE: create_fluxnet_lh_filename
! \label{create_fluxnet_lh_filename}
!
! !INTERFACE: 
subroutine create_fluxnet_lh_filename(odir,yr,filename)
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
! This routine creates a timestamped filename for FLUXNETmte_LH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the FLUXNETmte_LH file
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
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  
  write(unit=fyr, fmt='(i4.4)') yr
  
  filename = trim(odir)//'/EnsembleLEcor_May09_'//trim(fyr)//'.nc'
  
end subroutine create_fluxnet_lh_filename

!BOP
! 
! !ROUTINE: create_fluxnet_sh_filename
! \label{create_fluxnet_sh_filename}
!
! !INTERFACE: 
subroutine create_fluxnet_sh_filename(odir,yr,filename)
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
! This routine creates a timestamped filename for FLUXNETmte_SH data files 
! based on the given date (year, month, day)
!
!  The arguments are: 
!  \begin{description}
!   \item[odir]      GlboSnow base directory
!   \item[yr]        year of data
!   \item[filename]  Name of the FLUXNETmte_SH file
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
  character(len=*)             :: filename
!EOP

  character*4             :: fyr
  
  write(unit=fyr, fmt='(i4.4)') yr
  
  filename = trim(odir)//'/EnsembleHcor_Aug09_'//trim(fyr)//'.nc'
  
end subroutine create_fluxnet_sh_filename

