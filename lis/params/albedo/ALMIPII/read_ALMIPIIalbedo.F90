!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
!
! !ROUTINE: read_ALMIPIIalbedo
! \label{read_ALMIPIIalbedo}
!
! !REVISION HISTORY:
!  25 Jul 2005: Sujay Kumar; Initial Specification
!  20 Feb 2006: Sujay Kumar; Modified to support nesting
!
! !INTERFACE:
subroutine read_ALMIPIIalbedo(n, wt1,wt2,array1,array2)
! !USES:
  use ESMF
  use LIS_coreMod
  use LIS_logMod,         only : LIS_logunit, LIS_verify, LIS_endrun
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use LIS_albedoMod,  only : LIS_alb
  use LIS_timeMgrMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 

  integer, intent(in) :: n
  real                :: wt1
  real                :: wt2
  real, intent(inout) :: array1(LIS_rc%ntiles(n))
  real, intent(inout) :: array2(LIS_rc%ntiles(n))

! !DESCRIPTION:
!  This subroutine retrieves the greenness fraction climatology for the 
!  specified month and returns the values in the latlon projection
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!   index of the nest
!  \item[mo]
!   time index (month or quarter)
!  \item[array]
!   output field with the retrieved greenness fraction
!  \end{description}
!
!EOP      
  character(len=LIS_CONST_PATH_LEN) :: filename
  character*100               :: temp
  logical                     :: file_exists
  integer                     :: yr
  integer                     :: rc
  character*1                 :: fyr(4)
  integer                     :: ftn
  integer                     :: t1,t2
  logical                     :: albAlarmCheck
  integer                     :: alb1Id,alb2Id
  real, allocatable                       :: albvis(:,:)
  real, allocatable                       :: albnir(:,:)
  real, allocatable                       :: alb_t(:,:)
  real, allocatable                       :: localalb(:,:)
  integer                     :: i,c,r,t


  allocate(albvis(LIS_rc%gnc(n),LIS_rc%gnr(n)))
  allocate(albnir(LIS_rc%gnc(n),LIS_rc%gnr(n)))
  allocate(alb_t(LIS_rc%gnc(n),LIS_rc%gnr(n)))
  allocate(localalb(LIS_rc%lnc(n),LIS_rc%lnr(n)))

  if(LIS_alb(n)%firstInstance) then 
     albAlarmCheck = .true. 
  else
     albAlarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "LIS alb read alarm",&
          LIS_alb(n)%albIntervalType)
  endif

  call LIS_computeTemporalWeights(LIS_rc,LIS_alb(n)%albIntervalType, &
       t1,t2,wt1,wt2)

  if(albAlarmCheck) then 
     if(LIS_alb(n)%firstInstance) &
          LIS_alb(n)%firstInstance = .false. 
     
     array1 = LIS_rc%udef
     array2 = LIS_rc%udef

#if (defined USE_NETCDF3 ||  defined USE_NETCDF4) 
     write(unit=temp,fmt='(I4)') LIS_rc%yr
     read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
     
     filename = trim(LIS_alb(n)%albfile)//&
          "/ALMIP2_ECOCLIMAP2_"//trim(fyr(3))//trim(fyr(4))&
          //'.nc'
     
     inquire(file=trim(filename), exist=file_exists)
     if(.not.file_exists) then 
        write(LIS_logunit,*) 'albedo map ',trim(filename),' not found'
        write(LIS_logunit,*) 'Program stopping ...'
        call LIS_endrun
     endif
     
     write(LIS_logunit,*) 'opening ALB file ',trim(filename)
     call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
          ncid=ftn), 'nf90_open failed in read_ALMIPIIalb')
     
     call LIS_verify(nf90_inq_varid(ftn,'alb_vis',alb1Id),&
          'nf90_inq_varid failed in read_ALMIPIIalb')
     call LIS_verify(nf90_get_var(ftn,alb1Id,albvis,&
          start=(/1,1,t1/),&
          count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)),&
          'nf90_get_var failed in read_ALMIPIIalb')

     call LIS_verify(nf90_inq_varid(ftn,'alb_nir',alb2Id),&
          'nf90_inq_varid failed in read_ALMIPIIalb')
     call LIS_verify(nf90_get_var(ftn,alb2Id,albnir,&
          start=(/1,1,t1/),&
          count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)),&
          'nf90_get_var failed in read_ALMIPIIalb')
          

     do r=1,LIS_rc%gnr(n)
        do c=1,LIS_rc%gnc(n)
!           alb_t(c,r) = alb(c,LIS_rc%gnr(n)-r+1)
           alb_t(c,r) = 0.5*(albvis(c,r)+albnir(c,r))
        enddo
     enddo
     
     localalb(:,:) = &
          alb_t(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))
     
     do t=1,LIS_rc%ntiles(n)
        c = LIS_domain(n)%tile(t)%col
        r = LIS_domain(n)%tile(t)%row
        array1(t) = localalb(c,r)
     enddo
     
     call LIS_verify(nf90_close(ftn),'nf90_close failed in read_ALMIPIIalb')

     write(unit=temp,fmt='(I4)') LIS_rc%yr
     read(unit=temp,fmt='(4a1)') (fyr(i),i=1,4)
     
     filename = trim(LIS_alb(n)%albfile)//&
          "/ALMIP2_ECOCLIMAP2_"//trim(fyr(3))//trim(fyr(4))&
          //'.nc'
     
     inquire(file=trim(filename), exist=file_exists)
     if(.not.file_exists) then 
        write(LIS_logunit,*) 'albedo map ',trim(filename),' not found'
        write(LIS_logunit,*) 'Program stopping ...'
        call LIS_endrun
     endif
     
     write(LIS_logunit,*) 'opening ALB file ',trim(filename)

     call LIS_verify(nf90_open(path=trim(filename),mode=NF90_NOWRITE,&
          ncid=ftn), 'nf90_open failed in read_ALMIPIIalb')
     
     call LIS_verify(nf90_inq_varid(ftn,'alb_vis',alb1Id),&
          'nf90_inq_varid failed in read_ALMIPIIalb')
     call LIS_verify(nf90_get_var(ftn,alb1Id,albvis,&
          start=(/1,1,t2/),&
          count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)),&
          'nf90_get_var failed in read_ALMIPIIalb')
     call LIS_verify(nf90_inq_varid(ftn,'alb_nir',alb2Id),&
          'nf90_inq_varid failed in read_ALMIPIIalb')
     call LIS_verify(nf90_get_var(ftn,alb2Id,albnir,&
          start=(/1,1,t2/),&
          count=(/LIS_rc%gnc(n),LIS_rc%gnr(n),1/)),&
          'nf90_get_var failed in read_ALMIPIIalb')
     
     do r=1,LIS_rc%gnr(n)
        do c=1,LIS_rc%gnc(n)
!           alb_t(c,r) = alb(c,LIS_rc%gnr(n)-r+1)
           alb_t(c,r) = 0.5*(albvis(c,r) + albnir(c,r))
        enddo
     enddo
     
     localalb(:,:) = &
          alb_t(LIS_ews_halo_ind(n,LIS_localPet+1):&         
          LIS_ewe_halo_ind(n,LIS_localPet+1), &
          LIS_nss_halo_ind(n,LIS_localPet+1): &
          LIS_nse_halo_ind(n,LIS_localPet+1))
     
     do t=1,LIS_rc%ntiles(n)
        c = LIS_domain(n)%tile(t)%col
        r = LIS_domain(n)%tile(t)%row
        array2(t) = localalb(c,r)
     enddo
     
     call LIS_verify(nf90_close(ftn),'nf90_close failed in read_ALMIPIIalb')
  
#endif
  endif

  deallocate(albvis)
  deallocate(albnir)
  deallocate(alb_t)
  deallocate(localalb)
end subroutine read_ALMIPIIalbedo
