!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Data Toolkit (LDT) v7.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! !ROUTINE: read_nldas1
!  \label{read_nldas1}
!
! !REVISION HISTORY:
!  11 Apr 2000: Brian Cosgrove; changed code to use Forcing Mask (With inland
!               water filled in).  Deleteted unused variables.
!  27 Apr 2000: Brian Cosgrove; changed code to use the original 
!               mask again since that is the
!               mask which  NCEP has already applied to the forcing data
!               by the time NASA gets it......not possible to use the 
!               expanded NASA forcing mask
!  1  May 2000: Brian Cosgrove; changed code so that if parameter 11 (sw)
!               is not found in hourly ncep data, it will just use
!               edas-based shortwave from the hourly ncep files
!  20 Jun 2000: Brian Cosgrove; changed code so that it uses  LDAS%UDEF and
!                not a hard-wired undefined value of -999.9 and -999.0
!  18 Aug 2000: Brian Cosgrove; changed code so that FMASK and not MASK
!               is used when ungribbing.  NCEP data already has a mask applied
!               to it and so may not be able to supply forcing data to
!               all LDAS land forcing points.  In areas where LDAS
!               forcing mask states that land exists, but where NCEP forcing
!               data is non-existant, assign undefined value to forcing data.
!  22 Aug 2000: Brian Cosgrove; Altered code for US/Mexico/Canada Mask
!  05 Sep 2001: Brian Cosgrove; Removed dirnom and infile variables, changed
!               call to ungribncep to match removal.  Added code to make use
!               of precip weighting mask
!  02 Feb 2004: Sujay Kumar; Initial Specification in LDT
!  20 Oct 2007: Kristi Arsenault; Updated NLDAS-1 Code to Read EDAS Height Field
!  25 Jan 2012: Sujay Kumar; Switched to the use of grib-api library
!  15 Feb 2012: Kristi Arsenault; Accommodate GES DISC filename conventions
!  14 Mar 2014: David Mocko: Fixed the use of NLDAS-1 GOES SWdown radiation.
!                            Added NLDAS-1 precipitation and radiation
!                            field choices to the ldt.config file.
! 
! !INTERFACE:
subroutine read_nldas1(n, findex, order, name, ferror)
! !USES:
  use LDT_coreMod, only       : LDT_rc, LDT_domain
  use LDT_metforcingMod, only : LDT_forc
  use LDT_logMod, only        : LDT_logunit, LDT_verify, LDT_warning, LDT_endrun
  use nldas1_forcingMod, only : nldas1_struc
#if (defined USE_GRIBAPI)
  use grib_api
#endif

  implicit none
! !ARGUMENTS: 
  integer, intent(in)      :: n
  integer, intent(in)      :: findex
  integer, intent(in)      :: order
  character*80, intent(in) :: name
  integer, intent(out)     :: ferror
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  NLDAS-1 data, transforms into 9 LDT forcing 
!  parameters and interpolates to the LDT spatial domain.
!  The routine reads the gauge-corrected precip and GOES-based
!  radiation fields by default. If these fields are not present, 
!  then the code defaults to the corresponding EDAS fields. 
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[findex]
!    sequential index of the NLDAS-1 forcing among various forcings
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    hourly instance, order=2, read the next hourly instance)
!  \item[name]
!    name of the NLDAS-1 forcing file
!  \item[ferror]
!    flag to indicate success of the call (=1 indicates success)
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[interp\_nldas1](\ref{interp_nldas1}) \newline
!    spatially interpolates a NLDAS-1 variable
!  \item[nldas1\_ec\_removal](\ref{nldas1_ec_removal}) \newline
!    removes the native elevation correction applied to the 
!    NLDAS fields. This is performed so that the elevation
!    correction can be reapplied at the spatial resolution of
!    the LDT domain
!   \item[apply\_conusmask](\ref{apply_conusmask}) \newline
!    applies the CONUS landmask to exclude NLDAS-1 data over 
!    Canada and Mexico. 
!  \end{description}
!EOP
  integer                :: iv,ftn, iv_total
  logical                :: file_exists
  integer                :: iret,rc
  integer                :: nldas1,c,r,t,j
  integer                :: pds5_val, pds7_val, pds2_val
  real                   :: missingValue 
  integer                :: nvars
  integer                :: igrib
  logical                :: pcp_flag
  logical*1, allocatable :: lb(:)
  real, allocatable      :: f(:,:)
  real, allocatable      :: f2(:,:)
  real                   :: varfield(LDT_rc%lnc(n),LDT_rc%lnr(n))

  ferror = 1
  iv_total = 15
  nldas1 = (nldas1_struc(n)%nc*nldas1_struc(n)%nr)
! Original NLDAS-1 GRIB parameter table
!  pds5 = (/ 011,051,204,205,033,034,001,061,063,157/) !parameter
!  pds7 = (/ 002,002,000,000,010,010,000,000,000,000/) !htlev2
!  pds2 = (/ 84, 84, 84, 84, 84, 84, 84, 155, 84,84/)
! Updated NLDAS-1 GRIB parameter table (not valid for all hours in record)
!  pds5 = (/ 011,051,001,033,034,204,205,061,063,157,204,118,101,061,061/) !parameter
!  pds7 = (/ 002,002,000,010,010,000,000,000,000,000,000,000,000,000,000/) !htlev2
!  pds2 = (/ 084,084,084,084,084,084,084,084,084,084,154,154,154,155,156/)
  allocate(lb(nldas1))
  allocate(f(nldas1, iv_total))
  allocate(f2(nldas1,9))

  f = LDT_rc%udef
  f2 = LDT_rc%udef

#if(defined USE_GRIBAPI) 

  inquire (file=trim(name), exist=file_exists)
  if (file_exists) then   

     call grib_open_file(ftn,trim(name),'r',iret)
     if(iret.ne.0) then 
        write(LDT_logunit,*) &
             'Could not open file: ',trim(name)
        ferror = 0
        return
     endif
     call grib_count_in_file(ftn,nvars,iret)
     call LDT_verify(iret, 'error in grib_count_in_file in read_nldas1')
     
     do iv=1,nvars

        call grib_new_from_file(ftn, igrib, iret)
        call LDT_warning(iret,&
             'error in grib_new_from_file in read_nldas1')
        if(iret.ne.0) then 
           write(LDT_logunit,*) &
                'Failed with grib_new_from_file in read_nldas: ',trim(name)
           ferror = 0
           deallocate(f)
           deallocate(f2)
           deallocate(lb)
           return           
        endif
        
        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LDT_verify(rc,&
               'error in grib_get: indicatorOfParameter in read_nldas1')
        
        call grib_get(igrib,'level',pds7_val,rc)
        call LDT_verify(rc, 'error in grib_get: level in read_nldas1')
        
        call grib_get(igrib,'generatingProcessIdentifier',pds2_val,rc)
        call LDT_verify(rc,&
        'error in grib_get: generatingProcessIdentifier in read_nldas1')
        
        call grib_get(igrib,'values',f(:,iv),rc)
        call LDT_warning(rc,'error in grib_get:values in read_nldas1')
        
        if(rc.ne.0) then 
           write(LDT_logunit,*) &
                'Could not retrieve entries in file: ',trim(name),&
                'for variable ',iv
           ferror = 0
           deallocate(f)
           deallocate(f2)
           deallocate(lb)
           return
        endif
           
        call grib_get(igrib,'missingValue',missingValue,rc) 
        call LDT_verify(rc,&
                        'error in grib_get:missingValue in read_nldas1')
        
        call grib_release(igrib,rc)
        call LDT_verify(rc, 'error in grib_release in read_nldas1')       
     enddo

     call grib_close_file(ftn)

  else
     write(LDT_logunit,*) 'Could not find file: ',trim(name)
     ferror = 0
     deallocate(f)
     deallocate(f2)
     deallocate(lb)
     return
  endif
  
  if ( trim(LDT_rc%met_ecor(findex)) .ne. "none") then  
      do c = 1, nldas1
         if(lb(c)) then 
            call nldas1_ec_removal( n, c, f(c,1), f(c,2), f(c,7), f(c,3) )
         endif
      end do
   end if

! Map NLDAS-1 order of fields to the order expected in timeinterp - dmm
   if ((nldas1_struc(n)%prec_field.eq."NLDAS1").and.(nvars.le.13)) then
      write(LDT_logunit,*) 'NLDAS1 precipitation field not '//&
             'available for this hour, so using EDAS field '
   endif
   if ((nldas1_struc(n)%prec_field.eq."STAGEII").and.(nvars.le.14)) then
      write(LDT_logunit,*) 'STAGEII precipitation field not '//&
              'available for this hour, so using EDAS field '
   endif
   do t = 1,nldas1
      f2(t,1) = f(t,1)
      f2(t,2) = f(t,2)
      if (nldas1_struc(n)%swdn_field.eq."NLDAS1") then
!         if ( f(t,6) .ge. -9998.0 ) then
         if ( f(t,6) .ne. missingValue ) then
            f2(t,3) = f(t,6)
         else
            f2(t,3) = 0.0
         endif
!         if (f(t,11).ge.-9998.0) then
         ! I have seen values of -9999 and -9999000
         ! in message 11.  Since this is for swdown,
         ! exclude any negative values.
         if (f(t,11).ne.missingValue .and. &
             f(t,11).ge.0.0) then
            f2(t,3) = f(t,11)
         endif
      elseif (nldas1_struc(n)%swdn_field.eq."EDAS") then
!         if ( f(t,6) .ge. -9998.0 ) then
         if ( f(t,6) .ne. missingValue ) then
            f2(t,3) = f(t,6)
         else
            f2(t,3) = 0.0
         endif
      else
         write(LDT_logunit,*) 'NLDAS1 shortwave radiation field:'
         write(LDT_logunit,*) ' must be equal to either "NLDAS1"'
         write(LDT_logunit,*) ' or to "EDAS".  Please check your'
         write(LDT_logunit,*) ' config file.  Stopping....'
         call LDT_endrun()
      endif
      f2(t,4) = f(t,7)
      f2(t,5) = f(t,4)
      f2(t,6) = f(t,5)
      f2(t,7) = f(t,3)
      if (nldas1_struc(n)%prec_field.eq."NLDAS1") then
         if (nvars.ge.14) then
            f2(t,8) = f(t,14)
         else
            f2(t,8) = f(t,8)
         endif
         f2(t,9) = f(t,9)
      elseif (nldas1_struc(n)%prec_field.eq."EDAS") then
         f2(t,8) = f(t,8)
      elseif (nldas1_struc(n)%prec_field.eq."STAGEII") then
         if (nvars.ge.14) then
            f2(t,8) = f(t,14)
         else
            f2(t,8) = f(t,8)
         endif
         if ((nvars.ge.15).and.(f(t,15).ge.-9998.0)) then
            f2(t,8) = f(t,15)
         endif
      else
         write(LDT_logunit,*) 'NLDAS1 precipitation field:'
         write(LDT_logunit,*) ' must be equal to "NLDAS1", "EDAS",'
         write(LDT_logunit,*) ' or to "STAGEII".  Please check your'
         write(LDT_logunit,*) ' config file.  Stopping....'
         call LDT_endrun()
      endif
      f2(t,9) = max(min(f(t,9),f2(t,8)),0.0)
   enddo

!- Reinterpolate original NLDAS-1 forcing to LDT grid:
  do iv=1,9
     lb = .false. 
     do t=1,nldas1
        if(f2(t,iv).ne.missingValue) then 
           lb(t) = .true. 
        endif
     enddo
           
     pcp_flag = .false. 
     if(iv.eq.8.or.iv.eq.9) pcp_flag = .true. 
     
     call interp_nldas1(n, findex, pcp_flag,nldas1,f2(:,iv),&
          lb,LDT_rc%gridDesc(n,:),& 
          LDT_rc%lnc(n),LDT_rc%lnr(n),varfield)
     
     if(nldas1_struc(n)%applymask.eq.1) then 
        call apply_conusmask(n, LDT_rc%lnc(n), LDT_rc%lnr(n), &
             varfield)
     endif

      do r=1,LDT_rc%lnr(n)
         do c=1,LDT_rc%lnc(n)
            if(LDT_domain(n)%gindex(c,r).ne.-1) then 
               if(order.eq.1) then 
                  LDT_forc(n,findex)%metdata1(iv,LDT_domain(n)%gindex(c,r)) =&
                       varfield(c,r)
               elseif(order.eq.2) then 
                  LDT_forc(n,findex)%metdata2(iv,LDT_domain(n)%gindex(c,r)) =&
                       varfield(c,r)
               endif
            endif
         end do
      enddo
   enddo
   deallocate(f)
   deallocate(f2)
   deallocate(lb)
#endif     
      
   return
 end subroutine read_nldas1


