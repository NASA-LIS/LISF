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
! !ROUTINE: readprecip_RFE2Daily
! \label{readprecip_RFE2Daily}
!
! !REVISION HISTORY:
!  30 May 2010; Soni Yatheendradas, Initial LDT version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and documentation
!
! !INTERFACE:
subroutine readprecip_RFE2Daily( n, fname, findex, order, ferror_RFE2Daily )
! !USES:
  use LDT_coreMod,        only : LDT_rc, LDT_domain, LDT_masterproc
  use LDT_logMod,         only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_endrun, LDT_verify
  use LDT_metforcingMod,  only : LDT_forc
  use RFE2Daily_forcingMod, only : RFE2Daily_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  character(len=*)  :: fname
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_RFE2Daily
! 
! !DESCRIPTION:
!  For the given time, reads the RFE2Daily data
!  and interpolates to the LDT domain.
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the daily RFE2.0 file
!  \item[findex]
!    index of the supplemental forcing source
!  \item[ferror\_RFE2Daily]
!    flag to indicate success of the call (=0 indicates success)
!  \end{description}
!  
!  The routines invoked are:
!  \begin{description}
!  \item[reproject\_RFE2Daily](\ref{reproject_RFE2Daily}) \newline
!    Upscales/interpolates the RFE2Daily data
!  \end{description}
!EOP

!==== Local Variables=======================
  integer               :: c,r,t
  logical*1,allocatable :: lb1d(:)
  logical*1,allocatable :: lb2d(:,:)
  real, allocatable     :: rain1d(:)
  REAL, ALLOCATABLE     :: rain2d(:,:)
  integer               :: ftn, ios, ftn2, ftn3
  character(len=124)    :: fnametemp
  LOGICAL               :: file_exists
  real, dimension(LDT_rc%lnc(n), LDT_rc%lnr(n)) :: varfield ! reprojected arrray
!=== End Variable Definition =======================

  if(order.eq.1) then 
     LDT_forc(n,findex)%metdata1 = LDT_rc%udef
  elseif(order.eq.2) then 
     LDT_forc(n,findex)%metdata2 = LDT_rc%udef
  endif

  varfield = LDT_rc%udef
 
  allocate( lb1d(RFE2Daily_struc(n)%mi) )
  allocate( lb2d(RFE2Daily_struc(n)%nc,RFE2Daily_struc(n)%nr) )
  allocate( rain1d(RFE2Daily_struc(n)%mi) )
  allocate( rain2d(RFE2Daily_struc(n)%nc,RFE2Daily_struc(n)%nr) ) 

  INQUIRE(FILE=fname, EXIST=file_exists)   ! file_exists will be TRUE if the file
                                           ! exists and FALSE otherwise

  IF( file_exists .EQV. .TRUE. ) THEN

    ftn = LDT_getNextUnitNumber()
  
    open(unit=ftn,file=fname, access='direct', &
        recl=RFE2Daily_struc(n)%mi*4, &
        iostat=ios)
   
    if(ios .eq. 0) then 
     
     ! Read 2-d RFE2 precip values in:
      rain2d = LDT_rc%udef
      read (ftn,rec=1)  rain2d
      close(ftn)

      !- Upscaling/averaging 2D to 1D assignment:
      select case( LDT_rc%met_gridtransform(findex) )

        case( "average" )   ! Upscaling 
          lb1d = .false.
          do r=1, RFE2Daily_struc(n)%nr
            do c=1, RFE2Daily_struc(n)%nc
              if(rain2d(c,r).GT.(-999.0+1)) then
                t = c+(r-1)*RFE2Daily_struc(n)%nc
                lb1d(t) = .true.
                rain1d(t) = rain2d(c,r) 
              else
                rain1d(t)=LDT_rc%udef 
              endif
            enddo
          enddo
    
       case( "bilinear", "budget-bilinear", "neighbor"  )  
          !lb1d = .false. ! SY: See next line
          lb1d = .true. ! SY. NOTE. Any NODATA rain values also correspondingly fixed to valid 0 below  
          do r=1, RFE2Daily_struc(n)%nr
             do c=1, RFE2Daily_struc(n)%nc
                t = c+(r-1)*RFE2Daily_struc(n)%nc
                rain1d(t) = rain2d(c,r) 
                if( rain1d(t) .LT. (-999.0+1) ) then ! SY
                  rain1d(t) = 0.0
                  !lb1d(t) = .true. ! SY
                endif
             enddo
          enddo

       case default
          write(LDT_logunit,*)"RFE2Daily reprojection choice should be: "
          write(LDT_logunit,*) " average, neighbor, bilinear, or budget-bilinear"
          write(LDT_logunit,*) "Program stopping ... "
          call LDT_endrun()
      end select
    
    ! Spatially reproject RFE2 Daily files:
      call reproject_RFE2Daily(n,findex,&
                 RFE2Daily_struc(n)%mi, &
                 rain1d, lb1d, LDT_rc%gridDesc(n,:), &
                 LDT_rc%lnc(n),LDT_rc%lnr(n), varfield )
  
      do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
          if ( (LDT_domain(n)%gindex(c,r).ne.-1) .AND. &
               (varfield(c,r) .ge. 0) ) then       ! SY: leaves out -1 and NODATA values 
            if (varfield(c,r) .GT. 1000.0) THEN    ! 1000 mm threshold
              varfield(c,r) = 1000.0 
            endif
            if(order.eq.1) then 
               LDT_forc(n,findex)%metdata1(1,LDT_domain(n)%gindex(c,r)) = &
                    varfield(c,r)
            elseif(order.eq.2) then 
               LDT_forc(n,findex)%metdata2(1,LDT_domain(n)%gindex(c,r)) = &
                    varfield(c,r)
            endif
          endif 
        end do
      end do
  
      ferror_RFE2Daily = 0
      if(LDT_masterproc) write(LDT_logunit,*) &
                         "Read RFE2Daily data file: ", trim(fname)
    else
      if(LDT_masterproc) write(LDT_logunit,*) &
                     "Cannot open RFE2Daily precipitation file ", trim(fname)
      ferror_RFE2Daily = 1
      LDT_forc(n,findex)%metdata1 = LDT_rc%udef
      LDT_forc(n,findex)%metdata2 = LDT_rc%udef

    endif
    
    call LDT_releaseUnitNumber(ftn) 

  ELSE  ! For IF (file_exists .NEQV. .TRUE.) THEN

    if(LDT_masterproc) write(LDT_logunit,*) &
                   "Missing RFE2Daily precipitation file ", trim(fname)
    ferror_RFE2Daily = 1
    LDT_forc(n,findex)%metdata1 = LDT_rc%udef
    LDT_forc(n,findex)%metdata2 = LDT_rc%udef
    
  ENDIF 

  deallocate(lb1d)
  deallocate(lb2d)
  deallocate(rain1d)
  deallocate(rain2d)

end subroutine readprecip_RFE2Daily
