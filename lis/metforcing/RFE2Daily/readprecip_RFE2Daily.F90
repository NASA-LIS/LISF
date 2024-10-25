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
!  30 May 2010; Soni Yatheendradas, Initial LIS version for FEWSNET
!  20 Mar 2013; KR Arsenault, Cleaned up code and documentation
!
! !INTERFACE:
subroutine readprecip_RFE2Daily( n, kk, findex, fname, order, ferror_RFE2Daily)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
          LIS_releaseUnitNumber, LIS_endrun, LIS_verify
  use LIS_metforcingMod,  only : LIS_forc
  use RFE2Daily_forcingMod, only : RFE2Daily_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: kk
  integer, intent(in) :: findex
  character(len=*)   :: fname
  integer, intent(in) :: order
  integer             :: ferror_RFE2Daily
! 
! !DESCRIPTION:
!  For the given time, reads the RFE2Daily data
!  and interpolates to the LIS domain.
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
  LOGICAL               :: file_exists
  real, dimension(LIS_rc%lnc(n), LIS_rc%lnr(n)) :: varfield ! reprojected arrray

!=== End Variable Definition =======================

  if(order.eq.1) then 
     RFE2Daily_struc(n)%metdata1 = LIS_rc%udef
  elseif(order.eq.2) then 
     RFE2Daily_struc(n)%metdata2 = LIS_rc%udef
  endif

  varfield = LIS_rc%udef
 
  allocate(lb1d(NINT(RFE2Daily_struc(n)%gridDesci(2))*&
                NINT(RFE2Daily_struc(n)%gridDesci(3))))
  allocate(lb2d(NINT(RFE2Daily_struc(n)%gridDesci(2)),&
                NINT(RFE2Daily_struc(n)%gridDesci(3))))
  allocate(rain1d(NINT(RFE2Daily_struc(n)%gridDesci(2))*&
                  NINT(RFE2Daily_struc(n)%gridDesci(3))))
  allocate(rain2d(NINT(RFE2Daily_struc(n)%gridDesci(2)),&
                  NINT(RFE2Daily_struc(n)%gridDesci(3))))

  INQUIRE(FILE=fname, EXIST=file_exists)   ! file_exists will be TRUE if the file
                                           ! exists and FALSE otherwise

  IF (file_exists .EQV. .TRUE.) THEN

    ftn = LIS_getNextUnitNumber()
  
    open(unit=ftn,file=fname, access='direct', &
        recl=NINT(RFE2Daily_struc(n)%gridDesci(2))* &
             NINT(RFE2Daily_struc(n)%gridDesci(3))*4, &
        iostat=ios)
   
    if(ios .eq. 0) then 
     
      read (ftn,rec=1)  rain2d
      close(ftn)
   

   !- Upscaling/averaging 2D to 1D assignment:
      select case( LIS_rc%met_interp(findex) )

        case( "average" )   ! Upscaling 
          lb1d = .false.
          do r=1,NINT(RFE2Daily_struc(n)%gridDesci(3))
            do c=1,NINT(RFE2Daily_struc(n)%gridDesci(2))
              if(rain2d(c,r).GT.(-999.0+1)) then
                t = c+(r-1)*NINT(RFE2Daily_struc(n)%gridDesci(2))
                lb1d(t) = .true.
                rain1d(t) = rain2d(c,r) 
              else
                rain1d(t)=LIS_rc%udef 
              endif
            enddo
          enddo
    
       case( "bilinear", "budget-bilinear", "neighbor"  )  
          !lb1d = .false. ! SY: See next line
          lb1d = .true. ! SY. NOTE. Any NODATA rain values also correspondingly fixed to valid 0 below  
          do r=1,NINT(RFE2Daily_struc(n)%gridDesci(3))
             do c=1,NINT(RFE2Daily_struc(n)%gridDesci(2))
                t = c+(r-1)*NINT(RFE2Daily_struc(n)%gridDesci(2))
                rain1d(t) = rain2d(c,r) 
                if( rain1d(t) .LT. (-999.0+1) ) then ! SY
                  rain1d(t) = 0.0
                  !lb1d(t) = .true. ! SY
                endif
             enddo
          enddo

       case default
          write(LIS_logunit,*)"[ERR] RFE2Daily reprojection choice should be: "
          write(LIS_logunit,*) "  average, neighbor, bilinear,or budget-bilinear"
          write(LIS_logunit,*) "Program stopping ... "
          call LIS_endrun()
      end select
    
      call reproject_RFE2Daily(n,findex,&
                 NINT(RFE2Daily_struc(n)%gridDesci(2))* &
                 NINT(RFE2Daily_struc(n)%gridDesci(3)), &
                 rain1d, lb1d, LIS_rc%gridDesc(n,:),    &
                 LIS_rc%lnc(n),LIS_rc%lnr(n), varfield)
  
      do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
          if ( (LIS_domain(n)%gindex(c,r).ne.-1) .AND. &
               (varfield(c,r) .ge. 0) ) then       ! SY: leaves out -1 and NODATA values 
            IF (varfield(c,r) .GT. 1000.0) THEN    ! 1000 mm threshold
              varfield(c,r) = 1000.0 
            ENDIF
            if(order.eq.1) then 
               RFE2Daily_struc(n)%metdata1(kk,1,LIS_domain(n)%gindex(c,r)) = &
                    varfield(c,r)
            elseif(order.eq.2) then 
               RFE2Daily_struc(n)%metdata2(kk,1,LIS_domain(n)%gindex(c,r)) = &
                    varfield(c,r)
            endif
          endif 
        end do
      end do
  
      ferror_RFE2Daily = 0
      if(LIS_masterproc) write(LIS_logunit,*) &
                   "[INFO] Read RFE2Daily data file: ", trim(fname)
    else
      if(LIS_masterproc) write(LIS_logunit,*) &
                   "[WARN] Cannot open RFE2Daily precipitation file ", trim(fname)
      ferror_RFE2Daily = 1
    endif
    
    call LIS_releaseUnitNumber(ftn) 

  ELSE ! SY: For IF (file_exists .EQ. .TRUE.) THEN

    if(LIS_masterproc) write(LIS_logunit,*) &
              "[WARN] Missing RFE2Daily precipitation file ", fname
    ferror_RFE2Daily = 1
    
  ENDIF ! SY: For IF (file_exists .EQ. .TRUE.) THEN

  deallocate(lb1d)
  deallocate(lb2d)
  deallocate(rain1d)
  deallocate(rain2d)

end subroutine readprecip_RFE2Daily
