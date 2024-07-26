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
! !ROUTINE: readprecip_RFE2gdas
! \label{readprecip_RFE2gdas}
!
! !REVISION HISTORY:
!  30 May 2010; Soni Yatheendradas, Initial LIS version for FEWSNET
!
! !INTERFACE:
subroutine readprecip_RFE2gdas( n, kk,fname, month, findex, order, ferror_RFE2gdas, filehr)
! !USES:
  use LIS_coreMod,        only : LIS_rc, LIS_domain, LIS_masterproc
  use LIS_logMod,         only : LIS_logunit, LIS_getNextUnitNumber, &
          LIS_releaseUnitNumber, LIS_endrun
  use RFE2gdas_forcingMod, only : RFE2gdas_struc

  implicit none

! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: kk
  character(len=*)    :: fname
  integer, intent(in) :: month
  integer, intent(in) :: findex
  integer, intent(in) :: order
  integer             :: ferror_RFE2gdas
  integer             :: filehr
! 
! !DESCRIPTION:
!  For the given time, reads the RFE2gdas data
!  and interpolates to the LIS domain.
!  
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[fname]
!    name of the 6 hour RFE2gdas.0 file
!  \item[findex]
!    index of the supplemental forcing source
!  \item[ferror\_RFE2gdas]
!    flag to indicate success of the call (=0 indicates success)
!  \item[filehr]
!    current file hour
!  \end{description}
!  
!  The routines invoked are:
!  \begin{description}
!  \item[reproject\_RFE2gdas](\ref{reproject_RFE2gdas}) \newline
!    Upscales/interpolates the RFE2gdas data
!  \end{description}
!EOP
!
!==== Local Variables=======================

  integer               :: c,r,t
  logical*1,allocatable :: lb1d(:)
  logical*1,allocatable :: lb2d(:,:)
  real, allocatable     :: rain1d(:)
  REAL, ALLOCATABLE     :: rain2d(:,:)
  real, dimension(LIS_rc%lnc(n), LIS_rc%lnr(n)) :: varfield ! reprojected arrray
  integer               ::  ftn, ios, ftn2, ftn3

!=== End Variable Definition =======================

  varfield = LIS_rc%udef

  allocate(lb1d(NINT(RFE2gdas_struc(n)%gridDesci(2))*&
       NINT(RFE2gdas_struc(n)%gridDesci(3))))
  allocate(lb2d(NINT(RFE2gdas_struc(n)%gridDesci(2)),&
       NINT(RFE2gdas_struc(n)%gridDesci(3))))
  allocate(rain1d(NINT(RFE2gdas_struc(n)%gridDesci(2))*&
       NINT(RFE2gdas_struc(n)%gridDesci(3))))
  allocate(rain2d(NINT(RFE2gdas_struc(n)%gridDesci(2)),&
       NINT(RFE2gdas_struc(n)%gridDesci(3))))

  ftn = LIS_getNextUnitNumber()

  open(unit=ftn,file=fname, access='direct', &
       recl=NINT(RFE2gdas_struc(n)%gridDesci(2))* &
       NINT(RFE2gdas_struc(n)%gridDesci(3))*4, &
       iostat=ios)
  
  ! Read in precip field, if file opened properly:
  if( ios .eq. 0 ) then 

     read (ftn,rec=1)  rain2d
     close(ftn)

  !- Upscaling/averaging 2D to 1D assignment:
     select case( LIS_rc%met_interp(findex) )

       case( "average" )   ! Upscaling 
         lb2d = .false.
         lb1d = .false.
#if 0
! Soni's old code
         do r=1,NINT(RFE2gdas_struc(n)%gridDesci(3))
            do c=1,NINT(RFE2gdas_struc(n)%gridDesci(2))
               if(rain2d(c,r).GT.(-999.0+1)) then
                  lb2d(c,r) = .true.
               else
                  rain2d(c,r)=LIS_rc%udef 
               endif
            enddo
         enddo
#endif
         do r=1,NINT(RFE2gdas_struc(n)%gridDesci(3))
            do c=1,NINT(RFE2gdas_struc(n)%gridDesci(2))
               t = c+(r-1)*NINT(RFE2gdas_struc(n)%gridDesci(2))
               rain1d(t) = rain2d(c,r)
               if( rain1d(t) .GT. (-999.0+1) ) then ! SY
                 lb1d(t) = .true.
               else
                 rain1d(t)=LIS_rc%udef
               endif
            enddo
         enddo


      case( "bilinear", "budget-bilinear", "neighbor"  )

         !lb1d = .false. ! SY: See next line
         lb1d = .true. ! SY. NOTE. Any NODATA rain values also correspondingly fixed to valid 0 below  
         do r=1,NINT(RFE2gdas_struc(n)%gridDesci(3))
            do c=1,NINT(RFE2gdas_struc(n)%gridDesci(2))
               t = c+(r-1)*NINT(RFE2gdas_struc(n)%gridDesci(2))
               rain1d(t) = rain2d(c,r) 
               if( rain1d(t) .LT. (-999.0+1) ) then ! SY
                  rain1d(t) = 0.0
                  !lb1d(t) = .true. ! SY
               endif
            enddo
         enddo

     end select

     CALL reproject_RFE2gdas(n,findex,month,&
          NINT(RFE2gdas_struc(n)%gridDesci(2))* &
          NINT(RFE2gdas_struc(n)%gridDesci(3)),&
          rain1d,rain2d,lb1d,lb2d, LIS_rc%gridDesc(n,:),&
          LIS_rc%lnc(n),LIS_rc%lnr(n), varfield)

     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if ( (LIS_domain(n)%gindex(c,r).ne.-1) .AND. &
                (varfield(c,r) .ge. 0) ) then ! SY: leaves out -1 and NODATA values 
              if(order.eq.2) then 
                 RFE2gdas_struc(n)%metdata2(kk,1,LIS_domain(n)%gindex(c,r)) = &
                      varfield(c,r)
              endif
           endif
        end do
     end do

     ferror_RFE2gdas = 0
     if(LIS_masterproc) write(LIS_logunit,*) &
        "[INFO] Reading RFE2gdas data file: ", fname
  else
     if(LIS_masterproc) write(LIS_logunit,*) &
        "[WARN] Missing RFE2gdas data file: ", fname
     ferror_RFE2gdas = 1
  endif

  call LIS_releaseUnitNumber(ftn) 

  deallocate(lb1d)
  deallocate(lb2d)
  deallocate(rain1d)
  deallocate(rain2d)

end subroutine readprecip_RFE2gdas
