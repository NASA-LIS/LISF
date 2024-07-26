!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: read_GLIMS_glacierfraction
!  \label{read_GLIMS_glacierfraction}
!
! !REVISION HISTORY:
!  30 Mar 2018: Sujay Kumar; Initial Specification
!  23 Jun 2020: Mahdi Navari; Computing glacier fraction
!
! !INTERFACE:
subroutine read_GLIMS_glacierfraction(n, glacier_frac )

! !USES:
  use LDT_coreMod,   only : LDT_rc, LDT_localPet
  use LDT_gridmappingMod
  use LDT_logMod
  use LDT_glacierMod
  use LDT_fileIOMod
  use LDT_glacierFractionMod

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
  use netcdf
#endif

  implicit none

! !ARGUMENTS: 
  integer, intent(in)  :: n
!  real,    intent(out) :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
  real, intent(out)     :: glacier_frac(LDT_rc%lnc(n),LDT_rc%lnr(n),1)
! !DESCRIPTION:
!  This subroutine reads the landmask data and compute the 
!   glacier fraction.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of nest
!   \item[glacier_frac]
!    glacier fraction for the region of interest
!   \end{description}
!
!EOP      
  integer, parameter :: nr = 18000
  integer :: ftn, ios1,maskid
  logical :: file_exists
  integer :: c, r, t, i, line
  integer :: glpnc, glpnr             ! Parameter (global) total columns and rows
  integer :: subpnc, subpnr           ! Parameter subsetted columns and rows
  real    :: subparam_gridDesc(20)    ! Input parameter grid desc array
  real    :: cnt_0mask_1lc
  real    :: cnt_1mask_0lc
  integer, allocatable :: lat_line(:,:)
  integer, allocatable :: lon_line(:,:)
  real,    allocatable :: read_inputparm(:,:)  ! Read input parameter
  integer :: mi                       ! Total number of input param grid array points
  integer :: mo                       ! Total number of output LIS grid array points
  real,    allocatable  :: gi1(:)      ! Input parameter 1d grid
  logical*1,allocatable :: li1(:)      ! Input logical mask (to match gi)

  real                  :: go1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output lis 1d grid
  logical*1             :: lo1(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output logical mask (to match go)
  real                  :: param_gridDesc(20)       ! Input parameter grid desc array

! MN                                                                                                                
  real                  :: localmask(LDT_rc%lnc(n),LDT_rc%lnr(n),1)                                        
  integer, allocatable  :: n11(:)     ! Array that maps the location of each input grid                             
                                                      ! point in the output grid.                                 
  integer               :: pixels_pergrid(LDT_rc%lnc(n)*LDT_rc%lnr(n))                    
!  integer, parameter    :: num_bins = 1   
  integer               :: num_bins
  real                  :: go2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output lis 1d grid
  logical*1             :: lo2(LDT_rc%lnc(n)*LDT_rc%lnr(n))  ! Output logical mask (to match go)                                                  
  character(50)         :: gridtransform_opt
!  real      :: glacierfrac_cnt(LDT_rc%lnc(n),LDT_rc%lnr(n), num_bins)
!  real      :: isum
! ______________________________________________________________

!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = 36000 !subpnc
   param_gridDesc(3)  = 18000 !subpnr
   param_gridDesc(4)  = -89.995 !  + lat_line(1,1)*0.01   ! LL lat
   param_gridDesc(5)  = -179.995 ! + lon_line(1,1)*0.01    ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = 89.995  !+ (lat_line(1,1) + subpnc)*0.01   ! UR lat
   param_gridDesc(8)  = 179.856 !+ (lon_line(1,1) + subpnr)*0.01   ! UR lon
   param_gridDesc(9)  = 0.01     ! dy
   param_gridDesc(10) = 0.01      ! dx
   param_gridDesc(20) = 64

! ______________________________________________________________

!- Check for and open landmask file:
   inquire(file=trim(LDT_glacierfrac_struc(n)%glacierfracfile), exist=file_exists)
   if( file_exists ) then 
      write(LDT_logunit,*)"[INFO] GLIMS mask -- Reading ",trim(LDT_glacierfrac_struc(n)%glacierfracfile)
   else
      write(LDT_logunit,*) "[ERR] GLIMS map: ",trim(LDT_rc%glaciermask(n))," does not exist."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

! -------------------------------------------------------------------
!    PREPARE SUBSETTED PARAMETER GRID FOR READING IN NEEDED DATA
! -------------------------------------------------------------------
!- Map Parameter Grid Info to LIS Target Grid/Projection Info -- 
   subparam_gridDesc = 0.
   call LDT_RunDomainPts( n, LDT_glacierfrac_struc(n)%glacierfrac_proj, param_gridDesc(:), &
            glpnc, glpnr, subpnc, subpnr, subparam_gridDesc, lat_line, lon_line )

   allocate( read_inputparm(subpnc, subpnr) )
   read_inputparm = 0.
   
   ! -------------------------------------------------------------------

#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
   call LDT_verify(nf90_open(path=LDT_glacierfrac_struc(n)%glacierfracfile,&
        mode=nf90_nowrite,ncid=ftn),&
        'nf90_open file failed in read_GLIMS_glacierfraction')
   
  ! call LDT_verify(nf90_open(path=LDT_rc%glaciermask(n),&
  !      mode=nf90_nowrite,ncid=ftn),&
  !      'nf90_open file failed in read_GLIMS_glacierfraction')

   call LDT_verify(nf90_inq_varid(ftn,"glaciermask",maskid),&
        'nf90_inq_varid failed in read_GLIMS_glacierfraction')

   call LDT_verify(nf90_get_var(ftn,maskid,read_inputparm,&
        start=(/lon_line(1,1),nr-(lat_line(1,1)+subpnr)+1/),&
        count=(/subpnc,subpnr/)),&
        'nf90_get_var failed in read_GLIMS_glacierfraction')
   call LDT_verify(nf90_close(ftn))

#endif


!- Enter spatial downscaling/upscaling options to bring the data
!  to the target domain ...
    mi = subpnc*subpnr
    allocate( gi1(mi), li1(mi), n11(mi) )
    gi1 = LDT_rc%udef
    li1 = .false.
    mo = LDT_rc%lnc(n)*LDT_rc%lnr(n)
    lo1 = .false.

!- Assign 2-D array to 1-D for aggregation routines:
   do r = 1, subpnr
      do c = 1, subpnc
         i = c+(r-1)*subpnc
         gi1(i) = read_inputparm(c,r)
         if( gi1(i) .ge. 0. )  li1(i) = .true. 
      enddo
   enddo
#if 0 
   param_gridDesc = 0
!- Set parameter grid array inputs:
   param_gridDesc(1)  = 0.         ! Latlon
   param_gridDesc(2)  = subpnc
   param_gridDesc(3)  = subpnr
   param_gridDesc(4)  = -89.995 + lat_line(1,1)*0.01   ! LL lat
   param_gridDesc(5)  = -179.995 + lon_line(1,1)*0.01    ! LL lon
   param_gridDesc(6)  = 128
   param_gridDesc(7)  = -89.995 + (lat_line(1,1) + subpnc)*0.01   ! UR lat
   param_gridDesc(8)  = 179.856 + (lon_line(1,1) + subpnr)*0.01   ! UR lon
   param_gridDesc(9)  = 0.01     ! dy
   param_gridDesc(10) = 0.01      ! dx
   param_gridDesc(20) = 64
#endif
   deallocate( lat_line, lon_line )


# if 0 
!- Transform parameter from original grid to LIS output grid:
   call LDT_transform_paramgrid(n, LDT_glacierfrac_struc(n)%glacierfrac_gridtransform, &
        param_gridDesc, mi, 1, gi1, li1, mo, go1, lo1 )

!- Convert 1D to 2D grid output arrays:
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n); i = i + 1
         if( go1(i) <=  LDT_rc%gridcell_glacier_frac(n)) then  ! LDT_rc%gridcell_glacier_frac(n) this is cutoff value   
           localmask(c,LDT_rc%lnr(n)-r+1,1) = LDT_rc%udef
        else
           localmask(c,LDT_rc%lnr(n)-r+1,1) = 1.0
        end if
      enddo                                                                                                         
   enddo                                                                                                            
                                                                                                                     
# endif

!-------------------------------------------------------------------------------  
! Third method 
!-------------------------------------------------------------------------------
  ! Compute glacier fraction ! MN                                                                                    

 gridtransform_opt = LDT_glacierfrac_struc(n)%glacierfrac_gridtransform

 if( gridtransform_opt == "average" )  then


  !- Create mapping between parameter domain and LIS grid domain: --> n11                                                  
     call upscaleByAveraging_input( subparam_gridDesc, &
                             LDT_rc%gridDesc(n,:), mi, mo, n11 )


     call upscaleByAveraging( mi, mo, LDT_rc%udef, n11,  &
                                 li1, gi1, lo2, go2 )


  ! Compute glacier fraction ! MN        
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n); i = i + 1
         if( go2(i) <=  0.0001) then  ! LDT_rc%gridcell_glacier_frac(n) this is cutoff value   
           glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = 0
        else
          ! print *, 'i, go2(i)' , i, go2(i)
           glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = go2(i)
        end if
      enddo
   enddo
 else 
         write(LDT_logunit,*) "[ERR] This spatial transformation option ("//trim(gridtransform_opt)//") "
         write(LDT_logunit,*) " is not currently supported."
         write(LDT_logunit,*) " Program stopping ..."
         call LDT_endrun
endif
!-------------------------------------------------------------------------------- 



# if 0 
!-------------------------------------------------------------------------------
! Secod method 
! NOTE: this more likely to work if we add this to the readGlaciermask.F90
!-------------------------------------------------------------------------------
  ! Compute glacier fraction ! MN        
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n); i = i + 1
         if( go1(i) <=  LDT_rc%gridcell_glacier_frac(n)) then  ! LDT_rc%gridcell_glacier_frac(n) this is cutoff value   
           glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = 0
        else
!print *, 'i, go1(i)' , i, go1(i)
           glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = go1(i)
        end if
      enddo
   enddo

!-------------------------------------------------------------------------------
# endif 

# if 0
!------------------------------------------------------------------------------- 
! First method  , this works but need to be verified
!-------------------------------------------------------------------------------
  ! Compute glacier fraction ! MN                                                                                    
                                                                                                                     
  !- Create mapping between parameter domain and LIS grid domain: --> n11                                                  
     call upscaleByAveraging_input( subparam_gridDesc, &                                                             
                             LDT_rc%gridDesc(n,:), mi, mo, n11 )                                                     
  

  num_bins =  LDT_glacierfrac_struc(n)%glacierfrac%num_bins
  !- Calculate total counts of fine grid cells in each coarse gridcell: 
     call upscaleByCnt( mi, mo, num_bins, LDT_rc%udef, n11, li1, gi1, &
                      lo2, go2 )


  !- Estimate number of pixels per gridcell (coarse domain):    
   pixels_pergrid = 0. 
   do i = 1, mi
      if(li1(i)) then
         if( n11(i) .ne. 0 .and. gi1(i) .ne. LDT_rc%udef ) then
              pixels_pergrid(n11(i)) = pixels_pergrid(n11(i)) + 1.0
         endif
      endif
   enddo

     glacier_frac = 0.
     i = 0
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n);
            i = i + 1
            if( go2(i) > 0 ) then                                                     
             glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = go2(i)  / pixels_pergrid(i)
            else                                                                                                      
             glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = 0.0                                                                 
            end if                                                                  
        enddo
     enddo
!-------------------------------------------------------------------------------- 
# endif

! this method results in 0 ,1 
# if 0 

   glacierfrac_cnt = 0.
   i = 0
   do r = 1, LDT_rc%lnr(n)
      do c = 1, LDT_rc%lnc(n)
         i = i + 1
         do t = 1, num_bins
            glacierfrac_cnt(c,LDT_rc%lnr(n)-r+1,t) = go2(i)
         end do
      enddo
   enddo
   

     glacier_frac = 0.
     do r = 1, LDT_rc%lnr(n)
        do c = 1, LDT_rc%lnc(n);
       ! Calculate gridpoint column total:
         isum = sum(glacierfrac_cnt(c,r,1:num_bins))
       ! Estimate 2-D fraction for just irrigation:
         if( isum > 0. ) then           
             glacier_frac(c,r,1) = glacierfrac_cnt(c,r,1) / isum
            else
             glacier_frac(c,r,1) = 0.0
            end if
        enddo
     enddo
# endif 




#if 0                                                                                                                
  !- Estimate number of pixels per gridcell (coarse domain):                                                         
     do i = 1,  mi ! inpts                                                                                           
        if( n11(i) .ne. 0 ) then                                                                                     
           pixels_pergrid(n11(i)) = pixels_pergrid(n11(i)) + 1                                                       
        endif                                                                                                        
     enddo                                                                                                           
                                                                                                                     
     i = 0                                                                                                           
     do r = 1, LDT_rc%lnr(n)                                                                                         
        do c = 1, LDT_rc%lnc(n);                                                                                     
            i = i + 1                                                                                                
  !         if( go1(i) <=  LDT_rc%gridcell_glacier_frac(n)) then                                                     
             glacier_frac(c,LDT_rc%lnr(n)-r+1,1) = go1(i) / pixels_pergrid(n11(i))                                   
  !        else                                                                                                      
  !           localmask(c,LDT_rc%lnr(n)-r+1,1) = 1.0                                                                 
  !        end if                                                                  
     enddo
  enddo
#endif

  deallocate( li1, gi1, n11 )

end subroutine read_GLIMS_glacierfraction
