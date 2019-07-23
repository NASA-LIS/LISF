!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
#include "LVT_misc.h"
!BOP
! 
! !MODULE: LVT_domainMod
! \label(LVT_domanMod)
!
! !INTERFACE:
module LVT_domainMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_domain_pluginMod
  use LVT_LMLCMod
  use LVT_topoMod
  use LVT_soilsMod
  use LVT_mpiMod
  use LVT_logMod
  use LVT_statsDataMod
  use map_utils
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif


  implicit none
  PRIVATE
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! The code in this file provides interfaces to manages different running
! domain implementations
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_domainInit          !initialize specified domains
  public :: LVT_domain_setup
  public :: LVT_compute_domainIndices
  public :: LVT_quilt_domain         !generate quilted domains
  public :: LVT_readDataMask

!EOP

contains
!BOP
! 
! !ROUTINE: LVT_domainInit
! \label{LVT_domainInit}
!
! !INTERFACE: 
  subroutine LVT_domainInit()
! 
! !USES:

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine issues calls to set up the domain (grid and map projection)
!  information. The domain over which LVT analysis is run is generated first
!  If LIS output is one of the datastreams being compared, then the domain(s) of 
!  LIS simulation(s) are generated. 
! 
!   The routines invoked are: 
!   \begin{description}
!    \item[LVT\_domain\_plugin] (\ref{LVT_domain_plugin}) \newline
!      invokes the registry that includes the supported domains (map projections)
!    \item[readinput] (\ref{readinput}) \newline
!      reads the configuration options for the specified map projection
!    \item[create_LVT_gridspace] (\ref{create_LVT_gridspace}) \newline
!      creates the LVT gridspace based on the provided landmask
!    \item[create_LISoutput_domain] (\ref{create_LISoutput_domain}) \newline
!      creates the LIS output domain (if 'LIS output' is one of the datastreams)
!    \item[LVT\_readDataMask] (\ref{LVT_readDataMask}) \newline
!      reads the external data to be used for masking    
!   \end{description}
! 
!EOP

    integer :: c,r
    real    :: locallat, locallon
    integer :: status
    
    call LVT_domain_plugin
    call readinput(trim(LVT_rc%domain)//char(0))   
    call create_LVT_gridspace()

    if(LVT_rc%lis_output_obs) then
       write(LVT_logunit,*) '[INFO] '
       write(LVT_logunit,*) '[INFO] Setting up the LIS domain '
       call create_LISoutput_domain()
       write(LVT_logunit,*) '[INFO] Finished the LIS domain setup '
       write(LVT_logunit,*) '[INFO] '
    else

       LVT_domain%minLat = 200.0
       LVT_domain%minLon = 200.00
       LVT_domain%maxLat = -200.0
       LVT_domain%maxLon = -200.00
       
       do r=1,LVT_rc%lnr
          do c=1,LVT_rc%lnc
             call ij_to_latlon(LVT_domain%lvtproj,float(c),float(r),&
                  locallat,locallon)
             
             if(localLat.lt.LVT_domain%minLat) LVT_domain%minLat = localLat
             if(localLat.gt.LVT_domain%maxLat) LVT_domain%maxLat = localLat
             if(localLon.lt.LVT_domain%minLon) LVT_domain%minLon = localLon
             if(localLon.gt.LVT_domain%maxLon) LVT_domain%maxLon = localLon
             
          enddo
       enddo
       
    endif

    allocate(LVT_stats%datamask(LVT_rc%lnc, LVT_rc%lnr))    
    LVT_stats%datamask = 0.0  !initially set all to false.    
    call LVT_readDataMask()

  end subroutine LVT_domainInit


!BOP
! 
! !ROUTINE: LVT_quilt_domain
! \label{LVT_quilt_domain}
!
! !INTERFACE: 
  subroutine LVT_quilt_domain(nc, nr )
! 
! !USES:     

!
! !INPUT PARAMETERS: 
    integer, intent(in)  :: nc
    integer, intent(in)  :: nr
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine generates the quilted domain extents and sizes based on the 
! processor layout and the size of the specified halos. 
! 
!  The arguments are: 
!  \begin{description}
!    \item[nc]
!      size of the east-west dimension to be decomposed
!    \item[nr]
!       size of the north-south dimension to be decomposed
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[LVT\_mpDecomp] (\ref{LVT_mpDecomp}) \newline
!      decomposes the domain based on the processor layout
!   \end{description}
!EOP

    integer              :: ips, ipe, jps, jpe
    integer              :: Px, Py, P
    integer              :: mytask_x, mytask_y
    integer              :: i, j,ierr

     mytask_x = mod( LVT_localPet , LVT_rc%npesx)
     mytask_y = LVT_localPet / LVT_rc%npesx

     j = 1
     ips = -1
     do i=1, nc+1
        call LVT_mpDecomp(i,j,1,nc+1,1,nr+1,LVT_rc%npesx, LVT_rc%npesy,Px,Py,P)
        
        if(Px.eq. mytask_x) then 
           ipe = i 
           if(ips.eq.-1) ips = i 
        endif
     enddo
     
     i = 1
     jps = -1
     do j=1, nr+1
        call LVT_mpDecomp(i,j,1,nc+1,1,nr+1,LVT_rc%npesx,LVT_rc%npesy, Px, Py, P)
        if(Py.eq.mytask_y) then 
           jpe = j
           if(jps.eq.-1) jps = j 
        endif
     enddo
     LVT_ews_halo_ind(LVT_localPet+1) = max(ips-LVT_rc%haloy, 1)
     LVT_ewe_halo_ind(LVT_localPet+1) = min(ipe+LVT_rc%haloy, nc)
     LVT_nss_halo_ind(LVT_localPet+1) = max(jps-LVT_rc%halox, 1)
     LVT_nse_halo_ind(LVT_localPet+1) = min(jpe+LVT_rc%halox, nr)

     LVT_ews_ind(LVT_localPet+1) = ips
     LVT_ewe_ind(LVT_localPet+1) = min(ipe, nc)
     LVT_nss_ind(LVT_localPet+1) = jps
     LVT_nse_ind(LVT_localPet+1) = min(jpe, nr)

     LVT_rc%lnc = LVT_ewe_halo_ind(LVT_localPet+1)-&
          LVT_ews_halo_ind(LVT_localPet+1)+1
     LVT_rc%lnr = LVT_nse_halo_ind(LVT_localPet+1)-&
          LVT_nss_halo_ind(LVT_localPet+1)+1

     LVT_rc%lnc_red= (LVT_ewe_ind(LVT_localPet+1)-&
          LVT_ews_ind(LVT_localPet+1))+1
     LVT_rc%lnr_red= (LVT_nse_ind(LVT_localPet+1)-&
          LVT_nss_ind(LVT_localPet+1))+1

     write(LVT_logunit,*) '[INFO] Local domain :(',LVT_rc%lnc,LVT_rc%lnr,')'
     write(LVT_logunit,*) '[INFO] Local domain without halo',':(',&
          LVT_rc%lnc_red,LVT_rc%lnr_red,')'
     write(LVT_logunit,*) '[INFO] Running domain',':(',nc,nr,')'

     LVT_deltas(LVT_localPet) = LVT_rc%lnc_red*LVT_rc%lnr_red

#if (defined SPMD)
     call MPI_ALLGATHER(LVT_deltas(LVT_localPet),1,MPI_INTEGER,LVT_deltas(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_ews_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_ews_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_nss_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_nss_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_ewe_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_ewe_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_nse_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_nse_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)

     call MPI_ALLGATHER(LVT_ews_halo_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_ews_halo_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_nss_halo_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_nss_halo_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_ewe_halo_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_ewe_halo_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_nse_halo_ind(LVT_localPet+1),1,MPI_INTEGER,LVT_nse_halo_ind(:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
#endif   
     if(LVT_masterproc) then 
        LVT_offsets(0) = 0 
        do i=1,LVT_npes-1
           LVT_offsets(i) = LVT_offsets(i-1)+LVT_deltas(i-1)
        enddo
     end if
#if (defined SPMD)
     call MPI_BCAST(LVT_offsets(:), LVT_npes, MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
#endif

   end subroutine LVT_quilt_domain

!BOP
! 
! !ROUTINE: LVT_quilt_lis_domain
! \label{LVT_quilt_lis_domain}
!
! !INTERFACE: 
  subroutine LVT_quilt_lis_domain(source, nc, nr )
! 
! !USES:     

!
! !INPUT PARAMETERS: 
    integer, intent(in)  :: source
    integer, intent(in)  :: nc
    integer, intent(in)  :: nr
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! This routine generates the quilted domain extents and sizes based on the 
! processor layout and the size of the specified halos. 
! 
!  The arguments are: 
!  \begin{description}
!    \item[source]
!      index of the 'LIS output' data stream
!    \item[nc]
!      size of the LIS east-west dimension to be decomposed
!    \item[nr]
!       size of the LIS north-south dimension to be decomposed
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!    \item[LVT\_mpDecomp] (\ref{LVT_mpDecomp}) \newline
!      decomposes the domain based on the processor layout
!   \end{description}
!
!EOP

    integer              :: ips, ipe, jps, jpe
    integer              :: Px, Py, P
    integer              :: mytask_x, mytask_y
    integer              :: i, j,ierr

     mytask_x = mod( LVT_localPet , LVT_rc%npesx)
     mytask_y = LVT_localPet / LVT_rc%npesx

     j = 1
     ips = -1
     do i=1, nc+1
        call LVT_mpDecomp(i,j,1,nc+1,1,nr+1,LVT_rc%npesx, LVT_rc%npesy,Px,Py,P)
        
        if(Px.eq. mytask_x) then 
           ipe = i 
           if(ips.eq.-1) ips = i 
        endif
     enddo
     
     i = 1
     jps = -1
     do j=1, nr+1
        call LVT_mpDecomp(i,j,1,nc+1,1,nr+1,LVT_rc%npesx,LVT_rc%npesy, Px, Py, P)
        if(Py.eq.mytask_y) then 
           jpe = j
           if(jps.eq.-1) jps = j 
        endif
     enddo
     LVT_lis_ews_halo_ind(source,LVT_localPet+1) = max(ips-LVT_rc%haloy, 1)
     LVT_lis_ewe_halo_ind(source,LVT_localPet+1) = min(ipe+LVT_rc%haloy, nc)
     LVT_lis_nss_halo_ind(source,LVT_localPet+1) = max(jps-LVT_rc%halox, 1)
     LVT_lis_nse_halo_ind(source,LVT_localPet+1) = min(jpe+LVT_rc%haloy, nr)

     LVT_lis_ews_ind(source,LVT_localPet+1) = ips
     LVT_lis_ewe_ind(source,LVT_localPet+1) = min(ipe, nc)
     LVT_lis_nss_ind(source,LVT_localPet+1) = jps
     LVT_lis_nse_ind(source,LVT_localPet+1) = min(jpe, nr)

     LVT_LIS_rc(source)%lnc = LVT_lis_ewe_halo_ind(source,LVT_localPet+1)-&
          LVT_lis_ews_halo_ind(source,LVT_localPet+1)+1
     LVT_LIS_rc(source)%lnr = LVT_lis_nse_halo_ind(source,LVT_localPet+1)-&
          LVT_lis_nss_halo_ind(source,LVT_localPet+1)+1

     LVT_LIS_rc(source)%lnc_red= (LVT_lis_ewe_ind(source,LVT_localPet+1)-&
          LVT_lis_ews_ind(source,LVT_localPet+1))+1
     LVT_LIS_rc(source)%lnr_red= (LVT_lis_nse_ind(source,LVT_localPet+1)-&
          LVT_lis_nss_ind(source,LVT_localPet+1))+1

     write(LVT_logunit,*) '[INFO] LIS domain',':(',LVT_LIS_rc(source)%lnc,LVT_LIS_rc(source)%lnr,')'
     write(LVT_logunit,*) '[INFO] LIS domain without halo',':(',&
          LVT_LIS_rc(source)%lnc_red,LVT_LIS_rc(source)%lnr_red,')'
     write(LVT_logunit,*) '[INFO] running domain',':(',nc,nr,')'

     LVT_lis_deltas(source,LVT_localPet) = LVT_LIS_rc(source)%lnc_red*LVT_LIS_rc(source)%lnr_red

#if (defined SPMD)
     call MPI_ALLGATHER(LVT_lis_deltas(source,LVT_localPet),1,MPI_INTEGER,LVT_lis_deltas(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_ews_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_ews_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_nss_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_nss_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_ewe_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_ewe_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_nse_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_nse_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)

     call MPI_ALLGATHER(LVT_lis_ews_halo_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_ews_halo_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_nss_halo_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_nss_halo_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_ewe_halo_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_ewe_halo_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
     call MPI_ALLGATHER(LVT_lis_nse_halo_ind(source,LVT_localPet+1),1,MPI_INTEGER,LVT_lis_nse_halo_ind(source,:),1,MPI_INTEGER,&
          MPI_COMM_WORLD,ierr)
#endif   
     if(LVT_masterproc) then 
        LVT_lis_offsets(source,0) = 0 
        do i=1,LVT_npes-1
           LVT_lis_offsets(source,i) = LVT_lis_offsets(source,i-1)+LVT_lis_deltas(source,i-1)
        enddo
     end if
#if (defined SPMD)
     call MPI_BCAST(LVT_lis_offsets(source,:), LVT_npes, MPI_INTEGER,0, MPI_COMM_WORLD, ierr)
#endif

   end subroutine LVT_quilt_lis_domain

!BOP
! 
! !ROUTINE: create_LVT_gridspace
! \label{create_LVT_gridspace}
! 
! !INTERFACE: 
   subroutine create_LVT_gridspace()
!
! !DESCRIPTION:
! 
!  This routine reads an input landmask and creates the LVT gridspace
!  (number of grid points and the mapping of the LVT 1-d vector space to 
!  the 2-d domain) 
!EOP
    real          :: locallat, locallon
    integer       :: c, r, t, m
    integer       :: iv, it, ie, is, ia
    integer       :: vegt
    integer       :: soilt
    real          :: sand, clay, silt
    real          :: elev, slope, aspect
    integer       :: sf_index
    integer       :: gnc, gnr
    integer       :: kk
    integer       :: kk_sf(LVT_rc%max_model_types)
    real, allocatable :: sumv(:)
    integer       :: ntiles_surface
    integer       :: npatch_surface
    integer       :: ntiles_soilt
    integer       :: ntiles_soilf
    integer       :: ntiles_soil
    integer       :: ntiles_elev
    integer       :: ntiles_slope
    integer       :: ntiles_aspect
    integer       :: gid
    real          :: temp
    integer       :: soilf_index
    integer       :: elev_index
    integer       :: slope_index
    integer       :: aspect_index
    real          :: mask(LVT_rc%input_lnc,LVT_rc%input_lnr)
    real          :: mask_in(LVT_rc%input_lnc*LVT_rc%input_lnr)
    real          :: mask_out(LVT_rc%lnc*LVT_rc%lnr)
    logical*1     :: li(LVT_rc%input_lnc*LVT_rc%input_lnr)
    logical*1     :: lo(LVT_rc%lnc*LVT_rc%lnr)
    integer       :: count1
    integer       :: ios, ftn,maskid
    logical       :: file_exists, interp_flag


#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    inquire(file=LVT_rc%paramfile, exist=file_exists)
    if(file_exists) then 
       
       ios = nf90_open(path=LVT_rc%paramfile, &
            mode=NF90_NOWRITE,ncid=ftn)
       call LVT_verify(ios,'Error in nf90_open in LIS_domainMod')
              
       ios = nf90_inq_varid(ftn,'LANDMASK',maskid)
       call LVT_verify(ios,'LANDMASK field not found in the LVT param file')
       
       ios = nf90_get_var(ftn,maskid,mask)
       call LVT_verify(ios,'Error in nf90_get_var for LANDMASK in readinput_latlon')

       ios = nf90_close(ftn)
       call LVT_verify(ios,'Error in nf90_close')
    endif
#endif    
    
    do r=1,LVT_rc%input_lnr
       do c=1,LVT_rc%input_lnc
          mask_in(c+(r-1)*LVT_rc%input_lnc) = mask(c,r)
          if(mask(c,r).gt.0) then 
             li(c+(r-1)*LVT_rc%input_lnc) = .true. 
          else
             li(c+(r-1)*LVT_rc%input_lnc) = .false.                 
          endif
       enddo
    enddo

    if(LVT_rc%domain.eq."latlon") then 
       if(LVT_rc%input_gridDesc(10).gt.LVT_rc%gridDesc(10)) then !interpolate
          interp_flag = .true. 
       else
          interp_flag = .false. 
       endif
    elseif(LVT_rc%domain.eq."lambert") then 
       interp_flag = .false. 
    endif

    if(LVT_rc%spTransformAnlysDomain.eq."none") then 
       if(interp_flag) then 
          allocate(LVT_rc%rlat_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%rlon_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%n11_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%n12_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%n21_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%n22_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%w11_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%w12_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%w21_dn(LVT_rc%lnc*LVT_rc%lnr))
          allocate(LVT_rc%w22_dn(LVT_rc%lnc*LVT_rc%lnr))
          
          call bilinear_interp_input(LVT_rc%input_gridDesc,&
               LVT_rc%gridDesc,LVT_rc%input_lnc*LVT_rc%input_lnr,&
               LVT_rc%rlat_dn,LVT_rc%rlon_dn, &
               LVT_rc%n11_dn, LVT_rc%n12_dn, &
               LVT_rc%n21_dn, LVT_rc%n22_dn, &
               LVT_rc%w11_dn, LVT_rc%w12_dn, &
               LVT_rc%w21_dn, LVT_rc%w22_dn)
          
          write(LVT_logunit,*) &
               '[ERR] currently interpolation option for analysis grid is not supported'
          call LVT_endrun()
          
       else !upscale/average
          
          allocate(LVT_rc%n11_up(LVT_rc%input_lnc*LVT_rc%input_lnr))
          
          call upscaleByAveraging_input(LVT_rc%input_gridDesc,&
               LVT_rc%gridDesc,LVT_rc%input_lnc*LVT_rc%input_lnr, &
               LVT_rc%lnc*LVT_rc%lnr, LVT_rc%n11_up)
          
          call upscaleByAveraging(LVT_rc%input_lnc*LVT_rc%input_lnr,&
               LVT_rc%lnc*LVT_rc%lnr,LVT_rc%udef,&
               LVT_rc%n11_up,li,mask_in,lo,mask_out)
          
       endif
    elseif(LVT_rc%spTransformAnlysDomain.eq."neighbor") then 
       

       allocate(LVT_rc%rlat_dn(LVT_rc%lnc*LVT_rc%lnr))
       allocate(LVT_rc%rlon_dn(LVT_rc%lnc*LVT_rc%lnr))
       allocate(LVT_rc%n11_dn(LVT_rc%lnc*LVT_rc%lnr))
       
       call neighbor_interp_input (&
            LVT_rc%input_gridDesc,&
            LVT_rc%gridDesc,&
            LVT_rc%lnc*LVT_rc%lnr,&
            LVT_rc%rlat_dn,LVT_rc%rlon_dn, &
            LVT_rc%n11_dn)

       call neighbor_interp(LVT_rc%gridDesc,&
            li,mask_in,lo,mask_out,&
            LVT_rc%input_lnc*LVT_rc%input_lnr,&
            LVT_rc%lnc*LVT_rc%lnr, & 
            LVT_rc%rlat_dn,LVT_rc%rlon_dn,&
            LVT_rc%n11_dn,LVT_rc%udef, ios)

    endif
       
    LVT_rc%ngrid=0
    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          if(mask_out(c+(r-1)*LVT_rc%lnc).gt.0.99 .and. & 
               mask_out(c+(r-1)*LVT_rc%lnc).lt.3.01) then
             LVT_rc%ngrid=LVT_rc%ngrid+1
          endif
       enddo
    enddo
    write(LVT_logunit,*) '[INFO] Total number of grid points in the LVT domain',&
         LVT_rc%ngrid     
    
    count1 = 1
    allocate(LVT_domain%grid(LVT_rc%ngrid))
    allocate(LVT_domain%gindex(LVT_rc%lnc, LVT_rc%lnr))
    
    do r=1,LVT_rc%lnr
       do c=1,LVT_rc%lnc
          LVT_domain%gindex(c,r) = -1

          ! EMK...Bug fix to gridDesc indices
!          locallat = LVT_rc%gridDesc(4)+(r-1)*LVT_rc%gridDesc(9)
!          locallon = LVT_rc%gridDesc(5)+(c-1)*LVT_rc%gridDesc(10)
          locallat = LVT_rc%gridDesc(4)+(r-1)*LVT_rc%gridDesc(10)
          locallon = LVT_rc%gridDesc(5)+(c-1)*LVT_rc%gridDesc(9)
          
          if(mask_out(c+(r-1)*LVT_rc%lnc).gt.0.99 .and. & 
               mask_out(c+(r-1)*LVT_rc%lnc).lt.3.01) then
             LVT_domain%grid(count1)%lat = locallat
             LVT_domain%grid(count1)%lon = locallon
             LVT_domain%grid(count1)%col = c
             LVT_domain%grid(count1)%row = r
             LVT_domain%gindex(c,r) = count1
             count1 = count1+1
          endif
       enddo
    enddo
  
    call LVT_compute_domainIndices
  end subroutine create_LVT_gridspace




!BOP
! 
! !ROUTINE: LVT_compute_domainIndices
! \label{LVT_compute_domainIndices}
!
! !INTERFACE: 
  subroutine LVT_compute_domainIndices()
! 
! !USES:

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: i 
    integer :: k
    integer :: status
    
    LVT_ngrids(LVT_localPet) = LVT_rc%ngrid
    LVT_gdeltas(LVT_localPet) = LVT_rc%ngrid

!-----------------------------------------------------------------------------
!  The grid, tile space sizes of the decomposed domains are gathered 
!  from each processor to compute the total sizes of the entire domain. 
!-----------------------------------------------------------------------------
#if (defined SPMD)
    call MPI_ALLGATHER(LVT_ngrids(LVT_localPet),1,MPI_INTEGER,&
         LVT_ngrids(:),1,MPI_INTEGER,&
         MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(LVT_gdeltas(LVT_localPet),1,MPI_INTEGER,&
         LVT_gdeltas(:),1,MPI_INTEGER,&
         MPI_COMM_WORLD,ierr)
#endif

    LVT_rc%glbngrid = 0 
    do i=0,LVT_npes-1
       LVT_rc%glbngrid = LVT_rc%glbngrid + LVT_ngrids(i)
    enddo
    
    if(LVT_masterproc) then 
       LVT_goffsets(0) = 0 
       do i=1,LVT_npes-1
          LVT_goffsets(i) = LVT_goffsets(i-1)+LVT_gdeltas(i-1)
       enddo

    end if
    
#if (defined SPMD)
    call MPI_BCAST(LVT_goffsets(:), LVT_npes, MPI_INTEGER,0, &
         MPI_COMM_WORLD, ierr)
#endif
    
  end subroutine LVT_compute_domainIndices

!BOP
! 
! !ROUTINE: LVT_domain_setup
! \label{LVT_domain_setup}
!
! !INTERFACE: 
  subroutine LVT_domain_setup(source)
! 
! !USES: 

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routines computes the variables used to map between the tile, grid, 
!  spaces in LVT. The global indices including the halo regions are also 
!  computed. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    implicit none
    integer              :: source

    integer, allocatable :: ntiles_pergrid(:)
    integer, allocatable :: ntiles_pergrid_red(:)
    integer, allocatable :: npatch_pergrid(:,:)
    integer, allocatable :: npatch_pergrid_red(:,:)
    integer, allocatable :: gtmp(:,:)
    integer, allocatable :: gtmp1(:)
    integer              :: count, index
    integer              :: t, gid, gid_red, ierr
    integer              :: c,r, l,m, c1,c2, r1,r2
    integer              :: ntiles, npatch,stid

!gathering the number of tiles per grid to the master processor, in the
!right order
    allocate(LVT_LIS_domain(source)%ntiles_pergrid(&
         LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr))
    allocate(LVT_LIS_domain(source)%str_tind(&
         LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr))
    allocate(ntiles_pergrid(LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr))
    allocate(ntiles_pergrid_red(LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr))

    allocate(npatch_pergrid(&
         LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr,LVT_rc%max_model_types))
    allocate(npatch_pergrid_red(&
         LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr,LVT_rc%max_model_types))
    do m=1,LVT_rc%max_model_types
       allocate(LVT_surface(source, m)%npatch_pergrid(&
            LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr))
       allocate(LVT_surface(source,m)%str_patch_ind(&
            LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr))
    enddo

    do t=1,LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr
       ntiles_pergrid(t) = 0 
    enddo
    do t=1,LVT_LIS_rc(source)%ntiles
       index = LVT_LIS_domain(source)%gindex(&
            LVT_LIS_domain(source)%tile(t)%col,&
            LVT_LIS_domain(source)%tile(t)%row)
       if(index.ne.-1) then 
          gid = LVT_LIS_domain(source)%tile(t)%col+&              
                (LVT_LIS_domain(source)%tile(t)%row-1)*LVT_LIS_rc(source)%lnc
          ntiles_pergrid(gid) = ntiles_pergrid(gid)+1
       endif
    enddo

    do r=LVT_lis_nss_halo_ind(source,LVT_localPet+1), LVT_lis_nse_halo_ind(source,LVT_localPet+1)
       do c=LVT_lis_ews_halo_ind(source,LVT_localPet+1), LVT_lis_ewe_halo_ind(source,LVT_localPet+1)
          c1 = c-LVT_lis_ews_halo_ind(source,LVT_localPet+1)+1
          r1 = r-LVT_lis_nss_halo_ind(source,LVT_localPet+1)+1
          gid = c1+(r1-1)*LVT_LIS_rc(source)%lnc
          
          if(r.ge.LVT_lis_nss_ind(source,LVT_localPet+1).and. &
               r.le.LVT_lis_nse_ind(source,LVT_localPet+1).and.&
               c.ge.LVT_lis_ews_ind(source,LVT_localPet+1).and.&
               c.le.LVT_lis_ewe_ind(source,LVT_localPet+1))then !points not in halo
             c2 = c-LVT_lis_ews_ind(source,LVT_localPet+1)+1
             r2 = r-LVT_lis_nss_ind(source,LVT_localPet+1)+1
             gid_red = c2+(r2-1)*LVT_LIS_rc(source)%lnc_red
             ntiles_pergrid_red(gid_red) = ntiles_pergrid(gid)
          endif
       enddo
    enddo

    deallocate(ntiles_pergrid)


    if(LVT_masterproc) then 
       allocate(gtmp(LVT_LIS_rc(source)%gnc,LVT_LIS_rc(source)%gnr))
       allocate(gtmp1(LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr))
    else
       allocate(gtmp1(1))
    endif

#if (defined SPMD)
    call MPI_GATHERV(ntiles_pergrid_red,LVT_lis_deltas(source,LVT_localPet),MPI_INTEGER,gtmp1,&
         LVT_lis_deltas(source,:),LVT_lis_offsets(source,:),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
    gtmp1 = ntiles_pergrid_red
#endif
    deallocate(ntiles_pergrid_red)

    if(LVT_masterproc) then 
       count=1
       do l=1,LVT_npes
          do r=LVT_lis_nss_ind(source,l), LVT_lis_nse_ind(source,l)
             do c=LVT_lis_ews_ind(source,l), LVT_lis_ewe_ind(source,l)
                gtmp(c,r) = gtmp1(count)
                count = count+1
             enddo
          enddo
       enddo

       count=1
       do r=1,LVT_LIS_rc(source)%gnr
          do c=1,LVT_LIS_rc(source)%gnc
             LVT_LIS_domain(source)%ntiles_pergrid(count) = gtmp(c,r)
             count = count+1
          enddo
       enddo
       
       deallocate(gtmp)
       deallocate(gtmp1)
       
      
       if(LVT_LIS_domain(source)%ntiles_pergrid(1).ge.0) then 
          LVT_LIS_domain(source)%str_tind(1) = 1
       else
          LVT_LIS_domain(Source)%str_tind(1) = 0 
       endif
       
       do t=2,LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr
          LVT_LIS_domain(source)%str_tind(t) = LVT_LIS_domain(source)%str_tind(t-1)+&
               LVT_LIS_domain(source)%ntiles_pergrid(t-1)
       enddo       
    else
       deallocate(gtmp1)
    endif
#if (defined SPMD)
    call MPI_Bcast(LVT_LIS_domain(source)%ntiles_pergrid,LVT_rc%lis_gnc*LVT_rc%lis_gnr,&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(LVT_LIS_domain(source)%str_tind,LVT_rc%lis_gnc*LVT_rc%lis_gnr,&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif   

    if(LVT_masterproc) then 
       LVT_LIS_rc(source)%glbngrid_red=0
       LVT_LIS_rc(source)%glbntiles_red=0
       do l=1,LVT_npes
          do r=LVT_lis_nss_halo_ind(source,l),LVT_lis_nse_halo_ind(source,l)
             do c=LVT_lis_ews_halo_ind(source,l),LVT_lis_ewe_halo_ind(source,l)
                if(r.ge.LVT_lis_nss_ind(source,l).and.&
                     r.le.LVT_lis_nse_ind(source,l).and.&
                     c.ge.LVT_lis_ews_ind(source,l).and.&
                     c.le.LVT_lis_ewe_ind(source,l))then !points not in halo
                   gid = c+(r-1)*LVT_LIS_rc(source)%gnc
                   ntiles = LVT_LIS_domain(source)%ntiles_pergrid(gid)
                   stid = LVT_LIS_domain(source)%str_tind(gid)
                   if ( ntiles .ne. 0 ) then
                      LVT_LIS_rc(source)%glbngrid_red = LVT_LIS_rc(source)%glbngrid_red + 1
                   endif
                   do t=1,ntiles
                      LVT_LIS_rc(source)%glbntiles_red = LVT_LIS_rc(source)%glbntiles_red + 1
                   enddo
                endif
             enddo
          enddo
       enddo
    endif    
#if (defined SPMD)
    call MPI_Bcast(LVT_LIS_rc(source)%glbntiles_red,1,&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    call MPI_Bcast(LVT_LIS_rc(source)%glbngrid_red,1,&
         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

!new stuff
    do t=1,LVT_LIS_rc(source)%lnc*LVT_LIS_rc(source)%lnr
       npatch_pergrid(t,:) = 0
    enddo

    do m=1,LVT_rc%max_model_types
       do t=1,LVT_LIS_rc(source)%npatch(m)
          index =  LVT_LIS_domain(source)%gindex(LVT_surface(source,m)%lis_tile(t)%col,&
               LVT_surface(source,m)%lis_tile(t)%row)
          if(index.ne.-1) then 
             gid = LVT_surface(source,m)%lis_tile(t)%col + & 
                  (LVT_surface(source, m)%lis_tile(t)%row-1)*LVT_LIS_rc(source)%lnc
             npatch_pergrid(gid,m) = npatch_pergrid(gid,m) + 1
          endif
       enddo
    enddo

    do r=LVT_lis_nss_halo_ind(source,LVT_localPet+1), LVT_lis_nse_halo_ind(source,LVT_localPet+1)
       do c=LVT_lis_ews_halo_ind(source,LVT_localPet+1), LVT_lis_ewe_halo_ind(source,LVT_localPet+1)
          c1 = c-LVT_lis_ews_halo_ind(source,LVT_localPet+1)+1
          r1 = r-LVT_lis_nss_halo_ind(source,LVT_localPet+1)+1
          gid = c1+(r1-1)*LVT_LIS_rc(source)%lnc
          
          if(r.ge.LVT_lis_nss_ind(source,LVT_localPet+1).and. &
               r.le.LVT_lis_nse_ind(source,LVT_localPet+1).and.&
               c.ge.LVT_lis_ews_ind(source,LVT_localPet+1).and.&
               c.le.LVT_lis_ewe_ind(source,LVT_localPet+1))then !points not in halo
             c2 = c-LVT_lis_ews_ind(source,LVT_localPet+1)+1
             r2 = r-LVT_lis_nss_ind(source,LVT_localPet+1)+1
             gid_red = c2+(r2-1)*LVT_LIS_rc(source)%lnc_red
             npatch_pergrid_red(gid_red,:) = npatch_pergrid(gid,:)
          endif
       enddo
    enddo

    do m=1,LVT_rc%max_model_types
       if(LVT_masterproc) then 
          allocate(gtmp(LVT_LIS_rc(source)%gnc,LVT_LIS_rc(source)%gnr))
          allocate(gtmp1(LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr))
       else
          allocate(gtmp1(1))
       endif
       
#if (defined SPMD)
       call MPI_GATHERV(npatch_pergrid_red(:,m),&
            LVT_lis_deltas(source,LVT_localPet),MPI_INTEGER,gtmp1,&
            LVT_lis_deltas(source,n,:),LVT_lis_offsets(source,n,:),&
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#else
       gtmp1 = npatch_pergrid_red(:,m)
#endif

       if(LVT_masterproc) then 
          count=1
          do l=1,LVT_npes
             do r=LVT_lis_nss_ind(source,l), LVT_lis_nse_ind(source,l)
                do c=LVT_lis_ews_ind(source,l), LVT_lis_ewe_ind(source,l)
                   gtmp(c,r) = gtmp1(count)
                   count = count+1
                enddo
             enddo
          enddo
          
          count=1
          do r=1,LVT_LIS_rc(source)%gnr
             do c=1,LVT_LIS_rc(source)%gnc
                LVT_surface(source,m)%npatch_pergrid(count) = gtmp(c,r)
                count = count+1
             enddo
          enddo
          
          deallocate(gtmp)
          deallocate(gtmp1)
       
          
          if(LVT_surface(source,m)%npatch_pergrid(1).ge.0) then 
             LVT_surface(source,m)%str_patch_ind(1) = 1
          else
             LVT_surface(source,m)%str_patch_ind(1) = 0 
          endif
          
          do t=2,LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr
             LVT_surface(source,m)%str_patch_ind(t) = &
                  LVT_surface(source,m)%str_patch_ind(t-1)+&
                  LVT_surface(source,m)%npatch_pergrid(t-1)
          enddo
       else
          deallocate(gtmp1)
       endif
#if (defined SPMD)
       call MPI_Bcast(LVT_surface(source,m)%npatch_pergrid,&
            LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr,& 
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(LVT_surface(source,m)%str_patch_ind,&
            LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr,&
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif   

       if(LVT_masterproc) then 
          LVT_LIS_rc(source)%glbnpatch_red(m)=0
          do l=1,LVT_npes
             do r=LVT_lis_nss_halo_ind(source,l),LVT_lis_nse_halo_ind(source,l)
                do c=LVT_lis_ews_halo_ind(source,l),LVT_lis_ewe_halo_ind(source,l)
                   if(r.ge.LVT_lis_nss_ind(source,l).and.&
                        r.le.LVT_lis_nse_ind(source,l).and.&
                        c.ge.LVT_lis_ews_ind(source,l).and.&
                        c.le.LVT_lis_ewe_ind(source,l))then !points not in halo
                      gid = c+(r-1)*LVT_LIS_rc(source)%gnc
                      npatch = LVT_surface(source,m)%npatch_pergrid(gid)
                      stid = LVT_surface(source,m)%str_patch_ind(gid)
                      do t=1,npatch
                         LVT_LIS_rc(source)%glbnpatch_red(m) = &
                              LVT_LIS_rc(source)%glbnpatch_red(m) + 1
                      enddo
                   endif
                enddo
             enddo
          enddo
       endif
#if (defined SPMD)
       call MPI_Bcast(LVT_LIS_rc(source)%glbnpatch_red(m),1,&
            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

    enddo

  end subroutine LVT_domain_setup

  subroutine create_LISoutput_domain()

    use LVT_timeMgrMod

    integer              :: source
    integer              :: k
    integer              :: i,j
    integer              :: rc
    integer              :: ios
    integer              :: ftn
    integer              :: ncId, nrId
    logical              :: file_exists
    real                 :: lis_run_dd(8)
    character*10         :: time
!    logical              :: vic_flag(2)
!    character*100        :: vic_d1file(2)
!    character*100        :: vic_d2file(2)
!    character*100        :: vic_d3file(2)
!    real, allocatable    :: vic_depth(:,:)
!    integer              :: c,r,lis_gid

    if(LVT_rc%obs_duplicate) then 
       source =2 
    else
       source = 1
    endif

    call ESMF_ConfigFindLabel(LVT_config,"LIS output number of surface model types:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%nsf_model_types,rc=rc)
       call LVT_verify(rc,'LIS output number of surface model types: option not specified')
    enddo

    do k=1,source 
       allocate(LVT_LIS_rc(k)%sf_model_type_name_select(LVT_LIS_rc(k)%nsf_model_types))
       allocate(LVT_LIS_rc(k)%sf_model_type_select(LVT_LIS_rc(k)%nsf_model_types))
       
       allocate(LVT_LIS_rc(k)%sf_model_type_name(LVT_rc%max_model_types))
       allocate(LVT_LIS_rc(k)%sf_model_type(LVT_rc%max_model_types))

       allocate(LVT_LIS_rc(k)%npatch(LVT_rc%max_model_types))

       LVT_LIS_rc(k)%sf_model_type_name(1) = "LSM"
       LVT_LIS_rc(k)%sf_model_type(1)      = LVT_rc%lsm_index
       LVT_LIS_rc(k)%sf_model_type_name(2) = "Lake"
       LVT_LIS_rc(k)%sf_model_type(2)      = LVT_rc%lake_index
       LVT_LIS_rc(k)%sf_model_type_name(3) = "Glacier"
       LVT_LIS_rc(k)%sf_model_type(3)      = LVT_rc%glacier_index
       LVT_LIS_rc(k)%sf_model_type_name(4) = "Wetland"
       LVT_LIS_rc(k)%sf_model_type(4)      = LVT_rc%wetland_index
       LVT_LIS_rc(k)%sf_model_type_name(5) = "Openwater"
       LVT_LIS_rc(k)%sf_model_type(5)      = LVT_rc%openwater_index

    enddo
    
    call ESMF_ConfigFindLabel(LVT_config,"LIS output interval:",rc=rc)

    do k=1,source       
       call ESMF_ConfigGetAttribute(LVT_config,time,rc=rc)
       call LVT_verify(rc,"LIS output interval: not defined")
       call LVT_parseTimeString(time,LVT_LIS_rc(k)%ts)

       call LVT_update_timestep(LVT_rc, LVT_LIS_rc(k)%ts)
    enddo

    call ESMF_ConfigFindLabel(LVT_config,"LIS output surface model types:",rc=rc)

    do k=1,source
       do i=1,LVT_LIS_rc(k)%nsf_model_types
          call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%sf_model_type_name_select(i),rc=rc)
          call LVT_verify(rc,"LIS output surface model types: not defined")
          
          call LVT_mapSurfaceModelType(LVT_LIS_rc(k)%sf_model_type_name_select(i), &
               LVT_LIS_rc(k)%sf_model_type_select(i))
       enddo
    enddo

    call ESMF_ConfigFindLabel(LVT_config,"LIS output analysis data class:",rc=rc)

    do k=1,source       
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%anlys_data_class,rc=rc)
       call LVT_verify(rc,"LIS output analysis data class: not defined")
    enddo

    call ESMF_ConfigFindLabel(LVT_config,label="LIS output number of ensembles per tile:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%nensem,rc=rc)
       call LVT_verify(rc,'LIS output number of ensembles per tile: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config,label="LIS output model name:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%model_name,rc=rc)
       call LVT_verify(rc,'LIS output model name: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output maximum number of surface type tiles per grid:",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%surface_maxt, rc=rc)
       call LVT_verify(rc,'LIS output maximum number of surface type tiles per grid: not defined')
    enddo
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output minimum cutoff percentage (surface type tiles):",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%surface_minp, rc=rc)
       call LVT_verify(rc,'LIS output minimum cutoff percentage (surface type tiles): not defined')
    enddo
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output maximum number of soil texture tiles per grid:",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%soilt_maxt, rc=rc)
       call LVT_verify(rc,'LIS output maximum number of soil texture tiles per grid: not defined')
    enddo
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output minimum cutoff percentage (soil texture tiles):",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%soilt_minp, rc=rc)
       call LVT_verify(rc,'LIS output minimum cutoff percentage (soil texture tiles): not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output maximum number of soil fraction tiles per grid:",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%soilf_maxt, rc=rc)
       call LVT_verify(rc,'LIS output maximum number of soil fraction tiles per grid: not defined')
    enddo
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output minimum cutoff percentage (soil fraction tiles):",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%soilf_minp, rc=rc)
       call LVT_verify(rc,'LIS output minimum cutoff percentage (soil fraction tiles): not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output maximum number of elevation bands per grid:",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%elev_maxt, rc=rc)
       call LVT_verify(rc,'LIS output maximum number of elevation bands per grid: not defined')
    enddo
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output minimum cutoff percentage (elevation bands):",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%elev_minp, rc=rc)
       call LVT_verify(rc,'LIS output minimum cutoff percentage (elevation bands): not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output maximum number of slope bands per grid:",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%slope_maxt, rc=rc)
       call LVT_verify(rc,'LIS output maximum number of slope bands per grid: not defined')
    enddo
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output minimum cutoff percentage (slope bands):",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%slope_minp, rc=rc)
       call LVT_verify(rc,'LIS output minimum cutoff percentage (slope bands): not defined')
    enddo
    
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output maximum number of aspect bands per grid:",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%aspect_maxt, rc=rc)
       call LVT_verify(rc,'LIS output maximum number of aspect bands per grid: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output minimum cutoff percentage (aspect bands):",rc=rc)   
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%aspect_minp, rc=rc)
       call LVT_verify(rc,'LIS output minimum cutoff percentage (aspect bands): not defined')
    enddo

    
    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output domain and parameter file:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%domfile,rc=rc)
       call LVT_verify(rc,'LIS output domain and parameter file: not defined')
    enddo

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
    
    do k=1,source
       inquire(file=LVT_LIS_rc(k)%domfile, exist=file_exists)
       if(file_exists) then 
       
          ios = nf90_open(path=LVT_LIS_rc(k)%domfile,&
               mode=NF90_NOWRITE,ncid=ftn)
          call LVT_verify(ios,'Error in nf90_open in LISoutputInit')
          
          ios = nf90_get_att(ftn, NF90_GLOBAL, 'MAP_PROJECTION',&
               LVT_LIS_rc(k)%map_proj)
          call LVT_verify(ios, 'Error in nf90_get_att: MAP_PROJECTION')

          if(LVT_LIS_rc(k)%map_proj.eq."EQUIDISTANT CYLINDRICAL") then 

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
                  lis_run_dd(1))
             call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
                  lis_run_dd(2))
             call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',lis_run_dd(5))
             call LVT_verify(ios, 'Error in nf90_get_att: DX')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',lis_run_dd(6))
             call LVT_verify(ios, 'Error in nf90_get_att: DY')
             
             ios = nf90_inq_dimid(ftn,"east_west",ncId)
             call LVT_verify(ios,'Error in nf90_inq_dimid ')
             
             ios = nf90_inq_dimid(ftn,"north_south",nrId)
             call LVT_verify(ios,'Error in nf90_inq_dimid ')
             
             ios = nf90_inquire_dimension(ftn,ncId, len=LVT_LIS_rc(k)%gnc)
             call LVT_verify(ios,'Error in nf90_inquire_dimension')
             
             ios = nf90_inquire_dimension(ftn,nrId, len=LVT_LIS_rc(k)%gnr)
             call LVT_verify(ios,'Error in nf90_inquire_dimension ')
             
             LVT_LIS_rc(k)%lnc = LVT_LIS_rc(k)%gnc
             LVT_LIS_rc(k)%lnr = LVT_LIS_rc(k)%gnr

             lis_run_dd(4) = lis_run_dd(2) + (LVT_LIS_rc(k)%gnc-1)*lis_run_dd(5)
             lis_run_dd(3) = lis_run_dd(1) + (LVT_LIS_rc(k)%gnr-1)*lis_run_dd(6)
             LVT_LIS_rc(k)%gridDesc(:) = 0 
             LVT_LIS_rc(k)%gridDesc(1) = 0
             LVT_LIS_rc(k)%gridDesc(2) = LVT_LIS_rc(k)%gnc
             LVT_LIS_rc(k)%gridDesc(3) = LVT_LIS_rc(k)%gnr
             LVT_LIS_rc(k)%gridDesc(4) = lis_run_dd(1)
             LVT_LIS_rc(k)%gridDesc(5) = lis_run_dd(2)
             LVT_LIS_rc(k)%gridDesc(6) = 128
             LVT_LIS_rc(k)%gridDesc(7) = lis_run_dd(3)
             LVT_LIS_rc(k)%gridDesc(8) = lis_run_dd(4)
             LVT_LIS_rc(k)%gridDesc(9) = lis_run_dd(5)
             LVT_LIS_rc(k)%gridDesc(10) = lis_run_dd(6)
             LVT_LIS_rc(k)%gridDesc(11) = 64
             LVT_LIS_rc(k)%gridDesc(20) = 64
             
             if(LVT_LIS_rc(k)%gridDesc(10).gt.LVT_rc%gridDesc(10)) then !interpolate. 
                allocate(LVT_LIS_rc(k)%rlat_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%rlon_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%n11_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%n12_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%n21_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%n22_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%w11_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%w12_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%w21_dn(LVT_rc%lnc*LVT_rc%lnr))
                allocate(LVT_LIS_rc(k)%w22_dn(LVT_rc%lnc*LVT_rc%lnr))
                
                call bilinear_interp_input(LVT_LIS_rc(k)%gridDesc(:),&
                     LVT_rc%gridDesc,LVT_rc%lnc*LVT_rc%lnr,&
                     LVT_LIS_rc(k)%rlat_dn,LVT_LIS_rc(k)%rlon_dn, &
                     LVT_LIS_rc(k)%n11_dn, LVT_LIS_rc(k)%n12_dn, &
                     LVT_LIS_rc(k)%n21_dn, LVT_LIS_rc(k)%n22_dn, &
                     LVT_LIS_rc(k)%w11_dn, LVT_LIS_rc(k)%w12_dn, &
                     LVT_LIS_rc(k)%w21_dn, LVT_LIS_rc(k)%w22_dn)
             else
                allocate(LVT_LIS_rc(k)%n11_up(&
                     LVT_LIS_rc(k)%lnc*LVT_LIS_rc(k)%lnr))
                
                call upscaleByAveraging_input(LVT_LIS_rc(k)%gridDesc(:),&
                     LVT_rc%gridDesc,LVT_LIS_rc(k)%lnc*LVT_LIS_rc(k)%lnr,&
                     LVT_rc%lnc*LVT_rc%lnr, LVT_LIS_rc(k)%n11_up(:))
             endif
             
          elseif(LVT_LIS_rc(k)%map_proj.eq."LAMBERT CONFORMAL") then 
             
             LVT_LIS_rc(k)%lnc = LVT_LIS_rc(k)%gnc
             LVT_LIS_rc(k)%lnr = LVT_LIS_rc(k)%gnr

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LAT',&
                  lis_run_dd(1))
             call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LAT')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'SOUTH_WEST_CORNER_LON',&
                  lis_run_dd(2))
             call LVT_verify(ios, 'Error in nf90_get_att: SOUTH_WEST_CORNER_LON')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT2',&
                  lis_run_dd(3))
             call LVT_verify(ios, 'Error in nf90_get_att: TRUELAT2')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'TRUELAT1',&
                  lis_run_dd(4))
             call LVT_verify(ios, 'Error in nf90_get_att: TRUELAT1')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'DX',lis_run_dd(5))
             call LVT_verify(ios, 'Error in nf90_get_att: DX')
             
             ios = nf90_get_att(ftn, NF90_GLOBAL, 'DY',lis_run_dd(6))
             call LVT_verify(ios, 'Error in nf90_get_att: DY')

             ios = nf90_get_att(ftn, NF90_GLOBAL, 'STANDARD_LON',lis_run_dd(7))
             call LVT_verify(ios, 'Error in nf90_get_att: STANDARD_LON')

             ios = nf90_inq_dimid(ftn,"east_west",ncId)
             call LVT_verify(ios,'Error in nf90_inq_dimid ')
             
             ios = nf90_inq_dimid(ftn,"north_south",nrId)
             call LVT_verify(ios,'Error in nf90_inq_dimid ')
             
             ios = nf90_inquire_dimension(ftn,ncId, len=LVT_LIS_rc(k)%gnc)
             call LVT_verify(ios,'Error in nf90_inquire_dimension')
             
             ios = nf90_inquire_dimension(ftn,nrId, len=LVT_LIS_rc(k)%gnr)
             call LVT_verify(ios,'Error in nf90_inquire_dimension ')
             
             LVT_LIS_rc(k)%lnc = LVT_LIS_rc(k)%gnc
             LVT_LIS_rc(k)%lnr = LVT_LIS_rc(k)%gnr


             LVT_LIS_rc(k)%gridDesc(:) = 0 
             LVT_LIS_rc(k)%gridDesc(1) = 3
             LVT_LIS_rc(k)%gridDesc(2) = LVT_LIS_rc(k)%gnc
             LVT_LIS_rc(k)%gridDesc(3) = LVT_LIS_rc(k)%gnr
             LVT_LIS_rc(k)%gridDesc(4) = lis_run_dd(1)
             LVT_LIS_rc(k)%gridDesc(5) = lis_run_dd(2)
             LVT_LIS_rc(k)%gridDesc(6) = 8
             LVT_LIS_rc(k)%gridDesc(7) = lis_run_dd(3)
             LVT_LIS_rc(k)%gridDesc(8) = lis_run_dd(5)
             LVT_LIS_rc(k)%gridDesc(9) = lis_run_dd(6)
             LVT_LIS_rc(k)%gridDesc(10) = lis_run_dd(4)
             LVT_LIS_rc(k)%gridDesc(11) = lis_run_dd(7)

             if(LVT_LIS_rc(k)%gridDesc(8)/100.0.gt.LVT_rc%gridDesc(10)) then !interpolate. 
                print*, 'downscaling for LISlsmObs not yet supported '
                print*, 'stopping in LISlsm_obsMod'
                stop
             else
                allocate(LVT_LIS_rc(k)%n11_up(&
                     LVT_LIS_rc(k)%lnc*LVT_LIS_rc(k)%lnr))
                
                call upscaleByAveraging_input(LVT_LIS_rc(k)%gridDesc(:),&
                     LVT_rc%gridDesc,LVT_LIS_rc(k)%lnc*LVT_LIS_rc(k)%lnr,&
                     LVT_rc%lnc*LVT_rc%lnr, LVT_LIS_rc(k)%n11_up(:))
             endif

          else
             write(LVT_logunit,*) '[ERR] support for '//trim(LVT_LIS_rc(k)%map_proj)
             write(LVT_logunit,*) '[ERR] not currently implemented'
             call LVT_endrun()
          endif
          

          ios = nf90_close(ftn)
          call LVT_verify(ios,'Error in nf90_close')
       endif
    enddo
#endif

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output nest index:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%nest,rc=rc)
       call LVT_verify(rc,'LIS output nest index: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output directory:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%odir,rc=rc)
       call LVT_verify(rc,'LIS output directory: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output naming style:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%style,rc=rc)
       call LVT_verify(rc,'LIS output naming style: not defined')
    enddo


    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output methodology:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%wopt,rc=rc)
       call LVT_verify(rc,'LIS output methodology: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output format:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%format,rc=rc)
       call LVT_verify(rc,'LIS output format: not defined')
    enddo

    call ESMF_ConfigFindLabel(LVT_config, &
         label="LIS output attributes file:",rc=rc)
    do k=1,source
       call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%outputSpecFile,rc=rc)
       call LVT_verify(rc,'LIS output attributes file: not defined')
    end do

#if 0 
    vic_flag = .false. 
    do k=1,source
       if((LVT_LIS_rc(k)%model_name.eq."VIC").or.&
            (LVT_LIS_rc(k)%model_name.eq."VIC412L")) then 
          vic_flag(k) = .true. 
       endif
    enddo
    if(vic_flag(1).and.vic_flag(2)) then 
       call ESMF_ConfigFindLabel(LVT_config, &
            label="LIS output VIC soil depth1 file:",rc=rc)
       do k=1,source
          call ESMF_ConfigGetAttribute(LVT_config,vic_d1file(k),rc=rc)
          call LVT_verify(rc,'LIS output VIC soil depth1 file: not defined')
       enddo
       call ESMF_ConfigFindLabel(LVT_config, &
            label="LIS output VIC soil depth2 file:",rc=rc)
       do k=1,source
          call ESMF_ConfigGetAttribute(LVT_config,vic_d2file(k),rc=rc)
          call LVT_verify(rc,'LIS output VIC soil depth2 file: not defined')
       enddo

       call ESMF_ConfigFindLabel(LVT_config, &
            label="LIS output VIC soil depth3 file:",rc=rc)
       do k=1,source
          call ESMF_ConfigGetAttribute(LVT_config,vic_d3file(k),rc=rc)
          call LVT_verify(rc,'LIS output VIC soil depth3 file: not defined')
          
          allocate(LVT_LIS_rc(k)%vic_depth(3,LVT_rc%npts))
          allocate(vic_depth(LVT_LIS_rc(k)%lnc, LVT_LIS_rc(k)%lnr))

          ftn = LVT_getNextUnitNumber()
          open(ftn,file=vic_d1file(k),form='unformatted',access='direct',  &
               recl= LVT_LIS_rc(k)%lnc* LVT_LIS_rc(k)%lnr*4)
          read(ftn,rec=1) vic_depth
          call LVT_releaseUnitNumber(ftn)
          
          do r=1,LVT_LIS_rc(k)%lnr
             do c=1,LVT_LIS_rc(k)%lnc
                lis_gid = LVT_LIS_domain(k)%gindex(c,r)
                if(lis_gid.ne.-1) then 
                   LVT_LIS_rc(k)%vic_depth(1,lis_gid) = vic_depth(c,r)
                endif
             enddo
          enddo

          ftn = LVT_getNextUnitNumber()
          open(ftn,file=vic_d2file(k),form='unformatted',access='direct',  &
               recl= LVT_LIS_rc(k)%lnc* LVT_LIS_rc(k)%lnr*4)
          read(ftn,rec=1) vic_depth
          call LVT_releaseUnitNumber(ftn)

          do r=1,LVT_LIS_rc(k)%lnr
             do c=1,LVT_LIS_rc(k)%lnc
                lis_gid = LVT_LIS_domain(k)%gindex(c,r)
                if(lis_gid.ne.-1) then 
                   LVT_LIS_rc(k)%vic_depth(2,lis_gid) = vic_depth(c,r)
                endif
             enddo
          enddo

          ftn = LVT_getNextUnitNumber()
          open(ftn,file=vic_d3file(k),form='unformatted',access='direct',  &
               recl= LVT_LIS_rc(k)%lnc* LVT_LIS_rc(k)%lnr*4)
          read(ftn,rec=1) vic_depth
          call LVT_releaseUnitNumber(ftn)
          deallocate(vic_depth)

          do r=1,LVT_LIS_rc(k)%lnr
             do c=1,LVT_LIS_rc(k)%lnc
                lis_gid = LVT_LIS_domain(k)%gindex(c,r)
                if(lis_gid.ne.-1) then 
                   LVT_LIS_rc(k)%vic_depth(3,lis_gid) = vic_depth(c,r)
                endif
             enddo
          enddo
       enddo
    elseif(vic_flag(1)) then 
       call ESMF_ConfigGetAttribute(LVT_config,vic_d1file(1),&
            label="LIS output VIC soil depth1 file:",rc=rc)
       call LVT_verify(rc,'LIS output VIC soil depth1 file: not defined')

       call ESMF_ConfigGetAttribute(LVT_config,vic_d2file(1),&
            label="LIS output VIC soil depth2 file:",rc=rc)
       call LVT_verify(rc,'LIS output VIC soil depth2 file: not defined')

       call ESMF_ConfigGetAttribute(LVT_config,vic_d3file(1),&
            label="LIS output VIC soil depth3 file:",rc=rc)
       call LVT_verify(rc,'LIS output VIC soil depth3 file: not defined')
       
       allocate(LVT_LIS_rc(1)%vic_depth(3,LVT_rc%npts))
       allocate(vic_depth(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr))

       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d1file(1),form='unformatted',access='direct',  &
            recl= LVT_LIS_rc(1)%lnc* LVT_LIS_rc(1)%lnr*4)
       read(ftn,rec=1) vic_depth
       call LVT_releaseUnitNumber(ftn)

       do r=1,LVT_LIS_rc(1)%lnr
          do c=1,LVT_LIS_rc(1)%lnc
             lis_gid = LVT_LIS_domain(1)%gindex(c,r)
             if(lis_gid.ne.-1) then 
                LVT_LIS_rc(k)%vic_depth(1,lis_gid) = vic_depth(c,r)
             endif
          enddo
       enddo
       
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d2file(1),form='unformatted',access='direct',  &
            recl= LVT_LIS_rc(1)%lnc* LVT_LIS_rc(1)%lnr*4)
       read(ftn,rec=1) vic_depth
       call LVT_releaseUnitNumber(ftn)
       
       do r=1,LVT_LIS_rc(1)%lnr
          do c=1,LVT_LIS_rc(1)%lnc
             lis_gid = LVT_LIS_domain(1)%gindex(c,r)
             if(lis_gid.ne.-1) then 
                LVT_LIS_rc(k)%vic_depth(2,lis_gid) = vic_depth(c,r)
             endif
          enddo
       enddo
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d3file(1),form='unformatted',access='direct',  &
            recl= LVT_LIS_rc(1)%lnc* LVT_LIS_rc(1)%lnr*4)
       read(ftn,rec=1) vic_depth
       call LVT_releaseUnitNumber(ftn)

       do r=1,LVT_LIS_rc(1)%lnr
             do c=1,LVT_LIS_rc(1)%lnc
                lis_gid = LVT_LIS_domain(1)%gindex(c,r)
                if(lis_gid.ne.-1) then 
                   LVT_LIS_rc(k)%vic_depth(3,lis_gid) = vic_depth(c,r)
                endif
             enddo
          enddo
       deallocate(vic_depth)
    elseif(vic_flag(2)) then 
       call ESMF_ConfigGetAttribute(LVT_config,vic_d1file(2),&
            label="LIS output VIC soil depth1 file:",rc=rc)
       call LVT_verify(rc,'LIS output VIC soil depth1 file: not defined')

       call ESMF_ConfigGetAttribute(LVT_config,vic_d2file(2),&
            label="LIS output VIC soil depth2 file:",rc=rc)
       call LVT_verify(rc,'LIS output VIC soil depth2 file: not defined')

       call ESMF_ConfigGetAttribute(LVT_config,vic_d3file(2),&
            label="LIS output VIC soil depth3 file:",rc=rc)
       call LVT_verify(rc,'LIS output VIC soil depth3 file: not defined')
       
       allocate(LVT_LIS_rc(2)%vic_depth(3,LVT_rc%npts))
       allocate(vic_depth(LVT_LIS_rc(2)%lnc, LVT_LIS_rc(2)%lnr))
       
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d1file(2),form='unformatted',access='direct',  &
            recl= LVT_LIS_rc(2)%lnc* LVT_LIS_rc(2)%lnr*4)
       read(ftn,rec=1) vic_depth
       call LVT_releaseUnitNumber(ftn)
       
       do r=1,LVT_LIS_rc(2)%lnr
          do c=1,LVT_LIS_rc(2)%lnc
             lis_gid = LVT_LIS_domain(2)%gindex(c,r)
             if(lis_gid.ne.-1) then 
                LVT_LIS_rc(k)%vic_depth(1,lis_gid) = vic_depth(c,r)
             endif
          enddo
       enddo
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d2file(2),form='unformatted',access='direct',  &
            recl= LVT_LIS_rc(2)%lnc* LVT_LIS_rc(2)%lnr*4)
       read(ftn,rec=1)  vic_depth
       call LVT_releaseUnitNumber(ftn)
       
       do r=1,LVT_LIS_rc(2)%lnr
          do c=1,LVT_LIS_rc(2)%lnc
             lis_gid = LVT_LIS_domain(2)%gindex(c,r)
             if(lis_gid.ne.-1) then 
                LVT_LIS_rc(k)%vic_depth(2,lis_gid) = vic_depth(c,r)
             endif
          enddo
       enddo
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=vic_d3file(2),form='unformatted',access='direct',  &
            recl= LVT_LIS_rc(2)%lnc* LVT_LIS_rc(2)%lnr*4)
       read(ftn,rec=1) vic_depth
       call LVT_releaseUnitNumber(ftn)
       
       do r=1,LVT_LIS_rc(2)%lnr
          do c=1,LVT_LIS_rc(2)%lnc
             lis_gid = LVT_LIS_domain(2)%gindex(c,r)
             if(lis_gid.ne.-1) then 
                LVT_LIS_rc(k)%vic_depth(3,lis_gid) = vic_depth(c,r)
             endif
          enddo
       enddo
    else    
       call ESMF_ConfigFindLabel(LVT_config, &
            label="LIS output number of soil moisture layers:",rc=rc)
       do k=1,source
          call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%nsmlayers,rc=rc)
          call LVT_verify(rc,'LIS output number of soil moisture layers: not defined')
       enddo
       
       call ESMF_ConfigFindLabel(LVT_config, &
            label="LIS output number of soil temperature layers:",rc=rc)
       do k=1,source
          call ESMF_ConfigGetAttribute(LVT_config,LVT_LIS_rc(k)%nstlayers,rc=rc)
          call LVT_verify(rc,'LIS output number of soil temperature layers: not defined')
       enddo
       
       do k=1,source
          
          allocate(LVT_LIS_rc(k)%smthick(LVT_LIS_rc(k)%nsmlayers))
          allocate(LVT_LIS_rc(k)%smdepth(LVT_LIS_rc(k)%nsmlayers))
       
          allocate(LVT_LIS_rc(k)%stthick(LVT_LIS_rc(k)%nstlayers))
          allocate(LVT_LIS_rc(k)%stdepth(LVT_LIS_rc(k)%nstlayers))
          
       enddo
       call ESMF_ConfigFindLabel(LVT_config, 'LIS output soil moisture layer thickness:',rc=rc)
       call LVT_verify(rc,'LIS output soil moisture layer thickness: not defined')
       do k=1,source
          do j=1,LVT_LIS_rc(k)%nsmlayers
             call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%smthick(j),rc=rc)
          enddo
       enddo
       
       call ESMF_ConfigFindLabel(LVT_config, 'LIS output soil temperature layer thickness:',rc=rc)
       call LVT_verify(rc,'LIS output soil temperature layer thickness: not defined')
       do k=1,source
          do j=1,LVT_LIS_rc(k)%nstlayers
             call ESMF_ConfigGetAttribute(LVT_config, LVT_LIS_rc(k)%stthick(j),rc=rc)
          enddo
       enddo
       
       do k=1,source
          LVT_LIS_rc(k)%smdepth(1) = LVT_LIS_rc(k)%smthick(1)
          do j=2,LVT_LIS_rc(k)%nsmlayers
             LVT_LIS_rc(k)%smdepth(j) = LVT_LIS_rc(k)%smdepth(j-1)+LVT_LIS_rc(k)%smthick(j)
          enddo
          
          LVT_LIS_rc(k)%stdepth(1) = LVT_LIS_rc(k)%stthick(1)
          do j=2,LVT_LIS_rc(k)%nstlayers
             LVT_LIS_rc(k)%stdepth(j) = LVT_LIS_rc(k)%stdepth(j-1)+LVT_LIS_rc(k)%stthick(j)
          enddo
       enddo
    endif
#endif

    do k=1,source
       call LVT_quilt_lis_domain(k, LVT_LIS_rc(k)%gnc,LVT_LIS_rc(k)%gnr)
       call create_LIS_tilespace(k)
    end do
  end subroutine create_LISoutput_domain

  subroutine create_LIS_tilespace(k)

    implicit none

    integer      :: k 
    integer      :: rc
    integer      :: i
    integer      :: m
    integer, allocatable :: deblklist(:,:,:)
    integer :: stid, enid
    integer :: status
    type(ESMF_DistGrid) :: tileDG
    type(ESMF_DistGrid) :: gridDG
    type(ESMF_DistGrid) :: patchDG(LVT_rc%max_model_types)
    type(ESMF_DistGrid) :: gridEnsDG


!    if(LVT_LIS_rc(k)%wopt.eq."1d tilespace") then 
       
       call ESMF_ConfigFindLabel(LVT_config,&
            "LIS output elevation data source:",rc=rc)
       do i=1,k
          call ESMF_ConfigGetAttribute(LVT_config,&
               LVT_LIS_rc(k)%useelevationmap,rc=rc)
          call LVT_verify(rc,'LIS output elevation data source: not specified')
       enddo

       call ESMF_ConfigFindLabel(LVT_config,&
            "LIS output slope data source:",rc=rc)
       do i=1,k
          call ESMF_ConfigGetAttribute(LVT_config,&
               LVT_LIS_rc(k)%useslopemap,rc=rc)
          call LVT_verify(rc,'LIS output slope data source: not specified')
       enddo

       call ESMF_ConfigFindLabel(LVT_config,&
            "LIS output aspect data source:",rc=rc)
       do i=1,k
          call ESMF_ConfigGetAttribute(LVT_config,&
               LVT_LIS_rc(k)%useaspectmap,rc=rc)
          call LVT_verify(rc,'LIS output aspect data source: not specified')
       enddo

       call ESMF_ConfigFindLabel(LVT_config,&
            "LIS output soil texture data source:",rc=rc)
       do i=1,k
          call ESMF_ConfigGetAttribute(LVT_config,&
               LVT_LIS_rc(k)%usetexturemap,rc=rc)
          call LVT_verify(rc,'LIS output soil texture data source: not specified')
       enddo

       call ESMF_ConfigFindLabel(LVT_config,&
            "LIS output soil fraction data source:",rc=rc)
       do i=1,k
          call ESMF_ConfigGetAttribute(LVT_config,&
               LVT_LIS_rc(k)%usesoilfractionmap,rc=rc)
          call LVT_verify(rc,'LIS output soil fraction data source: not specified')
       enddo


       call LVT_LMLC_init(k,LVT_LIS_rc(k)%domfile)
       call LVT_topo_init(k,LVT_LIS_rc(k)%domfile)
       call LVT_soils_init(k,LVT_LIS_rc(k)%domfile)
       call make_LIS_domain(k)


       LVT_lis_ntiless(k,LVT_localPet) = LVT_LIS_rc(k)%ntiles
       LVT_lis_tdeltas(k,LVT_localPet) = LVT_LIS_rc(k)%ntiles
       
       do m=1,LVT_rc%max_model_types
          LVT_lis_npatches(k,m,LVT_localPet)     = LVT_LIS_rc(k)%npatch(m)
          LVT_lis_patch_deltas(k,m,LVT_localPet) = LVT_LIS_rc(k)%npatch(m)
       enddo
       
       LVT_lis_ngrids(k,LVT_localPet) = LVT_LIS_rc(k)%ngrid
       LVT_lis_gdeltas(k,LVT_localPet) = LVT_LIS_rc(k)%ngrid
       
       !-----------------------------------------------------------------------------
       !  The grid, tile space sizes of the decomposed domains are gathered 
       !  from each processor to compute the total sizes of the entire domain. 
       !-----------------------------------------------------------------------------
#if (defined SPMD)
       call MPI_ALLGATHER(LVT_lis_ntiless(k,LVT_localPet),1,MPI_INTEGER,&
            LVT_lis_ntiless(k,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(LVT_lis_ngrids(k,LVT_localPet),1,MPI_INTEGER,&
            LVT_lis_ngrids(k,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(LVT_lis_tdeltas(k,LVT_localPet),1,MPI_INTEGER,&
            LVT_lis_tdeltas(k,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       call MPI_ALLGATHER(LVT_lis_gdeltas(k,LVT_localPet),1,MPI_INTEGER,&
            LVT_lis_gdeltas(k,:),1,MPI_INTEGER,&
            MPI_COMM_WORLD,ierr)
       
       do m=1,LVT_rc%max_model_types
          call MPI_ALLGATHER(LVT_lis_npatches(k,m,LVT_localPet),1,MPI_INTEGER,&
               LVT_npatches(k,m,:),1,MPI_INTEGER, &
               MPI_COMM_WORLD,ierr)
       enddo
       do m=1,LVT_rc%max_model_types
          call MPI_ALLGATHER(LVT_lis_patch_deltas(k,m,LVT_localPet),1,MPI_INTEGER,&
               LVT_lis_patch_deltas(k,m,:),1,MPI_INTEGER, &
               MPI_COMM_WORLD,ierr)
       enddo
#endif
       
       LVT_LIS_rc(k)%glbntiles = 0 
       do i=0,LVT_npes-1
          LVT_LIS_rc(k)%glbntiles = LVT_LIS_rc(k)%glbntiles + LVT_lis_ntiless(k,i)
       enddo
       
       LVT_LIS_rc(k)%glbngrid = 0 
       do i=0,LVT_npes-1
          LVT_LIS_rc(k)%glbngrid = LVT_LIS_rc(k)%glbngrid + LVT_lis_ngrids(k,i)
       enddo
       
       do m=1,LVT_rc%max_model_types
          LVT_LIS_rc(k)%glbnpatch(m) = 0 
          do i=0,LVT_npes-1
             LVT_LIS_rc(k)%glbnpatch(m) = LVT_LIS_rc(k)%glbnpatch(m)+ &
                  LVT_lis_npatches(k,m,i)
          enddo
       enddo
       
       
       if(LVT_masterproc) then 
          LVT_lis_toffsets(k,0) = 0 
          do i=1,LVT_npes-1
             LVT_lis_toffsets(k,i) = LVT_lis_toffsets(k,i-1)+LVT_lis_tdeltas(k,i-1)
          enddo
          LVT_lis_goffsets(k,0) = 0 
          do i=1,LVT_npes-1
             LVT_lis_goffsets(k,i) = LVT_lis_goffsets(k,i-1)+LVT_lis_gdeltas(k,i-1)
          enddo
          
          LVT_lis_patch_offsets(k,:,0) =0
          do i=1,LVT_npes-1
             do m=1,LVT_rc%max_model_types
                LVT_lis_patch_offsets(k,m,i) = LVT_lis_patch_offsets(k,m,i-1)+&
                     LVT_lis_patch_deltas(k,m,i-1)
             enddo
          enddo
          
       end if
       
#if (defined SPMD)
       call MPI_BCAST(LVT_lis_toffsets(k,:), LVT_npes, MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)
       call MPI_BCAST(LVT_lis_goffsets(k,:), LVT_npes, MPI_INTEGER,0, &
            MPI_COMM_WORLD, ierr)
       do m=1,LVT_rc%max_model_types
          call MPI_BCAST(LVT_lis_patch_offsets(k,m,:), LVT_npes, MPI_INTEGER,0, &
               MPI_COMM_WORLD, ierr)
       enddo
       
#endif
       
       allocate(deblklist(1,2,LVT_npes))      
       do i=0,LVT_npes-1
          stid = LVT_lis_toffsets(k,i)+1
          enid = stid + LVT_lis_ntiless(k,i)-1
          
          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)          
       enddo
       
       tileDG = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/LVT_LIS_rc(k)%glbntiles/),&
            deBlockList=deblklist,rc=status)
       call LVT_verify(status, 'ESMF_DistGridCreate failed')
       
       do i=0,LVT_npes-1
          stid = LVT_lis_goffsets(k,i)+1
          enid = stid + LVT_lis_ngrids(k,i)-1
          
          deblklist(:,1,i+1) = (/stid/)
          deblklist(:,2,i+1) = (/enid/)
       enddo
       
       gridDG = ESMF_DistGridCreate(minIndex=(/1/), maxIndex=(/LVT_LIS_rc(k)%glbngrid/),&
            deBlockList=deblklist,rc=status)
       call LVT_verify(status,'ESMF_DistGridCreate failed')
       
       do m=1,LVT_rc%max_model_types
          
          do i=0,LVT_npes-1
             stid = LVT_lis_patch_offsets(k,m,i)+1
             enid = stid + LVT_lis_npatches(k,m,i)-1
             
             deblklist(:,1,i+1) = (/stid/)
             deblklist(:,2,i+1) = (/enid/)
          enddo
          
          patchDG(m) = ESMF_DistGridCreate(minIndex=(/1/),&
               maxIndex=(/LVT_LIS_rc(k)%glbnpatch(m)/),&
               deBlockList=deblklist,rc=status)
          call LVT_verify(status,'ESMF_DistGridCreate failed')
       enddo
       
       deallocate(deblklist)
       
       allocate(deblklist(2,2,LVT_npes))
       do i=0,LVT_npes-1
          stid = LVT_lis_goffsets(k,i)+1
          enid = stid+LVT_lis_ngrids(k,i)-1
          deblklist(:,1,i+1) = (/stid,1/)
          deblklist(:,2,i+1) = (/enid,LVT_LIS_rc(k)%nensem/)
       enddo
       
       gridEnsDG = ESMF_DistGridCreate(minIndex=(/1,1/),&
            maxIndex=(/LVT_LIS_rc(k)%glbngrid,LVT_LIS_rc(k)%nensem/),&
            deBlockList=deblklist,rc=status)
       call LVT_verify(status,'ESMF_DistGridCreate failed')
    
       LVT_vecTile(k) = ESMF_GridCreate(name = "LVT Tile Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = tileDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LVT_verify(status,'ESMF_DistGridCreate failed')
       
       LVT_vecGrid(k) = ESMF_GridCreate(name = "LVT Grid Space",&
            coordTypeKind=ESMF_TYPEKIND_R4, distGrid = gridDG,&
            gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
       call LVT_verify(status,'ESMF_DistGridCreate failed')
       
       LVT_ensOnGrid(k) = ESMF_GridCreate(name = "LVT Grid Ensemble Space",&
            coordTypeKind = ESMF_TYPEKIND_R4, distGrid = gridEnsDG,&
            gridEdgeLWidth=(/0,0/), gridEdgeUWidth=(/0,0/),rc=status)
       call LVT_verify(status,'ESMF_DistGridCreate failed')
       
       do m=1,LVT_rc%max_model_types
          LVT_vecPatch(k,m) = ESMF_GridCreate(name="LVT Patch Space",&
               coordTypeKind=ESMF_TYPEKIND_R4, distGrid = patchDG(m),&
               gridEdgeLWidth=(/0/), gridEdgeUWidth=(/0/),rc=status)
          call LVT_verify(status,'ESMF_DistGridCreate failed')
       enddo
       
       deallocate(deblklist)
       
!    endif

  end subroutine create_LIS_tilespace

!BOP
! 
! !ROUTINE: make_LIS_domain
! \label{make_LIS_domain}
!
! !INTERFACE: 
  subroutine make_LIS_domain(source)
!
! !DESCRIPTION: 
!  This subroutine generates the tile, patch and grid data structures
!  based on the vegetation, soils, and topography datasets. The
!  routine first computes the dominant distributions for vegetation
!  and soils and then generates the multi-dimensional tile structure.  
!
!EOP
    implicit none
    integer :: source

    integer :: c,r
    logical :: soilt_selected
    logical :: soilf_selected
    logical :: elev_selected
    logical :: slope_selected
    logical :: aspect_selected
    
    if(LVT_LIS_rc(source)%usetexturemap.ne."none") then 
       soilt_selected = .true. 
    else
       soilt_selected = .false.
    endif
    
    if(LVT_LIS_rc(source)%usesoilfractionmap.ne."none") then 
       soilf_selected = .true. 
    else
       soilf_selected = .false.
    endif
    
    if(LVT_LIS_rc(source)%useelevationmap.ne."none") then 
       elev_selected = .true. 
    else
       elev_selected = .false. 
    endif
    
    if(LVT_LIS_rc(source)%useslopemap.ne."none") then 
       slope_selected = .true. 
    else
       slope_selected = .false. 
    endif
    
    if(LVT_LIS_rc(source)%useaspectmap.ne."none") then 
       aspect_selected = .true. 
    else
       aspect_selected = .false. 
    endif
    !-----------------------------------------------------------------------
    ! normalize the vegetation distributions
    !-----------------------------------------------------------------------
    ! EMK...Special handling of single surface tile case.  First find
    ! dominant surface model type in each grid point (LSM, Openwater, etc),
    ! then find the dominant tile for that model type.
    if (LVT_LIS_rc(source)%surface_maxt == 1) then          
       call calculate_dominant_sfctile(source, &
            LVT_LMLC(source)%dommask, &
            LVT_LIS_rc(source)%nsurfacetypes, &
            LVT_LIS_rc(source)%surface_minp, &
            LVT_LIS_rc(source)%surface_maxt, &
            LVT_LMLC(source)%surfacetype,&
            LVT_LMLC(source)%landcover)
    else
       call calculate_domdistribution(source,LVT_LIS_rc(source)%nsurfacetypes, &
            LVT_LIS_rc(source)%surface_minp, LVT_LIS_rc(source)%surface_maxt, &
            LVT_LMLC(source)%landcover)
    endif

    if(soilt_selected) then 
       call calculate_domdistribution(source,LVT_LIS_rc(source)%nsoiltypes,&
            LVT_LIS_rc(source)%soilt_minp, LVT_LIS_rc(source)%soilt_maxt,&
            LVT_soils(source)%texture)
    endif
    
    if(soilf_selected) then 
       call calculate_domdistribution(source,LVT_LIS_rc(source)%nsoilfbands,&
            LVT_LIS_rc(source)%soilf_minp, LVT_LIS_rc(source)%soilf_maxt,&
            LVT_soils(source)%soilffgrd)
    endif
    
    if(elev_selected) then 
       call calculate_domdistribution(source,LVT_LIS_rc(source)%nelevbands, &
            LVT_LIS_rc(source)%elev_minp, LVT_LIS_rc(source)%elev_maxt,&
            LVT_topo(source)%elevfgrd)
    endif
    
    if(slope_selected) then 
       call calculate_domdistribution(source,LVT_LIS_rc(source)%nslopebands, &
            LVT_LIS_rc(source)%slope_minp, LVT_LIS_rc(source)%slope_maxt, &
            LVT_topo(source)%slopefgrd)
    endif
    
    if(aspect_selected) then 
       call calculate_domdistribution(source,LVT_LIS_rc(source)%naspectbands, &
            LVT_LIS_rc(source)%aspect_minp, LVT_LIS_rc(source)%aspect_maxt, &
            LVT_topo(source)%aspectfgrd)
    endif
    
    call create_tilespace(source,soilt_selected, soilf_selected, &
         elev_selected, slope_selected, aspect_selected)

  end subroutine make_LIS_domain

  ! EMK...Special subroutine for finding single dominant surface tile in each
  ! grid point. This first considers the single dominant surface model type 
  ! (LSM, Openwater, etc), and then identifies the dominant tile for that 
  ! model type.
  !
  ! The multiple surface tile case is more complicated and is not addressed
  ! here.  
  subroutine calculate_dominant_sfctile(n, dommask, &
       ntypes, minp, maxt, surfacetypes, fgrd)

     ! Defaults
     implicit none
     
     ! Arguments
     integer, intent(in)  :: n
     real, intent(in)     :: dommask(LVT_LIS_rc(n)%gnc,&
          LVT_LIS_rc(n)%gnr)
     integer, intent(in)  :: ntypes
     real, intent(in)     :: minp
     integer, intent(in)  :: maxt
     real, intent(in)     :: surfacetypes(LVT_LIS_rc(n)%gnc, &
         LVT_LIS_rc(n)%gnr, &
         ntypes)
     real                 :: fgrd(LVT_LIS_rc(n)%gnc,&
         LVT_LIS_rc(n)%gnr, &
         ntypes)
    
     ! Locals
     integer :: c, r, t, m, mm, tt
     real    :: maxv
     real :: model_type_fractions(LVT_rc%max_model_types)

     ! Sanity check
     if (maxt > 1) then
        write(LVT_logunit,*) &
             '[ERR], calculate_dominant_sfctile only works for maxt = 1!'
        call LVT_endrun()
     end if
     
     ! Loop through each patch, figure out the dominant surface model type,
     ! and then find the dominant land cover category for that model type.
     do r=1,LVT_LIS_rc(n)%gnr
        do c=1,LVT_LIS_rc(n)%gnc

           ! Skip points that are outside of the domain--these have no
           ! land cover values defined.
           if (dommask(c,r) .le. 0) cycle

           ! Consolidate the tile fractions into fractions by surface model 
           ! type (LSM, Openwater, etc)
           model_type_fractions(:) = 0
           do t=1,ntypes
              m = nint(surfacetypes(c,r,t))
              if (m .eq. 0) cycle ! Empty tile
              model_type_fractions(m) = &
                   model_type_fractions(m) + fgrd(c,r,t)
           end do ! t
           
           ! Now find which surface model type has the highest fraction.
           maxv = 0
           mm = 0
           do m = 1, LVT_rc%max_model_types
              ! Skip model type if not used
              if (.not. isSurfaceTypeSelected(n, LVT_LIS_rc(n)%sf_model_type(m))) cycle
              if (model_type_fractions(m) > maxv) then
                 maxv = model_type_fractions(m)
                 mm = m
              end if
           end do ! m
           
           ! Sanity check
           if (mm .eq. 0) then
              write(LVT_logunit,*) &
                   '[ERR] No dominant surface model type found!'
              write(LVT_logunit,*) &
                   'c,r,fgrd: ',c,r,fgrd(c,r,:)
              call LVT_flush(LVT_logunit)
              call LVT_endrun()
           end if

           ! We now know which surface model type is dominant.
           ! Identify and exclude tiles not used by that model type.
           do t = 1,ntypes
              if (nint(surfacetypes(c,r,t)) .eq. mm) cycle
              fgrd(c,r,t) = 0.0
           end do ! t

           ! At this point, remaining tiles are associated with only a single 
           ! surface model type. Exclude tiles below minimum tile grid area). 
           do t=1,ntypes
              if (fgrd(c,r,t).lt.minp) then
                 fgrd(c,r,t)=0.0 
              endif
           enddo ! t

           ! Now find the tile with the largest fraction.
           maxv = 0
           tt = 0
           do t=1,ntypes
              if (fgrd(c,r,t) > maxv) then
                 maxv = fgrd(c,r,t)
                 tt = t
              endif
           enddo ! t
           
           ! Sanity check: Make sure we still have one tile left!
           if (tt .eq. 0) then
              write(LVT_logunit,*) &
                   '[ERR] No surface tiles remain!'
              write(LVT_logunit,*) &
                   'c,r,fgrd: ',c,r,fgrd(c,r,:)
              call LVT_flush(LVT_logunit)
              call LVT_endrun()
           end if

           ! Exclude the smaller tiles, and set the dominant tile to 100%
           do t = 1,ntypes
              if (t .eq. tt) then
                 fgrd(c,r,t) = 1.0
              else
                 fgrd(c,r,t) = 0.0
              end if
           end do ! t

        end do ! c
     end do ! r
     
  end subroutine calculate_dominant_sfctile

!BOP
!
! !ROUTINE: calculate_domdistribution
! \label{calculate_domdistribution}
!
! !REVISION HISTORY:
!  09 Sept 2004: Sujay Kumar ; Initial version
!
! !INTERFACE:
  subroutine calculate_domdistribution(source,ntypes, minp, maxt, fgrd)
! !USES:
    implicit none

! !ARGUMENTS: 
    integer              :: source
    integer              :: ntypes
    real                 :: minp
    integer              :: maxt
    real                 :: fgrd(LVT_LIS_rc(source)%gnc,&
         LVT_LIS_rc(source)%gnr, &
         ntypes)
    !
! !DESCRIPTION:
!  This routine determines the percentages of 
!  dominant types for subgrid tiling.
!  The routine uses the distribution of data within a 
!  grid cell and imposes the specified cutoff percentage and 
!  the maximum number of data types allowed to create the 
!  number of tiles. 
!
!  The arguments are: 
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ntypes]
!    number of categories or bins in the input data. 
!   \item[minp]
!    minimum cutoff percentage to be used in the renormalization
!   \item[maxt]
!    maximum number of tiles per grid cell to be used in the renormalization
!   \item[fgrd]
!    fraction of grid covered by a data type. 
!  \end{description}
!EOP

    real,    allocatable :: tsum(:,:)  !Temporary processing variable
    integer, allocatable :: pfrac(:,:,:)
    integer :: c, r, t, ierr, i, j
    real    :: rsum         
    real    :: fvt(ntypes)  
    real    :: maxv

    allocate(tsum(LVT_LIS_rc(source)%gnc, LVT_LIS_rc(source)%gnr), stat=ierr)
    call LVT_verify(ierr,'Error allocating tsum.')
    tsum = 0.0
!----------------------------------------------------------------------      
! Exclude tiles with (minimum tile grid area),  
! normalize remaining tiles to 100%
!----------------------------------------------------------------------      
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc    
          rsum=0.0
          do t=1,ntypes
             if(fgrd(c,r,t).lt.minp)then
                fgrd(c,r,t)=0.0 
             endif
             rsum=rsum+fgrd(c,r,t)
          enddo
!----------------------------------------------------------------------      
! renormalize veg fractions within a grid to 1
!----------------------------------------------------------------------      
          if(rsum.gt.0.0) then  
             do t=1,ntypes
                if(rsum.gt.0.0)fgrd(c,r,t)=fgrd(c,r,t)/rsum
             enddo
             
             rsum=0.0
             do t=1,ntypes
                rsum=rsum+fgrd(c,r,t)
             enddo
             
             if(rsum.lt.0.9999.or.rsum.gt.1.0001)then 
                write(LVT_logunit,*) '[ERR] The tile distribution do not sum to 100%'
                call LVT_endrun()
             endif
          endif
       enddo
    enddo

    allocate(pfrac(LVT_LIS_rc(source)%gnc,LVT_LIS_rc(source)%gnr,ntypes))
!----------------------------------------------------------------------      
! Exclude tiles with SURFACE_MAXT (Maximum Tiles per grid), 
!   normalize remaining tiles to 100%
! Determine the grid predominance order of the tiles
!  PFRAC(NT) will contain the predominance order of tiles
!----------------------------------------------------------------------      
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc 
          do t=1,ntypes
             fvt(t)=fgrd(c,r,t)
             pfrac(c,r,t)=0
          enddo
          do i=1,ntypes
             maxv=0.0
             t=0
             do j=1,ntypes
                if(fvt(j).gt.maxv)then
                   if(fgrd(c,r,j).gt.0) then
                      maxv=fvt(j)
                      t=j
                   endif
                endif
             enddo
             if(t.gt.0) then
                pfrac(c,r,t)=i
                fvt(t)=-999.0       
             endif
          enddo
       enddo
    enddo
!----------------------------------------------------------------------      
! Impose MAXT Cutoff
!----------------------------------------------------------------------
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc
          rsum=0.0
          do t=1,ntypes
             if(pfrac(c,r,t).lt.1) then
                fgrd(c,r,t)=0.0    
                pfrac(c,r,t)=0  
             endif
             if(pfrac(c,r,t).gt.maxt) then
                fgrd(c,r,t)=0.0            
                pfrac(c,r,t)=0  
             endif
             rsum=rsum+fgrd(c,r,t)
          enddo
!----------------------------------------------------------------------
! renormalize veg fractions within a grid to 1
!----------------------------------------------------------------------
          if(rsum.gt.0.0) then  
             do t=1,ntypes
                if(rsum.gt.0.0)fgrd(c,r,t)= fgrd(c,r,t)/rsum
             enddo
             
             rsum=0.0
             do t=1,ntypes
                rsum=rsum+ fgrd(c,r,t)  !recalculate rsum to check 
             enddo
             tsum(c,r)=rsum
           
             if(rsum.lt.0.9999.or.rsum.gt.1.0001)then  !check renormalization
                write(LVT_logunit,*) &
                     '[ERR] Renormalization failed in calculate_domdistribution'
                write(LVT_logunit,*) '[ERR] ',c,r,fgrd(c,r,:)
                call LVT_endrun()
             endif
          endif
       enddo
    enddo
    deallocate(pfrac)
    deallocate(tsum)

  end subroutine calculate_domdistribution


!BOP
!
! !ROUTINE: create_tilespace
!  \label{create_tilespace}
!
! !REVISION HISTORY:
!  09 Sept 2004: Sujay Kumar ; Initial version
!
! !INTERFACE:
  subroutine create_tilespace(source, soilt_selected, soilf_selected, &
       elev_selected,slope_selected, aspect_selected)

    implicit none
! !ARGUMENTS: 
    integer, intent(in)   :: source
    logical, intent(in)   :: soilt_selected
    logical, intent(in)   :: soilf_selected
    logical, intent(in)   :: elev_selected
    logical, intent(in)   :: slope_selected
    logical, intent(in)   :: aspect_selected

! !DESCRIPTION: 
!  This routine creates the tilespace based on the specified dimensions
!  of tiling. Vegetation is used by default. Soil texture, soil fractions, 
!  elevation, slope, aspect can be specified as additional dimensions for
!  tiling. The routine also generates the patch data structures which are
!  the tile data structures in each surface model class. 
! 
!  The arguments are: 
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[soilt\_selected]
!    flag to indicate if soil texture based tiling is used
!   \item[soilf\_selected]
!    flag to indicate if soil fraction based tiling is used
!   \item[soilf\_selected]
!    flag to indicate if soil fraction based tiling is used
!   \item[elev\_selected]
!    flag to indicate if elevation based tiling is used
!   \item[slope\_selected]
!    flag to indicate if slope based tiling is used
!   \item[aspect\_selected]
!    flag to indicate if aspect based tiling is used
!  \end{description}
!EOP
    real          :: locallat, locallon
    integer       :: c, r, t, m
    integer       :: iv, it, ie, is, ia
    integer       :: vegt
    integer       :: soilt
    real          :: sand, clay, silt
    real          :: elev, slope, aspect
    integer       :: sf_index
    integer       :: gnc, gnr
    integer       :: kk
    integer       :: kk_sf(LVT_rc%max_model_types)
    real, allocatable :: sumv(:)
    integer       :: ntiles_surface
    integer       :: npatch_surface
    integer       :: ntiles_soilt
    integer       :: ntiles_soilf
    integer       :: ntiles_soil
    integer       :: ntiles_elev
    integer       :: ntiles_slope
    integer       :: ntiles_aspect
    integer       :: gid
    real          :: temp
    integer       :: soilf_index
    integer       :: elev_index
    integer       :: slope_index
    integer       :: aspect_index
    real          :: mask_in(LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr)
    real          :: mask_out(LVT_rc%lnc*LVT_rc%lnr)
    logical*1     :: li(LVT_LIS_rc(source)%gnc*LVT_LIS_rc(source)%gnr)
    logical*1     :: lo(LVT_rc%lnc*LVT_rc%lnr)
    integer       :: count1

    gnc = LVT_LIS_rc(source)%gnc
    gnr = LVT_LIS_rc(source)%gnr
    
    LVT_LIS_rc(source)%ntiles=0
    LVT_LIS_rc(source)%npatch(:) = 0 

    do r=1,LVT_LIS_rc(source)%gnr     
       do c=1,LVT_LIS_rc(source)%gnc             
          do m=1,LVT_LIS_rc(source)%nensem
             if(LVT_LMLC(source)%dommask(c,r).gt.0) then
                
                ntiles_surface = compute_ntiles_surface(source,c,r)
                ntiles_soilt = compute_ntiles_soilt(source,c,r,soilt_selected)
                ntiles_soilf = compute_ntiles_soilf(source,c,r,soilf_selected)
                ntiles_elev = compute_ntiles_elev(source,c,r,elev_selected)
                ntiles_slope = compute_ntiles_slope(source,c,r,slope_selected)
                ntiles_aspect = compute_ntiles_aspect(source,c,r,aspect_selected)

                if(soilf_selected) then 
                   LVT_LIS_rc(source)%ntiles = &
                        LVT_LIS_rc(source)%ntiles+&
                        ntiles_surface* & 
                        ntiles_soilf* &
                        ntiles_elev*&
                        ntiles_slope*& 
                        ntiles_aspect
                else
                   LVT_LIS_rc(source)%ntiles = &
                        LVT_LIS_rc(source)%ntiles+&
                        ntiles_surface* & 
                        ntiles_soilt* &
                        ntiles_elev*&
                        ntiles_slope*& 
                        ntiles_aspect
                endif
                   
                do t=1,LVT_LIS_rc(source)%nsurfacetypes
                   sf_index = nint(LVT_LMLC(source)%surfacetype(c,r,t))
                   npatch_surface = compute_npatches_surface(source,c,r,t,sf_index)
                   if(sf_index.ne.0) then
                      if(soilf_selected) then 
                         LVT_LIS_rc(source)%npatch(sf_index) = & 
                              LVT_LIS_rc(source)%npatch(sf_index) + &
                              npatch_surface*& 
                              ntiles_soilf* & 
                              ntiles_elev*& 
                              ntiles_slope*& 
                              ntiles_aspect                         
                      else
                         LVT_LIS_rc(source)%npatch(sf_index) = & 
                              LVT_LIS_rc(source)%npatch(sf_index) + &
                              npatch_surface*& 
                              ntiles_soilt* & 
                              ntiles_elev*& 
                              ntiles_slope*& 
                              ntiles_aspect
                      endif
                   endif
                enddo
                
             endif
          enddo
       enddo
    enddo
    write(LVT_logunit,*) &
         '[INFO] Total number of LIS surface tiles',LVT_LIS_rc(source)%ntiles
    do m=1,LVT_rc%max_model_types
       if(isSurfaceTypeSelected(source, LVT_LIS_rc(source)%sf_model_type(m))) then 
          write(LVT_logunit,*) &
               '[INFO] Total number of LIS tiles for ',&
               trim(LVT_LIS_rc(source)%sf_model_type_name(m)),&
               LVT_LIS_rc(source)%npatch(m)
       endif
    enddo

    LVT_LIS_rc(source)%ngrid=0
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc
          if(LVT_LMLC(source)%dommask(c,r).gt.0.and.&
               abs(sum(LVT_LMLC(source)%landcover(c,r,:))-1.0).lt.0.001.and.&
               ifContainsSelectedSurfaceType(source,&
               nint(LVT_LMLC(source)%surfacetype(c,r,:)))) then 
             LVT_LIS_rc(source)%ngrid = LVT_LIS_rc(source)%ngrid+1             
          elseif(LVT_LMLC(source)%dommask(c,r).gt.0.and.&
               sum(LVT_LMLC(source)%landcover(c,r,:)).ne.0.and.&
               ifContainsSelectedSurfaceType(source,&
               nint(LVT_LMLC(source)%surfacetype(c,r,:)))) then 
             write(LVT_logunit,*) '[ERR] The distribution of surface types in '
             write(LVT_logunit,*) '[ERR] the LANDCOVER field does not sum to 100%'
             write(LVT_logunit,*) '[ERR] mask ',LVT_LMLC(source)%dommask(c,r)
             write(LVT_logunit,*) '[ERR] veg ',LVT_LMLC(source)%landcover(c,r,:)
             call LVT_endrun()
          endif
       enddo
    enddo
    
    write(LVT_logunit,*) '[INFO] Total number of grid points',LVT_LIS_rc(source)%ngrid

    allocate(LVT_LIS_domain(source)%tile(LVT_LIS_rc(source)%ntiles))
    do m=1,LVT_rc%max_model_types
       allocate(LVT_surface(source,m)%lis_tile(LVT_LIS_rc(source)%npatch(m)))
    enddo
    allocate(LVT_LIS_domain(source)%grid(LVT_LIS_rc(source)%ngrid))
    allocate(LVT_LIS_domain(source)%gindex(LVT_LIS_rc(source)%gnc, LVT_LIS_rc(source)%gnr))
    
    LVT_domain%minLat = 200.0
    LVT_domain%minLon = 200.00
    LVT_domain%maxLat = -200.0
    LVT_domain%maxLon = -200.00
    
    kk = 1
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc
          LVT_LIS_domain(source)%gindex(c,r) = -1

          call ij_to_latlon(LVT_domain%lvtproj,float(c),float(r),&
               locallat,locallon)
          
          if(localLat.lt.LVT_domain%minLat) LVT_domain%minLat = localLat
          if(localLat.gt.LVT_domain%maxLat) LVT_domain%maxLat = localLat
          if(localLon.lt.LVT_domain%minLon) LVT_domain%minLon = localLon
          if(localLon.gt.LVT_domain%maxLon) LVT_domain%maxLon = localLon

          if(LVT_LMLC(source)%dommask(c,r).gt.0.and.&
               abs(sum(LVT_LMLC(source)%landcover(c,r,:))-1.0).lt.0.001.and.&
               ifContainsSelectedSurfaceType(source,&
               nint(LVT_LMLC(source)%surfacetype(c,r,:)))) then 
             LVT_LIS_domain(source)%grid(kk)%lat = locallat
             LVT_LIS_domain(source)%grid(kk)%lon = locallon
             LVT_LIS_domain(source)%grid(kk)%col = c
             LVT_LIS_domain(source)%grid(kk)%row = r
             LVT_LIS_domain(source)%gindex(c,r) = kk
             kk = kk+1
          endif
       enddo
    enddo

    kk = 0
    kk_sf = 0 
    do r=1,LVT_LIS_rc(source)%gnr
       do c=1,LVT_LIS_rc(source)%gnc
          do m=1,LVT_LIS_rc(source)%nensem
             if(LVT_LMLC(source)%dommask(c,r).gt.0) then
                
                ntiles_surface = compute_ntiles_surface(source,c,r)
                ntiles_soilt = compute_ntiles_soilt(source,c,r,soilt_selected)
                ntiles_soilf = compute_ntiles_soilf(source,c,r,soilf_selected)
                ntiles_elev = compute_ntiles_elev(source,c,r,elev_selected)
                ntiles_slope = compute_ntiles_slope(source,c,r,slope_selected)
                ntiles_aspect = compute_ntiles_aspect(source,c,r,aspect_selected)

                if(soilf_selected) then 
                   ntiles_soil = ntiles_soilf
                else
                   ntiles_soil = ntiles_soilt
                endif

                do iv = 1, ntiles_surface
                   call get_vegt_value(source,c,r,iv,vegt)
                   call get_surface_value(source,c,r,iv,sf_index)
                   do it = 1, ntiles_soil

                      sand = -1.0
                      clay = -1.0
                      silt = -1.0
                      soilt = -1.0

                      if(soilf_selected) then 
                         call get_soilf_value(source,c,r,it,soilf_selected,&
                              sand,clay,silt,soilf_index)
                      else
                         call get_soilt_value(source,c,r,it,soilt_selected,soilt)
                      endif

                      do ie = 1, ntiles_elev
                         call get_elev_value(source,c,r,ie,elev_selected,&
                              elev, elev_index)
                         do is=1, ntiles_slope
                            call get_slope_value(source,c,r,is,slope_selected,&
                                 slope, slope_index)
                            do ia=1, ntiles_aspect
                               call get_aspect_value(source,c,r,ia,aspect_selected,&
                                    aspect, aspect_index)
                            
                               kk = kk+1

                               LVT_LIS_domain(source)%tile(kk)%ensem = m
                               LVT_LIS_domain(source)%tile(kk)%row=r    
                               LVT_LIS_domain(source)%tile(kk)%col=c    
                               LVT_LIS_domain(source)%tile(kk)%index = &
                                    LVT_LIS_domain(source)%gindex(c,r)
                               LVT_LIS_domain(source)%tile(kk)%vegt=vegt
                               LVT_LIS_domain(source)%tile(kk)%soilt=soilt
                               LVT_LIS_domain(source)%tile(kk)%sand=sand
                               LVT_LIS_domain(source)%tile(kk)%clay=clay
                               LVT_LIS_domain(source)%tile(kk)%silt=silt
                               LVT_LIS_domain(source)%tile(kk)%sftype = &
                                    sf_index
                               LVT_LIS_domain(source)%tile(kk)%elev = &
                                    elev
                               LVT_LIS_domain(source)%tile(kk)%slope = &
                                    slope
                               LVT_LIS_domain(source)%tile(kk)%aspect = &
                                    aspect
                               LVT_LIS_domain(source)%tile(kk)%fgrd = &
                                    LVT_LMLC(source)%landcover(c,r,vegt)
                               
                               if(soilt_selected) then 
                                  LVT_LIS_domain(source)%tile(kk)%fgrd = & 
                                       LVT_LIS_domain(source)%tile(kk)%fgrd * &
                                       LVT_soils(source)%texture(c,r,soilt)
                               elseif(soilf_selected) then 
                                  LVT_LIS_domain(source)%tile(kk)%fgrd = & 
                                       LVT_LIS_domain(source)%tile(kk)%fgrd * &
                                       LVT_soils(source)%soilffgrd(c,r,soilf_index)
                               endif
                               
                               if(elev_selected) then 
                                  LVT_LIS_domain(source)%tile(kk)%fgrd=& 
                                       LVT_LIS_domain(source)%tile(kk)%fgrd*&
                                       LVT_topo(source)%elevfgrd(c,r,elev_index)
                               endif

                               if(slope_selected) then 
                                  LVT_LIS_domain(source)%tile(kk)%fgrd=& 
                                       LVT_LIS_domain(source)%tile(kk)%fgrd*&
                                       LVT_topo(source)%slopefgrd(c,r,slope_index)
                               endif
                            
                               if(aspect_selected) then 
                                  LVT_LIS_domain(source)%tile(kk)%fgrd=& 
                                       LVT_LIS_domain(source)%tile(kk)%fgrd*&
                                       LVT_topo(source)%aspectfgrd(c,r,aspect_index)
                               endif
                               
                               if(sf_index.gt.0) then 
                                  kk_sf(sf_index) = kk_sf(sf_index) + 1
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%ensem = m
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%tile_id = kk
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%row=r    
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%col=c    
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%index = &
                                       LVT_LIS_domain(source)%gindex(c,r)
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%vegt=vegt
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%soilt=soilt
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%sand=sand
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%clay=clay
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%silt=silt

                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%elev=elev
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%slope=slope
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%aspect=aspect
                                  LVT_surface(source,sf_index)%lis_tile(&
                                       kk_sf(sf_index))%fgrd=&
                                       LVT_LMLC(source)%landcover(c,r,vegt)
                                  if(soilt_selected) then 
                                     LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LVT_soils(source)%texture(c,r,soilt)
                                  elseif(soilf_selected) then 
                                     LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LVT_soils(source)%soilffgrd(c,r,soilf_index)
                                  endif

                                  if(elev_selected) then 
                                     LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LVT_topo(source)%elevfgrd(c,r,elev_index)
                                  endif
                                  if(slope_selected) then 
                                     LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LVT_topo(source)%slopefgrd(c,r,slope_index)
                                  endif
                                  if(aspect_selected) then 
                                     LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd = & 
                                          LVT_surface(source,sf_index)%lis_tile(&
                                          kk_sf(sf_index))%fgrd * & 
                                          LVT_topo(source)%aspectfgrd(c,r,aspect_index)
                                  endif                                  
                               endif

                            enddo
                         end do
                      end do
                   end do
                enddo
             end if
          end do
       end do
    enddo

    call LVT_domain_setup(source)

  end subroutine create_tilespace

!BOP
! 
! !ROUTINE: isSurfaceTypeSelected
! \label{isSurfaceTypeSelected}
! 
! !INTERFACE: 
  function isSurfaceTypeSelected(source, surface_type)
! !ARGUMENTS:     
    integer       :: source
    integer       :: surface_type
! 
! This function determines if the input surface type is among 
! the surface model types that are currently chosen in the 
! LVT run
!
!EOP
    logical       :: isSurfaceTypeSelected
    integer       :: i 

    isSurfaceTypeSelected = .false. 
    
    do i=1,LVT_LIS_rc(source)%nsf_model_types
       if(LVT_LIS_rc(source)%sf_model_type_select(i).eq.surface_type) then 
          isSurfaceTypeSelected = .true. 
          exit
       endif
    enddo

  end function isSurfaceTypeSelected

!BOP
! 
! !ROUTINE: ifContainsSelectedSurfaceType
! \label{ifContainsSelectedSurfaceType}
! 
! !INTERFACE: 
  function ifContainsSelectedSurfaceType(source,surface_types)
! !ARGUMENTS:     
    integer       :: source
    integer       :: surface_types(LVT_LIS_rc(source)%nsurfacetypes)
! 
! !DESCRIPTION: 
!  This function determines if any of the input list of surface types
!  is among the surface model types that are currently chosen in the 
!  LVT run
!EOP
    logical       :: ifContainsSelectedSurfaceType
    integer       :: i,k

    ifContainsSelectedSurfaceType = .false. 
    
    do k=1,LVT_LIS_rc(source)%nsurfacetypes       
       do i=1,LVT_LIS_rc(source)%nsf_model_types
          if(LVT_LIS_rc(source)%sf_model_type_select(i).eq.surface_types(k)) then 
             ifContainsSelectedSurfaceType = .true. 
             exit
          endif
       enddo
    enddo

  end function ifContainsSelectedSurfaceType

!BOP
! 
! !ROUTINE: compute_ntiles_surface
! \label{compute_ntiles_surface}
! 
! !INTERFACE: 
  function compute_ntiles_surface(source,c,r)
! !ARGUMENTS: 
    integer :: n 
    integer :: source
    integer :: c
    integer :: r
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of landcover types in a given grid cell
!EOP

    integer :: compute_ntiles_surface
    integer :: t
    
    compute_ntiles_surface = 0 
    do t=1,LVT_LIS_rc(source)%nsurfacetypes        
       if(LVT_LMLC(source)%dommask(c,r).gt.0.and.&
            LVT_LMLC(source)%landcover(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(source, &
            nint(LVT_LMLC(source)%surfacetype(c,r,t)))) then
          compute_ntiles_surface = compute_ntiles_surface + 1
       endif
    enddo
  end function compute_ntiles_surface

!BOP
! !ROUTINE: compute_npatches_surface
! \label{compute_npatches_surface}
! 
! !INTERFACE: 
  function compute_npatches_surface(source,c,r,t,sf_index)
! !ARGUMENTS: 
    integer :: n 
    integer :: source
    integer :: c
    integer :: r
    integer :: t
    integer :: sf_index
!
! !DESCRIPTION: 
!  This function computes the number of patches (for the given 
!  surface type) based on the distribution of landcover types 
!  in a given grid cell
!EOP
    integer :: compute_npatches_surface
    
    compute_npatches_surface = 0 
    if(LVT_LMLC(source)%dommask(c,r).gt.0.and.&
         LVT_LMLC(source)%landcover(c,r,t).gt.0.0.and.&
         isSurfaceTypeSelected(source, &
         nint(LVT_LMLC(source)%surfacetype(c,r,t)))) then
       compute_npatches_surface = compute_npatches_surface + 1
    endif
  end function compute_npatches_surface

!BOP
! !ROUTINE: compute_ntiles_soilt
! \label{compute_ntiles_soilt}
! 
! !INTERFACE: 
  function compute_ntiles_soilt(source,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: source
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of soil texture types in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_soilt
    integer :: t

    kk = 0 
    if(flag) then 
       do t=1,LVT_LIS_rc(source)%nsoiltypes        
          if(LVT_soils(source)%texture(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_soilt = max(1,kk)

  end function compute_ntiles_soilt

!BOP
! !ROUTINE: compute_ntiles_soilf
! \label{compute_ntiles_soilf}
! 
! !INTERFACE: 
  function compute_ntiles_soilf(source,c,r,flag)
! !ARGUMENTS: 
    integer :: source
    integer :: n 
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of soil fraction data in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_soilf
    integer :: t

    kk = 0 
    if(flag) then 
       do t=1,LVT_LIS_rc(source)%nsoilfbands
          if(LVT_soils(source)%soilffgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_soilf = max(1,kk)

  end function compute_ntiles_soilf

!BOP
!
! !ROUTINE: compute_ntiles_elev
! \label{compute_ntiles_elev}
! 
! !INTERFACE: 
  function compute_ntiles_elev(source,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: source
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of elevation bands in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_elev
    integer :: t

    kk = 0 
    
    if(flag) then 
       do t=1,LVT_LIS_rc(source)%nelevbands        
          if(LVT_topo(source)%elevfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_elev = max(1,kk)

  end function compute_ntiles_elev

!BOP
!
! !ROUTINE: compute_ntiles_slope
! \label{compute_ntiles_slope}
! 
! !INTERFACE: 
  function compute_ntiles_slope(source,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: source
    integer :: c
    integer :: r
    logical :: flag
!
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of slope bands in a given grid cell
!EOP 
    integer :: kk 
    integer :: compute_ntiles_slope
    integer :: t

    kk = 0 
    
    if(flag) then 
       do t=1,LVT_LIS_rc(source)%nslopebands        
          if(LVT_topo(source)%slopefgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_slope = max(1,kk)

  end function compute_ntiles_slope

!BOP
!
! !ROUTINE: compute_ntiles_aspect
! \label{compute_ntiles_aspect}
! 
! !INTERFACE: 
  function compute_ntiles_aspect(source,c,r,flag)
! !ARGUMENTS: 
    integer :: n 
    integer :: source
    integer :: c
    integer :: r
    logical :: flag
!
! !DESCRIPTION: 
!  This function computes the number of tiles based on the 
!  distribution of aspect bands in a given grid cell
!EOP
    integer :: kk 
    integer :: compute_ntiles_aspect
    integer :: t

    kk = 0 
    
    if(flag) then 
       do t=1,LVT_LIS_rc(source)%naspectbands        
          if(LVT_topo(source)%aspectfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
          endif
       enddo
    endif
    compute_ntiles_aspect = max(1,kk)

  end function compute_ntiles_aspect

!BOP
! 
! !ROUTINE: get_vegt_value
! \label{get_vegt_value}
! 
! !INTERFACE: 
  subroutine get_vegt_value(source,c,r,i,vegt)
! !ARGUMENTS: 
    integer  :: n 
    integer :: source
    integer  :: c
    integer  :: r
    integer  :: i
    integer  :: vegt
! 
! !DESCRIPTION: 
!  This routine computes the vegetation type value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk

    kk = 0 
    do t=1, LVT_LIS_rc(source)%nsurfacetypes
       if(LVT_LMLC(source)%landcover(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(source, &
            nint(LVT_LMLC(source)%surfacetype(c,r,t)))) then 
          kk = kk + 1
          if(kk.eq.i) then 
             vegt = t
             return
          endif
       endif
    enddo
  end subroutine get_vegt_value

!BOP
! !ROUTINE: get_surface_value
! \label{get_surface_value}
! 
! !INTERFACE:
  subroutine get_surface_value(source,c,r,i,sf_index)
    integer  :: source
    integer  :: n 
    integer  :: c
    integer  :: r
    integer  :: i
    integer  :: sf_index
! 
! !DESCRIPTION: 
!  This routine computes the surface type value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk

    kk = 0 
    do t=1, LVT_LIS_rc(source)%nsurfacetypes
       if(LVT_LMLC(source)%landcover(c,r,t).gt.0.0.and.&
            isSurfaceTypeSelected(source, &
            nint(LVT_LMLC(source)%surfacetype(c,r,t)))) then 
          kk = kk + 1
          if(kk.eq.i) then 
             sf_index = nint(LVT_LMLC(source)%surfacetype(c,r,t))
             return
          endif
       endif
    end do
  end subroutine get_surface_value

!BOP
! !ROUTINE: get_soilt_value
! \label{get_soilt_value}
!
! !INTERFACE:
  subroutine get_soilt_value(source,c,r,i,soilt_selected,soilt)
    integer  :: n 
    integer :: source
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: soilt_selected
    integer  :: soilt
! 
! !DESCRIPTION: 
!  This routine computes the soil texture value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    soilt = -1
    if(soilt_selected) then 
       kk = 0 
       do t=1, LVT_LIS_rc(source)%nsoiltypes
          if(LVT_soils(source)%texture(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                soilt = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_soilt_value

!BOP
! !ROUTINE: get_soilf_value
! \label{get_soilf_value}
!
! !INTERFACE:
  subroutine get_soilf_value(source,c,r,i,soilf_selected,sand,clay,silt,soilf_index)
    integer  :: n 
    integer :: source
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: soilf_selected
    real     :: sand
    real     :: clay
    real     :: silt
    integer  :: soilf_index
! 
! !DESCRIPTION: 
!  This routine computes the elevation value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    sand = -1
    clay = -1
    silt = -1
    soilf_index = -1
    if(soilf_selected) then 
       kk = 0 
       do t=1, LVT_LIS_rc(source)%nsoilfbands
          if(LVT_soils(source)%soilffgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                sand = LVT_soils(source)%sand(c,r,t)
                clay = LVT_soils(source)%clay(c,r,t)
                silt = LVT_soils(source)%silt(c,r,t)
                soilf_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_soilf_value

!BOP
! !ROUTINE: get_elev_value
! \label{get_elev_value}
!
! !INTERFACE:
  subroutine get_elev_value(source,c,r,i,elev_selected,elev,elev_index)
     use LVT_logMod, only: LVT_logunit ! EMK

    integer  :: n 
    integer :: source
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: elev_selected
    real     :: elev
    integer  :: elev_index
! 
! !DESCRIPTION: 
!  This routine computes the elevation value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk

    ! EMK BEGIN
!    write(LVT_logunit,*) 'EMK: shape(LVT_topo(source)%elevfgrd) = ', &
!         shape(LVT_topo(source)%elevfgrd)
!    write(LVT_logunit,*) 'EMK: shape(LVT_topo(source)%elevation) = ', &
!         shape(LVT_topo(source)%elevation)
    ! EMK END

    elev = -1
    elev_index = -1
    if(elev_selected) then 
       kk = 0 
       do t=1, LVT_LIS_rc(source)%nelevbands
          if(LVT_topo(source)%elevfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                elev = LVT_topo(source)%elevation(c,r,t)
                elev_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_elev_value

!BOP
! 
! !ROUTINE: get_slope_value
! \label{get_slope_value}
!
! !INTERFACE:
  subroutine get_slope_value(source,c,r,i,slope_selected,slope,slope_index)
    integer  :: n
    integer :: source 
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: slope_selected
    real     :: slope
    integer  :: slope_index
! 
! !DESCRIPTION: 
!  This routine computes the slope value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    slope = -1
    slope_index = -1
    if(slope_selected) then 
       kk = 0 
       do t=1, LVT_LIS_rc(source)%nslopebands
          if(LVT_topo(source)%slopefgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                slope = LVT_topo(source)%slope(c,r,t)
                slope_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_slope_value

!BOP
! !ROUTINE: get_aspect_value
! \label{get_aspect_value}
! 
! !INTERFACE: 
  subroutine get_aspect_value(source,c,r,i,aspect_selected,aspect,aspect_index)
! !ARGUMENTS: 
    integer  :: n 
    integer :: source
    integer  :: c
    integer  :: r
    integer  :: i
    logical  :: aspect_selected
    real     :: aspect
    integer  :: aspect_index
! 
! !DESCRIPTION: 
!  This routine computes the aspect value for given tile number
!  and grid cell. 
!EOP
    integer  :: t
    integer  :: kk
    
    aspect = -1
    aspect_index = -1
    if(aspect_selected) then 
       kk = 0 
       do t=1, LVT_LIS_rc(source)%naspectbands
          if(LVT_topo(source)%aspectfgrd(c,r,t).gt.0.0) then 
             kk = kk + 1
             if(kk.eq.i) then 
                aspect = LVT_topo(source)%aspect(c,r,t)
                aspect_index = t
                return
             endif
          end if
       enddo
    endif

  end subroutine get_aspect_value

!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: LVT_readDataMask
! \label{LVT_readDataMask}
! r
! !INTERFACE: 
  subroutine LVT_readDataMask
! !USES:     

    implicit none

! 
! !DESCRIPTION: 
!   This routine reads the external data mask, which is used to 
!   screen grid points, both spatially and temporally. 
!EOP

    character*100 :: maskfile
    logical       :: file_exists
    real          :: datamask(LVT_rc%lnc, LVT_rc%lnr)
    integer       :: ftn
    integer       :: c,r


!read the external mask, timestamped file      
    if(LVT_rc%datamask.eq.0) then 
       LVT_stats%datamask = 1
    elseif(LVT_rc%dataMask.eq.1) then  
       call create_datamask_filename(maskfile)
       inquire(file=(maskfile), exist=file_exists)
       if(file_exists) then 
          write(LVT_logunit,*) '[INFO] Reading mask ',(maskfile)
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=(maskfile), form='unformatted')

          read(ftn) datamask

          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(datamask(c,r).ne.LVT_rc%udef) then 
                   LVT_stats%datamask(c,r) = 1
                else
                   LVT_stats%datamask(c,r) = 0
                endif
             enddo
          enddo
          call LVT_releaseUnitNumber(ftn)          
       else
          LVT_stats%datamask = 0
       endif
    elseif(LVT_rc%dataMask.eq.2.and.LVT_rc%maskflag.eq.0) then  
!read the static mask file
       LVT_rc%maskflag = 1
       inquire(file=(LVT_rc%maskdir), exist=file_exists)

       if(file_exists) then 
          write(LVT_logunit,*) '[INFO] Reading mask ',(LVT_rc%maskdir)
          ftn = LVT_getNextUnitNumber()
          open(ftn,file=(LVT_rc%maskdir), form='unformatted')
 
          read(ftn) datamask
          
          do r=1,LVT_rc%lnr
             do c=1,LVT_rc%lnc
                if(datamask(c,r).ne.LVT_rc%udef) then 
                   LVT_stats%datamask(c,r) = 1
                else
                   LVT_stats%datamask(c,r) = 0
                endif
             enddo
          enddo
          call LVT_releaseUnitNumber(ftn)
       else
          write(LVT_logunit,*) '[ERR] Mask file not found ',(LVT_rc%maskdir)
          stop
          LVT_stats%datamask = 0 
       endif
    elseif(LVT_rc%dataMask.eq.3) then 
       if(LVT_rc%monthly_mask(LVT_rc%mo).eq.1) then 
          LVT_stats%datamask = 1
       else
          LVT_stats%datamask = 0 
       endif
    endif
    
  end subroutine LVT_readDataMask


!BOP
! 
! !ROUTINE:
!
! !INTERFACE:
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ROUTINE: create_datamask_filename
! \label{create_datamask_filename}
! 
! !INTERFACE: 
  subroutine create_datamask_filename(maskfile)
! !USES: 

    implicit none

! 
! !DESCRIPTION: 
!  This subroutine creates a timestamped filename for the external 
!  mask data
!EOP

    character(len=*)      :: maskfile

    character*12          :: fdate

    write(unit=fdate, fmt='(i4.4,i2.2,i2.2,i2.2,i2.2)') LVT_rc%yr, &
         LVT_rc%mo, LVT_rc%da, LVT_rc%hr, LVT_rc%mn

    maskfile = trim(LVT_rc%maskdir)//&
         '/'//trim(fdate)//'.1gs4r'

  end subroutine create_datamask_filename

end module LVT_domainMod

