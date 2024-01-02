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
! !ROUTINE: setmosp
! \label{setmosp}
!
! !REVISION HISTORY:
!  31 Jul 2001: Matt Rodell; Initial Specification
!  14 Feb 2002: Jon Gottschalck; Added allocated space for AVHRR LAI/DSAI     
!  07 Mar 2002: Brian Cosgrove; Corrected declaration of TEX1 var from real to int
!  14 Jan 2003: Urszula Jambor; Added conditional to check if need exists 
!               to allocate for AVHRR LAI/DSAI variables.
!  25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
!
! !INTERFACE:
subroutine setmosp(n)
! !USES:          
  use LIS_coreMod,  only : LIS_rc, LIS_surface
  use LIS_soilsMod, only : LIS_soils
  use LIS_logMod,   only : LIS_logunit, LIS_endrun
  use LIS_fileIOMod, only : LIS_read_param
  use mos_lsmMod      

   implicit none
! !ARGUMENTS: 
   integer, intent(in) :: n 
!
! !DESCRIPTION:
!  This subroutine sets the static mosaic parameters related to vegetation, 
!  soils and topography. 
!
!EOP
!----------------------------------------------------------------------------
!  Minimum values of sin(theta) based on vegetation type.
!----------------------------------------------------------------------------
   REAL, PARAMETER :: S8 = 0.57787      ! Closed shrubland
   REAL, PARAMETER :: S9 = 0.95504      ! Open shrubland
   REAL, PARAMETER :: S12 = 0.1736      ! Bare soil
   REAL, PARAMETER :: S0 = 0.05         ! All others 

!----------------------------------------------------------------------------
!  Maximum allowable porosity.
!----------------------------------------------------------------------------
   REAL, PARAMETER :: PORMAX=0.70
   integer, parameter   :: nt = 13
   integer              :: c,r,i,j,k,jj    ! loop counters
   real                 :: value(nt,mos_struc(n)%mos_nvegp)
   real                 :: basicset(nt,mos_struc(n)%mos_nsoilp)
   real, allocatable    :: soilpset(:,:)
   integer, allocatable :: tex1(:)
   real                 :: por1a,por2a,por3a
   real                 :: slope(LIS_rc%lnc(n), LIS_rc%lnr(n))
   real                 :: dsoil(LIS_rc%lnc(n), LIS_rc%lnr(n))
   integer              :: tid
   integer              :: mapVegToUMD
!----------------------------------------------------------------------------
! Get Vegetation Parameters for Mosaic Model in Tile Space
! Read in the Mosaic Static and Monthly Vegetation Parameter Files
!----------------------------------------------------------------------------

   allocate(tex1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
   open(unit=15,file=mos_struc(n)%mos_vfile,status='old')
   
   do j=1,mos_struc(n)%mos_nvegp
      read(15,*)
      read(15,*)(value(i,j),i=1,nt)
   enddo
   close(15)
   
!----------------------------------------------------------------------------
!  Assign STATIC vegetation parameters to each tile based on the
!  type of vegetation present in that tile.
!  These parameters will be stored in one long array--structured
!  as follows: Tile 1, all the parameters (1 through numparam)
!  then Tile 2, all the parameters. 
!  Then Tile 3, all the parameters etc.
!----------------------------------------------------------------------------
   do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      mos_struc(n)%mos(i)%vegt = mapVegToUMD(&
           LIS_surface(n,LIS_rc%lsm_index)%tile(i)%vegt,&
           LIS_rc%lcscheme)
   enddo
   do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
      do j=1,mos_struc(n)%mos_nvegp
         mos_struc(n)%mos(i)%vegp(j)=value(&
              mos_struc(n)%mos(i)%vegt,j)
      enddo !j
   enddo !i
   
!---------------------------------------------------------------------------
!  set the soils
!---------------------------------------------------------------------------
   if(LIS_rc%usetexturemap(n).eq."none".and.&
        LIS_rc%usesoilfractionmap(n).eq."none") then       
!----------------------------------------------------------------------------
!   Get Soil Parameters (Based on Vegetation) for Mosaic Model in Tile Space
!   Read in Soil Parameter Data
!----------------------------------------------------------------------------
      open(10,file=mos_struc(n)%mos_sfile,status='old', &
           access='sequential')
      
      do i=1,mos_struc(n)%mos_nsoilp
         read(10,*)(basicset(jj,i),jj=1,LIS_rc%nvegtypes)
      enddo
      close(10)
      
      do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
!         k=LIS_surface(n,LIS_rc%lsm_index)%tile(i)%vegt
         k=mos_struc(n)%mos(i)%vegt
         do j=1,mos_struc(n)%mos_nsoilp
            mos_struc(n)%mos(i)%soilp(j)=basicset(k,j)                  
         enddo !j 
      enddo !i 
   end if  !soil=0

   if(LIS_rc%usetexturemap(n).ne."none".or.&
        LIS_rc%usesoilfractionmap(n).ne."none") then     
     allocate(soilpset(mos_struc(n)%mos_nstxts,10))
     write(unit=LIS_logunit,fmt=*)' setmosp -- reading soil files'

     if(LIS_rc%usetexturemap(n).eq."none") then
        tex1 = -9999
        do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           call texture (&
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%sand, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%clay, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%silt, &
                tex1(i))
        enddo
     endif

     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)

         if(LIS_rc%usetexturemap(n).eq."none") then                
            mos_struc(n)%mos(i)%soiltype = tex1(i)
         else
            mos_struc(n)%mos(i)%soiltype = &
                 LIS_surface(n,LIS_rc%lsm_index)%tile(i)%soilt
         endif
           
         if (mos_struc(n)%mos(i)%soiltype .eq. -9999 ) then
            mos_struc(n)%mos(i)%soiltype = 6 
         endif

         if ( mos_struc(n)%mos(i)%soiltype .eq. 14 .and. &
              mos_struc(n)%mos(i)%vegt .ne. LIS_rc%waterclass ) then
            mos_struc(n)%mos(i)%soiltype = 7
         endif
     enddo

     !----------------------------------------
     ! Read in the Mosaic Soil Parameter File
     !----------------------------------------
     write(LIS_logunit,*) 'Reading soil parameter table ',&
          mos_struc(n)%mos_sfile

     open(unit=18,file=mos_struc(n)%mos_sfile,status='old',access='sequential')
     do i=1,10
        read(18,*)
        read(18,*)(soilpset(jj,i),jj=1,mos_struc(n)%mos_nstxts)
     enddo
     close(18)

     !-----------------------------------------------------------------------
     ! Assign SOIL Parameters to each tile based on the
     ! type of soil class present in that tile.
     !-----------------------------------------------------------------------

     ! Set tile soil b-parameter
     if (LIS_rc%usebexpmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           mos_struc(n)%mos(i)%soilp(1) = soilpset(mos_struc(n)%mos(i)%soiltype,1)
        end do
     else
        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           mos_struc(n)%mos(i)%soilp(1) = &
              lis_soils(n)%bexp(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
        end do
     end if

     ! Set tile soil Psi-sat
     if (LIS_rc%usepsisatmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           mos_struc(n)%mos(i)%soilp(2) = soilpset(mos_struc(n)%mos(i)%soiltype,2)
        end do
     else
        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           mos_struc(n)%mos(i)%soilp(2) = &
              lis_soils(n)%psisat(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
        end do
     end if

     ! Set tile soil K-sat
     if (LIS_rc%useksatmap(n).eq."none") then ! default, from look-up table
        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           mos_struc(n)%mos(i)%soilp(3) = soilpset(mos_struc(n)%mos(i)%soiltype,3)
        end do
     else
        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           mos_struc(n)%mos(i)%soilp(3) = &
              lis_soils(n)%ksat(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)
        end do
     end if

     ! Set tile soil depths
     do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        mos_struc(n)%mos(i)%soilp(4) = soilpset(mos_struc(n)%mos(i)%soiltype,4)
        mos_struc(n)%mos(i)%soilp(5) = soilpset(mos_struc(n)%mos(i)%soiltype,5)
        mos_struc(n)%mos(i)%soilp(6) = soilpset(mos_struc(n)%mos(i)%soiltype,6)
     end do
     if(mos_struc(n)%usedsoilmap.ne.0) then 

        call LIS_read_param(n,"SOILDEPTH",dsoil)
        
        do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
           c = LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col
           r = LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row

           if((dsoil(c,r)- mos_struc(n)%mos(i)%soilp(5)).gt.0) then 
              mos_struc(n)%mos(i)%soilp(6) = dsoil(c,r)
           endif
        enddo
     endif
     ! Set tile slope values
     call LIS_read_param(n,"SLOPE",slope)

     do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
        select case (mos_struc(n)%mos(i)%vegt)
        case (8)        ! closed shrubland
           mos_struc(n)%mos(i)%soilp(7) = &
                max(s8, sin(slope(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)))
        case (9)        ! open shrubland
           mos_struc(n)%mos(i)%soilp(7) = &
                max(s9, sin(slope(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)))
        case (12)       ! bare ground
           mos_struc(n)%mos(i)%soilp(7) = &
                max(s12, sin(slope(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)))
        case default    ! all other vegetation classes
           mos_struc(n)%mos(i)%soilp(7) = &
                max(s0, sin(slope(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col, &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row)))
        end select

        if (mos_struc(n)%mos(i)%soilp(7) .le. 0.0) then 
           write(LIS_logunit,*) 'col,row,LIS_domain(n)%tile,sin(theta)', &
                LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row, &
                i,mos_struc(n)%mos(i)%soilp(7)
        end if
     enddo

     ! Set tile soil maximum water content
     if (LIS_rc%useporositymap(n).eq."none") then ! default, from look-up table
        write(LIS_logunit,*) " setmosp -- no porosity datasource selected"
        write(LIS_logunit,*) " setmosp -- calling endrun"
        call LIS_endrun()
     else

        do i = 1,LIS_rc%npatch(n,LIS_rc%lsm_index)
            por1a = &
               lis_soils(n)%porosity(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row,1)
            por2a = &
               lis_soils(n)%porosity(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row,2)
            por3a = &
               lis_soils(n)%porosity(LIS_surface(n,LIS_rc%lsm_index)%tile(i)%col,LIS_surface(n,LIS_rc%lsm_index)%tile(i)%row,3)
            por1a=amin1(por1a,pormax)
            por2a=amin1(por2a,pormax)
            por3a=amin1(por3a,pormax)
            mos_struc(n)%mos(i)%soilp(8)=por1a*mos_struc(n)%mos(i)%soilp(4)*1000.00
            mos_struc(n)%mos(i)%soilp(9)=por2a*mos_struc(n)%mos(i)%soilp(5)*1000.00
            mos_struc(n)%mos(i)%soilp(10)=por3a*mos_struc(n)%mos(i)%soilp(6)*1000.00
        end do
     end if

     deallocate(soilpset)

  endif

  deallocate(tex1)
  return

end subroutine setmosp

function mapVegToUMD(vegt_in, classification)
  
  use LIS_logMod

  implicit none

  integer           :: vegt_in
  character(len=*)  :: classification
  integer           :: mapVegToUMD

  if(classification.eq."UMD") then !UMD, do nothing
     mapVegToUMD = vegt_in
  elseif(classification.eq."ECOCLIMAP2") then 
     if(vegt_in.eq.1) then
        mapVegToUMD = 12
     elseif(vegt_in.eq.2) then 
        mapVegToUMD = 12
     elseif(vegt_in.eq.3) then 
        mapVegToUMD = 12
     elseif(vegt_in.eq.4) then 
        mapVegToUMD = 4
     elseif(vegt_in.eq.5) then 
        mapVegToUMD = 3
     elseif(vegt_in.eq.6) then 
        mapVegToUMD = 1
     elseif(vegt_in.eq.7) then 
        mapVegToUMD = 11
     elseif(vegt_in.eq.8) then 
        mapVegToUMD = 11
     elseif(vegt_in.eq.9) then 
        mapVegToUMD = 11
     elseif(vegt_in.eq.10) then 
        mapVegToUMD = 10
     elseif(vegt_in.eq.11) then 
        mapVegToUMD = 10
     elseif(vegt_in.eq.12) then 
        mapVegToUMD = 7
     endif
  elseif(classification.eq."USGS") then 
     write(LIS_logunit,*) 'USGS classification not supported for Mosaic'
     call LIS_endrun()
  elseif(classification.eq."MODIS") then 
     write(LIS_logunit,*) 'MODIS classification not supported for Mosaic'
     call LIS_endrun()
  elseif(classification.eq."IGBPNCEP") then 
     if(vegt_in.eq.1) then
        mapVegToUMD = 1
     elseif(vegt_in.eq.2) then 
        mapVegToUMD = 2
     elseif(vegt_in.eq.3) then 
        mapVegToUMD = 3
     elseif(vegt_in.eq.4) then 
        mapVegToUMD = 4
     elseif(vegt_in.eq.5) then 
        mapVegToUMD = 5
     elseif(vegt_in.eq.6) then 
        mapVegToUMD = 8
     elseif(vegt_in.eq.7) then 
        mapVegToUMD = 9
     elseif(vegt_in.eq.8) then 
        mapVegToUMD = 7
     elseif(vegt_in.eq.9) then 
        mapVegToUMD = 10
     elseif(vegt_in.eq.10) then 
        mapVegToUMD = 10
     elseif(vegt_in.eq.11) then 
        mapVegToUMD = 10
     elseif(vegt_in.eq.12) then 
        mapVegToUMD = 11
     elseif(vegt_in.eq.13) then 
        mapVegToUMD = 13
     elseif(vegt_in.eq.14) then 
        mapVegToUMD = 11
     elseif(vegt_in.eq.15) then 
        mapVegToUMD = 12
     elseif(vegt_in.eq.16) then 
        mapVegToUMD = 12
     elseif(vegt_in.eq.17) then 
        mapVegToUMD = 12
     elseif(vegt_in.eq.18) then 
        mapVegToUMD = 12
     elseif(vegt_in.eq.19) then 
        mapVegToUMD = 12
     endif
  endif
  
end function mapVegToUMD
