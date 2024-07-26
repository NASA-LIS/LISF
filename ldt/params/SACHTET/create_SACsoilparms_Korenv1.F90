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
! !ROUTINE: create_SACsoilparms_Korenv1
! \label{create_SACsoilparms_Korenv1}
!
! !REVISION HISTORY:
!  10 Mar 2015: KR Arsenault; Modified to support SACHTET parameters
!
! !INTERFACE:
 subroutine create_SACsoilparms_Korenv1( n,     &
               soiltext, hsg, cosbysoils_table, & 
               uzfwm, uztwm, lztwm, lzfsm, lzfpm, &
               uzk, lzsk, lzpk, zperc, rexp,    &
               pfree, frz_stxt, dzup1, dzlw, txtlw, txtot )      

! !USES:
  use ESMF
  use LDT_coreMod,  only : LDT_rc, LDT_domain
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
                           LDT_releaseUnitNumber, LDT_endrun
  use LDT_paramDataMod
  use SACHTET_parmsMod

  implicit none

! !ARGUMENTS: 
! Inputs:
  integer,        intent(in) :: n
  character(140), intent(in) :: cosbysoils_table
  type(LDT_paramEntry), intent(in) :: soiltext
  type(LDT_paramEntry), intent(in) :: hsg
! Outputs:
  type(LDT_paramEntry), intent(inout) :: uzfwm
  type(LDT_paramEntry), intent(inout) :: uztwm
  type(LDT_paramEntry), intent(inout) :: lztwm
  type(LDT_paramEntry), intent(inout) :: lzfsm
  type(LDT_paramEntry), intent(inout) :: lzfpm
  type(LDT_paramEntry), intent(inout) :: uzk
  type(LDT_paramEntry), intent(inout) :: lzsk
  type(LDT_paramEntry), intent(inout) :: lzpk
  type(LDT_paramEntry), intent(inout) :: zperc
  type(LDT_paramEntry), intent(inout) :: rexp
  type(LDT_paramEntry), intent(inout) :: pfree
  type(LDT_paramEntry), intent(inout) :: frz_stxt
!
  type(LDT_paramEntry), intent(inout) :: dzup1
  type(LDT_paramEntry), intent(inout) :: dzlw
  type(LDT_paramEntry), intent(inout) :: txtlw
  type(LDT_paramEntry), intent(inout) :: txtot
  
!
! !DESCRIPTION:
!  This subroutine calculates the soil parameters required by
!   the SAC-HTET land/hydrology model.  The input datasets and
!   equations are described mostly in Koren et al. (1999), 
!   Koren et al. (2000), Koren et al. (2003) and Zhang et al. (2011).
!  
!  The arguments are:
! \begin{description}
!  \item[n]
!   index of the nest
!  \item[soiltext]
!   input layered (depth) soil texture fields   
!  \item[hsg]
!   input layered (index) hydrologic soil groups (HSG) field
!  \item[cosbysoils_table]
!   input Cosby et al. soil parameters based on textures
! \end{description}
!
!EOP      

  integer :: i, j, c, r
  integer :: ftn
  logical :: file_exists
  integer :: icls
  integer :: num_rows, num_cols

  integer, parameter :: maxclass = 12  ! The number of USDA soil classes

  real    :: SWMAX(maxclass)     ! Porosity (soil saturation) for soil classes
  real    :: SWFLD(maxclass)     ! Field capacity for soil classes
  real    :: SWWLT(maxclass)     ! Wilting point for soil classes
  real    :: KS(maxclass)        ! Hydraulic conductivity for soil classes (Ksat)
  real    :: SINIT(maxclass)     ! "Iinit" 
  real    :: hhed(maxclass)      ! Capillary water potential or pressure height (mm),
                                 !  based on Rawl's/Maidment(1993) - See Koren(2003), Table 1)
  real    :: hcapil(maxclass)    ! Capillary pressure height for soil classes (mm),
                                 !  (Rawls or Cosby et al. (1984) -- ? )

! ZPERC option:  'old' | 'new'      
!  'new' option used to calculate ZPERC using capillary head presssure
  character*3, parameter :: izperc = 'new'

! ----------------------------------------------------------
!  INITIAL LOSSES FOR A, B, C, D SOIL GROUPS WERE ESTIMATED 
!  FOR 'PASTURE OR RANGE' UNDER 'FAIR' HYDROLOGIC CONDITIONS
!  AND DRY SOIL CONDITIONS. CN VALUES IN THIS CASE ARE:
!    70, 63, 48, 28 (for D, C, B, A HSG)
!  See equation 1 in Koren et al. (2000)
!  original:  REAL SINHSG(5)/80.,33.,18.,13.,0./

  real SINHSG(5)/134.5,55.2,30.2,21.6,0./
! ----------------------------------------------------------

  real :: SINX, SX  ! Initial total losses (in mm)
  real :: shsgx     ! Accumulated total losses associated with HSG (in mm)
  real :: ZMXD      ! Depth to bottom of the soil column (mm)
  real :: UZFM      ! Scalar upper zone free-water depth
  real :: ZUPX      ! Depth of upper zone = sum(depths of DZUP)
  real :: IPRZ      ! Indicates if impermeable/bedrock layer already encountered 
  real :: DWUP      ! Depth of UZFWM
  real :: DZUP      ! Incremental depth, summed to compute total depth of lower zone
  integer :: IVALX     ! Soil texture for current soil layer
  integer :: IVAL_PRV  ! Soil texture in some layer above the current being used.
  real :: DZ        ! Total thickness (mm) of current layer.
  integer :: ILR    ! Loop paramteer for STATSGO soil layers 
  integer :: N16    ! Counter for texture layers with "other" in gridcell
  integer :: IMPRB  ! Counter for number of impermeable layers in gridcell
  real :: ZLWX      ! Thickness of the lower zone 
  real :: SMU       ! Vol. soil moisture at saturation for the upper zone
  real :: SFU       ! Vol. soil moisture at field capacity for the upper zone
  real :: SWU       ! Vol. soil moisture at wilting point for the upper zone
  real :: SML       ! Vol. soil moisture at saturation for the lower zone
  real :: SFL       ! Vol. soil moisture at field capacity for the lower zone
  real :: SWL       ! Vol. soil moisture at wilting point for the lower zone
  real :: SKL       ! Saturated hydraulic conductivity for the lower zone
  real :: SKU       ! Saturated hydraulic conductivity for the upper zone
  real :: shed      ! Capillary pressure height for each texture class
  integer :: NU        ! Number of STATSGO soil layers in upper zone
  integer :: NL        ! Number of STATSGO soil layers in the lower zone
  integer :: NE     ! Bottom-most STATSGO layer in the lower zone
  real :: DZU       ! Thickness of the current layer in upper zone.  Used to
                    !  compute depth-weighted average of porosity, field cap., WP, Ksat
  real :: DZL       ! Thickness of the current layer in the lower zone.  Used to 
                    !  compute depth-weighted average of poros,field cap, WP, Ksat
  real :: ZSK       ! Withdrawal coef. for the lower zone free supplemental baseflow (aka, LZSK)
  real :: ZPK       ! Withdrawal coef. for the lower zone free primary baseflow (aka, LZPK)
  real :: ZFSM      ! 
  real :: ZFLM      !
  real :: ZFPM      !
  real :: SFREE     ! aka, PFREE; % of percolated water that bypasses the 
                    !  Lower Zone Free Primary baseflow reservoir and goes directly
                    !  to Lower Zone Free and Supplemental baseflow reservoirs.
  real :: PBASE     ! Maximum amount of water that can exit the base flow reservoirs
  real :: pbase_hr  ! PBASE on hour increment (divided by 24)
  real :: skmin     ! Min. hydraulic conductivity (usually with lower zone)
  real :: pcap      ! Capillary pressure gradient (see Eq. 13 in Koren, 2008?)
  real :: tcurve    ! Time curve (used in Zperc calculation; See Eq 13 in Koren, 2008?)
  real :: zperc1    ! Multiplier used to specify the percolation curve
  integer :: ITXTUP(15) ! Counter for soil texture frequency for upper layer of gridcell
  integer :: ITXTLW(15) ! Counter for soil texture frequency for lower layer
  integer :: ITXTOT(15) ! Counter for soil texture frequency for total
  
  integer :: NSUM

! Depth to the TOP of the 11 STATSGO layers
  integer LAYER(12)/0,50,100,200,300,400,600,800,1000,1500,2000,2500/

! DNS is a stream channel density in 1/KM
  real, PARAMETER :: DNS = 2.5
! __________________________________________________________


!  print *, "Soiltext: ", LDT_rc%soiltext_gridDesc(n,1:10)
!  print *, soiltext%vlevels
!  print *, soiltext%value(1000,1000,1), soiltext%value(1000,1000,10)

!  print *, "HSG: ", LDT_rc%hsg_gridDesc(n,1:10)
!  print *, hsg%vlevels
!  print *, hsg%value(1000,1000,1), hsg%value(1000,1000,2)

!- Establish number of rows and columns for area to estimate parameters:
!-  Use soil texture field as default (and assume HSG matches this):
 
  select case( LDT_rc%soiltext_proj )  ! Base selection on parameter projection 
 
    case( "latlon" )
      num_cols = LDT_rc%soiltext_gridDesc(n,2)
      num_rows = LDT_rc%soiltext_gridDesc(n,3)

    case default
      write(LDT_logunit,*)"[ERR] SAC-HTET: No other soil parameter grid projections "
      write(LDT_logunit,*) " are supported at this time, except:  'latlon', "
      write(LDT_logunit,*) " or read in parameters (which can support HRAP grids)."
      call LDT_endrun
  end select

!- Cosby Soils Table step:
   inquire(file=trim(cosbysoils_table), exist=file_exists)
   if(.not.file_exists) then 
      write(LDT_logunit,*) "[ERR] SAC-HTET: No Cosby soils table file not found."
      write(LDT_logunit,*) " This table is required to calculate SAC soil parameters"
      write(LDT_logunit,*) " using the 'Koren_v1' method."
      write(LDT_logunit,*) "Program stopping ..."
      call LDT_endrun
   endif

   ftn = LDT_getNextUnitNumber()
   open( ftn, file=cosbysoils_table, form="formatted" )
   write(LDT_logunit,*) " -- Opening Cosby soils table ... ",trim(cosbysoils_table)

    read(ftn,*)
    read(ftn,*)
    do i = 1, maxclass
       read(ftn,*,END=999) icls, SWMAX(i), SWFLD(i), SWWLT(i),&
                           KS(i), sinit(i), hhed(i)
!       print *, icls, SWMAX(i), SWFLD(i), SWWLT(i),&
!                           KS(i), sinit(i), hhed(i)
    end do
999 if( icls .ne. maxclass ) THEN
        WRITE(*,*) "[ERR] NUMBER DEFINED CLASSES ARE LESS MAXCLASS "
        WRITE(*,*) " MAXCLASS = ",maxclass, ", ICLASS = ",icls
        call LDT_endrun
    endif
   call LDT_releaseUnitNumber(ftn)

!- Initialize SAC-HTET soil parameter fields:
   uzfwm%value = -1.  ! 0.
   uztwm%value = -1.  ! 0. 
   lztwm%value = -1.  ! 0. 
   lzfsm%value = -1.  ! 0. 
   lzfpm%value = -1.  ! 0.
   uzk%value   = -1.  ! 0. 
   lzsk%value  = -1.  ! 0.
   lzpk%value  = -1.  ! 0.  
   zperc%value = -1.  ! 0.
   rexp%value  = -1.  ! 0.
   pfree%value = -1.  ! 0. 
   frz_stxt%value = -1. ! 0. 

   txtlw%value = -1   ! 0.
   txtot%value = -10  ! 0.
   dzup1%value = -1.  ! 0.
   dzlw%value  = -1.  ! 0.

!- NEW STATEMENTS TO GET DOMINANT SOIL TEXTURE
   DO J=1,15
!      ITXTUP(J)=0
!      ITXTLW(J)=0
!      ITXTOT(J)=0
      ITXTUP(J)=-1
      ITXTLW(J)=-1
      ITXTOT(J)=-10
   ENDDO

   NSUM = 0
   IMPRB = 0 

!- Loop over each input col and row point in input grid(s):
!   print *, "rows, cols:", num_rows, num_cols
   do r = 1, num_rows
      do c = 1, num_cols

     !- Initialize the total losses (SINX) in mm and the sum of the 
     !  percentages of the HSGs:
        SINX=0.
        shsgx=0.0
        do I=1,5
           if( hsg%value(c,r,I) .ge. 0.0 ) then
!             SINX=SINX+SINHSG(I)*ihsgx(I)*0.01    ! Original
             SINX=SINX+SINHSG(I)*hsg%value(c,r,I)*0.01
!             shsgx=shsgx+ihsgx(i)
             shsgx=shsgx+hsg%value(c,r,I)
           endif
        end do

      ! There is a section of decision rules here in the annotated code
        if(shsgx .lt. 99.9) then
!          numbhsg=numbhsg+1
          if(shsgx .lt. 40.0) then
            sinx=0.0
            do i=1,5
!              sinx=sinx+0.25*sinhsg(i)
              sinx=sinx+(0.25*SINHSG(I))
            enddo
          else
            sinx=sinx/shsgx
          endif
        endif

     !- CALCULATE UPPER ZONE THICKNESS (ZUP):
     !  Assumptions based on Koren et al. (2000)
        ZMXD = 1750.   ! Middle of 10th STATSGO layer (in mm)
        UZFM = SINX
        SX = SINX
        IF(UZFM .EQ. -99.) GOTO 77
        ZUPX = 0.
        IPRZ = 0
        DWUP = 0.
        N16 = 0
!        DO ILR=1,NBAND-1
        DO ILR = 1, (soiltext%vlevels-1)
!           IVALX=itxtx(ILR)
           IVALX = int(soiltext%value(c,r,ILR))
         ! Assign any soil textures>16 to 10 (sandy clay) if soil_layer>1 (?):
           if(ivalx .gt. 16) then  
              if(ilr .gt. 1) then
                ivalx = 10
              else
                goto 77
              endif
           endif

         ! Handle missing or organic (texture=13) soil types: 
           IF(IVALX .LE. 0 .OR. IVALX .EQ. 13) THEN
           ! Recall texture for beneath the surface (1st) layer and
           !  set to previous layer texture:
             IF(ILR .GT. 1) THEN
               IVALX=IVAL_PRV
             ELSE
               ZUPX=-99.
               GOTO 77
             ENDIF

           ELSE
           ! If "other" soil type reached in layer ...
             IF(IVALX .EQ. 16) THEN
             ! Top layer: set parameters to 0, go to next cell...
               IF(ILR .EQ. 1) THEN
                 ZUPX=-1
                 IPRZ=-1
                 UZFM=-1
                 GOTO 77
               ELSE
               ! If there are < 6 'other' texture layers so far, 
               !  set the texture class to previous layer's texture:
                 IF(N16 .LE. 6) THEN
                   IVALX=IVAL_PRV
                   N16=N16+1
                 ELSE
                 ! If > 6 layers, declare gridcell to "missing" and parms = 0:
                   ZUPX=-1
                   IPRZ=-1
                   GOTO 77
                 ENDIF
               ENDIF
             ENDIF
 
           ! Compute the thickness, dz, of the current soil layer, ILR:
             DZ=LAYER(ILR+1)-LAYER(ILR)
 
           ! If current layer is water (14) or bedrock (15), then
           !  increment the number of impermeable layers (IMPRB):
             IF(IVALX .EQ. 14 .OR. IVALX .EQ. 15) THEN
               IF(ILR .EQ. 1) THEN
                 IMPRB=IMPRB+1
                 GOTO 33
               ENDIF
             ! If already have water or bedrock in deeper layers, 
             !  go to section and calculate the parameters:
               IF(IPRZ .EQ. 1) GOTO 777
   
             ! IF layer is first water or bedrock encountered, continue
             !  computing depth  to bottom of soil column (ZMXD):
               ZMXD=0.5*(LAYER(ILR)+LAYER(ILR+1))
               IPRZ=1
               DZ=ZMXD-LAYER(ILR)
               IVALX=IVAL_PRV
             ! Compute incremental Upper Zone depth containing space between
             !  porosity and field capacity and used to satifsy HSG losses(SX):
               DZUP=MIN(SX/(SWMAX(IVALX)-SWFLD(IVALX)),DZ)
               DWUP=DWUP+DZUP*(SWMAX(IVALX)-SWFLD(IVALX))
               SX=0.
             ! If incremental storage < 0.5*ilyr_dz, computer upper zone variables:
               IF(DZUP .LT. DZ) THEN
                 ZUPX=ZUPX+DZUP
                 DWUP=DWUP+DZUP*(SWMAX(IVALX)-SWFLD(IVALX))
                 UZFM=DWUP
               ELSE
                 ZUPX=ZMXD
                 UZFM=UZFM-SX+(SWMAX(IVALX)-SWFLD(IVALX))*DZ
                 DWUP=DWUP+DZ*(SWMAX(IVALX)-SWFLD(IVALX))
               ENDIF

             ELSE
             ! If remaining storage is to be satisified < 0, go to next layer
               if( sx .le. 0.) goto 888

             ! If current layer is 4+> and texture = not water/bedrock but > 8,
             !  set the UZFMW and assume all storages are 0:
               IF(ILR .GT. 3 .AND. IVALX .GT. 8) THEN
                 UZFM=DWUP
                 sx=0.
                 goto 888    ! Go to next layer
               ENDIF

             ! If none of the previous texture conditions have been met, then
             !  assume current layer is "normal" case layer and compute 
             !  remaining amount of total abstractions (SX) for subsequent layers:
               DZUP=MIN(SX/(SWMAX(IVALX)-SWFLD(IVALX)),DZ)
               ZUPX=ZUPX+DZUP
               SX=SX-(SWMAX(IVALX)-SWFLD(IVALX))*DZUP
             ! Compute current estimate of DWUP (final value is set to UZFWM??):
               DWUP=DWUP+DZUP*(SWMAX(IVALX)-SWFLD(IVALX))
             ENDIF
           ENDIF

 888       IVAL_PRV=IVALX

        END DO   ! End loop over total number of soil layers 
777     continue
 
!- Depth for the Upper soil zone (ZUPX or Zup) has now been calculated.
!  Next the thickness for the lower zone (ZLWX) gets computed, and then
!   the soil layers are iterated through in upper and lower zones to
!   accumulate the soil hydraulic properties.

     !- CALCULATE AVERAGE SOIL CHARACTERISTICS OVER EACH PROFILE
        ZLWX = ZMXD-ZUPX
        IPRZ=0
        N16=0
        SMU=0.
        SFU=0.
        SWU=0.
        SML=0.
        SFL=0.
        SWL=0.
        SKL=0.
        shed=0.
        sku=0.

     !- UPPER ZONE SOIL PARAMETERS AVERAGED OVER PROFILE:
     !  Iterate over 11 layers and determine number of layers in the
     !   upper (NU) and lower (NL) zones. Note: NL is the topmost layer
     !   in the lower zone.
!        DO ILR=1,NBAND-1
        DO ILR=1,(soiltext%vlevels-1)
           NU=ILR
           NL=ILR
           IF(LAYER(ILR+1) .GT. ZUPX) GOTO 200
           IF(LAYER(ILR+1) .EQ. ZUPX) GOTO 201
        ENDDO
201     NL=ILR+1
200     CONTINUE

     !- Loop over soil layers in the upper zone: 
        DO ILR=1,NU
!           IVALX=ITXTX(ILR)
           IVALX = int(soiltext%value(c,r,ILR))
           if(ivalx .gt. 16) then  ! if 'other' texture
             if(ilr .gt. 1) then   ! reassign 'other' to sandy-clay (?)
               ivalx = 10
             else
               goto 77
             endif
           endif
         ! If texture missing or "organic", go to next layer:
           IF(IVALX .LE. 0 .OR. IVALX .EQ. 13) THEN
             GOTO 77
           ENDIF
           IF(ILR .EQ. 1) IVAL_PRV=IVALX

         ! If 'other' texture, adopt previuous layer:
           IF(IVALX .EQ. 16) THEN
              IF(ILR .EQ. 1) THEN
                ZUPX=-1
                IPRZ=-1
                UZFM=-1
                GOTO 77
              ELSE
              ! As long as 6 or fewer occurrences of 'other':
                IF(N16 .LE. 6) THEN
                   IVALX=IVAL_PRV
                   N16=N16+1
                ELSE
                   ZUPX=-1
                   IPRZ=-1
                   GOTO 77
                ENDIF
              ENDIF
           ENDIF

         ! NEW CALCULATIONS OF UPPER LAYER DOMINANT SOIL TEXTURE CLASS
           ITXTUP(IVALX)=ITXTUP(IVALX)+1
           ITXTOT(IVALX)=ITXTOT(IVALX)+1
         !     NTXTUP(J)=NTXTUP(J)+1

         ! Compute the upper layer soil moisture and sat. hyd. cond. properties:
           DZU=LAYER(ILR+1)-LAYER(ILR)
           IF(IVALX .EQ. 15) IVALX=IVAL_PRV
           IF(ILR .EQ. NU) DZU=ZUPX-LAYER(ILR)
           IVAL_PRV=IVALX
           SMU=SMU+SWMAX(IVALX)*DZU/ZUPX
           SFU=SFU+SWFLD(IVALX)*DZU/ZUPX
           SWU=SWU+SWWLT(IVALX)*DZU/ZUPX
           sku=sku+ks(ivalx)*dzu/zupx

      ! End upper layers loop:
        ENDDO

    !- Examine STATSGO layers in the lower zone -!

      ! Lower zone soil parameters averaged over profile:

      ! If upper zone depth = max soil depth, then NO lower zone ...
        IF(ZMXD .EQ. ZUPX) GOTO 55
      ! (Then no lower zone parametes ... issue??
!        DO ILR=NL,NBAND-1
        DO ILR=NL,(soiltext%vlevels-1)
           IF(LAYER(ILR+1) .GT. ZMXD) GOTO 202
        ENDDO
202     NE=ILR
        N16=0
        IVAL_PRV=IVALX

      ! Loop over bottom or lower layers:
        DO ILR=NL,NE
!           IVALX=ITXTX(ILR)
           IVALX=int(soiltext%value(c,r,ILR))
         ! If texture>16 in top layer, go to next soil layer:
           if(ivalx .gt. 16) then
             if(ilr .gt. 1) then
               ivalx = 10
             else
               goto 77
             endif
           endif

         ! If texture is missing or organic (13) and layer<9, go to next layer:
           IF(IVALX .LE. 0 .OR. IVALX .EQ. 13) THEN
             IF(ILR .LT. 9) GOTO 77
             IVALX=IVAL_PRV
           ENDIF
           DZL=LAYER(ILR+1)-LAYER(ILR)
         ! If layer is bedrock, assign texture from previous layer:
           IF(IVALX .EQ. 15) IVALX=IVAL_PRV

         ! If texture is 'other' for 2nd or deeper layer and <6 layers of 'other',
         !  assign texture class from previous layer:
           IF(IVALX .EQ. 16 .AND. ILR .NE. 1) THEN
             IF(N16 .LE. 6) THEN
                IVALX=IVAL_PRV
                ivalx = 8
                N16=N16+1
             ELSE
                ZLWX=-99
                IPRZ=-1
                GOTO 77
             ENDIF
           ENDIF

        !- NEW CALCULATIONS OF LOWER ZONE DOMINANT SOIL TEXTURE CLASS
           IF(ILR .NE. NL .OR. NL .NE. NU) THEN
             ITXTLW(IVALX)=ITXTLW(IVALX)+1
             ITXTOT(IVALX)=ITXTOT(IVALX)+1
        !!     NTXTLW(J)=NTXTLW(J)+1
           ENDIF
           IF(NL .EQ. NE) THEN
             DZL=ZMXD-ZUPX
           ELSE
             IF(ILR .EQ. NL) DZL=LAYER(ILR+1)-ZUPX
             IF(ILR .EQ. NE) DZL=ZMXD-LAYER(ILR)
           ENDIF
         ! Compute the depth-weighted average of the lower zone porosity (SML),
         !  Field Capacity (SFL), Wilting pt (SWL) and sat. hydraulic cond. (SKL).
           SML=SML+SWMAX(IVALX)*DZL/ZLWX    ! DZL - thickness of layer
           SFL=SFL+SWFLD(IVALX)*DZL/ZLWX
           SWL=SWL+SWWLT(IVALX)*DZL/ZLWX
           SKL=SKL+KS(IVALX)*DZL/ZLWX
         ! Compute the capillary water potential shed using layers in lower zone.
           shed=shed+hhed(ivalx)*DZL/ZLWX  
           IVAL_PRV=IVALX

      ! End lower zone layers loop
        ENDDO

 55  NSUM=NSUM+1

! - At this point, we have computed the depths/thicknesses of the upper and
!    lower zones and the depth-weighted averages of the soil hydraulic
!    properties in the upper and lower zones.
! - Next compute the SAC-SMA parameters for each gridcell.
!    See Koren et al.(2000,2003) and Dry Region Parameters report for new ZPERC.

!      SAC_PAR(1)=UZFM
      uzfwm%value(c,r,1) = UZFM
!      SAC_PAR(2)=(SFU-SWU)*ZUPX
      uztwm%value(c,r,1) = (SFU-SWU)*ZUPX
!      SAC_PAR(3)=(SFL-SWL)*ZLWX
      lztwm%value(c,r,1) = (SFL-SWL)*ZLWX
!      SAC_PAR(6)=1.-(SFU/SMU)**1.6
      uzk%value(c,r,1) = 1.-(SFU/SMU)**1.6

      ZSK=(1.-(SFU/SMU)**1.6)/(1.+2.*(1.-SWL))

      IF(ZLWX .GT. 0.) THEN
        ZFSM=(SML-SFL)*ZLWX*(SWL/SML)**1.6
!        SAC_PAR(4)=SAC_PAR(4)+ZFSM
        lzfsm%value(c,r,1) = lzfsm%value(c,r,1)+ZFSM
        ZFLM=(SML-SFL)*ZLWX
        ZFPM=ZFLM-ZFSM
!        SAC_PAR(5)=ZFPM
        lzfpm%value(c,r,1) = ZFPM
        ZPK=1.-EXP(-0.00406*0.001*SKL*ZLWX*0.001*DNS**2/&
                  (3.5*(SML-SFL)**1.66))
      ! Note: ZPK = LZPK (Refer to Eq. A8.6 of Rev. Koren et al.,2003)
        PBASE=ZFSM*ZSK+ZFPM*ZPK

      ! 3/08  Added new option for ZPERC (capillary preasure):
        if(izperc .eq. 'old' .or. shed .gt. ZLWX) then
!          SAC_PAR(9)=((SFL-SWL)*ZLWX+ZFSM+ZFPM-&
!                       PBASE)/PBASE
          zperc%value(c,r,1) = ((SFL-SWL)*ZLWX+ZFSM+ZFPM-&
                                 PBASE)/PBASE
        else
     !     if(izperc .eq. 'new') then
        ! vk 01/09 select minimum hydraulic conductivity:
          skmin=skl
          if(sku .lt. skl) skmin=(sku+skl)*0.5
          pbase_hr=pbase/24.
          pcap=(shed+zupx*(smu-swu))*(sml-swl)
        ! Note: pcap is part of Eqt. 13 in Koren (2009)

        ! grid now in newzperc    tcurve=(0.1*skl*pcap/24.)**0.33333
        ! vk 01/09     tcurve=(0.1*skl*pcap)**0.33333
        ! vk 01/09     zperc=(skl*tcurve+sqrt(2*pcap*skl*tcurve)-&
        ! vk 01/09            pbase_hr*tcurve)/(pbase_hr*tcurve)

        ! From Equation 13 in Koren (2009)
          tcurve=(0.1*skmin*pcap)**0.33333
!          zperc=(skmin*tcurve+sqrt(2*pcap*skmin*tcurve)-&
          zperc1=(skmin*tcurve+sqrt(2*pcap*skmin*tcurve)-&
                  pbase_hr*tcurve)/(pbase_hr*tcurve)

        ! vk 01/2009 added limit on zperc maximum
!          if( zperc .gt. 2000.) then
          if( zperc1 .gt. 2000.) then
!             write(*,*) 'ZPERC above limit:',icol,irow,zperc1
             write(*,*) 'ZPERC above limit:',c,r,zperc1
!             zperc=2000.
             zperc1=2000.
          endif

!          SAC_PAR(9)=SAC_PAR(9)+zperc
          zperc%value(c,r,1) = zperc%value(c,r,1) + zperc1

!          if(zperc .le. 0.) then
          if(zperc1 .le. 0.) then
!             write(*,*) 'newZPERC=',zperc,shed,skl,zupx,pmax
             write(*,*) 'newZPERC=',zperc1,shed,skl,zupx !,pmax
             write(*,*) '         ',smu,swu,sml,swl
             write(*,*) '        =',zfsm,zsk,zfpm,zpk,pbase
             stop
          endif
        endif

        SFREE=(SWL/SML)**1.6

      ELSE
        ZPK=0.005
        SFREE=0.15
      ENDIF

!      SAC_PAR(7)=ZSK
      lzsk%value(c,r,1) = ZSK
!      SAC_PAR(8)=ZPK
      lzpk%value(c,r,1) = ZPK
!      SAC_PAR(10)=(SWL/(SWWLT(1)-0.001))**0.5
      rexp%value(c,r,1) = (SWL/(SWWLT(1)-0.001))**0.5
!      SAC_PAR(11)=SFREE
      pfree%value(c,r,1) = SFREE
!      SAC_PAR(15)=ZUPX
      dzup1%value(c,r,1) = ZUPX
!      SAC_PAR(16)=ZLWX
      dzlw%value(c,r,1) = ZLWX

!       write(*,*) '  ',uzfm,(SFU-SWU)*ZUPX,(SFL-SWL)*ZLWX,zfsm,zfpm, &
!        1.-(SFU/SMU)**1.6,zsk,zpk,((SFL-SWL)*ZLWX+ZFSM+ZFPM-PBASE)/PBASE, &
!           (SWL/(SWWLT(1)-0.001))**0.5,SFREE,zupx,zlwx

      GOTO 33

! 77  MISS=MISS+1
 77   continue

 33   continue

  7   continue

#if 0
   !-  NEW LOOP TO GET DOMINANT SOIL CLASS
      DTXTUP=0
      DTXTLW=0
      DTXTOT=0
      DO I=1,15
         IF(ITXTUP(I) .GT. DTXTUP) THEN
           DTXTUP=ITXTUP(I)
           ITXTUP(1)=I
         ENDIF
         IF(ITXTLW(I) .GT. DTXTLW) THEN
           DTXTLW=ITXTLW(I)
           ITXTLW(1)=I
         ENDIF
         IF(ITXTOT(I) .GT. DTXTOT) THEN
           DTXTOT=ITXTOT(I)
           ITXTOT(1)=I
         ENDIF
      ENDDO

   !- Fill zone dominant texture grids for upper, lower, and total zones
!      SAC_PAR(12)=ITXTUP(1)
      frz_stxt%value(c,r,1) = ITXTUP(1)
!      SAC_PAR(13)=ITXTLW(1)
      txtlw%value(c,r,1) = ITXTLW(1)
!      SAC_PAR(14)=ITXTOT(1)
      txtot%value(c,r,1) = ITXTOT(1)
#endif

      end do   ! End total column loop
   end do      ! End total row loop

   write(LDT_logunit,*) " ... Done ... Creating SAC-HTET Soil Parameters ... "

end subroutine create_SACsoilparms_Korenv1
