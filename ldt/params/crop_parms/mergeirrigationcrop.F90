!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA GSFC Land Data Toolkit (LDT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
!
! !ROUTINE: mergeirrigationcrop
!  \label{mergeirrigationcrop}
!
! !REVISION HISTORY:
!  20 Sep 2019: HK Beaudoing; Initial implementation for merging crop type
!                             and irrigation fraction information to align.
!                             The code is based on preprocess_irrig.F90 
!                             written by Sarith Mahanama.
!
! !INTERFACE:
subroutine mergeirrigationcrop( crop_classification, nest, num_types, &
           fgrd, mergesource, mergefile, merge_gridtransform, cropcaflag)

! !USES:
  use LDT_coreMod,  only : LDT_rc
  use LDT_logMod,   only : LDT_logunit, LDT_getNextUnitNumber, &
          LDT_releaseUnitNumber, LDT_verify, LDT_endrun
  use LDT_paramDataMod
  use LDT_LSMCropModifier_Mod
  use LDT_irrigationMod

  implicit none

! !ARGUMENTS: 
  character(len=*), intent(in)  :: crop_classification
  integer,          intent(in)  :: nest
  integer,          intent(in)  :: num_types
  real,          intent(inout)  :: fgrd(LDT_rc%lnc(nest),LDT_rc%lnr(nest),num_types)
  character(len=*), intent(in)  :: mergesource
  character(len=*), intent(in)  :: mergefile
  character(len=*), intent(in)  :: merge_gridtransform
  logical, intent(in)           :: cropcaflag
!
! !DESCRIPTION:
!  This subroutine returns crop fractions where irrigated crop fractions
!  are aligned with the choice of irrigation fraction data.
!  Additional source of irrigation fraction can be read in and merged to
!  the crop and irrigation datasets. 
!  Originally designed for merging GIA, GRIPC, and MIRCA datasets.
!  In this case, 
!    1. GIA provides the primary irrigation fraction (intensity). 
!    2. GRIPC provides fraction of paddy irrigation, irrigation and rainfed.
!    3. MIRCA provides irrigated crop and rainfed fractions.
!  The grid fraction returned from this routine for "MIRCA" and "MIRCA52" 
!  crop type source options are as follows:
!    1. num_types = 26, fgrd is a sum of irrigated and rainfed crops, or
!       rainfed crops where no irrigated crops.
!    2. num_types = 52, fgrd(1:26)= irrigated and fgrd(27:52)=rainfed
!
!  If cropcalender option is true, this routine adjusts the planting/hearvesting
!  dates where we plant wheat due to mismatch b/w MIRCA and GIA/GRIPC.
!  
!  The arguments are:
!  \begin{description}
!   \item[crop_classification]
!     Crop classification
!   \item[nest]
!     number of nesting
!   \item[num_types]
!     number of crop types
!   \item[fgrd]
!     the crop type fractions
!   \item[mergesource]
!     additional source of irrigation data to be merged
!   \item[mergefile]
!     irrigation data filename to be merged
!   \item[merge_gridtransform]
!     grid transformation for the merging data
!   \item[cropcaflag]
!     flag for computing crop calendar
!   \end{description}
!EOP      
!- Local:
   integer, parameter :: NCROPS = 26   !MIRCA
   integer            :: I,J,N
   real               :: GIA          ! GIA irrigfrac in fraction
   real               :: RFRAC(LDT_rc%lnc(nest),LDT_rc%lnr(nest)) ! GRIPC
   real               :: IFRAC(LDT_rc%lnc(nest),LDT_rc%lnr(nest)) ! GRIPC
   real               :: PFRAC(LDT_rc%lnc(nest),LDT_rc%lnr(nest)) ! GRIPC
   real               :: NONCROP(LDT_rc%lnc(nest),LDT_rc%lnr(nest)) ! GRIPC
   real               :: GRIPCTOT     ! GRIPC total fraction
   real               :: MICROP, MRCROP     ! MIRCA irr, rainfed crop
   real               :: MICROPA, MRCROPA   ! sum of non rice crops
   real               :: MIRICEA, MRRICEA   ! rice crop fraction
   integer            :: mc     ! multiple cropping seasons
   real, allocatable, dimension (:,:,:,:)   :: plantday,harvestday
   real, allocatable, dimension (:,:) :: frac_gridp, frac_gridh
   integer            :: imn,imx,jmn,jmx,l,c,r
   logical            :: foundPt
   
! _____________________________________

  write(LDT_logunit,*)"[INFO] Merging the crop type and irrigation for given classification: ",trim(crop_classification)
                      
! Read in additional merging dataset if selected
  if ( trim(mergesource) .ne. "none" ) then
    write(LDT_logunit,*)"[INFO] Merging additional irrigation dataset: ", & 
                        trim(mergesource)
    if ( trim(mergesource) == "GRIPC" )then
     LDT_irrig_struc(nest)%cropwatsrc%num_bins = 4
     allocate(LDT_irrig_struc(nest)%cropwatsrc%value(&
                    LDT_rc%lnc(nest),LDT_rc%lnr(nest),&
                    LDT_irrig_struc(nest)%cropwatsrc%num_bins))
     ! overwrite irrig_struc variables
     LDT_irrig_struc(nest)%irrigtypefile = trim(mergefile)
     LDT_irrig_struc(nest)%irrigtype_gridtransform = trim(merge_gridtransform)
     write(LDT_logunit,*)"[INFO] Overwritten irrigtypefile: ",trim(LDT_irrig_struc(nest)%irrigtypefile)  
     call read_GRIPC_irrigtype(nest, LDT_irrig_struc(nest)%cropwatsrc%value, &
          LDT_irrig_struc(nest)%cropwatsrc%num_bins )
     RFRAC = LDT_irrig_struc(nest)%cropwatsrc%value(:,:,1) !Rain-fed croplands
     IFRAC = LDT_irrig_struc(nest)%cropwatsrc%value(:,:,2) !Irrigated
     PFRAC = LDT_irrig_struc(nest)%cropwatsrc%value(:,:,3) !Paddy
     NONCROP = LDT_irrig_struc(nest)%cropwatsrc%value(:,:,4) !noncrop
    else
     write(LDT_logunit,*) "[ERR]  Merging additional irrigation for this "
     write(LDT_logunit,*) " dataset is not yet supported ",trim(mergesource)
     call LDT_endrun
    endif
  endif    ! mergesource
! Prepare planting/harvesting dates local arrays if selected
  if ( cropcaflag ) then
    allocate(plantday(LDT_rc%lnc(nest),LDT_rc%lnr(nest),num_types,LDT_LSMCrop_struc(nest)%multicroppingmax))
    allocate(harvestday(LDT_rc%lnc(nest),LDT_rc%lnr(nest),num_types,LDT_LSMCrop_struc(nest)%multicroppingmax))

    do mc = 1, LDT_LSMCrop_struc(nest)%multicroppingmax
     do N = 1, num_types
        plantday(:,:,N,mc) = LDT_LSMCrop_struc(nest)%plantday%value4d(:,:,N,mc) 
        harvestday(:,:,N,mc) = LDT_LSMCrop_struc(nest)%harvestday%value4d(:,:,N,mc) 
     enddo  !  N
    enddo  !mc
  endif    ! cropcaflag
                       
! ----------------------------------------------------------------------
! irrigcrop(I,J,CROPS) : irrigated crop fractions with rice is the 3rd slice
! rainfedcrop(I,J,CROPS) : rainfed crop fractions with rice is the 3rd slice
! 1. GIA-GRIPC 
!    Merge LDT_irrig_struc(n)%irrigfrac (irrigation fraction derived in 
!    LDT_irrigationMod, i.e. GIA) with the second irrigation dataset (i.e.GRIPC)
!    to separate into IRRIGFRAC and PADDYFRAC 
!     GIA yes -> scale by GRIPC IFRAC and PFRAC 
!         no  -> IFRAC = PFRAC = 0
! 2. GIA-GRIPC-MIRCA irrigated non-rice(a) and rice(b) CROPS
!    Adjust irrigation crop fractions of MIRCA to match
!     IFRAC = sum of irrigated crops except for rice
!     PFRAC = irrigated rice, scaled by GIA
! 3. GIA-GRIPC-MIRCA Rainfed CROPS
!    Adjust rainfed crop fractions of MIRCA to match
!     RFRAC = sum of all rainfed crops including rice
! 4. Align and fill in CROP CALENDAR dates for adjusted crop fractions accordingly 
! ----------------------------------------------------------------------

    DO J = 1, LDT_rc%lnr(nest)
       DO I = 1, LDT_rc%lnc(nest)
!          IF(LDT_LSMparam_struc(nest)%landmask%value(I,J,1) > 0.) THEN
             IF((maxval(LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,:)) > 0.).OR. &
                (maxval(LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,:)) > 0.).OR. &
                (IFRAC (I,J) > 0.).OR.(PFRAC (I,J) > 0.).OR.(RFRAC (I,J) > 0.)) THEN

                MICROP = 0.
                MRCROP = 0.
                MIRICEA= 0.
                MRRICEA= 0.
                MICROPA= 0.
                MRCROPA= 0.
                GIA = LDT_irrig_struc(nest)%irrigfrac%value(I,J,1) * 0.01

                DO N = 1, NCROPS
                   IF (N == 3) THEN
                      IF(MIRICEA < LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N)) &
                       MIRICEA = LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N)
                      IF(MRRICEA < LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N))  &
                       MRRICEA = LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N)
                   ELSE
                      MICROP = MICROP + LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N)
                      MRCROP = MRCROP + LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N)
                   ENDIF
                END DO

                IF(MICROPA < MICROP) MICROPA = MICROP !sum of icrops no rice
                IF(MRCROPA < MRCROP) MRCROPA = MRCROP !sum of rcrops no rice

                ! 1. GIA-GRIPC
                ! .........

                IF(GIA <= 0 ) THEN
                   ! MASK OUT non-irrigated per GIA
!HKB reassign them to RFRAC values    FRAC (I,J) = 0 
                   RFRAC (I,J) = RFRAC (I,J) + IFRAC (I,J) + PFRAC (I,J)
                   IFRAC (I,J) = 0.
                   PFRAC (I,J) = 0.
                ELSE
                   IF ((IFRAC (I,J)  + PFRAC (I,J)) <= 0.) THEN
                      ! GRIPC does not have data
                      PFRAC (I,J)  = 0.
                      IFRAC (I,J)  = GIA
                   ELSE
                      ! GRIPC has data
                      MICROP = PFRAC (I,J) + IFRAC (I,J) 
                         PFRAC (I,J)  = PFRAC (I,J) * GIA / MICROP
                         IFRAC (I,J)  = IFRAC (I,J) * GIA / MICROP
                   ENDIF
                ENDIF
!HKB            reset RFRAC and NONCROP in GRIPC according to new IFRAC/PFRAC
                GRIPCTOT = RFRAC(I,J)+PFRAC(I,J)+IFRAC(I,J)+NONCROP(I,J)
                if ( GRIPCTOT .GT. 0 ) then
                  RFRAC (I,J) = RFRAC(I,J) / GRIPCTOT
                  NONCROP (I,J) = NONCROP(I,J) / GRIPCTOT
                endif

                ! 2a. GIA-GRIPC-MIRCA IRRIGATED CROPS
                ! ...............................

                IF (IFRAC(I,J)  == 0) THEN
                   DO N = 1, NCROPS
                      IF(N /= 3) &
                       LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N) = 0.
                   END DO
                ELSE

                   IF(MICROPA > 0.) THEN

                      ! MIRCA irrigation crop area is present, thus scale 
                      ! crop fractions to match GIA-GRIPC (i.e. GIA)
                      DO N = 1, NCROPS
                         IF(N /= 3) &
                          LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N) = &
                           LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N) * &
                           IFRAC (I,J) /  MICROPA
                      END DO

                   ELSE

                      ! MIRCA irrigation crop area is not present, but GRIPC 
                      ! irrigation area is present
                      IF(MRCROPA > 0.) THEN

                         ! MIRCA rainfed area is present
                         DO N = 1, NCROPS
                            IF(N /= 3) THEN
                              LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N)= &
                               LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N) * IFRAC (I,J) /  MRCROPA
                               LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,N) = 0.
                            ENDIF
                         END DO
                         MRCROPA = 0.
                      ELSE

                         ! MIRCA irrigated and rainfed do not have data, thus
                         ! plant some wheat
                         LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,1) =  IFRAC (I,J)
                         ! Reset Planting/Harvesting dates to be filled later
                         if ( cropcaflag .and.  &
                            LDT_LSMCrop_struc(nest)%multicroppingmax.eq.2 ) then
                            plantday(I,J,1,1) = 999
                            plantday(I,J,1,2) = 0
                            harvestday(I,J,1,1) = 999
                            harvestday(I,J,1,2) = 0
                         endif     ! cropcaflag
                      ENDIF
                   ENDIF
                ENDIF


                ! 2b. GIA-GRIPC-MIRCA PADDY
                ! .....................
                IF(PFRAC (I,J) == 0.) THEN
                   LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,3) = 0.
                ELSE
                   IF(MIRICEA > 0.) THEN

                      ! MIRCA irrigated rice is present too, thus scale 
                      ! crop fractions to match GRIPC
                      LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,3) = &
                        LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,3) * &
                        PFRAC (I,J) /  MIRICEA

                   ELSE

                      ! MIRCA irrigated rice area is not present, but GRIPC has
                      ! paddy area present
                      IF(MRRICEA > 0.) THEN

                         ! Looks like MIRCA rainfed rice
                         LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,3) =  &
                          LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,3) * &
                          PFRAC (I,J) /  MRRICEA
                         LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,3) = 0.
                         MRRICEA = 0.
                      ELSE

                         ! MIRCA irrigated and rainfed do not have rice, thus
                         ! plant rice for PFRAC
                         LDT_LSMCrop_struc(nest)%irrigcrop%value (I,J,3) =  PFRAC (I,J)
                         ! Reset Planting/Harvesting dates to be filled later
                         if ( cropcaflag .and.  &
                            LDT_LSMCrop_struc(nest)%multicroppingmax.eq.2 ) then
                             plantday(I,J,3,1) = 999
                             plantday(I,J,3,2) = 0
                             harvestday(I,J,3,1) = 999
                             harvestday(I,J,3,2) = 0
                         endif   ! cropcaflag
                      ENDIF
                   ENDIF
                ENDIF

                ! 3. GIA-GRIPC-MIRCA Rainfed CROPS
                ! .............................

                IF( RFRAC (I,J) == 0.) THEN
                   LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,:) = 0.
                ELSE

                   IF(MRCROPA + MRRICEA > 0.) THEN

                      ! MIRCA rainfed crop area is present too, thus scale crop
                      ! fractions to match GRIPC to be used for aligning land cover crop types
                      DO N = 1, NCROPS
                        IF(N /= 3) &
                         LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,N) = &
                          LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,N) * &
                          RFRAC (I,J) / (MRCROPA + MRRICEA)
                      END DO

                      IF (MRRICEA > 0.) &
                       LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,3) = &
                        LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,3) * &
                        RFRAC (I,J) / (MRCROPA + MRRICEA)

                   ELSE

                      ! MIRCA rainfed crop area is not present but GRIPC has
                      ! rainfed area, thus plant some wheat
                      LDT_LSMCrop_struc(nest)%rainfedcrop%value (I,J,1) =  RFRAC (I,J)
                      if ( cropcaflag .and.  &
                          LDT_LSMCrop_struc(nest)%multicroppingmax.eq.2 ) then
                          plantday(I,J,1,1) = 999
                          plantday(I,J,1,2) = 0
                          harvestday(I,J,1,1) = 999
                          harvestday(I,J,1,2) = 0
                      endif   ! cropcaflag
                   ENDIF
                ENDIF

! Reset croptype fractions
! num_types=26: croptype fraction is a sum of irrigcrop + rainfedcrop, but 
!               the sum of irrigcrop maybe zero if GIA=0, then we still need
!               croptype fraction to fill cropland land cover if it exists
!               in this case, croptype is set to rainfedcrop
                DO N = 1, NCROPS
                 IF ( num_types .eq. 26 ) then  ! irrig croptypes 
                  IF (sum(LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,:)).eq.0& 
                      .and. RFRAC(I,J).gt.0 ) THEN
                  fgrd(I,J,N) = LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N)
                  ELSE
                  fgrd(I,J,N) = LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N)+&
                                LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N)
                  ENDIF

                 ELSE IF ( num_types .eq. 52 ) then  ! include rainfed 
                  fgrd(I,J,N) = LDT_LSMCrop_struc(nest)%irrigcrop%value(I,J,N)
                  fgrd(I,J,N+26) = LDT_LSMCrop_struc(nest)%rainfedcrop%value(I,J,N)
                 ELSE
                  write(LDT_logunit,*)"[ERR] Merge invalid num_types"
                  call LDT_endrun
                 ENDIF
                ENDDO
                
             ENDIF   ! GIA/GRIPC/MIRCA valid pixel
!          ENDIF   ! MASK > 0 
       END DO  ! I
    END DO  ! J

! 4. Fill in CROP CALENDAR dates for newly planted wheat and rice (N=1, 3)
! .............................
    if ( cropcaflag ) then
     allocate(frac_gridp(LDT_rc%lnc(nest),LDT_rc%lnr(nest)))
     allocate(frac_gridh(LDT_rc%lnc(nest),LDT_rc%lnr(nest)))
     DO J = 1, LDT_rc%lnr(nest)
        DO I = 1, LDT_rc%lnc(nest)
           IF(LDT_LSMparam_struc(nest)%landmask%value(I,J,1) > 0.) THEN
             DO N = 1,3,2
               frac_gridp (:, :)  = plantday(:,:,N,1) 
               frac_gridh (:, :)  = harvestday(:,:,N,1) 

               ! fill missing crop plant/harvest DOYs in irrigated and
               ! rainfed crops
               ! .......................................................
               IF(plantday(I,J,N,1) == 999) THEN
 
                IF(fgrd(I,J,N) > 0.) THEN
                 l = 1
                 foundPt = .false.
                 do while (.not. foundPt)
                    imn=MAX(I-l,1)
                    imx=MIN(I+l,LDT_rc%lnc(nest))
                    jmn=MAX(J-l,1)
                    jmx=MIN(J+l,LDT_rc%lnr(nest))
                    do r=jmn,jmx
                       do c=imn,imx
                       if(frac_gridp(c,r) < 998 .and.  &
                          frac_gridp(c,r) .gt. 0 ) then
                          plantday(I,J,N,1) = frac_gridp(c,r)
                          harvestday(I,J,N,1) = frac_gridh(c,r)
                          foundPt = .true.
                          exit
                       endif
                       end do
                    end do
                    l = l + 1
                 end do
                ELSE   ! fgrd = 0
                 plantday(I,J,N,1) = 998.
                 harvestday(I,J,N,1) = 998.
                ENDIF

               ENDIF
             END DO  !N
           ENDIF   ! MASK > 0 
        END DO  ! I
     END DO  ! J

     ! FINAL CROPCALENDAR VARIABLES  remove 998 and 999....
     do mc = 1, LDT_LSMCrop_struc(nest)%multicroppingmax
       do N = 1, LDT_rc%numcrop(nest)
         do J = 1, LDT_rc%lnr(nest)
           do I = 1, LDT_rc%lnc(nest)
              if ( plantday(I,J,N,mc) == 998. ) then
                LDT_LSMCrop_struc(nest)%plantday%value4d(I,J,N,mc) = LDT_rc%udef
              elseif ( LDT_LSMCrop_struc(nest)%plantday%value4d(I,J,N,mc) == 999. ) then
                LDT_LSMCrop_struc(nest)%plantday%value4d(I,J,N,mc) = plantday(I,J,N,mc)
              else
                LDT_LSMCrop_struc(nest)%plantday%value4d(I,J,N,mc) = plantday(I,J,N,mc)
              endif
              if ( harvestday(I,J,N,mc) == 998. ) then
                LDT_LSMCrop_struc(nest)%harvestday%value4d(I,J,N,mc) = LDT_rc%udef
              elseif ( harvestday(I,J,N,mc) == 999. ) then
                LDT_LSMCrop_struc(nest)%harvestday%value4d(I,J,N,mc) = harvestday(I,J,N,mc)
              else
                LDT_LSMCrop_struc(nest)%harvestday%value4d(I,J,N,mc) = harvestday(I,J,N,mc)
              endif
           enddo
         enddo
       enddo
     enddo
     deallocate(frac_gridp)
     deallocate(frac_gridh)
     deallocate(plantday)
     deallocate(harvestday)
    endif  ! cropcaflag

end subroutine mergeirrigationcrop
