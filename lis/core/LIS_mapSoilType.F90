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
! !ROUTINE: LIS_mapSoilType
!  \label{LIS_mapSoilType}
!
! !DESCRIPTION:
!  3 Jul 2001: Matt Rodell Initial specifiation
! 13 Dec 2004: Sujay Kumar Updated with support for STATSGO classification
! 16 Jan 2014: Greg Fall   Adjusted STATSGO classification logic to remove
!                          discontinuities and fix missing clay/silt loam
!                          regions of the USDA soil triangle. The logic is
!                          stated explicitly in the USDA training module
!                          "Soil Mechanics Level I, Module 3: USDA Textural
!                          Soil Classification," pp. 15-16, retrieved from
!                          ftp://ftp.wcc.nrcs.usda.gov/wntsc/H\&H/training/soilsOther/soil-USDA-textural-class.pdf
!
! !INTERFACE:
subroutine LIS_mapSoilType(class,sand,clay,silt,soiltyp)
! !USES:

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: class
  real, intent(in)       :: sand, clay, silt
  integer, intent(inout) :: soiltyp
!
! !DESCRIPTION:
!  This subroutine uses the percentages of sand, silt, and clay   
!  to convert to soil texture data. The transformation is based on the
!  type of classification used. This routine supports the transformation
!  to a Zobler (9 classes) or the STATSGO (19 classes) scheme. 
!
!  The arguments are: 
!  \begin{description}
!   \item[class]
!    soil classification scheme (1-zobler, 2-statsgo)
!   \item[sand]
!    array containing the sand fraction data
!   \item[clay]
!    array containing the clay fraction data
!   \item[silt]
!    array containing the silt fraction data
!   \item[soiltyp]
!    array containing the derived soil texture 
!   \end{description}
!EOP
  real :: sa,cl,si

  if(class .eq. 1) then !zobler
     if (clay .lt. 0.00) then
        soiltyp = -9999
        return
     else
        cl = clay
        sa = sand
     endif
     if (cl .lt. 0.23) then
        if (sa .lt. 0.50) then
           soiltyp = 8          ! loam
        else
           if (sa .lt. 0.75) then
              soiltyp = 4        ! sandy loam
           else
              soiltyp = 1        ! loamy sand
           end if
        end if
     else 
        if (cl .lt. 0.28) then
           if (sa .lt. 0.45) then
              soiltyp = 8        ! loam
           else
              soiltyp = 7        ! sandy clay loam
           endif
        else
           if (cl .lt. 0.37) then
              if (sa .lt. 0.2) then
                 soiltyp = 2      ! silty clay loam
              else
                 if (sa .lt. 0.43) then
                    soiltyp = 6    ! clay loam
                 else
                    soiltyp = 7    ! sandy clay loam
                 end if
              end if
           else
              if (cl .lt. 0.41) then
                 if (sa .lt. 0.2) then
                    soiltyp = 2   ! silty clay loam
                 else
                    if (sa .lt. 0.43) then
                       soiltyp = 6    ! clay loam
                    else
                       soiltyp = 5    ! sandy clay
                    end if
                 end if
              else
                 if (sa .lt. 0.43) then
                    soiltyp = 3      ! light clay
                 else
                    soiltyp = 5      ! sandy clay
                 end if
              end if
           end if
        end if
     end if
  elseif(class.eq.2) then !STATSGO texture
     if (clay .lt. 0.00) then
        soiltyp = -9999
        return
     else
        cl = clay
        sa = sand
        si = silt
     endif
     if((sa.ge.0.85).and.((cl*1.5+si).le.0.15)) then
        soiltyp = 1 !sand
     elseif((cl*2+si).le.0.30) then 
        soiltyp = 2 !loamy sand
     !Logan change begin
     elseif((cl.lt.0.07).and.&
            (si.lt.0.50).and.&
            (sa.gt.0.43).and.&
            (sa.lt.0.52)) then
        soiltyp = 3 !sandy loam
     elseif(((cl*2+si).gt.0.30).and.&
          (cl.le.0.20).and.&
          (sa.ge.0.52)) then
        soiltyp = 3 !sandy loam
     !elseif(((cl*2+si).gt.0.30).and.&
     !     (cl.lt.0.07).and.(si.lt.0.50)) then
     !   soiltyp = 3 !sandy loam
     !elseif(((cl*2+si).gt.0.30).and.&
     !     (cl.gt.0.07).and.(cl.le.0.20).and.&
     !     (sa.ge.0.52)) then
     !   soiltyp = 3 !sandy loam
     !Logan change end
     elseif((cl.ge.0.07).and.&
          (cl.le.0.27).and.&
          (si.ge.0.28).and.&
          (si.lt.0.50).and.&
          (sa.lt.0.52)) then
        soiltyp = 6 !loam
     elseif((si.ge.0.50).and.&
          (cl.ge.0.12).and.&
          (cl.le.0.27)) then
        soiltyp = 4 !silt loam
     !Logan add begin
     elseif((si.ge.0.50).and.&
            (si.lt.0.80).and.&
            (cl.lt.0.12)) then
        soiltyp = 4 !silt loam
     !Logan add end
     elseif((si.ge.0.80).and.&
          (cl.lt.0.12)) then 
        soiltyp = 5 !silt
     elseif((cl.gt.0.20).and.&
          (cl.lt.0.35).and.&
          (si.lt.0.28).and.&
          (sa.ge.0.45)) then
        soiltyp = 7 !sandy clay loam
     elseif((cl.gt.0.27).and.&
          (cl.lt.0.40).and.&
          (sa.gt.0.20).and.&
          (sa.lt.0.45)) then 
        soiltyp = 9 !clay loam
     elseif((cl.gt.0.27).and.&
          (cl.lt.0.40).and.&
          (sa.le.0.20)) then 
        soiltyp = 8 !silty clay loam
     elseif((cl.ge.0.35).and.&
          (sa.ge.0.45)) then
        soiltyp = 10 !sandy clay
     elseif((cl.ge.0.40).and.&
          !Logan change begin
          !(cl.le.0.60).and.&
          !Logan change end
          (sa.lt.0.45).and.&
          (si.lt.0.40)) then
        soiltyp = 12 !clay
     elseif((cl.ge.0.40).and.&
          (si.ge.0.40)) then 
        soiltyp = 11 !silty clay
!    elseif(gindex.eq.-1) then 
!       soiltyp = 14 !water
     else
        soiltyp = 15 !bedrock
     endif
  endif
end subroutine LIS_mapSoilType

