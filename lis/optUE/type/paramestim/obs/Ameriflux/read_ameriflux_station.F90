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
! !INTERFACE:
subroutine read_ameriflux_station(n, stn_index)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc
  use LIS_logMod,  only : LIS_logunit, LIS_getNextUnitNumber, & 
       LIS_releaseUnitNumber, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use LIS_timeMgrMod, only : LIS_calendar, LIS_localtime2gmt
  use AmerifluxobsMod, only : AmerifluxObs_struc
!
! !DESCRIPTION:
! Opens the file and reads line by line, inputting each variable into its array
! Further information on reading is provided at the code
! The array will then be processed in the above subroutine
!
! !REVISION HISTORY:
! Teodor Georgiev, initial version
!
!EOP
!----------------------------------------------------------------------
  implicit none

  integer,         intent(in)      :: stn_index
  integer :: n


  real :: Month     ! Month of current data
  real :: Day       ! Day of current data
  real :: Hour      ! Hour, from 0 to 23.5
  real :: min
  real :: sec
  real :: DoY       ! Decimal day of the year
  real :: Le        ! temporary latent heat flux 
  real :: H         ! temporary sensible heat flux
  real :: Precip    ! Precipitation in mm
  real :: Rg        ! global radiation W/m2
  real :: qf_Rg     ! quality flag for Rg:
                                                ! 0: ok, 1: possible spike
  real :: qf_Rad    ! quality flag for Rg:
                                                ! 0: ok, 1: possible bad values
  real :: Rn        ! net radiation W/m2
  real :: Rd        ! diffuse radiation W/m2
  real :: Ta        ! air temperature C
  real :: Ts1       ! soil temperature depth 1 C
  real :: Ts2       ! soil temperature depth 2 C
  real :: SWC1      ! soil water content depth 1 %vol
  real :: SWC2      ! soil water content depth 2 %vol
  real :: G1	! Soil heat flux
  real :: G2	! Soil heat flux 2
  real :: Rh        ! relative humidity as a %
  real :: WS        ! Wind horizontal speed m/s
  ! Initializes other variables. iCurrent is the location of the comma after
  ! the field currently being read. 
  ! i is the index of the arrays being accessed.
  Integer :: iCurrent, i          
  Integer :: year, arrayLen, stat
  integer :: ftn 
  integer :: status

  !This initializes the current line, first line, and filename variables. 
  Character (len = 300) :: currentLine
  Character (len=LIS_CONST_PATH_LEN) :: filename 
  character*4           :: fyr
  logical               :: file_exists

  type(ESMF_Time)       :: amerifluxtime
  real                  :: gmt
  integer               :: zone

  write(unit=fyr,fmt='(i4.4)') LIS_rc%yr
  filename = trim(AmerifluxObs_struc(n)%odir)//'/'//trim(AmerifluxObs_struc(n)%site_name(stn_index))&
       //'/'//trim(adjustl(AmerifluxObs_struc(n)%stn_name(stn_index)))//trim(fyr)//'_L3.txt'
  
  inquire(file=trim(filename),exist=file_exists)

  if(file_exists) then 
      write(LIS_logunit,*) 'Reading ',trim(filename)
      year = LIS_rc%yr
  
      !Calculates whether or not the current year is a leap year following the
      !method outlined by Microsoft. (Since the data is all recent, this 
      !method could simply be replaced with Mod(year, 4) == 0 but this is far
      !more accurate)
      !The value 17520 is 356 days * 24 hours * 2 measurements per hour
      !Likewise, 17568 is 366 days * 24 hours * 2 measurements per hour
      if (Mod(year, 4) == 0) Then            !Step1
         if (Mod(year, 100) == 0) Then       !Step2
            if (Mod(year, 400) == 0) Then    !Step3
               arrayLen = 17568              !Step4, leap year
            else                          
               arrayLen = 17520              !Step5, not leap year
            end if
         else
            arrayLen = 17568                 !Step4, leap year
         end if
      else
         arrayLen = 17520                    !Step5, not leap year
      End if
    
      !Opens the file and promptly discards the header
      ftn = LIS_getNextUnitNumber()
      Open (ftn, file = filename)
      READ (ftn, *)
  
      i = 1

      ! Repeatedly reads a line and processes it
      DO         
         READ (ftn, "(A300)", iostat = stat) currentLine
         if (stat /= 0) exit       
	
         ! looks for the first available comma and indexes it
         ! reads what is between the start of the line and the
         ! comma as the next value sets the current line to the
         ! current line minus that value and its comma. Values 
         ! that were not used are not read into any array but 
         ! are simply discarded. Basically the entire idea 
         ! revolves around the format of the data:
         ! value1,value2,value3,etc


         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *) Month
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *) Day
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         READ(currentLine(1: iCurrent - 1), *) Hour
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) DoY
         currentLine = currentLine(iCurrent + 1: Len(currentLine)) 
     
         iCurrent = Index(currentLine, ",") !CO2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !H2O
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !ZL
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !FC
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_Fc
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !NEE_st
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_NEE_st
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !NEE_or
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_NEE_or
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) H !AmerifluxObs_struc(n)%Qh(stn_index,i)
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Le !AmerifluxObs_struc(n)%Qle(stn_index,i)
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !ustar
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !qf_ust
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Precip
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rg
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !PPFD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !R_pot
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) qf_Rg
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) qf_Rad
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !Rr
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rn
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rd
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !APAR
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Ta
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Ts1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Ts2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) SWC1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) SWC2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) G1
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) G2
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",")
         Read(currentLine(1: iCurrent - 1), *) Rh
         currentLine = currentLine(iCurrent + 1: Len(currentLine))
     
         iCurrent = Index(currentLine, ",") !WD
         currentLine = currentLine(iCurrent + 1: Len(currentLine))

         ! by the last line the current line consists of only a value,
         ! no comma indexing is required
         Read(currentLine, *) WS
         
         min = (hour - floor(hour))*60.0
         sec = (min- floor(min))*60.0
         
         call LIS_localtime2gmt(gmt, AmerifluxObs_struc(n)%stnlon(stn_index), hour, zone)         
         call ESMF_TimeSet(amerifluxtime, yy=year, &
              mm = nint(month), &
              dd = nint(day), &
              h = floor(gmt), &
              m = nint(min), &
              s = nint(sec), & 
              calendar = LIS_calendar, &
              rc=status)
         call LIS_verify(status)
         
         n = nint((amerifluxtime - amerifluxobs_struc(n)%starttime)/amerifluxobs_struc(n)%timestep) + 1
!         print*, year, month, day, floor(gmt), min, sec, n
         
         if (H.ne.-9999.0000) Then
            AmerifluxObs_struc(n)%Qh(stn_index, n) = H
         end if
         
         if (Le.ne.-9999.0000) Then
            AmerifluxObs_struc(n)%Qle(stn_index, n) = Le
         end if
        
         !The current method only logs G if both values are present
         !This is sometimes ineffective since one of the two is
         !usually inaccurate, so I have included the code to accept
         !either one or the other.
         if ((G1.ne.-9999.0000).and.(G2.ne.-9999.0000)) Then
            AmerifluxObs_struc(n)%Qg(stn_index, n) = (G1 + G2)/2
        ! else if (G1.ne.-9999.0000) Then
        !    AmerifluxObs_struc(n)%Qg(stn_index, n) = G1   
        ! else if (G2.ne.-9999.0000) Then
        !    AmerifluxObs_struc(n)%Qg(stn_index, n) = G2            
         end if

!         if(Ta.ne.-9999.0000) then 
!            AmerifluxObs_struc(n)%Ta(stn_index,n) = Ta + 273.15
!         endif

         if(AmerifluxObs_struc(n)%nstlayers.eq.1) then             
            if (ts1.ne.-9999.0000) then 
               AmerifluxObs_struc(n)%sfst(stn_index,n) = &
                    (AmerifluxObs_struc(n)%sfst_wt(stn_index,1)*(ts1+273.15))
            endif
         elseif(AmerifluxObs_struc(n)%nstlayers.eq.2) then                
            if ((ts1.ne.-9999.0000.and.ts2.ne.-9999.0)) then 
               AmerifluxObs_struc(n)%sfst(stn_index,n) = &
                    (AmerifluxObs_struc(n)%sfst_wt(stn_index,1)*(ts1+273.15)+&
                    AmerifluxObs_struc(n)%sfst_wt(stn_index,2)*(ts2+273.15) )
            endif
         endif

         if(AmerifluxObs_struc(n)%nsmlayers.eq.1) then             
            if (swc1.ne.-9999.0000.and.swc1/100.0.le.0.5) then 
               AmerifluxObs_struc(n)%sfsm(stn_index,n) = &
                    (AmerifluxObs_struc(n)%sfsm_wt(stn_index,1)*swc1/100.0)
               if(stn_index.eq.40) print*, 'sm1',amerifluxobs_struc(n)%sfsm(stn_index,n)
               if(AmerifluxObs_struc(n)%sfsm(stn_index,n).gt.0.5) then 
                  print*, 'Note that soil moisture > 0.5 '
                  print*, ' -- likely due to the reporting being '
                  print*, 'in percentages. We need to multiply it '
                  print*, 'with the porosity for this location'
                  print*, stn_index,AmerifluxObs_struc(n)%stn_name(stn_index), &
                       AmerifluxObs_struc(n)%sfsm(stn_index,n), swc1                  
                  stop
               endif

            endif
         elseif(AmerifluxObs_struc(n)%nsmlayers.eq.2) then                
            if ((swc1.ne.-9999.0000.and.swc2.ne.-9999.0).and.&
                 (swc1/100.0.le.0.5.and.swc2/100.0.le.0.5))then 
               AmerifluxObs_struc(n)%sfsm(stn_index,n) = &
                    (AmerifluxObs_struc(n)%sfsm_wt(stn_index,1)*swc1/100.0+&
                    AmerifluxObs_struc(n)%sfsm_wt(stn_index,2)*swc2/100.0 )
               if(stn_index.eq.40) print*, 'sm2',amerifluxobs_struc(n)%sfsm(stn_index,n)
               if(AmerifluxObs_struc(n)%sfsm(stn_index,n).gt.0.5) then 
                  print*, 'Note that soil moisture > 0.5 '
                  print*, ' -- likely due to the reporting being '
                  print*, 'in percentages. We need to multiply it '
                  print*, 'with the porosity for this location'
                  print*, stn_index,n,AmerifluxObs_struc(n)%sfsm(stn_index,n), &
                       swc1, swc2, AmerifluxObs_struc(n)%sfsm_wt(stn_index,1), &
                       AmerifluxObs_struc(n)%sfsm_wt(stn_index,2)
               endif
            endif
         endif
         
!         if (precip.ne.-9999.0000) Then
!            AmerifluxObs_struc(n)%precip(stn_index, n) = precip/1800.0 !mm/s
!         end if
         
         i= i + 1
      END DO

      call LIS_releaseUnitNumber(ftn)

   end if

 end subroutine read_ameriflux_station
