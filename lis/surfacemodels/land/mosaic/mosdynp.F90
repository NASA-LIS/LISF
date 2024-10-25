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
! !ROUTINE: mosdynp
! \label{mosdynp}
!
! !REVISION HISTORY:
!  6  Apr 2001: Matt Rodell; Initial Specification
!  11 Feb 2002: Jon Gottschalck; Added use of AVHRR derived LAI/Greenness
!  01 Oct 2002: Jon Gottschalck; Modified to allow for MODIS LAI      
!  25 Sep 2007: Sujay Kumar; Upgraded for LIS5.0
!
! !INTERFACE:
subroutine mosdynp(n)
! !USES:      
  use LIS_coreMod,    only : LIS_rc,LIS_domain,LIS_surface
  use LIS_timeMgrMod, only : LIS_date2time
  use LIS_vegDataMod, only : LIS_gfrac,LIS_lai, LIS_sai, LIS_roughness
  use mos_lsmMod      
  
  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
!
! !DESCRIPTION:
!  This subroutine take all the monthly varying parameters
!  and the date and determine the actual value of the parameter for that date
!  this actual value is returned to the main program
!  The assumption is that the data point is valid for the 16th
!  of the given month at 00hr
!EOP

  integer :: index,i,j,k
  integer :: p,t                 ! loop counters
  real*8  :: time1,time2         ! temporary time variables
  integer :: yr1,mo1,yr2,mo2     ! temporary time variables
  integer :: nhmo1,nhmo2         ! temp var - n. hemisphere equiv month
  integer :: doy1,doy2           ! temporary time variables
  real    :: wt1,wt2,gmt1,gmt2   ! interpolation weights
  integer :: zeroi,numi          ! integer number holders
  real    :: valuemon(LIS_rc%nvegtypes,mos_struc(n)%mos_nmvegp,12)
  real    :: vegmp(LIS_rc%npatch(n,LIS_rc%lsm_index),mos_struc(n)%mos_nmvegp,12)   
  integer :: mapVegToUMD

  zeroi=0
  numi=16

  if(LIS_rc%da.lt.16)then
     mo1=LIS_rc%mo-1
     yr1=LIS_rc%yr 
     if(mo1.eq.0)then
        mo1=12
        yr1=LIS_rc%yr-1
     endif
     mo2=LIS_rc%mo
     yr2=LIS_rc%yr
  else
     mo1=LIS_rc%mo
     yr1=LIS_rc%yr
     mo2=LIS_rc%mo+1
     yr2=LIS_rc%yr
     if(mo2.eq.13)then
        mo2=1
        yr2=LIS_rc%yr+1
     endif
  endif
  call LIS_date2time(time1,doy1,gmt1,yr1,mo1, &
       numi,zeroi,zeroi,zeroi)
  call LIS_date2time(time2,doy2,gmt2,yr2,mo2, &
       numi,zeroi,zeroi,zeroi)
  wt1= (time2-LIS_rc%time)/(time2-time1)
  wt2= (LIS_rc%time-time1)/(time2-time1)


  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     if (mos_struc(n)%mos(t)%vegt .eq. LIS_rc%bareclass) then
        if(LIS_rc%uselaimap(n).ne."none".and.&
             LIS_rc%uselaimap(n).ne."LDT") then 
           mos_struc(n)%mos(t)%lai = LIS_lai(n)%tlai(t)
        else
           mos_struc(n)%mos(t)%lai = 0.001
           mos_struc(n)%mos(t)%dsai = 0.001
        endif

        if(LIS_rc%usegreennessmap(n).ne."none".and.&
             LIS_rc%usegreennessmap(n).ne."LDT") then 
           mos_struc(n)%mos(t)%green = LIS_gfrac(n)%greenness(t)
        else
           mos_struc(n)%mos(t)%green = 0.001
        endif
     else
        mos_struc(n)%mos(t)%lai = LIS_lai(n)%tlai(t)
! if gfrac is climo, then Mosaic calculated green. So we need to read
! SAI. 
        if(LIS_rc%usegreennessmap(n).ne."none") then 
           mos_struc(n)%mos(t)%green = LIS_gfrac(n)%greenness(t)
        else
           mos_struc(n)%mos(t)%dsai = lis_sai(n)%tsai(t)
           if ((mos_struc(n)%mos(t)%lai) .ne. 0.0.or.&
                mos_struc(n)%mos(t)%dsai.ne.0.0) then 
              mos_struc(n)%mos(t)%green = mos_struc(n)%mos(t)%lai / &
                   (mos_struc(n)%mos(t)%lai + mos_struc(n)%mos(t)%dsai)
           else
              mos_struc(n)%mos(t)%green = 0.001
           endif
        endif

        if(LIS_rc%uselaimap(n).ne."none".and.&
             LIS_rc%usesaimap(n).ne."none") then !climatology
           mos_struc(n)%mos(t)%lai = lis_lai(n)%tlai(t) + lis_sai(n)%tsai(t)
        elseif(LIS_rc%uselaimap(n).ne."none") then 
           mos_struc(n)%mos(t)%lai = lis_lai(n)%tlai(t)
        endif
    endif
 enddo

 open(unit=16,file=mos_struc(n)%mos_mvfile,status='old')
 do j=1,mos_struc(n)%mos_nmvegp
    read(16,*)
    do k=1,12
       read(16,*)(valuemon(i,j,k),i=1,LIS_rc%nvegtypes)
    enddo
 enddo
 close(16)
 
 do i=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     do j=1,mos_struc(n)%mos_nmvegp
        do k=1,12
           vegmp(i,j,k)=valuemon(mos_struc(n)%mos(i)%vegt,j,k)
        enddo 
     enddo 
  enddo 
  
  
  do t=1,LIS_rc%npatch(n,LIS_rc%lsm_index)
     
!      if tile is in the southern hemisphere, use veg parameters from the 
!         opposite (6 months away) time of year.
!         index = gindex(tile(t)%col, tile(t)%row)
     index = LIS_surface(n,LIS_rc%lsm_index)%tile(t)%index
     if (LIS_domain(n)%grid(index)%lat .lt. 0.0) then
        nhmo1 = mod((mo1+6),12)
        if (nhmo1 .eq. 0) nhmo1=12
        nhmo2 = mod((mo2+6),12)
        if (nhmo2 .eq. 0) nhmo2=12
     else
        nhmo1 = mo1
        nhmo2 = mo2
     end if

     if (LIS_rc%uselaimap(n).ne."none") then
        
        mos_struc(n)%mos(t)%vegip(1) = mos_struc(n)%mos(t)%green
        mos_struc(n)%mos(t)%vegip(2) = mos_struc(n)%mos(t)%lai
        if (mos_struc(n)%mos(t)%vegt.eq. 12) then
           mos_struc(n)%mos(t)%vegip(1) = 0.001
           mos_struc(n)%mos(t)%vegip(2) = 0.001
        endif
        
        do p=3,mos_struc(n)%mos_nmvegp
           mos_struc(n)%mos(t)%vegip(p) = wt1*vegmp(t,p,nhmo1) + &
                wt2*vegmp(t,p,nhmo2)
        enddo
        
        if(LIS_rc%useroughnessmap(n).ne."none") then 
           mos_struc(n)%mos(t)%vegip(3) = LIS_roughness(n)%roughness(t)
        endif
     else
        
        do p=1,mos_struc(n)%mos_nmvegp
           mos_struc(n)%mos(t)%vegip(p) = wt1*vegmp(t,p,nhmo1) +  &
                wt2*vegmp(t,p,nhmo2)
        enddo
        
     endif

  enddo

  return
end subroutine mosdynp
 
