!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! !REVISION HISTORY: 
! 27 Apr 2020: Augusto Getirana, Initial implementation
!                                
module HYMAP2_urbanMod
  !module defining urban drainage variables
contains
  ! ================================================               
  subroutine HYMAP2_calc_urban_drain_stonxt(zdt,zdrvel,zdrtotwth,zdrstomax,zdrout,zdrinf,rivdph,zdrsto,zrivsto)
    implicit none
    real, intent(in) :: zdt
    real, intent(in) :: zdrvel
    real, intent(in) :: zdrtotwth
    real, intent(in) :: zdrstomax
    real, intent(in) :: zdrout
    real, intent(in) :: zdrinf
    real, intent(in) :: rivdph
    real, intent(out) :: zrivsto
    real, intent(out) :: zdrsto
    
    real :: zdrcapnow
    real :: zdrflwpot
    real :: zsrfstoin
    
    !update urban drainage water storage with inflow and outflow
    zdrsto=max(0.,zdrsto-zdrout*zdt)
    zdrsto=zdrsto+zdrinf*zdt
    
    !interaction with surface
    if(zdrsto>zdrstomax)then
      !transfer water excess to surface
      zrivsto=zrivsto+max(0.,zdrsto-zdrstomax)
      zdrsto=zdrstomax
    else
      !water intake through gutters
      zdrcapnow=max(0.,zdrstomax-zdrsto)
      zdrflwpot=min(zrivsto,rivdph*zdrvel*zdrtotwth*zdt)
      zsrfstoin=min(zdrcapnow,zdrflwpot)
      zrivsto=zrivsto-zsrfstoin
      zdrsto=zdrsto+zsrfstoin
    endif    
  end subroutine HYMAP2_calc_urban_drain_stonxt
  ! ================================================               
  subroutine  HYMAP2_calc_urb_drain_out(zdt,zdrrad,zdrman,zdrslp,zdrstomax,zdrtotlgh,zdrnoutlet,zdrsto,zdrout)
    implicit none
    real, intent(in) :: zdt
    real, intent(in) :: zdrrad
    real, intent(in) :: zdrman
    real, intent(in) :: zdrslp
    real, intent(in) :: zdrstomax
    real, intent(in) :: zdrtotlgh
    real, intent(in) :: zdrnoutlet
    real, intent(in) :: zdrsto
    !real, intent(in) :: zdrsto_down
    real, intent(out) :: zdrout

    real :: zdph ! channel water depth [m]
    real :: zhrad !channel hydraulic radius [m]
    real :: za !channel wet area [m2]
 
    if(zdrsto>0)then
      zdph=zdrsto/zdrrad/zdrtotlgh
      za=zdrsto/zdrtotlgh
      zhrad=za/(2*zdph+zdrrad)
      !outflow in m3/s from a single pipe using Manning formula
      zdrout=(1./zdrman)*za*(zhrad**(2./3))*sqrt(zdrslp)
      !constrain outflow by available water storage
      zdrout=min(zdrsto/zdt,zdrout*zdrnoutlet)
      !constrain outflow by space availability (in m3) downstream
      !zdrout=min(zdrstomax_down-zdrsto_down,zdrout)
    else
      zdrout=0.
    endif
    
  end subroutine  HYMAP2_calc_urb_drain_out
    !============================================= 
  subroutine HYMAP2_get_urban_parameters(drfile,drwth,drhgt,drden,drvel,drblk,drrad,drlgh,drman,drslp) 
    
    use LIS_logMod
    implicit none
    character(*), intent(in)    :: drfile
    real, intent(out) :: drwth,drhgt,drden,drvel,drblk,drrad,drlgh,drman,drslp
    logical                     :: file_exists
    integer                     :: ftn

    inquire(file=drfile,exist=file_exists)
   
    if(file_exists)then
      ftn = LIS_getNextUnitNumber()
      open(ftn,file=trim(drfile),status='old')
      read(ftn,*)
      read(ftn,*)drwth
      read(ftn,*)drhgt
      read(ftn,*)drden
      read(ftn,*)drvel
      read(ftn,*)
      read(ftn,*)
      read(ftn,*)drblk
      read(ftn,*)drrad
      read(ftn,*)drlgh
      read(ftn,*)drman
      read(ftn,*)drslp
      close(ftn)
      call LIS_releaseUnitNumber(ftn)
    else
      write(LIS_logunit,*) '[ERR] urban drainage paramater file '//trim(drfile)
      write(LIS_logunit,*) '[ERR] does not exist'
      call LIS_endrun()
    endif 
  end subroutine HYMAP2_get_urban_parameters
    ! ================================================              
  subroutine HYMAP2_gen_urban_drain_maps(nseqall,drrad,drlgh,drden,drwth,drblk,&
                        grarea,next,flowmap,drstomax,drtotwth,drnoutlet,drtotlgh)
    implicit none
    integer, intent(in)   :: nseqall
    real,    intent(in)   :: drrad,drlgh !drainage width [m] and length density [m/m2]
    real,    intent(in)   :: drden,drwth,drblk  ! gutter density [units/m2], width [m] and street block length [m]
    real,    intent(in)   :: grarea(nseqall) !grid cell area
    integer, intent(in)   :: next(nseqall)   !downstream grid cell
    real,    intent(in)   :: flowmap(nseqall)   !downstream grid cell
    real,    intent(out)  :: drstomax(nseqall) !max. urban drainage storage [m3]
    real,    intent(out) :: drtotwth(nseqall) !sum of gutter width within HyMAP grid cell [m]
    real,    intent(out)  :: drnoutlet(nseqall) ! average number of drainage outlets within a grid cell [-]
    real,    intent(out)  :: drtotlgh(nseqall) ! total urban drainage network length within a grid cell [m]
    integer               :: i,j

  do i=1,nseqall
    if(flowmap(i)==4)then
      !get max. drainage storage
      drstomax(i)=(drrad**2)*drlgh*grarea(i)
      
      !sum of gutter width within grid cell
      drtotwth(i)=drden*drwth*grarea(i)
      
      !number of outlets within grid cell
      drnoutlet(i)=2*sqrt(grarea(i))/drblk

      !total drainage length within a grid cell
      drtotlgh(i)=drlgh*grarea(i)
    else
      drstomax(i)=0.
      drtotwth(i)=0.
      drnoutlet(i)=0.
      drtotlgh(i)=0.
    endif
   enddo

 end subroutine HYMAP2_gen_urban_drain_maps
  ! ================================================              
end module HYMAP2_urbanMod
