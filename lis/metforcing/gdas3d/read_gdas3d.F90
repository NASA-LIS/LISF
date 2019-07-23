!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: read_gdas3d
! \label{read_gdas3d}
!
! !REVISION HISTORY:
!  19 Mar 2009: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine read_gdas3d(n, m, order, sanlfile)
! !USES:  
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod,  only : LIS_getNextUnitNumber, LIS_releaseUnitNumber, &
       LIS_logunit, LIS_endrun
  use gdas3d_forcingMod, only : gdas3d_struc

  implicit none

  
! !ARGUMENTS:
!
  integer,   intent(in) :: n 
  integer,   intent(in) :: m
  integer,   intent(in) :: order
  character(len=*)      :: sanlfile

! !DESCRIPTION:
!
!EOP
  
  real,      allocatable    :: lat(:,:), lon(:,:)
  real,      allocatable    :: p_level(:,:,:)
  real,      allocatable    :: p_levels(:,:,:)
  real,      allocatable    :: l_temp(:,:,:)
  real,      allocatable    :: l_q2(:,:,:)
  real,      allocatable    :: l_o3(:,:,:)
  real,      allocatable    :: u_wind(:,:)
  real,      allocatable    :: v_wind(:,:)

  integer               :: ftn
  logical               :: file_exists
  integer               :: irec,ivar,iret,k,c,r
  logical*1, allocatable    :: lb(:)
  real,      allocatable    :: gi(:)
  real                  :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  logical*1             :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))

! read sanl file (atm)
  
  inquire(file=trim(sanlfile),exist=file_exists)

  if(.not.file_exists) then 
     write(LIS_logunit,*) 'GDAS sanl file',trim(sanlfile), ' does not exist'
     write(LIS_logunit,*) 'Program stopping ..'
     call LIS_endrun()
  endif

  allocate(lat(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr))
  allocate(lon(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr))
  allocate(p_level(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr,gdas3d_struc(n)%nlayer+1))
  allocate(p_levels(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr,gdas3d_struc(n)%nlayer))
  allocate(l_temp(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr,gdas3d_struc(n)%nlayer))
  allocate(l_q2(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr,gdas3d_struc(n)%nlayer))
  allocate(l_o3(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr,gdas3d_struc(n)%nlayer))
  allocate(u_wind(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr))
  allocate(v_wind(gdas3d_struc(n)%nc,gdas3d_struc(n)%nr))

  allocate(lb(gdas3d_struc(n)%nc*gdas3d_struc(n)%nr))
  allocate(gi(gdas3d_struc(n)%nc*gdas3d_struc(n)%nr))

  ftn = LIS_getNextUnitNumber()
  
  open(ftn,file=trim(sanlfile),status='old',form='unformatted',&
       access='direct',recl=gdas3d_struc(n)%nc*gdas3d_struc(n)%nr*4)

  irec = 1
  read(ftn,rec=irec) lon
  irec=irec+1
  read(ftn,rec=irec) lat

! level pressure (gdas3d_struc(n)%nlayer+1 levels) (mb)
  
  do k=1,gdas3d_struc(n)%nlayer+1
     irec=irec+1
     read(ftn,rec=irec) p_level(:,:,k)
  enddo
  p_level = p_level/100.0 ! ??

! p_levs is the delta pressure. Computing level pressure from 
! delta pressure. 

!Level Pressure
  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     lb = .true. 
     ivar =0
     do k=1,gdas3d_struc(n)%nlayer+1
        ivar = ivar+1
        do r=1,gdas3d_struc(n)%nr
           do c=1,gdas3d_struc(n)%nc
              gi(c+(r-1)*gdas3d_struc(n)%nc) = p_level(c,r,k)
           enddo
        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
             lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
             gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
             gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
             gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
             LIS_rc%udef, iret)
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then 
                    gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = &
                         go(c+(r-1)*LIS_rc%lnc(n))
                 elseif(order.eq.2) then 
                    gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = &
                         go(c+(r-1)*LIS_rc%lnc(n))
                 endif
              endif
           enddo
        enddo
               
     enddo
  endif

  do k=1,gdas3d_struc(n)%nlayer
     p_levels(:,:,k) = 0.5*(p_level(:,:,k)+p_level(:,:,k+1))
  enddo
  
!Atmospheric pressure. 
  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     lb = .true. 
     
     do k=1,gdas3d_struc(n)%nlayer
        ivar = ivar+1
        do r=1,gdas3d_struc(n)%nr
           do c=1,gdas3d_struc(n)%nc
              gi(c+(r-1)*gdas3d_struc(n)%nc) = p_levels(c,r,k)
           enddo
        enddo
        
!        do c=1,LIS_rc%lnc(n)*LIS_rc%lnr(n)
!           print*, 'lat/lon ',c,gdas3d_struc(n)%rlat1(c),&
!                gdas3d_struc(n)%rlon1(c),&
!                gdas3d_struc(n)%w111(c),  gdas3d_struc(n)%w121(c), &
!                gdas3d_struc(n)%w211(c),  gdas3d_struc(n)%w221(c), &
!                gdas3d_struc(n)%n111(c),  gdas3d_struc(n)%n121(c), &
!                gdas3d_struc(n)%n211(c),  gdas3d_struc(n)%n221(c)
!        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
             lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
             gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
             gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
             gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
             LIS_rc%udef, iret)
        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then 
                    gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = &
                         go(c+(r-1)*LIS_rc%lnc(n))
                 elseif(order.eq.2) then 
                    gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = &
                         go(c+(r-1)*LIS_rc%lnc(n))
                 endif
              endif
           enddo
        enddo
               
     enddo
  endif

!  read other atmospheric profiles

! layer temperature (K)
  do k=1,gdas3d_struc(n)%nlayer
     irec=irec+1
     read(ftn,rec=irec) l_temp(:,:,k)
  enddo

  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     lb = .true. 
     do k=1,gdas3d_struc(n)%nlayer
        ivar = ivar+1

        do r=1,gdas3d_struc(n)%nr
           do c=1,gdas3d_struc(n)%nc
              gi(c+(r-1)*gdas3d_struc(n)%nc) = l_temp(c,r,k)
           enddo
        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
             lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
             gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
             gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
             gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
             LIS_rc%udef, iret)

        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then 
                    gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 elseif(order.eq.2) then 
                    gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 endif
              endif
           enddo
        enddo
               
     enddo
  endif

! layer specific humidity (kg/kg)
  do k=1,gdas3d_struc(n)%nlayer
     irec=irec+1
     read(ftn,rec=irec) l_q2(:,:,k)
  enddo

  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     lb = .true. 
     do k=1,gdas3d_struc(n)%nlayer

        ivar = ivar+1

        do r=1,gdas3d_struc(n)%nr
           do c=1,gdas3d_struc(n)%nc
              gi(c+(r-1)*gdas3d_struc(n)%nc) = l_q2(c,r,k)
           enddo
        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
             lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
             gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
             gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
             gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
             LIS_rc%udef, iret)

        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then 
                    gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 elseif(order.eq.2) then 
                    gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 endif
              endif
           enddo
        enddo
               
     enddo
  endif


! layer O3 mixing ratio (kg/kg)
  do k=1,gdas3d_struc(n)%nlayer
     irec=irec+1
     read(ftn,rec=irec) l_o3(:,:,k)
  enddo

  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     lb = .true. 
     do k=1,gdas3d_struc(n)%nlayer

        ivar = ivar+1

        do r=1,gdas3d_struc(n)%nr
           do c=1,gdas3d_struc(n)%nc
              gi(c+(r-1)*gdas3d_struc(n)%nc) = l_o3(c,r,k)
           enddo
        enddo
        
        call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
             lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
             LIS_rc%lnc(n)*LIS_rc%lnr(n), &
             LIS_domain(n)%lat,LIS_domain(n)%lon,&
             gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
             gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
             gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
             gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
             LIS_rc%udef, iret)

        do r=1,LIS_rc%lnr(n)
           do c=1,LIS_rc%lnc(n)
              if(LIS_domain(n)%gindex(c,r).ne.-1) then 
                 if(order.eq.1) then 
                    gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 elseif(order.eq.2) then 
                    gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
                 endif
              endif
           enddo
        enddo
               
     enddo
  endif

! skip the next variable
  irec = irec+gdas3d_struc(n)%nlayer
  
! near surface wind 
  irec=irec+1
  read(ftn,rec=irec) u_wind

  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     ivar = ivar+1
     lb = .true. 
     do r=1,gdas3d_struc(n)%nr
        do c=1,gdas3d_struc(n)%nc
           gi(c+(r-1)*gdas3d_struc(n)%nc) = u_wind(c,r)
        enddo
     enddo
     
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
          lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_domain(n)%lat,LIS_domain(n)%lon,&
          gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
          gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
          gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
          gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
          LIS_rc%udef, iret)
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
              elseif(order.eq.2) then 
                 gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
              endif
           endif
        enddo
     enddo
               
  endif


  irec=irec+1
  read(ftn,rec=irec) v_wind

  if(trim(LIS_rc%met_interp(m)).eq."bilinear") then 
     ivar = ivar+1
     lb = .true. 
     do r=1,gdas3d_struc(n)%nr
        do c=1,gdas3d_struc(n)%nc
           gi(c+(r-1)*gdas3d_struc(n)%nc) = v_wind(c,r)
        enddo
     enddo
     
     call bilinear_interp(LIS_rc%gridDesc(n,:),lb,gi,&
          lo, go, gdas3d_struc(n)%nc*gdas3d_struc(n)%nr,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_domain(n)%lat,LIS_domain(n)%lon,&
          gdas3d_struc(n)%w111,  gdas3d_struc(n)%w121, &
          gdas3d_struc(n)%w211,  gdas3d_struc(n)%w221, &
          gdas3d_struc(n)%n111,  gdas3d_struc(n)%n121, &
          gdas3d_struc(n)%n211,  gdas3d_struc(n)%n221, &
          LIS_rc%udef, iret)
     
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(LIS_domain(n)%gindex(c,r).ne.-1) then 
              if(order.eq.1) then 
                 gdas3d_struc(n)%metdata1(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
              elseif(order.eq.2) then 
                 gdas3d_struc(n)%metdata2(ivar,LIS_domain(n)%gindex(c,r)) = go(c+(r-1)*LIS_rc%lnc(n))
              endif
           endif
        enddo
     enddo
               
  endif
!  print*, 'total ivar ',ivar
  call LIS_releaseUnitNumber(ftn)

  write(LIS_logunit,*) 'Read ',trim(sanlfile), 'successfully '
 

  deallocate(lb)
  deallocate(gi)

  deallocate(lat)
  deallocate(lon)
  deallocate(p_level)
  deallocate(p_levels)
  deallocate(l_temp)
  deallocate(l_q2)
  deallocate(l_o3)
  deallocate(u_wind)
  deallocate(v_wind)

end subroutine read_gdas3d

