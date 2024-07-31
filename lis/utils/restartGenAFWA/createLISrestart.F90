!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
! This program generates a global 1/2 degree restart file for LIS from 
! AGRMET outputx

program createLISrestart
  
  implicit none
  integer, parameter :: imax=512, jmax=512
  integer, parameter :: nc = 720, nr=300
  real, parameter :: xmeshl = 47.625
  real, parameter :: xpnmcaf = 257
  real, parameter :: ypnmcaf = 257
  character*100 :: filename1
  integer       :: ip
  real       :: udef1
  real       :: gi(2,imax,jmax)
  real       :: t1(nc,nr)
  real       :: cmc(nc,nr)
  real       :: snowh(nc,nr)
  real       :: sneqv(nc,nr)
  real       :: stc1(nc,nr)
  real       :: stc2(nc,nr)
  real       :: stc3(nc,nr)
  real       :: stc4(nc,nr)
  real       :: smc1(nc,nr)
  real       :: smc2(nc,nr)
  real       :: smc3(nc,nr)
  real       :: smc4(nc,nr)
  real       :: sh2o1(nc,nr)
  real       :: sh2o2(nc,nr)
  real       :: sh2o3(nc,nr)
  real       :: sh2o4(nc,nr)
  real       :: ch(nc,nr)
  real       :: cm(nc,nr)
  real, allocatable :: tt1(:)
  real, allocatable       :: tcmc(:)
  real, allocatable       :: tsnowh(:)
  real, allocatable       :: tsneqv(:)
  real, allocatable       :: tstc1(:)
  real, allocatable       :: tstc2(:)
  real, allocatable       :: tstc3(:)
  real, allocatable       :: tstc4(:)
  real, allocatable       :: tsmc1(:)
  real, allocatable       :: tsmc2(:)
  real, allocatable       :: tsmc3(:)
  real, allocatable       :: tsmc4(:)
  real, allocatable       :: tsh2o1(:)
  real, allocatable       :: tsh2o2(:)
  real, allocatable       :: tsh2o3(:)
  real, allocatable       :: tsh2o4(:)
  real, allocatable       :: tch(:)
  real, allocatable       :: tcm(:)
  real       :: mask(nc,nr)
  character*10 :: time
  character*8 :: date
  integer :: count,c,r,ntiles

  date = '20070101'
  time = '2007010100'

! Data read
  filename1 = "./AGRMET/FLUX3/"//date//'/skntmp_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/skntmp_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA tmp data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,t1)
  filename1 = "./AGRMET/FLUX3/"//date//'/cnpymc_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/cnpymc_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA cmc data'

  ip = 2
  udef1 = 0
  
  call interp_noahfield(ip,udef1,gi,cmc)

  filename1 = "./AGRMET/FLUX3/"//date//'/snodep_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/snodep_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA snowh data'

  ip = 2
  udef1 = 0
  
  call interp_noahfield(ip,udef1,gi,snowh)

  filename1 = "./AGRMET/FLUX3/"//date//'/snoeqv_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/snoeqv_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA sneqv data'

  ip = 2
  udef1 = 0
  
  call interp_noahfield(ip,udef1,gi,sneqv)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly1_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly1_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA stc1 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,stc1)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly2_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly2_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA stc2 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,stc2)
  
  filename1 = "./AGRMET/FLUX3/"//date//'/stcly3_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly3_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA stc3 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,stc3)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly4_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/stcly4_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA stc4 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,stc4)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly1_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly1_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA smc1 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,smc1)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly2_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly2_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA smc2 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,smc2)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly3_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly3_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA smc3 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,smc3)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly4_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smtly4_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA smc4 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,smc4)


  filename1 = "./AGRMET/FLUX3/"//date//'/smlly1_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smlly1_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA sh2o1 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,sh2o1)

  filename1 = "./AGRMET/FLUX3/"//date//'/smlly2_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smlly2_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA sh2o2 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,sh2o2)


  filename1 = "./AGRMET/FLUX3/"//date//'/smlly3_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smlly3_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA sh2o3 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,sh2o3)

  filename1 = "./AGRMET/FLUX3/"//date//'/smlly4_nh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/smlly4_sh.03hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA sh2o4 data'

  ip = 2
  udef1 = -1
  
  call interp_noahfield(ip,udef1,gi,sh2o4)

  filename1 = "./AGRMET/FLUX3/"//date//'/coeffh_nh.06hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/coeffh_sh.06hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA ch data'

  ip = 2
  udef1 = 0
  ch = 0 
  call interp_noahfield(ip,udef1,gi,ch)


  filename1 = "./AGRMET/FLUX3/"//date//'/coeffm_nh.06hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(1,:,:)
  close(100)

  filename1 = "./AGRMET/FLUX3/"//date//'/coeffm_sh.06hr.'//time
  open(100,file=filename1,form='unformatted',access='direct',recl=imax*jmax*4)
  read(100,rec=1) gi(2,:,:)
  close(100)
  
  print*,'Read the AFWA cm data'

  ip = 2
  udef1 = 0
  
  call interp_noahfield(ip,udef1,gi,cm)


  open(100,file='landmask.1gd4r',access='direct',recl=512*512*4)
  read(100,rec=1) mask
  close(100)

  count = 0 
  do r=1,nr
     do c=1,nc
        if(mask(c,r).ne.0) then 
           count = count+1
        endif
     enddo
  enddo
  print*,'Number of tiles ',count
  ntiles = count

  allocate(tt1(ntiles))
  allocate(tcmc(ntiles))
  allocate(tsnowh(ntiles))
  allocate(tsneqv(ntiles))
  allocate(tstc1(ntiles))
  allocate(tstc2(ntiles))
  allocate(tstc3(ntiles))
  allocate(tstc4(ntiles))
  allocate(tsmc1(ntiles))
  allocate(tsmc2(ntiles))
  allocate(tsmc3(ntiles))
  allocate(tsmc4(ntiles))
  allocate(tsh2o1(ntiles))
  allocate(tsh2o2(ntiles))
  allocate(tsh2o3(ntiles))
  allocate(tsh2o4(ntiles))
  allocate(tch(ntiles))
  allocate(tcm(ntiles))  
  
  count = 0 
  do r=1,nr
     do c=1,nc
        if(mask(c,r).ne.0) then 
           count = count+1
           tt1(count) = t1(c,r)
           tcmc(count) = cmc(c,r)
           tsnowh(count) = snowh(c,r)
           tsneqv(count) = sneqv(c,r)*snowh(c,r)
           tstc1(count) = stc1(c,r)
           tstc2(count) = stc2(c,r)
           tstc3(count) = stc3(c,r)
           tstc4(count) = stc4(c,r)
           tsmc1(count) = smc1(c,r)
           tsmc2(count) = smc2(c,r)
           tsmc3(count) = smc3(c,r)
           tsmc4(count) = smc4(c,r)
           tsh2o1(count) = sh2o1(c,r)
           tsh2o2(count) = sh2o2(c,r)
           tsh2o3(count) = sh2o3(c,r)
           tsh2o4(count) = sh2o4(c,r)
           tch(count) = ch(c,r)
           tcm(count) = cm(c,r)
        endif
     enddo
  enddo


!  open(100,file='landmask.bin',form='unformatted')
!  write(100) mask
!  close(100)
  print*,'writing LIS restart ','noah'//date//'.rst'
  open(100,file='noah'//date//'.rst',form='unformatted')
  write(100) nc,nr,ntiles
  write(100) tt1
  write(100) tcmc
  write(100) tsnowh
  write(100) tsneqv
  write(100) tstc1
  write(100) tstc2
  write(100) tstc3
  write(100) tstc4
  write(100) tsmc1
  write(100) tsmc2
  write(100) tsmc3
  write(100) tsmc4
  write(100) tsh2o1
  write(100) tsh2o2
  write(100) tsh2o3
  write(100) tsh2o4
  write(100) tch
  write(100) tcm
  close(100) 

end program createLISrestart

subroutine interp_noahfield(ip,udef1,gi,varfield)
  
  implicit none
  integer, parameter :: imax=512, jmax=512
  integer, parameter :: nc = 720, nr=300
  real, parameter :: xmeshl = 47.625
  real, parameter :: xpnmcaf = 257
  real, parameter :: ypnmcaf = 257
  integer :: ip
  integer :: udef1
  real    :: gi(2,imax,jmax)
  real    :: varfield(nc,nr)
  real          :: gridDesci(50)
  real          :: gridDesco(50)
  integer       :: hemi
  real :: xmesh, orient,xi1,xj1
  real :: alat1,alon1
  real, allocatable      :: rlat1_nh(:)
  real, allocatable      :: rlat1_sh(:)
  real, allocatable      :: rlon1_nh(:)
  real, allocatable      :: rlon1_sh(:)
  integer, allocatable   :: n111_nh(:)
  integer, allocatable   :: n111_sh(:)
  integer, allocatable   :: n121_nh(:)
  integer, allocatable   :: n121_sh(:)
  integer, allocatable   :: n211_nh(:)
  integer, allocatable   :: n211_sh(:)
  integer, allocatable   :: n221_nh(:)
  integer, allocatable   :: n221_sh(:)
  real, allocatable      :: w111_nh(:)
  real, allocatable      :: w121_nh(:)
  real, allocatable      :: w111_sh(:)
  real, allocatable      :: w121_sh(:)
  real, allocatable      :: w211_nh(:)
  real, allocatable      :: w221_nh(:)
  real, allocatable      :: w211_sh(:)
  real, allocatable      :: w221_sh(:)     
  real, allocatable      :: rlat2_nh(:)
  real, allocatable      :: rlon2_nh(:)
  real, allocatable      :: rlat2_sh(:)
  real, allocatable      :: rlon2_sh(:)
  integer, allocatable   :: n112_nh(:)
  integer, allocatable   :: n112_sh(:)
  integer                :: mo1,mo2
  logical*1, allocatable :: lo_nh(:),lo_sh(:)
  real, allocatable      :: go_nh(:),go_sh(:)
  real, allocatable      :: comb_data(:)
  real      :: gi_temp(imax*jmax)
  logical*1      :: li(imax*jmax)
  integer :: iret
  integer :: i,j
  integer                :: ibi,ibo

  gridDesci =  0 
  gridDesco =  0 
  
! Setup interpolation weights
  do hemi = 1, 2
     gridDesco(1) = 0 
     gridDesco(2) = nc
     gridDesco(5) = -179.75
     gridDesco(8) = 179.75
     gridDesco(6) = 128
     gridDesco(9) = 0.50
     gridDesco(10) = 0.50
     gridDesco(20) = 255
     if(hemi.eq.1) then
        gridDesco(4) = 0.25
        gridDesco(7) = 89.75
        gridDesco(3) = 180
     else
        gridDesco(4) = -59.75
        gridDesco(7) = -0.25
        gridDesco(3) = 120             
     endif
     
     if(hemi.eq.1) then 
        xmesh = xmeshl
        orient = 100.0
     else
        xmesh = -xmeshl
        orient = 280.0
     endif
     xj1 = float(1)-ypnmcaf
     xi1 = float(1)-xpnmcaf
     
     call polarToLatLon(xi1,xj1,xmesh,orient,alat1,alon1)

     gridDesci = 0 
     gridDesci(1) = 5
     gridDesci(2) = imax
     gridDesci(3) = jmax
     gridDesci(4) = alat1
     gridDesci(5) = alon1
     gridDesci(6) = 8
     gridDesci(7) = orient
     gridDesci(8) = xmesh
     gridDesci(9) = xmesh
     gridDesci(10) = 0.0
     if(hemi .eq.2) then 
        gridDesci(20) = 128
        gridDesci(11) = 128
     endif
     gridDesci(20) = 255
             
     if(hemi.eq.1) then 
        mo1 = gridDesco(2)*gridDesco(3)
        allocate(rlat1_nh(mo1))
        allocate(rlon1_nh(mo1))
        allocate(n111_nh(mo1))
        allocate(n121_nh(mo1))
        allocate(n211_nh(mo1))
        allocate(n221_nh(mo1))
        allocate(w111_nh(mo1))
        allocate(w121_nh(mo1))
        allocate(w211_nh(mo1))
        allocate(w221_nh(mo1))
        call bilinear_interp_input(gridDesci,gridDesco,mo1,rlat1_nh,&
             rlon1_nh,n111_nh,n121_nh,n211_nh,n221_nh,w111_nh,w121_nh,&
             w211_nh,w221_nh)
        allocate(rlat2_nh(mo1))
        allocate(rlon2_nh(mo1))
        allocate(n112_nh(mo1))
        call neighbor_interp_input(gridDesci,gridDesco,&
             mo1,rlat2_nh,&
             rlon2_nh,n112_nh)
     else
        mo2 = gridDesco(2)*gridDesco(3)
        allocate(rlat1_sh(mo2))
        allocate(rlon1_sh(mo2))
        allocate(n111_sh(mo2))
        allocate(n121_sh(mo2))
        allocate(n211_sh(mo2))
        allocate(n221_sh(mo2))
        allocate(w111_sh(mo2))
        allocate(w121_sh(mo2))
        allocate(w211_sh(mo2))
        allocate(w221_sh(mo2))
        call bilinear_interp_input(gridDesci,gridDesco,mo2,rlat1_sh,&
             rlon1_sh,n111_sh,n121_sh,n211_sh,n221_sh,w111_sh,w121_sh,&
             w211_sh,w221_sh)
        allocate(rlat2_sh(mo2))
        allocate(rlon2_sh(mo2))
        allocate(n112_sh(mo2))
        call neighbor_interp_input(gridDesci,gridDesco,&
             mo2,rlat2_sh,&
             rlon2_sh,n112_sh)
     endif
  enddo


!Interpolation section  
  
  if(ip.eq.1) then 
     ibi = 1
     gi_temp = udef1
     allocate(lo_nh(mo1))
     allocate(lo_sh(mo2))
     allocate(go_nh(mo1))
     allocate(go_sh(mo2))
     allocate(comb_data(mo1+mo2))
     do hemi = 1,2
        li = .false.
        do i=1,imax
           do j=1,jmax
              if(gi(hemi,i,j).ne.udef1) then 
                 li(i+(j-1)*imax) = .true.
                 gi_temp(i+(j-1)*jmax) = gi(hemi,i,j)
              endif
           enddo
        enddo
        gridDesco(1) = 0 
        gridDesco(2) = nc
        gridDesco(5) = -179.75
        gridDesco(8) = 179.75
        gridDesco(6) = 128
        gridDesco(9) = 0.50
        gridDesco(10) = 0.50
        gridDesco(20) = 255
        if(hemi.eq.1) then
           gridDesco(4) = 0.25
           gridDesco(7) = 89.75
           gridDesco(3) = 180
        else
           gridDesco(4) = -59.75
           gridDesco(7) = -0.25
           gridDesco(3) = 120
        endif
        if(hemi.eq.1) then 
           call bilinear_interp(gridDesco,ibi,li,gi_temp,&
                ibo,lo_nh,go_nh,imax*jmax,mo1,&
                rlat1_nh,rlon1_nh,&
                w111_nh,w121_nh,&
                w211_nh,w221_nh,&
                n111_nh,n121_nh,&
                n211_nh,n221_nh,-9999,iret)
        elseif(hemi.eq.2) then 
           call bilinear_interp(gridDesco,ibi,li,gi_temp,&
                ibo,lo_sh,go_sh,imax*jmax,mo2,&
                rlat1_sh, rlon1_sh,&
                w111_sh,w121_sh,&
                w211_sh,w221_sh,&
                n111_sh,n121_sh,&
                n211_sh,n221_sh,-9999,iret)
        endif        
     enddo     
     comb_data(1:mo2) = go_sh(:)
     comb_data(mo2+1:mo1+mo2) = go_nh(:)
     varfield = RESHAPE(comb_data(1:mo1+mo2),(/nc,nr/))
     deallocate(lo_nh)
     deallocate(lo_sh)
     deallocate(go_nh)
     deallocate(go_sh)
  elseif(ip.eq.2 ) then !neighbor search
     ibi = 1
     li = .false.
     gi_temp = udef1
     allocate(lo_nh(mo1))
     allocate(lo_sh(mo2))
     allocate(go_nh(mo1))
     allocate(go_sh(mo2))
     allocate(comb_data(mo1+mo2))
     do hemi = 1,2
        li = .false.
        gi_temp = udef1
        do i=1,imax
           do j=1,jmax              
              if(gi(hemi,i,j).ne.udef1) then 
                 li(i+(j-1)*imax) = .true.
                 gi_temp(i+(j-1)*jmax) = gi(hemi,i,j)
              endif
           enddo
        enddo
        gridDesco(1) = 0 
        gridDesco(2) = nc
        gridDesco(3) = nr/2
        gridDesco(5) = -179.75
        gridDesco(8) = 179.75
        gridDesco(6) = 128
        gridDesco(9) = 0.50
        gridDesco(10) = 0.50
        gridDesco(20) = 255
        if(hemi.eq.1) then
           gridDesco(4) = 0.25
           gridDesco(7) = 89.75
           gridDesco(3) = 180
        else
           gridDesco(4) = -59.75
           gridDesco(7) = -0.25
           gridDesco(3) = 120
        endif
        if(hemi.eq.1) then 
           call neighbor_interp(gridDesco,ibi,li,gi_temp,&
                ibo,lo_nh,go_nh,imax*jmax,mo1,rlat2_nh, rlon2_nh,&
                n112_nh,0.0,iret)
        elseif(hemi.eq.2) then 
           call neighbor_interp(gridDesco,ibi,li,gi_temp,&
                ibo,lo_sh,go_sh,imax*jmax,mo2,rlat2_sh, rlon2_sh,&
                n112_sh,0.0,iret)
        endif
     enddo
     comb_data(1:mo2) = go_sh(:)
     comb_data(mo2+1:mo1+mo2) = go_nh(:)
     varfield = RESHAPE(comb_data(1:mo1+mo2),(/nc,nr/))

     deallocate(lo_nh)
     deallocate(lo_sh)
     deallocate(go_nh)
     deallocate(go_sh)
  endif
  open(110,file='landmask.bin',form='unformatted')
  write(110) varfield
  close(110)
end subroutine interp_noahfield
