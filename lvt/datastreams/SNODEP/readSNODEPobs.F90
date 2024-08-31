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
! 
! !ROUTINE: readSNODEPobs
! \label{readSNODEPobs}
! 
! !REVISION HISTORY: 
!  23 APR 2010: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readSNODEPobs(source)
! !USES:   
  use ESMF
  use grib_api
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod
  use LVT_timeMgrMod
  use SNODEP_obsMod,    only : snodepobs
  use map_utils

  implicit none
! !ARGUMENTS: 

! 
! !DESCRIPTION: 
! 
! 
!EOP
  integer                :: source
  real                   :: varfield(LVT_rc%lnc,LVT_rc%lnr)
  integer                :: status,rc,iret
  real                   :: value
  real                   :: missingValue
  real                   :: gi(2,snodepobs(source)%pmax,snodepobs(source)%pmax)
  real                   :: f1(snodepobs(source)%pmax*snodepobs(source)%pmax)
  real                   :: f2(snodepobs(source)%pmax*snodepobs(source)%pmax)
  integer                :: hemi,c,r
  integer                :: igrib
  integer                :: npts
  integer                :: ftn1,ftn2
  integer                :: pds5,pds7,pds5_val,pds7_val
  character*100          :: name_nh, name_sh
  type(ESMF_Time)        :: snodeptime, snodeptime1
  logical                :: file_exists1,file_exists2,read_flag
  real                   :: timenow
  logical                :: alarmCheck

  varfield = LVT_rc%udef
  timenow = float(LVT_rc%dhr(source))*3600 + 60*LVT_rc%dmn(source) + &
       LVT_rc%dss(source)
  alarmcheck = (mod(timenow, 86400.0).eq.0)

  if(snodepobs(source)%startflag.or.alarmCheck.or.&
       LVT_rc%resetFlag(source)) then 
     
     LVT_rc%resetFlag(source) = .false. 
     snodepobs(source)%startflag = .false. 
  
     hemi = 1
     call create_SNODEP_filename(name_nh, &
          snodepobs(source)%mesh, hemi, &
          snodepobs(source)%odir, LVT_rc%dyr(source), &
          LVT_rc%dmo(source), LVT_rc%dda(source), LVT_rc%dhr(source), &
          LVT_rc%dmn(source))

     hemi = 2
     call create_SNODEP_filename(name_sh, &
          snodepobs(source)%mesh, hemi, &
          snodepobs(source)%odir, LVT_rc%dyr(source), &
          LVT_rc%dmo(source), LVT_rc%dda(source), LVT_rc%dhr(source), &
          LVT_rc%dmn(source))
     
     
     inquire(file=trim(name_nh),exist=file_exists1)
     inquire(file=trim(name_sh),exist=file_exists2)

     if(file_exists1.and.file_exists2) then 
        read_flag = .true. 
     else
        read_flag = .false.
     endif
     if(file_exists1) then 
        write(LVT_logunit,*) '[INFO] Reading SNODEP data ',&
             trim(name_nh)
        
        npts = snodepobs(source)%pmax*snodepobs(source)%pmax
        if(snodepobs(source)%mesh.eq.8) then 
           pds5 = 174 
           pds7 = 0
        else
           pds5 = 66 
           pds7 = 0
        endif
        call grib_open_file(ftn1,trim(name_nh),'r',iret)
        call LVT_verify(iret, 'error in grib_open_file in read_SNODEPobs')


        call grib_new_from_file(ftn1,igrib,iret)
        call LVT_verify(iret,'error in grib_new_from_file in read_SNODEPobs')

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LVT_verify(rc,'error in grib_get:indicatorOfParameter in read_SNODEPobs')

        call grib_get(igrib,'level',pds7_val,rc)
        call LVT_verify(rc,'error in grib_get:level in read_SNODEPobs')

        f1 = -9999.0
        if((pds5.eq.pds5_val).and.(pds7.eq.pds7_val)) then

           call grib_get(igrib,'values',f1,rc)
           call LVT_verify(rc,'error in grib_get:values in read_SNODEPobs')

           call grib_get(igrib,'missingValue',missingValue,rc)
           call LVT_verify(rc, 'error in grib_get:missingValue in read_SNODEPobs')
           
           call grib_release(igrib,rc)
           call LVT_verify(rc, 'error in grib_release in read_SNODEPobs')
        else
           write(LVT_logunit,*) '[ERR] Unable to retrieve values from ',trim(name_nh)
           call LVT_endrun()

        endif
        call grib_close_file(ftn1)

        do r=1,snodepobs(source)%pmax
           do c=1,snodepobs(source)%pmax
              value = f1(c+(r-1)*snodepobs(source)%pmax)
              gi(1,c,r) = value
              if(value.gt.408.9) gi(1,c,r) = LVT_rc%udef
           enddo
        enddo
     else
        gi(1,:,:) = LVT_rc%udef
     endif
     
     if(file_exists2) then 
        npts = snodepobs(source)%pmax*snodepobs(source)%pmax
        if(snodepobs(source)%mesh.eq.8) then 
           pds5 = 174 
           pds7 = 0
        else
           pds5 = 66 
           pds7 = 0
        endif
        write(LVT_logunit,*) '[INFO] reading SNODEP data ',name_sh

        call grib_open_file(ftn2,trim(name_sh),'r',iret)
        call LVT_verify(iret,'error in grib_open_file in read_SNODEPobs')

        call grib_new_from_file(ftn2,igrib,iret)
        call LVT_verify(iret,'error in grib_new_from_file in read_SNODEPobs')

        call grib_get(igrib,'indicatorOfParameter',pds5_val,rc)
        call LVT_verify(rc,'error in grib_get:indicatorOfParameter in read_SNODEPobs')

        call grib_get(igrib,'level',pds7_val,rc)
        call LVT_verify(rc,'error in grib_get:level in read_SNODEPobs')

        f2 = -9999.0
        if((pds5.eq.pds5_val).and.(pds7.eq.pds7_val)) then

           call grib_get(igrib,'values',f2,rc)
           call LVT_verify(rc,'error in grib_get:values in read_SNODEPobs')

           call grib_get(igrib,'missingValue',missingValue,rc)
           call LVT_verify(rc, 'error in grib_get:missingValue in read_SNODEPobs')
           
           call grib_release(igrib,rc)
           call LVT_verify(rc, 'error in grib_release in read_SNODEPobs')
        else
           write(LVT_logunit,*) '[ERR] Unable to retrieve values from ',trim(name_sh)
           call LVT_endrun()

        endif
        call grib_close_file(ftn2)

        do r=1,snodepobs(source)%pmax
           do c=1,snodepobs(source)%pmax
              value = f2(c+(r-1)*snodepobs(source)%pmax)
              gi(2,c,r) = value
              if(value.gt.408.9) gi(2,c,r) = LVT_rc%udef
           enddo
        enddo
     else
        gi(2,:,:)=LVT_rc%udef
     endif
     
     if(file_exists1.or.file_exists2) then 
        call interp_SNODEPfield(source, 1,gi,LVT_rc%udef,varfield)
     else
        varfield = LVT_rc%udef
     endif
  endif

  call LVT_logSingleDataStreamVar(LVT_MOC_snowdepth,source,&
       varfield,vlevel=1,units="m")

end subroutine readSNODEPobs


subroutine create_SNODEP_filename(name, mesh, hemi, ndir, yr, mo,da,hr,mn)

  implicit none

  character*100     :: name
  integer           :: hemi
  integer           :: mesh
  integer           :: yr, mo, da, hr,mn
  character (len=*) :: ndir
  character (len=4) :: fyr
  character (len=2) :: fmo,fda,fhr,fmn
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
  write(unit=fhr, fmt='(i2.2)') hr
  write(unit=fmn, fmt='(i2.2)') mn  

  if(mesh.eq.8) then 
     if(hemi.eq.1) then 
       name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
            '/depth_nh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
!        name = trim(ndir)//'/depth_nh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
     elseif(hemi.eq.2) then 
       name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//trim(fda)//&
            '/depth_sh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
!        name = trim(ndir)//'/depth_sh.'//trim(fyr)//trim(fmo)//trim(fda)//trim(fhr)
     endif
  elseif(mesh.eq.16) then 
     if(hemi.eq.1) then 
       name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//&
            '/SNODEP_16_NH_'//trim(fyr)//trim(fmo)//&
            trim(fda)//'12.GR1'
     elseif(hemi.eq.2) then 
       name = trim(ndir)//'/'//trim(fyr)//trim(fmo)//&
            '/SNODEP_16_SH_'//trim(fyr)//trim(fmo)//&
            trim(fda)//'12.GR1'
     endif
  endif
end subroutine Create_SNODEP_filename

!BOP
! !ROUTINE: interp_SNODEPfield
!  \label{interp_SNODEPfield}
! 
! !INTERFACE:
subroutine interp_SNODEPfield(source,ip,gi,udef,varfield)
! !USES:
  use LVT_coreMod,       only : LVT_rc
  use SNODEP_obsMod,     only : snodepobs

  implicit none
! !ARGUMENTS: 
  integer, intent(in)    :: source
  integer, intent(in)    :: ip
  real, intent(in)       :: gi(2,snodepobs(source)%pmax,snodepobs(source)%pmax)
  real, intent(in)       :: udef
  real, intent(inout)    :: varfield(LVT_rc%lnc,LVT_rc%lnr)
!
! !DESCRIPTION:
!   This subroutine interpolates a given AFWA field 
!   to the LIS domain (from polar stereographic grid to the LIS grid)
!
!  The arguments are: 
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ip]
!    interpolation option
!  \item[gi]
!    input AGRMET field (both hemispheres)
!  \item[udef]
!    undefined value in the input field
!  \item[varfield]
!    interpolated field in the LIS grid
!  \end{description}
!
!  The routines invoked are: 
!  \begin{description}
!  \item[bilinear\_interp](\ref{bilinear_interp}) \newline
!    spatially interpolate the forcing data using bilinear interpolation
!  \item[neighbor\_interp](\ref{neighbor_interp}) \newline
!    spatially interpolate the forcing data using neighbor interpolation
! \end{description}
!EOP
  real                   :: gi_temp(snodepobs(source)%mi)
  logical*1              :: li(snodepobs(source)%mi)
  integer                :: ihemi,i,j,iret
  real                   :: gridDesco(200)
  logical*1, allocatable :: lo_nh(:),lo_sh(:)
  real, allocatable      :: go_nh(:),go_sh(:)
  real, allocatable      :: comb_data(:)

  if(ip.eq.1) then 
     if(snodepobs(source)%gridspan.eq.1) then 
        allocate(lo_nh(snodepobs(source)%mo1))
        allocate(go_nh(snodepobs(source)%mo1))
     elseif(snodepobs(source)%gridspan.eq.2) then 
        allocate(lo_sh(snodepobs(source)%mo2))
        allocate(go_sh(snodepobs(source)%mo2))
     else
        allocate(lo_nh(snodepobs(source)%mo1))
        allocate(lo_sh(snodepobs(source)%mo2))
        allocate(go_nh(snodepobs(source)%mo1))
        allocate(go_sh(snodepobs(source)%mo2))
     endif
     allocate(comb_data(snodepobs(source)%mo1+snodepobs(source)%mo2))
     
     do ihemi = snodepobs(source)%shemi,snodepobs(source)%nhemi
        li = .false.
        do i=1,snodepobs(source)%pmax
           do j=1,snodepobs(source)%pmax
              if(gi(ihemi,i,j).ne.udef) then 
                 li(i+(j-1)*snodepobs(source)%pmax) = .true.
                 gi_temp(i+(j-1)*snodepobs(source)%pmax) = gi(ihemi,i,j)
!                 if(gi(ihemi,i,j).gt.30) then 
!                    print*, ihemi,i,j,gi(ihemi,i,j)
!                 endif
              endif
           enddo
        enddo
        gridDesco(1) = 0 
        gridDesco(2) = snodepobs(source)%hemi_nc(ihemi)
        gridDesco(3) = snodepobs(source)%hemi_nr(ihemi)
        gridDesco(5) = LVT_rc%gridDesc(5)
        gridDesco(8) = LVT_rc%gridDesc(8)
        gridDesco(6) = LVT_rc%gridDesc(6)
        gridDesco(9) = LVT_rc%gridDesc(9)
        gridDesco(10) = LVT_rc%gridDesc(10)
        gridDesco(20) = 0 
        if(snodepobs(source)%gridspan.eq.1.or.snodepobs(source)%gridspan.eq.2) then 
           gridDesco(4) = LVT_rc%gridDesc(4)
           gridDesco(7) = LVT_rc%gridDesc(7)
        else
           if(ihemi.eq.1) then 
              gridDesco(4) = LVT_rc%gridDesc(9)/2
              gridDesco(7) = LVT_rc%gridDesc(7)
           else
              gridDesco(4) = LVT_rc%gridDesc(4)
              gridDesco(7) = -LVT_rc%gridDesc(9)/2
           endif
        endif
        if(ihemi.eq.1) then 
           call bilinear_interp(gridDesco,li,gi_temp,&
                lo_nh,go_nh,snodepobs(source)%mi,snodepobs(source)%mo1,&
                snodepobs(source)%rlat1_nh,snodepobs(source)%rlon1_nh,&
                snodepobs(source)%w111_nh,snodepobs(source)%w121_nh,&
                snodepobs(source)%w211_nh,snodepobs(source)%w221_nh,&
                snodepobs(source)%n111_nh,snodepobs(source)%n121_nh,&
                snodepobs(source)%n211_nh,snodepobs(source)%n221_nh,LVT_rc%udef,&
                iret)
        elseif(ihemi.eq.2) then 
           call bilinear_interp(gridDesco,li,gi_temp,&
                lo_sh,go_sh,snodepobs(source)%mi,snodepobs(source)%mo2,&
                snodepobs(source)%rlat1_sh, snodepobs(source)%rlon1_sh,&
                snodepobs(source)%w111_sh,snodepobs(source)%w121_sh,&
                snodepobs(source)%w211_sh,snodepobs(source)%w221_sh,&
                snodepobs(source)%n111_sh,snodepobs(source)%n121_sh,&
                snodepobs(source)%n211_sh,snodepobs(source)%n221_sh,&
                LVT_rc%udef,iret)
        endif
     end do
     if(snodepobs(source)%gridspan.eq.1) then
        comb_data(1:snodepobs(source)%mo2) = LVT_rc%udef
        comb_data(snodepobs(source)%mo2+1:snodepobs(source)%mo1+snodepobs(source)%mo2) = go_nh(:)
     elseif(snodepobs(source)%gridspan.eq.2) then 
        comb_data(1:snodepobs(source)%mo2) = go_sh(:)
        comb_data(snodepobs(source)%mo2+1:snodepobs(source)%mo1+snodepobs(source)%mo2) = LVT_rc%udef
     else
        comb_data(1:snodepobs(source)%mo2) = go_sh(:)
        comb_data(snodepobs(source)%mo2+1:snodepobs(source)%mo1+snodepobs(source)%mo2) = go_nh(:)
     endif
     varfield = RESHAPE(comb_data(1:snodepobs(source)%mo1+snodepobs(source)%mo2),(/LVT_rc%lnc,LVT_rc%lnr/))
     
     if(snodepobs(source)%gridspan.eq.1) then 
        deallocate(lo_nh)
        deallocate(go_nh)
     elseif(snodepobs(source)%gridspan.eq.2) then 
        deallocate(lo_sh)
        deallocate(go_sh)
     else
        deallocate(lo_nh)
        deallocate(lo_sh)
        deallocate(go_nh)
        deallocate(go_sh)
     endif
     deallocate(comb_data)     
  endif

end subroutine interp_SNODEPfield
