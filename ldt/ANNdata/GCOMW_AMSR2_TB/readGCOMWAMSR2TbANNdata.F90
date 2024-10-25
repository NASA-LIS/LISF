!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LDT_misc.h"
!BOP
! 
! !ROUTINE: readGCOMWAMSR2TbANNdata
! \label{readGCOMWAMSR2TbANNdata}
! 
! !REVISION HISTORY: 
!  20 Jan 2018: Sujay Kumar, Initial Specification
! 
! !INTERFACE: 
subroutine readGCOMWAMSR2TBANNdata(n,iomode, p_s, p_e)
! !USES:   
  use ESMF
  use LDT_coreMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN
  use LDT_ANNMod
  use GCOMWAMSR2TB_ANNdataMod
  use map_utils

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n
  integer, intent(in) :: iomode
  integer, intent(in) :: p_s
  integer, intent(in) :: p_e
! 
! !DESCRIPTION: 
! 
! This subroutine provides the data reader for the synthetic
! soil moisture retrieval product. 
!
!EOP

  real              :: timenow
  logical           :: alarmCheck
  logical           :: file_exists
  integer           :: c,r,i,j,k
  integer           :: ftn
  character(len=LDT_CONST_PATH_LEN) :: fnames(2) ! 2 channels
  character*3       :: fnest
  real              :: TB_H(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: TB_V(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real              :: dtb_ah(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real              :: dtb_dh(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real              :: dtb_av(LDT_rc%lnc(n), LDT_rc%lnr(n))
  real              :: dtb_dv(LDT_rc%lnc(n), LDT_rc%lnr(n))


  dtb_ah = -9999.0
  dtb_dh = -9999.0
  dtb_av = -9999.0
  dtb_dv = -9999.0

  write(fnest,'(i3.3)') n
  alarmCheck = LDT_isAlarmRinging(LDT_rc,"GCOMW AMSR2 Tb alarm "//trim(fnest))

  if(alarmCheck) then
     GCOMWAMSR2TBobs(n)%TB_A_H = LDT_rc%udef
     GCOMWAMSR2TBobs(n)%TB_A_V = LDT_rc%udef

     GCOMWAMSR2TBobs(n)%TB_D_H = LDT_rc%udef
     GCOMWAMSR2TBobs(n)%TB_D_V = LDT_rc%udef
        
     GCOMWAMSR2TBobs(n)%startmode = .false. 
     
     call create_ANNGCOMWAMSR2TB_filenames(GCOMWAMSR2TBobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, 'A',fnames)
     
     do k=1,2
        inquire(file=trim(fnames(k)),exist=file_exists)
        if(file_exists) then
           write(LDT_logunit,*) '[INFO] Reading ..',trim(fnames(k))

           TB_H = LDT_rc%udef
           TB_V = LDT_rc%udef
           call readGCOMWAMSR2Tb_SC_data(n, fnames(k),TB_H, TB_V)
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if(TB_H(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
                    GCOMWAMSR2TBobs(n)%TB_A_H(c,r,k) = TB_H(c+(r-1)*LDT_rc%lnc(n))
                 endif
                 if(TB_V(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
                    GCOMWAMSR2TBobs(n)%TB_A_V(c,r,k) = TB_V(c+(r-1)*LDT_rc%lnc(n))
                 endif
              enddo
           enddo
        endif
     enddo

     call create_ANNGCOMWAMSR2TB_filenames(GCOMWAMSR2TBobs(n)%odir, &
          LDT_rc%yr, LDT_rc%mo, LDT_rc%da, 'D',fnames)
     
     do k=1,2
        inquire(file=trim(fnames(k)),exist=file_exists)
        if(file_exists) then
           write(LDT_logunit,*) '[INFO] Reading ..',trim(fnames(k))

           TB_H = LDT_rc%udef
           TB_V = LDT_rc%udef
           call readGCOMWAMSR2Tb_SC_data(n, fnames(k),TB_H, TB_V)
           
           do r=1,LDT_rc%lnr(n)
              do c=1,LDT_rc%lnc(n)
                 if(TB_H(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
                    GCOMWAMSR2TBobs(n)%TB_D_H(c,r,k) = TB_H(c+(r-1)*LDT_rc%lnc(n))
                 endif
                 if(TB_V(c+(r-1)*LDT_rc%lnc(n)).ne.-9999.0) then 
                    GCOMWAMSR2TBobs(n)%TB_D_V(c,r,k) = TB_V(c+(r-1)*LDT_rc%lnc(n))
                 endif
              enddo
           enddo
        endif
     enddo

#if 0 
     do k=1,2
        call LDT_logSingleANNdata(n,&
             GCOMWAMSR2TBobs(n)%TB_A_H(:,:,k),  &
             pindex=k+p_s-1, &
             iomode = iomode, &
             name = GCOMWAMSR2TBobs(n)%TB_A_H_name(k),&
             units="K")
     enddo

     do k=3,4
        call LDT_logSingleANNdata(n,&
             GCOMWAMSR2TBobs(n)%TB_A_V(:,:,k-2),  &
             pindex=k+p_s-1, &
             iomode = iomode, &
             name = GCOMWAMSR2TBobs(n)%TB_A_V_name(k-2),&
             units="K")
     enddo

     do k=5,6
        call LDT_logSingleANNdata(n,&
             GCOMWAMSR2TBobs(n)%TB_D_H(:,:,k-4),  &
             pindex=k+p_s-1, &
             iomode = iomode, &
             name = GCOMWAMSR2TBobs(n)%TB_A_H_name(k-4),&
             units="K")
     enddo

     do k=7,8
        call LDT_logSingleANNdata(n,&
             GCOMWAMSR2TBobs(n)%TB_D_V(:,:,k-6),  &
             pindex=k+p_s-1, &
             iomode = iomode, &
             name = GCOMWAMSR2TBobs(n)%TB_A_H_name(k-6),&
             units="K")
     enddo
#endif

     do r=1,LDT_rc%lnr(n)
        do c=1,LDT_rc%lnc(n)
           if(GCOMWAMSR2TBobs(n)%TB_A_H(c,r,1).ne.-9999.0.and.&
                GCOMWAMSR2TBobs(n)%TB_A_H(c,r,2).ne.-9999.0) then 
              dtb_ah(c,r) = GCOMWAMSR2TBobs(n)%TB_A_H(c,r,1) - &
                   GCOMWAMSR2TBobs(n)%TB_A_H(c,r,2)              
           endif
           if(GCOMWAMSR2TBobs(n)%TB_D_H(c,r,1).ne.-9999.0.and.&
                GCOMWAMSR2TBobs(n)%TB_D_H(c,r,2).ne.-9999.0) then 
              dtb_dh(c,r) = GCOMWAMSR2TBobs(n)%TB_D_H(c,r,1) - &
                   GCOMWAMSR2TBobs(n)%TB_D_H(c,r,2)              
           endif
           if(GCOMWAMSR2TBobs(n)%TB_A_V(c,r,1).ne.-9999.0.and.&
                GCOMWAMSR2TBobs(n)%TB_A_V(c,r,2).ne.-9999.0) then 
              dtb_av(c,r) = GCOMWAMSR2TBobs(n)%TB_A_V(c,r,1) - &
                   GCOMWAMSR2TBobs(n)%TB_A_V(c,r,2)              
           endif
           if(GCOMWAMSR2TBobs(n)%TB_D_V(c,r,1).ne.-9999.0.and.&
                GCOMWAMSR2TBobs(n)%TB_D_V(c,r,2).ne.-9999.0) then 
              dtb_dv(c,r) = GCOMWAMSR2TBobs(n)%TB_D_V(c,r,1) - &
                   GCOMWAMSR2TBobs(n)%TB_D_V(c,r,2)              
           endif
        enddo
     enddo
     call LDT_logSingleANNdata(n,&
          DTB_AH,&
          pindex=p_s, &
          iomode = iomode, &
          name = "DTB_A_H",&
          units="K")
     call LDT_logSingleANNdata(n,&
          DTB_DH,&
          pindex=1+p_s, &
          iomode = iomode, &
          name = "DTB_D_H",&
          units="K")

     call LDT_logSingleANNdata(n,&
          DTB_DH,&
          pindex=2+p_s, &
          iomode = iomode, &
          name = "DTB_A_V",&
          units="K")
     call LDT_logSingleANNdata(n,&
          DTB_DH,&
          pindex=3+p_s, &
          iomode = iomode, &
          name = "DTB_D_V",&
          units="K")


  endif

end subroutine readGCOMWAMSR2TBANNdata

subroutine readGCOMWAMSR2Tb_SC_data(n, fname,TB_H, TB_V)

  use LDT_coreMod
  use LDT_logMod
  use GCOMWAMSR2TB_ANNdataMod

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none
  
  integer              :: n 
  character(len=*)     :: fname
  real                 :: TB_H(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  real                 :: TB_V(LDT_rc%lnc(n)*LDT_rc%lnr(n))
 
  integer              :: ftn,status
  integer              :: c,r,c1,r1,t
  integer              :: tb_h_id, tb_v_id

  real                 :: tb_h_in(GCOMWAMSR2TBobs(n)%nc, GCOMWAMSR2TBobs(n)%nr)
  real                 :: tb_v_in(GCOMWAMSR2TBobs(n)%nc, GCOMWAMSR2TBobs(n)%nr)

  logical*1            :: tb_h_b(GCOMWAMSR2TBobs(n)%nc*GCOMWAMSR2TBobs(n)%nr)
  real                 :: tb_h_1d(GCOMWAMSR2TBobs(n)%nc*GCOMWAMSR2TBobs(n)%nr)
  logical*1            :: tb_h_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))

  logical*1            :: tb_v_b(GCOMWAMSR2TBobs(n)%nc*GCOMWAMSR2TBobs(n)%nr)
  real                 :: tb_v_1d(GCOMWAMSR2TBobs(n)%nc*GCOMWAMSR2TBobs(n)%nr)
  logical*1            :: tb_v_b_ip(LDT_rc%lnc(n)*LDT_rc%lnr(n))
  integer              :: ios
  

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

  ios = nf90_open(path=trim(fname),mode=NF90_NOWRITE, ncid=ftn)
  call LDT_verify(ios, 'Error opening file '//trim(fname))
  
  ios = nf90_inq_varid(ftn,"Brightness Temperature (H)",tb_h_id)
  call LDT_verify(ios,"Error in nf90_inq_varid: Tb(H)")

  ios = nf90_inq_varid(ftn,"Brightness Temperature (V)",tb_v_id)
  call LDT_verify(ios,"Error in nf90_inq_varid: Tb(V)")

  ios = nf90_get_var(ftn,tb_h_id, TB_H_IN)
  call LDT_verify(ios,"Error in nf90_get_var: TB_H")

  ios = nf90_get_var(ftn,tb_v_id, TB_V_IN)
  call LDT_verify(ios,"Error in nf90_get_var: TB_V")

  ios = nf90_close(ftn)
  call LDT_verify(ios, "Error in closing file "//trim(fname))

  tb_h_b = .false. 
  tb_v_b = .false. 

  tb_h_1d = LDT_rc%udef
  tb_v_1d = LDT_rc%udef

  do r=1,GCOMWAMSR2TBobs(n)%nr
     do c=1,GCOMWAMSR2TBobs(n)%nc        
        r1 = GCOMWAMSR2TBobs(n)%nr-r+1

        c1 = c+ GCOMWAMSR2TBobs(n)%nc/2 -1
        if(c1.gt.GCOMWAMSR2TBobs(n)%nc) then 
           c1 = c-GCOMWAMSR2TBobs(n)%nc/2 - 1
        endif
        if(tb_h_in(c1,r1).ne.65534) then 
           tb_h_1d(c+(r-1)*GCOMWAMSR2TBobs(n)%nc) = & 
                tb_h_in(c1,r1)*0.01
           tb_v_1d(c+(r-1)*GCOMWAMSR2TBobs(n)%nc) = & 
                tb_v_in(c1,r1)*0.01

           tb_h_b(c+(r-1)*GCOMWAMSR2TBobs(n)%nc) = .true. 
           tb_v_b(c+(r-1)*GCOMWAMSR2TBobs(n)%nc) = .true. 
        endif
     enddo
  enddo

        
!  open(100,file='test_in.bin',form='unformatted')
!  write(100) tb_h_1d
!  write(100) tb_v_1d
!  close(100)
  
!--------------------------------------------------------------------------
! Interpolate to the LDT running domain
!-------------------------------------------------------------------------- 
  call neighbor_interp(LDT_rc%gridDesc(n,:),&
       tb_h_b, tb_h_1d, tb_h_b_ip, tb_h, &
       GCOMWAMSR2TBobs(n)%nc*GCOMWAMSR2TBobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       GCOMWAMSR2TBobs(n)%n11,&
       LDT_rc%udef, status)

  call neighbor_interp(LDT_rc%gridDesc(n,:),&
       tb_v_b, tb_v_1d, tb_v_b_ip, tb_v, &
       GCOMWAMSR2TBobs(n)%nc*GCOMWAMSR2TBobs(n)%nr, &
       LDT_rc%lnc(n)*LDT_rc%lnr(n), &
       LDT_domain(n)%lat, LDT_domain(n)%lon,&
       GCOMWAMSR2TBobs(n)%n11,&
       LDT_rc%udef, status)

!  open(100,file='test_out.bin',form='unformatted')
!  write(100) tb_h
!  close(100)
!  stop
#endif

end subroutine readGCOMWAMSR2Tb_SC_data
!BOP
! !ROUTINE: create_ANNGCOMWAMSR2TB_filename
! \label{create_ANNGCOMWAMSR2TB_filename}
! 
! !INTERFACE: 
subroutine create_ANNGCOMWAMSR2TB_filenames(ndir, yr, mo,da, pass, filenames)
! !USES:   

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: filenames(2)
  integer           :: yr, mo, da
  character (len=*) :: pass
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
!  This subroutine creates the synthetic filename based on the time and date 
! 
!  AMSR2 filename convention: GW1AM2_YYYYMMDD_ttt_PPWX_LLxxKKKrdvaaappp
!   GW1 : Satellite(Fixed Value) AM2 : Sensor(Fixed Value)
!   YYYYMMDD : Start date of observation (UTC), "DD"="00" means monthly data
!   ttt : Calculating period (01D:Daily, 01M:Monthly)
!   PP : Map projection (EQ: Equirectangular, PN: Polar stereo(North hemisphere), PS: Polar stereo(South hemisphere))
!   W : Calculating method (M: Mean, O: Overwrite) X : Orbit direction(A: Ascending, D: Descending)
!   LL : Process Level (L3:Level 3) xx : Process Kind(SG: Standard operation product)
!   KKK : Product ID
!
!  The arguments are: 
!  \begin{description}
!  \item[ndir] name of the synthetic soil moisture directory
!  \item[yr]  current year
!  \item[mo]  current month
!  \item[da]  current day
!  \item[filenames] Generated synthetic filename
! \end{description}
!EOP

  character (len=4) :: fyr
  character (len=2) :: fmo,fda
  
  write(unit=fyr, fmt='(i4.4)') yr
  write(unit=fmo, fmt='(i2.2)') mo
  write(unit=fda, fmt='(i2.2)') da
 
  filenames(1) = trim(ndir)//'/TB18GHz_25/'//trim(fyr)//'/GW1AM2_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'_01D_EQM'//trim(pass)//&
       '_L3SGT18LA2220220.h5'
  filenames(2) = trim(ndir)//'/TB36GHz_25/'//trim(fyr)//'/GW1AM2_' &
       //trim(fyr)//trim(fmo)//trim(fda)//'_01D_EQM'//trim(pass)//&
       '_L3SGT36LA2220220.h5'

      
end subroutine create_ANNGCOMWAMSR2TB_filenames
