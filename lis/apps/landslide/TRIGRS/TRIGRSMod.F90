!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module TRIGRSMod
!BOP
!
! !MODULE:
! 
! !DESCRIPTION:
!  Module for program to compute pore-pressure response and
!  factor of safety for saturated and unsaturated infiltration.
!  by Rex L. Baum and W.Z. Savage, USGS
! 
! !REVISION HISTORY: 
! ! 27 April 2013: Soni Yatheendradas, Recoded from earlier bugged version to current current rain-looping one, all physics looping in TRIGRS_final; native code TRIGRS 2.0.09z 
  use ESMF

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: TRIGRS_init
  PUBLIC :: TRIGRS_run
  PUBLIC :: TRIGRS_output
  PUBLIC :: TRIGRS_final
!EOP

contains

!BOP
! 
! !ROUTINE: TRIGRS_init
! \label{TRIGRS_init}
! 
! !INTERFACE:
  subroutine TRIGRS_init()
! !USES:   
    use LIS_coreMod, only : LIS_rc, LIS_config
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    use LIS_constantsMod,  only : LIS_CONST_PI
    use LIS_timeMgrMod, only : LIS_parseTimeString, &
         LIS_update_timestep, LIS_registerAlarm ! SY

    use input_file_defs
    use input_vars
    use grids
    use model_vars

!
! !DESCRIPTION:
!EOP
!   SY: Subroutine indentation performed for integration into LIS

    implicit none

    character*100       :: trigrs_initfile
    integer             :: rc
    integer             :: n

    integer:: i,j,k!,imx1,mnd !,m ! "m" removed 1 Feb 2013, RLB ! SY
    integer:: patlen !added umax, 12/6/2010, RLB ! SY
    integer:: maxzo ! SY
    real::mnzmx,mndep ! SY
    character (len=1):: tb
    character (len=255):: outfil,infil
    character (len=8):: date
    character (len=10):: time
    character (len=7):: vrsn
    character (len=11):: bldate
    character (len=2)::pid(3)
    integer:: ftn ! SY

! first executable statement ............
    call date_and_time(date,time)
    test=-9999.D0; test1=-9999.
    pid=(/'TI','GM','TR'/)
    pi=3.141592653589793
    dg2rad=pi/180.D0
    vrsn='2.0.09z'; bldate='06 May 2013'
    smt=0.01d0; lard=12.d0 ! test values for early-time (moved 2/15/2012, RLB)
    write (*,*) ''
    write (*,*) 'TRIGRS: Transient Rainfall Infiltration'
    write (*,*) 'and Grid-based Regional Slope-Stability'
    write (*,*) '               Analysis'
    write (*,*) '       Version ', vrsn,', ',bldate
    write (*,*) '  By Rex L. Baum and William Z. Savage'
    write (*,*) '       U.S. Geological Survey'
    write (*,*) '-----------------------------------------'
    write (*,*) 'Portions of this program include material'
    write (*,*) '           copyrighted (C) by'
    write (*,*) '      Absoft Corporation 1988-2010.'
    write (*,*) ''
    tb=char(9)
!  open log file ! 14 Feb 2013 Added adjustl() statements for compatability with other compilers
    !outfil='TrigrsLog.txt'; outfil=adjustl(outfil) !SY
    !open (u(19),file=trim(outfil),status='unknown',err=410) ! SY
    write (LIS_logunit,*) '' ! SY
    write (LIS_logunit,*) 'Starting TRIGRS ', vrsn,' ',bldate ! SY
    write (LIS_logunit,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4) ! SY
    write (LIS_logunit,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6) ! SY

    call ESMF_ConfigGetAttribute(LIS_config,trigrs_initfile,&
         label="TRIGRS initialization file:",rc=rc)
    call LIS_verify(rc,'TRIGRS initialization file: option not specified in the config file')

    ! SY: Start for obtaining intended TRIGRS model time step
    call ESMF_ConfigFindLabel(LIS_config,"TRIGRS app timestep:",rc=rc)
    do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,time,rc=rc)
     call LIS_verify(rc,'TRIGRS app timestep: not defined')

     call LIS_parseTimeString(time,TRIGRS_timestep)

     call LIS_update_timestep(LIS_rc, n, TRIGRS_timestep)
  
     call LIS_registerAlarm("TRIGRS app alarm",&
          TRIGRS_timestep,&
          TRIGRS_timestep)
     
    enddo
    ! SY: End for obtaining intended TRIGRS model time step

!read initialization file
    call trini(trigrs_initfile,dg2rad)

! SY: Begin reading choice of forcings source: LIS interface or native TRIGRS input files-style
    call ESMF_ConfigGetAttribute(LIS_config,RainReadSource,&
         label="TRIGRS rain input source style:",rc=rc)
    call LIS_verify(rc,'TRIGRS rain input source style: option not specified in the config file')
! SY: End reading choice of forcings source: LIS interface or native TRIGRS input files-style

! determine grid size parameters RLB 4/18/2011
    patlen=scan(elevfil,'/\',.true.) ! find end of folder name
    elfoldr=elevfil(1:patlen) ! path to elevation grid
    ans=.false.
    do i=1,3
      elfoldr=adjustl(elfoldr)
      infil=trim(elfoldr)//pid(i)//'grid_size.txt'
      infil=adjustl(infil)
      inquire (file=trim(infil),exist=ans)
      write(*,*) trim(infil), ans
      if(ans) exit
    end do
    if(ans) then
      ftn = LIS_getNextUnitNumber() ! SY 
      !open (u(22),file=trim(infil),status='unknown',err=420) ! SY
      open (ftn,file=trim(infil),status='unknown',err=420) ! SY
      read (ftn,*) heading ! SY
      read (ftn,*) imax,row,col,nwf ! SY
      close (ftn) ! SY
      call LIS_releaseUnitNumber(ftn) ! SY
    else
      infil=elevfil; infil=adjustl(infil)
      call ssizgrd(row,col,celsiz,nodat,imax,infil,&
        &trigrs_header) ! SY
      outfil=trim(elfoldr)//pid(3)//'grid_size.txt'
      outfil=adjustl(outfil)
      ftn = LIS_getNextUnitNumber() ! SY 
      open (ftn,file=trim(outfil),status='unknown',err=410) ! SY
      write (ftn,*) 'imax      row      col      nwf' ! SY
      nwf=1 ! dsctr is computed by TopoIndex; dsctr=1 is default value for no runoff routing. ! SY
      write (ftn,*) imax,row,col,nwf ! SY
      write (ftn,*) '' ! SY
      close (ftn) ! SY
      call LIS_releaseUnitNumber(ftn) ! SY
    end if
    write(LIS_logunit,*) 'Grid size parameters from ', trim(infil) ! SY
    write (LIS_logunit,*) heading ! SY
    write (LIS_logunit,*) imax,row,col,nwf ! SY

    ! SY: Begin domain size checking    
    n = 1
    if ((row .NE. LIS_rc%lnr(n)) .OR. &
       & (col .NE. LIS_rc%lnc(n))) then
     write(LIS_logunit,*) &
       & 'Domain size mismatch between lis.config & tr_in.txt'
     write(*,*) &
       & 'Domain size mismatch between lis.config & tr_in.txt'
     call LIS_endrun() ! SY
    end if 
    ! SY: End domain size checking    

! Allocate & initialize arrays needed for runoff routing
    trigrs_grd=row*col ! SY
    trigrs_imx1=imax
    allocate (pf2(trigrs_grd),indx(imax),nxt(imax))
    allocate (dsctr(imax+1),slo(imax))
    allocate (pf1(trigrs_grd),rizero(imax))
    allocate (pf1_Slope(trigrs_grd)) ! SY
    allocate (pf1_Mat(nper,trigrs_grd)) ! SY
    allocate (ri(imax),rik(imax*nper),ro(imax))
    allocate (ri_Mat(nper,imax)) ! SY
    allocate (rikzero(imax),temp(col),itemp(col))
    allocate (depth(imax),zmax(imax))
    allocate (zo(imax),ir(imax),tfg(imax))
    allocate (elev(imax)) ! added 4/21/2010
    elev=0.
    pf2=0
    indx=0
    nxt=0
    dsctr=0
    zo=0
    pf1=0.
    pf1_Slope=0. ! SY
    pf1_Mat=0. ! SY
    slo=0.
    rizero=0.
    ri=0.
    ri_Mat=0. ! SY
    rik=0.
    ro=0.
    rikzero=0.
    temp=0.;itemp=0
    ir=0.
    depth=0.
    zmax=0.
! Choose file extension for grid files, Added 4/14/2010
    grxt='.txt'
    if(lasc) grxt='.asc'
! *****************************************************************
!read gridded data from GIS
    write (*,*) 'Reading input grids'
    write(LIS_logunit,*) 'Input file name,            Cell count' ! SY
!  read slope angles
    call srdgrd(trigrs_grd,col,&
      &   trigrs_ncol,trigrs_nrow,&
      &   celsiz,nodat,&
      &   slo,pf1,trigrs_sctr,imax,temp,slofil,&
      &   param,trigrs_header) ! SY
    ! SY: Begin for assisting rain assignment from 2-d LIS grid to 1-d TRIGRS array
    pf1_Slope=pf1
    nodat_Slope=nodat
    ! SY: End for assisting rain assignment from 2-d LIS grid to 1-d TRIGRS array
    write(LIS_logunit,*) 'Slope angle grid' ! SY
    write(LIS_logunit,*) trim(slofil),trigrs_sctr,&
      &  ' data cells' ! SY
    if(trigrs_sctr/=imax &
      & .or. trigrs_ncol/=col &
      & .or. trigrs_nrow/=row) then
     write(*,*) 'Grid mismatch: ', trim(slofil)
     write(*,*) 'Check slope grid and (or) initialization file against elevation grid.'
     write(LIS_logunit,*) 'Grid mismatch: ', trim(slofil) ! SY
     write(LIS_logunit,*) 'Check slope grid and (or) initialization file against elevation grid.' ! SY
    end if
    slo=slo*dg2rad ! convert slope angles to radians
!  read property zone numbers, zo
    if(nzon==1) then
      zo=1 ! if only one zone, all values of zone grid equal 1.
      write(*,*) 'One property zone, no grid required!'
      write(LIS_logunit,*) 'One property zone, no grid required!' ! SY
      parami=param ! added 7/29/2008 RLB
    else
      call irdgrd(trigrs_grd,col, &
       &trigrs_ncol,trigrs_nrow,celsiz,&
       &trigrs_nodata,trigrs_mnd,&
       !&zo,pf2,trigrs_sctr,imax,temp,u(15),zonfil,& ! SY : see next line
       &zo,pf2,trigrs_sctr,imax,& 
       &itemp,zonfil,& 
       &parami,trigrs_header) ! SY
      write(LIS_logunit,*) 'Property zone grid' ! SY
      write(LIS_logunit,*) trim(zonfil),trigrs_sctr,&
        &' data cells' ! SY
      if(trigrs_sctr/=imax .or. &
        & trigrs_ncol/=col .or. &
        & trigrs_nrow/=row) then
        write (*,*) 'Grid mismatch ',trim(zonfil)
        write (*,*) 'Correct property-zone grid and/or initializtion file.'
        write (LIS_logunit,*) 'Grid mismatch ',trim(zonfil) ! SY
        write (LIS_logunit,*) 'Correct property-zone grid and/or initializtion file.' ! SY
        !close(LIS_logunit) ! SY: DO NOT close log file
        !pause 'Press return/enter to quit' ! SY
        write (LIS_logunit,*) '-1' ! SY
        write (*,*) '-1' ! SY
        !stop '-1' ! SY
        call LIS_endrun() ! SY
      end if
      maxzo=maxval(zo)
! SY: Begin possibly for parallelized domain
!      if (maxzo/=nzon) then
      if (maxzo>nzon) then ! SY
!        write (*,*) 'Maximum zone number does not equal number of property zones!'
        write (*,*) 'Maximum zone number greater than number of property zones!'
        write (*,*) 'Correct property-zone grid and/or initializtion file.'
!        write (LIS_logunit,*) 'Maximum zone number does not equal number of property zones!'
        write (LIS_logunit,*) 'Maximum zone number greater than number of property zones!'
        write (LIS_logunit,*) 'Correct property-zone grid and/or initializtion file.'
        !close(LIS_logunit) ! SY: DO NOT close log file
        !pause 'Press return/enter to quit' ! SY
        write (LIS_logunit,*) '-1' ! SY
        write (*,*) '-1' ! SY
        !stop '-1' ! SY
        call LIS_endrun() ! SY
      end if
! SY: End possibly for parallelized domain
    end if
! *********************
!  read background infiltration rate, Isub0
    if (crizero.lt.0) then
     call srdgrd(trigrs_grd,col,&
       &trigrs_ncol,trigrs_nrow,&
       &celsiz,nodat,&
       &rizero,pf1,trigrs_sctr,imax,temp,&
       &rizerofil,param,trigrs_header) ! SY
     write(LIS_logunit,*) 'Background infiltration rate grid'
     write(LIS_logunit,*) trim(rizerofil),trigrs_sctr,&
       &' data cells'
     if(trigrs_sctr/=imax .or. &
       &trigrs_ncol/=col .or. &
       & trigrs_nrow/=row)&
       & write (LIS_logunit,*) 'Grid mismatch ',trim(rizerofil)
    else
     rizero=crizero
    end if
!  read initial depth to water table,
    if (dep.lt.0) then
     call srdgrd(trigrs_grd,col,trigrs_ncol,&
       &trigrs_nrow,celsiz,nodat,&
       &depth,pf1,trigrs_sctr,&
       &imax,temp,depfil,param,trigrs_header) ! SY
     write(LIS_logunit,*) 'Initial water-table depth grid'
     write(LIS_logunit,*) trim(depfil),trigrs_sctr,&
       &' data cells'
     if(trigrs_sctr/=imax .or. &
       &trigrs_ncol/=col .or. &
       &trigrs_nrow/=row)&
       &write (LIS_logunit,*) 'Grid mismatch ',trim(depfil)
    else
     depth=dep
    end if
!  read depth to base of potential slide, zmax
    if (czmax.lt.0) then
     call srdgrd(trigrs_grd,col,trigrs_ncol,&
       &trigrs_nrow,celsiz,nodat,&
       &zmax,pf1,trigrs_sctr,imax,temp,&
       &zfil,param,trigrs_header) ! SY
     write(LIS_logunit,*) 'Maximum depth grid'
     write(LIS_logunit,*) trim(zfil),trigrs_sctr,&
       &' data cells'
     if(trigrs_sctr/=imax .or. &
       &trigrs_ncol/=col .or. &
       &trigrs_nrow/=row) &
       &write (LIS_logunit,*) 'Grid mismatch ',trim(zfil)
    else
     zmax=czmax
    end if
! Trap error conditions for zmin values.
    mndep=minval(depth) !Added 28 Jan 2013, RLB
    mnzmx=minval(zmax) !Added 28 Jan 2013, RLB
    if(zmin>mnzmx .or. zmin>mndep) zmin=0.
!  read digital elevations, elev !! added 2/24/ 2010, unit number changed to 12 12/6/2010
    call srdgrd(trigrs_grd,col,trigrs_ncol,&
      &trigrs_nrow,celsiz,nodat,&
      &  elev,pf1,trigrs_sctr,imax,temp,elevfil,&
      &  param,trigrs_header) ! SY
    write(LIS_logunit,*) 'Elevation grid'
    write(LIS_logunit,*) trim(elevfil),trigrs_sctr,&
      &' data cells'
    if(trigrs_sctr/=imax .or. &
      &trigrs_ncol/=col .or. &
      &trigrs_nrow/=row) &
      &write (LIS_logunit,*) 'Grid mismatch ',trim(zfil)
! *****************************************************************
    write(LIS_logunit,*) '---------------******---------------'
! test and adjust (if necessary) steady background infiltration rates
    call steady(sumex,trigrs_imx1) ! SY

    i_per = 0 ! SY

    return
! Error reporting
410 continue
    write (*,*) 'Error opening output file'
    write (*,*) '--> ',outfil
    write (*,*) 'Check file path and status'
    write (LIS_logunit,*) 'Error opening output file'
    write (LIS_logunit,*) '--> ',outfil
    write (LIS_logunit,*) 'Check file path and status'
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '410' ! SY
    write (*,*) '410' ! SY
    !stop '410'  ! SY
    call LIS_endrun() ! SY
420 continue
    write (*,*) 'Error opening input file'
    write (*,*) '--> ',infil
    write (*,*) 'Check file path and status'
    write (LIS_logunit,*) 'Error opening input file'
    write (LIS_logunit,*) '--> ',infil
    write (LIS_logunit,*) 'Check file path and status'
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '420' ! SY
    write (*,*) '420' ! SY
    !stop '420'  ! SY
    call LIS_endrun() ! SY

  end subroutine TRIGRS_init

!BOP
! 
! !ROUTINE: TRIGRS_run
! \label{TRIGRS_run}
! 
! !INTERFACE:
  subroutine TRIGRS_run(n)
! !USES:
    use LIS_coreMod,  only : LIS_rc, LIS_domain
    use LIS_logMod,         only : LIS_endrun, LIS_verify, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use LIS_metforcingMod, only : LIS_FORC_State
    use LIS_FORC_AttributesMod
    use LIS_timeMgrMod,    only : LIS_isAlarmRinging, LIS_is_last_step ! SY

    use input_file_defs
    use input_vars
    use grids
    use model_vars
!
! !DESCRIPTION:
!EOP
!   SY: Subroutine indentation performed for integration into LIS
    implicit none

    integer, intent(IN) :: n
    integer           :: cc,rr
    type(ESMF_Field)  :: pcpField
    real,     pointer :: pcp(:)
    integer           :: status

    integer:: mm,ii ! SY
    logical             :: alarmCheck ! SY 
    logical :: IsLastStep ! SY
    integer:: i,j,k!,imx1,mnd !,m ! "m" removed 1 Feb 2013, RLB ! SY
    !integer:: ncol,nrow,u(ulen),maxzo,ncc,nccs ! SY
    integer:: ncc,nccs ! SY
    !real::x1,mnzmx,mndep ! SY
    real::x1 !,mnzmx,mndep ! SY
    !character (len=1):: tb ! SY
    character (len=255):: outfil,infil
    character (len=14):: fminfil='TRfs_min_'
    character (len=14):: zfminfil='TRz_at_fs_min_'
    character (len=14):: pminfil='TRp_at_fs_min_'
    character (len=8):: wtabfil='TRwater_'
    character (len=18):: profil='TRlist_z_p_fs_'
    !character (len=14):: header(6) ! SY
    character (len=13):: ncvfil='TRnon_convrg_'
    character (len=8):: date
    character (len=10):: time
    character (len=4):: stp
    character (len=31):: scratch,irfil
    character (len=7):: vrsn
    !character (len=11):: bldate ! SY
    !character (len=2)::pid(3) ! SY
    logical :: lwarn
    integer:: ftn ! SY
    logical :: unit_open ! SY

    alarmCheck = LIS_isAlarmRinging(LIS_rc,"TRIGRS app alarm") ! SY
    if(alarmCheck) then

     call ESMF_StateGet(LIS_FORC_State(n),trim(LIS_FORC_Rainf%varname(1)),pcpField,&
          rc=status)
     call LIS_verify(status,'TRIGRS_run: error getting Rainf')
 
     call ESMF_FieldGet(pcpField,localDE=0, farrayPtr=pcp,rc=status)
     call LIS_verify(status,'TRIGRS_run: error retrieving pcp')

     i_per=i_per+1 ! SY
 
     ! SY: Begin check if LIS timestep beyond that in tr_in.txt
     if (i_per .GT. nper) then
      write (LIS_logunit,*) 'i_per became > nper!!!' ! SY
      write (*,*) 'i_per became > nper!!!' ! SY
      call LIS_endrun() ! SY
     endif
     ! SY: End check if LIS timestep beyond that in tr_in.txt
 
     j=i_per ! SY
 
 !  read precipitation intensity, I
 ! SY: Begin extra code
     if (RainReadSource .EQ. 1) then ! SY: Forcings coming through LIS interface
 
      do rr=1,LIS_rc%lnr(n)
       do cc=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(cc,rr).ne.-1) then
         !pf1(cc+(rr-1)*LIS_rc%lnc(n))=& ! SY: Needs to be FLIPPED, see next line
         pf1(cc+(LIS_rc%lnr(n)-rr)*LIS_rc%lnc(n))=&
           pcp(LIS_domain(n)%gindex(cc,rr))/1000.0 !convert to m/s
         !print*,'blk1',cc,rr,ri(cc+(rr-1)*LIS_rc%lnc(n)) ! SY
        else
         !pf1(cc+(rr-1)*LIS_rc%lnc(n))= LIS_rc%udef ! SY: Needs to be FLIPPED, see next line
         pf1(cc+(LIS_rc%lnr(n)-rr)*LIS_rc%lnc(n))= LIS_rc%udef
        endif
       enddo
      end do
     ! SY: Begin from srdgrd, guided by nodata values of slofil
     trigrs_sctr=0
     do mm=1,trigrs_nrow
      do ii=1,trigrs_ncol
       temp(ii)=pf1_Slope(ii+(mm-1)*trigrs_ncol)
       if(temp(ii).ne.nodat_Slope) then
        trigrs_sctr=trigrs_sctr+1
        if (trigrs_sctr>imax) then
         write (LIS_logunit,*) 'trigrs_sctr>imax in TRIGRS_run'
         write (*,*) 'trigrs_sctr>imax in TRIGRS_run'
         call LIS_endrun()
        endif
        if (pf1(ii+(mm-1)*trigrs_ncol) .ne. LIS_rc%udef) then ! SY: make sure valid TRIGRS pixels have valid rain values
         ri(trigrs_sctr)=pf1(ii+(mm-1)*trigrs_ncol)
         temp(ii)=pf1(ii+(mm-1)*trigrs_ncol)
        else
         write (LIS_logunit,*) 'valid TRIGRS pixel not having valid rain value: TRIGRS row, col are ', mm, ii
         write (*,*) 'valid TRIGRS pixel not having valid rain value: TRIGRS row, col are ', mm, ii
         call LIS_endrun()
        endif
       else ! SY: required to assign nodata values for pixels where LIS sees data
        pf1(ii+(mm-1)*trigrs_ncol)=LIS_rc%udef
       endif
      enddo
     enddo
     ! SY: End from srdgrd, guided by nodata values of slofil
! SY: Begin extra code
     do  i=1,trigrs_grd
      pf1_Mat(j,i)=pf1(i)
     end do
! SY: End extra code
     if(trigrs_sctr/=trigrs_imx1) &
       &write (LIS_logunit,*) 'Grid mismatch ',trim(rifil(j)) ! SY
 
     else if (RainReadSource .EQ. 2) then ! SY: Forcing files that native TRIGRS inputs; for testing purposes primarily
 ! SY: End extra code
 
      if (cri(j).lt.0) then
       call srdgrd(trigrs_grd,col,&
         &trigrs_ncol,trigrs_nrow,&
         &celsiz,nodat,ri,&
         &pf1,trigrs_sctr,imax,temp,&
         &rifil(j),param,trigrs_header) ! SY
 ! SY: Begin extra code
       do  i=1,trigrs_grd
        pf1_Mat(j,i)=pf1(i)
       end do
 ! SY: End extra code
       write(LIS_logunit,*) 'Precipitation intensity grid ',j ! SY
       write(LIS_logunit,*) trim(rifil(j)),trigrs_sctr,&
         &' data cells' ! SY
       if(trigrs_sctr/=trigrs_imx1) &
         &write (LIS_logunit,*) 'Grid mismatch ',trim(rifil(j)) ! SY
      else
       do i=1,trigrs_imx1
        ri(i)=cri(j)
       end do
       ! SY: Rex Baum needs to be contacted about pf1 not being assigned in this case (and hence pf1_Mat not assigned in LIS)
      end if
 
 ! SY: Begin extra code
     else
 
      write (LIS_logunit,*) 'RainReadSource should be 1 or 2' ! SY
      write (*,*) 'RainReadSource should be 1 or 2' ! SY
      call LIS_endrun() ! SY
 
     end if
 ! SY: End extra code
 
 ! SY: Begin extra code
     do  i=1,imax
      ri_Mat(j,i)=ri(i)
     end do
 ! SY: End extra code

! SY: Begin running physics if last time step
     IsLastStep = LIS_is_last_step(LIS_rc)
     if (IsLastStep) then

      ! SY: Begin check if TRIGRS_run was correctly called nper times
      !if (RainReadSource .EQ. 2) then
       if (i_per .NE. nper) then
        write (LIS_logunit,*) 'TRIGRS_run was not called nper times!!!' ! SY
        write (*,*) 'TRIGRS_run was not called nper times!!!' ! SY
        call LIS_endrun() ! SY
       endif
      !endif
      ! SY: End check if TRIGRS_run was correctly called nper times
  
      !pid=(/'TI','GM','TR'/) ! SY
      pi=3.141592653589793
      dg2rad=pi/180.D0
      vrsn='2.0.09z'!; bldate='06 May 2013' ! SY
      smt=0.01d0; lard=12.d0 ! test values for early-time (moved 2/15/2012, RLB)
      fminfil=adjustl(fminfil);zfminfil=adjustl(zfminfil)
      profil=adjustl(profil); pminfil=adjustl(pminfil);
  
  ! conduct runoff routing and adjust transient infiltration rates
      call rnoff(trigrs_grd,sumex,&
      & trigrs_imx1,celsiz,param,parami,nodat,&
      & trigrs_nodata,trigrs_mnd,&
      & trigrs_sctr,trigrs_ncol,&
      & trigrs_nrow,trigrs_header,&
      & test1) ! SY
  ! Deallocate arrays that are no longer needed
      deallocate (ri,ro,wf)
      deallocate (ri_Mat) ! SY
      deallocate (pf2,indx,nxt,dsctr,dsc)
      write(LIS_logunit,*) '---------------******---------------' ! SY
      write(LIS_logunit,*) 'Input file name,          Cell count' ! SY
  ! *****************************************************************
  ! compute pore pressure distributions for either fully saturated or
  ! partially saturated conditions.
  ! Partially saturated zone overlies saturated zone
  ! Allocate and initialize new arrays
      allocate (fsmin(imax*nout),pmin(imax*nout),zfmin(imax*nout))
      allocate (p(nzs+1),ptran(nzs+1),pzero(nzs+1),bline(nzs+1))
      allocate (fc(nzs+1),fw(nzs+1),thz(nzs+1),kz(nzs+1),trz(nzs+1))
      allocate (nvu(imax),nv(imax),uwsp(nzs+1),gs(nzon))
      allocate (chi(nzs+1))
      if(outp(1)) allocate(wtab(imax*nout))
      if (flag<= -4) then ! moved 5/3/2010, flag =-4,...,-9 produces ijz or xyz output 12/22/2010 RLB
       allocate(ix(imax),jy(imax))
       ix=0;jy=0
      end if
      fsmin=0.
      zfmin=0.
      pmin=0.; if(outp(1)) wtab=0.
      p=0.
      ptran=0.
      pzero=0.
      bline=0.
      fc=0.
      fw=0.
      nv=0
      nvu=0
  ! determine number of time steps needed
      kper=nper
  ! SY: Note that the following 'if' part of the decision will not happen because I already enforced in trini that t=capt(nper+1)
      if (t>capt(nper+1)) then
       kper=nper+1
      else
       do k=1,nper ! find the period that contains t
        if(t>=capt(k) .and. t<=capt(k+1)) kper=k
       end do
      end if
      if (tx<1) tx=1
      nts=kper*tx ! number of time-steps from 0 to t
      tns=float(nts)
      tmin=0.
      tmax=t
      tinc=(tmax-tmin)/tns
  ! compute output pointers
      allocate(jsav(nts+1))
      jsav=0
      write (LIS_logunit,*) '******** Output times ********' ! SY
      write (LIS_logunit,*) 'number, timestep #,  time' ! SY
      write (*,*) '******** Output times ********'
      write (*,*) 'number, timestep #,  time'
      lwarn=.false.
      do k=1,nout
       ts=tmin
       do j=1,nts
        if(tsav(k)>=ts .and. tsav(k)<ts+tinc) then
         if(tsav(k)/=ts) lwarn=.true.
         jsav(j)=k
         ksav(k)=j
         tsav(k)=ts
         exit
        else if(tsav(k)>=tmax) then
         jsav(nts+1)=k
         ksav(k)=nts+1
         tsav(k)=tmax
         exit
        end if
        ts=ts+tinc
       end do
       if(lwarn) then
        write(LIS_logunit,*) 'One or more specified output times unavailable, ' ! SY
        write(LIS_logunit,*) 'Nearest available time substituted.' ! SY
        write(*,*) 'One or more specified output times unavailable, '
        write(*,*) 'Nearest available time substituted.'
       end if
       write(LIS_logunit,*) k,ksav(k),tsav(k) ! SY
       write(*,*) k,ksav(k),tsav(k)
      end do
  ! allocate and initialize additonal model arrays
      allocate (r(nmax),q(kper),qtime(2*nts+1),qb(nts+1),tcap(nts+2)) ! qb dimension changed to nts+1 29 Jan 2013, RLB
      allocate(qts(nts+1)) ! qts assigns q, transient surface flux, to timesteps
      if(outp(7)) then
       allocate(rik1(imax*(nts+1))) ! corrected 6/30/2008, RLB
       rik1=0.
      end if
      eps=1.0e-18
      nmax2=0
      nmin=1+nmax
      ncc=0;nccs=0
      tis=tiny(x1)
      r=0.; tcap=0; qts=0. ! Added 21 Feb 2013, RLB
      write(*,*) 'Starting computations of pressure head &
              &and factor of safety'
  ! prepare file, header, and arrays for generating list files, Added 4/21/2010, revised 5/3/2010, 12/6/2010
      if(flag < 0) then
       !umax=maxval(u) ! added 12/6/2010, RLB ! SY
       !uijz(1)=umax+1 ! SY
       uijz(1)=LIS_getNextUnitNumber() ! SY
       zmn(1)=minval(elev);zmx(1)=maxval(elev) ! added 9/14/2011, RLB
       write(LIS_logunit,*) 'DEM minimum & maximum elevations:', zmn(1), zmx(1) ! added 2/14/2012 RLB
       ftn = LIS_getNextUnitNumber() ! SY: for native code's u(2) that is closed at the end of this subroutine  
       call prpijz(ftn,profil,col,row,trigrs_header,&
        &vrsn) !Revised 12/23/2010 ! SY
      end if
  ! ---------------------------------------------------------------
      if(unsat0) then
       do j=1,nzon ! compute specific gravity of solids (gs(j)) from saturated unit weight of soil (uws(j))
        gs(j)=((uws(j)/uww)-ths(j))/(1-ths(j))
       end do
       if(mmax.lt.0) then ! infinite depth model
        mmax=20; nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
        write(*,*) 'Calling unsaturated infinite-depth model'
  !  4/14/2010 RLB added logging of main subroutines handling infiltration
        write(LIS_logunit,*) 'Unsaturated infinite-depth model, unsinf()' ! SY
        call unsinf(trigrs_imx1,&
          &ftn,ncc,nccs) ! SY
       else ! finite depth model
        write(*,*) 'Calling unsaturated finite-depth model'
        write(LIS_logunit,*) 'Unsaturated finite-depth model, unsfin()' ! SY
        call unsfin(trigrs_imx1,&
          &ftn,ncc,nccs)
       end if
      else ! Saturated zone extends to ground surface
       write(*,*) 'Ignoring unsaturated zone'
       if(tx==1 .and. nout==1) then
        deallocate (rizero,ks,ir)
  !  compute pore-pressure distributions and factor of safety
        outfil=trim(folder)//trim(profil)//trim(suffix)//'.txt'
        nccs=0
        if(mmax.lt.0) then ! infinite depth model
         mmax=20; nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
         write(*,*) 'Calling saturated infinite-depth model'
         write(LIS_logunit,*) 'Saturated infinite-depth model, iverson()' ! SY
         call iverson(trigrs_imx1,ftn,&
           &outfil) ! SY
        else ! finite depth model
          write(*,*) 'Calling saturated finite-depth model'
          write(LIS_logunit,*) 'Saturated finite-depth model, savage()' ! SY
          nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
          call savage(trigrs_imx1,ftn,&
            &outfil,nccs) ! SY
        end if
       else
        if(mmax.gt.0) then
         write(*,*) 'Calling multistep saturated finite-depth model'
         write(LIS_logunit,*) 'Multistep saturated finite-depth model, satfin()' ! SY
         nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
         call satfin(trigrs_imx1,&
           &ftn,nccs) ! SY
        else
         write(*,*) 'Calling multistep saturated infinite-depth model'
         write(LIS_logunit,*) 'Multistep saturated infinite-depth model, satinf()' ! SY
         mmax=20; nmn=1+mmax ! initialization of nmn added 1/3/2012 RLB
         call satinf(trigrs_imx1,&
          &ftn,nccs) ! SY
        end if
       end if
      end if
      !if(flag<=-1) close(u(2)) ! moved from subroutines 1 Dec 2011 RLB ! SY: see next few following lines 
      ! SY: Begin replacement code
      if(flag<=-1) then
       close (ftn) ! SY
       call LIS_releaseUnitNumber(ftn) ! SY
      end if
      ! SY: End replacement code
  ! *****************************************************************
  !  write output grid files
  ! 4/14/2010 added option to let file extension be either ".asc" or ".txt"
  101 continue
      write(*,*) 'Saving results'
      ti=tiny(param(1)) ! Changed from param(m) to param(1), 28 Jan 2013, RLB
      do j=1,nout
       write(stp,'(i4)') j
       stp=adjustl(stp)
       if (outp(3)) then ! minimum factor of safety
        tfg=0.
        do i=1,trigrs_imx1
         tfg(i)=fsmin(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(fminfil)//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,test1,param,&
        &outfil,ti,trigrs_header) ! SY
       end if
       if (outp(4)) then ! depth of minimum factor of safety
        tfg=0.
        do i=1,trigrs_imx1
         tfg(i)=zfmin(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(zfminfil)//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,test1,param,&
        &outfil,ti,trigrs_header) ! SY
       end if
       if (outp(5)) then ! pressure head at depth of minimum factor of safety
        tfg=0.
        do i=1,trigrs_imx1
         tfg(i)=pmin(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(pminfil)//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,test1,param,&
        &outfil,ti,trigrs_header) ! SY
       end if
       if (outp(1)) then ! computed water table
        tfg=0.
        do i=1,trigrs_imx1
         tfg(i)=wtab(i+(j-1)*imax)
        end do
        outfil=trim(folder)//trim(wtabfil)//el_or_dep//'_'//trim(suffix)//'_'//trim(stp)//grxt
        call ssvgrd(tfg,imax,pf1,row,col,test1,param,&
        &outfil,ti,trigrs_header) ! SY
       end if
      end do
      if (outp(7) .and. unsat0) then ! incremental basal flux, unsaturated zone
       do j=1,nts
        ir=0.
        do i=1,trigrs_imx1
         ir(i)=rik1(i+(j-1)*imax)
        end do
        ti=tiny(param(1))
        write(scratch,'(i6)') j
        scratch=adjustl(scratch)
        irfil='TRunszfluxTS'//trim(scratch)//trim(suffix)//grxt
        irfil=adjustl(irfil)
        outfil=trim(folder)//trim(irfil)
        call ssvgrd(ir,imax,pf1,row,col,test1,&
  &        param,outfil,ti,trigrs_header) ! SY
       end do
      end if
      if (ncc>0) then ! non-convergent cells, unsaturated zone (12/6/2010, changed unit # from 7 to 14)
       outfil=trim(folder)//ncvfil//'UZ_'//trim(suffix)//grxt
       call isvgrd(nvu,imax,pf1,row,col,test,test,&
       &trigrs_mnd,&
       & parami,outfil,ti,trigrs_header) ! SY
      end if
      if (nccs>0) then ! non-convergent cells, saturated zone
       outfil=trim(folder)//ncvfil//'SZ_'//trim(suffix)//grxt
       call isvgrd(nv,imax,pf1,row,col,test,test,&
       & trigrs_mnd,&
       & parami,outfil,ti,trigrs_header) ! SY
      end if
  
  ! SY: Begin extra code to release file unit numbers for uijz's
      do j=1,nout
       inquire (uijz(j),opened=unit_open)
       if(unit_open) then
        close (uijz(j)) ! SY
        call LIS_releaseUnitNumber(uijz(j)) ! SY
       end if
      end do
  ! SY: End extra code to release file unit numbers for uijz's
  
      write (*,*) 'TRIGRS finished!'
      write (LIS_logunit,*) 'TRIGRS finished normally' ! SY
      call date_and_time(date,time)
      write (LIS_logunit,*) 'Date ',date(5:6),'/',date(7:8),'/',date(1:4) ! SY
      write (LIS_logunit,*) 'Time ',time(1:2),':',time(3:4),':',time(5:6) ! SY
      !close(LIS_logunit) ! SY: DO NOT close log file
      write (LIS_logunit,*) '0' ! SY
      write (*,*) '0' ! SY
      !stop '0' ! SY
      !call LIS_endrun() ! SY
      return ! SY
  ! Error reporting
  410 continue
      write (*,*) 'Error opening output file'
      write (*,*) '--> ',outfil
      write (*,*) 'Check file path and status'
      write (LIS_logunit,*) 'Error opening output file' ! SY
      write (LIS_logunit,*) '--> ',outfil ! SY
      write (LIS_logunit,*) 'Check file path and status' ! SY
      !pause 'Press RETURN to exit' ! SY
      write (LIS_logunit,*) '410' ! SY
      write (*,*) '410' ! SY
      !stop '410' ! SY
      call LIS_endrun() ! SY

     end if
! SY: End running physics if last time step

    end if ! SY: if(alarmCheck) then

  end subroutine TRIGRS_run


!BOP
! !ROUTINE: TRIGRS_output
! 
! !INTERFACE: 
  subroutine TRIGRS_output(n)
! !USES: 
    use LIS_coreMod, only : LIS_rc, LIS_domain, LIS_masterproc
    use LIS_logMod, only : LIS_getNextUnitNumber, LIS_releaseUnitNumber
    use LIS_fileIOMod,  only : LIS_create_output_directory
    use LIS_historyMod, only : LIS_writevar_bin
! !ARGUMENTS: 
    integer, intent(in) :: n 
! 
! !DESCRIPTION: 
! 
!EOP
  end subroutine TRIGRS_output


!BOP
! 
! !ROUTINE: TRIGRS_final
! \label{TRIGRS_final}
! 
! !INTERFACE:
  subroutine TRIGRS_final()
! !USES:
!
! !DESCRIPTION:
!EOP

  end subroutine TRIGRS_final


  subroutine trini(filename,dg2rad)
! reads initialization file, R.L. Baum, USGS, latest revision 29 Mar 2013 
! SY: Subroutine indentation performed for integration into LIS
!
    use LIS_constantsMod, only : LIS_CONST_PI
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use LIS_coreMod,  only : LIS_rc ! SY

    use input_file_defs
    use input_vars
    use model_vars, only: smt

    implicit none
    
    character(len=*)   :: filename
    integer            :: ftn 
    integer            :: status
    logical            :: file_exists

    integer:: i,j,iz,linct
    real, intent(in):: dg2rad
    real :: tstar,tdif
    logical :: ltdif
    character (len=31):: scratch

    filename=adjustl(filename)
    inquire (file=trim(filename),exist=ans)
    ftn = LIS_getNextUnitNumber() ! SY
    if(ans) then
      open (ftn,file=trim(filename),status='old',err=201) ! SY
      write (*,*) 'Opening default initialization file'
    else
      write (LIS_logunit,*) 'Cannot locate initialization file' ! SY
      write (*,*) 'Cannot locate initialization file' ! SY
!      write (*,*) 'Cannot locate default initialization file, <tr_in.txt>' ! SY
!      write (*,*) 'Type name of initialization file and' ! SY
!      write (*,*) 'press RETURN to continue' ! SY
!      read (*,'(a)') init ! SY
!      init=adjustl(init) ! SY
!      open (uini,file=trim(init),status='old',err=201) ! SY
      call LIS_endrun() ! SY
    end if
    write (LIS_logunit,*) 'initialization file -->',trim(filename) ! SY
    write (LIS_logunit,*) '-- LISTING OF INITIALIZATION FILE --' ! SY
! write copy of data to log file
    linct=1
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) title; linct=linct+1 ! SY
    title=adjustl(title)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(title) ! SY
    write (*,*) title
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) tx,nmax,mmax,nzon; linct=linct+1 ! SY
    write (LIS_logunit,*) trim(heading) ! SY
    if(nmax<2) nmax=2 ! set minimum value for nmax
    write (LIS_logunit,*) tx,nmax,mmax,nzon ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) nzs,zmin,uww,nper,t ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) nzs,zmin,uww,nper,t ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) czmax,dep,crizero,slomin; linct=linct+1 ! SY
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) czmax,dep,crizero,slomin ! SY
    slomin=slomin*dg2rad ! convert minimum slope angle to radians
! allocate & read arrays for zone properties and initial conditions  	
    allocate (c(nzon),phi(nzon),uws(nzon),dif(nzon),&
      & ks(nzon),ths(nzon),thr(nzon),alp(nzon),unsat(nzon))
    allocate (igcap(nzon)) ! Added 2/15/2012 RLB
    c=0;phi=0
    uws=0;dif=0;ks=0;ths=0;thr=0;alp=0
    unsat=.true.;unsat0=.false.; igcap=.false.; igcapf=.true.
    do i=1,nzon
      read (ftn,*,err=420) scratch,iz ! property zone number ! SY
      linct=linct+1
      scratch=adjustl(scratch)
      read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
      heading=adjustl(heading)
      read (ftn,*,err=420) c(iz),phi(iz),&
      & uws(iz),dif(iz),ks(iz),ths(iz),thr(iz),alp(iz) ! zone parameters ! SY
      linct=linct+1
      write (LIS_logunit,*) trim(scratch),': ',iz ! SY
      write (LIS_logunit,*) trim(heading) ! SY
      write (LIS_logunit,*) c(iz),phi(iz),&
      & uws(iz),dif(iz),ks(iz),ths(iz),thr(iz),alp(iz) ! SY 
      if (c(iz)<0. .or. phi(iz)<0. .or. uws(iz)<0. .or.&
      & dif(iz)<0. .or. ks(iz)<0.) then ! check neg. vaules RLB 12/7/2010
        write(*,*) 'Error, negative property value in line ',linct
        write(*,*) 'Edit tr_in.txt and restart program TRIGRS'
        write(LIS_logunit,*) 'Error, negative property value in line ',linct ! SY
        write(LIS_logunit,*) 'Edit tr_in.txt and restart program TRIGRS' ! SY
        write (LIS_logunit,*) 'Negative property in trini()' ! SY
        write (*,*) 'Negative property in trini()' ! SY
        !stop 'Negative property in trini()' ! SY
        call LIS_endrun() ! SY
      end if
      if(ths(iz)<thr(iz)) then 
        write(*,*)'Error, Theta-resid. > Theta-sat. for property zone',iz
        write(*,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
        write(LIS_logunit,*)'Error, Theta-resid. > Theta-sat. for property zone',iz ! SY
        write(LIS_logunit,*)'Saturated infiltration model will be used for cells in zone', iz,'.' ! SY
        unsat(iz)=.false.
      end if
      if(alp(iz)<=0) then
        unsat(iz)=.false.
        write(*,*)'Negative or zero value of Alpha for property zone',iz
        write(*,*)'Saturated infiltration model will be used for cells in zone', iz,'.'
        write(LIS_logunit,*)'Negative or zero value of Alpha for property zone',iz ! SY
        write(LIS_logunit,*)'Saturated infiltration model will be used for cells in zone', iz,'.' ! SY
      end if
      if(unsat(iz)) then
        unsat0=.true. ! tracks whether any property zones are unsat.
        write(*,*)'Unsaturated infiltration model selected for cells in zone', iz,'.'
        write(LIS_logunit,*)'Unsaturated infiltration model selected for cells in zone', iz,'.' ! SY
      end if
    end do
    phi=phi*dg2rad ! convert phi angles to radians
! Allocate & read arrays for storm period data	
    allocate (cri(nper),capt(nper+2),rifil(nper)) 
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) (cri(j), j=1,nper) ! List of rainfall rates ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) (cri(j), j=1,nper) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) (capt(j), j=1,nper+1) ! List of times corresponding to change of rate ! SY
    !SY : Begin extra code to enforce requirement that t=capt(nper+1)
    if (t .NE. capt(nper+1)) then
      write (LIS_logunit,*) 'In LIS, t needs to be equal to final capt value in tr_in.txt' ! SY
      write (*,*) 'In LIS, t needs to be equal to final capt value in tr_in.txt' ! SY
      call LIS_endrun() ! SY
    end if
    !SY : End extra code to enforce requirement that t=capt(nper+1)
    !SY : Begin extra code to check if TRIGRS time steps are same as LIS time step
    do i=1,nper
      if ((capt(i+1)-capt(i)) .NE. TRIGRS_timestep) then
        write (LIS_logunit,*) 'tr_in.txt time steps are not same as LIS time step!!' ! SY
        write (*,*) 'tr_in.txt time steps are not same as LIS time step!!' ! SY
        call LIS_endrun() ! SY
      end if
    end do
    !SY : End extra code to check if TRIGRS time step is same as LIS time step
    capt(nper+2)=t ! for cases where t>capt(nper+1)
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) (capt(j), j=1,nper+1) ! SY
    ltdif=.false. ! Test time-step order and error message added 29 Jan 2013, RLB 
    do j=1,nper
      tdif=capt(j+1)-capt(j)
      if(tdif<0.) ltdif=.true. 
    end do
    if(ltdif) goto 424
!  path names of input files
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) slofil ! File name of slope angle grid (slofil) ! SY
    linct=linct+1
    slofil=adjustl(slofil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(slofil) ! SY
! added 4/21/2010 ******************************************
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) elevfil ! File name of digital elevation grid (elevfil) ! SY
    linct=linct+1
    elevfil=adjustl(elevfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(elevfil) ! SY
! **********************************************************
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) zonfil ! File name of property zone grid (zonfil) ! SY
    linct=linct+1
    zonfil=adjustl(zonfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(zonfil) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) zfil ! File name of depth grid (zfil) ! SY
    linct=linct+1
    zfil=adjustl(zfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(zfil) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) depfil ! File name of initial depth of water table grid (depfil) ! SY
    linct=linct+1
    depfil=adjustl(depfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(depfil) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=420) rizerofil ! File name of initial infiltration rate grid (rizerofil) ! SY
    linct=linct+1
    rizerofil=adjustl(rizerofil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(rizerofil) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    write (LIS_logunit,*) trim(heading) ! SY
    do j=1,nper
      read (ftn,'(a)',err=421) rifil(j) ! List of file names of rainfall intensity for each period, (rifil()) ! SY
      linct=linct+1
      rifil(j)=adjustl(rifil(j))
      write (LIS_logunit,*) trim(rifil(j)) ! SY
    end do
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=421) nxtfil ! File name of grid of D8 runoff receptor cell numbers (nxtfil) ! SY
    linct=linct+1
    nxtfil=adjustl(nxtfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(nxtfil) ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=421) ndxfil ! File name of list of defining runoff computation order (ndxfil) ! SY
    linct=linct+1
    ndxfil=adjustl(ndxfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(ndxfil) ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=421) dscfil ! File name of list of all runoff receptor cells  (dscfil) ! SY
    linct=linct+1
    dscfil=adjustl(dscfil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(dscfil) ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=421) wffil ! File name of list of runoff weighting factors  (wffil) ! SY
    linct=linct+1
    wffil=adjustl(wffil)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(wffil) ! SY
!  location of output files	
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=421) folder ! Folder where output grid files will be stored  (folder) ! SY
    linct=linct+1
    folder=adjustl(folder)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(folder) ! SY
!  output-file ID code
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,'(a)',err=421) suffix ! Identification code to be added to names of output files (suffix) ! SY
    linct=linct+1
    suffix=adjustl(suffix)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) trim(suffix) ! SY
!  output file selections	
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=421) rodoc ! Save grid files of runoff? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) rodoc ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=421) outp(3) ! Save grid of minimum factor of safety? ! SY
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(3) ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=421) outp(4) ! Save grid of depth of minimum factor of safety? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(4) ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=421) outp(5) ! Save grid of pore pressure at depth of minimum factor of safety? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(5) ! SY
    read (ftn,'(a)',err=421) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=422) outp(1), el_or_dep ! Save grid of computed water table depth or elevation for specified time steps? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(1), el_or_dep ! SY
    read (ftn,'(a)',err=422) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=422) outp(6) ! Save grid files of actual infiltration rate? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(6) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) outp(7) ! Save grid files of unsaturated zone basal flux? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(7) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=423) flag,spcg! Save listing of pressure head and factor of safety ("flag") & increment (spcg)? ! Added spcg 2/14/2012 RLB ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) flag,spcg! Added spcg 2/14/2012 RLB ! SY
    read (ftn,'(a)',err=423) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=423) nout ! Number of times to save output grids ! SY
    linct=linct+1
    if(nout<1) nout=1 ! must save at least one time; negative number not allowed
! Add code to limit nout to less than some specific number of values?
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) nout ! SY
    allocate (tsav(nout),ksav(nout),uijz(nout)) !Added uijz 12/07/2010 RLB
    tsav=0.;ksav=0
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) (tsav(j), j=1,nout) ! Times of output grids ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) (tsav(j), j=1,nout) ! SY
! user options  	
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) lskip ! Skip other timesteps? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) lskip ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) lany ! Use analytic solution for fillable porosity? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) lany ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) llus ! Estimate positive pressure head in rising water table zone ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) llus ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) lps0 ! Use psi0=-1/alpha? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) lps0 ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) outp(8) ! Log mass balance results? ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) outp(8) ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) flowdir ! Specify flow direction ! SY
    linct=linct+1! 29 Mar 2013, RLB, Added afer each read statement to correct error in counting line numbers from this point to end of input.
    flowdir=adjustl(flowdir)
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) flowdir ! SY
! added 19 Aug 2009 RLB-------------------
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) bkgrof ! Specify background flux offset  ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) bkgrof ! SY
! added 14-15 Apr 2010 RLB-------------------
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) lasc ! Specify grid file extension (.asc if true, default is .txt) ! SY 
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) lasc ! SY
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) lpge0 ! Ignore negative pore pressures in FS (p=0 if true) ! SY 
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) lpge0 ! SY
! added 2-15/2012 RLB ----
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) igcapf ! Ignore height of capilary fringe greater than depth to water table in selecting infiltration model? ! SY 
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) igcapf ! SY
!Parameters for deep pore-pressure estimate in SCOOPS ijz output (Added 20 Apr 2011, RLB): 
    read (ftn,'(a)',err=420) heading; linct=linct+1 ! SY
    heading=adjustl(heading)
    read (ftn,*,err=420) deepz,deepwat ! Deep point (floating point) depth below ground surface,  pressure option ('zero' or 'flow') ! SY
    linct=linct+1
    write (LIS_logunit,*) trim(heading) ! SY
    write (LIS_logunit,*) deepz,deepwat ! SY
!-----------------------------------------  	
    close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    write (LIS_logunit,*) '-- END OF INITIALIZATION DATA --' ! SY
    write (LIS_logunit,*) '' ! SY
    write (LIS_logunit,*) trim(title) ! SY
    write (LIS_logunit,*) '' ! SY
    do iz=1,nzon
      if(alp(iz)>=0 .and. unsat(iz)) then ! Revised 2/15/2012 RLB
        if(igcapf) igcap(iz)=.true. ! Use unsaturated model even for large capillary fringe
 ! Use saturated model for situations where much of interval is small time
        tstar=t*ks(iz)*alp(iz)/(ths(iz)-thr(iz))
        if(tstar < 4.d0*smt) igcap(iz)=.false.
      end if
      if (igcap(iz)) then
        write(*,*) '******** Zone ',iz, ' *********'
        write(*,*)'Using unsaturated infiltration model.'
        write(LIS_logunit,*) '******** Zone ',iz, ' *********' ! SY
        write(LIS_logunit,*)'Using unsaturated infiltration model.' ! SY
      else if (igcapf) then ! ................................................................ RLB 2/21/2012
        write(*,*) '******** Zone ',iz, ' *********'
        write(*,*)'Using saturated infiltration model to avoid'
        write(*,*)'early-time errors in unsaturated infiltration model.'
        write(LIS_logunit,*) '******** Zone ',iz, ' *********' ! SY
        write(LIS_logunit,*)'Using saturated infiltration model to avoid' ! SY
        write(LIS_logunit,*)'early-time errors in unsaturated infiltration model.' ! SY
      else if (alp(iz)<0) then ! ................................................................ RLB 1/30/2013 
        write(*,*) '******** Zone ',iz, ' *********'
        write(*,*)'Using saturated infiltration model; Alpha<0.'
        write(LIS_logunit,*) '******** Zone ',iz, ' *********' ! SY
        write(LIS_logunit,*)'Using saturated infiltration model; Alpha<0.' ! SY
      else
! ................................................................ RLB 2/15/2012          
        write(*,*) '******** Zone ',iz, ' *********'
        write(*,*)'Using unsaturated infiltration model.'
        write(*,*) 'Cells where water table is shallower than '
        write(*,*) '           ', 1./alp(iz)
        write(*,*) 'treated as tension saturated--Saturated infiltration model used.'
        write(LIS_logunit,*) '******** Zone ',iz, ' *********' ! SY
        write(LIS_logunit,*)'Using unsaturated infiltration model.' ! SY
        write(LIS_logunit,*) 'Cells where water table is shallower than ' ! SY
        write(LIS_logunit,*) '           ', 1./alp(iz) ! SY
        write(LIS_logunit,*) 'treated as tension saturated--Saturated infiltration model used.' ! SY
      end if
    end do
    write(*,*) '********  ********  ********  *********'
    write(LIS_logunit,*) '********  ********  ********  *********' ! SY
    return
201 continue
    write (*,*) '*** Error opening intialization file ***'
    write (*,*) '--> ',trim(filename)
    write (*,*) 'Check file location and name'
    write (LIS_logunit,*) '*** Error opening intialization file ***' ! SY
    write (LIS_logunit,*) '--> ',trim(filename) ! SY
    write (LIS_logunit,*) 'Check file location and name' ! SY
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '201 in trini()' ! SY
    write (*,*) '201 in trini()' ! SY
    !stop '201 in trini()' ! SY
    call LIS_endrun() ! SY
420 continue
    write (*,*) 'Error reading initialization file'
    write (*,*) '--> ',trim(filename), ' at line ',linct
    write (*,*) 'Check file contents and organization'
    write (LIS_logunit,*) 'Error reading initialization file' ! SY
    write (LIS_logunit,*) '--> ',trim(filename), ' at line ',linct ! SY
    write (LIS_logunit,*) 'Check file contents and organization' ! SY
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '420 in trini()' ! SY
    write (*,*) '420 in trini()' ! SY
    !stop '420 in trini()' ! SY
    call LIS_endrun() ! SY
421 continue
    write (*,*) 'Error reading initialization file'
    write (*,*) '--> ',trim(filename), ' at line ',linct
    write (*,*) 'Check file contents and organization'
    write (*,*) 'Number of file names/place holders for rainfall data'
    write (*,*) 'must equal nper.  List each on a separate line.'
    write (LIS_logunit,*) 'Error reading initialization file' ! SY
    write (LIS_logunit,*) '--> ',trim(filename), ' at line ',linct ! SY
    write (LIS_logunit,*) 'Check file contents and organization' ! SY
    write (LIS_logunit,*) 'Number of file names/place holders for rainfall data' ! SY
    write (LIS_logunit,*) 'must equal nper.  List each on a separate line.' ! SY
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '421 in trini()' ! SY
    write (*,*) '421 in trini()' ! SY
    !stop '421 in trini()' ! SY
    call LIS_endrun() ! SY
422 continue
    write (*,*) 'Error reading initialization file'
    write (*,*) '--> ',trim(filename), ' at line ',linct
    write (*,*) 'Check file contents and organization'
    write (*,*) 'Be sure to specify water table elevation or depth'
    write (*,*) 'along with choice or whether or not to save to file.'
    write (LIS_logunit,*) 'Error reading initialization file' ! SY
    write (LIS_logunit,*) '--> ',trim(filename), ' at line ',linct ! SY
    write (LIS_logunit,*) 'Check file contents and organization' ! SY
    write (LIS_logunit,*) 'Be sure to specify water table elevation or depth' ! SY
    write (LIS_logunit,*) 'along with choice or whether or not to save to file.' ! SY
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '422 in trini()' ! SY
    write (*,*) '422 in trini()' ! SY
    !stop '422 in trini()' ! SY
    call LIS_endrun() ! SY
423 continue
    write (*,*) 'Error reading initialization file'
    write (*,*) '--> ',trim(filename), ' at line ',linct
    write (*,*) 'Check file contents and organization'
    write (*,*) 'Be sure to specify vertical spacing increment for list file'
    write (*,*) 'on same line as output flag.'
    write (LIS_logunit,*) 'Error reading initialization file' ! SY
    write (LIS_logunit,*) '--> ',trim(filename), ' at line ',linct ! SY
    write (LIS_logunit,*) 'Check file contents and organization' ! SY
    write (LIS_logunit,*) 'Be sure to specify vertical spacing increment for list file' ! SY
    write (LIS_logunit,*) 'on same line as output flag.' ! SY
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '423 in trini()' ! SY
    write (*,*) '423 in trini()' ! SY
    !stop '423 in trini()' ! SY
    call LIS_endrun() ! SY
424 continue ! Error message added 29 Jan 2013, RLB 
    write (*,*) 'Time steps out of order in file'
    write (*,*) '--> ',trim(filename), ' at line ',linct
    write (*,*) 'List time steps in increasing order'
    write (LIS_logunit,*) 'Time steps out of order in file4' ! SY
    write (LIS_logunit,*) '--> ',trim(filename), ' at line ',linct ! SY
    write (LIS_logunit,*) 'List time steps in increasing order' ! SY
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '424 in trini()' ! SY
    write (*,*) '424 in trini()' ! SY
    !stop '424 in trini()' ! SY
    call LIS_endrun() ! SY
    
  end subroutine trini


!   SY: Subroutine indentation performed for integration into LIS
 !  subroutine to read an ascii grid file and store in a 1-d array
 !  by Rex L. Baum, USGS May 2001 latest revison 13 Mar 2013, RLB
 !  single precision
 !
  subroutine srdgrd(grd,pth,ncol,nrow,celsiz,nodat, &
     pf,pf1,ctr,imax,temp,infil,param,header) ! SY
     !pf,pf1,ctr,imax,temp,u,infil,param,header,u1) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    implicit none
    integer grd,pth,i,m,ncol,nrow,ctr,imax!,u,u1 ! SY
    double precision param(6),nodat,celsiz,cns,cew
    double precision east,west,north,south
    real pf(imax),pf1(grd),temp(pth),nodats
    character*14 header(6)
    character*255 infil
    integer ftn ! SY
 !  
    ftn = LIS_getNextUnitNumber() ! SY  
    !open(u,file=trim(infil),status='old',err=23) ! SY
    open(ftn,file=trim(infil),status='old',err=23) ! SY
    do 200, m=1,6
     read(ftn,*) header(m),param(m)
     header(m)=adjustl(header(m))
200 continue
 ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
    nodat=-9999.d0
    celsiz=-10.d0
    do 210, m=1,6
     if (trim(header(m)).eq.'ncols') ncol=int(param(m))
     if (trim(header(m)).eq.'nrows') nrow=int(param(m))
     if (trim(header(m)).eq.'cellsize') celsiz=param(m)
     if (trim(header(m)).eq.'NODATA_value') nodat=param(m)
     if (trim(header(m)).eq.'nodata_value') nodat=param(m)
     if (trim(header(m)).eq.'cols:') ncol=int(param(m))
     if (trim(header(m)).eq.'rows:') nrow=int(param(m))
     if (trim(header(m)).eq.'east:') east=param(m)
     if (trim(header(m)).eq.'west:') west=param(m)
     if (trim(header(m)).eq.'north:') north=param(m)
     if (trim(header(m)).eq.'south:') south=param(m)
210 continue
    if (celsiz.le.0) then
      cew=abs(east-west)/ncol
      cns=abs(north-south)/nrow
      if (cew.eq.cns) then
        celsiz=cew
      else
        celsiz=sqrt(cew*cns)
        write(*,*) 'Rectangular cells ',cew, ' X ', cns
        write(LIS_logunit,*) 'Rectangular cells ',cew, ' X ', cns
      end if
    end if
    if (ncol*nrow .gt. grd) then
     write(*,*) 'Grid file exceeds array size'
     write (*,*) '--> ',trim(infil)
     write(*,*) 'Check intialization file row and column values.'
     write(LIS_logunit,*) 'Grid file exceeds array size'
     write (LIS_logunit,*) '--> ',trim(infil)
     write(LIS_logunit,*) 'Check intialization file row and column values.'
     !pause 'Press RETURN to exit' ! SY
     close(ftn) ! SY
     call LIS_releaseUnitNumber(ftn) ! SY
     !close(LIS_logunit) ! SY: DO NOT close log file
     write(LIS_logunit,*) 'stopping in srdgrd'
     write(*,*) 'stopping in srdgrd'
     call LIS_endrun() ! SY
     !stop ! SY
    end if
    nodats=nodat
    ctr=0
    do 120, m=1,nrow
 !  next sequence of lines read data in but skips no_data values
 !  count maintained by ctr should coincide with node numbers from GIS
 !  pf1() keeps track of positions of nodata values so that results
 !  can be written out in grid format.
     read(ftn,*,end=125) (temp(i), i=1,ncol) ! SY
     do 250, i=1,ncol
      pf1(i+(m-1)*ncol)=temp(i)
      if(temp(i).ne.nodats) then
        ctr=ctr+1
        if (ctr>imax) then
          write(*,*) 'Number of data cells exceeds array size'
          write (*,*) '--> ',trim(infil)
          write(*,*) 'Check imax value in intialization file.'
          write(LIS_logunit,*) 'Number of data cells exceeds array size'
          write (LIS_logunit,*) '--> ',trim(infil)
          write(LIS_logunit,*) 'Check imax value in intialization file.'
          !pause 'Press RETURN to exit' ! SY
          close(ftn) ! SY
          call LIS_releaseUnitNumber(ftn) ! SY
          !close(LIS_logunit) ! SY: DO NOT close log file
        end if
        pf(ctr)=temp(i)
      end if
250  continue
120 continue
125 close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    return
23  continue
    write (*,*) '*** Error opening input file ***'
    write (*,*) '--> ',trim(infil)
    write (*,*) 'Check file name and location'
    write (LIS_logunit,*) '*** Error opening input file ***'
    write (LIS_logunit,*) '--> ',trim(infil)
    write (LIS_logunit,*) 'Check file name and location'
    !pause 'Press RETURN to exit' ! SY
    close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    !close(LIS_logunit) ! SY: DO NOT close log file
    write (LIS_logunit,*) '-23 in srdgrd()' ! SY
    write (*,*) '-23 in srdgrd()' ! SY
    call LIS_endrun() ! SY
    !stop '-23 in srdgrd()' ! SY
  end subroutine srdgrd ! SY


 !   SY: Subroutine indentation performed for integration into LIS 
 !  subroutine to read an ascii grid file of integer values
 !  and store in a 1-d array
 ! Rex L. Baum, USGS, spring 2001
 !  Dec 2003, added code to handle GRASS GIS ASCII grid files
 ! 
  subroutine irdgrd(grd,pth,ncol,nrow,celsiz,inodat,mnd, &! SY
      y,y1,ctr,imax,itemp,infil,param,header) ! SY
      !y,y1,ctr,imax,itemp,u,infil,param,header,u1) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    implicit none
    !integer grd,pth,i,m,ncol,nrow,ctr,imax,u,inodat,u1 ! SY
    integer grd,pth,i,m,ncol,nrow,ctr,imax,inodat ! SY
    integer y(imax),y1(grd),itemp(pth),mnd
    double precision param(6),celsiz,cns,cew,nodat
    double precision east,west,north,south
    character*14 header(6)
    character*255 infil
    integer ftn ! SY
 !  
    infil=adjustl(infil)
    ftn = LIS_getNextUnitNumber() ! SY
    !open(u,file=trim(infil),status='old',err=10) ! SY
    open(ftn,file=trim(infil),status='old',err=10) ! SY
    do 200, m=1,6
     read(ftn,*) header(m),param(m) ! SY
     header(m)=adjustl(header(m))
200 continue
 ! set default value of nodat & celsiz for use with GRASS GIS ascii files   
    inodat=-9999
    celsiz=-10.
    do 210, m=1,6
     if (trim(header(m)).eq.'ncols') ncol=int(param(m))
     if (trim(header(m)).eq.'nrows') nrow=int(param(m))
     if (trim(header(m)).eq.'cellsize') celsiz=param(m)
     if (trim(header(m)).eq.'cols:') ncol=int(param(m))
     if (trim(header(m)).eq.'rows:') nrow=int(param(m))
     if (trim(header(m)).eq.'east:') east=param(m)
     if (trim(header(m)).eq.'west:') west=param(m)
     if (trim(header(m)).eq.'north:') north=param(m)
     if (trim(header(m)).eq.'south:') south=param(m)
 ! 7/19/06 added allocatable for nodata postion in param array
     if (trim(header(m)).eq.'NODATA_value') then
       nodat=(param(m))
       mnd=m
       inodat=int(nodat)
     end if
     if (trim(header(m)).eq.'nodata_value') then
       nodat=(param(m))
       mnd=m 
       inodat=int(nodat)
     end if
210 continue
    if (celsiz.le.0) then
      cew=abs(east-west)/ncol
      cns=abs(north-south)/nrow
      if (cew.eq.cns) then
        celsiz=cew
      else
        celsiz=sqrt(cew*cns)
        write(*,*) 'Rectangular cells ',cew, ' X ', cns
        write(LIS_logunit,*) 'Rectangular cells ',cew, ' X ', cns
      end if
    end if
    if (ncol*nrow .gt. grd) then
     write(*,*) 'Grid file exceeds array size'
     write (*,*) '--> ',trim(infil)
     write(*,*) 'Check intialization file row and column values.'
     write(LIS_logunit,*) 'Grid file exceeds array size'
     write (LIS_logunit,*) '--> ',trim(infil)
     write(LIS_logunit,*) 'Check intialization file row and column values.'
     !pause 'Press RETURN to exit' ! SY
     close(ftn) ! SY
     call LIS_releaseUnitNumber(ftn) ! SY
     !close(LIS_logunit) ! SY: DO NOT close log file
     write(LIS_logunit,*) '-11 in subroutine irdgrd' ! SY
     write(*,*) '-11 in subroutine irdgrd' ! SY
     call LIS_endrun() ! SY
     !stop '-11 in subroutine irdgrd' ! SY
    end if
    ctr=0
    do 120, m=1,nrow
 !  next sequence of lines read data in but skips no_data values
 !  count maintained by ctr should coincide with node numbers from GIS
 !  y1() keeps track of positions of nodata values so that results
 !  can be written out in grid format.
     read(ftn,*,end=125) (itemp(i), i=1,ncol) ! SY
     do 250, i=1,ncol
      y1(i+(m-1)*ncol)=itemp(i)
      if(itemp(i).ne.inodat) then
       ctr=ctr+1
       if (ctr>imax) then
         write(*,*) 'Number of data cells exceeds array size'
         write (*,*) '--> ',trim(infil)
         write(*,*) 'Check imax value in intialization file.'
         write(LIS_logunit,*) 'Number of data cells exceeds array size'
         write (LIS_logunit,*) '--> ',trim(infil)
         write(LIS_logunit,*) 'Check imax value in intialization file.'
         !pause 'Press RETURN to exit' ! SY
         close(ftn) ! SY
         call LIS_releaseUnitNumber(ftn) ! SY
         !close(LIS_logunit) ! SY: DO NOT close log file
         write(LIS_logunit,*) '-12 in subroutine irdgrd' ! SY
         write(*,*) '-12 in subroutine irdgrd' ! SY
         call LIS_endrun() ! SY
         !stop '-12 in subroutine irdgrd' ! SY
       end if
       y(ctr)=itemp(i)
      end if
250  continue
120 continue
125 close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    return
10  continue
    write (*,*) '*** Error opening input file ***'
    write (*,*) '--> ',trim(infil)
    write (*,*) 'Check file name and location'
    write (LIS_logunit,*) '*** Error opening input file ***'
    write (LIS_logunit,*) '--> ',trim(infil)
    write (LIS_logunit,*) 'Check file name and location'
    !pause 'Press RETURN to exit' ! SY
    close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    !close(LIS_logunit) ! SY: DO NOT close log file
    write(LIS_logunit,*) '-10 in subroutine irdgrd' ! SY
    write(*,*) '-10 in subroutine irdgrd' ! SY
    call LIS_endrun() ! SY
    !stop '-10 in subroutine irdgrd'
  end subroutine irdgrd ! SY


  subroutine steady(sumex,imx1)
    use LIS_logMod, only : LIS_logunit
!   SY: Subroutine indentation performed for integration into LIS
    use grids
    use input_vars
    implicit none
    !integer:: i,acnt,ulog,imx1 ! SY
    integer:: i,acnt,imx1
    integer,parameter:: double=kind(1d0)
    real (double) :: sumex,b,rslo
    !  By Rex L. Baum, 1 April 2004
    !  Compute initial estimate of Isteady/Ks, rikzero(), and test values.
    !  Isteady must be < Ks.  If cos(slo) *cos(slo)<(Isteady/Ks), then an 
    !  inverted water table results.
    write(*,*) 'Testing and adjusting steady infiltration rates'
    acnt=0
    do i=1,imx1
       if(ks(zo(i))==0.) then ! prevent division by zero errors
          rikzero(i)=1.
       else
          rikzero(i)=rizero(i)/ks(zo(i))
       end if
    end do
    sumex=0
    do i=1,imx1
       rslo=slo(i)
       b=cos(rslo)*cos(rslo)
       if(rikzero(i)>=b) then
          rikzero(i)=cos(rslo) 
          acnt=acnt+1
          write(LIS_logunit,*) 'Adjusted steady infiltration rate, cell ',i
          write(LIS_logunit,*) 'Corrected rizero/Ks = ',rikzero(i)
          write(LIS_logunit,*) 'Original rizero, ks',rizero(i),&
               ks(zo(i))
          write(LIS_logunit,*) 'Set pore presssures to zero'
       end if
       if (depth(i)==0 .and. rizero(i)<0) then
          sumex=sumex-rizero(i)
       end if
    end do
    write(*,*) 'Adjusted steady infiltration rate at '&
         &,acnt,' cells'
    write(LIS_logunit,*) 'Adjusted steady infiltration rate at '&
         &,acnt,' cells'
    return
  end subroutine steady


! SY: Subroutine indentation performed for integration into LIS
  subroutine rnoff(grd,sumex,imx1,celsiz,param,parami,nodat,&
       &nodata,mnd,sctr,ncol,nrow,header,test1) ! SY

! Runoff routing for TRIGRS        
! By Rex L. Baum, 1 April 2004, latest revision 1 Jul 2011
    use LIS_logMod, only : LIS_endrun, LIS_logunit, LIS_getNextUnitNumber, &
         LIS_releaseUnitNumber, LIS_verify

    use input_file_defs
    use input_vars
    use grids

    implicit none
    integer,parameter:: double=kind(1d0)
    integer:: nodata,sctr,id,next 
    integer:: i,j,l,imx1,roflg,mnd !,m "m" removed 1 Feb 2013, RLB 
    integer:: grd!,unum ! SY
    integer:: ncol,nrow!,u(unum) ! SY
    logical:: logi(5)
    character (len=255):: infil,outfil
    character (len=31):: rofil, scratch,irfil
    character (len=14):: header(6)
    real:: rnof,inflx,test1
    real (double):: nodat,sumin,sumro,sumrf,sumex
    real (double):: celsiz,param(6),parami(6),ti
    integer:: ftn ! SY

! verify that the runoff routing input files are named in the initialization file
    logi(5)=.true.
    if(trim(nxtfil)=='no_input') logi(5)=.false.
    if(trim(ndxfil)=='no_input') logi(5)=.false.
    if(trim(dscfil)=='no_input') logi(5)=.false.
    if(trim(wffil)=='no_input') logi(5)=.false.
    ans=logi(5)
    if(logi(5)) then
! .... and that the named files acually exist
       inquire (file=trim(nxtfil),exist=logi(1))
       inquire (file=trim(ndxfil),exist=logi(2))
       inquire (file=trim(dscfil),exist=logi(3))
       inquire (file=trim(wffil),exist=logi(4))
       infil=trim(elfoldr)//'TIgrid_size.txt'
       infil=adjustl(infil)
       inquire (file=trim(infil),exist=logi(5)) ! added 2/24/2011, Revised to add folder 21 Feb 2013, RLB
       do i=1,5
        if(logi(i)) then
         continue
        else
         ans=.false.
        end if
       end do
    end if
    if(ans) then
     write(*,*) 'Starting runoff-routing computations'
! Added 2/2/2011, revised 4/1/2011 RLB
     if(nwf<imax) then 
      nwf=imax*2
      write(LIS_logunit,*) 'Error, nwf <imax, reset nwf to ',nwf *2
      write(LIS_logunit,*) 'If errors occur in runoff routing, check "TI_grid_size.txt" for correct value of nwf!'
      write(*,*) 'Error, nwf <imax, reset nwf to ',nwf*2 
      write(*,*) 'If errors occur in runoff routing, check "TI_grid_size.txt" for correct value of nwf!'
     end if
     allocate (dsc(nwf),wf(nwf))
     wf=0.; dsc=0
!  read numbers of subjacent cells
     call irdgrd(grd,col,ncol,nrow,celsiz,nodata,mnd,&
          nxt,pf2,sctr,&
          imax,itemp,&
          nxtfil,parami,header) ! SY
     write(LIS_logunit,*) 'Downslope cell grid'
     write(LIS_logunit,*) trim(nxtfil),sctr,' data cells'
     if(sctr/=imx1) write (LIS_logunit,*) 'Grid mismatch ',trim(nxtfil)
! read cell number index to determine order of computation for 
! runoff routing
     infil=ndxfil
     ftn = LIS_getNextUnitNumber() ! SY
     !open (u(18),file=infil,status='old',err=400) ! SY
     open (ftn,file=infil,status='old',err=400) ! SY
     do i=1,imax
      read (ftn,*,end=430) j,indx(j) ! SY
     end do
     close (ftn) ! SY
     call LIS_releaseUnitNumber(ftn) ! SY
! read cell numbers and weighting factors for runoff routing
     infil=dscfil
     ftn = LIS_getNextUnitNumber() ! SY
     !open (u(20),file=infil,status='old',err=400) ! SY
     open (ftn,file=infil,status='old',err=400) ! SY
     call irdswm(nwf,imax,ftn,nodata,dsc,dsctr) ! SY
     close (ftn) ! SY
     call LIS_releaseUnitNumber(ftn) ! SY
     infil=wffil
     ftn = LIS_getNextUnitNumber() ! SY
     !open (u(21),file=infil,status='old',err=400) ! SY
     open (ftn,file=infil,status='old',err=400) ! SY
     call srdswm(nwf,imax,ftn,test1,wf,dsctr) ! SY
     close (ftn) ! SY
     call LIS_releaseUnitNumber(ftn) ! SY
!  compute normalized infiltration intensities, Itransient/Ks, rik,
     write(LIS_logunit,*) ''
     write(LIS_logunit,*) '*** Runoff-routing computations ***'
     do j=1,nper
      sumin=0.d0
      sumrf=0.d0
      sumro=0.d0
      roflg=0
! SY: Begin extra code
      do  i=1,imax
       ri(i)=ri_Mat(j,i)
      end do
! SY: End extra code
! initialize ro and ir at beginning of each period
      do i=1,imx1
       ir(i)=ks(zo(i))
       ro(i)=0.
       sumrf=sumrf+ri(i)
      end do
! compute infiltration and runoff at each cell for each time period  
      do i=1,imx1
       id=indx(i)
       next=nxt(id)
       inflx=ro(id)+ri(id)
! case 1. exfiltration at cells where water table is initially at
! at the ground surface; does not track outflow from cells where
! the water table was initially below the surface and later filled up.
       if(depth(id)==0.0 .and. rizero(id)<0.0) then
        ir(id)=0
        rik(id+(j-1)*imax)=0.
        rnof=inflx-rizero(id)
        ro(id)=rnof
        do l=dsctr(id), dsctr(id+1)-1
         if (dsc(l).eq.id) then
          ro(dsc(l))=ro(dsc(l))+rnof*(wf(l)-1.)   
         else
          if (dsc(l)<1) write (*,*) dsc(l)
          if (dsc(l)>imx1) write (*,*) dsc(l)
          ro(dsc(l))=ro(dsc(l))+rnof*wf(l)
         end if
        end do    
        sumin=sumin+0
! case 2. infiltration, but available water > Ks    
       else if (ks(zo(id)).lt.inflx) then
        rik(id+(j-1)*imax)=1.d0
        rnof=inflx-ks(zo(id))
        ro(id)=rnof
        do l=dsctr(id), dsctr(id+1)-1
         if (dsc(l).eq.id) then
          ro(dsc(l))=ro(dsc(l))+rnof*(wf(l)-1.)   
         else
          if (dsc(l)<1) write (*,*) dsc(l)
          if (dsc(l)>imx1) write (*,*) dsc(l)
          ro(dsc(l))=ro(dsc(l))+rnof*wf(l)
         end if
        end do    
        sumin=sumin+ks(zo(id))
       else
! case 3. available water < Ks
        ir(id)=inflx
        rik(id+(j-1)*imax)=inflx/ks(zo(id))
        rnof=0.
        ro(id)=rnof
        ro(next)=ro(next)+rnof
        sumin=sumin+inflx
       end if
      end do
! compute total runoff  
      do i=1,imx1
       if (ro(i)>0.) roflg=1
       if(i==nxt(i)) sumro=sumro+ro(i)
      end do
      write(LIS_logunit,*) ''
      write (LIS_logunit,*) 'Mass Balance Totals for period ',j
      write (LIS_logunit,*) 'Precipitation + Exfiltration = &
       &Infiltration + Runoff'
      write (LIS_logunit,*) sumrf+sumex,' : ',sumin+sumro
      write (LIS_logunit,*) 'Infiltration   Runoff'
      write (LIS_logunit,*) sumin,sumro
      write (LIS_logunit,*) 'Precipitation   Exfiltration'
      write (LIS_logunit,*) sumrf,sumex
      write(LIS_logunit,*) '------------------------------------------'
! output ro() arrays for grids that have non-zero runoff
      if(roflg==1 .and. rodoc) then
! write only grids that are at or before the current time.
       if(capt(j)<=t) then
        ti=tiny(param(1))  ! Changed from param(m) to param(1), 28 Jan 2013, RLB
        write(scratch,'(i6)') j
        scratch=adjustl(scratch)
        rofil='TRrunoffPer'//trim(scratch)//trim(suffix)//grxt !4/26/2010 replaced "'.txt'" w/ "grxt"
        outfil=trim(folder)//trim(rofil)
! SY: Begin extra code
        do  i=1,grd
         pf1(i)=pf1_Mat(j,i)
        end do
! SY: End extra code
        !call ssvgrd(ro,imax,pf1,nrow,ncol,u(3),test1,param,u(19),& ! SY
        call ssvgrd(ro,imax,pf1,nrow,ncol,test1,param,&
          &    outfil,ti,header) ! SY
! save infiltration rate for grids that have non-zero runoff (changed unit number, from 6 to 8, 12/6/2010, RLB)
        if (outp(6)) then
         irfil='TRinfilratPer'//trim(scratch)//trim(suffix)//grxt
         outfil=trim(folder)//trim(irfil)
! SY: Begin extra code
         do  i=1,grd
          pf1(i)=pf1_Mat(j,i)
         end do
! SY: End extra code
         !call ssvgrd(ir,imax,pf1,nrow,ncol,u(8),test1,param,u(19),& ! SY
         call ssvgrd(ir,imax,pf1,nrow,ncol,test1,param,&
         &      outfil,ti,header) ! SY
        end if
       end if
! RLB--27 June 2011  
      else
       write(scratch,'(i6)') j
       scratch=adjustl(scratch)
       write(LIS_logunit,*) 'Runoff was zero at all grid cells for period ', scratch 
      end if  
     end do
    else
! skip runoff routing if any of the input files are missing
!  compute normalized infiltration intensities, Itransient/Ks, rik,
     write (*,*) 'Skipped runoff-routing computations;'
     write (*,*) 'Runoff routing input data did not exist.'
     write(LIS_logunit,*) ''
     write (LIS_logunit,*) 'Skipped runoff-routing computations;'
     write (LIS_logunit,*) 'Runoff routing input data did not exist.'
     write(LIS_logunit,*) '*******************************************'
! Added 2/2/2011 RLB
     nwf=1 
     allocate (dsc(nwf),wf(nwf))
     do j=1,nper
! SY: Begin extra code
      do  i=1,imax
       ri(i)=ri_Mat(j,i)
      end do
! SY: End extra code
! compute infiltration at each cell for each time period  
      do i=1,imx1
       if (ks(zo(i)).lt.ri(i)) then
        ir(i)=ks(zo(i))
        rik(i+(j-1)*imax)=1.d0
       else
        ir(i)=ri(i)
        rik(i+(j-1)*imax)=ri(i)/ks(zo(i))
       end if
      end do
! save actual infiltration rates to files  
! write only grids that are at or before the current time.
      if (outp(6)) then
       if(capt(j)<=t) then
        ti=tiny(param(1))
        write(scratch,'(i6)') j
        scratch=adjustl(scratch)
        irfil='TRinfilratPer'//trim(scratch)//trim(suffix)//grxt
        outfil=trim(folder)//trim(irfil)
! SY: Begin extra code
        do  i=1,grd
         pf1(i)=pf1_Mat(j,i)
        end do
! SY: End extra code
        !call ssvgrd(ir,imax,pf1,nrow,ncol,u(8),test1,param,u(19),& ! SY
        call ssvgrd(ir,imax,pf1,nrow,ncol,test1,param,&
          &      outfil,ti,header) ! SY
       end if
      end if
     end do
    end if ! end runouff routing and infiltration rates 
    return
400 continue
    write (*,*) '*** Error opening file ***'
    write (*,*) '--> ',infil
    write (*,*) 'Check file name and location'
    write (LIS_logunit,*) '*** Error opening file ***'
    write (LIS_logunit,*) '--> ',infil
    write (LIS_logunit,*) 'Check file name and location'
    !pause 'Press RETURN to exit' ! SY
    !stop '400' ! SY
    write (LIS_logunit,*) '400' ! SY
    write (*,*) '400' ! SY
    call LIS_endrun() ! SY
430 continue
    write (*,*) 'Attempt to read past end of file'
    write (*,*) '--> ',infil
    write (*,*) 'Check file contents and value of Imax'
    write (LIS_logunit,*) 'Attempt to read past end of file'
    write (LIS_logunit,*) '--> ',infil
    write (LIS_logunit,*) 'Check file contents and value of Imax'
    !pause 'Press RETURN to exit' ! SY
    !stop '430' ! SY
    write (LIS_logunit,*) '430' ! SY
    write (*,*) '430' ! SY
    call LIS_endrun() ! SY
  end subroutine rnoff


!  
!  Reads a list of integer values that represent a 2-D matrix into 
!  two 1-D arrays.  One array is a allocatable that tracks the starting
!  location of each row in the matrix.  This scheme is advantageous 
!  for storing sparse "sawtooth" arrays. Such arrays have rows of 
!  variable length and all of the non-zero values are at the left 
!  ends of the rows. 
!  by Rex L. Baum, spring 2001 
!  
! Implementation:
!  	open (u,file=infil,status='old')
!  	call irdswm(len,jmax,u,test,x,u2,ctr)
!  	close(u)
!  
! 
  subroutine irdswm(len,jmax,u,test,x,ctr)
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
        LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    implicit none
    integer j,k,jmax,u,val,test,len,x(len)
    integer ctr(jmax+1)
    j=0
    ctr(1)=0
    k=0
100 continue
    read (u,*,end=110) val
    !  Check for marker value at end of a row
    if (val.eq.test) then
       read (u,*,end=110) j
       ctr(j)=k+1
       if (j.gt.jmax) go to 200
       go to 100
    end if
    k=k+1
    x(k)=val
    go to 100
110 continue
    ctr(jmax+1)=k+1
    return
200 continue
    write(*,*) 'Array bound <jmax> exceeded'
    !pause 'Press return/enter to quit.' ! SY
    write(LIS_logunit,*) 'Stopping in irdswm' ! SY
    write(*,*) 'Stopping in irdswm' ! SY
    call LIS_endrun() ! SY 
    !stop ! SY
  end subroutine irdswm


!  
!  Reads a list of real values that represent a 2-D matrix into 
!  two 1-D arrays.  One array is a allocatable that tracks the starting
!  location of each row in the matrix.  This scheme is advantageous 
!  for storing sparse "sawtooth" arrays. Such arrays have rows of 
!  variable length and all of the non-zero values are at the left 
!  ends of the rows. 
!   by Rex L. Baum, USGS, Spring 2001
!  
! 
  subroutine srdswm(len,jmax,u,test,x,ctr)
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    implicit none
    integer j,k,jmax,u,len
    integer ctr(jmax+1)
    real x(len),val,test
    j=0
    ctr(1)=0
    k=0
100 continue
    read (u,*,end=110) val
    !  Check for marker value at end of a row
    if (val.eq.test) then
       read (u,*,end=110) j
       ctr(j)=k+1
       if (j.gt.jmax) go to 200
       go to 100
    end if
    k=k+1
    x(k)=val
    go to 100
110 continue
    ctr(jmax+1)=k+1
    return
200 continue
    write(*,*) 'Array bound <jmax> exceeded'
    !pause 'Press return/enter to quit.' ! SY
    write(LIS_logunit,*) 'stop in srdswm' ! SY
    write(*,*) 'stop in srdswm' ! SY
    call LIS_endrun() ! SY
    !stop ! SY
  end subroutine srdswm


! SY: Subroutine indentation performed for integration into LIS
 !  subroutine to save a 1-d array as a 2-d grid, including nodata values
 !  By Rex L. Baum USGS, December 2000, latest revison 4/19/2010
 !
  !subroutine ssvgrd(z3,grd,z2,nrow,ncol,u,nodata,param,u1, & ! SY
  subroutine ssvgrd(z3,grd,z2,nrow,ncol,nodata,param, &
       outfil,ti,header) ! SY
! !USES:
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify

    implicit none
    !integer i,j,nrow,ncol,u,grd,ctr,u1,ia,m ! SY
    integer i,j,nrow,ncol,grd,ctr,ia,m ! SY
 !  Map 2-d array, z2, onto 1-d array, z3.
 !  Reverse rows and columns, because Fortran stores arrays
 !  by column.  z2 includes nodata values and is used to map
 !  z3 to the grid format
    real z2(ncol,nrow),z3(grd),nodata
    double precision test, param(6),ti,a
    character*14 header(6)
    character*1 sp
    character*255 outfil
    character*31 scratch
    integer ftn ! SY
    sp=char(32)
    ctr=0
 !  open and initialize output grid files 
    outfil=adjustl(outfil)
    ftn = LIS_getNextUnitNumber() ! SY
    !open (u,file=trim(outfil), status='unknown',err=20) !SY
    open (ftn,file=trim(outfil), status='unknown',err=20) ! SY
    do 100, m=1,6
     a=param(m)
     ia=int(param(m))
 ! 19 Apr 2010 added special cases for m=1,2, and 6
     if(m.le.2) then
       write(scratch,'(i16)') ia
     else if(m.eq.6) then
       write(scratch,'(f15.0)') a      
     else if(abs(a-ia).le.ti) then
       write(scratch,'(i16)') ia
     else
       write(scratch,*) a
     end if
     scratch=adjustl(scratch)
     write(ftn,1012) trim(header(m)),trim(scratch) ! SY
100 continue
!  write grid 
    do 120, i=1,nrow
     do 110, j=1,ncol
      test=abs(z2(j,i)-nodata)
      if (test.le.0.1)then
       write(scratch,1004) z2(j,i)
       scratch=adjustl(scratch)
       if(j.ne.ncol) then
        write(ftn,1010) trim(scratch),sp ! SY 
       else
        write(ftn,1011) trim(scratch) ! SY
       end if
      else
       ctr=ctr+1
       if (z3(ctr).eq.0)  then
        write(scratch,1003) z3(ctr)
       else
        write(scratch,1004) z3(ctr)
       end if
       scratch=adjustl(scratch)
       if(j.ne.ncol) then
        write(ftn,1010) trim(scratch),sp 
       else
        write(ftn,1011) trim(scratch) 
       end if
      end if
110  continue
120 continue
    close (ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    return
20  continue
    write (*,*) 'Error opening output file'
    write (*,*) '--> ',outfil
    write (*,*) 'Check file name and status'
    write (LIS_logunit,*) 'Error opening output file'
    write (LIS_logunit,*) '--> ',outfil
    write (LIS_logunit,*) 'Check file name and status'
    !pause 'Press RETURN to exit' ! SY
    write (LIS_logunit,*) '-20 in ssvgrd()' ! SY
    write (*,*) '-20 in ssvgrd()' ! SY
    !stop '-20 in ssvgrd()' ! SY
    call LIS_endrun() ! SY
1003 format(f2.0)
1004 format(g12.4)
1010 format(a,a,$)
1011 format(a)
1012 format(t1,a,t15,a)
  end subroutine ssvgrd ! SY 


!   SY: Subroutine indentation performed for integration into LIS
  !subroutine unsinf(imx1,ulog,u1,ncc,nccs) ! SY
  subroutine unsinf(imx1,u1,ncc,nccs) ! SY
! By Rex L. Baum, USGS, Latest revision, 29 Jan 2013 
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    use grids
    use input_vars 
    use model_vars
    use input_file_defs
    implicit none
    !integer::i,j,jf,k,ulog,u1,imx1,ncc,nccs,nmax0 ! ,nmx ! SY
    integer::i,j,jf,k,u1,imx1,ncc,nccs,nmax0 ! ,nmx ! SY
    integer::nmn1,nmin1,nmax3,nmns ! ,nmxp,nmnp,nmxs
    integer:: ivctr
    logical:: lcv,lwt ! ,lcvs
    real:: delwt,dwt,zwt,qbij(nts+1) 
    real::testqk,tolqk ! Added 2/2/2011 RLB
    real (double)::rf(nzs+1),finf,vqt,qta,al,qzmax 
    real (double)::ddwt,sqin,intq(nts+1),b,dhwt(nts+1),delh 
    real (double)::qtn(2*nts+1),intq1(nts+1),vqtn,cd
    nmax3=0;nmax0=0
    nmn1=nmax+1;nmin1=nmax+1
    ivctr=0 !nmxp=0;nmnp=0;   ! Added 29 Jan 2013, RLB 
    intq=0.d0; intq1=0.d0 ! Added 29 Jan 2013, RLB 
    write(LIS_logunit,*) 'Starting coupled saturated & unsaturated zone'
    write(LIS_logunit,*) 'computations for infinite-depth saturated zone'
    write(*,*) 'Starting coupled saturated & unsaturated zone'
    write(*,*) 'computations for infinite-depth saturated zone'
    write(*,*) 'Cells completed: '
! loop over all grid cells
    finf=10.
    grid_loop: do i=1,imx1 
     !if (mod(i-1,2000)==0) write (*,*) i-1 ! cells completed ! SY
     ! SY: Begin
     if (mod(i-1,2000)==0) then
       write (*,*) i-1 ! cells completed
       write (LIS_logunit,*) i-1 ! cells completed
     end if
     ! SY: End
     if(slo(i)<slomin) then ! default values for gently sloping cells 
       do jf=1,nout
         fsmin(i+(jf-1)*imax)=finf+1.
         zfmin(i+(jf-1)*imax)=zmax(i)
         pmin(i+(jf-1)*imax)=0.
       end do
       cycle
     end if
     lcv=.true. !;lcvs=.true.
     q=0.;qb=0 ! qb initialization added 29 Jan 2013, RLB 
     tolqk=ks(zo(i))*5.e-07 ! Moved 29 Jan 2013, RLB 
     do j=1,kper
       if(j>nper) then
         q(j)=0.
       else
! 8/18/2009 RLB added optional offset of background flux to prevent excessive drying during periods of zero infiltration.      
        if(bkgrof) then
          q(j)=ks(zo(i))*(rik(i+(j-1)*imax)+rikzero(i))
! RLB 2/2/2001 revised test
          testqk=q(j)-(ks(zo(i))+rizero(i))
!          tolqk=ks(zo(i))*5.e-07
          !if(testqk>tolqk) write (LIS_logunit,*) '*q>Ks+ri!', i,j,q(j),ks(zo(i))+rizero(i) ! Added 2/2/2011 RLB ! SY: see next line
          if(testqk>tolqk) write (LIS_logunit,'(A,2I3,2ES14.6E2)') '*q>Ks+ri!', i,j,q(j),ks(zo(i))+rizero(i) ! Added 2/2/2011 RLB ! SY
        else
          q(j)=ks(zo(i))*rik(i+(j-1)*imax)
          testqk=q(j)-ks(zo(i))
          tolqk=ks(zo(i))*5.e-07          
          !if(testqk>tolqk) write (LIS_logunit,*) '*q>Ks!', i,j,q(j),ks(zo(i)) ! Moved 2/2/2011 RLB ! SY: see next line
          if(testqk>tolqk) write (LIS_logunit,'(A,2I3,2ES14.6E2)') '*q>Ks!', i,j,q(j),ks(zo(i)) ! Moved 2/2/2011 RLB ! SY
        end if
       end if
     end do
     qmax=maxval(q)
     b=cos(slo(i))
! next lines compute depth to top of capilarly fringe.
     if(unsat(zo(i))) then
       dcf=depth(i)-1.d0/(alp(zo(i)))
     else
       dcf=0.
     end if 
     if(lps0 .and. unsat(zo(i))) then
       dusz=depth(i)-1.d0/(alp(zo(i)))
     else 
       dusz=depth(i)
     end if
! set value of beta (Iverson's beta line)
! and maximum drainage rate at water table, qzmax    
     cd=1.0d0 ! complete drainage for infinite depth basal flow boundary         
     select case (flowdir)
      case ('slope')
       beta=b*b
       qzmax=(1.d0-beta)*cd*ks(zo(i)) 
      case ('hydro')
       beta=1.d0
       qzmax=0.d0 
      case default
       beta=b*b-rikzero(i) ! 2/12/09 corrected formula for beta & qzmax
       qzmax=(1.d0-beta)*cd*ks(zo(i))-cd*rizero(i) 
     end select
     delwt=0.;dwt=depth(i);zwt=depth(i);ddwt=depth(i)
     lwt=.false.
     ts=tmin
     if(dcf>0. .and. (unsat(zo(i)) .or. igcap(zo(i)))) then ! updated to enforce non-zero depth, 4/18/2013, RLB
! compute flux and pore-pressure rise at each time step   
      vqt=0.;vqtn=0.;qta=rizero(i);sqin=0.
      al=alp(zo(i))*b*b
      if(outp(8)) write (LIS_logunit,*) 'ts,    qt '! times and basal flux to log file
      call roots(nmax,r,dusz,al,eps,pi)
      flux_loop: do j=1,2*nts+1
       call flux(i,kper,ts,j,lcv,ncc,nvu(i),lwt) 
       if(nmax1>nmax2) nmax2=nmax1
       if(nmn<nmin) nmin=nmn
       if(qt<rizero(i)) qt=rizero(i)
! RLB 2/3/2011, Added case for bkgrof=.true.               
       if(bkgrof) then
         if(qt>ks(zo(i))+rizero(i)+tolqk) then
           write(LIS_logunit,*) 'Error! Basal flux exceeds Ks at'
           write(LIS_logunit,*) 'cell ',i, ', timestep ',j 
           write(LIS_logunit,*) 'flux ',qt, ', Ks ',ks(zo(i))+rizero(i) 
           qt=ks(zo(i))
         end if
       else
         if(qt>ks(zo(i))+tolqk) then
           write(LIS_logunit,*) 'Error! Basal flux exceeds Ks at'
           write(LIS_logunit,*) 'cell ',i, ', timestep ',j 
           write(LIS_logunit,*) 'flux ',qt, ', Ks ',ks(zo(i))
           qt=ks(zo(i))
         end if
       end if
       if(outp(8)) write (LIS_logunit,*) ts,qt ! times and basal flux to log file
! drain off excess basal flux
       if(qt>qzmax) then
         qtime(j)=qt-qzmax
       else
         qtime(j)=0.d0
       end if
       qtn(j)=qt
       qta=qt
       ts=ts+tinc/2.d0
      end do flux_loop
      call dsimps(nts,tinc/2.d0,qtime,intq)
      call dsimps(nts,tinc/2.d0,qtn,intq1)
      !if(outp(8)) write (LIS_logunit,*) 'Time, Cumulative volume in, Cumulative background flux,& ! SY: Changed to line below
      if(outp(8)) write (LIS_logunit,'(500A)') 'Time, Cumulative volume in, Cumulative background flux,& 
              & Cumul. volume out,Cuml. absorbed, Cuml. qin-qout, qout not drained, Water table rise'
      ts=tmin
      wt_rise_loop: do j=1,nts+1
        jf=jsav(j)
        rf=0.0 
        vqt=intq(j)-ts*rizero(i)
        if(vqt<0.) vqt=0.d0
        vqtn=intq1(j)-ts*rizero(i)
        sqin=0.
        sum_q_in_loop: do k=1,nper
          if(ts>capt(k) .and. ts<=capt(k+1)) then
                qts(j)=q(k) !; write(*,*) 'j,k,qts(j),q(k),ts,capt(k) ', j,k,qts(j),q(k),ts,capt(k) ! Added 21 Feb 2013, RLB 
          endif
          if(ts>=capt(k+1)) then
            sqin=sqin+(capt(k+1)-capt(k))*q(k)
          end if
          if(ts>capt(k) .and. ts<capt(k+1)) then
            sqin=sqin+(ts-capt(k))*q(k)
          end if
        end do sum_q_in_loop
        if(jf>0 .and. lskip) then
! compute usaturated zone pressure & water content
          call unsth(i,j,ncc,kper,ts,nmax0,& 
          &lcv,vqt,delh,nmn1,sqin,vqtn) ! SY
          !&lcv,LIS_logunit,vqt,delh,nmn1,sqin,vqtn) ! SY
          dhwt(j)=delh
          if(nmax0>nmax3) nmax3=nmax0
          if(nmn1<nmin1) nmin1=nmn1
        else if(lskip) then
          continue
        else
          call unsth(i,j,ncc,kper,ts,nmax0,& 
          &lcv,vqt,delh,nmn1,sqin,vqtn) ! SY
          !&lcv,LIS_logunit,vqt,delh,nmn1,sqin,vqtn) ! SY
          dhwt(j)=delh
          if(nmax0>nmax3) nmax3=nmax0
          if(nmn1<nmin1) nmin1=nmn1
        end if        
        tcap(j)=ts ! pass to diffusion subroutine
        tcap(j+1)=ts+tinc
        ts=ts+tinc
        if(jf>0 .and. lskip) then
! compute pressure diffusion in saturated zone
!          call pstpi(u1,dhwt,dwt,outfil,&
          call pstpi(u1,dhwt,dwt,&
                  & i,j-1,rf,tcap(j),jf) ! Revised 29 Jan 2013, RLB ! SY
                  !& LIS_logunit,i,j-1,rf,tcap(j),jf) ! Revised 29 Jan 2013, RLB ! SY 
!                  & LIS_logunit,i,j-1,rf,nccs,lcvs,nmx,tcap(j),jf)
!                  nmxp=nmx;nmnp=nmn
          if(j>1) then
! Check change in water table depth and adjust dusz if needed
            delwt=abs(dwt-zwt)*1000.
            if(delwt>dwt) then
              lwt=.true.
              dwt=zwt
            end if
          end if   
        else if(lskip) then
          continue
        else
! compute pressure diffusion in saturated zone
          call pstpi(u1,dhwt,dwt,&
                  & i,j-1,rf,tcap(j),jf) ! Revised 29 Jan 2013, RLB ! SY
                  !& LIS_logunit,i,j-1,rf,tcap(j),jf) ! Revised 29 Jan 2013, RLB ! SY 
!                  & LIS_logunit,i,j-1,rf,nccs,lcvs,tcap(j),jf)
!                  nmxp=nmx;nmnp=nmn
          if(j>1) then
! Check change in water table depth and adjust if needed
            delwt=abs(dwt-zwt)*1000.
            if(delwt>dwt) then
              lwt=.true.
              dwt=zwt
            end if
          end if   
        end if        
      end do wt_rise_loop
! map unsaturated zone outflux to grid
! there are 2*nts+1 increments in qtime()
      do k=1,nts  
       qb(k)=qtime(2*k+1)
       if(outp(7)) rik1(i+(k-1)*imax)=qb(k)/ks(zo(i))
      end do
     else ! top of capillary fringe at ground surface, so use surface flux
      qb=0. ! initialize qb for case where ts>capt(nper+1)
      dwt=depth(i)
      delh=0.;rf=0.;zwt=depth(i)
      do j=1,nts+1
       do k=1,kper
        if(ts>=capt(k) .and. ts<=capt(k+1)) qb(j)=q(k)
       end do
       if(outp(7)) rik1(i+(j-1)*imax)=qb(j)/ks(zo(i))
       tcap(j)=ts ! pass to diffusion subroutine
       ts=ts+tinc
      end do
      do j=1,nts+1
       qbij(j)=qb(j)/ks(zo(i))
      end do
      rf=0.
      ivctr=ivctr+1
      call ivestp(u1,qbij,&
           & i,rf) ! Revised 29 Jan 2013, RLB ! SY
           !& LIS_logunit,i,rf) ! Revised 29 Jan 2013, RLB ! SY 
!                & LIS_logunit,i,rf,nccs,lcvs,nmxs)
      nmns=nmn
     end if
    end do grid_loop
    write (*,*) imx1, ' cells completed' 
    write (LIS_logunit,*) imx1, ' cells completed' 
!
    if(nmin>nmax2) nmin=nmax2; if(nmin1>nmax3) nmin1=nmax3
    write(LIS_logunit,*) 'Convergence data for unsaturated zone:'
    write(LIS_logunit,*) 'Maximum terms used by Fourier series', nmax2,nmax3
    write(LIS_logunit,*) 'Minimum terms used by Fourier series', nmin,nmin1
    write(LIS_logunit,*) 'Unsaturated zone nonconvergent cells: '
    write(LIS_logunit,*) ncc
!if(nmnp>nmxp) nmnp=nmxp; if(nmns>nmxs) nmns=nmxs
    write(LIS_logunit,*) 'Convergence data for saturated zone:'
! ivestp() and pstpi() use single-term solutions to compute pressure head--Convergence is not an issue 
    write(LIS_logunit,*) 'Terms used by infinite-depth '
    write(LIS_logunit,*) 'solutions ivestp() and pstpi():', 1, 1
    write(LIS_logunit,*) 'Saturated-zone nonconvergent cells: '
    write(LIS_logunit,*) nccs !initialized to zero in trigrs main
    write(LIS_logunit,*) 'Cells using ivestp() = ',ivctr
    return
  end subroutine unsinf


!   SY: Subroutine indentation performed for integration into LIS
  subroutine roots(nmax,r,d,beta,eps,pi)
 ! by Rex L. Baum and W.Z. Savage, May 2004
! Latest revision 28 Jan 2013, RLB 	
    implicit none
    integer i,n,nmax,imax,mmax
 !real beta promoted to double precision 12/6/2006, df1 & dfsgn added 29 Jan 2013, RLB 
    double precision a,b,x0,f,df,x1,d,df1,dfsgn
    double precision r(*),pi,eps,beta
    double precision rlb,rub,rtb,fl,fu,dx
    a=beta*d
    b=2.0
    x0=eps
    do 15 n=1,nmax
 !  set upper and lower bound for each interval
     rlb=(pi/a)*float(n-1)+eps
     rub=(pi/a)*float(n)
     fl= sin(a*rlb)+b*rlb*cos(a*rlb)
     fu= sin(a*rub)+b*rub*cos(a*rub)
     if(fl.lt.0) then
       rtb=rlb
       dx=rub-rlb
     else
       rtb=rub
       dx=rlb-rub
     end if
 ! find initial estimate of roots by bisection
     mmax=5
     call dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
     if(imax.ge.1) go to 11
     do 10 i=1,100
      f= sin(a*x0)+b*x0*cos(a*x0)
 !    df=(a-b)*cos(a*x0)-a*b*x0*sin(a*x0) !formula appears to contain sign error 28 Jan 2013, RLB
      df=(a+b)*cos(a*x0)-a*b*x0*sin(a*x0)
 ! Next lines added to prevent division by zero, 28-30 Jan 2013, RLB.    
      df1=abs(df) 
      if(df1 .lt. eps) then
       dfsgn=1.d0
       if(df .lt. 0.d0) dfsgn=-1.d0
       df=eps*dfsgn
      end if
      x1=x0-f/df
      if(abs((x1-x0)/x0) .lt. eps) then
        if(x1.gt.rlb .and. x1.lt.rub) then
          imax=i
          x0=x1
          goto 11
        else
 ! if Newton-Raphson method converges on the wrong root, then
 ! revert to bisection
         mmax=200
         call dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
         imax=imax+5
         if(imax.ge.6) go to 11
        end if
      else
        x0=x1
      end if
10   continue
 ! if Newton-Raphson method did not converge, then
 ! revert to bisection
     mmax=200
     call dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
     imax=imax+5
11   continue
     r(n)=x0
     x0=x0+pi/a
15  continue
    return
  end subroutine roots ! SY


!   SY: Subroutine indentation performed for integration into LIS
  subroutine dbsct(mmax,eps,a,b,rtb,dx,x0,imax)
 ! By Rex L. Baum, May 14, 2004
 ! Uses bisection method as described in Press and others, 1986, p. 246-247.
    implicit none
    integer imax,m,mmax
    double precision a,b,x0,eps
    double precision rtb,fm,dx
    imax=0
    do 10,m=1,mmax
     dx=dx/2.
     x0=dx+rtb
     fm=sin(a*x0)+b*x0*cos(a*x0)
     if(fm.le.0.) rtb=x0
     if(abs(dx).le.eps .or. fm.eq.0.) then
       imax=m
       return
     end if 
10  continue 
    return
  end subroutine dbsct ! SY


!   SY: Subroutine indentation performed for integration into LIS
! by W.Z. Savage with modifications by R.L. Baum, Latest revision 29 Jan 2013, RLB
  subroutine flux(ic,iper,t0,j1,lcv,ncc,nv0,lwt)
    use grids
    use input_vars
    use model_vars
    implicit none
    integer:: ic
    integer:: j1,iper,i,k,k1,ncc,nv0 !,i1 removed 12/28/2010
    real (double):: al,qa 
    real (double):: t0,tf,q1a,q1b,q2a,q2b,qta,qtb,qtop !,bot
    real (double):: q2old,delq2,tdif1,tdif2,tstar
    real (double):: q0a,q0b 
    real (double):: z,psih,qlb,ck,th,b1
    real (double):: bot(nmax),qtops(nmax),dseep ! added qtops & made arrays 12/28/2010
!
    logical:: lcv,lwt
    b1=cos(slo(ic))
! coordinate transformation, corrected 9/5/06
    al=alp(zo(ic))*b1*b1
    qa=rizero(ic)
    tf=al*ks(zo(ic))/(ths(zo(ic))-thr(zo(ic)))
    qt=0.0; qta=0.0; qtb=0.0 
    nmax1=0
    nmn=nmax+1
    do k=1,nmax
     bot(k)=1.+al*dusz/2.+2.*al*dusz*r(k)**2
     qtops(k)=r(k)*sin(r(k)*al*dusz)
    end do
    flux_time_loop: do  i=1,iper
!  i1=1+(i-1)*tx
     tdif1=t0-capt(i)
     if (tdif1 < 0.) exit ! jump out of loop rather than compute extra zeros.
     if(tdif1 > 0.0) then
      tstar=tf*tdif1
      if(tstar < smt) then
! Early-time solution (ETS) ...............
       z=dusz
       qlb=ks(zo(ic))
! Corrected ths, thr to include index 12/20/2010, RLB
       call smallt(tstar,al,qa,q(i),dusz,&
                 &qlb,ths(zo(ic)),thr(zo(ic)),qta,psih,ck,th,z)  
       nmn=1
       if(outp(8)) write(*,*) 'cell, time, t* ', ic,t0,tstar, ' Using ETS for basal flux' !Revised 2/2/2011 RLB
      else    
! later-time solution ...............    
       q2a=0.0
       do k=1,nmax
        qtop=qtops(k)*exp(-(r(k)**2)*tstar) !12/28/2010
        q2old=q2a
        q2a=q2a+qtop/bot(k)  !12/28/2010
        delq2=abs((q2a-q2old)/q2a)
        k1=k
        if(delq2<=1.0e-06) exit 
        if(abs(q2a)<=eps .and. k>3) exit 
       end do
       if(lcv) then
         if((delq2>1.0e-06) .and. (abs(q2a)>eps .and. k>3)) then
           ncc=ncc+1
           nv0=1
           lcv=.false.
         end if
       end if
       if(k1>nmax1) nmax1=k1
       if(k1<nmn) nmn=k1
       q0a=q(i)
       q1a=4.d0*(q(i)-qa)*exp(al*dusz/2.)
       qta=q0a-q1a*q2a*exp(-tstar/4.)
      end if
     else 
      qta=0.0
     end if
     tdif2=t0-capt(i+1)
     if(tdif2 > 0.0) then
      tstar=tf*tdif2
      if(tstar < smt) then
! Early-time solution ...............    
       z=dusz  ! inititalizations added 29 Jan 2013, RLB 
       qlb=ks(zo(ic))
! Corrected ths, thr to include index 12/20/2010, RLB
       call smallt(tstar,al,qa,q(i),dusz,&
                &qlb,ths(zo(ic)),thr(zo(ic)),qtb,psih,ck,th,z)  
       nmn=1
       if(outp(8)) write(*,*) 'cell, time, t* ', ic,t0,tstar, ' Using ETS for basal flux' !Revised 2/2/2011 RLB
      else
! later-time solution ...............    
       q2b=0.0
       do  k=1,nmax
        qtop=qtops(k)*exp(-(r(k)**2)*tstar) !12/28/2010
        q2old=q2b
        q2b=q2b+qtop/bot(k) !12/28/2010
        delq2=abs((q2b-q2old)/q2b)
        k1=k
        if(delq2<=1.0e-06) exit 
        if(abs(q2b)<=eps .and. k>3) exit 
       end do
       if(lcv) then
         if((delq2>1.0e-06) .and. (abs(q2b)>eps .and. k>3)) then
           ncc=ncc+1
           nv0=1
           lcv=.false.
         end if
       end if
       if(k1>nmax1) nmax1=k1
       if(k1<nmn) nmn=k1
       q0b=q(i)
       q1b=4.d0*(q(i)-qa)*exp(al*dusz/2.)
       qtb=q0b-q1b*q2b*exp(-tstar/4.)
      end if
     else
      qtb=0.0
     end if
     qt=qt+qta-qtb
    end do flux_time_loop
    if(qmax<qt .or. qt<0) then ! moved after end of loop, 7 Jan 2013
      dseep=ks(zo(ic))/(ths(zo(ic))*t0) !Added check of distance traversed at saturated seepage velocity 12/28/2010 RLB
      if(dseep < dusz) then
        qt=0.d0
      else
        write(*,*) 'Error computing basal flux!' !Revised 12/28/2010, RLB
        write(*,*) 'Cell, Depth, Time, t*, Max. Input flux, Basal flux: ',ic,dusz,t0,tstar,qmax,qt
        write(*,*) ''
      end if
    end if
    return
  end subroutine flux


!   SY: Subroutine indentation performed for integration into LIS
  subroutine smallt(tstar,alfa,qa,qb,d,&
        &qlb,ths,thr,qt,psih,ck,th,z)
! By Rex L. Baum and W.Z. Savage, USGS, October 1, 2004        
    integer,parameter:: double=kind(1d0)
    real (double), intent(in):: tstar,alfa,qa,d,qlb,z ! added intent 12/22/2010, RLB ,ths,thr,qb
    real, intent(in):: ths,thr,qb ! corrected to single 1/21/2011, 2/13/2013 RLB
    real (double), intent(out):: qt,psih,ck,th ! added intent 12/22/2010, RLB
    real (double):: ck0,psi0,th0,f0,ft,gt,arg1,arg2 
    real (double):: f1a,f1b,f1c,f1,f2a,f2b,f2c,f2
    real (double):: gt1,gt2,derfc
    ck0=qa-(qa-qlb)*exp(-alfa*(d-z))
    psi0=(log(ck0))/alfa
    th0=thr+(ths-thr)*exp(alfa*psi0)
    f0=2.*(qb-qa)*exp(alfa*z/2.)
    ft=0.0
    gt=0.0
    arg1=(alfa*z)/(2.*sqrt(tstar))
    f1a=derfc(arg1)
    f1b=exp(arg1*sqrt(tstar)+tstar/4.)
    f1c=derfc(arg1+sqrt(tstar/4.))
    f1=f1a-f1b*f1c
    arg2=alfa*(2.*d-z)/(2.*sqrt(tstar))
    f2a=derfc(arg2)
    f2b=exp(arg2*sqrt(tstar)+tstar/4.)
    f2c=derfc(arg2+sqrt(tstar/4.))
    f2=-f2a+f2b*f2c
    ft=f1+f2 
    gt1=f1b*f1c
    gt2=f2b*f2c
    gt=gt1-gt2
    ck=ck0+f0*ft
    if(ck < 0.0) ck=-ck
    psih=log(ck)/alfa
    th=thr+(ths-thr)*exp(alfa*psih)
    qt=f0*gt/2.0
    return
  end subroutine smallt


!   SY: Subroutine indentation performed for integration into LIS
  subroutine dsimps(n,h,f,intf)
! by Rex L. Baum 7 July 2004,	Latest revision 29 Jan 2013, RLB
! uses Simpson's 3-point rule to numerically integrate a function at n successive points
! on an interval that has been subdivided into 2n equal segments.
    integer, parameter:: double=kind(1d0)
    integer:: i,n
    real (double):: f(2*n+1)
    real (double):: h,intf(n+1),sumf,fl,fmid,fr ! n changed to n+1, 29 Jan 2013, RLB
    sumf=0.d0
    do i=1,n
      fl=f(2*i-1)
      fmid=f(2*i)
      fr=f(2*i+1)
      sumf=sumf+fl+4.d0*fmid+fr
      intf(i+1)=sumf*h/3.d0 ! value of integral over the interval f(1) and f(2*i+1)
    end do
    return
  end subroutine dsimps


!   SY: Subroutine indentation performed for integration into LIS
! compute K, Pressure head and volumetric water content in the unsaturated zone
! uses analytic solution for normalized transmissivity
! By Rex L. Baum and W.Z. Savage, USGS, latest revison 28 Jan 2013, RLB
  !subroutine unsth(i,j1,ncc,iper,t1,nmx1,lcv,ulog,vqt,delh,nmn1,sqin,vqtn) ! SY
  subroutine unsth(i,j1,ncc,iper,t1,nmx1,lcv,vqt,delh,nmn1,sqin,vqtn) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
        LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use model_vars
    use input_vars 
    use grids
    implicit none
    !integer::i,j,j1,k,k1,m,ncc,nmx1,nmn1,iper,ulog ! SY
    integer::i,j,j1,k,k1,m,ncc,nmx1,nmn1,iper!,ulog ! SY 
    logical lcv,lnu
    real (double):: t1,d,sqin 
    real (double):: al,cks,cths,cthr,tf,zns,qa,rslo,b1 
    real (double):: z,zinc,tstar,tdif1,tdif2,kas,kbs,ka,kb 
    real (double):: k2a,k2b,k2old,bot,top,delk,tol 
    real (double):: derfc,ck0,f0,ft,f1,f1a,f1b,f1c 
    real (double):: f2,f2a,f2b,f2c,arg1,arg2,cn,psi0,th0,tstar1
    real (double):: vfil,trvty,vqt,vqtn,vdif1 
    real (double):: vdif,vhi,zhi,vlo,zlo,vf(nzs+1),dusznew,delh
    real (double):: vfa(nzs+1),trsa,trsb,trva,trvb
    real (double):: trold,trbot,tra,trb,tr1,tr2,deltr,qlb
    real (double):: trtop,trtop1,trtop2,trtop3,ck1  
    tol=1.0e-06
    rslo=slo(i)
    b1=cos(rslo)
    al=alp(zo(i))
    al=al*b1*b1 ! for coordinate transformation
    cks=ks(zo(i))
    qlb=cks
    qa=cks*rikzero(i)
    cths=ths(zo(i))
    cthr=thr(zo(i))
    tf=al*cks/(cths-cthr)
    zns=float(nzs)
    zinc=(zmax(i)-zmin)/zns
    z=zmin
    pzero=0.d0
    d=dusz
    kz=0.d0;trz=0.d0
    Z_loop: do j=1,nzs+1
    !write(*,*) 'iper, z', iper,z  !d
     kas=0.d0;kbs=0.d0
     trsa=0.d0;trsb=0.d0
     if(z >= d) then 
       kz(j)=cks
       thz(j)=cths
       trz(j)=(zmax(i)-z)*cks
       z=z+zinc
       cycle
     end if
     ck0=qa-(qa-cks)*exp(-al*(d-z))  ! based on psi* = psi - psi0 (see Savage et al. 2004)
     psi0=log(ck0/cks)/al
     th0=cthr+(cths-cthr)*exp(al*psi0)
     f0=-(qts(j1)-qa)*exp(al*(d-z)/2.) ! changed q to qts() 21 Feb 2013, to solve array boundary error RLB 
     ft=0.0
     t_loop: do m=1,iper
       delk=0.d0;deltr=0.d0
       tdif1=t1-capt(m)
       if (tdif1 < 0.d0) exit ! jump out of loop rather than compute extra zeros.
       tstar=0.d0; tstar1=0.d0 ! added 2 Jan 2012 RLB
       if(tdif1 > 1.0e-10) then ! (1.0e-10 is effectively zero)
         tstar=tf*tdif1
         tstar1=tstar
         if(tstar < smt) then ! early-time solution
          cn=float(2*m-1)
          arg1=al*(cn*d-z)/(2.*sqrt(tstar))
          f1a=exp(-tstar/4.)*derfc(arg1)
          f1b=exp(arg1*sqrt(tstar))
          f1c=derfc(arg1+sqrt(tstar/4.))
          f1=f1a-f1b*f1c
          ka=ck0+2.*f1*f0 
          trva=ka*zinc ! approximate formula for  normalized transmissivity
         else ! later-time series solution
          k2a=0.d0; tra=0.d0
          k_series_a: do k=1,nmax
           top=sin(r(k)*al*(d-z))*sin(r(k)*al*d)*&
            &exp(-(r(k)**2)*tstar)
           bot=1.+al*d/2.+2.*al*d*r(k)**2
           k2old=k2a
           k2a=k2a+top/bot
           delk=abs((k2a-k2old)/k2a)
           trtop1=sin(al*r(k)*(d-z))+2.*r(k)*&
            &cos(r(k)*al*(d-z))
           trtop2=2.*r(k)*exp(al*d/2.)-exp(al*z/2.)*trtop1
           trtop3=sin(r(k)*al*d)*exp(-(r(k)**2)*tstar)
           trtop=trtop2*trtop3
           trbot=(r(k)**2+1./4.)*bot
           trold=tra
           tra=tra+trtop/trbot 
           deltr=abs((tra-trold)/tra)
           k1=k
           if(delk<=tol .and. deltr<=tol) exit k_series_a
          end do k_series_a
          if(lcv .and. delk>tol) then
            ncc=ncc+1
            nvu(i)=1
            lcv=.false.
            write(*,*) 'noncv_th1!', k1,t1,tdif1,delk,k2a
          end if
          f1=2.*(q(m)-qa)*exp(al*z/2.)*exp(-tstar/4.)
          ck1=q(m)-(q(m)-qlb)*exp(-al*(d-z)) !  based on psi* = psi - psi0 (see Savage et al. 2004)
          ka=ck1-2.*f1*k2a
          tr1=q(m)*(d-z)-(q(m)-qlb)*(1.0-exp(-al*(d-z)))/al
          tr2=2.*(q(m)-qa)*exp(-tstar/4.)/al
          trva=tr1-tr2*tra
         end if
         ! SY: Begin making sure k1 is used only when defined
         if(tstar >= smt) then ! later-time series solution
          if(k1>nmx1) nmx1=k1
          if(k1<nmn1) nmn1=k1
         end if
         ! SY: End making sure k1 is used only when defined
       else
         ka=0.0
         f1=0.0
         trva=0.0
       end if
       kas=kas+ka
       trsa=trsa+trva
!    write(*,*) 'j1,m,t1,tdif1,tstar,tstar1,kas,trsa', j1,m,t1,tdif1,tstar,tstar1,kas,trsa !d
       tdif2=t1-capt(m+1)
       tstar=0.d0!; tstar1=0.d0 ! added 2 Jan 2012 RLB
       if(tdif2 > 1.0e-10) then ! (1.0e-10 is effectively zero)
         tstar=tf*tdif2
         if(tstar < smt) then ! early-time solution
          cn=float(2*m-1)! added 28 Jan 2013, RLB 
          arg2=al*(cn*d+z)/(2.*sqrt(tstar))
          f2a=exp(-tstar/4.)*derfc(arg2)
          f2b=exp(arg2*sqrt(tstar))
          f2c=derfc(arg2+sqrt(tstar/4.))
          f2=f2a-f2b*f2c
          kb=ck0+2.*f2*f0 
          trvb=ka*zinc ! approximate formula for normalized transmissivity
         else ! later-time series solution
          k2b=0.d0; trb=0.d0
          k_series_b: do k=1,nmax
           top=sin(r(k)*al*(d-z))*sin(r(k)*al*d)*&
            &exp(-(r(k)**2)*tstar)
           bot=1.+al*d/2.+2.*al*d*r(k)**2
           k2old=k2b
           k2b=k2b+top/bot
           delk=abs((k2b-k2old)/k2b)
           trtop1=sin(al*r(k)*(d-z))+2.*r(k)*&
            &cos(r(k)*al*(d-z))
           trtop2=2.*r(k)*exp(al*d/2.)-exp(al*z/2.)*trtop1
           trtop3=sin(r(k)*al*d)*exp(-(r(k)**2)*tstar)
           trtop=trtop2*trtop3
           trbot=(r(k)**2+1./4.)*bot
           trold=trb
           trb=trb+trtop/trbot 
           deltr=abs((trb-trold)/trb)
           k1=k
           if(delk<=tol .and. deltr<=tol) exit k_series_b
          end do k_series_b
          if(lcv .and. delk>tol) then
           ncc=ncc+1
           nvu(i)=1
           lcv=.false.
           write(*,*) 'noncv_th2!', k1,t1,tdif2,delk,k2b
          end if
          f1=2.*(q(m)-qa)*exp(al*z/2.)*exp(-tstar/4.)
          ck1=q(m)-(q(m)-qlb)*exp(-al*(d-z)) ! based on psi* = psi - psi0 (see Savage et al. 2004)
          kb=ck1-2.*f1*k2b
          tr1=q(m)*(d-z)-(q(m)-qlb)*(1.0-exp(-al*(d-z)))/al
          tr2=2.*(q(m)-qa)*exp(-tstar/4.)/al
          trvb=tr1-tr2*trb
         end if
         ! SY: Begin making sure k1 is used only when defined
         if(tstar >= smt) then ! later-time series solution
          if(k1>nmx1) nmx1=k1
          if(k1<nmn1) nmn1=k1
         end if
         ! SY: End making sure k1 is used only when defined
       else
         kb=0.0
         f2=0.0
         trvb=0.0
       end if
       kbs=kbs+kb
       trsb=trsb+trvb
!        write(*,*) 'j1,m,t1,tdif2,tstar,tstar1,kbs,trsb', j1,m,t1,tdif2,tstar,tstar1,kbs,trsb !d
     end do t_loop  
     if(t1==0.d0) then
       kz(j)=ck0
       ptran(j)=psi0*b1
! Analytic formulas for normalized transmissivity, trz()          
       if(lps0) then 
         trz(j)=qa*(d-z)+(qa-cks)*(exp(-al*(d-z))-1)/al ! use if psi0=-1/al,   see Savage et al. 2004
       else
         trz(j)=qa*(d-z)+(qa-cks)*(exp(-al*(d-z))-1)/al ! use if psi0=0
       end if
     else
       kz(j)=kz(j)+kas-kbs
       trz(j)=trz(j)+trsa-trsb
     end if
! confirm that ptran(j)<0
     if(kz(j)>0.d0) then
       if(lps0) then ! multipliying by b1* corrects for slope & preserves water content profile
         ptran(j)=b1*(log(kz(j)/cks)-1.)/al ! use if psi0=-1/al
       else
         ptran(j)=b1*(log(kz(j)/cks))/al ! use if psi0=0
       end if
     else
       ptran(j)=-b1/al 
     end if  
     thz(j)=cthr+(cths-cthr)*exp(al*ptran(j)) 
     if (thz(j)>cths) thz(j)=cths
     z=z+zinc
    end do Z_loop
    z=zmax(i)-zinc
    vfil=0.d0;vf(nzs+1)=0.d0;vfa(nzs+1)=0.d0
! Numerical integration of normalized transmissivity (trvty)
    trvty=-(zmax(i)-d)*cks 
    do j=nzs,1,-1
      trvty=trvty+zinc*(kz(j)+kz(j+1))/2. 
      if(tstar1 < smt) then
        trz(j)=trvty
      end if
      if(z>=d) then
        vf(j)=0.d0
      else
        vf(j)=(cths-cthr)*((d-z)-trvty/cks) 
        if(vf(j)<0.d0) vf(j)=-vf(j)
      end if
      if(z>d) then
        vfa(j)=0.d0
      else
        vfa(j)=(cths-cthr)*((d-z)-trz(j)/cks) 
        if(vfa(j)<0.d0) vfa(j)=-vfa(j)
      end if
      z=z-zinc
    end do
    lnu=.true.
    zhi=0.d0; zlo=zmax(i) ! initalization to default values added 28 Jan 2013, RLB
    if(lany) vf=vfa
    if(t1==0.d0) vf0=vf(1)
    if(vqt>0.d0) then
      z=zmin; vlo=0.d0; vhi=cths*zmax(i)
      do j=1,nzs 
        if(vf(j)>0.d0) then
          vdif=vf(j)-vqt
          if(vdif==0.d0) then
            dusznew=z
            exit
          else if(vdif>0.d0) then
            vdif1=vf(j+1)-vqt
            if(vdif1<0) then
              vhi=vf(j)
              zhi=z
              vlo=vf(j+1)
              zlo=z+zinc
            end if
          end if  
        end if
        z=z+zinc
      end do
      if(vqt>vf(1)) then
        zlo=zmax(i)
        zhi=0.d0
        vhi=vqt
        vlo=0.d0
      end if
    else
      lnu=.false.
    end if
    delh=0.d0
    if(llus) then ! estimate water-table rise
      if(lnu) then
        if(zlo<d) then
          dusznew=zhi+zinc*(vhi-vqt)/(vhi-vlo)
        else
          dusznew=d-zinc*(vqt)/(vhi-vlo)
        end if
        delh=d-dusznew
        if(depth(i)>zmax(i)) delh=delh-(depth(i)-zmax(i)) ! accounts for initial wt below zmax(i)
        if(vqt>vf(1)) then
          dusznew=0.d0
          delh=depth(i)-dusznew
        end if
      else 
        delh=0.d0
        dusznew=d-delh
      end if
    end if
!if(outp(8)) write (ulog,*) t1,sqin,t1*rizero(i),vqtn,vf0-vf(1),sqin-vqtn-t1*rizero(i),vqt,delh ! SY: See next line
!if(outp(8)) write (ulog,'(1G,7ES24.15E3)') t1,sqin,t1*rizero(i),vqtn,vf0-vf(1),sqin-vqtn-t1*rizero(i),vqt,delh ! SY: Change again to below
        !SY : start formatting
    if (t1 .eq. 0) then
     if(outp(8)) write (LIS_logunit,'(3X,1F17.15,7ES24.15E3)') t1,sqin,t1*rizero(i),vqtn,vf0-vf(1),sqin-vqtn-t1*rizero(i),vqt,delh ! SY
    else
     if(outp(8)) write (LIS_logunit,'(3X,1F16.9,7ES24.15E3)') t1,sqin,t1*rizero(i),vqtn,vf0-vf(1),sqin-vqtn-t1*rizero(i),vqt,delh ! SY
    end if
        !SY : stop formatting
    return
  end subroutine unsth


!   SY: Subroutine indentation performed for integration into LIS
!  diffusion model for pressure head from applied load of rising water table
!	By R.L. Baum and W.Z. Savage, USGS, May 2004, latest revision 29 Jan 2013, RLB 
  subroutine pstpi(u1,dhwt,dwt,&
        & i,iper,rf,t0,jf) ! 1/29/2013 RLB ! SY
        !& ulog,i,iper,rf,t0,jf) ! 1/29/2013 RLB ! SY
    use grids
    use input_vars
    use model_vars
    implicit none
    !integer:: i,j,jf,u1,ulog,m,iper,jmark,nstp,jtop !,n1,n,nccs,nmx ! SY
    integer:: i,j,jf,u1,m,iper,jmark,nstp,jtop !,n1,n,nccs,nmx ! SY
    real:: dwt!,chi 12/22/2010 made chi a 1-d array & moved to model_vars
    real (double):: finf,t1,t2,term1,term2,znew,t0 !,rn
    real (double):: ferfc1,ferfc3,dhwt(nts+1),dh
    real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
    real (double):: ar1,ar3,fs,rf(nzs+1) !,ar2
    real (double):: rfa,rfb,derfc,ff,rslo,rphi,fmn,ptest
    real (double):: tol !,delt1,delt2,t1old,t2old
    real (double):: ddg2rad,dusz0,dusz1,dlz,d1,zm
    real (double):: uwt1,uwsum !,flt1,flt2,tstar
    pi=3.141592653589793
    ddg2rad=pi/180.D0
    !  maximum value of Factor of Safety
    finf=10.
    nmn=1+mmax
    tol=1.e-06
    dh=dhwt(iper+1)
    rslo=slo(i)
    if (flag<=-1 .and. flag>=-3) then !12/14/10 added flag>=-3
      write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
      & 'cell ',i,rslo/ddg2rad,iper+1,t0 ! Added label 4/22/2010 RLB
    end if
    rphi=phi(zo(i))
    a1=sin(rslo)
    b1=cos(rslo)
    select case (flowdir) ! set value of beta (Iverson's beta line)
      case ('slope')
       beta=b1*b1
      case ('hydro')
       beta=1.d0
      case default
       beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
    end select
    if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
    if (abs(rslo)>1.e-5) then
      ff=tan(rphi)/tan(rslo)
    else !  set factor of safety to fixed value for flat slopes   
      ff=finf
    end if
    zns=float(nzs)
    zinc=(zmax(i)-zmin)/zns
    nstp=int((1./alp(zo(i))/zinc))
    z=zmin
    fmn=1.e25
  ! compute depth of unsaturated zone
    if(unsat(zo(i))) then
      if(lps0) then
        dusz0=depth(i)-1./alp(zo(i)) ! use when psi0=-1/alp(zo(i))
      else
        dusz0=depth(i) ! use if psi0=0
      end if
      dusz1=dwt-1./alp(zo(i))
    else
      dusz0=0.  
      dusz1=0.  
    end if
    dlz=zmax(i)-dusz0
    zm=zmax(i)
    uwsum=0. ! for computing depth-averaged unit weight
    uwsp=uws(zo(i))
    rf=0.0 ! initialize rf
    jmark=1 ! Added in case depth(i)==0., 28 Jan 2013, RLB 
    Z_loop: do j=1,nzs+1
     znew=z
     if(z<depth(i))then 
       jmark=j
     end if
     if(znew < 1.0e-30) znew =1.0e-30
     pzero(j)=beta*(z-depth(i))
     if(z <= depth(i))then 
       pzero(j)=0.
       if(lps0) then ! psi0=-1/alpha
         if(llus) then ! compute water table rise  
           if(z > dusz0-dh .and. z < depth(i)-dh) then  
             rf(j)=(-1./alp(zo(i))+(z-dusz0+dh))
           else if(z >= depth(i)-dh) then  
             rf(j)=beta*(z-depth(i)+dh)
           else 
             rf(j)=ptran(j)
           end if
         else ! static water table
           if(z > dusz0 .and. z < depth(i)) then  
             rf(j)=(-1./alp(zo(i))+(z-dusz0))
           else if(z >= depth(i)) then  
             rf(j)=beta*(z-depth(i))
           else
             rf(j)=ptran(j) 
           end if
         end if
       else ! psi0=0
         if(z >= dusz0-dh .and. llus) then ! water table rise
           rf(j)=beta*(z+dh-dusz0)
         else ! static water table
           rf(j)=ptran(j)
         end if
       end if
     else if(z>=zm .and. dlz<0.001) then
       rf(j)=beta*(z+dh-dusz0)      
     else if(llus) then ! compute pressure rise below wt only if llus=.true.
       temporal_loop: do m=1,iper
         tdif1=t0-tcap(m)
         if(tdif1 > 0.0) then
 ! corrected diffusivity term in next line (divide by b1*b1) 
           d1=dif(zo(i))/(b1*b1)
           t1=sqrt(d1*tdif1)
           if (t1<1.0e-29) t1=1.0e-29
           term1=0.0
 ! Solution for all times
           ar1=(z-dusz0)/(2.*t1)
           ferfc1=derfc(ar1) 
           term1=ferfc1
           rfa=term1
         else
           rfa=0.0
         end if
         tdif2=t0-tcap(m+1)
         if(tdif2 > 0.0) then
 ! corrected diffusivity term in next line (divide by b1*b1)       
           d1=dif(zo(i))/(b1*b1)
           t2=sqrt(d1*tdif2)
           if (t2<1.0e-29) t2=1.0e-29
           term2=0.0 
 ! Solution for all times
           ar3=(z-dusz0)/(2.*t2)
           ferfc3=derfc(ar3)
           term2=ferfc3
           rfb=term2
         else
           rfb=0.0
         end if
         rf(j)=rf(j)+dh*beta*(rfa-rfb) ! slope & Hyd. grad. correction
       end do temporal_loop 
     end if
     bline(j)=z*beta
     ptran(j)=rf(j)
     p(j)=pzero(j)+ptran(j)
     ptest=p(j)-bline(j)
     if(ptest > 0.0) then
       p(j)=bline(j)
     end if
 ! estimate partially saturated unit weight      
     if(p(j)<-1.d0/alp(zo(i))) then
       uwt1=(gs(zo(i))*(1-ths(zo(i)))+thz(j))*uww
     else
       uwt1=uws(zo(i))
     end if
     uwsum=uwsum+uwt1
     uwsp(j)=uwsum/float(j)
     z=z+zinc
    end do Z_loop  
! find new height of rising water table in zones of upward seepage   
    if(rikzero(i)<0.0) then
     zinc=(zmax(i)-zmin)/zns
     z=zmin
     newdep=0.0
     do j=1,nzs+1
      if(p(j)<0.0) newdep=z
      z=z+zinc
     end do
! adjust presures 
     z=zmin
     do j=1,nzs+1
      if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
      if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
      z=z+zinc    
     end do
    end if
! Smooth data near initial water table
    if(llus .and. zm>depth(i)) then
      jtop=jmark-nstp
      if(jtop<1) jtop=1
      do j=jmark,jmark-5,-1
        if(ptran(j)>ptran(j+1)) then 
          ptran(j)=ptran(j+1)
          p(j)=ptran(j)+pzero(j)
        end if
      end do
    end if  
! Compute factor of safety & save results    
    z=zmin
    Z_FS_loop: do j=1,nzs+1    
     chi(j)=1.0 
 ! Approximate suction stress computation, moved 4/22/2010 RLB
     if(p(j)<0.) then ! for suction 2/25/2009 (formerly for suctions exceeding the air-entry value, -1/alpha)
       chi(j)=(thz(j)-thr(zo(i)))/(ths(zo(i))-thr(zo(i)))
     else
       chi(j)=1.0 
     end if
     if (abs(a1)>1.e-5 .and. z>0.) then
 ! Approximate suction stress computation    
       fw(j)=-(chi(j)*p(j)*uww*tan(rphi))/(uwsp(j)*z*a1*b1) ! same as saturated formula with factor Chi
       fc(j)=c(zo(i))/(uwsp(j)*z*a1*b1) ! uses depth averaged unit weight
     else
       fw(j)=0.d0
       fc(j)=0.d0
     end if
     fs=ff+fw(j)+fc(j)
 ! frictional strength cannot be less than zero (correction added 7/12/02)
     if ((ff+fw(j))<0.) fs=fc(j)
     if (fs>finf) fs=finf
     if (z<=1.e-02) fs=finf 
     if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
     if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),&
        & pzero(j),ptran(j),bline(j),fs
     if (flag==-3) write(u1,'(6(g12.5,1x):)') z,p(j),fs,chi(j) ! added 4/14/2010 RLB
     if (fs<fmn) then
       fmn=fs
       if (jf>0) then
         zfmin(i+(jf-1)*imax)=z
         pmin(i+(jf-1)*imax)=p(j)
       end if
     end if
     z=z+zinc
    end do Z_FS_loop
    if (jf>0) then
      fsmin(i+(jf-1)*imax)=fmn
      if(fmn==finf) then ! Added 30 Jan 2013, RLB 
        pmin(i+(jf-1)*imax)=p(nzs+1)
        zfmin(i+(jf-1)*imax)=zmax(i)
      end if 
      if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
        !if (flag>=-6) call svijz(i,jf,dh,newdep,ulog) ! SY
        if (flag>=-6) call svijz(i,jf,dh,newdep) ! SY
        !if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,dh,newdep,ulog)  ! Added 2/10/2012 ! SY
        if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,dh,newdep)  ! Added 2/10/2012 ! SY
      end if
    end if
    return
  end subroutine pstpi


!   SY: Subroutine indentation performed for integration into LIS
!  Implementation of Iverson's (2000) method of computing pore pressure
!  for rain infiltration.
!  by W.Z. Savage, spring 2001, with modifications by R.L. Baum
!  (both) USGS, Latest revsion 29 Jan 2013, RLB 
!  
  !subroutine ivestp(u1,rikf,ulog,i,rf) ! SY
  subroutine ivestp(u1,rikf,i,rf) ! SY
    use grids
    use input_vars
    use model_vars
    implicit none
    !integer:: i,j,jf,u1,ulog,nmx ! ,nccs ! SY
    integer:: i,j,jf,u1,nmx ! ,nccs ! SY
    integer:: n,nn
    real:: rikf(nts+1),finf 
    real (double):: derfc,a1,b1,ff,zns,zinc,z,t0,znew
    real (double):: zstar,tstar,x1,x2,x3,x4
    real (double):: rf1,rf2,rf3,rf4,rfa,rfb,rf(nzs+1)
    real (double):: fs,rslo,rphi,fmn,ptest,pmn,dhat ! pmn added 4/15/2010, dhat added 8 Jan 2013
    real (double):: captstar1,captstar2,tdif1,tdif2
    real (double):: dusz1,ddg2rad,newdep
    !logical:: lcvs 
    pi=3.141592653589793
    ddg2rad=pi/180.D0
    !  maximum value of Factor of Safety
    finf=10.
    nmx=0
! nmn=1+mmax ! initialization moved to main program 07 Jan 2013 RLB 
    rslo=slo(i) 
    rphi=phi(zo(i))
    a1=sin(rslo)
    b1=cos(rslo)
    dhat=4.*dif(zo(i))/(b1*b1) ! added 8 Jan 2013, RLB
    p0zmx=0. ! Added 6 May 2013, RLB 
    select case (flowdir) ! set value of beta (Iverson's beta line)
      case ('slope')
        beta=b1*b1
      case ('hydro')
        beta=1.d0
      case default
        beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
    end select
    if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
    if (abs(rslo)>1.e-5) then
      ff=tan(rphi)/tan(rslo)
    else !  set factor of safety to fixed value for flat slopes
      ff=finf
    end if
    zns=float(nzs)
    zinc=(zmax(i)-zmin)/zns
    dusz1=0.
    ptran=0. ! added 17Aug2009 RLB
    temporal_loop: do n=1,nts+1
      t0=tcap(n)
      fmn=1.e25
      jf=jsav(n)
      if (flag<=-1 .and. flag>=-3) then !12/14/10 added flag>=-3
        write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
        'cell ',i,rslo/ddg2rad,n,t0 ! added label 4/22/2010 RLB
      end if
      z=zmin
      Z_loop: do j=1,nzs+1
        znew=z
        if(znew < 1.0e-30) znew =1.0e-30
        if (abs(a1)>1.e-5) then
          fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
        else
          fc(j)=0.d0
        end if
        pzero(j)=beta*(z-depth(i))
        if(z<dusz1 .or. rikf(n)==0.) then
          rf(j)=0.0
        else
          rf(j)=0.0
          if (abs(z)>0.) then ! Formulas in next 2 lines apply only for z>0, test added 21 Feb 2013, RLB
            zstar=z**2/dhat ! formula simplified 8 Jan 2013, RLB
            tstar=t0/zstar
          end if
          temporal_loop_1: do nn=1,nper
            if(z==0.) then ! exact formula added for case of z=0, 8 Jan 2013, RLB
              tdif1=t0-capt(nn)
              if(tdif1>0.) then
                rfa=sqrt(tdif1*dhat/pi)
              else
                rfa=0.0
              end if
              tdif2=t0-capt(nn+1)
              if(tdif2>0.) then
                rfb=sqrt(tdif1*dhat/pi)
              else
                rfb=0.0
              end if
            else ! z>0
              captstar1=capt(nn)/zstar 
              tdif1=tstar-captstar1
              if(tdif1 > 0.0) then 
                x1=1./tdif1
                x2=1./(sqrt(tdif1))
                rf1=sqrt(1./(x1*pi))*exp(-x1)
                rf2=derfc(x2)
                rfa=rf1-rf2
              else
                rfa=0.0
              end if
              captstar2=capt(nn+1)/zstar
              tdif2=tstar-captstar2
              if(tdif2 > 0.0) then
                x3=1./tdif2
                x4=1./(sqrt(tdif2))
                rf3=sqrt(1./(x3*pi))*exp(-x3)
                rf4=derfc(x4)
                rfb=rf3-rf4
              else
                rfb=0.0
              endif
            endif
            rf(j)=rf(j)+rik(i+(nn-1)*imax)*(rfa-rfb)
    !         if (z==0) write(*,*) 'i, j, rfa, rfb, rf(j), tdif1, tdif2',  i, j, rfa, rfb, rf(j), tdif1, tdif2
          end do temporal_loop_1
        end if
        bline(j)=z*beta
        if(abs(rf(j))>0.0) then ! added 17AUG2009 RLB 
          ptran(j)=z*rf(j) ! correction, added "z*" 4/27/2012, RLB
          if(z==0) ptran(j)=rf(j) ! formula for ptran at z=0 not normalized, 9 Jan 2013, RLB 
        end if
        p(j)=pzero(j)+ptran(j)
        ptest=p(j)-bline(j)
        if(ptest > 0.0) then
          p(j)=bline(j)
        end if
        z=z+zinc
      end do   Z_loop
      if(n==1) p0zmx=p(nzs+1) ! Added 6 May 2013, RLB 
    ! find new height of rising water table in zones of upward seepage   
      if(rikzero(i)<0.0) then
        zinc=(zmax(i)-zmin)/zns
        z=zmin
        newdep=0.0
        do j=1,nzs+1
          if(p(j)<0.0) newdep=z
          z=z+zinc
        end do
    ! adjust presures 
        z=zmin
        do j=1,nzs+1
          if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
          if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
          z=z+zinc    
        end do
      end if
    ! Compute factor of safety & save results    
      z=zmin
      Z_FS_loop: do j=1,nzs+1 
        if (abs(a1)>1.e-5 .and. z>1.e-30) then
          if(lpge0 .and. p(j)<0.) then !option added 4/15/2010
            fw(j)=0.d0
          else
            fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
          end if
        else
          fw(j)=0.d0
        end if
        fs=ff+fw(j)+fc(j)
    ! frictional strength cannot be less than zero 
        if ((ff+fw(j))<0.) fs=fc(j)
        if (fs>finf) fs=finf
        if (z<=1.e-02) fs=finf 
        if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
        if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j)&
               & ,pzero(j),ptran(j),bline(j),fs
        if (flag==-3 .and. unsat0) then ! revised 12/23/2010
          write(u1,'(6(g12.5,1x):)') z,p(j),fs,1. ! added 4/14/2010 RLB
        else if (flag==-3) then
          write(u1,'(6(g12.5,1x):)') z,p(j),fs
        end if
        if(jf>0) then
          if (fs<fmn) then
            fmn=fs
            zfmin(i+(jf-1)*imax)=z
            pmn=p(j) ! revised 4/15/2010
          end if
        end if
        z=z+zinc
      end do Z_FS_loop
      if (jf>0) then  ! revised 4/29/2010 to include pmin() RLB
        fsmin(i+(jf-1)*imax)=fmn
        if(fmn==finf) then ! Added 30 Jan 2013, RLB 
          pmn=p(nzs+1)
          zfmin(i+(jf-1)*imax)=zmax(i)
        end if 
        if(lpge0 .and. pmn<0.) then !option added 4/15/2010
          pmin(i+(jf-1)*imax)=0.
        else
          pmin(i+(jf-1)*imax)=pmn
        end if
        if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
          dcf=0. ! this line used only for saturated infiltration model
          chi=1.d0 !added 12/23/2010
          !if (flag>=-6) call svijz(i,jf,0.d0,newdep,ulog) ! SY
          if (flag>=-6) call svijz(i,jf,0.d0,newdep) ! SY
          !if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,0.d0,newdep,ulog)  ! Added 2/10/2012 ! SY
          if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,0.d0,newdep)  ! Added 2/10/2012 ! SY
        end if
      end if
    end do temporal_loop
    return
  end subroutine ivestp


!   SY: Subroutine indentation performed for integration into LIS
  !subroutine unsfin(imx1,ulog,u1,ncc,nccs) ! SY
  subroutine unsfin(imx1,u1,ncc,nccs) ! SY
! By Rex L. Baum, USGS	Latest revision, 29 Jan 2013 
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    use grids
    use input_vars 
    use model_vars
    use input_file_defs
    implicit none
    !integer::i,j,jf,k,ulog,u1,imx1,ncc,nccs,nmx,nmax0 ! SY
    integer::i,j,jf,k,u1,imx1,ncc,nccs,nmx,nmax0 ! SY
    integer::nmn1,nmin1,nmax3,nmxp,nmxs,nmnp,nmns
    integer:: svgctr
    logical:: lcv,lcvs,lwt
    real:: delwt,dwt,zwt,qbij(nts+1) 
    real::testqk,tolqk ! Added 2/2/2011 RLB
    real (double)::rf(nzs+1),finf,vqt,qta,al,qzmax 
    real (double)::ddwt,sqin,intq(nts+1),b,dhwt(nts+1),delh 
    real (double)::qtn(2*nts+1),intq1(nts+1),vqtn,cd
    nmax3=0;nmax0=0
    nmn1=nmax+1;nmin1=nmax+1
    nmxp=0; nmnp=0; svgctr=0 ! Added 29 Jan 2013, RLB 
    intq=0.d0; intq1=0.d0 ! Added 29 Jan 2013, RLB 
    write(LIS_logunit,*) 'Starting coupled saturated & unsaturated zone'
    write(LIS_logunit,*) 'computations for finite-depth saturated zone'
    write(*,*) 'Starting coupled saturated & unsaturated zone'
    write(*,*) 'computations for finite-depth saturated zone'
    write(*,*) 'Cells completed: '
    ! loop over all grid cells
    finf=10.
    grid_loop: do i=1,imx1 
     if (mod(i-1,2000)==0) write (*,*) i-1 ! cells completed
     if(slo(i)<slomin) then ! default values for gently sloping cells 
       do jf=1,nout
         fsmin(i+(jf-1)*imax)=finf+1.
         zfmin(i+(jf-1)*imax)=zmax(i)
         pmin(i+(jf-1)*imax)=0.
       end do
       cycle
     end if
     lcv=.true.;lcvs=.true.
     q=0.;qb=0 ! qb initialization added 29 Jan 2013, RLB 
     tolqk=ks(zo(i))*5.e-07 ! Moved 29 Jan 2013, RLB 
     do j=1,kper
       if(j>nper) then
         q(j)=0.
       else
! 8/18/2009 RLB added optional offset of background flux to prevent excessive drying during periods of zero infiltration.      
        if(bkgrof) then
          q(j)=ks(zo(i))*(rik(i+(j-1)*imax)+rikzero(i))
! RLB 2/2/2001 revised test
          testqk=q(j)-(ks(zo(i))+rizero(i))
!          tolqk=ks(zo(i))*5.e-07
          !if(testqk>tolqk) write (LIS_logunit,*) '*q>Ks+ri!', i,j,q(j),ks(zo(i))+rizero(i) ! Added 2/2/2011 RLB ! SY: See next line
          if(testqk>tolqk) write (LIS_logunit,'(A,2I3,2ES14.6E2)') '*q>Ks+ri!', i,j,q(j),ks(zo(i))+rizero(i) ! Added 2/2/2011 RLB ! SY
        else
          q(j)=ks(zo(i))*rik(i+(j-1)*imax)
          testqk=q(j)-ks(zo(i))
          tolqk=ks(zo(i))*5.e-07          
          !if(testqk>tolqk) write (LIS_logunit,*) '*q>Ks!', i,j,q(j),ks(zo(i)) ! Moved 2/2/2011 RLB ! SY: See next line
          if(testqk>tolqk) write (LIS_logunit,'(A,2I3,2ES14.6E2)') '*q>Ks!', i,j,q(j),ks(zo(i)) ! Moved 2/2/2011 RLB ! SY
  
        end if
       end if
     end do
     qmax=maxval(q)
     b=cos(slo(i))
     if(unsat(zo(i))) then
       dcf=depth(i)-1.d0/(alp(zo(i)))
     else
       dcf=0.
     end if 
     if(lps0 .and. unsat(zo(i))) then
       dusz=depth(i)-1.d0/(alp(zo(i)))
     else 
       dusz=depth(i)
     end if
! set value of beta (Iverson's beta line)
! and maximum drainage rate at water table, qzmax    
     cd=0.1d0 !partial drainage for low-permeability basal boundary    
     select case (flowdir)         
       case ('slope')
        beta=b*b
        qzmax=(1.d0-beta)*cd*ks(zo(i)) 
       case ('hydro')
        beta=1.d0
        qzmax=0.d0 
       case default
        beta=b*b-rikzero(i)  ! 2/12/09 corrected formula for beta & qzmax
        qzmax=(1.d0-beta)*cd*ks(zo(i))-cd*rizero(i) 
     end select
     delwt=0.;dwt=depth(i);zwt=depth(i);ddwt=depth(i)
     lwt=.false.
     ts=tmin
     if(dcf>0. .and. (unsat(zo(i)) .or. igcap(zo(i)))) then ! updated to enforce non-zero depth, 4/18/2013, RLB
! compute flux and pore-pressure rise at each time step   
      vqt=0.;vqtn=0.;qta=rizero(i);sqin=0.
      al=alp(zo(i))*b*b
      if(outp(8)) write (LIS_logunit,*) 'ts,    qt '! times and basal flux to log file
      call roots(nmax,r,dusz,al,eps,pi)
      flux_loop: do j=1,2*nts+1
       call flux(i,kper,ts,j,lcv,ncc,nvu(i),lwt) 
       if(nmax1>nmax2) nmax2=nmax1
       if(nmn<nmin) nmin=nmn
       if(qt<rizero(i)) qt=rizero(i)
! RLB 2/3/2011, Added case for bkgrof=.true.               
       if(bkgrof) then
         if(qt>ks(zo(i))+rizero(i)+tolqk) then
           write(LIS_logunit,*) 'Error! Basal flux exceeds Ks at'
           write(LIS_logunit,*) 'cell ',i, ', timestep ',j 
           write(LIS_logunit,*) 'flux ',qt, ', Ks ',ks(zo(i))+rizero(i) 
           qt=ks(zo(i))
         end if
       else
         if(qt>ks(zo(i))+tolqk) then
           write(LIS_logunit,*) 'Error! Basal flux exceeds Ks at'
           write(LIS_logunit,*) 'cell ',i, ', timestep ',j 
           write(LIS_logunit,*) 'flux ',qt, ', Ks ',ks(zo(i))
           qt=ks(zo(i))
         end if
       end if
! RLB
       if(outp(8)) write (LIS_logunit,*) ts,qt ! times and basal flux to log file
! drain off excess basal flux 
       if(qt>qzmax) then
         qtime(j)=qt-qzmax
       else
         qtime(j)=0.d0
       end if
       qtn(j)=qt
       qta=qt
       ts=ts+tinc/2.d0
      end do flux_loop
      call dsimps(nts,tinc/2.d0,qtime,intq)
      call dsimps(nts,tinc/2.d0,qtn,intq1)
      !if(outp(8)) write (LIS_logunit,*) 'Time, Cumulative volume in, Cumulative background flux,& ! SY: Changed to line below
      if(outp(8)) write (LIS_logunit,'(500A)') 'Time, Cumulative volume in, Cumulative background flux,& 
              & Cumul. volume out,Cuml. absorbed, Cuml. qin-qout, qout not drained, Water table rise'
      ts=tmin
      wt_rise_loop: do j=1,nts+1
        jf=jsav(j)
        rf=0.0 
        vqt=intq(j)-ts*rizero(i)
        if(vqt<0.) vqt=0.d0
        vqtn=intq1(j)-ts*rizero(i)
        sqin=0.
        sum_q_in_loop: do k=1,nper
          if(ts>capt(k) .and. ts<=capt(k+1)) then
                qts(j)=q(k) !; write(*,*) 'j,k,qts(j),q(k),ts,capt(k) ', j,k,qts(j),q(k),ts,capt(k) ! Added 21 Feb 2013, RLB 
          endif
          if(ts>=capt(k+1)) then
            sqin=sqin+(capt(k+1)-capt(k))*q(k)
          end if
          if(ts>capt(k) .and. ts<capt(k+1)) then
            sqin=sqin+(ts-capt(k))*q(k)
          end if
        end do sum_q_in_loop
        if(jf>0 .and. lskip) then
          call unsth(i,j,ncc,kper,ts,nmax0,&
          &lcv,vqt,delh,nmn1,sqin,vqtn) ! SY
          !&lcv,LIS_logunit,vqt,delh,nmn1,sqin,vqtn) ! SY
          dhwt(j)=delh
          if(nmax0>nmax3) nmax3=nmax0
          if(nmn1<nmin1) nmin1=nmn1
        else if(lskip) then
          continue
        else
          call unsth(i,j,ncc,kper,ts,nmax0,&
          &lcv,vqt,delh,nmn1,sqin,vqtn) ! SY
          !&lcv,LIS_logunit,vqt,delh,nmn1,sqin,vqtn) ! SY
          dhwt(j)=delh
          if(nmax0>nmax3) nmax3=nmax0
          if(nmn1<nmin1) nmin1=nmn1
        end if        
        tcap(j)=ts ! pass to diffusion subroutine
        tcap(j+1)=ts+tinc
        ts=ts+tinc
        if(jf>0 .and. lskip) then
! compute pressure diffusion in saturated zone
          call pstpf(u1,dhwt,dwt,&
                  & i,j-1,rf,nccs,lcvs,nmx,tcap(j),jf) ! SY
                  !& LIS_logunit,i,j-1,rf,nccs,lcvs,nmx,tcap(j),jf) ! SY
          nmxp=nmx;nmnp=nmn
          if(j>1) then
! Check change in water table depth and adjust dusz if needed
            delwt=abs(dwt-zwt)*1000.
            if(delwt>dwt) then
              lwt=.true.
              dwt=zwt
            end if
          end if   
        else if(lskip) then
          continue
        else
! compute pressure diffusion in saturated zone
          call pstpf(u1,dhwt,dwt,&
                  & i,j-1,rf,nccs,lcvs,nmx,tcap(j),jf) ! SY
                  !& LIS_logunit,i,j-1,rf,nccs,lcvs,nmx,tcap(j),jf) ! SY
          nmxp=nmx;nmnp=nmn
          if(j>1) then
! Check change in water table depth and adjust if needed
            delwt=abs(dwt-zwt)*1000.
            if(delwt>dwt) then
              lwt=.true.
              dwt=zwt
            end if
          end if   
        end if        
      end do wt_rise_loop
! map unsaturated zone outflux to grid
! there are 2*nts+1 increments in qtime()
      do k=1,nts  
       qb(k)=qtime(2*k+1)
       if(outp(7)) rik1(i+(k-1)*imax)=qb(k)/ks(zo(i))
      end do
     else ! top of capillary fringe at ground surface, so use surface flux
      qb=0. ! initialize qb for case where ts>capt(nper+1)
      dwt=depth(i)
      delh=0.;rf=0.;zwt=depth(i)
      do j=1,nts+1
       do k=1,kper
        if(ts>=capt(k) .and. ts<=capt(k+1)) qb(j)=q(k)
       end do
       if(outp(7)) rik1(i+(j-1)*imax)=qb(j)/ks(zo(i))
       tcap(j)=ts ! pass to diffusion subroutine
       ts=ts+tinc
      end do
      do j=1,nts+1
       qbij(j)=qb(j)/ks(zo(i))
      end do
      rf=0.
      svgctr=svgctr+1
      !call svgstp(u1,qbij,LIS_logunit,i,rf,nccs,lcvs,nmxs) ! SY
      call svgstp(u1,qbij,i,rf,nccs,lcvs,nmxs) ! SY
      nmns=nmn
     end if
    end do grid_loop
    write (*,*) imx1, ' cells completed' 
    write (LIS_logunit,*) imx1, ' cells completed' 
!
    if(nmin>nmax2) nmin=nmax2; if(nmin1>nmax3) nmin1=nmax3
    write(LIS_logunit,*) 'Convergence data for unsaturated zone:'
    write(LIS_logunit,*) 'Maximum terms used by Fourier series', nmax2,nmax3
    write(LIS_logunit,*) 'Minimum terms used by Fourier series', nmin,nmin1
    write(LIS_logunit,*) 'Unsaturated zone nonconvergent cells: '
    write(LIS_logunit,*) ncc
    if(nmnp>nmxp) nmnp=nmxp; if(nmns>nmxs) nmns=nmxs
    write(LIS_logunit,*) 'Convergence data for saturated zone:'
    write(LIS_logunit,*) 'Max. terms used by error-function series', nmxp,nmxs
    write(LIS_logunit,*) 'Min. terms used by error-function series', nmnp,nmns
    write(LIS_logunit,*) 'Saturated-zone nonconvergent cells: '
    write(LIS_logunit,*) nccs
    write(LIS_logunit,*) 'Cells using svgstp() = ',svgctr
    return
  end subroutine unsfin


!   SY: Subroutine indentation performed for integration into LIS
!  diffusion model for pressure head from applied load of rising water table
!By R.L. Baum and W.Z. Savage, USGS, May 2004, latest revision 21 Jan 2011
  subroutine pstpf(u1,dhwt,dwt,&
       & i,iper,rf,nccs,lcvs,nmx,t0,jf) ! SY
       !& ulog,i,iper,rf,nccs,lcvs,nmx,t0,jf) ! SY
    use grids
    use input_vars
    use model_vars
    implicit none
    !integer:: i,j,jf,u1,nmx,ulog,n1,nccs,n,m,iper,jmark,nstp,jtop ! SY
    integer:: i,j,jf,u1,nmx,n1,nccs,n,m,iper,jmark,nstp,jtop ! SY
    real:: dwt !,chi 12/22/2010 made chi a 1-d array & moved to model_vars 
    real (double):: finf,t1,t2,term1,term2,rn,znew,t0
    real (double):: ferfc1,ferfc2,ferfc3,ferfc4,dhwt(nts+1),dh
    real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
    real (double):: ar1,ar2,ar3,ar4,fs,rf(nzs+1)
    real (double):: rfa,rfb,derfc,ff,rslo,rphi,fmn,ptest
    real (double):: tol,delt1,delt2,t1old,t2old
    real (double):: ddg2rad,dusz0,dusz1,dlz,d1,zm
    real (double):: tstar,flt1,flt2,uwt1,uwsum
    logical:: lcvs 
    pi=3.141592653589793
    ddg2rad=pi/180.D0
    !  maximum value of Factor of Safety
    finf=10.
    nmx=0
    nmn=1+mmax
    tol=1.e-06
    dh=dhwt(iper+1)
    rslo=slo(i)
    if (flag<=-1 .and. flag>=-3) then !12/14/10 added flag>=-3
      write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
      & 'cell ',i,rslo/ddg2rad,iper+1,t0 ! Added label 4/22/2010 RLB
    end if
    rphi=phi(zo(i))
    a1=sin(rslo)
    b1=cos(rslo)
    select case (flowdir) ! set value of beta (Iverson's beta line)
      case ('slope')
       beta=b1*b1
      case ('hydro')
       beta=1.d0
      case default
       beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
    end select
    if(abs(b1-rikzero(i))<1.e-6) beta=0.d0 
    if (abs(rslo)>1.e-5) then
      ff=tan(rphi)/tan(rslo)
    else !  set factor of safety to fixed value for flat slopes   
      ff=finf
    end if
    zns=float(nzs)
    zinc=(zmax(i)-zmin)/zns  
    nstp=int((1./alp(zo(i))/zinc))
    z=zmin
    fmn=1.e25
  ! compute depth of unsaturated zone, 
    if(unsat(zo(i))) then
      if(lps0) then
        dusz0=depth(i)-1./alp(zo(i)) ! use when psi0=-1/alp(zo(i))
      else
        dusz0=depth(i) ! use if psi0=0
      end if
      dusz1=dwt-1./alp(zo(i))
    else
      dusz0=0.  
      dusz1=0.  
    end if
    dlz=zmax(i)-dusz0
    zm=zmax(i)
    uwsum=0. ! for computing depth-averaged unit weight
    uwsp=uws(zo(i))
    rf=0.0 
    jmark=1 ! Added in case depth(i)==0., 28 Jan 2013, RLB 
    Z_loop: do j=1,nzs+1
     znew=z
     if(z<depth(i))then 
       jmark=j
     end if
     if(znew < 1.0e-30) znew =1.0e-30
     pzero(j)=beta*(z-depth(i))
     if(z <= depth(i))then 
       pzero(j)=0.
       if(lps0) then ! psi0=-1/alpha
         if(llus) then ! compute water table rise
           if(z > dusz0-dh .and. z < depth(i)-dh) then  
             rf(j)=(-1./alp(zo(i))+(z-dusz0+dh))
           else if(z >= depth(i)-dh) then  
             rf(j)=beta*(z-depth(i)+dh)
           else 
             rf(j)=ptran(j)
           end if
         else ! static water table
           if(z > dusz0 .and. z < depth(i)) then  
             rf(j)=(-1./alp(zo(i))+(z-dusz0))
           else if(z >= depth(i)) then  
             rf(j)=beta*(z-depth(i))
           else
             rf(j)=ptran(j) 
           end if
         end if
       else ! psi0=0
         if(z >= dusz0-dh .and. llus) then ! water table rise
           rf(j)=beta*(z+dh-dusz0)
         else ! static water table
           rf(j)=ptran(j)
         end if
       end if
     else if(z>=zm .and. dlz<0.001) then
       rf(j)=beta*(z+dh-dusz0)      
     else if(llus) then  ! compute pressure rise below wt only if llus=.true.
       temporal_loop: do m=1,iper
         tdif1=t0-tcap(m)
         if(tdif1 > 0.0) then
          d1=dif(zo(i))/(b1*b1)
          t1=sqrt(d1*tdif1)
          tstar=t1/(dlz*dlz)
          if (t1<1.0e-29) t1=1.0e-29
          term1=0.0
          if(tstar>5.) then ! solution for later time...
            series_a_lt: do n=1,mmax
              rn=float(n)
              ar1=(2.d0*rn-1.d0)*pi*(z/dlz-1.d0)/2.d0
              ar2=(2.d0*rn-1.d0)*(2.d0*rn-1.d0)*pi*pi*tstar/4.d0
              flt1=cos((rn-1.d0)*pi)*exp(-ar2)*cos(ar1)/(2.d0*rn-1.d0)
! test for convergence of series       
              t1old=term1
              term1=term1+flt1
              delt1=abs(term1-t1old)
              if (term1/=0.) then
               delt1=delt1/abs(term1)
              end if
              n1=n
              if(delt1<=tol) exit
              if(abs(term1)<=tis .and. n>=4) then
                term1=0.d0
                delt1=0.d0
                exit
              end if
            end do series_a_lt
            term1=1.d0-term1*4.d0/pi
          else ! solution for early time
            series_a_et: do n=1,mmax
              rn=float(n)
              ar1=((2.*rn-1.)*dlz-(zm-z))/(2.*t1)
              ar2=((2.*rn-1.)*dlz+(zm-z))/(2.*t1)
              ferfc1=cos((rn+1.0)*pi)*derfc(ar1)
              ferfc2=cos((rn+1.0)*pi)*derfc(ar2)
! test for convergence of series       
              t1old=term1
              term1=term1+ferfc1+ferfc2
              delt1=abs(term1-t1old)
              if (term1/=0.) then
               delt1=delt1/abs(term1)
              end if
              n1=n
              if(delt1<=tol) exit
            end do series_a_et
          end if
          if(lcvs .and. delt1>tol) then
            nccs=nccs+1
            nv(i)=1
            lcvs=.false.
            write(*,*) 'SZ noncvg a',i,m,n1,t0,term1,t1old
            write(*,*) 'delt1,tol',delt1,tol
          end if
          if(n1>nmx) nmx=n1
          if(n1<nmn) nmn=n1
          rfa=term1
         else
          rfa=0.0
         end if
         tdif2=t0-tcap(m+1)
         if(tdif2 > 0.0) then
          d1=dif(zo(i))/(b1*b1)
          t2=sqrt(d1*tdif2)
          tstar=t2/(dlz*dlz)
          if (t2<1.0e-29) t2=1.0e-29
          term2=0.0 
          if(tstar>5.) then ! solution for later time...
            series_b_lt: do n=1,mmax
              rn=float(n) ! Added 14Sep2009, RLB, previously omitted in error.
              ar1=(2.d0*rn-1.d0)*pi*(z/dlz-1.d0)/2.d0
              ar2=(2.d0*rn-1.d0)*(2.d0*rn-1.d0)*pi*pi*tstar/4.d0
              flt2=cos((rn-1.d0)*pi)*exp(-ar2)*cos(ar1)/(2.d0*rn-1.d0)
! test for convergence of series      
              t2old=term2
              term2=term2+flt2
              delt2=abs(term2-t2old)
              if (term2/=0.) then
                delt2=delt2/abs(term2)
              end if
              n1=n
              if(n>10) write (*,*) 'b_lt: n,term2,flt1,tstar',n,term2,flt2,tstar
              if(delt2<=tol) exit
              if(abs(term2)<=tis .and. n>=4) then
                term2=0.d0
                delt2=0.d0
                exit
              end if
            end do series_b_lt
            term2=1.d0-term2*4.d0/pi
          else ! solution for early time
            series_b_et: do n=1,mmax 
              rn=float(n)
              ar3=((2.*rn-1.)*dlz-(zm-z))/(2.*t2)
              ar4=((2.*rn-1.)*dlz+(zm-z))/(2.*t2)
              ferfc3=cos((rn+1.0)*pi)*derfc(ar3)
              ferfc4=cos((rn+1.0)*pi)*derfc(ar4)
! test for convergence of series       
              t2old=term2
              term2=term2+ferfc3+ferfc4
              delt2=abs(term2-t2old)
              if (term2/=0.) then
                delt2=delt2/abs(term2)
              end if
              n1=n
              if(delt2<=tol) exit
            end do series_b_et
          end if
          if(lcvs .and. delt2>tol) then
            nccs=nccs+1
            nv(i)=1
            lcvs=.false.
            write(*,*)  'SZ noncvg b',i,m,n1,t0,term2,t2old
            write(*,*) 'delt2,tol',delt2,tol
          end if
          if(n1>nmx) nmx=n1
          if(n1<nmn) nmn=n1
          rfb=term2
         else
          rfb=0.0
         end if
         rf(j)=rf(j)+dh*beta*(rfa-rfb)
       end do temporal_loop 
     end if
     bline(j)=z*beta
     ptran(j)=rf(j)
     p(j)=pzero(j)+ptran(j)
     ptest=p(j)-bline(j)
     if(ptest > 0.0) then
       p(j)=bline(j)
     end if
 ! estimate partially saturated unit weight      
     if(p(j)<-1.d0/alp(zo(i))) then
       uwt1=(gs(zo(i))*(1-ths(zo(i)))+thz(j))*uww
     else
       uwt1=uws(zo(i))
     end if
     uwsum=uwsum+uwt1
     uwsp(j)=uwsum/float(j)
     z=z+zinc
    end do Z_loop  
! find new height of rising water table in zones of upward seepage   
    if(rikzero(i)<0.0) then
     zinc=(zmax(i)-zmin)/zns
     z=zmin
     newdep=0.0
     do j=1,nzs+1
      if(p(j)<0.0) newdep=z
      z=z+zinc
     end do
! adjust presures
     z=zmin
     do j=1,nzs+1
      if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
      if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
      z=z+zinc    
     end do
    end if
! Smooth data near initial water table
    if(llus .and. zm>depth(i)) then
      jtop=jmark-nstp
      if(jtop<1) jtop=1
      do j=jmark,jmark-5,-1
        if(ptran(j)>ptran(j+1)) then 
          ptran(j)=ptran(j+1)
          p(j)=ptran(j)+pzero(j)
        end if
      end do
    end if  
! Compute factor of safety & save results    
    z=zmin
    Z_FS_loop: do j=1,nzs+1    
      chi(j)=1.0 
  ! Approximate suction stress computation, moved 4/22/2010 RLB
      if(p(j)<0.) then ! for suction 2/25/2009 (formerly for suctions exceeding the air-entry value, -1/alpha)
        chi(j)=(thz(j)-thr(zo(i)))/(ths(zo(i))-thr(zo(i)))
      else
        chi(j)=1.0 
      end if
      if (abs(a1)>1.e-5 .and. z>0.) then
        fw(j)=-(chi(j)*p(j)*uww*tan(rphi))/(uwsp(j)*z*a1*b1) ! same as saturated formula with factor Chi
        fc(j)=c(zo(i))/(uwsp(j)*z*a1*b1) ! uses depth averaged unit weight
      else
        fw(j)=0.d0
        fc(j)=0.d0
      end if
      fs=ff+fw(j)+fc(j)
  ! frictional strength cannot be less than zero 
      if ((ff+fw(j))<0.) fs=fc(j)
      if (fs>finf) fs=finf
      if (z<=1.e-02) fs=finf 
      if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
      if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),&
      & pzero(j),ptran(j),bline(j),fs
      if (flag==-3) write(u1,'(6(g12.5,1x):)') z,p(j),fs,chi(j) ! added 4/14/2010 RLB
      if (fs<fmn) then
       fmn=fs
       if (jf>0) then
         zfmin(i+(jf-1)*imax)=z
         pmin(i+(jf-1)*imax)=p(j)
       end if
      end if
      z=z+zinc
    end do Z_FS_loop
    if (jf>0) then ! Revised 12/22/2010 RLB
      fsmin(i+(jf-1)*imax)=fmn
      if(fmn==finf) then ! Added 30 Jan 2013, RLB 
        pmin(i+(jf-1)*imax)=p(nzs+1)
        zfmin(i+(jf-1)*imax)=zmax(i)
      end if 
      if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
        !if (flag>=-6) call svijz(i,jf,dh,newdep,ulog) ! SY
        if (flag>=-6) call svijz(i,jf,dh,newdep) ! SY
        !if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,dh,newdep,ulog)  ! Added 2/10/2012 ! SY
        if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,dh,newdep)  ! Added 2/10/2012 ! SY
      end if
    end if
    return
  end subroutine pstpf


!   SY: Subroutine indentation performed for integration into LIS
!  finite depth diffusion model for rain infiltration.
!  by W.Z. Savage, Spring 2001, with modifications by R.L. Baum
!  (both) USGS, Latest revision 20 Jan 2011 RLB
  !subroutine svgstp(u1,rikf,ulog,i,rf,nccs,lcvs,nmx) ! SY
  subroutine svgstp(u1,rikf,i,rf,nccs,lcvs,nmx) ! SY
    use grids
    use input_vars
    use model_vars
    implicit none
    !integer:: i,j,u1,nmx,ulog,n1,nccs,n,m,mm,jf ! SY
    integer:: i,j,u1,nmx,n1,nccs,n,m,mm,jf ! SY
    real:: rikf(nts+1) 
    real (double):: finf,t1,t2,term1,term2,rn,znew,t0
    real (double):: fierfc1,fierfc2,fierfc3,fierfc4
    real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
    real (double):: ar1,ar2,ar3,ar4,fs,rf(nzs+1)
    real (double):: rfa,rfb,derfc,ff,rslo,rphi,fmn,ptest,pmn ! pmn added 4/15/2010
    real (double):: tol,delt1,delt2,t1old,t2old
    real (double):: ddg2rad,dusz1,dlz
    logical:: lcvs 
    pi=3.141592653589793
    ddg2rad=pi/180.D0
    finf=10.
    nmx=0
    !nmn=1+mmax ! initialization moved to main program 07 Jan 2013 RLB 
    rslo=slo(i)
    rphi=phi(zo(i))
    a1=sin(rslo)
    b1=cos(rslo)
    p0zmx=0. ! Added 6 May 2013, RLB 
    select case (flowdir) ! set value of beta (Iverson's beta line)
      case ('slope')
        beta=b1*b1
      case ('hydro')
        beta=1.d0
      case default
        beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
    end select
    if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
    if (abs(rslo)>1.e-5) then
      ff=tan(rphi)/tan(rslo)
    else !  set factor of safety to fixed value for flat slopes   
      ff=finf
    end if
    zns=float(nzs)
    zinc=(zmax(i)-zmin)/zns
    dusz1=0.  
    dlz=zmax(i)-dusz1
    ptran=0. ! added 17Aug2009 RLB
    temporal_loop: do m=1,nts+1 
      t0=tcap(m) 
      fmn=1.e25
      jf=jsav(m) 
      if (flag<=-1 .and. flag>=-3) then !12/14/10 added flag>=-3
        write (u1,'(a5,i12,f6.1,i12,2x,g14.8)')&
        'cell ',i,rslo/ddg2rad,m,t0 ! added label 4/22/2010 RLB
      end if
      z=zmin
      Z_loop: do j=1,nzs+1
        znew=z
        if(znew < 1.0e-30) znew =1.0e-30
        if (abs(a1)>1.e-5) then
          fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
        else
          fc(j)=0.d0
        end if
        pzero(j)=beta*(z-depth(i))
        if(z<dusz1 .or. rikf(m)==0.)then
          rf(j)=0.0 
        else
          rf(j)=0.0 
          temporal_loop_1: do mm=1,nper
            tdif1=t0-capt(mm)
            if(tdif1 > 0.0) then
    !  corrected diffusivity term in next line (divide by b1*b1)       
              t1=sqrt(tdif1*dif(zo(i))/(b1*b1))
              if (t1<1.0e-29) t1=1.0e-29
              term1=0.0
              series_a: do n=1,mmax
                rn=float(n)
                ar1=((2.*rn-1.)*dlz-(zmax(i)-z))/(2.*t1)
                ar2=((2.*rn-1.)*dlz+(zmax(i)-z))/(2.*t1)
                fierfc1=exp(-ar1**2)/sqrt(pi)-ar1*derfc(ar1)
                fierfc2=exp(-ar2**2)/sqrt(pi)-ar2*derfc(ar2)
    ! test for convergence of series to within 1/10000 of previous value      
                t1old=term1
                tol=term1/10000.
                term1=term1+fierfc1+fierfc2
                delt1=abs(term1-t1old)
                n1=n
                if(delt1<=tol) exit
              end do series_a
              if(lcvs .and. delt1>tol) then
                nccs=nccs+1
                nv(i)=1
                lcvs=.false.
              end if
              if(n1>nmx) nmx=n1
              if(n1<nmn) nmn=n1
              rfa=2.*t1*term1
            else
              rfa=0.0
            end if
            tdif2=t0-capt(mm+1)
            if(tdif2 > 0.0) then
    !  corrected diffusivity term in next line (divide by b1*b1)       
              t2=sqrt(tdif2*dif(zo(i))/(b1*b1))
              if (t2<1.0e-29) t2=1.0e-29
              term2=0.0 
              series_b: do n=1,mmax
                rn=float(n)
                ar3=((2.*rn-1.)*dlz-(zmax(i)-z))/(2.*t2)
                ar4=((2.*rn-1.)*dlz+(zmax(i)-z))/(2.*t2)
                fierfc3=exp(-ar3**2)/sqrt(pi)-ar3*derfc(ar3)
                fierfc4=exp(-ar4**2)/sqrt(pi)-ar4*derfc(ar4)
    ! test for convergence of series to within 1/10000 of previous value      
                t2old=term2
                tol=term2/10000.
                term2=term2+fierfc3+fierfc4
                delt2=abs(term2-t2old)
                n1=n
                if(delt2<=tol) exit
              end do series_b
              if(lcvs .and. delt2>tol) then
                nccs=nccs+1
                nv(i)=1
                lcvs=.false.
              end if
              if(n1>nmx) nmx=n1
              if(n1<nmn) nmn=n1
              rfb=2.*t2*term2
            else
              rfb=0.0
            end if
            rf(j)=rf(j)+rik(i+(mm-1)*imax)*(rfa-rfb)
          end do temporal_loop_1 
        end if
        bline(j)=z*beta
        if(abs(rf(j))>0.0) then ! added 17AUG2009 RLB 
          ptran(j)=rf(j)
        end if
        p(j)=pzero(j)+ptran(j)
        ptest=p(j)-bline(j)
        if(ptest > 0.0) then
          p(j)=bline(j)
        end if
        z=z+zinc
      end do Z_loop  
      if(m==1) p0zmx=p(nzs+1) ! Added 6 May 2013, RLB 
    ! find new height of rising water table in zones of upward seepage   
      if(rikzero(i)<0.0) then
        zinc=(zmax(i)-zmin)/zns
        z=zmin
        newdep=0.0
        do j=1,nzs+1
          if(p(j)<0.0) newdep=z
          z=z+zinc
        end do
    ! adjust presures 
        z=zmin
        do j=1,nzs+1
          if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
          if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
          z=z+zinc    
        end do
      end if
    ! Compute factor of safety & save results    
      z=zmin
      Z_FS_loop: do j=1,nzs+1    
        if (abs(a1)>1.e-5 .and. z>1.e-30) then
          if(lpge0 .and. p(j)<0.) then !option added 4/15/2010
            fw(j)=0.d0
          else
            fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
          end if
        else
          fw(j)=0.d0
        end if
        fs=ff+fw(j)+fc(j)
    ! frictional strength cannot be less than zero
        if ((ff+fw(j))<0.) fs=fc(j)
        if (fs>finf) fs=finf
        if (z<=1.e-02) fs=finf 
        if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
        if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),&
        & pzero(j),ptran(j),bline(j),fs
        if (flag==-3 .and. unsat0) then ! revised 12/23/2010
          write(u1,'(6(g12.5,1x):)') z,p(j),fs,1. ! added 4/14/2010 RLB
        else if (flag==-3) then
          write(u1,'(6(g12.5,1x):)') z,p(j),fs
        end if
        if (jf>0) then
          if (fs<fmn) then
            fmn=fs 
            zfmin(i+(jf-1)*imax)=z
            pmn=p(j) ! revised 4/15/2010
          end if
        end if
        z=z+zinc
      end do Z_FS_loop
    !  next statement assumes that computations begin at surface and work downward   
      if (jf>0) then  ! revised 4/29/2010 to include pmin() RLB
        fsmin(i+(jf-1)*imax)=fmn
        if(fmn==finf) then ! Added 30 Jan 2013, RLB 
          pmn=p(nzs+1)
          zfmin(i+(jf-1)*imax)=zmax(i)
        end if 
        if(lpge0 .and. pmn<0.) then !option added 4/15/2010
          pmin(i+(jf-1)*imax)=0.
        else
          pmin(i+(jf-1)*imax)=pmn
        end if
        if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
            dcf=0. ! this line used only for saturated infiltration model
            chi=1.d0 !added 12/23/2010
          !if (flag>=-6) call svijz(i,jf,0.d0,newdep,ulog) ! SY
          if (flag>=-6) call svijz(i,jf,0.d0,newdep) ! SY
          !if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,0.d0,newdep,ulog)  ! Added 2/10/2012 ! SY
          if (flag<=-7 .and. flag>=-9) call svxmdv(i,jf,0.d0,newdep)  ! Added 2/10/2012 ! SY
        end if
      end if
    end do temporal_loop 
    return
  end subroutine svgstp


!   SY: Subroutine indentation performed for integration into LIS
!  Implementation of Iverson's (2000) method of computing pore pressure
!  for rain infiltration.
!  by W.Z. Savage, spring 2001, with modifications by R.L. Baum
  !subroutine iverson(imx1,u1,profil,ulog) ! SY
  subroutine iverson(imx1,u1,profil) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
        LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use grids
    use input_vars
    use model_vars
    implicit none
    !integer:: j,i,u1,ulog,imx1,n ! SY
    integer:: j,i,u1,imx1,n ! SY
    real:: finf 
    real (double) :: derfc,a1,b1,ff,zns,zinc,z,znew
    real (double) :: zstar,tstar,x1,x2,x3,x4
    real (double) :: rf1,rf2,rf3,rf4,rfa,rfb,rf
    real (double) :: fs,rslo,rphi,fmn,ptest,pmn,dhat ! pmn added 4/15/2010, dhat added 8 Jan 2013
    real (double) :: newdep,captstar1,captstar2,tdif1,tdif2
    character (len=255) profil
    write(LIS_logunit,*) 'Starting saturated-zone' ! SY
    write(LIS_logunit,*) 'computations for infinite-depth' ! SY
    write(*,*) 'Starting saturated-zone'
    write(*,*) 'computations for infinite-depth'
    pi=3.141592653589793
    dg2rad=pi/180.D0
    !  maximum value of Factor of Safety
    finf=10.
    !  loop steps through all grid cells
    write(*,*) 'Cells completed: '
    grid_loop: do i=1,imx1
     rslo=slo(i) 
     if(rslo<slomin) then
       fsmin(i)=finf+1.
       zfmin(i)=zmax(i)
       pmin(i)=0.
       if (mod(i,2000)==0) write (*,*) i
       cycle grid_loop
     end if
     if (flag<=-1 .and. flag>=-3) write (u1,'(a5,i12,f6.1,i2,2x,g14.8)')&
     & 'cell ',i,rslo/dg2rad,1,t ! added label, t 4/22/2010, 12/14/10 added flag>=-3
     rphi=phi(zo(i))
     a1=sin(rslo)
     b1=cos(rslo)
     dhat=4.*dif(zo(i))/(b1*b1) ! added 8 Jan 2013, RLB
     select case (flowdir) ! set value of beta (Iverson's beta line)
       case ('slope')
        beta=b1*b1
       case ('hydro')
        beta=1.d0
       case default
        beta=b1*b1-rikzero(i)
     end select
     if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
     if (abs(rslo)>1.e-5) then
       ff=tan(rphi)/tan(rslo)
     else
   !  set factor of safety to fixed value for flat slopes   
       ff=finf
     end if
     zns=float(nzs)
     zinc=(zmax(i)-zmin)/zns
     z=zmin
     fmn=1.e25
     rf=0.0
     z_loop: do j=1,nzs+1
       znew=z
       if(znew < 1.0e-30) znew =1.0e-30
       if (abs(a1)>1.e-5) then
         fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
       else
         fc(j)=0.d0
       end if
       pzero(j)=beta*(z-depth(i))
       if (abs(z)>0.) then ! Formulas in next 2 lines apply only for z>0, test added 21 Feb 2013, RLB
         zstar=z**2/dhat ! formula simplified 8 Jan 2013, RLB
         tstar=t/zstar
       end if
       rf=0.0
       temporal_loop: do n=1,nper
         if(z==0.) then ! exact formula added for case of z=0, 8 Jan 2013, RLB
           tdif1=t-capt(n)
           if(tdif1>0.) then
             rfa=sqrt(tdif1*dhat/pi)
           else
             rfa=0.0
           end if
           tdif2=t-capt(n+1)
           if(tdif2>0.) then
             rfb=sqrt(tdif1*dhat/pi)
           else
             rfb=0.0
           end if
         else ! z>0
           captstar1=capt(n)/zstar
           tdif1=tstar-captstar1
           if(tdif1 > 0.0) then 
             x1=1./tdif1
             x2=1./(sqrt(tdif1))
             rf1=sqrt(1./(x1*pi))*exp(-x1)
             rf2=derfc(x2)
             rfa=rf1-rf2
           else
             rfa=0.0
           end if
           captstar2=capt(n+1)/zstar
           tdif2=tstar-captstar2
           if(tdif2 > 0.0) then
             x3=1./tdif2
             x4=1./(sqrt(tdif2))
             rf3=sqrt(1./(x3*pi))*exp(-x3)
             rf4=derfc(x4)
             rfb=rf3-rf4
           else
             rfb=0.0
            end if
         end if
         rf=rf+rik(i+(n-1)*imax)*(rfa-rfb)
       end do temporal_loop
       ptran(j)=z*rf
       if(z==0) ptran(j)=rf ! formula for ptran at z=0 not normalized, 9 Jan 2013, RLB 
       p(j)=pzero(j)+ptran(j)
       bline(j)=z*beta
       ptest=p(j)-bline(j)
       if(ptest > 0.0) then
         p(j)=bline(j)
       end if
       if (abs(a1)>1.e-5) then
         if(lpge0 .and. p(j)<0.) then !option added 4/15/2010
           fw(j)=0.d0
         else if (z>0.) then ! Added z>0 condition 2/12/2013, RLB
           fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
         end if
       else
         fw(j)=0.d0
       end if
       z=z+zinc
     end do z_loop
   ! find new height of rising water table in zones of upward seepage   
     if(rikzero(i)<0.0) then
       zinc=(zmax(i)-zmin)/zns
       z=zmin
       newdep=0.0
       z_loop_a: do j=1,nzs+1
        if(p(j)<0.0) newdep=z
        z=z+zinc
       end do z_loop_a
! adjust presures 
       z=zmin
       z_loop_b: do j=1,nzs+1
         if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
         if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
         z=z+zinc    
       end do z_loop_b
     end if
     z=zmin
     fs_loop: do j=1,nzs+1
       fs=ff+fw(j)+fc(j)
   ! frictional strength cannot be less than zero 
       if ((ff+fw(j))<0.) fs=fc(j)
       if (fs>finf) fs=finf
       if (z<=1.e-02) fs=finf 
       if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
       if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),pzero(j),ptran(j),&
              & bline(j),fs
       if (flag==-3) write(u1,'(6(g12.5,1x):)') z,p(j),fs ! added 4/14/2010 RLB, Revised 12/23/2010
       if (fs<fmn) then
         fmn=fs
         zfmin(i)=z
         pmn=p(j) ! revised 4/15/2010
       end if
       z=z+zinc
     end do fs_loop
   !  next statement assumes that computations begin at surface and work downward   
     fsmin(i)=fmn
     if(fmn==finf) then ! Added 30 Jan 2013, RLB 
      pmn=p(nzs+1)
      zfmin(i)=zmax(i)
     end if 
     if(lpge0 .and. pmn<0.) then !option added 4/15/2010
       pmin(i)=0.
     else
       pmin(i)=pmn
     end if
     if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
       dcf=0. ! this line used only for saturated infiltration model
       chi=1.d0 !added 12/23/2010
       !if (flag>=-6) call svijz(i,1,0.d0,newdep,ulog) ! SY
       if (flag>=-6) call svijz(i,1,0.d0,newdep) ! SY
       !if (flag<=-7 .and. flag>=-9) call svxmdv(i,1,0.d0,newdep,ulog)  ! Added 2/10/2012 ! SY
       if (flag<=-7 .and. flag>=-9) call svxmdv(i,1,0.d0,newdep)  ! Added 2/10/2012 ! SY
     end if
     if (mod(i,2000)==0) write (*,*) i
    end do grid_loop
    write(*,*) imax, ' cells completed' 
    !write(ulog,*) imax, ' cells completed' ! SY
    write(LIS_logunit,*) imax, ' cells completed' ! SY
    return
  end subroutine iverson


!   SY: Subroutine indentation performed for integration into LIS
!  finite depth diffusion solution for rain infiltration.
!  by W.Z. Savage, Spring 2001, with modifications by R.L. Baum
!  (both) USGS
  !subroutine savage(imx1,u1,profil,ulog,nccs) ! SY
  subroutine savage(imx1,u1,profil,nccs) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
        LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use grids
    use input_vars
    use model_vars
    implicit none
    !integer :: j,i,imx1,u1,nmx,ulog,n1,nccs,n,m ! SY
    integer :: j,i,imx1,u1,nmx,n1,nccs,n,m ! SY
    real (double):: finf,t1,t2,term1,term2,rn,znew
    real (double):: fierfc1,fierfc2,fierfc3,fierfc4
    real (double):: a1,b1,zns,zinc,tdif1,tdif2,z,newdep
    real (double):: ar1,ar2,ar3,ar4,fs
    real (double):: rfa,rfb,rf,derfc,ff,rslo,rphi,fmn,ptest,pmn ! pmn added 4/15/2010
    real (double):: tol,delt1,delt2,t1old,t2old
    real (double):: ddg2rad,dlz 
    logical :: lcv 
    character (len=255)::  profil
    write(LIS_logunit,*) 'Starting saturated-zone' ! SY
    write(LIS_logunit,*) 'computations for finite-depth' ! SY
    write(*,*) 'Starting saturated-zone'
    write(*,*) 'computations for finite-depth'
    pi=3.141592653589793
    ddg2rad=pi/180.D0
    !  maximum value of Factor of Safety
    finf=10.
    !  loop steps through all grid cells
    write(*,*) 'Cells completed: '
    nmx=0
    !nmn=1+mmax ! initalization of nmn moved to main program 07 Jan 2013, RLB 
    grid_loop: do i=1,imx1
      rslo=slo(i)
      if(rslo<slomin) then
        fsmin(i)=finf+1.
        zfmin(i)=zmax(i)
        pmin(i)=0.
        if (mod(i,2000)==0) write (*,*) i
        cycle grid_loop
      end if
      if (flag<=-1 .and. flag>=-3) write (u1,'(a5,i12,f6.1,i2,2x,g14.8)')&
      & 'cell ',i,rslo/dg2rad,1,t ! added label, t 4/22/2010, 12/14/10 added flag>=-3
      rphi=phi(zo(i))
      a1=sin(rslo)
      b1=cos(rslo)
      select case (flowdir) ! set value of beta (Iverson's beta line)
        case ('slope')
         beta=b1*b1
        case ('hydro')
         beta=1.d0
        case default
         beta=b1*b1-rikzero(i)  ! 2/12/09 corrected formula for beta
      end select
      if(abs(b1-rikzero(i))<1.e-6) beta=0.d0
      if (abs(rslo)>1.e-5) then
        ff=tan(rphi)/tan(rslo)
      else
    !  set factor of safety to fixed value for flat slopes   
        ff=finf
      end if
      zns=float(nzs)
      zinc=(zmax(i)-zmin)/zns
      z=zmin
      lcv=.true.
    ! compute depth of saturated zone, dlz
      dlz=zmax(i) 
      Z_loop: do j=1,nzs+1
        znew=z
        if(znew < 1.0e-30) znew =1.0e-30
        if (abs(a1)>1.e-5) then
          fc(j)=c(zo(i))/(uws(zo(i))*znew*a1*b1)
        else
          fc(j)=0.d0
        end if
        pzero(j)=beta*(z-depth(i))
        rf=0.0 
        temporal_loop: do m=1,nper
          tdif1=t-capt(m)
    ! the next if block assumes that no change in transient pore pressure
    ! occurs if no flux at the top of the saturated zone.
          if(rik(i+(m-1)*imax)==0.) then
            cycle temporal_loop      
          end if
          if(tdif1 > 0.0) then
    ! corrected diffusivity term in next line (divide by b1*b1)       
            t1=sqrt(tdif1*dif(zo(i))/(b1*b1))
            if (t1<1.0e-29) t1=1.0e-29
            term1=0.0
            series_a: do n=1,mmax
              rn=float(n)
              ar1=((2.*rn-1.)*dlz-(zmax(i)-z))/(2.*t1)
              ar2=((2.*rn-1.)*dlz+(zmax(i)-z))/(2.*t1)
              fierfc1=exp(-ar1**2)/sqrt(pi)-ar1*derfc(ar1)
              fierfc2=exp(-ar2**2)/sqrt(pi)-ar2*derfc(ar2)
    ! test for convergence of series to within 1/10000 of previous value      
              t1old=term1
              tol=term1/10000.
              term1=term1+fierfc1+fierfc2
              delt1=abs(term1-t1old)
              n1=n
              if(delt1<=tol) exit
            end do series_a
            if(lcv .and. delt1>tol) then
              nccs=nccs+1
              nv(i)=1
              lcv=.false.
            end if
            if(n1>nmx) nmx=n1
            if(n1<nmn) nmn=n1
            rfa=2.*t1*term1
          else
            rfa=0.0
          end if
          tdif2=t-capt(m+1)
          if(tdif2 > 0.0) then
    !  corrected diffusivity term in next line (divide by b1*b1)       
            t2=sqrt(tdif2*dif(zo(i))/(b1*b1))
            if (t2<1.0e-29) t2=1.0e-29
            term2=0.0 
            series_b: do n=1,mmax
              rn=float(n)
              ar3=((2.*rn-1.)*dlz-(zmax(i)-z))/(2.*t2)
              ar4=((2.*rn-1.)*dlz+(zmax(i)-z))/(2.*t2)
              fierfc3=exp(-ar3**2)/sqrt(pi)-ar3*derfc(ar3)
              fierfc4=exp(-ar4**2)/sqrt(pi)-ar4*derfc(ar4)
    ! test for convergence of series to within 1/10000 of previous value      
              t2old=term2
              tol=term2/10000.
              term2=term2+fierfc3+fierfc4
              delt2=abs(term2-t2old)
              n1=n
              if(delt2<=tol) exit
            end do series_b
            if(lcv .and. delt2>tol) then
              nccs=nccs+1
              nv(i)=1
              lcv=.false.
            end if
            if(n1>nmx) nmx=n1
            if(n1<nmn) nmn=n1
            rfb=2.*t2*term2
          else
            rfb=0.0
          end if
          rf=rf+rik(i+(m-1)*imax)*(rfa-rfb)
        end do temporal_loop 
        ptran(j)=rf
        p(j)=pzero(j)+ptran(j)
        bline(j)=z*beta
        ptest=p(j)-bline(j)
        if(ptest > 0.0) then
          p(j)=bline(j)
        end if
        if (abs(a1)>1.e-5) then
          if(lpge0 .and. p(j)<0.) then !option added 4/15/2010
            fw(j)=0.d0
          else if (z>0.) then ! Added z>0 condition 2/12/2013, RLB
            fw(j)=-(p(j)*uww*tan(rphi))/(uws(zo(i))*z*a1*b1)
          end if
        else
          fw(j)=0.d0
        end if
        z=z+zinc
      end do Z_loop  
    ! find new height of rising water table in zones of upward seepage   
      if(rikzero(i)<0.0) then
        zinc=(zmax(i)-zmin)/zns
        z=zmin
        newdep=0.0
        do j=1,nzs+1
          if(p(j)<0.0) newdep=z
          z=z+zinc
        end do
    ! adjust presures 
        z=zmin
        do j=1,nzs+1
          if(p(j)>0.0 .and. z<newdep) p(j)=0.d0
          if(p(j)>=0.0 .and. z>=newdep) p(j)=beta*(z-newdep)
          z=z+zinc    
        end do
      end if
      z=zmin
      fmn=1.e25
      Z_FS_loop: do j=1,nzs+1
        fs=ff+fw(j)+fc(j)
    ! frictional strength cannot be less than zero
        if ((ff+fw(j))<0.) fs=fc(j)
        if (fs>finf) fs=finf
        if (z<=1.e-02) fs=finf 
        if (flag==-1) write(u1,'(6(g12.5,1x):)') z,p(j),fs
        if (flag==-2) write(u1,'(6(g12.5,1x):)') z,p(j),&
        & pzero(j),ptran(j),bline(j),fs
        if (flag==-3) write(u1,'(6(g12.5,1x):)') z,p(j),fs ! added 4/14/2010 RLB, Revised 12/23/2010
        if (fs<fmn) then
          fmn=fs
          zfmin(i)=z
          pmn=p(j) ! revised 4/15/2010
        end if
        z=z+zinc
      end do Z_FS_loop
    !  next statement assumes that computations begin at surface and work downward   
      fsmin(i)=fmn
      if(fmn==finf) then ! Added 30 Jan 2013, RLB 
        pmn=p(nzs+1)
        zfmin(i)=zmax(i)
      end if 
      if(lpge0 .and. pmn<0.) then !option added 4/15/2010
        pmin(i)=0.
      else
        pmin(i)=pmn
      end if
      if (flag<=-4 .or. outp(1)) then ! Added 12/22/2010 RLB, Rev. 8/5/2011, 2/10/2012
        dcf=0. ! this line used only for saturated infiltration model
        chi=1.d0 !added 12/23/2010
        !if (flag>=-6) call svijz(i,1,0.d0,newdep,ulog) ! SY
        if (flag>=-6) call svijz(i,1,0.d0,newdep) ! SY
        !if (flag<=-7 .and. flag>=-9) call svxmdv(i,1,0.d0,newdep,ulog)  ! Added 2/10/2012 ! SY
        if (flag<=-7 .and. flag>=-9) call svxmdv(i,1,0.d0,newdep)  ! Added 2/10/2012 ! SY
      end if
      if (mod(i,2000)==0) write (*,*) i
    end do grid_loop
    write(*,*) imx1, ' cells completed'
    write(LIS_logunit,*) imx1, ' cells completed' ! SY
    if(t==0 .and. nper==1) then ! added 4/14/2010
      nmx=0; nmn=0
    end if
    write(LIS_logunit,*) 'Max. terms used by error-function series', nmx ! SY
    write(LIS_logunit,*) 'Min. terms used by error-function series', nmn ! SY
    write(LIS_logunit,*) 'Saturated-zone nonconvergent cells: ' ! SY
    write(LIS_logunit,*) nccs ! SY
    return
  end subroutine savage


!   SY: Subroutine indentation performed for integration into LIS
  !subroutine satfin(imx1,ulog,u1,nccs) ! SY
  subroutine satfin(imx1,u1,nccs) ! SY
! 7/19/2006 Rex L. Baum, USGS 
!  computations for time series in the fully saturated case.	
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    use grids
    use input_vars 
    use model_vars
    use input_file_defs
    implicit none
    !integer::i,j,jf,k,ulog,u1,imx1,nccs ! SY
    integer::i,j,jf,k,u1,imx1,nccs ! SY
    integer::nmn1,nmin1,nmax3,nmxs,nmns,nmax0
    logical:: lcvs 
    real:: qbij(nts+1)
    real (double)::rf(nzs+1),finf
    nmax3=0;nmax0=0
    nmn1=nmax+1;nmin1=nmax+1
    write(LIS_logunit,*) 'Starting computations'
    write(LIS_logunit,*) 'for finite-depth saturated zone'
    write(*,*) 'Starting computations'
    write(*,*) 'for finite-depth saturated zone'
    write(*,*) 'Cells completed: '
! loop over all grid cells
    finf=10.
    do i=1,imx1 
     if (mod(i-1,2000)==0) write (*,*) i-1 ! cells completed
     if(slo(i)<slomin) then ! default values for gently sloping cells 
       do jf=1,nout
         fsmin(i+(jf-1)*imax)=finf+1.
         zfmin(i+(jf-1)*imax)=zmax(i)
         pmin(i+(jf-1)*imax)=0.
       end do
       cycle
     end if
     lcvs=.true.
     q=0.
     do j=1,kper
       if(j>nper) then
         q(j)=0.
       else
         q(j)=ks(zo(i))*rik(i+(j-1)*imax)
         !if(q(j)>ks(zo(i))) write (LIS_logunit,*) '*q>Ks!', i,j,q(j),ks(zo(i)) ! SY: See next line
         if(q(j)>ks(zo(i))) write (LIS_logunit,'(A,2I3,2ES14.6E2)') '*q>Ks!', i,j,q(j),ks(zo(i)) ! SY
       end if
     end do
!  use surface flux in infiltration computations
     qb=0. ! initialize qb for case where ts>capt(nper+1)
     ts=0.
     do j=1,nts+1
      do k=1,kper
       if(ts>=capt(k) .and. ts<=capt(k+1)) qb(j)=q(k)
      end do
      if(outp(7)) rik1(i+(j-1)*imax)=qb(j)/ks(zo(i))
      tcap(j)=ts ! pass to diffusion subroutine
      ts=ts+tinc
     end do
     do j=1,nts+1
      qbij(j)=qb(j)/ks(zo(i))
     end do
     rf=0. 
     call svgstp(u1,qbij,&
             & i,rf,nccs,lcvs,nmxs) ! SY
             !& LIS_logunit,i,rf,nccs,lcvs,nmxs) ! SY
     nmns=nmn
    end do
    write (*,*) imx1, ' cells completed' 
    write (LIS_logunit,*) imx1, ' cells completed' 
!
    if(nmns>nmxs) nmns=nmxs
    write(LIS_logunit,*) 'Convergence data for saturated zone:'
    write(LIS_logunit,*) 'Max. terms used by error-function series', nmxs
    write(LIS_logunit,*) 'Min. terms used by error-function series', nmns
    write(LIS_logunit,*) 'Saturated-zone nonconvergent cells: '
    write(LIS_logunit,*) nccs
    return
  end subroutine satfin


!   SY: Subroutine indentation performed for integration into LIS
  !subroutine satinf(imx1,ulog,u1,nccs) ! SY
  subroutine satinf(imx1,u1,nccs) ! SY
! 7/20/2006 Rex L. Baum, USGS, Latest revision 29 Jan 2013 
!  computations for time series in the fully saturated case.
!  calls ivestp() for infinite depth solution.	
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    use grids
    use input_vars 
    use model_vars
    use input_file_defs
    implicit none
    !integer::i,j,jf,k,ulog,u1,imx1,nccs ! SY
    integer::i,j,jf,k,u1,imx1,nccs ! SY
!integer::nmn1,nmin1,nmax3,nmax0,nmns,nmxs
!logical:: lcvs 
    real:: qbij(nts+1)
    real (double)::rf(nzs+1),finf
!nmax3=0;nmax0=0
!nmn1=nmax+1;nmin1=nmax+1
    write(LIS_logunit,*) 'Starting computations'
    write(LIS_logunit,*) 'for infinite-depth saturated zone'
    write(*,*) 'Starting computations'
    write(*,*) 'for infinite-depth saturated zone'
    write(*,*) 'Cells completed: '
! loop over all grid cells
    finf=10.
    do i=1,imx1 
     if (mod(i-1,2000)==0) write (*,*) i-1 ! cells completed
     if(slo(i)<slomin) then ! default values for gently sloping cells 
       do jf=1,nout
         fsmin(i+(jf-1)*imax)=finf+1.
         zfmin(i+(jf-1)*imax)=zmax(i)
         pmin(i+(jf-1)*imax)=0.
       end do
       cycle
     end if
!    lcvs=.true.
     q=0.
     do j=1,kper
       if(j>nper) then
         q(j)=0.
       else
         q(j)=ks(zo(i))*rik(i+(j-1)*imax)
         !if(q(j)>ks(zo(i))) write (LIS_logunit,*) '*q>Ks!', i,j,q(j),ks(zo(i)) ! SY: See next line
         if(q(j)>ks(zo(i))) write (LIS_logunit,'(A,2I3,2ES14.6E2)') '*q>Ks!', i,j,q(j),ks(zo(i)) ! SY
       end if
     end do
!  use surface flux in infiltration computations
     qb=0. ! initialize qb for case where ts>capt(nper+1)
     ts=0.
     do j=1,nts+1
      do k=1,kper
       if(ts>=capt(k) .and. ts<=capt(k+1)) qb(j)=q(k)
      end do
      if(outp(7)) rik1(i+(j-1)*imax)=qb(j)/ks(zo(i))
      tcap(j)=ts ! pass to diffusion subroutine
      ts=ts+tinc
     end do
     do j=1,nts+1
      qbij(j)=qb(j)/ks(zo(i))
     end do
     rf=0. 
     call ivestp(u1,qbij,&
            & i,rf) ! Revised 28 Jan 2013, RLB ! SY 
            !& LIS_logunit,i,rf) ! Revised 28 Jan 2013, RLB ! SY 
!                & LIS_logunit,i,rf,nccs,lcvs,nmxs)
!                nmns=nmn
    end do
    write (*,*) imx1, ' cells completed' 
    write (LIS_logunit,*) imx1, ' cells completed' 
    write(LIS_logunit,*) 'Saturated-zone nonconvergent cells: '
    write(LIS_logunit,*) nccs !initialized in trigrs main
    return
  end subroutine satinf


! SY: Subroutine indentation performed for integration into LIS
 !  subroutine to save a 1-d array as a 2-d grid, including nodata values
 !  By Rex L. Baum USGS, December 2000, latest revision 1/20/2011
   !subroutine isvgrd(z3,grd,z2,nrow,ncol,u,nodat,nodata, & ! SY
   !    mnd,param,u1,outfil,ti,header) ! SY
  subroutine isvgrd(z3,grd,z2,nrow,ncol,nodat,nodata, &
       mnd,param,outfil,ti,header) ! SY
! !USES:
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify

    implicit none
    !integer i,j,nrow,ncol,u,grd,ctr,u1,ia,m,mnd ! SY
    integer i,j,nrow,ncol,grd,ctr,ia,m,mnd ! SY
 !  Map 2-d array, z2, onto 1-d array, z3.
 !  Reverse rows and columns, because Fortran stores arrays
 !  by column.  z2 includes nodata values and is used to map
 !  z3 to the grid format
    real z2(ncol,nrow),nodats,test
    integer z3(grd)
    double precision nodat,nodata,param(6),ti,a
    character*14 header(6)
    character*1 sp
    character*255 outfil
    character*31 scratch
    integer ftn ! SY
    sp=char(32)
    ctr=0
    nodats=nodat
 !  open and initialize output grid files 
    outfil=adjustl(outfil)
    ftn = LIS_getNextUnitNumber() ! SY
    !open (u,file=trim(outfil), status='unknown',err=20) ! SY
    open (ftn,file=trim(outfil), status='unknown',err=20) ! SY
! 7/18/2006 added allocatable to nodata value in param array
    param(mnd)=nodata
    do 100, m=1,6
     a=param(m)
     ia=int(param(m))
 ! 19 Apr 2010 added special cases for m=1,2, and 6
     if(m.le.2 .or. m.eq.6) then
       write(scratch,'(i16)') ia
     else if(abs(a-ia).le.ti) then
       write(scratch,'(i16)') ia
     else
       write(scratch,*) a
     end if
     scratch=adjustl(scratch)
     write(ftn,1012) trim(header(m)),trim(scratch) ! SY
100 continue
 !  write grid 
    do 120, i=1,nrow
     do 110, j=1,ncol
      test=abs(z2(j,i)-nodats)
      if (test.le.0.1)then
       write(scratch,1005) int(z2(j,i))
       scratch=adjustl(scratch)
       if(j.ne.ncol) then
        write(ftn,1010) trim(scratch),sp ! SY
       else
        write(ftn,1011) trim(scratch) ! SY
       end if
      else
       ctr=ctr+1
       write(scratch,1005) z3(ctr)
       scratch=adjustl(scratch)
       if(j.ne.ncol) then
         write(ftn,1010) trim(scratch),sp ! SY 
       else
         write(ftn,1011) trim(scratch) ! SY 
       end if
      end if
110  continue
120 continue
    close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    return
20  continue
    write (*,*) 'Error opening output file'
    write (*,*) '--> ',trim(outfil)
    write (*,*) 'Check file name and status'
    write (*,*) 'Subroutine isvgrd()'
    write (LIS_logunit,*) 'Error opening output file'
    write (LIS_logunit,*) '--> ',trim(outfil)
    write (LIS_logunit,*) 'Check file name and status'
    write (LIS_logunit,*) 'Subroutine isvgrd()'
    write (LIS_logunit,*) '20 in isvgrd()' ! SY
    write (*,*) '20 in isvgrd()' ! SY
    !stop '20 in isvgrd()' ! SY
    call LIS_endrun() ! SY
1005 format(i12)
1010 format(a,a,$)
1011 format(a)
1012 format(t1,a,t15,a)
  end subroutine isvgrd ! SY


!  Routine to create file headers for various list file formats (1/20/2011)
!
!  Rex L. Baum, 21 Apr 2010, latest revision 07 Feb 2012
  subroutine prpijz(u1,profil,ncol,nrow,header,vrsn)
! !USES:
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify

    use grids
    use input_vars
    use model_vars
    use input_file_defs
!
!   SY: Subroutine indentation performed for integration into LIS
    implicit none
    !integer::i,j,jj,ctr,ncol,nrow,ulog,u1,ctr1,swt,m,nvar,outrow,zval ! SY
    integer::i,j,jj,ctr,ncol,nrow,u1,ctr1,swt,m,nvar,outrow,zval ! SY
    real (double):: test0
    real (double):: xmin,xmax,ymin,ymax,zlo,zhi,plo,pmax,chimin,chimax
    character (len=255)::outfil
    character (len=18):: profil
    character (len=14):: header(6)
    character (len=8)::fini(2)
    character (len=7):: vrsn ! added 12/28/2010
    character (len=4)::stp
    fini=(/'Infinite','Finite  '/)
    swt=1
    if(mmax>0) swt=2
    fini(swt)=adjustl(fini(swt))
! Open file and write header
! list file for storing depth profile on each cell
    if (flag<=-1 .and. flag>=-3) then
      outfil=trim(folder)//trim(profil)//trim(suffix)//'.txt'
      outfil=adjustl(outfil)
      open (u1,file=trim(outfil),status='unknown',err=10)
      write(u1,'(a,a)',err=11) 'Depth profiles at each cell, TRIGRS, v. ',vrsn
      write(u1,'(a,a)',err=11) trim(fini(swt)),'-depth no-flow boundary' ! format corrected 5/3/10 RLB
      write(u1,'(a)',err=11) 'Cell Number, Slope angle, Step#, Time'
      if (flag==-1) write(u1,'(a)',err=11) 'Z         P         FS'
      if (flag==-2) write(u1,'(a)',err=11) &
      &'Z         P      Pzero     Ptran      Pbeta       FS'
      if (flag==-3 .and. unsat0) then
        write(u1,'(a)',err=11) 'Z         P         FS      Chi' !added 4/14/2010 RLB
      else if (flag==-3) then
        write(u1,'(a)',err=11) 'Z         P         FS'
      end if
      return ! Reorganized if block and added return statement 9/6/2011, RLB
    end if
! ijz file for depth profile of each cell
    if(flag==-4 .or. flag==-5 .or. flag==-6) then ! revised 12/23/2010
      if(deepz>0) then
        zval=nzs+3
      else
        zval=nzs+2
      end if
      profil='TR_ijz_p_ch_'
      profil=adjustl(profil)
      do j=1,nout ! added loop to create series of files 15 Nov 2010 RLB
        !uijz(j)=j-1+uijz(1) !12/6/2010, added consecutive list of file unit numbers. ! SY
        uijz(j)=LIS_getNextUnitNumber() !12/6/2010, added consecutive list of file unit numbers. ! SY: Note redundancy for j = 1 since LIS already set unit number for uijz(1) in TRIGRS_final() 
        write(stp,'(i4)') j
        stp=adjustl(stp)
        outfil=trim(folder)//trim(profil)//trim(suffix)//'_'//trim(stp)//'.txt'
        outfil=adjustl(outfil)
        open(uijz(j),file=trim(outfil),status='unknown',err=10) ! Revised formating in succeeding write statements, RLB 9/2/2011
        write(uijz(j),'(a,a)',err=11) '# Comments: Quasi-3-D pressure head data from TRIGRS, v. ',vrsn
        write(uijz(j),'(a,a)',err=11) '# ', title ! revised format 1/6/2011 RLB
        if (unsat0) then
          write(uijz(j),'(a)',err=11) '# i  j  z  p  chi'
        else
          write(uijz(j),'(a)',err=11) '# i  j  z  p'
        end if
        write(uijz(j),'(a12,1x,i12)',err=11) '# timestep= ', ksav(j)
        write(uijz(j),'(a8,1x,g20.7)',err=11) '# time= ', tsav(j)
! New SCOOPS header format, updated 7/27/2011 RLB, Revised 2/28/2012 RLB
        write(uijz(j),'(a)',err=11) 'coords'
        write(uijz(j),'(a)',err=11) 'ijz'
      end do
      ctr=0;ctr1=0
    end if
    if(flag<=-4) then
!  Map i-j values to cell numbers for flag=-4, -5, -6
      do j=1,nrow
        do i=1,ncol
          ctr1=ctr1+1
          test0=abs(pf1(ctr1)-test1) ! check against default no-data value
          if (test0.le.0.1) then
            cycle
          else
            ctr=ctr+1
            jj=nrow-j+1 ! transform from origin at upper left to lower left corner of grid
            ix(ctr)=i; jy(ctr)=jj
          end if
        end do
      end do
    end if
! Copy coordinates of SW corner of grid for xyz file formats
    if (flag<=-7 .and. flag>=-9) then
      do m=1,6
        if (trim(header(m)).eq.'xllcorner') xllc=param(m)
        if (trim(header(m)).eq.'yllcorner') yllc=(param(m))
        if (trim(header(m)).eq.'west:') xllc=param(m)
        if (trim(header(m)).eq.'south:') yllc=param(m)
      end do
    end if
!
    if(flag<=-7 .and. flag>=-9) then ! Added 7/27/2011, RLB, Xmdv format, readable by VisIt 
      profil='TR_xyz_p_chi_'
      profil=adjustl(profil)
      nvar=5
! xmin, ymin at lower left corner of lower left cell, xmax, ymax at upper right corner of upper right grid cell
      xmin=xllc
      xmax=xllc+dfloat(ncol)*celsiz
      ymin=yllc
      ymax=yllc+dfloat(nrow)*celsiz
      pmax=maxval(zmax)
      plo=-pmax
      if(deepz > pmax) pmax=deepz
      zlo=minval(elev)-pmax
      zhi=maxval(elev)
      chimin=0
      chimax=1
      do j=1,nout
        !uijz(j)=j-1+uijz(1) ! SY
        uijz(j)=LIS_getNextUnitNumber() ! SY: Note redundancy for j = 1 since LIS already set unit number for uijz(1) in TRIGRS_final() 
        write(stp,'(i4)') j
        stp=adjustl(stp)
        outfil=trim(folder)//trim(profil)//trim(suffix)//'_'//trim(stp)//'.okc'
        outfil=adjustl(outfil)
        open(uijz(j),file=trim(outfil),status='unknown',err=10)
        if(deepz>0) then
          outrow=imax*(nzs+3)
        else
          outrow=imax*(nzs+2)
        end if
        write(uijz(j),*,err=11) nvar, outrow, 12 ! last entry on this line not used in all XMDV implemetations, but needed for VisIt to read database correctly RLB 2/21/2012
        write(uijz(j),*,err=11) 'x'
        write(uijz(j),*,err=11) 'y'
        write(uijz(j),*,err=11) 'z'
        write(uijz(j),*,err=11) 'p'
        write(uijz(j),*,err=11) 'chi'
        write(uijz(j),fmt='(2(g18.9,1x),i3)',err=11) xmin,xmax,10
        write(uijz(j),fmt='(2(g18.9,1x),i3)',err=11) ymin,ymax,10
        write(uijz(j),fmt='(2(g15.5,1x),i3)',err=11) zlo,zhi,10
        write(uijz(j),fmt='(2(g11.4,1x),i3)',err=11) plo,pmax,10
        write(uijz(j),fmt='(2(g11.4,1x),i3)',err=11) chimin,chimax,10
      end do
    end if
!
    return
10  continue
    write(*,*) 'Error opening output file ',outfil
    write(LIS_logunit,*) 'Error opening output file ',outfil
    close (LIS_logunit)
    stop '10 in subroutine prpijz()'
11  continue
    write(*,*) 'Error writing output file ',outfil
    write(LIS_logunit,*) 'Error writing output file ',outfil
    close (LIS_logunit)
    stop '11 in subroutine prpijz()'
  end subroutine prpijz


! SY: Subroutine indentation performed for integration into LIS
!  subroutine to read an ascii grid file (elevations) and determine its size (rows, columns, & data cells)
!  by Rex L. Baum, USGS Feb 2011 latest revison 13 Mar 2013
!  single precision
!
        !subroutine ssizgrd(row,col,celsiz,nodat,ctr,u,infil,header,ulog) ! SY
  subroutine ssizgrd(row,col,celsiz,nodat,ctr,infil,header) ! SY
! !USES:
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
         LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify
    implicit none
    integer,parameter:: double=kind(1d0)
    integer ::i,m,ctr!,u,ulog ! SY
    integer :: col,row,ctall
    real (double):: nodat,celsiz,cns,cew,param(6)
    real (double):: east,west,north,south
    real:: nodats
    real, allocatable:: temp(:)
    character (len=14):: header(6)
    character (len=255):: infil
    integer:: ftn ! SY
    !  
    infil=adjustl(infil) ! adjustl statements added to improve compatibility with other compilers 14 Feb 2013 RLB
    ftn = LIS_getNextUnitNumber() ! SY
    !open(u,file=trim(infil),status='old',err=23) ! SY
    open(ftn,file=trim(infil),status='old',err=23) ! SY
    do m=1,6
      read(ftn,*) header(m),param(m)
      header(m)=adjustl(header(m))
    end do
! set default value of nodat & celsiz for use with GRASS GIS ascii files   
    nodat=-9999.d0
    celsiz=-10.d0
    do m=1,6
     if (trim(header(m)).eq.'ncols') col=int(param(m))
     if (trim(header(m)).eq.'nrows') row=int(param(m))
     if (trim(header(m)).eq.'cellsize') celsiz=param(m)
     if (trim(header(m)).eq.'NODATA_value') nodat=param(m)
     if (trim(header(m)).eq.'nodata_value') nodat=param(m)
     if (trim(header(m)).eq.'cols:') col=int(param(m))
     if (trim(header(m)).eq.'rows:') row=int(param(m))
     if (trim(header(m)).eq.'east:') east=param(m)
     if (trim(header(m)).eq.'west:') west=param(m)
     if (trim(header(m)).eq.'north:') north=param(m)
     if (trim(header(m)).eq.'south:') south=param(m)
    end do
    if (celsiz.le.0) then
      cew=abs(east-west)/col
      cns=abs(north-south)/row
      if (cew.eq.cns) then
       celsiz=cew
      else
       celsiz=sqrt(cew*cns)
       write(*,*) 'Rectangular cells ',cew, ' X ', cns
       write(LIS_logunit,*) 'Rectangular cells ',cew, ' X ', cns
      end if
    end if
    nodats=nodat
    allocate(temp(col))
    ctr=0; ctall=0
    row_loop: do m=1,row
!  next sequence of lines read data in but skips no_data values
!  count maintained by ctr should coincide with node numbers from GIS
     read(ftn,*,end=125) (temp(i), i=1,col) 
     col_loop: do i=1,col
      ctall=ctall+1
      if(temp(i) /= nodats) then
       ctr=ctr+1
      end if
     end do col_loop
    end do row_loop
125 close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    write(*,*) ctr,' = number of data cells'
    write(*,*) ctall,' = total number of cells'
    return
23  continue
    write (*,*) '*** Error opening input file ***'
    write (*,*) '--> ',trim(infil)
    write (*,*) 'Check file name and location'
    write (LIS_logunit,*) '*** Error opening input file ***'
    write (LIS_logunit,*) '--> ',trim(infil)
    write (LIS_logunit,*) 'Check file name and location'
    !pause 'Press RETURN to exit' ! SY
    close(ftn) ! SY
    call LIS_releaseUnitNumber(ftn) ! SY
    !close(LIS_logunit) ! SY: DO NOT close log file
    write (LIS_logunit,*) '- Error in ssizgrd()' ! SY
    write (*,*) '- Error in ssizgrd()' ! SY
    !stop '- Error in ssizgrd()' ! SY
    call LIS_endrun() ! SY
  end subroutine ssizgrd


!   SY: Subroutine indentation performed for integration into LIS
!CS    REAL FUNCTION ERFC(X)
  DOUBLE PRECISION FUNCTION DERFC(X)
!C--------------------------------------------------------------------
!C
!C This subprogram computes approximate values for erfc(x).
!C   (see comments heading CALERF).
!C
!C   Author/date: W. J. Cody, January 8, 1985
!C
!C--------------------------------------------------------------------
    INTEGER JINT
!CS    REAL             X, RESULT
    DOUBLE PRECISION X, RESULT
!C------------------------------------------------------------------
    JINT = 1
    CALL CALERF(X,RESULT,JINT)
!CS    ERFC = RESULT
    DERFC = RESULT
    RETURN
!C---------- Last card of DERFC ----------
  END FUNCTION DERFC ! SY


!   SY: Subroutine indentation performed for integration into LIS
  SUBROUTINE CALERF(ARG,RESULT,JINT)
!C------------------------------------------------------------------
!C
!C This packet evaluates  erf(x),  erfc(x),  and  exp(x*x)*erfc(x)
!C   for a real argument  x.  It contains three FUNCTION type
!C   subprograms: ERF, ERFC, and ERFCX (or DERF, DERFC, and DERFCX),
!C   and one SUBROUTINE type subprogram, CALERF.  The calling
!C   statements for the primary entries are:
!C
!C                   Y=ERF(X)     (or   Y=DERF(X)),
!C
!C                   Y=ERFC(X)    (or   Y=DERFC(X)),
!C   and
!C                   Y=ERFCX(X)   (or   Y=DERFCX(X)).
!C
!C   The routine  CALERF  is intended for internal packet use only,
!C   all computations within the packet being concentrated in this
!C   routine.  The function subprograms invoke  CALERF  with the
!C   statement
!C
!C          CALL CALERF(ARG,RESULT,JINT)
!C
!C   where the parameter usage is as follows
!C
!C      Function                     Parameters for CALERF
!C       call              ARG                  Result          JINT
!C
!C     ERF(ARG)      ANY REAL ARGUMENT         ERF(ARG)          0
!C     ERFC(ARG)     ABS(ARG) .LT. XBIG        ERFC(ARG)         1
!C     ERFCX(ARG)    XNEG .LT. ARG .LT. XMAX   ERFCX(ARG)        2
!C
!C   The main computation evaluates near-minimax approximations
!C   from "Rational Chebyshev approximations for the error function"
!C   by W. J. Cody, Math. Comp., 1969, PP. 631-638.  This
!C   transportable program uses rational functions that theoretically
!C   approximate  erf(x)  and  erfc(x)  to at least 18 significant
!C   decimal digits.  The accuracy achieved depends on the arithmetic
!C   system, the compiler, the intrinsic functions, and proper
!C   selection of the machine-dependent constants.
!C
!C*******************************************************************
!C*******************************************************************
!C
!C Explanation of machine-dependent constants
!C
!C   XMIN   = the smallest positive floating-point number.
!C   XINF   = the largest positive finite floating-point number.
!C   XNEG   = the largest negative argument acceptable to ERFCX;
!C            the negative of the solution to the equation
!C            2*exp(x*x) = XINF.
!C   XSMALL = argument below which erf(x) may be represented by
!C            2*x/sqrt(pi)  and above which  x*x  will not underflow.
!C            A conservative value is the largest machine number X
!C            such that   1.0 + X = 1.0   to machine precision.
!C   XBIG   = largest argument acceptable to ERFC;  solution to
!C            the equation:  W(x) * (1-0.5/x**2) = XMIN,  where
!C            W(x) = exp(-x*x)/[x*sqrt(pi)].
!C   XHUGE  = argument above which  1.0 - 1/(2*x*x) = 1.0  to
!C            machine precision.  A conservative value is
!C            1/[2*sqrt(XSMALL)]
!C   XMAX   = largest acceptable argument to ERFCX; the minimum
!C            of XINF and 1/[sqrt(pi)*XMIN].
!C
!C   Approximate values for some important machines are:
!C
!C                          XMIN       XINF        XNEG     XSMALL
!C
!C  CDC 7600      (S.P.)  3.13E-294   1.26E+322   -27.220  7.11E-15
!C  CRAY-1        (S.P.)  4.58E-2467  5.45E+2465  -75.345  7.11E-15
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (S.P.)  1.18E-38    3.40E+38     -9.382  5.96E-8
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (D.P.)  2.23D-308   1.79D+308   -26.628  1.11D-16
!C  IBM 195       (D.P.)  5.40D-79    7.23E+75    -13.190  1.39D-17
!C  UNIVAC 1108   (D.P.)  2.78D-309   8.98D+307   -26.615  1.73D-18
!C  VAX D-Format  (D.P.)  2.94D-39    1.70D+38     -9.345  1.39D-17
!C  VAX G-Format  (D.P.)  5.56D-309   8.98D+307   -26.615  1.11D-16
!C
!C
!C                          XBIG       XHUGE       XMAX
!C
!C  CDC 7600      (S.P.)  25.922      8.39E+6     1.80X+293
!C  CRAY-1        (S.P.)  75.326      8.39E+6     5.45E+2465
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (S.P.)   9.194      2.90E+3     4.79E+37
!C  IEEE (IBM/XT,
!C    SUN, etc.)  (D.P.)  26.543      6.71D+7     2.53D+307
!C  IBM 195       (D.P.)  13.306      1.90D+8     7.23E+75
!C  UNIVAC 1108   (D.P.)  26.582      5.37D+8     8.98D+307
!C  VAX D-Format  (D.P.)   9.269      1.90D+8     1.70D+38
!C  VAX G-Format  (D.P.)  26.569      6.71D+7     8.98D+307
!C
!C*******************************************************************
!C*******************************************************************
!C
!C Error returns
!C
!C  The program returns  ERFC = 0      for  ARG .GE. XBIG;
!C
!C                       ERFCX = XINF  for  ARG .LT. XNEG;
!C      and
!C                       ERFCX = 0     for  ARG .GE. XMAX.
!C
!C
!C Intrinsic functions required are:
!C
!C     ABS, AINT, EXP
!C
!C
!C  Author: W. J. Cody
!C          Mathematics and Computer Science Division
!C          Argonne National Laboratory
!C          Argonne, IL 60439
!C
!C  Latest modification: March 19, 1990
!C
!C------------------------------------------------------------------
    INTEGER I,JINT
!CS    REAL
    DOUBLE PRECISION &
        A,ARG,B,C,D,DEL,FOUR,HALF,P,ONE,Q,RESULT,SIXTEN,SQRPI, &
        TWO,THRESH,X,XBIG,XDEN,XHUGE,XINF,XMAX,XNEG,XNUM,XSMALL, &
        Y,YSQ,ZERO
    DIMENSION A(5),B(4),C(9),D(8),P(6),Q(5)
!C------------------------------------------------------------------
!C  Mathematical constants
!C------------------------------------------------------------------
!CS    DATA FOUR,ONE,HALF,TWO,ZERO/4.0E0,1.0E0,0.5E0,2.0E0,0.0E0/,
!CS   1     SQRPI/5.6418958354775628695E-1/,THRESH/0.46875E0/,
!CS   2     SIXTEN/16.0E0/
    DATA FOUR,ONE,HALF,TWO,ZERO/4.0D0,1.0D0,0.5D0,2.0D0,0.0D0/, &
        SQRPI/5.6418958354775628695D-1/,THRESH/0.46875D0/, &
        SIXTEN/16.0D0/
!C------------------------------------------------------------------
!C  Machine-dependent constants
!C------------------------------------------------------------------
!CS    DATA XINF,XNEG,XSMALL/3.40E+38,-9.382E0,5.96E-8/,
!CS   1     XBIG,XHUGE,XMAX/9.194E0,2.90E3,4.79E37/
    DATA XINF,XNEG,XSMALL/1.79D308,-26.628D0,1.11D-16/, &
         XBIG,XHUGE,XMAX/26.543D0,6.71D7,2.53D307/
!C------------------------------------------------------------------
!C  Coefficients for approximation to  erf  in first interval
!C------------------------------------------------------------------
!CS    DATA A/3.16112374387056560E00,1.13864154151050156E02,
!CS   1       3.77485237685302021E02,3.20937758913846947E03,
!CS   2       1.85777706184603153E-1/
!CS    DATA B/2.36012909523441209E01,2.44024637934444173E02,
!CS   1       1.28261652607737228E03,2.84423683343917062E03/
    DATA A/3.16112374387056560D00,1.13864154151050156D02, &
           3.77485237685302021D02,3.20937758913846947D03, &
           1.85777706184603153D-1/
    DATA B/2.36012909523441209D01,2.44024637934444173D02, &
           1.28261652607737228D03,2.84423683343917062D03/
!C------------------------------------------------------------------
!C  Coefficients for approximation to  erfc  in second interval
!C------------------------------------------------------------------
!CS    DATA C/5.64188496988670089E-1,8.88314979438837594E0,
!CS   1       6.61191906371416295E01,2.98635138197400131E02,
!CS   2       8.81952221241769090E02,1.71204761263407058E03,
!CS   3       2.05107837782607147E03,1.23033935479799725E03,
!CS   4       2.15311535474403846E-8/
!CS    DATA D/1.57449261107098347E01,1.17693950891312499E02,
!CS   1       5.37181101862009858E02,1.62138957456669019E03,
!CS   2       3.29079923573345963E03,4.36261909014324716E03,
!CS   3       3.43936767414372164E03,1.23033935480374942E03/
    DATA C/5.64188496988670089D-1,8.88314979438837594D0, &
           6.61191906371416295D01,2.98635138197400131D02, &
           8.81952221241769090D02,1.71204761263407058D03, &
           2.05107837782607147D03,1.23033935479799725D03, &
           2.15311535474403846D-8/
    DATA D/1.57449261107098347D01,1.17693950891312499D02, &
           5.37181101862009858D02,1.62138957456669019D03, &
           3.29079923573345963D03,4.36261909014324716D03, &
           3.43936767414372164D03,1.23033935480374942D03/
!C------------------------------------------------------------------
!C  Coefficients for approximation to  erfc  in third interval
!C------------------------------------------------------------------
!CS    DATA P/3.05326634961232344E-1,3.60344899949804439E-1,
!CS   1       1.25781726111229246E-1,1.60837851487422766E-2,
!CS   2       6.58749161529837803E-4,1.63153871373020978E-2/
!CS    DATA Q/2.56852019228982242E00,1.87295284992346047E00,
!CS   1       5.27905102951428412E-1,6.05183413124413191E-2,
!CS   2       2.33520497626869185E-3/
    DATA P/3.05326634961232344D-1,3.60344899949804439D-1, &
           1.25781726111229246D-1,1.60837851487422766D-2, &
           6.58749161529837803D-4,1.63153871373020978D-2/
    DATA Q/2.56852019228982242D00,1.87295284992346047D00, &
           5.27905102951428412D-1,6.05183413124413191D-2, &
           2.33520497626869185D-3/
!C------------------------------------------------------------------
    X = ARG
    Y = ABS(X)
    IF (Y .LE. THRESH) THEN
!C------------------------------------------------------------------
!C  Evaluate  erf  for  |X| <= 0.46875
!C------------------------------------------------------------------
       YSQ = ZERO
       IF (Y .GT. XSMALL) YSQ = Y * Y
       XNUM = A(5)*YSQ
       XDEN = YSQ
       DO 20 I = 1, 3
          XNUM = (XNUM + A(I)) * YSQ
          XDEN = (XDEN + B(I)) * YSQ
20     CONTINUE
       RESULT = X * (XNUM + A(4)) / (XDEN + B(4))
       IF (JINT .NE. 0) RESULT = ONE - RESULT
       IF (JINT .EQ. 2) RESULT = EXP(YSQ) * RESULT
       GO TO 800
!C------------------------------------------------------------------
!C  Evaluate  erfc  for 0.46875 <= |X| <= 4.0
!C------------------------------------------------------------------
    ELSE IF (Y .LE. FOUR) THEN
       XNUM = C(9)*Y
       XDEN = Y
       DO 120 I = 1, 7
          XNUM = (XNUM + C(I)) * Y
          XDEN = (XDEN + D(I)) * Y
120    CONTINUE
       RESULT = (XNUM + C(8)) / (XDEN + D(8))
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
!C------------------------------------------------------------------
!C  Evaluate  erfc  for |X| > 4.0
!C------------------------------------------------------------------
    ELSE
       RESULT = ZERO
       IF (Y .GE. XBIG) THEN
          IF ((JINT .NE. 2) .OR. (Y .GE. XMAX)) GO TO 300
          IF (Y .GE. XHUGE) THEN
             RESULT = SQRPI / Y
             GO TO 300
          END IF
       END IF
       YSQ = ONE / (Y * Y)
       XNUM = P(6)*YSQ
       XDEN = YSQ
       DO 240 I = 1, 4
          XNUM = (XNUM + P(I)) * YSQ
          XDEN = (XDEN + Q(I)) * YSQ
240    CONTINUE
       RESULT = YSQ *(XNUM + P(5)) / (XDEN + Q(5))
       RESULT = (SQRPI -  RESULT) / Y
       IF (JINT .NE. 2) THEN
          YSQ = AINT(Y*SIXTEN)/SIXTEN
          DEL = (Y-YSQ)*(Y+YSQ)
          RESULT = EXP(-YSQ*YSQ) * EXP(-DEL) * RESULT
       END IF
    END IF
!C------------------------------------------------------------------
!C  Fix up for negative argument, erf, etc.
!C------------------------------------------------------------------
300 IF (JINT .EQ. 0) THEN
       RESULT = (HALF - RESULT) + HALF
       IF (X .LT. ZERO) RESULT = -RESULT
    ELSE IF (JINT .EQ. 1) THEN
       IF (X .LT. ZERO) RESULT = TWO - RESULT
    ELSE
       IF (X .LT. ZERO) THEN
          IF (X .LT. XNEG) THEN
                RESULT = XINF
          ELSE
                YSQ = AINT(X*SIXTEN)/SIXTEN
                DEL = (X-YSQ)*(X+YSQ)
                Y = EXP(YSQ*YSQ) * EXP(DEL)
                RESULT = (Y+Y) - RESULT
          END IF
       END IF
    END IF
800 RETURN
!C---------- Last card of CALERF ----------
  END SUBROUTINE CALERF ! SY


!   SY: Subroutine indentation performed for integration into LIS
!  Routine to prepare pressure head & relative saturation (chi) data for export to ijz format
!  Calling routine should set value of dcf if it is not initialized in main.
!  Rex L. Baum, 25 Feb 2010, latest revision 27 Apr 2012
  !subroutine svijz(i,jf,delh,newdep,ulog) ! SY
  subroutine svijz(i,jf,delh,newdep) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
        LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use grids
    use input_vars 
    use model_vars
    use input_file_defs
    implicit none
    !integer::i,jf,ulog,m,m1,nw,mfrst,mlast!,mnp(1),mctr,mmid ! SY
    integer::i,jf,m,m1,nw,mfrst,mlast!,mnp(1),mctr,mmid ! SY
    integer::mrpx(1),mrpn(1),dinc,wtctr,wctr,wptr(nzs+1),wtptr(nzs+1)
    real (double)::delh,newdep,zns,zinc,fdinc !,pab(nzs+1),pmn
    real (double):: dwat,z,ptop,pbot,zwat(nzs+1),zbot,ztop,chtop,chbot,zwat0
    real (double):: zdeep,pdeep,chdeep,zrpn,zrpx,relhgt!,x,y,zmid,zmn(1),zmx(1)
    real (double):: tol1,tol2,mp,plin(nzs+1),rp(nzs+1),rpmin,rpmax
    logical ::ldeep, lopmx,lopmn,loxlo
    !  compute water-table elevations from surface elevation and depths
    zbot=elev(i)-zmax(i)
    ztop=elev(i)-zmin
    zwat0=elev(i)-depth(i) ! Initial water table depth
    zwat=zwat0 ! initialize array values to initial water table depth
    ! Check for deepz <= 0 to skip deep node
    ldeep=.true.
    zdeep=elev(i)-deepz 
    if(zdeep>ztop) ldeep=.false.
    if(zdeep>zbot .and. ldeep) then
      zdeep=zbot-zmax(i) 
    end if
    zns=float(nzs) 
    zinc=(zmax(i)-zmin)/zns
    ptop=p(1)
    pbot=p(nzs+1)
    chtop=chi(1)
    chbot=chi(nzs+1)
    chdeep=1.
    ! estimate height of water table from results of unsaturated infiltration model
    if(dcf>0. .and. (unsat(zo(i)) .or. igcap(zo(i)))) then
      if(rikzero(i)>=0.0) then
        zwat(1)=elev(i)-(dusz-delh) ! statement in unsth() accounts for initial wt below zmax(i)
      else
        zwat(1)=elev(i)-newdep
      end if
      wtctr=1;wtptr(wtctr)=1;wptr(wtctr)=1  ! changed wtptr(wctr) to wtptr(wtctr) 18 April 2013, RLB
    else if (dcf<=0.) then
    ! find water table(s) from results of saturated infiltration model
      call dzero_brac(nzs,p,wtctr,wtptr) ! wtptr(j)=1 for intervals containing a water table
      chi=1.d0
      wctr=0;dwat=0.d0;wptr=0
      do nw=1,nzs+1
        if(wtptr(nw)>0) then
         if(wtptr(nw)==2) then ! water table at node
           wctr=wctr+1
           z=zmin+zinc*float(nw-1)
           dwat=z 
         end if
         if(wtptr(nw)==1) then ! water table between nodes
           wctr=wctr+1
           z=zmin+zinc*float(nw-1)
           dwat=z+zinc*(0.d0-p(nw))/(p(nw+1)-p(nw))
         end if
         if(rikzero(i)>=0.0) then
           zwat(nw)=elev(i)-dwat 
           wptr(wctr)=nw
         else
           zwat(nw)=elev(i)-newdep
           wptr(wctr)=nw
         end if
        end if
        if(wctr==wtctr) exit  ! all zero points accounted for
      end do
      if(p(nzs+1)<0.d0) then ! water table below bottom node (below zmax(i))
        wtptr(nzs+1)=1
        wtctr=wtctr+1;wctr=wctr+1
        zwat(nzs+1)=elev(i)-depth(i)+beta*(p(nzs+1)-p0zmx)  ! Approximate formula RLB 5/3/2013.
        if(zwat(nzs+1)<(elev(i)-depth(i))) zwat(nzs+1)=elev(i)-depth(i) ! prevents water table from falling below initial value
        wptr(wtctr)=nzs+1
      end if
    end if
    if(outp(1)) then
      if(wtctr==0)then
        if(p(1)>0. .and. zmin>0.)then ! water table at ground surface and zmin is below surface.  Added 5/3/2013, RLB 
          wtctr=wtctr+1;wctr=wctr+1
          wptr(wtctr)=1
          zwat(1)=elev(i);p(1)=0.
        else 
          write(*,*) 'Error in svijz(): wtctr not initialized!'
          write(*,*) 'wptr(wtctr), wtctr', wptr(wtctr), wtctr
        endif
      endif
      select case (el_or_dep) ! water table elevation or depth
        case('eleva')
          wtab(i+(jf-1)*imax)=zwat(wptr(wtctr)) ! use lowest computed water table
        case('depth')
          wtab(i+(jf-1)*imax)=elev(i)-zwat(wptr(wtctr))
        case default
          wtab(i+(jf-1)*imax)=zwat(wptr(wtctr)) ! Use lowest computed water table
      end select
    end if
    if(flag>-4 .or. flag<-7) return 
    select case (deepwat) ! pore pressure at deep node for SCOOPS 
      case('zero') ! zero pressure
        pdeep=0.
      case('flow') ! pressure consistent with user specified flow direction, use lowest computed water table
  !      pdeep=beta*(zwat(wptr(wtctr))-zdeep)
        pdeep=beta*(elev(i)-depth(i)-zdeep) ! use initial depth to water table
      case('hydr') ! hydrostatic pressure option
  !      pdeep=(zwat(wptr(wtctr))-zdeep)
        pdeep=(elev(i)-depth(i)-zdeep) ! use initial depth to water table
      case('relh') ! hydrostatic pressure reduced by relative height and constant factor
        relhgt=(elev(i)-zmn(1))/(zmx(1)-zmn(1))
  !      pdeep=(zwat(wptr(wtctr))-zdeep)*(1.-relhgt/3.)
        pdeep=(elev(i)-depth(i)-zdeep)*(1.-relhgt/3.) ! use initial depth to water table
      case default
        pdeep=0.
    end select
! Save data in ijz format, full (flag==-4) or downsampled (flag==-5) depth profile
    if(flag==-4 .or. flag==-5) then
      if(spcg<1) spcg=1
      dinc=1; if(flag==-5) dinc=spcg !dinc used for downsampled ijz file.  !Added spcg 2/14/2012 RLB
      fdinc=float(dinc)
      z=ztop
      if(dcf>0. .and. unsat(zo(i))) then !Unsaturated infiltration model
        do m=1,nzs+1,dinc
          write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),z,p(m),chi(m)
          if(zwat(1) < z .and. zwat(1) > (z-fdinc*zinc) .and. p(m)<0.) then ! eliminate repeated entries for water table
              if(zwat(1) > (zbot*1.00001)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & ix(i),jy(i),zwat(1),0.d0,1.d0
          end if
          z=z-(zinc*fdinc)
        end do
        if(zwat(1) < (zbot/1.00001)) then ! water table below zmax(i)
           write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(1),0.d0,1.d0 
        end if
        if(ldeep) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zdeep,pdeep,chdeep
      else if (dcf<=0.) then ! Saturated infiltration model
        do m=1,nzs+1,dinc
          write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),z,p(m),chi(m)
          if(dinc==1) then
            if(wtptr(m)==1) then
              write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & ix(i),jy(i),zwat(m),0.d0,chi(m)
            end if
          else if(dinc>1) then
            do m1=m, m+dinc-1
              if(wtptr(m1)==1) then
                write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
                & ix(i),jy(i),zwat(m1),0.d0,1.d0
              end if
            end do
          end if
          z=z-(zinc*fdinc)
        end do
    !    if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i)
    !       write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
    !        & ix(i),jy(i),zwat(nzs+1),0.d0,1.d0
    !    end if
        if(ldeep) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zdeep,pdeep,chdeep
      end if
    end if  
    if(flag==-6) then ! Sparce output option
    ! find points at maximum (+) and minimum (-) excursions from linear trend between ground surface and basal water table
      lopmn=.false.; lopmx=.false.; loxlo=.false.
      plin=0.d0; rp=0.d0
      mp=ptop/(ztop-zwat(nzs+1))
      z=ztop
      plin(1)=ptop
      do m=2,nzs+1
        z=z-zinc
        if(z<zwat(1) .and. dcf>0. .and. unsat(zo(i))) exit
        plin(m)=ptop-mp*(ztop-z)
        rp(m)=p(m)-plin(m)
      end do
      rpmin=minval(rp); mrpn=minloc(rp)
      rpmax=maxval(rp); mrpx=maxloc(rp)
      tol1=abs(p(mrpn(1))/100.)
      tol2=abs(p(mrpx(1))/100.)
      zrpn=ztop-zinc*float(mrpn(1)-1)
      zrpx=ztop-zinc*float(mrpx(1)-1)
      if(abs(rpmin) >= tol1 .and. rpmin<0.d0 .and. zrpn>zwat(1)) lopmn=.true.
      if(abs(rpmax) >= tol2 .and. rpmax>0.d0 .and. zrpx>zwat(1)) lopmx=.true.
      mfrst=mrpn(1);mlast=mrpx(1)
      if(mrpn(1)>mrpx(1)) then
        loxlo=.true.
        mfrst=mrpx(1);mlast=mrpn(1)
      end if
    ! Save data in ijz format
      if(dcf>0. .and. unsat(zo(i))) then 
        write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),ztop,ptop,chtop !,'top1' 
        if(loxlo) then
            if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zrpx,p(mrpx(1)),chi(mrpx(1)) !,'max1'
            if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zrpn,p(mrpn(1)),chi(mrpn(1)) !,'min1'
        else
            if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zrpn,p(mrpn(1)),chi(mrpn(1)) !,'min2'
            if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zrpx,p(mrpx(1)),chi(mrpx(1)) !,'max2'
        end if
        if(zwat(1) < ztop .and. ptop < 0.) then ! eliminate repeated entries for water table Unsaturated model has only one water table
            if(zwat(1) > (zbot*1.00001)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(1),0.d0,1.d0 !,'watu1'
        end if
        write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zbot,pbot,chbot !,'bot1'
        if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i)
            write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(nzs+1),0.d0,1.d0 !,'watd1'
        end if
        if(ldeep ) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zdeep,pdeep,chdeep !,'dp1'
      else if (dcf<=0.d0) then ! Saturated infiltration model
        write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),ztop,ptop,1.d0 !,'top2'
        if(loxlo) then
          if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & ix(i),jy(i),ztop-zinc*float(mrpx(1)),p(mrpx(1)),chi(mrpx(1)) !,'max3'
          do m=mfrst, mlast-1 ! Check for wt between min & max 4/25/2012 RLB
            if(wtptr(m)==1 .or. wtptr(m)==2) then
              write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & ix(i),jy(i),zwat(m),0.d0,chi(m)
            end if
          end do
          if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & ix(i),jy(i),ztop-zinc*float(mrpn(1)),p(mrpn(1)),chi(mrpn(1)) !,'min3'
        else
          if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & ix(i),jy(i),ztop-zinc*float(mrpn(1)),p(mrpn(1)),chi(mrpn(1)) !,'min4'
            do m=mfrst, mlast-1 ! Check for wt between min & max 4/25/2012 RLB
              if(wtptr(m)==1 .or. wtptr(m)==2) then
                write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
                & ix(i),jy(i),zwat(m),0.d0,1.d0
              end if
            end do
            if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),ztop-zinc*float(mrpx(1)),p(mrpx(1)),chi(mrpx(1)) !,'max4'
        end if
        do m=mlast, nzs ! find lower water table
          if(wtptr(m)==1 .or. wtptr(m)==2) then !
            if (zwat(m) > (zbot*1.00001)) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & ix(i),jy(i),zwat(m),0.d0,1.d0 !,'watu2'
          end if    
        end do
        write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zbot,pbot,1.d0 !,'bot2'
        if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i) 
          write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & ix(i),jy(i),zwat(nzs+1),0.d0,1.d0 !,'watd2'
        end if
        if(ldeep) write(uijz(jf),fmt='(2(i6,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & ix(i),jy(i),zdeep,pdeep,chdeep !,'dp2'
      end if
    end if
    return
11  continue
    write(*,*) 'Error writing ijz output file, at step ',jf
    write(LIS_logunit,*) 'Error writing ijz output file, at step ',jf ! SY
    !close (LIS_logunit) ! SY: DO NOT close log file
    write(LIS_logunit,*) '11 in subroutine svijz()' ! SY
    write(*,*) '11 in subroutine svijz()' ! SY
    !stop '11 in subroutine svijz()'
    call LIS_endrun() ! SY
  end subroutine svijz


!   SY: Subroutine indentation performed for integration into LIS
!  Routine to prepare pressure head & relative saturation (chi) data for export to xmdv formats
!  Calling routine should set value of dcf if it is not initialized in main.
!  Rex L. Baum, 10 Feb 2010, based on svijz(), Lastest revision 4/27/2012
  !subroutine svxmdv(i,jf,delh,newdep,ulog) ! SY
  subroutine svxmdv(i,jf,delh,newdep) ! SY
    use LIS_logMod, only : LIS_endrun, LIS_logunit, &
        LIS_getNextUnitNumber, LIS_releaseUnitNumber, LIS_verify ! SY
    use grids
    use input_vars 
    use model_vars
    use input_file_defs
    implicit none
    !integer::i,jf,ulog,m,m1,nw,mfrst,mlast ! SY
    integer::i,jf,m,m1,nw,mfrst,mlast ! SY
    integer::mrpx(1),mrpn(1),dinc,wtctr,wctr,wptr(nzs+1),wtptr(nzs+1)
    real (double)::delh,newdep,zns,zinc,fdinc 
    real (double):: dwat,z,ptop,pbot,zwat(nzs+1),zbot,ztop,chtop,chbot,zwat0
    real (double):: zdeep,pdeep,chdeep,zrpn,zrpx,relhgt,x,y 
    real (double):: tol1,tol2,mp,plin(nzs+1),rp(nzs+1),rpmin,rpmax
    logical ::ldeep, lopmx,lopmn,loxlo
    ! compute water-table elevations from surface elevation and depths
    zbot=elev(i)-zmax(i)
    ztop=elev(i)-zmin
    zwat0=elev(i)-depth(i) ! Initial water table depth
    zwat=zwat0
    ! Check for deepz <= 0 to skip deep node
    ldeep=.true.
    zdeep=elev(i)-deepz 
    if(zdeep>ztop) ldeep=.false.
    if(zdeep>zbot .and. ldeep) then
      zdeep=zbot-zmax(i) 
    end if
    zns=float(nzs) 
    zinc=(zmax(i)-zmin)/zns
    ptop=p(1)
    pbot=p(nzs+1)
    chtop=chi(1)
    chbot=chi(nzs+1)
    chdeep=1.
    ! estimate height of water table from results of unsaturated infiltration model
    if(dcf>0. .and. (unsat(zo(i)) .or. igcap(zo(i)))) then
      if(rikzero(i)>=0.0) then
        zwat(1)=elev(i)-(dusz-delh) ! statement in unsth() accounts for initial wt below zmax(i)
      else
        zwat(1)=elev(i)-newdep
      end if
      wtctr=1;wtptr(wtctr)=1;wptr(wtctr)=1  ! changed wtptr(wctr) to wtptr(wtctr) 18 April 2013, RLB
    else if (dcf<=0.) then
    ! find water table(s) from results of saturated infiltration model
      call dzero_brac(nzs,p,wtctr,wtptr) ! wtptr(j)=1 for intervals containing a water table
      chi=1.d0
      wctr=0;dwat=0.d0;wptr=0
      do nw=1,nzs+1
        if(wtptr(nw)>0) then
         if(wtptr(nw)==2) then ! water table at node
           wctr=wctr+1
           z=zmin+zinc*float(nw-1)
           dwat=z 
         end if
         if(wtptr(nw)==1) then ! water table between nodes
           wctr=wctr+1
           z=zmin+zinc*float(nw-1)
           dwat=z+zinc*(0.d0-p(nw))/(p(nw+1)-p(nw))
 !              write(*,*) 'i,nw,dwat ', i,nw,dwat
         end if
         if(rikzero(i)>=0.0) then
           zwat(nw)=elev(i)-dwat 
           wptr(wctr)=nw
         else
           zwat(nw)=elev(i)-newdep
           wptr(wctr)=nw
         end if
        end if
        if(wctr==wtctr) exit  ! all zero points accounted for
      end do
      if(p(nzs+1)<0.d0) then ! water table below bottom node (below zmax(i))
        wtptr(nzs+1)=1
        wtctr=wtctr+1;wctr=wctr+1
        zwat(nzs+1)=elev(i)-depth(i)+beta*(p(nzs+1)-p0zmx)  ! Approximate formula RLB 5/3/2013.
        if(zwat(nzs+1)<(elev(i)-depth(i))) zwat(nzs+1)=elev(i)-depth(i) ! prevents water table from falling below initial value
        wptr(wtctr)=nzs+1
      end if
    end if
    if(outp(1)) then
      if(wtctr==0)then
        if(p(1)>0. .and. zmin>0.)then ! water table at ground surface and zmin is below surface.  Added 5/3/2013, RLB 
          wtctr=wtctr+1;wctr=wctr+1
          wptr(wtctr)=1
          zwat(1)=elev(i);p(1)=0.
        else 
          write(*,*) 'Error in svxmdv(): wtctr not initialized!'
          write(*,*) 'wptr(wtctr), wtctr', wptr(wtctr), wtctr
        endif
      endif
      select case (el_or_dep) ! water table elevation or depth
        case('eleva')
         wtab(i+(jf-1)*imax)=zwat(wptr(wtctr)) ! use lowest computed water table
        case('depth')
         wtab(i+(jf-1)*imax)=elev(i)-zwat(wptr(wtctr))
        case default
         wtab(i+(jf-1)*imax)=zwat(wptr(wtctr)) ! Use lowest computed water table
      end select
    end if
    if(flag>-7 .or. flag<-9) return 
    select case (deepwat) ! pore pressure at deep node for SCOOPS 
      case('zero') ! zero pressure
        pdeep=0.
      case('flow') ! pressure consistent with user specified flow direction, use lowest computed water table
  !      pdeep=beta*(zwat(wptr(wtctr))-zdeep)
        pdeep=beta*(elev(i)-depth(i)-zdeep) ! use initial depth to water table
      case('hydr') ! hydrostatic pressure option
  !      pdeep=(zwat(wptr(wtctr))-zdeep)
        pdeep=(elev(i)-depth(i)-zdeep) ! use initial depth to water table
      case('relh') ! hydrostatic pressure reduced by relative height and constant factor
        relhgt=(elev(i)-zmn(1))/(zmx(1)-zmn(1))
  !      pdeep=(zwat(wptr(wtctr))-zdeep)*(1.-relhgt/3.)
        pdeep=(elev(i)-depth(i)-zdeep)*(1.-relhgt/3.) ! use initial depth to water table
      case default
        pdeep=0.
    end select
! x,y positions at centers of grid cells
    x=xllc+(0.5d0+dfloat(ix(i)-1))*celsiz
    y=yllc+(0.5d0+dfloat(jy(i)-1))*celsiz
! save data in XMDV file format, full depth profile (flag==-7) or downsampled (flag==-8)   
    if(flag==-7 .or. flag==-8) then
      if(spcg<1) spcg=1
      dinc=1; if(flag==-8) dinc=spcg ! dinc used for downsampled xmdv file ! Added spcg 2/14/2012 RLB
      fdinc=float(dinc)
      z=ztop
      if(dcf>0. .and. unsat(zo(i))) then !Unsaturated infiltration model
        do m=1,nzs+1,dinc
          write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & x,y,z,p(m),chi(m)
          if(zwat(1) < z .and. zwat(1) > (z-fdinc*zinc) .and. p(m)<0.) then ! eliminate repeated entries for water table
            if(zwat(1) > (zbot*1.00001)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & x,y,zwat(1),0.d0,1.d0
          end if
          z=z-(zinc*fdinc)
        end do
        if(zwat(1) < (zbot/1.00001)) then ! water table below zmax(i)
          write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
           & x,y,zwat(1),0.d0,1.d0 
        end if
        if(ldeep) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,zdeep,pdeep,chdeep
      else if (dcf<=0.) then ! Saturated infiltration model
        do m=1,nzs+1,dinc
          write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & x,y,z,p(m),chi(m)
          if(dinc==1) then
            if(wtptr(m)==1) then
              write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & x,y,zwat(m),0.d0,chi(m)
            end if
          else if(dinc>1) then
            do m1=m, m+dinc-1
              if(wtptr(m1)==1) then
                write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
                & x,y,zwat(m1),0.d0,1.d0
              end if
            end do
          end if
          z=z-(zinc*fdinc)
        end do
!    if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i)
!       write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
!        & x,y,zwat(nzs+1),0.d0,1.d0
!    end if
        if(ldeep) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,zdeep,pdeep,chdeep
      end if
    end if  
    if(flag==-9) then ! save data in "sparse" XMDV file format.
! find points at maximum (+) and minimum (-) excursions from linear trend between ground surface and water table
      lopmn=.false.; lopmx=.false.; ;loxlo=.false.
      plin=0.d0; rp=0.d0
      mp=ptop/(ztop-zwat(nzs+1))
      z=ztop
      plin(1)=ptop
      do m=2,nzs+1
        z=z-zinc
        if(z<zwat(1) .and. dcf>0. .and. unsat(zo(i))) exit
        plin(m)=ptop-mp*(ztop-z)
        rp(m)=p(m)-plin(m)
      end do
      rpmin=minval(rp); mrpn=minloc(rp)
      rpmax=maxval(rp); mrpx=maxloc(rp)
      tol1=abs(p(mrpn(1))/100.)
      tol2=abs(p(mrpx(1))/100.)
      zrpn=ztop-zinc*float(mrpn(1)-1)
      zrpx=ztop-zinc*float(mrpx(1)-1)
      if(abs(rpmin) >= tol1 .and. rpmin<0.d0 .and. zrpn>zwat(1)) lopmn=.true.
      if(abs(rpmax) >= tol2 .and. rpmax>0.d0 .and. zrpx>zwat(1)) lopmx=.true.
      mfrst=mrpn(1);mlast=mrpx(1)
      if(mrpn(1)>mrpx(1)) then
        loxlo=.true.
        mfrst=mrpx(1);mlast=mrpn(1)
      end if
! Save data in XMDV format
      if(dcf>0. .and. unsat(zo(i))) then 
        write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,ztop,ptop,chtop !,'top1' 
        if(loxlo) then
          if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zrpx,p(mrpx(1)),chi(mrpx(1)) !,'max1'
          if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zrpn,p(mrpn(1)),chi(mrpn(1)) !,'min1'
        else
          if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zrpn,p(mrpn(1)),chi(mrpn(1)) !,'min2'
          if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zrpx,p(mrpx(1)),chi(mrpx(1)) !,'max2'
        end if
        if(zwat(1) < ztop .and. ptop < 0.) then ! eliminate repeated entries for water table Unsaturated model has only one water table
          if(zwat(1) > (zbot*1.00001)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zwat(1),0.d0,1.d0 !,'watu1'
        end if
        write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,zbot,pbot,chbot !,'bot1'
        if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i)
          write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zwat(nzs+1),0.d0,1.d0 !,'watd1'
        end if
        if(ldeep ) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,zdeep,pdeep,chdeep !,'dp1'
      else if (dcf<=0.d0) then ! Saturated infiltration model
        write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,ztop,ptop,1.d0 !,'top2'
        if(loxlo) then
          if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,ztop-zinc*float(mrpx(1)),p(mrpx(1)),chi(mrpx(1)) !,'max3'
          do m=mfrst, mlast-1 ! Check for wt between min & max 4/25/2012 RLB
            if(wtptr(m)==1 .or. wtptr(m)==2) then
              write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & x,y,zwat(m),0.d0,chi(m)
            end if
          end do
          if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,ztop-zinc*float(mrpn(1)),p(mrpn(1)),chi(mrpn(1)) !,'min3'
        else
          if(lopmn .and. mrpn(1)>1 .and. mrpn(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,ztop-zinc*float(mrpn(1)),p(mrpn(1)),chi(mrpn(1)) !,'min4'
          do m=mfrst, mlast-1 ! Check for wt between min & max 4/25/2012 RLB
            if(wtptr(m)==1 .or. wtptr(m)==2) then
              write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
              & x,y,zwat(m),0.d0,1.d0
            end if
          end do
          if(lopmx .and. mrpx(1)>1 .and. mrpx(1)<(nzs+1)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,ztop-zinc*float(mrpx(1)),p(mrpx(1)),chi(mrpx(1)) !,'max4'
        end if
        do m=mlast, nzs ! find lower water table
          if(wtptr(m)==1 .or. wtptr(m)==2) then !
            if (zwat(m) > (zbot*1.00001)) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
            & x,y,zwat(m),0.d0,1.d0 !,'watu2'
          end if    
        end do
        write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,zbot,pbot,1.d0 !,'bot2'
        if(zwat(nzs+1) < (zbot/1.00001)) then ! water table below zmax(i) 
          write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
          & x,y,zwat(nzs+1),0.d0,1.d0 !,'watd2'
        end if
        if(ldeep) write(uijz(jf),fmt='(2(g18.9,1x),g15.5,1x,2(g11.4,1x))',err=11)&
        & x,y,zdeep,pdeep,chdeep!,'dp2'
      end if
    end if
    return
11  continue
    write(*,*) 'Error writing XMDV output file, at step ',jf
    write(LIS_logunit,*) 'Error writing XMDV output file, at step ',jf ! SY
    !close (LIS_logunit) ! SY: DO NOT close log file
    write(LIS_logunit,*) '11 in subroutine svxmdv()' ! SY
    write(*,*) '11 in subroutine svxmdv()' ! SY
    !stop '11 in subroutine svxmdv()' ! SY
    call LIS_endrun() ! SY
  end subroutine svxmdv


!   SY: Subroutine indentation performed for integration into LIS
  subroutine dzero_brac(nstp,x,zroctr,zrptr)
! Procedure to bracket locations of zero in a list of values using change of sign between values
! Rex L Baum, USGS, 2 March 2012
    implicit none
    integer, parameter:: double=kind(1d0)
    integer:: n,nstp,zroctr,zrptr(nstp+1)
    real (double):: x(nstp+1)
    zroctr=0; zrptr=0
    do n=1,nstp+1 ! Check first for zero values in list
      if(x(n)==0.d0) then
        zroctr=zroctr+1
        zrptr(n)=2
      end if
    end do
    do n=1,nstp  ! Now check for zero values between points in list
      if(x(n)*x(n+1)<0.d0) then
        zroctr=zroctr+1
        zrptr(n)=1
      else
        cycle
      end if
    end do
  end subroutine dzero_brac


end module TRIGRSMod
