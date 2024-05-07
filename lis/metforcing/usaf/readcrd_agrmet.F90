!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

#include "LIS_misc.h"

!BOP
!
! !ROUTINE: readcrd_agrmet
! \label{readcrd_agrmet}
!
! !REVISION HISTORY:
! 29Jul2005; Sujay Kumar, Initial Code
! 14Aug2008: Chris Franks, Add code to get precipitaion ob blacklist file
! 04Mar2010: Ryan Ruhge, Added option to put timestamp only on gfs directory
! 31 MAR 2010 added various resolution options - "8th", "16th", "64th",
!             latlon and associated masks.................Michael Shaw/WXE 
! 10 JUL 2010 added max_sfcobs and max_pcpobs.....Chris Franks/16WS/WXE/SEMS
! 15 NOV 2010 removed retreival of blacklist filename and change references 
!             from CDMS to JMOBS..................Chris Franks/16WS/WXE/SEMS
! 23 JUN 2011 Added the root of the input LIS files for retrospective mode
!             ....................................Chris Franks/16WS/WXE/SEMS
! 13 May 2013 added options for reading CMORPH data
!             ......................................Ryan Ruhge/16WS/WXE/SEMS
! 28 Aug 2018 Added IMERG...........................Eric Kemp/NASA/SSAI
! 21 Feb 2020 added support for 10-km GALWEM........Eric Kemp/NASA/SSAI
! 05 Mar 2020 added support for new GFS filename version...Eric Kemp/NASA/SSAI
!
! !INTERFACE:    
subroutine readcrd_agrmet()
! !USES:
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_config, LIS_masterproc
  use LIS_logMod,     only : LIS_logunit, LIS_verify, LIS_abort, &
       LIS_endrun
#if (defined SPMD)
  use LIS_mpiMod
#endif
  use LIS_pluginIndices, only : LIS_agrmetrunId
  use AGRMET_forcingMod, only : agrmet_struc

  implicit none
!
! !DESCRIPTION:
!
!  This routine reads the options specific to AGRMET algorithms from 
!  the LIS configuration file. 
!  
!EOP
  integer:: n,rc
  character(len=10)       :: cdate
  character(len=255) :: message(20) ! EMK
  real :: tmp_max_dist ! EMK
  character(len=201) :: c_string ! EMK
  integer :: ios ! EMK
  integer, external :: LIS_create_subdirs
  integer :: tmp_imerg_plp_thresh
  integer :: ierr
  logical :: use_nrt_bias_files ! EMK

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET forcing directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%agrmetdir,rc=rc)
     write(LIS_logunit,*)'[INFO] Using AGRMET forcing'
     write(LIS_logunit,*) '[INFO] AGRMET forcing directory: ', trim(agrmet_struc(n)%agrmetdir)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET first guess source:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%first_guess_source,rc=rc)
     call LIS_verify(rc,'AGRMET first guess source: option not specified in the config file')
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET analysis directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%analysisdir,rc=rc)
      write(unit=cdate, fmt='(a2,i2.2)') '.d',n      
     agrmet_struc(n)%analysisdir = trim(agrmet_struc(n)%analysisdir)&
          //trim(cdate)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET surface fields directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%sfcalcdir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET merged precip directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mrgpcpdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET cloud data directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%clouddir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GFS data directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%gfsdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GALWEM data directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%galwemdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET SSMI data directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%ssmidir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GEOPRECIP data directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%geodir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET JMOBS data directory:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cdmsdir,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use timestamp on directories:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%use_timestamp,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use timestamp on gfs:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%gfs_timestamp,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 8th polar mask file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%maskfile8,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 8th polar terrain file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%terrainfile8,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 16th polar mask file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%maskfile16,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 64th polar mask file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%maskfile64,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 16th polar terrain file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%terrainfile16,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 64th polar terrain file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%terrainfile64,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET latlon mask file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%maskfilell,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET sfcalc cntm file:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%sfcntmfile,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET precip climatology:",rc=rc)

  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%climodir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET nogaps wind weight:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%wndwgt,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET minimum wind speed:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%minwnd,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use present/past weather estimate:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%pwswch,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use precip observations:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%pcpobswch,rc=rc)
  enddo

  ! EMK...Precip observation file formats
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET precip obs file format:", &
       rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%pcpobsfmt, rc=rc)
     if (agrmet_struc(n)%pcpobsfmt .ne. 1 .and. &
          agrmet_struc(n)%pcpobsfmt .ne. 2) then
        write(LIS_logunit,*) &
             "[ERR] Bad 'AGRMET precip obs file format:' option"
        write(LIS_logunit,*) &
             '[ERR] Expected 1 or 2, found ', agrmet_struc(n)%pcpobsfmt
        write(LIS_logunit,*) '[ERR] Aborting...'

        flush(LIS_logunit)
        message(1) = &
             '[ERR] Illegal value for AGRMET precip obs file format'
#if (defined SPMD)
        call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
        if (LIS_masterproc) then
           call LIS_abort(message)
        else
           call sleep(10)
           call LIS_endrun()
        end if
     end if
  enddo

  ! EMK...Sfc observation file formats
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET sfc obs file format:", &
       rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%sfcobsfmt, rc=rc)
     if (agrmet_struc(n)%sfcobsfmt .ne. 1 .and. &
          agrmet_struc(n)%sfcobsfmt .ne. 2) then
        write(LIS_logunit,*) &
             "[ERR] Bad 'AGRMET sfc obs file format:' option"
        write(LIS_logunit,*) &
             '[ERR] Expected 1 or 2, found ', agrmet_struc(n)%sfcobsfmt
        write(LIS_logunit,*) '[ERR] Aborting...'

        flush(LIS_logunit)
        message(1) = &
             '[ERR] Illegal value for AGRMET sfc obs file format'
#if (defined SPMD)
        call MPI_Barrier(LIS_MPI_COMM, ierr)
#endif
        if (LIS_masterproc) then
           call LIS_abort(message)
        else
           call sleep(10)
           call LIS_endrun()
        end if
     end if
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET native imax:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imaxnative,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET native jmax:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%jmaxnative,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use CMORPH data:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmorswch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH minimum temperature threshold:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmorminthresh,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH maximum temperature threshold:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmormaxthresh,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GEO_PRECIP minimum temperature threshold:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%geominthresh,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GEO_PRECIP maximum temperature threshold:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%geomaxthresh,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH data directory:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmordir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH imax:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imaxcmor,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH jmax:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%jmaxcmor,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH min lat:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmorminlat,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH max lat:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmormaxlat,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH min lon:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmorminlon,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH max lon:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmormaxlon,rc=rc)
  enddo
    call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH dx:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmordx,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CMORPH dy:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cmordy,rc=rc)
  enddo

!EMK Add IMERG
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use IMERG data:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imerg_swch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET IMERG temperature threshold:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imerg_t_thresh,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET IMERG data directory:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imerg_dir,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET IMERG product:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imerg_product, &
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET IMERG version:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imerg_version, &
          rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config, &
       "AGRMET IMERG Probability Liquid Precip Threshold:",rc=rc)
  do n=1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          tmp_imerg_plp_thresh, &
          rc=rc)
     if (tmp_imerg_plp_thresh .gt. 100 .or. &
          tmp_imerg_plp_thresh .lt. 0) then
        write(LIS_logunit,*) &
             '[ERR] Illegal value for IMERG Probability Liquid Precip Threshold'
        write(LIS_logunit,*)'ABORTING!'
        flush(LIS_logunit)
        message(1) = &             
             '[ERR] Illegal value for IMERG Probability Liquid Precip Threshold'           
#if (defined SPMD)
        call MPI_Barrier(LIS_MPI_COMM, ierr)           
#endif    
        if (LIS_masterproc) then
           call LIS_abort(message)          
        else
           call sleep(10)
           call LIS_endrun()
        end if
     end if
     ! 4-byte integer to 2-byte integer
     agrmet_struc(n)%imerg_plp_thresh = tmp_imerg_plp_thresh
  enddo

! EMK Add OBA
  agrmet_struc(:)%oba_switch = 0
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET output OBA data:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%oba_switch,rc=rc)
  end do
  agrmet_struc(:)%skip_backqc = 0
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET skip backQC:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%skip_backqc,rc=rc)
  end do
  agrmet_struc(:)%skip_superstatqc = 0
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET skip superstatQC:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%skip_superstatqc,rc=rc)
  end do


  agrmet_struc(:)%gfsprecswch = 0
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use GFS precip:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%gfsprecswch,rc=rc)
  enddo
  agrmet_struc(:)%galwemprecswch = 0
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use GALWEM precip:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%galwemprecswch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use SSMI data:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%raswch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET SSMI imax:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imaxsmi,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET SSMI jmax:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%jmaxsmi,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use CDFSII-based estimate:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cdfs2swch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use GEOPRECIP estimate:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%geoswch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GEOPRECIP imax:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%imaxgp,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET GEOPRECIP jmax:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%jmaxgp,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET CDFSII time interval:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cdfs2int,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET use precip climatology:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%clswch,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET SSMI zero use switch:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%razero,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET snow distribution shape parameter:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%salp,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET alternate monthly weighting factor:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%clmult,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET minimum 3hr climo value:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mnpr,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET maximum 3hr climo value:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mxpr,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET minimum precip-per-precip day multiplier:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mnpd,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET maximum precip-per-precip day multiplier:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mxpd,rc=rc)
  enddo
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET cloud threshold to generate CDFSII estimate:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%cldth,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET median cloud cover percentage1:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mdmv,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET median cloud cover percentage2:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%mdpe,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET overcast percentage:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%ovpe,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET 3hr maximum precip ceiling:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%pcap,rc=rc)
  enddo
  
  call ESMF_ConfigFindLabel(LIS_config,"AGRMET maximum surface obs:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%max_sfcobs,&
          default=25000,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET maximum precip obs:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%max_pcpobs,&
          default=50000,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET retrospective root filename:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%retroFileRoot,rc=rc)
  enddo

  call ESMF_ConfigFindLabel(LIS_config,"AGRMET radiation derived from:",rc=rc)
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,agrmet_struc(n)%compute_radiation,default="cloud types", rc=rc)
  enddo

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%security_class,&
       label="AGRMET security classification:",rc=rc)

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%distribution_class,&
       label="AGRMET distribution classification:",rc=rc)

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%data_category,&
       label="AGRMET data category:",rc=rc)

  call ESMF_ConfigGetAttribute(LIS_config,LIS_rc%area_of_data,&
       label="AGRMET area of data:",rc=rc)

  !EMK NEW...Add map info of AGRMET forcing data if not running in ops mode
  if (LIS_rc%runmode .ne. LIS_agrmetrunId) then

     call ESMF_ConfigFindLabel(LIS_config,"AGRMET forcing map projection:", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing map projection: option not specified in the config file')
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingMapProj,rc=rc)
        call LIS_verify(rc, &
          '[ERR] AGRMET forcing map projection: value not specified in the config file')
        if ( trim(agrmet_struc(n)%forcingMapProj) .ne. "latlon") then
           message(1) = &
                '[ERR], AGRMET forcing map projection must be latlon!'
           write(LIS_logunit,*) message(1)
           call LIS_abort(message)          
        end if
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET forcing domain lower left lat:", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain lower left lat: option not specified in the config file')

     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingLowerLeftLat,rc=rc)
        call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain lower left lat: value not specified in the config file')
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET forcing domain lower left lon:", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain lower left lon: option not specified in the config file')
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingLowerLeftLon,rc=rc)
        call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain lower left lon: value not specified in the config file')
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET forcing domain upper right lat:", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain upper right lat: option not specified in the config file')
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingUpperRightLat,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain upper right lat: value not specified in the config file')
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET forcing domain upper right lon:", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain upper right lon: option not specified in the config file')
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingUpperRightLon,rc=rc)
        call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain upper right lon: value not specified in the config file')
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET forcing domain resolution (dx):", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain resolution (dx): option not specified in the config file')
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingResDx,rc=rc)
        call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain resolution (dx): value not specified in the config file')
     enddo

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET forcing domain resolution (dy):", &
          rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain resolution (dy): option not specified in the config file')
     do n=1,LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%forcingResDy,rc=rc)
        call LIS_verify(rc, &
          '[ERR] AGRMET forcing domain resolution (dy): value not specified in the config file')
     enddo

  end if ! EMK NEW

  ! EMK NEW...Error statistics for Bratseth precip analysis using GALWEM
  ! background.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip background error scale length (m): option not specified in the config file')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM Precip background error scale length (m): value not specified in the config file')
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip background error variance: option not specified in the config file')
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM Precip background error variance: value not specified in the config file')
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip Gauge observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip Gauge observation error variance: option not specified in the config file') 
 do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_gauge_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM Precip Gauge observation error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip GEOPRECIP observation error scale length (m):",rc=rc)
    call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM Precip GEOPRECIP observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_geoprecip_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM Precip GEOPRECIP observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip GEOPRECIP observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip GEOPRECIP observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_geoprecip_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM Precip GEOPRECIP observation error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip SSMI observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip SSMI observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_ssmi_err_scale_length,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip SSMI observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip SSMI observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip SSMI observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_ssmi_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip SSMI observation error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip CMORPH observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip CMORPH observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_cmorph_err_scale_length,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip CMORPH observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip CMORPH observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip CMORPH observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_cmorph_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip CMORPH observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK IMERG support
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip IMERG observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip IMERG observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_imerg_err_scale_length,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip IMERG observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM Precip IMERG observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip IMERG observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_precip_imerg_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM Precip IMERG observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK...Need to calculate maximum distance to compare or spread precip
  ! differences.  When using Gaussian function, this is 2 * largest scale 
  ! length (for correlation of ~0.02)
  do n=1,LIS_rc%nnest
     tmp_max_dist = max( &
          agrmet_struc(n)%galwem_precip_back_err_scale_length, &
          agrmet_struc(n)%galwem_precip_geoprecip_err_scale_length, &
          agrmet_struc(n)%galwem_precip_ssmi_err_scale_length, &
          agrmet_struc(n)%galwem_precip_cmorph_err_scale_length, &
          agrmet_struc(n)%galwem_precip_imerg_err_scale_length)
     agrmet_struc(n)%galwem_precip_max_dist = 2*tmp_max_dist
  end do ! n

  ! EMK NEW...Error statistics for Bratseth 2-m temperature analysis using 
  ! GALWEM background.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM T2M background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM T2M background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_t2m_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM T2M background error scale length (m): value not specified in the config file') 
     ! Maximum distance for spreading data for Gaussian function correlation
     ! is 2 * scale length (for ~0.02 correlation).
     agrmet_struc(n)%galwem_t2m_max_dist = &
          2*agrmet_struc(n)%galwem_t2m_back_err_scale_length
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM T2M background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM T2M background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_t2m_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM T2M background error variance: value not specified in the config file')      
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM T2M station observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM T2M station observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_t2m_stn_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM T2M station observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK NEW...Error statistics for Bratseth 2-m relative humidity analysis 
  ! using GALWEM background.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM RH2M background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM RH2M background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_rh2m_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM RH2M background error scale length (m): value not specified in the config file') 
     ! Maximum distance for spreading data for Gaussian function correlation
     ! is 2 * scale length (for ~0.02 correlation).
     agrmet_struc(n)%galwem_rh2m_max_dist = &
          2*agrmet_struc(n)%galwem_rh2m_back_err_scale_length
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM RH2M background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM RH2M background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_rh2m_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM RH2M background error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM RH2M station observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM RH2M station observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_rh2m_stn_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM RH2M station observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK NEW...Error statistics for Bratseth 10-m wind speed analysis 
  ! using GALWEM background.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM SPD10M background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM SPD10M background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_spd10m_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM SPD10M background error scale length (m): value not specified in the config file') 
     ! Maximum distance for spreading data for Gaussian function correlation
     ! is 2 * scale length (for ~0.02 correlation).
     agrmet_struc(n)%galwem_spd10m_max_dist = &
          2*agrmet_struc(n)%galwem_spd10m_back_err_scale_length
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM SPD10M background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET SPD10M background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_spd10m_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET SPD10M background error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM SPD10M station observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM SPD10M station observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_spd10m_stn_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM SPD10M station observation error variance: value not specified in the config file') 
  enddo ! n
  
  ! EMK NEW...Error statistics for Bratseth precip analysis using GFS
  ! background. Note that GFS is used as emergency backup if GALWEM is
  ! not available.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip background error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip background error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip Gauge observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip Gauge observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_gauge_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip Gauge observation error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip GEOPRECIP observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip GEOPRECIP observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_geoprecip_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip GEOPRECIP observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip GEOPRECIP observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip GEOPRECIP observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_geoprecip_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip GEOPRECIP observation error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip SSMI observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip SSMI observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_ssmi_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip SSMI observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip SSMI observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip SSMI observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_ssmi_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip SSMI observation error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip CMORPH observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip CMORPH observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_cmorph_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip CMORPH observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip CMORPH observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip CMORPH observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_cmorph_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS Precip CMORPH observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK IMERG support
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip IMERG observation error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip IMERG observation error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_imerg_err_scale_length,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip IMERG observation error scale length (m): value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS Precip IMERG observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip IMERG observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_precip_imerg_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
       '[ERR] AGRMET GFS Precip IMERG observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK...Need to calculate maximum distance to compare or spread precip
  ! differences.  When using Gaussian function, this is 2 * largest scale 
  ! length (for correlation of ~0.02)
  do n=1,LIS_rc%nnest
     tmp_max_dist = max( &
          agrmet_struc(n)%gfs_precip_back_err_scale_length, &
          agrmet_struc(n)%gfs_precip_geoprecip_err_scale_length, &
          agrmet_struc(n)%gfs_precip_ssmi_err_scale_length, &
          agrmet_struc(n)%gfs_precip_cmorph_err_scale_length, &
          agrmet_struc(n)%gfs_precip_imerg_err_scale_length)
     agrmet_struc(n)%gfs_precip_max_dist = 2*tmp_max_dist
  end do ! n

  ! EMK NEW...Error statistics for Bratseth 2-m temperature analysis using 
  ! GFS background. Note that GFS is used as emergency backup if GALWEM is
  ! not available.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS T2M background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS T2M background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_t2m_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS T2M background error scale length (m): value not specified in the config file') 
     ! Maximum distance for spreading data for Gaussian function correlation
     ! is 2 * scale length (for ~0.02 correlation).
     agrmet_struc(n)%gfs_t2m_max_dist = &
          2*agrmet_struc(n)%gfs_t2m_back_err_scale_length
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS T2M background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS T2M background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_t2m_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS T2M background error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS T2M station observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS T2M station observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_t2m_stn_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS T2M station observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK NEW...Error statistics for Bratseth 2-m relative humidity analysis 
  ! using GFS background.  Note that GFS is used as emergency backup if
  ! GALWEM is not available.
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS RH2M background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS RH2M background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_rh2m_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS RH2M background error scale length (m): value not specified in the config file') 
     ! Maximum distance for spreading data for Gaussian function correlation
     ! is 2 * scale length (for ~0.02 correlation).
     agrmet_struc(n)%gfs_rh2m_max_dist = &
          2*agrmet_struc(n)%gfs_rh2m_back_err_scale_length
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS RH2M background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS RH2M background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_rh2m_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS RH2M background error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS RH2M station observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS RH2M station observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_rh2m_stn_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS RH2M station observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK NEW...Error statistics for Bratseth 10-m wind speed analysis 
  ! using GFS background.  Note that GFS is used as emergency backup if
  ! GALWEM is not available
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS SPD10M background error scale length (m):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS SPD10M background error scale length (m): option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_spd10m_back_err_scale_length,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS SPD10M background error scale length (m): value not specified in the config file') 
     ! Maximum distance for spreading data for Gaussian function correlation
     ! is 2 * scale length (for ~0.02 correlation).
     agrmet_struc(n)%gfs_spd10m_max_dist = &
          2*agrmet_struc(n)%gfs_spd10m_back_err_scale_length
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS SPD10M background error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS SPD10M background error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_spd10m_back_sigma_b_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS SPD10M background error variance: value not specified in the config file') 
  enddo ! n
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GFS SPD10M station observation error variance:",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GFS SPD10M station observation error variance: option not specified in the config file') 
  do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%gfs_spd10m_stn_sigma_o_sqr,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GFS SPD10M station observation error variance: value not specified in the config file') 
  enddo ! n

  ! EMK END Bratseth settings

  ! EMK...Set the GALWEM resolution.  Currently either 17-km or 10-km     
  call ESMF_ConfigFindLabel(LIS_config,&
       "AGRMET GALWEM nominal resolution (km):",rc=rc)
  call LIS_verify(rc, &
       '[ERR] AGRMET GALWEM nominal resolution (km): option not specified in the config file') 
   do n=1,LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config,&
          agrmet_struc(n)%galwem_res,rc=rc)
     call LIS_verify(rc, &
          '[ERR] AGRMET GALWEM nominal resolution (km): value not specified in the config file') 
  enddo ! n

  !EMK...Set the GFS filename version.  
  call ESMF_ConfigFindLabel(LIS_config, &
       "AGRMET GFS filename version:",rc=rc)
  call LIS_verify(rc, &
       "[ERR] AGRMET GFS filename version: not specified in config file")
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%gfs_filename_version, rc=rc)
     call LIS_verify(rc, &
          "[ERR] AGRMET GFS filename version: not specified in config file")
  enddo ! n

  ! EMK Add WWMCA GRIB1 option
  call ESMF_ConfigFindLabel(LIS_config, &
       "AGRMET WWMCA GRIB1 read option:",rc=rc)
  call LIS_verify(rc, &
       "[ERR] AGRMET WWMCA GRIB1 read option: not specified in config file")
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%read_wwmca_grib1, rc=rc)
     call LIS_verify(rc, &
          "[ERR] AGRMET WWMCA GRIB1 read option: not specified in config file")
  enddo ! n

  ! EMK Add background bias correction option
  call ESMF_ConfigFindLabel(LIS_config, &
       "AGRMET PPT Background bias correction option:",rc=rc)
  call LIS_verify(rc, &
       "[ERR] AGRMET PPT Background bias correction option: not specified in config file")
  do n = 1, LIS_rc%nnest
     call ESMF_ConfigGetAttribute(LIS_config, &
          agrmet_struc(n)%back_bias_corr, rc=rc)
     call LIS_verify(rc, &
          "[ERR] AGRMET PPT Background bias correction option: not specified in config file")
     if (agrmet_struc(n)%back_bias_corr .lt. 0 .or. &
          agrmet_struc(n)%back_bias_corr .gt. 2) then
        call LIS_verify(rc, &
             "[ERR] AGRMET PPT Background bias correction option: bad value in config file, set 0, 1, or 2")
     end if
     if (agrmet_struc(n)%back_bias_corr .eq. 1) then
        allocate(agrmet_struc(n)%pcp_back_bias_ratio(LIS_rc%gnc(n),LIS_rc%gnr(n)))
        agrmet_struc(n)%pcp_back_bias_ratio = 1.
        agrmet_struc(n)%pcp_back_bias_ratio_month = 0
     else if (agrmet_struc(n)%back_bias_corr .eq. 2) then

        allocate(agrmet_struc(n)%gfs_nrt_bias_ratio( &
             LIS_rc%gnc(n),LIS_rc%gnr(n)))
        allocate(agrmet_struc(n)%galwem_nrt_bias_ratio( &
             LIS_rc%gnc(n),LIS_rc%gnr(n)))
        agrmet_struc(n)%gfs_nrt_bias_ratio = 1.
        agrmet_struc(n)%galwem_nrt_bias_ratio = 1.
        agrmet_struc(n)%pcp_back_bias_ratio_month = 0
     end if
  end do

  ! EMK Add support for NRT bias files for GFS and GALWEM
  use_nrt_bias_files = .false.
  do n = 1, LIS_rc%nnest
     if (agrmet_struc(n)%back_bias_corr .eq. 2) then
        use_nrt_bias_files = .true.
        exit
     end if
  end do

  if (use_nrt_bias_files) then
     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET PPT GFS NRT bias file:", rc=rc)
     call LIS_verify(rc, &
          "[ERR] AGRMET PPT GFS NRT bias file: not specified in config file")

     do n = 1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%gfs_nrt_bias_ratio_file, rc=rc)
        call LIS_verify(rc, &
             "[ERR] AGRMET PPT GFS NRT bias file: not specified in config file")
     end do

     call ESMF_ConfigFindLabel(LIS_config, &
          "AGRMET PPT GALWEM NRT bias file:", rc=rc)
     call LIS_verify(rc, &
          "[ERR] AGRMET PPT GALWEM NRT bias file: not specified in config file")
     do n = 1, LIS_rc%nnest
        call ESMF_ConfigGetAttribute(LIS_config, &
             agrmet_struc(n)%galwem_nrt_bias_ratio_file, rc=rc)
        call LIS_verify(rc, &
             "[ERR] AGRMET PPT GALWEM NRT bias file: not specified in config file")
     end do
  end if

  do n=1,LIS_rc%nnest
     agrmet_struc(n)%radProcessInterval = 1
     agrmet_struc(n)%radProcessAlarmTime = 0.0
     
     agrmet_struc(n)%agrmettime1 = 3000.0
     agrmet_struc(n)%agrmettime2 = 0.0
          
     agrmet_struc(n)%pcpProcessInterval = 3
     agrmet_struc(n)%pcpProcessAlarmTime = 0.0

     agrmet_struc(n)%agrmetpcptime1 = 3000.0
     agrmet_struc(n)%agrmetpcptime2 = 0.0
     
     agrmet_struc(n)%readgfsInterval = 6
     agrmet_struc(n)%gfsAlarmTime = 0.0
! create the analysis directories. 
     
     !write(LIS_logunit,*)'EMK: Calling system in readcrd_agrmet...'
     !call system('mkdir -p '//trim(agrmet_struc(n)%analysisdir))
     c_string = trim(agrmet_struc(n)%analysisdir)
     ! EMK: Avoid race condition
     if (LIS_masterproc) then
        ios = LIS_create_subdirs(len_trim(c_string),trim(c_string))
        if (ios .ne. 0) then
           write(LIS_logunit,*)'[ERR] Cannot create directory ', &
                trim(agrmet_struc(n)%analysisdir)
           flush(LIS_logunit)
        end if
     end if

!EMK...Make sure directory has been created before continuing.
#if (defined SPMD)
     call mpi_barrier(LIS_mpi_comm, ierr)
#endif

  enddo
end subroutine readcrd_agrmet
