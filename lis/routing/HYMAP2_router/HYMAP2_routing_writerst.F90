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
#include "LIS_NetCDF_inc.h"

!BOP
! !ROUTINE: HYMAP2_routing_writerst
! \label{HYMAP2_routing_writerst} 
!
! !REVISION HISTORY:
! 15 Nov 2011: Augusto Getirana;  Initial implementation
! 19 Jan 2016: Augusto Getirana;  Inclusion of four Local Inertia variables
! 10 Mar 2019: Sujay Kumar;       Added support for NetCDF and parallel 
!                                 processing. 
! 27 Apr 2020: Augusto Getirana;  Added support for urban drainage
!  
! !INTERFACE: 
subroutine HYMAP2_routing_writerst(n)

!
! !DESCRIPTION: 
!  This routine writes restart files for HYMAP2. The restart files
!  are in NetCDF format. 
!
!EOP

  use ESMF
  use LIS_coreMod
  use LIS_logMod
  use LIS_fileIOMod
  use LIS_timeMgrMod
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use HYMAP2_routingMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
  use netcdf
#endif

  implicit none
  
  integer, intent(in)   :: n 
  
  character(len=LIS_CONST_PATH_LEN) :: filename
  integer               :: ftn
  integer               :: status
  logical               :: alarmCheck
  
  alarmCheck = LIS_isAlarmRinging(LIS_rc,&
          "HYMAP2 router restart alarm")
  if(alarmCheck.or.(LIS_rc%endtime ==1)) then 
     if(LIS_masterproc) then 
        call LIS_create_output_directory('ROUTING')
        call LIS_create_restart_filename(n,filename,&
             'ROUTING','HYMAP2_router',&
             wformat="netcdf")
        write(LIS_logunit,*) '[INFO] Writing routing restart ',trim(filename)
       
#if (defined USE_NETCDF4)
        status = nf90_create(path=filename, cmode=nf90_hdf5, ncid = ftn)
        call LIS_verify(status,"Error in nf90_open in HYMAP2_routing_writerst")
#endif
#if (defined USE_NETCDF3)
        status = nf90_create(Path = filename, cmode = nf90_clobber, ncid = ftn)
        call LIS_verify(status, "Error in nf90_open in HYMAP2_routing_writerst")
#endif
     endif
     call HYMAP2_dump_restart(n, ftn)
          
     if (LIS_masterproc) then
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
        status = nf90_close(ftn)
        call LIS_verify(status, "Error in nf90_close in HYMAP2_routing_writerst")
#endif
     endif
     
  endif

end subroutine HYMAP2_routing_writerst


!BOP
!
! !ROUTINE: HYMAP2_dump_restart
! \label{HYMAP2_dump_restart}
!
! !REVISION HISTORY:
!
! !INTERFACE:
subroutine HYMAP2_dump_restart(n, ftn)

! !USES:
    use LIS_coreMod, only : LIS_rc, LIS_masterproc
    use LIS_logMod, only  : LIS_logunit
    use LIS_historyMod
    use HYMAP2_routingMod

    implicit none

    integer, intent(in) :: ftn
    integer, intent(in) :: n

!
! !DESCRIPTION:
!  This routine gathers the necessary restart variables and performs
!  the actual write statements to create the restart files.
!
!  The arguments are:
!  \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    unit number for the restart file
!  \end{description}
!
! 
!EOP 
        
    integer :: dimID(1)
    integer :: rivsto_ID
    integer :: fldsto_ID
    integer :: rnfsto_ID
    integer :: bsfsto_ID
    integer :: rivout_pre_ID
    integer :: rivdph_pre_ID
    integer :: fldout_pre_ID
    integer :: flddph_pre_ID
    !ag (27Apr2020)
    integer :: drsto_ID
    integer :: drout_ID
    ! write the header of the restart file
    call HYMAP2_writeGlobalHeader_restart(ftn, n, &
         "HYMAP2", &
         dimID)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, rivsto_ID, "RIVSTO", &
         "river storage", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, fldsto_ID, "FLDSTO", &
         "flood storage", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, rnfsto_ID, "RNFSTO", &
         "runoff storage", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, bsfsto_ID, "BSFSTO", &
         "baseflow storage", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, rivout_pre_ID, "RIVOUT_PRE", &
         "rivoutpre", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, rivdph_pre_ID, "RIVDPH_PRE", &
         "river depth_pre", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, fldout_pre_ID, "FLDOUT_PRE", &
         "flood discharge pre", &
         "-", 1, -99999.0, 99999.0)
    call HYMAP2_writeHeader_restart(ftn, n, dimID, flddph_pre_ID, "FLDDPH_PRE", &
         "flood depth pre", &
         "-", 1, -99999.0, 99999.0)
!ag (27Apr2020)
    if(HYMAP2_routing_struc(n)%flowtype==4)then    
      call HYMAP2_writeHeader_restart(ftn, n, dimID, drsto_ID, "DRSTO", &
           "urban drainage storage", &
           "-", 1, -99999.0, 99999.0)
      call HYMAP2_writeHeader_restart(ftn, n, dimID, drout_ID, "DROUT", &
           "urban drainage discharge", &
           "-", 1, -99999.0, 99999.0)
    endif
    call HYMAP2_closeHeader_restart(ftn)
    
    ! write state variables into restart file
    if(HYMAP2_routing_struc(n)%useens.eq.0) then 
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%rivsto, &
            rivsto_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%fldsto, &
            fldsto_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%rnfsto, &
            rnfsto_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%bsfsto, &
            bsfsto_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%rivout_pre, &
            rivout_pre_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%rivdph_pre, &
            rivdph_pre_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%fldout_pre, &
            fldout_pre_ID)    
       call HYMAP2_writevar_restart(ftn,n, &
            HYMAP2_routing_struc(n)%flddph_pre, &
            flddph_pre_ID)    
       !ag (27Apr2020)
       !urban drainage storage   
       if(HYMAP2_routing_struc(n)%flowtype==4)then    
         call HYMAP2_writevar_restart(ftn,n, &
              HYMAP2_routing_struc(n)%drsto, &
              drsto_ID)    
         call HYMAP2_writevar_restart(ftn,n, &
              HYMAP2_routing_struc(n)%drout, &
              drout_ID)
       endif
    else

       call HYMAP2_writevar_restart_ens(ftn,n,&
            HYMAP2_routing_struc(n)%rivsto, &
            rivsto_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n,&
            HYMAP2_routing_struc(n)%fldsto, &
            fldsto_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n,&
            HYMAP2_routing_struc(n)%rnfsto, &
            rnfsto_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n,&
            HYMAP2_routing_struc(n)%bsfsto, &
            bsfsto_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n, &
            HYMAP2_routing_struc(n)%rivout_pre, &
            rivout_pre_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n, &
            HYMAP2_routing_struc(n)%rivdph_pre, &
            rivdph_pre_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n, &
            HYMAP2_routing_struc(n)%fldout_pre, &
            fldout_pre_ID)    
       call HYMAP2_writevar_restart_ens(ftn,n, &
            HYMAP2_routing_struc(n)%flddph_pre, &
            flddph_pre_ID)  
       !ag (27Apr2020)
       !urban drainage storage       
       if(HYMAP2_routing_struc(n)%flowtype==4)then    
         call HYMAP2_writevar_restart_ens(ftn,n, &
              HYMAP2_routing_struc(n)%drsto, &
              drsto_ID)    
         call HYMAP2_writevar_restart_ens(ftn,n, &
              HYMAP2_routing_struc(n)%drout, &
              drout_ID)
       endif  
    endif

  end subroutine HYMAP2_dump_restart

!BOP
! !ROUTINE: HYMAP2_writeGlobalHeader_restart
! \label{HYMAP2_writeGlobalHeader_restart}
! 
! !INTERFACE: HYMAP2_writeGlobalHeader_restart
  subroutine HYMAP2_writeGlobalHeader_restart(ftn,n,&
       model_name, dimID)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
    use HYMAP2_routingMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
    use netcdf
#endif
    implicit none

! !ARGUMENTS: 
    integer,   intent(in)     :: n 
    integer,   intent(in)     :: ftn
    character(len=*), intent(in) :: model_name
    integer                   :: dimID(1)

! 
! !DESCRIPTION: 
!  This routine writes an output file in the NETCDF format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!EOP

    integer, dimension(8) :: values
    character(len=8)  :: date
    character(len=10) :: time
    character(len=5)  :: zone
    character*20      :: wout

#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
    call date_and_time(date,time,zone,values)       
    if(LIS_masterproc) then 
       
       if(HYMAP2_routing_struc(n)%useens.eq.0) then 
          call LIS_verify(nf90_def_dim(ftn,'ntiles',&
               LIS_rc%glbnroutinggrid(n),&
               dimID(1)),&
               'nf90_def_dim failed for ntiles in HYMAP2_writeGlobalHeader_restart')
       else
          call LIS_verify(nf90_def_dim(ftn,'ntiles',&
               LIS_rc%glbnroutinggrid(n)*LIS_rc%nensem(n),&
               dimID(1)),&
               'nf90_def_dim failed for ntiles in HYMAP2_writeGlobalHeader_restart')
       endif
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", &
            LIS_rc%udef),'nf90_put_att failed for missing_value')
       
       
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
            "LIS land surface model restart"),&
            'nf90_put_att failed for title')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
            trim(LIS_rc%institution)),&
            'nf90_put_att failed for institution')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
            trim(model_name)),&
            'nf90_put_att failed for source')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
            'nf90_put_att failed for history')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
            "Kumar_etal_EMS_2006, Peters-Lidard_etal_ISSE_2007"),&
            'nf90_put_att failed for references')
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
            "CF-1.6"),'nf90_put_att failed for conventions') !CF version 1.6
       call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"),&
            'nf90_put_att failed for comment')
       
       if(LIS_rc%lis_map_proj.eq."latlon") then !latlon
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "EQUIDISTANT CYLINDRICAL"),&
               'nf90_put_att failed for MAP_PROJECTION')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,9)),&
               'nf90_put_att failed for DX')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,10)),&
               'nf90_put_att failed for DY')
       elseif(LIS_rc%lis_map_proj.eq."mercator") then 
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "MERCATOR"),&
               'nf90_put_att failed for MAP_PROJECTION')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LON') 
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)),&
               'nf90_put_att failed for TRUELAT1')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)),&
               'nf90_put_att failed for STANDARD_LON')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)),&
               'nf90_put_att failed for DX')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)),&
               'nf90_put_att failed for DY')
       elseif(LIS_rc%lis_map_proj.eq."lambert") then !lambert conformal
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "LAMBERT CONFORMAL"),&
               'nf90_put_att failed for MAP_PROJECTION')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)),&
               'nf90_put_att failed for TRUELAT1')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
               LIS_rc%gridDesc(n,7)),&
               'nf90_put_att failed for TRUELAT2')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)),&
               'nf90_put_att failed for STANDARD_LON')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)),&
               'nf90_put_att failed for DX')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)),&
               'nf90_put_att failed for DY')
          
       elseif(LIS_rc%lis_map_proj.eq."polar") then ! polar stereographic
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"),&
               'nf90_put_att failed for MAP_PROJECTION')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LAT", &
               LIS_rc%gridDesc(n,4)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,&
               "SOUTH_WEST_CORNER_LON", &
               LIS_rc%gridDesc(n,5)),&
               'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LIS_rc%gridDesc(n,10)),&
               'nf90_put_att failed for TRUELAT1')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
               LIS_rc%gridDesc(n,7)),&
               'nf90_put_att failed for ORIENT')
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LIS_rc%gridDesc(n,11)),&
               'nf90_put_att failed for STANDARD_LON')                  
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LIS_rc%gridDesc(n,8)),&
               'nf90_put_att failed for DX')                  
          call LIS_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LIS_rc%gridDesc(n,9)),&
               'nf90_put_att failed for DY')                  
       endif
    endif
#endif    
  end subroutine HYMAP2_writeGlobalHeader_restart

!BOP
! !ROUTINE: HYMAP2_writeHeader_restart
! \label{HYMAP2_writeHeader_restart}
! 
! !INTERFACE: 
  subroutine HYMAP2_writeHeader_restart(ftn,n,dimID, vid, standard_name, &
       long_name, units, vlevels, valid_min, valid_max)
! !USES: 
    use LIS_coreMod
    use LIS_logMod
#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    use netcdf
#endif
    implicit none

! !ARGUMENTS:     
    integer                    :: ftn
    integer                    :: n
    integer                    :: dimID(1)
    integer                    :: vid
    character(len=*)           :: standard_name
    character(len=*)           :: long_name
    character(len=*)           :: units
    integer                    :: vlevels
    real                       :: valid_min
    real                       :: valid_max

! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single HYMAP2 variable
! 
!   The arguments are: 
!   \begin{description}
!   \item[n]
!    index of the nest
!   \item[ftn]
!    NETCDF file unit handle
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!   \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LIS\_endrun](\ref{LIS_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LIS\_verify](\ref{LIS_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
    integer :: data_index
    integer :: shuffle, deflate, deflate_level
    integer :: dimID_t(2)
    integer :: fill_value

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    dimID_t(1) = dimID(1)
    fill_value = LIS_rc%udef

    if(LIS_masterproc) then 
       call LIS_verify(nf90_def_var(ftn,trim(standard_name),&
            nf90_float,dimids = dimID_t(1:1), varID=vid),&
            'nf90_def_var(1d) failed in HYMAP2_writeHeader_restart')
       
#if(defined USE_NETCDF4)
       call LIS_verify(nf90_def_var_fill(ftn,&
            vid, 1,fill_value), 'nf90_def_var_fill failed for '//&
            standard_name)
       
       call LIS_verify(nf90_def_var_deflate(ftn,&
            vid, shuffle, deflate, deflate_level),&
            'nf90_def_var_deflate(1d) failed in HYMAP2_writeHeader_restart')
#endif

       call LIS_verify(nf90_put_att(ftn,vid,&
            "units",trim(units)),&
            'nf90_put_att failed for units')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "standard_name",trim(standard_name)))
       call LIS_verify(nf90_put_att(ftn,vid,&
            "long_name",trim(long_name)),&
            'nf90_put_att failed for long_name')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "scale_factor",1.0),&
            'nf90_put_att failed for scale_factor')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "add_offset",0.0),&
            'nf90_put_att failed for add_offset')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "missing_value",LIS_rc%udef),&
            'nf90_put_att failed for missing_value')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "_FillValue",LIS_rc%udef),&
            'nf90_put_att failed for _FillValue')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "vmin",valid_min),&
            'nf90_put_att failed for vmin')
       call LIS_verify(nf90_put_att(ftn,vid,&
            "vmax",valid_max),&
            'nf90_put_att failed for vmax')
#endif
    endif
  end subroutine HYMAP2_writeHeader_restart
  

!BOP
! !ROUTINE: HYMAP2_closeHeader_restart
! \label{HYMAP2_closeHeader_restart}
! 
! !INTERFACE: 
  subroutine HYMAP2_closeHeader_restart(ftn)

    use LIS_coreMod
    use LIS_logMod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
    use netcdf
#endif

    implicit none

! !ARGUMENTS:     
    integer            :: ftn
! 
! !DESCRIPTION: 
!    This routine closes the required NETCDF header.
! 
!   The arguments are: 
!   \begin{description}
!   \item[ftn]
!    NETCDF file unit handle
!   \item[n]
!    index of the nest
!   \item[m]
!    index of the surface type
!   \item[dimID]
!    NETCDF dimension ID corresponding to the variable
!   \item[rstInterval]
!    the restart interval
!   \end{description}

!EOP

#if(defined USE_NETCDF3 || defined USE_NETCDF4)
       
    if(LIS_masterproc) then 
       call LIS_verify(nf90_enddef(ftn),&
            'nf90_enddef failed')       
    endif
#endif
  end subroutine HYMAP2_closeHeader_restart

!BOP
! !ROUTINE: HYMAP2_writevar_restart
! \label{HYMAP2_writevar_restart}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_restart
  subroutine HYMAP2_writevar_restart(ftn, n, var, varid)
! !USES:
    use LIS_coreMod
    use LIS_routingMod
    use LIS_mpiMod
    use LIS_logMod
    use HYMAP2_routingMod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
    use netcdf
#endif

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    real                :: var(LIS_rc%nroutinggrid(n))
    integer             :: varid

! !DESCRIPTION:
!  Writes a real variable to a NetCDF restart file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variable being written, dimensioned in the tile space
!   \item [varid]
!     index of the variable being written
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp1(:)
    integer           :: l,i,ix,iy,ix1,iy1,m
    integer :: ierr

    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbnroutinggrid(n)))
       allocate(gtmp1(LIS_rc%glbnroutinggrid(n)))
    else
       allocate(gtmp1(1))
    endif
#if (defined SPMD)     
    call MPI_GATHERV(var,LIS_rc%nroutinggrid(n),&
         MPI_REAL,gtmp1,&
         LIS_routing_gdeltas(n,:),&
         LIS_routing_goffsets(n,:),&
         MPI_REAL,0,LIS_mpi_comm,ierr)
#else 
       gtmp1 = var
#endif
    if(LIS_masterproc) then
       do l=1,LIS_npes
          do i=1,LIS_routing_gdeltas(n,l-1)
             ix = HYMAP2_routing_struc(n)%seqx_glb(i+&
                  LIS_routing_goffsets(n,l-1))
             iy = HYMAP2_routing_struc(n)%seqy_glb(i+&
                  LIS_routing_goffsets(n,l-1))
             ix1 = ix + LIS_ews_halo_ind(n,l) - 1
             iy1 = iy + LIS_nss_halo_ind(n,l)-1
             gtmp(HYMAP2_routing_struc(n)%sindex(ix1,iy1)) = &
                  gtmp1(i+LIS_routing_goffsets(n,l-1))
          enddo
       enddo
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
       ierr = nf90_put_var(ftn,varid,gtmp,(/1/),&
               (/LIS_rc%glbnroutinggrid(n)/))
       call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
#endif   
       deallocate(gtmp)       
    endif
    deallocate(gtmp1)
  end subroutine HYMAP2_writevar_restart


!BOP
! !ROUTINE: HYMAP2_writevar_restart_ens
! \label{HYMAP2_writevar_restart_ens}
! 
! !INTERFACE:
! Private name: call using LIS_writevar_restart_ens
  subroutine HYMAP2_writevar_restart_ens(ftn, n, var, varid)
! !USES:
    use LIS_coreMod
    use LIS_routingMod
    use LIS_mpiMod
    use LIS_logMod
    use HYMAP2_routingMod
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
    use netcdf
#endif

    implicit none
! !ARGUMENTS: 
    integer, intent(in) :: ftn
    integer, intent(in) :: n
    real                :: var(LIS_rc%nroutinggrid(n),LIS_rc%nensem(n))
    integer             :: varid

! !DESCRIPTION:
!  Writes a real variable to a NetCDF restart_ens file. 
!
!  The arguments are: 
!  \begin{description}
!   \item [n]
!     index of the domain or nest.
!   \item [ftn]
!     unit number of the binary output file
!   \item [var]
!     variable being written, dimensioned in the tile space
!   \item [varid]
!     index of the variable being written
!  \end{description}
!EOP
    real, allocatable :: gtmp(:)
    real, allocatable :: gtmp2(:)
    real, allocatable :: gtmp1(:,:)
    integer           :: m
    integer           :: l,i,ix,iy,ix1,iy1
    integer :: ierr

    if(LIS_masterproc) then 
       allocate(gtmp(LIS_rc%glbnroutinggrid(n)*LIS_rc%nensem(n)))
       allocate(gtmp1(LIS_rc%glbnroutinggrid(n),LIS_rc%nensem(n)))
    else
       allocate(gtmp1(1,LIS_rc%nensem(n)))
    endif
    do m=1,LIS_rc%nensem(n)  
#if (defined SPMD)    
       call MPI_GATHERV(var(:,m),LIS_rc%nroutinggrid(n),&
            MPI_REAL,gtmp1(:,m),&
            LIS_routing_gdeltas(n,:),&
            LIS_routing_goffsets(n,:),&
            MPI_REAL,0,LIS_mpi_comm,ierr)
#else
       gtmp1(:,m) = var(:,m)
#endif
    enddo
    if(LIS_masterproc) then
       do l=1,LIS_npes
          do i=1,LIS_routing_gdeltas(n,l-1)
             do m=1,LIS_rc%nensem(n)  

                ix = HYMAP2_routing_struc(n)%seqx_glb(i+&
                     LIS_routing_goffsets(n,l-1))
                iy = HYMAP2_routing_struc(n)%seqy_glb(i+&
                     LIS_routing_goffsets(n,l-1))
                ix1 = ix + LIS_ews_halo_ind(n,l) - 1
                iy1 = iy + LIS_nss_halo_ind(n,l)-1

                gtmp(m + (HYMAP2_routing_struc(n)%sindex(ix1,iy1) -1)* &
                     LIS_rc%nensem(n)) = &
                     gtmp1(i+LIS_routing_goffsets(n,l-1),m)
             enddo
          enddo
       enddo    
#if ( defined USE_NETCDF3 || defined USE_NETCDF4 )
       ierr = nf90_put_var(ftn,varid,gtmp,(/1/),&
            (/LIS_rc%glbnroutinggrid(n)*LIS_rc%nensem(n)/))
       call LIS_verify(ierr,'nf90_put_var failed in LIS_historyMod')
#endif   
       deallocate(gtmp)       
    endif
    deallocate(gtmp1)

  end subroutine HYMAP2_writevar_restart_ens
