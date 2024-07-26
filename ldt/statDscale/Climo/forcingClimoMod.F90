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
#include "LDT_NetCDF_inc.h"
module forcingClimoMod
!BOP
!
! !MODULE: forcingClimoMod
!
! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  2 Mar 2016: Sujay Kumar; Initial implementation
!
  use ESMF
  use LDT_coreMod
  use LDT_metforcingMod
  use LDT_timeMgrMod
  use LDT_logMod
  use LDT_fileIOMod
  use LDT_FORC_AttributesMod

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
  use netcdf
#endif

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: forcingClimo_init
  public :: forcingClimo_diagnose
  public :: forcingClimo_compute
  public :: forcingClimo_output
!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  type, public :: forcingClimoEntry

     real,    allocatable :: sx_mu(:,:,:)
     integer, allocatable :: sx_count(:,:,:)

  end type forcingClimoEntry

  type, public :: forcingClimo_struc
     logical                               :: startFlag
     integer                               :: ntimes
     integer                               :: nDataHours
     integer                               :: nVars
     integer,                 allocatable  :: varId(:)
     character*20,            allocatable  :: units(:)
     character*100,           allocatable  :: varName(:) 
     character*100,           allocatable  :: shortName(:) 
     type(forcingClimoEntry), allocatable  :: dataEntry(:)
  end type forcingClimo_struc

  type(forcingClimo_struc) :: LDT_forcingClimo_struc

!EOP

contains

  subroutine forcingClimo_init()

    integer             :: rc
    integer             :: n
    integer             :: i,k
    character*20        :: stime
    real                :: nsecs

    n = 1
    
    LDT_forcingClimo_struc%startFlag = .true. 

    
    call ESMF_ConfigGetAttribute(LDT_config,stime, &
         label="Forcing climatology temporal frequency of data:",&
         rc=rc)
    call LDT_verify(rc,&
         'Forcing climatology temporal frequency of data: not defined')
    call LDT_parseTimeString(stime, nsecs)

    LDT_forcingClimo_struc%nDataHours = nsecs/3600
    ! In number of hours in a day
    LDT_forcingClimo_struc%ntimes = nint(24*3600.0/nsecs)


  end subroutine forcingClimo_init


  subroutine forcingClimo_diagnose(n, pass)

    integer :: n
    integer :: pass

    integer :: rc
    integer                    :: i,t,l,m,fobjcount
    character*100, allocatable :: forcobjs(:)
    real,          pointer     :: forcdata1(:)
    type(ESMF_Field)           :: f1Field
    integer                    :: varCount,f1flag

    call ESMF_StateGet(LDT_FORC_Base_State(n,1),itemCount=fobjcount,rc=rc)
    call LDT_verify(rc,'ESMF_StateGet failed for objcount in forcingClimo_diagnose')
    
    allocate(forcobjs(fobjcount))

    call ESMF_StateGet(LDT_FORC_Base_State(n,1),itemNameList=forcobjs,rc=rc)
    call LDT_verify(rc,'ESMF_StateGet failed for forcobjs1 in forcingClimo_diagnose')
    
    varCount = 0

    l = LDT_rc%doy

    if((mod(LDT_rc%yr,4) .eq. 0 .and. mod(LDT_rc%yr, 100).ne.0) & !leap year
         .or.(mod(LDT_rc%yr,400) .eq.0)) then 
       ! decrease the doy by a day after Feb 28:
       if(LDT_rc%doy.gt.59) then 
          l = l - 1
       endif      
    endif

    m = LDT_rc%hr/LDT_forcingClimo_struc%nDataHours + 1

    LDT_forcingClimo_struc%nvars = fobjCount

    if(LDT_forcingClimo_struc%startFlag) then 
       allocate(LDT_forcingClimo_struc%dataEntry(LDT_forcingClimo_struc%nvars))
       allocate(LDT_forcingClimo_struc%varName(LDT_forcingClimo_struc%nvars))
       allocate(LDT_forcingClimo_struc%units(LDT_forcingClimo_struc%nvars))
       allocate(LDT_forcingClimo_struc%shortName(LDT_forcingClimo_struc%nvars))
       allocate(LDT_forcingClimo_struc%varId(LDT_forcingClimo_struc%nvars))
       
       do i=1,LDT_forcingClimo_struc%nvars
          allocate(LDT_forcingClimo_struc%dataEntry(i)%sx_mu(LDT_rc%ntiles(n),&
               365,LDT_forcingClimo_struc%ntimes))
          allocate(LDT_forcingClimo_struc%dataEntry(i)%sx_count(LDT_rc%ntiles(n),&
               365,LDT_forcingClimo_struc%ntimes))
          LDT_forcingClimo_struc%dataEntry(i)%sx_mu = 0
          LDT_forcingClimo_struc%dataEntry(i)%sx_count = 0
       end do
    endif

    do i=1,fobjcount

       if(LDT_forcingClimo_struc%startFlag) then 
          LDT_forcingClimo_struc%varName(i) = forcobjs(i)
          
          call mapShortVarname(forcobjs(i),LDT_forcingClimo_struc%shortName(i))

          LDT_forcingClimo_struc%units(i) = "TBD" !for now

       endif

       call ESMF_StateGet(LDT_FORC_Base_State(n,1),forcobjs(i),f1Field,&
            rc=rc)
       
       call ESMF_AttributeGet(f1Field, "Enabled", f1flag,rc=rc)

       if(f1flag.eq.1) then 
          varCount = varCount + 1

          call ESMF_FieldGet(f1Field,localDE=0,farrayPtr=forcdata1, &
               rc=rc)
          call LDT_verify(rc,'ESMF_FieldGet (ffield) failed in forcingClimo_diagnose')
          
          if(pass.eq.1) then 
             !read forcing datasets, save into data structures to compute PDFs.
             do t=1,LDT_rc%ntiles(n)
                if(forcdata1(t).ne.-9999.0) then 
                   
                   LDT_forcingClimo_struc%dataEntry(varCount)%sx_mu(t,l,m) = &
                        LDT_forcingClimo_struc%dataEntry(varCount)%sx_mu(t,l,m) &
                        + forcdata1(t)

                   LDT_forcingClimo_struc%dataEntry(varCount)%sx_count(t,l,m) = &
                        LDT_forcingClimo_struc%dataEntry(varCount)%sx_count(t,l,m) + 1
                endif
             enddo

          else

          end if
       endif
    enddo
    LDT_forcingClimo_struc%startFlag = .false.

  end subroutine forcingClimo_diagnose


  subroutine forcingClimo_compute(pass)

    integer :: n
    integer :: pass

    integer :: rc
    integer :: i,t,l,m

    n = 1

    do i=1,LDT_forcingClimo_struc%nvars
       do t=1,LDT_rc%ntiles(n)
          do l=1,365
             do m=1,LDT_forcingClimo_struc%ntimes
                if(LDT_forcingClimo_struc%dataEntry(i)%sx_count(t,l,m).gt.0) then 
                   LDT_forcingClimo_struc%dataEntry(i)%sx_mu(t,l,m) = &
                        LDT_forcingClimo_struc%dataEntry(i)%sx_mu(t,l,m)/ & 
                        LDT_forcingClimo_struc%dataEntry(i)%sx_count(t,l,m)

                else
                   LDT_forcingClimo_struc%dataEntry(i)%sx_mu(t,l,m) = LDT_rc%udef
                endif
             enddo
          enddo
       enddo
    enddo

  end subroutine forcingClimo_compute


  subroutine forcingClimo_output(pass)

    integer :: n
    integer :: pass

    integer           :: l
    character*3       :: cdoy
    integer           :: ftn
    integer           :: iret
    character*100     :: outfile

    n = 1
    if(LDT_masterproc) then
       call LDT_create_output_directory('FORCING')

       do l=1,365
          write(unit=cdoy, fmt='(i3.3)') l
          outfile  = trim(LDT_rc%odir)//'/FORCING/LDT_FORC_CLIMO_'//&
               trim(cdoy)//'.nc'

#if (defined USE_NETCDF3 || defined USE_NETCDF4)
#if (defined USE_NETCDF4)
          iret = nf90_create(path=outfile,cmode=nf90_netcdf4,&
               ncid = ftn)
          call LDT_verify(iret,'creating netcdf file failed')
#endif
#if (defined USE_NETCDF3)
          iret = nf90_create(path=outfile,cmode=nf90_clobber,&
               ncid = ftn)
          call LDT_verify(iret,'creating netcdf file failed')
#endif
#endif
         ! Write out meteorological forcing output:
          write(LDT_logunit,*)"[INFO] Writing output file: ",trim(outfile) 
          call writeNetcdfClimoOutput(n,l,ftn)
       
#if (defined USE_NETCDF3 || defined USE_NETCDF4)
          if(LDT_masterproc) then 
             iret = nf90_close(ftn)
          endif
#endif
       enddo
    endif
     
  end subroutine forcingClimo_output

!BOP
! !ROUTINE: writeNetcdfClimoOutput
! \label{writeNetcdfClimoOutput}
! 
! !INTERFACE: writeNetcdfClimoOutput
  subroutine writeNetcdfClimoOutput(n, l, ftn)
! !USES: 

! !ARGUMENTS: 
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: l
    integer,   intent(in)   :: ftn
! 
! !DESCRIPTION: 
!  This routine writes an output file in the NETCDF format based on the 
!  list of selected output variables. 
!  The arguments are: 
!  \begin{description}
!    \item[n] index of the nest
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[outInterval]   history output frequency
!    \item[nsoillayers]  Number of soil layers
!    \item[lyrthk]   Thickness of soil layers
!    \item[model\_name] Name of the model that generates the output
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[defineNETCDFheadervar](\ref{defineNetCDFClimoHeaderVar})
!     writes the required headers for a single variable
!   \item[writeSingleNETCDFvar](\ref{writeSingleNETCDFvar})
!     writes a single variable into a netcdf formatted file. 
!   \item[LDT\_verify](\ref{LDT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP
    integer               :: dimID(3)
    integer               :: xlatId, xlonId
    integer               :: tdimID
    integer               :: t,c,r,index1,i,m
    type(LDT_forcdataEntry), pointer :: xlat, xlong
    character*8           :: xtime_begin_date
    character*6           :: xtime_begin_time
    character*50          :: xtime_units
    integer               :: iret
! Note that the fix to add lat/lon to the NETCDF output will output
! undefined values for the water points. 
    character(len=8)      :: date
    character(len=10)     :: time
    character(len=5)      :: zone
    integer, dimension(8) :: values
    type(LDT_forcdataEntry), pointer :: dataEntry

    character(100)        :: zterp_flag
! __________________________________________________________

#if (defined USE_NETCDF3 || defined USE_NETCDF4)           
    call date_and_time(date,time,zone,values)
    
    allocate(xlat)
    allocate(xlong)

    xlat%short_name = "lat"
    xlat%standard_name = "latitude"
    xlat%units = "degree_north"
    xlat%nunits = 1
    xlat%vlevels = 1
    xlat%timeAvgOpt = 0 
    xlat%selectOpt = 1
    allocate(xlat%modelOutput(1,LDT_rc%ntiles(n),xlat%vlevels))
    allocate(xlat%count(1,xlat%vlevels))
    xlat%count = 1
    allocate(xlat%unittypes(1))
    xlat%unittypes(1) = "degree_north"
    xlat%valid_min = 0.0
    xlat%valid_max = 0.0

    xlong%short_name = "lon"
    xlong%standard_name = "longitude"
    xlong%units = "degree_east"
    xlong%nunits = 1
    xlong%vlevels = 1
    xlong%timeAvgOpt = 0 
    xlong%selectOpt = 1
    allocate(xlong%modelOutput(1,LDT_rc%ntiles(n),xlong%vlevels))
    allocate(xlong%count(1,xlong%vlevels))
    xlong%count = 1
    allocate(xlong%unittypes(1))
    xlong%unittypes(1) = "degree_east"
    xlong%valid_min = 0.0
    xlong%valid_max = 0.0

    if(LDT_masterproc) then 
       if(LDT_rc%wopt.eq."1d tilespace") then 
          call LDT_verify(nf90_def_dim(ftn,'ntiles',LDT_rc%glbntiles(n),&
               dimID(1)),&
               'nf90_def_dim for ntiles failed in LDT_metforcingMod')
          call LDT_verify(nf90_def_dim(ftn,'ntimes',LDT_forcingClimo_struc%ntimes,&
               dimID(2)),&
               'nf90_def_dim for ntimes failed in LDT_metforcingMod')
       elseif(LDT_rc%wopt.eq."2d gridspace") then 
          call LDT_verify(nf90_def_dim(ftn,'east_west',LDT_rc%gnc(n),&
               dimID(1)),&
               'nf90_def_dim for east_west failed in LDT_metforcingMod')
          call LDT_verify(nf90_def_dim(ftn,'north_south',LDT_rc%gnr(n),&
               dimID(2)),&
               'nf90_def_dim for north_south failed in LDT_metforcingMod')
          call LDT_verify(nf90_def_dim(ftn,'ntimes',LDT_forcingClimo_struc%ntimes,&
               dimID(3)),&
               'nf90_def_dim for ntimes failed in LDT_metforcingMod')
       endif

!LDT output is always writing output for a single time
       call LDT_verify(nf90_def_dim(ftn,'time',1,tdimID),&
            'nf90_def_dim for time failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", LDT_rc%udef),&
            'nf90_put_att for missing_value failed in LDT_metforcingMod')
       
       call defineNetCDFClimoHeaderVar(n,ftn,dimID, "lat",&
            "latitude", xlatid,"-",non_model_field=.true.)
       call defineNetCDFClimoHeaderVar(n,ftn,dimID, "lon",&
            "longitude", xlonid,"-",non_model_field=.true.)
       
       do i=1,LDT_forcingClimo_struc%nvars

          call defineNetCDFClimoHeaderVar(n,ftn,dimId,&
               LDT_forcingClimo_struc%shortName(i),&
               LDT_forcingClimo_struc%varName(i),&
               LDT_forcingClimo_struc%varID(i),&
               LDT_forcingClimo_struc%units(i))
       enddo

       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
            "LDT metforcing output"),&
            'nf90_put_att for title failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
            trim(LDT_rc%institution)),&
            'nf90_put_att for institution failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
            "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
            date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
            'nf90_put_att for history failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
            "Arsenault_etal_GMD_2018, Kumar_etal_EMS_2006"),&
            'nf90_put_att for references failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
            "CF-1.6"),'nf90_put_att for conventions failed in LDT_metforcingMod')
       call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
            "website: http://lis.gsfc.nasa.gov/"),&
            'nf90_put_att for comment failed in LDT_metforcingMod')

     ! -- Grid information --
       select case ( LDT_rc%lis_map_proj(n) )

        case( "latlon" )
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "EQUIDISTANT CYLINDRICAL"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,9)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,10)))       

        case( "mercator" )
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "MERCATOR"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

        case( "lambert" )   ! Lambert conformal
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "LAMBERT CONFORMAL"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

        case( "polar" )    ! polar stereographic
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
               "POLAR STEREOGRAPHIC"))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LAT", &
               LDT_rc%gridDesc(n,4)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"SOUTH_WEST_CORNER_LON", &
               LDT_rc%gridDesc(n,5)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
               LDT_rc%gridDesc(n,10)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
               LDT_rc%gridDesc(n,7)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
               LDT_rc%gridDesc(n,11)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
               LDT_rc%gridDesc(n,8)))
          call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
               LDT_rc%gridDesc(n,9)))

       case default 
       end select
       call LDT_verify(nf90_enddef(ftn))
    endif

    do t=1,LDT_rc%ntiles(n)
       c = LDT_domain(n)%tile(t)%col
       r = LDT_domain(n)%tile(t)%row
       index1 = LDT_domain(n)%gindex(c,r)
       xlat%modelOutput(1,t,1) = LDT_domain(n)%grid(index1)%lat
       xlong%modelOutput(1,t,1) = LDT_domain(n)%grid(index1)%lon
    enddo
    call writeSingleNETCDFnonmodelfield(ftn,n,xlatid,xlat%modelOutput(1,:,1))
    call writeSingleNETCDFnonmodelfield(ftn,n,xlonid,xlong%modelOutput(1,:,1))

    do i=1,LDT_forcingClimo_struc%nvars
       call writeSingleNetCDFforcClimovar(ftn,n,&
            LDT_forcingClimo_struc%varID(i),&
            LDT_forcingClimo_struc%dataEntry(i)%sx_mu(:,l,:))
    enddo

    deallocate(xlat%modelOutput)
    deallocate(xlat%count)
    deallocate(xlat%unittypes)
    deallocate(xlat)

    deallocate(xlong%modelOutput)
    deallocate(xlong%count)
    deallocate(xlong%unittypes)
    deallocate(xlong)
#endif

  end subroutine writeNetcdfClimoOutput


!BOP
! !ROUTINE: defineNetCDFClimoHeaderVar
! \label{defineNetCDFClimoHeaderVar}
! 
! !INTERFACE: 
  subroutine defineNetCDFClimoHeaderVar(n,ftn,dimID, shortName, &
       standardName, varID,units,non_model_field)
! !USES: 

! !ARGUMENTS:     
    integer               :: n
    integer               :: ftn
    character(len=*)      :: shortName
    character(len=*)      :: standardName
    character(len=*)      :: units
    integer               :: varId
    integer               :: dimID(3)
    logical, optional     :: non_model_field
! 
! !DESCRIPTION: 
!    This routine writes the required NETCDF header for a single variable
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
!   \item[LDT\_endrun](\ref{LDT_endrun})
!     call to abort program when a fatal error is detected. 
!   \item[LDT\_verify](\ref{LDT_verify})
!     call to check if the return value is valid or not.
!   \end{description}
!EOP

    integer       :: shuffle, deflate, deflate_level
    integer       :: fill_value
    logical       :: n_model

    if(present(non_model_field)) then 
       n_model = non_model_field 
    else
       n_model = .false. 
    endif

#if(defined USE_NETCDF3 || defined USE_NETCDF4)

    shuffle = NETCDF_shuffle
    deflate = NETCDF_deflate
    deflate_level =NETCDF_deflate_level

    if(n_model) then 
       if(LDT_rc%wopt.eq."1d tilespace") then              
          call LDT_verify(nf90_def_var(ftn,trim(shortName),&
               nf90_float,&
               dimids = dimID(1:1), varID=varId),&
               'nf90_def_var for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')
          call LDT_verify(nf90_def_var_fill(ftn,&
               varID, 1,fill_value), 'nf90_def_var_fill failed for '//&
               shortName)
#if(defined USE_NETCDF4)                
          call LDT_verify(nf90_def_var_deflate(ftn,&
               varId,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')                     
#endif                
       elseif(LDT_rc%wopt.eq."2d gridspace") then 
          call LDT_verify(nf90_def_var(ftn,trim(shortName),&
               nf90_float,&
               dimids = dimID(1:2), varID=varID),&
               'nf90_def_var for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')                     
          call LDT_verify(nf90_def_var_fill(ftn,&
               varID, &
               1,fill_value), 'nf90_def_var_fill failed for '//&
               shortName)
#if(defined USE_NETCDF4)
          call LDT_verify(nf90_def_var_deflate(ftn,&
               varId,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')  
#endif
       endif
    else
       if(LDT_rc%wopt.eq."1d tilespace") then              
          call LDT_verify(nf90_def_var(ftn,trim(shortName),&
               nf90_float,&
               dimids = dimID(1:2), varID=varId),&
               'nf90_def_var for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')
          call LDT_verify(nf90_def_var_fill(ftn,&
               varID, 1,fill_value), 'nf90_def_var_fill failed for '//&
               shortName)
#if(defined USE_NETCDF4)                
          call LDT_verify(nf90_def_var_deflate(ftn,&
               varId,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')                     
#endif                
       elseif(LDT_rc%wopt.eq."2d gridspace") then 
          call LDT_verify(nf90_def_var(ftn,trim(shortName),&
               nf90_float,&
               dimids = dimID(1:3), varID=varID),&
               'nf90_def_var for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')                     
          call LDT_verify(nf90_def_var_fill(ftn,&
               varID, &
               1,fill_value), 'nf90_def_var_fill failed for '//&
               shortName)
#if(defined USE_NETCDF4)
          call LDT_verify(nf90_def_var_deflate(ftn,&
               varId,&
               shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate for '//trim(shortName)//&
               'failed in defineNETCDFheadervar')                     
#endif                
       endif
    endif
    call LDT_verify(nf90_put_att(ftn,varID,&
         "units",trim(units)),&
         'nf90_put_att for units failed in defineNetCDFClimoHeaderVar')
    call LDT_verify(nf90_put_att(ftn,varID,&
         "standard_name",trim(standardname)),&
         'nf90_put_att for standard_name failed in defineNetCDFClimoHeaderVar')
    call LDT_verify(nf90_put_att(ftn,varID,&
         "scale_factor",1.0),&
         'nf90_put_att for scale_factor failed in defineNetCDFClimoHeaderVar')
    call LDT_verify(nf90_put_att(ftn,varID,&
         "add_offset",0.0),&
         'nf90_put_att for add_offset failed in defineNetCDFClimoHeaderVar')
    call LDT_verify(nf90_put_att(ftn,varID,&
         "missing_value",LDT_rc%udef),&
         'nf90_put_att for missing_value failed in defineNetCDFClimoHeaderVar')
    call LDT_verify(nf90_put_att(ftn,varID,&
         "_FillValue",LDT_rc%udef),&
         'nf90_put_att for _FillValue failed in defineNetCDFClimoHeaderVar')
#endif

  end subroutine defineNetCDFClimoHeaderVar


!BOP
! !ROUTINE: writeSingleNETCDFforcClimoVar
! \label{writeSingleNETCDFforcClimoVar}
!
! !INTERFACE: 
  subroutine writeSingleNETCDFforcClimoVar(ftn, n, varid, varData)
! !USES: 
    use LDT_coreMod,    only : LDT_rc
    use LDT_historyMod, only : LDT_writevar_netcdf

    implicit none

    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: varid
    real                    :: varData(LDT_rc%ntiles(n),LDT_forcingClimo_struc%ntimes)


! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LDT\_writevar\_netcdf](\ref{LDT_writevar_netcdf})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!EOP    
    integer           :: m 

    do m=1,LDT_forcingClimo_struc%ntimes
       call LDT_writevar_netcdf(ftn, n,&
            varData(:,m),&
            varId, &
            dim1=m)
    enddo
  end subroutine writeSingleNETCDFforcClimoVar


!BOP
! !ROUTINE: writeSingleNETCDFnonmodelField
! \label{writeSingleNETCDFnonmodelField}
!
! !INTERFACE: 
  subroutine writeSingleNETCDFnonmodelfield(ftn, n, varid, varData)
! !USES: 
    use LDT_coreMod,    only : LDT_rc
    use LDT_historyMod, only : LDT_writevar_netcdf

    implicit none

    integer,   intent(in)   :: ftn
    integer,   intent(in)   :: n 
    integer,   intent(in)   :: varid
    real                    :: varData(LDT_rc%ntiles(n))


! 
! !DESCRIPTION: 
!  This routine writes a single variable to a NETCDF file
!  The arguments are: 
!  \begin{description}
!    \item[ftn] file unit for the output file
!    \item[ftn\_stats] file unit for the output statistics file
!    \item[n] index of the nest
!   \item[dataEntry]
!    object containing the values and attributes of the variable to be 
!    written
!  \end{description}
!
!   The routines invoked are: 
!   \begin{description}
!   \item[LDT\_writevar\_netcdf](\ref{LDT_writevar_netcdf})
!     writes a variable into a netcdf formatted file. 
!   \end{description}
!EOP    
    
    call LDT_writevar_netcdf(ftn, n,&
         varData,&
         varId, &
         dim1=1)
  end subroutine writeSingleNETCDFnonmodelfield

  subroutine mapShortVarname(fullname,shortName)

    character(len=*)   :: fullName
    character(len=*)   :: shortName

    if(fullName.eq."Eastward Wind Level 001") then 
       shortName = "Uwind"
    elseif(fullName.eq."Incident Longwave Radiation Level 001") then 
       shortName = "LWdown"
    elseif(fullName.eq."Incident Shortwave Radiation Level 001") then 
       shortName = "SWdown"
    elseif(fullName.eq."Near Surface Air Temperature Level 001") then 
       shortName = "Tair"
    elseif(fullName.eq."Near Surface Specific Humidity Level 001") then 
       shortName = "Qair"
    elseif(fullName.eq."Northward Wind Level 001") then 
       shortName = "Vwind"
    elseif(fullName.eq."Rainfall Rate Level 001") then 
       shortName = "Rainf"
    elseif(fullName.eq."Surface Pressure Level 001") then 
       shortName = "Psurf"
    elseif(fullName.eq."Convective Rainfall Rate Level 001") then 
       shortName = "CRainf"
    endif

  end subroutine mapShortVarname
end module forcingClimoMod
