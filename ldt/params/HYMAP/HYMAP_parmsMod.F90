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
module HYMAP_parmsMod
!BOP
!
! !MODULE: HYMAP_parmsMod
!
! !DESCRIPTION:
!  The code in this file implements routines to process HYMAP 
!  parameter data. 
!
! !REVISION HISTORY:
!
!  29 Mar 2013: Sujay Kumar; Initial implementation
!   2 Dec 2015: Augusto Getirana: Included drainage area and basin maps
!   1 Nov 2017: Augusto Getirana: Included flow type maps, baseflow and surface runoff dwi maps
!   9 Jun 2020: Yeosang Yoon: Support flexible grid setting (dx~=dy)
!  24 Aug 2021: Hiroko Beaudoing: Fix boundary for global domain
!
  use ESMF
  use LDT_coreMod
  use LDT_historyMod
  use LDT_paramDataMod
  use LDT_logMod
  use LDT_constantsMod, only : LDT_CONST_PATH_LEN

  implicit none

  PRIVATE
!------------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!------------------------------------------------------------------------------
  public :: HYMAPparms_init    !allocates memory for required structures
  public :: HYMAPparms_writeHeader
  public :: HYMAPparms_writeData

!------------------------------------------------------------------------------
! !PUBLIC TYPES:
!------------------------------------------------------------------------------
  public :: HYMAP_struc

  type, public :: hymap_type_dec
     real          :: hymapparms_gridDesc(20)
     character*50  :: hymap_proj
     character*50  :: hymap_gridtransform

     integer       :: bfdwicode, rundwicode, rivflocode

     character(len=LDT_CONST_PATH_LEN) :: riverwidthfile
     character(len=LDT_CONST_PATH_LEN) :: riverheightfile
     character(len=LDT_CONST_PATH_LEN) :: riverlengthfile
     character(len=LDT_CONST_PATH_LEN) :: riverzfile
     character(len=LDT_CONST_PATH_LEN) :: fldheightfile
     character(len=LDT_CONST_PATH_LEN) :: fldzfile
     character(len=LDT_CONST_PATH_LEN) :: flowdirxfile
     character(len=LDT_CONST_PATH_LEN) :: flowdiryfile
     character(len=LDT_CONST_PATH_LEN) :: gridelevfile
     character(len=LDT_CONST_PATH_LEN) :: griddistfile
     character(len=LDT_CONST_PATH_LEN) :: gridareafile
     character(len=LDT_CONST_PATH_LEN) :: drainareafile
     character(len=LDT_CONST_PATH_LEN) :: basinfile
     character(len=LDT_CONST_PATH_LEN) :: runoffdelayfile
     character(len=LDT_CONST_PATH_LEN) :: runoffdelaymfile
     character(len=LDT_CONST_PATH_LEN) :: baseflowdelayfile
     character(len=LDT_CONST_PATH_LEN) :: refqfile
     character(len=LDT_CONST_PATH_LEN) :: basinmaskfile
     character(len=LDT_CONST_PATH_LEN) :: flowtypefile 
     character(len=LDT_CONST_PATH_LEN) :: baseflowdwiratiofile 
     character(len=LDT_CONST_PATH_LEN) :: runoffdwiratiofile 

     type(LDT_paramEntry) :: hymap_river_width
     type(LDT_paramEntry) :: hymap_river_height
     type(LDT_paramEntry) :: hymap_river_length
     type(LDT_paramEntry) :: hymap_river_z
     type(LDT_paramEntry) :: hymap_fld_z
     type(LDT_paramEntry) :: hymap_fld_height
     type(LDT_paramEntry) :: hymap_flow_dir_x
     type(LDT_paramEntry) :: hymap_flow_dir_y
     type(LDT_paramEntry) :: hymap_grid_elev
     type(LDT_paramEntry) :: hymap_grid_dist
     type(LDT_paramEntry) :: hymap_grid_area
     type(LDT_paramEntry) :: hymap_drain_area
     type(LDT_paramEntry) :: hymap_basin
     type(LDT_paramEntry) :: hymap_runoff_delay
     type(LDT_paramEntry) :: hymap_runoff_delay_m
     type(LDT_paramEntry) :: hymap_baseflow_delay
     type(LDT_paramEntry) :: hymap_mask
     type(LDT_paramEntry) :: hymap_flow_type
     type(LDT_paramEntry) :: hymap_baseflow_dwi_ratio
     type(LDT_paramEntry) :: hymap_runoff_dwi_ratio
  end type hymap_type_dec

  type(hymap_type_dec), allocatable :: HYMAP_struc(:)

contains

!BOP
! 
! !ROUTINE: HYMAPparms_init
! \label{HYMAPparms_init}
! 
! !INTERFACE:
  subroutine HYMAPparms_init
! !USES:
    use LDT_fileIOMod, only : LDT_readDomainConfigSpecs
    use LDT_logMod,    only : LDT_verify
    use LDT_paramOptCheckMod, only: LDT_noahparmsOptChecks, &
                       LDT_gridOptChecks
!
! !DESCRIPTION:
!
! Allocates memory for data structures for reading 
! the HYMAP datasets 
! 
!  The routines invoked are: 
!  \begin{description}
!   \item[HYMAPparmssetup](\ref{HYMAPparmssetup}) \newline
!    calls the registry to invoke the HYMAPparms setup methods. 
!  \end{description}
!
!EOP
   implicit none
   integer       :: n
   integer       :: rc
   integer       :: bfdwicode, rundwicode, rivflocode
   integer, allocatable :: nextx(:,:)
   integer, allocatable :: nexty(:,:)
   integer, allocatable :: mask(:,:)
   logical       :: hymap_params_selected
   character*50  :: hymap_proj
   real,    allocatable :: hymapparms_gridDesc(:,:)

   hymap_params_selected = .false. 

   allocate(hymapparms_gridDesc(LDT_rc%nnest,20))
   allocate(HYMAP_struc(LDT_rc%nnest))

   do n=1,LDT_rc%nnest

      call set_param_attribs(HYMAP_struc(n)%hymap_river_width,&
           "HYMAP_river_width")
      call set_param_attribs(HYMAP_struc(n)%hymap_river_height,&
           "HYMAP_river_height")
      call set_param_attribs(HYMAP_struc(n)%hymap_river_z,&
           "HYMAP_river_roughness")
      call set_param_attribs(HYMAP_struc(n)%hymap_fld_z,&
           "HYMAP_floodplain_roughness")
      call set_param_attribs(HYMAP_struc(n)%hymap_river_length,&
           "HYMAP_river_length")
      call set_param_attribs(HYMAP_struc(n)%hymap_fld_height,&
           "HYMAP_floodplain_height")
      call set_param_attribs(HYMAP_struc(n)%hymap_flow_dir_x,&
           "HYMAP_flow_direction_x")
      call set_param_attribs(HYMAP_struc(n)%hymap_flow_dir_y,&
           "HYMAP_flow_direction_y")
      call set_param_attribs(HYMAP_struc(n)%hymap_grid_elev,&
           "HYMAP_grid_elevation")
      call set_param_attribs(HYMAP_struc(n)%hymap_grid_dist,&
           "HYMAP_grid_distance")
      call set_param_attribs(HYMAP_struc(n)%hymap_grid_area,&
           "HYMAP_grid_area")
      call set_param_attribs(HYMAP_struc(n)%hymap_drain_area,&
           "HYMAP_drain_area")
      call set_param_attribs(HYMAP_struc(n)%hymap_basin,&
           "HYMAP_basin")
      call set_param_attribs(HYMAP_struc(n)%hymap_runoff_delay,&
           "HYMAP_runoff_time_delay")
      call set_param_attribs(HYMAP_struc(n)%hymap_runoff_delay_m,&
           "HYMAP_runoff_time_delay_multiplier")
      call set_param_attribs(HYMAP_struc(n)%hymap_baseflow_delay,&
           "HYMAP_baseflow_time_delay")
      call set_param_attribs(HYMAP_struc(n)%hymap_mask,&
           "HYMAP_basin_mask")
      call set_param_attribs(HYMAP_struc(n)%hymap_flow_type,&
           "HYMAP_river_flow_type")
      call set_param_attribs(HYMAP_struc(n)%hymap_baseflow_dwi_ratio,&
           "HYMAP_baseflow_dwi_ratio")
      call set_param_attribs(HYMAP_struc(n)%hymap_runoff_dwi_ratio,&
           "HYMAP_runoff_dwi_ratio")
      
   end do

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP floodplain height levels:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%hymap_fld_height%num_bins,rc=rc)
      call LDT_verify(rc,'HYMAP floodplain height levels: not specified')
   enddo
 
   do n=1,LDT_rc%nnest

      allocate(HYMAP_struc(n)%hymap_river_width%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_river_width%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_river_height%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_river_height%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_river_length%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_river_length%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_river_z%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_river_z%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_fld_z%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_fld_z%vlevels))       
      
      HYMAP_struc(n)%hymap_fld_height%vlevels = &
           HYMAP_struc(n)%hymap_fld_height%num_bins
      allocate(HYMAP_struc(n)%hymap_fld_height%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_fld_height%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_flow_dir_x%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_flow_dir_x%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_flow_dir_y%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_flow_dir_y%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_grid_elev%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_grid_elev%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_grid_dist%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_grid_dist%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_grid_area%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_grid_area%vlevels))       

      allocate(HYMAP_struc(n)%hymap_drain_area%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_drain_area%vlevels))       

      allocate(HYMAP_struc(n)%hymap_basin%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_basin%vlevels))       

      allocate(HYMAP_struc(n)%hymap_runoff_delay%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_runoff_delay%vlevels))       
      allocate(HYMAP_struc(n)%hymap_runoff_delay_m%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_runoff_delay_m%vlevels))       
      
      allocate(HYMAP_struc(n)%hymap_baseflow_delay%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_baseflow_delay%vlevels))       
!      allocate(HYMAP_struc(n)%hymap_ref_q%value(&
!           LDT_rc%lnc(n),LDT_rc%lnr(n),&
!           HYMAP_struc(n)%hymap_ref_q%vlevels))       
      allocate(HYMAP_struc(n)%hymap_mask%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_mask%vlevels))                        

      allocate(HYMAP_struc(n)%hymap_flow_type%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_flow_type%vlevels))                        

      allocate(HYMAP_struc(n)%hymap_baseflow_dwi_ratio%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_baseflow_dwi_ratio%vlevels))                        

      allocate(HYMAP_struc(n)%hymap_runoff_dwi_ratio%value(&
           LDT_rc%lnc(n),LDT_rc%lnr(n),&
           HYMAP_struc(n)%hymap_runoff_dwi_ratio%vlevels))                        
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP river width map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%riverwidthfile,rc=rc)
      call LDT_verify(rc,'HYMAP river width map: not specified')
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP river height map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%riverheightfile,rc=rc)
      call LDT_verify(rc,'HYMAP river height map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP river length map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%riverlengthfile,rc=rc)
      call LDT_verify(rc,'HYMAP river length map: not specified')
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP river roughness map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%riverzfile,rc=rc)
      call LDT_verify(rc,'HYMAP river roughness map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP floodplain height map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%fldheightfile,rc=rc)
      call LDT_verify(rc,'HYMAP floodplain height map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP floodplain roughness map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%fldzfile,rc=rc)
      call LDT_verify(rc,'HYMAP floodplain roughness map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP flow direction x map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%flowdirxfile,rc=rc)
      call LDT_verify(rc,'HYMAP flow direction x map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP flow direction y map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%flowdiryfile,rc=rc)
      call LDT_verify(rc,'HYMAP flow direction y map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP grid elevation map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%gridelevfile,rc=rc)
      call LDT_verify(rc,'HYMAP grid elevation map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP grid distance map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%griddistfile,rc=rc)
      call LDT_verify(rc,'HYMAP grid distance map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP grid area map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%gridareafile,rc=rc)
      call LDT_verify(rc,'HYMAP grid area map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP drainage area map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%drainareafile,rc=rc)
      call LDT_verify(rc,'HYMAP drainage area map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP basin map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%basinfile,rc=rc)
      call LDT_verify(rc,'HYMAP basin map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP runoff time delay map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%runoffdelayfile,rc=rc)
      call LDT_verify(rc,'HYMAP runoff time delay map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP runoff time delay multiplier map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%runoffdelaymfile,rc=rc)
      call LDT_verify(rc,'HYMAP runoff time delay multiplier map: not specified')
   enddo
   call ESMF_ConfigFindLabel(LDT_config,"HYMAP baseflow time delay map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%baseflowdelayfile,rc=rc)
      call LDT_verify(rc,'HYMAP baseflow time delay map: not specified')
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP basin mask map:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%basinmaskfile,rc=rc)
      call LDT_verify(rc,'HYMAP basin mask map: not specified')
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP river flow type map:",rc=rivflocode)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%flowtypefile,rc=rc)
      if (rivflocode.eq.0) then
         call LDT_verify(rc,'HYMAP river flow type map: not specified')
      endif
      HYMAP_struc(n)%rivflocode = rivflocode
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP baseflow dwi ratio map:",rc=bfdwicode)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%baseflowdwiratiofile,rc=rc) 
      if (bfdwicode.eq.0) then
         call LDT_verify(rc,'HYMAP baseflow dwi ratio map: not specified')
      endif
      HYMAP_struc(n)%bfdwicode = bfdwicode
   enddo

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP runoff dwi ratio map:",rc=rundwicode)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%runoffdwiratiofile,rc=rc) 
      if (rundwicode.eq.0) then
         call LDT_verify(rc,'HYMAP runoff dwi ratio map: not specified')
      endif
      HYMAP_struc(n)%rundwicode = rundwicode
   enddo

   write(LDT_logunit,*)" - - - - - - - - HYMAP Router Parameters - - - - - - - - - -"
   
   call ESMF_ConfigGetAttribute(LDT_config,hymap_proj,&
        label="HYMAP params map projection:",rc=rc)
   call LDT_verify(rc,'HYMAP params map projection: option not specified in the config file')
   
   HYMAP_struc(:)%hymap_proj = hymap_proj

   call ESMF_ConfigFindLabel(LDT_config,"HYMAP params spatial transform:",rc=rc)
   do n=1,LDT_rc%nnest
      call ESMF_ConfigGetAttribute(LDT_config,HYMAP_struc(n)%hymap_gridtransform,&
           rc=rc)
      call LDT_verify(rc,'HYMAP params spatial transform: option not specified in the config file')
   enddo
   
   call LDT_readDomainConfigSpecs("HYMAP params",hymap_proj,&
        hymapparms_gridDesc)
   
    do n=1,LDT_rc%nnest
       HYMAP_struc(n)%hymapparms_gridDesc(:) = hymapparms_gridDesc(n,:)

       call LDT_gridOptChecks( n, "HYMAP river width", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%riverwidthfile)
       call read_HYMAP_river_width(n, &
            HYMAP_struc(n)%hymap_river_width%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%riverwidthfile)
       
       call LDT_gridOptChecks( n, "HYMAP river height", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%riverheightfile)
       call read_HYMAP_river_height(n,&
            HYMAP_struc(n)%hymap_river_height%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%riverheightfile)
       
       call LDT_gridOptChecks( n, "HYMAP river length", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, & 
              hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%riverlengthfile)
       call read_HYMAP_river_length(n,&
            HYMAP_struc(n)%hymap_river_length%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%riverlengthfile)
       
       call LDT_gridOptChecks( n, "HYMAP river roughness", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%riverzfile)
       call read_HYMAP_river_z(n,&
            HYMAP_struc(n)%hymap_river_z%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%riverzfile)

       call LDT_gridOptChecks( n, "HYMAP floodplain roughness", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%fldzfile)
       call read_HYMAP_fld_z(n,&
            HYMAP_struc(n)%hymap_fld_z%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%fldzfile)
       
       call LDT_gridOptChecks( n, "HYMAP floodplain height", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%fldheightfile)
       call read_HYMAP_fld_height(n,&
            HYMAP_struc(n)%hymap_fld_height%value(:,:,:))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%fldheightfile)
       
       call LDT_gridOptChecks( n, "HYMAP flow direction x", &
              HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
              hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%flowdirxfile)
       call read_HYMAP_flow_dir_x(n,&
            HYMAP_struc(n)%hymap_flow_dir_x%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%flowdirxfile)

       call LDT_gridOptChecks( n, "HYMAP flow direction y", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%flowdiryfile)
       call read_HYMAP_flow_dir_y(n,&
            HYMAP_struc(n)%hymap_flow_dir_y%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%flowdiryfile)
       
       call LDT_gridOptChecks( n, "HYMAP grid elevation", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%gridelevfile)
       call read_HYMAP_grid_elev(n,&
            HYMAP_struc(n)%hymap_grid_elev%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%gridelevfile)
       
       call LDT_gridOptChecks( n, "HYMAP grid distance", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
              hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%griddistfile)
       call read_HYMAP_grid_dist(n,&
            HYMAP_struc(n)%hymap_grid_dist%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%griddistfile)

       call LDT_gridOptChecks( n, "HYMAP grid area", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )       
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%gridareafile)
       call read_HYMAP_grid_area(n,&
            HYMAP_struc(n)%hymap_grid_area%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%gridareafile)

       call LDT_gridOptChecks( n, "HYMAP drainage area", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )       
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%drainareafile)
       call read_HYMAP_drain_area(n,&
            HYMAP_struc(n)%hymap_drain_area%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%drainareafile)

       call LDT_gridOptChecks( n, "HYMAP basin", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )       
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%basinfile)
       call read_HYMAP_basin(n,&
            HYMAP_struc(n)%hymap_basin%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%basinfile)

       call LDT_gridOptChecks( n, "HYMAP runoff delay", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )       
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%runoffdelayfile)
       call read_HYMAP_runoff_delay(n,&
            HYMAP_struc(n)%hymap_runoff_delay%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%runoffdelayfile)

       call LDT_gridOptChecks( n, "HYMAP runoff delay multiplier", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%runoffdelaymfile)
       call read_HYMAP_runoff_delaym(n,&
            HYMAP_struc(n)%hymap_runoff_delay_m%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%runoffdelaymfile)
       
       call LDT_gridOptChecks( n, "HYMAP baseflow delay", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj,  &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%baseflowdelayfile)
       call read_HYMAP_baseflow_delay(&
            n,HYMAP_struc(n)%hymap_baseflow_delay%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%baseflowdelayfile)
       
       call LDT_gridOptChecks( n, "HYMAP basin mask", &
            HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
            hymapparms_gridDesc(n,9) )
       write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%basinmaskfile)
       call read_HYMAP_basin_mask(&
            n,HYMAP_struc(n)%hymap_mask%value(:,:,1))
       write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%basinmaskfile)
       
       if (HYMAP_struc(n)%rivflocode.eq.0) then
          call LDT_gridOptChecks( n, "HYMAP river flow type map", &
               HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
               hymapparms_gridDesc(n,9) )

          write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%flowtypefile)
          call read_HYMAP_flow_type(&
               n,HYMAP_struc(n)%hymap_flow_type%value(:,:,1))
          write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%flowtypefile)
       endif

       if (HYMAP_struc(n)%bfdwicode.eq.0) then
          call LDT_gridOptChecks( n, "HYMAP baseflow DWI ratio map", &
               HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
               hymapparms_gridDesc(n,9) )

          write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%baseflowdwiratiofile)
          call read_HYMAP_baseflow_dwi_ratio(&
               n,HYMAP_struc(n)%hymap_baseflow_dwi_ratio%value(:,:,1))
          write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%baseflowdwiratiofile)
       endif

       if (HYMAP_struc(n)%rundwicode.eq.0) then
          call LDT_gridOptChecks( n, "HYMAP runoff DWI ratio map", &
               HYMAP_struc(n)%hymap_gridtransform, hymap_proj, &
               hymapparms_gridDesc(n,9) )

          write(LDT_logunit,*) 'reading '//trim(HYMAP_struc(n)%runoffdwiratiofile)
          call read_HYMAP_runoff_dwi_ratio(&
               n,HYMAP_struc(n)%hymap_runoff_dwi_ratio%value(:,:,1))
          write(LDT_logunit,*) 'Done reading '//trim(HYMAP_struc(n)%runoffdwiratiofile)
       endif

       allocate(nextx(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       allocate(nexty(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       allocate(mask(LDT_rc%lnc(n),LDT_rc%lnr(n)))
       
       nextx = nint(HYMAP_struc(n)%hymap_flow_dir_x%value(:,:,1))
       nexty = nint(HYMAP_struc(n)%hymap_flow_dir_y%value(:,:,1))
       mask = nint(HYMAP_struc(n)%hymap_mask%value(:,:,1))

       ! Yeosang Yoon: support flexible grid setting
       call adjust_nextxy(&
            LDT_rc%gnc(n),&
            LDT_rc%gnr(n),&
            -9999, &
            nextx, &
            nexty, &
            mask, &
            LDT_rc%gridDesc(n,5),&   ! lon min of smaller domain
            LDT_rc%gridDesc(n,4),&   ! lat min of smaller domain
            hymapparms_gridDesc(n,5),&   ! lon min of larger domain
            hymapparms_gridDesc(n,4),&   ! lat min of larger domain
            hymapparms_gridDesc(n,9),&   ! dx spatial resolution
            hymapparms_gridDesc(n,10))   ! dy spatial resolution
 
       HYMAP_struc(n)%hymap_flow_dir_x%value(:,:,1) = real(nextx(:,:))
       HYMAP_struc(n)%hymap_flow_dir_y%value(:,:,1) = real(nexty(:,:))
       
       LDT_rc%routing_grid_count=count(nextx/=-9999.and.mask>0)

       deallocate(nextx)
       deallocate(nexty)
       deallocate(mask)
    enddo
    
  end subroutine HYMAPparms_init

  subroutine HYMAPparms_writeHeader(n,ftn,dimID,monthID)

#if(defined USE_NETCDF3 || defined USE_NETCDF4) 
    use netcdf
#endif

    integer   :: n 
    integer   :: ftn
    integer   :: dimID(3)
    integer   :: monthID
    integer   :: tdimID(3)

    logical   :: hymap_params_selected
    
    tdimID(1) = dimID(1)
    tdimID(2) = dimID(2)

#if(defined USE_NETCDF3 || defined USE_NETCDF4)     
    call LDT_verify(nf90_def_dim(ftn,'fld_height_levels',&
         HYMAP_struc(n)%hymap_fld_height%vlevels,tdimID(3)))
#endif

    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_river_width)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_river_height)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_river_length)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_river_z)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_fld_z)
    call LDT_writeNETCDFdataHeader(n,ftn,tdimID,&
         HYMAP_struc(n)%hymap_fld_height)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_flow_dir_x)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_flow_dir_y)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_grid_elev)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_grid_dist)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_grid_area)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_drain_area)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_basin)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_runoff_delay)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_runoff_delay_m)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_baseflow_delay)
!    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
!         HYMAP_struc(n)%hymap_ref_q)
    call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
         HYMAP_struc(n)%hymap_mask)
    if (HYMAP_struc(n)%rivflocode.eq.0) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            HYMAP_struc(n)%hymap_flow_type)
    endif
    if (HYMAP_struc(n)%bfdwicode.eq.0) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            HYMAP_struc(n)%hymap_baseflow_dwi_ratio)
    endif
    if (HYMAP_struc(n)%rundwicode.eq.0) then
       call LDT_writeNETCDFdataHeader(n,ftn,dimID,&
            HYMAP_struc(n)%hymap_runoff_dwi_ratio)
    endif
  end subroutine HYMAPparms_writeHeader
  
  subroutine HYMAPparms_writeData(n,ftn)

    integer   :: n 
    integer   :: ftn
    logical   :: hymap_params_selected
    

    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_river_width)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_river_height)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_river_length)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_river_z)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_fld_z)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_fld_height)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_flow_dir_x)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_flow_dir_y)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_grid_elev)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_grid_dist)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_grid_area)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_drain_area)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_basin)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_runoff_delay)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_runoff_delay_m)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_baseflow_delay)
!    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_ref_q)
    call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_mask)
    if (HYMAP_struc(n)%rivflocode.eq.0) then
       call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_flow_type)
    endif
    if (HYMAP_struc(n)%bfdwicode.eq.0) then
       call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_baseflow_dwi_ratio)
    endif
    if (HYMAP_struc(n)%rundwicode.eq.0) then
       call LDT_writeNETCDFdata(n,ftn,HYMAP_struc(n)%hymap_runoff_dwi_ratio)
    endif

  end subroutine HYMAPparms_writeData

  subroutine adjust_nextxy(nx,ny,imis,i2nextx,i2nexty,i2mask,zgx,zgy,zpx,zpy,xres,yres)
    ! ================================================
    ! to adjust flow direction matrixes for smaller domain
    ! by Augusto GETIRANA
    ! on 13th Mar 2012
    ! at HSL/GSFC/NASA
    ! ================================================   
    
    use LDT_logmod, only : LDT_logunit

    implicit none       
  
    integer, intent(in)    :: nx                  ! number of grids in horizontal
    integer, intent(in)    :: ny                  ! number of grids in vertical
    integer, intent(in)    :: imis                ! integer undefined value
    integer, intent(inout) :: i2nextx(nx,ny)      ! point downstream horizontal
    integer, intent(inout) :: i2nexty(nx,ny)      ! point downstream horizontal
    integer, intent(in)    :: i2mask(nx,ny)       ! mask limiting modeled region (0: out; >=1: in)
    real*4,  intent(in)    :: zgx                 ! lon min of smaller domain
    real*4,  intent(in)    :: zgy                 ! lat min of smaller domain
    real*4,  intent(in)    :: zpx                 ! lon min of larger domain
    real*4,  intent(in)    :: zpy                 ! lat min of larger domain
    real*4,  intent(in)    :: xres,yres           ! spatial resolution
    
    integer, parameter :: ibound = -9
    
    integer             ::  idx,idy

    idx=int((zgx-zpx)/xres)
    idy=int((zgy-zpy)/yres)
    
    where(i2nextx>0)i2nextx=i2nextx-idx
    where(i2nexty>0)i2nexty=i2nexty-idy

!Hiroko: do not insert boundary if global domain
!        this fix only works on single processor run
    if ( idx.eq.0 .and. idy.eq.0 ) then
     write(LDT_logunit,*) '[INFO] HYMAP parameter global'
     where(i2nextx<1.and.i2nextx/=imis.and.i2mask>0)
        i2nextx=ibound
        i2nexty=ibound
     endwhere
     where(i2nextx>nx.and.i2mask>0)
        i2nextx=ibound
        i2nexty=ibound
     endwhere
    else    ! local domain, insert boundary
     i2nexty(1,:)=imis
     i2nexty(nx,:)=imis
     i2nexty(:,1)=imis
     i2nexty(:,ny)=imis
    
     i2nextx(1,:)=imis
     i2nextx(nx,:)=imis
     i2nextx(:,1)=imis
     i2nextx(:,ny)=imis
 
     where(i2nextx<=1.and.i2nextx/=imis.and.i2mask>0)
        i2nextx=ibound
        i2nexty=ibound
     endwhere
    
     where(i2nextx>=nx.and.i2mask>0)
        i2nextx=ibound
        i2nexty=ibound
     endwhere
    endif    ! global
    
    where(i2nexty<=1.and.i2nexty/=imis.and.i2mask>0)
       i2nextx=ibound
       i2nexty=ibound
    endwhere
    
    where(i2nexty>=ny.and.i2mask>0)
       i2nextx=ibound
       i2nexty=ibound
    endwhere
         
    print*, ''
  end subroutine adjust_nextxy

!BOP
! !ROUTINE:  set_param_attribs
! \label{set_param_attribs}
!
! !INTERFACE:
  subroutine set_param_attribs(paramEntry, short_name)

! !DESCRIPTION:
!   This routine reads over the parameter attribute entries
!   in the param_attribs.txt file.
!
! !USES:
   type(LDT_paramEntry),intent(inout) :: paramEntry
   character(len=*),    intent(in)    :: short_name

! ____________________________________________________
    
   
   paramEntry%short_name = trim(short_name)
   paramEntry%vlevels = 1
   paramEntry%selectOpt = 1
   paramEntry%source = "HYMAP"
   paramEntry%units ="none"
   paramEntry%num_times = 1
   paramEntry%num_bins = 1
   paramEntry%standard_name = trim(short_name)

  end subroutine set_param_attribs

end module HYMAP_parmsMod
