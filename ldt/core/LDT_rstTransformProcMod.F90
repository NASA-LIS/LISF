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
module LDT_rstTransformProcMod
!BOP
!
! !MODULE: LDT_rstTransformProcMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage 
!   the processing of climatological restart files for land 
!   surface models and routing models. 
!
! !REVISION HISTORY: 
!  26 Jan 2016    Sujay Kumar  Initial Specification
! 
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_fileIOMod
  use LDT_timeMgrMod
  use LDT_paramDataMod
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_rstTransformProcInit

!-----------------------------------------------------------------------------
! !PUBLIC TYPES:
!-----------------------------------------------------------------------------
  public :: LDT_rstTransform_struc

  type, public :: ldtrsttransdec
     integer :: tmp
  end type ldtrsttransdec
  
  type(ldtrsttransdec)   :: LDT_rstTransform_struc

contains

!BOP
! !ROUTINE: LDT_rst
! label{LDT_rstTransformProcInit}
! 
! !INTERFACE: 
  subroutine LDT_rstTransformProcInit
! 
! !DESCRIPTION: 
! 
!  This routine performs initialization in the Restart transformation processing mode by 
!  reading the relevant config file entries. 
!  
!  Example entries in the config file are as follows: 
!
! \begin{verabtim}
! LIS restart source:        "LSM"
! Input restart filename:    LIS_RST_NOAHMP401_201907010000.d01.coarse.nc
! Output restart filename:   LIS_RST_NOAHMP401_201503312345.d01.fine.nc
!  \end{verbatim}
!
!EOP  
  
    integer                   :: n,i 
    integer                   :: status
    character*20              :: stime
    character*100             :: model_name      
    
    n = 1
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%rstsource,&
         label="LIS restart source:",&
         rc=status)
    call LDT_verify(status,'LIS restart source: not defined')
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%inputrst,&
         label="Input restart filename:",&
         rc=status)
    call LDT_verify(status,'Input restart filename: not defined')
    
    call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%outputrst,&
         label="Output restart filename:",&
         rc=status)
    call LDT_verify(status,'Output restart filename: not defined')
    
    if(LDT_rc%rstsource.eq."LSM") then 
       if(LDT_rc%lsm.eq."Noah.3.2") then 
          model_name = "NOAH32"
       elseif(LDT_rc%lsm.eq."Noah.3.3") then 
          model_name = "NOAH33"
       elseif(LDT_rc%lsm.eq."Noah.3.6") then 
          model_name = "NOAH36"
       elseif(LDT_rc%lsm.eq."Noah.3.9") then 
          model_name = "NOAH39"
       elseif(LDT_rc%lsm.eq."Noah.2.7.1") then 
          model_name = "NOAH271"
       elseif(LDT_rc%lsm.eq."Noah-MP.3.6") then 
          model_name = "NOAHMP36"
       elseif(LDT_rc%lsm.eq."Noah-MP.4.0.1") then 
          model_name = "NOAHMP401"
       elseif(LDT_rc%lsm.eq."CLSMF2.5") then 
          model_name = "CLSMF25"
       elseif(LDT_rc%lsm.eq."RUC.3.7") then 
          model_name = "RUC37"
       elseif(LDT_rc%lsm.eq."VIC.4.1.1") then 
          model_name = "VIC411"
       elseif(LDT_rc%lsm.eq."VIC.4.1.2") then 
          model_name = "VIC412"
       elseif(LDT_rc%lsm .eq. "JULES.5.0") then
          model_name = "JULES50"
       else
          write(LDT_logunit,*) "[INFO] Restart File transform - LSMs supported: "
          write(LDT_logunit,*) "  -- CLSMF2.5, Noah.3.2, Noah.3.3, Noah.3.6, Noah.3.9, "
          write(LDT_logunit,*) "  -- Noah-MP.3.6, Noah-MP.4.0.1, "
          write(LDT_logunit,*) "     Noah.2.7.1, RUC.3.7, VIC.4.1.1, VIC.4.1.2 "
          write(LDT_logunit,*) "[ERR] No other LSMs supported at this time ... stopping."
          call LDT_endrun() 
       endif
    else
       write(LDT_logunit,*) "[INFO] Restart File transform - Only support the following LSM: "
       write(LDT_logunit,*) "  -- CLSMF2.5, Noah.3.2, Noah.3.3, Noah.3.6, Noah.3.9, "
       write(LDT_logunit,*) "  -- Noah-MP.3.6, Noah-MP.4.0.1, "
       write(LDT_logunit,*) "     Noah.2.7.1, RUC.3.7, VIC.4.1.1, VIC.4.1.2 "
       write(LDT_logunit,*) "[ERR] No other surface models supported at this time ... stopping."
       call LDT_endrun()
    endif
    
    call convertCoarseRSTtoFineRST()

  end subroutine LDT_rstTransformProcInit
  

  subroutine convertCoarseRSTtoFineRST()

    implicit none

    integer               :: ftn, ftn2
    integer               :: k,i,t,m,c,r,kk
    integer               :: iret
    integer               :: nDims
    integer               :: nVars
    integer               :: nGlobalAtts
    integer               :: unlimdimID
    integer               :: nvardims
    integer               :: nLoop
    character*50          :: varName
    integer, allocatable  :: n_dimids(:)
    character*50          :: tIncr
    character*8           :: beg_date
    character*6           :: beg_time
    character*2           :: fd
    character*100         :: model_name   
    character*100         :: standard_name
    character*33          :: units
    real                  :: scale_factor
    real                  :: offset
    real                  :: vmin
    real                  :: vmax
    integer               :: xtimeId
    integer               :: tdimId
    integer               :: ntimes
    integer                   :: c1,c2,r1,r2,cc,rr
    integer                   :: rad, maxrad
    logical                   :: found
    character*50, allocatable :: dimName(:)
    integer,     allocatable  :: nvardimIds(:)
    real   ,     allocatable  :: var(:,:)
    real   ,     allocatable  :: var_new(:,:)
    real   ,     allocatable  :: var3d(:,:,:)
    integer,     allocatable  :: dims(:)
    integer,     allocatable  :: dimID(:),dimID2(:)
    real   ,     allocatable  :: var1_2d(:)
    logical*1,   allocatable  :: var1_b_2d(:)
    real   ,     allocatable  :: var2_2d(:)
    logical*1,   allocatable  :: var2_b_2d(:)
    real   ,     allocatable  :: rlat(:)
    real   ,     allocatable  :: rlon(:)
    real   ,     allocatable  :: w11(:)
    real   ,     allocatable  :: w12(:)
    real   ,     allocatable  :: w21(:)
    real   ,     allocatable  :: w22(:)
    real   ,     allocatable  :: n11(:)
    real   ,     allocatable  :: n12(:)
    real   ,     allocatable  :: n21(:)
    real   ,     allocatable  :: n22(:)

      ! Generate router model ensemble restart file:
    if(LDT_rc%rstsource.eq."LSM") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
        
       call LDT_verify(nf90_open(path=LDT_rc%inputrst,& 
            mode=nf90_NOWRITE,ncid=ftn),&
            'failed to open '//trim(LDT_rc%inputrst))
       
#if (defined USE_NETCDF4)
       call LDT_verify(nf90_create(path=LDT_rc%outputrst,cmode=nf90_netcdf4,&
            ncid = ftn2),&
            'creating netcdf file failed in LDT_ensRstMod')
#endif
#if (defined USE_NETCDF3)
       call LDT_verify(nf90_create(path=LDT_rc%outputrst,cmode=nf90_clobber,&
            ncid = ftn2),&
            'creating netcdf file failed in LDT_ensRstMod')
#endif
       call LDT_verify(nf90_inquire(ftn,nDims,nVars,nGlobalAtts,unlimdimId),&
            'nf90_inquire failed in LDT_ensRstMod')
       
       allocate(dimID(nDims))
       allocate(dimID2(nDims))
       allocate(dims(nDims))
       allocate(dimName(nDims))   
       
       call LDT_verify(nf90_inq_dimId(ftn,"ntiles",dimId(1)),&
            'nf90_inq_dimId failed for ntiles in LDT_ensRstMod')
       call LDT_verify(nf90_inquire_dimension(ftn,dimId(1),len=dims(1)),&
            'nf90_inquire_dimension failed in LDT_ensRstMod')
       
       allocate(rlat(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(rlon(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(n11(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(n12(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(n21(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(n22(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(w11(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(w12(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(w21(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
       allocate(w22(LDT_rc%lnc(2)*LDT_rc%lnr(2)))

       call neighbor_interp_input_withgrid (LDT_rc%gridDesc(1,:),&
            LDT_rc%gridDesc(2,:), LDT_rc%lnc(2)*LDT_rc%lnr(2),&
            rlat,rlon,n11)

       if(dims(1).ne.LDT_rc%npatch(1,1)) then 
          write(LDT_logunit,*) "[ERR] The tile dimension in the input restart file"
          write(LDT_logunit,*) "   does not match the tile dimension computed in LDT."
          call LDT_endrun()
       endif
       
       do k=2,nDims - 1
          ! EMK Fix formatting for JULES
          if (k-1 .lt. 10) then
             write(unit=fd,fmt='(I1)') k-1
          else
             write(unit=fd,fmt='(I2)') k-1
          end if
          dimName(k) = 'dim'//trim(fd)
          call LDT_verify(nf90_inq_dimId(ftn,dimName(k),dimId(k)),&
               'nf90_inq_dimId failed for '//dimName(k)//' in LDT_ensRstMod')
          call LDT_verify(nf90_inquire_dimension(ftn,dimID(k),len=dims(k)),&
               'nf90_inquire failed in LDT_ensRstMod')
       enddo
       
       call LDT_verify(nf90_inq_dimId(ftn,"time",tdimId),&
            'nf90_inq_dimId failed for time LDT_ensRstMod')

       call LDT_verify(nf90_inquire_dimension(ftn,tdimID,len=ntimes),&
            'nf90_inquire_dimension failed in LDT_ensRstMod')

       if(LDT_rc%lsm .eq. "Noah.3.3") then 
          model_name = "Noah version 3.3"
       elseif(LDT_rc%lsm .eq. "CLSMF2.5") then 
          model_name = "Catchment"
       elseif(LDT_rc%lsm .eq. "NoahMP.3.6") then
          model_name = "NOAHMP36"
       elseif(LDT_rc%lsm .eq. "NoahMP.3.9") then
          model_name = "NOAHMP39"       
       elseif(LDT_rc%lsm .eq. "JULES.5.0") then 
          model_name = "JULES50"
       endif
       
       dims(1) = LDT_rc%npatch(2,1)

       call writeglobalheader(ftn2, model_name, nDims, dims, &
            LDT_rc%nensem(2), dimID2)
       
       call LDT_verify(nf90_enddef(ftn2))
       
       do k=1,nVars
          call LDT_verify(nf90_inquire_variable(ftn,k,varName,ndims=nvardims),&
               'nf90_inquire_variable failed in LDT_ensRstMod')
          allocate(n_dimids(nvardims))
          call LDT_verify(nf90_inquire_variable(ftn,k,varName,dimids=n_dimids),&
               'nf90_inquire_variable failed in LDT_ensRstMod')
          if(varname.eq."time") then 
             call LDT_verify(nf90_get_att(ftn,k,"units",units),&
                  'nf90_get_att failed for time in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"time_increment",tincr),&
                  'nf90_get_att failed for time_increment in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"begin_date",beg_date),&
                  'nf90_get_att failed for begin_date in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"begin_time",beg_time),&
                  'nf90_get_att failed for begin_time in LDT_ensRstMod')
          else
             call LDT_verify(nf90_get_att(ftn,k,"units",units),&
                  'nf90_get_att failed in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"standard_name",standard_name),&
                  'nf90_get_att failed in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"scale_factor",scale_factor),&
                  'nf90_get_att failed in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"add_offset",offset),&
                  'nf90_get_att failed in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"vmin",vmin),&
                  'nf90_get_att failed in LDT_ensRstMod')
             call LDT_verify(nf90_get_att(ftn,k,"vmax",vmax),&
                  'nf90_get_att failed in LDT_ensRstMod')
          endif
          
          if(varname.eq."time") then

             call LDT_verify(nf90_def_var(ftn2,'time',&
                  nf90_float, dimids=ntimes, varID=xtimeId))

             call LDT_verify(nf90_put_att(ftn2,k,&
                  "units",trim(units)),&
                  'nf90_put_att failed for time units')
             call LDT_verify(nf90_put_att(ftn2,k,&
                  "time_increment",trim(tincr)),&
                  'nf90_put_att failed for time_increment')
             call LDT_verify(nf90_put_att(ftn2,k,&
                  "begin_date",trim(tincr)),&
                  'nf90_put_att failed for begin_date')
             call LDT_verify(nf90_put_att(ftn2,k,&
                  "begin_time",trim(tincr)),&
                  'nf90_put_att failed for begin_time')

             call LDT_verify(nf90_put_var(ftn2,xtimeID,0.0),&
                  'nf90_put_var failed for xtimeID')            
          else

             call writeheader_restart(ftn2,nvardims,&
                  n_dimIds,&
                  k,&
                  standard_name,units,&
                  scale_factor, offset,&
                  vmin,vmax)
          endif
          deallocate( n_dimids )
       enddo
       
       do k=1,nVars

          call LDT_verify(nf90_inquire_variable(ftn,k,ndims=nvardims),&
               'nf90_inquire_variable failed in LDT_ensRstMod')
          allocate(nvardimIds(nvardims))
          call LDT_verify(nf90_inquire_variable(ftn,k,varName,&
               dimIDs=nvardimIDs),&
               'nf90_inquire_variable failed in LDT_ensRstMod')
          
          if(varname.ne."time") then 
!----------------------------------------------------------------------------
! Assume that restart files do not have variables with more than 2 dimensions. 
!----------------------------------------------------------------------------
             if(nvardims.gt.1) then
                allocate(var(LDT_rc%npatch(1,1),dims(nvarDimIDs(2))))
                allocate(var_new(LDT_rc%npatch(2,1), dims(nvarDimIDs(2))))
             else
                allocate(var(LDT_rc%npatch(1,1),1))
                allocate(var_new(LDT_rc%npatch(2,1),1))
             endif
             
             call LDT_verify(nf90_get_var(ftn,k,var),&
                  'nf90_get_var failed in LDT_ensRstMod')

             allocate(var1_2d(LDT_rc%lnc(1)*LDT_rc%lnr(1)))
             allocate(var1_b_2d(LDT_rc%lnc(1)*LDT_rc%lnr(1)))
             
             allocate(var2_2d(LDT_rc%lnc(2)*LDT_rc%lnr(2)))
             allocate(var2_b_2d(LDT_rc%lnc(2)*LDT_rc%lnr(2)))


             nLoop = 1
             if(nvardims.gt.1) then 
                nLoop = dims(nvarDimIDs(2))
             endif
             do kk=1,nLoop

                var1_2d(:) = -9999.0
                var1_b_2d(:) = .false. 
                do t=1,LDT_rc%npatch(1,1)
                   r = LDT_surface(1,1)%tile(t)%row
                   c = LDT_surface(1,1)%tile(t)%col
                   var1_2d(c+(r-1)*LDT_rc%lnc(1)) = var(t,kk)
                   var1_b_2d(c+(r-1)*LDT_rc%lnc(1)) = .true.
                enddo


                call neighbor_interp(LDT_rc%gridDesc(2,:),&
                     var1_b_2d,var1_2d(:),var2_b_2d,var2_2d(:),&
                     LDT_rc%lnc(1)*LDT_rc%lnr(1),&
                     LDT_rc%lnc(2)*LDT_rc%lnr(2),&
                     rlat,rlon,n11,&
                     LDT_rc%udef,iret)


                !gapfilling
             
                do r=1,LDT_rc%lnr(2)
                   do c=1,LDT_rc%lnc(2)
                      if(LDT_LSMparam_struc(2)%landmask%value(c,r,1).gt.0.and.&
                           var2_2d(c+(r-1)*LDT_rc%lnc(2)).eq.LDT_rc%udef) then 
                         
                         found = .false. 
                         maxrad = 50
                         rad = 1
                         do while (.not.found) 
                            c1 = max(1,c-rad)
                            c2 = min(LDT_rc%lnc(2),c+rad)
                            r1 = max(1,r-rad)
                            r2 = min(LDT_rc%lnr(2),r+rad)
                            
                            do rr=r1,r2
                               do cc=c1,c2
                                  if(var2_2d(cc+(rr-1)*LDT_rc%lnc(2)).ne.&
                                       LDT_rc%udef) then 
                                     var2_2d(c+(r-1)*LDT_rc%lnc(2)) = &
                                          var2_2d(cc+(rr-1)*LDT_rc%lnc(2))
                                     found = .true. 
                                     exit;
                                  endif
                               enddo
                            enddo
                            rad = rad+1
                            if(rad.gt.maxrad) then 
                               write(LDT_logunit,*) "[ERR] neighbor search failed at c,r "  ,c,' ',r  
                               call LDT_endrun()
                            endif
                            
                         enddo
                      endif
                   enddo
                enddo
             
                do t=1,LDT_rc%npatch(2,1)
                   r = LDT_surface(2,1)%tile(t)%row
                   c = LDT_surface(2,1)%tile(t)%col
                   var_new(t,kk) = var2_2d(c+(r-1)*LDT_rc%lnc(2)) 
                   if(var_new(t,kk).eq.-9999.0) then
                      write(LDT_logunit,*) "[ERR] problem at c,r "  ,c,' ',r      
                      call LDT_endrun()
                   endif
                enddo
             enddo

             deallocate(var1_2d)
             deallocate(var1_b_2d)
             
             deallocate(var2_2d)
             deallocate(var2_b_2d)


             if(nvardims.gt.1) then 
                call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                     (/dims(1)*LDT_rc%nensem(2),dims(nvarDimIDs(2))/)),&
                     'nf90_put_var failed in LDT_ensRstMod for '//& 
                     trim(varName))
             else
                call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                     (/dims(1)*LDT_rc%nensem(2),1/)),&
                     'nf90_put_var failed in LDT_ensRstMod for '//&
                     trim(varName))
             endif
             
             deallocate(var)
             deallocate(var_new)
          endif
          deallocate(nvardimIDs)
       enddo
       
       deallocate(dimID)
       deallocate(dimID2)
       deallocate(dims)
       deallocate(dimName)
#endif
    endif
  end subroutine convertCoarseRSTtoFineRST
  
  subroutine writeglobalheader(ftn,model_name, &
       ndims, dims, nens, dimID)
    
      integer           :: ftn
      character(len=*)  :: model_name
      integer           :: nens
      integer           :: ndims
      integer           :: dims(ndims)
      integer           :: dimID(ndims)
      integer           :: k 
      character*50      :: dimName
      character*2       :: fd
      character(len=8)  :: date
      character(len=10) :: time
      character(len=5)  :: zone
      integer, dimension(8) :: values
! _________________________________________________


#if (defined USE_NETCDF3 || defined USE_NETCDF4)       
      call date_and_time(date,time,zone,values)


      call LDT_verify(nf90_def_dim(ftn,'ntiles',dims(1),&
           dimID(1)),&
           'nf90_def_dim failed for ntiles in LDT_ensRstMod')

      do k=2,nDims-1 
         ! EMK Fix format for JULES
         if (k-1 .lt. 10) then
            write(unit=fd,fmt='(I1)') k-1
         else
            write(unit=fd,fmt='(I2)') k-1
         end if
         dimName = 'dim'//trim(fd)
         call LDT_verify(nf90_def_dim(ftn,dimName,dims(k),&
              dimID(k)),&
              'nf90_def_dim failed for ntiles in LDT_ensRstMod')
      enddo
      call LDT_verify(nf90_def_dim(ftn,'time',1,&
           dimID(nDims)),&
           'nf90_def_dim failed for ntiles in LDT_ensRstMod')

      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"missing_value", &
           LDT_rc%udef),'nf90_put_att failed for missing_value')
      
      
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"title", &
           "LIS land surface model restart"),&
           'nf90_put_att failed for title')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"institution", &
           trim(LDT_rc%institution)),&
           'nf90_put_att failed for institution')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"source",&
           trim(model_name)),&
           'nf90_put_att failed for source')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"history", &
           "created on date: "//date(1:4)//"-"//date(5:6)//"-"//&
           date(7:8)//"T"//time(1:2)//":"//time(3:4)//":"//time(5:10)),&
           'nf90_put_att failed for history')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"references", &
           "Arsenault_etal_GMD_2018, Kumar_etal_EMS_2006"),&
           'nf90_put_att failed for references')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"conventions", &
           "CF-1.6"),'nf90_put_att failed for conventions') !CF version 1.6
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"comment", &
           "website: http://lis.gsfc.nasa.gov/"),&
           'nf90_put_att failed for comment')
      
!grid information
      if(LDT_rc%lis_map_proj(2).eq."latlon") then !latlon
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "EQUIDISTANT CYLINDRICAL"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(2,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(2,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(2,9)),&
              'nf90_put_att failed for DX')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(2,10)),&
                  'nf90_put_att failed for DY')
      elseif(LDT_rc%lis_map_proj(2).eq."mercator") then 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "MERCATOR"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(2,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(2,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON') 
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
              LDT_rc%gridDesc(2,10)),&
              'nf90_put_att failed for TRUELAT1')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
              LDT_rc%gridDesc(2,11)),&
              'nf90_put_att failed for STANDARD_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(2,8)),&
              'nf90_put_att failed for DX')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(2,9)),&
              'nf90_put_att failed for DY')
      elseif(LDT_rc%lis_map_proj(2).eq."lambert") then !lambert conformal
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "LAMBERT CONFORMAL"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(2,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(2,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
              LDT_rc%gridDesc(2,10)),&
              'nf90_put_att failed for TRUELAT1')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
              LDT_rc%gridDesc(2,7)),&
              'nf90_put_att failed for TRUELAT2')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
              LDT_rc%gridDesc(2,11)),&
              'nf90_put_att failed for STANDARD_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(2,8)),&
              'nf90_put_att failed for DX')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(2,9)),&
              'nf90_put_att failed for DY')
         
      elseif(LDT_rc%lis_map_proj(2).eq."polar") then ! polar stereographic
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
              "POLAR STEREOGRAPHIC"),&
              'nf90_put_att failed for MAP_PROJECTION')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LAT", &
              LDT_rc%gridDesc(2,4)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
              "SOUTH_WEST_CORNER_LON", &
              LDT_rc%gridDesc(2,5)),&
              'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
              LDT_rc%gridDesc(2,10)),&
              'nf90_put_att failed for TRUELAT1')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
              LDT_rc%gridDesc(2,7)),&
              'nf90_put_att failed for ORIENT')
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
              LDT_rc%gridDesc(2,11)),&
              'nf90_put_att failed for STANDARD_LON')                  
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
              LDT_rc%gridDesc(2,8)),&
              'nf90_put_att failed for DX')                  
         call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
              LDT_rc%gridDesc(2,9)),&
              'nf90_put_att failed for DY')                  
      endif
      
#endif
    end subroutine writeglobalheader
    
    subroutine writeheader_restart(ftn,nvardims,&
         dimID,varId,&
         standardName,units,&
         scale_factor, offset,&
         vmin,vmax)


      integer            :: ftn
      integer            :: nvardims
      integer            :: dimID(nvardims)
      integer            :: varID
      character(len=*)   :: standardName
      character(len=*)   :: units
      real               :: scale_factor
      real               :: offset
      real               :: vmin
      real               :: vmax

      integer :: shuffle, deflate, deflate_level
      integer :: dimID_t(2)

      shuffle = NETCDF_shuffle
      deflate = NETCDF_deflate
      deflate_level =NETCDF_deflate_level

      dimID_t(1) = dimID(1)
      
      if(nvardims.gt.1) then 
         dimID_t(2) = dimID(2)

         call LDT_verify(nf90_def_var(ftn,trim(standardName),&
              nf90_float, dimids = dimID_t(1:2), varID=varID),&
              'nf90_def_var(2d) failed in LDT_ensRstMod')
         
#if(defined USE_NETCDF4)
         call LDT_verify(nf90_def_var_deflate(ftn,&
               varID,shuffle, deflate, deflate_level),&
               'nf90_def_var_deflate(2d) failed in LDT_ensRstMod')
#endif
      else
         call LDT_verify(nf90_def_var(ftn,trim(standardName),&
              nf90_float,dimids = dimID_t(1:1), varID=varID),&
              'nf90_def_var(1d) failed in LDT_ensRstMod')
#if(defined USE_NETCDF4)                
         call LDT_verify(nf90_def_var_deflate(ftn,&
              varID, shuffle, deflate, deflate_level),&
              'nf90_def_var_deflate(1d) failed in LDT_ensRstMod')
#endif
      endif
      
      call LDT_verify(nf90_put_att(ftn,varID,&
           "units",trim(units)),&
           'nf90_put_att failed for units')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "standard_name",trim(standardName)))
      call LDT_verify(nf90_put_att(ftn,varID,&
           "long_name",trim(standardName)),&
           'nf90_put_att failed for long_name')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "scale_factor",scale_factor),&
           'nf90_put_att failed for scale_factor')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "add_offset",offset),&
           'nf90_put_att failed for add_offset')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "missing_value",LDT_rc%udef),&
           'nf90_put_att failed for missing_value')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "_FillValue",LDT_rc%udef),&
           'nf90_put_att failed for _FillValue')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "vmin",vmin),&
           'nf90_put_att failed for vmin')
      call LDT_verify(nf90_put_att(ftn,varID,&
           "vmax",vmax),&
           'nf90_put_att failed for vmax')
    end subroutine writeheader_restart



end module LDT_rstTransformProcMod
