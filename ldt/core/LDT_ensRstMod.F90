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
module LDT_ensRstMod
!BOP
!
! !MODULE: LDT_ensRstMod
! 
! !DESCRIPTION: 
!   The code in this file provides interfaces to manage 
!   the upscaling and downscaling of ensemble restarts
!
! !REVISION HISTORY: 
!  16 Aug 2012    Sujay Kumar  Initial Specification
!   1 Nov 2017 Augusto Getirana: Add HyMAP2 parameters
! 
  use ESMF
  use LDT_coreMod
  use LDT_logMod
  use LDT_ran2_gasdev
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
  use netcdf
#endif

  implicit none

  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LDT_ensRstInit
  public :: LDT_transformEnsRst
  
  contains

!BOP
! !ROUTINE: LDT_ensRstInit
! label{LDT_ensRstInit}
! 
! !INTERFACE: 
    subroutine LDT_ensRstInit
! 
! !DESCRIPTION: 
! 
!EOP  
  
      integer             :: n 
      integer             :: status
      n = 1
      
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%rstsource,&
           label="LIS restart source:",&
           rc=status)
      call LDT_verify(status,'LIS restart source: not defined')

      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%ensrstmode,&
           label="Ensemble restart generation mode:",&
           rc=status)
      call LDT_verify(status,'Ensemble restart generation mode: not defined')
      

      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%ensrstsampling,&
           label="Ensemble restart generation sampling strategy:",&
           rc=status)
      if(status.ne.0) then 
         write(LDT_logunit,*)'[ERR] Ensemble restart generation sampling strategy: not defined'
         write(LDT_logunit,*)"[ERR] options are ..'none' and 'random sampling'"

         call LDT_endrun()
      endif
      
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%inputrst,&
           label="Input restart filename:",&
           rc=status)
      call LDT_verify(status,'Input restart filename: not defined')
      
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%outputrst,&
           label="Output restart filename:",&
           rc=status)
      call LDT_verify(status,'Output restart filename: not defined')
      
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%nens_in,&
           label="Number of ensembles per tile (input restart):",&
           rc=status)
      call LDT_verify(status,&
           'Number of ensembles per tile (input restart): not defined')
      
      call ESMF_ConfigGetAttribute(LDT_config,LDT_rc%nens_out,&
           label="Number of ensembles per tile (output restart):",&
           rc=status)
      call LDT_verify(status,&
           'Number of ensembles per tile (output restart): not defined')
            
      if(LDT_rc%nensem(n).ne.LDT_rc%nens_in) then 
         write(LDT_logunit,*) '[ERR] Number of ensembles in the LDT domain is not '
         write(LDT_logunit,*) '  equal to the Number of ensembles per tile specified'
         write(LDT_logunit,*) '  for the input domain'
         call LDT_endrun()
      endif
      
      if(LDT_rc%ensrstmode.eq."upscale") then !going from 1 to many
         if(LDT_rc%nensem(n).gt.LDT_rc%nens_out) then 
            write(LDT_logunit,*) '[ERR] In the upscale mode, the number of ensemble '
            write(LDT_logunit,*) '  members in the target domain must be greater than the'
            write(LDT_logunit,*) '  number of members in the input domain.'
            call LDT_endrun()
         endif
      elseif(LDT_rc%ensrstmode.eq."downscale") then !going from 1 to many
         if(LDT_rc%nensem(n).lt.LDT_rc%nens_out) then 
            write(LDT_logunit,*) '[ERR] In the downscale mode, the number of ensemble '
            write(LDT_logunit,*) '  members in the target domain must be greater than the'
            write(LDT_logunit,*) '  number of members the input domain.'
            call LDT_endrun()
         endif
      endif
      
    end subroutine LDT_ensRstInit

!BOP
! 
! !ROUTINE: LDT_transformEnsRst
!  \label{LDT_transformEnsRst}
! 
! !INTERFACE: 
    subroutine LDT_transformEnsRst
!
! !DESCRIPTION: 
! 
!EOP
      if(LDT_rc%ensrstmode.eq."upscale") then !going from 1 to many
         call upscale_ensembleRst
      elseif(LDT_rc%ensrstmode.eq."downscale") then !going from many to 1
         call downscale_ensembleRst
      endif

    end subroutine LDT_transformEnsRst

!BOP
! 
! !ROUTINE: upscale_ensembleRst
!  \label{upscale_ensembleRst}
! 
! !INTERFACE: 
    subroutine upscale_ensembleRst()
!
! !DESCRIPTION:
!  This routine generates a restart file with larger number of ensembles
!  from an input file with smaller number of ensembles.
!EOP
      
      integer               :: n 
      integer               :: ios
      integer               :: ftn,ftn2
      integer               :: nDims
      integer               :: nVars
      integer               :: nGlobalAtts
      integer               :: unlimdimID
      integer               :: dimID
      integer,     allocatable  :: dimID2(:)
      character*50          :: dimName
      character*50          :: varName
      integer,     allocatable  :: dims(:)
      integer,     allocatable  :: n_dimids(:)
      integer               :: nvardims
      integer,     allocatable  :: nvardimIds(:)
      real   ,     allocatable  :: var(:,:)
      real   ,     allocatable  :: var_new(:,:)
      real   ,     allocatable  :: var3d(:,:,:)
      character*2           :: fd
      character*100         :: model_name   
      character*100         :: standard_name
      character*100         :: units
      real                  :: scale_factor
      real                  :: offset
      real                  :: vmin
      real                  :: vmax
      integer               :: tdimId
      integer               :: xtimeId
      character*8           :: beg_date
      character*6           :: beg_time
      character*50          :: tIncr
      integer               :: k,i,t,m
      integer               :: t1,t2,m1
      integer               :: numvars
      logical               :: file_exists
      integer               :: seed(NRANDSEED)
      real                  :: rand
      !ag (1Nov2017)
      real   ,     allocatable  :: var1d(:)
      integer,     allocatable  :: var_map(:)
! __________________________________________________

      n = 1
      seed = -1000
      
      write(LDT_logunit,*) " Generating restart file: ",trim(LDT_rc%outputrst)

      ! Generate router model ensemble restart file:
      if(LDT_rc%rstsource.eq."LSM") then 
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
      
         call LDT_verify(nf90_open(path=LDT_rc%inputrst,& 
              mode=nf90_NOWRITE,ncid=ftn),'failed to open '//trim(LDT_rc%inputrst))
         
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
         
         allocate(dimID2(nDims))
         allocate(dims(nDims))
         
         call LDT_verify(nf90_inq_dimId(ftn,"ntiles",dimId),&
              'nf90_inq_dimId failed for ntiles in LDT_ensRstMod')
         call LDT_verify(nf90_inquire_dimension(ftn,dimId,len=dims(1)),&
              'nf90_inquire_dimension failed in LDT_ensRstMod')
         
         if(dims(1).ne.LDT_rc%npatch(n,1)) then 
            write(LDT_logunit,*) "[ERR] The tile dimension in the input restart file"
            write(LDT_logunit,*) "  does not match the tile dimension computed in LDT."
            call LDT_endrun()
         endif
                  
         do k=2,nDims - 1
            ! EMK Fix formatting for JULES
            if (k-1 .lt. 10) then
               write(unit=fd,fmt='(I1)') k-1
            else
               write(unit=fd,fmt='(I2)') k-1
            end if
            dimName = 'dim'//trim(fd)
            call LDT_verify(nf90_inq_dimId(ftn,dimName,dimId),&
                 'nf90_inq_dimId failed for '//dimName//' in LDT_ensRstMod')
            call LDT_verify(nf90_inquire_dimension(ftn,dimID,len=dims(k)),&
                 'nf90_inquire failed in LDT_ensRstMod')
         enddo
         
         call LDT_verify(nf90_inq_dimId(ftn,"time",tdimId),&
              'nf90_inq_dimId failed for time LDT_ensRstMod')
         
         if(LDT_rc%lsm .eq. "Noah.3.3") then 
            model_name = "Noah version 3.3"
         elseif(LDT_rc%lsm .eq. "CLSMF2.5") then 
            model_name = "Catchment"
         elseif(LDT_rc%lsm .eq. "NoahMP.3.6") then
            model_name = "NOAHMP36"
         elseif(LDT_rc%lsm .eq. "JULES.5.0") then 
            model_name = "JULES50"
         endif
         
         call writeglobalheader(ftn2, model_name, nDims, dims, &
              LDT_rc%nens_in,LDT_rc%nens_out, dimID2)
         
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
                    nf90_float, dimids=tdimId, varID=xtimeId))
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
               
               call LDT_verify(nf90_put_var(ftn2,xtimeID,0.0))            
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

         if(LDT_rc%ensrstsampling.eq."random sampling") then
            allocate(var_map(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in))
            do i=1,dims(1)/LDT_rc%nens_in
               do m=1,LDT_rc%nens_out
                  t2 = (i-1)*LDT_rc%nens_out+m
                  if(m.gt.LDT_rc%nens_in) then 
                     !randomly select an ensemble member
                     call nr_ran2(seed, rand)
                     m1 = 1+nint(rand*(LDT_rc%nens_in-1))
                     t1 = (i-1)*LDT_rc%nens_in+m1
                     var_map(t2) = t1
                  else
                     t1 = (i-1)*LDT_rc%nens_in+m
                     var_map(t2) = t1
                  endif
               enddo
            enddo
         endif
         
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
                  allocate(var(dims(1),dims(nvarDimIDs(2))))
                  allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in, &
                       dims(nvarDimIDs(2))))
               else
                  allocate(var(dims(1),1))
                  allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1))
               endif
               
               call LDT_verify(nf90_get_var(ftn,k,var),&
                    'nf90_get_var failed in LDT_ensRstMod')

               if(LDT_rc%ensrstsampling.eq."none") then 
                  do t=1,dims(1)
                     do m=1,LDT_rc%nens_out
                        var_new((t-1)*LDT_rc%nens_out+m,:) = var(t,:)
                     enddo
                  enddo
               elseif(LDT_rc%ensrstsampling.eq."random sampling") then
                  do i=1,dims(1)/LDT_rc%nens_in
                     do m=1,LDT_rc%nens_out
                        t2 = (i-1)*LDT_rc%nens_out+m
                        t1 = var_map(t2)
                        var_new(t2,:) = var(t1,:)
                     enddo
                  enddo                  
               endif
               if(nvardims.gt.1) then 
                  call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                       (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,dims(nvarDimIDs(2))/)),&
                       'nf90_put_var failed in LDT_ensRstMod for '//& 
                       trim(varName))
               else
                  call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                       (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1/)),&
                       'nf90_put_var failed in LDT_ensRstMod for '//&
                       trim(varName))
               endif
               
               deallocate(var)
               deallocate(var_new)
            endif
            deallocate(nvardimIDs)
         enddo
         
         deallocate(dimID2)
         deallocate(dims)
         if(LDT_rc%ensrstsampling.eq."random sampling") then
            deallocate(var_map)
         endif
         call LDT_verify(nf90_close(ftn))
         call LDT_verify(nf90_close(ftn2))
#endif
         
      ! Generate router model ensemble restart file:
      elseif(LDT_rc%rstsource.eq."Routing") then
    
         ! HYMAP Router:
         !ag (1Nov2017)
         if( LDT_rc%routingmodel .eq. "HYMAP") then 

            write(LDT_logunit,*)"[INFO] 'Inflating' ensemble restart for routing model: "&
                  //trim(LDT_rc%routingmodel)

            allocate(var(LDT_rc%gnc(n),LDT_rc%gnr(n)))
            allocate(var3d(LDT_rc%gnc(n),LDT_rc%gnr(n),LDT_rc%nens_out))

            ftn = LDT_getNextUnitNumber()
            ftn2 = LDT_getNextUnitNumber()

            ! Check if input restart file is present:
            inquire( file=trim(LDT_rc%inputrst), exist=file_exists )
            if( file_exists ) then
               write(LDT_logunit,*) "[INFO] Opening HYMAP input restart file, "
               write(LDT_logunit,*)  trim(LDT_rc%inputrst)
               open(ftn,file=LDT_rc%inputrst,status='old',&
                    form='unformatted',access='sequential',iostat=ios)
            else
               write(LDT_logunit,*) "[ERR] HYMAP input binary restart file, "
               write(LDT_logunit,*)  trim(LDT_rc%inputrst)
               write(LDT_logunit,*) " is missing. Stopping run ..."
               call LDT_endrun
            endif

            ! If output ensemble restart file exists, provide warning ...
            inquire( file=trim(LDT_rc%outputrst), exist=file_exists )
            if(file_exists) then
               write(LDT_logunit,*) "[WARN] HYMAP binary restart file, "
               write(LDT_logunit,*)  trim(LDT_rc%outputrst)
               write(LDT_logunit,*) "  already exists! Overwriting the file ..." 
!               call LDT_endrun
            endif

            open(ftn2,file=LDT_rc%outputrst,status='new',&
                 form='unformatted',access='sequential',iostat=ios)

            numvars = 4  

            do i=1,numvars

               read(ftn) var

               do m=1,LDT_rc%nens_out
                  var3d(:,:,m) = var(:,:)
               enddo
            
               write(ftn2) var3d
            enddo

            call LDT_releaseUnitNumber(ftn)
            call LDT_releaseUnitNumber(ftn2)

            deallocate(var)
            deallocate(var3d)

         ! HYMAP2 Router
         elseif(LDT_rc%routingmodel .eq. "HYMAP2") then 

            write(LDT_logunit,*)"[INFO] 'Inflating' ensemble restart for routing model: "&
                  //trim(LDT_rc%routingmodel)

#if (defined USE_NETCDF3 || defined USE_NETCDF4)

         call LDT_verify(nf90_open(path=LDT_rc%inputrst,&
              mode=nf90_NOWRITE,ncid=ftn),'failed to open '//trim(LDT_rc%inputrst))

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
         
         allocate(dimID2(nDims))
         allocate(dims(nDims))

         call LDT_verify(nf90_inq_dimId(ftn,"ntiles",dimId),&
              'nf90_inq_dimId failed for ntiles in LDT_ensRstMod')
         call LDT_verify(nf90_inquire_dimension(ftn,dimId,len=dims(1)),&
              'nf90_inquire_dimension failed in LDT_ensRstMod')

         model_name = "HYMAP2"
         call writeglobalheader(ftn2, model_name, nDims, dims, &
              LDT_rc%nens_in,LDT_rc%nens_out, dimID2)

         call LDT_verify(nf90_enddef(ftn2))
         
         do k=1,nVars
            call LDT_verify(nf90_inquire_variable(ftn,k,varName,ndims=nvardims),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            allocate(n_dimids(nvardims))
            call LDT_verify(nf90_inquire_variable(ftn,k,varName,dimids=n_dimids),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            
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
            
            call writeheader_restart(ftn2,nvardims,&
                 n_dimIds,&
                 k,&
                 standard_name,units,&
                 scale_factor, offset,&
                 vmin,vmax)
            deallocate(n_dimids)
         enddo

         if(LDT_rc%ensrstsampling.eq."random sampling") then
           allocate(var_map(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in))
            do i=1,dims(1)/LDT_rc%nens_in
               do m=1,LDT_rc%nens_out
                  t2 = (i-1)*LDT_rc%nens_out+m
                  if(m.gt.LDT_rc%nens_in) then
                     !randomly select an ensemble member
                     call nr_ran2(seed, rand)
                     m1 = 1+nint(rand*(LDT_rc%nens_in-1))
                     t1 = (i-1)*LDT_rc%nens_in+m1
                     var_map(t2) = t1
                  else
                     t1 = (i-1)*LDT_rc%nens_in+m
                     var_map(t2) = t1
                  endif
               enddo
            enddo            
         endif
         
         do k=1,nVars
            call LDT_verify(nf90_inquire_variable(ftn,k,ndims=nvardims),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            allocate(nvardimIds(nvardims))
            call LDT_verify(nf90_inquire_variable(ftn,k,varName,&
                 dimIDs=nvardimIDs),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            
            if(nvardims.gt.1) then
               allocate(var(dims(1),dims(nvarDimIDs(2))))
               allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,&
                    dims(nvarDimIDs(2))))
            else
               allocate(var(dims(1),1))
               allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1))
            endif
            
            call LDT_verify(nf90_get_var(ftn,k,var),&
                 'nf90_get_var failed in LDT_ensRstMod')

            if(LDT_rc%ensrstsampling.eq."none") then 
               do t=1,dims(1)
                  do m=1,LDT_rc%nens_out
                     var_new((t-1)*LDT_rc%nens_out+m,:) = var(t,:)
                  enddo
               enddo
            elseif(LDT_rc%ensrstsampling.eq."random sampling") then 
               do i=1,dims(1)/LDT_rc%nens_in
                  do m=1,LDT_rc%nens_out
                     t2 = (i-1)*LDT_rc%nens_out+m
                     t1 = var_map(t2)
                     var_new(t2,:) = var(t1,:)
                  enddo
               enddo
            endif
            
            if(nvardims.gt.1) then
               call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                    (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,&
                    dims(nvarDimIDs(2))/)),&
                    'nf90_put_var failed in LDT_ensRstMod for '//&
                    trim(varName))
            else
               call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                    (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1/)),&
                    'nf90_put_var failed in LDT_ensRstMod for '//&
                    trim(varName))
            endif
            
            deallocate(var)
            deallocate(var_new)
            deallocate(nvardimIDs)
         enddo     
         deallocate(dimID2)
         deallocate(dims)

         if(LDT_rc%ensrstsampling.eq."random sampling") then
           deallocate(var_map)
         endif
         
         call LDT_verify(nf90_close(ftn))
         call LDT_verify(nf90_close(ftn2))

#endif

! The following code below is support for past HYMAP binary restart files
#if 0 
         write(LDT_logunit,*)"[INFO] 'Inflating' ensemble restart for routing model: "&
              //trim(LDT_rc%routingmodel)
         
         allocate(var1d(LDT_rc%routing_grid_count))
         allocate(var(LDT_rc%routing_grid_count,LDT_rc%nens_out))
         
         ftn = LDT_getNextUnitNumber()
         ftn2 = LDT_getNextUnitNumber()
         
         ! Check if input restart file is present:
         inquire( file=trim(LDT_rc%inputrst), exist=file_exists )
         if( file_exists ) then
            write(LDT_logunit,*) "[INFO] Opening HYMAP input restart file, "
            write(LDT_logunit,*)  trim(LDT_rc%inputrst)
            open(ftn,file=LDT_rc%inputrst,status='old',&
                 form='unformatted',iostat=ios)
            !open(ftn,file=LDT_rc%inputrst,status='old',&
            !     form='unformatted',access='sequential',iostat=ios)
         else
            write(LDT_logunit,*) "[ERR] HYMAP input binary restart file, "
            write(LDT_logunit,*)  trim(LDT_rc%inputrst)
            write(LDT_logunit,*) " is missing. Stopping run ..."
            call LDT_endrun
         endif
         
         ! If output ensemble restart file exists, provide warning...
         inquire( file=trim(LDT_rc%outputrst), exist=file_exists )
         if(file_exists) then
            write(LDT_logunit,*) "[WARN] If the HYMAP2 binary restart file, "
            write(LDT_logunit,*)  trim(LDT_rc%outputrst)
            write(LDT_logunit,*) "  already exists! Overwriting the file ..."
            !               call LDT_endrun
         endif
         
         !open(ftn2,file=LDT_rc%outputrst,status='new',&
         !     form='unformatted',access='sequential',iostat=ios)
         open(ftn2,file=LDT_rc%outputrst,status='new',&
              form='unformatted',iostat=ios)
         numvars = 8  
         
         do i=1,numvars
            
            read(ftn) var1d
            
            do m=1,LDT_rc%nens_out
               var(:,m) = var1d(:)
            enddo
            
            write(ftn2) var
         enddo
         
         call LDT_releaseUnitNumber(ftn)
         call LDT_releaseUnitNumber(ftn2)
         
         deallocate(var)
         deallocate(var1d)
! End of HYMAP restart binary file set of code
#endif

      else
         write(LDT_logunit,*) '[ERR] Ensemble restart for '//trim(LDT_rc%routingmodel)
         write(LDT_logunit,*) '   is not currently supported.'
         call LDT_endrun()
      endif

      ! Other ensemble restart sources are not currently supported
   else
      write(LDT_logunit,*) '[ERR] Ensemble restart source for '//trim(LDT_rc%rstsource)
      write(LDT_logunit,*) 'is not currently supported.'
      call LDT_endrun()
   endif
   
   write(LDT_logunit,*) " Successfully generated restart file: ",&
        trim(LDT_rc%outputrst)
   
 end subroutine upscale_ensembleRst
       
!BOP
! 
! !ROUTINE: downscale_ensembleRst
!  \label{downscale_ensembleRst}
! 
 ! !INTERFACE: 
 subroutine downscale_ensembleRst()
!
! !DESCRIPTION:
!  This routine generates a restart file with a smaller number of ensembles
!  from an input file with larger number of ensembles.
!EOP

   
   integer               :: n 
   integer               :: ios
   integer               :: ftn,ftn2
   integer               :: nDims
   integer               :: nVars, numvars
   integer               :: nGlobalAtts
   integer               :: unlimdimID
   integer               :: dimID
   integer,     allocatable  :: dimID2(:)
   character*50          :: dimName
   character*50          :: varName
   integer,     allocatable  :: dims(:)
   integer,     allocatable  :: n_dimids(:)
   integer               :: nvardims
   integer,     allocatable  :: nvardimIDs(:)
   real   ,     allocatable  :: var1d(:)
   real   ,     allocatable  :: var(:,:)
   real   ,     allocatable  :: var3d(:,:,:)
   real   ,     allocatable  :: var_new(:,:)
   character*2           :: fd
   character*100         :: model_name   
   character*100         :: standard_name
   character*100         :: units
   real                  :: scale_factor
   real                  :: offset
   real                  :: vmin
   real                  :: vmax
   integer               :: tdimId
   integer               :: xtimeId
   character*8           :: beg_date
   character*6           :: beg_time
   character*50          :: tIncr
   integer               :: k,t,m,i,p
   integer               :: m1, t1,t2
   integer               :: nc, nr
   integer               :: st, en
   logical               :: file_exists
   integer,     allocatable  :: var_map(:)
   logical                   :: cycl_check
   logical                   :: dupl_check
   integer               :: seed(NRANDSEED)
   real                  :: rand

   ! _______________________________________________________

   n = 1
   seed = -1000

   ! Generate single LSM model member (average):
   if(LDT_rc%rstsource.eq."LSM") then

      write(LDT_logunit,*)"[INFO] Downscaling ensemble restart for: "&
           //trim(LDT_rc%rstsource)

#if (defined USE_NETCDF3 || defined USE_NETCDF4) 

      call LDT_verify(nf90_open(path=LDT_rc%inputrst,& 
           mode=nf90_NOWRITE,ncid=ftn),'failed to open '//trim(LDT_rc%inputrst))

      write(LDT_logunit,*) " Generating restart file: ",trim(LDT_rc%outputrst)
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


      allocate(dimID2(nDims))
      allocate(dims(nDims))

      call LDT_verify(nf90_inq_dimId(ftn,"ntiles",dimId),&
           'nf90_inq_dimId failed for ntiles in LDT_ensRstMod')
      call LDT_verify(nf90_inquire_dimension(ftn,dimId,len=dims(1)),&
           'nf90_inquire_dimension failed in LDT_ensRstMod')

      if(dims(1).ne.LDT_rc%npatch(n,1)) then 
         write(LDT_logunit,*) "[ERR] The tile dimension in the input restart file"
         write(LDT_logunit,*) "  does not match the tile dimension computed in LDT."
         call LDT_endrun()
      endif

      do k= 2, nDims-1
         ! EMK Format fix for JULES
         if (k-1 .lt. 10) then
            write(unit=fd,fmt='(I1)') k-1
         else
            write(unit=fd,fmt='(I2)') k-1
         end if
         dimName = 'dim'//trim(fd)
         call LDT_verify(nf90_inq_dimId(ftn,dimName,dimId),&
              'nf90_inq_dimId failed for '//dimName//' in LDT_ensRstMod')
         call LDT_verify(nf90_inquire_dimension(ftn,dimID,len=dims(k)),&
              'nf90_inquire failed in LDT_ensRstMod')
      enddo

      call LDT_verify(nf90_inq_dimId(ftn,"time",tdimId),&
           'nf90_inq_dimId failed for time for LDT_ensRstMod')

      ! Does this need to be expanded?
      if(LDT_rc%lsm .eq. "Noah.3.3") then 
         model_name = "Noah version 3.3"
      elseif(LDT_rc%lsm .eq. "CLSMF2.5") then 
         model_name = "Catchment"
      elseif(LDT_rc%lsm .eq. "NoahMP.3.6") then 
         model_name = "NOAHMP36"
      elseif(LDT_rc%lsm .eq. "JULES.5.0") then 
         model_name = "JULES50"
      endif

      call writeglobalheader(ftn2, model_name, nDims, dims, &
           LDT_rc%nens_in,LDT_rc%nens_out, dimID2)
         
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
                 nf90_float, dimids=tdimId, varID=xtimeId))
            call LDT_verify(nf90_put_att(ftn2,k,&
                 "units",trim(units)),&
                 'nf90_put_att failed for time units in LDT_ensRstMod')
            call LDT_verify(nf90_put_att(ftn2,k,&
                 "time_increment",trim(tincr)),&
                 'nf90_put_att failed for time_increment in LDT_ensRstMod')
            call LDT_verify(nf90_put_att(ftn2,k,&
                 "begin_date",trim(tincr)),&
                 'nf90_put_att failed for begin_date in LDT_ensRstMod')
            call LDT_verify(nf90_put_att(ftn2,k,&
                 "begin_time",trim(tincr)),&
                 'nf90_put_att failed for begin_time in LDT_ensRstMod')

            call LDT_verify(nf90_put_var(ftn2,xtimeID,0.0))            
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

  
      if(LDT_rc%ensrstsampling.eq."random sampling") then
         allocate(var_map(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in))
         var_map = -1
         do i=1,dims(1)/LDT_rc%nens_in

            st = (i-1)*LDT_rc%nens_out+1
            en = (i-1)*LDT_rc%nens_out+LDT_rc%nens_out

            do m=1,LDT_rc%nens_out
               t2 = (i-1)*LDT_rc%nens_out+m
               cycl_check = .true. 

               do while(cycl_check) 
                  !randomly select an ensemble member                  
                  call nr_ran2(seed, rand)
                  m1 = 1+nint(rand*(LDT_rc%nens_in-1))
                  t1 = (i-1)*LDT_rc%nens_in+m1

                  dupl_check =.false. 
                  do p=st,en
                     if(t1.eq.var_map(p)) then
                        dupl_check = .true.
                     endif
                  enddo
                  if(dupl_check) then
                     cycl_check = .true.
                  else
                     cycl_check = .false.
                  endif
               enddo
               
               var_map(t2) = t1
            enddo
         enddo
      endif


      do k=1,nVars
         call LDT_verify(nf90_inquire_variable(ftn,k,ndims=nvardims),&
              'nf90_inquire_variable failed in LDT_ensRstMod')
         allocate(nvardimIDs(nvardims))
         call LDT_verify(nf90_inquire_variable(ftn,k,varName,&
              dimIDs=nvardimIDs),&
              'nf90_inquire_variable failed in LDT_ensRstMod')

         if(varname.ne."time") then 
!----------------------------------------------------------------------------
! Assume that restart files do not have variables with more than 2 dimensions. 
!----------------------------------------------------------------------------
            if(nvardims.gt.1) then
               allocate(var(dims(1),dims(nvarDimIDs(2))))
               allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in, &
                    dims(nvarDimIDs(2))))
            else
               allocate(var(dims(1),1))
               allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1))
            endif

            call LDT_verify(nf90_get_var(ftn,k,var),&
                 'nf90_get_var failed in LDT_ensRstMod')

            if(LDT_rc%ensrstsampling.eq."none") then
               do t=1,dims(1)/LDT_rc%nens_in
                  var_new(t,:) = 0.0
                  do m=1,LDT_rc%nens_in
                     var_new(t,:) = var_new(t,:) + var((t-1)*LDT_rc%nens_in+m,:)
                  enddo
                  var_new(t,:) = var_new(t,:)/LDT_rc%nens_in
               enddo
            elseif(LDT_rc%ensrstsampling.eq."random sampling") then
               do i=1,dims(1)/LDT_rc%nens_in
                  do m=1,LDT_rc%nens_out
                     t2 = (i-1)*LDT_rc%nens_out+m
                     t1 = var_map(t2)
                     var_new(t2,:) = var(t1,:)
                  enddo
               enddo
            endif

            if(nvardims.gt.1) then 
               call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                    (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,&
                    dims(nvarDimIDs(2))/)),&
                    'nf90_put_var failed in LDT_ensRstMod for '//& 
                    trim(varName))
            else

               call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                    (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1/)),&
                    'nf90_put_var failed in LDT_ensRstMod for '//&
                    trim(varName))
            endif

            deallocate(var)
            deallocate(var_new)
         endif
         deallocate(nvardimIDs)
      enddo

      deallocate(dimID2)
      deallocate(dims)

      if(LDT_rc%ensrstsampling.eq."random sampling") then
         deallocate(var_map)
      endif

#endif

      ! -- Routing Model Option --

      ! Generate router model ensemble restart file:
   elseif(LDT_rc%rstsource.eq."Routing") then

      ! HYMAP Router:
      !ag (1Nov2017)
      if( LDT_rc%routingmodel .eq. "HYMAP") then

         write(LDT_logunit,*)"[INFO] Downscaling ensemble restart for routing model: "&
              //trim(LDT_rc%routingmodel)

         ! Open input binary file:
         ftn = LDT_getNextUnitNumber()

         ! Check if input restart file is present:
         inquire( file=trim(LDT_rc%inputrst), exist=file_exists )
         if( file_exists ) then
            write(LDT_logunit,*) "[INFO] Opening HYMAP input restart file, "
            write(LDT_logunit,*)  trim(LDT_rc%inputrst)
            open(ftn,file=LDT_rc%inputrst,status='old',&
                 form='unformatted',access='sequential',iostat=ios)
         else
            write(LDT_logunit,*) "[ERR] HYMAP input binary restart file, "
            write(LDT_logunit,*)  trim(LDT_rc%inputrst)
            write(LDT_logunit,*) " is missing. Stopping run ..."
            call LDT_endrun
         endif

         ! If output ensemble restart file exists, provide warning ...
         inquire( file=trim(LDT_rc%outputrst), exist=file_exists )
         if(file_exists) then
            write(LDT_logunit,*) "[WARN] HYMAP binary restart file, "
            write(LDT_logunit,*)  trim(LDT_rc%outputrst)
            write(LDT_logunit,*) "  already exists! Overwriting the file ..."
            !               call LDT_endrun
         endif

         ! Create output binary file:
         write(LDT_logunit,*) " Generating restart file: ",trim(LDT_rc%outputrst)
         ftn2 = LDT_getNextUnitNumber()
         open(ftn2,file=LDT_rc%outputrst,status='new',&
              form='unformatted',access='sequential',iostat=ios)

         numvars = 4

         ! Loop over, read and average each HYMAP variable ensemble: 
         do i=1,numvars

            allocate(var3d(LDT_rc%gnc(n),LDT_rc%gnr(n),LDT_rc%nens_in))
            allocate(var_new(LDT_rc%gnc(n),LDT_rc%gnr(n)))

            read(ftn) var3d

            var_new = 0.0
            do nr = 1, LDT_rc%gnr(n)
               do nc = 1, LDT_rc%gnc(n)
                  if( count( mask = var3d(nc,nr,:) .ne. LDT_rc%udef ) == 0 ) then
                     var_new(nc,nr) = LDT_rc%udef
                  else 
                     var_new(nc,nr) = sum(var3d(nc,nr,:), mask = var3d(nc,nr,:) .ne. LDT_rc%udef) &
                          / count( mask = var3d(nc,nr,:) .ne. LDT_rc%udef )
                  endif
               enddo
            enddo

            write(ftn2) var_new

            deallocate(var3d)
            deallocate(var_new)

         enddo

         call LDT_releaseUnitNumber(ftn)
         call LDT_releaseUnitNumber(ftn2)

      ! HYMAP2 Router:
      elseif(LDT_rc%routingmodel .eq. "HYMAP2") then

         write(LDT_logunit,*)"[INFO] Downscaling ensemble restart for: "&
              //trim(LDT_rc%rstsource)
         
#if (defined USE_NETCDF3 || defined USE_NETCDF4) 
         
         call LDT_verify(nf90_open(path=LDT_rc%inputrst,& 
              mode=nf90_NOWRITE,ncid=ftn),'failed to open '//trim(LDT_rc%inputrst))
         
         write(LDT_logunit,*) " Generating restart file: ",trim(LDT_rc%outputrst)
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
         
         
         allocate(dimID2(nDims))
         allocate(dims(nDims))
         
         call LDT_verify(nf90_inq_dimId(ftn,"ntiles",dimId),&
              'nf90_inq_dimId failed for ntiles in LDT_ensRstMod')
         call LDT_verify(nf90_inquire_dimension(ftn,dimId,len=dims(1)),&
              'nf90_inquire_dimension failed in LDT_ensRstMod')
                  
         do k= 2, nDims-1
            if (k-1 .lt. 10) then
               write(unit=fd,fmt='(I1)') k-1
            else
               write(unit=fd,fmt='(I2)') k-1
            end if
            dimName = 'dim'//trim(fd)
            call LDT_verify(nf90_inq_dimId(ftn,dimName,dimId),&
                 'nf90_inq_dimId failed for '//dimName//' in LDT_ensRstMod')
            call LDT_verify(nf90_inquire_dimension(ftn,dimID,len=dims(k)),&
                 'nf90_inquire failed in LDT_ensRstMod')
         enddo
         
         ! Does this need to be expanded?
         model_name = "HYMAP2"

         call writeglobalheader(ftn2, model_name, nDims, dims, &
              LDT_rc%nens_in,LDT_rc%nens_out, dimID2)
         
         call LDT_verify(nf90_enddef(ftn2))
         
         do k=1,nVars
            
            call LDT_verify(nf90_inquire_variable(ftn,k,varName,ndims=nvardims),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            
            allocate(n_dimids(nvardims))
            call LDT_verify(nf90_inquire_variable(ftn,k,varName,dimids=n_dimids),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')

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

            call writeheader_restart(ftn2,nvardims,&
                 n_dimIds,&
                 k,&
                 standard_name,units,&
                 scale_factor, offset,&
                 vmin,vmax)

            deallocate( n_dimids )
         enddo

  
         if(LDT_rc%ensrstsampling.eq."random sampling") then
            allocate(var_map(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in))
            var_map = -1
            do i=1,dims(1)/LDT_rc%nens_in
               
               st = (i-1)*LDT_rc%nens_out+1
               en = (i-1)*LDT_rc%nens_out+LDT_rc%nens_out
               
               do m=1,LDT_rc%nens_out
                  t2 = (i-1)*LDT_rc%nens_out+m
                  cycl_check = .true. 
                  
                  do while(cycl_check) 
                     !randomly select an ensemble member                  
                     call nr_ran2(seed, rand)
                     m1 = 1+nint(rand*(LDT_rc%nens_in-1))
                     t1 = (i-1)*LDT_rc%nens_in+m1
                     
                     dupl_check =.false. 
                     do p=st,en
                        if(t1.eq.var_map(p)) then
                           dupl_check = .true.
                        endif
                     enddo
                     if(dupl_check) then
                        cycl_check = .true.
                     else
                        cycl_check = .false.
                     endif
                  enddo
                  
                  var_map(t2) = t1
               enddo
            enddo
         endif
         
         
         do k=1,nVars
            call LDT_verify(nf90_inquire_variable(ftn,k,ndims=nvardims),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            allocate(nvardimIDs(nvardims))
            call LDT_verify(nf90_inquire_variable(ftn,k,varName,&
                 dimIDs=nvardimIDs),&
                 'nf90_inquire_variable failed in LDT_ensRstMod')
            
!----------------------------------------------------------------------------
! Assume that restart files do not have variables with more than 2 dimensions. 
!----------------------------------------------------------------------------
            if(nvardims.gt.1) then
               allocate(var(dims(1),dims(nvarDimIDs(2))))
               allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in, &
                    dims(nvarDimIDs(2))))
            else
               allocate(var(dims(1),1))
               allocate(var_new(dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1))
            endif
            
            call LDT_verify(nf90_get_var(ftn,k,var),&
                 'nf90_get_var failed in LDT_ensRstMod')
            
            if(LDT_rc%ensrstsampling.eq."none") then
               do t=1,dims(1)/LDT_rc%nens_in
                  var_new(t,:) = 0.0
                  do m=1,LDT_rc%nens_in
                     var_new(t,:) = var_new(t,:) + var((t-1)*LDT_rc%nens_in+m,:)
                  enddo
                  var_new(t,:) = var_new(t,:)/LDT_rc%nens_in
               enddo
            elseif(LDT_rc%ensrstsampling.eq."random sampling") then
               do i=1,dims(1)/LDT_rc%nens_in
                  do m=1,LDT_rc%nens_out
                     t2 = (i-1)*LDT_rc%nens_out+m
                     t1 = var_map(t2)
                     var_new(t2,:) = var(t1,:)
                  enddo
               enddo
            endif
            
            if(nvardims.gt.1) then 
               call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                    (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,&
                    dims(nvarDimIDs(2))/)),&
                    'nf90_put_var failed in LDT_ensRstMod for '//& 
                    trim(varName))
            else
               
               call LDT_verify(nf90_put_var(ftn2,k,var_new,(/1,1/),&
                    (/dims(1)*LDT_rc%nens_out/LDT_rc%nens_in,1/)),&
                    'nf90_put_var failed in LDT_ensRstMod for '//&
                    trim(varName))
            endif
            
            deallocate(var)
            deallocate(var_new)
            deallocate(nvardimIDs)
         enddo
         deallocate(dimID2)
         deallocate(dims)

         if(LDT_rc%ensrstsampling.eq."random sampling") then
            deallocate(var_map)
         endif

#endif
      endif

      ! OTHER downscaled ensemble restart sources/methods are currently
      !  not supported/recognized at this time:
   else
      write(LDT_logunit,*) "[ERR] Ensemble restart source for "//trim(LDT_rc%rstsource)
      write(LDT_logunit,*) "  is not currently supported or recognized."
      call LDT_endrun()
   endif

   write(LDT_logunit,*) " Successfully generated restart file: ",&
        trim(LDT_rc%outputrst)

 end subroutine downscale_ensembleRst


 subroutine writeglobalheader(ftn,model_name, &
      ndims, dims, nens_in,nens_out, dimID)
   
   integer           :: n 
   integer           :: ftn
   character(len=*)  :: model_name
   integer           :: nens_in
   integer           :: nens_out
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

   n = 1
#if (defined USE_NETCDF3 || defined USE_NETCDF4)       
   call date_and_time(date,time,zone,values)
   if(LDT_rc%ensrstmode.eq."upscale") then 
      call LDT_verify(nf90_def_dim(ftn,'ntiles',dims(1)*nens_out/nens_in,&
           dimID(1)),&
           'nf90_def_dim failed for ntiles in LDT_ensRstMod')
   elseif(LDT_rc%ensrstmode.eq."downscale") then 
      call LDT_verify(nf90_def_dim(ftn,'ntiles',dims(1)*nens_out/nens_in,&
           dimID(1)),&
           'nf90_def_dim failed for ntiles in LDT_ensRstMod')
   endif

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
   if(LDT_rc%lis_map_proj(n).eq."latlon") then !latlon
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
           "EQUIDISTANT CYLINDRICAL"),&
           'nf90_put_att failed for MAP_PROJECTION')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LAT", &
           LDT_rc%gridDesc(n,4)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LON", &
           LDT_rc%gridDesc(n,5)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
      !  Add NORTH_EAST CORNER POINTS??
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
           LDT_rc%gridDesc(n,9)),&
           'nf90_put_att failed for DX')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
           LDT_rc%gridDesc(n,10)),&
           'nf90_put_att failed for DY')
   elseif(LDT_rc%lis_map_proj(n).eq."mercator") then 
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
           "MERCATOR"),&
           'nf90_put_att failed for MAP_PROJECTION')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LAT", &
           LDT_rc%gridDesc(n,4)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LON", &
           LDT_rc%gridDesc(n,5)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LON') 
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
           LDT_rc%gridDesc(n,10)),&
           'nf90_put_att failed for TRUELAT1')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
           LDT_rc%gridDesc(n,11)),&
           'nf90_put_att failed for STANDARD_LON')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
           LDT_rc%gridDesc(n,8)),&
           'nf90_put_att failed for DX')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
           LDT_rc%gridDesc(n,9)),&
           'nf90_put_att failed for DY')
   elseif(LDT_rc%lis_map_proj(n).eq."lambert") then !lambert conformal
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
           "LAMBERT CONFORMAL"),&
           'nf90_put_att failed for MAP_PROJECTION')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LAT", &
           LDT_rc%gridDesc(n,4)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LON", &
           LDT_rc%gridDesc(n,5)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
           LDT_rc%gridDesc(n,10)),&
           'nf90_put_att failed for TRUELAT1')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT2", &
           LDT_rc%gridDesc(n,7)),&
           'nf90_put_att failed for TRUELAT2')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
           LDT_rc%gridDesc(n,11)),&
           'nf90_put_att failed for STANDARD_LON')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
           LDT_rc%gridDesc(n,8)),&
           'nf90_put_att failed for DX')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
           LDT_rc%gridDesc(n,9)),&
           'nf90_put_att failed for DY')

   elseif(LDT_rc%lis_map_proj(n).eq."polar") then ! polar stereographic
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"MAP_PROJECTION", &
           "POLAR STEREOGRAPHIC"),&
           'nf90_put_att failed for MAP_PROJECTION')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LAT", &
           LDT_rc%gridDesc(n,4)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LAT')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,&
           "SOUTH_WEST_CORNER_LON", &
           LDT_rc%gridDesc(n,5)),&
           'nf90_put_att failed for SOUTH_WEST_CORNER_LON')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"TRUELAT1", &
           LDT_rc%gridDesc(n,10)),&
           'nf90_put_att failed for TRUELAT1')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"ORIENT", &
           LDT_rc%gridDesc(n,7)),&
           'nf90_put_att failed for ORIENT')
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"STANDARD_LON", &
           LDT_rc%gridDesc(n,11)),&
           'nf90_put_att failed for STANDARD_LON')                  
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DX", &
           LDT_rc%gridDesc(n,8)),&
           'nf90_put_att failed for DX')                  
      call LDT_verify(nf90_put_att(ftn,NF90_GLOBAL,"DY", &
           LDT_rc%gridDesc(n,9)),&
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

end module LDT_ensRstMod
