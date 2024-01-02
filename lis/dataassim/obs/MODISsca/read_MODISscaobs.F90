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
! !ROUTINE: read_MODISscaobs
! \label{read_MODISscaobs}
!
! !REVISION HISTORY:
!  31 Mar 08    Jiarui Dong; Initial Specification
!  10 Sep 08    Sujay Kumar; Adopted in LIS 6.0, included upscaling algorithm
!  17 Apr 09    Sujay Kumar; Added support for gap filled product
!
! !INTERFACE: 
subroutine read_MODISscaobs(n, OBS_State, OBS_Pert_State)
! !USES: 
  use ESMF
  use LIS_historyMod, only : LIS_readvar_gridded
  use LIS_coreMod,  only : LIS_rc, LIS_domain
  use LIS_logMod,     only : LIS_logunit, LIS_verify
  use LIS_constantsMod, only : LIS_CONST_PATH_LEN
  use MODISscaobs_module, only : MODISsca_obs_obj
  use LIS_mpiMod

  implicit none
! !ARGUMENTS: 
  integer, intent(in) :: n 
  type(ESMF_State)    :: OBS_State
  type(ESMF_State)    :: OBS_Pert_State
!
! !DESCRIPTION:
!  
!  reads the MODIS snow cover area observations 
!  from HDF-EOS files and packages it 
!  into an ESMF State with certain predefined 
!  attributes
!  
!  The arguments are: 
!  \begin{description}
!  \item[n] index of the nest
!  \item[OBS\_State] observations state
!  \end{description}
!
!EOP
  type(ESMF_Field)    :: scafield

  real,    pointer    :: obsl(:)
  real                :: tmp_obsl(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  integer             :: gid(LIS_rc%ngrid(n))
  integer             :: assimflag(LIS_rc%ngrid(n))

  character(len=LIS_CONST_PATH_LEN) :: obsdir, name
  logical             :: data_update
  logical             :: file_exists

  logical             :: readflag
  integer             :: status
  integer             :: t,c,r


  call ESMF_AttributeGet(OBS_State,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(OBS_State,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  file_exists = .false.
  if(LIS_rc%hr .eq. 10 .and. LIS_rc%mn .eq. 30) then  ! Terra
     if(MODISsca_obs_obj(n)%use_fill.eq.1) then 
        call MODISscaGF_filename(name,obsdir)     
     else
        call MODISsca_filename(name,obsdir)     
     endif
     inquire(file=name,exist=file_exists)
  endif

  if(file_exists) then 
     readflag = .true. 
  else 
     readflag = .false.
  endif

  if (readflag) then 
     write(LIS_logunit,*) 'Reading MODIS file ',trim(name)
     call ESMF_StateGet(OBS_State,"Observation01",scafield,&
          rc=status)
     call LIS_verify(status)

     call ESMF_FieldGet(scafield,localDE=0,farrayPtr=obsl,rc=status)
     call LIS_verify(status)
     
     call getMODISsca(n,name,tmp_obsl)

!     open(100,file='obs.bin',form='unformatted')
!     write(100) tmp_obsl
!     close(100)
!     stop
     obsl = LIS_rc%udef 

     do r =1,LIS_rc%lnr(n)
        do c =1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r) .ne. -1)then
              obsl(LIS_domain(n)%gindex(c,r))=tmp_obsl(c+LIS_rc%lnc(n)*(r-1))
           end if
        end do
     end do

     readflag = .false.

     do t=1,LIS_rc%ngrid(n)
        gid(t) = t
        if(obsl(t).ne.-9999.0) then 
           assimflag(t) = 1
        else
           assimflag(t) = 0
        endif
     enddo

     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .true. , rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(scafield,"Grid Number",&
          gid,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)

     call ESMF_AttributeSet(scafield,"Assimilation Flag",&
          assimflag,itemCount=LIS_rc%ngrid(n),rc=status)
     call LIS_verify(status)     
  else
     call ESMF_AttributeSet(OBS_State,"Data Update Status",&
          .false., rc=status)
     call LIS_verify(status)
     return
  endif
  
end subroutine read_MODISscaobs

!BOP
! !ROUTINE: MODISscaGF_filename
! \label{MODISscaGF_filename}
!
! !REVISION HISTORY:
!  17 Apr 09    Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine MODISscaGF_filename(name, ndir)
! !USES:   
  use LIS_coreMod,only : LIS_rc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
! 
! This routine creates the name of the MODIS sca observation 
! based on the current model date. 
! 
! The arguments are:
! \begin{description}
!  \item[name] Name of the MODIS sca observation file
!  \item[ndir] Name of the MODIS sca root directory
! \end{description}
!
!EOP

  integer            :: yr_lp, doy

  character (len=4)  :: fyr
  character (len=3)  :: fdoy

  integer :: DaysOfPrev(12)
  integer :: DaysOfPrev1(12)=(/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
  integer :: DaysOfPrev2(12)=(/ 0,31,60,91,121,152,182,213,244,274,305,335 /)
  
  write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
  
  do yr_lp = 1904, 2096, 4
    if ( LIS_rc%yr == yr_lp ) then;  DaysOfPrev=DaysOfPrev2
    else;   DaysOfPrev=DaysOfPrev1
    end if
  end do
    
  doy=LIS_rc%da+DaysOfPrev(LIS_rc%mo)
  write(unit=fdoy, fmt='(i3.3)') doy

  name = trim(ndir)//'/'//trim(fyr)//'/'//trim(fyr)//trim(fdoy)//'gapfill.hdf'
end subroutine MODISscaGF_filename


!BOP
! !ROUTINE: MODISsca_filename
! \label{MODISsca_filename}
!
! !REVISION HISTORY:
!  10 Sep 08    Sujay Kumar; Initial Specification
!
! !INTERFACE: 
subroutine MODISsca_filename(name, ndir)
! !USES:   
  use LIS_coreMod,only : LIS_rc
  use LIS_logMod, only : LIS_logunit

  implicit none
! !ARGUMENTS: 
  character(len=*)  :: name
  character (len=*) :: ndir
! 
! !DESCRIPTION: 
! 
! This routine creates the name of the MODIS sca observation 
! based on the current model date. 
! 
! The arguments are:
! \begin{description}
!  \item[name] Name of the MODIS sca observation file
!  \item[ndir] Name of the MODIS sca root directory
! \end{description}
!
!EOP
 
  integer           :: yr_lp, doy

  character (len=4) :: fyr
  character (len=3) :: fdoy

  integer :: DaysOfPrev(12)
  integer :: DaysOfPrev1(12)=(/ 0,31,59,90,120,151,181,212,243,273,304,334 /)
  integer :: DaysOfPrev2(12)=(/ 0,31,60,91,121,152,182,213,244,274,305,335 /)
  
  write(unit=fyr, fmt='(i4.4)') LIS_rc%yr
  
  do yr_lp = 1904, 2096, 4
    if ( LIS_rc%yr == yr_lp ) then;  DaysOfPrev=DaysOfPrev2
    else;   DaysOfPrev=DaysOfPrev1
    end if
  end do
    
  doy=LIS_rc%da+DaysOfPrev(LIS_rc%mo)
  write(unit=fdoy, fmt='(i3.3)') doy

  name = trim(ndir)//'/'//trim(fyr)//'/'//'MOD10C1.A'&
            //trim(fyr)//trim(fdoy)//'.005.hdf'
end subroutine MODISsca_filename

!BOP
!
! !ROUTINE: getMODISsca
! 
! !INTERFACE: 
subroutine getMODISsca(n,name,tmp_obsl)
! !USES: 
  use LIS_coreMod,only : LIS_rc,LIS_domain
  use LIS_logMod, only : LIS_logunit, LIS_verify
  use MODISscaobs_module, only : MODISsca_obs_obj
  
  implicit none
! 
! !DESCRIPTION: 
!  This routine retrievs the MODIS sca observations. The data is read
!  from a file in HDF-EOS format followed by the application of the 
!  confidence index and cloud cover flags. Finally the data is upscaled
!  or downscaled to the LIS domain
!EOP

#if (defined USE_HDF4) 
#include "hdf.f90"
#endif
  
  integer              :: n 
  character(len=*)     :: name
  real                 :: tmp_obsl(LIS_rc%lnc(n)*LIS_rc%lnr(n))

#if (defined USE_HDF4)
  integer, parameter   :: modis_nc=7200, modis_nr=3600
  integer              :: local_nc, local_nr

  !declare the hdf-eos library functions 
  integer              :: gdopen,gdattach,gddefboxreg,gdreginfo
  integer              :: gdextreg,gddetach,gdclose
  character*50         :: grid_name,ps_name,ci_name,pc_name,qa_name
  integer              :: ntype,rank,dims(2),size
  integer              :: file_id,grid_id,region_id,ret
  integer*1,allocatable    :: ps(:),ci(:),pc(:),qa(:)
  real*8               :: upleftpt(2),lowrightpt(2)
  real*8               :: cornerlon(2),cornerlat(2)

  integer              :: sfstart, sfselect, sfrdata, sfend
  integer              :: s_id1, s_id2, s_id3, sd_index, istat
  integer              :: stride(2), edges(2), start(2)
  integer*1, allocatable   :: data_map(:,:),cloud_map(:,:)
  integer*1, allocatable   :: cloud_pers_map(:,:)
  integer*1, allocatable   :: gf_snow_map(:,:)
  real,      allocatable   :: sca1(:)
  integer              :: c1,r1,c2,r2,c3,r3,c,r,gid
  
  if(MODISsca_obs_obj(n)%use_fill.ne.1) then 
     !Grid and field names
     grid_name ="MOD_CMG_Snow_5km"
     ps_name   ="Day_CMG_Snow_Cover"
     ci_name   ="Day_CMG_Confidence_Index"
     pc_name   ="Day_CMG_Cloud_Obscured"
     qa_name   ="Snow_Spatial_QA"

     !open the hdf file

     file_id = gdopen(trim(name),DFACC_READ)
     if (file_id.eq.-1)then
        write(LIS_logunit,*)"Failed to open hdf file",trim(name)
        return
     end if
     !  write(LIS_logunit,*) 'opened file',file_id

     !get the grid id
     grid_id = gdattach(file_id,grid_name)
     if (grid_id.eq.-1)then
        write(LIS_logunit,*)"Failed to attach grid: ",grid_name,trim(name)
        ret = gdclose(file_id)
        return
     end if
     !  write(LIS_logunit,*) 'gdattach',grid_id

     !get the LIS domain
     cornerlat(1)=MODISsca_obs_obj(n)%gridDesci(4)
     cornerlon(1)=MODISsca_obs_obj(n)%gridDesci(5)
     cornerlat(2)=MODISsca_obs_obj(n)%gridDesci(7)
     cornerlon(2)=MODISsca_obs_obj(n)%gridDesci(8)

     region_id = gddefboxreg(grid_id,cornerlon,cornerlat)
     if (region_id <0)then
        write(LIS_logunit,*)"Failed to obtain region id: gddefboxreg"
        ret = gdclose(file_id)
        return
     end if
     !  write(LIS_logunit,*) 'gddefboxreg',region_id

     ! find the dimensions of the subregion:dims(2)
     ret = gdreginfo(grid_id,region_id,ps_name,ntype,rank,dims,size,upleftpt,lowrightpt)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get region info: gdreginfo"
        ret = gdclose(file_id)
        return
     end if

     !    write(LIS_logunit,*) 'gdreginfo',ret

     ! get the percent snow (ps) cover field
     allocate(ps(dims(1)*dims(2)))
     ret = gdextreg(grid_id,region_id,ps_name,ps)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the ps field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        return
     end if

     ! get the confidence index (ci) field
     allocate(ci(dims(1)*dims(2)))
     ret = gdextreg(grid_id,region_id,ci_name,ci)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the CI field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        return
     end if

     ! get the percent cloud (pc) cover field
     allocate(pc(dims(1)*dims(2)))
     ret = gdextreg(grid_id,region_id,pc_name,pc)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the cloud field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        return
     end if

     ! get the qa flag of the grid
     allocate(qa(dims(1)*dims(2)))
     ret = gdextreg(grid_id,region_id,qa_name,qa)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to get the QA field"
        ret=gddetach(grid_id)
        ret=gdclose(file_id)
        return
     end if

     call calculate_sca(n, dims(1), dims(2), ps,ci,pc, tmp_obsl)

     deallocate(ps)
     deallocate(ci)
     deallocate(pc)
     deallocate(qa)

     ret=gddetach(grid_id)
     if (ret <0)then
        write(LIS_logunit,*)"Failed to detach grid_id: ",grid_id
     end if
     ret=gdclose(file_id)
     ! write(LIS_logunit,*) 'gdclose',file_id
     if (ret <0)then
        write(LIS_logunit,*)"Failed to close file: ",file_id
     end if
  else
     !open the hdf file
     
     file_id = sfstart(trim(name),DFACC_READ)
     if (file_id.eq.-1)then
        call LIS_verify(file_id,'Failed to open file '//trim(name))
     end if

     sd_index = 0 ! first record
     s_id1 = sfselect(file_id,sd_index)

     if(s_id1.eq.-1) then 
        call LIS_verify(s_id1,'Failed to select record')
     end if

     start = 0 
     edges(1) = modis_nc
     edges(2) = modis_nr
     stride(1) = 1
     stride(2) = 1
     
     allocate(data_map(modis_nc,modis_nr))
     allocate(gf_snow_map(modis_nc,modis_nr))
     allocate(cloud_map(modis_nc,modis_nr))
     allocate(cloud_pers_map(modis_nc,modis_nr))

     istat = sfrdata(s_id1,start,stride,edges,data_map)

     do r=1,modis_nr
        do c=1,modis_nc
           gf_snow_map(c,r) = data_map(c,modis_nr-r+1)
        enddo
     enddo

     if(istat.ne.0) then 
        call LIS_verify(istat,'Data Extraction of gap filled snow failed')
     endif

     sd_index = 1
     s_id2 = sfselect(file_id,sd_index)

     istat = sfrdata(s_id2,start,stride,edges,data_map)

     do r=1,modis_nr
        do c=1,modis_nc
           cloud_map(c,r) = data_map(c,modis_nr-r+1)
        enddo
     enddo

     if(istat.ne.0) then 
        call LIS_verify(istat,'Data Extraction of cloud map failed')
     endif

     sd_index = 3
     s_id3 = sfselect(file_id,sd_index)

     istat = sfrdata(s_id3,start,stride,edges,data_map)

     do r=1,modis_nr
        do c=1,modis_nc
           cloud_pers_map(c,r) = data_map(c,modis_nr-r+1)
        enddo
     enddo
     deallocate(data_map)

     if(istat.ne.0) then 
        call LIS_verify(istat,'Data Extraction of cloud persistence map failed')
     endif
     
     istat = sfend(file_id)

     !subset the snow map based on the thresholds

     c1 = nint((MODISsca_obs_obj(n)%gridDesci(5)+179.975)/0.05)+1
     r1 = nint((MODISsca_obs_obj(n)%gridDesci(4)+89.975)/0.05)+1

     c2 = nint((MODISsca_obs_obj(n)%gridDesci(8)+179.975)/0.05)+1
     r2 = nint((MODISsca_obs_obj(n)%gridDesci(7)+89.975)/0.05)+1


     local_nc = MODISsca_obs_obj(n)%gridDesci(2)
     local_nr = MODISsca_obs_obj(n)%gridDesci(3)

!test
!     allocate(sca1(modis_nc*modis_nr))
!     do r=1,modis_nr
!        do c=1,modis_nc
!           sca1(c+modis_nr*(r-1)) = gf_snow_map(c,r)
!        enddo
!     enddo
!     open(100,file='sca.bin',form='unformatted')
!     write(100) gf_snow_map
!     close(100)
!     stop

     allocate(sca1(local_nc*local_nr))
     sca1 = LIS_rc%udef

     do r=r1,r2
        do c=c1,c2
!           if(cloud_map(c,r).le.MODISsca_obs_obj(n)%cloud_thres.and.&
           if(cloud_pers_map(c,r).le.MODISsca_obs_obj(n)%cloud_pers_thres) then 
              c3 = c-c1+1
              r3 = r-r1+1
!              print*, c,r,cloud_map(c,r), cloud_pers_map(c,r)
              gid = c3+local_nc*(r3-1)
              sca1(gid) = gf_snow_map(c,r)
           endif
        enddo
     enddo
!     stop
    
 !    open(100,file='sca.bin',form='unformatted')
 !    write(100) sca1
 !    close(100)

     deallocate(gf_snow_map)
     deallocate(cloud_map)
     deallocate(cloud_pers_map)
     
     call calculate_sca_gf(n, local_nc, local_nr, sca1, tmp_obsl)
     
!     open(100,file='sca.bin',form='unformatted')
!     write(100) tmp_obsl
!     close(100)
!     stop
  endif

#endif
end subroutine getMODISsca


subroutine calculate_sca_gf(n, dim1, dim2, sca, tmp_obsl)

  use MODISscaobs_module, only: MODISsca_obs_obj
  use LIS_coreMod, only : LIS_rc, LIS_domain
  use LIS_logMod,  only : LIS_verify

  implicit none

  integer             :: n 
  integer             :: dim1
  integer             :: dim2
  real                :: sca(dim1*dim2)
  real                :: tmp_obsl(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  integer            :: c,r
  integer            :: iret
  logical*1, allocatable :: li(:)
  logical*1          :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  if(LIS_rc%gridDesc(n,10).le.0.05) then 
     allocate(li(dim1*dim2))
     li = .false.
     do r=1,dim2
        do c=1,dim1
           if(sca(c+(r-1)*dim1).ne.-9999.0) then 
              li(c+(r-1)*dim1) = .true.
           endif
        enddo
     enddo
     call bilinear_interp(LIS_rc%gridDesc(n,:), li, sca, lo, tmp_obsl,&
          MODISsca_obs_obj(n)%mi, LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          MODISsca_obs_obj(n)%w11,MODISsca_obs_obj(n)%w12,&
          MODISsca_obs_obj(n)%w21,MODISsca_obs_obj(n)%w22,&
          MODISsca_obs_obj(n)%n11,MODISsca_obs_obj(n)%n12,&
          MODISsca_obs_obj(n)%n21,MODISsca_obs_obj(n)%n22,LIS_rc%udef,iret)
     deallocate(li)
  else
     call LIS_verify(1,'Currently only downscaling of MODIS  data is supported')
  endif
end subroutine calculate_sca_gf

!BOP
! 
! !ROUTINE: calculate_sca
! \label{calculate_sca}
!
! !INTERFACE: 
subroutine calculate_sca(n, dim1, dim2, ps,ci,pc,go)
! !USES: 
  use MODISscaobs_module, only: MODISsca_obs_obj
  use LIS_coreMod, only : LIS_rc, LIS_domain

  implicit none
! !ARGURMENTS: 
  integer,  intent(IN) :: n 
  integer,  intent(IN) :: dim1
  integer,  intent(IN) :: dim2
  integer*1            :: ps(dim1*dim2)
  integer*1            :: ci(dim1*dim2)
  integer*1            :: pc(dim1*dim2)
  real                 :: go(LIS_rc%lnc(n)*LIS_rc%lnr(n))
! 
! !DESCRIPTION: 
!  This routine processes the raw sca observations by eliminating
!  the cloud covered points and points of low CI. This routine also
!  upscales/downscales the data MODIS sca data to the LIS resolution
! 
!EOP

  logical*1      :: lo(LIS_rc%lnc(n)*LIS_rc%lnr(n))

  logical*1, allocatable :: li(:)

  integer        :: i
  real           :: sca(dim1*dim2)
  real, allocatable  :: sca1(:)

  integer        :: iret
  integer        :: c,r

  do i = 1,dim1*dim2
    if ((ps(i) >= 0.0001 .and. ps(i)<=100).and.&
         (pc(i).le.10.and.pc(i).ge.0001))then
       sca(i)=real(ps(i))*100.0/real(ci(i))
    else
      sca(i)=-9999.0
    end if
  end do

!  open(120,file='obs1.bin',form='unformatted')
!  write(120) sca
!  close(120)
 
  if(LIS_rc%gridDesc(n,10).le.0.05) then 
     allocate(sca1(dim1*dim2))
     allocate(li(dim1*dim2))
     li = .false.
     do r=1,dim2
        do c=1,dim1
           sca1(c+(r-1)*dim1) = sca(c+(dim2-r+1-1)*dim1)
           if(sca1(c+(r-1)*dim1).ne.-9999.0) then 
              li(c+(r-1)*dim1) = .true.
           endif
        enddo
     enddo
     call bilinear_interp(LIS_rc%gridDesc(n,:), li, sca1, lo, go,&
          MODISsca_obs_obj(n)%mi, LIS_rc%lnc(n)*LIS_rc%lnr(n), &
          LIS_domain(n)%lat, LIS_domain(n)%lon,&
          MODISsca_obs_obj(n)%w11,MODISsca_obs_obj(n)%w12,&
          MODISsca_obs_obj(n)%w21,MODISsca_obs_obj(n)%w22,&
          MODISsca_obs_obj(n)%n11,MODISsca_obs_obj(n)%n12,&
          MODISsca_obs_obj(n)%n21,MODISsca_obs_obj(n)%n22,LIS_rc%udef,iret)
     deallocate(sca1)
     deallocate(li)
  else
     allocate(sca1(dim1*dim2))
     allocate(li(dim1*dim2))

     li = .false. 
     do r=1,dim2
        do c=1,dim1
           sca1(c+(r-1)*dim1) = sca(c+(dim2-r+1-1)*dim1)
           if(sca1(c+(r-1)*dim1).ne.-9999.0) then 
              li(c+(r-1)*dim1) = .true.              
!              print*, c+(r-1)*dim1, sca(c+(dim2-r+1-1)*dim1), li(c+(r-1)*dim1)

           endif
        enddo
     enddo
!     print*, dim1, dim2
!     open(100,file='sca.bin',form='unformatted')
!     write(100) sca1
!     close(100)
!     stop
     call upscaleByAveraging(MODISsca_obs_obj(n)%mi,&
          LIS_rc%lnc(n)*LIS_rc%lnr(n),LIS_rc%udef, &
          MODISsca_obs_obj(n)%n11,li, sca1, lo,go)
     deallocate(sca1)
     deallocate(li)

!     open(100,file='snow.bin',form='unformatted')
!     write(100) go
!     close(100)
!     stop
  endif
end subroutine calculate_sca
