!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !ROUTINE: read_ISMNsmobs
! \label{read_ISMNsmobs}
!
! !REVISION HISTORY:
!  21 Sep 2018   Sujay Kumar;   Initial Specification
!
! !INTERFACE: 
subroutine read_ISMNsmobs(Obj_Space)
! !USES: 
  use ESMF
  use LIS_coreMod,    only : LIS_rc, LIS_domain
  use LIS_timeMgrMod, only : LIS_calendar, LIS_tick
  use LIS_logMod,     only : LIS_logunit, LIS_verify, &
       LIS_getNextUnitNumber, LIS_releaseUnitNumber
  use LIS_fileIOMod,  only : LIS_readData
  use ISMNsm_obsMod,  only : ISMNsm_obs_struc
  use map_utils

#if (defined USE_NETCDF3 || defined USE_NETCDF4)  
  use netcdf
#endif

  implicit none
! !ARGUMENTS: 
  type(ESMF_State)    :: Obj_Space
!
! !DESCRIPTION:
!  
!  The arguments are: 
!  \begin{description}
!  \item[n]    index of the nest
!  \item[Obj\_Space] Objective Space
!  \end{description}
!
!EOP
  character*400         :: cline
  real,    pointer      :: smc(:)
  integer, pointer      :: nsmc(:)
  type(ESMF_Field)      :: smcField
  character*100         :: obsdir
  logical               :: data_update
  integer               :: i
  character*200         :: filename
  integer               :: ios
  integer               :: yr,doy,mo,da,hr,mn,ss
  logical               :: file_exists
  integer               :: status
  type(ESMF_Time)       :: obstime,obstime1
  integer               :: ftn
  integer               :: k,c,r,gid
  integer               :: tind
  integer               :: num_files 
  integer               :: n,v,iloc
  real                        :: lat, lon, elev
  real                        :: depthfrom, depthto
  real                        :: sm_value
  character*1                 :: sm_flag
  real                        :: sfsm, rzsm
  logical                     :: dataCheck
  real                        :: col,row
  integer                     :: stn_col,stn_row
  real,       allocatable     :: depth(:)
  real,       allocatable     :: sf_wt(:),rz_wt(:)
  real,       allocatable     :: obs_mask(:,:)
  integer                     :: dimID(2),maskid,iret
  integer                     :: syr, smo, sda, shr, smn, sss
  n = 1

  call ESMF_AttributeGet(Obj_Space,"Data Directory",&
       obsdir, rc=status)
  call LIS_verify(status)
  call ESMF_AttributeGet(Obj_Space,"Data Update Status",&
       data_update, rc=status)
  call LIS_verify(status)

  if(ISMNsm_obs_struc(n)%yr.ne.LIS_rc%yr) then 

     ISMNsm_obs_struc(n)%yr = LIS_rc%yr

     call getNumberOfISMNfiles(ISMNsm_obs_struc(n)%odir,&
          LIS_rc%yr, num_files)

     call generateISMNstationInfo(n, &
          ISMNsm_obs_struc(n)%odir, LIS_rc%yr, num_files)

     call ESMF_TimeSet(ISMNsm_obs_struc(n)%startTime, yy=LIS_rc%yr, &
          mm=1, dd=1, h=0, m=0,s = 0, &
          calendar=LIS_calendar, rc=status)
     call LIS_verify(status, 'ISMNsm starttime set failed')

     do k = 1,ISMNsm_obs_struc(n)%n_stns

        allocate(depth(ISMNsm_obs_struc(n)%stn(k)%vlevels))
        allocate(sf_wt(ISMNsm_obs_struc(n)%stn(k)%vlevels))
        allocate(rz_wt(ISMNsm_obs_struc(n)%stn(k)%vlevels))

        depth = LIS_rc%udef
        
        do v=1, ISMNsm_obs_struc(n)%stn(k)%vlevels
           
           inquire(file=trim(ISMNsm_obs_struc(n)%stn(k)%fname(v)),&
                exist=file_exists)

           if(file_exists) then 
              write(LIS_logunit,*) '[INFO] Reading ',&
                   trim(ISMNsm_obs_struc(n)%stn(k)%fname(v))
              ftn = LIS_getNextUnitNumber()
              open(ftn,file=trim(ISMNsm_obs_struc(n)%stn(k)%fname(v)),&
                   form='formatted')
              ios = 0 

              do while(ios.eq.0) 
                 read(ftn,'(a)',iostat=ios) cline
                 if(ios.ne.0.or.len(trim(adjustl(cline))).eq.0) exit

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),&
                      '(I4.4,1X,I2.2,1X,I2.2)') yr,mo,da
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),&
                      '(I2.2,1X,I2.2)') hr,mn
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) lat
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) lon
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) elev
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) depthfrom
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) depthto
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) sm_value
                 cline = adjustl(cline(iloc+1:len(cline)))

                 iloc = index(cline," ")
                 read(cline(1:iloc-1),*) sm_flag

                 call ESMF_TimeSet(obstime, yy=yr,mm=mo,dd=da,h=hr,m=mn,&
                      calendar=LIS_calendar,rc=status)
                 call LIS_verify(status, 'ESMF_TimeSet in readISMNObs')
                 
                 
                 tind = nint((obstime-ISMNsm_obs_struc(n)%startTime)/&
                      ISMNsm_obs_struc(n)%timestep)+1
                 ISMNsm_obs_struc(n)%stn(k)%lat = lat
                 ISMNsm_obs_struc(n)%stn(k)%lon = lon

                 if(tind.gt.0.and.tind.le.8784) then 
                    
!if we choose only good ('G') data, then we have no data. So including the
!unchecked ones as well. 
                    
                    if(sm_flag.eq.'G'.or.sm_flag.eq.'U') then 
                       
                       ISMNsm_obs_struc(n)%stn(k)%sm(tind,v) = sm_value                 
                       depth(v) = (depthfrom+depthto)/2.0
                    endif
                 endif
              enddo
              call LIS_releaseUnitNumber(ftn)
           endif
        enddo

        call compute_vinterp_weights(&
             ISMNsm_obs_struc(n)%stn(k)%vlevels,&
             ISMNsm_obs_struc(n)%lis_sf_d,&
             ISMNsm_obs_struc(n)%lis_rz_d,&
             depth,sf_wt,rz_wt)
        
        sfsm = 0 
        rzsm = 0 

        do tind = 1, ISMNsm_obs_struc(n)%nts
           dataCheck = .true. 
           do v=1,ISMNsm_obs_struc(n)%stn(k)%vlevels
              if(ISMNsm_obs_struc(n)%stn(k)%sm(tind,v).eq.LIS_rc%udef) then 
                 dataCheck = .false. 
              endif
           enddo
           if(dataCheck) then 
              sfsm = 0 
              rzsm = 0         
              do v=1,ISMNsm_obs_struc(n)%stn(k)%vlevels
                 if(ISMNsm_obs_struc(n)%stn(k)%sm(tind,v).ne.LIS_rc%udef) then 
                    sfsm = sfsm + sf_wt(v)*ISMNsm_obs_struc(n)%stn(k)%sm(tind,v)
                    rzsm = rzsm + rz_wt(v)*ISMNsm_obs_struc(n)%stn(k)%sm(tind,v)
                 endif
              enddo
              if(sfsm.gt.0.001) then 
                 ISMNsm_obs_struc(n)%stn(k)%sfsm(tind) = sfsm
              endif
              if(rzsm.gt.0.001) then 
                 ISMNsm_obs_struc(n)%stn(k)%rzsm(tind) = rzsm
              endif
           endif
        enddo
        deallocate(depth)
        deallocate(sf_wt)
        deallocate(rz_wt)
     end do

#if 0 
!computing the "obs mask" 
     allocate(obs_mask(LIS_rc%lnc(n),LIS_rc%lnr(n)))
     obs_mask = 0.0

     do i=1,ISMNsm_obs_struc(n)%n_stns
        call latlon_to_ij(LIS_domain(n)%lisproj,&
             ISMNsm_obs_struc(n)%stn(i)%lat, &
             ISMNsm_obs_struc(n)%stn(i)%lon,col,row)
        stn_col = nint(col)
        stn_row = nint(row)

        if(stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n).and.&
             stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n)) then 
           gid = LIS_domain(n)%gindex(stn_col, stn_row)
           
           if(gid.ne.-1) then 
              do k=1,ISMNsm_obs_struc(n)%nts
                 if(ISMNsm_obs_struc(n)%stn(i)%sfsm(k).gt.0) then 
                    obs_mask(stn_col,stn_row) = obs_mask(stn_col,stn_row) +&
                         ISMNsm_obs_struc(n)%stn(i)%sfsm(k)
                 endif
              enddo
           endif
        endif
     enddo
     do r=1,LIS_rc%lnr(n)
        do c=1,LIS_rc%lnc(n)
           if(obs_mask(c,r).le.0) then 
              obs_mask(c,r) = -9999.0
           else
              obs_mask(c,r) = 1.0
           endif
        enddo
     enddo

     iret = nf90_create('obsmask.nc',NF90_CLOBBER, &
          ftn)
     iret = nf90_def_dim(ftn,"east_west",LIS_rc%lnc(n),dimID(1))
     iret = nf90_def_dim(ftn,"north_south",LIS_rc%lnr(n),dimID(2))
     iret = nf90_def_var(ftn,'LANDMASK',NF90_FLOAT, dimID, maskid)
     iret = nf90_enddef(ftn)
     
     iret = nf90_put_var(ftn,maskid,obs_mask)
     iret = nf90_close(ftn)
     stop
#endif

     
  end if


  call ESMF_TimeSet(obstime1, yy=LIS_rc%yr, &
       mm=LIS_rc%mo, dd=LIS_rc%da, h=LIS_rc%hr, m=LIS_rc%mn, &
       s = LIS_rc%ss, calendar=LIS_calendar, rc=status)
  call LIS_verify(status, 'obstime1 set failed')


  call ESMF_StateGet(Obj_Space,"ISMN_sm",smcField,&
       rc=status)
  call LIS_verify(status)
  
  call ESMF_FieldGet(smcField,localDE=0,farrayPtr=smc,rc=status)
  call LIS_verify(status)


  allocate(nsmc(LIS_rc%ngrid(n)))
  smc = 0.0
  nsmc = 0

  do i=1,ISMNsm_obs_struc(n)%n_stns
     call latlon_to_ij(LIS_domain(n)%lisproj,&
          ISMNsm_obs_struc(n)%stn(i)%lat, &
          ISMNsm_obs_struc(n)%stn(i)%lon,col,row)
     stn_col = nint(col)
     stn_row = nint(row)

     if(stn_col.ge.1.and.stn_col.le.LIS_rc%lnc(n).and.&
          stn_row.ge.1.and.stn_row.le.LIS_rc%lnr(n)) then 
        gid = LIS_domain(n)%gindex(stn_col, stn_row)
     
        if(gid.ne.-1) then 
           
           tind = nint((obstime1 - ISMNsm_obs_struc(n)%starttime)/&
                ISMNsm_obs_struc(n)%timestep)+1
           
           if((tind.gt.0.and.tind.lt.8784).and.&
                LIS_rc%mn.eq.0.and.LIS_rc%ss.eq.0) then !only use at hourly intervals 
              smc(gid) = smc(gid) + ISMNsm_obs_struc(n)%stn(i)%sfsm(tind)
              nsmc(gid) = nsmc(gid) + 1
           endif
        endif
     endif
  end do

  do i=1,LIS_rc%ngrid(n)
     if(nsmc(i).gt.0) then 
        smc(i) = smc(i)/nsmc(i)
     else
        smc(i) = LIS_rc%udef
     endif
  enddo
  deallocate(nsmc)

  call ESMF_AttributeSet(Obj_Space,"Data Update Status",&
       .true., rc=status)
  call LIS_verify(status)

end subroutine read_ISMNsmobs

!BOP
!
! !ROUTINE: create_ISMNsmobs_filename
! \label(create_ISMNsmobs_filename)
!
! !INTERFACE:
subroutine create_ISMNsmobs_filename(odir, stn, yr, filename)
!
! !USES:
  implicit none
!
! !INPUT PARAMETERS:
  character(len=*), intent(in) :: odir
  character(len=*), intent(in) :: stn
  integer  ,        intent(in) :: yr
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
  character(len=*)             :: filename

  character*4 :: fyr

  write(unit=fyr, fmt='(i4.4)') yr

  filename = trim(odir)//'/'//trim(stn)//'_'//trim(fyr)//'.txt'

end subroutine create_ISMNsmobs_filename


!BOP
! 
! !ROUTINE: getNumberOfISMNfiles
! \label(getNumberOfISMNfiles)
!
! !INTERFACE:
subroutine getNumberOfISMNfiles(odir, yr, num_files)
! 
! !USES:  
  use LIS_coreMod
  use LIS_logMod

  implicit none
!
! !INPUT PARAMETERS: 
  character(len=*), intent(in) :: odir
  integer  ,        intent(in) :: yr
  integer                      :: num_files
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
  character*100                :: temp1
  character*500                :: ls_comm
  character*500                :: cmd2
  character*4                  :: fyr
  integer                      :: ftn
  character*1                  :: fproc(4) 

  write(unit=temp1,fmt='(i4.4)') LIS_localPet
  read(unit=temp1,fmt='(4a1)')fproc

  write(unit=fyr, fmt='(i4.4)') yr  
  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/*sm*stm > file_list.'//&
       fproc(1)//fproc(2)//fproc(3)//fproc(4)
  cmd2 = 'wc -w file_list.'//fproc(1)//fproc(2)//fproc(3)//fproc(4)//&
       ' > file_list.wc.'//fproc(1)//fproc(2)//fproc(3)//fproc(4)

  call system(ls_comm)
  call system(cmd2)

  ftn = LIS_getNextUnitNumber()
  open(ftn,file='file_list.wc.'//fproc(1)//fproc(2)//fproc(3)//fproc(4),&
       form='formatted',action='read')
  read(ftn,*) num_files
  call LIS_releaseUnitNumber(ftn)

end subroutine getNumberOfISMNfiles

subroutine generateISMNstationInfo(source, odir, yr, num_files)

  use LIS_coreMod
  use LIS_logMod
  use ISMNsm_obsMod

  implicit none

  integer  ,        intent(in) :: source
  character(len=*), intent(in) :: odir
  integer  ,        intent(in) :: yr
  integer  ,        intent(in) :: num_files


  character*100                :: temp1
  real                         :: depthfrom
  real                         :: depthto
  
  character*500, allocatable   :: filenames(:)
  character*500                :: filename2
  character*500                :: name1, name2, name3
  character*500                :: stnname
  character*500                :: checkString

  
  integer                      :: iloc
  integer                      :: n_stns
  character*500                :: ls_comm
  character*4                  :: fyr
  integer                      :: k 
  integer                      :: ftn
  character*1                  :: fproc(4) 

  write(unit=temp1,fmt='(i4.4)') LIS_localPet
  read(unit=temp1,fmt='(4a1)')fproc

  write(unit=fyr, fmt='(i4.4)') yr  

  ls_comm = 'ls '//trim(odir)//'/'//trim(fyr)//'/*sm*stm > file_list.'//&
       fproc(1)//fproc(2)//fproc(3)//fproc(4)

  call system(ls_comm)

  allocate(filenames(num_files))

  if(num_files.gt.0) then 
     ftn = LIS_getNextUnitNumber()
     do k=1,num_files
        open(ftn,file='file_list.'//fproc(1)//fproc(2)//fproc(3)//fproc(4),&
             form='formatted',action='read')
        read(ftn,'(a)') filenames(k)
     enddo
     call LIS_releaseUnitNumber(ftn)
  endif

  n_stns = 0 
  checkString = "dummy"

  do k=1, num_files
     filename2 = filenames(k)

     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name1
     filename2 = filename2(iloc+1:len(filename2))

     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name2
     filename2 = filename2(iloc+1:len(filename2))

     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name3
     filename2 = filename2(iloc+1:len(filename2))

     stnname = trim(name1)//"_"//trim(name2)//"_"//trim(name3)
     if(checkString.ne.stnname) then 
        n_stns = n_stns + 1
        checkString = stnname
     endif
  enddo

  ISMNsm_obs_struc(source)%n_stns = n_stns

  if(allocated(ISMNsm_obs_struc(source)%stn)) then 
     deallocate(ISMNsm_obs_struc(source)%stn)
  endif
  allocate(ISMNsm_obs_struc(source)%stn(ISMNsm_obs_struc(source)%n_stns))

  do k=1,ISMNsm_obs_struc(source)%n_stns

     ISMNsm_obs_struc(source)%stn(k)%vlevels = 0 

  enddo

  n_stns = 0 
  checkString = "dummy"
  do k=1, num_files
     filename2 = filenames(k)
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name1
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name2
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name3
     filename2 = filename2(iloc+1:len(filename2))
     
     stnname = trim(name1)//"_"//trim(name2)//"_"//trim(name3)
     
     if(checkString.ne.stnname) then 
        n_stns = n_stns + 1
        checkString = stnname
     endif         
   
     ISMNsm_obs_struc(source)%stn(n_stns)%vlevels = ISMNsm_obs_struc(source)%stn(n_stns)%vlevels +1

  enddo

  do k=1,ISMNsm_obs_struc(source)%n_stns 
     if(allocated(ISMNsm_obs_struc(source)%stn(k)%fname)) then 
        deallocate(ISMNsm_obs_struc(source)%stn(k)%fname)
        deallocate(ISMNsm_obs_struc(source)%stn(k)%sm)
        deallocate(ISMNsm_obs_struc(source)%stn(k)%sfsm)
        deallocate(ISMNsm_obs_struc(source)%stn(k)%rzsm)
     endif

     allocate(ISMNsm_obs_struc(source)%stn(k)%fname(ISMNsm_obs_struc(source)%stn(k)%vlevels))
     allocate(ISMNsm_obs_struc(source)%stn(k)%sm(ISMNsm_obs_struc(source)%nts,ISMNsm_obs_struc(source)%stn(k)%vlevels))
     allocate(ISMNsm_obs_struc(source)%stn(k)%sfsm(ISMNsm_obs_struc(source)%nts))
     allocate(ISMNsm_obs_struc(source)%stn(k)%rzsm(ISMNsm_obs_struc(source)%nts))
     ISMNsm_obs_struc(source)%stn(k)%vlevels = 0 
     ISMNsm_obs_struc(source)%stn(k)%sm = LIS_rc%udef
     ISMNsm_obs_struc(source)%stn(k)%sfsm = LIS_rc%udef
     ISMNsm_obs_struc(source)%stn(k)%rzsm = LIS_rc%udef
  enddo

  n_stns = 0
  checkString = "dummy"
  do k=1, num_files
     filename2 = filenames(k)
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name1
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name2
     filename2 = filename2(iloc+1:len(filename2))
     
     iloc = index(filename2,"_")
     read(filename2(1:iloc-1),'(a)') name3
     filename2 = filename2(iloc+1:len(filename2))
     
     stnname = trim(name1)//"_"//trim(name2)//"_"//trim(name3)
     
     if(checkString.ne.stnname) then 
        n_stns = n_stns + 1
        checkString = stnname
     endif         

     ISMNsm_obs_struc(source)%stn(n_stns)%vlevels = ISMNsm_obs_struc(source)%stn(n_stns)%vlevels +1
     ISMNsm_obs_struc(source)%stn(n_stns)%fname(ISMNsm_obs_struc(source)%stn(n_stns)%vlevels) = filenames(k)
  enddo

end subroutine generateISMNstationInfo
