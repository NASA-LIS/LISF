!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_snotel
! \label{read_snotel}
! 
! !REVISION HISTORY:
!
!  08 Jun 2011: Yuqiong Liu; Initial Specification, based on LIS SNOTEL obs reader by S. Kumar 
!
! !INTERFACE:      
subroutine read_snotel(n,ftn,findex,order) 
! !USES:
  use LIS_logMod, only         : LIS_logunit, LIS_endrun
  use LIS_coreMod, only        : LIS_rc,LIS_domain
  use LIS_metforcingMod, only : LIS_forc
  use LIS_constantsMod,  only : LIS_CONST_PATH_LEN
  use snotel_forcingMod,    only : snotel_struc
  use map_utils,    only : latlon_to_ij

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: n
  integer, intent(in)           :: ftn
  integer, intent(in)           :: findex
  integer, intent(in)           :: order
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  the correct SNOTEL station data (ASCII), transforms into LIS forcing 
!  parameters and interpolates to the LIS domain.
!
!  The arguments are: 
!  \begin{description}
!  \item[order]
!    flag indicating which data to be read (order=1, read the previous 
!    instance, order=2, read the next instance)
!  \item[n]
!    index of the nest
!  \item[ftn]
!    unit number for the SNOTEL station data 
!  \end{description}
! 
!  The routines invoked are: 
!  \begin{description}
!  \item[normalize\_stnwts](\ref{normalize_stnwts}) \newline
!    renormalizes the station weights accounting for
!    missing data
!  \item[interp\_stndata](\ref{interp_stndata}) \newline
!    spatially interpolates the station data onto the 
!    LIS grid.
!  \end{description}
!EOP  

  integer :: c,r,count1,i,j
  real :: pcp(snotel_struc(n)%nstns),tmppcp
!  real :: varfield(npts)
  real :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer :: npcp(LIS_rc%lnc(n),LIS_rc%lnr(n))
  character(len=LIS_CONST_PATH_LEN) :: snotel_filename
  character(len=500) :: line
  integer :: yr,num,dum,hr,mo,da,mint,sec
  logical :: file_exists, readflag
  integer :: ios,kk,ios1,iloc
  real    :: swe_data, prcp_data
  character*10  :: datestring
  real    :: col, row
  integer :: stn_col, stn_row

  pcp = snotel_struc(n)%undef
  varfield1 = 0
  npcp = 0

  do i=1,snotel_struc(n)%nstns
!generate the snotel filename
     yr = LIS_rc%yr
     if (LIS_rc%mo .ge. 10) yr = yr+1
    
     call create_SNOTEL_filename(trim(snotel_struc(n)%snoteldir),snotel_struc(n)%statename(i), &
          snotel_struc(n)%stnid(i), yr, snotel_filename)

     inquire(file=snotel_filename,exist=file_exists)
     if(file_exists) then 
        write(LIS_logunit,*) 'Reading SNOTEL file ',trim(snotel_filename)
        open(ftn,file=snotel_filename,form='formatted',status='old')
        read(ftn,*) ! skip the header line

        readflag = .true.
        do while (readflag)
           read(ftn,'(a)',iostat=ios)line
           if(ios.ne.0) then
              readflag = .false.
              write(LIS_logunit,*) '****** ios=0'
              exit
           endif
           do kk=1,len(line)
              if(line(kk:kk)==achar(9)) line(kk:kk) =';'
           enddo
           iloc = index(line,";")
           read(line(1:iloc-1),*) datestring
           line = line(iloc+1:len(line))

          ! read(datestring,'(i2.2,i2.2,i2.2)') mo,da,yr
           read(datestring,'(i2,i2,i2)') mo,da,yr
           if(yr.gt.60) then
               yr = yr + 1900
           else
               yr = yr + 2000
           endif

           if ( LIS_rc%da== da .and. LIS_rc%mo == mo .and. LIS_rc%yr == yr ) then  
              
              iloc = index(line,";")
              read(line(1:iloc-1),*,iostat=ios1) swe_data
              line = line(iloc+1:len(line))
              if(ios1.ne.0) swe_data = -9999.0

              do kk=1,4 !skip
                 iloc = index(line,";")
                 line = line(iloc+1:len(line))
              enddo

              read(line,*,iostat=ios1) prcp_data
              if(ios1.ne.0) prcp_data = -9999.0
              if(prcp_data.ge.0) then
                 pcp(i) = prcp_data*25.4/3600/24 !convert from inches/day to mm/s
              else
                 pcp(i) = LIS_rc%udef
              endif
              write(LIS_logunit,*) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, pcp(i) 
              readflag = .false.
              exit
           endif
        enddo
        close(ftn)
     endif
  enddo

!  call normalize_stnwts(pcp,snotel_struc(n)%nstns,&
!       LIS_rc%lnc(n)*LIS_rc%lnr(n),snotel_struc(n)%undef,snotel_struc(n)%stnwt)

!  call interp_stndata(snotel_struc(n)%stnwt,snotel_struc(n)%undef,pcp,varfield(:),&
!       LIS_rc%lnc(n)*LIS_rc%lnr(n),snotel_struc(n)%nstns)
 
  do i=1,snotel_struc(n)%nstns

      call latlon_to_ij(LIS_domain(n)%lisproj, snotel_struc(n)%stnlat(i),& 
           snotel_struc(n)%stnlon(i),col,row)
        stn_col = nint(col)
        stn_row = nint(row)

      if(stn_col.gt.0.and.stn_row.gt.0.and. &
             stn_col.le.LIS_rc%lnc(n).and.stn_row.le.LIS_rc%lnr(n)) then
          if(pcp(i).ne.snotel_struc(n)%undef) then
              varfield1(stn_col, stn_row) = varfield1(stn_col, stn_row) + pcp(i)
              npcp(stn_col, stn_row) = npcp(stn_col, stn_row) + 1 !# of stations with pcp obs in local domain
          endif
      endif

  enddo

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(npcp(c,r).gt.0) then
           varfield1(c,r) = varfield1(c,r)/npcp(c,r) !assign station average to grid cell
        else
           varfield1(c,r) = snotel_struc(n)%undef
        endif
     enddo
  enddo

 
!  count1 = 0 
!  do r=1,LIS_rc%lnr(n)
!     do c=1,LIS_rc%lnc(n)
!        varfield1(c,r) = varfield(c+count1)
!     enddo
!     count1 = count1 + LIS_rc%lnc(n)
!  enddo

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              snotel_struc(n)%metdata1(1,LIS_domain(n)%gindex(c,r)) =&
                   varfield1(c,r)
           elseif(order.eq.2) then 
              snotel_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) =&
                   varfield1(c,r)
           endif
        endif
     enddo
  enddo

 end subroutine read_snotel

subroutine create_SNOTEL_filename(odir, stateid, stnid, yr,snotelname)

  use LIS_String_Utility
  implicit none

! !ARGUMRENTS:
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: stateid
  character(len=*), intent(in)  :: stnid
  integer,          intent(in)  :: yr
  character(len=*), intent(out) :: snotelname
! !DESCRIPTION:
!
! This routine creates a filename for the SNOTEL station
!
!  The arguments are:
!  \begin{description}
!   \item[stnid] Station ID
!  \end{description}
!EOP
  character*4             :: fyr

  write(fyr, '(i4.4)' ) yr

  snotelname = trim(odir)//'/'//trim(stateid)//'/'//LIS_StrLowCase(trim(stnid))//'_'&
       //trim(fyr)//'.tab'

end subroutine create_SNOTEL_filename

