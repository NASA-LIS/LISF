!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_coop
! \label{read_coop}
! 
! !REVISION HISTORY:
!
!  13 Jul 2011: Yuqiong Liu; Initial Specification, based on LVT COOP obs reader by S. Kumar 
!
! !INTERFACE:      
subroutine read_coop(n,ftn,findex,order) 
! !USES:
  use LIS_logMod, only         : LIS_logunit, LIS_endrun
  use LIS_coreMod, only        : LIS_rc,LIS_domain
  use LIS_metforcingMod,  only : LIS_forc
  use coop_forcingMod,    only : coop_struc
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
!  the correct COOP station data (ASCII), transforms into LIS forcing 
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
!    unit number for the COOP station data 
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

  integer :: c,r,count1,i,j, jj
  real :: pcp(coop_struc(n)%nstns)
!  real :: varfield(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer :: npcp(LIS_rc%lnc(n),LIS_rc%lnr(n))
  character*80 :: coop_filename
!  character(len=500) :: line
  integer :: yr,mo
  logical :: file_exists, readflag
  integer :: ios
 ! real    :: swe_data, prcp_data
 ! character*10  :: datestring
  integer                :: dorig, stnid, wban,div
  character*30           :: stnname
  character*4            :: mettype
  character*2            :: units
  character*1            :: flag1(31),flag2(31)
  integer                :: da(31),hr(31)
  integer                :: prcp(31)
  integer                :: nstns
  integer                :: stn_col, stn_row
  real                   :: col,row


  pcp = coop_struc(n)%undef
  varfield1 = 0
  npcp = 0

  do i=1,coop_struc(n)%nstates
!generate the coop filename
     yr = LIS_rc%yr
     mo = LIS_rc%mo
    
     call create_COOP_filename(trim(coop_struc(n)%coopdir),coop_struc(n)%statename(i), &
        yr, mo, coop_filename)

     inquire(file=coop_filename,exist=file_exists)
     if(file_exists) then 
        write(LIS_logunit,*) 'Reading COOP file ',coop_filename
        open(ftn,file=coop_filename,form='formatted',status='old')
        read(ftn,*)
        read(ftn,*) ! skip the two header line

        readflag = .true.
        nstns = 0
        do while (readflag)
!        do while (nstns .lt. coop_struc(n)%nstns)
           read(ftn,fmt=300, iostat=ios) dorig,stnid,wban,stnname,div, &
                mettype, units,yr,mo,(da(j),hr(j), prcp(j),flag1(j),flag2(j),j=1,31)
           do j=1,coop_struc(n)%nstns
               !write(LIS_logunit,*) coop_struc(n)%stnid(j), stnid
              if(coop_struc(n)%stnid(j) .eq. stnid) then
                 jj = LIS_rc%da 
                 if ((trim(flag1(jj)).eq.'').and. (trim(flag2(jj)).eq.'0') .and. &
                     (prcp(jj).ge.0 ) ) then

                     pcp(j) = prcp(jj)

                     if(trim(adjustl(units)).eq.'I') then
                         pcp(j)= pcp(j)*25.4/24/3600 !convert to mm/s.
                     elseif(trim(adjustl(units)).eq.'TI') then
                         pcp(j)= pcp(j)*2.54/24/3600
                     elseif(trim(adjustl(units)).eq.'HI') then
                         pcp(j)= pcp(j)*0.254/24/3600
                     else
                         write(LIS_logunit,*) 'precip Units not supported at ',&
                           stnname, 'at ',yr,mo
                         write(LIS_logunit,*) 'Units are ',trim(units)
                         call LIS_endrun
                     endif
                     if (pcp(j).gt. 100.0/24/3600) pcp(j) = coop_struc(n)%undef
                 endif 
                 nstns = nstns + 1
                 write(LIS_logunit,*) LIS_rc%yr, LIS_rc%mo, LIS_rc%da, pcp(j)
 
                 exit
              endif
           enddo
           if(ios.ne.0 .or. nstns.eq.coop_struc(n)%nstns) then
              readflag = .false.
           endif

        enddo
!        if (nstns .ne. coop_struc(n)%nstns) then
!            write(LIS_logunit,*) '# of coop stations not consistent!', coop_struc(n)%nstns, nstns
!            call LIS_endrun
!        endif 

        close(ftn)
     endif
  enddo

300 format(I4,1X,I6,1X,I5,1X,A30,1X,I2,1X,A4,1X,A2,1X,I4.4,I2.2,1X,&
         31(I2.2,I2.2,1X,I6,1X,A1,1X,A1,1X))

!  call normalize_stnwts(pcp,coop_struc(n)%nstns,&
!       LIS_rc%lnc(n)*LIS_rc%lnr(n),coop_struc(n)%undef,coop_struc(n)%stnwt)

!  call interp_stndata(coop_struc(n)%stnwt,coop_struc(n)%undef,pcp,varfield(:),&
!       LIS_rc%lnc(n)*LIS_rc%lnr(n),coop_struc(n)%nstns)
  
  do i=1,coop_struc(n)%nstns

      call latlon_to_ij(LIS_domain(n)%lisproj, coop_struc(n)%stnlat(i),&
           coop_struc(n)%stnlon(i),col,row)
        stn_col = nint(col)
        stn_row = nint(row)

      if(stn_col.gt.0.and.stn_row.gt.0.and. &
             stn_col.le.LIS_rc%lnc(n).and.stn_row.le.LIS_rc%lnr(n)) then
          if(pcp(i).ne.coop_struc(n)%undef) then
              varfield1(stn_col, stn_row) = varfield1(stn_col, stn_row) + pcp(i)
              npcp(stn_col, stn_row) = npcp(stn_col, stn_row) + 1
          endif
      endif

  enddo

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(npcp(c,r).gt.0) then
           varfield1(c,r) = varfield1(c,r)/npcp(c,r) !assign station average to grid cell
        else
           varfield1(c,r) = coop_struc(n)%undef
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
              coop_struc(n)%metdata1(1,LIS_domain(n)%gindex(c,r)) =&
                   varfield1(c,r)
           elseif(order.eq.2) then
              coop_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) =&
                   varfield1(c,r)
           endif
           !write(LIS_logunit,*) 'pcp=',pcp, '  suppl=',varfield1(c,r)
        endif
     enddo
  enddo

 end subroutine read_coop

subroutine create_COOP_filename(odir, stateid, yr, mo, coopname)

  use LIS_String_Utility
  implicit none

! !ARGUMRENTS:
  character(len=*), intent(in)  :: odir
  character(len=*), intent(in)  :: stateid
  integer,          intent(in)  :: yr
  integer,          intent(in)  :: mo 
  character(len=*), intent(out) :: coopname
! !DESCRIPTION:
!
! This routine creates a filename for the COOP station
!
!  The arguments are:
!  \begin{description}
!   \item[stnid] Station ID
!  \end{description}
!EOP

  character*4             :: fyr
  character*2             :: fmo

  write(fyr, '(i4.4)' ) yr
  write(fmo, '(i2.2)' ) mo

  coopname = trim(odir)//'/'//trim(stateid)//'/'//trim(stateid)//'_'&
       //trim(fyr)//trim(fmo)//'_PRCPdat.txt'

end subroutine create_COOP_filename

