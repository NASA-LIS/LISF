!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center Land Information System (LIS) v7.2
!
! Copyright (c) 2015 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
!
! !ROUTINE: read_FASSTsingle
! \label{read_FASSTsingle}
! 
! !REVISION HISTORY:
! 13 Apr 2007: Bailing Li, Initial Specification 
! 06 Oct 2010: David Mocko, Updated for FASST single-point test case
!
! !INTERFACE:
subroutine read_FASSTsingle(n,ftn,findex,metdata,itime)
! !USES:
  use LIS_logMod, only            : LIS_logunit,LIS_endrun
  use LIS_coreMod, only           : LIS_rc,LIS_domain
  use LIS_timeMgrMod, only        : LIS_date2time,LIS_doy2date
  use FASSTsingle_forcingMod, only : FASSTsingle_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in) :: n
  integer, intent(in) :: ftn
  integer, intent(in) :: itime
  integer, intent(in) :: findex
  real(kind=8) :: metdata(LIS_rc%met_nf(findex),LIS_rc%ngrid(n))
!
! !DESCRIPTION:
!  For the given time, reads parameters from the correct FASSTsingle
!  station data (ASCII), transforms into LIS forcing parameters, and
!  interpolates to the LIS domain.
!
!  The arguments are:
!  \begin{description}
!  \item[n]
!    index of the nest
!  \item[ftn]
!    unit number for the FASSTsingle station data
!  \item[metdata]
!    supplemental data (output as read from station file)
!  \item[itime]
!    index for reading current time or reading a future time
!  \end{description}
!
!  The routines invoked are:
!  \begin{description}
!  \item[normalize\_stnwts](\ref{normalize_stnwts}) \newline
!    renormalizes the station weights accounting for
!    missing data
!  \item[interp\_stndata](\ref{interp_stndata}) \newline
!    spatially interpolates the station data onto the LIS grid.
!  \end{description}
!EOP

  integer :: i,c,r,f,count1
  real    :: pcp(FASSTsingle_struc(n)%nstns)
  real    :: psurf(FASSTsingle_struc(n)%nstns)
  real    :: tair(FASSTsingle_struc(n)%nstns)
  real    :: qair(FASSTsingle_struc(n)%nstns)
  real    :: swdown(FASSTsingle_struc(n)%nstns)
  real    :: lwdown(FASSTsingle_struc(n)%nstns)
  real    :: u(FASSTsingle_struc(n)%nstns)
  real    :: v(FASSTsingle_struc(n)%nstns)
  real(kind=8) :: varfield(35,LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real(kind=8) :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  real*8  :: listime,fassttime
  real    :: lisgmt,fasstgmt
  integer :: lisdoy,fasstdoy
  integer :: fasstyr,fasstmon,fasstday,fassthr,fasstmin,fasstsec
  logical :: file_exists
  character*80       :: FASSTsingle_filename
  character(len=500) :: line

  !      write(LIS_logunit,*) 'starting read_FASSTsingle'
  do i = 1,FASSTsingle_struc(n)%nstns
     ! Generate the FASSTsingle filename and see if it exists
     FASSTsingle_filename = trim(FASSTsingle_struc(n)%FASSTsinglefile)
     write(LIS_logunit,*) 'Reading FASSTsingle file: ',            &
          trim(FASSTsingle_filename)
     inquire(file=FASSTsingle_filename,exist=file_exists)
     if (file_exists) then
        !            write(LIS_logunit,*) 'File is open!!'
        open(ftn,file=FASSTsingle_filename,form='formatted',       &
             status='old')
        ! segment to skip the header
        do f = 1,3
           read(ftn,'(a)') line
           !               write(LIS_logunit,*) line(1:30)
        enddo
        !            write(LIS_logunit,*) 'Starting data read'
        ! Actual data reading section
        do
           read(ftn,*) (FASSTsingle_struc(n)%metm(f),f=1,32),      &
                FASSTsingle_struc(n)%metm(34),             &
                FASSTsingle_struc(n)%metm(35)
           fasstyr = int(FASSTsingle_struc(n)%metm(1))
           fasstdoy = int(FASSTsingle_struc(n)%metm(2))
           fassthr = int(FASSTsingle_struc(n)%metm(3))
           fasstmin = int(FASSTsingle_struc(n)%metm(4))
           call LIS_doy2date(fasstyr,fasstdoy,fasstmon,fasstday)
           !               write(LIS_logunit,*) fasstyr,fasstmon,fasstday,fassthr,fasstmin
           fasstsec = 0
           call LIS_date2time(fassttime,fasstdoy,fasstgmt,fasstyr,fasstmon,  &
                fasstday,fassthr,fasstmin,fasstsec)
           call LIS_date2time(listime,lisdoy,lisgmt,LIS_rc%yr,     &
                LIS_rc%mo,LIS_rc%da,LIS_rc%hr,       &
                LIS_rc%mn,LIS_rc%ss)
           !               write(LIS_logunit,*) 'read_time: ',fassttime,fasstyr,fasstmon,&
           !                                                   fasstday,fassthr,fasstmin
           !               write(LIS_logunit,*) 'LIS_time: ',listime,LIS_rc%yr,    &
           !                               LIS_rc%mo,LIS_rc%da,LIS_rc%hr,LIS_rc%mn
           if ((fassttime.ge.listime).and.(itime.eq.1)) exit
           if ((fassttime.gt.listime).and.(itime.eq.2)) exit
        enddo
        close(ftn)
     else
        write(LIS_logunit,*) 'Filename not found; stopping program'
        call LIS_endrun
     endif
  enddo

  do f = 1,35
     call normalize_stnwts_fasst(FASSTsingle_struc(n)%metm(f),     &
          FASSTsingle_struc(n)%nstns,             &
          LIS_rc%lnc(n)*LIS_rc%lnr(n),            &
          FASSTsingle_struc(n)%undef,             &
          FASSTsingle_struc(n)%stnwt)
     call interp_stndata_fasst(FASSTsingle_struc(n)%stnwt,         &
          FASSTsingle_struc(n)%undef,               &
          FASSTsingle_struc(n)%metm(f),             &
          varfield(f,:),LIS_rc%lnc(n)*LIS_rc%lnr(n),&
          FASSTsingle_struc(n)%nstns)
  enddo

  do f = 1,35
     count1 = 0
     do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
           varfield1(c,r) = varfield(f,c+count1)
        enddo
        count1 = count1 + LIS_rc%lnc(n)
     enddo

     do r = 1,LIS_rc%lnr(n)
        do c = 1,LIS_rc%lnc(n)
           if (LIS_domain(n)%gindex(c,r).ne.-1) then
              metdata(f,LIS_domain(n)%gindex(c,r)) = varfield1(c,r)
           endif
        enddo
     enddo
  enddo

end subroutine read_FASSTsingle

!-----------------------------------------------------------------------

subroutine normalize_stnwts_fasst(stndata,nstns,npts,undef,wt)
  implicit none
! !ARGUMENTS: 
  integer    :: nstns
  integer    :: npts
  real(kind=8) :: stndata(nstns)
  real       :: wt(npts,nstns)
  real       :: undef
!
! !DESCRIPTION: 
!  This subroutine normalizes the interpolation weights to convert 
!  station data into a gridded set. 
!
!  The arguments are: 
!  \begin{description}
!    \item[nstns]
!     number of stations used in interpolation
!    \item[npts]
!     integer maximum number of coordinates
!    \item[stndata]
!     input data (station observations)
!    \item[undef]
!     undefined value used in observations
!    \item[wt] 
!     interpolation weights
!   \end{description}
!
!EOP
  integer :: i,j
  real :: norm(npts)
  norm = 0 
  do i=1,npts
     do j=1,nstns
        if(stndata(j).ne.undef) then
           norm(i) = norm(i)+wt(i,j)
        endif
     enddo
  enddo
  do i=1,npts
     do j=1,nstns
        if(norm(i).eq.0) norm(i) = 1
        wt(i,j) = wt(i,j)/norm(i)
     enddo
  enddo
end subroutine normalize_stnwts_fasst

!-----------------------------------------------------------------------

subroutine interp_stndata_fasst(wt,undef,vari,varo,npts,nstns)

  use LIS_coreMod,  only : LIS_rc

  implicit none
! !ARGUMENTS: 

  integer :: nstns,npts
  real    :: undef
  real    :: wt(npts,nstns)
  real(kind=8) :: vari(nstns)
  real(kind=8) :: varo(npts)  

! !DESCRIPTION: 
!  This subroutine interpolates station data into a gridded set. 
!  Currently works only on a lat/lon grid
!
!  The arguments are:
!  \begin{description}
!  \item[nstns]
!    number of stations
!  \item[npts]
!    number of points in the output field
!  \item[vari]
!    input variable (array of observations from the stations)
!  \item[varo]
!    output variable (gridded data on the output field)
!  \item[wt]
!    interpolation weights 
!  \end{description}
!EOP

  integer:: i,j
  integer :: nundefs

  varo = 0.0
  nundefs = 0

  do i=1,npts
     do j=1,nstns
        if(vari(j).ne.undef) then 
           varo(i) = varo(i) + wt(i,j)*vari(j)
        else
           nundefs = nundefs + 1
        endif
        !        if(i==1.and.j==1) print*, i,j,varo(i),wt(i,j),vari(j),npts
     enddo
     if(nundefs .eq. nstns) then !all stns have undef values. 
        varo(i) = LIS_rc%udef
     endif
     nundefs = 0 
     !     print*, i,varo(i)
  enddo
end subroutine interp_stndata_fasst
