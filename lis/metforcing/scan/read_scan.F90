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
! !ROUTINE: read_scan
! \label{read_scan}
! 
! !REVISION HISTORY:
!
!  13 Apr 2007: Bailing Li; Initial Specification 
!
! !INTERFACE:      
subroutine read_scan(n,ftn,findex,order)
! !USES:
  use LIS_logMod, only         : LIS_logunit
  use LIS_metforcingMod,  only : LIS_forc
  use LIS_coreMod, only        : LIS_rc,LIS_domain
  use LIS_constantsMod,   only : LIS_CONST_PATH_LEN
  use scan_forcingMod,    only : scan_struc

  implicit none
! !ARGUMENTS:
  integer, intent(in)           :: n
  integer, intent(in)           :: ftn
  integer, intent(in)           :: findex
  integer, intent(in)           :: order
!
! !DESCRIPTION:
!  For the given time, reads parameters from
!  the correct SCAN station data (ASCII), transforms into LIS forcing 
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
!    unit number for the SCAN station data 
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

  integer :: c,r,count1,i
  real :: pcp(scan_struc(n)%nstns),tmppcp,dum
  real :: varfield(LIS_rc%lnc(n)*LIS_rc%lnr(n))
  real :: varfield1(LIS_rc%lnc(n),LIS_rc%lnr(n))
  character(len=LIS_CONST_PATH_LEN) :: scan_filename
  character(len=500) :: line
  integer :: yr,num,hr,mon,day,mint,sec
  logical :: file_exists

  pcp = scan_struc(n)%undef

  do i=1,scan_struc(n)%nstns
!generate the scan filename
     write(scan_filename,'(i4,a,i4,i2.2,a)')scan_struc(n)%stnid(i),'_',&
          LIS_rc%yr,LIS_rc%mo,'.txt'
     scan_filename = trim(scan_struc(n)%scandir)//'/'//trim(scan_filename)
     write(LIS_logunit,*) 'Reading SCAN file ',trim(scan_filename)
     inquire(file=scan_filename,exist=file_exists)
     if(file_exists) then 
        open(ftn,file=scan_filename,form='formatted',status='old')

! segment to skip the header
        do 
           read(ftn,'(a)')line
           if (line(1:4) .eq. 'Date') exit
        enddo
        
!  actual data reading section. 
        do
           read(ftn,'(a)')line
           !The empty space line signals the end of the first data section
           if (line(1:1) .lt. '0' .or. line(1:1) .gt. '9') exit 
           read(line,40)yr,mon,day,hr,mint,sec,num,dum,tmppcp
           if ( LIS_rc%da== day .and. LIS_rc%hr == hr ) then
              pcp(i) = tmppcp*25.4/3600                 !rainfall rate mm/s
              exit
           endif
        enddo
        close(ftn)
     endif
  enddo

40  format(7i2,2x,f5.2,2x,f5.2)

  call normalize_stnwts(pcp,scan_struc(n)%nstns,&
       LIS_rc%lnc(n)*LIS_rc%lnr(n),scan_struc(n)%undef,scan_struc(n)%stnwt)

  call interp_stndata(scan_struc(n)%stnwt,scan_struc(n)%undef,pcp,varfield(:),&
       LIS_rc%lnc(n)*LIS_rc%lnr(n),scan_struc(n)%nstns)
  
  count1 = 0 
  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        varfield1(c,r) = varfield(c+count1)
     enddo
     count1 = count1 + LIS_rc%lnc(n)
  enddo

  do r=1,LIS_rc%lnr(n)
     do c=1,LIS_rc%lnc(n)
        if(LIS_domain(n)%gindex(c,r).ne.-1) then 
           if(order.eq.1) then 
              scan_struc(n)%metdata1(1,LIS_domain(n)%gindex(c,r)) =&
                   varfield1(c,r)
           elseif(order.eq.2) then 
              scan_struc(n)%metdata2(1,LIS_domain(n)%gindex(c,r)) =&
                   varfield1(c,r)
           endif
        endif
     enddo
  enddo
 end subroutine read_scan
