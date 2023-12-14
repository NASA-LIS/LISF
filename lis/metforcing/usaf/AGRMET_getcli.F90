!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: AGRMET_getcli
! \label{AGRMET_getcli}
!
! !REVISION HISTORY:
!
!    04 Dec 97  Initial version............................Capt Andrus/AFWA/DNXM
!    31 Mar 99  Changed data retrieval from files to new UNIX based system.
!               Increased array and loop sizes to full hemisphere..Mr Moore/DNXM
!    08 Oct 99  Ported to IBM SP-2, incorporated FORTRAN 90....Capt Hidalgo/DNXM
!    21 Feb 01  Reformatted diagnostic prints.  
!               Added intent attribute to arguments................Mr Gayno/DNXM
!    10 Jun 02  Changed all references from rtneph to cdfs2........Mr Gayno/DNXM
!    01 Nov 05  Adopted in LIS, split out read_pcpclimodata into
!               separate routine........................Mr Sujay Kumar/NASA/GSFC
!    11 May 07  Made tempdata array allocatable.
!               Updated putget argument list....................Mr Lewiston/A8TM
!    29 May 15  Updated read in of files................Ryan Ruhge/SEMS/16WS/WXE
!
! !INTERFACE:
subroutine AGRMET_getcli(n, filename,rtn,clidat)
! !USES: 
  use LIS_coreMod, only       : LIS_rc, LIS_ews_halo_ind, LIS_ewe_halo_ind,&
                                LIS_nss_halo_ind, LIS_nse_halo_ind, LIS_localPet
  use LIS_logMod, only        : LIS_logunit, LIS_getNextUnitNumber, &
       LIS_releaseUnitNumber

  implicit none
! !ARGUMENTS: 
  integer, intent(in)     :: n 
  character*100           :: filename
  real,    intent(out)    :: clidat(LIS_rc%lnc(n),LIS_rc%lnr(n))
  integer                 :: rtn
! 
! !DESCRIPTION: 
! 
!   This routine retrieves monthly climo data from file. The same routine is 
!   called to retrieve the three types of precip climo data. 
!    climatological rtneph percent cloud cover \newline
!    3 hour precip climo \newline
!    precip-per-precip day amount \newline
!   
!    The arguments and variables are: 
!  \begin{description}
!   \item[n] 
!     index of the nest
!   \item[filename]
!     name of the climo file
!   \item[rtn]
!     flag to indicate that the data should be read as integers
!   \item[clidat]
!     output climatological values
!   \item[exists]
!     a logical that indicates whether or not a file exists
!   \item[tempdata]
!     a temp array used to read in climo files that are stored
!     as integers
!  \end{description}
!EOP

  logical               :: exists
  integer               :: ftn
  real                  :: data_in(LIS_rc%gnc(n), LIS_rc%gnr(n))
  integer               :: istat

  inquire( file = trim(filename), exist = exists)
  if ( exists )  then 
     write(LIS_logunit,*) '[INFO] READING ', trim(filename)
     
!     ------------------------------------------------------------------
!     the rtneph climo files are real valued files, while the precip
!     per precip day files and precip files are integer file.
!     ------------------------------------------------------------------
     if ( rtn.eq.1) then

        ftn= LIS_getNextUnitNumber()
        open(ftn, file=trim(filename), access='direct',&
             status='old', form="unformatted", recl=LIS_rc%gnr(n)*LIS_rc%gnc(n)*4)
     
        read(ftn, rec=1, iostat=istat) data_in
        clidat(:,:) =  data_in(&
           LIS_ews_halo_ind(n,LIS_localPet+1):&
           LIS_ewe_halo_ind(n,LIS_localPet+1), &
           LIS_nss_halo_ind(n,LIS_localPet+1): &
           LIS_nse_halo_ind(n,LIS_localPet+1))
        call LIS_releaseUnitNumber(ftn)
     else
!     ------------------------------------------------------------------
!     otherwise, use putget_int and then convert to real values
!     to get real-values
!     RR - the above comment seems out of date and does not match the code
!     ------------------------------------------------------------------
       print*, 'Get CLI being called with rtn ',rtn
       print*, 'Stopping ..'
       stop
!        allocate (tempdata (agrmet_struc(n)%imax, agrmet_struc(n)%jmax))
!        call LIS_putget( tempdata, 'r', filename, routine_name, &
!                     agrmet_struc(n)%imax, agrmet_struc(n)%jmax )
!        clidat = float(tempdata)
!        deallocate (tempdata)
     endif
     
  endif
end subroutine AGRMET_getcli
