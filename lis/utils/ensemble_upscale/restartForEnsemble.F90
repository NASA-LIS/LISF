!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! !ROUTINE: restartForEnsemble
! \label{restartForEnsemble}
! 
! !DESCRIPTION: 
!  This program generates a restart file for a LIS run that uses 
!  ensembles greater than 1 per tile, from a restart file generated
!  without ensembles. An example would be to generate a restart 
!  file for an EnKF run (that uses ensembles), but the spinup
!  using only a single ensemble per tile. 
!  
! !REMARKS:
!  Note that this routines assumes that 
!   * the fields do not contain any integer fields \newline
!  
! !REVISION HISTORY:
! 29 Jan 2007 : Sujay Kumar , Initial Specification
! 
! !ROUINE: 
program restartForEnsemble

! need : filename, number of fields, 

  integer       :: i , j, k, m
  character*80  :: inputf, outputf
  integer       :: nens, nfields, ntiles
  integer       :: lsm
  character*80  :: ctmp
  integer       :: gnc,gnr,glbntiles
  real, allocatable :: temp1(:),temp2(:)
  integer, allocatable :: itemp1(:),itemp2(:)

  i = iargc()
  if(i.ne.5) then 
     write(*,*) "Usage:"
     write(*, *)"restartForEnsemble lsm input_file output_file nfields nens"
     write(*, *)"Usage example: "
     write(*, *)"restartForEnsemble 1 noah.rst noah12.rst 18 12"
     write(*, *)" where 18 is the number "
     write(*, *)" of fields and 12 is the number of ensembles."
     write(*, *)" Hint: nfields =18 for noah 2.7.1, 22 for noah 3.1 23 for noah 3.2, 102 for CLM2, 25 for CLSM, 8 for Mosaic"
     stop
  end if
  
  call getarg(1,ctmp)
  read(ctmp,*) lsm
  call getarg(2,inputf)
  call getarg(3,outputf)
  call getarg(4,ctmp)
  read(ctmp,*) nfields
  call getarg(5,ctmp)
  read(ctmp,*) nens
  
  print*, lsm,inputf,outputf,nfields,nens

  open(40,file=trim(inputf),form='unformatted',status='old')
  open(41,file=trim(outputf),form='unformatted',status='new')
  read(40) gnc, gnr, glbntiles
  write(41) gnc, gnr, glbntiles*nens
  print*,'dims ',gnc,gnr,glbntiles
  allocate(temp1(glbntiles))
  allocate(temp2(glbntiles*nens))
  
  if(lsm.ne.2) then 
     do j=1,nfields
        print*,'reading field ',j
        read(40) temp1
        do k=1,glbntiles
           do m=1,nens
              temp2((k-1)*nens+m) = temp1(k)
           enddo
        enddo
        write(41) temp2 
     enddo
  elseif(lsm.eq.2) then 
     allocate(itemp1(glbntiles))
     allocate(itemp2(glbntiles*nens))
     do j=1,nfields
        print*,'reading field ',j
        if(j.ge.3) then 
           read(40) temp1
           do k=1,glbntiles
              do m=1,nens
                 temp2((k-1)*nens+m) = temp1(k)
              enddo
           enddo
           write(41) temp2 
        else
           read(40) itemp1
           do k=1,glbntiles
              do m=1,nens
                 itemp2((k-1)*nens+m) = itemp1(k)
              enddo
           enddo
           write(41) itemp2 
        endif
     enddo
     deallocate(itemp1)
     deallocate(itemp2)
  endif
  deallocate(temp1)
  deallocate(temp2)
  print*,'program finished successfully'

end program restartForEnsemble
