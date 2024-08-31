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
! !ROUTINE: restartFromEnsemble
! \label{restartFromEnsemble}
! 
! !DESCRIPTION: 
!  This program generates a restart file for a LIS run that uses 
!  a single ensemble member, from a restart file generated
!  with ensembles. An example would be to generate a restart 
!  file for a single member coupled WRF-LIS run, but the spinup
!  used a multi-member ensemble simulation. 
!  
! !REMARKS:
!  Note that this routines assumes that 
!   * the fields do not contain any integer fields \newline
!  
! !REVISION HISTORY:
! 16 Jul 2010 : Sujay Kumar , Initial Specification
! 
! !ROUINE: 
program restartFromEnsemble

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
     write(*, *)"restartFromEnsemble lsm input_file output_file nfields nens"
     write(*, *)"Usage example: "
     write(*, *)"restartFromEnsemble 1 noah12.rst noah1.rst 18 12"
     write(*, *)" where 18 is the number "
     write(*, *)" of fields and 12 is the number of ensembles."
     write(*, *)" Hint: nfields =18 for noah, 102 for CLM2, 25 for CLSM, 8 for Mosaic"
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
  open(41,file=trim(outputf),form='unformatted')
  read(40) gnc, gnr, glbntiles
  write(41) gnc, gnr, glbntiles/nens
  print*,'dims ',gnc,gnr,glbntiles/nens
  allocate(temp1(glbntiles))
  allocate(temp2(glbntiles/nens))
  
  if(lsm.ne.2) then 
     do j=1,nfields
        print*,'reading field ',j
        read(40) temp1
        do k=1,glbntiles/nens
           temp2(k) = 0 
           do m=1,nens
              temp2(k) = temp2(k) + temp1((k-1)*nens+m)
           enddo
           temp2(k) = temp2(k)/nens
        enddo
        write(41) temp2 
     enddo
  elseif(lsm.eq.2) then 
     print*, 'not tested for lsm 2 '
     stop
     allocate(itemp1(glbntiles/nens))
     allocate(itemp2(glbntiles))
     do j=1,nfields
        print*,'reading field ',j
        if(j.ge.3) then 
           read(40) temp1
           do k=1,glbntiles/nens
              temp2(k) = 0 
              do m=1,nens
                 temp2(k) = temp2(k) + temp1((k-1)*nens+m)
              enddo
              temp2(k) = temp2(k)/nens
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

end program restartFromEnsemble
