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
! !ROUTINE: generateUEinitFile
! \label{generateUEinitFile}
! 
! !DESCRIPTION: 
!  This program generates an initialization file for UE algorithms
!  (MCMC, DEMC) from an optimization algorithm output. 
!  
! !REMARKS:
!  Note that this routines assumes that 
!   * the fields do not contain any integer fields \newline
!  
! !REVISION HISTORY:
! 3 Aug 2010 : Sujay Kumar, Ken Harrison , Initial Specification
! 
! !ROUINE: 
program generateUEinitFile

  integer       :: i,k,m,c,r
  character*80  :: opt_outfile, ue_initfile
  integer       :: nens_ue, glbngrid, gnc, gnr, nparam
  integer       :: opt, iterno, count, obj_func
  character*80  :: ctmp
  real, allocatable :: temp1(:,:),fitness(:,:),temp2(:),temp3(:)

  i = iargc()
  if(i.ne.9) then 
     write(*,*) "Usage:"
     write(*, *)"generateUEinitFile <opt> <opt_outfile> <gnc> <gnr> <glbngrid> <nparam> <obj_func> <ue_initfile> <nens_ue>"
     write(*, *)"Usage example: "
     write(*, *)"generateUEinitFile 2 GA100.1gd4r 5 5 23 4 1 DEMC_init.dat 1000"
     write(*, *)" where 1000 is the number "
     write(*, *)" of ensembles in DEMC"
     write(*, *)" Hint: <opt> = 1 for LM, 2 for GA, 3 for SCEUA"
     stop
  end if
  
  call getarg(1,ctmp)
  read(ctmp,*) opt
  call getarg(2,opt_outfile)
  call getarg(3,ctmp)
  read(ctmp,*) gnc
  call getarg(4,ctmp)
  read(ctmp,*) gnr
  call getarg(5,ctmp)
  read(ctmp,*) glbngrid
  call getarg(6,ctmp)
  read(ctmp,*) nparam
  call getarg(7,ctmp)
  read(ctmp,*) obj_func
  call getarg(8,ue_initfile)
  call getarg(9,ctmp)
  read(ctmp,*) nens_ue
  
  open(40,file=trim(opt_outfile),form='unformatted',status='old')
  open(41,file=trim(ue_initfile),form='unformatted')

  allocate(temp1(gnc,gnr))
  allocate(fitness(gnc,gnr))
  allocate(temp2(glbngrid))
  allocate(temp3(glbngrid*nens_ue))

  if(opt.eq.2) then 
     read(40) fitness  !fitness
     read(40) temp1  !avgfitness

     do i=1,nparam
        read(40) temp1 !parameter
        count = 0
        do c=1,gnc
           do r=1,gnr              
              if(temp1(c,r).ne.-9999.0) then 
                 count = count + 1
                 temp2(count) = temp1(c,r)
              endif
           enddo
        enddo
        if(count.ne.glbngrid) then 
           print*, 'Error: the domain sizes do not match '
           print*, 'gnc,gnr, glbngrid,count = ',gnc,gnr,glbngrid, count
           stop
        endif
        do k=1,glbngrid
           do m=1,nens_ue
              temp3((k-1)*nens_ue+m) = temp2(k)
              print*,               temp3((k-1)*nens_ue+m), (k-1)*nens_ue+m
           enddo
        enddo
        write(41) temp3
     enddo
!write sigma, if objective function is least squares. 
!here we assume that the fitness =1/standard_deviation
     if(obj_func.eq.1) then 
        count = 0 
        do c=1,gnc
           do r=1,gnr
              if(fitness(c,r).ne.-9999.0) then 
                 count = count + 1
                 temp2(count) = 1/fitness(c,r)
                 print*, 'sigma ',temp2(count)
              endif
           enddo
        enddo
        do k=1,glbngrid
           do m=1,nens_ue
              temp3((k-1)*nens_ue+m) = temp2(k)
           enddo
        enddo
        write(41) temp3
     endif

  endif
  deallocate(temp1)
  deallocate(temp2)
  deallocate(temp3)
  deallocate(fitness)

  close(41)
  print*,'program finished successfully'

end program generateUEinitFile
