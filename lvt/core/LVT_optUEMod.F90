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
! !MODULE: LVT_optUEMod
! \label(LVT_optUEMod)
!
! !INTERFACE:
module LVT_optUEMod
! 
! !USES:   
  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   The code in this file provides interfaces to manages data extraction 
!   from the optimization (Genetic Algorithm) output
! 
!  NOTE : Currently the LIS domain and LVT domain are assumed to be 
!  identical
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  02 Oct 2008    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_optUEInit          !initialize specified domains
  public :: LVT_readoptUEdata
  public :: LVT_computeoptUEstats
!EOP
  public :: LVT_optuectl

  type, public ::  optuectl
     integer                :: optueAlg
     character*10           :: algname
     character*100          :: decspaceAttribsFile
     integer                :: maxIter
     integer                :: nparam
     integer                :: nparam_total
     integer                :: computeTS
     integer                :: ntslocs
     integer                :: tsspecstyle
     character*100          :: tsspecfile
     character*40,  allocatable :: vname(:)
     real,   allocatable        :: fitness(:,:)
     real,   allocatable        :: avgfitness(:,:)
     character*100, allocatable :: param_name(:)
     real,   allocatable        :: param_val(:,:)
     real,   allocatable        :: param_val_ga(:,:,:)
  end type optuectl

  type, public :: optue_ts_struc
     real         :: tslat1
     real         :: tslon1
!     real         :: tslat2
!     real         :: tslon2
     integer      :: ts_cindex1
     integer      :: ts_rindex1
!     integer      :: ts_cindex2
!     integer      :: ts_rindex3
     integer      :: ts_tindex
     integer,       allocatable :: ftn_ts_loc(:)
     character*40           :: tslocname
     character*100, allocatable :: tslocfile(:)
     integer                :: ftn_ts_fitloc
     character*100          :: tsloc_fitfile
  end type optue_ts_struc

  type(optuectl)                   :: LVT_optuectl
  type(optue_ts_struc), allocatable    :: LVT_optuets(:)

contains
!BOP
! 
! !ROUTINE: LVT_optUEInit
! \label{LVT_optUEInit}
!
! !INTERFACE: 
  subroutine LVT_optUEInit()
! 
! !USES:   
    use ESMF
    use LVT_coreMod
    use LVT_logMod
    use map_utils

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine performs the initialization steps for optimization or 
!  uncertainty estimation data extraction. 
!  It reads the decision space attributes file (from the LIS run), 
!  and other configurable options. The time series file that specifies the 
!  location of the points to be extracted are also read.
!    
!   The locations can be specified
!   in three different formats: (1) using the lat/lon values (2) using the 
!   column/row indices and (3) using the tile indices. A sample file is 
!   shown below: \newline
!   
!   \begin{verbatim}
!    #Number of locations
!    1
!    #Location style (1-lat/lon, 2-col/row, 3-tile)
!    2 
!    #site name 
!    Site1
!    244  236 
!   \end{verbatim}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer :: i,j,k
    integer :: status
    integer :: ftn,ftn2

    real          :: parmax
    real          :: parmin
    integer       :: selectOpt
    real          :: col, row
    integer       :: count
    character*40  :: algname
    character*40  :: vname

!read decision space file, number of generations, 
    call ESMF_ConfigGetAttribute(LVT_config,LVT_rc%odir,&
         label="OptUE output data directory:", rc=status)
    call LVT_verify(status, 'OptUE output data directory: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,algname,&
         label="OptUE algorithm used:", rc=status)
    call LVT_verify(status, 'OptUE algorithm used: not defined')

    if(algname.eq."Genetic algorithm") then 
       LVT_optuectl%algname = 'GA'
       LVT_optuectl%optuealg = 2
    elseif(algname.eq."Monte carlo sampling") then 
       LVT_optuectl%algname = 'MCSIM'
       LVT_optuectl%optuealg = 4
    elseif(algname.eq. & 
         "Random walk markov chain monte carlo") then 
       LVT_optuectl%algname = 'RWMCMC'
       LVT_optuectl%optuealg = 5
    elseif(algname.eq.&
         "Differential evolution monte carlo") then 
       LVT_optuectl%algname = 'DEMC'
       LVT_optuectl%optuealg = 6
    endif

    call ESMF_ConfigGetAttribute(LVT_config,LVT_optuectl%decspaceAttribsFile,&
         label="OptUE decision space attributes file:", rc=status)
    call LVT_verify(status, 'OptUE decision space attributes file: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_optuectl%maxIter,&
         label="OptUE number of iterations:",rc=status)
    call LVT_verify(status, 'OptUE number of iterations: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_optuectl%computeTS,&
         label="OptUE compute time series:",rc=status)
    call LVT_verify(status, 'OptUE compute time series: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_optuectl%tsspecfile,&
         label="OptUE time series location file:",rc=status)
    call LVT_verify(status, 'OptUE time series location file: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_optuectl%nparam_total,&
         label="OptUE total number of parameters:",rc=status)
    call LVT_verify(status, 'OptUE total number of parameters: not defined')

    call ESMF_ConfigGetAttribute(LVT_config,LVT_optuectl%nparam,&
         label="OptUE total number of selected parameters:",rc=status)
    call LVT_verify(status, 'OptUE total number of selected parameters: not defined')

    allocate(LVT_optuectl%vname(LVT_optuectl%nparam))
    

    ftn = LVT_getNextUnitNumber()
    write(LVT_logunit,*) '[INFO] Reading decision space attributes ...', &
         LVT_optuectl%decspaceAttribsFile
    open(ftn,file=(LVT_optuectl%decspaceAttribsFile),status='old')   
    read(ftn,*)
   
    count = 0 
    do i=1,LVT_optuectl%nparam_total
       read(ftn,*) selectOpt, vname, parmin,parmax
       write(LVT_logunit,*) '[INFO] vname ',vname,&
            parmin,parmax
       if(selectOpt.eq.1) then 
          count = count + 1
          LVT_optuectl%vname(count) = vname
       endif
    enddo
    
    call LVT_releaseUnitNumber(ftn)
    write(LVT_logunit,*) '[INFO] Finished reading decision space attributes ..'

    if(LVT_optuectl%computeTS.eq.1) then 
       ftn = LVT_getNextUnitNumber()
       open(ftn,file=(LVT_optuectl%tsspecfile),status='old')
       read(ftn,*) 
       read(ftn,*) LVT_optuectl%ntslocs
       read(ftn,*) 
       read(ftn,*) LVT_optuectl%tsspecstyle
       read(ftn,*) 
       
       allocate(LVT_optuets(LVT_optuectl%ntslocs))
       do k=1,LVT_optuectl%ntslocs
          allocate(LVT_optuets(k)%tslocfile(1))
          allocate(LVT_optuets(k)%ftn_ts_loc(1))
       enddo
       
       write(LVT_logunit,*) '[INFO] Reading location data ....'
       do k=1,LVT_optuectl%ntslocs
          read(ftn,fmt='(a40)') LVT_optuets(k)%tslocname
          write(LVT_logunit,*) '[INFO] Station name: ',LVT_optuets(k)%tslocname
          
          if(LVT_optuectl%tsspecstyle.eq.1) then 
             read(ftn,*)  LVT_optuets(k)%tslat1, LVT_optuets(k)%tslon1
             write(LVT_logunit,*) '[INFO] Station location: ',&
                  LVT_optuets(k)%tslat1, LVT_optuets(k)%tslon1
          elseif(LVT_optuectl%tsspecstyle.eq.2) then 
             read(ftn,*)  LVT_optuets(k)%ts_cindex1, LVT_optuets(k)%ts_rindex1
             write(LVT_logunit,*) '[INFO] Station location: ',&
                  LVT_optuets(k)%ts_cindex1, LVT_optuets(k)%ts_rindex1
          elseif(LVT_optuectl%tsspecstyle.eq.3) then 
             read(ftn,*) LVT_optuets(k)%ts_tindex
             write(LVT_logunit,*) '[INFO] Station location: ',&
                  LVT_optuets(k)%ts_tindex
          endif
       enddo
       call LVT_releaseUnitNumber(ftn)
    endif
       
    if(LVT_optuectl%computeTS.eq.1) then       
       if(LVT_optuectl%tsspecstyle.eq.1) then 
          do i=1,LVT_optuectl%ntslocs
             call latlon_to_ij(LVT_domain%lvtproj, &
                  LVT_optuets(i)%tslat1, LVT_optuets(i)%tslon1, col,row)

             LVT_optuets(i)%ts_cindex1 = nint(col)
             LVT_optuets(i)%ts_rindex1 = nint(row)
             
             LVT_optuets(i)%ts_tindex = LVT_domain%gindex(nint(col),nint(row))
          enddo
       endif
    endif
     if(LVT_optuectl%optUEAlg.eq.2) then !GA 
       allocate(LVT_optuectl%fitness(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr))
       allocate(LVT_optuectl%avgfitness(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr))
       allocate(LVT_optuectl%param_val_ga(LVT_LIS_rc(1)%lnc,LVT_LIS_rc(1)%lnr,LVT_optuectl%nparam))
       allocate(LVT_optuectl%param_name(LVT_optuectl%nparam))
       if(LVT_optuectl%computeTS.eq.1) then 
          call system('mkdir -p '//(LVT_rc%statsodir))
          do i=1,LVT_optuectl%ntslocs
             LVT_optuets(i)%tslocfile(1) = trim(LVT_rc%statsodir)//'/LIS_GA_'&
                  //trim(LVT_optuets(i)%tslocname)//&
                  '.dat'
             
             LVT_optuets(i)%ftn_ts_loc(1) = LVT_getNextUnitNumber()
             open(LVT_optuets(i)%ftn_ts_loc(1),file=(LVT_optuets(i)%tslocfile(1)),&
                  form='formatted')
          enddo                    
       endif
    elseif(LVT_optuectl%optuealg.eq.4.or.&
         LVT_optuectl%optuealg.eq.5.or.LVT_optuectl%optuealg.eq.6) then 
       
       allocate(LVT_optuectl%fitness(LVT_LIS_rc(1)%ntiles,1))
       allocate(LVT_optuectl%param_val(LVT_LIS_rc(1)%ntiles,LVT_optuectl%nparam))
       allocate(LVT_optuectl%param_name(LVT_optuectl%nparam))

       if(LVT_optuectl%computeTS.eq.1) then 
          do i=1,LVT_optuectl%ntslocs
             allocate(LVT_optuets(i)%tslocfile(LVT_optuectl%nparam))
             allocate(LVT_optuets(i)%ftn_ts_loc(LVT_optuectl%nparam))
          enddo
       endif          
       
       if(LVT_optuectl%computeTS.eq.1) then 
          
          call system('mkdir -p '//(LVT_rc%statsodir))
          do i=1,LVT_optuectl%ntslocs
             do j=1,LVT_optuectl%nparam
                LVT_optuets(i)%tslocfile(j) = &
                     trim(LVT_rc%statsodir)//'/LIS_OPTUE_'&
                     //trim(LVT_optuets(i)%tslocname)//&
                     '_'//trim(LVT_optuectl%vname(j))//'.dat'
                
                
                LVT_optuets(i)%ftn_ts_loc(j) = LVT_getNextUnitNumber()
                open(LVT_optuets(i)%ftn_ts_loc(j),&
                     file=(LVT_optuets(i)%tslocfile(j)),&
                     form='formatted')
                
             enddo
             LVT_optuets(i)%tsloc_fitfile = &
                  trim(LVT_rc%statsodir)//'/LIS_OPTUE_'&
                  //trim(LVT_optuets(i)%tslocname)//&
                  '_'//'fitness.dat'
             
             LVT_optuets(i)%ftn_ts_fitloc = LVT_getNextUnitNumber()
             open(LVT_optuets(i)%ftn_ts_fitloc,&
                  file=(LVT_optuets(i)%tsloc_fitfile),&
                  form='formatted')
          enddo
       endif
    endif

  end subroutine LVT_optUEInit

!BOP
! 
! !ROUTINE: LVT_readoptUEdata
!  \label{LVT_readoptUEdata}
!
! !INTERFACE: 
  subroutine LVT_readoptUEdata(iterNo)
! 
! !USES:   
    use LVT_coreMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine reads the optUE output for the specified iteration number. 
!  For the optimization algorithms (GA), the output file is read, and for the 
!  uncertainty estimation algorithms (MCMC, DEMC), the restart file is read. 
!
!  The arguments are: 
!  \begin{description}
!   \item[iterNo] generation number 
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    integer :: iterNo
!EOP
    integer :: ftn
    integer :: t, iter_f
    character*100 :: filen
    real    :: dummy(LVT_LIS_rc(1)%lnc, LVT_LIS_rc(1)%lnr)


    ftn = LVT_getNextUnitNumber()
    
    call create_optue_filename(iterNo,LVT_optuectl%optUEAlg, filen)

    write(LVT_logunit,*) '[INFO] Reading ',filen
    open(ftn,file=(filen), status='unknown',form='unformatted')
    
    if(LVT_optuectl%optUEalg.eq.2) then        
       read(ftn) LVT_optuectl%nparam
       do t=1,LVT_optuectl%nparam+2
          if(t==LVT_optuectl%nparam+1) then 
             read(ftn) LVT_optuectl%fitness
          elseif(t==LVT_optuectl%nparam+2) then 
             read(ftn) LVT_optuectl%avgfitness
          else
             read(ftn) LVT_optuectl%param_name(t)
             read(ftn) LVT_optuectl%param_val_ga(:,:,t)
          endif
       enddo
    elseif(LVT_optuectl%optUEalg.eq.4.or.&
         LVT_optuectl%optUEalg.eq.5.or.LVT_optuectl%optUEalg.eq.6) then 

       read(ftn) iter_f
       do t=1,LVT_optuectl%nparam+1
          if(t==LVT_optuectl%nparam+1) then 
             read(ftn) LVT_optuectl%fitness(:,1)          
          else
             read(ftn) LVT_optuectl%param_name(t)
             read(ftn) LVT_optuectl%param_val(:,t)             
          endif
       enddo
    endif
    call LVT_releaseUnitNumber(ftn)

    
  end subroutine LVT_readoptUEdata

!BOP
! 
! !ROUTINE: LVT_computeoptUEstats
! \label{LVT_computeoptUEstats}
!
! !INTERFACE: 
  subroutine LVT_computeoptUEstats(iterNo)
! 
! !USES:   
    use LVT_coreMod

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This subroutine writes the optUE fitness information (at the specified 
!  locations) to an ASCII text file. For the optimization algorithm, 
!  the fitness values for the best solution and the average fitness are written.  
!  For the uncertainty estimation algorithms, the fitness and the parameter
!  values for the entire ensemble is written (into separate files). 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS:  
    integer :: iterNo
!EOP
    integer       :: i,j  
    character*10  :: fn 
    character*30  :: fmt_line
    integer       :: tid, stid, enid
    integer       :: kk

    if(LVT_optuectl%optuealg.eq.2) then 

       write(fn, '(i10)' ) LVT_optuectl%nparam + 2
       fmt_line = '(I3.3,'//trim(fn)//'E14.6)'

       do i=1,LVT_optuectl%ntslocs
          write(LVT_optuets(i)%ftn_ts_loc(1),fmt_line) &
               iterNo, &
               LVT_optuectl%fitness(LVT_optuets(i)%ts_cindex1, &
               LVT_optuets(i)%ts_rindex1), &
               LVT_optuectl%avgfitness(LVT_optuets(i)%ts_cindex1, &
               LVT_optuets(i)%ts_rindex1), &
               (LVT_optuectl%param_val_ga(LVT_optuets(i)%ts_cindex1, &
               LVT_optuets(i)%ts_rindex1,j),j=1,LVT_optuectl%nparam)
       enddo
    elseif(LVT_optuectl%optuealg.eq.4.or.&
         LVT_optuectl%optuealg.eq.5.or.LVT_optuectl%optuealg.eq.6) then 

       write(fn, '(i10)' ) LVT_LIS_rc(1)%nensem
       fmt_line = '(I3.3,'//trim(fn)//'E14.6)'
       do i=1,LVT_optuectl%ntslocs
          
          tid =  LVT_optuets(i)%ts_tindex
          stid = (tid-1)*LVT_LIS_rc(1)%nensem+1
          enid = stid+LVT_LIS_rc(1)%nensem-1
          
          do j=1,LVT_optuectl%nparam
             do kk=1,LVT_optuectl%nparam
                if(trim(LVT_optuectl%param_name(j)).eq.LVT_optuectl%vname(kk)) then 
                   exit; 
                endif
             enddo
            
             write(LVT_optuets(i)%ftn_ts_loc(kk),fmt_line) &
                  iterNo, &
                  LVT_optuectl%param_val(stid:enid,j)
          enddo
          
          write(LVT_optuets(i)%ftn_ts_fitloc,fmt_line) &
               iterNo, &
               LVT_optuectl%fitness(stid:enid,1)
       enddo
    endif

111 format(I3.3,2E14.6)

  end subroutine LVT_computeoptUEstats

!BOP
! 
! !ROUTINE: create_optue_filename
! \label{create_optue_filename}
!
! !INTERFACE: 
  subroutine create_optue_filename(no, alg, filen)
! 
! !USES:   
    use LVT_coreMod, only : LVT_rc
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine generates the name of the optUE output file (written by the 
!  LIS simulation). Note that for the uncertainty estimation algorithms 
!  (MCMC, DEMC), the restart file is read. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer          :: no

    character(len=*) :: filen
    integer          :: alg
    character(len=4) :: fiter 
    
    write(unit=fiter, fmt='(i4.4)') no
    if(alg.eq.2) then ! GA output file 
       filen = trim(LVT_rc%odir)//'/GA/GA.'&
            //trim(fiter)//'.1gd4r'
    elseif(alg.eq.4) then !MCSIM restart file
       filen = trim(LVT_rc%odir)//'/MCSIM/MCSIM.'&
            //trim(fiter)//'.MCSIMrst'
    elseif(alg.eq.5) then !MCMC restart file
       filen = trim(LVT_rc%odir)//'/MCMC/MCMC.'&
            //trim(fiter)//'.MCMCrst'
    elseif(alg.eq.6) then !DEMC restart file
       filen = trim(LVT_rc%odir)//'/DEMC/DEMC.'&
            //trim(fiter)//'.DEMCrst'
    endif
  end subroutine create_optue_filename

end module LVT_optUEMod

