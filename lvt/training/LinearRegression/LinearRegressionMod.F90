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
! !MODULE: LinearRegressionMod
! \label{LinearRegressionMod}
!
! !INTERFACE:
module LinearRegressionMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_logMod

  implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   The code in this file provides the computations for single variable
!   and two-variable linear regression. During the first pass, the 
!   the regression coefficients are computed and during the second pass, 
!   the output values are computed using the regression formulation. 
!   
!   For single variable regression, the equation is assumed to be of the   
!   form y = a + b1x 
!   
!   For the one variable regression: 
!         sum(xy)
!     b1 = --------
!         sum(x2)
!     a = mean(y) - b1 * mean (x)
! 
!   For two variable case (y = a + b1 x1 + b2 x2
!    
!         (sum(x_2^2))*(sum(x_1*y)) - (sum(x_1*x_2))*(sum(x_2*y))
!    b1 =  ------------------------------------------------------
!               (sum(x_1^2)) * (sum(x_2^2)) - sum(x_1 * x_2)^2
! 
!         (sum(x_1^2))*(sum(x_2*y)) - (sum(x_1*x_2))*(sum(x_1*y))
!    b2 =  ------------------------------------------------------
!               (sum(x_1^2)) * (sum(x_2^2)) - sum(x_1 * x_2)^2
!
!    a  = mean(y) - b1 * mean(x1) - b2 * mean(x2)
!
! !FILES USED:
!
! !NOTES: Currently do not support ensemble mode
!
! !REVISION HISTORY: 
!  14 Jul 2015    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
  PRIVATE
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: linearRegression_init
  public :: linearRegression_run
!EOP

  type, public :: linearRegDec
     character*100        :: mode
     integer              :: out_of_sample
     real, allocatable    :: sy(:,:)
     real, allocatable    :: sx1(:,:)
     real, allocatable    :: sx1x1(:,:)
     real, allocatable    :: sx1y(:,:)
     
     real, allocatable    :: sx2(:,:)
     real, allocatable    :: sx2x2(:,:)      
     real, allocatable    :: sx1x2(:,:)      
     real, allocatable    :: sx2y(:,:)      

     integer, allocatable :: nval(:,:)

     real,    allocatable :: a_coef(:,:)
     real,    allocatable :: b1_coef(:,:)
     real,    allocatable :: b2_coef(:,:)
  end type linearRegDec
  
  type (linearRegDec)       :: linearReg_struc
contains
!BOP
! 
! !ROUTINE: linearRegression_init
! \label{linearRegression_init}
!
! !INTERFACE: 
  subroutine linearRegression_init()
! 
! !USES:   
    implicit none
!
! !DESCRIPTION: 
!  This routine performs the initialization of variables
!  required for linear regression fitting. The routine
!  also reads the user specified options. 
! 
!EOP
    integer             :: nsize
    integer             :: rc

    if(LVT_rc%computeEnsMetrics.eq.1) then
       nsize = LVT_LIS_rc(1)%ntiles
    else
       nsize = LVT_rc%ngrid
    endif

    call ESMF_ConfigGetAttribute(LVT_config,linearReg_struc%mode,&
         label="Linear regression mode:",rc=rc)
    if(rc.ne.0) then 
       write(LVT_logunit,*) '[ERR] Linear regression mode: not defined'
       write(LVT_logunit,*) "[ERR] Options are: "
       write(LVT_logunit,*) "[ERR] 'single-variable', 'two-variable' "
       call LVT_endrun()
    endif

    call ESMF_ConfigGetAttribute(LVT_config,linearReg_struc%out_of_sample,&
         label="Linear regression use out of sample method:",rc=rc)
    call LVT_verify(rc,'Linear regression use out of sample method: not defined')


    if(linearReg_struc%mode.eq. "single-variable") then 
       allocate(linearReg_struc%sx1(nsize,1))
       allocate(linearReg_struc%sy(nsize,1))
       allocate(linearReg_struc%sx1x1(nsize,1))
       allocate(linearReg_struc%sx1y(nsize,1))
       allocate(linearReg_struc%nval(nsize,1))
       
       linearReg_struc%sx1   = 0 
       linearReg_struc%sy   = 0 
       linearReg_struc%sx1x1  = 0 
       linearReg_struc%sx1y  = 0 
       linearReg_struc%nval = 0 

       
       allocate(linearReg_struc%a_coef(nsize,1))
       allocate(linearReg_struc%b1_coef(nsize,1))
       linearReg_struc%a_coef  = LVT_rc%udef
       linearReg_struc%b1_coef = LVT_rc%udef

    elseif(linearReg_struc%mode.eq. "two-variable") then 

       allocate(linearReg_struc%sx1y(nsize,1))
       allocate(linearReg_struc%sx1x1(nsize,1))
       allocate(linearReg_struc%sx2x2(nsize,1))
       allocate(linearReg_struc%sx1x2(nsize,1))
       allocate(linearReg_struc%sx2y(nsize,1))
       allocate(linearReg_struc%sx1(nsize,1))
       allocate(linearReg_struc%sx2(nsize,1))
       allocate(linearReg_struc%sy(nsize,1))
       allocate(linearReg_struc%nval(nsize,1))
       
       linearReg_struc%sx1   = 0 
       linearReg_struc%sx2   = 0 
       linearReg_struc%sy    = 0        
       linearReg_struc%sx1y  = 0 
       linearReg_struc%sx1x2 = 0 
       linearReg_struc%sx2x2 = 0 
       linearReg_struc%sx1x2 = 0 
       linearREg_struc%sx2y  = 0 
       linearREg_struc%nval  = 0 

       allocate(linearReg_struc%a_coef(nsize,1))
       allocate(linearReg_struc%b1_coef(nsize,1))
       allocate(linearReg_struc%b2_coef(nsize,1))

       linearReg_struc%a_coef  = LVT_rc%udef
       linearReg_struc%b1_coef = LVT_rc%udef
       linearReg_struc%b2_coef = LVT_rc%udef
       
    endif

  end subroutine LinearRegression_init

!BOP
! 
! !ROUTINE: linearRegression_run
!  \label{linearRegression_run}
!
! !INTERFACE: 
  subroutine linearRegression_run(pass)
! 
! !USES:   
    use LVT_coreMod
    use LVT_logMod

    implicit none
!
! !INPUT PARAMETERS: 
    integer, intent(in)            :: pass
! 

! !DESCRIPTION: 
! 
!  This routine performs two tasks. During the first 
!  pass through the data, this routine computes the 
!  the linear regression coefficients. During the second
!  pass, the routine calculates the outputs using the
!  regression formulation that was developed in the 
!  first pass. 
! 
!EOP

    integer               :: i, j,m
    integer               :: nsize

    type(LVT_metadataEntry), pointer :: ds1_inp1
    type(LVT_metadataEntry), pointer :: ds1_inp2
    type(LVT_metadataEntry), pointer :: ds2

    real                             :: mean_x1, mean_x2, mean_y
    real                             :: temp_val
    logical                          :: data_flag

    if(LVT_rc%computeEnsMetrics.eq.1) then
       nsize = LVT_LIS_rc(1)%ntiles
    else
       nsize = LVT_rc%ngrid
    endif
    
    if(pass.eq.1) then 
       if(LVT_rc%computeFlag) then 
          call LVT_getDataStream1Ptr(ds1_inp1)
          call LVT_getDataStream2Ptr(ds2)

          if(linearReg_struc%mode.eq. "two-variable") then 
             ds1_inp2 => ds1_inp1%next
          endif
          
          if(linearReg_struc%mode.eq. "single-variable") then 
             if(ds1_inp1%selectNlevs.ge.1.and.ds2%selectNlevs.ge.1) then 
                
                do m=1,LVT_rc%nensem
                   do i=1,nsize
                      if(ds1_inp1%count(i,m,1).ne.0.and.ds2%count(i,m,1).ne.0) then 
                         do j=1,nsize
                            if(ds1_inp1%count(j,m,1).ne.0.and.ds2%count(j,m,1).ne.0) then 
                               data_flag = .true. 
                               if(linearReg_struc%out_of_sample.eq.1) then 
                                  if(i.eq.j) then 
                                     data_flag = .false. 
                                  endif
                               endif
                               if(data_flag) then 
                                  linearReg_struc%sx1(j,1) = linearReg_struc%sx1(j,1) + & 
                                       ds1_inp1%value(j,m,1)
                                  linearReg_struc%sx1x1(j,1) = linearReg_struc%sx1x1(j,1)+ & 
                                       ds1_inp1%value(j,m,1)*ds1_inp1%value(j,m,1)
                                  
                                  linearReg_struc%sy(j,1) = linearReg_struc%sy(j,1) + & 
                                       ds2%value(j,m,1)
                                  
                                  linearReg_struc%sx1y(j,1) = linearReg_struc%sx1y(j,1)+ & 
                                       ds1_inp1%value(j,m,1)*ds2%value(j,m,1)
                                  
                                  linearReg_struc%nval(j,1) = linearReg_struc%nval(j,1) + 1
                               endif
                            endif
                         enddo
                      endif
                   enddo
                enddo
             endif
          elseif(linearReg_struc%mode.eq. "two-variable") then 
             
             if(ds1_inp1%selectNlevs.ge.1.and.ds1_inp2%selectNlevs.ge.1&
                  .and.ds2%selectNlevs.ge.1) then 
                do m=1,LVT_rc%nensem
                   do i=1,nsize
                      if(ds1_inp1%count(i,m,1).ne.0.and.ds1_inp2%count(i,m,1).ne.0&
                           .and.ds2%count(i,m,1).ne.0) then 
                         do j=1,nsize
                            if(ds1_inp1%count(j,m,1).ne.0.and.ds1_inp2%count(j,m,1).ne.0&
                                 .and.ds2%count(j,m,1).ne.0) then 
                               data_flag = .true. 
                               if(linearReg_struc%out_of_sample.eq.1) then 
                                  if(i.eq.j) then 
                                     data_flag = .false. 
                                  endif
                               endif
                               
                               if(data_flag) then 
                                  
                                  linearReg_struc%sy(j,1) = linearReg_struc%sy(j,1)+ & 
                                       ds2%value(j,m,1)
                                  
                                  linearReg_struc%sx1(j,1) = linearReg_struc%sx1(j,1)+ & 
                                       ds1_inp1%value(j,m,1)
                                  
                                  linearReg_struc%sx2(j,1) = linearReg_struc%sx2(j,1)+ & 
                                       ds1_inp2%value(j,m,1)
                                  
                                  linearReg_struc%sx1x1(j,1) = linearReg_struc%sx1x1(j,1)+ & 
                                       ds1_inp1%value(j,m,1)*ds1_inp1%value(j,m,1)
                                  
                                  linearReg_struc%sx2x2(j,1) = linearReg_struc%sx2x2(j,1)+ & 
                                       ds1_inp2%value(j,m,1)*ds1_inp2%value(j,m,1)
                                  
                                  linearReg_struc%sx1y(j,1) = linearReg_struc%sx1y(j,1)+ & 
                                       ds1_inp1%value(j,m,1)*ds2%value(j,m,1)
                                  
                                  linearReg_struc%sx2y(j,1) = linearReg_struc%sx2y(j,1)+ & 
                                       ds1_inp2%value(j,m,1)*ds2%value(j,m,1)
                                  
                                  
                                  linearReg_struc%nval(j,1) = linearReg_struc%nval(j,1) + 1
                                  
                               endif
                            endif
                         enddo
                      endif
                   enddo
                enddo
             endif
          endif

       endif
       
       if(LVT_rc%endtime.eq.1) then   
          if(linearReg_struc%mode.eq. "single-variable") then 
             do i=1,nsize
                if( linearReg_struc%nval(i,1).gt.0) then 
                   !---------------------------------------------------
                   !         sum(xy)
                   !     b = --------
                   !         sum(x2)
                   !---------------------------------------------------

                   linearReg_struc%b1_coef(i,1) = linearReg_struc%sx1y(i,1) / & 
                        linearReg_struc%sx1x1(i,1)
                   !---------------------------------------------------
                   !          a = mean(y) - b * mean (x)
                   !---------------------------------------------------
                   mean_y = linearReg_struc%sy(i,1)/linearReg_struc%nval(i,1)
                   mean_x1 = linearReg_struc%sx1(i,1)/linearReg_struc%nval(i,1)
                   
                   linearReg_struc%a_coef(i,1) = mean_y - &
                        linearReg_struc%b1_coef(i,1)*mean_x1
                else                                
                   linearReg_struc%a_coef(i,1) = LVT_rc%udef
                   linearReg_struc%b1_coef(i,1) = LVT_rc%udef
                endif
             enddo
          elseif(linearReg_struc%mode.eq. "two-variable") then 
             do i=1,nsize
                if( linearReg_struc%nval(i,1).gt.0) then 
                   !------------------------------------------------------------------
                   !         (sum(x_2^2))*(sum(x_1*y)) - (sum(x_1*x_2))*(sum(x_2*y))
                   !    b1 =  ------------------------------------------------------
                   !               (sum(x_1^2)) * (sum(x_2^2)) - sum(x_1 * x_2)^2
                   ! 
                   !         (sum(x_1^2))*(sum(x_2*y)) - (sum(x_1*x_2))*(sum(x_1*y))
                   !    b2 =  ------------------------------------------------------
                   !               (sum(x_1^2)) * (sum(x_2^2)) - sum(x_1 * x_2)^2
                   !-----------------------------------------------------------------

                   linearReg_struc%b1_coef(i,1) = (linearReg_struc%sx2x2(i,1)* & 
                        linearReg_struc%sx1y(i,1) - linearReg_struc%sx1x2(i,1)*& 
                        linearReg_struc%sx2y(i,1)) / &
                        (linearReg_struc%sx1x1(i,1) * linearReg_struc%sx2x2(i,1) - & 
                        linearReg_struc%sx1x2(i,1) )

                   linearReg_struc%b2_coef(i,1) = (linearReg_struc%sx1x1(i,1)* & 
                        linearReg_struc%sx2y(i,1) - linearReg_struc%sx1x2(i,1)*& 
                        linearReg_struc%sx1y(i,1)) / &
                        (linearReg_struc%sx1x1(i,1) * linearReg_struc%sx2x2(i,1) - & 
                        linearReg_struc%sx1x2(i,1) )

                   mean_y = linearReg_struc%sy(i,1)/linearReg_struc%nval(i,1)
                   mean_x1 = linearReg_struc%sx1(i,1)/linearReg_struc%nval(i,1)
                   mean_x2 = linearReg_struc%sx2(i,1)/linearReg_struc%nval(i,1)
                   !-----------------------------------------------------------------
                   !    a  = mean(y) - b1 * mean(x1) - b2 * mean(x2)                 
                   !-----------------------------------------------------------------  
                   linearReg_struc%a_coef(i,1) = mean_y - &
                        linearReg_struc%b1_coef(i,1)*mean_x1 - & 
                        linearReg_struc%b2_coef(i,1)*mean_x2 
                else
                   linearReg_struc%a_coef(i,1) = LVT_rc%udef
                   linearReg_struc%b1_coef(i,1) = LVT_rc%udef
                   linearReg_struc%b2_coef(i,1) = LVT_rc%udef
                endif
             enddo

          endif
       endif
    else
!----------------------------------------------------------------------
! During the second pass, apply the regression model and store them 
! in second data stream
!----------------------------------------------------------------------
       if(LVT_rc%computeEnsMetrics.eq.1) then
          nsize = LVT_LIS_rc(1)%ntiles
       else
          nsize = LVT_rc%ngrid
       endif

       if(linearReg_struc%mode.eq. "single-variable") then        
          call LVT_getDataStream1Ptr(ds1_inp1)
          call LVT_getDataStream2Ptr(ds2)
          
          if(ds1_inp1%selectNlevs.ge.1.and.ds2%selectNlevs.ge.1) then 
             do i=1,nsize
                if(ds1_inp1%count(i,m,1).ne.0.and.ds2%count(i,m,1).ne.0) then                 
                   temp_val = linearReg_struc%a_coef(i,1) + &
                        linearReg_struc%b1_coef(i,1) * ds1_inp1%value(i,m,1)
                   ds2%value(i,m,1) = temp_val
                else
                   ds2%value(i,m,1) = LVT_rc%udef
                endif
             enddo
          endif

       elseif(linearReg_struc%mode.eq. "two-variable") then        
          call LVT_getDataStream1Ptr(ds1_inp1)
          call LVT_getDataStream2Ptr(ds2)

          ds1_inp2 => ds1_inp1%next

          if(ds1_inp1%selectNlevs.ge.1.and.ds1_inp2%selectNlevs.ge.1&
               .and.ds2%selectNlevs.ge.1) then 
             do i=1,nsize
                if(ds1_inp1%count(i,m,1).ne.0.and.ds1_inp2%count(i,m,1).ne.0.and.&
                     ds2%count(i,m,1).ne.0) then                 
                   temp_val = linearReg_struc%a_coef(i,1) + &
                        linearReg_struc%b1_coef(i,1) * ds1_inp1%value(i,m,1) + & 
                        linearReg_struc%b2_coef(i,1) * ds1_inp2%value(i,m,1)
                   ds2%value(i,m,1) = temp_val
                else
                   ds2%value(i,m,1) = LVT_rc%udef
                endif
             enddo
          endif
       endif

    endif
  end subroutine linearRegression_run

end module LinearRegressionMod

