!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

! YX: This module is used to contain subroutines related to Tb prediction using trained SVMs
! The interface is used to call C function

!!**********************************************************************************************
module mwSVM_routines 

!use clsm_ensdrv_glob_param, ONLY:    &    !original
!    logit,                           &
!    logunit,                         &
!    nodata_generic

!use ldas_exceptions,        ONLY:    &    !original
!    ldas_abort,                      &
!    LDAS_GENERIC_ERROR

implicit none

private

!public :: output_svmTb, catch2mwSVM_vars    !original
public :: output_svmTb                      !kyh20170420


interface
     subroutine predict_func(input_file,output_file,model_file) bind(c, name="predict_func")  !original 
        use,intrinsic :: ISO_C_BINDING
        implicit none
        character(kind=c_char) :: input_file(*), output_file(*), model_file(*)
        !real(c_float) :: swe_lis_c      !kyh20170609
     end subroutine predict_func
end interface


contains

!**********************************************************************************************
subroutine predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                                fstr_filenamein, fstr_filenameout,fstr_modelname)   

       use LIS_logMod,   only: LIS_logunit         !kyh20170602
       use,intrinsic :: ISO_C_BINDING
       implicit none

       integer, intent(in)            :: len_filenamein, len_filenameout, len_modelname
       character(len=len_filenamein)  :: fstr_filenamein
       character(len=1,kind=c_char)   :: cstr_filenamein(len_filenamein+1)
       character(len=len_filenameout) :: fstr_filenameout
       character(len=1,kind=c_char)   :: cstr_filenameout(len_filenameout+1)
       character(len=len_modelname)   :: fstr_modelname
       character(len=1,kind=c_char)   :: cstr_modelname(len_modelname+1)
       integer::n,i
       integer::a,b
       integer::c,d

       real                           :: swe_lis               !kyh20170602
       !real(kind=c_float)             :: swe_lis_c             !kyh20170609

       !write(LIS_logunit,*) 'swe_lis_predict_via_svm=', swe_lis                    !kyh20170602

       n=len_trim(fstr_filenamein)
       do i=1,n
          cstr_filenamein(i)=fstr_filenamein(i:i)
       end do
          cstr_filenamein(n+1)=C_NULL_CHAR

       a=len_trim(fstr_filenameout)
       do b=1,a
          cstr_filenameout(b)=fstr_filenameout(b:b)
       end do
          cstr_filenameout(a+1)=C_NULL_CHAR

       c=len_trim(fstr_modelname)
       do d=1,c
          cstr_modelname(d)=fstr_modelname(d:d)
       end do
          cstr_modelname(c+1)=C_NULL_CHAR

       !swe_lis_c = swe_lis                  !kyh20170609

       !write(LIS_logunit,*) 'swe_lis_c_predict_via_svm=', swe_lis_c    !kyh20170609

!svk commenting out for now 
! need a better strategy for interfacing
!       call predict_func(cstr_filenamein, cstr_filenameout, cstr_modelname)  

end subroutine predict_via_svm 
!*********************************************************************************************
!subroutine output_svmTb(year, dateofyear, account_name, num_of_inputs,     &
!                        training_period, training_target, domain,          & 
!                        scaling_option, forest_decoup_option,              &
!                        atm_decoup_option, myid, ease_rowind, ease_colind, &
!                        predictTb_10H_36H, predictTb_10V_36V,              &
!                        predictTb_18H_36H, predictTb_18V_36V,              &
!                        predictTb_10H,                                     &
!                        predictTb_10V,                                     &
!                        predictTb_18H,                                     &
!                        predictTb_18V,                                     &
!                        predictTb_36H,                                     &
!                        predictTb_36V)                                            !original

subroutine output_svmTb(year, dateofyear, case_directory, svm_directory,   &
                        num_of_inputs,training_period, training_target,    &
                        domain, scaling_option,                            &
                        forest_decoup_option,                              &
                        atm_decoup_option, myid, ind_obs,ind_svm,      &
                        predictTb_10H_36H, predictTb_10V_36V,              &
                        predictTb_18H_36H, predictTb_18V_36V,              &
                        predictTb_10H,                                     &
                        predictTb_10V,                                     &
                        predictTb_18H,                                     &
                        predictTb_18V,                                     &
                        predictTb_36H,                                     &
                        predictTb_36V)                                          !kyh20170801  

! The definition of start_period depends on training period
! num_of_inputs == 1,2,3,4,5,6,7,8,9,10,11
! training_period == 'fortnight', 'month', 'seasonal'
! training_target == 'Tb','db'
! scaling_option == 'none','linear','standardization','unitvector'
! YX: The subroutine is used to output Tb predictions at certain frequency/frequency combinations based on training target and number of inputs
! YX: trained SVM parameters were placed in a well-structured txt file

!*********************************************************************************************!

use LIS_logMod,   only: LIS_logunit      !kyh20170423

implicit none

integer                     :: len_filenamein, len_filenameout, len_modelname
integer                     :: myid, suffix_id
logical                     :: forest_decoup_option
logical                     :: atm_decoup_option

character(:), allocatable   :: case_directory, svm_directory               !kyh201770801

character(:), allocatable   :: fstr_filenamein

character(:), allocatable   :: fstr_filenameout_10H_36H, fstr_filenameout_10V_36V, &
                               fstr_filenameout_18H_36H, fstr_filenameout_18V_36V
character(:), allocatable   :: fstr_filenameout_10H, fstr_filenameout_10V, &
                               fstr_filenameout_18H, fstr_filenameout_18V, &
                               fstr_filenameout_36H, fstr_filenameout_36V

character(:), allocatable   :: fstr_modelname_10H_36H, fstr_modelname_10V_36V, &
                               fstr_modelname_18H_36H, fstr_modelname_18V_36V
character(:), allocatable   :: fstr_modelname_10H, fstr_modelname_10V, &
                               fstr_modelname_18H, fstr_modelname_18V, &
                               fstr_modelname_36H, fstr_modelname_36V

integer                     :: year, start_period, num_of_inputs, dateofyear
integer                     :: i, ierr
!integer                     :: row_of_interest, col_of_interest              !original


real, intent(out),optional           :: predictTb_10H_36H, predictTb_18H_36H 
real, intent(out),optional           :: predictTb_10V_36V, predictTb_18V_36V
real, intent(out),optional           :: predictTb_10H, predictTb_18H, predictTb_36H
real, intent(out),optional           :: predictTb_10V, predictTb_18V, predictTb_36V


!integer                     :: ease_rowind           !original
!integer                     :: ease_colind           !original

logical                        :: OK
integer                        :: isize
!character(len=10)              :: account_name          !original
character(len=9)               :: training_period
character(len=2)               :: FOY
character(len=4)               :: YYYY
character(len=3)               :: DOY
!character(len=3)               :: str_row_of_interest    !original
!character(len=4)               :: str_col_of_interest    !original
character(len=6)               :: str_ind_svm             !kyh20170425
character(len=6)               :: str_ind_obs             !kyh20171016
character(len=2)               :: str_num_of_inputs, str_myid 
character(len=2)               :: training_target
character(len=4)               :: scaling_option
!character(len=3)               :: domain                 !original
character(len=9)               :: domain                 !kyh20170421
character(:), allocatable      :: filepath
character(len=30),parameter    :: Iam = 'output_svmTb'
character(len=100)             :: err_msg

integer                        :: ind_svm                !kyh20170421
integer                        :: ind_obs                !kyh20171016
real                           :: swe_lis                !kyh20170602

integer                        :: IOstatus               !kyh20170912

!write(LIS_logunit, *) 'case_directory_output_svm=', case_directory    !kyh20170801
!write(LIS_logunit, *) 'svm_directory_output_svm=', svm_directory      !kyh20170801

write (str_num_of_inputs, '(I2.2)') num_of_inputs
write (DOY, '(I3.3)') dateofyear

! Define training period
select case (training_period)
  case ('fortnight')
      start_period = ceiling(dateofyear / 14.0)
      if (start_period .GE. 27) start_period = 26
      write (FOY, '(I2.2)') start_period
      write (YYYY, '(I4.4)') year 
end select

suffix_id = myid + 10
write(str_myid, '(I2.2)') suffix_id

! Find the right folder
if ( (forest_decoup_option) .and. (atm_decoup_option) ) then

  select case (training_target)

    case('db')
       filepath = svm_directory// 'Decouple_DTB_SVM_Param_'//trim(domain)//'_'//str_num_of_inputs

    case('Tb')
       filepath = svm_directory// 'Decouple_TB_SVM_Param_'//trim(domain)//'_'//str_num_of_inputs

  end select


else

  select case (training_target)

    case('db')
       filepath = svm_directory//'ALLYEARS_DTB_SVM_Param_'//trim(domain)//'_'//str_num_of_inputs//'/SVM_txt'
 
    case('Tb')
       filepath = svm_directory//'ALLYEARS_TB_SVM_Param_'//trim(domain)//'_'//str_num_of_inputs//'/SVM_txt'

  end select

end if
 
! reference table for specific number of inputs
!select case (num_of_inputs)
!    case(1)
!      compare_str = 'SWE_m'
!    case(2)
!      compare_str = 'SWE_m' // 'SLWCT'
!    case(3)
!      compare_str = 'Tp1_K' // 'SLWCT' // 'SWE_m'
!    case(4)
!      compare_str = 'Tp1_K' // 'SLWCT' // 'SWE_m' // 'Tsurf'
!    case(5)
!      compare_str = 'Tp1_K' // 'SLWCT' // 'SWE_m' // 'Tsurf' // 'rhsn3'
!    case(6)
!      compare_str = 'Tp1_K' // 'SLWCT' // 'SWE_m' // 'Tsurf' // 'rhsn2' // 'rhsn3'
!    case(7)
!      compare_str = 'Tp1_K' // 'TpsnN' // 'SLWCT' // 'SWE_m' // 'Tsurf' // 'rhsn2' // 'rhsn3'
!    case(8)
!      compare_str = 'Tp1_K' // 'TpsnN' // 'SLWCT' // 'SWE_m' // 'Tsurf' // 'rhsn1' // 'rhsn2' // 'rhsn3'
!    case(9)
!      compare_str = 'Tp1_K' // 'TpsnN' // 'SLWCT' // 'SWE_m' // 'Tsurf' // 'TairK' // 'rhsn1' // 'rhsn2' // 'rhsn3'
!    case(10)
!      compare_str = 'Tp1_K' // 'Tpsn1' // 'TpsnN' // 'SLWCT' // 'SWE_m' // 'Tsurf' // 'TairK' // 'rhsn1' // 'rhsn2' // 'rhsn3'
!end select


!row_of_interest = ease_rowind
!col_of_interest = ease_colind

!write(str_row_of_interest,'(I3.3)') row_of_interest
!write(str_col_of_interest,'(I4.4)') col_of_interest 

write(str_ind_svm,'(I6.6)') ind_svm             !kyh20170425
write(str_ind_obs,'(I6.6)') ind_obs           !kyh20171016

!write(LIS_logunit,*) 'ind_obs=', ind_obs         !kyh20171017

! Concatenate right svm parameter name
select case (training_target)

  case('db')
     !-------------------------------------original
     !fstr_modelname_10H_36H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10H_36H'
     !fstr_modelname_10V_36V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10V_36V'
     !fstr_modelname_18H_36H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18H_36H'
     !fstr_modelname_18V_36V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18V_36V'

     !-------------------------------------------kyh20170425
     fstr_modelname_10H_36H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10H_36H'
     fstr_modelname_10V_36V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10V_36V'
     fstr_modelname_18H_36H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18H_36H'
     fstr_modelname_18V_36V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18V_36V'
     !------------------------------------------

   case('Tb')
     !----------------------------------------------original
     !fstr_modelname_10H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10H'
     !fstr_modelname_10V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10V'
     !fstr_modelname_18H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18H'
     !fstr_modelname_18V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18V' 
     !fstr_modelname_36H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_36H'
     !fstr_modelname_36V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
     !       '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
     !       '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_36V' 

     !----------------------------------------------kyh20170425
     fstr_modelname_10H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10H'
     fstr_modelname_10V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_10V'
     fstr_modelname_18H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18H'
     fstr_modelname_18V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_18V'
     fstr_modelname_36H = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_36H'
     fstr_modelname_36V = filepath // '/svm_fn' // FOY // 'ALLYEARS' // &
            '_svmgid'// str_ind_svm // '_' // trim(training_target) // '_global_' // str_num_of_inputs // '_36V'
     !----------------------------------------------

end select

!Inquire the existence of parameter txt file
!only need to check one of the frequency
!E.G. if the parameter file for 10H exists, then parameter file for 10V, 18V, 18H, 36H, 36V also exist

select case (training_target)

   case('db') 
      inquire(file=fstr_modelname_10H_36H, exist=OK)

      !kyh_debug
      !write(LIS_logunit,*) 'fstr_modelname_10H_36H=', fstr_modelname_10H_36H
      !write(LIS_logunit,*) 'OK=', OK

   case('Tb')
      inquire(file=fstr_modelname_10H, exist=OK)

      !kyh_debug
      !write(LIS_logunit,*) 'fstr_modelname_10H=', fstr_modelname_10H
      !write(LIS_logunit,*) 'OK=', OK

end select


if(OK) then 

  !--------------------------------------------kyh20170912
  suffix_id = myid + 10
  write(str_myid, '(I2.2)') suffix_id
  !--------------------------------------------kyh20170912

  ! Give proper names to output files

    select case (training_target)

      case('db')
        !-------------------------------------------------original
     	!fstr_filenameout_10H_36H = '/lustre/' // trim(account_name) // '/db_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_db_' // domain // '_' // str_num_of_inputs // '_10H_36H'
     	!fstr_filenameout_18H_36H = '/lustre/' // trim(account_name) // '/db_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_db_' // domain // '_' // str_num_of_inputs // '_18H_36H'
     	!fstr_filenameout_10V_36V = '/lustre/' // trim(account_name) // '/db_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_db_' // domain // '_' // str_num_of_inputs // '_10V_36V'
    	!fstr_filenameout_18V_36V = '/lustre/' // trim(account_name) // '/db_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_db_' // domain // '_' // str_num_of_inputs // '_18V_36V'
      
        !-------------------------------------------------kyh20170425
        fstr_filenameout_10H_36H = case_directory// &
                                   'OUTPUT/db_output/PS_SVM_fn' &
                                   // FOY // 'vy' // YYYY // &
                                   '_svmgid'// str_ind_svm // &
                                   '_db_' // trim(domain) // '_' // str_num_of_inputs // '_10H_36H_'// str_myid
        fstr_filenameout_18H_36H = case_directory// &
                                   'OUTPUT/db_output/PS_SVM_fn' &
                                   // FOY // 'vy' // YYYY // &
                                   '_svmgid'// str_ind_svm // &
                                   '_db_' // trim(domain) // '_' // str_num_of_inputs // '_18H_36H_'// str_myid   
        fstr_filenameout_10V_36V = case_directory// &
                                   'OUTPUT/db_output/PS_SVM_fn' &
                                   // FOY // 'vy' // YYYY // &
                                   '_svmgid'// str_ind_svm // &
                                   '_db_' // trim(domain) // '_' // str_num_of_inputs // '_10V_36V_'// str_myid
        fstr_filenameout_18V_36V = case_directory// &
                                   'OUTPUT/db_output/PS_SVM_fn' &
                                   // FOY // 'vy' // YYYY // &
                                   '_svmgid'// str_ind_svm // &
                                   '_db_' // trim(domain) // '_' // str_num_of_inputs // '_18V_36V_'// str_myid


        !write(LIS_logunit,*) 'fstr_filenameout_10H_36H =', fstr_filenameout_10H_36H
        !write(LIS_logunit,*) 'fstr_filenameout_18H_36H =', fstr_filenameout_18H_36H
        !write(LIS_logunit,*) 'fstr_filenameout_10V_36V =', fstr_filenameout_10V_36V
        !write(LIS_logunit,*) 'fstr_filenameout_18V_36V =', fstr_filenameout_18V_36V
        !------------------------------------------------------------

      case('Tb')
        !--------------------------------------------------original
     	!fstr_filenameout_10H = '/lustre/' // trim(account_name) // '/Tb_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_Tb_' // domain // '_' // str_num_of_inputs // '_10H'
     	!fstr_filenameout_18H = '/lustre/' // trim(account_name) // '/Tb_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_Tb_' // domain // '_' // str_num_of_inputs // '_18H'
     	!fstr_filenameout_10V = '/lustre/' // trim(account_name) // '/Tb_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_Tb_' // domain // '_' // str_num_of_inputs // '_10V'
     	!fstr_filenameout_18V = '/lustre/' // trim(account_name) // '/Tb_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_Tb_' // domain // '_' // str_num_of_inputs // '_18V'
     	!fstr_filenameout_36H = '/lustre/' // trim(account_name) // '/Tb_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_Tb_' // domain // '_' // str_num_of_inputs // '_36H'
     	!fstr_filenameout_36V = '/lustre/' // trim(account_name) // '/Tb_output/PS_SVM_fn' // FOY // 'vy' // YYYY // &
        !    '_r'// str_row_of_interest // '_c' // str_col_of_interest // &
        !    '_Tb_' // domain // '_' // str_num_of_inputs // '_36V'

        !----------------------------------------------------kyh20170425
        fstr_filenameout_10H = case_directory// &
                               'OUTPUT/Tb_output/PS_SVM_fn' &
                               // FOY // 'vy' // YYYY // &
                               '_svmgid'// str_ind_svm // &
                               '_Tb_' // trim(domain) // '_' // str_num_of_inputs // '_10H_'// str_myid
        fstr_filenameout_18H = case_directory// &
                               'OUTPUT/Tb_output/PS_SVM_fn' &
                               // FOY // 'vy' // YYYY // &
                               '_svmgid'// str_ind_svm // &
                               '_Tb_' // trim(domain) // '_' // str_num_of_inputs // '_18H_'// str_myid
        fstr_filenameout_10V = case_directory// &
                               'OUTPUT/Tb_output/PS_SVM_fn' &
                               // FOY // 'vy' // YYYY // &
                               '_svmgid'// str_ind_svm // &
                               '_Tb_' // trim(domain) // '_' // str_num_of_inputs // '_10V_'// str_myid
        fstr_filenameout_18V = case_directory// &
                               'OUTPUT/Tb_output/PS_SVM_fn' &
                               // FOY // 'vy' // YYYY // &
                               '_svmgid'// str_ind_svm // &
                               '_Tb_' // trim(domain) // '_' // str_num_of_inputs // '_18V_'// str_myid
        fstr_filenameout_36H = case_directory// &
                               'OUTPUT/Tb_output/PS_SVM_fn' &
                               // FOY // 'vy' // YYYY // &
                               '_svmgid'// str_ind_svm // &
                               '_Tb_' // trim(domain) // '_' // str_num_of_inputs // '_36H_'// str_myid
        fstr_filenameout_36V = case_directory// &
                               'OUTPUT/Tb_output/PS_SVM_fn' &
                               // FOY // 'vy' // YYYY // &
                               '_svmgid'// str_ind_svm // &
                               '_Tb_' // trim(domain) // '_' // str_num_of_inputs // '_36V_'// str_myid
 
        !write(LIS_logunit,*) 'fstr_filenameout_10H =', fstr_filenameout_10H
        !write(LIS_logunit,*) 'fstr_filenameout_18H =', fstr_filenameout_18H
        !write(LIS_logunit,*) 'fstr_filenameout_10V =', fstr_filenameout_10V
        !write(LIS_logunit,*) 'fstr_filenameout_18V =', fstr_filenameout_18V  
        !write(LIS_logunit,*) 'fstr_filenameout_36H =', fstr_filenameout_36H
        !write(LIS_logunit,*) 'fstr_filenameout_36V =', fstr_filenameout_36V
        !---------------------------------------------------------------

   end select

     !-------------------------------------------original
     ! Concatenate right LDASsa input file names
     !fstr_filenamein = '/lustre/' // trim(account_name) // '/LDASsa_input/PS_LDASsa_'//DOY//'vy'//YYYY//'_'//str_num_of_inputs//'_r'//str_row_of_interest//'_c'//str_col_of_interest//'_'//str_myid
     
     !-------------------------------------------kyh20170425
     ! Concatenate right LIS input file names
     fstr_filenamein = case_directory// &
                       'OUTPUT/SVM_inputs/LIS_svm_' &
                       //DOY//'vy'//YYYY//'_'//str_num_of_inputs//'_obsgid'//str_ind_obs//'_'//str_myid

     !write(LIS_logunit,*) 'fstr_filenamein=', fstr_filenamein
     !------------------------------------------------------     
                       

    select case (training_target)

       case('db')
          len_filenamein = len(trim(fstr_filenamein))
          len_filenameout = len(trim(fstr_filenameout_10H_36H))
          len_modelname = len(trim(fstr_modelname_10H_36H))
       case('Tb')
          len_filenamein = len(trim(fstr_filenamein))
          len_filenameout = len(trim(fstr_filenameout_10H))
          len_modelname = len(trim(fstr_modelname_10H))

    end select
        
    
    select case (training_target)

      case('db')

     	! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_10H_36H, fstr_modelname_10H_36H)   
   
     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_10H_36H, size=isize)
        !write(LIS_logunit,*) 'isize=', isize                   !kyh_debug
     	if (isize>0) then
        	open(unit=70+myid, file=fstr_filenameout_10H_36H,status='old',action='read')
        	!read(70+myid,*) predictTb_10H_36H        !original
                read(70+myid,*,IOSTAT=IOstatus) predictTb_10H_36H     !kyh20170912
        	if ((predictTb_10H_36H  .LT. -100.) .or. (predictTb_10H_36H  .GT. 100.)) predictTb_10H_36H  = -9999.              
        	close(70+myid,status='delete')    !original
                !close(70+myid)                     !kyh20170427
        else
        	predictTb_10H_36H = -9999.
        end if

        !kyh check state variables
        !write(LIS_logunit,*) 'swe_lis_output_svm_Tb=', swe_lis

        ! Call C-function based svm predictors
        call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_10V_36V, fstr_modelname_10V_36V)   
 
        ! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_10V_36V, size=isize)
        !write(LIS_logunit,*) 'isize=', isize               !kyh_debug
     	if (isize>0) then
        	open(unit=71+myid, file=fstr_filenameout_10V_36V,status='old',action='read')
        	!read(71+myid,*) predictTb_10V_36V     !original
                read(71+myid,*,IOSTAT=IOstatus) predictTb_10V_36V     !kyh20170912
        	if ((predictTb_10V_36V  .LT. -100.) .or. (predictTb_10V_36V  .GT. 100.)) predictTb_10V_36V  = -9999.
        	close(71+myid,status='delete')      !original
                !close(71+myid)                       !kyh20170427
     	else
       		predictTb_10V_36V = -9999.
     	end if


        ! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_18H_36H, fstr_modelname_18H_36H)  

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_18H_36H, size=isize)
     	if (isize>0) then
        	open(unit=72+myid, file=fstr_filenameout_18H_36H,status='old',action='read')
        	!read(72+myid,*) predictTb_18H_36H            !original
                read(72+myid,*,IOSTAT=IOstatus) predictTb_18H_36H     !kyh20170912
        	if ((predictTb_18H_36H  .LT. -100.) .or. (predictTb_18H_36H  .GT. 100.)) predictTb_18H_36H  = -9999.
        	close(72+myid,status='delete')       !original
                !close(72+myid)                        !kyh20170427 
     	else
        	predictTb_18H_36H = -9999.
     	end if

    
     	! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_18V_36V, fstr_modelname_18V_36V)   

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_18V_36V, size=isize)
     	if (isize>0) then
        	open(unit=73+myid, file=fstr_filenameout_18V_36V,status='old',action='read')
        	!read(73+myid,*) predictTb_18V_36V             !original
                read(73+myid,*,IOSTAT=IOstatus) predictTb_18V_36V     !kyh20170912
        	if ((predictTb_18V_36V  .LT. -100.) .or. (predictTb_18V_36V  .GT. 100.)) predictTb_18V_36V  = -9999.
        	close(73+myid,status='delete')         !original
                !close(73+myid)                          !kyh20170427
     	else
       		predictTb_18V_36V = -9999.
     	end if

     !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     case('Tb')

        ! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_10H, fstr_modelname_10H)       

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_10H, size=isize)

        !kyh_debug
        !write(LIS_logunit,*) 'filename=', fstr_filenameout_10H
        !write(LIS_logunit,*) 'isize=', isize

     	if (isize>0) then
        	open(unit=70+myid, file=fstr_filenameout_10H,status='old',action='read')
        	!read(70+myid,*) predictTb_10H    !original
                read(70+myid,*,IOSTAT=IOstatus) predictTb_10H     !kyh20170912 
                !write(*,*) 'kyh_IOSTAT_10H=', IOstatus                
        	if ((predictTb_10H  .LT. 100.) .or. (predictTb_10H  .GT. 400.)) predictTb_10H  = -9999.
        	close(70+myid,status='delete')      !original
                !close(70+myid)                       !kyh20170427
     	else
        	predictTb_10H = -9999.
     	end if


     	! Call C-function based svm predictors
    	 call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_18H, fstr_modelname_18H)      

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_18H, size=isize)
     	if (isize>0) then
        	open(unit=71+myid, file=fstr_filenameout_18H,status='old',action='read')
        	!read(71+myid,*) predictTb_18H      !original
                read(71+myid,*,IOSTAT=IOstatus) predictTb_18H     !kyh20170912
                !write(*,*) 'kyh_IOSTAT_18H=', IOstatus
        	if ((predictTb_18H  .LT. 100.) .or. (predictTb_18H  .GT. 400.)) predictTb_18H  = -9999.
        	close(71+myid,status='delete')        !original
                !close(71+myid)                         !kyh20170427
     	else
        	predictTb_18H = -9999.
     	end if

     	! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_36H, fstr_modelname_36H)     

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_36H, size=isize)
     	if (isize>0) then
        	open(unit=72+myid, file=fstr_filenameout_36H,status='old',action='read')
        	!read(72+myid,*) predictTb_36H           !original
                read(72+myid,*,IOSTAT=IOstatus) predictTb_36H     !kyh20170912
                !write(*,*) 'kyh_IOSTAT_36H=', IOstatus
        	if ((predictTb_36H  .LT. 100.) .or. (predictTb_36H  .GT. 400.)) predictTb_36H  = -9999.
        	close(72+myid,status='delete')           !original
                !close(72+myid)                            !kyh20170427
     	else
        	predictTb_36H = -9999.
     	end if

     	! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_10V, fstr_modelname_10V)        
       
     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_10V, size=isize)
     	if (isize>0) then
        	open(unit=73+myid, file=fstr_filenameout_10V,status='old',action='read')
        	!read(73+myid,*) predictTb_10V              !original
                read(73+myid,*,IOSTAT=IOstatus) predictTb_10V     !kyh20170912
                !write(*,*) 'kyh_IOSTAT_10V=', IOstatus
        	if ((predictTb_10V  .LT. 100.) .or. (predictTb_10V  .GT. 400.)) predictTb_10V  = -9999.
        	close(73+myid,status='delete')         !original
                !close(73+myid)                          !kyh20170427
     	else
        	predictTb_10V = -9999.
     	end if

     	! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_18V, fstr_modelname_18V)  

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_18V, size=isize)
     	if (isize>0) then
        	open(unit=74+myid, file=fstr_filenameout_18V,status='old',action='read')
        	!read(74+myid,*) predictTb_18V                          !original
                read(74+myid,*,IOSTAT=IOstatus) predictTb_18V     !kyh20170912
                !write(*,*) 'kyh_IOSTAT_18V=', IOstatus
        	if ((predictTb_18V  .LT. 100.) .or. (predictTb_18V  .GT. 400.)) predictTb_18V  = -9999.
        	close(74+myid,status='delete')          !original
                !close(74+myid)                           !kyh20170427
     	else
        	predictTb_18V = -9999.
     	end if

     	! Call C-function based svm predictors
     	call predict_via_svm(len_filenamein, len_filenameout, len_modelname, &
                     fstr_filenamein, fstr_filenameout_36V, fstr_modelname_36V)         

     	! Organize prediction into one dimension, following the i/j of tile_data
     	inquire(file=fstr_filenameout_36V, size=isize)
     	if (isize>0) then
        	open(unit=75+myid, file=fstr_filenameout_36V,status='old',action='read')
        	!read(75+myid,*) predictTb_36V                  !original
                read(75+myid,*,IOSTAT=IOstatus) predictTb_36V     !kyh20170912
                !write(*,*) 'kyh_IOSTAT_36V=', IOstatus
        	if ((predictTb_36V  .LT. 100.) .or. (predictTb_36V  .GT. 400.)) predictTb_36V  = -9999.
        	close(75+myid,status='delete')            !origianl
                !close(75+myid)                             !kyh20170427
     	else
        	predictTb_36V = -9999.
     	end if

    end select

   !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     else
       
!          if(logit) write (logunit,*) 'No file exist for point #' // ii_str // ' move on to the next'

     select case (training_target)

        case('db')
           predictTb_10H_36H = -9999.
           predictTb_10V_36V = -9999.
           predictTb_18H_36H = -9999.
           predictTb_18V_36V = -9999.
        case('Tb')
           predictTb_10H = -9999.
           predictTb_10V = -9999.
           predictTb_18H = -9999.
           predictTb_18V = -9999.
           predictTb_36H = -9999.
           predictTb_36V = -9999.

      end select


end if


end subroutine output_svmTb 

end module mwSVM_routines 
