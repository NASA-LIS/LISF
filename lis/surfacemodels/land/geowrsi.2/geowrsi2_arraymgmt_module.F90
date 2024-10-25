!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module geowrsi2_arraymgmt_module

! Basic, generic fortran 90 memory management library for nullifying,
! allocating, and deleting/accessing in a controlled fashion
! for WRSI pointers, arrays, etc.
!
! Translated by Brad Wind, ESSIC, 2011
! 25 Oct 2013: KR Arsenault;  Added GeoWRSI2.0 model to LIS-7
!
   use fbil_module

  implicit none

Contains

!
! The following set of routines utilize optional arguments but are intended to 
! be called with certain subsets of those potential argument lists.
! The dtypePtr#d variables are ideally mutually-exclusive within one call 
! and must be accompanied by the correct number of dimension-size specifiers.
!

 subroutine nullify_ptr(ptrIsNull, lgclptr1d, lgclptr2d,   &
                       int2ptr2d, int2ptr3d,   &
                       int4ptr1d, int4ptr2d,   &
                       real4ptr2d, real8ptr2d, &
                       charNptr1d  )

! This subroutine nullfies arrays.  It does not deallocate to free up memory. 
   logical, intent(inout)                        :: ptrIsNull
   logical, pointer, dimension(:),      optional :: lgclptr1d
   logical, pointer, dimension(:,:),    optional :: lgclptr2d
   integer*2, pointer, dimension(:,:),  optional :: int2ptr2d
   integer*2, pointer, dimension(:,:,:),optional :: int2ptr3d
   integer*4, pointer, dimension(:),    optional :: int4ptr1d
   integer*4, pointer, dimension(:,:),  optional :: int4ptr2d
   real*4, pointer, dimension(:,:),     optional :: real4ptr2d
   real*8, pointer, dimension(:,:),     optional :: real8ptr2d
   type(charN), pointer, dimension(:),  optional :: charNptr1d

   if (present(lgclptr1d)) nullify(lgclptr1d)
   if (present(lgclptr2d)) nullify(lgclptr2d)
   if (present(int2ptr2d)) nullify(int2ptr2d)
   if (present(int2ptr3d)) nullify(int2ptr3d)
   if (present(int4ptr1d)) nullify(int4ptr1d)
   if (present(int4ptr2d)) nullify(int4ptr2d)
   if (present(real4ptr2d)) nullify(real4ptr2d)
   if (present(real8ptr2d)) nullify(real8ptr2d)
   if (present(charNptr1d)) nullify(charNptr1d)

   ptrIsNull = .true.

 end subroutine nullify_ptr


 subroutine dealloc_arr(ptrIsNull, dim1PrevSz,  &
                       lgclptr1d, lgclptr2d,   &
                       int2ptr2d, int2ptr3d,   &
                       int4ptr1d, int4ptr2d,   &
                       real4ptr2d, real8ptr2d, &
                       charNptr1d  )

! This subroutine deallocates to free up memory and nullifies the arrays.
   logical, intent(inout)                         :: ptrIsNull
   integer*4,                            optional :: dim1PrevSz
   logical, pointer, dimension(:),       optional :: lgclptr1d
   logical, pointer, dimension(:,:),     optional :: lgclptr2d
   integer*2, pointer, dimension(:,:),   optional :: int2ptr2d
   integer*2, pointer, dimension(:,:,:), optional :: int2ptr3d
   integer*4, pointer, dimension(:),     optional :: int4ptr1d
   integer*4, pointer, dimension(:,:),   optional :: int4ptr2d
   real*4, pointer, dimension(:,:),      optional :: real4ptr2d
   real*8, pointer, dimension(:,:),      optional :: real8ptr2d
   type(charN), pointer, dimension(:),   optional :: charNptr1d

   integer*4 :: i
!_______________________________________________________________

   if (.not. ptrIsNull) then

      if (present(lgclptr1d)) then
         deallocate(lgclptr1d)
         call nullify_ptr(ptrIsNull, lgclptr1d=lgclptr1d)
      endif
      if (present(lgclptr2d)) then
         deallocate(lgclptr2d)
         call nullify_ptr(ptrIsNull, lgclptr2d=lgclptr2d)
      endif
      if (present(int2ptr2d)) then
         deallocate(int2ptr2d)
         call nullify_ptr(ptrIsNull, int2ptr2d=int2ptr2d)
      endif
      if (present(int2ptr3d)) then
         deallocate(int2ptr3d)
         call nullify_ptr(ptrIsNull, int2ptr3d=int2ptr3d)
      endif
      if (present(int4ptr1d)) then
         deallocate(int4ptr1d)
         call nullify_ptr(ptrIsNull, int4ptr1d=int4ptr1d)
      endif
      if (present(int4ptr2d)) then
         deallocate(int4ptr2d)
         call nullify_ptr(ptrIsNull, int4ptr2d=int4ptr2d)
      endif
      if (present(real4ptr2d)) then
         deallocate(real4ptr2d)
         call nullify_ptr(ptrIsNull, real4ptr2d=real4ptr2d)
      endif
      if (present(real8ptr2d)) then
         deallocate(real8ptr2d)
         call nullify_ptr(ptrIsNull, real8ptr2d=real8ptr2d)
      endif
      if (present(charNptr1d)) then

         if(present(dim1PrevSz)) then
            do i=1, dim1PrevSz, 1
               write(charNptr1d(i)%str, '(a)') ''
               deallocate(charNptr1d(i)%str)
               nullify(charNptr1d(i)%str)
            end do
         endif

         deallocate(charNptr1d)
         call nullify_ptr(ptrIsNull, charNptr1d=charNptr1d)
      endif

   endif

 end subroutine dealloc_arr

 function alloc_arr( ptrIsNull, dim1PrevSz,  &
                    dim1Sz, dim2Sz, dim3Sz, &
                    lgclptr1d, lgclptr2d,   &
                    int2ptr2d, int2ptr3d,   &
                    int4ptr1d, int4ptr2d,   &
                    real4ptr2d, real8ptr2d, charNptr1d )

! This function deallocates to free up memory as necessary and allocates memory.
! Upon successful allocation returns .false. (i.e. to assign 'not null' to ptrIsNull)
! Otherwise, returns .true.

   logical                                        :: alloc_arr  ! Function result

   logical, intent(inout)                         :: ptrIsNull
   integer*4, intent(in),                optional :: dim1PrevSz ! This is for dealloc purposes
   integer*4, intent(in)                          :: dim1Sz     ! The new size to allocate
   integer*4, intent(in),                optional :: dim2Sz
   integer*4, intent(in),                optional :: dim3Sz
   logical, pointer, dimension(:),       optional :: lgclptr1d
   logical, pointer, dimension(:,:),     optional :: lgclptr2d
   integer*2, pointer, dimension(:,:),   optional :: int2ptr2d
   integer*2, pointer, dimension(:,:,:), optional :: int2ptr3d
   integer*4, pointer, dimension(:),     optional :: int4ptr1d
   integer*4, pointer, dimension(:,:),   optional :: int4ptr2d
   real*4, pointer, dimension(:,:),      optional :: real4ptr2d
   real*8, pointer, dimension(:,:),      optional :: real8ptr2d
   type(charN), pointer, dimension(:),   optional :: charNptr1d

   integer    :: allocSuccess0else1
   integer*4  :: i, j
! ________________________________________________________________

   alloc_arr = .true.
   allocSuccess0else1 = 1

   if (.not. ptrIsNull) then
      if (present(lgclptr1d))  call dealloc_arr(ptrIsNull, lgclptr1d=lgclptr1d) 
      if (present(lgclptr2d))  call dealloc_arr(ptrIsNull, lgclptr2d=lgclptr2d) 
      if (present(int2ptr2d))  call dealloc_arr(ptrIsNull, int2ptr2d=int2ptr2d)
      if (present(int2ptr3d))  call dealloc_arr(ptrIsNull, int2ptr3d=int2ptr3d)
      if (present(int4ptr1d))  call dealloc_arr(ptrIsNull, int4ptr1d=int4ptr1d)
      if (present(int4ptr2d))  call dealloc_arr(ptrIsNull, int4ptr2d=int4ptr2d)
      if (present(real4ptr2d)) call dealloc_arr(ptrIsNull, real4ptr2d=real4ptr2d)
      if (present(real8ptr2d)) call dealloc_arr(ptrIsNull, real8ptr2d=real8ptr2d)
      if (present(charNptr1d)) then
         if(.not.(present(dim1PrevSz))) then
            call dealloc_arr(ptrIsNull, charNptr1d=charNptr1d)
         else
            call dealloc_arr(ptrIsNull, dim1PrevSz=dim1PrevSz, charNptr1d=charNptr1d) 
         endif
      endif
   endif

   if (present(lgclptr1d)) then
      allocate(lgclptr1d(dim1Sz), STAT=allocSuccess0else1) 
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, lgclptr1d=lgclptr1d)
   endif 
   if (present(lgclptr2d)) then
      allocate(lgclptr2d(dim1Sz, dim2Sz), STAT=allocSuccess0else1) 
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, lgclptr2d=lgclptr2d)
   endif 
   if (present(int2ptr2d)) then
      allocate(int2ptr2d(dim1Sz,dim2Sz), STAT=allocSuccess0else1)
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, int2ptr2d=int2ptr2d)
   endif
   if (present(int2ptr3d)) then
      allocate(int2ptr3d(dim1Sz,dim2Sz,dim3Sz), STAT=allocSuccess0else1)
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, int2ptr3d=int2ptr3d)
   endif
   if (present(int4ptr1d)) then
      allocate(int4ptr1d(dim1Sz), STAT=allocSuccess0else1)
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, int4ptr1d=int4ptr1d)
   endif
   if (present(int4ptr2d)) then
      allocate(int4ptr2d(dim1Sz,dim2Sz), STAT=allocSuccess0else1)
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, int4ptr2d=int4ptr2d)
   endif
   if (present(real4ptr2d)) then
      allocate(real4ptr2d(dim1Sz,dim2Sz), STAT=allocSuccess0else1)
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, real4ptr2d=real4ptr2d)
   endif
   if (present(real8ptr2d)) then
      allocate(real8ptr2d(dim1Sz,dim2Sz), STAT=allocSuccess0else1)
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) &
         call nullify_ptr(ptrIsNull, real8ptr2d=real8ptr2d)
   endif
   if (present(charNptr1d)) then
      allocate(charNptr1d(dim1Sz), STAT=allocSuccess0else1) 
      if(allocSuccess0else1 /= gALLOCATE_SUCCESS) then
         call nullify_ptr(ptrIsNull, charNptr1d=charNptr1d)
      else
         do i=1, dim1Sz, 1
            allocate(charNptr1d(i)%str, STAT=allocSuccess0else1)

            if(allocSuccess0else1 /= gALLOCATE_SUCCESS) then
            ! This prevents potential partial memory allocations for charNptr1d
               if(i > 1) then
                  do j=i, 1, (0-1)
                     write(charNptr1d(j)%str, '(a)') ''
                     deallocate(charNptr1d(j)%str)
                     nullify(charNptr1d(j)%str)
                  end do
               endif
               call dealloc_arr(ptrIsNull, charNptr1d=charNptr1d) 
               exit
            endif

         end do
      endif
   endif

   ! This reports only success or failure of the last attempted allocation
   ! Therefore this function ideally should not be used for more than one allocation.
   if (allocSuccess0else1 /= gALLOCATE_SUCCESS) then
      alloc_arr = .true.
   else
      alloc_arr = .false.
   endif

 end function alloc_arr


!!! Begin controlled array access section !!!

 function int2Array2dLookup(isNull, nulVal, int2ptr2d, x, y)

   integer*2,    pointer :: int2Array2dLookup

   logical, intent(in)   :: isNull
   integer*2, target     :: nulVal
   integer*2, pointer, dimension(:,:) :: int2ptr2d
   integer*4, intent(in) :: x
   integer*4, intent(in) :: y

   if(isNull .eqv. .true.) then 

      ! Assume: upon attempt to access the data, the data array is null
      ! then we are instead using the user-setting for one number everywhere
      ! i.e., the configuration file -specified value
      ! or whatever value this number has to prior to assigning it the user-setting
      ! i.e., the value given in subroutine startup
      int2Array2dLookup => nulVal ! as-read from the configuration file, for instance

      ! this same pattern follows for other map -or- one-number-everywhere-set-by user
      ! variables such as WHC, LGP, and SOS, etc.  Anywhere we may desire to intercept 
      ! data-lookup upon accessing a quantity.  Pattern applies for other data types as well.

      return
   endif

   int2Array2dLookup => int2ptr2d(x, y)

 end function int2Array2dLookup 


 function int4Array2dLookup(isNull, nulVal, int4ptr2d, x, y)

   integer*4, pointer :: int4Array2dLookup

   logical, intent(in)   :: isNull
   integer*4, target     :: nulVal
   integer*4, pointer, dimension(:,:) :: int4ptr2d
   integer*4, intent(in) :: x
   integer*4, intent(in) :: y

   if(isNull .eqv. .true.) then 
      int4Array2dLookup => nulVal; return
   endif

   int4Array2dLookup => int4ptr2d(x, y)

 end function int4Array2dLookup 


 function lgclArray1dLookup(isNull, nulVal, lgclptr1d, i)

   logical, pointer :: lgclArray1dLookup

   logical, intent(in)   :: isNull
   logical, target       :: nulVal
   logical, pointer, dimension(:) :: lgclptr1d
   integer*4, intent(in) :: i

   if(isNull .eqv. .true.) then
      lgclArray1dLookup => nulVal; return
   endif

   lgclArray1dLookup => lgclptr1d(i)

 end function lgclArray1dLookup


 function int2Array3dLookup(isNull, nulVal, int2ptr3d, x, y, z)

   integer*2,    pointer :: int2Array3dLookup

   logical, intent(in)   :: isNull
   integer*2, target     :: nulVal
   integer*2, pointer, dimension(:,:,:) :: int2ptr3d
   integer*4, intent(in) :: x
   integer*4, intent(in) :: y
   integer*4, intent(in) :: z

   if(isNull .eqv. .true.) then 
      int2Array3dLookup => nulVal; return
   endif

   int2Array3dLookup => int2ptr3d(x, y, z)

 end function int2Array3dLookup 


 function charNarray1dLookup(isNull, nulCrN, nulVal, charNptr1d, i, offset_arg)

   type(charN), pointer  :: charNarray1dLookup

   logical, intent(in)   :: isNull
   type(charN), target   :: nulCrN
   character*1           :: nulVal
   type(charN), pointer, dimension(:) :: charNptr1d
   integer*4, intent(in) :: i
   integer, optional     :: offset_arg

   integer :: offset

   offset = 0
   if(present(offset_arg)) offset = offset_arg
   if(isNull .eqv. .true.) then
      write(nulCrN%str, '(a)') nulVal
      charNarray1dLookup => nulCrN; return
   endif

   charNarray1dLookup => charNptr1d(i + offset)

 end function charNarray1dLookup


 function real4Array2dLookup(isNull, nulVal, real4ptr2d, x, y)

   real*4, pointer :: real4Array2dLookup

   logical, intent(in)   :: isNull
   real*4, target        :: nulVal
   real*4, pointer, dimension(:,:) :: real4ptr2d
   integer*4, intent(in) :: x
   integer*4, intent(in) :: y

   if(isNull .eqv. .true.) then 
      real4Array2dLookup => nulVal; return
   endif

   real4Array2dLookup => real4ptr2d(x, y)

 end function real4Array2dLookup 


 function real8Array2dLookup(isNull, nulVal, real8ptr2d, x, y)

   real*8,       pointer :: real8Array2dLookup

   logical, intent(in)   :: isNull
   real*8, target        :: nulVal
   real*8, pointer, dimension(:,:) :: real8ptr2d
   integer*4, intent(in) :: x
   integer*4, intent(in) :: y

   if(isNull .eqv. .true.) then 
      real8Array2dLookup => nulVal; return
   endif

   real8Array2dLookup => real8ptr2d(x, y)

 end function real8Array2dLookup 

end module geowrsi2_arraymgmt_module
