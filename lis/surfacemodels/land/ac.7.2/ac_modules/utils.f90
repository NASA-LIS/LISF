module ac_utils
!! A place for various small utilities.

use ac_kinds, only: dp, &
                    int8, &
                    int32
use, intrinsic :: iso_c_binding, only: c_f_pointer, &
                                       c_loc, &
                                       c_null_char, &
                                       c_ptr
implicit none


interface roundc
    module procedure roundc_int8
    module procedure roundc_int32
end interface roundc


contains


subroutine assert(condition, message)
    !! Prints an error message if the condition is not met,
    !! and then shuts down the whole program.
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
        print *, 'ABORT: ', message
        stop 1
    end if
end subroutine assert


function GetAquaCropDescription() result(str)
    !! Returns a string with the aquacrop version and release date.
    character(len=:), allocatable :: str

    str = 'AquaCrop ' // GetVersionString() // ' (' // GetReleaseDate() // ')'
end function GetAquaCropDescription


function GetAquaCropDescriptionWithTimeStamp() result(str)
    !! Same as GetAquaCropDescription(), but with a time stamp.
    character(len=:), allocatable :: str
    character(len=10) :: datestr
    character(len=8)  :: timestr

    integer, dimension(8) :: d

    call date_and_time(values=d)

    WRITE(datestr,10)d(3),d(2),d(1)
    WRITE(timestr,8)d(5),d(6),d(7)
 10 FORMAT(I2.2, '-', I2.2, '-', I4.4)
 8  FORMAT(I2.2, ':', I2.2, ':', I2.2)
    str = GetAquaCropDescription() // ' - Output created on (date) : ' // &
          datestr // '   at (time) : ' // timestr
    
end function GetAquaCropDescriptionWithTimeStamp


function GetReleaseDate() result(str)
    !! Returns a string containing the month and year of the release.
    character(len=:), allocatable :: str

    str = 'August 2024'
end function GetReleaseDate


function GetVersionString() result(str)
    !! Returns a string containing the version number (major+minor).
    character(len=:), allocatable :: str

    str = '7.2'
end function GetVersionString


function int2str(num) result(str)
    !! Converts an integer into a string.
    integer, intent(in) :: num
    character(len=:), allocatable :: str

    character(len=128) :: tmp

    write(tmp, '(i0)') num
    str = trim(tmp)
end function int2str


function roundc_int32(x, mold) result(y)
    !! Returns commercial rounds, following Pascal's banker's rules for rounding
    real(dp), intent(in) :: x
        !! Value to be rounded to an integer
    integer(int32), intent(in) :: mold
        !! Integer determining the kind of the integer result
    integer(int32) :: y

   if (abs(x - floor(x, kind=int32) - 0.5_dp) < epsilon(0._dp)) then
       if (x > 0) then
          if (mod(abs(trunc(x)),2) == 0) then
              y = floor(x, kind=int32)
          else
              y = ceiling(x, kind=int32)
          end if
       else
          if (mod(abs(trunc(x)),2) == 0) then
              y = ceiling(x, kind=int32)
          else
              y = floor(x, kind=int32)
          end if
       end if
    else !standard round for values not ending on 0.5
       y = nint(x, kind=int32)
    end if
end function roundc_int32


function roundc_int8(x, mold) result(y)
    !! Returns commercial rounds, following Pascal's banker's rules for rounding
    real(dp), intent(in) :: x
        !! Value to be rounded to an integer
    integer(int8), intent(in) :: mold
        !! Integer determining the kind of the integer result
    integer(int8) :: y

    if (abs(x - floor(x, kind=int32) - 0.5_dp) < epsilon(0._dp)) then
       if (x > 0) then
          if (mod(abs(trunc(x)),2) == 0) then
             y = floor(x, kind=int8)
          else
              y = ceiling(x, kind=int8)
          end if
       else
          if (mod(abs(trunc(x)),2) == 0) then
              y = ceiling(x, kind=int8)
          else
              y = floor(x, kind=int8)
          end if
       end if
    else !standard round for values not ending on 0.5
       y = nint(x, kind=int8)
    end if
end function roundc_int8


function trunc(x) result(y)
    !! Returns the integer part of x, which is always smaller than (or equal to) x
    !! in absolute value.
    real(dp), intent(in) :: x
    integer(int32) :: y

    if (x > 0) then
        y = floor(x, kind=int32)
    else
        y = ceiling(x, kind=int32)
    end if
end function trunc


subroutine upper_case(str)
    !! Converts a string to upper case.
    !!
    !! Adapted from https://fortranwiki.org/fortran/show/m_strings.
    character(len=*), intent(inout) :: str

    integer :: i

    do i = 1, len(str)
        select case(str(i:i))
        case("a":"z")
            str(i:i) = achar(iachar(str(i:i)) - 32)
        end select
    end do
end subroutine upper_case


function pointer2string(c_pointer, strlen) result(string)
    !! Returns a Fortran string from a C-pointer plus the string length.
    type(c_ptr), intent(in) :: c_pointer
        !! C-style pointer
    integer(int32), intent(in) :: strlen
        !! Length of the string
    character(len=strlen) :: string

    character, pointer, dimension(:) :: f_pointer
    integer :: i

    call c_f_pointer(c_pointer, f_pointer, [strlen])

    do i = 1, strlen
        string(i:i) = f_pointer(i)
    end do
end function pointer2string


function string2pointer(string) result(c_pointer)
    !! Returns a C-pointer from a Fortran string.
    character(len=*), intent(in) :: string
    type(c_ptr) :: c_pointer

    character(len=:), allocatable, target, save :: f_string

    f_string = string // c_null_char
    c_pointer = c_loc(f_string)
end function string2pointer


subroutine twostrings2twopointers(string1, string2, c_pointer1, c_pointer2)
    !! Returns two C-pointers for two Fortran strings.
    character(len=*), intent(in) :: string1
    character(len=*), intent(in) :: string2
    type(c_ptr), intent(inout) :: c_pointer1
    type(c_ptr), intent(inout) :: c_pointer2

    character(len=:), allocatable, target, save :: f_string1
    character(len=:), allocatable, target, save :: f_string2

    f_string1 = string1 // c_null_char
    c_pointer1 = c_loc(f_string1)

    f_string2 = string2 // c_null_char
    c_pointer2 = c_loc(f_string2)
end subroutine twostrings2twopointers


subroutine threestrings2threepointers(string1, string2, string3, &
                                    c_pointer1, c_pointer2, c_pointer3)
    !! Returns three C-pointers for three Fortran strings.
    character(len=*), intent(in) :: string1
    character(len=*), intent(in) :: string2
    character(len=*), intent(in) :: string3
    type(c_ptr), intent(inout) :: c_pointer1
    type(c_ptr), intent(inout) :: c_pointer2
    type(c_ptr), intent(inout) :: c_pointer3

    character(len=:), allocatable, target, save :: f_string1
    character(len=:), allocatable, target, save :: f_string2
    character(len=:), allocatable, target, save :: f_string3

    f_string1 = string1 // c_null_char
    c_pointer1 = c_loc(f_string1)

    f_string2 = string2 // c_null_char
    c_pointer2 = c_loc(f_string2)

    f_string3 = string3 // c_null_char
    c_pointer3 = c_loc(f_string3)
end subroutine threestrings2threepointers


subroutine write_file(fhandle, line, advance, iostat)
    !! Writes one line to a file.
    integer, intent(in) :: fhandle
        !! file handle of an already-opened file
    character(len=*), intent(in) :: line
        !! line to write to the file
    logical, intent(in) :: advance
        !! whether or not to append a newline character
    integer, intent(out) :: iostat
        !! IO status returned by write()

    character(len=:), allocatable :: advance_str

    if (advance) then
        advance_str = 'yes'
    else
        advance_str = 'no'
    end if

    write(fhandle, '(a)', advance=advance_str, iostat=iostat) line
end subroutine write_file


subroutine open_file(fhandle, filename, mode, iostat)
    !! Opens a file in the given mode.
    integer, intent(out) :: fhandle
        !! file handle to be used for the open file
    character(len=*), intent(in) :: filename
        !! name of the file to assign the file handle to
    character, intent(in) :: mode
        !! open the file for reading ('r'), writing ('w') or appending ('a')
    integer, intent(out) :: iostat
        !! IO status returned by open()

    logical :: file_exists

    inquire(file=filename, exist=file_exists)

    if (mode == 'r') then
        open(newunit=fhandle, file=trim(filename), status='old', &
             action='read', iostat=iostat)
    elseif (mode == 'a') then
        if (file_exists) then
            open(newunit=fhandle, file=trim(filename), status='old', &
                 position='append', action='write', iostat=iostat)
        else
            open(newunit=fhandle, file=trim(filename), status='new', &
                 action='write', iostat=iostat)
        end if
    elseif (mode == 'w') then
        open(newunit=fhandle, file=trim(filename), status='replace', &
             action='write', iostat=iostat)
    end if
end subroutine open_file


end module ac_utils
