module ac_kinds

use iso_fortran_env, only: int8, &
                           int16, &
                           int32, &
                           real32, &
                           real64
implicit none


integer, parameter :: sp = real32
    !! single precision real kind
integer, parameter :: dp = real64
    !! double precision real kind
integer, parameter :: intEnum = int8
    !! integer kind for enumerated types

end module ac_kinds
