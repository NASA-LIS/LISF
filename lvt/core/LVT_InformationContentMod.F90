!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LVT_InformationContentMod

  implicit  none
  
  PRIVATE   
  
  PUBLIC :: LVT_initInformationContent
  PUBLIC :: LVT_findWordIndex
  PUBLIC :: LVT_findWordTransitionIndex

  PUBLIC :: LVT_ICwordLength

  integer, parameter   :: LVT_ICwordlength=3
  character*3, allocatable :: LVT_ICword(:)
  
contains

  subroutine LVT_initInformationContent
    
    allocate(LVT_ICword(2**LVT_ICwordlength))
  
    if(LVT_ICwordlength==3) then 
!8 possibilities
       LVT_ICword(1) = '000'
       LVT_ICword(2) = '001'
       LVT_ICword(3) = '010'
       LVT_ICword(4) = '011'
       LVT_ICword(5) = '100'
       LVT_ICword(6) = '101'
       LVT_ICword(7) = '110'
       LVT_ICword(8) = '111'
    endif

  end subroutine LVT_initInformationContent


  subroutine LVT_findWordIndex(word,binval)

    character*1  :: word(LVT_ICwordlength)
    integer      :: binval

    character*3  :: word1
    integer      :: l
    
    word1 = ''
    do l=1,LVT_ICwordlength
       word1 = trim(word1)//trim(word(l))
    enddo

    binval = -1
    do l=1,2**LVT_ICwordlength
       if(trim(word1).eq.LVT_ICword(l)) then 
          binval = l
       endif
    enddo

  end subroutine LVT_findWordIndex

  subroutine LVT_findWordTransitionIndex(word,word_prev, col,row)

    character*1  :: word(LVT_ICwordlength)
    character*1  :: word_prev(LVT_ICwordlength)
    integer      :: col,row

    integer      :: binval_curr
    integer      :: binval_prev

    character*3  :: word1
    integer      :: l
    
    word1 = ''
    do l=1,LVT_ICwordlength
       word1 = trim(word1)//trim(word(l))
    enddo

    row = -1
    do l=1,2**LVT_ICwordlength
       if(trim(word1).eq.LVT_ICword(l)) then 
          row = l
       endif
    enddo

    word1 = ''
    do l=1,LVT_ICwordlength
       word1 = trim(word1)//trim(word_prev(l))
    enddo
    
    col = -1
    do l=1,2**LVT_ICwordlength
       if(trim(word1).eq.LVT_ICword(l)) then 
          col = l
       endif
    enddo

!the position of the transition is col, row (prev -> curr)




  end subroutine LVT_findWordTransitionIndex

end module LVT_InformationContentMod
