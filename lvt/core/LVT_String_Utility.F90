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
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
MODULE LVT_String_Utility 
   IMPLICIT NONE 
   PRIVATE 
   PUBLIC :: LVT_StrUpCase 
   PUBLIC :: LVT_StrLowCase 
   CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz' 
   CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
CONTAINS 
   FUNCTION LVT_StrUpCase ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN ) :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n 
     ! -- Copy input string 
     Output_String = Input_String 
     ! -- Loop over string elements 
     DO i = 1, LEN( Output_String ) 
       ! -- Find location of letter in lower case constant string 
       n = INDEX( LOWER_CASE, Output_String( i:i ) ) 
       ! -- If current substring is a lower case letter, make it upper case 
       IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n ) 
     END DO 
   END FUNCTION LVT_StrUpCase
   FUNCTION LVT_StrLowCase ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN ) :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n 
     ! -- Copy input string 
     Output_String = Input_String 
     ! -- Loop over string elements 
     DO i = 1, LEN( Output_String ) 
       ! -- Find location of letter in upper case constant string 
       n = INDEX( UPPER_CASE, Output_String( i:i ) ) 
       ! -- If current substring is an upper case letter, make it lower case 
       IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n ) 
     END DO 
   END FUNCTION LVT_StrLowCase
 END MODULE LVT_String_Utility
  
