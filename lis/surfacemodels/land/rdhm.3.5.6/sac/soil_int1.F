C  INTERPOLATE SOIL VARIABLES FROM A NON-REGULAR PROFILE INTO
C  USER DEFINED PROFILE
C
      subroutine SOIL_INT1(ymod,nmod,sdmod,sdint,nint,yint)

#include "config.h"

      real ymod(nmod),yint(nint),sdmod(nmod),sdint(nint)
C
C  YMOD - SIMULATED SOIL PROPERTIES AT MODEL DEFINED PROFILE
C  YINT - SOIL PROPERTY AT DESIRED PROFILE
C  SDMOD - DEFINITION OF MODELING SOIL PROFILE
C  SDINT - DEFINITION OF DESIRED SOIL PROFILE
C  NMOD - NUMBER OF MODELING SOIL LAYERS INCLUDING SURFACE AND BOTTOM DEPTHS
C  NINT - NUMBER OF DESIRED SOIL LAYERS
C
        
C      WRITE(*,*) 'nmod = ', nmod, 'nint = ', nint
      do i=1,nint
C       WRITE(*,*) 'yint(', i, ')=', yint(i)
C       WRITE(*,*) 'sdint(', i, ')=', sdint(i)
       do j=1,nmod
C         WRITE(*,*) 'ymod(', j , ')=', ymod(j)
C         WRITE(*,*) 'sdmod(', j, ')=', sdmod(j)
C   Changed .lt. to .le. to be capable to interpolate if sdint=sdmod
C VK 8/2012        if(sdint(i) .lt. sdmod(j)) goto 1
        if(sdint(i) .le. sdmod(j)) goto 1
       enddo
#ifdef ENABLE_DEBUGINFO
       write(*,*) ' WARNING:: Desired soil layer out of model range:',
     +            ' sdmodmax=',sdmod(nmod),', sdint(i)=',sdint(i)
       write(*,*) ' Number of desired soil layers was changed from',
     +             nint,' to', nint-1
#endif
       nint=nint-1
       if(nint .eq. 0) then
        write(*,*) ' ERROR: No layers selected, nint=0'
        stop 777
       endif
       goto 2
    1  continue
       yint(i)=(ymod(j-1)*(sdmod(j)-sdint(i))+ymod(j)*(sdint(i)-
     +         sdmod(j-1)))/(sdmod(j)-sdmod(j-1))
      enddo
    
    2 continue
      return
      end
