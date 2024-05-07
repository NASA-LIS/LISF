!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
module LIS_numerRecipesMod

!BOP
!
! !MODULE: LIS_numerRecipesMod
!
! !DESCRIPTION:
!  The code in this file contains a few borrowed routines from 
!  the Fortran Numerical Recipes book
!   
! !REVISION HISTORY:
  implicit none
  
  PRIVATE

  public :: LIS_gasdev
  public :: LIS_rand_func
!EOP
  contains

!BOP
! 
! !ROUTINE: LIS_gasdev
! \label{LIS_gasdev}
! 
! !DESCRIPTION:
! 
! Returns a normally distributed deviate with zero mean and unit
! variance, using ran2(idum) as the source of uniform deviates. 
! This routine is adapted from the Numerical Recipies for Fortran 
!
! !REVISION HISTORY: 
!  27 Feb 2005: Sujay Kumar : Specification in LIS.
!
! !INTERFACE:
    FUNCTION LIS_gasdev(idum)
!EOP
      INTEGER idum
      REAL LIS_gasdev
      INTEGER iset
      REAL fac,gset,rsq,v1,v2
      SAVE iset,gset
      DATA iset/0/
      if (idum.lt.0) iset=0
      if (iset.eq.0) then
1        v1=2.*ran2(idum)-1.
         v2=2.*ran2(idum)-1.
         rsq=v1**2+v2**2
         if(rsq.ge.1..or.rsq.eq.0.)goto 1
         fac=sqrt(-2.*log(rsq)/rsq)
         gset=v1*fac
         LIS_gasdev=v2*fac
         iset=1
      else
         LIS_gasdev=gset
         iset=0
      endif
      return
    END FUNCTION LIS_gasdev


!BOP
! 
! !ROUTINE: ran2
! \label{ran2}
!
! !DESCRIPTION: 
! Long period (>2!1e18) random number generator of L'Ecuyer with
! Bays-Durham shuffle and added safeguards. Returns a uniform
! random deviate between 0.0 and 1.0 (exclusive of the endpoint
! values). Call with ``idum'' a negative integer to initialize;
! thereafter, do not alter ``idum'' between successive deviates
! in a sequence. RNMX should approximate the largest floating 
! value that is less than 1. This function is adapted from 
! Numerical Recipies for Fortran. 
! 
! !REVISION HISTORY: 
!  27 Feb 2005: Sujay Kumar : Specification in LIS.
!
! !INTERFACE:
    FUNCTION ran2(idum)
!EOP

      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
           IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
           NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
         idum=max(-idum,1)
         idum2=idum
         do j=NTAB+8,1,-1
            k=idum/IQ1
            idum=IA1*(idum-k*IQ1)-k*IR1
            if (idum.lt.0) idum=idum+IM1
            if (j.le.NTAB) iv(j)=idum
         end do
         iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
    END function ran2
   

  subroutine LIS_rand_func(idum,rand)

!  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to
!  any negative value to initialize or reinitialize the sequence.
!  This function is taken from W.H. Press', ``Numerical Recipes'' p. 199.

    implicit none
    integer :: seed
    integer :: count, count_rate, count_max
    integer :: idum
    real   :: rand
      
    call system_clock(count,count_rate,count_max)
    seed = -count-count_max
    call random_seed(seed)
    call random_number(rand)
#if 0 
    real, parameter :: mbig=4000000.,mseed=1618033.,mz=0.
    real, parameter :: fac=1./4000000.
    
!  According to Knuth, any large mbig, and any smaller (but still large)
!  mseed can be substituted for the above values.
    integer         ::  ma(55)
    integer         ::  iff
    integer         ::  idum
    integer         ::  i,j,k,ii,jj
    real            ::  mk
    integer         ::  inext
    integer         ::  inextp
    real            ::  mj
    real            ::  rand
    
    iff = 0 
    ma = 0 
    
    if (idum.lt.0 .or. iff.eq.0) then
       iff=1
       mj=mseed-float(abs(idum))
       mj=mod(mj,mbig)      
       ma(55)=mj
       mk=1
       do i=1,54
          ii=mod(21*i,55)
          ma(ii)=mk
          mk=mj-mk
          if(mk.lt.mz) mk=mk+mbig
          mj=ma(ii)
       enddo
       do k=1,4
          do i=1,55
             ma(i)=ma(i)-ma(1+mod(i+30,55))
             if(ma(i).lt.mz) ma(i)=ma(i)+mbig
          enddo
       enddo
       inext=0
       inextp=31
       idum=1
    endif
    inext=inext+1
    if(inext.eq.56) inext=1
    inextp=inextp+1
    if(inextp.eq.56) inextp=1
    mj=ma(inext)-ma(inextp)
    if(mj.lt.mz) mj=mj+mbig
    ma(inext)=mj
    rand=mj*fac
    return
#endif
  end subroutine LIS_rand_func


  end module LIS_numerRecipesMod
