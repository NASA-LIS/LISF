!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!
!BOP
! 
! !MODULE: LVT_NumericalRecipesMod
! \label(LVT_NumericalRecipesMod)
!
! !INTERFACE:
module LVT_NumericalRecipesMod
! 
! !USES:   
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the Confidence Interval computations
! 
! !FILES USED:
!
!  !REVISION HISTORY: 
!  27 Feb 2011    Sujay Kumar  Initial Specification
! 
!EOP
!
  implicit none

  private
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  PUBLIC :: ks2d2s
  PUBLIC :: log2

contains

  SUBROUTINE ks2d2s(x1,y1,n1,x2,y2,n2,d,prob) 
    INTEGER n1,n2
    REAL d,prob,x1(n1),x2(n2),y1(n1),y2(n2)
    
! USESpearsn,probks,quadct

!Two-dimensional Kolmogorov-Smirnov test on two samples.
! Given the x and y coordinates of the first sample as n1 values 
! in arrays x1(1:n1) and y1(1:n1), and likewise for the second sample,
! n2 values in arrays x2 and y2, this routine returns the two-dimensional, 
! two- sample K-S statistic as d, and its significance level as prob. 
! Small values of prob show that the two samples are significantly different.
! Note that the test is slightly distribution- dependent, so prob is only 
! an estimate.
    INTEGER j
    REAL d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen!,probks 
    d1=0.0
    
    do j=1,n1 !First, use points in the first sample as origins. 
       call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
       call quadct(x1(j),y1(j),x2,y2,n2,ga,gb,gc,gd) 
       d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
    enddo
    d2=0.0
    do j=1,n2 !Then, use points in the second sample as origins.
       call quadct(x2(j),y2(j),x1,y1,n1,fa,fb,fc,fd)
       call quadct(x2(j),y2(j),x2,y2,n2,ga,gb,gc,gd) 
       d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
    enddo
    d=0.5*(d1+d2) !Average the K-S statistics. 
    sqen=sqrt(float(n1)*float(n2)/float(n1+n2))
    call pearsn(x1,y1,n1,r1,dum,dumm) !Get the linear correlation coefficient for each sample. 
    call pearsn(x2,y2,n2,r2,dum,dumm)
    rr=sqrt(1.0-0.5*(r1**2+r2**2))
!  Estimate the probability using the K-S probability function probks.
    prob=probks(d*sqen/(1.0+rr*(0.25-0.75/sqen)))
    return
  END SUBROUTINE ks2d2s


  SUBROUTINE quadct(x,y,xx,yy,nn,fa,fb,fc,fd) 
    INTEGER nn
    REAL fa,fb,fc,fd,x,y,xx(nn),yy(nn)
    !  Given an origin (x, y), and an array of nn points with 
    !  coordinates xx and yy, count how many of them are in each
    !  quadrant around the origin, and return the normalized frac- tions. 
    !  Quadrants are labeled alphabetically, counterclockwise from the upper 
    !  right. Used by ks2d1s and ks2d2s.
    INTEGER k,na,nb,nc,nd 
    REAL ff
    
    na=0
    nb=0
    nc=0
    nd=0 
    do k=1,nn
       if(yy(k).gt.y) then 
          if(xx(k).gt.x) then 
             na = na + 1
          else
             nb = nb + 1
          endif
       else
          if(xx(k).gt.x) then 
             nd = nd + 1
          else
             nc = nc + 1
          endif
       endif
    enddo
    ff = 1.0/nn
    fa = ff*na
    fb = ff*nb
    fc = ff*nc
    fd = ff*nd
    return
  end SUBROUTINE quadct

  SUBROUTINE pearsn(x,y,n,r,prob,z) 
    INTEGER n
    REAL prob,r,z,x(n),y(n),TINY 
    PARAMETER (TINY=1.e-20)
    ! USES betai
    !Will regularize the unusual case of com- plete correlation.
    !Given two arrays x(1:n) and y(1:n), this routine computes their 
    ! correlation coefficient r (returned as r), the significance level
    ! at which the null hypothesis of zero correlation is disproved 
    !(prob whose small value indicates a significant correlation), 
    !and Fisher’s z (returned as z), whose value can be used in further 
    !statistical tests as described above.
    INTEGER j
    REAL ax,ay,df,sxx,sxy,syy,t,xt,yt!,betai
    ax=0.
    ay=0.
    do j=1,n 
       ax=ax+x(j)
       ay=ay+y(j) 
    enddo
    ax=ax/n
    ay=ay/n
    sxx=0.
    syy=0.
    sxy=0.

    do j=1,n
       xt = x(j) - ax
       yt = y(j) - ay
       sxx = sxx + xt**2
       syy = syy + yt**2
       sxy = sxy + xt*yt
    enddo
    r=sxy/(sqrt(sxx*syy)+TINY) 
    z=0.5*log(((1.+r)+TINY)/((1.-r)+TINY)) 
    df=n-2 
    t=r*sqrt(df/(((1.-r)+TINY)*((1.+r)+TINY))) 
    prob=betai(0.5*df,0.5,df/(df+t**2))
    !prob=erfcc(abs(z*sqrt(n-1.))/1.4142136) 
    return
  END SUBROUTINE pearsn


  FUNCTION probks(alam)

    implicit none

    REAL probks,alam,EPS1,EPS2 
    PARAMETER (EPS1=0.001, EPS2=1.e-8)
    !  Kolmogorov-Smirnov probability function.
    INTEGER j
    REAL a2,fac,term,termbf 
    a2=-2.*alam**2
    fac=2.
    probks=0.
    termbf=0.
    !  Previous term in sum.
    do j=1,100
       term=fac*exp(a2*j**2)
       probks=probks+term 
       if(abs(term).le.EPS1*termbf.or.abs(term).le.EPS2*probks)return
       fac=-fac
       !     Alternating signs in sum.
       !Get here only by failing to converge.
       termbf=abs(term) 
    enddo
    probks=1. 
    return 
  END FUNCTION probks


  FUNCTION betai(a,b,x)
    REAL betai,a,b,x
    ! USESbetacf,gammln
    !Returns the incomplete beta function Ix(a,b).
    REAL bt!,betacf,gammln
    if(x.lt.0..or.x.gt.1.) then
      print*, 'bad argument x in betai'
      stop
    endif
    if(x.eq.0..or.x.eq.1.) then
       bt=0.
    else !Factors in front of the continued fraction.
       bt=exp(gammln(a+b)-gammln(a)-gammln(b) +a*log(x)+b*log(1.-x))
    endif
    if(x.lt.(a+1.)/(a+b+2.))then
       betai=bt*betacf(a,b,x)/a
       return 
    else
       betai=1.-bt*betacf(b,a,1.-x)/b
       return 
    endif
  END FUNCTION betai

  FUNCTION gammln(xx) 
    integer :: j 
    REAL gammln,xx
    !Returns the value ln[Γ(xx)] for xx > 0. INTEGER j
    DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
    !Internal arithmetic will be done in double precision, a nicety that you can omit if five-figure
    !accuracy is good enough.
    SAVE cof,stp
    DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,&
         24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,&
         -.5395239384953d-5,2.5066282746310005d0/ 
    x=xx
    y=x
    tmp=x+5.5d0 
    tmp=(x+0.5d0)*log(tmp)-tmp 
    ser=1.000000000190015d0
    do j=1,6
       y=y+1.d0
       ser=ser+cof(j)/y 
    enddo
    gammln=tmp+log(stp*ser/x) 
    return
  END FUNCTION gammln

  FUNCTION betacf(a,b,x)
    INTEGER MAXIT
    REAL betacf,a,b,x,EPS,FPMIN
    PARAMETER (MAXIT=100,EPS=3.e-7,FPMIN=1.e-30)
    !  Used by betai: Evaluates continued fraction for incomplete beta function by modified
    !Lentz’s method (§5.2). 
    INTEGER m,m2
    REAL aa,c,d,del,h,qab,qam,qap
    qab=a+b !These q’s will be used in factors that occur in the qap=a+1. coefficients (6.4.6).
    qap =a+1
    qam=a-1.
    c=1. !First step of Lentz’s method.
    d=1.-qab*x/qap
    if(abs(d).lt.FPMIN)d=FPMIN
    d=1./d
    h=d
    do m=1,MAXIT
       m2=2*m
       aa=m*(b-m)*x/((qam+m2)*(a+m2)) 
       d=1.+aa*d
       if(abs(d).lt.FPMIN)d=FPMIN
       c=1.+aa/c
       if(abs(c).lt.FPMIN)c=FPMIN
       d=1./d
       h=h*d*c 
       aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
       d=1.+aa*d 
       if(abs(d).lt.FPMIN)d=FPMIN
       c=1.+aa/c 
       if(abs(c).lt.FPMIN)c=FPMIN
       d=1./d
       del=d*c
       h=h*del 
       if(abs(del-1.).lt.EPS)goto 1
    enddo
    print*, 'a or b too big, or MAXIT too small in betacf'
    stop
1   betacf=h
    return 
  END FUNCTION betacf

!BOP 
! 
! !INTERFACE: 
  function log2(x)
! !DESCRIPTION: 
!  computes log with a base of 2
!EOP
    real  :: log2
    real  :: x

    log2= log(x)/log(2.0)

  end function log2

end module LVT_NumericalRecipesMod
