!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------

!  $Id: DragCoefficients.F90,v 1.7 2007-10-04 20:38:23 dasilva Exp $

module clsmf25_DragCoefficientsMod

!----------------------------------------------------------------------
!BOP

! !MODULE: 

!    DragCoefficientsMod -- Container for the GEOS drag coefficient utility.

! !USES:

use clsmf25_MAPL_Constants

  implicit none
  private

! !PUBLIC MEMBER FUNCTIONS:

  public GETCDS
  public GETCDM
  public GETCDH

! !INTERFACE:
!
  interface GETCDS ! (TVA,UA,DZ,TVS,Z0,ZT,ZQ,  CT,CM,CQ,CN,RI,DCT,MSK)
    module procedure GETCDS_1D
    module procedure GETCDS_2D
  end interface

  interface GETCDM ! (TVA,UA,DZ,TVS,Z0, CM,CN,RI)
    module procedure GETCDM_1D
    module procedure GETCDM_2D
  end interface

  interface GETCDH ! (TVA,UA,DZ,TVS,ZH,ZQ,CN,RI,  CT,CQ,DCT,DCQ)
    module procedure GETCDH_1D
    module procedure GETCDH_2D
  end interface

!EOP

!BOP

! !IROUTINE: GETCDS -- Computes drag coefficients at the surface

! !DESCRIPTION: A utility component used to compute drag coefficients
!               using the Louis scheme. It computes two coefficients,
!  one for momentum , $C_D^u$, and one for heat and moisture, $C_D^h$.
!  It also, optionally, computes the derivative of the drag coefficients
!  with respect to surface temperature and surface humidity.
!
!  Surface-air quantities are subscripted $a$ and surface quantities are
!  sucscripted $s$. The difference of the two is, for example, 
!  $\delta T = T_a - T_s$.
! 
!  The drag coefficients have the form:
!  $$
!  C_D^m = C_D^n f_m({\rm Ri})
!  $$
!  and
!  $$
!  C_D^h = C_D^n f_h({\rm Ri}),
!  $$
!  where
!  $$
!  C_D^n = (\frac{\kappa}{\log{(\zeta)}})^2
!  $$ 
!  is the neutral drag coefficient and $\zeta=\frac{\delta z}{z_0} + 1$.
!  The Von Karman constant, $\kappa$, is 
!  taken as 0.40 and $z_0$ is the surface roughness. The height corresponding
!  to surface-air quantities is taken as one-half the thickness of the 
!  lowest model layer.
!
!  The surface bulk Richardson number, Ri, is defined as 
!  $$
!  {\rm Ri} = \delta T_v \frac{\frac{g}{\delta z}}{T_o*\frac{u_a}{\delta z}^2}
!  $$
!  $$
!  {\rm \,} = \delta T_v \alpha,
!  $$
!  where $T_v=T(1+\epsilon q)$ is the virtual temperature, 
!  $\epsilon=\frac{M_a}{M_w}-1$, $M_a$ and $M_w$ are the molecular weights of
!  dry air and water, $u_a$ is the surface-air wind speed, and $\alpha$ is
!  defined to simplify the calculation. In the code we use $T_v$ in place of
!  $T_o$, but assume that the factor $\alpha$ is constatnt when differentiating
!  the Ricahrdson number with respect to $\delta T_v$.
!  
!  The two universal functions of the Richardson number,  $f_m$ and $f_h$,
!  are taken from Louis et al (1979). For unstable conditions (Ri$\le 0$),
!  they are:
!  $$
!  f_m = (1 - 2b \psi)
!  $$
!  and
!  $$
!  f_h = (1 - 3b \psi),
!  $$
!  where
!  $$
!  \psi = \frac{ {\rm Ri} }{ 1 + 3bcC_D^n\sqrt{-{\rm Ri}\zeta} }.
!  $$
!  
!  For stable condition (Ri$\ge 0$), they are
!  $$
!  f_m = \frac{1}{1.0 + 2b\frac{{\rm Ri}}{\psi}}
!  $$
!  and
!  $$
!  f_h = \frac{1}{1.0 + 3b{\rm Ri}\psi},
!  $$
!  where
!  $$
!  \psi =  \sqrt{1+d{\rm Ri}}.
!  $$
!  As in Louis et al (1979) the parameters appearing in these are taken  
!  as $b = c = d = 5$. 
!  
!  We also require the derivative of $ C_D^h$ with respect to $\delta T_v$.
!  $$
!  \frac{\partial C_D^h}{\partial \delta T_v} = 
!    C_D^n \frac{{\rm d} f_h}{{\rm d} {\rm Ri}} \alpha
!  $$
!  
!  The derivatives of $f_h$ are as follows.
!
!  For unstable conditions:
!  $$
!  \frac{{\rm d} f_h}{{\rm d} {\rm Ri}} = - 3b \frac{{rm d} \psi}{{\rm d} {\rm Ri}},
!  $$
!  where
!  $$
!  \frac{{\rm d} \psi}{{\rm d} {\rm Ri}} = \frac{\psi}{{\rm Ri}} 
!     (1 + \frac{3bcC_D^n\psi\zeta}{2\sqrt{-{\rm Ri}\zeta} })
!  $$
!    
!  For stable conditions:
!  $$
!  \frac{{\rm d} f_h}{{\rm d} {\rm Ri}} = -f_h^2  3b  (\psi + {\rm Ri}\frac{d}{2\psi}),
!  $$
!  
!  
! \bigskip
! \hrulefill
! \bigskip

! !INTERFACE:
!
!  subroutine GETCDS (TVA,UA,DZ,TVS,Z0,LOUIS,  CT,CM,CN,RI,DCT)

!
! !INPUT ARGUMENTS:
!
! \ev
!
! \begin{itemize}
!   \item[]
!\makebox[1in][l]{\bf TVA} \makebox[1in][l]{ (K)}
!                   \parbox[t]{4in}{Virtual surface air temperature}
!   \item[]
!\makebox[1in][l]{\bf UA} \makebox[1in][l]{ (m s$^{-1}$)}
!                   \parbox[t]{4in}{Surface wind speed}
!   \item[]
!\makebox[1in][l]{\bf UU} \makebox[1in][l]{ (m)}
!                   \parbox[t]{4in}{Altitude of UA}
!   \item[]
!\makebox[1in][l]{\bf TVS} \makebox[1in][l]{ (K)}
!                   \parbox[t]{4in}{Virtual surface air temperature}
!   \item[]
!\makebox[1in][l]{\bf Z0} \makebox[1in][l]{ (m)}
!                   \parbox[t]{4in}{Surface roughness}
!   \item[]
!\makebox[1in][l]{\bf LOUIS} \makebox[1in][l]{}
!                   \parbox[t]{4in}{Louis scheme parameter (usually 5)}
! \end{itemize}

! \bv

! !OUTPUT ARGUMENTS:

! \ev

! \begin{itemize}
!   \item[]
!\makebox[1in][l]{\bf CT} \makebox[1in][l]{}
!                   \parbox[t]{4in}{Drag coefficient for heat and scalars}
!   \item[]
!\makebox[1in][l]{\bf CM} \makebox[1in][l]{}
!                   \parbox[t]{4in}{Drag coefficient for momentum}
!   \item[]
!\makebox[1in][l]{\bf CN} \makebox[1in][l]{}
!                   \parbox[t]{4in}{Neutral drag coefficient}
!   \item[]
!\makebox[1in][l]{\bf RI} \makebox[1in][l]{}
!                   \parbox[t]{4in}{Surface Bulk Richardson number }
!   \item[]
!\makebox[1in][l]{\bf DCT} \makebox[1in][l]{(K$^{-1}$)}
!                   \parbox[t]{4in}{$ \frac{\partial C_D^h}{\partial \delta T_v}$}
! \end{itemize}
!

! \bv

!EOP

real, parameter :: LOUIS = 5.
real, parameter :: MAXRI = 2.

contains

#undef DIMS
#define DIMS 1
#include "getcds.code"


#undef DIMS
#define DIMS 2
#include "getcds.code"

#undef DIMS
#define DIMS 1
#include "getcdm.code"


#undef DIMS
#define DIMS 2
#include "getcdm.code"

#undef DIMS
#define DIMS 1
#include "getcdh.code"


#undef DIMS
#define DIMS 2
#include "getcdh.code"



end module clsmf25_DragCoefficientsMod
