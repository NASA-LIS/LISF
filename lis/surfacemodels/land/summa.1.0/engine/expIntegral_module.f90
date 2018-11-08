module expIntegral_module
USE nrtype
implicit none
private
public::expint
contains

 ! Numerical recipes routines removed; use code from UEB-Veg

 ! ****************** EXPONENTIAL INTEGRAL FUNCTION *****************************************
 ! From UEB-Veg
 ! Computes the exponential integral function for the given value
 FUNCTION EXPINT (LAI)
 REAL(DP) LAI
 REAL(DP) EXPINT
 REAL(DP) a0,a1,a2,a3,a4,a5,b1,b2,b3,b4
 real(dp),parameter :: verySmall=tiny(1.0_dp)     ! a very small number
 IF (LAI < verySmall)THEN
  EXPINT=1._dp

 ELSEIF (LAI.LE.1.0) THEN
  a0=-.57721566_dp
  a1=.99999193_dp
  a2=-.24991055_dp
  a3=.05519968_dp
  a4=-.00976004_dp
  a5=.00107857_dp

  EXPINT = a0+a1*LAI+a2*LAI**2+a3*LAI**3+a4*LAI**4+a5*LAI**5 - log(LAI)

 ELSE
  a1=8.5733287401_dp
  a2=18.0590169730_dp
  a3=8.6347637343_dp
  a4=.2677737343_dp
  b1=9.5733223454_dp
  b2=25.6329561486_dp
  b3=21.0996530827_dp
  b4=3.9584969228_dp

  EXPINT=(LAI**4+a1*LAI**3+a2*LAI**2+a3*LAI+a4)/ &
      ((LAI**4+b1*LAI**3+b2*LAI**2+b3*LAI+b4)*LAI*exp(LAI))

 END IF
 RETURN
 END FUNCTION EXPINT

END MODULE expIntegral_module
