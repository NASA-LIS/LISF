!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.5
!
! Copyright (c) 2024 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
MODULE algo_vpol_m
        IMPLICIT NONE
 
      CONTAINS
       FUNCTION algo_vpol_function(X) RESULT(simtbv)
        USE varsio_m
        USE mironov_m

        IMPLICIT NONE

        REAL*4    :: roh, rov, rsh, rsv, exptauh, exptauv, Ah, Av, X
        COMPLEX(4) :: er_r, c_er, er
        REAL*4    :: algo_vpol_output, simtbv, simtbh

        CALL mironov (freq,X,clay,er_r) !freq: frequency in GHZ

        c_er = er_r
        er = SQRT (c_er - SIN (inc*d2r)**2)
        roh = ABS ( (COS (inc*d2r) - er) / (COS (inc*d2r) + er) )**2
        rov = ABS ( (c_er * COS (inc*d2r) - er) / (c_er * COS (inc*d2r) + er) )**2
        rsh = ( (1-Q) * roh + Q * rov) * EXP (-h * COS (inc*d2r) * COS (inc*d2r))
        rsv = ( (1-Q) * rov + Q * roh) * EXP (-h * COS (inc*d2r) * COS (inc*d2r))
        exptauh = EXP (-tau)
        exptauv = EXP (-tau)
        Ah = Ts * (1 - omega) * (1 - exptauh)
        Av = Ts * (1 - omega) * (1 - exptauv)
        simtbh = Ts * (1 - rsh) * exptauh + Ah * (1 + rsh * exptauh)
        simtbv = Ts * (1 - rsv) * exptauv + Av * (1 + rsv * exptauv)
        algo_vpol_output = simtbv

      END FUNCTION algo_vpol_function

        SUBROUTINE algo_vpol (ii,jj,x1,x2,exitstate)
          USE varsio_m
          IMPLICIT NONE
 
          REAL(4), INTENT(IN)                      :: ii, jj
          REAL(4), INTENT(OUT)                     :: x1, x2
          INTEGER(1), INTENT(OUT)                  :: exitstate
          REAL(4)                                  :: NEDT = 2.0
          REAL(4)                                  :: lowerbound = 0.02
          REAL(4)                                  :: upperbound
          REAL(4)                                  :: x
          REAL(4)                                  :: incvsm
          INTEGER(4)                               :: numvsm
          INTEGER(4)                               :: vv, opt
          REAL(4), DIMENSION(:), ALLOCATABLE, SAVE :: vsmvec, tbvvec
          incvsm = 0.01
          upperbound = 1 - bulkdensity/2.65
          numvsm = FLOOR ((upperbound - lowerbound)/incvsm)
          ALLOCATE (vsmvec(numvsm))
          ALLOCATE (tbvvec(numvsm))
          DO vv = 1,numvsm
             vsmvec(vv) = lowerbound + (numvsm-1)*incvsm - (vv-1)*incvsm
             tbvvec(vv) = algo_vpol_function(vsmvec(vv))
          ENDDO
          IF (tbv >= tbvvec(1) - NEDT .AND. tbv <= tbvvec(numvsm) + NEDT) THEN
              IF (tbv < tbvvec(1)) THEN
                  tbv = tbvvec(1)
              ENDIF
              IF (tbv > tbvvec(numvsm)) THEN
                  tbv = tbvvec(numvsm)
              ENDIF
               IF (tbvvec(numvsm) - tbvvec(1) > NEDT) THEN
                   opt=MINLOC(ABS(tbv-tbvvec),1)
                   x = vsmvec(opt)
                   exitstate = 0
               ELSE
                   x = FillValue_float32
                   exitstate = 1
               ENDIF
          ELSEIF (tbv > tbvvec(numvsm) + NEDT) THEN
              IF (topigbptype >= 1 .AND. topigbptype <= 5) THEN
                  x = upperbound
              ELSE
                  x = lowerbound
              ENDIF
              exitstate = 1
          ELSEIF (tbv < tbvvec(1) - NEDT) THEN
              x = upperbound
              exitstate = 1
          ELSE
              x = FillValue_float32
              exitstate = 1
          ENDIF
          x1 = x
          x2 = tau
          DEALLOCATE (vsmvec)
          DEALLOCATE (tbvvec)
 
        END SUBROUTINE algo_vpol
 
      END MODULE algo_vpol_m
