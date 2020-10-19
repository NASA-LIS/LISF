SUBROUTINE DR_HOOK_UTIL(LDHOOK,CDNAME,KCASE,PKEY,CDFILENAME,KSIZEINFO)
USE PARKIND1  ,ONLY : JPIM     ,JPRB    ,JPRD
USE OML_MOD,ONLY : OML_MAX_THREADS,OML_MY_THREAD,OML_INIT
#ifdef SFX_MPI
USE MPL_INIT_MOD, ONLY : MPL_INIT
USE MPL_ARG_MOD, ONLY : MPL_GETARG
#endif
USE YOMGSTATS, ONLY : LAST_KNUM,LAST_KSWITCH,LDETAILED_STATS,MYPROC_STATS, &
                      NHOOK_MESSAGES,TIME_LAST_CALL
USE YOMHOOKSTACK, ONLY : LL_THREAD_FIRST,ISAVE,IMAXSTACK,CSTACK   ! For monitoring thread stack usage
#ifdef NAG
USE F90_UNIX_ENV, ONLY: GETARG
#endif

IMPLICIT NONE
LOGICAL,INTENT(INOUT)       :: LDHOOK
CHARACTER(LEN=*),INTENT(IN) :: CDNAME,CDFILENAME
INTEGER(KIND=JPIM),INTENT(IN) :: KCASE,KSIZEINFO
REAL(KIND=JPRB),INTENT(INOUT) :: PKEY
#ifdef RS6K
INTEGER(KIND=JPIM) :: INEWMASK, IOLDMASK, UMASK
#endif
#ifdef CRAYXT
INTEGER(KIND=JPIM) :: IRET, SETVBUF3F
#endif
LOGICAL,SAVE :: LL_FIRST_TIME = .TRUE.
CHARACTER(LEN=512) :: CLENV
INTEGER(KIND=JPIM) INUMTIDS, IMYTID
LOGICAL :: LLMPI
INTEGER*8 :: MAXMEM=0
INTEGER*8 :: GETMAXMEM
INTEGER*8 GETMAXLOC
LOGICAL :: LLFINDSUMB=.FALSE.
REAL(KIND=JPRD) :: ZCLOCK
REAL(KIND=JPRB) :: ZDIFF
CHARACTER(LEN=7) CLSTR

INTEGER*8 ILOC         ! For monitoring thread stack usage
CHARACTER(LEN=3) CHEAP ! For monitoring heap usage
INTEGER          JHEAP ! For monitoring heap usage
DATA JHEAP/0/

#include "user_clock.h"


! -----------------------------------------------------------------

IF (.NOT.LDHOOK) RETURN

IMYTID = OML_MY_THREAD()
INUMTIDS = OML_MAX_THREADS()
IF (LL_FIRST_TIME) THEN
  LL_FIRST_TIME = .FALSE.
#ifdef CRAYXT
  IRET = SETVBUF3F(0, 1, 0) ! Set unit#0 into line-buffering mode to avoid messy output
#endif
  CALL OML_INIT()
  CALL GET_ENVIRONMENT_VARIABLE('DR_HOOK_NOT_MPI',CLENV)
  IF (CLENV == ' ' .OR. CLENV == '0' .OR. &
    & CLENV == 'false' .OR. CLENV == 'FALSE') THEN
    LLMPI=.TRUE.
#ifdef SFX_MPI
    CALL MPL_INIT(LDINFO=.FALSE.) ! Do not produce any output
#endif
  ELSE
    LLMPI=.FALSE.
  ENDIF
  CALL GET_ENVIRONMENT_VARIABLE('DR_HOOK',CLENV)
  IF (CLENV == ' ' .OR. CLENV == '0' .OR. &
    & CLENV == 'false' .OR. CLENV == 'FALSE') THEN
    LDHOOK = .FALSE.
    CALL C_DRHOOK_SET_LHOOK(0)
  ENDIF
  IF (LLMPI) THEN
#ifdef SFX_MPI
    CALL MPL_GETARG(0, CLENV)  ! Get executable name & also propagate args
#endif
  ELSE
    CALL GETARG(0, CLENV)
  ENDIF
  IF (.NOT.LDHOOK) RETURN
  
  CALL C_DRHOOK_INIT(CLENV, INUMTIDS)

!JFH---Initialisation to monitor stack usage by threads-------------
  CALL GET_ENVIRONMENT_VARIABLE('DR_HOOK_STACKCHECK',CSTACK)
  IF (CSTACK == 'yes' .OR. CSTACK == 'YES' ) THEN
    IF(IMYTID == 1 ) THEN
      ALLOCATE(LL_THREAD_FIRST(INUMTIDS))
      ALLOCATE(ISAVE(INUMTIDS))
      ALLOCATE(IMAXSTACK(INUMTIDS))
      LL_THREAD_FIRST=.TRUE.
      ISAVE=0
      IMAXSTACK=0
    ENDIF
  ENDIF
!JFH------------ End ---------------------------------------------
!JFH---Initialisation to monitor heap usage-----------------------
  JHEAP=0
  CALL GET_ENVIRONMENT_VARIABLE('DR_HOOK_HEAPCHECK',CHEAP)
  IF (CHEAP == 'yes' .OR. CHEAP == 'YES' ) JHEAP=1
  IF (CHEAP == 'trb' .OR. CHEAP == 'TRB' ) JHEAP=2
  IF(IMYTID == 1) THEN
    IF(JHEAP>0) THEN
!     write(0,*) "HEAPCHECK=",CHEAP,JHEAP
      CALL SETHEAPCHECK()
    ENDIF
  ENDIF
!JFH------------ End ---------------------------------------------

ENDIF

!JFH---Code to monitor stack usage by threads---------------------
#ifndef NAG
IF (CSTACK == 'yes' .or. CSTACK == 'YES' ) THEN
  IF(IMYTID > 1) THEN
    IF(LL_THREAD_FIRST(IMYTID))THEN 
      LL_THREAD_FIRST(IMYTID)=.FALSE.
      ISAVE(IMYTID)=LOC(LLMPI)
    ENDIF
    ILOC=LOC(LLMPI)
    IF(ISAVE(IMYTID)-ILOC > IMAXSTACK(IMYTID)) THEN
      IMAXSTACK(IMYTID)=ISAVE(IMYTID)-ILOC
      WRITE(0,'(A,I3,A,I12,2X,A)')"STACKCHECK Max stack usage by thread",imytid," =",IMAXSTACK(IMYTID),CDNAME
    ENDIF
  ENDIF
ENDIF
#endif
!JFH------------ End ---------------------------------------------

IF (KCASE == 0) THEN
  CALL C_DRHOOK_START(CDNAME, IMYTID, PKEY, CDFILENAME, KSIZEINFO)
!JFH---Code to monitor heap usage -------------------------
  IF(IMYTID == 1 .AND. MYPROC_STATS == 1 .AND. JHEAP>0) THEN
    GETMAXMEM=GETMAXLOC()
    IF(GETMAXMEM .GT. MAXMEM) THEN
      MAXMEM = GETMAXMEM
      WRITE(0,*) "HEAPCHECK Max heap at beg of routine =",MAXMEM," ",CDNAME
#ifdef RS6K
      IF(JHEAP == 2) CALL XL__TRBK()
#endif
    ENDIF
  ENDIF
!JFH------------ End ---------------------------------------------
ELSE IF (KCASE == 1) THEN
!JFH---Code to monitor heap usage -------------------------
  IF(IMYTID == 1 .AND. MYPROC_STATS == 1 .AND. JHEAP>0) THEN
    GETMAXMEM=GETMAXLOC()
    IF(GETMAXMEM .GT. MAXMEM) THEN
      MAXMEM = GETMAXMEM
      WRITE(0,*) "HEAPCHECK Max heap at end of routine =",MAXMEM," ",CDNAME
#ifdef RS6K
      IF(JHEAP == 2) CALL XL__TRBK()
#endif
    ENDIF
  ENDIF
!JFH------------ End ---------------------------------------------
  CALL C_DRHOOK_END  (CDNAME, IMYTID, PKEY, CDFILENAME, KSIZEINFO)
ENDIF

!GM---Code to find gstats SUMB time-------------------------------
IF( LDETAILED_STATS .AND. LLFINDSUMB )THEN
  IF( IMYTID==1 .AND. LAST_KNUM>=500 .AND. MYPROC_STATS <= 2 )THEN
    IF( LAST_KSWITCH==1 .OR. LAST_KSWITCH==2 )THEN
      CALL USER_CLOCK(PELAPSED_TIME=ZCLOCK)
      ZDIFF=ZCLOCK-TIME_LAST_CALL
      IF( ZDIFF > 0.1_JPRB )THEN
        IF( KCASE == 0 )THEN
          CLSTR='ENTERED'
        ELSE
          CLSTR='EXITED'
        ENDIF
        IF( NHOOK_MESSAGES < 100000 )THEN
          WRITE(0,'("DR_HOOK_UTIL: ",A,2X,A," TIMESUMB=",F10.6)')CDNAME,CLSTR,ZDIFF
          NHOOK_MESSAGES=NHOOK_MESSAGES+1
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDIF
!GM------------ End ---------------------------------------------

END SUBROUTINE DR_HOOK_UTIL
