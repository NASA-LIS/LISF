C  ****   New calculation of average temperature (TAVG)   **********
C  ****   in freezing/thawing layer using UP, DOWN, and MIDDLE   ***
c  ****   layer temperatures (TUP, TDN, TM)               **********
c!!!  Ken's version and early my version was wrong because all IF
c     statements compared to 0. not T0 temperature. However, results were
c     reasonably good due to a reasonable assumption that TAVG in all 
c     combinations of TUP, TM & TDN equals (TUP+2*TM+TDN)/4.
c     In this version I replaced all IF statements to compare to T0,
c     and corrected one old statement (marked c!!!) that was wrong.

      FUNCTION ST_AVG1(TUP,TM,TDN,DZ)
      
      PARAMETER (T0=273.16)

      DZH=DZ*0.5
      IF(TUP .LT. T0) THEN

        IF(TM .LT. T0) THEN

          IF(TDN .LT. T0) THEN
C1  ******   TUP, TM, TDN < 0.  CASE  ************          
            TAVG=(TUP+2*TM+TDN)/4.
            GOTO 777 
C       C111  ****   END TDN_111  THEN

          ELSE
C2  ******   TUP & TM < 0.,  TDN > 0. CASE   *****
            X0=(T0-TM)*DZH/(TDN-TM)
            TAVG=0.5*(TUP*DZH+TM*(DZH+X0)+T0*(2.*DZH-X0))/DZ
            GOTO 777           
C       C112   ****  END TDN_111  ELSE   ****
          ENDIF      

C   C11  ***********  END TM_1 THEN   ***********
        ELSE
        
          IF(TDN .LT. T0) THEN
C3  *******  TUP < 0.  TM > 0.  TDN < 0. CASE   ****
            XUP=(T0-TUP)*DZH/(TM-TUP)
            XDN=DZH-(T0-TM)*DZH/(TDN-TM)
            TAVG=0.5*(TUP*XUP+T0*(2.*DZ-XUP-XDN)+TDN*XDN)/DZ
            GOTO 777          
C       C121   ****   END TDN_121  THEN  *****

          ELSE
C4   ******  TUP < 0  TM > 0  TDN > 0  CASE   ******
            XUP=(T0-TUP)*DZH/(TM-TUP)
            TAVG=0.5*(TUP*XUP+T0*(2.*DZ-XUP))/DZ
            GOTO 777          

C       C122   ***   END TDN_121  ELSE   ***
          ENDIF   
        
C   C12  ***********  END TM_1 ELSE   ********** 
        ENDIF

C1    ********************  END TUP THEN   ***********    
      ELSE

        IF(TM .LT. T0) THEN

          IF(TDN .LT. T0) THEN
C5   *****  TUP > 0  TM < 0 TDN < 0 CASE   *********
            XUP=DZH-(T0-TUP)*DZH/(TM-TUP)
            TAVG=0.5*(T0*(DZ-XUP)+TM*(DZH+XUP)+TDN*DZH)/DZ
            GOTO 777          
C       C211   ****   END TDN_211  THEN  *****

          ELSE
C6   *****  TUP > 0  TM < 0  TDN > 0 CASE  *********
            XUP=DZH-(T0-TUP)*DZH/(TM-TUP)
            XDN=(T0-TM)*DZH/(TDN-TM)
c!!!  It was wrong statement
c!!!            TAVG=0.5*(T0*(2.*DZ-XUP+XDN)+TM*(XUP+XDN))/DZ
            TAVG=0.5*(T0*(2.*DZ-XUP-XDN)+TM*(XUP+XDN))/DZ
            GOTO 777                       

C       C212   ***   END TDN_211  ELSE   ***
          ENDIF   

C   C21   ************   END TM_2  THEN   ************
        ELSE

          IF(TDN .LT. T0) THEN
C7   *****  TUP > 0  TM > 0  TDN < 0 CASE   ********
            XDN=DZH-(T0-TM)*DZH/(TDN-TM)
            TAVG=(T0*(DZ-XDN)+0.5*(T0+TDN)*XDN)/DZ
            GOTO 777      
C       C221   ****   END TDN_221  THEN  *****

          ELSE
C8   *****  TUP > 0  TM > 0  TDN > 0 CASE   ********
            TAVG=(TUP+2.*TM+TDN)/4.
            GOTO 777          
C       C222   ***   END TDN_221  ELSE   ***
          ENDIF           

C   C22   ************   END TM_2  ELSE   ************
        ENDIF

C   *********************  END TUP ELSE  **************
      ENDIF                      
777   ST_AVG1=TAVG

      RETURN
      END
      