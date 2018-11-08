      subroutine penmanadj(grn_ind,year,month,day,penadj)
      
      parameter (pi = 3.14159)
      
      integer year,month,day,ndays,mdays(12) 
      real grn_ind,penkmax,penkmin,penadj

      data mdays/286,317,346,11,41,72,102,133,164,194,225,255/ 
      
c------------------------------------------------------------
c Estimation of daily adjustment factors to Penman-based PET
c Sin curve is used to interpolate between penkmax & penkmin
c values. Relationships for penkmax & penkmin derived from 
c analysis of OHD climate PET and Empirical Penman PET. 
c Climate index (grn_ind) of available data was 0.26 - 0.63. 
c Therefore, estimates beyond this range should be used with 
c caution.
c------------------------------------------------------------

      if(grn_ind .gt. 0.51) then
       penkmax = 0.725
      else
       penkmax = -0.6071*grn_ind + 0.8718
      endif
      
      if(grn_ind .gt. 0.70) then
       penkmin = 0.35
      else
       penkmin = -0.7372*grn_ind + 0.8659
      endif

cvk just for sensitivity tests
cvk 7148400 basin kmin and kmax parameters:
c      penkmax=0.723
c      penkmin=0.626
cvk 7186000
cvk      penkmax=0.975
cvk      penkmin=0.573
cvk 7231000
cvk      penkmax=0.645
cvk      penkmin=0.566      
cvk ELDO2
cvk      penkmax=0.815
cvk      penkmin=0.389
cvk BLUO2
cvk      penkmax=0.759
cvk      penkmin=0.624
cvk 7147070
cvk      penkmax=0.723
cvk      penkmin=0.492
                   
      ndays = mdays(month)+day
      if(month .eq. 3 .and. day .ge. 21) ndays = day - 20
      penadj = 0.5*(penkmax+penkmin) + sin(2*ndays*pi/366)*
     +                                0.5*(penkmax-penkmin) 
      
      return
      end
        
