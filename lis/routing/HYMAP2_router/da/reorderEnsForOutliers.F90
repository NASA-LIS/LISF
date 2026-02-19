subroutine reorderEnsForOutliers(nensem, statevec, minvalue)

  implicit none

  integer              :: nensem
  real                 :: statevec(nensem)
  real                 :: minvalue

  real                 :: minvT, maxvT, minvG, maxvG
  integer              :: k
  real                 :: spread_total, spread_good, spread_ratio

  !Ensemble spread (total and with 'good' ensemble members
  minvT = 1E10
  maxvT = -1E10
  minvG = 1E10
  maxvG = -1E10

  do k=1,nensem

     if(statevec(k).lt.minvT) then
        minvT = statevec(k)
     endif
     if(statevec(k).gt.maxvT) then
        maxvT = statevec(k)
     endif

     if(statevec(k).gt.minvalue) then
        if(statevec(k).lt.minvG) then
           minvG = statevec(k)
        endif
        if(statevec(k).gt.maxvG) then
           maxvG = statevec(k)
        endif
     endif
  enddo

  if(minvG.eq.1E10.and.maxvG.eq.-1E10) then
     statevec = minvalue
  else
     spread_total = (maxvT - minvT)
     spread_good  = (maxvG - minvG)

     spread_ratio = spread_good/spread_total

     !rescale the ensemble

     do k=1,nensem-1
        statevec(k) = statevec(nensem) + &
             (statevec(k) - statevec(nensem))*spread_ratio
     enddo
  endif

end subroutine reorderEnsForOutliers

