!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------------
! NASA GSFC Land surface Verification Toolkit (LVT) V1.0
!-------------------------END NOTICE -- DO NOT EDIT-----------------------------
!
!BOP
! 
! !MODULE: LVT_KmeansMod
! \label(LVT_KmeansMod)
!
! !INTERFACE:
module LVT_KmeansMod
! 
! !USES:   
  use ESMF
  use LVT_coreMod
  use LVT_histDataMod
  use LVT_statsDataMod
  use LVT_historyMod
  use LVT_TSMod
  use LVT_logMod
  use LVT_CIMod

  implicit none

!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This module handles the computations required to perform K-means
!  clustering on a given time series. K-means clustering is a type of
!  unsupervised learning, used to find groups within a dataset, with 
!  the number of groups represented by the variable K. The algorithm 
!  works iteratively to assign each data point to one of the K groups
!  based on the features that are provide. LVT computes the clusters based
!  on the mean and standard deviation of each variable. 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
!  1 Nov 2018    Sujay Kumar  Initial Specification
! 
!EOP
!BOP
!-----------------------------------------------------------------------------
! !PUBLIC MEMBER FUNCTIONS:
!-----------------------------------------------------------------------------
  public :: LVT_initKmeans
  public :: LVT_diagnoseKmeans
  public :: LVT_computeKmeans
  public :: LVT_writeMetric_Kmeans
  public :: LVT_resetMetric_Kmeans
  public :: LVT_writerestart_Kmeans
  public :: LVT_readrestart_Kmeans
!EOP


  type, public :: kmeansdec
     integer          :: cluster_num
     integer          :: it_max
     integer          :: it_num
  end type kmeansdec

  type(kmeansdec), save :: LVT_kmeans_struc
  
  private

contains
  subroutine LVT_initKmeans(selectNlevs, stats,metric)
! !ARGUMENTS: 
    integer                 :: selectNlevs(LVT_rc%nDataStreams)
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!
! !DESCRIPTION: 
!  This subroutine initializes the variables required for the 
!  K-means clustering calculations. 
! 
!EOP
    integer                 :: m
    integer                 :: cluster_num
    integer                 :: rc


    allocate(stats%kmeans(LVT_rc%nensem))

    call ESMF_ConfigGetAttribute(LVT_config,&
         cluster_num,&
         label="Number of clusters to use in K-means calculation:",rc=rc)
    call LVT_verify(rc,"Number of clusters to use in K-means calculation: option not specified in the config file")
    
    LVT_kmeans_struc%cluster_num = cluster_num
    LVT_kmeans_struc%it_max = 20

    do m=1,LVT_rc%nensem
       if(metric%selectOpt.eq.1) then 
          allocate(stats%kmeans(m)%model_value_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%kmeans(m)%model_value_total = 0.0
          allocate(stats%kmeans(m)%count_model_value_total(LVT_rc%ngrid,&
               selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%kmeans(m)%count_model_value_total = 0
          
          allocate(stats%kmeans(m)%model_value_total_sxsx(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))
          stats%kmeans(m)%model_value_total_sxsx = 0.0

          allocate(stats%kmeans(m)%model_stdev_total(LVT_rc%ngrid, selectNlevs(1), &
               LVT_rc%strat_nlevels))

          allocate(stats%kmeans(m)%model_mean_cluster(LVT_rc%ngrid,&
               selectNlevs(1), &
               LVT_rc%strat_nlevels))

          allocate(stats%kmeans(m)%model_stdev_cluster(LVT_rc%ngrid,&
               selectNlevs(1), &
               LVT_rc%strat_nlevels))

          allocate(stats%kmeans(m)%model_mean_cluster_se(LVT_rc%ngrid,&
               selectNlevs(1),&
               LVT_rc%strat_nlevels))

          allocate(stats%kmeans(m)%model_stdev_cluster_se(LVT_rc%ngrid,&
               selectNlevs(1),&
               LVT_rc%strat_nlevels))
          
          stats%kmeans(m)%model_mean_cluster_se = 0 
          stats%kmeans(m)%model_stdev_cluster_se = 0 

          if(LVT_rc%obssource(2).ne."none") then 
             if(selectNlevs(2).ge.1) then 
                allocate(stats%kmeans(m)%obs_value_total(LVT_rc%ngrid, &
                     selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%kmeans(m)%obs_value_total = 0.0
                allocate(stats%kmeans(m)%count_obs_value_total(LVT_rc%ngrid, &
                     selectNlevs(2), &
                     LVT_rc%strat_nlevels))
                stats%kmeans(m)%count_obs_value_total = 0
                
                allocate(stats%kmeans(m)%obs_value_total_sxsx(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                stats%kmeans(m)%obs_value_total_sxsx = 0.0

                allocate(stats%kmeans(m)%obs_stdev_total(LVT_rc%ngrid, selectNlevs(1), &
                     LVT_rc%strat_nlevels))
                
                allocate(stats%kmeans(m)%obs_mean_cluster(LVT_rc%ngrid,&
                     selectNlevs(1), &
                     LVT_rc%strat_nlevels))

                allocate(stats%kmeans(m)%obs_stdev_cluster(LVT_rc%ngrid,&
                     selectNlevs(1), &
                     LVT_rc%strat_nlevels))

                allocate(stats%kmeans(M)%obs_mean_cluster_se(LVT_rc%ngrid,&
                     selectNlevs(1),&
                     LVT_rc%strat_nlevels))
                
                allocate(stats%kmeans(M)%obs_stdev_cluster_se(LVT_rc%ngrid,&
                     selectNlevs(1),&
                     LVT_rc%strat_nlevels))

                stats%kmeans(m)%obs_mean_cluster_se = 0 
                stats%kmeans(m)%obs_stdev_cluster_se = 0 
             endif
          endif
       end if
    enddo
!-------------------------------------------------------------------------
! Number of passes required to compute the metric
!-------------------------------------------------------------------------

    metric%npass = 1    
    if(LVT_rc%obssource(2).ne."none") then 
       metric%obsData = .true. 
    else
       metric%obsData = .false. 
    endif
    metric%stdevFlag = .false. 

    metric%customNames = .true. 
    metric%nFields     = 4 

    allocate(metric%mName(metric%nFields))
    metric%mName(1) = "kmeans_mean"
    metric%mName(2) = "kmeans_stdev"
    metric%mName(3) = "kmeans_mean_se" !squared error mean
    metric%mName(4) = "kmeans_stdev_se" !squared error mean

  end subroutine LVT_initKmeans

!BOP
! 
! !ROUTINE: LVT_diagnoseKmeans
! \label{LVT_diagnoseKmeans}
!
! !INTERFACE: 
  subroutine LVT_diagnoseKmeans(pass)
! 
! !USES:     

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to update the computations for 
!   the k-means calculations
!
!   The methods invoked are: 
!   \begin{description}
!    \item[diagnoseSingleModelKmeans](\ref{diagnoseSingleModelKmeans})
!     updates the k-means computation for a single variable 
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer       :: pass
!EOP
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats

    if(pass.eq.LVT_metrics%kmeans%npass) then
       if(LVT_metrics%kmeans%selectOpt.eq.1.or.&
            LVT_metrics%kmeans%timeOpt.eq.1) then 

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             call diagnoseSingleKmeans(model,obs,stats,&
                  LVT_metrics%kmeans)
             model => model%next
             obs => obs%next
             stats => stats%next
          enddo
       endif
    endif
  end subroutine LVT_diagnoseKmeans

!BOP
! 
! !ROUTINE: diagnoseSingleKmeans
! \label{diagnoseSingleKmeans}
!
! !INTERFACE: 
  subroutine diagnoseSingleKmeans(model, obs, stats,metric)
! 
! !USES:   
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This routine updates the k-means computation of the 
!   specified variable. 
!
!  The arguments are: 
!
!  \begin{description}
!   \item[model] model variable object
!   \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    type(LVT_metaDataEntry) :: model
    type(LVT_metadataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP
    integer    :: t,k,n,m,m_k,o_k,tind

    if(stats%selectOpt.eq.1.and.&
         model%selectNlevs.ge.1) then   
       do t=1,LVT_rc%ngrid
          do k=1,model%selectNlevs
             do m=1,LVT_rc%nensem
                m_k = k+model%startNlevs -1
                if(model%count(t,m,m_k).gt.0) then
                   if(metric%selectOpt.eq.1) then 
                      if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                         stats%kmeans(m)%model_value_total(t,k,1) = &
                              stats%kmeans(m)%model_value_total(t,k,1) + &
                              model%value(t,m,m_k)

                         stats%kmeans(m)%model_value_total_sxsx(t,k,1) = &
                              stats%kmeans(m)%model_value_total_sxsx(t,k,1) + &
                              model%value(t,m,m_k)*model%value(t,m,m_k)

                         stats%kmeans(m)%count_model_value_total(t,k,1) = &
                              stats%kmeans(m)%count_model_value_total(t,k,1) + 1
                      endif
                   endif
                   
                   if(LVT_rc%strat_nlevels.gt.1) then 
                      if(LVT_stats%strat_var(t,m,k).gt.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then  
                               stats%kmeans(m)%model_value_total(t,k,2) = &
                                    stats%kmeans(m)%model_value_total(t,k,2) + &
                                    model%value(t,m,m_k)

                               stats%kmeans(m)%model_value_total_sxsx(t,k,2) = &
                                    stats%kmeans(m)%model_value_total_sxsx(t,k,2) + &
                                    model%value(t,m,m_k)*model%value(t,m,m_k)

                               stats%kmeans(m)%count_model_value_total(t,k,2) = &
                                    stats%kmeans(m)%count_model_value_total(t,k,2) + 1
                            endif
                         endif
                      elseif(LVT_stats%strat_var(t,m,k).le.&
                           LVT_rc%strat_var_threshold) then
                         if(metric%selectOpt.eq.1) then 
                            if(model%value(t,m,m_k).ne.LVT_rc%udef) then 
                               stats%kmeans(m)%model_value_total(t,k,3) = &
                                    stats%kmeans(m)%model_value_total(t,k,3) + &
                                    model%value(t,m,m_k)

                               stats%kmeans(m)%model_value_total_sxsx(t,k,3) = &
                                    stats%kmeans(m)%model_value_total_sxsx(t,k,3) + &
                                    model%value(t,m,m_k)*model%value(t,m,m_k)

                               stats%kmeans(m)%count_model_value_total(t,k,3) = &
                                    stats%kmeans(m)%count_model_value_total(t,k,3) + 1
                            endif
                         endif
                      endif
                   endif
                endif
             end do
          enddo
          do k=1,obs%selectNlevs
             do m=1,LVT_rc%nensem
                o_k = k+obs%startNlevs -1
                if(LVT_rc%obssource(2).ne."none") then
                   if(obs%selectNlevs.ge.1) then
                      if(obs%count(t,m,o_k).gt.0) then 
                         if(metric%selectOpt.eq.1.and.&
                              obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                            stats%kmeans(m)%obs_value_total(t,k,1) = &
                                 stats%kmeans(m)%obs_value_total(t,k,1) + &
                                 obs%value(t,m,o_k)

                            stats%kmeans(m)%obs_value_total_sxsx(t,k,1) = &
                                 stats%kmeans(m)%obs_value_total_sxsx(t,k,1) + &
                                 obs%value(t,m,o_k)*obs%value(t,m,o_k)

                            stats%kmeans(m)%count_obs_value_total(t,k,1) = &
                                 stats%kmeans(m)%count_obs_value_total(t,k,1) + 1
                         endif
                         
                         if(LVT_rc%strat_nlevels.gt.1) then 
                            if(LVT_stats%strat_var(t,m,k).gt.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1 &
                                    .and.obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%kmeans(m)%obs_value_total(t,k,2) = &
                                       stats%kmeans(m)%obs_value_total(t,k,2) + &
                                       obs%value(t,m,o_k)


                                  stats%kmeans(m)%obs_value_total_sxsx(t,k,2) = &
                                       stats%kmeans(m)%obs_value_total_sxsx(t,k,2) + &
                                       obs%value(t,m,o_k)*obs%value(t,m,o_k)
                                  

                                  stats%kmeans(m)%count_obs_value_total(t,k,2) = &
                                       stats%kmeans(m)%count_obs_value_total(t,k,2) + 1
                               endif
                            elseif(LVT_stats%strat_var(t,m,k).le.&
                                 LVT_rc%strat_var_threshold) then
                               if(metric%selectOpt.eq.1.and.&
                                    obs%value(t,m,o_k).ne.LVT_rc%udef) then 
                                  stats%kmeans(m)%obs_value_total(t,k,3) = &
                                       stats%kmeans(m)%obs_value_total(t,k,3) + &
                                       obs%value(t,m,o_k)


                                  stats%kmeans(m)%obs_value_total_sxsx(t,k,3) = &
                                       stats%kmeans(m)%obs_value_total_sxsx(t,k,3) + &
                                       obs%value(t,m,o_k)*obs%value(t,m,o_k)

                                  stats%kmeans(m)%count_obs_value_total(t,k,3) = &
                                       stats%kmeans(m)%count_obs_value_total(t,k,3) + 1
                               endif
                               
                               
                            endif
                         endif
                      endif
                   endif
                endif
             enddo
          enddo
       enddo
    endif

  end subroutine diagnoseSingleKmeans

!BOP
! 
! !ROUTINE: LVT_computeKmeans
! \label{LVT_computeKmeans}
!
! !INTERFACE: 
  subroutine LVT_computeKmeans(pass,alarm)
! 
! !USES:   

    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!   This subroutine issues the calls to perform the k-means calculation of
!   desired variables.
!
!   The methods invoked are: 
!   \begin{description}
!    \item[computeSingleModelKmeans](\ref{computeSingleModelKmeans})
!     computes the k-means calculations for a single variable
!   \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
    integer               :: pass
    logical               :: alarm
!EOP
    integer     :: i,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats
    
    if(pass.eq.LVT_metrics%kmeans%npass) then
       if(LVT_metrics%kmeans%selectOpt.eq.1.or.&
            LVT_metrics%kmeans%timeOpt.eq.1) then 
          if(alarm) then 
             if(LVT_metrics%kmeans%timeOpt.eq.1.and.&
                  LVT_metrics%kmeans%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%kmeans%ftn_ts_loc(i,m),200,advance='no') &
                              LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                              LVT_rc%hr,'',LVT_rc%mn, '' 
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%kmeans%ftn_ts_loc(i,1),200,advance='no') &
                           LVT_rc%yr, '',LVT_rc%mo, '', LVT_rc%da, '', &
                           LVT_rc%hr,'',LVT_rc%mn, '' 
                   enddo
                endif
             endif
          endif
             
200       format(I4, a1, I2.2, a1, I2.2, a1, I2.2, a1, I2.2,a1)

          call LVT_getDataStream1Ptr(model)
          call LVT_getDataStream2Ptr(obs)
          call LVT_getstatsEntryPtr(stats)

          do while(associated(model))
             
             call computeSingleKmeans(alarm,model,obs,stats,&
                  LVT_metrics%kmeans)
             
             model => model%next
             obs   => obs%next
             stats => stats%next
          enddo
          
          if(alarm) then 
             if(LVT_metrics%kmeans%timeOpt.eq.1.and.&
                  LVT_metrics%kmeans%extractTS.eq.1) then 
                if(LVT_rc%lvt_wopt.eq."2d ensemble gridspace") then 
                   do m=1,LVT_rc%nensem
                      do i=1,LVT_rc%ntslocs
                         write(LVT_metrics%kmeans%ftn_ts_loc(i,m),fmt='(a1)') ''
                      enddo
                   enddo
                else
                   do i=1,LVT_rc%ntslocs
                      write(LVT_metrics%kmeans%ftn_ts_loc(i,1),fmt='(a1)') ''
                   enddo
                endif
             endif
          endif
       endif
    endif
  end subroutine LVT_computeKmeans
  

!BOP
! 
! !ROUTINE: computeSingleKmeans
! \label{computeSingleKmeans}
!
! !INTERFACE: 
  subroutine computeSingleKmeans(alarm,model,obs,stats,metric)
! 
! !USES:   
        
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
!  This routine performs the k-means calculations. 
!
!  The arguments are: 
!
!  \begin{description}
!    \item[model] model variable object
!    \item[stats] object to hold the updated statistics
!  \end{description}
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP
!BOP
! !ARGUMENTS: 
    logical                 :: alarm
    type(LVT_metaDataEntry) :: model
    type(LVT_metaDataEntry) :: obs
    type(LVT_statsEntry)    :: stats
    type(LVT_metricEntry)   :: metric
!EOP

    integer  :: t,l,k,m,i
    real,    allocatable :: model_value_avg(:,:,:)
    real,    allocatable :: obs_value_avg(:,:,:)   
    integer              :: seed
    real                 :: cluster_center(1,LVT_kmeans_struc%cluster_num)
    real                 :: cluster_energy(LVT_kmeans_struc%cluster_num)
    integer              :: cluster_population(LVT_kmeans_struc%cluster_num)

    real                 :: term1, term2
    real,    allocatable :: obs_value(:)
    integer, allocatable :: obs_cluster(:)
    integer              :: count_obs(LVT_rc%nensem, obs%selectNlevs)

    if(LVT_rc%endtime.eq.1.and.metric%selectOpt.eq.1) then 
       if(stats%selectOpt.eq.1.and.&
            model%selectNlevs.ge.1) then 
          
          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%kmeans(m)%count_model_value_total(t,k,l).gt.&
                           LVT_rc%obsCountThreshold) then 

                         stats%kmeans(m)%model_value_total(t,k,l) = &
                              stats%kmeans(m)%model_value_total(t,k,l)/&
                              stats%kmeans(m)%count_model_value_total(t,k,l)

                         term1 = (stats%kmeans(m)%model_value_total_sxsx(t,k,l)/&
                              stats%kmeans(m)%count_model_value_total(t,k,l)) - & 
                              (stats%kmeans(m)%model_value_total(t,k,l)/&
                              stats%kmeans(m)%count_model_value_total(t,k,l))**2

                         if(term1.gt.0) then 
                            stats%kmeans(m)%model_stdev_total(t,k,l) = &
                                 sqrt(term1)
                         else
                            stats%kmeans(m)%model_stdev_total(t,k,l) = LVT_rc%udef
                         endif

                      else
                         stats%kmeans(m)%model_value_total(t,k,l) = LVT_rc%udef
                         stats%kmeans(m)%model_stdev_total(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
                do k=1,obs%selectNlevs
                   if(LVT_rc%obssource(2).ne."none") then 
                      if(obs%selectNlevs.ge.1) then
                         do l=1,LVT_rc%strat_nlevels                      
                            if(stats%kmeans(m)%count_obs_value_total(t,k,l).gt.&
                                 LVT_rc%obsCountThreshold) then 
                               stats%kmeans(m)%obs_value_total(t,k,l) = &
                                    stats%kmeans(m)%obs_value_total(t,k,l)/&
                                    stats%kmeans(m)%count_obs_value_total(t,k,l)       

                               term1 = (stats%kmeans(m)%obs_value_total_sxsx(t,k,l)/&
                                    stats%kmeans(m)%count_obs_value_total(t,k,l)) - & 
                                    (stats%kmeans(m)%obs_value_total(t,k,l)/&
                                    stats%kmeans(m)%count_obs_value_total(t,k,l))**2
                               
                               if(term1.gt.0) then 
                                  stats%kmeans(m)%obs_stdev_total(t,k,l) = &
                                       sqrt(term1)
                               else
                                  stats%kmeans(m)%obs_stdev_total(t,k,l) = LVT_rc%udef
                               endif
                               
                            else
                               stats%kmeans(m)%obs_value_total(t,k,l) = LVT_rc%udef
                               stats%kmeans(m)%obs_stdev_total(t,k,l) = LVT_rc%udef
                            endif
                         enddo
                      endif
                   endif
                enddo
             enddo
          enddo

          do m=1,LVT_rc%nensem
             do k=1,model%selectNlevs
                do l=1,LVT_rc%strat_nlevels

                   seed = 123456789

                   call cluster_initialize( 1, &
                        LVT_rc%ngrid,&
                        LVT_kmeans_struc%cluster_num,&
                        stats%kmeans(m)%model_value_total(:,k,l),&
                        seed, &
                        cluster_center )

                   call kmeans_01(1,&
                        LVT_rc%ngrid,&
                        LVT_kmeans_struc%cluster_num,&
                        LVT_kmeans_struc%it_max,&
                        LVT_kmeans_struc%it_num,&
                        stats%kmeans(m)%model_value_total(:,k,l),&
                        stats%kmeans(m)%model_mean_cluster(:,k,l),&
                        cluster_center, &
                        cluster_population,&
                        cluster_energy)

                   do t=1,LVT_rc%ngrid
                      if(stats%kmeans(m)%count_model_value_total(t,k,l).ne.0.and.&
                           stats%kmeans(m)%model_mean_cluster(t,k,l).ne.0) then
                         stats%kmeans(m)%model_mean_cluster_se(t,k,l) = &
                              (stats%kmeans(m)%model_value_total(t,k,l) - &
                              cluster_center(1,&
                              stats%kmeans(m)%model_mean_cluster(t,k,l)))**2
                      else
                         stats%kmeans(m)%model_mean_cluster_se(t,k,l) = LVT_rc%udef
                      endif
                   enddo

                   seed = 123456789

                   call cluster_initialize( 1, &
                        LVT_rc%ngrid,&
                        LVT_kmeans_struc%cluster_num,&
                        stats%kmeans(m)%model_stdev_total(:,k,l),&
                        seed, &
                        cluster_center )

                   call kmeans_01(1,&
                        LVT_rc%ngrid,&
                        LVT_kmeans_struc%cluster_num,&
                        LVT_kmeans_struc%it_max,&
                        LVT_kmeans_struc%it_num,&
                        stats%kmeans(m)%model_stdev_total(:,k,l),&
                        stats%kmeans(m)%model_stdev_cluster(:,k,l),&
                        cluster_center, &
                        cluster_population,&
                        cluster_energy)

                   do t=1,LVT_rc%ngrid
                      if(stats%kmeans(m)%count_model_value_total(t,k,l).ne.0.and.&
                           stats%kmeans(m)%model_stdev_cluster(t,k,l).ne.0) then
                         stats%kmeans(m)%model_stdev_cluster_se(t,k,l) = &
                              (stats%kmeans(m)%model_stdev_total(t,k,l) - &
                              cluster_center(1,&
                              stats%kmeans(m)%model_stdev_cluster(t,k,l)))**2
                      else
                         stats%kmeans(m)%model_stdev_cluster_se(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo
          enddo
          if(LVT_rc%obssource(2).ne."none") then 
             do m=1,LVT_rc%nensem
                do k=1,obs%selectNlevs
                   count_obs(m,k) = 0
                   do t=1,LVT_rc%ngrid
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%kmeans(m)%obs_value_total(t,k,l).ne.LVT_rc%udef) then 
                            count_obs(m,k) = count_obs(m,k) + 1
                         endif
                      enddo
                   enddo
                enddo
             enddo

             do m=1,LVT_rc%nensem
                do k=1,obs%selectNlevs

                   allocate(obs_value(count_obs(m,k)))
                   allocate(obs_cluster(count_obs(m,k)))

                   count_obs(m,k) = 1
                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         if(stats%kmeans(m)%obs_value_total(t,k,l).ne.LVT_rc%udef) then 
                            obs_value(count_obs(m,k)) = &
                                 stats%kmeans(m)%obs_value_total(t,k,l)
                            count_obs(m,k) = count_obs(m,k) + 1
                         endif
                      enddo
                   enddo

                   seed = 123456789
                   
                   call cluster_initialize( 1, &
                        count_obs(m,k),&
                        LVT_kmeans_struc%cluster_num,&
                        obs_value, &
                        seed, &
                        cluster_center )
                   
                   call kmeans_01(1,&
                        count_obs(m,k),&
                        LVT_kmeans_struc%cluster_num,&
                        LVT_kmeans_struc%it_max,&
                        LVT_kmeans_struc%it_num,&
                        obs_value,&
                        obs_cluster,&
                        cluster_center, &
                        cluster_population,&
                        cluster_energy)
                   
                   count_obs(m,k) = 1
                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         if(stats%kmeans(m)%obs_value_total(t,k,l).ne.LVT_rc%udef) then 
                            stats%kmeans(m)%obs_mean_cluster(t,k,l) = & 
                                 obs_cluster(count_obs(m,k))
                            count_obs(m,k) = count_obs(m,k) + 1  
                         endif
                      enddo
                   enddo

                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         if(stats%kmeans(m)%count_obs_value_total(t,k,l).ne.0.and.&
                              stats%kmeans(m)%obs_mean_cluster(t,k,l).ne.0) then
                            stats%kmeans(m)%obs_mean_cluster_se(t,k,l) = &
                                 (stats%kmeans(m)%obs_value_total(t,k,l) - &
                                 cluster_center(1,stats%kmeans(m)%obs_mean_cluster(t,k,l)))**2
                         else
                            stats%kmeans(m)%obs_mean_cluster_se(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   enddo
!stdev
                   count_obs(m,k) = 1
                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         if(stats%kmeans(m)%obs_stdev_total(t,k,l).ne.LVT_rc%udef) then 
                            obs_value(count_obs(m,k)) = &
                                 stats%kmeans(m)%obs_stdev_total(t,k,l)
                            count_obs(m,k) = count_obs(m,k) + 1
                         endif
                      enddo
                   enddo
                   
                   seed = 123456789
                   
                   call cluster_initialize( 1, &
                        count_obs(m,k),&
                        LVT_kmeans_struc%cluster_num,&
                        obs_value,&
                        seed, &
                        cluster_center )
                   
                   call kmeans_01(1,&
                        count_obs(m,k),&
                        LVT_kmeans_struc%cluster_num,&
                        LVT_kmeans_struc%it_max,&
                        LVT_kmeans_struc%it_num,&
                        obs_value,&
                        obs_cluster,&
                        cluster_center, &
                        cluster_population,&
                        cluster_energy)

                   count_obs(m,k) = 1
                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         if(stats%kmeans(m)%obs_stdev_total(t,k,l).ne.LVT_rc%udef) then 
                            stats%kmeans(m)%obs_stdev_cluster(t,k,l) = & 
                                 obs_cluster(count_obs(m,k))
                            count_obs(m,k) = count_obs(m,k) + 1  
                         endif
                      enddo
                   enddo

                   do l=1,LVT_rc%strat_nlevels
                      do t=1,LVT_rc%ngrid
                         if(stats%kmeans(m)%count_obs_value_total(t,k,l).ne.0.and.&
                              stats%kmeans(m)%obs_stdev_cluster(t,k,l).ne.0) then
                            stats%kmeans(m)%obs_stdev_cluster_se(t,k,l) = &
                                 (stats%kmeans(m)%obs_stdev_total(t,k,l) - &
                                 cluster_center(1,stats%kmeans(m)%obs_stdev_cluster(t,k,l)))**2
                         else
                            stats%kmeans(m)%obs_stdev_cluster_se(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   enddo

                   deallocate(obs_value)
                   deallocate(obs_cluster)
                enddo
             enddo
          endif

          do t=1,LVT_rc%ngrid
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels
                      if(stats%kmeans(m)%count_model_value_total(t,k,l).eq.0) then
                         stats%kmeans(m)%model_mean_cluster(t,k,l) = LVT_rc%udef
                         stats%kmeans(m)%model_stdev_cluster(t,k,l) = LVT_rc%udef
                      endif
                   enddo
                enddo
             enddo
          enddo

          if(LVT_rc%obssource(2).ne."none") then 
             do t=1,LVT_rc%ngrid
                do m=1,LVT_rc%nensem
                   do k=1,obs%selectNlevs
                      do l=1,LVT_rc%strat_nlevels
                         if(stats%kmeans(m)%count_obs_value_total(t,k,l).eq.0) then
                            stats%kmeans(m)%obs_mean_cluster(t,k,l) = LVT_rc%udef
                            stats%kmeans(m)%obs_stdev_cluster(t,k,l) = LVT_rc%udef
                         endif
                      enddo
                   enddo
                enddo
             enddo
          endif
       endif
    endif
  end subroutine computeSingleKmeans

!BOP
! 
! !ROUTINE: LVT_writeMetric_Kmeans
! \label(LVT_writeMetric_Kmeans)
!
! !INTERFACE:
  subroutine LVT_writeMetric_Kmeans(pass,final,vlevels,stats,model,obs)
! 
! !USES:   
    use LVT_statsMod, only : LVT_writeSummaryStats
    use LVT_pluginIndices
    implicit none
!
! !INPUT PARAMETERS: 
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                 :: pass
    integer                 :: final
    integer                 :: vlevels
    type(LVT_metaDataEntry) :: model
    type(LVT_statsEntry)    :: stats
    type(LVT_metaDataEntry) :: obs

    integer                 :: k,l,m,tind

    real,    allocatable    :: model_cluster(:,:)
    real,    allocatable    :: stdev_model_cluster(:,:)
    real,    allocatable    :: model_mean_cluster_se(:,:)
    real,    allocatable    :: model_stdev_cluster_se(:,:)
    real,    allocatable    :: model_value_ci(:)

    real,    allocatable    :: obs_cluster(:,:)
    real,    allocatable    :: stdev_obs_cluster(:,:)
    real,    allocatable    :: obs_mean_cluster_se(:,:)
    real,    allocatable    :: obs_stdev_cluster_se(:,:)
    real,    allocatable    :: obs_value_ci(:)

    integer, allocatable    :: count_cluster(:,:)

    if(pass.eq.LVT_metrics%kmeans%npass) then
       if(pass.eq.LVT_metrics%kmeans%npass) then
          if(stats%selectOpt.eq.1) then
             
             allocate(model_cluster(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(stdev_model_cluster(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(model_value_ci(LVT_rc%nensem))
             allocate(model_mean_cluster_se(LVT_rc%ngrid,LVT_rc%nensem))
             allocate(model_stdev_cluster_se(LVT_rc%ngrid,LVT_rc%nensem))

             allocate(count_cluster(LVT_rc%ngrid,LVT_rc%nensem))
             
             count_cluster = 1

             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                allocate(obs_cluster(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(stdev_obs_cluster(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(obs_mean_cluster_se(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(obs_stdev_cluster_se(LVT_rc%ngrid,LVT_rc%nensem))
                allocate(obs_value_ci(LVT_rc%nensem))
             endif
             
             do k=1,vlevels
                do l=1,LVT_rc%strat_nlevels
                   do m=1,LVT_rc%nensem
                      model_cluster(:,m) = &
                           real(stats%kmeans(m)%model_mean_cluster(:,k,l))
!                      stdev_model_cluster(:,m) = & 
!                           stats%kmeans(m)%stdev_model_value_total(:,k,l)
                      stdev_model_cluster(:,m) = &
                           real(stats%kmeans(m)%model_stdev_cluster(:,k,l))
                      model_mean_cluster_se(:,m) = &
                           stats%kmeans(m)%model_mean_cluster_se(:,k,l)
                      model_stdev_cluster_se(:,m) = &
                           stats%kmeans(m)%model_stdev_cluster_se(:,k,l)
                      model_value_ci(m) = LVT_rc%udef
                      
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         obs_cluster(:,m) = &
                              real(stats%kmeans(m)%obs_mean_cluster(:,k,l))
                         stdev_obs_cluster(:,m) = &
                              real(stats%kmeans(m)%obs_stdev_cluster(:,k,l))
                         obs_mean_cluster_se(:,m) = &
                              stats%kmeans(m)%obs_mean_cluster_se(:,k,l)
                         obs_stdev_cluster_se(:,m) = &
                              stats%kmeans(m)%obs_stdev_cluster_se(:,k,l)

                         obs_value_ci(m) = LVT_rc%udef
                         
                      endif
                   enddo
                   if(LVT_metrics%kmeans%selectOpt.eq.1) then 

                      call LVT_writevar_gridded(LVT_metrics%kmeans%ftn_total, &
                           model_cluster(:,:),&
                           stats%vid_total(LVT_Kmeansid,1),k)
                      !counts are not output
                      call LVT_writevar_gridded(LVT_metrics%kmeans%ftn_total, &
                           real(stdev_model_cluster(:,:)),&
                           stats%vid_total(LVT_Kmeansid,3),k)

                      call LVT_writevar_gridded(LVT_metrics%kmeans%ftn_total, &
                           model_mean_cluster_se(:,:),&
                           stats%vid_total(LVT_Kmeansid,5),k)

                      call LVT_writevar_gridded(LVT_metrics%kmeans%ftn_total, &
                           model_stdev_cluster_se(:,:),&
                           stats%vid_total(LVT_Kmeansid,7),k)
                      
                      if(LVT_rc%obssource(2).ne."none".and.&
                           obs%selectNlevs.ge.1) then 
                         call LVT_writevar_gridded(&
                              LVT_metrics%kmeans%ftn_total, &
                              obs_cluster(:,:),&
                              stats%vid_total(LVT_Kmeansid,2),k)
                         call LVT_writevar_gridded(&
                              LVT_metrics%kmeans%ftn_total, &
                              real(stdev_obs_cluster(:,:)),&
                              stats%vid_total(LVT_Kmeansid,4),k)

                         call LVT_writevar_gridded(LVT_metrics%kmeans%ftn_total, &
                              obs_mean_cluster_se(:,:),&
                              stats%vid_total(LVT_Kmeansid,6),k)
                         
                         call LVT_writevar_gridded(LVT_metrics%kmeans%ftn_total, &
                              obs_stdev_cluster_se(:,:),&
                              stats%vid_total(LVT_Kmeansid,6),k)
                         
                      endif
                      
                      call LVT_writeSummaryStats(&
                           LVT_metrics%kmeans%ftn_summ,&
                           l,&
                           LVT_metrics%kmeans%short_name,&
                           LVT_rc%ngrid,&
                           model_cluster(:,:), &
                           count_cluster(:,:),&
                           stats%standard_name,&
                           model_value_ci(:))
                      if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                         call LVT_writeSummaryStats(&
                              LVT_metrics%kmeans%ftn_summ,&
                              l,&
                              LVT_metrics%kmeans%short_name,&
                              LVT_rc%ngrid,&
                              obs_cluster(:,:), &
                              count_cluster(:,:),&
                              "DS2_"//trim(stats%standard_name),&
                              obs_value_ci(:))
                      endif
                   endif
                enddo
             enddo
             
             deallocate(model_cluster)
             deallocate(count_cluster)
             deallocate(model_value_ci)
             deallocate(model_mean_cluster_se)
             deallocate(model_stdev_cluster_se)

             
             if(LVT_rc%obssource(2).ne."none".and.obs%selectNlevs.ge.1) then 
                deallocate(obs_cluster)
                deallocate(obs_value_ci)
                deallocate(obs_mean_cluster_se)
                deallocate(obs_stdev_cluster_se)
             endif
             
          endif
       endif
    endif
  end subroutine LVT_writeMetric_Kmeans

!BOP
! 
! !ROUTINE: LVT_resetMetric_Kmeans
! \label(LVT_resetMetric_Kmeans)
!
! !INTERFACE:
  subroutine LVT_resetMetric_Kmeans(alarm)
! 
! !USES:   
!
! !INPUT PARAMETERS: 
    logical                :: alarm
! 
! !OUTPUT PARAMETERS:
!
! !DESCRIPTION: 
! 
! !FILES USED:
!
! !REVISION HISTORY: 
! 
!EOP

    integer                :: i,k,l,m
    type(LVT_metadataEntry), pointer :: model
    type(LVT_metadataEntry), pointer :: obs
    type(LVT_statsEntry)   , pointer :: stats


  end subroutine LVT_resetMetric_Kmeans

!BOP
! 
! !ROUTINE: LVT_writerestart_Kmeans
! 
! !INTERFACE:
  subroutine LVT_writerestart_Kmeans(ftn,pass)
! !USES: 

! 
! !ARGUMENTS: 
    integer                 :: ftn
    integer                 :: pass

! !DESCRIPTION: 
!  This routine writes the restart file for Kmeans metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP
    integer              :: k,l,m
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

#if 0 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%kmeans%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_writevar_restart(ftn,&
                           stats%kmeans(m)%model_value_total(:,k,l))
                      call LVT_writevar_restart(ftn,&
                           stats%kmeans(m)%count_model_value_total(:,k,l))
                   enddo
                enddo
                
                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_writevar_restart(ftn,&
                                 stats%kmeans(m)%obs_value_total(:,k,l))
                            call LVT_writevar_restart(ftn,&
                                 stats%kmeans(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo
                   
                   endif
                endif
             enddo
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
#endif
  end subroutine LVT_writerestart_Kmeans


!BOP
! 
! !ROUTINE: LVT_readrestart_Kmeans
! 
! !INTERFACE:
  subroutine LVT_readrestart_Kmeans(ftn)
! !USES: 
! 
! !ARGUMENTS: 
    integer                 :: ftn

! !DESCRIPTION: 
!  This routine reads the restart file for Kmeans metric computations
! 
!EOP
    
!
! !DESCRIPTION: 
! 
!EOP

    integer              :: k,l,m,index
    type(LVT_metaDataEntry), pointer :: model
    type(LVT_metaDataEntry), pointer :: obs
    type(LVT_statsEntry),    pointer :: stats

#if 0 
    call LVT_getDataStream1Ptr(model)
    call LVT_getDataStream2Ptr(obs)
    call LVT_getstatsEntryPtr(stats)
    
    do while(associated(model))
       if(LVT_metrics%kmeans%selectOpt.eq.1) then 
          if(stats%selectOpt.eq.1.and.&
               model%selectNlevs.ge.1) then 
             do m=1,LVT_rc%nensem
                do k=1,model%selectNlevs
                   do l=1,LVT_rc%strat_nlevels         
                      call LVT_readvar_restart(ftn,&
                           stats%kmeans(m)%model_value_total(:,k,l))
                      call LVT_readvar_restart(ftn,&
                           stats%kmeans(m)%count_model_value_total(:,k,l))
                   enddo
                enddo

                if(LVT_rc%obssource(2).ne."none") then 
                   if(obs%selectNlevs.ge.1) then 
                      do k=1,obs%selectNlevs
                         do l=1,LVT_rc%strat_nlevels
                            call LVT_readvar_restart(ftn,&
                                 stats%kmeans(m)%obs_value_total(:,k,l))
                            call LVT_readvar_restart(ftn,&
                                 stats%kmeans(m)%count_obs_value_total(:,k,l))
                         enddo
                      enddo

                   endif
                endif
             end do
          endif
       endif
       model => model%next
       obs   => obs%next
       stats => stats%next
    end do
#endif
  end subroutine LVT_readrestart_Kmeans

  subroutine cluster_initialize( dim_num, point_num, cluster_num, point, &
       seed, cluster_center )
    
!*****************************************************************************
!
!! CLUSTER_INITIALIZE initializes the cluster centers to random values.
!
!  Discussion:
!
!    In this case, each cluster center is a random convex combination 
!    of the data points.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    06 June 2006
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the coordinates 
!    of the points.
!
!    Input/output, integer ( kind = 4 ) SEED, a seed for the random 
!    number generator.
!
!    Output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the coordinates of the cluster centers.
!
    implicit none
    
    integer ( kind = 4 ) cluster_num
    integer ( kind = 4 ) dim_num
    integer ( kind = 4 ) point_num
    
    real ( kind = 4 ) cluster_center(dim_num,cluster_num)
    real ( kind = 4 ) column_sum
    real ( kind = 4 ) factor(point_num,cluster_num)
    integer ( kind = 4 ) j
    real ( kind = 4 ) point(dim_num,point_num)
    integer ( kind = 4 ) seed
    !
    !  Get a PxC block of random factors.
    !
    call r8mat_uniform_01 ( point_num, cluster_num, seed, factor )
    !
    !  Make each column of factors have unit sum.
    !
    do j = 1, cluster_num
       column_sum = sum ( factor(1:point_num,j) )
       factor(1:point_num,j) = factor(1:point_num,j) / column_sum
    end do
    !
    !  Set centers = points * factors.
    !
    cluster_center(1:dim_num,1:cluster_num) = &
         matmul ( point(1:dim_num,1:point_num), &
         factor(1:point_num,1:cluster_num) )
    
    return
  end subroutine cluster_initialize

  subroutine r8mat_uniform_01 ( m, n, seed, r )

!*****************************************************************************80
!
!! R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
!
!  Discussion:
!
!    An R8MAT is an array of real ( kind = 8 ) values.
!
!    For now, the input quantity SEED is an integer ( kind = 4 ) variable.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    05 July 2006
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Paul Bratley, Bennett Fox, Linus Schrage,
!    A Guide to Simulation,
!    Springer Verlag, pages 201-202, 1983.
!
!    Bennett Fox,
!    Algorithm 647:
!    Implementation and Relative Efficiency of Quasirandom
!    Sequence Generators,
!    ACM Transactions on Mathematical Software,
!    Volume 12, Number 4, pages 362-376, 1986.
!
!    Peter Lewis, Allen Goodman, James Miller
!    A Pseudo-Random Number Generator for the System/360,
!    IBM Systems Journal,
!    Volume 8, pages 136-143, 1969.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns 
!    in the array.
!
!    Input/output, integer ( kind = 4 ) SEED, the "seed" value, which 
!    should NOT be 0.  On output, SEED has been updated.
!
!    Output, real ( kind = 8 ) R(M,N), the array of pseudorandom values.
!
    implicit none
    
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n
    
    integer ( kind = 4 ) i
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) seed
    real ( kind = 4 ) r(m,n)
    
    if ( seed == 0 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
       write ( *, '(a)' ) '  Input value of SEED = 0.'
       stop 1
    end if
    
    do j = 1, n
       
       do i = 1, m
          
          k = seed / 127773
          
          seed = 16807 * ( seed - k * 127773 ) - k * 2836
          
          if ( seed < 0 ) then
             seed = seed + 2147483647
          end if
          
          r(i,j) = real ( seed, kind = 8 ) * 4.656612875D-10
          
       end do
    end do
    
    return
  end subroutine r8mat_uniform_01


  subroutine kmeans_01 ( dim_num, point_num, cluster_num, it_max, it_num, point, &
       cluster, cluster_center, cluster_population, cluster_energy )
    
!*****************************************************************************80
!
!! KMEANS_01 applies the K-Means algorithm.
!
!  Discussion:
!
!    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
!    observations are to be allocated to CLUSTER_NUM clusters in such 
!    a way that the within-cluster sum of squares is minimized.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    18 November 2004
!
!  Author:
!
!    Original FORTRAN77 version by David Sparks.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    David Sparks,
!    Algorithm AS 58: 
!    Euclidean Cluster Analysis,
!    Applied Statistics,
!    Volume 22, Number 1, 1973, pages 126-130.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) DIM_NUM, the number of spatial dimensions.
!
!    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) CLUSTER_NUM, the number of clusters.
!
!    Input, integer ( kind = 4 ) IT_MAX, the maximum number of iterations.
!
!    Output, integer ( kind = 4 ) IT_NUM, the number of iterations taken.
!
!    Input, real ( kind = 8 ) POINT(DIM_NUM,POINT_NUM), the points.
!
!    Output, integer ( kind = 4 ) CLUSTER(POINT_NUM), indicates which cluster
!    each point belongs to.
!
!    Input/output, real ( kind = 8 ) CLUSTER_CENTER(DIM_NUM,CLUSTER_NUM),
!    the cluster centers.
!
!    Output, integer ( kind = 4 ) CLUSTER_POPULATION(CLUSTER_NUM), the number 
!    of points in each cluster.
!
!    Output, real ( kind = 8 ) CLUSTER_ENERGY(CLUSTER_NUM), the 
!    cluster energies.
!
    implicit none

    integer ( kind = 4 ) cluster_num
    integer ( kind = 4 ) dim_num
    integer ( kind = 4 ) point_num
    
    integer ( kind = 4 ) cluster(point_num)
    real ( kind = 4 ) cluster_center(dim_num,cluster_num)
    real ( kind = 4 ) cluster_energy(cluster_num)
    integer ( kind = 4 ) cluster_population(cluster_num)
    real ( kind = 4 ) dc
    real ( kind = 4 ) de
    real ( kind = 4 ) f(point_num)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) il
    integer ( kind = 4 ) ir
    integer ( kind = 4 ) it_max
    integer ( kind = 4 ) it_num
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) list(1)
    real ( kind = 4 ) point(dim_num,point_num)
    integer ( kind = 4 ) swap

    it_num = 0
    !
    !  Idiot checks.
    !
    if ( cluster_num < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       write ( *, '(a)' ) '  CLUSTER_NUM < 1.'
       stop 1
    end if

    if ( dim_num < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       write ( *, '(a)' ) '  DIM_NUM < 1.'
       stop 1
    end if

    if ( point_num < 1 ) then
       write ( *, '(a)' ) ' '
       write ( *, '(a)' ) 'KMEANS_01 - Fatal error!'
       write ( *, '(a)' ) '  POINT_NUM < 1.'
       stop 1
    end if
    !
    !  For each observation, calculate the distance from each cluster
    !  center, and assign to the nearest.
    !
    do i = 1, point_num

       do j = 1, cluster_num
          cluster_energy(j) = sum ( &
               ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
       end do

       list = minloc ( cluster_energy(1:cluster_num) )
       cluster(i) = list(1)

    end do
    !
    !  Determine the cluster population counts.
    !
    cluster_population(1:cluster_num) = 0

    do i = 1, point_num
       j = cluster(i)
       cluster_population(j) = cluster_population(j) + 1
    end do
    !
    !  Calculate the mean and sum of squares for each cluster.
    !
    cluster_center(1:dim_num,1:cluster_num) = 0.0D+00

    do i = 1, point_num
       j = cluster(i)
       cluster_center(1:dim_num,j) = cluster_center(1:dim_num,j) &
            + point(1:dim_num,i)
    end do

    do i = 1, cluster_num
       if ( 0 < cluster_population(i) ) then
          cluster_center(1:dim_num,i) = cluster_center(1:dim_num,i) &
               / real ( cluster_population(i), kind = 4 )
       end if
    end do
    !
    !  Set the point energies.
    !
    f(1:point_num) = 0.0D+00

    do i = 1, point_num
       j = cluster(i)
       f(i) = sum ( ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )
    end do
    !
    !  Set the cluster energies.
    !
    cluster_energy(1:cluster_num) = 0.0D+00

    do i = 1, point_num
       j = cluster(i)
       cluster_energy(j) = cluster_energy(j) + f(i)
    end do
    !
    !  Adjust the point energies by a weight factor.
    !
    do i = 1, point_num
       j = cluster(i)
       if ( 1 < cluster_population(j) ) then
          f(i) = f(i) * real ( cluster_population(j), kind = 4 ) &
               / real ( cluster_population(j) - 1, kind = 4 )
       end if
    end do
    !
    !  Examine each observation in turn to see if it should be
    !  reassigned to a different cluster.
    !
    it_num = 0

    do while ( it_num < it_max )

       it_num = it_num + 1

       swap = 0

       do i = 1, point_num

          il = cluster(i)
          ir = il

          if ( cluster_population(il) <= 1 ) then
             cycle
          end if

          dc = f(i)

          do j = 1, cluster_num

             if ( j /= il ) then

                de = sum ( &
                     ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 ) &
                     * real ( cluster_population(j), kind = 4 ) &
                     / real ( cluster_population(j) + 1, kind = 4 )

                if ( de < dc ) then
                   dc = de
                   ir = j
                end if

             end if

          end do
          !
          !  If the lowest value was obtained by staying in the current cluster,
          !  then cycle.
          !
          if ( ir == il ) then
             cycle
          end if
          !
          !  Reassign the point from cluster IL to cluster IR.
          !
          cluster_center(1:dim_num,il) = &
               ( cluster_center(1:dim_num,il) &
               * real ( cluster_population(il), kind = 4 ) &
               - point(1:dim_num,i) ) / real ( cluster_population(il) - 1, kind = 4 )

          cluster_center(1:dim_num,ir) = &
               ( cluster_center(1:dim_num,ir) &
               * real ( cluster_population(ir), kind = 4 ) &
               + point(1:dim_num,i) ) / real ( cluster_population(ir) + 1, kind = 4 )

          cluster_energy(il) = cluster_energy(il) - f(i)
          cluster_energy(ir) = cluster_energy(ir) + dc
          cluster_population(ir) = cluster_population(ir) + 1
          cluster_population(il) = cluster_population(il) - 1

          cluster(i) = ir
          !
          !  Adjust the value of F for points in clusters IL and IR.
          !
          do j = 1, point_num

             k = cluster(j)

             if ( k == il .or. k == ir ) then

                f(j) = sum ( &
                     ( point(1:dim_num,j) - cluster_center(1:dim_num,k) )**2 )

                if ( 1 < cluster_population(k) ) then
                   f(j) = f(j) * real ( cluster_population(k), kind = 4 ) &
                        / ( real ( cluster_population(k) - 1, kind = 4 ) )
                end if

             end if

          end do

          swap = swap + 1

       end do
       !
       !  Exit if no reassignments were made during this iteration.
       !
       if ( swap == 0 ) then
          exit
       end if

    end do
    !
    !  Compute the cluster energies.
    !
    cluster_energy(1:cluster_num) = 0.0D+00

    do i = 1, point_num

       j = cluster(i)

       cluster_energy(j) = cluster_energy(j) + sum ( &
            ( point(1:dim_num,i) - cluster_center(1:dim_num,j) )**2 )

    end do

    return
  end subroutine kmeans_01

end module LVT_KmeansMod
