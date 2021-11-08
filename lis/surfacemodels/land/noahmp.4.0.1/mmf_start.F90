subroutine mmf_start(n)
    use LIS_coreMod
    use NoahMP401_lsmMod
    use module_sf_noahmpdrv_401
    implicit none 
    integer :: n, row, col, t, ridx, cidx 
    real :: wtddt 
   ! SW, MMF 
    integer, allocatable,dimension(:,:) :: isltyp, ivgtyp
    !real, allocatable :: fdepth(:,:)
    real, allocatable,dimension(:,:) ::  fdepth, topo , area, rechclim, rivercond, &
                                        wtd, riverbed, eqwtd, pexp, smcwtdxy, &
                                        deeprechxy, rechxy, qslatxy, qrfsxy, qspringsxy  
    real, allocatable,dimension(:,:,:) :: smois, sh2o, smoiseq
    wtddt = int(LIS_rc%ts/60) ! in minutes? 

    allocate(isltyp(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(ivgtyp(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(fdepth(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(topo(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(area(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(rechclim(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(rivercond(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(wtd(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(riverbed(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(eqwtd(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(pexp(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(smcwtdxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(deeprechxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(rechxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(qslatxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(qrfsxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(qspringsxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(smois(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(sh2o(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(smoiseq(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    
    do row=NOAHMP401_struc(n)%row_min, NOAHMP401_struc(n)%row_max
        do col=NOAHMP401_struc(n)%col_min, NOAHMP401_struc(n)%col_max
            t = NOAHMP401_struc(n)%rct_idx(row,col)  ! rct_idx is row by col, i.e. row major, Shugong Wang
            if(t .ne. LIS_rc%udef) then
                isltyp(col,row) = NOAHMP401_struc(n)%noahmp401(t)%soiltype
                ivgtyp(col,row) = NOAHMP401_struc(n)%noahmp401(t)%vegetype
            else
                isltyp(col,row) = NOAHMP401_struc(n)%soil2d(row,col) ! soil2d is row major, row by col, Shugong Wang
                ivgtyp(col,row) = NOAHMP401_struc(n)%vege2d(row,col) ! vege2d is row major, row by col, Shugong Wang
            endif
            smois(row,:,col) = NOAHMP401_struc(n)%init_smc(:)
            sh2o(row,:,col) = NOAHMP401_struc(n)%init_smc(:)
            smoiseq(row,:,col) = 0.0
        enddo
    enddo
    !!! 2-D, MMF, SW
    do row=NOAHMP401_struc(n)%row_min, NOAHMP401_struc(n)%row_max
        do col=NOAHMP401_struc(n)%col_min, NOAHMP401_struc(n)%col_max
            ridx = row - NOAHMP401_struc(n)%row_min + 1
            cidx = col - NOAHMP401_struc(n)%col_min + 1
            fdepth(col,row)    = NOAHMP401_struc(n)%fdepth(ridx, cidx)
            topo(col,row)      = NOAHMP401_struc(n)%topo(ridx, cidx)
            area(col,row)      = NOAHMP401_struc(n)%area(ridx, cidx)
            riverbed(col,row)  = NOAHMP401_struc(n)%riverbed(ridx, cidx)
            eqwtd(col,row)     = NOAHMP401_struc(n)%eqwtd(ridx, cidx)
            rivercond(col,row) = NOAHMP401_struc(n)%rivercond(ridx, cidx)
            rechclim(col,row)  = NOAHMP401_struc(n)%rechclim(ridx, cidx)
            pexp(col,row)      = 1.0
        enddo
    enddo

    call  groundwater_init (noahmp401_struc(n)%nsoil,  & !nsoil ,
                            noahmp401_struc(n)%sldpth, & !dzs, 
                            isltyp, ivgtyp, wtddt ,    &
                            fdepth, topo, riverbed, eqwtd, rivercond, pexp , area ,wtd ,  &
                            smois,sh2o, smoiseq, smcwtdxy, deeprechxy, rechxy ,  &
                            qslatxy, qrfsxy, qspringsxy,                  &
                            rechclim  ,                                   &
                            NOAHMP401_struc(n)%row_min, & !ids,
                            NOAHMP401_struc(n)%row_max, & !ide, +1 for test
                            NOAHMP401_struc(n)%col_min, & !jds,
                            NOAHMP401_struc(n)%col_max, & !jde, 
                            1, 1,                       & !kds,kde,
                            NOAHMP401_struc(n)%row_min, & !ims,
                            NOAHMP401_struc(n)%row_max, & !ime, 
                            NOAHMP401_struc(n)%col_min, & !jms,
                            NOAHMP401_struc(n)%col_max, & !jme, 
                            1, 1,                       & !kms,kme,
                            NOAHMP401_struc(n)%row_min, & !ips,
                            NOAHMP401_struc(n)%row_max, & !ipe, 
                            NOAHMP401_struc(n)%col_min, & !jps,
                            NOAHMP401_struc(n)%col_max, & !jpe, 
                            1,1,                        & !kps,kpe,
                            NOAHMP401_struc(n)%row_min, & !its,
                            NOAHMP401_struc(n)%row_max, & !ite, 
                            NOAHMP401_struc(n)%col_min, & !jts,
                            NOAHMP401_struc(n)%col_max, & !jte, 
                            1,1)                          !kts,kte

    do row=NOAHMP401_struc(n)%row_min, NOAHMP401_struc(n)%row_max
        do col=NOAHMP401_struc(n)%col_min, NOAHMP401_struc(n)%col_max
            t = NOAHMP401_struc(n)%rct_idx(row,col)
            NOAHMP401_struc(n)%noahmp401(t)%wtd       = wtd(col,row) 
            NOAHMP401_struc(n)%noahmp401(t)%zwt       = wtd(col,row)  !!!! zwt should be the same as wtd 
            NOAHMP401_struc(n)%noahmp401(t)%rivercond = rivercond(col,row) 
            NOAHMP401_struc(n)%noahmp401(t)%smc(:)    = smois(col,:,row) ! smois
            NOAHMP401_struc(n)%noahmp401(t)%sh2o(:)   = sh2o(col,:,row)
            NOAHMP401_struc(n)%noahmp401(t)%smoiseq(:)= smoiseq(col,:,row) 
            if(isltyp(col,row) .eq. 14) then
                NOAHMP401_struc(n)%noahmp401(t)%smcwtd = 1.0
            else
                NOAHMP401_struc(n)%noahmp401(t)%smcwtd = smcwtdxy(col,row)
            endif
            NOAHMP401_struc(n)%noahmp401(t)%deeprech  = deeprechxy(col,row)
            NOAHMP401_struc(n)%noahmp401(t)%rech      = rechxy(col,row)
            NOAHMP401_struc(n)%noahmp401(t)%qslat     = qslatxy(col,row) 
            NOAHMP401_struc(n)%noahmp401(t)%qrfs      = qrfsxy(col,row) 
            NOAHMP401_struc(n)%noahmp401(t)%qsprings  = qspringsxy(col,row) 
            NOAHMP401_struc(n)%rivercond(row, col)    = rivercond(col,row) !!! make a copy to the 2D paramter data structure 
            NOAHMP401_struc(n)%riverbed(row, col)     = riverbed(col,row)  !!! make a copy to the 2D paramter data structure 
            NOAHMP401_struc(n)%eqwtd(row, col)        = eqwtd(col,row)     !!! make a copy 
        enddo
    enddo 
    deallocate(isltyp)
    deallocate(ivgtyp)
    deallocate(fdepth)
    deallocate(topo)
    deallocate(area)
    deallocate(rechclim)
    deallocate(rivercond)
    deallocate(wtd)
    deallocate(riverbed)
    deallocate(eqwtd)
    deallocate(pexp)
    deallocate(smcwtdxy)
    deallocate(deeprechxy)
    deallocate(rechxy)
    deallocate(qslatxy)
    deallocate(qrfsxy)
    deallocate(qspringsxy)
    deallocate(smois)
    deallocate(sh2o)
    deallocate(smoiseq)
end subroutine mmf_start
