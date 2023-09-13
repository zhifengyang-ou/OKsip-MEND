MODULE MOD_MEND  
! File:   TMEND.F90
! Author: GANGSHENG WANG @ ORNL
! Updated: April 28, 2015
! Created on February 26, 2013, 11:02 AM
! Updated on 7/15/2015
    
    USE MOD_MEND_TYPE
    USE MOD_USRFS
    IMPLICIT NONE

    PRIVATE:: subMEND_INI
    PRIVATE:: subMEND_output
    PRIVATE:: subMEND_output_rate
    PRIVATE:: subMEND_RUN
    PRIVATE:: subMEND
    PRIVATE:: subMEND_PAR
    PRIVATE:: subMEND_CPOOL_UPDATE1
    PRIVATE:: subMEND_NPOOL_UPDATE1
    PRIVATE:: subMEND_CN_UPDATE
    PRIVATE:: sOUT_OPT_h
    
    PUBLIC :: sINP_Read
    PUBLIC :: sOUT_Day2Mon
    PUBLIC :: sOUT_OPT
    PUBLIC :: sOUT_ALL_tscale
    
    PUBLIC:: fMEND_OBJ
    PUBLIC:: fV_MM     !!Michaelis-Menten Kinetics
    PUBLIC:: fAds      !!adsorption-desorption
    
    PUBLIC:: fTArrhenius
    PUBLIC:: fTArh0
    PUBLIC:: fTArh
    PUBLIC:: fT_Linear
    PUBLIC:: fT_CUE
    PUBLIC:: fTAG
    PUBLIC:: fTQ10
    PUBLIC:: fTQ10_MEND
    
    PUBLIC:: fSWP2SWC
    PUBLIC:: fSWC2SWP
    PUBLIC:: fSWP0
    PUBLIC:: fSWP
    PUBLIC:: fSWP_Death
    PUBLIC:: fSWP_OPT
    PUBLIC:: fSWP_CLM
    PUBLIC:: fSWPsat
    PUBLIC:: fSWP_A2D
    PUBLIC:: fSWP_D2A
    PUBLIC:: fWFP_PieceWise0
    PUBLIC:: fWFP_PieceWise
    PUBLIC:: fO2_CONC
    PUBLIC:: fO2_scalar
    PUBLIC:: fNLimit_MB
    PUBLIC:: fpH
    PUBLIC:: fpH0
    PUBLIC:: fPermil
    
    
    CONTAINS
!!---------------------------------------------------------------------------- 
!! [1] MEND MODEL: INITILIZATION, MODEL, RUN CONTROL, COMPUTATION OF OBF
!! [2] MODEL KINETICS
!! [3] SOIL TEMPERATURE SCALAR
!! [4] SOIL WATER SCALAR
!! [5] SOIL pH SCALAR
!! [6] ISOTOPES
!! [7] OUTPUT PROCESSING: VARIABLE EXTRACTION, TIME SCALE CONVERSION
!!----------------------------------------------------------------------------

!=============================================================================
!MEND: Microbial-Enzyme-mediated Nitrification Denitrification & Decomposition
!-----------------------------------------------------------------------------
SUBROUTINE subMEND_INI(sINI)
    !model initialization

    !!ARGUMENTS:
    TYPE(sMEND_INI), intent(inout) :: sINI
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP):: sINP
    INTEGER nObs_time, nObs_var, nObs_col != nSoil*nSubstrate*nObs_var  !# of columns in observation data file

    INTEGER i, j, k, lp
    INTEGER ifobs  
    REAL(8) frISO(const_nISO), frISOadd(const_nISO)
    REAL(8) frac(3), CN_DOM, rCN_LIG2LAB  !!rCN_LIG2LAB = CN_LIGNIN/CN_LABILE 

    CHARACTER(len = 100) sRead!!propName(nProp) !name of 10 properties 
    INTEGER iRead  !!iTemp
    REAL(8) rRead  !!rTemp
    
    
    ifobs = 2
!
!    !!READ observations for model calibration
    lp = 0
    do i=1,sINI%nVARopt
        j = sINI%VARopt_int(i,1)
        sRead = trim(sINI%dirinp)//trim(sINI%VARfile(j))
        open(unit = ifobs, file = sRead, status = 'old')
        read(ifobs,*)
        read(ifobs,*)
        do k=1,sINI%VARopt_int(i,2)
            lp = lp + 1
            read(ifobs,*)sINI%dOBS_opt(lp,1:2)
            sINI%dOBS_opt(lp,3) = i
        end do
        close(ifobs)
    end do
    
!    do i=1,sINI%nOBS_tot
!        print*,sINI%dOBS_opt(i,1:2)
!    end do

!
    frISO(2) = const_Rstd(1)/(1d0 + const_Rstd(1)) !standard ratio of C14/C12 - 1e-12
    frISO(1) = 1d0 - frISO(2) !e.g., C12 

!
!! sINP%CPOOL        
!    sINI % soilDepth    = sINI%dINI(1)  !![cm], soil depth
    sINI % fTexture     = (/sINI%dINI(2),1.d0-sINI%dINI(2)-sINI%dINI(3),sINI%dINI(3)/) !!fraction of SAND, SILT, CLAY
    sINI % CN_MB_avg    = sINI%dINI(4)
    sINI % CN_MB_min    = sINI%dINI(5)
    sINI % CN_MB_max    = sINI%dINI(6)
    sINI % CN_ENZP      = (/sINI%dINI(7), sINI%dINI(8)/)
    sINI % CN_ENZM      = sINI%dINI(9)
    !!fQOM = sINI%dINI(10)
    
    !carbon pools
    sINP % CPOOL % POM = sINI%dINI(12) * (/sINI%LCI0, 1.d0 - sINI%LCI0/) ![mg C/g soil],Particulate Organic Carbon, size of {} determined by nPOM, !(/DBLE(0.958), DBLE(2.290)/)  !reshape((/DBLE(0.958), DBLE(2.290)/),shape(sINP%CPOOL%POM))
    sINP % CPOOL % MOM = sINI%dINI(13)*(1.d0 - sINI%dINI(10)) ![mg C/g soil],Mineral Associate Organic Carbon  !DBLE(27.868) 
    sINP % CPOOL % QOM = sINI%dINI(13) * sINI%dINI(10) ![mg C/g soil],adsorbed phase of DOM  !DBLE(0.100)
    sINP % CPOOL % DOM = sINI%dINI(14) ![mg C/g soil],Dissolved Organic Carbon     !DBLE(0.210)  
    sINP % CPOOL % MB  = sINI%dINI(15) ![mg C/g soil],Microbial Biomass Carbon  !DBLE(0.640)    
    sINP % CPOOL % MBA = sINP % CPOOL % MB * sINI % r0 ![mg C/g soil],Active Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % MBD = sINP % CPOOL % MB * (1.d0 - sINI % r0) ![mg C/g soil],Active Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % ENZP= (/sINI%dINI(16),sINI%dINI(17)/)!!(/1.1d-3, 1.1d-3/) ![mg C/g soil],ENZyme for POM        !reshape((/DBLE(1.1d-3),DBLE(1.1d-3)/),shape(sINP%CPOOL%ENZP))   
    sINP % CPOOL % ENZM= sINI%dINI(18) !!1.4d-3 ![mg C/g soil],Enzyme for MAOC  
    
    sINP % MNPOOL % CO2 = 0.d0 ![mg C/g soil],CO2 in the soil
    CALL subMEND_CPOOL_UPDATE1(sINP % CPOOL)
    
    IF(.not.sINI%Carbon_only) THEN
        !!CN_SOM = sINI%dINI(19)
        !!CN_POM = sINI%dINI(20)
        !!CN_MOM = sINI%dINI(21)
        CN_DOM = sINI%dINI(22)
        !!CN_MB = sINI%dINI(23)
        !!COMPUTE CN_POM1
        frac = (/sINP % CPOOL % POM(1)/(sINI%dINI(13)+sINI%dINI(12)),&
                sINP % CPOOL % POM(2)/(sINI%dINI(13)+sINI%dINI(12)), &
                sINI%dINI(13)/(sINI%dINI(13)+sINI%dINI(12))/)
        rRead = frac(1)/(1.D0/sINI%dINI(19) - frac(3)/sINI%dINI(21) - frac(2)/const_CN_Cellulose)
        
        sINP % NPOOL % POM(1)= sINP % CPOOL % POM(1)/rRead !!sINI%dINI(20)  !!lignin
        sINP % NPOOL % POM(2)= sINP % CPOOL % POM(2)/const_CN_Cellulose  !!cellulose
        sINP % NPOOL % MOM   = sINP % CPOOL % MOM/sINI%dINI(21)
        sINP % NPOOL % QOM   = sINP % CPOOL % QOM/sINI%dINI(21)
        sINP % NPOOL % DOM   = sINP % CPOOL % DOM/sINI%dINI(22)
        sINP % NPOOL % MBA   = sINP % CPOOL % MBA/sINI%dINI(23)
        sINP % NPOOL % MBD   = sINP % CPOOL % MBD/sINI%dINI(23)
        sINP % NPOOL % ENZP  = sINP % CPOOL % ENZP/sINI % CN_ENZP
        sINP % NPOOL % ENZM  = sINP % CPOOL % ENZM/sINI % CN_ENZM

        sINP % MNPOOL % NH4  = sINI%dINI(24)
        sINP % MNPOOL % NO3  = sINI%dINI(25)
        
        CALL subMEND_NPOOL_UPDATE1(sINP % NPOOL, sINP % MNPOOL)
        CALL subMEND_CN_UPDATE(sINP % CPOOL, sINP % NPOOL, sINP % CN)
        !!COMPUTE CN_LITT(3), CN_WOOD(2), CN_ROOT(3)
        rCN_LIG2LAB = sINI%dINI(26)

        sINI % CN_LITT(2) = const_CN_Cellulose
        sINI % CN_WOOD(2) = const_CN_Cellulose
        sINI % CN_ROOT(2) = const_CN_Cellulose

!        sINI % CN_LITT(1) = frac(1)/(1.D0/sINI%CN_LITT_avg - frac(3)/sINI % CN_LITT(3) - frac(2)/sINI % CN_LITT(2))
        frac = sINI%SIN_frac
        sINI % CN_LITT(3) = (frac(1)/rCN_LIG2LAB + frac(3))/(1.D0/sINI%CN_LITT_avg - frac(2)/sINI % CN_LITT(2))
        sINI % CN_LITT(1) = sINI % CN_LITT(3)*rCN_LIG2LAB

        rRead = sum(sINI%SIN_other(1,:))
        if(rRead.lt.1.D-20) then
            frac = 1.D0/3
        else 
            frac = sINI%SIN_other(1,:)/rRead
        end if
        sINI % CN_WOOD(3) = (frac(1)/rCN_LIG2LAB + frac(3))/(1.D0/sINI%CN_WOOD_avg - frac(2)/sINI % CN_WOOD(2))
        sINI % CN_WOOD(1) = sINI % CN_WOOD(3)*rCN_LIG2LAB

        rRead = sum(sINI%SIN_other(2,:))
        if(rRead.lt.1.D-20) then
            frac = 1.D0/3
        else 
            frac = sINI%SIN_other(2,:)/rRead
        end if
!        frac = sINI%SIN_other(2,:)/sum(sINI%SIN_other(2,:))
        sINI % CN_ROOT(3) = (frac(1)/rCN_LIG2LAB + frac(3))/(1.D0/sINI%CN_ROOT_avg - frac(2)/sINI % CN_ROOT(2))
        sINI % CN_ROOT(1) = sINI % CN_ROOT(3)*rCN_LIG2LAB
        
        sINI%SIN_NH4 = sINI%SIN_NH4/sINI%soilDepth  !!convert [mg N/cm2/h] to [mg N/cm3/h]
        sINI%SIN_NO3 = sINI%SIN_NO3/sINI%soilDepth  !!convert [mg N/cm2/h] to [mg N/cm3/h]
        
!        print*,sINI%SIN_NH4(1:5)
!        print*,sINI%SIN_NO3(1:5)
!        print*,sINI % CN_LITT
!        print*,sINI % CN_WOOD
!        print*,sINI % CN_ROOT
!        print*,"debug"
    END IF !!IF(.not.sINI%Carbon_only)
!!     END TYPE sMEND_INP

!!ISOTOPE:
!    do j = 1, const_nPOM
!        sINP % CPOOLIFR % POM(j) = frISO
!        sINP % CPOOLIFR % ENZP(j) = frISO
!
!        sINP % CPOOLI % POM(j) = sINP % CPOOLIFR % POM(j) * sINP % CPOOL % POM(j)
!        sINP % CPOOLI % ENZP(j) = sINP % CPOOLIFR % ENZP(j) * sINP % CPOOL % ENZP(j)
!
!!!        sINP % CADDI % POMadd(j) = sINP % CADD % POMadd(j) * frISOadd
!    end do
!
!!!    sINP % CADDI % DOMadd = sINP % CADD % DOMadd * frISOadd
!
!    sINP % CPOOLIFR % MOM = frISO
!    sINP % CPOOLIFR % ENZM = frISO
!    sINP % CPOOLIFR % QOM = frISO
!    sINP % CPOOLIFR % DOM = frISO
!    sINP % CPOOLIFR % MB = frISO
!    sINP % CPOOLIFR % MBA = frISO
!    sINP % CPOOLIFR % MBD = frISO
!    sINP % CPOOLIFR % CO2 = frISO

    !     print*, "C12 added = ",sINP%CADDI(1)
    !     print*, "C14 added = ",sINP%CADDI(2) 

!    sINP % CPOOLI % MOM = sINP % CPOOLIFR % MOM * sINP % CPOOL % MOM
!    sINP % CPOOLI % ENZM = sINP % CPOOLIFR % ENZM * sINP % CPOOL % ENZM
!    sINP % CPOOLI % QOM = sINP % CPOOLIFR % QOM * sINP % CPOOL % QOM
!    sINP % CPOOLI % DOM = sINP % CPOOLIFR % DOM * sINP % CPOOL % DOM
!    sINP % CPOOLI % MB = sINP % CPOOLIFR % MB * sINP % CPOOL % MB
!    sINP % CPOOLI % MBA = sINP % CPOOLIFR % MBA * sINP % CPOOL % MBA
!    sINP % CPOOLI % MBD = sINP % CPOOLIFR % MBD * sINP % CPOOL % MBD
!    sINP % CPOOLI % CO2 = sINP % CPOOLIFR % CO2 * sINP % CPOOL % CO2


!    do j = 1, const_nISO
!        sINP % CPOOLI(j) % SOM  = sum(sINP % CPOOLI(j) % POM) &
!        &                       + sINP % CPOOLI(j) % MOM + sINP % CPOOLI(j) % QOM
!    end do
!    sINP % CPOOLIFR % SOM = sINP % CPOOLI % SOM/sINP % CPOOL % SOM
!
!    do j = 1, const_nISO - 1
!        do k = 1, const_nPOM
!            sINP % CPOOLI_SIG(j) % POM(k) = &
!            & fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % POM(k), sINP % CPOOLI(j + 1) % POM(k))
!            sINP % CPOOLI_SIG(j) % ENZP(k) = &
!            & fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % ENZP(k), sINP % CPOOLI(j + 1) % ENZP(k))
!        end do
!
!        sINP % CPOOLI_SIG(j) % MOM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MOM, sINP % CPOOLI(j + 1) % MOM)
!        sINP % CPOOLI_SIG(j) % QOM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % QOM, sINP % CPOOLI(j + 1) % QOM)
!        sINP % CPOOLI_SIG(j) % MB = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MB, sINP % CPOOLI(j + 1) % MB)
!        sINP % CPOOLI_SIG(j) % MBA = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBA, sINP % CPOOLI(j + 1) % MBA)
!        sINP % CPOOLI_SIG(j) % MBD = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBD, sINP % CPOOLI(j + 1) % MBD)
!        sINP % CPOOLI_SIG(j) % DOM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % DOM, sINP % CPOOLI(j + 1) % DOM)
!        sINP % CPOOLI_SIG(j) % ENZM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % ENZM, sINP % CPOOLI(j + 1) % ENZM)
!        sINP % CPOOLI_SIG(j) % CO2 = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % CO2, sINP % CPOOLI(j + 1) % CO2)
!        sINP % CPOOLI_SIG(j) % SOM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % SOM, sINP % CPOOLI(j + 1) % SOM)
!    end do !j = 1, const_nISO - 1

    if (sINI % iModel .eq. 0) then !output results for model simulation
!        write(sINI%iFout_VAR_hour, '(i10,13e20.3,13e20.3,13e20.3)') &
!                    0, sINP % CPOOL,sINP % CPOOLI(1),sINP % CPOOLI(2)
        write(sINI%iFout_VAR_hour, '(i10,100e20.3)') &
                    0, sINP % CPOOL, sINP%MNPOOL, sINP%NPOOL, sINP%CN !!,sINP % CPOOLI(1),sINP % CPOOLI(2)
    end if


    sINI % sINP = sINP
!    DEALLOCATE(dINI)
END !!subroutine subMEND_INI

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_INI_Read(dINI,sFile_INI)
    !model initialization

    !!ARGUMENTS:
    REAL(8), ALLOCATABLE :: dINI(:)
    CHARACTER(LEN=*) sFile_INI
    
    !!LOCAL VARIABLES:
    INTEGER nProp !# of soil properties
    INTEGER j !!i, j, k, lp
    INTEGER ifini !!, ifobs  

    CHARACTER(len = 100) sRead!!propName(nProp) !name of 10 properties 
    INTEGER iRead  !!iTemp

    ifini = 1
    open(unit = ifini, file = sFile_INI, status = 'old')
    read(ifini, *) sRead, nProp
    read(ifini, *) sRead !!Line-2
    read(ifini, *) sRead !!Line-3: head: !!ID	Property	 Value
    
    ALLOCATE(dINI(nProp))
   
    do j = 1, nProp
        read(ifini, *) iRead, sRead, dINI(j)
    end do
    close(ifini)
END !!SUBROUTINE subMEND_INI_READ
!-----------------------------------------------------------------------------
SUBROUTINE subMEND_CPOOL_UPDATE1(CPOOL)
    !!ARGUMENTS:
    TYPE(sMEND_POOL), intent(inout) :: CPOOL

    CPOOL % MB  = CPOOL % MBA + CPOOL % MBD
    CPOOL % ENZ = sum(CPOOL%ENZP) + CPOOL%ENZM
    CPOOL % SOM = sum(CPOOL%POM) + CPOOL%MOM + CPOOL%QOM
    CPOOL % TOM = CPOOL%SOM + CPOOL%DOM + CPOOL%MB + CPOOL%ENZ
    CPOOL % TM  = CPOOL % TOM
END !!SUBROUTINE subMEND_POOL_UPDATE1(CPOOL)
!-----------------------------------------------------------------------------

SUBROUTINE subMEND_NPOOL_UPDATE1(NPOOL, MNPOOL)
    !!ARGUMENTS:
    TYPE(sMEND_POOL)    , intent(inout) :: NPOOL
    TYPE(sMEND_POOL_MN) , intent(inout) :: MNPOOL

    NPOOL % MB  = NPOOL % MBA + NPOOL % MBD
    NPOOL % ENZ = sum(NPOOL%ENZP) + NPOOL%ENZM
    NPOOL % SOM = sum(NPOOL%POM) + NPOOL%MOM + NPOOL%QOM
    NPOOL % TOM = NPOOL%SOM + NPOOL%DOM + NPOOL%MB + NPOOL%ENZ
    
    MNPOOL%Nmine= MNPOOL%NH4 + MNPOOL%NO3
    NPOOL % TM  = NPOOL % TOM + MNPOOL%Nmine
END !!SUBROUTINE subMEND_POOL_UPDATE1(NPOOL)
!-----------------------------------------------------------------------------
SUBROUTINE subMEND_CN_UPDATE(CPOOL, NPOOL, CN)
    !!ARGUMENTS:
    TYPE(sMEND_POOL), intent(in)    :: CPOOL
    TYPE(sMEND_POOL), intent(in)    :: NPOOL
    TYPE(sMEND_POOL), intent(OUT)   :: CN

!    if(isnan(NPOOL % POM(1)).or.(NPOOL % POM(1).lt.1.0D-20)) then
!        write(*,*)"Error:POM1:",NPOOL % POM
!    end if
!    if(isnan(NPOOL % POM(2)).or.(NPOOL % POM(2).lt.1.0D-20)) then
!        write(*,*)"Error:POM2:",NPOOL % POM
!    end if
!    
!    if(isnan(NPOOL % QOM).or.(NPOOL % QOM.lt.1.0D-20)) then
!        write(*,*)"Error:QOM:",NPOOL % QOM
!    end if
    
    CN % TM     = CPOOL % TM/NPOOL % TM
    CN % TOM    = CPOOL % TOM/NPOOL % TOM
    CN % SOM    = CPOOL % SOM/NPOOL % SOM
    CN % POM    = CPOOL % POM/NPOOL % POM
    CN % MOM    = CPOOL % MOM/NPOOL % MOM
    CN % QOM    = CPOOL % QOM/NPOOL % QOM
    CN % DOM    = CPOOL % DOM/NPOOL % DOM
    CN % MB     = CPOOL % MB/NPOOL % MB
    CN % MBA    = CPOOL % MBA/NPOOL % MBA
    CN % MBD    = CPOOL % MBD/NPOOL % MBD
    CN % ENZ    = CPOOL % ENZ/NPOOL % ENZ
    CN % ENZP   = CPOOL % ENZP/NPOOL % ENZP
    CN % ENZM   = CPOOL % ENZM/NPOOL % ENZM
    CN % TM_err = CN % TM

END !!SUBROUTINE subMEND_POOL_UPDATE1(CPOOL)
!-----------------------------------------------------------------------------
    
REAL(8) function fMEND_OBJ(xx, sPAR, sINI, sOUT)
!    USE MOD_MEND
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(inout) :: sPAR
    TYPE(sMEND_INI), intent(inOut) :: sINI
    TYPE(sMEND_OUT) sOUT
    REAL(8) xx(sINI%nPar)
    
    !!LCOAL VARIABLES:
    TYPE(sMEND_INP) sINP
    REAL(8) sum1, sum2, stddev
    INTEGER nObj
    INTEGER ibeg, iend,i,j
    
!    INTEGER, PARAMETER :: nSub = 2 !# of substrates for model calibration, i.e., Glucose, Starch
    
!    do i = 1, sPAR%nPar
!        write(*,*)'xx=',i,xx(i)
!    end do
    
    !!Assignment of Parameter Values are moved to SUBROUTINE subMEND_PAR() 
    sINI%LCI0       = xx(1) !initial LCI
    sINI%r0         = xx(2) !fraction of active biomass
!    print*,"ok:",sINI%VARopt_int
    CALL subMEND_INI(sINI) !initialization: initial pool sizes
    CALL subMEND_RUN(xx, sPAR, sINI, sOUT)

    do i = 1, sINI%nVARopt 
        j = sINI%VARopt_int(i,1)
        if(i.eq.1) then
            ibeg = 1
        else
            ibeg = sum(sINI%VARopt_int(1:(i-1),2)) + 1
        end if
        iend = sum(sINI%VARopt_int(1:i,2))
        if(sINI % iModel.eq.4) then !! MCMC
            if(trim(sINI%VARobj(j)).eq."NSEn".or.trim(sINI%VARobj(j)).eq."MARn".or.trim(sINI%VARobj(j)).eq."CORR") then
                sINI%rOBJ(i) = fSSE_Normalize(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), 1, &
                                const_FillValue)
            else
                sINI%rOBJ(i) = fSSE_Normalize(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), 0, &
                                const_FillValue)
            end if
        else  !!((sINI % iModel.ne.4) then !! MCMC
            if(trim(sINI%VARobj(j)).eq."NSEC") then
                sINI%rOBJ(i) = f1NSE(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), &
				sINI%VARobj_tol(j), const_FillValue)
            else if (trim(sINI%VARobj(j)).eq."MARE") then !! MARE 
                sINI%rOBJ(i) = fMARE(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), const_FillValue)
            else if (trim(sINI%VARobj(j)).eq."MARt") then !! MARE with tolerance
                sINI%rOBJ(i) = fMARE_tolerance(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), &
                            sINI%VARobj_tol(j), const_FillValue)
            else if(trim(sINI%VARobj(j)).eq."NSEn") then  !! NSEC with normalized data
                sINI%rOBJ(i) = f1NSE_norm(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), &
                                const_FillValue)
            else if(trim(sINI%VARobj(j)).eq."MARn") then  !! MARE with normalized data
                sINI%rOBJ(i) = fMARE_norm(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), &
                                const_FillValue)
            else if(trim(sINI%VARobj(j)).eq."CORR") then  !! Correlation Coefficient
                sINI%rOBJ(i) = f1CORR(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), 0, &
                                sINI%VARobj_tol(j),const_FillValue)
            else if(trim(sINI%VARobj(j)).eq."CORl") then  !! Correlation Coefficient: log10 transformation of data
                sINI%rOBJ(i) = f1CORR(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), 1, &
                                sINI%VARobj_tol(j),const_FillValue)
            else if (trim(sINI%VARobj(j)).eq."AVGr") then !! dabs(ratio-fAVG_sim/fAVG_obs)
                sINI%rOBJ(i) = f1RAVG_ratio(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), &
                            sINI%VARobj_tol(j), const_FillValue)
            else !! 
                Write(*,'(a11,i2,a2,a6,a50,10a6)')"ERROR! obj-",j," (",trim(sINI%VARobj(j)), &
                        ") is NOT defined in the Objective Functions: ",const_obj_name
                STOP
!            sINI%rOBJ(i) = fMARE(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), const_FillValue)
            end if
        end if
    end do
    
    if(sINI % iModel.eq.4) then !! MCMC
        fMEND_OBJ = fSUM(sINI%nVARopt, sINI%rOBJ, const_FillValue)
!        fMEND_OBJ = sum(sINI%VARopt_int(1:sINI%nVARopt,2))*fSUM(sINI%nVARopt, sINI%rOBJ, const_FillValue)
    else
        fMEND_OBJ = fWAVG(sINI%nVARopt,sINI%rOBJw,sINI%rOBJ,const_FillValue)
    end if
    if(isnan(fMEND_OBJ)) then
        write(*,'(A10,2x,50A10)')"Obj",sINI%Name_PAR
        write(*,'(A10,50E10.2)')"NaN",sPAR
        STOP
    end if 
!    write(*,*)sINI%rOBJ
!    write(*,*)sINI%rOBJw
!    write(*,*)fMEND_OBJ

END !!function fMEND_OBJ

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_output(sDate,ihr, sPAR, sINI, sOUT)
!    USE MOD_MEND
!!Hourly output for all state variables and fluxes
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(in) :: sPAR
    TYPE(sMEND_INI), intent(in) :: sINI
    TYPE(sMEND_OUT), intent(in) :: sOUT
    CHARACTER(LEN=8),intent(in) :: sDate
    INTEGER,         intent(in) :: ihr
    
    !!LOCAL VARIABLES:
    CHARACTER(LEN=2)  str2
    CHARACTER(LEN=10) sDateHr
    
    CALL sInt2Str(ihr,2,str2)
    sDateHr = sDate//str2
!    if (sINI%iModel .eq. 0) then
!    write(sINI%iFout_VAR_hour, '(A10,13e20.3,13e20.3,13e20.3)') &
!        sDateHr, sOUT % CPOOL,sOUT % CPOOLI(1),sOUT % CPOOLI(2)
    write(sINI%iFout_VAR_hour, '(A10,200e20.6)') &
        sDateHr, sOUT % CPOOL, sOUT%MNPOOL, sOUT%NPOOL, sOUT%CN !!,sOUT % CPOOLI(1),sOUT % CPOOLI(2)
    write(sINI%iFout_FLX_hour, '(A10,200e20.6)') sDateHr, sOUT % CFLUX, sOUT%MNFLUX, sOUT%NFLUX
    write(sINI%iFout_PAR_hour, '(A10,100e20.6)') sDateHr, sPAR
!    end if
END !!subroutine subMEND_output

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_output_rate(sDate,ihr, sINI, sPAR, sINP, sOUT)
!    USE MOD_MEND
!!Hourly output for derived parameters
!!"Hour","kPOM1","kPOM2","kMOM","kDOM","kMB","phi","rMBA"
    !!ARGUMENTS:
    TYPE(sMEND_INI), intent(in) :: sINI
    TYPE(sMEND_PAR), intent(in) :: sPAR
    TYPE(sMEND_INP), intent(in) :: sINP
    TYPE(sMEND_OUT), intent(in) :: sOUT   
    CHARACTER(LEN=8),intent(in) :: sDate
    INTEGER,         intent(in) :: ihr
    
    !!LOCAL VARIABLES 
    REAL(8) kPOM1           !!Equivalent 1st-order decomposition rate; k=Vd*ENZP/(POM + Ks)
    REAL(8) kPOM2           !!Equivalent 1st-order decomposition rate; k=Vd*ENZP/(POM + Ks)
    REAL(8) kMOM            !!Equivalent 1st-order decomposition rate; k=Vd*ENZM/(MOM + Ks)
    REAL(8) kDOM            !!Equivalent 1st-order turnover rate of DOM; k=[(Vg+Vm)/Yg]*MB/(DOM + Ks)
    REAL(8) kMBa            !!Equivalent 1st-order turnover rate of Active Mcirobes; k=respiration_rate + mortality_rate+Enzyme_production_rate
    REAL(8) kMBa_in         !!Equivalent 1st-order assimilation rate of active microbes 
    REAL(8) kMBd            !!Equivalent 1st-order turnover rate of dormant microbes
    REAL(8) kMBd_in         !!Equivalent 1st-order microbial dormancy rate 
    REAL(8) kMB             !!Equivalent 1st-order turnover rate of microbes
    REAL(8) kMB_in          !!Equivalent 1st-order assimilation rate of total microbes
    REAL(8) phi             !!DOM saturation level; phi=DOM/(DOM+Ks)
    REAL(8) rMBa            !!Active Fraction of Microbes, r=MBA/MB
    REAL(8) CUE             !!Apparent Microbial Carbon Use Efficiency, CUE=[DOM_to_MBA - CO2_gmo]/DOM_to_MBA
    REAL(8) NUE             !!Apparent Nitrogen Use Efficiency, NUE=YgN
    
    CHARACTER(LEN=2)  str2
    CHARACTER(LEN=10) sDateHr
    
    CALL sInt2Str(ihr,2,str2)
    sDateHr = sDate//str2
    
    phi     = sINP%CPOOL%DOM/(sPAR%KsDOM+sINP%CPOOL%DOM)
    rMBa    = sINP%CPOOL%MBA/sINP%CPOOL%MB
    kDOM    = (sPAR%Vg + sPAR%Vm)/sPAR%Yg*sINP%CPOOL%MBA/(sPAR%KsDOM+sINP%CPOOL%DOM)
    kMBa    = (sPAR%Vg + sPAR%Vm)*(1.D0/sPAR%Yg - 1.D0)*phi &
            + sPAR%rMORT + (sPAR%pENZP+sPAR%pENZM) * sPAR%Vm 
    kMBa_in = (sPAR%Vg + sPAR%Vm)/sPAR%Yg*phi + sOUT%CFLUX%MBD_to_MBA/max(1.0D-20,sINP%CPOOL%MBA)
    kMBd_in = sOUT%CFLUX%MBA_to_MBD/sINP%CPOOL%MBD
    kMBd    = (sOUT%CFLUX%MBD_to_MBA + sOUT%MNFLUX%CO2_maintn_MBD)/sINP%CPOOL%MBD
    kMB     = (sOUT%MNFLUX%CO2_gm + sOUT%CFLUX%MBA_PM)/sINP%CPOOL%MB
    kMB_in  = sOUT%CFLUX%DOM_to_MBA/sINP%CPOOL%MB  
    if (sOUT%CFLUX%DOM_to_MBA > 1.0D-20) then
        CUE = (sOUT%CFLUX%DOM_to_MBA - sOUT%MNFLUX%CO2_gmo)/sOUT%CFLUX%DOM_to_MBA
    else
        CUE = 0.D0
    end if
    NUE     = sPAR%YgN
    
    SELECT CASE (sINI%iKinetics)
        CASE (1,11)  !!First Order
            kPOM1 = sPAR%VdPOM(1)
            kPOM2 = sPAR%VdPOM(2) 
            kMOM = sPAR%VdMOM
        CASE (2)  !!Second Order
            kPOM1 = sPAR%VdPOM(1)*sINP%CPOOL%ENZP(1)
            kPOM2 = sPAR%VdPOM(2)*sINP%CPOOL%ENZP(2)
            kMOM = sPAR%VdMOM*sINP%CPOOL%ENZM
        CASE DEFAULT !!M-M Kinetics, iType = 0
            kPOM1 = sPAR%VdPOM(1)*sINP%CPOOL%ENZP(1)/(sPAR%KsPOM(1)+sINP%CPOOL%POM(1)) 
            kPOM2 = sPAR%VdPOM(2)*sINP%CPOOL%ENZP(2)/(sPAR%KsPOM(2)+sINP%CPOOL%POM(2)) 
            kMOM = sPAR%VdMOM*sINP%CPOOL%ENZM/(sPAR%KsMOM+sINP%CPOOL%MOM)     
    END SELECT
!    if (sINI%iModel .eq. 0) then
    write(sINI%iFout_rate_hour, '(A10,30e20.3)') &
            sDateHr, kPOM1,kPOM2,kMOM,kDOM,kMBa,kMBa_in,kMBd,kMBd_in,kMB,kMB_in,phi,rMBa, CUE, NUE,&
            sOUT%CPOOL%TM_err,sINP%CPOOL%TM,sOUT%CPOOL%TM,sOUT%CFLUX%TOTinp,sOUT%CFLUX%TOTout,&
            sINP%tmp, sINP%SWC, sINP%SWP, sINP%pH
!    end if
END !!subroutine subMEND_output_rate
!-----------------------------------------------------------------------------
SUBROUTINE subMEND_RUN(xx, sPAR, sINI, sOUT)
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(inout) :: sPAR
    TYPE(sMEND_INI), intent(inout) :: sINI
    TYPE(sMEND_OUT), intent(inout) :: sOUT
    REAL(8)        , intent(in)    :: xx(sINI%nPar)
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP) sINP  
    INTEGER*4 i, k 
    INTEGER j, lp, jbeg,jend, tstep, ibeg, iend,iday
    INTEGER iFunc, iObs_time, iHour
    REAL(8) frISOadd(const_nISO)
    INTEGER nHour,nday,nmon
    REAL(8), ALLOCATABLE:: dSIM_d(:,:)  !!daily
    REAL(8), ALLOCATABLE:: dSIM_m(:,:)  !!monthly
    INTEGER nOBS
    
    CHARACTER(len=6) sYM
    CHARACTER(len=8) sDate
    INTEGER iYYYYMM, iyr,imo,ida,ihr, iyr_end_all
    REAL(8) dSIM_h(sINI%nHour_sim,sINI%nVARopt)  !!store hourly output
    REAL(8) dSIM(24,sINI%nVARopt)  !!store hourly output
    
    character(len=200) sFile_inp
    character(len=200) sRead
    INTEGER iRead
    REAL(8) rRead(5)
    
    INTEGER iJDay
    REAL(8) SIN0,STP0,SWC0,SWP0,SpH0
    REAL(8) SIN_mean(const_nHourInYear)
    REAL(8) STP_mean(const_nHourInYear)
    REAL(8) SWC_mean(const_nHourInYear)
    REAL(8) SWP_mean(const_nHourInYear)
    REAL(8) SpH_mean(const_nHourInYear)
    
    
    nday = nDaysbwDates(sINI%sDate_end_all,sINI%sDate_end_sim)
    if(nday.gt.0) then
        !!projection period beyond observation period (sINI%sDate_end_all)
        !!only used for sINI%iScenario == 1, i.e., 1-yr mean hourly data to drive projection run
        !!when nday < 1, do NOT need these mean hourly values 
        sFile_inp = trim(sINI%dirinp)//'ITW_hour_mean.dat'
        open(unit=9001,file=sFile_inp,status='old')
        read(9001,*)sRead
        do j = 1, const_nHourInYear
            read(9001,*)sRead,SIN_mean(j),STP_mean(j),SWC_mean(j),SWP_mean(j),SpH_mean(j)
        end do
        close(9001)
    end if
!!ATTENTION: use ASSOMIATE will result in subroutine invisible in Navigator
!    ASSOMIATE(sINP => sINI%sINP) 
    
    !!replace sDate_beg_all & sDate_end_all with sDate_beg_sim & sDate_end_sim 
    nday = nDaysbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
    nmon = nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
    ALLOCATE(dSIM_d(nday,sINI%nVARopt*2)) !!mean & sd
    ALLOCATE(dSIM_m(nmon,sINI%nVARopt*2)) !!mean & sd
    
!    sINI%LCI0       = xx(1) !initial LCI
!    sINI%r0         = xx(2) !fraction of active biomass
!     
!    CALL subMEND_INI(sINI) !initialization: initial pool sizes
    
    frISOadd(2) = 1d0/(1d0 + sINI%SIN_C12_C14) !C14
    frISOadd(1) = 1d0 - frISOadd(2)            !C12
    
    sINP = sINI%sINP

    !    print*, sINP%CPOOL
    !    print*, sINP%CPOOLI(1)
    !    print*, sINP%CPOOLI(2)
    jbeg = 24*(nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_beg_inp2)-1)+1
    jend = 24*nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_end_inp2)
    
    !!Initial Step:
    sOUT % CPOOL    = sINP % CPOOL
    sOUT % NPOOL    = sINP % NPOOL
    sOUT % MNPOOL   = sINP % MNPOOL
    sOUT % CN       = sINP % CN
!    sOUT % CPOOLI   = sINP % CPOOLI
!    sOUT % CPOOLIFR = sINP % CPOOLIFR
    
    !    nHour = sINI % nHour
    nHour = sINI % nHour_sim
    ibeg = 24*(nDaysbwDates(sINI%sDate_beg_all,sINI%sDate_beg_sim)-1)  !!t=0 for simulation (sDate_beg_sim>=sDate_beg_all)

!!wgs[7/12/2015]: REVERT to use ARRAY instead of READING ITW each step; 
!!wgs[7/12/2015]: REASON: improve computational efficiency since INPUTs are usually available in a couple of years
!    sFile_inp = trim(sINI%dirout)//trim("ITW_hour.dat")
!    open(unit = sINI%iFout_ITW_hour, file = sFile_inp, status = "old")
!    !!skip the 1st two lines and hours before sINI%sDate_beg_sim
!    do lp = 1, ibeg + 2
!        read(sINI%iFout_ITW_hour,*)sRead
!    end do
!!wgs: END
    
    do iday = 1,nday
        CALL sDate_After(iday,sINI%sDate_beg_sim,sDate)
        if(sINI % iModel.eq.0) then
!            write(*,'(a6,I10,a6,I10,a10,a,$)')'ndays=',nday,'iday=',iday,sDate,cBackspace
            write(*,'(a6,a10,f16.2,a,a,$)')'Date=',sDate,iday*100.0/nday,'%',cBackspace
        end if
 
        do lp = 1, 24  !!
            i = (iday - 1)*24 + lp  !!hour
            k = ibeg + i
!            
            sINP % CPOOL    = sOUT % CPOOL
            sINP % NPOOL    = sOUT % NPOOL
            sINP % MNPOOL   = sOUT % MNPOOL
            sINP % CN       = sOUT % CN
            if(isnan(sOUT%CN%MBA)) then
                write(*,*)"iday; hour=", iday, lp, sOUT%CPOOL%MBA, sOUT%NPOOL%MBA, sOUT%CPOOL%DOM, sOUT%NPOOL%DOM
                write(*,"(A10,100E10.2)")"CFLUX=",sOUT%CFLUX
                write(*,"(A10,100E10.2)")"NFLUX=",sOUT%NFLUX
                write(*,"(A10,100E10.2)")"MNFLUX=",sOUT%MNFLUX
                STOP
            end if
!            sINP % CPOOLI = sOUT % CPOOLI
!            sINP % CPOOLIFR = sOUT % CPOOLIFR
    !        end if
            
            if(k.le.sINI%nHour) then
!!wgs[7/12/2015]: REVERT: use ARRAY instead of READING inputs at each step    
!                read(sINI%iFout_ITW_hour,*)iRead,rRead(1:5)
!                sINP%SIN = rRead(1)
!                sINP%tmp = rRead(2)
!                sINP%SWP = rRead(4)
!                sINP%pH  = rRead(5)
!!wgs: END
                
                sINP%SIN = sINI%SIN(k)*sPAR%fINP
                sINP%tmp = sINI%STP(k)
                sINP%SWC = sINI%SWC(k)
                sINP%SWP = sINI%SWP(k)
                sINP%pH  = sINI%SpH(k)
                
                !!External INPUT: BEGIN
                sOUT % CFLUX % POMadd(1) = sINP%SIN * sINI%SIN_frac(1) + sINI%SIN_other(1,1)  !![mg POM/g soil/h], inputs to POM
                sOUT % CFLUX % POMadd(2) = sINP%SIN * sINI%SIN_frac(2) + sINI%SIN_other(1,2)
                sOUT % CFLUX % DOMadd    = sINP%SIN * sINI%SIN_frac(3) + sINI%SIN_other(1,3)  !![mg DOM/g soil/h], inputs to DOM

                if(i.ge.jbeg.and.i.le.jend) then
                    sOUT%CFLUX%POMadd(1) = sOUT%CFLUX%POMadd(1)+sINI%SIN_other(2,1)/DBLE(jend-jbeg+1)
                    sOUT%CFLUX%POMadd(2) = sOUT%CFLUX%POMadd(2)+sINI%SIN_other(2,2)/DBLE(jend-jbeg+1)
                    sOUT%CFLUX%DOMadd    = sOUT%CFLUX%DOMadd   +sINI%SIN_other(2,3)/DBLE(jend-jbeg+1)
                end if
        !        write(*,*)i,sINP % CADD % POMadd(1:2),sINP % CADD % DOMadd
                
                IF(.NOT.sINI%Carbon_only) THEN
                    sOUT % NFLUX % POMadd(1) = sINP%SIN * sINI%SIN_frac(1)/sINI%CN_LITT(1) + sINI%SIN_other(1,1)/sINI%CN_WOOD(1)  !![mg POM/g soil/h], inputs to POM
                    sOUT % NFLUX % POMadd(2) = sINP%SIN * sINI%SIN_frac(2)/sINI%CN_LITT(2) + sINI%SIN_other(1,2)/sINI%CN_WOOD(2)
                    sOUT % NFLUX % DOMadd    = sINP%SIN * sINI%SIN_frac(3)/sINI%CN_LITT(3) + sINI%SIN_other(1,3)/sINI%CN_WOOD(3)  !![mg DOM/g soil/h], inputs to DOM
                    
                    if(i.ge.jbeg.and.i.le.jend) then
                        sOUT%NFLUX%POMadd(1) = sOUT%NFLUX%POMadd(1)+sINI%SIN_other(2,1)/DBLE(jend-jbeg+1)/sINI%CN_ROOT(1)
                        sOUT%NFLUX%POMadd(2) = sOUT%NFLUX%POMadd(2)+sINI%SIN_other(2,2)/DBLE(jend-jbeg+1)/sINI%CN_ROOT(2)
                        sOUT%NFLUX%DOMadd    = sOUT%NFLUX%DOMadd   +sINI%SIN_other(2,3)/DBLE(jend-jbeg+1)/sINI%CN_ROOT(3)
                    end if
                    
                    sOUT%MNFLUX%NH4_dep = sINI%SIN_NH4(k)
                    sOUT%MNFLUX%NO3_dep = sINI%SIN_NO3(k)
                END IF
            else !! POST-DATA PERIOD: simulation period w/o input data, sINI%nHour<i<=sINI%nHour_sim
                !!use available data to fill the missing data after the available period (>sDate_end_all)
                CALL sDate2YMD(sINI%sDate_end_all,iyr_end_all,imo,ida)
                CALL sDate2YMD(sDate,iyr,imo,ida)
                !!---------------------------------------------
                if(sINI%iScenario == 1) then
                    !!METHOD1: use 1-year mean hourly data
                    iJday = iJulianDay(iyr,imo,ida)
                    j = (iJday - 1)*24 + lp
                    SIN0 = SIN_mean(j)*sPAR%fINP
                    STP0 = STP_mean(j)
                    SWC0 = SWC_mean(j)
                    SWP0 = SWP_mean(j)
                    SpH0 = SpH_mean(j)
                !!---------------------------------------------
                else
                    !!METHOD2: use multiple-year available data
                    j = mod(k,sINI%nHour)   !! mod(k,sINI%nHour)+1  
                    if(j==0) j=sINI%nHour
                    SIN0 = sINI%SIN(j)*sPAR%fINP
                    STP0 = sINI%STP(j)
                    SWC0 = sINI%SWC(j)
                    SWP0 = sINI%SWP(j)
                    SpH0 = sINI%SpH(j)
                !!---------------------------------------------
                end if
                
                
                sINP%tmp = STP0
                sINP%SWC = SWC0
                sINP%SWP = SWP0
                sINP%pH  = SpH0
                  
                
!                write(*,*)sINI%delta_STP
!                write(*,*)sINI%delta_SWC
!                write(*,*)SINI % waterRetention_vg
                
                sINP%tmp = STP0 + (iyr - iyr_end_all)*sINI%STP_delta
!                sINP%SWC = sINI%SWC(j) + (iyr - iyr_end_all)*sINI%delta_SWC
!                sINP%SWC = sINI%SWC(j)*(sINI%delta_SWC**(iyr - iyr_end_all))
!                rRead(1) = (0.5d0 + 1d0/(1d0 + dexp(sINI%SWC_logis*(iyr - iyr_end_all))))
                
!                rRead(1) =  (sINI%SWC_logis(1)+(1d0-sINI%SWC_logis(1))/ &
!                            (sINI%SWC_logis(1)+(1d0-sINI%SWC_logis(1))*dexp(sINI%SWC_logis(2)*(iyr - sINI%SIN_logis(3)))))/ &
!                            (sINI%SWC_logis(1)+(1d0-sINI%SWC_logis(1))/ &
!                            (sINI%SWC_logis(1)+(1d0-sINI%SWC_logis(1))*dexp(sINI%SWC_logis(2)*(iyr_end_all - sINI%SIN_logis(3)))))
!               
                rRead(1) = (1d0-(1d0-sINI%SWC_logis(1))*dexp(-sINI%SWC_logis(2)*(iyr_end_all - sINI%SWC_logis(3))))/ &
                           (1d0-(1d0-sINI%SWC_logis(1))*dexp(-sINI%SWC_logis(2)*(iyr - sINI%SWC_logis(3))))
                           
                 !! 8/21/2016, When const Soil Moisture is used, only SWP is set in mend.ini, thus we need to get const SWC first as follows
!                !! 12/6/2016, When const Soil Moisture is used, only SWC is set in mend.ini
!                           sINI%SWC(j) = fSWP2SWC(sINI%SWP(j),"MPa", &
!                                            SINI % waterRetention_vg(1),SINI % waterRetention_vg(2),&
!                                            SINI % waterRetention_vg(3),SINI % waterRetention_vg(4))
                
                sINP%SWC = SWC0*rRead(1) 
!                write(*,*)'SWC_logis = ',rRead(1),iyr,sINI%SWC(j),sINP%SWC
                
                sINP%SWP = fSWC2SWP(sINP%SWC,SINI % waterRetention_vg(1),SINI % waterRetention_vg(2),&
                                             SINI % waterRetention_vg(3),SINI % waterRetention_vg(4),const_SWPmin)
                
                !!External INPUT: BEGIN; assume root input = leaf litter input
                
!               sOUT%CFLUX%POMadd(1) = sINI%SIN(j) * sINI%SIN_frac(1) * sINI%SIN_Multiplier !!+ sINI%SIN_other(1,1)  !![mg POM/g soil/h], inputs to POM
!               sOUT%CFLUX%POMadd(2) = sINI%SIN(j) * sINI%SIN_frac(2) * sINI%SIN_Multiplier !!+ sINI%SIN_other(1,2)
!               sOUT%CFLUX%DOMadd    = sINI%SIN(j) * sINI%SIN_frac(3) * sINI%SIN_Multiplier
                !! update sINI%SIN_Multiplier according to the logistic equation
                rRead(5) = (1+dexp(sINI%SIN_logis(1)-sINI%SIN_logis(2)*(iyr_end_all - sINI%SIN_logis(3))))/ &
                           (1+dexp(sINI%SIN_logis(1)-sINI%SIN_logis(2)*(iyr - sINI%SIN_logis(3))))
!                write(*,*)'SIN_logis = ',rRead(1), iyr, sINI%SIN_logis
                
                !!wgs[9/27/2016]: currently only assumes linear increase in DOM fraction, needs to be modified 
                rRead(3) = min(sINI%SIN_frac(3)*2d0, sINI%SIN_frac(3) + sINI%SIN_logis(4)*(iyr - iyr_end_all))  !! new DOM fraction
                rRead(1) = sINI%SIN_frac(1)/(sINI%SIN_frac(1) + sINI%SIN_frac(2)) * (1d0 - rRead(3))             !! new POM1 fraction
                rRead(2) = sINI%SIN_frac(2)/(sINI%SIN_frac(1) + sINI%SIN_frac(2)) * (1d0 - rRead(3))             !! new POM2 fraction
                
                sOUT%CFLUX%POMadd(1) = SIN0 * rRead(1) * sINI%SIN_Multiplier * rRead(5) !!+ sINI%SIN_other(1,1)  !![mg POM/g soil/h], inputs to POM
                sOUT%CFLUX%POMadd(2) = SIN0 * rRead(2) * sINI%SIN_Multiplier * rRead(5) !!+ sINI%SIN_other(1,2)
                sOUT%CFLUX%DOMadd    = SIN0 * rRead(3) * sINI%SIN_Multiplier * rRead(5) !!+ sINI%SIN_other(1,3)  !![mg DOM/g soil/h], inputs to DOM
                
                
                IF(.NOT.sINI%Carbon_only) THEN
                    sOUT%NFLUX%POMadd(1) = SIN0 * sINI%SIN_frac(1) * sINI%SIN_Multiplier/sINI%CN_LITT(1) !!+ sINI%SIN_other(1,1)  !![mg POM/g soil/h], inputs to POM
                    sOUT%NFLUX%POMadd(2) = SIN0 * sINI%SIN_frac(2) * sINI%SIN_Multiplier/sINI%CN_LITT(2) !!+ sINI%SIN_other(1,2)
                    sOUT%NFLUX%DOMadd    = SIN0 * sINI%SIN_frac(3) * sINI%SIN_Multiplier/sINI%CN_LITT(3)
                    
                    sOUT%MNFLUX%NH4_dep = sINI%SIN_NH4(j)   !!TODO
                    sOUT%MNFLUX%NO3_dep = sINI%SIN_NO3(j)   !!TODO
                END IF
            end if

!            do j = 1, const_nPOM
!                sINP % CADDI % POMadd(j) = sINP%CADD%POMadd(j)*frISOadd
!            end do
!            sINP % CADDI % DOMadd = sINP%CADD%DOMadd*frISOadd
            !!External INPUT: END

            sINI%sINP = sINP
            CALL subMEND_PAR(xx, sPAR, sINI)   !!modify parameter values by soil temperature, water potential, pH
            CALL subMEND(xx, sPAR, sINI, sOUT) !!MEND model 
            if (sINI % iModel .eq. 0) then !output results for model simulation
                CALL subMEND_output_rate(sDate,lp, sINI, sPAR, sINP, sOUT)
                CALL subMEND_output(sDate,lp, sPAR, sINI, sOUT)
            end if

            !!Extract output variables to compare with observations 
            CALL sOUT_OPT_h(sINI%nVARopt,24,lp,dSIM,sPAR,sOUT,sINI%VARopt_int)
            dSIM_h((iday-1)*24+lp,:) = dSIM(lp,:)   !!continuous hourly output
        end do !!lp = 1, 24  
        
        !!convert hourly data to daily average data
        do j = 1,sINI%nVARopt
            dSIM_d(iday,j*2-1) = fAVG(24,dSIM(:,j),const_FillValue)
            dSIM_d(iday,j*2)   = fSTDDEV(24,dSIM(:,j),const_FillValue)
        end do
    
    end do !!iday = 1,nday
!    close(sINI%iFout_ITW_hour)
    
!    do i=1,nHour
!        write(*,*)i,dSIM(i,:)
!    end do
    
    !!convert daily to monthly average data
    do j = 1,sINI%nVARopt
        CALL sOUT_Day2Mon(nday,dSIM_d(:,j*2-1),nmon,dSIM_m(:,j*2-1:j*2),sINI%sDate_beg_sim,sINI%sDate_end_sim,2,1)
    end do
    
    if (sINI % iModel .eq. 0) then !output results for model simulation       
        do i = 1,nday
            CALL sDate_After(i,sINI%sDate_beg_sim,sDate)
            write(sINI%iFout_SIM_day,'(A10,20f20.6)')sDate,dSIM_d(i,:)
        end do

        do i = 1,nmon
            CALL sYM_After(i,sINI%sDate_beg_sim(1:6),sYM)
            write(sINI%iFout_SIM_mon,'(A10,20f20.6)')sYM,dSIM_m(i,:)
        end do
        write(*,*)">>>Extracting data matching the observations on the same days or months"
    end if
!    
!    do i = 1,nmon
!        write(*,*)i,dSIM_m(i,:)
!    end do
    
    !!Extract those data matching the observations on the same days or months
    CALL sDate2YMD(sINI%sDate_beg_sim,iyr,imo,ida)
    do j = 1,sINI%nVARopt
        nOBS = sINI%VARopt_int(j,2)
        tstep = sINI%VARopt_int(j,3)
        if(j.eq.1) then
            jbeg = 1
        else
            jbeg = sum(sINI%VARopt_int(1:(j-1),2)) + 1
        end if
        do k =1,nOBS 
            lp = jbeg+k-1
            sINI%dSIM_opt(lp,1) = sINI%dOBS_opt(lp,1)
            if(tstep.eq.0) then  !!hourly
                iRead = int(sINI%dOBS_opt(lp,1)/100)    !!int(YYYYMMDDHH/100) = YYYYMMDD
                write(sDate,'(I8)')iRead                !!convert to CHARACTER(len=8)
                ihr = int(sINI%dOBS_opt(lp,1) - iRead*100)
                i = nDaysbwDates(sINI%sDate_beg_sim, sDate)*24 - (24 - ihr)
                sINI%dSIM_opt(lp,2) = dSIM_h(i,j)
                sINI%dSIM_opt(lp,3) = 0.D0        !!no STANDARD DEVIATION for hourly data
            else if(tstep.eq.1) then  !!daily
                write(sDate,'(I8)')int(sINI%dOBS_opt(lp,1))  !!convert REAL(8) to CHARACTER(len=8)
                i = nDaysbwDates(sINI%sDate_beg_sim, sDate)
!                sINI%dSIM_opt(lp,1) = sINI%dOBS_opt(lp,1)
                if(i.GE.1.and.i.LE.nday) then
                    sINI%dSIM_opt(lp,2) = dSIM_d(i,j*2-1)
                    sINI%dSIM_opt(lp,3) = dSIM_d(i,j*2)
                end if
            else if(tstep.eq.2) then  !!monthly
                iYYYYMM = int(sINI%dOBS_opt(lp,1))
                i = nMonths(iyr*100+imo,iYYYYMM)
!                sINI%dSIM_opt(lp,1) = sINI%dOBS_opt(lp,1)
                if(i.GE.1.and.i.LE.nmon) then
                    sINI%dSIM_opt(lp,2) = dSIM_m(i,j*2-1)
                    sINI%dSIM_opt(lp,3) = dSIM_m(i,j*2)
                end if
            else if(tstep.eq.4) then  !!yearly
                iend = int(sINI%dOBS_opt(lp,1)) !!year
                ibeg = max(1,nMonths(iyr*100+imo,iend*100+1))       !!in case of: imo<>1
                iend = min(nmon,nMonths(iyr*100+imo,iend*100+12))   !!in case of: YYYYMM_end < iend*100 + 12
                sINI%dSIM_opt(lp,2) = fAVG((iend-ibeg+1),dSIM_m(ibeg:iend,j*2-1),const_FillValue)  
                sINI%dSIM_opt(lp,3) = fSTDDEV((iend-ibeg+1),dSIM_m(ibeg:iend,j*2-1),const_FillValue)
            else if(tstep.eq.5) then  !!overall mean value, k=1
!                sINI%dSIM_opt(lp,2) = fAVG(nmon,dSIM_m(:,j*2-1),const_FillValue) 
!                sINI%dSIM_opt(lp,3) = fSTDDEV(nmon,dSIM_m(:,j*2-1),const_FillValue)
!                sINI%dSIM_opt(lp,3) = fSTDDEV(nday,dSIM_d(:,j*2-1),const_FillValue)
                ibeg = max(1,int(nday*0.1)) 
                iend = nday
                sINI%dSIM_opt(lp,2) = fAVG(nday,dSIM_d(ibeg:iend,j*2-1),const_FillValue) 
                sINI%dSIM_opt(lp,3) = fSTDDEV(nday,dSIM_d(ibeg:iend,j*2-1),const_FillValue)
            end if
        end do !!do k =1,nOBS
    end do !!do j = 1,sINI%nVARopt
    
    DEALLOCATE(dSIM_d)
    DEALLOCATE(dSIM_m)
!    END ASSOMIATE
!    do i=1,sINI%nOBS_tot 
!        write(*,*)i,sINI%dOBS_opt(i,:),sINI%dSIM_opt(i,:)
!    end do
END !!subroutine subMEND_RUN
!-----------------------------------------------------------------------------

SUBROUTINE subMEND(xx, sPAR, sINI, sOUT)
    !model components 
!    USE MOD_MEND
!    IMPLICIT NONE
    !!ARGUMENTS:
    
    TYPE(sMEND_PAR), intent(inout)  :: sPAR
    TYPE(sMEND_INI), intent(in)     :: sINI
    TYPE(sMEND_OUT), intent(out)    :: sOUT
    REAL(8)        , intent(in)     :: xx(sINI%nPar)  !!used in subMEND_PAR
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP) sINP
    TYPE(sMM_PAR) smmPAR
    TYPE(sMM_INP) smmINP

    TYPE(sSORP_PAR) sSorpPAR
    TYPE(sSORP_INP) sSorpINP
    TYPE(sSORP_OUT) sSorpOUT

    INTEGER i, j, k, iFunc
    REAL(8) sumPOM, sumENZP, phi  !!OC1, OC2,
    REAL(8) frPOM(const_nPOM), frENZP(const_nPOM), ENZP_loss(const_nPOM) !!POMdec(const_nPOM),
    REAL(8) rtp1, rtp2, rtp3, rtp4
    
!    ASSOCIATE(sINP => sINI%sINP) 
    sINP = sINI%sINP

!    write(*,'(a15,50d10.3)')"sINP%CPOOL=",sINP%CPOOL
    !compute initial total mass of the system 
    sOUT%MNFLUX%Nmine_dep = sOUT%MNFLUX%NH4_dep + sOUT%MNFLUX%NO3_dep
    sOUT%CFLUX%TOTinp = sum(sOUT%CFLUX%POMadd) + sOUT%CFLUX%DOMadd
    sOUT%NFLUX%TOTinp = sum(sOUT%NFLUX%POMadd) + sOUT%NFLUX%DOMadd + sOUT%MNFLUX%Nmine_dep
    
    sOUT % CPOOL    = sINP % CPOOL
    sOUT % NPOOL    = sINP % NPOOL
    sOUT % MNPOOL   = sINP % MNPOOL
    !Flux 2: POM decomposition
    DO i = 1, const_nPOM
        smmINP % substrate = sINP % CPOOL % POM(i)
        smmINP % enzyme = sINP % CPOOL % ENZP(i)
        smmPAR % vm = sPAR % VdPOM(i)
        smmPAR % km = sPAR % KsPOM(i)
        sOUT%CFLUX%POMdec(i) = fV_MM(sINI%iKinetics, smmPAR, smmINP)
        
        sOUT%CFLUX%POMdec_to_MOM(i) = (1.d0 - sPAR%frPOM2DOM) * sOUT%CFLUX%POMdec(i) 
        
        if(sINI%iKinetics.lt.10) then
            sOUT%CFLUX%POMdec_to_DOM(i) = sPAR%frPOM2DOM*sOUT % CFLUX % POMdec(i) 
        end if
    end DO!for i = 1:sINP%nPOM
    
    if(sINI%iKinetics.eq.11) then  !!First-Order w/o Microbes/Enzymes
        sOUT % MNFLUX % CO2_gm = sPAR%frPOM2DOM*sum(sOUT%CFLUX%POMdec)     
    end if
    
    !Flux 3: MOM decomposition
    smmINP % substrate = sINP % CPOOL % MOM
    smmINP % enzyme = sINP % CPOOL % ENZM
    smmPAR % vm = sPAR % VdMOM
    smmPAR % km = sPAR % KsMOM
    
    sOUT%CFLUX%MOMdec = fV_MM(sINI%iKinetics, smmPAR, smmINP) !MOM decomposition
    
    if(sINI%iKinetics.lt.10) then
        sOUT % CFLUX % MOM_to_DOM = sOUT%CFLUX%MOMdec
    else if(sINI%iKinetics.eq.11) then  !!First-Order w/o Microbes/Enzymes
        sOUT % MNFLUX % CO2_gm = sOUT % MNFLUX % CO2_gm + sOUT%CFLUX%MOMdec     
    end if
    
    !!------------------------------------------------------------------------
    !!MEND with MICROBES/ENZYMES: BEG
    IF(sINI%iKinetics.lt.10) then
       !Flux 1: DOM uptake by MB
        smmINP % substrate = sINP % CPOOL % DOM
        smmINP % enzyme = sINP % CPOOL % MBA
        smmPAR % vm = (sPAR % Vg + sPAR % Vm)/sPAR % Yg !GROWTH + MAINTENANCE    
        smmPAR % km = sPAR % KsDOM
        sOUT % CFLUX % DOM_to_MBA = fV_MM(0, smmPAR, smmINP)

        !!FIRST UPDATE, used for sorption-desorption calculation
        sOUT % CPOOL % DOM = sINP % CPOOL % DOM - sOUT % CFLUX % DOM_to_MBA  

        !Flux 4 & 5: Growth and Maintenance Respiration
        !emission of CO2 from soil to air NOT considered  
        sOUT % MNFLUX % CO2_maintn_MBD = sPAR%VmD * sINP%CPOOL%MBD  !!Microbial Maintenance of Dormant Microbes 
        
        !---------------------------------------------
        !! Calculation of CO2_growth & CO2_maintn_MBA
        !! 2 methods may result in different results under DOM-limited condition (DOC_to_MBA = DOC).
        if (sINI%iHR .eq. 0) then
            !! fraction of potential (growth/maintenance) uptake
            smmINP % substrate = sINP % CPOOL % DOM
            smmINP % enzyme = sINP % CPOOL % MBA
            smmPAR % vm = sPAR % Vg * (1.D0/sPAR % Yg - 1.D0)
            smmPAR % km = sPAR % KsDOM
            sOUT % MNFLUX % CO2_growth = fV_MM(0, smmPAR, smmINP)

            smmPAR % vm = sPAR % Vm * (1.D0/sPAR % Yg - 1.D0)
            sOUT % MNFLUX % CO2_maintn_MBA = fV_MM(0, smmPAR, smmINP) 
        else
            !! fraction of actual total uptake
            phi = sPAR % Vm/(sPAR % Vg + sPAR % Vm)
            sOUT % MNFLUX % CO2_growth      = sOUT % CFLUX % DOM_to_MBA * (1.D0 - sPAR % Yg) * (1.D0 - phi)
            sOUT % MNFLUX % CO2_maintn_MBA  = sOUT % CFLUX % DOM_to_MBA * (1.D0 - sPAR % Yg) * phi
        end if

        
        !---------------------------------------------

        sOUT % MNFLUX % CO2_maintn = sOUT % MNFLUX % CO2_maintn_MBA + sOUT % MNFLUX % CO2_maintn_MBD
        sOUT % MNFLUX % CO2_gm = sOUT%MNFLUX%CO2_growth + sOUT%MNFLUX%CO2_maintn 

        !Flux 6 & 7: DOM adsorption-desorption
        sSorpINP % adsorbate = sOUT % CPOOL % DOM !!DOM has been update;  = sINP%CPOOL%DOM - sOUT%CFLUX%DOM_to_MB
        sSorpINP % adsorbent = sINP % CPOOL % QOM
        sSorpPAR % Qmax = sPAR % Qmax
        sSorpPAR % Kads = sPAR % Kads
        sSorpPAR % Kdes = sPAR % Kdes
        iFunc = fAds(sSorpPAR, sSorpINP, sSorpOUT)
        sOUT % CFLUX % QOM_to_DOM = sSorpOUT % des
        sOUT % CFLUX % DOM_to_QOM = sSorpOUT % ads
        sOUT % CFLUX % DOM_to_QOM_net = sSorpOUT % ads_net

        !Flux 8: microbial mortality
    !    sOUT % CFLUX % MB_PM = sPAR % Vm * sINP % CPOOL % MBA
    !    sOUT % CFLUX % MB_mortality = (1.d0 - sPAR % pENZP - sPAR % pENZM) * sOUT % CFLUX % MB_PM
        sOUT % CFLUX % MBA_mortality= sPAR%rMORT * sINP%CPOOL%MBA
        sOUT % CFLUX % MBA_to_DOM   = sPAR % frMB2DOM * sOUT % CFLUX % MBA_mortality
        sOUT % CFLUX % MBA_to_POM   = (1.d0 - sPAR % frMB2DOM) * sOUT % CFLUX % MBA_mortality
        !Flux 9: enzyme synthesis
        sumPOM = sum(sINP % CPOOL % POM)
        frPOM = sINP % CPOOL % POM/sumPOM !compute the fraction of each type of POM, e%g%, lignin, cellulose 
    !    sOUT % CFLUX % MB_to_ENZP = frPOM * sPAR % pENZP * sOUT % CFLUX % MB_PM
    !    sOUT % CFLUX % MB_to_ENZM = sPAR % pENZM * sOUT % CFLUX % MB_PM
        sOUT % CFLUX % MBA_to_ENZP = frPOM * sPAR % pENZP * sPAR % Vm * sINP % CPOOL % MBA
        sOUT % CFLUX % MBA_to_ENZM =         sPAR % pENZM * sPAR % Vm * sINP % CPOOL % MBA

        sOUT % CFLUX % MBA_PM = sOUT % CFLUX % MBA_mortality + sum(sOUT % CFLUX % MBA_to_ENZP) + sOUT % CFLUX % MBA_to_ENZM

        !Flux 10: enzyme turnover
        sOUT % CFLUX % ENZP_to_DOM = sPAR % rENZP * sINP % CPOOL % ENZP
        sOUT % CFLUX % ENZM_to_DOM = sPAR % rENZM * sINP % CPOOL % ENZM

        !Internal Flux: microbial dormancy and reactivation
!        IF(sINI%Carbon_only) THEN
        phi = sINP % CPOOL % DOM/(sINP % CPOOL % DOM + sPAR % KsDOM)
!        ELSE
!            phi = (sINP % CPOOL % DOM/(sINP % CPOOL % DOM + sPAR % KsDOM) &
!                  +sINP%MNPOOL%NH4/(sINP%MNPOOL%NH4 + sPAR%KsNH4_MB) &
!                  +sINP%MNPOOL%NO3/(sINP%MNPOOL%NO3 + sPAR%KsNO3_MB))/3.D0
!        END IF
        sOUT % CFLUX % MBA_to_MBD = (1.d0 - phi) * sPAR % VmA2D * sINP % CPOOL % MBA
        sOUT % CFLUX % MBD_to_MBA = phi * sPAR % VmD2A * sINP % CPOOL % MBD
    
    END IF !!IF(sINI%iKinetics.lt.10)
    sOUT % MNFLUX % CO2_gmo = sOUT % MNFLUX % CO2_gm  !!no C-overflow for Carbon-only
    !!MEND with MICROBES/ENZYMES: END
    !!------------------------------------------------------------------------
    
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !!C MASS BALANCE: BEG   
    sOUT % CPOOL % POM      = sINP % CPOOL % POM + sOUT % CFLUX % POMadd - sOUT % CFLUX % POMdec !Lignin and Cellulose
    sOUT % CPOOL % POM(1)   = sOUT % CPOOL % POM(1) + sOUT % CFLUX % MBA_to_POM !ONLY Lignin!     
    sOUT % CPOOL % MOM      = sINP % CPOOL % MOM + sum(sOUT % CFLUX % POMdec_to_MOM) - sOUT%CFLUX%MOMdec!!sOUT % CFLUX % MOM_to_DOM
    
    IF(sINI%iKinetics.lt.10) THEN
      
        sOUT % CPOOL % ENZP     = sINP % CPOOL % ENZP + sOUT % CFLUX % MBA_to_ENZP - sOUT % CFLUX % ENZP_to_DOM
        sOUT % CPOOL % ENZM     = sINP % CPOOL % ENZM + sOUT % CFLUX % MBA_to_ENZM - sOUT % CFLUX % ENZM_to_DOM
        
        sOUT % CPOOL % QOM      = sINP % CPOOL % QOM + sOUT % CFLUX % DOM_to_QOM_net !QOM mass balance
        sOUT % CPOOL % MB       = sINP % CPOOL % MB + sOUT % CFLUX % DOM_to_MBA &
        &                       - sOUT % MNFLUX % CO2_gm - sOUT % CFLUX % MBA_PM
        
        !! sOUT % CPOOL % MBA = rtp1
        sOUT % CPOOL % MBA      = sINP % CPOOL % MBA + sOUT % CFLUX % DOM_to_MBA &
        &                       - (sOUT % MNFLUX % CO2_growth + sOUT % MNFLUX % CO2_maintn_MBA) &
        &                       - sOUT % CFLUX % MBA_PM &
        &                       - sOUT % CFLUX % MBA_to_MBD + sOUT % CFLUX % MBD_to_MBA
        
        if(sOUT % CPOOL % MBA.lt.0.D0) then  !!11/20/2018
            sOUT % CFLUX % MBA_to_MBD = sOUT % CFLUX % MBA_to_MBD + sOUT % CPOOL % MBA
            sOUT % CPOOL % MBA = 0.D0
        end if
        
        sOUT % CPOOL % MBD      = sINP % CPOOL % MBD - sOUT % MNFLUX % CO2_maintn_MBD &
        &                       + sOUT % CFLUX % MBA_to_MBD - sOUT % CFLUX % MBD_to_MBA
        
        sOUT % CPOOL % DOM      = sINP % CPOOL % DOM + sOUT % CFLUX % DOMadd &
        &                       + sum(sOUT % CFLUX % POMdec_to_DOM) &
        &                       + sOUT % CFLUX % MBA_to_DOM &
        &                       + sOUT % CFLUX % MOM_to_DOM &
        &                       + sum(sOUT % CFLUX % ENZP_to_DOM) &
        &                       + sOUT % CFLUX % ENZM_to_DOM &
        &                       - sOUT % CFLUX % DOM_to_MBA &
        &                       - sOUT % CFLUX % DOM_to_QOM_net
        

    ELSE IF(sINI%iKinetics.eq.11) THEN
        sOUT % CPOOL % MOM  = sOUT % CPOOL % MOM + sOUT % CFLUX % DOMadd
    END IF!!IF(sINI%iKinetics.lt.10)
    
    sOUT % CPOOL % ENZ      = sum(sOUT % CPOOL % ENZP) + sOUT % CPOOL % ENZM
    sOUT % CPOOL % SOM      = sum(sOUT % CPOOL % POM) + sOUT % CPOOL % MOM + sOUT % CPOOL % QOM
    sOUT % CPOOL % TOM      = sOUT % CPOOL % SOM + sOUT % CPOOL % DOM + sOUT % CPOOL % MB + sOUT % CPOOL % ENZ
    sOUT % CPOOL % TM       = sOUT % CPOOL % TOM
    
    sOUT % MNPOOL % CO2     = sINP % MNPOOL % CO2 + sOUT % MNFLUX % CO2_gmo
    
    !!C MASS BALANCE CHECK
    sOUT%CFLUX%TOTout = sOUT%MNFLUX%CO2_gmo
    sOUT%CPOOL%TM_err = (sOUT%CPOOL%TM - sINP%CPOOL%TM) - (sOUT%CFLUX%TOTinp - sOUT%CFLUX%TOTout)
    
!    if(isnan(sOUT % CPOOL%MBA)) then
!        write(*,'(a15,50d10.3)')"sOUT%CPOOL=",sOUT%CPOOL
!        write(*,'(a15,50d10.3)')"sOUT%CFLUX=",sOUT%CFLUX
!        stop
!    end if
    
    if(dabs(sOUT%CPOOL%TM_err).gt.1.D-10.or.sOUT%CPOOL%MBA.lt.0) then
        print*,"Negative Microbial Biomass or Balance Check ERROR: sOUT%CPOOL%TM_err=",sOUT%CPOOL%TM_err
        write(sINI%iFout_ini,'(/,a,e15.3)')"Balance Check ERROR: sOUT%CPOOL%TM_err =",sOUT%CPOOL%TM_err
        write(sINI%iFout_ini,'(/,a10,6x,100a15)')"NAME=",sINI%Name_PAR0
        write(sINI%iFout_ini,'(a10,100e15.3)')"PARxx=",xx
        write(sINI%iFout_ini,'(/,a10,6x,100a15)')"NAME=",sINI%Name_PAR
        write(sINI%iFout_ini,'(a10,100e15.3)')"PAR=",sPAR
        write(sINI%iFout_ini,'(/,a10,6x,100a15)')"NAME=",sINI%Name_POOL
        write(sINI%iFout_ini,'(a10,100e15.3)')"INP_CPOOL=",sINP%CPOOL
        write(sINI%iFout_ini,'(a10,100e15.3)')"OUT_CPOOL=",sOUT%CPOOL
        write(sINI%iFout_ini,'(/,a10,6x,100a15)')"NAME=",sINI%Name_FLUX
        write(sINI%iFout_ini,'(a10,100e15.3)')"CFLUX=",sOUT%CFLUX
        write(sINI%iFout_ini,'(/,a10,6x,100a15)')"NAME=",sINI%Name_MNFLUX
        write(sINI%iFout_ini,'(a10,100e15.3)')"MNFLUX=",sOUT%MNFLUX
        close(sINI%iFout_ini)
        stop
    end if
    !!C MASS BALANCE: END
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !!------------------------------------------------------------------------
    IF(.NOT.sINI%Carbon_only) THEN
    !!N FLUXES:BEG
    sOUT%NFLUX%POMdec           = sOUT%CFLUX%POMdec/sINP%CN%POM
    sOUT%NFLUX%POMdec_to_MOM    = sOUT%CFLUX%POMdec_to_MOM/sINP%CN%POM
    sOUT%NFLUX%POMdec_to_DOM    = sOUT%CFLUX%POMdec_to_DOM/sINP%CN%POM
    sOUT%NFLUX%MOMdec           = sOUT%CFLUX%MOMdec/sINP%CN%MOM
    sOUT%NFLUX%MOM_to_DOM       = sOUT%CFLUX%MOM_to_DOM/sINP%CN%MOM
    sOUT%NFLUX%QOM_to_DOM       = sOUT%CFLUX%QOM_to_DOM/sINP%CN%QOM
    sOUT%NFLUX%DOM_to_QOM       = sOUT%CFLUX%DOM_to_QOM/sINP%CN%DOM
    sOUT%NFLUX%DOM_to_QOM_net   = sOUT%NFLUX%DOM_to_QOM - sOUT%NFLUX%QOM_to_DOM    
    sOUT%NFLUX%MBA_to_DOM       = sOUT%CFLUX%MBA_to_DOM/sINP%CN%MBA
    sOUT%NFLUX%MBA_to_POM       = sOUT%CFLUX%MBA_to_POM/sINP%CN%MBA
    sOUT%NFLUX%MBA_mortality    = sOUT%CFLUX%MBA_mortality/sINP%CN%MBA
    sOUT%NFLUX%MBA_to_ENZP      = sOUT%CFLUX%MBA_to_ENZP/sINI%CN_ENZP
    sOUT%NFLUX%MBA_to_ENZM      = sOUT%CFLUX%MBA_to_ENZM/sINI%CN_ENZM
    sOUT%NFLUX%MBA_PM           = sOUT%NFLUX%MBA_mortality + sum(sOUT%NFLUX%MBA_to_ENZP) + sOUT%NFLUX%MBA_to_ENZM
    sOUT%NFLUX%ENZP_to_DOM      = sOUT%CFLUX%ENZP_to_DOM/sINI%CN_ENZP
    sOUT%NFLUX%ENZM_to_DOM      = sOUT%CFLUX%ENZM_to_DOM/sINI%CN_ENZM
    sOUT%NFLUX%MBA_to_MBD       = sOUT%CFLUX%MBA_to_MBD/sINP%CN%MBA
    sOUT%NFLUX%MBD_to_MBA       = sOUT%CFLUX%MBD_to_MBA/sINP%CN%MBD
    
    !!N FLUXES: MICROBIAL UPTAKE  !!PLANT UPTAKE is IGNORED CURRENTLY
!    rtp1 = 1.d0 + sINP%CPOOL%MBA/sPAR%KsDOM + sINP%CPOOL%DOM/sPAR%KsDOM &
!                + sINP%MNPOOL%NH4/sPAR%KsNH4_MB + sINP%MNPOOL%NO3/sPAR%KsNO3_MB
    !!sOUT%NFLUX%DOM_to_MBA  = (sPAR%Vg + sPAR%Vm)/sPAR%Yg * sINP%CPOOL%DOM * sINP%NPOOL%MBA/(sPAR%KsDOM*rtp1)
    sOUT%NFLUX%DOM_to_MBA   = sOUT%CFLUX%DOM_to_MBA/SINP%CN%DOM 
    
    !!MINERALIZATION: FIRST TIME
!    phi = 0.5D0*(sINP%MNPOOL%NH4/(sINP%MNPOOL%NH4 + sPAR%KsNH4_MB)+sINP%MNPOOL%NO3/(sINP%MNPOOL%NO3 + sPAR%KsNO3_MB))
!    phi = ((sINP%CN%MBA - sINI%CN_MB_min)/(sINI%CN_MB_max - sINI%CN_MB_min))**sPAR%wdie  !!higher phi means N limited
    phi = fNLimit_MB(sINP%CN%MBA,sINI%CN_MB_min,sINI%CN_MB_max,sPAR%wdorm)  !!wdorm or wdie
    
!    phi = fNLimit_MB(sINP%CN%DOM,sINI%CN_MB_min,sINI%CN_MB_max,sPAR%wdorm)
    sPAR%YgN        = phi  !![5/25/2017]  sPAR%YgN=sPAR%YgN * phi 
    sPAR%VNup_MB    = sPAR%VNup_MB * phi  !!5/25/2017
    sOUT%MNFLUX%Nmn_MBA = (1.D0 - sPAR%YgN)*sOUT%NFLUX%DOM_to_MBA
    !!IMMOBILIZATION: FIRST TIME
    rtp1 = 1.d0 + sINP%MNPOOL%NH4/sPAR%KsNH4_MB + sINP%MNPOOL%NO3/sPAR%KsNO3_MB
    sOUT%MNFLUX%Nim_NH4 = sPAR%VNup_MB * sINP%MNPOOL%NH4 &
                        * sINP%NPOOL%MBA/(sPAR%KsNH4_MB*(rtp1 + sINP%NPOOL%MBA/sPAR%KsNH4_MB))
    sOUT%MNFLUX%Nim_NO3 = sPAR%VNup_MB * sINP%MNPOOL%NO3 &
                        * sINP%NPOOL%MBA/(sPAR%KsNO3_MB*(rtp1 + sINP%NPOOL%MBA/sPAR%KsNO3_MB))
    sOUT%MNFLUX%Nim     = sOUT%MNFLUX%Nim_NH4 + sOUT%MNFLUX%Nim_NO3
    !!N FLUXES:END
    !!------------------------------------------------------------------------
    !!N MASS BALANCE: BEG
    sOUT % NPOOL % POM      = sINP % NPOOL % POM + sOUT % NFLUX % POMadd - sOUT % NFLUX % POMdec !Lignin and Cellulose
    sOUT % NPOOL % POM(1)   = sOUT % NPOOL % POM(1) + sOUT % NFLUX % MBA_to_POM !ONLY Lignin!     
    sOUT % NPOOL % MOM      = sINP % NPOOL % MOM + sum(sOUT % NFLUX % POMdec_to_MOM) - sOUT%NFLUX%MOMdec!!sOUT % CFLUX % MOM_to_DOM
    
    IF(sINI%iKinetics.lt.10) THEN
        sOUT % NPOOL % QOM      = sINP % NPOOL % QOM + sOUT % NFLUX % DOM_to_QOM_net !QOM mass balance
        sOUT % NPOOL % DOM      = sINP % NPOOL % DOM + sOUT % NFLUX % DOMadd &
        &                       + sum(sOUT % NFLUX % POMdec_to_DOM) &
        &                       + sOUT % NFLUX % MOM_to_DOM &
        &                       + sOUT % NFLUX % MBA_to_DOM &
        &                       + sum(sOUT % NFLUX % ENZP_to_DOM) &
        &                       + sOUT % NFLUX % ENZM_to_DOM &
        &                       - sOUT % NFLUX % DOM_to_MBA &
        &                       - sOUT % NFLUX % DOM_to_QOM_net
        sOUT % NPOOL % ENZP     = sINP % NPOOL % ENZP + sOUT % NFLUX % MBA_to_ENZP - sOUT % NFLUX % ENZP_to_DOM
        sOUT % NPOOL % ENZM     = sINP % NPOOL % ENZM + sOUT % NFLUX % MBA_to_ENZM - sOUT % NFLUX % ENZM_to_DOM
        sOUT % NPOOL % MBA      = sINP % NPOOL % MBA + sOUT % NFLUX % DOM_to_MBA &
        &                       - sOUT % NFLUX % MBA_PM &
        &                       - sOUT % NFLUX % MBA_to_MBD + sOUT % NFLUX % MBD_to_MBA &
        &                       + sOUT % MNFLUX % Nim &
        &                       - sOUT%MNFLUX%Nmn_MBA
        sOUT % NPOOL % MBD      = sINP % NPOOL % MBD &
        &                       + sOUT % NFLUX % MBA_to_MBD - sOUT % NFLUX % MBD_to_MBA
        
!        if(isnan(sOUT % NPOOL % MBA)) then
!        write(*,*)"NaN: ",sINP % NPOOL % MBA, sOUT % NFLUX % DOM_to_MBA, &
!        &                       sOUT % NFLUX % MBA_PM,&
!        &                       sOUT % NFLUX % MBA_to_MBD, sOUT % NFLUX % MBD_to_MBA,&
!        &                       sOUT % MNFLUX % Nim,&
!        &                       sOUT%MNFLUX%Nmn_MBA
!        STOP
!        end if
    ELSE IF(sINI%iKinetics.eq.11) THEN
        sOUT % NPOOL % MOM  = sOUT % NPOOL % MOM + sOUT % NFLUX % DOMadd
    END IF!!IF(sINI%iKinetics.lt.10)
    
    !!FIRST BALANCE CHECK: N: BEG
    sOUT % MNPOOL % NH4 = sINP%MNPOOL%NH4 + sOUT%MNFLUX%NH4_dep - sOUT%MNFLUX%Nim_NH4 + sOUT%MNFLUX%Nmn_MBA !!+ sOUT%MNFLUX%Nmn_MBD
                        
    sOUT % MNPOOL % NO3 = sINP%MNPOOL%NO3 + sOUT%MNFLUX%NO3_dep - sOUT%MNFLUX%Nim_NO3
    CALL subMEND_NPOOL_UPDATE1(sOUT%NPOOL, sOUT%MNPOOL)
    
    sOUT%NFLUX%TOTout = 0.D0
    sOUT%NPOOL%TM_err = (sOUT%NPOOL%TM - sINP%NPOOL%TM) - (sOUT%NFLUX%TOTinp - sOUT%NFLUX%TOTout)
    if(dabs(sOUT%NPOOL%TM_err).gt.1.D-8) then
        print*,"1ST Balance Check ERROR: sOUT%NPOOL%TM_err=",sOUT%NPOOL%TM_err
    end if
    !!FIRST BALANCE CHECK: N: END
    !!N MASS BALANCE: END
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    !!CHECK MBC:MBN
    !!(1) MBD
    sOUT%MNFLUX%CO2_ovflow_MBD = 0.d0
    sOUT%MNFLUX%Nmn_MBD = 0.d0
    rtp1 = sOUT % CPOOL % MBD/sOUT % NPOOL % MBD  !!CN ratio of MBD
    if(rtp1.lt.sINI%CN_MB_min) then  !!C limited, N mineralization
        rtp2 = sOUT%CPOOL%MBD/sINI%CN_MB_min 
        sOUT%MNFLUX%Nmn_MBD = sOUT%NPOOL%MBD - rtp2
        sOUT%NPOOL%MBD = rtp2
    else if(rtp1.gt.sINI%CN_MB_max) then !!N limited, overflow
        rtp2 = sOUT%NPOOL%MBD*sINI%CN_MB_max
        sOUT%MNFLUX%CO2_ovflow_MBD = sOUT%CPOOL%MBD - rtp2
        sOUT%CPOOL%MBD = rtp2    
    end if
    
    !!(2) MBA: STEP1
    sOUT%MNFLUX%CO2_ovflow_MBA = 0.d0
!    sOUT%MNFLUX%Nmn_MBA = 0.d0
    rtp1 = sOUT % CPOOL % MBA/sOUT % NPOOL % MBA  !!CN ratio of MBA
    if(rtp1.lt.sINI%CN_MB_min) then  !!C limited, N mineralization
        rtp2 = sOUT%CPOOL%MBA/sINI%CN_MB_min 
!        sOUT%MNFLUX%Nmn_MBA = sOUT%NPOOL%MBA - rtp2
        sOUT%MNFLUX%Nmn_MBA = sOUT%MNFLUX%Nmn_MBA + (sOUT%NPOOL%MBA - rtp2)
        sOUT%NPOOL%MBA = rtp2
    else if(rtp1.gt.sINI%CN_MB_max) then !!N limited, N IMMOBILIZATION + overflow
        rtp2 = sOUT%CPOOL%MBA/sINI%CN_MB_max
        rtp3 = rtp2 - sOUT%NPOOL%MBA  !!N deficit  
        
        !!(2) MBA: STEP2
        sOUT%MNPOOL%NH4 = sINP%MNPOOL%NH4 + sOUT%MNFLUX%NH4_dep - sOUT%MNFLUX%Nim_NH4 &
                        + sOUT%MNFLUX%Nmn_MBA + sOUT%MNFLUX%Nmn_MBD
        sOUT%MNPOOL%NO3 = sINP%MNPOOL%NO3 + sOUT%MNFLUX%NO3_dep - sOUT%MNFLUX%Nim_NO3
        rtp4 = sOUT%MNPOOL%NH4 + sOUT%MNPOOL%NO3 !!N Availability
        !!IMMOBILIZATION: SECOND TIME
        if(rtp3.le.rtp4) then !!deficit <= availability
            sOUT%MNFLUX%Nim_NH4 = sOUT%MNFLUX%Nim_NH4 + rtp3*sOUT%MNPOOL%NH4/rtp4
            sOUT%MNFLUX%Nim_NO3 = sOUT%MNFLUX%Nim_NO3 + rtp3*sOUT%MNPOOL%NO3/rtp4
            sOUT%MNFLUX%Nim     = sOUT%MNFLUX%Nim_NH4 + sOUT%MNFLUX%Nim_NO3
            sOUT%NPOOL%MBA = rtp2
        else  !!deficit > availability
            sOUT%MNFLUX%Nim_NH4 = sOUT%MNFLUX%Nim_NH4 + sOUT%MNPOOL%NH4
            sOUT%MNFLUX%Nim_NO3 = sOUT%MNFLUX%Nim_NO3 + sOUT%MNPOOL%NO3
            sOUT%MNFLUX%Nim     = sOUT%MNFLUX%Nim_NH4 + sOUT%MNFLUX%Nim_NO3
            sOUT%NPOOL%MBA = sOUT%NPOOL%MBA + rtp4  !!uptake all Mineral N
            sOUT%MNFLUX%CO2_ovflow_MBA = sOUT%CPOOL%MBA - sOUT%NPOOL%MBA*sINI%CN_MB_max
            sOUT%CPOOL%MBA = sOUT%NPOOL%MBA*sINI%CN_MB_max
        end if
    end if
    sOUT%MNFLUX%Nmn     = sOUT%MNFLUX%Nmn_MBD + sOUT%MNFLUX%Nmn_MBA
    sOUT%MNFLUX%Nmn_net = sOUT%MNFLUX%Nmn - sOUT%MNFLUX%Nim 
    !!revisit C BALANCE
    CALL subMEND_CPOOL_UPDATE1(sOUT%CPOOL)
    
    sOUT % MNFLUX % CO2_ovflow  = sOUT % MNFLUX % CO2_ovflow_MBA + sOUT % MNFLUX % CO2_ovflow_MBD
    sOUT % MNFLUX % CO2_gmo = sOUT % MNFLUX % CO2_gm + sOUT % MNFLUX % CO2_ovflow
    sOUT % MNPOOL % CO2     = sINP % MNPOOL % CO2 + sOUT % MNFLUX % CO2_gmo
    
    !!N BALANCE
    sOUT % MNPOOL % NH4 = sINP%MNPOOL%NH4 + sOUT%MNFLUX%NH4_dep + sOUT%MNFLUX%Nmn_MBA + sOUT%MNFLUX%Nmn_MBD &
                        - sOUT%MNFLUX%Nim_NH4
    sOUT % MNPOOL % NO3 = sINP%MNPOOL%NO3 + sOUT%MNFLUX%NO3_dep - sOUT%MNFLUX%Nim_NO3
    
    !!NITRIFICATION & DENITRIFICATION
    sOUT%MNFLUX%Nitrif      = sPAR%Vnit * sOUT % MNPOOL % NH4
    sOUT%MNFLUX%Denitrif    = sPAR%Vdenit * sOUT % MNPOOL % NO3
    phi = 1.D0 - fO2_scalar(sINI%sINP%SWC,sINI%porosity) !!higher SWC, lower fO2_scalar, higher phi
    sOUT%MNFLUX%Denitrif    = sOUT%MNFLUX%Denitrif &                                        !!DeNitrification
                            + sOUT%MNFLUX%Nitrif*phi    !!Nitrifier DeNitrification !!sINI%sINP%SWC/sINI%porosity
    sOUT % MNPOOL % NH4     = sOUT % MNPOOL % NH4 - sOUT%MNFLUX%Nitrif
    sOUT % MNPOOL % NO3     = sOUT % MNPOOL % NO3 + sOUT%MNFLUX%Nitrif - sOUT%MNFLUX%Denitrif

    CALL subMEND_NPOOL_UPDATE1(sOUT%NPOOL, sOUT%MNPOOL)
    CALL subMEND_CN_UPDATE(sOUT % CPOOL, sOUT % NPOOL, sOUT % CN)
    
    !!C MASS BALANCE CHECK
    sOUT%CFLUX%TOTout = sOUT%MNFLUX%CO2_gmo
    sOUT%CPOOL%TM_err = (sOUT%CPOOL%TM - sINP%CPOOL%TM) - (sOUT%CFLUX%TOTinp - sOUT%CFLUX%TOTout)
    
    if(dabs(sOUT%CPOOL%TM_err).gt.1.D-8) then
        print*,"Balance Check ERROR: sOUT%CPOOL%TM_err=",sOUT%CPOOL%TM_err
    end if
    
    !!N MASS BALANCE CHECK
    sOUT%NFLUX%TOTout = sOUT%MNFLUX%Denitrif
    sOUT%NPOOL%TM_err = (sOUT%NPOOL%TM - sINP%NPOOL%TM) - (sOUT%NFLUX%TOTinp - sOUT%NFLUX%TOTout)
    if(dabs(sOUT%NPOOL%TM_err).gt.1.D-8) then
        print*,"2ND Balance Check ERROR: sOUT%NPOOL%TM_err=",sOUT%NPOOL%TM_err
    end if
    !!N MASS BALANCE: END
    END IF !!IF(.NOT.sINI%Carbon_only) 
    !!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    
    
!    do j = 2, const_nISO !j = 1 for lightest isotope, e.g., C12
!        !Isotopic components
!
!        sOUT % CPOOLI(j) % POM      = sINP % CPOOLI(j) % POM + sINP % CADDI(j) % POMadd &
!        &                           - sOUT % CFLUX % POMdec * sINP % CPOOLIFR(j) % POM !Lignin and Cellulose                             
!        sOUT % CPOOLI(j) % POM(1)   = sOUT % CPOOLI(j) % POM(1) & 
!        &                           + sOUT % CFLUX % MB_to_POM * sINP % CPOOLIFR(j) % MBA !ONLY Lignin!    
!        sOUT % CPOOLI(j) % MOM      = sINP % CPOOLI(j) % MOM &
!        &                           + sum(sOUT % CFLUX % POMdec_to_MOM * sINP % CPOOLIFR(j) % POM) &
!        &                           - sOUT % CFLUX % MOM_to_DOM * sINP % CPOOLIFR(j) % MOM
!        sOUT % CPOOLI(j) % QOM      = sINP % CPOOLI(j) % QOM &
!        &                           + sOUT % CFLUX % DOM_to_QOM * sINP % CPOOLIFR(j) % DOM &
!        &                           - sOUT % CFLUX % QOM_to_DOM * sINP % CPOOLIFR(j) % QOM!QOM mass balance
!        sOUT % CPOOLI(j) % MB      = sINP % CPOOLI(j) % MB &
!        &                           + sOUT % CFLUX % DOM_to_MB * sINP % CPOOLIFR(j) % DOM &
!        &                           - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn + sOUT % CFLUX % MB_PM) &
!        &                           * sINP % CPOOLIFR(j) % MBA &
!        &                           - sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBD
!        sOUT % CPOOLI(j) % MBA     = sINP % CPOOLI(j) % MBA &
!        &                           + sOUT % CFLUX % DOM_to_MB * sINP % CPOOLIFR(j) % DOM &
!        &                           - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn + sOUT % CFLUX % MB_PM) &
!        &                           * sINP % CPOOLIFR(j) % MBA &
!        &                           - sOUT % CFLUX % MBA_to_MBD * sINP % CPOOLIFR(j) % MBA &
!        &                           + sOUT % CFLUX % MBD_to_MBA * sINP % CPOOLIFR(j) % MBD
!        sOUT % CPOOLI(j) % MBD     = sINP % CPOOLI(j) % MBD &
!        &                           - sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBD &
!        &                           + sOUT % CFLUX % MBA_to_MBD * sINP % CPOOLIFR(j) % MBA &
!        &                           - sOUT % CFLUX % MBD_to_MBA * sINP % CPOOLIFR(j) % MBD
!        sOUT % CPOOLI(j) % DOM      = sINP % CPOOLI(j) % DOM + sINP % CADDI(j) % DOMadd &
!        &                           + sum(sOUT % CFLUX % POMdec_to_DOM * sINP % CPOOLIFR(j) % POM) &
!        &                           + sOUT % CFLUX % MB_to_DOM * sINP % CPOOLIFR(j) % MBA &
!        &                           + sOUT % CFLUX % MOM_to_DOM * sINP % CPOOLIFR(j) % MOM &
!        &                           + sum(sOUT % CFLUX % ENZP_to_DOM * sINP % CPOOLIFR(j) % ENZP) &
!        &                           + sOUT % CFLUX % ENZM_to_DOM * sINP % CPOOLIFR(j) % ENZM &
!        &                           - sOUT % CFLUX % DOM_to_MB * sINP % CPOOLIFR(j) % DOM &
!        &                           - sOUT % CFLUX % DOM_to_QOM * sINP % CPOOLIFR(j) % DOM &
!        &                           + sOUT % CFLUX % QOM_to_DOM * sINP % CPOOLIFR(j) % QOM
!        sOUT % CPOOLI(j) % ENZP     = sINP % CPOOLI(j) % ENZP &
!        &                           + sOUT % CFLUX % MB_to_ENZP * sINP % CPOOLIFR(j) % MBA &
!        &                           - sOUT % CFLUX % ENZP_to_DOM * sINP % CPOOLIFR(j) % ENZP
!        sOUT % CPOOLI(j) % ENZM     = sINP % CPOOLI(j) % ENZM &
!        &                           + sOUT % CFLUX % MB_to_ENZM * sINP % CPOOLIFR(j) % MBA &
!        &                           - sOUT % CFLUX % ENZM_to_DOM * sINP % CPOOLIFR(j) % ENZM
!        sOUT % CPOOLI(j) % CO2      = sINP % CPOOLI(j) % CO2 &
!                                    + (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) * sINP % CPOOLIFR(j) % MBA &
!        &                           + sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBD
!
!        !!CO2_C14 or C13 flux
!        sOUT%CO2_gm_iso             = (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) * sINP%CPOOLIFR(2)%MBA &
!        &                           + sOUT % CFLUX % CO2_maintn_dorm * sINP%CPOOLIFR(2)%MBD
!        
!        sOUT % CPOOLI(j) % SOM      = sum(sOUT % CPOOLI(j) % POM) &
!        &                           + sOUT % CPOOLI(j) % MOM + sOUT % CPOOLI(j) % QOM
!
!    end do !j = 1, const_nISO
!
!    !j = 1, compute the lightest isotope
!    DO i = 1, const_nPOM
!        sOUT % CPOOLI(1) % POM(i) = sOUT % CPOOL % POM(i) - sum(sOUT % CPOOLI(2:const_nISO) % POM(i)) !Lignin and Cellulose
!        sOUT % CPOOLI(1) % ENZP(i) = sOUT % CPOOL % ENZP(i) - sum(sOUT % CPOOLI(2:const_nISO) % ENZP(i))
!    END DO
!    sOUT % CPOOLI(1) % MOM = sOUT % CPOOL % MOM - sum(sOUT % CPOOLI(2:const_nISO) % MOM)
!    sOUT % CPOOLI(1) % QOM = sOUT % CPOOL % QOM - sum(sOUT % CPOOLI(2:const_nISO) % QOM)
!    sOUT % CPOOLI(1) % MB = sOUT % CPOOL % MB - sum(sOUT % CPOOLI(2:const_nISO) % MB)
!    sOUT % CPOOLI(1) % MBA = sOUT % CPOOL % MBA - sum(sOUT % CPOOLI(2:const_nISO) % MBA)
!    sOUT % CPOOLI(1) % MBD = sOUT % CPOOL % MBD - sum(sOUT % CPOOLI(2:const_nISO) % MBD)
!    sOUT % CPOOLI(1) % DOM = sOUT % CPOOL % DOM - sum(sOUT % CPOOLI(2:const_nISO) % DOM)
!    sOUT % CPOOLI(1) % ENZM = sOUT % CPOOL % ENZM - sum(sOUT % CPOOLI(2:const_nISO) % ENZM)
!    sOUT % CPOOLI(1) % CO2 = sOUT % CPOOL % CO2 - sum(sOUT % CPOOLI(2:const_nISO) % CO2)
!    sOUT % CPOOLI(1) % SOM = sOUT % CPOOL % SOM - sum(sOUT % CPOOLI(2:const_nISO) % SOM)
!
!    do j = 1, const_nISO
!        !Isotopic fractions
!        sOUT % CPOOLIFR(j) % POM = sOUT % CPOOLI(j) % POM/sOUT % CPOOL % POM
!        sOUT % CPOOLIFR(j) % MOM = sOUT % CPOOLI(j) % MOM/sOUT % CPOOL % MOM
!        sOUT % CPOOLIFR(j) % QOM = sOUT % CPOOLI(j) % QOM/sOUT % CPOOL % QOM
!        sOUT % CPOOLIFR(j) % MB = sOUT % CPOOLI(j) % MB/sOUT % CPOOL % MB
!        sOUT % CPOOLIFR(j) % MBA = sOUT % CPOOLI(j) % MBA/sOUT % CPOOL % MBA
!        sOUT % CPOOLIFR(j) % MBD = sOUT % CPOOLI(j) % MBA/sOUT % CPOOL % MBD
!        sOUT % CPOOLIFR(j) % DOM = sOUT % CPOOLI(j) % DOM/sOUT % CPOOL % DOM
!        sOUT % CPOOLIFR(j) % ENZP = sOUT % CPOOLI(j) % ENZP/sOUT % CPOOL % ENZP
!        sOUT % CPOOLIFR(j) % ENZM = sOUT % CPOOLI(j) % ENZM/sOUT % CPOOL % ENZM
!        sOUT % CPOOLIFR(j) % CO2 = sOUT % CPOOLI(j) % CO2/sOUT % CPOOL % CO2
!        sOUT % CPOOLIFR(j) % SOM = sOUT % CPOOLI(j) % SOM/sOUT % CPOOL % SOM
!
!    end do !j = 1, const_nISO        
!
!    !    print*, sum(sOUT%CPOOLIFR%MB)!sum(sOUT%CPOOLI(2:const_nISO)%MB),sOUT%CPOOLI(1:const_nISO)%MB
!
!    !Convert Isotope Concentration to Signatures [Permil]
!    do j = 1, const_nISO - 1
!
!        do k = 1, const_nPOM
!            sOUT % CPOOLI_SIG(j) % POM(k) = &
!            & fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % POM(k), sOUT % CPOOLI(j + 1) % POM(k))
!            sOUT % CPOOLI_SIG(j) % ENZP(k) = &
!            & fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % ENZP(k), sOUT % CPOOLI(j + 1) % ENZP(k))
!        end do
!
!        sOUT % CPOOLI_SIG(j) % MOM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MOM, sOUT % CPOOLI(j + 1) % MOM)
!        sOUT % CPOOLI_SIG(j) % QOM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % QOM, sOUT % CPOOLI(j + 1) % QOM)
!        sOUT % CPOOLI_SIG(j) % MB = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MB, sOUT % CPOOLI(j + 1) % MB)
!        sOUT % CPOOLI_SIG(j) % MBA = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBA, sOUT % CPOOLI(j + 1) % MBA)
!        sOUT % CPOOLI_SIG(j) % MBD = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBD, sOUT % CPOOLI(j + 1) % MBD)
!        sOUT % CPOOLI_SIG(j) % DOM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % DOM, sOUT % CPOOLI(j + 1) % DOM)
!        sOUT % CPOOLI_SIG(j) % ENZM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % ENZM, sOUT % CPOOLI(j + 1) % ENZM)
!        sOUT % CPOOLI_SIG(j) % CO2 = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % CO2, sOUT % CPOOLI(j + 1) % CO2)
!        sOUT % CPOOLI_SIG(j) % SOM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % SOM, sOUT % CPOOLI(j + 1) % SOM)
!    end do !j = 1, const_nISO - 1

    
!    END ASSOMIATE
END !!subroutine subMEND!MEND model

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_PAR(xx, sPAR, sINI)
!    MEND MODEL CODE 
!    USE MOD_MEND
!    IMPLICIT NONE
!    REAL(8), PARAMETER :: Tref = 20.d0  !!moved to MOD_MEND_TYPE
!    REAL(8) fTArh, fT_CUE !!FUNCTIONS
!    REAL(8) fSWP, fSWP_OPT, fSWP_A2D, fSWP_D2A !!FUNCTIONS
!    REAL(8) fpH  !!FUNCTIONS
!    ATTENTION: currently the following parameters are NOT used: xx(6): KPC (=KPL), xx(16): pEM (=pEP)
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(inout) :: sPAR
    TYPE(sMEND_INI), intent(inout) :: sINI
    REAL(8)        , intent(in)    :: xx(sINI%nPar)
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_PAR) sPARbase  !!base PAR value = xx(:)
    REAL(8) CUE_slope,CUE_ref, SWPmin !!CUE_slope = (-0.005, -0.016)  
    REAL(8) Vm0  !!specific maintenance rate without any modification
    
    REAL(8) tp_scalar       !!temperature scalar
    REAL(8) wp_scalar       !!INCrease with increasing wp 
    REAL(8) wp_scalar_low   !!DECrease with increasing wp 
    REAL(8) wp_scalar_opt   !!with increasing wp, INcrease first to OPTimal condition, then decrease
    REAL(8) pH_scalar   !!pH scalar
    
    REAL(8) tp
    REAL(8) wp, wfp !!wp=water matric potential, wfp= water-filled porosity
    REAL(8) pH
    
    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    CHARACTER(len=20) sCase
    
    SWPmin = -13.86d0  !!for Microbial Mortality
    
!    sPAR%iKinetics = sINI%iKinetics
    BIOME = trim(sINI%BIOME)
    SOM   = trim(sINI%SOM)
    
    tp = sINI%sINP%tmp
    wp = sINI%sINP%SWP
    wfp = min(1.D0,sINI%sINP%SWC/sINI%porosity)
    pH = sINI%sINP%pH
    
    !!sPARbase
    sPARbase%fINP       = xx(3)
    sPARbase%VdPOM(1)   = xx(4)
    sPARbase%VdPOM(2)   = xx(4) !!xx(4)
    sPARbase%VdMOM      = xx(4) !!xx(5)
    sPARbase%KsPOM(1)   = xx(6)
    sPARbase%KsPOM(2)   = xx(6)*0.1  !!xx(7)
    sPARbase%KsMOM      = xx(6)*10.0 !!xx(8)
    sPARbase%Qmax       = xx(9)
    sPARbase%Kba        = xx(10)
    sPARbase%Kdes       = xx(11)
    sPARbase%Kads       = sPARbase%Kdes*sPARbase%Kba
    sPARbase%rENZM      = xx(12)
    sPARbase%rENZP      = sPARbase%rENZM
    sPARbase%pENZP      = xx(13)
    sPARbase%pENZM      = sPARbase%pENZP*xx(14)
    sPARbase%frPOM2DOM  = xx(15)
    sPARbase%frMB2DOM   = xx(16)
    sPARbase%Vg         = xx(17)
    sPARbase%alpha      = xx(18)
    sPARbase%Vm         = sPARbase%Vg * sPARbase%alpha/(1.D0 - sPARbase%alpha)
    sPARbase%KsDOM      = xx(19)
    sPARbase%Yg         = xx(20)
!    sPARbase%Ygsl       = xx(21)  !!no Ygsl in sPAR
    CUE_slope           = -1.D0*xx(21)
    sPARbase%Q10        = xx(22)  
    sPARbase%gamma      = xx(23)
    sPARbase%rMORT      = min(0.99D0,sPARbase%Vm*sPARbase%gamma)
    sPARbase%beta       = xx(24)
    sPARbase%VmD        = sPARbase%Vm*sPARbase%beta
    sPARbase%SWP_A2D    = xx(25)
    sPARbase%tau        = xx(26)
    sPARbase%wdorm      = xx(27)
    
    sPARbase%VNup_MB    = xx(28)
    sPARbase%KsNH4_MB   = xx(29)
    sPARbase%KsNO3_MB   = xx(30)
    sPARbase%Vnit       = xx(31)
    sPARbase%Vdenit     = xx(32)
    sPARbase%YgN        = xx(33)
    
    sPAR % fINP = sPARbase % fINP
    sPAR % Q10 = sPARbase % Q10  
    !![1] DECOMPOSITION of POM1 & POM2
    sCase = "LIG"
    if(sINI%iTmp_Func.eq.0) then
        tp_scalar = fTArh(sCase, tp, const_Tref)
    else
        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
    end if
    wp_scalar_opt   = fSWP_OPT(wp)
    pH_scalar       = fpH(sCase,pH)
    sPAR % VdPOM(1)     = sPARbase % VdPOM(1)*tp_scalar*wp_scalar_opt*pH_scalar !!(/xx(3)/) ![mg POM/mg ENZP/h], maximum reaction rate for conversion of POM by ENZP
    
    sCase = "CEL"
    if(sINI%iTmp_Func.eq.0) then
        tp_scalar = fTArh(sCase, tp, const_Tref)
    else
        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
    end if
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
    sPAR % VdPOM(2)     = sPARbase % VdPOM(2)*tp_scalar*wp_scalar*pH_scalar  !!xx(4) !!set Vd_cel = Vd_lig
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    
    sPAR % KsPOM        = sPARbase % KsPOM*tp_scalar ![mg POM/cm3], half-saturation constant for conversion of POM by ENZP (/xx(6), xx(7)/)

    !![2] DECOMPOSITION of MOM
    sCase = "ENZ"
    if(sINI%iTmp_Func.eq.0) then
        tp_scalar = fTArh(sCase, tp, const_Tref)
    else
        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
    end if
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
    sPAR % VdMOM        = sPARbase % VdMOM*tp_scalar*wp_scalar*pH_scalar ![mg MOM/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC xx(5)
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % KsMOM        = sPARbase % KsMOM*tp_scalar ![mg MOM/cm3], half-saturation constant for conversion of MAOC by ENZMAOC xx(8)

    !![3] ADSORPTION-DESORPTION
    sPAR % Qmax         = sPARbase % Qmax !!xx(9) ![mg C/g soil], adsorption capacity
    sPAR % Kba          = sPARbase % Kba !!xx(10) ![mg C/g soil/h], binding affinity
    
    sCase = "Kdes"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % Kdes         = sPARbase % Kdes*tp_scalar ![mg DOM/h],desorption rate constant  xx(11)
    
    sCase = "Kads"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % Kads         = sPARbase % Kads*tp_scalar ![mg DOM/mg DOM/h], adsorption rate constant = Kba*Kdes
    
    !![4] ENZYME TURNOVER & PRODUCTION
    sPAR % rENZP        = sPARbase % rENZP!!(/xx(12), xx(12)/) !![1/h],turnover rate of ENZP
    sPAR % rENZM        = sPARbase % rENZM!!xx(12) !![1/h], turnover rate of ENZMAOC
    sPAR % pENZP        = sPARbase % pENZP!!xx(13) ![mg ENZP/mg MB/h], production rate of ENZP 
    sPAR % pENZM        = sPARbase % pENZM!!xx(13)*xx(14) !!wgs:6/22/2015: set pENZM = xx(16)*pENZP, xx(16) ![mg ENZM/mg MB/h], production rate of ENZMAOC
    
    !![5] ALLOCATION to DOM
    sPAR % frPOM2DOM    = sPARbase % frPOM2DOM !!xx(15) ![-], fraction of decomposed POM allocated to DOM
    sPAR % frMB2DOM     = sPARbase % frMB2DOM !!xx(16) ![-], fraction of dead MB allocated to DOM
   
    !![6] MICROBIAL UPTAKE of DOM
    sCase = "DOM"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    wp_scalar = fSWP(wp, BIOME, SOM)  !!same as Cellulose-decomposition
    sPAR % Vg        = sPARbase % Vg*tp_scalar !!*wp_scalar !![mg DOM/mg MB/h], maximum uptake rate of DOM by MB  xx(17)
    sPAR % alpha     = sPARbase % alpha!!xx(18) ![-], alpha = Vm/(Vm + Vg) <= 0.5 
    
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
!    Vm0 = xx(17) * sPAR%alpha/(1.D0 - sPAR%alpha)
    sPAR % Vm         = sPARbase % Vm*tp_scalar !!*wp_scalar ![1/h], specific microbial maintenance rate
     
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % KsDOM        = sPARbase % KsDOM*tp_scalar ![mg DOM/cm3], half-saturation constant for uptake of DOM by MB  xx(19)
    
    sCase = "CUE"  !!label, no-use
    CUE_ref   = sPARbase%Yg !!xx(20)
!    CUE_slope = -1.D0*xx(21)
    sPAR % Yg = fT_CUE(tp, const_Tref, CUE_slope, CUE_ref) ![-], carbon use efficiency in uptake of DOM by MB 
    
    !![7] MICROBIAL MORTALITY
!    sPAR%wdie = sPARbase%wdorm
    sPAR%gamma= sPARbase%gamma !!xx(23) 
    !! use the same fSWP as fSWP_A2D
!    if(sINI%iSWP_die.eq.0) then
    wp_scalar_low = fSWP_A2D(wp, sPARbase % SWP_A2D, sPARbase%wdorm)  !!negative effect  !!wdie
!    else
!        wp_scalar_low =  fSWP_Death(wp,SWPmin,sPAR%wdie)  
!    end if
    
!    sPAR % rMORT = (sPAR%gamma*Vm0)*wp_scalar_low !!wgs: 6/19/2015,  !!tp_scalar
    sPAR % rMORT = sPARbase % rMORT*wp_scalar_low
       
    !![8] MICROBIAL DORMANCY & RESUSCITATION
    sPAR % beta       = sPARbase % beta !!xx(24) !beta = Vm_dormant/Vm
    sPAR % VmD        = sPARbase % VmD!!sPAR % beta*Vm0 !!*wp_scalar_low  !!<5/15/2015>
    sPAR % SWP_A2D    = sPARbase % SWP_A2D!!xx(25)
    sPAR % tau        = sPARbase % tau!!xx(26)
    sPAR % SWP_D2A    = sPARbase % tau * sPARbase % SWP_A2D 
    sPAR % wdorm      = sPARbase % wdorm!!xx(27)
    
    wp_scalar     = fSWP_D2A(wp, sPAR % SWP_D2A, sPAR % wdorm)  !!positive effect, increase with increasing SWC or SWP
    wp_scalar_low = fSWP_A2D(wp, sPAR % SWP_A2D, sPAR % wdorm)  !!negative effect
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % VmA2D = sPARbase%Vm * tp_scalar*wp_scalar_low 
    sPAR % VmD2A = sPARbase%Vm * tp_scalar*wp_scalar
    
    !![9] MICROBIAL UPTAKE of NH4 & NO3
    sCase = "NH4"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % VNup_MB      = sPARbase % VNup_MB*tp_scalar !![mg N/mg MBN/h], maximum uptake rate of NH4 or NO3 by MB  xx(28)
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    sPAR % KsNH4_MB     = sPARbase % KsNH4_MB*tp_scalar ![mg N/cm3], half-saturation constant xx(29)
    sPAR % KsNO3_MB     = sPARbase % KsNH4_MB*tp_scalar ![mg N/cm3], half-saturation constant xx(30)
    
    !![10] NITRIFICATION & DENITRIFICATION
    sCase = "NITRIF"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    wp_scalar_opt   = fWFP_PieceWise(sCase,wfp)
    sPAR % Vnit     = sPARbase % Vnit*tp_scalar*wp_scalar_opt !!xx(31)
    
    sCase = "DENITRIF_NO3"
    tp_scalar = fTArh(sCase, tp, const_Tref)
!    if(sINI%iTmp_Func.eq.0) then
!        tp_scalar = fTArh(sCase, tp, const_Tref)
!    else
!        tp_scalar = fTQ10_MEND(sCase,tp, const_Tref, sPAR % Q10)
!    end if
    wp_scalar_opt   = fWFP_PieceWise(sCase,wfp)
    sPAR % Vdenit   = sPARbase % Vdenit*tp_scalar*wp_scalar_opt  !!xx(32)

    sPAR % YgN      = sPARbase % YgN!!xx(33)  
END !!subroutine subMEND_PAR
!-----------------------------------------------------------------------------

!==================================================================
!! MODEL KINETICS
!------------------------------------------------------------------  
REAL(8) FUNCTION fV_MM(iType, sPAR, sINP)
    !Michaelis-Menten Kinetics (Equation)
    !!ARGUMENTS:
    INTEGER,       INTENT(IN) :: iType
    TYPE(sMM_INP), INTENT(IN) :: sINP
    TYPE(sMM_PAR), INTENT(IN) :: sPAR
    
    SELECT CASE (iType)
        CASE (1,11)  !!First Order
            fV_MM = sPAR%vm * sINP%substrate
        CASE (2)  !!Second Order
            fV_MM = sPAR%vm * sINP%enzyme * sINP%substrate
        CASE DEFAULT !!M-M Kinetics, iType = 0
            fV_MM = sPAR % vm * sINP % substrate * sINP % enzyme/(sPAR % km + sINP % substrate)  
    END SELECT
    fV_MM = min(fV_MM, sINP % substrate)
    return
END !!FUNCTION fV_MM
!------------------------------------------------------------------       
INTEGER FUNCTION fAds(sPAR, sINP, sOUT)
    !Adsorption and Desorption of DOM
    !eg: adsorbent: QOM, adsorbate: DOM
    !!ARGUMENTS:
    TYPE(sSORP_INP), INTENT(IN) :: sINP
    TYPE(sSORP_PAR), INTENT(IN) :: sPAR
    TYPE(sSORP_OUT), INTENT(OUT) :: sOUT
    
    sOUT % ads = sPAR % Kads * sINP % adsorbate * (1.d0 - sINP % adsorbent/sPAR % Qmax)
    sOUT % des = sPAR % Kdes * sINP % adsorbent/sPAR % Qmax
    if (sOUT % des > (sINP % adsorbent + sOUT % ads)) then
        sOUT % des = sINP % adsorbent + sOUT % ads
    elseif(sOUT % ads > (sINP % adsorbate + sOUT % des)) then
        sOUT % ads = sINP % adsorbate + sOUT % des
    end if
    sOUT % ads_net = sOUT % ads - sOUT % des
    return
END !!FUNCTION fAds
!------------------------------------------------------------------

!!=================================================================
!! TEMPERATURE SCALAR: BEGIN
!------------------------------------------------------------------
REAL(8) function fTArrhenius(tmp, Ea)
    !!ARGUMENTS:
    REAL(8)tmp,Ea
    !Ea [kJ/mol]: Arrhenius activation energy
    !temp (C): temperature
    fTArrhenius = dexp(-Ea * 1d3/(const_R * (tmp + const_tmp_C2K)))
    return
END !!function fTArrhenius
!------------------------------------------------------------------
REAL(8) function fTArh0(T, Tref, Ea)
!    USE MOD_MEND
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !Ea [KJ/mol]: Arrhenius activation energy
    !temp (C): temperature
    !!ARGUMENTS:
    REAL(8), intent(in) :: T, Tref, Ea
    REAL(8) TKref, TK
    TKref = Tref + const_tmp_C2K
    TK = T + const_tmp_C2K
    fTArh0 = dexp(Ea*1.D3/const_R * (1.d0/TKref - 1.d0/TK))
    return
END !!function fTArh0

!------------------------------------------------------------------
REAL(8) function fTArh(sCase, T, Tref)
!    USE MOD_MEND
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !Ea [J/mol]: Arrhenius activation energy
    !temp (C): temperature
    !!ARGUMENTS:  
    CHARACTER(*), intent(in) :: sCase
    REAL(8)     , intent(in) :: T, Tref
    
    !!LOCAL VARIABLES:
    REAL(8) TKref, TK, Ea
    
!    sCase = trim(sCase)
    SELECT CASE (trim(sCase))
        CASE ("BG")  !! Beta-glucosidase
            Ea = 42.2d0
        CASE ("CBH") !! Cellobiohydrolase
            Ea = 32.2d0
        CASE ("EG")  !! Endo-glucanase
            Ea = 34.4d0
        CASE ("PER") !! Peroxidase
            Ea = 52.9d0
        CASE ("POX") !! Phenol oxidase
            Ea = 53.1d0
        CASE ("LIG") !! Ligninases
            Ea = 53.0d0
        CASE ("CEL") !! Cellulases
            Ea = 36.3d0
        CASE ("NAG") !! N-acetylglutamate synthase
            Ea = 52.9d0
        CASE ("PAC") !! Acid phosphatases
            Ea = 36.3d0
        CASE ("PAL") !! Alkaline phosphatases
            Ea = 23.6d0
        CASE ("PHO") !! PHOSPHATASES
            Ea = 34.4d0
        CASE ("Km")  !! half-saturation constant
            Ea = 30d0
        CASE ("MR")  !! Microbial Maintenance
            Ea = 20d0
        CASE ("Kads")!! adsorption
            Ea = 5d0
        CASE ("Kdes")!! desorption
            Ea = 20d0
        CASE ("DOM") !! DOM uptake
            Ea = 47d0
        CASE DEFAULT
            Ea = 47d0
    END SELECT
    fTArh = fTArh0(T, Tref, Ea)
    return
END !!function fTArh
!------------------------------------------------------------------
REAL(8) function fT_Linear(T, Tref, slope, intercept)
    !fKmT: [mg/m3], half-saturation constant in M-M kinetics (Km)modified by Temperature
    !slope: [mg/m3/C]
    !intercept: [mg/m3]
    !!ARGUMENTS:
    REAL(8), intent(in) :: T, Tref, slope, intercept
    
    fT_Linear = slope * (T - Tref) + intercept
    return
END
!------------------------------------------------------------------
REAL(8) function fT_CUE(T,Tref,slope,CUEref)
    !fT_CUE: [-], half-saturation constant in M-M kinetics (Km)modified by Temperature
    !slope: [1/degree C]
    !intercept: [-]
    !! parameter values: Wang et al. (2015), ISMEJ 
    !!ARGUMENTS:
    REAL(8), intent(in):: T, Tref,slope,CUEref
!    REAL(8), PARAMETER :: Tref      = 0    !![degree C]
!    REAL(8), PARAMETER :: slope     = -0.01
!    REAL(8), PARAMETER :: intercept = 0.56  !! []
    
    !!LOCAL VARIABLES:
    REAL(8), PARAMETER:: CUEmax = 0.9D0
    REAL(8), PARAMETER:: CUEmin = 0.1D-1
    
    fT_CUE = fT_Linear(T, Tref, slope, CUEref)
    if(fT_CUE.gt.CUEmax) then
        fT_CUE = CUEmax
    elseif(fT_CUE.lt.CUEmin) then
        fT_CUE = CUEmin
    end if
    
    return
END
!------------------------------------------------------------------
REAL(8) function fTAG(tmp)
    !Arrhenius Equation used in ECOSYS
    !Grant, et al., 1993. SBB 25, 1317-1329.
    !!ARGUMENTS:
    REAL(8) tmp
    
    !!LOCAL VARIABLES:
    REAL(8) A, S, Ha, Hdl, Hdh, b1, b2, a1, a2, a3
    A = 17.1124D0 !a parameter selected such that fmt=1 at temp = 30 C
    S = 710D0 ![J/mol/K], the change in entropy
    Ha = 57500D0 ![J/mol], energy of activation
    Hdl = 192500D0![J/mol], energy of low temperature deactivation
    Hdh = 222500D0![J/mol], energy of high temperature deactivation
    b1 = const_tmp_C2K + tmp
    b2 = const_R * b1
    a1 = dexp(A - Ha/b2)
    a2 = dexp((Hdl - S * b1)/b2)
    a3 = dexp((-Hdh + S * b1)/b2)
    fTAG = b1 * a1/(1.0D0 + a2 + a3)
    return
END !!function fTAG

!------------------------------------------------------------------
REAL(8) function fTQ10(T, Tref, Q10)
!    USE MOD_MEND
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !temp (C): temperature
    !!ARGUMENTS:
    REAL(8), intent(in) :: T, Tref, Q10
    fTQ10 = Q10**((T-Tref)/10.D0)
    return
END !!function fTArh0

!------------------------------------------------------------------
REAL(8) function fTQ10_MEND(sCase, T, Tref,Q10_base)
!    USE MOD_MEND
    !fTvm,[-] temperature-adjusted factor for maximum reaction rate in M-M kinetics 
    !fT = 1 at T = Tref
    !Q10_base: base Q10 for Ea = 47 KJ/mol/K
    !temp (C): temperature
    !!ARGUMENTS:  
    CHARACTER(*), intent(in) :: sCase
    REAL(8)     , intent(in) :: T, Tref, Q10_base
    
    !!LOCAL VARIABLES:
!    REAL(8) TKref, TK, Ea
    REAL(8) Q10, Q10_multiplier
    
!    sCase = trim(sCase)
    SELECT CASE (trim(sCase))
        CASE ("BG")  !! Beta-glucosidase
            Q10_multiplier = 0.94d0
        CASE ("CBH") !! Cellobiohydrolase
            Q10_multiplier = 0.82d0
        CASE ("EG")  !! Endo-glucanase
            Q10_multiplier = 0.84d0
        CASE ("PER") !! Peroxidase
            Q10_multiplier = 1.08d0
        CASE ("POX") !! Phenol oxidase
            Q10_multiplier = 1.09d0
        CASE ("LIG") !! Ligninases
            Q10_multiplier = 1.08d0
        CASE ("CEL") !! Cellulases
           Q10_multiplier = 0.87d0
        CASE ("NAG") !! N-acetylglutamate synthase
            Q10_multiplier = 1.08d0
        CASE ("PAC") !! Acid phosphatases
            Q10_multiplier = 0.87d0
        CASE ("PAL") !! Alkaline phosphatases
            Q10_multiplier = 0.73d0
        CASE ("PHO") !! PHOSPHATASES
            Q10_multiplier = 0.84d0
        CASE ("Km")  !! half-saturation constant
            Q10_multiplier = 0.79d0
        CASE ("MR")  !! Microbial Maintenance
            Q10_multiplier = 0.69d0
        CASE ("Kads")!! adsorption
            Q10_multiplier = 0.57d0
        CASE ("Kdes")!! desorption
            Q10_multiplier = 0.69d0
        CASE ("DOM") !! DOM uptake
            Q10_multiplier = 1.00d0
        CASE DEFAULT
            Q10_multiplier = 1.00d0
    END SELECT
    Q10 = max(1.0d0,Q10_base*Q10_multiplier)
    fTQ10_MEND = fTQ10(T, Tref, Q10)
    return
END !!function fTArh
!!-----------------------------------------------------------------
!! TEMPERATURE SCALAR: END
!!=================================================================

!!=================================================================
!! pH Scalar: BEGIN
!------------------------------------------------------------------
Subroutine subVGparameter(SoilTexture,SWCres,SWCsat,alpha,n)
!! Soil Water Characteristic Curve
!! van-Genuchten parameters
!! Tuller M, Or D (2004) Retention of water in soil and the soil water characteristic curve. 
!! Encyclopedia of soils in the environment, 4, 278-289.

    !!ARGUMENTS:
    CHARACTER(len=*) , intent(in) :: SoilTexture
    REAL(8)          , intent(out) :: SWCres,SWCsat,alpha,n 
    
    !!LOCAL VARIABLES
    
!    sCase = trim(sCase)
    SELECT CASE (trim(SoilTexture))
        CASE ("Sand")  
            SWCres = 0.058
            SWCsat = 0.37
            alpha  = 0.035
            n      = 3.19
        CASE ("Loamy-Sand") 
            SWCres = 0.074
            SWCsat = 0.39
            alpha  = 0.035
            n      = 2.39
        CASE ("Sandy-Loam")  
            SWCres = 0.067
            SWCsat = 0.37
            alpha  = 0.021
            n      = 1.61
        CASE ("Loam") 
            SWCres = 0.083
            SWCsat = 0.46
            alpha  = 0.025
            n      = 1.31
        CASE ("Silt") 
            SWCres = 0.123
            SWCsat = 0.48
            alpha  = 0.006
            n      = 1.53
        CASE ("Silt-Loam") 
            SWCres = 0.061
            SWCsat = 0.43
            alpha  = 0.012
            n      = 1.39
        CASE ("Sandy-Clay-Loam") 
            SWCres = 0.086
            SWCsat = 0.40
            alpha  = 0.033
            n      = 1.49
        CASE ("Clay-Loam") 
            SWCres = 0.129
            SWCsat = 0.47
            alpha  = 0.030
            n      = 1.37
        CASE ("Silty-Clay-Loam") 
            SWCres = 0.098
            SWCsat = 0.55
            alpha  = 0.027
            n      = 1.41
        CASE ("Silty-Clay") 
            SWCres = 0.163
            SWCsat = 0.47
            alpha  = 0.023
            n      = 1.39
        CASE ("Clay") 
            SWCres = 0.102
            SWCsat = 0.51
            alpha  = 0.021
            n      = 1.20
        CASE DEFAULT
            SWCres = 0.095
            SWCsat = 0.45
            alpha  = 0.024
            n      = 1.66
    END SELECT
    
    return
END !!subVGparameter
!------------------------------------------------------------------
!!=================================================================
!! Soil Water Scalar: BEGIN
!!-----------------------------------------------------------------
REAL(8) function fSWP2SWC(SWP,SWP_units,SWCres,SWCsat,alpha,n)
!!van-Genuchten equation
!!convert SWP(cm) to SWC (0-1)
!!SWP will be converted to cm
    !!ARGUMENTS:
    CHARACTER(len=*) SWP_units  
    REAL(8),intent(in):: SWP,SWCres,SWCsat,alpha,n
    
    !!LOCAL VARIABLES:
    REAL(8) SWP_cm, m, eff_sat  !!effective saturation

    if (trim(SWP_units).eq."MPa") then
        SWP_cm = SWP/const_cm2MPa
    else if (trim(SWP_units).eq."bar") then  !!1 bar = 100 kPa = 0.1 MPa
        SWP_cm = SWP*0.1/const_cm2MPa
    else if (trim(SWP_units).eq."kPa") then  !!1 kPa = 1d-3 MPa
        SWP_cm = SWP*1d-3/const_cm2MPa
    else if (trim(SWP_units).eq."Pa") then  !!1 kPa = 1d-6 MPa
        SWP_cm = SWP*1d-6/const_cm2MPa
    else if (trim(SWP_units).eq."mm") then
        SWP_cm = SWP*0.1
    else if (trim(SWP_units).eq."m") then
        SWP_cm = SWP*100
    end if
    
    m = 1d0-1d0/n
    
    eff_sat = (1/(1+(alpha*dabs(SWP_cm))**n))**m
    fSWP2SWC = SWCres + (SWCsat - SWCres)*eff_sat
    return
   
END !!function fSWP2SWC
!!-----------------------------------------------------------------
REAL(8) function fSWC2SWP(SWC0,SWCres,SWCsat,alpha,n,SWPmin)
!!van-Genuchten equation, SWP in units of [cm], SWC if fraction (0-1)
!!this function converts SWP [cm] to SWP [MPa]
!!return fSWC2SWP), actually fSWC2SWP<0
!!convert SWC(0-1) to SWP(MPa)
!!SWC needs be converted to fraction first
!    USE MOD_MEND
!    IMPLICIT NONE
    REAL(8),PARAMETER::rlim = 1.01d0
    !!ARGUMENTS:
    REAL(8),intent(in):: SWC0
    REAL(8) SWCres,SWCsat,alpha,n
    
    !!LOCAL VARIABLES:
    REAL(8) SWC,SWPmin  !!min SWP <0
    REAL(8) m, eff_sat  !!effective saturation
    
    m = 1.d0-1.d0/n

    if (SWC0.le.SWCres*rlim) then
!        fSWC2SWP = SWPmin
        SWC = SWCres*rlim
    else
        SWC = SWC0
    end if

    if (SWC.lt.SWCsat) then
        eff_sat = (SWC - SWCres)/(SWCsat - SWCres)
        fSWC2SWP = (1.d0/(eff_sat**(1.d0/m)) - 1.d0)**(1.d0/n)/alpha
        fSWC2SWP = -1.d0*fSWC2SWP*const_cm2MPa
    else
        fSWC2SWP = 0
    end if

    return
END !!function fSWC2SWP
!------------------------------------------------------------------
REAL(8) function fSWP0(SWP,SWPmin,w)
    !Rate Scalar for Soil Water Potential (SWP)
    !Manzoni et al (2012), Ecology, 93: 930-938
!    USE MOD_MEND
!    IMPLICIT NONE
    
    REAL(8), PARAMETER :: SWP_FC = -0.033 ![MPa], field capacity SWP 
    !!ARGUMENTS:
    REAL(8) SWP, SWPmin ![MPa]
    REAL(8) w
    
    if (SWP.lt.SWPmin) then
        fSWP0 = 0.0d0  
    else if (SWP.lt.SWP_FC) then
        fSWP0 = 1-(dlog(SWP/SWP_FC)/dlog(SWPmin/SWP_FC))**w
    else 
        fSWP0 = 1.0d0
    end if
    return
END !!function fSWP0

!------------------------------------------------------------------
REAL(8) function fSWP(SWP,BIOME,SOM)
    !SWP Scalar for SOM (Cellulose) decomposition
    !Manzoni et al (2012), Ecology, 93: 930-938
    !!fSWP increases with increasing SWP (wetter condition)
    
    !REAL(8), PARAMETER :: SWP_FC = -0.033 ![MPa], field capacity SWP 

    !!ARGUMENTS:
    REAL(8)         , intent(in) :: SWP ![MPa]
    CHARACTER(len=3), intent(in) :: BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3), intent(in) :: SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    
    !!LOCAL VARIABLES:
    REAL(8) SWPmin, w ![MPa]
    
    if (trim(SOM).eq."SOD") then
        SWPmin = -1.71d0
        w      = 1.43d0
    else if (trim(SOM).eq."SOL") then
        SWPmin = -13.86d0
        w      = 1.20d0
    else if (trim(SOM).eq."LIT") then
        SWPmin = -36.49d0
        w      = 1.04d0
    end if
    
    SELECT CASE (trim(BIOME)) 
        CASE ("ASM")  !!Arid, Semiarid, & Mediterranean
            if (trim(SOM).eq."SOL") then
                SWPmin = -10.95d0
                w      = 1.26d0
            end if
        CASE ("MGC")  !!Mesic Grassland & Cropland
            if (trim(SOM).eq."SOD") then
                SWPmin = -1.71d0
                w      = 1.43d0
            else if (trim(SOM).eq."SOL") then
                SWPmin = -22.61d0
                w      = 1.11d0
            else if (trim(SOM).eq."LIT") then
                SWPmin = -39.73d0
                w      = 0.89d0
            end if
        CASE ("MDF")    !!Mesic Deciduous Forest
            if (trim(SOM).eq."SOL") then
                SWPmin = -4.97d0
                w      = 1.07d0
            else if (trim(SOM).eq."LIT") then
                SWPmin = -29.00d0
                w      = 1.27d0
            end if
        CASE ("MCF")    !!MEsic Conifer Forest
            if (trim(SOM).eq."SOL") then
                SWPmin = -8.24d0
                w      = 1.40d0
            else if (trim(SOM).eq."LIT") then
                SWPmin = -39.85d0
                w      = 1.06d0
            end if
        CASE DEFAULT    !!All Biome Average
            if (trim(SOM).eq."SOD") then
                SWPmin = -1.71d0
                w      = 1.43d0
            else if (trim(SOM).eq."SOL") then
                SWPmin = -13.86d0
                w      = 1.20d0
            else if (trim(SOM).eq."LIT") then
                SWPmin = -36.49d0
                w      = 1.04d0
            end if
    END SELECT
    fSWP = fSWP0(SWP,SWPmin,w)
    return
END !!function fSWP
!------------------------------------------------------------------
REAL(8) FUNCTION fSWP_Death(SWP,SWPmin,w)
    !!SWP Scalar for Microbial Mortality
    !!fSWP_Death DEcreases with increasing SWP (wetter condition)

    !!ARGUMENTS:
    REAL(8) SWP ![MPa]
    !!ARGUMENTS:
    REAL(8) SWPmin ![MPa]
    REAL(8) w
!    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
!    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    
    fSWP_Death = 1.d0 - fSWP0(SWP,SWPmin,w)
!    fSWP_Death = 1.d0 - fSWP(SWP,BIOME,SOM)
    
    return
END !!function fSWP_MicrobeMortality
!------------------------------------------------------------------
REAL(8) function fSWP_OPT(SWP)
    !SWP Scalar for SOM (lignin) decomposition
    !Hansen et al (1990), DAISY Model, page 105, Eq (6-16)
!    
!    REAL(8), PARAMETER :: SWPmin = -dexp(4.5*dlog(10d0))   !![MPa]
!    REAL(8), PARAMETER :: SWPlow = -dexp(0.5*dlog(10d0))
!    REAL(8), PARAMETER :: SWPhigh= -dexp(-0.5*dlog(10d0))
!    REAL(8), PARAMETER :: SWPmax = -dexp(-2.0*dlog(10d0))
!    !!ARGUMENTS:
!    REAL(8) SWP ![MPa]
!    
!    if (SWP.le.SWPmin) then
!        fSWP_OPT = 0.0d0  
!    else if (SWP.le.SWPlow) then
!        fSWP_OPT = 1.625-0.25*dlog10(1.0d2*dabs(SWP))
!    else if (SWP.le.SWPhigh) then
!        fSWP_OPT = 1.0d0
!    else if (SWP.le.SWPmax) then
!        fSWP_OPT = 0.6+0.4/1.5*dlog10(1.0d2*dabs(SWP))
!    else
!        fSWP_OPT = 0.6
!    end if
    
    REAL(8), PARAMETER :: SWPmin = -dexp(2.5*dlog(10d0))   !![MPa]
    REAL(8), PARAMETER :: SWPlow = -dexp(-1.5*dlog(10d0))
    REAL(8), PARAMETER :: SWPhigh= -dexp(-2.5*dlog(10d0))
    REAL(8), PARAMETER :: SWPmax = -dexp(-4.0*dlog(10d0))
    !!ARGUMENTS:
    REAL(8) SWP ![MPa]
    
    if (SWP.le.SWPmin) then
        fSWP_OPT = 0.0d0  
    else if (SWP.le.SWPlow) then
        fSWP_OPT = 0.625-0.25*dlog10(dabs(SWP))
    else if (SWP.le.SWPhigh) then
        fSWP_OPT = 1.0d0
    else if (SWP.le.SWPmax) then
        fSWP_OPT = (2.5+0.4*dlog10(dabs(SWP)))/1.5
    else
        fSWP_OPT = 0.6
    end if
    return
END !!function fSWP_OPT
!------------------------------------------------------------------
REAL(8) function fSWP_CLM(SWP,SWPsat)
    !Rate Scalar for Soil Water Potential (SWP)
    !CLM4.5 Technical Note, page 285, Eq (15.6)
    !Andren and Paustain (1987), Orchard and Cook (1983).
    
    REAL(8), PARAMETER :: SWPmin = -10 ![MPa], lower limit for SWP control on decomposition rate
    !!ARGUMENTS:
    REAL(8) SWP, SWPsat ![MPa], saturated SWP
    
    if (SWP.lt.SWPmin) then
        fSWP_CLM = 0.0d0
    else if (SWP.le.SWPsat) then
        fSWP_CLM = dlog(SWPmin/SWP)/dlog(SWPmin/SWPsat)
    else
        fSWP_CLM = 1.0d0
    end if
    return
END !!function fSWP_CLM
!------------------------------------------------------------------
REAL(8) function fWFP_PieceWise0(WFP,WFPcr,Slope,Intercept)
    !SWP Scalar for Nitrification/Denitrification
    !Muller (1999)
    
    !!ARGUMENTS:
    REAL(8) WFP ![-], water-filled porosity
    REAL(8) WFPcr(4), Slope(2), Intercept(2)
    
    if (WFP.le.WFPcr(1)) then
        fWFP_PieceWise0 = 0.0d0  
    else if (WFP.le.WFPcr(2)) then
        fWFP_PieceWise0 = Intercept(1) + Slope(1)*WFP
    else if (WFP.le.WFPcr(3)) then
        fWFP_PieceWise0 = 1.0d0
    else if (WFP.le.WFPcr(4)) then
        fWFP_PieceWise0 = Intercept(2) + Slope(2)*WFP
    else
        fWFP_PieceWise0 = 0.0d0  
    end if
    fWFP_PieceWise0 = max(0.d0, min(1.d0, fWFP_PieceWise0))
    return
END !!fWFP_PieceWise0
!------------------------------------------------------------------
REAL(8) function fWFP_PieceWise(sCase,WFP)
    !SWP Scalar for Nitrification/Denitrification
    !Muller (1999)
    
    !!ARGUMENTS:
    CHARACTER(len=*), intent(in) :: sCase
    REAL(8) WFP ![-], water-filled porosity
    
    !!LOCAL VARIABLES:
    REAL(8) WFPcr(4), Slope(2), Intercept(2)
    
    SELECT CASE (trim(sCase))
        CASE ("NITRIF")  !!
            WFPcr = (/0.09d0, 0.54d0,0.69d0,1.00d0/)
            Slope       = (/ 2.20d0, -3.23d0/)
            Intercept   = (/-0.19d0,  3.23d0 /)
        CASE ("DENITRIF_NO3")
            WFPcr       = (/0.36d0, 1.d0, 1.d0, 1.01d0 /)
            Slope       = (/ 1.56d0, 0.d0 /)
            Intercept   = (/-0.56d0, 1.d0/)
        CASE ("DENITRIF_NO2")
            WFPcr       = (/0.4d0, 0.6d0, 0.66d0, 0.7d0/)
            Slope       = (/5.d0, -20d0/)
            Intercept   = (/-2.d0, 14d0/)    
        CASE ("DENITRIF_NO")
            WFPcr       = (/ 0.1d0, 0.8d0, 0.9d0, 1.d0/)
            Slope       = (/1.43d0, -10d0/)
            Intercept   = (/-0.14d0, 10.d0/)
        CASE ("DENITRIF_N2O")
            WFPcr       = (/0.40d0, 0.85d0, 1.d0, 1.01d0/)
            Slope       = (/2.22d0, 0.d0/)
            Intercept   = (/-0.89d0,1.d0/)
        CASE DEFAULT
            
    END SELECT
    
    fWFP_PieceWise = fWFP_PieceWise0(WFP,WFPcr,Slope,Intercept)
    return
END !!fWFP_PieceWise
!------------------------------------------------------------------
REAL(8) function fSWPsat(p_sand,p_clay)
    !Saturated Soil Water Potential (SWP) [MPa]
    !CLM4.5 Technical Note, page 285, Eq (15.7)
    !Cosby et al. (1984)
    !p_sand (0-100): volume % of sand
    !p_clay (0-100): volume % of clay

    
    REAL(8), PARAMETER :: SWP_base = -9.8d-5 ![MPa], base SWP
    REAL(8), PARAMETER :: a = 1.54d0
    REAL(8), PARAMETER :: b = 9.5d-3
    REAL(8), PARAMETER :: c = 6.3d-3
    
    !!ARGUMENTS:
    REAL(8) p_sand, p_clay, p_silt
    
    p_silt = 100d0 - p_sand - p_clay
    fSWPsat = SWP_base*dexp((a-b*p_sand+c*p_silt)*dlog(10d0))
    return
END !!function fSWPsat

!------------------------------------------------------------------
REAL(8) function fSWP_A2D(SWP, SWP_A2D, w)
    !Soil Water Scalar for Microbial Dormancy
    !! Manzoni (2014) SBB, 73: 69-83

    !!ARGUMENTS:
    REAL(8), intent(in):: SWP_A2D !! = -0.4 [MPa], 
    REAL(8), intent(in):: SWP
    REAL(8), intent(in):: w       !! = 4 
    
    fSWP_A2D = dabs(SWP)**w/(dabs(SWP)**w + dabs(SWP_A2D)**w)
!    if(ISNAN(fSWP_A2D)) then
!        print*,"wp_scalar_low=",fSWP_A2D
!    end if
    return
END !!function fSWP_A2D

!------------------------------------------------------------------
REAL(8) function fSWP_D2A(SWP, SWP_D2A, w)
    !!Soil Water Scalar for Microbial reactivation
    !!ARGUMENTS:
    REAL(8), intent(in):: SWP_D2A !! = 1/4*SWP_A2D [MPa]
    REAL(8), intent(in):: SWP
    REAL(8), intent(in):: w       !! = 4 
    
    fSWP_D2A = dabs(SWP_D2A)**w/(dabs(SWP)**w + dabs(SWP_D2A)**w)
    return
END !!function fSWP_D2A

!------------------------------------------------------------------
REAL(8) function fO2_CONC(SWC,Porosity)
    !!fraction of O2 concentration in soil [m3 O2/cm3 air] 
    !!Davidson et al., 2012, GCB, 18: 371-384
    !!ARGUMENTS:
    REAL(8), intent(in):: SWC
    REAL(8), intent(in):: Porosity       
    
    !!LOCAL VARIABLES
    REAL(8), PARAMETER :: Dgas = 1.67D0 !!dimensionless diffusion coefficient
    REAL(8), PARAMETER :: frac_O2_in_air = 0.209D0
    REAL(8) AFP !!air-filled porosity
    
    AFP = max(0.D0,Porosity - SWC)
    fO2_CONC = Dgas*frac_O2_in_air*AFP**(4.D0/3.D0)
    return
END !!function fO2_CONC
!------------------------------------------------------------------
REAL(8) function fO2_scalar(SWC,Porosity)
    !!O2 saturation level in soil
    !!Davidson et al., 2012, GCB, 18: 371-384
    !!ARGUMENTS:
    REAL(8), intent(in):: SWC
    REAL(8), intent(in):: Porosity       
    
    !!LOCAL VARIABLES
    REAL(8) Ks_O2 !!half-saturation constant for O2
    REAL(8) O2_CONC !!O2 CONCENTRATION, [m3 O2/cm3 air]
!    
    Ks_O2   = fO2_CONC(0.5D0*Porosity,Porosity)
    O2_CONC = fO2_CONC(SWC,Porosity)
    
    fO2_scalar = O2_CONC/(O2_CONC + Ks_O2)
    return
END !!function fO2_scalar
!------------------------------------------------------------------
REAL(8) function fNLimit_MB(MB_CN,MB_CN_min,MB_CN_max,wexp)
    !!Nitrogen limitation level of microbes
    !!higher MB_CN, higher fNLimit_MB, higher N limitation 
    !!ARGUMENTS:
    REAL(8), intent(in):: MB_CN,MB_CN_min,MB_CN_max
    REAL(8), intent(in):: wexp !!exponential, >=1
    
    !!LOCAL VARIABLES
    REAL(8) CN
    if(MB_CN.gt.MB_CN_max) then
        CN = MB_CN_max
    else if(MB_CN.lt.MB_CN_min) then
        CN = MB_CN_min
    else
        CN = MB_CN
    end if
            
!    fNLimit_MB = ((MB_CN - MB_CN_min)/(MB_CN_max - MB_CN_min))**wexp
    fNLimit_MB = ((CN - MB_CN_min)/(MB_CN_max - MB_CN_min))**wexp
    
    return
END !!fNLimit_MB

!!-----------------------------------------------------------------
!! Soil Water Scalar: END
!!=================================================================


!!=================================================================
!! pH Scalar: BEGIN
!------------------------------------------------------------------
REAL(8) FUNCTION fpH(sCase,pH)
!    REAL(8) fpH0 !! function
    !!ARGUMENTS:
    CHARACTER(len=*) , intent(in) :: sCase
    REAL(8)          , intent(in) :: pH !pH value
    
    !!LOCAL VARIABLES
    REAL(8) pHopt !optimum pH
    REAL(8) pHsen !pH sensitivity
    
!    sCase = trim(sCase)
    SELECT CASE (trim(sCase))
        CASE ("BG")  !! Beta-glucosidase
            pHopt = 5.6
            pHsen = 1.7
        CASE ("CBH") !! Cellobiohydrolase
            pHopt = 5.1
            pHsen = 2.1
        CASE ("EG")  !! Endo-glucanase
            pHopt = 5.1
            pHsen = 1.6
        CASE ("PER") !! Peroxidase
            pHopt = 4.5
            pHsen = 1.5
        CASE ("POX") !! Phenol oxidase
            pHopt = 4.1
            pHsen = 1.4
        CASE ("LIG") !! Ligninases
            pHopt = 4.2
            pHsen = 1.4
        CASE ("CEL") !! Cellulases
            pHopt = 5.3
            pHsen = 1.7
        CASE ("PAC") !! Acid phosphatases
            pHopt = 5.2
            pHsen = 1.8
        CASE ("PAL") !! Alkaline phosphatases
            pHopt = 9.5
            pHsen = 2.6
        CASE ("PHO") !! PHOSPHATASES
            pHopt = 6.0
            pHsen = 2.0
        CASE ("ENZ") !! ENZ for MOM (Mineral-Associated Organic Matter)
            pHopt = 4.8
            pHsen = 1.6
        CASE DEFAULT
            pHopt = 6.0 !!mean pH of 763 soil samples
            pHsen = 2.0
    END SELECT
    
    fpH = fpH0(pH,pHopt,pHsen)
    return
END !!FUNCTION fpH
!------------------------------------------------------------------
REAL(8) FUNCTION fpH0(pH, pHopt, pHsen)
    !!ARGUMENTS:
    REAL(8), intent(in) :: pH !pH value
    REAL(8), intent(in) :: pHopt !optimum pH
    REAL(8), intent(in) :: pHsen !pH sensitivity
    fpH0 = dexp(-1d0 * ((pH - pHopt)/pHsen)**2d0)
    return
END !!FUNCTION fpH0
!!-----------------------------------------------------------------
!! pH Scalar: END
!!=================================================================

!!=================================================================
!! ISOTOPES: BEGIN
!!-----------------------------------------------------------------
REAL(8) FUNCTION fPermil(iOpt, Rstd, iso1, iso2)
    !Convert isotope concentration to signature/abundance [‰]
    !Rastetter et al. 2005. Ecological Applications 15, 1772-1782.
    !!ARGUMENTS:
    INTEGER, intent(in) :: iOpt !option, 0-input concentration of both iso1 and iso2, otherwise, iso2 = ratio of iso2/iso1
    REAL(8), intent(in) :: Rstd !standard ratio of iso2/iso1, e.g., C14/C12 = 1d-12, C13/C12 = 0.0112372
    REAL(8), intent(in) :: iso1 !concentration of iso1
    REAL(8), intent(in) :: iso2 !concentration of iso2, or iso2/iso1
    if (iOpt .eq. 0) then
        fPermil = ((iso2/iso1)/Rstd - 1d0) * 1000d0
    else
        fPermil = (iso2/Rstd - 1d0) * 1000d0
    end if
    return
END !!FUNCTION fPermil
!!-----------------------------------------------------------------
!! ISOTOPES: END
!!=================================================================

!!=================================================================
!! OUTPUT PROCESSING: BEGIN
!!-----------------------------------------------------------------
SUBROUTINE sOUT_OPT_h(nVAR,nHour,iHour,dSIM,sPAR,sOUT,VARopt_int)
    !!EXTRACT hourly output for opt RESPONSE VARIABLES
    !!ARGUMENTS:
    INTEGER, intent(in)     :: nVAR
    INTEGER, intent(in)     :: iHour,nHour
    INTEGER, intent(in)     :: VARopt_int(nVAR,3)
    REAL(8), intent(out)    :: dSIM(nHour,nVAR)
!    Real(8) CO2_ISO_inp  !!cumulative CO2 isotope (e.g., C14 or C13) at the beginning of the time-step
!    REAL(8), intent(in)     :: vENZ(2)      !!specific ENZ activity, ligninase & cellulase
    TYPE(sMEND_PAR),intent(in):: sPAR
    TYPE(sMEND_OUT),intent(in):: sOUT
    
    !!LOCAL VARIABLES:
    INTEGER j,k
    
    do j = 1,nVAR
        k = VARopt_int(j,1)
        select case (k)  !!see "MEND.ini"
            case (1)    !!CO2
                dSIM(iHour,j) = sOUT%MNFLUX%CO2_gmo
            case (2)    !!CO2_ISO,e.g., C14_CO2
                dSIM(iHour,j) = sOUT%CO2_gm_iso
!                dSIM(iHour,j) = sOUT%CPOOLI(2)%CO2 - CO2_ISO_inp!!sINP%CPOOLI(2)%CO2
            case (3)    !!MB
                dSIM(iHour,j) = sOUT%CPOOL%MB
!            case (4)    !!MB_ISO
!                dSIM(iHour,j) = sOUT%CPOOLI(2)%MB
            case (5)    !!DOM
                dSIM(iHour,j) = sOUT%CPOOL%MBA/sOUT%CPOOL%MB
            case (6)    !!SOM
                dSIM(iHour,j) = sOUT%CPOOL%SOM
!            case (7)    !!SOM_ISO
!                dSIM(iHour,j) = sOUT%CPOOLI(2)%SOM
            case (8)    !!ENZ_MOM, ENZ Activity = ENZ concentration*specific enzyme activity
                dSIM(iHour,j) = sOUT%CPOOL%ENZM*sPAR%VdMOM
            case (9)    !!ENZ_LIG, ENZ Activity = ENZ concentration*specific enzyme activity
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(1)*sPAR%VdPOM(1)
            case (10)   !!ENZ_CEL, ENZ Activity
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(2)*sPAR%VdPOM(2)
            case (11)    !!SOM_CN
                dSIM(iHour,j) = sOUT%CN%SOM
            case (12)    !!dPOC(1)
                dSIM(iHour,j) = sOUT%CFLUX%POMdec(1)
            case (13)    !!dPOC(2)
                dSIM(iHour,j) = sOUT%CFLUX%POMdec(2)
            case (14)    !!dMOC
                dSIM(iHour,j) = sOUT%CFLUX%MOMdec
            case (15)    !!MBC_grow
                dSIM(iHour,j) = sOUT%MNFLUX%CO2_growth
            case (16)    !! Active fraction
                dSIM(iHour,j) = sOUT%CPOOL%MBA/sOUT%CPOOL%MB
            case (17)    !!ENZC_LIG, ENZ concentration
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(1)
            case (18)   !!ENZC_CEL, 
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(2)
            case (19)   !!ENZC_MOM, 
                dSIM(iHour,j) = sOUT%CPOOL%ENZM
            case (20)   !!ENZC_SOM, 
                dSIM(iHour,j) = sOUT%CPOOL%ENZ
            case default
                dSIM(iHour,j) = sOUT%MNFLUX%CO2_gmo
        end select
    end do
    
END !!subroutine sOUT_OPT_h
!!-----------------------------------------------------------------
SUBROUTINE sOUT_OPT(nh,dSIM_h,nt,dSIM_t,sDate_beg,sDate_end,tstep,flag_avg)
    !!convert hourly data to any time step
    !!dSIM_h(nh):hourly data
    !!sDate_beg: YYYYMMDD, beginning date
    !!tstep = 1(daily),2(monthly),3(seasonal),4(yearly)
    !!flag_avg = 1(average), 0(hourly data at the end of tstep)
    !!ARGUMENTS:
    INTEGER nh, nt, tstep, flag_avg
    REAL(8) dSIM_h(nh)
    REAL(8) dSIM_t(nt,2)  !!mean & sd
    CHARACTER(len=8) sDate_beg,sDate_end
    
    !!LOCAL VARIABLES:
    INTEGER iyr0, imo0, ida0, iyr, imo
    INTEGER nda_in_ym
    INTEGER i,j, nhr,nda,nmo, ibeg, iend,nbe
    
    nda = nDaysbwDates(sDate_beg,sDate_end)
    nmo = nMonsbwDates(sDate_beg,sDate_end)
    nhr = nda*24
    if(nh.ne.nhr) then
        write(*,*)"Error: nRow<>nHour in 'sOUT_OPT'"
    end if
    
    CALL sDate2YMD(sDate_beg,iyr0,imo0,ida0)
    
    select case (tstep)
        case (1) !!daily
            do i = 1,nda
                ibeg=(i-1)*24+1
                iend=i*24
                nbe = iend - ibeg + 1
                if(flag_avg.eq.1) then
                    dSIM_t(i,1) = fAVG2(nh,dSIM_h,ibeg,iend,const_FillValue)
!                    write(*,*)i,dSIM_t(i)
                else
                    dSIM_t(i,1) = dSIM_h(iend)
                end if
                dSIM_t(i,2) = fSTDDEV(nbe,dSIM_h(ibeg:iend),const_FillValue)
            end do
        case (2) !!monthly
            iyr = iyr0
            imo = imo0
            nda = 0
            do j=1,nmo
!                write(*,*)'iyr,imo = ',iyr,imo
                nda_in_ym = nDaysofMon(iyr,imo)
                ibeg = nda*24+1
                iend = (nda + nda_in_ym)*24
                nbe = iend-ibeg+1
                nda = nda + nda_in_ym
                if(flag_avg.eq.1) then
                    dSIM_t(j,1) = fAVG2(nh,dSIM_h,ibeg,iend,const_FillValue)
                else
                    dSIM_t(j,1) = dSIM_h(iend)
                end if
                dSIM_t(j,2) = fSTDDEV(nbe,dSIM_h(ibeg:iend),const_FillValue)
                imo = imo + 1
                if(imo.le.12) then
                    iyr = iyr
                else
                    iyr = iyr+1
                    imo = 1
                end if
            end do
!        case(3) !!seasonal
!            
!        case(4) !!yearly
            
        case default
        
    end select
END !!subroutine sOUT_OPT

!!-----------------------------------------------------------------
SUBROUTINE sOUT_Day2Mon(nday,dSIM_d,nt,dSIM_t,sDate_beg,sDate_end,tstep,flag_avg)
    !!convert daily data to any time step
    !!dSIM_d(nday):daily data
    !!sDate_beg: YYYYMMDD, beginning date
    !!tstep = 1(daily),2(monthly),3(seasonal),4(yearly)
    !!flag_avg = 1(average), 0(hourly data at the end of tstep)
    !!ARGUMENTS:
    INTEGER nday, nt, tstep, flag_avg
    REAL(8) dSIM_d(nday)
    REAL(8) dSIM_t(nt,2)  !!mean & sd
    CHARACTER(len=8) sDate_beg,sDate_end
    
    !!LOCAL VARIABLES:
    INTEGER iyr0, imo0, ida0, iyr, imo
    INTEGER nda_in_ym
    INTEGER i,j,nda,nmo, ibeg, iend,nbe
    
    nda = nDaysbwDates(sDate_beg,sDate_end)
    nmo = nMonsbwDates(sDate_beg,sDate_end)
!    nhr = nda*24
    if(nday.ne.nda) then
        write(*,*)"Error: nRow<>nday in 'sOUT_Day2Mon'"
    end if
    
    CALL sDate2YMD(sDate_beg,iyr0,imo0,ida0)
    
    select case (tstep)
!        case (1) !!daily
!            do i = 1,nda
!                ibeg=(i-1)*24+1
!                iend=i*24
!                nbe = iend - ibeg + 1
!                if(flag_avg.eq.1) then
!                    dSIM_t(i,1) = fAVG2(nh,dSIM_h,ibeg,iend,const_FillValue)
!!                    write(*,*)i,dSIM_t(i)
!                else
!                    dSIM_t(i,1) = dSIM_h(iend)
!                end if
!                dSIM_t(i,2) = fSTDDEV(nbe,dSIM_h(ibeg:iend),const_FillValue)
!            end do
        case (2) !!monthly
            iyr = iyr0
            imo = imo0
            nda = 0
            do j=1,nmo
!                write(*,*)'iyr,imo = ',iyr,imo
                nda_in_ym = nDaysofMon(iyr,imo)
                ibeg = nda+1
                iend = nda + nda_in_ym
                nbe = nda_in_ym
                nda = nda + nda_in_ym
                if(flag_avg.eq.1) then
                    dSIM_t(j,1) = fAVG2(nday,dSIM_d,ibeg,iend,const_FillValue)
                else
                    dSIM_t(j,1) = dSIM_d(iend)
                end if
                dSIM_t(j,2) = fSTDDEV(nbe,dSIM_d(ibeg:iend),const_FillValue)
                imo = imo + 1
                if(imo.le.12) then
                    iyr = iyr
                else
                    iyr = iyr+1
                    imo = 1
                end if
            end do
!        case(3) !!seasonal
!            
!        case(4) !!yearly
            
        case default
        
    end select
END !!subroutine sOUT_OPT
!!-----------------------------------------------------------------
SUBROUTINE sOUT_ALL_tscale(sFile_hour,sFile_t,nRow_skip,nVAR,sDate_beg,sDate_end,tstep,flag_avg)
    !!convert hourly data to any time step
    !!sFile_hour: file name with hourly data
    !!nRow_skip: >=2, first nRow to be skipped: 1st-row: time period; 2nd-row: head; >=3rd-row: others, e.g., values at t=0
    !!nday: # of days
    !!nVAR: # of variables (columns) in datafile
    !!nh: # of hours
    !!dSIM_h(nh):hourly data
    !!sDate_beg: YYYYMMDD, beginning date
    !!tstep = 1(daily),2(monthly),3(seasonal),4(yearly)
    !!flag_avg = 1(average), 0(hourly data at the end of tstep)
    !!ARGUMENTS:
    INTEGER nRow_skip,nVAR,tstep, flag_avg !!nh, nt,
    CHARACTER(len=*)sFile_hour,sFile_t
    CHARACTER(len=8) sDate_beg,sDate_end
!    REAL(8) dSIM_h(nh)
!    REAL(8) dSIM_t(nt,2)  !!mean & sd
    
    !!LOCAL VARIABLES:
    INTEGER iyr0, imo0, ida0, iyr, imo, iyear,imonth,iday
    INTEGER nda_in_ym
    INTEGER i,j,k 
    INTEGER nhr,nda,nmo,nyr, ibeg, iend,nbe
    INTEGER iRead
    CHARACTER(LEN=10) sDateHr
    CHARACTER(len=8) sDate
    CHARACTER(len=6) sYM
    CHARACTER(len=2000)sRead1,sRead2
    CHARACTER(len=200)format1
    REAL(8), DIMENSION(24,nVAR)::dSIM_h1
    REAL(8), DIMENSION(31,nVAR)::dSIM_h2 !!days in 1 month, at most = 31, hours in 1 month, at most =24*31
    REAL(8), DIMENSION(366,nVAR)::dSIM_h3 !!days in 1 year, at most =366, hours in 1 year, at most =366*24
    REAL(8), DIMENSION(nVAR)   ::dSIM_t
    dSIM_h1 = const_FillValue
    dSIM_h2 = const_FillValue
    dSIM_h3 = const_FillValue
    dSIM_t = const_FillValue
    
!    write(format1,*)"(i10,",nVAR,"(e20.6))"
    if(tstep.eq.4) then !!yearly
        write(format1,*)"(I10,",nVAR,"(e20.6))"
    else
        write(format1,*)"(A10,",nVAR,"(e20.6))"
    end if
    
    open(2,file=sFile_t,status='unknown')
    
    nda = nDaysbwDates(sDate_beg,sDate_end)
    nmo = nMonsbwDates(sDate_beg,sDate_end)
    nyr = nYearsbwDates(sDate_beg,sDate_end)
!    nhr = nda*24
!    if(nh.ne.nhr) then
!        write(*,*)"Error: nRow<>nHour in 'sOUT_OPT'"
!    end if
!    
    CALL sDate2YMD(sDate_beg,iyr0,imo0,ida0)
    open(1,file=sFile_hour,status='old')
    read(1,'(a)')sRead1      !1st line: time period
    write(2,'(a)')sRead1
    read(1,'(a10,a)')sRead1,sRead2  !2nd line: head
    if(tstep.eq.1) then
        write(2,'(a10,a)')"Day",sRead2
    elseif(tstep.eq.2) then
        write(2,'(a10,a)')"Mon",sRead2
    elseif(tstep.eq.4) then
        write(2,'(a10,a)')"Year",sRead2
    end if
    
    do i=1,nRow_skip-2
        read(1,'(a)')sRead1
        write(2,'(a)')sRead1
    end do
    
    select case (tstep)
        case (1) !!convert hourly to daily
            if(nda > 0) then
                do i = 1,nda
                    ibeg = 1
                    iend = 24
                    nbe = iend - ibeg + 1
                    do j=1,24 !!24 hours in 1 day
                        read(1,*)sDateHr,dSIM_h1(j,1:nVAR)
                    end do
                    do j=1,nVAR
                        if(flag_avg.eq.1) then
                            dSIM_t(j) = fAVG2(iend,dSIM_h1(:,j),ibeg,iend,const_FillValue)
        !                    write(*,*)i,dSIM_t(i)
                        else
                            dSIM_t(j) = dSIM_h1(iend,j)
                        end if
                    end do
                    CALL sDate_After(i,sDate_beg,sDate)
                    write(2,format1)sDate,dSIM_t
                end do
            end if !! nda > 0
        case (2) !!convert daily to monthly
            if(nmo > 0) then
                iyr = iyr0
                imo = imo0
                nda = 0
                do i=1,nmo
                    nda_in_ym = nDaysofMon(iyr,imo)
    !                write(*,*)'iyr,imo,nday = ',iyr,imo,nda_in_ym
                    ibeg = 1
!                    iend = nda_in_ym*24
                    iend = nda_in_ym
                    do j = ibeg,iend
                        read(1,*)sDate,dSIM_h2(j,1:nVAR)
                    end do
                    do j=1,nVAR
                        if(flag_avg.eq.1) then
                            dSIM_t(j) = fAVG2(iend,dSIM_h2(:,j),ibeg,iend,const_FillValue)
                        else
                            dSIM_t(j) = dSIM_h2(iend,j)
                        end if
                    end do
                    imo = imo + 1
                    if(imo.le.12) then
                        iyr = iyr
                    else
                        iyr = iyr+1
                        imo = 1
                    end if
                    CALL sYM_After(i,sDate_beg(1:6),sYM)
                    write(2,format1)sYM,dSIM_t
                end do
            end if !!(nmo>0)
!        case(3) !!seasonal
!            
        case(4) !!convert daily to yearly
            if(nyr>0) then
                do i=1,nyr
                    iyr = iyr0 + i -1
                    ibeg = 1
                    if(i.eq.1) then
                        call sYMD2Date(iyr, 12, 31, sDate)
                        iend = nDaysbwDates(sDate_beg, sDate)
                    else if(i.eq.nyr) then
                        call sYMD2Date(iyr, 1, 1, sDate)
                        iend = nDaysbwDates(sDate, sDate_end)
                    else !! 1<i<nyr
                        iend = nDaysofYear(iyr)
                    end if
                    
                    do j = ibeg,iend
                        read(1,*)sDate,dSIM_h3(j,1:nVAR)
                    end do
                    
                    do j=1,nVAR
                        if(flag_avg.eq.1) then
                            dSIM_t(j) = fAVG2(iend,dSIM_h3(:,j),ibeg,iend,const_FillValue)
                        else
                            dSIM_t(j) = dSIM_h3(iend,j)
                        end if
                    end do
                    write(2,format1)iyr,dSIM_t
                end do !!do i=1,nyr
            end if !!if(nyr > 0)
            
        case default
        
    end select
    close(1)
    close(2)
END !!subroutine sOUT_OPT
!!-----------------------------------------------------------------
SUBROUTINE sOUT_tscale(dirout,sDate_beg,sDate_end)
!!Convert OUTPUTS from HOURLY to DAILY & MONTHLY for all state variables & fluxes
    !!ARGUMENTS:
    CHARACTER(LEN = *)dirout
    CHARACTER(LEN = 8)sDate_beg,sDate_end
    
    !!LOCAL VARIABLES:
    CHARACTER(LEN = 200) sFile_inp, sFile_out
    INTEGER nRow_skip, nVAR, tstep,flag_avg
    
    print*,""
    print*,">>>Convert OUTPUTs from Hourly to Daily-Monthly-Yearly: BEG"
    
    print*,">>>[1] STATE VARIABLEs:"
    sFile_inp = trim(dirout)//"VAR_hour.out"   
    sFile_out = trim(dirout)//"VAR_day.out"
!    nRow_skip=3; nVAR=const_nPOOL*(const_nISO+1); tstep=1; flag_avg=1
    nRow_skip=3; nVAR=const_nPOOL*3 + const_nPOOL_MN; tstep=1; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    CALL system('gzip -f '//sFile_inp)
    
    sFile_inp = trim(dirout)//"VAR_day.out" 
    sFile_out = trim(dirout)//"VAR_mon.out"
!    nRow_skip=3; nVAR=const_nVARc*(const_nISO+1); tstep=2; flag_avg=1
    nRow_skip=3; nVAR=const_nPOOL*3 + const_nPOOL_MN; tstep=2; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_inp = trim(dirout)//"VAR_day.out"   
    sFile_out = trim(dirout)//"VAR_year.out"
    nRow_skip=3; nVAR=const_nPOOL*3 + const_nPOOL_MN; tstep=4; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    print*,">>>[2] FLUXes:"
    sFile_inp = trim(dirout)//"FLX_hour.out"
    sFile_out = trim(dirout)//"FLX_day.out"
    nRow_skip=2; nVAR=const_nFLUX*2 + const_nFLUX_MN; tstep=1; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    CALL system('gzip -f '//sFile_inp)
    
    sFile_inp = trim(dirout)//"FLX_day.out"   
    sFile_out = trim(dirout)//"FLX_mon.out"
    nRow_skip=2; nVAR=const_nFLUX*2 + const_nFLUX_MN; tstep=2; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_inp = trim(dirout)//"FLX_day.out"   
    sFile_out = trim(dirout)//"FLX_year.out"
    nRow_skip=2; nVAR=const_nFLUX*2 + const_nFLUX_MN; tstep=4; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    print*,">>>[3] PARAMETERs:"
    sFile_inp = trim(dirout)//"PAR_hour.out"
    sFile_out = trim(dirout)//"PAR_day.out"
    nRow_skip=2; nVAR=const_nPAR; tstep=1; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    CALL system('gzip -f '//sFile_inp)
    
    sFile_inp = trim(dirout)//"PAR_day.out"
    sFile_out = trim(dirout)//"PAR_mon.out"
    nRow_skip=2; nVAR=const_nPAR; tstep=2; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_inp = trim(dirout)//"PAR_day.out"
    sFile_out = trim(dirout)//"PAR_year.out"
    nRow_skip=2; nVAR=const_nPAR; tstep=4; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    print*,">>>[4] DERIVED RATEs:"
    sFile_inp = trim(dirout)//"RATE_hour.out"
    sFile_out = trim(dirout)//"RATE_day.out"
    nRow_skip=2; nVAR=const_nRATE; tstep=1; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    CALL system('gzip -f '//sFile_inp)
    
    sFile_inp = trim(dirout)//"RATE_day.out"
    sFile_out = trim(dirout)//"RATE_mon.out"
    nRow_skip=2; nVAR=const_nRATE; tstep=2; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_inp = trim(dirout)//"RATE_day.out"
    sFile_out = trim(dirout)//"RATE_year.out"
    nRow_skip=2; nVAR=const_nRATE; tstep=4; flag_avg=1
    CALL sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    print*,">>>Convert OUTPUTs from Hourly to Daily-Monthly-Yearly: END"

END !!subroutine sOUT_tscale
!!-----------------------------------------------------------------


SUBROUTINE sINP_Read(nfile,sfilename,dirinp,ststep,is_total,nMon,nHour,rINP)
!!Read inputs: STP, SWC, SIN
    !!ARGUMENTS:
    INTEGER,         intent(in)     :: nfile, nHour, nMon
    INTEGER,         intent(in)     :: is_total  !!=1: need to convert to hourly rate; =0: directly assign value
    CHARACTER(LEN=*),intent(in)     :: sfilename(nfile)
    CHARACTER(LEN=*),intent(in)     :: dirinp, ststep
    REAL(8),         intent(inout)  :: rINP(nHour)
    
    !!LOCAL VARIABLES:
    INTEGER i,j,k,lp
    INTEGER ndays, nmons, iyr,imo, eof
    REAL(8) rRead
    CHARACTER(len=50)sRead
    CHARACTER(len=200) sfilename_full
    CHARACTER(len=8) sDate_beg,sDate_end
    
    k = 0 !!hour
      if(trim(ststep).eq."hourly") then             
          do i = 1,nfile
              sfilename_full = trim(dirinp)//trim(sfilename(i))
              open(101,file=sfilename_full,status='old')
              read(101,'(a)')sDate_beg
              read(101,'(a)')sDate_end
              ndays = nDaysbwDates(sDate_beg,sDate_end)
              do j=1,ndays*24
                  k = k + 1
                  read(101,*,iostat=eof)rRead
                  if(eof<0) exit
                  rINP(k) = rRead
              end do
              close(101)
          end do
      elseif(trim(ststep).eq."daily") then
          do i = 1,nfile
              sfilename_full = trim(dirinp)//trim(sfilename(i))
              open(101,file=sfilename_full,status='old')
              read(101,'(a)')sDate_beg
              read(101,'(a)')sDate_end
              ndays = nDaysbwDates(sDate_beg,sDate_end)
              do j=1,ndays
                  read(101,*,iostat=eof)rRead
                  if(eof<0) exit
                  do lp = 1,24 !!24 hours in 1 day
                    k = k + 1
                    if(is_total.eq.1) then
                      rINP(k) = rRead/DBLE(24)
                    else
                      rINP(k) = rRead
                    end if
                  end do
              end do
              close(101)
          end do
      elseif(trim(ststep).eq."monthly") then !!only 1 file is allowed
          sfilename_full = trim(dirinp)//trim(sfilename(1))
          open(101,file=sfilename_full,status='old')
          read(101,'(a)')sRead !!head
          nmons = nMon !!nMonsbwDates(sINI%sDate_beg_all,sINI%sDate_end_all)
          do j=1,nmons
              read(101,*,iostat=eof)iyr,imo,rRead  !!SWC or SWP
              if(eof<0) exit
              ndays = nDaysofMon(iyr,imo)
              do lp = 1,ndays*24
                k = k+1
                if(is_total.eq.1) then
                    rINP(k) = rRead/DBLE(ndays*24)
                else
                    rINP(k) = rRead
                end if
              end do
          end do   
          close(101)
      end if
    
END !!SUBROUTINE sINP_Read
!!-----------------------------------------------------------------
!! OUTPUT PROCESSING: END
!!=================================================================

!!---------------------------------------------------------    
    SUBROUTINE TEST()
!! TEST: BEGIN 
!        USE MOD_MEND, ONLY: fSWC2SWP,fSWP2SWC
!    REAL(8), dimension(8) :: v = (/ 2,4,4,4,5,5,7,9 /)
!    REAL(8) sd
!    sd = fSTDDEV(6,v(2:7),const_FillValue)
!    print *, "std dev = ", sd

!    REAL(8) fSWC2SWP
    INTEGER i,j, ni, nj,iRead(10)
    REAL(8) rRead(10)
    CHARACTER(len=20) sRead(10),units
    character(len=200):: sFile_inp,sFile_out
!    REAL(8) fSWPsat1,fSWP_dec1
    REAL(8) SWCres,SWCsat,alpha,rn
    REAL(8) SWC,SWP
!    fSWPsat1 = fSWPsat(20d0,30d0)
!    fSWP_dec1 = fSWP_dec(-1d0,fSWPsat1)
    
!    print*, fSWP_OPT(-dexp(3.5*dlog(10d0)))
!    print*, fSWP_OPT(-dexp(0.0*dlog(10d0)))
!    print*, fSWP_OPT(-dexp(-0.5*dlog(10d0)))
!    print*, fSWP_OPT(-dexp(-1*dlog(10d0)))
    
!    SWCres = 0.162603606
!    SWCsat = 0.579564599
!    alpha  = 0.021
!    rn     = 1.5630346
    
    SWCres = 0.083
    SWCsat = 0.46
    alpha  = 0.025
    rn     = 1.31
    
    units = "perc"
    sFile_inp = "/Users/wg4/Dropbox (ORNL)/Model/MEND/userio/inp/SWC_mon_Harvard.dat"
    sFile_out = "/Users/wg4/Dropbox (ORNL)/Model/MEND/userio/inp/SWP_mon_Harvard.dat"
    ni = 264
    nj = 3
    open(unit=1,file=sFile_inp,status='unknown')
    open(unit=2,file=sFile_out,status='unknown')
    read(1,*)sRead(1:6)
    write(2,*)sRead(1:6)
    do i=1,ni
        read(1,*)sRead(1:3),rRead(1:3)
        do j=1,nj
            SWC = rRead(j)
            SWP = fSWC2SWP(SWC,SWCres,SWCsat,alpha,rn,-1d3)
            rRead(nj+j) = SWP
        end do
        write(2,*)sRead(1:3),rRead((nj+1):(nj*2))
    end do
    
    close(1)
    close(2)
    
!    write(*,*)SWC,SWP
!    units = "MPa"
!    SWC = fSWP2SWC(SWP,units,SWCres,SWCsat,alpha,rn)
!    write(*,*)SWC,SWP

!    open(unit=1,file = './userio/inp/SWC2009.txt',status='old')
!    open(unit=2,file = './userio/inp/SWP2009.dat',status='unknown')
!    !!head
!    read(1,*) sRead(1:5)
!    write(*,*)sRead(1:5)
!    do i = 1, 3672
!        read(1,*)sRead(1:4),SWC
!        SWP = fSWC2SWP(SWC,units,SWCres,SWCsat,alpha,rn,-1d3)
!!        write(*,*)SWC,SWP
!        write(2,*),sRead(1:4),SWC,SWP
!    end do
!    
!    close(1)
!    close(2)
!! TEST: END
!!---------------------------------------------------------
    END

!------------------------------------------------------------------

!REAL(8) function fMEND_OBJ0(xx, sPAR, sINI, sOUT)
!    USE STRUCT_MEND
!    TYPE(sMEND_PAR), intent(inout) :: sPAR
!    TYPE(sMEND_INI), intent(inOut) :: sINI
!    TYPE(sMEND_OUT) sOUT
!    TYPE(sMEND_INP) sINP
!    REAL(8) sum1, fNSE !function
!    REAL(8) xx(sPAR % nPar)
!    INTEGER nObj
!
!    sPAR % VdPOM = (/xx(1), xx(2)/) ![mg POM/mg ENZP/h], maximum reaction rate for conversion of POM by ENZP
!    sPAR % KsPOM = (/xx(3), xx(4)/) ![mg POM/cm3], half-saturation constant for conversion of POM by ENZP
!    sPAR % rENZP = (/xx(5), xx(5)/) ![1/h],turnover rate of ENZP
!    sPAR % pENZP = xx(6) ![mg ENZP/mg MB/h], production rate of ENZP
!    sPAR % frPOM2DOM = xx(7) ![-], fraction of decomposed POM allocated to DOM
!    sPAR % frMB2DOM = xx(8) ![-], fraction of dead MB allocated to DOM
!    sPAR % Vg = xx(9) ![mg DOM/mg MB/h], maximum uptake rate of DOM by MB
!    sPAR % KsDOM = xx(10) ![mg DOM/cm3], half-saturation constant for uptake of DOM by MB
!    sPAR % Yg = xx(11) ![-], carbon use efficiency in uptake of DOM by MB
!    sPAR % Vm = xx(12) ![1/h], specific microbial maintenance rate
!    sPAR % rENZM = xx(13) ![1/h], turnover rate of ENZMAOC
!    sPAR % pENZM = xx(14) ![mg ENZM/mg MB/h], production rate of ENZMAOC
!    sPAR % VdMOM = xx(15) ![mg MOM/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC
!    sPAR % KsMOM = xx(16) ![mg MOM/cm3], half-saturation constant for conversion of MAOC by ENZMAOC
!    sPAR % Qmax = xx(17) ![mg C/g soil], adsorption capacity
!    sPAR % Kba = xx(18) ![mg C/g soil/h], binding affinity
!    sPAR % Kdes = xx(19) ![mg DOM/h],desorption rate constant
!    sPAR % Kads = sPAR % Kdes * sPAR % Kba ![mg DOM/mg DOM/h], adsorption rate constant = Kba*Kdes
!
!
!    !     print*, sINI%nHour
!    CALL subMEND_RUN(xx, sPAR, sINI, sOUT)
!
!    nObj = 0
!    sum1 = 0d0
!    sum2 = 0d0
!    do i = 1, sINI % nObs_var
!        sINI % rNSE(i) = fNSE(sINI % nObs_time, sINI % dObs_comp(i,:), sINI % dSim_comp(i,:), const_FillValue)
!        !         print*, i, sINI%rNSE(i)
!        if (sINI % rNSE(i) .ne. const_FillValue) then
!            nObj = nObj + 1
!            sum1 = sum1 + sINI % rNSE(i) * sINI % rNSE_weight(i)
!            sum2 = sum2 + sINI % rNSE_weight(i)
!        end if
!    end do
!
!    fMEND_OBJ = 1.0 - sum1/sum2 !sINI%nObs_var  !Minimization
!
!end function fMEND_OBJ0

!end module TMEND           
END MODULE MOD_MEND

