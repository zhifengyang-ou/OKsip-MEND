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
    INTEGER, PARAMETER :: nProp = 13 !# of soil properties
    REAL(8) dINI(nProp)
    INTEGER nObs_time, nObs_var, nObs_col != nSoil*nSubstrate*nObs_var  !# of columns in observation data file

    INTEGER i, j, k, lp
    INTEGER ifini, ifobs  
    REAL(8) frISO(const_nISO), frISOadd(const_nISO)
    REAL(8) rLig_Cel(2) !proportion of Lignin and Cellulose in POC
    REAL(8) rQOC !proportion of QOC in total MOC

    CHARACTER(len = 100) sRead!!propName(nProp) !name of 10 properties 
    INTEGER iRead  !!iTemp
    REAL(8) rRead  !!rTemp
    
    rLig_Cel = (/sINI%LCI0, 1.0 - sINI%LCI0/) !proportion of Lignin and Cellulose in POC
    rQOC = 0.01d0 !proportion of QOC in total MOC

    ifini = 1
    ifobs = 2
    sRead = trim(sINI%dirinp_case)//trim(sINI%SOIL_INI_file)
    open(unit = ifini, file = sRead, status = 'old')
    do i = 1, 3
        read(ifini, *) sRead !head: !!ID	Property	 Value
    end do
    do j = 1, nProp
        read(ifini, *) iRead, sRead, dINI(j)
    end do
    close(ifini)
!
!    !!READ observations for model calibration
    lp = 0
    do i=1,sINI%nVARopt
        j = sINI%VARopt_int(i,1)
        sRead = trim(sINI%dirinp_case)//trim(sINI%VARfile(j))
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
    frISO(2) = 0 !const_Rstd(1)/(1d0 + const_Rstd(1)) !standard ratio of C14/C12 - 1e-12
    frISO(1) = 1d0 - frISO(2) !e.g., C12 

!
!! sINP%CPOOL        !carbon pools
    sINI % soilDepth = dINI(13)  !![cm], soil depth
    sINP % CPOOL % POC = dINI(3) * rLig_Cel ![mg C/g soil],Particulate Organic Carbon, size of {} determined by nPOC, !(/DBLE(0.958), DBLE(2.290)/)  !reshape((/DBLE(0.958), DBLE(2.290)/),shape(sINP%CPOOL%POC))
    sINP % CPOOL % MOC = dINI(4)*(1.d0 - rQOC) ![mg C/g soil],Mineral Associate Organic Carbon  !DBLE(27.868) 
    sINP % CPOOL % QOC = dINI(4) * rQOC ![mg C/g soil],adsorbed phase of DOC  !DBLE(0.100)
    sINP % CPOOL % MBC = dINI(5) ![mg C/g soil],Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % DOC = dINI(6) ![mg C/g soil],Dissolved Organic Carbon     !DBLE(0.210)  
    sINP % CPOOL % MBCA = sINP % CPOOL % MBC * sINI % r0 ![mg C/g soil],Active Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % MBCD = sINP % CPOOL % MBC * (1.d0 - sINI % r0) ![mg C/g soil],Active Microbial Biomass Carbon  !DBLE(0.640) 
    sINP % CPOOL % ENZP = (/dINI(7),dINI(8)/)!!(/1.1d-3, 1.1d-3/) ![mg C/g soil],ENZyme for POC        !reshape((/DBLE(1.1d-3),DBLE(1.1d-3)/),shape(sINP%CPOOL%ENZP))   
    sINP % CPOOL % ENZM = dINI(9) !!1.4d-3 ![mg C/g soil],Enzyme for MAOC  
    sINP % CPOOL % CO2 = 0d0 ![mg C/g soil],CO2 in the soil
    sINP % CPOOL % SOC = sum(sINP % CPOOL % POC) + sINP % CPOOL % MOC + sINP % CPOOL % QOC
!!     END TYPE sMEND_INP
    do j = 1, const_nPOC
        sINP % CPOOLIFR % POC(j) = frISO
        sINP % CPOOLIFR % ENZP(j) = frISO

        sINP % CPOOLI % POC(j) = sINP % CPOOLIFR % POC(j) * sINP % CPOOL % POC(j)
        sINP % CPOOLI % ENZP(j) = sINP % CPOOLIFR % ENZP(j) * sINP % CPOOL % ENZP(j)

!!        sINP % CADDI % POCadd(j) = sINP % CADD % POCadd(j) * frISOadd
    end do

!!    sINP % CADDI % DOCadd = sINP % CADD % DOCadd * frISOadd

    sINP % CPOOLIFR % MOC = frISO
    sINP % CPOOLIFR % ENZM = frISO
    sINP % CPOOLIFR % QOC = frISO
    sINP % CPOOLIFR % DOC = frISO
    sINP % CPOOLIFR % MBC = frISO
    sINP % CPOOLIFR % MBCA = frISO
    sINP % CPOOLIFR % MBCD = frISO
    sINP % CPOOLIFR % CO2 = frISO

    !     print*, "C12 added = ",sINP%CADDI(1)
    !     print*, "C14 added = ",sINP%CADDI(2) 

    sINP % CPOOLI % MOC = sINP % CPOOLIFR % MOC * sINP % CPOOL % MOC
    sINP % CPOOLI % ENZM = sINP % CPOOLIFR % ENZM * sINP % CPOOL % ENZM
    sINP % CPOOLI % QOC = sINP % CPOOLIFR % QOC * sINP % CPOOL % QOC
    sINP % CPOOLI % DOC = sINP % CPOOLIFR % DOC * sINP % CPOOL % DOC
    sINP % CPOOLI % MBC = sINP % CPOOLIFR % MBC * sINP % CPOOL % MBC
    sINP % CPOOLI % MBCA = sINP % CPOOLIFR % MBCA * sINP % CPOOL % MBCA
    sINP % CPOOLI % MBCD = sINP % CPOOLIFR % MBCD * sINP % CPOOL % MBCD
    sINP % CPOOLI % CO2 = sINP % CPOOLIFR % CO2 * sINP % CPOOL % CO2


    do j = 1, const_nISO
        sINP % CPOOLI(j) % SOC  = sum(sINP % CPOOLI(j) % POC) &
        &                       + sINP % CPOOLI(j) % MOC + sINP % CPOOLI(j) % QOC
    end do
    sINP % CPOOLIFR % SOC = sINP % CPOOLI % SOC/sINP % CPOOL % SOC

    do j = 1, const_nISO - 1
        do k = 1, const_nPOC
            sINP % CPOOLI_SIG(j) % POC(k) = &
            & fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % POC(k), sINP % CPOOLI(j + 1) % POC(k))
            sINP % CPOOLI_SIG(j) % ENZP(k) = &
            & fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % ENZP(k), sINP % CPOOLI(j + 1) % ENZP(k))
        end do

        sINP % CPOOLI_SIG(j) % MOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MOC, sINP % CPOOLI(j + 1) % MOC)
        sINP % CPOOLI_SIG(j) % QOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % QOC, sINP % CPOOLI(j + 1) % QOC)
        sINP % CPOOLI_SIG(j) % MBC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBC, sINP % CPOOLI(j + 1) % MBC)
        sINP % CPOOLI_SIG(j) % MBCA = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBCA, sINP % CPOOLI(j + 1) % MBCA)
        sINP % CPOOLI_SIG(j) % MBCD = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % MBCD, sINP % CPOOLI(j + 1) % MBCD)
        sINP % CPOOLI_SIG(j) % DOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % DOC, sINP % CPOOLI(j + 1) % DOC)
        sINP % CPOOLI_SIG(j) % ENZM = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % ENZM, sINP % CPOOLI(j + 1) % ENZM)
        sINP % CPOOLI_SIG(j) % CO2 = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % CO2, sINP % CPOOLI(j + 1) % CO2)
        sINP % CPOOLI_SIG(j) % SOC = fPermil(0, const_Rstd(j), sINP % CPOOLI(1) % SOC, sINP % CPOOLI(j + 1) % SOC)
    end do !j = 1, const_nISO - 1

    if (sINI % iModel .eq. 0) then !output results for model simulation
        write(sINI%iFout_VAR_hour, '(i10,13e20.3,13e20.3,13e20.3)') &
                    0, sINP % CPOOL,sINP % CPOOLI(1),sINP % CPOOLI(2)
    end if


    sINI % sINP = sINP

END !!subroutine subMEND_INI

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
    REAL(8) sum1, sum2
    INTEGER iCase, nObj, eof
    INTEGER ibeg, iend,i,j, k
    
    INTEGER iRead
    REAL(8) rRead
    CHARACTER(LEN = 10) sRead
    CHARACTER(LEN = 200) format510,format521, dirout
    CHARACTER(LEN = 200) sFile_inp
    INTEGER, ALLOCATABLE :: VARid(:)
    INTEGER, ALLOCATABLE :: VARid_unique(:,:)
    REAL(8), ALLOCATABLE :: rVARid(:)
    REAL(8), ALLOCATABLE :: dOBS_SIM(:,:)
!    INTEGER, PARAMETER :: nSub = 2 !# of substrates for model calibration, i.e., Glucose, Starch
    
!    do i = 1, sPAR%nPar
!        write(*,*)'xx=',i,xx(i)
!    end do
    
    !!Assignment of Parameter Values are moved to SUBROUTINE subMEND_PAR() 
    sINI%LCI0       = xx(1) !initial LCI
    sINI%r0         = xx(2) !fraction of active biomass
    
      !!output files for response variables
      !!--------------------------------------------------------------------------
      sFile_inp = trim(sINI%dirout)//'SIM_obs.out'
      open(unit = sINI%iFout_SIM_obs_cases, file = sFile_inp, status = 'unknown')
      write(sINI%iFout_SIM_obs_cases,*)"SIMULATION vs. OBSERVATION for ALL CASES: "
      write(sINI%iFout_SIM_obs_cases,'(2a5,a15,3a20)')"ID","VAR","Date","OBS_avg","SIM_avg","SIM_sd"
      !!--------------------------------------------------------------------------
!    print*,"ok:",sINI%VARopt_int
    do iCase=1, sINI%nCase
        call MENDIN_CASE(iCase,sINI)
        call subMEND_INI(sINI) !initialization: initial pool sizes
        call subMEND_RUN(xx, sPAR, sINI, sOUT)
    
        do i = 1, sINI%nVARopt 
            j = sINI%VARopt_int(i,1)
            if(i.eq.1) then
                ibeg = 1
            else
                ibeg = sum(sINI%VARopt_int(1:(i-1),2)) + 1
            end if
            iend = sum(sINI%VARopt_int(1:i,2))
            if(trim(sINI%VARobj(j)).eq."NSEC") then
                sINI%rOBJ(i) = f1NSE(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), const_FillValue)
            else !! "MARE"
                sINI%rOBJ(i) = fMARE(sINI%VARopt_int(i,2), sINI%dOBS_opt(ibeg:iend,2), sINI%dSIM_opt(ibeg:iend,2), const_FillValue)
            end if
        end do

        fMEND_OBJ = fWAVG(sINI%nVARopt,sINI%rOBJw,sINI%rOBJ,const_FillValue)

        do i=1,sINI%nOBS_tot 
            write(sINI%iFout_SIM_obs_cases,'(2i5,i15,3f20.6)')i,int(sINI%dOBS_opt(i,3)),int(sINI%dOBS_opt(i,1)),&
                                    sINI%dOBS_opt(i,2),sINI%dSIM_opt(i,2),sINI%dSIM_opt(i,3)   !!obs_mean,sim_mean,sim_sd
        end do

        if(sINI%iModel.eq.0) then
            WRITE(*,'(/,A10,I5)')">>>iCase = ",iCase
            do i=1,sINI%nOBS_tot 
                write(sINI%iFout_SIM_obs,'(2i5,i15,3f20.6)')i,int(sINI%dOBS_opt(i,3)),int(sINI%dOBS_opt(i,1)),&
                                    sINI%dOBS_opt(i,2),sINI%dSIM_opt(i,2),sINI%dSIM_opt(i,3)   !!obs_mean,sim_mean,sim_sd
            end do

            write(sINI%iFout_SIM_obs,'(/,a)')"PARAMETERS & CRITERIA:"
            write(format510,*)"(/,",sINI%nPar,"(a12),","' |  CRITERION'",",a10, I10)"   
!            sRead = "    OBJ-1:"
    !        write(sINI%iFout_SIM_obs,format510)sPAR_SCE%parName,sRead,sINI%nVARopt
            write(format521,*)"(",sINI%nPar,"f12.6,","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
            write(sINI%iFout_SIM_obs,format521) xx,fMEND_OBJ,sINI%rOBJ

            write(*,'(a20,f10.4)')"fOBJ (best=0) = ", fMEND_OBJ
            write(*,'(a20,6f10.4)')"fOBJ[i] = ",sINI%rOBJ
            write(*,'(a20,6f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw

            !!close all output files
            close(sINI%iFout_SIM_obs)
            close(sINI%iFout_SIM_day)
            close(sINI%iFout_SIM_mon)
            close(sINI%iFout_VAR_hour)
            close(sINI%iFout_FLX_hour)
            close(sINI%iFout_RATE_hour)
            close(sINI%iFout_PAR_hour)

            dirout = trim(sINI%dirout)//"/"//trim(sINI%CASEname(iCase))//"/"//trim(sINI%SITE)//"_"
            !!Convert OUTPUTS from HOURLY to DAILY & MONTHLY for all STATE VARIABLEs | FLUXes | RATEs
            call sOUT_tscale(dirout,sINI%sDate_beg_sim,sINI%sDate_end_sim)
        end if
    
        DEALLOCATE(sINI%STP)
        DEALLOCATE(sINI%SWC)
        DEALLOCATE(sINI%SWP)
        DEALLOCATE(sINI%SpH)
        DEALLOCATE(sINI%SIN) !!SOC input
        DEALLOCATE(sINI%VARopt)
        DEALLOCATE(sINI%VARstep)
        DEALLOCATE(sINI%VARfile)
        DEALLOCATE(sINI%VARcol)
        DEALLOCATE(sINI%VARobj)
        DEALLOCATE(sINI%VARobjw)
        DEALLOCATE(sINI%VARopt_int)
        DEALLOCATE(sINI%rOBJ)
        DEALLOCATE(sINI%rOBJw)
        DEALLOCATE(sINI%dOBS_opt) !!date,obs,iVARopt
        DEALLOCATE(sINI%dSIM_opt) !!date,sim,sim_sd
    end do !!do iCase=1, sINI%nCase
    close(sINI%iFout_SIM_obs_cases)  !!save outputs
    
    !!(1) determine # of data-points
!    sFile_inp = trim(sINI%dirout)//'SIM_obs.out'
    open(unit = sINI%iFout_SIM_obs_cases, file = sFile_inp, status = 'unknown')
    REWIND(unit = sINI%iFout_SIM_obs_cases)
    read(sINI%iFout_SIM_obs_cases,*)sRead
    read(sINI%iFout_SIM_obs_cases,*)sRead
    k = 0
    read(sINI%iFout_SIM_obs_cases,*,iostat=eof)sRead !!Line-3: 1st data
    do while (eof.ge.0) 
        k = k + 1
        read(sINI%iFout_SIM_obs_cases,*,iostat=eof)sRead
    end do             
!    close(sINI%iFout_SIM_obs_cases)
    
    !!(2) assign values for dOBS, dSIM
    ALLOCATE(VARid(k))
    ALLOCATE(VARid_unique(k,2))
    ALLOCATE(rVARid(k))
    ALLOCATE(dOBS_SIM(k,2))

    REWIND(unit = sINI%iFout_SIM_obs_cases)
!    open(unit = sINI%iFout_SIM_obs_cases, file = sFile_inp, status = 'unknown')
!   ID  VAR           Date             OBS_avg             SIM_avg              SIM_sd
!    1    1       20090820            0.001228            0.001334            0.000107
    read(sINI%iFout_SIM_obs_cases,*)sRead
    read(sINI%iFout_SIM_obs_cases,*)sRead
    do i=1,k
        read(sINI%iFout_SIM_obs_cases,*)iRead, VARid(i), sRead, dOBS_SIM(i,1:2),rRead
    end do
    close(sINI%iFout_SIM_obs_cases)
    
    rVARid = DBLE(VARid)  !!convert to real(8) for sorting
    call sort(k,2,dOBS_SIM,rVARid)
!! count unique VARid, i.e., # of variable for calibration 
    call CountUniqueElement(k,VARid,sINI%nVARopt_cases,VARid_unique)
       
!    ALLOCATE(sINI%rOBJ(sINI%nVARopt_cases))
!    ALLOCATE(sINI%rOBJw(sINI%nVARopt_cases))
    
    sINI%rOBJw_cases(1:sINI%nVARopt_cases) =  sINI%VARobjw_cases(1:sINI%nVARopt_cases)
    call Array_Normalize(sINI%nVARopt_cases,sINI%rOBJw_cases,const_FillValue)
    
    do i = 1, sINI%nVARopt_cases 
        j = VARid_unique(i,1)
        if(i.eq.1) then
            ibeg = 1
        else
            ibeg = sum(VARid_unique(1:(i-1),2)) + 1
        end if
        iend = sum(VARid_unique(1:i,2))
!        if(VARid_unique(i,2).gt.10) then  !!NSEC
        if(trim(sINI%VARobj_cases(i)).eq."NSEC") then
            sINI%rOBJ_cases(i) = f1NSE(VARid_unique(i,2), dOBS_SIM(ibeg:iend,1), dOBS_SIM(ibeg:iend,2), const_FillValue)
        else !! "MARE"
            sINI%rOBJ_cases(i) = fMARE(VARid_unique(i,2), dOBS_SIM(ibeg:iend,1), dOBS_SIM(ibeg:iend,2), const_FillValue)
        end if
    end do

    fMEND_OBJ = fWAVG(sINI%nVARopt_cases,sINI%rOBJw_cases,sINI%rOBJ_cases,const_FillValue)
!    if(sINI%iModel.eq.0) then
!        write(*,'(/,A)')">>>FINAL OBF for ALL CASES:"
!        write(*,'(a20,f10.4)')"fOBJ (best=0) = ", fMEND_OBJ
!        write(*,'(a20,6f10.4)')"fOBJ[i] = ",sINI%rOBJ
!        write(*,'(a20,6f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw
!    end if
!    do i = 1, k
!        print*,VARid_unique(i,1:2)
!    end do 
!    do i=1,k
!        write(*,'(I10,f20.0,2f20.6)')i,rVARid(i),dOBS_SIM(i,1:2)
!    end do
    
    DEALLOCATE(VARid)
    DEALLOCATE(VARid_unique)
    DEALLOCATE(rVARid)
    DEALLOCATE(dOBS_SIM)

!    DEALLOCATE(sINI%rOBJ)
!    DEALLOCATE(sINI%rOBJw)
    
    if(isnan(fMEND_OBJ)) then
        write(*,*)"NaN",sPAR
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
    
    call sInt2Str(ihr,2,str2)
    sDateHr = sDate//str2
!    if (sINI%iModel .eq. 0) then
    write(sINI%iFout_VAR_hour, '(A10,13e20.3,13e20.3,13e20.3)') &
        sDateHr, sOUT % CPOOL,sOUT % CPOOLI(1),sOUT % CPOOLI(2)
    write(sINI%iFout_FLX_hour, '(A10,50e20.3)') sDateHr, sOUT % CFLUX
    write(sINI%iFout_PAR_hour, '(A10,50e20.3)') sDateHr, sPAR
!    end if
END !!subroutine subMEND_output

!-----------------------------------------------------------------------------
SUBROUTINE subMEND_output_rate(sDate,ihr, sINI, sPAR, sINP, sOUT)
!    USE MOD_MEND
!!Hourly output for derived parameters
!!"Hour","kPOC1","kPOC2","kMOC","kDOC","kMBC","phi","rMBA"
    !!ARGUMENTS:
    TYPE(sMEND_INI), intent(in) :: sINI
    TYPE(sMEND_PAR), intent(in) :: sPAR
    TYPE(sMEND_INP), intent(in) :: sINP
    TYPE(sMEND_OUT), intent(in) :: sOUT   
    CHARACTER(LEN=8),intent(in) :: sDate
    INTEGER,         intent(in) :: ihr
    
    !!LOCAL VARIABLES 
    REAL(8) kPOC1           !!Equivalent 1st-order decomposition rate; k=Vd*ENZP/(POC + Ks)
    REAL(8) kPOC2           !!Equivalent 1st-order decomposition rate; k=Vd*ENZP/(POC + Ks)
    REAL(8) kMOC            !!Equivalent 1st-order decomposition rate; k=Vd*ENZM/(MOC + Ks)
    REAL(8) kDOC            !!Equivalent 1st-order turnover rate of DOC; k=[(Vg+Vm)/Yg]*MBC/(DOC + Ks)
    REAL(8) kMBa            !!Equivalent 1st-order turnover rate of Active Mcirobes; k=respiration_rate + mortality_rate+Enzyme_production_rate
    REAL(8) kMBa_in         !!Equivalent 1st-order assimilation rate of active microbes 
    REAL(8) kMBd            !!Equivalent 1st-order turnover rate of dormant microbes
    REAL(8) kMBd_in         !!Equivalent 1st-order microbial dormancy rate 
    REAL(8) kMB             !!Equivalent 1st-order turnover rate of microbes
    REAL(8) kMB_in          !!Equivalent 1st-order assimilation rate of total microbes
    REAL(8) phi             !!DOC saturation level; phi=DOC/(DOC+Ks)
    REAL(8) rMBa            !!Active Fraction of Microbes, r=MBCA/MBC
    REAL(8) CUE             !!Apparent Microbial Carbon Use Efficiency, CUE=[DOM_to_MBA - CO2_gmo - death]/DOM_to_MBA
    
    CHARACTER(LEN=2)  str2
    CHARACTER(LEN=10) sDateHr
    
    call sInt2Str(ihr,2,str2)
    sDateHr = sDate//str2
    
    phi = sINP%CPOOL%DOC/(sPAR%KsDOC+sINP%CPOOL%DOC)
    rMBa = sINP%CPOOL%MBCA/sINP%CPOOL%MBC
    kPOC1 = sPAR%VdPOC(1)*sINP%CPOOL%ENZP(1)/(sPAR%KsPOC(1)+sINP%CPOOL%POC(1)) 
    kPOC2 = sPAR%VdPOC(2)*sINP%CPOOL%ENZP(2)/(sPAR%KsPOC(2)+sINP%CPOOL%POC(2)) 
    kMOC = sPAR%VdMOC*sINP%CPOOL%ENZM/(sPAR%KsMOC+sINP%CPOOL%MOC)
    kDOC = (sPAR%Vg + sPAR%Vm)/sPAR%Yg*sINP%CPOOL%MBCA/(sPAR%KsDOC+sINP%CPOOL%DOC)
    kMBa = (sPAR%Vg + sPAR%Vm)*(1.D0/sPAR%Yg - 1.D0)*phi &
            + sPAR%rMORT + (sPAR%pENZP+sPAR%pENZM) * sPAR%Vm 
    kMBa_in = (sPAR%Vg + sPAR%Vm)/sPAR%Yg*phi + sOUT%CFLUX%MBCD_to_MBCA/sINP%CPOOL%MBCA
    kMBd_in = sOUT%CFLUX%MBCA_to_MBCD/sINP%CPOOL%MBCD
    kMBd = (sOUT%CFLUX%MBCD_to_MBCA + sOUT%CFLUX%CO2_maintn_dorm)/sINP%CPOOL%MBCD
    kMB = (sOUT%CFLUX%CO2_gm + sOUT%CFLUX%MBC_PM)/sINP%CPOOL%MBC
    kMB_in = sOUT%CFLUX%DOC_to_MBC/sINP%CPOOL%MBC       
    if (sOUT%CFLUX%DOC_to_MBC > 0.D0) then
        CUE = (sOUT%CFLUX%DOC_to_MBC - sOUT%CFLUX%CO2_gm - sOUT%CFLUX%MBC_PM)/sOUT%CFLUX%DOC_to_MBC
    else
        CUE = 0.D0
    end if
!    if (sINI%iModel .eq. 0) then
    write(sINI%iFout_rate_hour, '(A10,30e20.3)') &
            sDateHr, kPOC1,kPOC2,kMOC,kDOC,kMBa,kMBa_in,kMBd,kMBd_in,kMB,kMB_in,phi,rMBa,CUE, &
            sOUT%RE,sOUT%TOCbeg,sOUT%TOCend,sOUT%TOCinp,sOUT%TOCout
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
    INTEGER j, lp, jbeg,jend, tstep, ibeg, iend, iday
    INTEGER iFunc, iObs_time, iHour
    REAL(8) frISOadd(const_nISO)
    INTEGER nHour,nday,nmon
    REAL(8), ALLOCATABLE:: dSIM_d(:,:)  !!daily
    REAL(8), ALLOCATABLE:: dSIM_m(:,:)  !!monthly
    INTEGER nOBS
    
    CHARACTER(len=6) sYM
    CHARACTER(len=8) sDate
    INTEGER iYYYYMM, iyr,imo,ida, ihr
    REAL(8) dSIM_h(sINI%nHour_sim,sINI%nVARopt)  !!store hourly output
    REAL(8) dSIM(24,sINI%nVARopt)  !!store hourly output
    
    character(len=200) sFile_inp
    character(len=200) sRead
    INTEGER iRead
    REAL(8) rRead(5)

!!ATTENTION: use ASSOCIATE will result in subroutine invisible in Navigator
!    ASSOCIATE(sINP => sINI%sINP) 
    
    !!replace sDate_beg_all & sDate_end_all with sDate_beg_sim & sDate_end_sim 
    nday = nDaysbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
    nmon = nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
    ALLOCATE(dSIM_d(nday,sINI%nVARopt*2)) !!mean & sd
    ALLOCATE(dSIM_m(nmon,sINI%nVARopt*2)) !!mean & sd
    
!    sINI%LCI0       = xx(1) !initial LCI
!    sINI%r0         = xx(2) !fraction of active biomass
!     
!    call subMEND_INI(sINI) !initialization: initial pool sizes
    
    frISOadd(2) = 1 !1d0/(1d0 + sINI%SIN_C12_C14) !C14
    frISOadd(1) = 1d0 - frISOadd(2)            !C12
    
    sINP = sINI%sINP

    !    print*, sINP%CPOOL
    !    print*, sINP%CPOOLI(1)
    !    print*, sINP%CPOOLI(2)
    jbeg = 24*(nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_beg_inp2)-1)+1
    jend = 24*nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_end_inp2)
    
    !!Initial Step:
    sOUT % CPOOL    = sINP % CPOOL
    sOUT % CPOOLI   = sINP % CPOOLI
    sOUT % CPOOLIFR = sINP % CPOOLIFR
    
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
            write(*,'(a6,a10,f16.2,a,a,$)')'Date=',sDate,iday*100.0/nday,'%',cBackspace
        end if
! print*, 'iday=',iday,sDate
 
        do lp = 1, 24  !!
            i = (iday - 1)*24 + lp  !!hour
            k = ibeg + i
!            
            sINP % CPOOL = sOUT % CPOOL
            sINP % CPOOLI = sOUT % CPOOLI
            sINP % CPOOLIFR = sOUT % CPOOLIFR
    !        end if

            if(k.le.sINI%nHour) then
!!wgs[7/12/2015]: REVERT: use ARRAY instead of READING inputs at each step    
!                read(sINI%iFout_ITW_hour,*)iRead,rRead(1:5)
!                sINP%SIN = rRead(1)
!                sINP%tmp = rRead(2)
!                sINP%SWP = rRead(4)
!                sINP%pH  = rRead(5)
!!wgs: END
                
                sINP%SIN = sINI%SIN(k)
                sINP%tmp = sINI%STP(k)
                sINP%SWP = sINI%SWP(k)
                sINP%pH  = sINI%SpH(k)
                
                !!External INPUT: BEGIN
                sINP % CADD % POCadd(1) = sINP%SIN * sINI%SIN_frac(1) + sINI%SIN_other(1,1)  !![mg POC/g soil/h], inputs to POC
                sINP % CADD % POCadd(2) = sINP%SIN * sINI%SIN_frac(2) + sINI%SIN_other(1,2)
                sINP % CADD % DOCadd    = sINP%SIN * sINI%SIN_frac(3) + sINI%SIN_other(1,3)  !![mg DOC/g soil/h], inputs to DOC

                if(i.ge.jbeg.and.i.le.jend) then
                    sINP%CADD%POCadd(1) = sINP%CADD%POCadd(1)+sINI%SIN_other(2,1)/DBLE(jend-jbeg+1)
                    sINP%CADD%POCadd(2) = sINP%CADD%POCadd(2)+sINI%SIN_other(2,2)/DBLE(jend-jbeg+1)
                    sINP%CADD%DOCadd    = sINP%CADD%DOCadd   +sINI%SIN_other(2,3)/DBLE(jend-jbeg+1)
                end if
        !        write(*,*)i,sINP % CADD % POCadd(1:2),sINP % CADD % DOCadd
            else !! POST-DATA PERIOD: simulation period w/o input data, sINI%nHour<i<=sINI%nHour_sim
                j = mod(k,sINI%nHour)+1  
                !!use available data to fill the missing data after the available period (>sDate_end_all)
                sINP%tmp = sINI%STP(j)
                sINP%SWP = sINI%SWP(j)
                sINP%pH  = sINI%SpH(j)

                !!External INPUT: BEGIN; assume root input = leaf litter input
                sINP % CADD % POCadd(1) = 0 !sINI%SIN(j) * sINI%SIN_frac(1) * sINI%SIN_Multiplier !!+ sINI%SIN_other(1,1)  !![mg POC/g soil/h], inputs to POC
                sINP % CADD % POCadd(2) = 0 !sINI%SIN(j) * sINI%SIN_frac(2) * sINI%SIN_Multiplier !!+ sINI%SIN_other(1,2)
                sINP % CADD % DOCadd    = 0 !sINI%SIN(j) * sINI%SIN_frac(3) * sINI%SIN_Multiplier !!+ sINI%SIN_other(1,3)  !![mg DOC/g soil/h], inputs to DOC
            end if

            do j = 1, const_nPOC
                sINP % CADDI % POCadd(j) = sINP%CADD%POCadd(j)*frISOadd
            end do
            sINP % CADDI % DOCadd = sINP%CADD%DOCadd*frISOadd
            !!External INPUT: END

            sINI%sINP = sINP
            call subMEND_PAR(xx, sPAR, sINI)   !!modify parameter values by soil temperature, water potential, pH
            call subMEND(sPAR, sINI, sOUT) !!MEND model 
            if (sINI % iModel .eq. 0) then !output results for model simulation
                call subMEND_output_rate(sDate,lp, sINI, sPAR, sINP, sOUT)
                call subMEND_output(sDate,lp, sPAR, sINI, sOUT)
            end if

            !!Extract output variables to compare with observations 
            call sOUT_OPT_h(sINI%nVARopt,24,lp,dSIM,sOUT,sINI%VARopt_int,sPAR%VdPOC(1:2))
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
            call sDate_After(i,sINI%sDate_beg_sim,sDate)
            write(sINI%iFout_SIM_day,'(A10,20f20.6)')sDate,dSIM_d(i,:)
        end do

        do i = 1,nmon
            call sYM_After(i,sINI%sDate_beg_sim(1:6),sYM)
            write(sINI%iFout_SIM_mon,'(A10,20f20.6)')sYM,dSIM_m(i,:)
        end do
    end if
!    
!    do i = 1,nmon
!        write(*,*)i,dSIM_m(i,:)
!    end do
    
    !!Extract those data matching the observations on the same days or months
    call sDate2YMD(sINI%sDate_beg_sim,iyr,imo,ida)
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
                sINI%dSIM_opt(lp,2) = dSIM_d(i,j*2-1)
                sINI%dSIM_opt(lp,3) = dSIM_d(i,j*2)
            else if(tstep.eq.2) then  !!monthly
                iYYYYMM = int(sINI%dOBS_opt(lp,1))
                i = nMonths(iyr*100+imo,iYYYYMM)
!                sINI%dSIM_opt(lp,1) = sINI%dOBS_opt(lp,1)
                sINI%dSIM_opt(lp,2) = dSIM_m(i,j*2-1)
                sINI%dSIM_opt(lp,3) = dSIM_m(i,j*2)
            else if(tstep.eq.4) then  !!yearly
                iend = int(sINI%dOBS_opt(lp,1)) !!year
                ibeg = max(1,nMonths(iyr*100+imo,iend*100+1))       !!in case of: imo<>1
                iend = min(nmon,nMonths(iyr*100+imo,iend*100+12))   !!in case of: YYYYMM_end < iend*100 + 12
                sINI%dSIM_opt(lp,2) = fAVG((iend-ibeg+1),dSIM_m(ibeg:iend,j*2-1),const_FillValue)  
                sINI%dSIM_opt(lp,3) = fSTDDEV((iend-ibeg+1),dSIM_m(ibeg:iend,j*2-1),const_FillValue)
            else if(tstep.eq.5) then  !!overall mean value, k=1
                sINI%dSIM_opt(lp,2) = fAVG(nmon,dSIM_m(:,j*2-1),const_FillValue)  
                sINI%dSIM_opt(lp,3) = fSTDDEV(nday,dSIM_d(:,j*2-1),const_FillValue)
            end if
        end do
    end do
    
    DEALLOCATE(dSIM_d)
    DEALLOCATE(dSIM_m)
!    END ASSOCIATE
!    do i=1,sINI%nOBS_tot 
!        write(*,*)i,sINI%dOBS_opt(i,:),sINI%dSIM_opt(i,:)
!    end do
END !!subroutine subMEND_RUN
!-----------------------------------------------------------------------------

SUBROUTINE subMEND(sPAR, sINI, sOUT)
    !model components 
!    USE MOD_MEND
!    IMPLICIT NONE
    !!ARGUMENTS:
    TYPE(sMEND_PAR), intent(in) :: sPAR
    TYPE(sMEND_INI), intent(in) :: sINI
    TYPE(sMEND_OUT), intent(out):: sOUT
!    REAL(8)        , intent(in) :: xx(sPAR % nPar)  !!used in subMEND_PAR
    
    !!LOCAL VARIABLES:
    TYPE(sMEND_INP) sINP
    TYPE(sMM_PAR) smmPAR
    TYPE(sMM_INP) smmINP

    TYPE(sSORP_PAR) sSorpPAR
    TYPE(sSORP_INP) sSorpINP
    TYPE(sSORP_OUT) sSorpOUT

    INTEGER i, j, k, iFunc
    REAL(8) OC1, OC2, sumPOC, sumENZP, phi
    REAL(8) frPOC(const_nPOC), frENZP(const_nPOC), POCdec(const_nPOC), ENZP_loss(const_nPOC)
    
!    ASSOCIATE(sINP => sINI%sINP) 
    sINP = sINI%sINP

!    write(*,'(a15,50d10.3)')"sINP%CPOOL=",sINP%CPOOL
    !compute initial total mass of the system 
    sOUT%TOCbeg = sINP % CPOOL % SOC + sINP % CPOOL % DOC + sINP % CPOOL % MBC &
                + sum(sINP % CPOOL % ENZP) + sINP % CPOOL % ENZM
    sOUT%TOCinp = sum(sINP % CADD % POCadd) + sINP % CADD % DOCadd
    !Flux 2: POC decomposition
    DO i = 1, const_nPOC
        smmINP % substrate = sINP % CPOOL % POC(i)
        smmINP % enzyme = sINP % CPOOL % ENZP(i)
        smmPAR % vm = sPAR % VdPOC(i)
        smmPAR % km = sPAR % KsPOC(i)
        POCdec(i) = fV_MM(sINI%iKinetics, smmPAR, smmINP)
        sOUT % CFLUX % POCdec(i) = min(POCdec(i), sINP % CPOOL % POC(i))
        sOUT%CFLUX%POCdec_to_DOC(i) = sPAR%frPOC2DOC*sOUT % CFLUX % POCdec(i) 
        sOUT%CFLUX%POCdec_to_MOC(i) = (1.d0 - sPAR%frPOC2DOC) * sOUT%CFLUX%POCdec(i) 
    end DO!for i = 1:sINP%nPOC

    !Flux 3: MOC decomposition
    smmINP % substrate = sINP % CPOOL % MOC
    smmINP % enzyme = sINP % CPOOL % ENZM
    smmPAR % vm = sPAR % VdMOC
    smmPAR % km = sPAR % KsMOC
    sOUT % CFLUX % MOC_to_DOC = fV_MM(sINI%iKinetics, smmPAR, smmINP) !MOC decomposition
    sOUT % CFLUX % MOC_to_DOC = min(sOUT % CFLUX % MOC_to_DOC, sINP % CPOOL % MOC)

   !Flux 1: DOC uptake by MB
    smmINP % substrate = sINP % CPOOL % DOC
    smmINP % enzyme = sINP % CPOOL % MBCA
    smmPAR % vm = (sPAR % Vg + sPAR % Vm)/sPAR % Yg !GROWTH + MAINTENANCE    
    smmPAR % km = sPAR % KsDOC
    sOUT % CFLUX % DOC_to_MBC = fV_MM(0, smmPAR, smmINP)
    
    if(sOUT % CFLUX % DOC_to_MBC.gt.sINP % CPOOL % DOC) then
!        print*,"DEBUG:",sOUT % CFLUX % DOC_to_MBC,sINP % CPOOL % DOC,min(sOUT % CFLUX % DOC_to_MBC,sINP % CPOOL % DOC)
        sOUT % CFLUX % DOC_to_MBC = sINP % CPOOL % DOC       
    end if
    !!FIRST UPDATE, used for sorption-desorption calculation
    sOUT % CPOOL % DOC = sINP % CPOOL % DOC - sOUT % CFLUX % DOC_to_MBC  
    
    !Flux 4 & 5: Growth and Maintenance Respiration
    !emission of CO2 from soil to air NOT considered  
    sOUT % CFLUX % CO2_maintn_dorm = sPAR%VmD * sINP%CPOOL%MBCD  !!Microbial Maintenance of Dormant Microbes 
    
!    sOUT % CFLUX % CO2_growth  = sOUT % CFLUX % DOC_to_MBC * (1.D0 - sPAR % Yg) * sPAR % Vg/(sPAR % Vg + sPAR % Vm)
!    sOUT % CFLUX % CO2_maintn  = sOUT % CFLUX % DOC_to_MBC * (1.D0 - sPAR % Yg) - sOUT % CFLUX % CO2_growth
    
    smmINP % substrate = sINP % CPOOL % DOC
    smmINP % enzyme = sINP % CPOOL % MBCA
    smmPAR % vm = sPAR % Vg * (1.D0/sPAR % Yg - 1.D0)
    smmPAR % km = sPAR % KsDOC
    sOUT % CFLUX % CO2_growth = fV_MM(0, smmPAR, smmINP)

    smmPAR % vm = sPAR % Vm * (1.D0/sPAR % Yg - 1.D0)
    sOUT % CFLUX % CO2_maintn = fV_MM(0, smmPAR, smmINP) 
                                    
    sOUT % CFLUX % CO2_gm = sOUT%CFLUX%CO2_growth + sOUT%CFLUX%CO2_maintn + sOUT%CFLUX%CO2_maintn_dorm

    !Flux 6 & 7: DOC adsorption-desorption
    sSorpINP % adsorbate = sOUT % CPOOL % DOC !!DOC has been update;  = sINP%CPOOL%DOC - sOUT%CFLUX%DOC_to_MBC
    sSorpINP % adsorbent = sINP % CPOOL % QOC
    sSorpPAR % Qmax = sPAR % Qmax
    sSorpPAR % Kads = sPAR % Kads
    sSorpPAR % Kdes = sPAR % Kdes
    iFunc = fAds(sSorpPAR, sSorpINP, sSorpOUT)
    sOUT % CFLUX % QOC_to_DOC = sSorpOUT % des
    sOUT % CFLUX % DOC_to_QOC = sSorpOUT % ads
    sOUT % CFLUX % DOC_to_QOC_net = sSorpOUT % ads_net

    !Flux 8: microbial mortality
!    sOUT % CFLUX % MBC_PM = sPAR % Vm * sINP % CPOOL % MBCA
!    sOUT % CFLUX % MBC_mortality = (1.d0 - sPAR % pENZP - sPAR % pENZM) * sOUT % CFLUX % MBC_PM
    sOUT % CFLUX % MBC_mortality = sPAR%rMORT * sINP%CPOOL%MBCA
    sOUT % CFLUX % MBC_to_DOC = sPAR % frMBC2DOC * sOUT % CFLUX % MBC_mortality
    sOUT % CFLUX % MBC_to_POC = (1.d0 - sPAR % frMBC2DOC) * sOUT % CFLUX % MBC_mortality
    !Flux 9: enzyme synthesis
    sumPOC = sum(sINP % CPOOL % POC)
    frPOC = sINP % CPOOL % POC/sumPOC !compute the fraction of each type of POC, e%g%, lignin, cellulose 
!    sOUT % CFLUX % MBC_to_ENZP = frPOC * sPAR % pENZP * sOUT % CFLUX % MBC_PM
!    sOUT % CFLUX % MBC_to_ENZM = sPAR % pENZM * sOUT % CFLUX % MBC_PM
    sOUT % CFLUX % MBC_to_ENZP = frPOC * sPAR % pENZP * sPAR % Vm * sINP % CPOOL % MBCA
    sOUT % CFLUX % MBC_to_ENZM =         sPAR % pENZM * sPAR % Vm * sINP % CPOOL % MBCA
    
    sOUT % CFLUX % MBC_PM = sOUT % CFLUX % MBC_mortality + sum(sOUT % CFLUX % MBC_to_ENZP) + sOUT % CFLUX % MBC_to_ENZM

    !Flux 10: enzyme turnover
    sOUT % CFLUX % ENZP_to_DOC = sPAR % rENZP * sINP % CPOOL % ENZP
    sOUT % CFLUX % ENZM_to_DOC = sPAR % rENZM * sINP % CPOOL % ENZM
    
    !Internal Flux: microbial dormancy and reactivation
    phi = sINP % CPOOL % DOC/(sINP % CPOOL % DOC + sPAR % KsDOC)
    sOUT % CFLUX % MBCA_to_MBCD = (1.d0 - phi) * sPAR % VmA2D * sINP % CPOOL % MBCA
    sOUT % CFLUX % MBCD_to_MBCA = phi * sPAR % VmD2A * sINP % CPOOL % MBCD

    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !Mass Balance
    !Total Carbon
!    write(1001,*)"fMEND = ",sum(sINP % CADD % POCadd)+sINP % CADD % DOCadd, sOUT % CFLUX % CO2_gm
    
    sOUT % CPOOL % POC      = sINP % CPOOL % POC + sINP % CADD % POCadd - sOUT % CFLUX % POCdec !Lignin and Cellulose
    sOUT % CPOOL % POC(1)   = sOUT % CPOOL % POC(1) + sOUT % CFLUX % MBC_to_POC !ONLY Lignin!     
    sOUT % CPOOL % MOC      = sINP % CPOOL % MOC + sum(sOUT % CFLUX % POCdec_to_MOC) - sOUT % CFLUX % MOC_to_DOC
    sOUT % CPOOL % QOC      = sINP % CPOOL % QOC + sOUT % CFLUX % DOC_to_QOC_net !QOC mass balance
    sOUT % CPOOL % MBC      = sINP % CPOOL % MBC + sOUT % CFLUX % DOC_to_MBC &
    &                       - sOUT % CFLUX % CO2_gm - sOUT % CFLUX % MBC_PM
    sOUT % CPOOL % MBCA     = sINP % CPOOL % MBCA + sOUT % CFLUX % DOC_to_MBC &
    &                       - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) &
    &                       - sOUT % CFLUX % MBC_PM &
    &                       - sOUT % CFLUX % MBCA_to_MBCD + sOUT % CFLUX % MBCD_to_MBCA
    sOUT % CPOOL % MBCD     = sINP % CPOOL % MBCD - sOUT % CFLUX % CO2_maintn_dorm &
    &                       + sOUT % CFLUX % MBCA_to_MBCD - sOUT % CFLUX % MBCD_to_MBCA
    sOUT % CPOOL % DOC      = sINP % CPOOL % DOC + sINP % CADD % DOCadd &
    &                       + sum(sOUT % CFLUX % POCdec_to_DOC) &
    &                       + sOUT % CFLUX % MBC_to_DOC &
    &                       + sOUT % CFLUX % MOC_to_DOC &
    &                       + sum(sOUT % CFLUX % ENZP_to_DOC) &
    &                       + sOUT % CFLUX % ENZM_to_DOC &
    &                       - sOUT % CFLUX % DOC_to_MBC &
    &                       - sOUT % CFLUX % DOC_to_QOC_net
    sOUT % CPOOL % ENZP     = sINP % CPOOL % ENZP + sOUT % CFLUX % MBC_to_ENZP - sOUT % CFLUX % ENZP_to_DOC
    sOUT % CPOOL % ENZM     = sINP % CPOOL % ENZM + sOUT % CFLUX % MBC_to_ENZM - sOUT % CFLUX % ENZM_to_DOC
    sOUT % CPOOL % CO2      = sINP % CPOOL % CO2 + sOUT % CFLUX % CO2_gm
    sOUT % CPOOL % SOC      = sum(sOUT % CPOOL % POC) + sOUT % CPOOL % MOC + sOUT % CPOOL % QOC

    do j = 2, const_nISO !j = 1 for lightest isotope, e.g., C12
        !Isotopic components

        sOUT % CPOOLI(j) % POC      = sINP % CPOOLI(j) % POC + sINP % CADDI(j) % POCadd &
        &                           - sOUT % CFLUX % POCdec * sINP % CPOOLIFR(j) % POC !Lignin and Cellulose                             
        sOUT % CPOOLI(j) % POC(1)   = sOUT % CPOOLI(j) % POC(1) & 
        &                           + sOUT % CFLUX % MBC_to_POC * sINP % CPOOLIFR(j) % MBCA !ONLY Lignin!    
        sOUT % CPOOLI(j) % MOC      = sINP % CPOOLI(j) % MOC &
        &                           + sum(sOUT % CFLUX % POCdec_to_MOC * sINP % CPOOLIFR(j) % POC) &
        &                           - sOUT % CFLUX % MOC_to_DOC * sINP % CPOOLIFR(j) % MOC
        sOUT % CPOOLI(j) % QOC      = sINP % CPOOLI(j) % QOC &
        &                           + sOUT % CFLUX % DOC_to_QOC * sINP % CPOOLIFR(j) % DOC &
        &                           - sOUT % CFLUX % QOC_to_DOC * sINP % CPOOLIFR(j) % QOC!QOC mass balance
        sOUT % CPOOLI(j) % MBC      = sINP % CPOOLI(j) % MBC &
        &                           + sOUT % CFLUX % DOC_to_MBC * sINP % CPOOLIFR(j) % DOC &
        &                           - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn + sOUT % CFLUX % MBC_PM) &
        &                           * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBCD
        sOUT % CPOOLI(j) % MBCA     = sINP % CPOOLI(j) % MBCA &
        &                           + sOUT % CFLUX % DOC_to_MBC * sINP % CPOOLIFR(j) % DOC &
        &                           - (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn + sOUT % CFLUX % MBC_PM) &
        &                           * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % MBCA_to_MBCD * sINP % CPOOLIFR(j) % MBCA &
        &                           + sOUT % CFLUX % MBCD_to_MBCA * sINP % CPOOLIFR(j) % MBCD
        sOUT % CPOOLI(j) % MBCD     = sINP % CPOOLI(j) % MBCD &
        &                           - sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBCD &
        &                           + sOUT % CFLUX % MBCA_to_MBCD * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % MBCD_to_MBCA * sINP % CPOOLIFR(j) % MBCD
        sOUT % CPOOLI(j) % DOC      = sINP % CPOOLI(j) % DOC + sINP % CADDI(j) % DOCadd &
        &                           + sum(sOUT % CFLUX % POCdec_to_DOC * sINP % CPOOLIFR(j) % POC) &
        &                           + sOUT % CFLUX % MBC_to_DOC * sINP % CPOOLIFR(j) % MBCA &
        &                           + sOUT % CFLUX % MOC_to_DOC * sINP % CPOOLIFR(j) % MOC &
        &                           + sum(sOUT % CFLUX % ENZP_to_DOC * sINP % CPOOLIFR(j) % ENZP) &
        &                           + sOUT % CFLUX % ENZM_to_DOC * sINP % CPOOLIFR(j) % ENZM &
        &                           - sOUT % CFLUX % DOC_to_MBC * sINP % CPOOLIFR(j) % DOC &
        &                           - sOUT % CFLUX % DOC_to_QOC * sINP % CPOOLIFR(j) % DOC &
        &                           + sOUT % CFLUX % QOC_to_DOC * sINP % CPOOLIFR(j) % QOC
        sOUT % CPOOLI(j) % ENZP     = sINP % CPOOLI(j) % ENZP &
        &                           + sOUT % CFLUX % MBC_to_ENZP * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % ENZP_to_DOC * sINP % CPOOLIFR(j) % ENZP
        sOUT % CPOOLI(j) % ENZM     = sINP % CPOOLI(j) % ENZM &
        &                           + sOUT % CFLUX % MBC_to_ENZM * sINP % CPOOLIFR(j) % MBCA &
        &                           - sOUT % CFLUX % ENZM_to_DOC * sINP % CPOOLIFR(j) % ENZM
        sOUT % CPOOLI(j) % CO2      = sINP % CPOOLI(j) % CO2 &
                                    + (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) * sINP % CPOOLIFR(j) % MBCA &
        &                           + sOUT % CFLUX % CO2_maintn_dorm * sINP % CPOOLIFR(j) % MBCD

        !!CO2_C14 or C13 flux
        sOUT%CO2_gm_iso             = (sOUT % CFLUX % CO2_growth + sOUT % CFLUX % CO2_maintn) * sINP%CPOOLIFR(2)%MBCA &
        &                           + sOUT % CFLUX % CO2_maintn_dorm * sINP%CPOOLIFR(2)%MBCD
        
        sOUT % CPOOLI(j) % SOC      = sum(sOUT % CPOOLI(j) % POC) &
        &                           + sOUT % CPOOLI(j) % MOC + sOUT % CPOOLI(j) % QOC

    end do !j = 1, const_nISO

    !j = 1, compute the lightest isotope
    DO i = 1, const_nPOC
        sOUT % CPOOLI(1) % POC(i) = sOUT % CPOOL % POC(i) - sum(sOUT % CPOOLI(2:const_nISO) % POC(i)) !Lignin and Cellulose
        sOUT % CPOOLI(1) % ENZP(i) = sOUT % CPOOL % ENZP(i) - sum(sOUT % CPOOLI(2:const_nISO) % ENZP(i))
    END DO
    sOUT % CPOOLI(1) % MOC = sOUT % CPOOL % MOC - sum(sOUT % CPOOLI(2:const_nISO) % MOC)
    sOUT % CPOOLI(1) % QOC = sOUT % CPOOL % QOC - sum(sOUT % CPOOLI(2:const_nISO) % QOC)
    sOUT % CPOOLI(1) % MBC = sOUT % CPOOL % MBC - sum(sOUT % CPOOLI(2:const_nISO) % MBC)
    sOUT % CPOOLI(1) % MBCA = sOUT % CPOOL % MBCA - sum(sOUT % CPOOLI(2:const_nISO) % MBCA)
    sOUT % CPOOLI(1) % MBCD = sOUT % CPOOL % MBCD - sum(sOUT % CPOOLI(2:const_nISO) % MBCD)
    sOUT % CPOOLI(1) % DOC = sOUT % CPOOL % DOC - sum(sOUT % CPOOLI(2:const_nISO) % DOC)
    sOUT % CPOOLI(1) % ENZM = sOUT % CPOOL % ENZM - sum(sOUT % CPOOLI(2:const_nISO) % ENZM)
    sOUT % CPOOLI(1) % CO2 = sOUT % CPOOL % CO2 - sum(sOUT % CPOOLI(2:const_nISO) % CO2)
    sOUT % CPOOLI(1) % SOC = sOUT % CPOOL % SOC - sum(sOUT % CPOOLI(2:const_nISO) % SOC)

    do j = 1, const_nISO
        !Isotopic fractions
        sOUT % CPOOLIFR(j) % POC = sOUT % CPOOLI(j) % POC/sOUT % CPOOL % POC
        sOUT % CPOOLIFR(j) % MOC = sOUT % CPOOLI(j) % MOC/sOUT % CPOOL % MOC
        sOUT % CPOOLIFR(j) % QOC = sOUT % CPOOLI(j) % QOC/sOUT % CPOOL % QOC
        sOUT % CPOOLIFR(j) % MBC = sOUT % CPOOLI(j) % MBC/sOUT % CPOOL % MBC
        sOUT % CPOOLIFR(j) % MBCA = sOUT % CPOOLI(j) % MBCA/sOUT % CPOOL % MBCA
        sOUT % CPOOLIFR(j) % MBCD = sOUT % CPOOLI(j) % MBCD/sOUT % CPOOL % MBCD
        sOUT % CPOOLIFR(j) % DOC = sOUT % CPOOLI(j) % DOC/sOUT % CPOOL % DOC
        sOUT % CPOOLIFR(j) % ENZP = sOUT % CPOOLI(j) % ENZP/sOUT % CPOOL % ENZP
        sOUT % CPOOLIFR(j) % ENZM = sOUT % CPOOLI(j) % ENZM/sOUT % CPOOL % ENZM
        sOUT % CPOOLIFR(j) % CO2 = sOUT % CPOOLI(j) % CO2/sOUT % CPOOL % CO2
        sOUT % CPOOLIFR(j) % SOC = sOUT % CPOOLI(j) % SOC/sOUT % CPOOL % SOC

    end do !j = 1, const_nISO        

    !    print*, sum(sOUT%CPOOLIFR%MBC)!sum(sOUT%CPOOLI(2:const_nISO)%MBC),sOUT%CPOOLI(1:const_nISO)%MBC

    !Convert Isotope Concentration to Signatures [Permil]
    do j = 1, const_nISO - 1

        do k = 1, const_nPOC
            sOUT % CPOOLI_SIG(j) % POC(k) = &
            & fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % POC(k), sOUT % CPOOLI(j + 1) % POC(k))
            sOUT % CPOOLI_SIG(j) % ENZP(k) = &
            & fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % ENZP(k), sOUT % CPOOLI(j + 1) % ENZP(k))
        end do

        sOUT % CPOOLI_SIG(j) % MOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MOC, sOUT % CPOOLI(j + 1) % MOC)
        sOUT % CPOOLI_SIG(j) % QOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % QOC, sOUT % CPOOLI(j + 1) % QOC)
        sOUT % CPOOLI_SIG(j) % MBC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBC, sOUT % CPOOLI(j + 1) % MBC)
        sOUT % CPOOLI_SIG(j) % MBCA = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBCA, sOUT % CPOOLI(j + 1) % MBCA)
        sOUT % CPOOLI_SIG(j) % MBCD = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % MBCD, sOUT % CPOOLI(j + 1) % MBCD)
        sOUT % CPOOLI_SIG(j) % DOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % DOC, sOUT % CPOOLI(j + 1) % DOC)
        sOUT % CPOOLI_SIG(j) % ENZM = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % ENZM, sOUT % CPOOLI(j + 1) % ENZM)
        sOUT % CPOOLI_SIG(j) % CO2 = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % CO2, sOUT % CPOOLI(j + 1) % CO2)
        sOUT % CPOOLI_SIG(j) % SOC = fPermil(0, const_Rstd(j), sOUT % CPOOLI(1) % SOC, sOUT % CPOOLI(j + 1) % SOC)
    end do !j = 1, const_nISO - 1

    !mass balance check for the whole system
    sOUT%TOCend = sOUT % CPOOL % SOC + sOUT % CPOOL % DOC + sOUT % CPOOL % MBC &
                + sum(sOUT % CPOOL % ENZP) + sOUT % CPOOL % ENZM 
    sOUT%TOCout = sOUT % CFLUX % CO2_gm
    sOUT%RE = (sOUT%TOCend - sOUT%TOCbeg) - (sOUT%TOCinp - sOUT%TOCout)
    if(isnan(sOUT % CPOOL%MBCA)) then
!        write(*,*)sPAR
        write(*,'(a15,50d10.3)')"sOUT%CPOOL=",sOUT%CPOOL
        write(*,'(a15,50d10.3)')"sOUT%CFLUX=",sOUT%CFLUX
        stop
    end if
    
    if(dabs(sOUT % RE).gt.1.D-8) then
        print*,"Balance Check ERROR: sOUT%RE=",sOUT % RE
    else
!        print*,"Balance Check: ",sOUT % RE
    end if
    
!    END ASSOCIATE
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
    REAL(8) CUE_slope,CUE_ref, SWPmin !!CUE_slope = (-0.005, -0.016)  
    REAL(8) Vm0  !!specific maintenance rate without any modification
    
    REAL(8) tp_scalar       !!temperature scalar
    REAL(8) wp_scalar       !!INCrease with increasing wp 
    REAL(8) wp_scalar_low   !!DECrease with increasing wp 
    REAL(8) wp_scalar_opt   !!with increasing wp, INcrease first to OPTimal condition, then decrease
    REAL(8) pH_scalar   !!pH scalar
    
    REAL(8) tp
    REAL(8) wp
    REAL(8) pH
    
    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    CHARACTER(len=10) sCase
    
    SWPmin = -13.86d0  !!for Microbial Mortality
    
!    sPAR%iKinetics = sINI%iKinetics
    BIOME = trim(sINI%BIOME)
    SOM   = trim(sINI%SOM)
    
    tp = sINI%sINP%tmp
    wp = sINI%sINP%SWP
    pH = sINI%sINP%pH
    
    !![1] DECOMPOSITION of POC1 & POC2
    sCase = "LIG"
    tp_scalar       = fTArh(sCase, tp, const_Tref)
    wp_scalar_opt   = fSWP_OPT(wp)
    pH_scalar       = fpH(sCase,pH)
    sPAR % VdPOC(1)     = xx(3)*tp_scalar*wp_scalar_opt*pH_scalar !!(/xx(1), xx(2)/) ![mg POC/mg ENZP/h], maximum reaction rate for conversion of POC by ENZP
    
    sCase = "CEL"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
!    sPAR % VdPOC(2)     = xx(4)*tp_scalar*wp_scalar*pH_scalar
    sPAR % VdPOC(2)     = xx(3)*tp_scalar*wp_scalar*pH_scalar  !!set Vd_cel = Vd_lig
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % KsPOC        = (/xx(6), xx(7)/)*tp_scalar ![mg POC/cm3], half-saturation constant for conversion of POC by ENZP

    !![2] DECOMPOSITION of MOC
    sCase = "ENZ"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    wp_scalar = fSWP(wp, BIOME, SOM)
    pH_scalar = fpH(sCase,pH)
!    sPAR % VdMOC        = xx(5)*tp_scalar*wp_scalar*pH_scalar ![mg MOC/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC
    sPAR % VdMOC        = xx(3)*tp_scalar*wp_scalar*pH_scalar 
    
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % KsMOC        = xx(8)*tp_scalar ![mg MOC/cm3], half-saturation constant for conversion of MAOC by ENZMAOC

    !![3] ADSORPTION-DESORPTION
    sPAR % Qmax         = xx(9) ![mg C/g soil], adsorption capacity
    sPAR % Kba          = xx(10) ![mg C/g soil/h], binding affinity
    
    sCase = "Kdes"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % Kdes         = xx(11)*tp_scalar ![mg DOC/h],desorption rate constant
    
    sCase = "Kads"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % Kads         = (xx(11) * sPAR % Kba)*tp_scalar ![mg DOC/mg DOC/h], adsorption rate constant = Kba*Kdes
    
    !![4] ENZYME TURNOVER & PRODUCTION
    sPAR % rENZP        = (/xx(12), xx(12)/) !![1/h],turnover rate of ENZP
    sPAR % rENZM        = xx(12) !![1/h], turnover rate of ENZMAOC
    sPAR % pENZP        = xx(13) ![mg ENZP/mg MBC/h], production rate of ENZP 
    sPAR % pENZM        = xx(13)*xx(14) !!wgs:6/22/2015: set pENZM = xx(16)*pENZP, xx(16) ![mg ENZM/mg MBC/h], production rate of ENZMAOC
    
    !![5] ALLOCATION to DOC
    sPAR % frPOC2DOC    = xx(15) ![-], fraction of decomposed POC allocated to DOC
    sPAR % frMBC2DOC    = xx(16) ![-], fraction of dead MBC allocated to DOC
   
    !![6] MICROBIAL UPTAKE of DOM
    sCase = "DOM"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % Vg        = xx(17)*tp_scalar !![mg DOC/mg MBC/h], maximum uptake rate of DOC by MB
    sPAR % alpha     = xx(18) ![-], alpha = Vm/(Vm + Vg) <= 0.5 
    
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    Vm0 = xx(17) * sPAR%alpha/(1.D0 - sPAR%alpha)
    sPAR % Vm         = Vm0*tp_scalar !!*wp_scalar ![1/h], specific microbial maintenance rate
     
    sCase = "Km"
    tp_scalar = fTArh(sCase, tp, const_Tref)  
    sPAR % KsDOC        = xx(19)*tp_scalar ![mg DOC/cm3], half-saturation constant for uptake of DOC by MB
    
    sCase = "CUE"  !!label, no-use
    CUE_ref   = xx(20)
    CUE_slope = -1.D0*xx(21)
    sPAR % Yg = fT_CUE(tp, const_Tref, CUE_slope, CUE_ref) ![-], carbon use efficiency in uptake of DOC by MB 
    
    !![7] MICROBIAL MORTALITY
!    sPAR%wdie = xx(22)
    sPAR%gamma= xx(23) 
!    wp_scalar_low =  fSWP_Death(wp,SWPmin,sPAR%wdie)  
!    sPAR % rMORT = (1.D0 - sPAR%pENZP - sPAR%pENZM)*Vm0*wp_scalar_low !!wgs: 6/19/2015,  !!tp_scalar
!    wp_scalar_low = fSWP_A2D(wp, sPAR % SWP_A2D, sPAR%wdorm)  !!negative effect  !!wdie
       
    !![8] MICROBIAL DORMANCY & RESUSCITATION
    sPAR % beta       = xx(24) !beta = Vm_dormant/Vm
    sPAR % VmD        = sPAR % beta*Vm0 !!*wp_scalar_low  !!<5/15/2015>
    sPAR % SWP_A2D      = xx(25)
    sPAR % tau          = xx(26)
    sPAR % SWP_D2A      = sPAR % tau * sPAR % SWP_A2D 
    sPAR % wdorm        = xx(27)
    
    wp_scalar     = fSWP_D2A(wp, sPAR % SWP_D2A, sPAR % wdorm)  !!positive effect, increase with increasing SWC/SWP
    wp_scalar_low = fSWP_A2D(wp, sPAR % SWP_A2D, sPAR % wdorm)  !!negative effect
    sCase = "MR"
    tp_scalar = fTArh(sCase, tp, const_Tref)
    sPAR % VmA2D = Vm0*tp_scalar*wp_scalar_low 
    sPAR % VmD2A = Vm0*tp_scalar*wp_scalar
    
    sPAR % rMORT = (sPAR%gamma*Vm0)*wp_scalar_low 

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
        CASE (1)  !!First Order
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
    !Adsorption and Desorption of DOC
    !eg: adsorbent: QOC, adsorbate: DOC
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
    CHARACTER(*), intent(inout) :: sCase
    REAL(8)          , intent(in) :: T, Tref
    
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
!!-----------------------------------------------------------------
!! TEMPERATURE SCALAR: END
!!=================================================================

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
    REAL(8) SWP ![MPa]
    CHARACTER(len=3) BIOME  !!'ASM' = arid/semiarid/mediterranean; 'MGC'=mesic grassland & cropland; 'MDF'=Mesic Deciduous Forest; 'MCF'=MEsic Conifer Forest
    CHARACTER(len=3) SOM    !!'SOD' = disturbed soil; 'SOL'=intact soil; 'LIT'=litter
    
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
        CASE ("ASM")  !!Arid, Semiarid, & MediterManean
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
SUBROUTINE sOUT_OPT_h(nVAR,nHour,iHour,dSIM,sOUT,VARopt_int,vENZ)
    !!EXTRACT hourly output for opt RESPONSE VARIABLES
    !!ARGUMENTS:
    INTEGER nVAR
    INTEGER iHour,nHour
    INTEGER VARopt_int(nVAR,3)
    REAL(8) dSIM(nHour,nVAR)
!    Real(8) CO2_ISO_inp  !!cumulative CO2 isotope (e.g., C14 or C13) at the beginning of the time-step
    REAL(8) vENZ(2)      !!specific ENZ activity, ligninase & cellulase
    TYPE(sMEND_OUT),intent(in):: sOUT
    
    INTEGER j,k
    do j = 1,nVAR
        k = VARopt_int(j,1)
        select case (k)  !!see "MEND.ini"
            case (1)    !!CO2
                dSIM(iHour,j) = sOUT%CFLUX%CO2_gm
            case (2)    !!CO2_ISO,e.g., C14_CO2
                dSIM(iHour,j) = sOUT%CO2_gm_iso
!                dSIM(iHour,j) = sOUT%CPOOLI(2)%CO2 - CO2_ISO_inp!!sINP%CPOOLI(2)%CO2
            case (3)    !!MBC
                dSIM(iHour,j) = sOUT%CPOOL%MBC
            case (4)    !!MBC_ISO
                dSIM(iHour,j) = sOUT%CPOOLI(2)%MBC
            case (5)    !!DOC
                dSIM(iHour,j) = sOUT%CPOOL%DOC
            case (6)    !!DOC_SIO
                dSIM(iHour,j) = sOUT%CPOOLI(2)%DOC
            case (7)    !!SOC
                dSIM(iHour,j) = sOUT%CPOOL%SOC
            case (8)    !!SOC_ISO
                dSIM(iHour,j) = sOUT%CPOOLI(2)%SOC
            case (9)    !!ENZC_LIG, ENZ Activity = ENZ concentration*specific enzyme activity
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(1)  !!*vENZ(1)  sPAR%VdPOC(1)
            case (10)   !!ENZC_CEL, ENZ Activity
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(2)  !!*vENZ(2)  !!sPAR%VdPOC(2)
            case (11)    !!CO2_TOT
                dSIM(iHour,j) = sOUT%CPOOL%CO2
            case (12)    !!CO2_TOT: 13C or 14C
                dSIM(iHour,j) = sOUT%CPOOLI(2)%CO2
            case (13)    !!CO2_TOT: 12C
                dSIM(iHour,j) = sOUT%CPOOLI(1)%CO2
            case (14)    !!RR_LIG, response ratio of Geochip oxidative gene abundance
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(1)  !!*vENZ(1)  sPAR%VdPOC(1)
            case (15)   !!RR_LIG, response ratio of Geochip hydrolytic gene abundance
                dSIM(iHour,j) = sOUT%CPOOL%ENZP(2)  !!*vENZ(2)  !!sPAR%VdPOC(2)    
            case default
                dSIM(iHour,j) = sOUT%CFLUX%CO2_gm
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
    
    call sDate2YMD(sDate_beg,iyr0,imo0,ida0)
    
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
    
    call sDate2YMD(sDate_beg,iyr0,imo0,ida0)
    
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
    INTEGER iyr0, imo0, ida0, iyr, imo
    INTEGER nda_in_ym
    INTEGER i,j 
    INTEGER nhr,nda,nmo, ibeg, iend,nbe
    INTEGER iRead
    CHARACTER(LEN=10) sDateHr
    CHARACTER(len=8) sDate
    CHARACTER(len=6) sYM
    CHARACTER(len=2000)sRead1,sRead2
    CHARACTER(len=200)format1
    REAL(8), DIMENSION(24,nVAR)::dSIM_h1
    REAL(8), DIMENSION(744,nVAR)::dSIM_h2 !!hours in 1 month, at most =24*31
    REAL(8), DIMENSION(nVAR)   ::dSIM_t
    dSIM_h1 = const_FillValue
    dSIM_h2 = const_FillValue
    dSIM_t = const_FillValue
    
!    write(format1,*)"(i10,",nVAR,"(e20.6))"
    write(format1,*)"(A10,",nVAR,"(e20.6))"
    open(2,file=sFile_t,status='unknown')
    
    nda = nDaysbwDates(sDate_beg,sDate_end)
    nmo = nMonsbwDates(sDate_beg,sDate_end)
!    nhr = nda*24
!    if(nh.ne.nhr) then
!        write(*,*)"Error: nRow<>nHour in 'sOUT_OPT'"
!    end if
!    
    call sDate2YMD(sDate_beg,iyr0,imo0,ida0)
    open(1,file=sFile_hour,status='old')
    read(1,'(a)')sRead1      !1st line: time period
    write(2,'(a)')sRead1
    read(1,'(a10,a)')sRead1,sRead2  !2nd line: head
    if(tstep.eq.1) then
        write(2,'(a10,a)')"Day",sRead2
    elseif(tstep.eq.2) then
        write(2,'(a10,a)')"Mon",sRead2
    end if
    do i=1,nRow_skip-2
        read(1,'(a)')sRead1
        write(2,'(a)')sRead1
    end do
    
    select case (tstep)
        case (1) !!daily
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
        case (2) !!monthly
            if(nmo > 0) then
                iyr = iyr0
                imo = imo0
                nda = 0
                do i=1,nmo
                    nda_in_ym = min(nda, nDaysofMon(iyr,imo))  !!10/2/2019: temporarily used when nmo = 1 and nda = a few days within a month; need further revision if any dates are used for beginning and ending 
    !                write(*,*)'iyr,imo,nday = ',iyr,imo,nda_in_ym
                    ibeg = 1
                    iend = nda_in_ym*24
                    do j = ibeg,iend
                        read(1,*)sDateHr,dSIM_h2(j,1:nVAR)
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
!        case(4) !!yearly
            
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
    print*,">>>Convert OUTPUTs from HOURLY to DAILY & MONTHLY: BEG"
    
    print*,">>>[1] STATE VARIABLEs:"
    sFile_inp = trim(dirout)//"VAR_hour.out"   
    sFile_out = trim(dirout)//"VAR_day.out"
    nRow_skip=3; nVAR=const_nVARc*(const_nISO+1); tstep=1; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_out = trim(dirout)//"VAR_mon.out"
    nRow_skip=3; nVAR=const_nVARc*(const_nISO+1); tstep=2; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    call system('gzip -f '//sFile_inp)
    
    print*,">>>[2] FLUXes:"
    sFile_inp = trim(dirout)//"FLX_hour.out"
    sFile_out = trim(dirout)//"FLX_day.out"
    nRow_skip=2; nVAR=const_nFLXc; tstep=1; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_out = trim(dirout)//"FLX_mon.out"
    nRow_skip=2; nVAR=const_nFLXc; tstep=2; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    call system('gzip -f '//sFile_inp)
    
    print*,">>>[3] PARAMETERs:"
    sFile_inp = trim(dirout)//"PAR_hour.out"
    sFile_out = trim(dirout)//"PAR_day.out"
    nRow_skip=2; nVAR=33; tstep=1; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_out = trim(dirout)//"PAR_mon.out"
    nRow_skip=2; nVAR=33; tstep=2; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    call system('gzip -f '//sFile_inp)
    
    print*,">>>[4] DERIVED RATEs:"
    sFile_inp = trim(dirout)//"RATE_hour.out"
    sFile_out = trim(dirout)//"RATE_day.out"
    nRow_skip=2; nVAR=18; tstep=1; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    
    sFile_out = trim(dirout)//"RATE_mon.out"
    nRow_skip=2; nVAR=18; tstep=2; flag_avg=1
    call sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sDate_beg, sDate_end,tstep,flag_avg)
    call system('gzip -f '//sFile_inp)
    print*,">>>Convert OUTPUTs from HOURLY to DAILY & MONTHLY: END"
    print*,"" 
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
    CHARACTER(len=20) sRead(10),units
!    REAL(8) fSWPsat1,fSWP_dec1
    REAL(8) SWCres,SWCsat,alpha,rn
    REAL(8) SWC,SWP
!    fSWPsat1 = fSWPsat(20d0,30d0)
!    fSWP_dec1 = fSWP_dec(-1d0,fSWPsat1)
    
!    print*, fSWP_OPT(-dexp(3.5*dlog(10d0)))
!    print*, fSWP_OPT(-dexp(0.0*dlog(10d0)))
!    print*, fSWP_OPT(-dexp(-0.5*dlog(10d0)))
!    print*, fSWP_OPT(-dexp(-1*dlog(10d0)))
    
    SWCres = 0.162603606
    SWCsat = 0.579564599
    alpha  = 0.021
    rn     = 1.5630346
    units = "perc"
    SWC = 0.195
    SWP = fSWC2SWP(SWC,SWCres,SWCsat,alpha,rn,-1d3)
    write(*,*)SWC,SWP
    units = "MPa"
    SWC = fSWP2SWC(SWP,units,SWCres,SWCsat,alpha,rn)
    write(*,*)SWC,SWP

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
!    sPAR % VdPOC = (/xx(1), xx(2)/) ![mg POC/mg ENZP/h], maximum reaction rate for conversion of POC by ENZP
!    sPAR % KsPOC = (/xx(3), xx(4)/) ![mg POC/cm3], half-saturation constant for conversion of POC by ENZP
!    sPAR % rENZP = (/xx(5), xx(5)/) ![1/h],turnover rate of ENZP
!    sPAR % pENZP = xx(6) ![mg ENZP/mg MBC/h], production rate of ENZP
!    sPAR % frPOC2DOC = xx(7) ![-], fraction of decomposed POC allocated to DOC
!    sPAR % frMBC2DOC = xx(8) ![-], fraction of dead MBC allocated to DOC
!    sPAR % Vg = xx(9) ![mg DOC/mg MBC/h], maximum uptake rate of DOC by MB
!    sPAR % KsDOC = xx(10) ![mg DOC/cm3], half-saturation constant for uptake of DOC by MB
!    sPAR % Yg = xx(11) ![-], carbon use efficiency in uptake of DOC by MB
!    sPAR % Vm = xx(12) ![1/h], specific microbial maintenance rate
!    sPAR % rENZM = xx(13) ![1/h], turnover rate of ENZMAOC
!    sPAR % pENZM = xx(14) ![mg ENZM/mg MBC/h], production rate of ENZMAOC
!    sPAR % VdMOC = xx(15) ![mg MOC/mg ENZMAOC/h], maximum reaction rate for conversion of MAOC by ENZMAOC
!    sPAR % KsMOC = xx(16) ![mg MOC/cm3], half-saturation constant for conversion of MAOC by ENZMAOC
!    sPAR % Qmax = xx(17) ![mg C/g soil], adsorption capacity
!    sPAR % Kba = xx(18) ![mg C/g soil/h], binding affinity
!    sPAR % Kdes = xx(19) ![mg DOC/h],desorption rate constant
!    sPAR % Kads = sPAR % Kdes * sPAR % Kba ![mg DOC/mg DOC/h], adsorption rate constant = Kba*Kdes
!
!
!    !     print*, sINI%nHour
!    call subMEND_RUN(xx, sPAR, sINI, sOUT)
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

