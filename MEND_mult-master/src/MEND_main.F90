PROGRAM MEND_main
!MEND: Microbial-ENzyme-mediated Decomposition model
!C  Author: GANGSHENG WANG @ ORNL
!C  Environmental Sciences Division
!C  Oak Ridge National Laboratory
!C  Oak Ridge, TN 37831-6301
!C  EMAIL: WANGG@ORNL.GOV
!C  March, 2014
!C  Updated: July 17, 2015

USE MOD_MEND_TYPE
USE MOD_MEND, ONLY: fMEND_OBJ, sOUT_tscale, TEST 
USE MOD_OPT_TYPE
USE MOD_OPT, ONLY: SCEUA
USE MOD_USRFS, ONLY: iRandSeedGen, indexx
USE MOD_USRFS, ONLY: Sec2HMS, nMonsbwDates

    IMPLICIT NONE

    TYPE(sMEND_PAR) sPAR
    TYPE(sMEND_INP) sINP
    TYPE(sMEND_OUT) sOUT
    TYPE(sMEND_INI) sINI
    TYPE(sSCE_PAR) sPAR_SCE


    INTEGER i, j, nRun, eof  !!, iDay
    INTEGER iFpar, iFvar !output file for model comparison: Sim vs. Obs; input file for parameter samples; output file for response variables  

    real(8) fObj
    real(8), allocatable :: xx(:), bestPar(:,:), bestObj(:)
    integer, allocatable :: iwk(:)

    CHARACTER(LEN = 10) soil_select, soil, sRead
    CHARACTER(LEN = 200) format100, format101,format510,format521
       
    INTEGER t_start,t_end,t_rate,t_elapse  !!time_start, time_end, time_elapsed
    INTEGER, DIMENSION(3):: tHMS  !!call subroutine "Sec2HMS" in "WGSfunc.f"

!!--------------------------------------------------------- 
!! TEST: BEGIN 
!    CALL TEST()
!! TEST: END
!!--------------------------------------------------------- 
   
    print*, "*---------------------------------------------*"
    print*, "* Microbial-ENzyme Decomposition (MEND) Model *"
    print*, "*----------- Multiple-Case Version -----------*"
    print*, "*---------wangg@ornl.gov; May 1, 2015---------*"
    print*, "*---------------------------------------------*"
    
    write (*,*) '>>MEND RUN BEGINS:'
    call system_clock(t_start,t_rate)
    
    sINI%nOutStep = 24    ![h], output interval, 24 h = daily

    call MENDIN(sPAR_SCE,sINI) !read Model parameters: initial value, lower and upper bounds 
    print*,"nCase=",sINI%nCase
    Allocate(xx(sPAR_SCE%nPar))
    xx = sPAR_SCE%a
    
    !Select Model Run (0: run model; 1: optimization; 2: uncertainty)   
    !!sINI%iModel = 1!see options above, moved to SUBROUTINE MENDIN()  
    SELECT CASE (sINI%iModel)
    CASE (1) !SCEUA optimization
        write(sPAR_SCE%iFout2,*)"SCE-UA Results from Multiple Runs:"
        sRead = "    OBJ-1:"
        write(format510,*)"(/,",sPAR_SCE%nPar,"(a12),","' |  CRITERION'",",a10, I10)"                
        write(sPAR_SCE%iFout2,format510)sPAR_SCE%parName,sRead,sINI%nVARopt  !!(1:sPAR_SCE%nPar)  !(xname(j),j=1,nopt2)
        
        if (sPAR_SCE%nRun .gt. 0) then
            nRun = min(sPAR_SCE%nRun, 200)
        else
            nRun = 1
        end if
        
        print*, ">>MODEL CALIBRATION/OPTIMIZATION..."
!        write(*,*)"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
!        write(*,*)"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
!        write(*,'(a10,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
        write(*,'(a10,I5)')"nRun= ", nRun
        write(*,'(a10,I5)')"nPar= ", sPAR_SCE%nPar
        write(*,'(a10,I5)')"nOpt= ", sPAR_SCE%nOpt
        write(format510,*)"(a12,",sPAR_SCE%nOpt,"(I3))"
        write(*,format510)"PAR_Opt[i]= ",sPAR_SCE%iOpt
        write(*,*)
        
        ALLOCATE(iwk(nRun))
        ALLOCATE(bestObj(nRun))
        ALLOCATE(bestPar(nRun, sPAR_SCE%nPar))

        do i = 1, nRun
            call iRandSeedGen(sPAR_SCE%iseed)  !!see MOD_USRFS
            call srand(sPAR_SCE%iseed)  !!reinitialize the random number generator
            write (*, '(A20,I5,A15,I15)') 'SCEUA Run Number = ', i, '; Random Seed = ', sPAR_SCE%iseed

            call SCEUA(sPAR_SCE, sPAR, sINI, sOUT)
            bestObj(i) = sPAR_SCE % bestObj
            bestPar(i,:) = sPAR_SCE % bestPar
        end do

        call indexx(nRun, bestObj, iwk)     !rank best objective function value (OBF)
        xx = bestPar(iwk(1),:)              !pick parameter values resulting in best OBF

    CASE (2) !read parameter samples to compute statistics of response variables
        write(*,*)">>SIMPLE UNCERTAINTY ANALYSIS, MULTIPLE PAR VALUES ARE PROVIDED in <MEpar.dat>"
        iFpar = 301             !File to store parameter values for sensitivity/uncertainty analysis
        iFvar = 302            !File to store response variable outputs from sensitivity/uncertainty analysis
        open(unit = iFpar, file = trim(sINI%dirinp)//'/SApar.dat', status = 'old')
        open(unit = iFvar, file = trim(sINI%dirout)//'/SAvar.out', status = 'unknown')
        write(format101, *) "(", sINI%nOBS_tot, "f10.4)"
        read(iFpar, *) !head
        DO
            read(iFpar, *, iostat = eof) soil, fobj, xx
            IF (eof < 0) THEN !end of file has reached
                EXIT
            ELSE IF (eof > 0) THEN !input error
                STOP
            ELSE
                IF (trim(soil) .eq. trim(soil_select)) THEN
                print*, soil
                fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
                write(iFvar, format101)sINI%dSIM_opt(:,2)  !!(sINI%dSim_comp2(i,:), i = 1, sINI % nObs_var)
                END IF
            END IF
        END DO

        close(iFpar)
        close(iFvar)
    END SELECT !!CASE (sINI%iModel)
    
!RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "SCEIN.DAT" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
    write (*,*)
    write (*,*) '>>FINAL RUN <fMEND_OBJ> with GIVEN or BEST PAR...'

    do i=1,sPAR_SCE%nPar
        write(*,'(I3,a10,f12.6)')i,sPAR_SCE%parName(i),xx(i)
    end do
    
!    write(*,*)"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
!    write(*,*)"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
!    write(*,'(a10,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
    write(*,'(/,A)')">>>MEND is RUNNING, PLEASE BE PATIENT..."

    sINI % iModel = 0 !RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "SCEIN.DAT" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
    fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
        

    
!    write(sINI%iFout_SIM_obs,'(/,a)')"PARAMETERS & CRITERIA:"
!    sRead = "    OBJ-1:"
!    write(format510,*)"(/,",sPAR_SCE%nPar,"(a12),","' |  CRITERION'",",a10, I10)"                
!    write(sINI%iFout_SIM_obs,format510)sPAR_SCE%parName,sRead,sINI%nVARopt
!    write(format521,*)"(",sPAR_SCE%nPar,"f12.6,","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
!    write(sINI%iFout_SIM_obs,format521) xx,fObj,sINI%rOBJ
    
    write(*,'(/,A)')">>>FINAL OBF for ALL CASES:"
    write(*,'(a20,f10.4)')"fOBJ (best=0) = ", fObj
    write(*,'(a20,6f10.4)')"fOBJ[i] = ",sINI%rOBJ_cases(1:sINI%nVARopt_cases)
    write(*,'(a20,6f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw_cases(1:sINI%nVARopt_cases)
    
    !!close all output files
!    close(sINI%iFout_SIM_obs)
!    close(sINI%iFout_SIM_day)
!    close(sINI%iFout_SIM_mon)
!    close(sINI%iFout_VAR_hour)
!    close(sINI%iFout_FLX_hour)
!    close(sINI%iFout_RATE_hour)
!    close(sINI%iFout_PAR_hour)
    
    close(sPAR_SCE%iFout1)
    close(sPAR_SCE%iFout2)
    close(sPAR_SCE%iFout3)
    
    !!Convert OUTPUTS from HOURLY to DAILY & MONTHLY for all STATE VARIABLEs | FLUXes | RATEs
!    call sOUT_tscale(sINI%dirout,sINI%sDate_beg_sim,sINI%sDate_end_sim)
      
    call system_clock(t_end)!!timer(t_end)
    t_elapse = (t_end - t_start)/real(t_rate)
    call Sec2HMS(t_elapse,tHMS)
    write(*,'(a16,3(I3,a8))')">>Elapsed Time = ",tHMS(1),"Hours", tHMS(2),"Minutes",tHMS(3),"Seconds"
    write (*,'(a23)') ">>MEND RUN COMPLETED :)"
    STOP
END PROGRAM MEND_main
!!==================================================================================================================================    

