PROGRAM MEND_main
!MEND: Microbial-ENzyme Decomposition model
!C  Author: GANGSHENG WANG @ ORNL
!C  Environmental Sciences Division
!C  Oak Ridge National Laboratory
!C  Oak Ridge, TN 37831-6301
!C  EMAIL: WANGG@ORNL.GOV
!C  March, 2014
!C  Updated: August 28, 2015

USE MOD_MEND_TYPE
USE MOD_MEND, ONLY: fMEND_OBJ, sOUT_tscale, TEST 
USE MOD_OPT_TYPE
USE MOD_OPT, ONLY: SCEUA
USE MOD_MCMC, ONLY: MCMC
USE MOD_USRFS, ONLY: iRandSeedGen, indexx, gasdev,gasdev0
USE MOD_USRFS, ONLY: Sec2HMS, nMonsbwDates

    IMPLICIT NONE

!    TYPE(sMEND_INP) sINP
    TYPE(sMEND_PAR) sPAR
    TYPE(sMEND_OUT) sOUT
    TYPE(sMEND_INI) sINI
    TYPE(sSCE_PAR) sPAR_SCE


    INTEGER i, j, nRun, eof  !!, iDay
    INTEGER iFpar,iFpar_UQ, iFvar, iFobj, nrow, ncol !output file for model comparison: Sim vs. Obs; input file for parameter samples; output file for response variables  

    INTEGER iRead(3)
    real(8) rRead, fObj, fObj_cr
    real(8), ALLOCATABLE :: xx(:), bestPar(:,:), bestObj(:)
    integer, ALLOCATABLE :: iwk(:)

    CHARACTER(LEN = 1000) sRead 
!    CHARACTER(LEN = 20) sRead1(5)
    CHARACTER(LEN = 200) format100, format101,format510,format521
       
    INTEGER t_start,t_end,t_rate,t_elapse  !!time_start, time_end, time_elapsed
    INTEGER, DIMENSION(3):: tHMS  !!CALL subroutine "Sec2HMS" in "WGSfunc.f"
    INTEGER time_end(8)

!!--------------------------------------------------------- 
!! TEST: BEGIN 
!!    CALL TEST()
!! TEST: END
!!--------------------------------------------------------- 
   
    print*, "*---------------------------------------------*"
    print*, "* Microbial-ENzyme Decomposition (MEND) Model *"
    print*, "*-------Carbon-Nitrogen Coupled Version-------*"
    print*, "*---------wangg@ornl.gov; Aug 1, 2015---------*"
    print*, "*---------------------------------------------*"
    
!    write (*,*) '>>MEND RUN BEGINS:'
    CALL system_clock(t_start,t_rate)
    
    sINI%nOutStep = 24    ![h], output interval, 24 h = daily

    CALL MENDIN(sPAR_SCE,sINI) !read Model parameters: initial value, lower and upper bounds 
    Allocate(xx(sPAR_SCE%nPar))
    xx = sPAR_SCE%a
    
    !Select Model Run (0: run model; 1: optimization; 2: uncertainty)    
    SELECT CASE (sINI%iModel)
    CASE (1) !SCEUA optimization
        write(sPAR_SCE%iFout_end,*)"SCE-UA Results from Multiple Runs:"
        sRead = "    OBJ-1:"
        write(format510,*)"(/,",sPAR_SCE%nPar,"(a16),","' |  CRITERION'",",a10, I10)"                
        write(sPAR_SCE%iFout_end,format510)sPAR_SCE%parName,sRead,sINI%nVARopt  
        
        if (sPAR_SCE%nRun .gt. 0) then
            nRun = min(sPAR_SCE%nRun, 200)
        else
            nRun = 1
        end if
        
        print*, ">>MODEL CALIBRATION/OPTIMIZATION..."
        write(*,*)"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
        write(*,*)"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
        write(*,'(a10,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
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
            CALL iRandSeedGen(sPAR_SCE%iseed)  !!see MOD_USRFS
!!            sPAR_SCE%iseed = -1919936880  !! for DEBUG: see MOD_USRFS.F90, very small variance<0 in fSTDDEV (caused by identical xx), bug fixed on 11/16/2018
!            sPAR_SCE%iseed = -1869941160  !!for DEBUG: Active Microbial Biomass (MBA <0), see MOD_MEND.F90 Line 1097, bug fixed on 11/20/2018
            CALL srand(sPAR_SCE%iseed)  !!reinitialize the random number generator
            write (*, '(A20,I5,A15,I15)') 'SCEUA Run Number = ', i, '; Random Seed = ', sPAR_SCE%iseed
            
            !!wgs: test random numbers
!            do j=1,1000
!                print*,j,gasdev0(sPAR_SCE%iseed),gasdev(sPAR_SCE%iseed)
!            end do
            
            CALL SCEUA(sPAR_SCE, sPAR, sINI, sOUT)
            bestObj(i) = sPAR_SCE % bestObj
            bestPar(i,:) = sPAR_SCE % bestPar
        end do

        CALL indexx(nRun, bestObj, iwk)     !rank best objective function value (OBF)
        xx = bestPar(iwk(1),:)              !pick parameter values resulting in best OBF
        
        DEALLOCATE(bestPar)
        DEALLOCATE(bestObj)
        DEALLOCATE(iwk)
        
    CASE (2) !read parameter from SCE samples to compute statistics of response variables
        write(*,*)">>UNCERTAINTY QUANTIFICATION, MULTIPLE PAR VALUES ARE PROVIDED in <UQpar.dat>"
        iFpar    = 301            !File to store parameter values for sensitivity/uncertainty analysis
        iFpar_UQ = 302            !PAR file with fObj < fObj_cr
        iFvar    = 303            !File to store response variable outputs from sensitivity/uncertainty analysis
        
        open(unit = iFpar, file = trim(sINI%dirinp)//'/UQpar.dat', status = 'old')
        open(unit = iFpar_UQ, file = trim(sINI%dirout)//'UQpar.out', status = 'unknown')
        open(unit = iFvar, file = trim(sINI%dirout)//'UQvar.out', status = 'unknown')
        write(format101, *) "(", sINI%nOBS_tot, "E15.3)"
        read(iFpar, *)sRead,fObj_cr  !!Critical Objective Function Value
        write(iFpar_UQ,'(A20,f10.4)')sRead,fObj_cr
        write(*,'(A20,f10.4)')sRead,fObj_cr
        read(iFpar, '(A)')sRead      !!field names
        write(iFpar_UQ,*)sRead
        !! LOOP TRIALS COMPLXS  BESTF   WORSTF   PAR-RNG             LF0  
        j = 0
        DO
            read(iFpar, *, iostat = eof) rRead,xx(1:sPAR_SCE%nPar)
            
            IF (eof < 0) THEN !end of file has reached
                EXIT
            ELSE IF (eof > 0) THEN !input error
                print*,">>>Skip this Line!"
            ELSE
                IF(rRead.le.fObj_cr) then  !!fObj = rRead(1)
                    j = j+1
    !                print*, sRead(1:3), rRead,sRead(4:5),xx
                    fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
                    write(*,'(A5,I6,A1,2f20.4)')"fObj[",j,"]",rRead,fObj
                    write(iFpar_UQ, '(f10.4,50f15.8)')rRead,xx
                    write(iFvar, format101)sINI%dSIM_opt(:,2)  !!(sINI%dSim_comp2(i,:), i = 1, sINI % nObs_var)
!                ELSE
!                    print*,">>>Skip this Line!"
                END IF
            END IF
        END DO

        close(iFpar)
        close(iFpar_UQ)
        close(iFvar)
    CASE (3) !read parameter samples to compute statistics of response variables
        write(*,*)">>SENSITIVITY ANALYSIS, MULTIPLE PAR SETS ARE PROVIDED in <SApar.txt>"
        iFpar    = 301            !File to store parameter values for sensitivity/uncertainty analysis
        iFobj    = 302            !File to store fObj
        iFvar    = 303            !File to store response variable outputs from sensitivity/uncertainty analysis
        
        open(unit = iFpar, file = trim(sINI%dirinp)//'/SApar.txt', status = 'old')
!        open(unit = iFpar, file = trim(sINI%dirinp)//'/SApar_less.dat', status = 'old')
        open(unit = iFobj, file = trim(sINI%dirout)//'SAobj.out', status = 'unknown')
        open(unit = iFvar, file = trim(sINI%dirout)//'SAvar.out', status = 'unknown')
        write(format101, *) "(", sINI%nOBS_tot, "E15.3)"
        read(iFpar, '(A)')sRead  !!Par names
!        print*, sRead
        j = 0
        DO
            read(iFpar, *, iostat = eof) xx(1:sPAR_SCE%nPar)            
            IF (eof < 0) THEN !end of file has reached
                print*,'eof=',eof
                EXIT
            ELSE IF (eof > 0) THEN !input error
                print*,">>>Skip this Line!"
            ELSE
                j = j+1
                if(j.ge.sINI%iSA_range(1) .and. j.le.sINI%iSA_range(2)) then
                    fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
                    write(*,'(A5,I6,A1,f20.4)')"fObj[",j,"]",fObj
                    write(iFobj, *)fobj
                    write(iFvar, format101)sINI%dSIM_opt(:,2)  !!(sINI%dSim_comp2(i,:), i = 1, sINI % nObs_var)
                end if
            END IF
        END DO

        close(iFpar)
        close(iFobj)
        close(iFvar)
    CASE (4) ! MCMC
!        write(sPAR_SCE%iFout_end,*)"SCE-UA Results from Multiple Runs:"
        CALL iRandSeedGen(sPAR_SCE%iseed)  !!see MOD_USRFS
        CALL srand(sPAR_SCE%iseed)  !!reinitialize the random number generator
        write (*, '(A20,I5,A15,I15)') 'SCEUA Run Number = ', i, '; Random Seed = ', sPAR_SCE%iseed
        
        CALL MCMC(sPAR_SCE, sPAR, sINI, sOUT)
        xx = sPAR_SCE % bestPar
    END SELECT !!CASE (sINI%iModel)
    
!RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "SCEIN.DAT" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
    IF(sINI % iModel.lt.2 .or.sINI % iModel.eq.4) then
        write (*,*)
        write (*,*) '>>FINAL RUN <fMEND_OBJ> with GIVEN or BEST PAR...'
        ncol = 11
        nrow = sPAR_SCE%nPar/ncol
        if (nrow.lt.1) then
            write(*,'(20a12)')sPAR_SCE%parName
            write(*,'(20f12.6)')xx
        else
            do i=1,nrow
                write(*,'(20a12)')sPAR_SCE%parName(((i-1)*ncol+1):(i*ncol))
                write(*,'(20f12.6)')xx(((i-1)*ncol+1):(i*ncol))
            end do
            
            if(sPAR_SCE%nPar.gt.(nrow*ncol)) then
                write(*,'(20a12)')sPAR_SCE%parName((nrow*ncol+1):sPAR_SCE%nPar)
                write(*,'(20f12.6)')xx((nrow*ncol+1):sPAR_SCE%nPar)
            end if
        end if
        
        
!        do i=1,sPAR_SCE%nPar
!            write(*,'(I3,a10,f12.6)')i,sPAR_SCE%parName(i),xx(i)
!        end do

        write(*,'(/,4A)')"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
        write(*,'(4A)')"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
        write(*,'(A,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)

        if(sINI%Carbon_only) then
            write(*,*)">>>MEND is RUNNING, PLEASE BE PATIENT: Carbon-Only"  
        else
            write(*,*)">>>MEND is RUNNING, PLEASE BE PATIENT: Carbon-Nitrogen Coupled"  
        end if

        sINI % iModel = 0 !RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "MEND_namelist.nml" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
        fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)

        !!write SIM vs. OBS
        Print*,">>>Write output for SIM vs. OBS: "
        Print*,"ATTENTION: It does NOT matter if there are void rOBJ values due to NO OBS data available."
        do i=1,sINI%nOBS_tot 
            write(sINI%iFout_SIM_obs,'(2i5,i15,3f20.8)')i,int(sINI%dOBS_opt(i,3)),int(sINI%dOBS_opt(i,1)),&
                                sINI%dOBS_opt(i,2),sINI%dSIM_opt(i,2),sINI%dSIM_opt(i,3)   !!obs_mean,sim_mean,sim_sd
        end do

        write(sINI%iFout_SIM_obs,'(/,a)')"OBJ-FUNCTIONS & PARAMETERS:"            
        write(sINI%iFout_SIM_obs,'(a20,f10.4)')"fOBJ (best=0) = ", fObj
        write(sINI%iFout_SIM_obs,'(a20,6f10.4)')"fOBJ[i] = ",sINI%rOBJ
        write(sINI%iFout_SIM_obs,'(a20,6f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw

        write(format510,*)"(/,",sPAR_SCE%nPar,"(a16))"                
        write(sINI%iFout_SIM_obs,format510)sPAR_SCE%parName

        write(format521,*)"(",sPAR_SCE%nPar,"(f15.8,','))"
!        write(format521,*)"(",sPAR_SCE%nPar,"(f30.20,','))"
        write(sINI%iFout_SIM_obs,format521) xx

        !!close all output files
        close(sINI%iFout_SIM_obs)
        close(sINI%iFout_SIM_day)
        close(sINI%iFout_SIM_mon)
        close(sINI%iFout_VAR_hour)
        close(sINI%iFout_FLX_hour)
        close(sINI%iFout_RATE_hour)
        close(sINI%iFout_PAR_hour)
        
        close(sINI%iFout_mcmc)

        !!Write to Screen:
        write(*,'(a20,f10.4)')"fOBJ (best=0) = ", fObj
        write(*,'(a20,6f10.4)')"fOBJ[i] = ",sINI%rOBJ
        write(*,'(a20,6f10.4)')"fOBJ_weight[i] = ",sINI%rOBJw
        
        !!Convert OUTPUTS from HOURLY to DAILY & MONTHLY for all STATE VARIABLEs | FLUXes | RATEs
        CALL sOUT_tscale(sINI%dirout,sINI%sDate_beg_sim,sINI%sDate_end_sim)
        
    END IF !!IF(sINI % iModel.lt.2)
    
    close(sPAR_SCE%iFout_all)
    close(sPAR_SCE%iFout_end)
    close(sPAR_SCE%iFout_ini)


    DEALLOCATE(xx)  
    DEALLOCATE(sINI%STP)
    DEALLOCATE(sINI%SWC)
    DEALLOCATE(sINI%SWP)
    DEALLOCATE(sINI%SpH)
    DEALLOCATE(sINI%SIN)
    
    IF(.NOT.sINI%Carbon_only) THEN
        DEALLOCATE(sINI%SIN_NH4)
        DEALLOCATE(sINI%SIN_NO3)
    END IF
    
    CALL system_clock(t_end) !!timer(t_end)
    t_elapse = (t_end - t_start)/real(t_rate)
    CALL Sec2HMS(t_elapse,tHMS)
    CALL DATE_AND_TIME(values = time_end) ! Get the current time

    write(*,'(a16,3(I3,a8))')">>Elapsed Time = ",tHMS(1),"Hours", tHMS(2),"Minutes",tHMS(3),"Seconds"
    write (*,'(a25,I4,a1,I2,a1,I2,I3,a1,I2,a1,I2)') ">>MEND RUN COMPLETED at: ", &
            time_end(1),'-',time_end(2),'-',time_end(3),time_end(5),':',time_end(6),':',time_end(7)
    STOP
END PROGRAM MEND_main
!!==================================================================================================================================    

