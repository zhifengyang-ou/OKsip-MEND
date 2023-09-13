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
    TYPE(sSCE_PAR) sSCE


    INTEGER i, j, nRun, eof  !!, iDay
    INTEGER iFpar, iFpar_UQ,iFvar !output file for model comparison: Sim vs. Obs; input file for parameter samples; output file for response variables  

    real(8) fObj,fObj0
    real(8) fObj_cr(2)  !!upper and lower bound for COFI
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

    call MENDIN(sSCE,sINI) !read Model parameters: initial value, lower and upper bounds 
    print*,"nCase=",sINI%nCase
    Allocate(xx(sSCE%nPar))
    xx = sSCE%a
    
    !Select Model Run (0: run model; 1: optimization; 2: uncertainty)   
    !!sINI%iModel = 1!see options above, moved to SUBROUTINE MENDIN()  
    SELECT CASE (sINI%iModel)
    CASE (1) !SCEUA optimization
        write(sSCE%iFout_end,*)"SCE-UA Results from Multiple Runs:"
        sRead = "    OBJ-1:"
        write(format510,*)"(/,",sSCE%nPar,"(a12),","' |  CRITERION'",",a10, I10)"                
        write(sSCE%iFout_end,format510)sSCE%parName,sRead,sINI%nVARopt  !!(1:sSCE%nPar)  !(xname(j),j=1,nopt2)
        
        if (sSCE%nRun .gt. 0) then
            nRun = min(sSCE%nRun, 200)
        else
            nRun = 1
        end if
        
        print*, ">>MODEL CALIBRATION/OPTIMIZATION..."
!        write(*,*)"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
!        write(*,*)"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
!        write(*,'(a10,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
        write(*,'(a10,I5)')"nRun= ", nRun
        write(*,'(a10,I5)')"nPar= ", sSCE%nPar
        write(*,'(a10,I5)')"nOpt= ", sSCE%nOpt
        write(format510,*)"(a12,",sSCE%nOpt,"(I3))"
        write(*,format510)"PAR_Opt[i]= ",sSCE%iOpt
        write(*,*)
        
        ALLOCATE(iwk(nRun))
        ALLOCATE(bestObj(nRun))
        ALLOCATE(bestPar(nRun, sSCE%nPar))

        do i = 1, nRun
            call iRandSeedGen(sSCE%iseed)  !!see MOD_USRFS
            call srand(sSCE%iseed)  !!reinitialize the random number generator
            write (*, '(A20,I5,A15,I15)') 'SCEUA Run Number = ', i, '; Random Seed = ', sSCE%iseed

            call SCEUA(sSCE, sPAR, sINI, sOUT)
            bestObj(i) = sSCE % bestObj
            bestPar(i,:) = sSCE % bestPar
        end do

        call indexx(nRun, bestObj, iwk)     !rank best objective function value (OBF)
        xx = bestPar(iwk(1),:)              !pick parameter values resulting in best OBF

    CASE (2,3) !read parameter samples to compute statistics of response variables
        write(*,*)">>UNCERTAINTY QUANTIFICATION, MULTIPLE PAR VALUES ARE PROVIDED in <UQpar.dat>"
        iFpar    = 301            !File to store parameter values for sensitivity/uncertainty analysis
        iFpar_UQ = 302            !PAR file with fObj < fObj_cr
        iFvar    = 303            !File to store response variable outputs from sensitivity/uncertainty analysis
        sINI%iFout_UQvar = iFvar
        
        open(unit = iFpar, file = trim(sINI%dirinp)//'/COFIpar.dat', status = 'old')
        open(unit = iFpar_UQ, file = trim(sINI%dirout)//'COFIpar.out', status = 'unknown')
        open(unit = iFvar, file = trim(sINI%dirout)//'COFIvar.out', status = 'unknown')
        
!        write(format101, *) "(", sINI%nOBS_tot, "E15.3)"
                        
        read(iFpar, *)sRead,fObj_cr(1:2)  !!Critical Objective Function Value
        write(iFpar_UQ,'(A20,2f10.4)')sRead,fObj_cr(1:2)
        write(*,'(A20,2f10.4)')sRead,fObj_cr(1:2)
        read(iFpar, '(A)')sRead      !!field names
        write(iFpar_UQ,*)sRead
        !! fOBJ            LF0             r0           fINP             VP 
        j = 0
        DO
            read(iFpar, *, iostat = eof) fObj0,xx(1:sSCE%nPar)
            
            IF (eof < 0) THEN !end of file has reached
                EXIT
            ELSE IF (eof > 0) THEN !input error
                print*,">>>Skip this Line!"
            ELSE
                IF(fObj0.le.fObj_cr(1).and.fObj0.gt.fObj_cr(2)) then  !!fObj = rRead(1)
                    j = j+1
                    fObj = fObj0
                    if (sINI%iModel.eq.3) then
                        fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
!                        write(iFvar, format101)sINI%dSIM_opt(:,2)  !! write output in fMEND_OBJ
                    end if
                    write(*,'(A5,I6,A1,2f20.4)')"fObj[",j,"]",fObj0,fObj
                    write(iFpar_UQ, '(f10.4,50f15.8)')fObj,xx
!                ELSE
!                    print*,">>>Skip this Line!"
                END IF
            END IF
        END DO

        close(iFpar)
        close(iFpar_UQ)
        close(iFvar)
    END SELECT !!CASE (sINI%iModel)
    
!RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "SCEIN.DAT" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
    write (*,*)
    write (*,*) '>>FINAL RUN <fMEND_OBJ> with GIVEN or BEST PAR...'

    do i=1,sSCE%nPar
        write(*,'(I3,a10,f12.6)')i,sSCE%parName(i),xx(i)
    end do
    
!    write(*,*)"Input-Date_Period = ",sINI%sDate_beg_all," -- ",sINI%sDate_end_all
!    write(*,*)"Simulation_Period = ",sINI%sDate_beg_sim," -- ",sINI%sDate_end_sim
!    write(*,'(a10,I5)')"nMon= ", nMonsbwDates(sINI%sDate_beg_sim, sINI%sDate_end_sim)
    write(*,'(/,A)')">>>MEND is RUNNING, PLEASE BE PATIENT..."

    sINI % iModel = 0 !RUN MEND MODEL USING (i)PARAMETER VALUES IN LAST LINE OF "SCEIN.DAT" or (ii)"BEST" PARAMETER VALUES FROM OPTIMIZATION
    fObj = fMEND_OBJ(xx, sPAR, sINI, sOUT)
        

    
!    write(sINI%iFout_SIM_obs,'(/,a)')"PARAMETERS & CRITERIA:"
!    sRead = "    OBJ-1:"
!    write(format510,*)"(/,",sSCE%nPar,"(a12),","' |  CRITERION'",",a10, I10)"                
!    write(sINI%iFout_SIM_obs,format510)sSCE%parName,sRead,sINI%nVARopt
!    write(format521,*)"(",sSCE%nPar,"f12.6,","' | '",",f10.4,",sINI%nVARopt,"f10.4)"
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
    
    close(sSCE%iFout_all)
    close(sSCE%iFout_end)
    close(sSCE%iFout_ini)
    
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

