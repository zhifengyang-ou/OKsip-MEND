SUBROUTINE MENDIN_CASE(iCase,sINI)
!c$debug
!c
!c   THIS SUBROUTINE READS 
!c    (1) MEND MODEL INITILIZATION
!c    (2) INPUT DATA/CALIBRATION DATA
!c    (3) INPUT VARIABLES FOR MODEL OPTIMIZATION
!c    (4) MODEL PARAMETERS 
!c  
!c  AUTHOR: GANGSHENG WANG
!C  Environmental Sciences Division
!C  Oak Ridge National Laboratory
!C  Oak Ridge, TN 37831-6301
!C  March, 2013 
!C  Updated: May 5, 2015

      USE MOD_OPT_TYPE
      USE MOD_MEND_TYPE
      USE MOD_MEND  , only: fSWC2SWP,sINP_Read,sOUT_ALL_tscale!!sOUT_OPT  !!function
      USE MOD_String, only: StrCompress
      USE MOD_USRFS,  only: nDaysbwDates,nMonsbwDates,nDaysofMon
      USE MOD_USRFS,  only: sDate_After,sInt2Str
      USE MOD_USRFS,  only: Array_Normalize !!function
      
      IMPLICIT NONE
      !!ARGUMENTS:
!      TYPE(sSCE_PAR), intent(inout):: sPAR_SCE
      INTEGER,        intent(in)   :: iCase
      TYPE(sMEND_INI),intent(inout):: sINI
      
      !!LOCAL VARIABLES:
      INTEGER iFin  !input file unit
      
!      real pcenta
!      integer, allocatable:: iPar(:), iPar_opt(:)
!      character*10 pcntrl,deflt,usrsp
!      character*4 reduc,initl,ysflg,noflg
      integer i, j, k, lp !!, ierror, iwarn
      integer eof !!end of file
      INTEGER is_total  !!=1: need to convert to hourly rate; =0: directly assign value

!      data deflt/' DEFAULT  '/
!      data usrsp/'USER SPEC.'/
!      data ysflg/'YES '/
!      data noflg/'NO  '/
      
      character(len=8) sDate_beg,sDate_end,sDate
      integer ifdata,nfile, ndays,nmons
      integer iyr,imo, nVAR
      CHARACTER(LEN=2) str2
      character(len=10) sUnits,ststep
      character(len=20) sfilename(10)
      character(len=200) sfilename_full,dirinp, dirout,sFile_inp,sFile_out
      character(len=50) sRead,sRead1
      character(len=200) sRead2
      real(8) rRead
      INTEGER iRead
      real(8) vg_SWCres,vg_SWCsat,vg_alpha,vg_n!!van-Genuchten equation
      
      character(len=20)   :: Name_POOL(const_nVARc)

      DATA Name_POOL /"SOC","POC1","POC2","MOC","QOC","DOC","MBC","MBCA","MBCD",&
                                "ENZ_POC1","ENZ_POC2","ENZ_MOC","CO2_TOT"/
      !      write (*,*) '>>ENTER SUBROUTINE <MENDIN_CASE>'

!c  INITIALIZE I/O VARIABLES
      iFin = 10
      dirinp = trim(sINI%dirinp)//"/"//trim(sINI%CASEname(iCase))//"/"
      sINI%dirinp_case = dirinp  !!used by subMEND_INI()
      
      sfilename_full = trim(dirinp)//trim(sINI%CASEname(iCase))//".ini"
      open(unit=iFin,file=sfilename_full,status='old')
      
      read(iFin,*)sRead                 !!Line 1
      read(iFin,*)sRead                 !!Line 2
      read(iFin,*)sRead                 !!Line 3
      read(iFin,*)sINI%SITE             !!Line 4
      read(iFin,*)sRead                 !!Line 5
      read(iFin,*)sINI%BIOME            !!Line 6
      read(iFin,*)sRead                 !!Line 7      
      read(iFin,*)sINI%SOM              !!Line 8
      read(iFin,*)sRead                 !!Line 9
      read(iFin,*)sRead                 !!Line 10
      read(iFin,*)sINI%sDate_beg_all,sINI%sDate_end_all !!Line 11
      read(iFin,*)sRead                 !!Line 12
      read(iFin,*)sINI%sDate_beg_sim,sINI%sDate_end_sim !!Line 13
      
      ndays = nDaysbwDates(sINI%sDate_beg_all,sINI%sDate_beg_sim)
      if(ndays.lt.1) then
          print*,"ERROR: Date_beg_sim = ",sINI%sDate_beg_sim, " < Date_beg_all = ",sINI%sDate_beg_all
          stop
      end if
      
      if(sINI%iModel.eq.1) then !!optimization
        ndays = nDaysbwDates(sINI%sDate_end_all,sINI%sDate_end_sim)
        if(ndays.gt.0) then
            sINI%sDate_end_sim = sINI%sDate_end_all
!        else
!            write(*,*)"Ignore the WARNING, Simulation Period is OK (within the Input-Data Period)"
        end if
      end if
           
      sINI%iFout_SIM_obs    = iCase*100 + 11
      sINI%iFout_SIM_day    = iCase*100 + 12
      sINI%iFout_SIM_mon    = iCase*100 + 13
      sINI%iFout_VAR_hour   = iCase*100 + 14
      sINI%iFout_FLX_hour   = iCase*100 + 15
      sINI%iFout_RATE_hour  = iCase*100 + 16
      sINI%iFout_PAR_hour   = iCase*100 + 17
      sINI%iFout_ITW_hour   = iCase*100 + 18
           
      if(sINI%iModel.eq.0) then
      !!output files for response variables
      !!--------------------------------------------------------------------------
      dirout = trim(sINI%dirout)//"/"//trim(sINI%CASEname(iCase))//"/"//trim(sINI%SITE)//"_"
      sfilename_full = trim(dirout)//'SIM_obs.out'
      open(unit = sINI%iFout_SIM_obs, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_SIM_obs,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_SIM_obs,'(2a5,a15,3a20)')"ID","VAR","Date","OBS_avg","SIM_avg","SIM_sd"
      
      sfilename_full = trim(dirout)//'SIM_day.out'
      open(unit = sINI%iFout_SIM_day, file = sfilename_full, status = 'unknown')  !!daily simulation output; see MOD_MEND::subMEND_RUN() 
      write(sINI%iFout_SIM_day,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_SIM_day,'(a10,a50)')"Day","AVG & STDDEV for VAR[1:n]"
      
      sfilename_full = trim(dirout)//'SIM_mon.out'
      open(unit = sINI%iFout_SIM_mon, file = sfilename_full, status = 'unknown')  !!daily simulation output; see MOD_MEND::subMEND_RUN() 
      write(sINI%iFout_SIM_mon,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_SIM_mon,'(a10,a50)')"Mon","AVG & STDDEV for VAR[1:n]"
      
      sfilename_full = trim(dirout)//"VAR_hour.out"
      open(unit = sINI%iFout_VAR_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_VAR_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_VAR_hour,'(a10,100a20)')"Hour",(trim(Name_POOL(i)),i=1,const_nVARc),&
                                                    (trim(Name_POOL(i))//"_I1",i=1,const_nVARc),&
                                                    (trim(Name_POOL(i))//"_I2",i=1,const_nVARc)

      sfilename_full = trim(dirout)//"FLX_hour.out"
      open(unit = sINI%iFout_FLX_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_FLX_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_FLX_hour,'(a10,50a20)')"Hour","POCdec1","POCdec2","POCdec2DOC1","POCdec2DOC2",&
            "POCdec2MOC1","POCdec2MOC2","MOC2DOC","QOC2DOC","DOC2QOC","DOC2QOCnet","DOC2MBC",&
            "CO2_growth","CO2_maintn","CO2_maintn_dorm","CO2_gm","MBC_mortality", "MBC2EP1","MBC2EP2","MBC2EM","MBC_PM",&
            "EP2DOC1","EP2DOC2","EM2DOC","MBC2DOC","MBC2POC","MBA2MBD","MBD2MBA"

      sfilename_full = trim(dirout)//"RATE_hour.out"
      open(unit = sINI%iFout_rate_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_RATE_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_RATE_hour,'(a10,20a20)')"Hour","kPOC1","kPOC2","kMOC","kDOC","kMBa","kMBa_in","kMBd",&
                                               "kMBd_in","kMB","kMB_in","phi","Active_Fraction","CUE",&
                                               "Balance_Error","TOCbeg","TOCend","TOCinp","TOCout"
      
      sfilename_full = trim(dirout)//"PAR_hour.out"
      open(unit = sINI%iFout_par_hour, file = sfilename_full, status = 'unknown')
      write(sINI%iFout_PAR_hour,*)"Simulation_Period = ",sINI%sDate_beg_sim, " -- ",sINI%sDate_end_sim
      write(sINI%iFout_PAR_hour,'(a10,33a20)')"Hour","VP1","VP2","VM","KP1","KP2","KM","Qmax","Kba","Kdes","Kads", &
                                              "rEP1","rEP2","rEM","pEP","pEM","fD","gD",&
                                              "Vg","alpha","Vm","KD","Yg","wdie","gamma","rMORT",&
                                              "beta","Vm_dorm","VmA2D","VmD2A","SWP_A2D","tau","SWP_D2A","wdorm"
      
      sfilename_full = trim(dirout)//trim("ITW_hour.dat")
      open(unit = sINI%iFout_ITW_hour, file = sfilename_full, status = "unknown")
      write(sINI%iFout_ITW_hour,*)"Data_Period = ",sINI%sDate_beg_all, " -- ",sINI%sDate_end_all
      write(sINI%iFout_ITW_hour,'(a10,5a20)')"Hour","SIN_mg/cm3/h","STP_oC","SWC","SWP_MPa","pH"
      !!--------------------------------------------------------------------------
      end if !!if(sINI%iModel.eq.0) 
      
      
      nmons = nMonsbwDates(sINI%sDate_beg_all,sINI%sDate_end_all)
      ndays = nDaysbwDates(sINI%sDate_beg_all,sINI%sDate_end_all)
      sINI%nHour = ndays*24
      sINI%nHour_sim = 24*nDaysbwDates(sINI%sDate_beg_sim,sINI%sDate_end_sim)
      ALLOCATE(sINI%STP(sINI%nHour))
      ALLOCATE(sINI%SWC(sINI%nHour))
      ALLOCATE(sINI%SWP(sINI%nHour))
      ALLOCATE(sINI%SpH(sINI%nHour))
      ALLOCATE(sINI%SIN(sINI%nHour)) !!SOC input
      
      !!INPUT INFO:
      read(iFin,*)sRead         !!Line 14
      read(iFin,*)sRead         !!Line 15
      !!Soil temperature info: Hourly or Daily
      read(iFin,*)sRead         !!Line 16
      read(iFin,*)ifdata,sUnits,ststep,nfile   !!Line-17
      read(iFin,*)sRead         !!Line 18
      if(ifdata.eq.1) then  !!read hourly data
          read(iFin,*)sfilename(1:nfile)     !!Line-19

          sUnits = StrCompress(sUnits)
          ststep = StrCompress(ststep)  !!TODO NEXT: convert data with time-step (ststep=monthly,daily) to hourly
          is_total = 0
          call sINP_Read(nfile,sfilename,dirinp,ststep,is_total,nmons,sINI%nHour,sINI%STP)
      else  !!constant temperature
          read(iFin,*)rRead  !![oC],Line-19
          sINI%STP = rRead    
      end if
      
      !!Soil water info: Hourly or Daily
      read(iFin,*)sRead                         !!Line 20
      read(iFin,*)ifdata,sUnits,ststep,nfile    !!Line-21
      read(iFin,*)sRead                         !!Line 22
      if(ifdata.eq.1) then  !!read hourly data
          read(iFin,*)sfilename(1:nfile)        !!Line-23
          read(iFin,*)sRead                     !!Line-24
          read(iFin,*)vg_SWCres,vg_SWCsat,vg_alpha,vg_n  !!Line-25: van-Genuchten equation

          sUnits = StrCompress(sUnits)
          ststep = StrCompress(ststep)
          is_total = 0
          call sINP_Read(nfile,sfilename,dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SWC)
          
          if (trim(sUnits).eq."perc") then  !!%, percent
            sINI%SWC = sINI%SWC/1.d2
          end if
          
          do k = 1, sINI%nHour
              sINI%SWP(k) = fSWC2SWP(sINI%SWC(k),vg_SWCres,vg_SWCsat,vg_alpha,vg_n,const_SWPmin)
          end do
          
      else  !!constant water potential (MPa), negative value
          read(iFin,*)rRead          !!Line-23
          sINI%SWC = const_FillValue  !!we only use SWP
          sINI%SWP = rRead
          !!skip Line-24:25
          read(iFin,*)sRead
          read(iFin,*)sRead
      end if
      
      !!External input, e.g., litter fall, fertilizer: Daily or Monthly
      read(iFin,*)sRead                         !!Line-26
      read(iFin,*)ifdata,sUnits,ststep,nfile    !!Line-27
      read(iFin,*)sRead                         !!Line-28
      if(ifdata.eq.1) then  !!nfile = 1 if ststep='monthly'
          read(iFin,*)sfilename(1:nfile)        !!Line-29

          sUnits = StrCompress(sUnits)
          ststep = StrCompress(ststep)  !!convert data with time-step (ststep=monthly,daily) to hourly
          is_total = 1  !!usually litter_input is the total amount during a period
          call sINP_Read(nfile,sfilename,dirinp,ststep,is_total,nmons,sINI%nHour,sINI%SIN)
      
      else !!constant litter input [mgC/cm3/h]
          read(iFin,*)rRead  !!Line-29
          sINI%SIN = rRead
      end if
      
      read(iFin,*)sRead                 !!Line-30
      read(iFin,*)sINI%SIN_frac(1:3)    !!Line-31, fraction of input to POC1, POC2, & DOC
      read(iFin,*)sRead                 !!Line-32
      read(iFin,*)sINI%SIN_other(1,1:3) !!Line-33: [mgC/cm3/h], other constant inputs, e.g., coarse wood, roots
      read(iFin,*)sRead                 !!Line-34
      read(iFin,*)sINI%SIN_other(2,1:3),sINI%sDate_beg_inp2,sINI%sDate_end_inp2  !!Line-35: [mgC/cm3/h], other constant inputs, e.g., coarse wood, roots
      read(iFin,*)sRead                 !!Line-36
      read(iFin,*)sINI%SIN_C12_C14, sINI%SIN_Multiplier    !!Line-37: ratio of C12 to C14 in SOC input; multiplier for litter input during post-data period 
      sUnits = StrCompress(sUnits)
      ststep = StrCompress(ststep)
      sINI%SIN_frac(3) = 1.d0 - sINI%SIN_frac(1) - sINI%SIN_frac(2) !!ensure total = 100%
      
      !!soil pH
      read(iFin,*)sRead                         !!Line-38
      read(iFin,*)ifdata,sUnits,ststep,nfile    !!Line-39
      read(iFin,*)sRead                         !!Line-40
      read(iFin,*)rRead                         !!Line-41
      sUnits = StrCompress(sUnits)
      ststep = StrCompress(ststep)
      if(trim(ststep).eq."const") then
          sINI%SpH(1:sINI%nHour) = rRead
!          write(*,*)"pH = ", rRead
      end if
      
      if(sINI%iModel.eq.0) then
        !!write hourly input data into 1 file
        do k = 1,ndays
            call sDate_After(k,sINI%sDate_beg_all,sDate)
            do j=1,24 !!sINI%nHour
                call sInt2Str(j,2,str2)
                i = (k - 1)*24 + j
                write(sINI%iFout_ITW_hour,'(A10,5f20.6)')sDate//str2,sINI%SIN(i),sINI%STP(i),sINI%SWC(i),sINI%SWP(i),sINI%SpH(i)
            end do
        end do
        close(sINI%iFout_ITW_hour)

        !!compute daily/monthly STP,SWP,SIN-----------------------------------------------------------BEG
        !    print*,">>>Inputs, Temperature, Water Content & Potential:"
          sFile_inp = trim(dirout)//"ITW_hour.dat"
          sFile_out = trim(dirout)//"ITW_day.dat"
  !        sOUT_ALL_tscale(sFile_inp,sFile_out,nRow_skip,nVAR, sINI%sDate_beg_all, sINI%sDate_beg_end,tstep,flag_avg)
          call sOUT_ALL_tscale(sFile_inp,sFile_out,2,5, sINI%sDate_beg_all, sINI%sDate_end_all,1,1)

          sFile_out = trim(dirout)//"ITW_mon.dat"
          call sOUT_ALL_tscale(sFile_inp,sFile_out,2,5, sINI%sDate_beg_all, sINI%sDate_end_all,2,1)
        !!compute daily/monthly STP,SWP,SIN-----------------------------------------------------------END  
      end if !!if(sINI%iModel.eq.0) 
      
      read(iFin,*)sRead                         !!Line-42
      read(iFin,*)sINI%SOIL_INI_file            !!Line-43: soil initialization file     
     
      !!Calibration Variables & Data
      read(iFin,*)sRead                         !!Line-44
      read(iFin,*)sRead                         !!Line-45
      read(iFin,*)sRead, nVAR                   !!Line-46 !!# of output variables  !!sINI%nObs_var
      read(iFin,'(a)')sRead2                    !!Line-47
!      write(sPAR_SCE%iFout3,*)
!      write(sPAR_SCE%iFout3,*)"MODEL CALIBRATION VARIABLES & DATA-FILE (see MEND.ini):"
!      write(sPAR_SCE%iFout3,*)sRead2
!      write(sPAR_SCE%iFout3,'(100(1h-))')
      ALLOCATE(sINI%VARopt(nVAR))
      ALLOCATE(sINI%VARstep(nVAR))
      ALLOCATE(sINI%VARfile(nVAR))
      ALLOCATE(sINI%VARcol(nVAR))
      ALLOCATE(sINI%VARobj(nVAR))
      ALLOCATE(sINI%VARobjw(nVAR))
      sINI%nVARopt = 0
      !!Line-48:57
      do i=1,nVAR
          read(iFin,*)iRead,sRead1,sRead2,sINI%VARopt(i),sINI%VARstep(i),sINI%VARfile(i),&
                                        sINI%VARcol(i),sINI%VARobj(i),sINI%VARobjw(i)
!          write(*,*)sINI%VARopt(i),sINI%VARstep(i),sINI%VARfile(i),&
!                                        sINI%VARcol(i),sINI%VARobj(i),sINI%VARobjw(i)
          if(sINI%VARopt(i).gt.0) then
              sINI%nVARopt = sINI%nVARopt + 1
!              write(sPAR_SCE%iFout3,'(I5,5x,a10,a15,2I5,5x,a15,I5,a10,f10.0)') &
!                                      iRead,sRead1,sRead2,sINI%VARopt(i),sINI%VARstep(i),sINI%VARfile(i),&
!                                      sINI%VARcol(i),sINI%VARobj(i),sINI%VARobjw(i)
!!1 CO2 mgC-cm3-h   1 	1   HR.obs 	2   NSEC	10
          end if
      end do
      if(sINI%nVARopt.lt.1) then
          write(*,*)'No Variables Available for Model Optimization!!!'
      end if
      close(iFin)
      
      ALLOCATE(sINI%VARopt_int(sINI%nVARopt, 3))
      ALLOCATE(sINI%rOBJ(sINI%nVARopt))
      ALLOCATE(sINI%rOBJw(sINI%nVARopt))
      sINI%VARopt_int = 0  !!initialization

      j=0
      do i=1,nVAR
          if(sINI%VARopt(i).gt.0) then
              j=j+1
              sINI%VARopt_int(j,1) = i  !!index of par_opt
          end if
      end do

      do i=1,sINI%nVARopt
          j = sINI%VARopt_int(i,1)
          sINI%rOBJw(i) = sINI%VARobjw(j)
          sINI%VARopt_int(i,3) = sINI%VARstep(j)  !!time-step
          sfilename_full = trim(dirinp)//trim(sINI%VARfile(j))
          open(101,file=sfilename_full,status='old')
          read(101,*)sRead,iRead    !!# of observations
          sINI%VARopt_int(i,2) = iRead
          close(101)
!          print*, sINI%VARopt_int(i,:)
      end do
      
      sINI%nOBS_tot = sum(sINI%VARopt_int(1:sINI%nVARopt,2))
      ALLOCATE(sINI%dOBS_opt(sINI%nOBS_tot,3)) !!date,obs,iVARopt
      ALLOCATE(sINI%dSIM_opt(sINI%nOBS_tot,3)) !!date,sim,sim_sd
      sINI%dOBS_opt = const_FillValue
      sINI%dSIM_opt = const_FillValue
      call Array_Normalize(sINI%nVARopt,sINI%rOBJw,const_FillValue)
      
!      write (*,*) '>>>EXIT SUBROUTINE <MENDIN_CASE>'

      return
END !!SUBROUTINE MENDIN

