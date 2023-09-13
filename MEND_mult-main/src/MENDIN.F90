SUBROUTINE MENDIN(sPAR_SCE,sINI)
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
      TYPE(sSCE_PAR), intent(inout):: sPAR_SCE
      TYPE(sMEND_INI),intent(inout):: sINI
      
      !!LOCAL VARIABLES:
      INTEGER iFin  !input file unit
      
      real pcenta
      integer, allocatable:: iPar(:), iPar_opt(:)
      character*10 pcntrl,deflt,usrsp
      character*4 reduc,initl,ysflg,noflg
      integer i, j, k, lp, ierror, iwarn
      integer eof !!end of file
      INTEGER is_total  !!=1: need to convert to hourly rate; =0: directly assign value

      data deflt/' DEFAULT  '/
      data usrsp/'USER SPEC.'/
      data ysflg/'YES '/
      data noflg/'NO  '/
      
      character(len=8) sDate_beg,sDate_end,sDate
      integer ifdata,nfile, ndays,nmons
      integer iyr,imo, nVAR
      CHARACTER(LEN=2) str2
      character(len=10) sUnits,ststep
      character(len=20) sfilename(10)
      character(len=200) sfilename_full,sFile_inp,sFile_out
      character(len=50) sRead,sRead1
      character(len=200) sRead2
      real(8) rRead
      INTEGER iRead
      real(8) vg_SWCres,vg_SWCsat,vg_alpha,vg_n!!van-Genuchten equation

      write (*,*) '>>ENTER SUBROUTINE <MEND_INI>'

!c  INITIALIZE I/O VARIABLES
      iFin = 10
      open(unit=iFin,file='MEND.ini',status='old')

      ierror = 0
      iwarn = 0

!c  READ THE SCE CONTROL PARAMETERS
      sPAR_SCE%ideflt = 0
      !!BEGIN to read 'MEND.ini'
      read(iFin,*)sRead                 !!Line 1
      read(iFin,*)sRead                 !!Line 2 
      read(iFin,*)sRead                 !!Line 3
      read(iFin,*)sINI%iModel           !!Line 4
      read(iFin,*)sRead                 !!Line 5
      read(iFin,*)sRead                 !!Line 6
      read(iFin,'(a)')sRead             !!Line 7
      sINI%dirinp = StrCompress(sRead)
      read(iFin,*)sRead                 !!Line 8
      read(iFin,'(a)')sRead             !!Line 9
      sINI%dirout = StrCompress(sRead)
      read(iFin,*)sRead                 !!Line 10
      read(iFin,*)sRead                 !!Line 11
      read(iFin,*)sRead, sINI%nCase     !!Line 12
      read(iFin,*)sRead                 !!Line 13
      ALLOCATE(sINI%CASEname(sINI%nCase)) 
      read(iFin,*)sINI%CASEname(1:sINI%nCase)   !!Line 14  
      !!-------------------------------------------------
      read(iFin,*)sRead                 !!Line 15
      read(iFin,*)sRead                 !!Line 16: "SCE Parameters"
      read(iFin,*)sRead                 !!Line 17
      !!Line 18:
      read(iFin,*)sPAR_SCE%maxn,sPAR_SCE%kstop,sPAR_SCE%pcento,sPAR_SCE%ngs,sPAR_SCE%nRun,sPAR_SCE%ideflt 

!c  IF ideflt IS EQUAL TO 1, READ THE SCE CONTROL PARAMETERS
      read(iFin,*)sRead                 !!Line 19
      !!Line 20:
      if (sPAR_SCE%ideflt .eq. 1) Then
        read(iFin,*)sPAR_SCE%npg,sPAR_SCE%nps,sPAR_SCE%nspl,sPAR_SCE%mings,sPAR_SCE%iniflg,sPAR_SCE%iprint
        pcntrl = usrsp
      else
        read(iFin,*)
        pcntrl = deflt
      end if
!      write (*,*)sRead
!      write (*,810)sPAR_SCE%npg,sPAR_SCE%nps,sPAR_SCE%nspl,sPAR_SCE%mings,sPAR_SCE%iniflg,sPAR_SCE%iprint
  800 format(2i5,f10.4,3i5)
  810 format(6i5)

!! READ Decomposition Kinetics: sPAR%iKinetics: 0-Michaelis-Menten, 1-First Order, 2-Second Order
      read(iFin,*)sRead                 !!Line 21: dash line
      read(iFin,*)sRead                 !!Line-22, Comments
      read(iFin,*)sINI%iKinetics        !!Line-23
      
!! READ OBJECTIVE FUNCTION WEIGHTING FACTORS: VARobjw_cases[1:10]
      read(iFin,*)sRead                 !!Line 24: dash line
      read(iFin,*)sRead                 !!Line-25, Comments
      read(iFin,*)sINI%VARobjw_cases(1:10) !!Line-26
      read(iFin,*)sINI%VARobj_cases(1:10)
      
!c  READ THE INITIAL PARAMETER VALUES AND THE PARAMETER BOUNDS
      read(iFin,*)sRead                 !!Line 27
      read(iFin,*)sRead                 !!Line 28: model parameters
      read(iFin,*)sRead, sPAR_SCE%nPar  !!Line 29
      sINI%nPar = sPAR_SCE%nPar
!      write(*,'(a10,I5)')sRead, sPAR_SCE%nPar
      
      read(iFin,*)sRead                 !!Line 30:column names
      
      ALLOCATE(iPar(sPAR_SCE%nPar))
      ALLOCATE(iPar_opt(sPAR_SCE%nPar))
      ALLOCATE(sPAR_SCE%parName(sPAR_SCE%nPar))
      ALLOCATE(sPAR_SCE%a(sPAR_SCE%nPar))
      ALLOCATE(sPAR_SCE%bl(sPAR_SCE%nPar))
      ALLOCATE(sPAR_SCE%bu(sPAR_SCE%nPar))
      ALLOCATE(sPAR_SCE%bestPar(sPAR_SCE%nPar))
      
      sPAR_SCE%nOpt = 0
      !!Line 31-57
      do j = 1, sPAR_SCE%nPar
        read(iFin,*)iPar(j), sPAR_SCE%parName(j), sPAR_SCE%a(j), sPAR_SCE%bl(j), sPAR_SCE%bu(j), iPar_opt(j)
!        print*, sPAR_SCE%parName(j)
        sPAR_SCE%parName(j) = ADJUSTR(sPAR_SCE%parName(j))
!        print*, sPAR_SCE%parName(j)
!        print*, iPar(j), sPAR_SCE%parName(j), sPAR_SCE%a(j), sPAR_SCE%bl(j), sPAR_SCE%bu(j), iPar_opt(j)
        
        if(iPar_opt(j).gt.0) then
            sPAR_SCE%nOpt = sPAR_SCE%nOpt + 1
        end if
      end do
  830 format(3f10.3)
  
      read(iFin,*)              !!Line 58:"------"
      read(iFin, *)             !!Line 59: parameter names
      read(iFin,*)              !!Line 60:"------"
      read(iFin,*)sPAR_SCE%a    !!Line 61
 
      close(iFin)
      
      
      ALLOCATE(sPAR_SCE%iOpt(sPAR_SCE%nOpt))
      i = 0
      do j = 1, sPAR_SCE%nPar
          if(iPar_opt(j).gt.0) then
              i = i + 1
              sPAR_SCE%iOpt(i) = j
          end if
      end do
!      print*, sPAR_SCE%iOpt

!c  IF ideflt IS EQUAL TO 0, SET THE SCE CONTROL PARAMETERS TO
!c  THE DEFAULT VALUES
      if (sPAR_SCE%ideflt .eq. 0) then
        sPAR_SCE%npg = 2*sPAR_SCE%nOpt + 1
        sPAR_SCE%nps = sPAR_SCE%nOpt + 1
        sPAR_SCE%nspl = sPAR_SCE%npg
        sPAR_SCE%mings = sPAR_SCE%ngs
        sPAR_SCE%iniflg = 0
        sPAR_SCE%iprint = 0
      end if

 
      !!define file units
      sPAR_SCE%iFout_all = 11
      sPAR_SCE%iFout_end = 12
      sPAR_SCE%iFout_ini = 13
      sINI%iFout_SIM_obs_cases    = 14
      
      !!3 output files for model optimization
      !!--------------------------------------------------------------------------
!      sINI%dirout = trim(sINI%dirout)//"/"//trim(sINI%SITE)//"_"
      CALL system('mkdir '//sINI%dirout)
      sINI%dirout = trim(sINI%dirout)//"/"
      sfilename_full = trim(sINI%dirout)//"OPT_all.out"
      open(unit=sPAR_SCE%iFout_all,file=sfilename_full,status='unknown')
      sfilename_full = trim(sINI%dirout)//"OPT_end.out"
      open(unit=sPAR_SCE%iFout_end,file=sfilename_full,status='unknown')
      sfilename_full = trim(sINI%dirout)//"OPT_ini.out"
      open(unit=sPAR_SCE%iFout_ini,file=sfilename_full,status='unknown')
      !!output files for response variables
      !!--------------------------------------------------------------------------
!      sfilename_full = trim(sINI%dirout)//'SIM_obs.out'
!      open(unit = sINI%iFout_SIM_obs_cases, file = sfilename_full, status = 'unknown')
!      write(sINI%iFout_SIM_obs_cases,*)"SIMULATION vs. OBSERVATION for ALL CASES: "
!      write(sINI%iFout_SIM_obs_cases,'(2a5,a15,3a20)')"ID","VAR","Date","OBS_avg","SIM_avg","SIM_sd"
      !!--------------------------------------------------------------------------
      write(sPAR_SCE%iFout_ini,700)
  700 format(10x,'SHUFFLED COMPLEX EVOLUTION GLOBAL OPTIMIZATION',&
     &       /,10x,46(1h=))
     
      do i=1,sINI%nCase
          sRead = "mkdir "//trim(sINI%dirout)//trim(sINI%CASEname(i))
          call system(sRead)
      end do
      !!--------------------------------------------------------------------------
      
!c  CHECK IF THE SCE CONTROL PARAMETERS ARE VALID
      if (sPAR_SCE%ngs .lt. 1 .or. sPAR_SCE%ngs .ge. 1320) then
        write(sPAR_SCE%iFout_ini,900) sPAR_SCE%ngs
  900   format(//,1x,'**ERROR** NUMBER OF COMPLEXES IN INITIAL ',&
     &         ' POPULATION ',i5,' IS NOT A VALID CHOICE')
        ierror = ierror + 1
      end if

      if (sPAR_SCE%kstop .lt. 0 .or. sPAR_SCE%kstop .ge. 20) then
        write(sPAR_SCE%iFout_ini,901) sPAR_SCE%kstop
  901   format(//,1x,'**WARNING** THE NUMBER OF SHUFFLING LOOPS IN',&
     &  ' WHICH THE CRITERION VALUE MUST CHANGE ',/,13x,'SHOULD BE',&
     &  ' GREATER THAN 0 AND LESS THAN 10.  ','kstop = ',i2,&
     &  ' WAS SPECIFIED.'/,13x,'BUT kstop = 5 WILL BE USED INSTEAD.')
        iwarn = iwarn + 1
        sPAR_SCE%kstop=5
      end if

      if (sPAR_SCE%mings .lt. 1 .or. sPAR_SCE%mings .gt. sPAR_SCE%ngs) then
        write(sPAR_SCE%iFout_ini,902) sPAR_SCE%mings
  902   format(//,1x,'**WARNING** THE MINIMUM NUMBER OF COMPLEXES ',&
     &         i2,' IS NOT A VALID CHOICE. SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%mings = sPAR_SCE%ngs
      end if

      if (sPAR_SCE%npg .lt. 2 .or. sPAR_SCE%npg .gt. 1320/max(sPAR_SCE%ngs,1)) then
        write(sPAR_SCE%iFout_ini,903) sPAR_SCE%npg
  903   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A COMPLEX ',&
     &         I4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%npg = 2*sPAR_SCE%nOpt+1
      end if

      if (sPAR_SCE%nps.lt.2 .or. sPAR_SCE%nps.gt.sPAR_SCE%npg .or. sPAR_SCE%nps.gt.50) then
        write(sPAR_SCE%iFout_ini,904) sPAR_SCE%nps
  904   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN A SUB-',&
     &  'COMPLEX ',i4,' IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%nps = sPAR_SCE%nOpt + 1
      end if

      if (sPAR_SCE%nspl .lt. 1) then
        write(sPAR_SCE%iFout_ini,905) sPAR_SCE%nspl
  905   format(//,1x,'**WARNING** THE NUMBER OF EVOLUTION STEPS ',&
     &         'TAKEN IN EACH COMPLEX BEFORE SHUFFLING ',I4,/,13x,&
     &         'IS NOT A VALID CHOICE, SET IT TO DEFAULT')
        iwarn = iwarn + 1
        sPAR_SCE%nspl = sPAR_SCE%npg
      end if

!c  COMPUTE THE TOTAL NUMBER OF POINTS IN INITIAL POPULATION
      sPAR_SCE%npt = sPAR_SCE%ngs * sPAR_SCE%npg

      if (sPAR_SCE%npt .gt. 1320) then
        write(sPAR_SCE%iFout_ini,906) sPAR_SCE%npt
  906   format(//,1x,'**WARNING** THE NUMBER OF POINTS IN INITIAL ',&
     &         'POPULATION ',i5,' EXCEED THE POPULATION LIMIT,',/,13x,&
     &         'SET NGS TO 2, AND NPG, NPS AND NSPL TO DEFAULTS')
        iwarn = iwarn + 1
        sPAR_SCE%ngs = 2
        sPAR_SCE%npg = 2*sPAR_SCE%nOpt + 1
        sPAR_SCE%nps = sPAR_SCE%nOpt + 1
        sPAR_SCE%nspl = sPAR_SCE%npg
      end if

!c  PRINT OUT THE TOTAL NUMBER OF ERROR AND WARNING MESSAGES
      if (ierror .ge. 1) write(sPAR_SCE%iFout_all,907) ierror
  907 format(//,1x,'*** TOTAL NUMBER OF ERROR MESSAGES IS ',i2)

      if (iwarn .ge. 1) write(sPAR_SCE%iFout_all,908) iwarn
  908 format(//,1x,'*** TOTAL NUMBER OF WARNING MESSAGES IS ',i2)

      if (sPAR_SCE%mings .lt. sPAR_SCE%ngs) then
        reduc = ysflg
      else
        reduc = noflg
      end if

      if (sPAR_SCE%iniflg .ne. 0) then
        initl = ysflg
      else
        initl = noflg
      end if


!c  PRINT SHUFFLED COMPLEX EVOLUTION OPTIMIZATION OPTIONS
  104 write(sPAR_SCE%iFout_ini,910)
  910 format(//,2x,'SCE CONTROL',5x,'MAX TRIALS',5x,&
     &'REQUIRED IMPROVEMENT',5x,'RANDOM',/,3x,'PARAMETER',8x,&
     &'ALLOWED',6x,'PERCENT',4x,'NO. LOOPS',6x,'SEED',/,&
     &2x,11(1h-),5x,10(1H-),5x,7(1h-),4x,9(1h-),5x,6(1h-))

      pcenta=sPAR_SCE%pcento*100.
      write(sPAR_SCE%iFout_ini,912) pcntrl,sPAR_SCE%maxn,pcenta,sPAR_SCE%kstop,sPAR_SCE%iseed
  912 format(3x,a10,7x,i5,10x,f3.1,9x,i2,9x,i5)
      write(sPAR_SCE%iFout_ini,914) sPAR_SCE%ngs,sPAR_SCE%npg,sPAR_SCE%npt,sPAR_SCE%nps,sPAR_SCE%nspl
  914 format(//,18x,'SCE ALGORITHM CONTROL PARAMETERS',/,18x,32(1H=),&
     &//,2x,'NUMBER OF',5x,'POINTS PER',5x,'POINTS IN',6x,'POINTS PER',&
     &4x,'EVOL. STEPS',/,2x,'COMPLEXES',6X,'COMPLEX',6x,'INI. POPUL.',&
     &5x,'SUB-COMPLX',4x,'PER COMPLEX',/,2x,9(1h-),5x,10(1h-),4x,&
     &11(1h-),5x,10(1h-),4x,11(1h-),5x,/,2x,5(i5,10x))
      write(sPAR_SCE%iFout_ini,915) reduc,sPAR_SCE%mings,initl
  915 format(//,15x,'COMPLX NO.',5x,'MIN COMPLEX',5x,'INI. POINT',/,&
     &15x,'REDUCTION',6x,'NO. ALLOWED',6x,'INCLUDED',/,&
     &15x,10(1h-),5x,11(1h-),5x,10(1h-),/,18x,a4,6x,i8,13x,a4)
      write(sPAR_SCE%iFout_ini,916)
  916 format(//,8x,'INITIAL PARAMETER VALUES AND PARAMETER BOUNDS',/,&
     &       8x,45(1h=),//,2x,'PARAMETER',5x,'INITIAL VALUE',5x,&
     &       'LOWER BOUND',5x,'UPPER BOUND',6x,'OPT-Y/N',/,2x,9(1h-),5x,13(1h-),5x,&
     &       11(1h-),5x,11(1h-),5x,11(1h-))
     
      do 920 i = 1, sPAR_SCE%nPar
        write(sPAR_SCE%iFout_ini,918) sPAR_SCE%parName(i),sPAR_SCE%a(i),sPAR_SCE%bl(i),sPAR_SCE%bu(i),iPar_opt(i)
  920 continue
  918   format(a10,4x,(3x,f12.6),2(4x,f12.6),5x,I10)
  
      if (ierror .ge. 1) then
      write(sPAR_SCE%iFout_ini,922)
  922 format(//,'*** THE OPTIMIZATION SEARCH IS NOT CONDUCTED BECAUSE',&
     &       ' OF INPUT DATA ERROR ***')
      stop
      end if
      
      
      write (*,*) '>>>EXIT SUBROUTINE <MEND_INI>'
      return
END !!SUBROUTINE MENDIN

