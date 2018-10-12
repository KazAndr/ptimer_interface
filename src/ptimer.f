*=========================================================================
*=========================================================================
* PTIMER v.1.1
* Pulsar timing program        ╧ЁюуЁрььр фы  їЁюэюьхЄЁшЁютрэш  яєы№ёрЁют
* ┬√ўшёыхэшх срЁшЎхэЄЁшўхёъшї ╠╧╚ яєы№ёрЁют ш єЄюўэхэшх ярЁрьхЄЁют
* ╬фшэюўэ√ї ш фтющэ√ї яєы№ёрЁют
*=========================================================================
*=========================================================================
*  PRAO ASC FIAN, 2018
*  
*=========================================================================
* last update 
*
* v 1.1  07.08.2018

     
*=========================================================================
        subroutine gtnarg(n)           ! Linux/Win
* (Returns number of arguments passed to program)
        integer*4 n
*(comment one of the two lines below)
        n=NARGS()-1  ! for Win
c        n=IARGC()    ! for Linux
        return
        end
*=========================================================================
        subroutine gtstr(i,sbuf)       ! Linux/Win
*(Returns i-th parameter of input string of parameters into sbuf)
        character sbuf*60
        integer*4 i,st
*(comment one of the lines)
        call GETARG(i,sbuf,st)       ! for Win
c        call GETARG(i,sbuf)          ! for Linux
        return
        end
*=========================================================================
        subroutine gtdate(ardt)        ! Linux/Win
*(Returns computer's day,month,year into array ardt(3))
	INTEGER*4 ardt(3)
*(comment one of the two below)
         CALL getdat(ardt(3),ardt(2),ardt(1))     ! for Win
c        CALL idate(ardt)                         ! for Linux
        return
        end
*=========================================================================
        subroutine gttime(ardt)        ! Linux/Win
*(Returns computer's hrs,min,secs into array ardt(3))
	INTEGER*4 ardt(3),i100
*(comment one of the two below)
        CALL gettim(ardt(1),ardt(2),ardt(3),i100)     ! for Win
c       CALL itime(ardt)                              ! for Linux
        return
        end

*=========================================================================


        PROGRAM PTIMER

*** (Main program PTIMER for analysis of the slow-down model of
*** (PSR rotation and LSQ fitting of observed data)

        implicit none
        INTEGER*4 NMI,IDEGS,MaxSit,MaxJUMP,MaxCLK 

* (Maximal number of TOAs NMI and degree of freedom IDEGS are here)
        PARAMETER(NMI=9999,MaxSit=100,IDEGS=224,MaxJUMP=41,MaxCLK=9999)

        REAL*4 Y(NMI),U(NMI,2),dDM(NMI),mfunct,age,
     &   FR0D(NMI),ERCN,ERP,ERPT,ERPTT,ERAL,ERDL,ERALM,
     &   ERDLM,EM,SEC,QM,DXPOLE(NMI),DYPOLE(NMI),DUT1(NMI),YMEAN,
     &   mc,minWEG,maxWEG,facWEG,eqdWEG,frqMIN,frqMAX,weg00,wegRMS,
     &   rndm
        REAL*8 JDD(NMI),RM(NMI),WEG(NMI),TDBi(NMI),SCTOPS(NMI),WGM,
     &	 WEGHT,ONE(IDEGS,IDEGS),COEF(IDEGS),ER(IDEGS),RHS0(IDEGS),
     &   DXP,DYP,DUTC1,TDBSAV(NMI),TTSSAV(NMI),TMP,TMP1,
     &   APS,DPS,P0,P1,P2,Pepoch,DM(11),PMRA,PMDEC,DIST,
     &   DISTMIN,DISTMAX,JDstart,JDend,
     &   FIPLR,LPLR,APSR,DPSR,R001,dTHET(3),XTOP,
     &	 YTOP,ZTOP,HEIGH,R0,SIG0,JD,UTC,TDT,TDB,JDB,XX,YY,ZZ,TAD,
     &	 AE,DE,RE,XP,YP,ZP,TA,TOS,TTS,SCTOP,FR,TFR,
     &	 f(11),f0(11),erf(11),APSR0,DPSR0,PMRA0,PMDEC0,
     &   R00,DTB,DNP,DRN,DTI,DT1,DT2,DTIS,
     &   FRPDAT,SIG,SIGu,SWEGH,DISP,RSG(3),RSG0(3),ErSG(3),
     &   Pbdot(3),Pbdot0(3),ErPbdot(3),JDD00,JDDi,FR0D00,
     &   Pbd2(3),Pbd20(3),ErPbd2(3),Pbd3(3),Pbd30(3),ErPbd3(3),
     &    Xd2(3),Xd20(3),ErXd2(3),
     &   SC,TPR,C,SUM,RHS,D0,BUF,BUF405,TABOUT,DNUT,SCS,SCJ,
     &    SSUN,SJUP,SSAT,XN,YN,ZN,ZNn(3),XB,YB,ZB,XR,YR,ZR,RNB,TSHAPS,
     &   TSHAPJ,VPR,SS,Si(3),Si0(3),ErSi(3),RG(3),RG0(3),TSHPSS,
     &   ErRG(3),DNP0k,SpdL,ErX(3),ErT0(3),ErOm(3),Gc3(3),
     &   ErGc3(3),Gc30(3),ScDM,e(3),X(3),Pb(3),T0(3),Om(3),
     &   Om1,RTP,EXC(3),eSQ(3),AF(3),BT(3),
     &   DNP0,PI2,RDA,ABSI,Omdot(3),CHI2,ErPb(3),Ere(3),Pb0(3),RigKT(3),
     &   Ex(3),dthet0(3),Edthet(3),X0(3),Om0(3),e0(3),CW,SW,Cu(3),Su(3),
     &   ArcRd,Sday,Omdot0(3),ErOmdot(3),Xdot(3),Xdot0(3),ErXdot(3),
     &   edot0(3),edot(3),Eredot(3),ShapT,ShapC,RoemEp,DTBj,TDBj,DNPj,
     &   PRLX,ErPrlx,PRLX0,DM0(11),ErDM(11),PR,NCIRC,JDBPOS,Ane(3),
     &   DTPER(3),RTOPO(3),TOAmin,TOAmax,R0N(MaxJUMP),erR0N(MaxJUMP),
     &   R0N0(MaxJUMP),DTOAS,phORB1,phORB2,X1,tstep,Pbar,Pobs,xvel,
     &   cFi,sFi,scl,Vpsr,Vobs,dclk
        REAL*8 XTOPi(MaxSit),YTOPi(MaxSit),ZTOPi(MaxSit),
     &   clkdate(MAXCLK),clockc(MAXCLK)
        real*8   am,  bp,bpp,dr,a0,b0                     ! not yet used
	  INTEGER*4 IH,IMM,IS,IH0,IMM0,IS0,iline,
     &   IO,NumPar,Niter,ISS0,ISS1,funit,fuinc,toaun,
     &   NOPTS,III,IDEGEQ,ISM,K,N,ardt(3),NFORM,
     &   IEIN,IPNT,IER,II,NN,IM,IYE,IE,IREQ,NZP,ICW,ICENT,I,J,
     &   KCNT,J1,I1,IT,KT,NLOOPS,NSITES,NORB,neph,LACK,JDINT,JCORR
        INTEGER*2 IndSIT(NMI),DNJUMP(NMI),NJUMP,NJUMPN,DNPHAS(NMI),
     &   DNPHAS0,IndORD(NMI)
        CHARACTER WW*5,PNM*7,WR(IDEGS)*7,PP*1,SBUF*60,WLET(IDEGS)*6, 
     &   KS*40,WrsA*17,WrsD*17,chmon(12)*3,SiteC(MaxSit)*2,ChSIT*2,
     &   ST80*80,head100*100,filpar*80,fileph*80,fileop*80,filobs*80,
     &   filtoa*80,filutc*80,psrname*9,SiteL(MaxSit)*1,codobs*1,
     &   clsite(MAXCLK)*1,confile*80,filleap*80,PW*1
        LOGICAL SHOWPR,DELMEAN,skipTOA,corrPHAS,skipFIT,skipORB,
     &   SHCOM,AKHILL,BINY,mDADE,mDADEP,mBTP,
     &   JUMPED,mBT,mBTPP,EOP,IWEG,WRITRS,FAKET,GENPER,EPHTM,MKINP,
     &   MTRCR,ORBSH,FREFM,found,clkcor,jdsegm,cde405
        
* (Common blocks for partial derivatives)
	  COMMON/CETBL7a/BINY,mDADE,mBTPP,mDADEP,mBT
          COMMON/CETBL7b/IDEGEQ
          COMMON/CETBL7c/DTB,DTPER,FRPDAT,RE,AE,DE,APSR,DPSR,DTIS,WEGHT,
     &        VPR,AF,BT,RigKT,Cu,Su,ex,Znn,X,Gc3,EsQ,ScDM,Ane,EXC,Pb
* (Common block for matrixes of LSQ)
	  COMMON/CETBL6/C(IDEGS),SUM(IDEGS,IDEGS),RHS(IDEGS),D0,IE(IDEGS)
* (Common blocks used by DE200/LE200)
	  COMMON/CETBL1/EM
	  COMMON/CETBL2/ICW,ICENT,IREQ(13)
	  COMMON/CETBL3/BUF(826)
	  COMMON/CETBL4/TABOUT(6,12),DNUT(4)
	  COMMON/NZP/NZP
	  COMMON/CETBL5/BUF405(1018)

	data chmon/'Jan','Feb','Mar','Apr','May','Jun',
     &             'Jul','Aug','Sep','Oct','Nov','Dec'/
	data WLET/'.Const','.F','.Fdot ',' F2dot',' F3dot',    !   1-5
     &  '.Al','.Dl','.mAl ','.mDl ',' pi',' DM ',              !  6-11
     &  ' x1',' Pb1',' e1',' To1',' OM1',' Gm1',' r1',' s1',   ! 12-19 << 1orb
     &  ' OMd1',' xd1',' Pbd1',' ed1',' thet1',                ! 20-24
     &  ' x2',' Pb2',' e2',' To2',' OM2',' Gm2',' r2',' s2',   ! 25-32 << 2 orb
     &  ' OMd2',' xd2',' Pbd2',' ed2',' thet2',                ! 33-37
     &  ' x3',' Pb3',' e3',' To3',' OM3',' Gm3',' r3',' s3',   ! 38-45 << 3 orb
     &  ' OMd3',' xd3',' Pbd3',' ed3',' thet3',                ! 46-50
     &  'x','x','x','x','x','x','x','x','x',                   ! 51-59
     &  ' F4dot',' F5dot',' F6dot',' F7dot',' F8dot',          ! 60-64
     &  ' Pbdd1',' Xdd1','x','x','x',                          ! 65-69
     &  ' DMd1 ',' DMd2 ',' DMd3 ',' DMd4 ',' DMd5 ',' DMd6 ', ! 70-75
     &   25*' Jump',                                           ! 76-100
     &  ' Pbdd1',' Xdd1',' Wdd1',                              ! 101-103
     &  ' Pbd3','x','x','x','x','x','x','x','x','x',           ! 104-113
     &  ' Pbdd2',' Xdd2',' Wdd2',                              ! 114-116
     &  'x','x','x','x','x','x','x','x','x','x',               ! 117-126
     &  ' Pbdd3',' Xdd3',' Wdd3',                              ! 127-129
     &   95*' Jump'/

* (Some constants)
        PI2  =6.28318530717958648d0
        ArcRd=206264.806247d0
        Sday =86400.0d0
        SpdL =299792458.d0       ! Speed of light, m/s
        TA  = 499.004783703d0    ! Seconds in 1-AU
        TAD = TA/Sday
        SSUN=2953.0D0/Spdl       ! Schwarzschild radius of Sun
        SJUP=2.82D0/Spdl         ! Schwarzschild radius of Jupiter
        SSAT=0.844D0/Spdl        ! Schwarzschild radius of Saturn
        
        NITER = 5                ! Maximal number of iterations (default).
        NJUMP=0                  ! JUMPs counter
        funit=11                 ! free-form format file unit
        fuinc=14                 ! free-form format file unit (included file)
        NFORM=5                  ! format number of TOA data
        norb=0                   ! number of orbits
        
*----------------------------------------------------------------------
* Configuration file path
 
c        confile='./ptimer.cfg'              
        confile='c:/ptimer/ptimer.cfg'       

        
* (Input/Output devices selection and flags setting)
*  -------------------------------------------------
        mDADE  =.FALSE.        ! Damour & Deruelle model (DD)
        mDADEP =.FALSE.        ! Damour & Deruelle model (DD+)
        mBT    =.FALSE.        ! Blandford & Teukolsky + (BT)
        mBTP   =.FALSE.        ! Blandford & Teukolsky + (BT+)
        mBTPP  =.FALSE.        ! Blandford & Teukolsky + (BT++)
        
        BINY  =.FALSE.        ! PSR is binary
        SHCOM =.FALSE.        ! Show comments
        AKHILL=.FALSE.        ! How to reach minima
        SHOWPR=.TRUE.         ! display used PSR parameters
        DELMEAN=.TRUE.        ! delete mean in res.
        JUMPED =.FALSE.       !
        corrPHAS=.FALSE.      ! correction to the pulsar phase
        skipFIT =.FALSE.      !
        skipORB =.FALSE.      !
        EOP   =.TRUE.         ! EOP corrections: dUT1,dX,dY
        IWEG  =.TRUE.         ! weighted fit
        WRITRS=.TRUE.         ! write residuals output to files
        FAKET =.FALSE.        ! generate fake TOAs
	GENPER=.FALSE.        ! generate apparent periods
        EPHTM =.FALSE.        ! if TOA are in TB scale
        MKINP =.FALSE.        ! create input file
        MTRCR =.FALSE.        ! print precise cor.matrix
        ORBSH =.FALSE.        ! show orb.phase in GENPER mode
        FREFM =.FALSE.        ! use free format
        CLKCOR=.FALSE.        ! clock correction to TOAs
        JDSEGM=.false.        ! selected JD
	CDE405=.false.        ! de405 ephemerides used/de200 
	
c=======================================================================
* (Set Constants & Parameters of DE200/LE200)
*
	DO 6 II=1,13
 6	IREQ(II)= 0
* (Select Barycentric Frame of reference)
	ICENT   = 0
* (Earth's Coord. & Velocities to be computed).
	IREQ(3) = 2
* (Jupiter's Coord.)
	IREQ(5) = 1
* (Saturn's Coord.)
	IREQ(6) = 1
        
* (Solar Coord.)
	IREQ(10)= 1
* (Nutation)
	IREQ(13)= 1
	EM	= 81.300587d0
	ICW	= 0
* (init variables)
	do i=1,nmi
	 dxpole(i)=0.0
	 dypole(i)=0.0
         dut1(i)=0.0
	enddo
        
* ( Read string of options )
* ( Functions NARGC(S), GETARG() etc. are not compatible)
* ( for PC/UNIX/VAX. Please modify calls of incompatible)
* ( subroutines - they are all above main source PROGRAM PTIMER)

        call GTNARG(NumPar)
        IF (NumPar.ge.1) then
		do i=1,NumPar
                call GTSTR(i,sbuf)
                IF (SBUF(1:2).eq.'-h') SHCOM =.TRUE.
			  IF (SBUF(1:2).eq.'-d') confile=SBUF(3:60)
                IF (SBUF(1:2).eq.'-z') AKHILL=.TRUE.
                IF (SBUF(1:2).eq.'-c') SHOWPR=.FALSE.
                IF (SBUF(1:2).eq.'-r') DELMEAN=.FALSE.
                IF (SBUF(1:2).eq.'-o') ORBSH =.TRUE. 
			  IF (SBUF(1:2).eq.'-p') THEN
				       PNM(1:7)=SBUF(3:9)
				       ENDIF
                IF (SBUF(1:2).eq.'-i') MKINP =.TRUE.
                IF (SBUF(1:2).eq.'-m') MTRCR =.TRUE.
			  IF (SBUF(1:2).eq.'-f') FAKET =.TRUE.


                enddo !(i=1,NumPar)
* (check if file name was last parameter        )
* (if the file exist, use free-format data input)
                filpar=SBUF
                open(10,file=filpar,status='OLD',Err=8001)
                close(10)
                FREFM=.TRUE.            ! use then free-form format
 8001           continue
	else
        WRITE(6,1) confile
        stop
	endif ! (NumPar.ge.1)
        
*(Check configuration file ptimer.cfg)
        
        open(30,file=confile,status='old',err=700)
        read(30,'(a)',err=701) fileph
        read(30,'(a)',err=702) filobs
        read(30,'(a)',err=703) fileop
        read(30,'(a)',err=704) filleap
        read(30,'(a)',err=705,end=705) filutc
        goto 802
 700    write(*,*) ' >> !! Configuration file not found >'//confile
        stop ' '
 701    fileph='de200.bin'
 702    filobs='obsys.dat'
 703    fileop='eop.dat'
 704    filleap='leap.sec'
 705    filutc='time.dat'
 802    close(30)

* DE200 / DE405 defined from file name of Ephemerides
	  do i = 1, len(fileph)-3
	  if (fileph(i:i+2).eq.'405') then
* CDE405 = .true. for DE405, = .false. for DE200
	  CDE405 = .true.
	  endif
	  enddo
        
* (Open output file)
	IO=8
        open(unit=IO, file='_tim.apr')

* (Display comments and short help)

        IF (SHCOM) then
        WRITE(6,1) confile
        stop
        endif
        
*(Replace +/- in PSR name and build string for input filename)
*(Check if it is free-format input - TEMPO2/TIMAPR2 or higher - used)
*(Or old format - TEMPO/TIMAPR1)
        
        if (.NOT.FREFM) then
	   if (PNM(5:5).EQ.'-') then
	    PP='m'
           else
	     if (PNM(5:5).EQ.'+') then
	      PP='p'
	     else                 ! check if free-format file used
         FREFM=.TRUE.           ! use then free-form format

 8002    continue
	endif
	endif !(PNM(5:5).EQ.'-'/'+')
        endif ! (.NOT.FREFM)
                
* (Open Pulsar's parameters - file PNM.p)
        
        if (.NOT.FREFM) then
        filpar=PNM(1:4)//PP//PNM(6:7)//'.p'
	  filtoa=PNM(1:4)//PP//PNM(6:7)
        open(unit=10,file=filpar,status='OLD',Err=9005)
        endif ! (.NOT.FREFM)
        
        psrname(1:7)=PNM(1:7)
        
*(If free-form format then read PSR parameters========)
        
        if (FREFM) then
        open(unit=funit,file=filpar,status='OLD',
     &      form='formatted',Err=9005)
       call readpar(funit,psrname,found,norb,neph,f,
     +     p0,p1,p2,DIST,Pepoch,aps,dps,pmra,pmdec,dm,
     +     X,E,EDOT,T0,pb,om,OMDOT,RG,PBDOT,SI,am,Gc3,
     +     xdot,
     +     bp,bpp,dr,dthet,a0,b0,pbd2,pbd3,xd2,IE,NITER)
       if(norb.gt.0) BINY=.TRUE.
       if (BINY) then
c convert units
       do i=1,norb
       Gc3(i)=Gc3(i)*4.92549    ! Msun->us
       Pbdot(i)=Pbdot(i)*365.25*Sday*1.d-06       ! mcs/yr <- s/s^-12
       enddo
       endif
        if (.not.found) then
           close(funit)
           write(*,*) ' > STOP: no TOAS line found in '//filpar
           stop
        endif
        endif ! (FRFM)
        
* (Show program version & date of calculation)
        
        CALL gtdate(ardt)
	ih =ardt(3)
	imm=ardt(2)
	is =ardt(1)
        WRITE(* ,'(2X,14HPTIMER (v1.00),2x,
     &	I2,1H-,a3,1H-,I4)') IS,chmon(IMM),IH
	WRITE(* ,'(2X,11HPulsar  PSR,1x,A9/)') psrname        
        WRITE(IO,'(2X,14HPTIMER (v1.00),2x,
     &	I2,1H-,a3,1H-,I4)') IS,chmon(IMM),IH
	WRITE(IO,'(2X,11HPulsar  PSR,1x,A9/)') psrname
        
        if(.not.FREFM) then
        
* (Read header of the parameter file in old format)
        READ(10,'(a100)') head100
        read(head100,'(1x,55I1)')  (IE(I),I=1,55)
        endif !(.not.FREFM)
        
*(Generate apparent periods)
        if(IE(55).eq.1) GENPER=.TRUE.
* (Generate artifitial TOAs)
        if(IE(55).eq.2) FAKET=.TRUE.
* (Read parameters if need P or fake TOAs)
        if(FAKET.OR.GENPER) then
        if(.not.FREFM) then
* (read parameters for P or fake TOA in old format)
        if(FAKET) read(10,*,err=2001)  FR0D00,JDD00,tstep,NOPTS,WegRMS
	  if(GENPER) read(10,*,err=2001) FR0D00,JDD00,tstep,NOPTS
        goto 2002
 2001   write(*,*)' ERROR: incorrect format of parameter string',filpar
        write(*,*)' EXAMPLE: 1410 2451344.5 0.1 100 10'
        stop
 2004   write(*,*)' ERROR: incorrect format of parameter string',filpar
        write(*,*)' EXAMPLE:EF 1410 2451344.5 0.1 100 10'
        stop
 2002   continue
        if((JDD00.gt.3.0d6).or.(JDD00.lt.2.4d6)) goto 2001
        else ! (not.FREFM)
* (read parameters for P or fake TOA in free-form format)
        if(FAKET) read(funit,*,err=2004) head100(57:58),
     &   FR0D00,JDD00,tstep,NOPTS,WegRMS
	   if(GENPER)read(funit,*,err=2004) head100(57:58),
     &   FR0D00,JDD00,tstep,NOPTS
        endif! (not.FREFM)
	endif ! (GENPER.FAKET)
        
        if(.not.FREFM) then        
* (f4)
        if(IE(5).eq.2) then
        IE(5)=1  ! f3
        IE(60)=1 ! f4
        endif
* (f5)
        if(IE(5).eq.3) then
        IE(5)=1  ! f3
        IE(60)=1 ! f4
        IE(61)=1 ! f5
        endif
* (f6)
        if(IE(5).eq.4) then
        IE(5)=1  ! f3
        IE(60)=1 ! f4
        IE(61)=1 ! f5
        IE(62)=1 ! f6
        endif
* (f7)
        if(IE(5).eq.5) then
        IE(5)=1  ! f3
        IE(60)=1 ! f4
        IE(61)=1 ! f5
        IE(62)=1 ! f6
        IE(63)=1 ! f7
        endif
* (f8)
        if(IE(5).ge.6) then
        IE(5)=1  ! f3
        IE(60)=1 ! f4
        IE(61)=1 ! f5
        IE(62)=1 ! f6
        IE(63)=1 ! f7
        IE(64)=1 ! f8
        endif
* (DMdot1)
        if(IE(11).eq.2) then
        IE(11)=1 ! DM(1)
        IE(70)=1 ! DM(2)
        endif
* (DMdot2)
        if(IE(11).eq.3) then
        IE(11)=1 ! DM(1)
        IE(70)=1 ! DM(2)
        IE(71)=1 ! DM(3)
        endif
* (DMdot3)
        if(IE(11).eq.4) then
        IE(11)=1 ! DM(1)
        IE(70)=1 ! DM(2)
        IE(71)=1 ! DM(3)
        IE(72)=1 ! DM(4)
        endif
* (DMdot4)
        if(IE(11).eq.5) then
        IE(11)=1 ! DM(1)
        IE(70)=1 ! DM(2)
        IE(71)=1 ! DM(3)
        IE(72)=1 ! DM(4)
        IE(73)=1 ! DM(5)
        endif
* (DMdot5)
        if(IE(11).eq.6) then
        IE(11)=1 ! DM(1)
        IE(70)=1 ! DM(2)
        IE(71)=1 ! DM(3)
        IE(72)=1 ! DM(4)
        IE(73)=1 ! DM(5)
        IE(74)=1 ! DM(6)
        endif
* (DMdot6)
        if(IE(11).ge.7) then
        IE(11)=1 ! DM(1)
        IE(70)=1 ! DM(2)
        IE(71)=1 ! DM(3)
        IE(72)=1 ! DM(4)
        IE(73)=1 ! DM(5)
        IE(74)=1 ! DM(6)
        IE(75)=1 ! DM(7)
        endif
        
* (Read a priory values of parameters of PSR from file)

	READ(10,*,err=44)
     &  APS,DPS,P0,P1,P2,Pepoch,DM(1),DIST,PMRA,PMDEC
        go to 45
 44     write(*,'(1x,23HERROR in data format of,3x,a50)') filpar
        stop
 45     continue
        
* (check if PSR is binary)
         read(head100,'(i1)') i
         norb=i
         if (norb.gt.3) norb=3     ! not more than 3 companions
         if (norb.ge.1) biny=.true.

* (Read orbital and PNn parameters if PSR is binary)

        if (BINY) then
         do k=1,norb
         READ(10,*,err=44) Pb(k),X(k),e(k),T0(k),Om(k),RG(k),
     &   Omdot(k),Pbdot(k),Si(k),Gc3(k)

* (Guess initial values for these small parameters to be zero)

         edot(k) =0.0d0
c         Xdot(k) =0.0d0
         dthet(k)=0.0d0
         RSG(k)  =0.0d0
         enddo                     ! (k=1,norb)
        endif !(BINY)

* (Close file of PSR parameters)
        CLOSE(unit=10)

* (Calculate spin frequency and derivatives)

	F(1) =  1.D0/P0
	F(2) = -P1/P0/P0
	F(3) = -P2/P0/P0+P1**2*2.D0/P0/P0/P0
        if(P2.EQ.0.0d0) F(3)=0.0d0
        
* (in old format these initial f(n) set by 0)
        do i=4,9
        F(i) = 0.0d0
        enddo
        do i=2,7
        DM(i) = 0.0d0
        enddo
        endif ! (NOT FREFM===============================)
        
* (Select method of LSQ - Weighted or not)
        if(IE(52).eq.1) IWEG=.FALSE.

* (Earth Orientation Parameters dUT1,dX,dY)
        if(IE(53).eq.1) EOP=.FALSE.

* (Don't write residuals to disk)
        if(IE(54).eq.1) WRITRS=.FALSE.

* (Convert ( h.m s ) & ( o.' " ) to radians)

	  CALL RADIAN(DPS,0,DPSR)
	  CALL RADIAN(APS,1,APSR)

* (Print assumed values of parameters)

	  CALL RADEG(APSR,1,APS,WrsA)
        CALL RADEG(DPSR,0,DPS,WrsD)

        WRITE(IO,76) WrsA,WrsD,PMRA,PMDEC,P0,P1,P2,
     &  Pepoch,DM(1),DIST
        if(SHOWPR) WRITE(*,76) WrsA,WrsD,PMRA,PMDEC,P0,P1,P2,
     &  Pepoch,DM(1),DIST

* (Print spin and spin derivatives)
        if(F(1).ne.0.0d0) write(io,'(8x,6Hf0   :,1x,(1pd24.16),1x,
     &  3H1/s)') f(1)
        do i=2,9
        write(pp,'(i1)') i-1
        write(pw,'(i1)') i
        if(F(i).ne.0.0d0) write(io,'(8x,1Hf,a1,4Hdot:,1x,
     &  (1pd24.12),1x,3H1/s,a1)') pp,f(i),pw
        enddo ! i=2,9
* (DM derivatives )
        do i=2,7
        write(pp,'(i1)') i-1
        if(DM(i).ne.0.0d0)write(io,'(8x,3HDMd,a1,2H :,1x,(1pd24.12),1x,
     &   9Hpc*cm3/yr,a1)') pp,dm(i)*Sday**(i-1)*365.25**(i-1),pp
        enddo ! i=2,7
        
        if (BINY) then
         do k=1,norb
	 if (norb.gt.1) write(IO,'(/1x,i3,1x,5HOrbit)') k
         if (SHOWPR) WRITE(*,77) Pb(k),X(k),e(k),Om(k),T0(k),
     &                      RG(k),Omdot(k),Pbdot(k),Gc3(k),Si(k)
         WRITE(IO,77) Pb(k),X(k),e(k),Om(k),T0(k),
     &                      RG(k),Omdot(k),Pbdot(k),Gc3(k),Si(k)


* (Convert units of parameters to those used by the program)
         Pb(k)=Pb(k)/PI2
         Om(k)=Om(k)/360.0*PI2                   ! deg -> rad
         RG(k)=RG(k)*1.d-06                      ! mcs -> sec
         Omdot(k)=Omdot(k)/360.0*PI2/365.25/Sday ! deg/yr -> rad/s
         Pbdot(k)=Pbdot(k)/365.25/Sday*1.d-06/PI2! mcs/yr -> s/s
c         edot(k) =edot(k)/365.25/Sday            ! 1/yr -> 1/s
         Gc3(k)  = Gc3(k)*1.d-06                 ! mcs -> s
        enddo !(k=1,norb)
        endif !(BINY)
        
        write(IO,'(80(1H-))')

        PRLX= 1.20962D-03/DIST                   ! pc  -> sec
	Pepoch=Pepoch+2400000.5D0
	JDBPOS=Pepoch              ! Epoch for PSR position same as for Period
        PMRA =PMRA/(ArcRd*Sday*365.25D0*1.D+03)  ! mas/yr -> rad/s
        PMDEC =PMDEC/(ArcRd*Sday*365.25D0*1.D+03)! mas/yr -> rad/s
        
* (Read Number of TOA point & Observatory coordinates)
*
* The format of data is written into header of TOAs-file of data.
* Observatories - in 1st symbol, defines specific format in data:

        write(IO,'(a)') ' Configuration file               >'//confile
        write(IO,'(a)') ' Leap second   file               >'//filleap
        
* (Coordinates of observatory. Read from file OBSYS.DAT)
        
        open(10,file=filobs,status='old',err=5)
	open(20,file='_site.dat')
        goto 9
 5      write(*,*) ' ERROR: Not found ',filobs
	write(*,*) ' Check config file ptimer.cfg'
        stop
 9      continue
        write(IO,'(a)') ' Observatory data (see _site.dat) >'//filobs
        NSITES=0
 7      read(10,'(3f15.0,25x,a1,2x,a2)',end=8)XTOP,YTOP,ZTOP,
     +  SiteL(NSITES+1),SiteC(NSITES+1)
        NSITES=NSITES+1
	if((XTOP+YTOP+ZTOP).NE.0.0) then
* (The site coordinates given as Geodetic (Lat,Long,Elevat))
        if (DABS(ZTOP).LT.10000.0d0) then
* (Convert (h.ms) & (o.'") to radians)
        CALL RADIAN(XTOP*1.0d-4,0,FIPLR)
        CALL RADIAN(YTOP*1.0d-4,0,LPLR)   ! FIPL,LPL in deg
        HEIGH=ZTOP
* (Projections of observatory [XTOP,YTOP,ZTOP] in GRS in sec.)
        CALL GEO(FIPLR,LPLR,HEIGH,XTOP,YTOP,ZTOP)
        else ! (ZTOP<10000)
* (The site coordinates given as (X,Y,Z) )
* (Geodetical coordinates of site)
        CALL XYZ2GEO(XTOP,YTOP,ZTOP,FIPLR,LPLR,HEIGH)
	XTOP=XTOP/SpdL
	YTOP=YTOP/SpdL
	ZTOP=ZTOP/SpdL
        endif ! (ZTOP<10000)
	endif ! (XTOP+YTOP+ZTOP)

	XTOPi(NSITES)=XTOP
	YTOPi(NSITES)=YTOP
	ZTOPi(NSITES)=ZTOP
* Write observatory coordinates and codes in _site.dat
        CALL RADEG(FIPLR,0,TMP,WrsA)
        CALL RADEG(LPLR,1,TMP,WrsD)
        WRITE(20,
     & '(1x,a1,1x,a2,1x,4HXYZ=,3f11.1,1x,3HFi=,a11,1x,3HLa=,a11)')
     & SiteL(NSITES),SiteC(NSITES),XTOP*SpDl,YTOP*SpDl,ZTOP*SpDl,
     & WrsA(1:11),WrsD(1:11)
        if (NSITES.GT.MAXSIT) then
        write(*,*) ' Max.No of Sites=',MAXSIT,' Increase MaxSIT.'
        stop
        else 
        goto 7
        endif !(NSITES>MaxSit)
 8      close(10)
	close(20)
        
*(Check if clock correction file 'filutc' exist. Unit 124)
         open(124,file=filutc,status='OLD',Err=8005)
         CLKCOR=.TRUE.
         close(124)
* READ CLOCK CORRECTION DATA if file=filutc exist
        if (CLKCOR) then
        write(IO,'(a)') ' Observatory clock corrections    >'//filutc
        call READCLK(io,filutc,clkdate,clsite,clockc)
        endif
 8005   continue        
        
* (Reading of TOA in UTC - scale or in TT (JPL format only))
* (Read observed TOAs using specific format of observatory)

        WGM=0.0D0
        TOAmin=1.0d10
        TOAmax=0.0d0
        i=0
        dTOAs=0.0d0
        minWEG=0.0d0
        maxWEG=1.0d12
        facWEG=1.0d0
        eqdWEG=0.0d0
        frqMIN=0.0
        frqMAX=1.0d12
        skipTOA=.false.
        dNPHAS0=0
        weg00=0.0d0
        phORB1=0.0d0
        phORB2=1.0d0
        JDstart=0.0d0
        JDend  =100000.0d0

        if(FAKET.OR.GENPER) then
	      R0=0.0d0
* (Get Observatory code)
            ChSIT(1:2)=head100(57:58)
            IndSIT(1)=-1
            do k=1,NSITES
            if(ChSIT.EQ.SiteC(k)) IndSIT(1)=k
            enddo!(k)
            if (IndSIT(1).EQ.-1) then
            write(*,*)' ERROR: Unknown Observatory code ',ChSIT
            stop
            endif!(IndSIT(1)=-1)

        if(FAKET)write(IO,'(/1x,24HGENERATION of fake TOAs:)')
        if(FAKET)write(*,'(/1x,24HGENERATION of fake TOAs:)')
	if(NOPTS.GT.NMI) then
        write(*,'(10X,29H! Too many points! Reduced to,i7)')NMI
        NOPTS=NMI
	endif ! (I.GT.NMI)
        
        if(GENPER)write(IO,'(/1x,32HCALCULATION of observed periods:)')
        if(GENPER)write(*,'(/1x,32HCALCULATION of observed periods:)')
        write(IO,'(1x,5HData:,f24.12,1H,f24.12)')JDD00,JDD00+NOPTS*tstep
        write(IO,'(1x,5HStep:,f24.12,2Hd=,f24.12,1Hs)')tstep,tstep*Sday
        write(IO,'(1x,10HFrequency:,f11.5,1x,3HMHz)') FR0D00
        write(IO,'(1x,10HMean RMS :,f11.5,1x,3Hmcs)') WegRMS
        write(IO,'(1x,11HObservatory,1x,a2,1x,
     &  41H(codes in _site.dat, data from obsys.dat))') ChSIT
        WEGRMS=WEGRMS*0.333
*(generate artifitial set of TOAs)
        
        if(FAKET) then
            do i=1,NOPTS
            JDDi  = JDD00+dfloat(i-1)*tstep+1.111111d-6
            JDD(i)= dint(JDDi)+0.5d0
            RM(i) = (JDDi-JDD(i))*Sday
            if (RM(i).LT.0.0d0) then
               JDD(i)=JDD(i)-1.0d0
               RM(i)=RM(i) + Sday
            endif
            DNJUMP(i)=1
            dNPHAS(i)=0.0d0
            Y(I)  =JDDi-2400000.5d0
c TOA error - simulate
            Weg(i)=rndm(i) * WegRMS * 2.0
	      if (Weg(i).lt.0.0 ) then
	        Weg(i) = - weg(i)
	      endif
	      if (Weg(i).le.0.01) then
		    Weg(i) = WegRMS
	      endif
c 	      Weg(i) = rand(i) * WegRMS
            IndSIT(i)=IndSIT(1)
            WGM   =WEG(i)+WGM
            FR0D(i)=FR0D00
            enddo! i
            TOAmin=JDD(1)
            TOAmax=JDD(NOPTS)
*(generate apparent periods)
        else ! (GENPER)
	    KCNT=1
            EOP=.FALSE.
            mBT=.TRUE.
            FR0D(1)=FR0D00
            open(123,file='_obsp.dat')
            goto 63
        endif ! (FAKET.GENPER)

c=======================================================================
        else !         if(FAKET.OR.GENPER) ??????
        
        write(*,*) ' Reading TOA data..'
        if(FREFM)write(IO,'(a)')
     &   ' Free format input data from      >'//filpar
        
* (Log unit number for TOA data)
* funit = 11 for free-format file (2145.tpo) <-filpar new format
* fuinc = 13 for free-format INCLUDE file )  <-filpar new format
* 
* toaun = fuinc/funit (switcher) in free-format
        
* (Open Times of Arrival TOAs - file PNM)
        if (.not.FREFM) then
        toaun=9                             ! <-filtoa old format
	if((.NOT.FAKET).AND.(.NOT.GENPER))
     &	open(toaun,file=filtoa,status='OLD',Err=9004)
        write(IO,'(a)') ' TOA data file                    >'//filtoa
        else
        toaun = funit            ! =11 (xxxx.tpo)        
        endif ! (.not.FREFM)
        
        i=0
        iline=0
c=======================================================================
        
 18     i=i+1              ! counter for lines containing TOAs
        iline=iline+1      ! counter of line number
        
            READ(toaun,'(a80)',end=21) ST80
            
*(check if there were a command in TOA file ---------------------)
* C,c,# are commentary lines
            
            if (ST80(1:1).eq.'C'.or.ST80(1:1).eq.'c'
     &       .or.ST80(1:1).eq.'#') then
            i=i-1
            goto 18
            endif!(St80='C')
*(SKIP) skip data block in TOA
            if (ST80(1:4).eq.'SKIP') then
            skipTOA=.true.
            i=i-1
            goto 18
            endif!(St80='SKIP')
*(NOSKIP) end of skip data block in TOA
            if (ST80(1:6).eq.'NOSKIP') then
            skipTOA=.false.
            i=i-1
            goto 18
            endif!(St80='NOSKIP')
*(END) line
            if (ST80(1:3).eq.'END') then
            i=i-1
            goto 19
            endif!(St80='END')
* check if there were SKIP line before then
* ignore string. Only SKIP and NOSKIP may appear before any command.
            if (skipTOA) then
            i=i-1
            goto 18
            endif ! (.not.skipTOA)
            
* skip blank lines if free form format file currently read
            
            if (FREFM) then                                   !***
            if ((toaun.eq.funit).and.(ST80(1:1).EQ.' ')) then !**
            i=i-1
            goto 18
            else  ! (toaun.eq.funit)                          !**
            
* open file followed by INCLUDE command
* if the data read from the parameter file
              if(ST80(1:7).eq.'INCLUDE') then
                 k=8                          ! delete blanks in filename
 6001            if (st80(k:k).eq.' ') k=k+1
                 if((k.lt.80).and.(st80(k:k).eq.' ')) goto 6001
                filtoa=ST80(k:80)
                open(fuinc,file=filtoa,status='OLD',Err=8003)
                toaun=fuinc
                write(IO,'(a)')
     &          ' TOA data file                    >'//filtoa
                goto 8004
 8003           continue
                write(*,*) ' !! No TOA file >'//filtoa
 8004           continue
                i=i-1
                goto 18
              endif ! ('INCLUDE')
            endif ! (toaun.eq.funit)                          !**
            endif ! (FREFM)                                   !***
            
*(check if there were a command in TOA file ---------------------)

*(JUMP)     
            if (ST80(1:4).eq.'JUMP') then
            JUMPED=.NOT.JUMPED
            if (JUMPED) then
            NJUMP=NJUMP+1
            endif!(JUMPED)
            i=i-1
            goto 18
            endif!(St80='JUMP')
*(DTOAS) add time corrections to the TOAs    
            if (ST80(1:4).eq.'TIME') then
            read(ST80,'(4x,f24.0)',err=26) dTOAs
            i=i-1
            goto 18
            endif!(St80='TIME')
*(EMIN) minimal uncertainty
            if (ST80(1:4).eq.'EMIN') then
            read(ST80,'(4x,f24.0)',err=26) minWEG
            i=i-1
            goto 18
            endif!(St80='EMIN')
*(EMAX) minimal uncertainty
            if (ST80(1:4).eq.'EMAX') then
            read(ST80,'(4x,f24.0)',err=26) maxWEG
            i=i-1
            goto 18
            endif!(St80='EMAX')
*(FMIN) minimal obs.frequency
            if (ST80(1:4).eq.'FMIN') then
            read(ST80,'(4x,f24.0)',err=26) frqMIN
            i=i-1
            goto 18
            endif!(St80='FMIN')
*(FMAX) maximal obs.frequency
            if (ST80(1:4).eq.'FMAX') then
            read(ST80,'(4x,f24.0)',err=26) frqMAX
            i=i-1
            goto 18
            endif!(St80='FMAX')
*(EFAC)     
            if (ST80(1:4).eq.'EFAC') then
            read(ST80,'(4x,f24.0)',err=26) facWEG
            i=i-1
            goto 18
            endif!(St80='EFAC')
*(EQUAD)    
            if (ST80(1:4).eq.'EQUA') then
            read(ST80,'(5x,f24.0)',err=26) eqdWEG
            i=i-1
            goto 18
            endif!(St80='EQUAD')
*(PHASE) add N phase turns to the TOAs
            if (ST80(1:5).eq.'PHASE') then
            read(ST80,'(5x,i24)',err=26) dNPHAS0
            corrPHAS=.true.
            i=i-1
            goto 18
            endif!(St80='PHASE')
*(PHA1) beginning of orbit phase which is ignored in fit
            if (ST80(1:4).eq.'PHA1') then
            read(ST80,'(4x,f24.0)',err=26) phORB1
            skipORB=.true.
            i=i-1
            goto 18
            endif!(St80='PHA1')
*(PHA2) end orbit phase which is ignored in fit
            if (ST80(1:4).eq.'PHA2') then
            read(ST80,'(4x,f24.0)',err=26) phORB2
            skipORB=.true.
            i=i-1
            goto 18
            endif!(St80='PHA1')
*(JDSTART) start JD, earlier data will be ignored in fit
            if (ST80(1:6).eq.'JDSTAR') then
            read(ST80,'(7x,f24.0)',err=26) JDstart
            JDSEGM=.true.
            if(JDstart.gt.2400000.0) JDstart=JDstart-2400000.5d0
            i=i-1
            goto 18
            endif!(St80='JDSTA')
*(JDEND) end JD, later data will be ignored in fit
            if (ST80(1:4).eq.'JDEN') then
            read(ST80,'(5x,f24.0)',err=26) JDend
            JDSEGM=.true.
            if(JDend.gt.2400000.0) JDend=JDend-2400000.5d0
            i=i-1
            goto 18
            endif!(St80='JDEN')
*(SIGMA) force uncertainties equal to sigma
            if (ST80(1:5).eq.'SIGMA') then
            read(ST80,'(5x,f24.0)',err=26) weg00
            i=i-1
            goto 18
            endif!(St80='PHASE')
*(MODE) not in ptimer
            if (ST80(1:4).eq.'MODE') then
            i=i-1
            goto 18
            endif!(St80='MODE')
*(END) not in ptimer
            if (ST80(1:3).eq.'END') then
            i=i-1
            goto 18
            endif!(St80='END')
            
*(unknown command)
 26         continue
            goto 27
            if (ST80(1:1).ne.' ') then
            write(*,'(a,i4,1x,3a)')
     & ' >> Unknown command ignored, line>',
     &   iline,' file: ', filpar(1:13), ST80(1:40)
            i=i-1
            goto 18
            endif! (St80=' ')
*(skip TOA block)
 27         if (skipTOA) then
            i=i-1
            goto 18
            endif!(skipTOA)
*(---------------------------------------------------------------)
	if(I.GT.NMI) then
	close(9)
        write(*,'(10X,34H! >> Too many points! Increase NMI)')
	stop
	endif ! (I.GT.NMI)
*(---------------------------------------------------------------)
* Find out which format is used and
* Read TOAs
* in nfmt=0,1 18x insteed of 15x because ddm read wrong
 1005   format(11x,i5,f14.13,f9.2,f12.5,f9.6,1x,a2)     ! nform=5 TIMAPR1
 1006   format(a1,14x,f9.0,1x,i5,f15.14,f8.2,18x,f10.0) ! nform=0 Princeton1
 1007   format(a1,14x,f9.0,i5,f15.14,f9.2,18x,f10.0)    ! nform=1 Princeton2
 1010   format(a1,14x,f9.0,i5,f15.14,3x,f6.2,18x,f10.0) ! nform=4 Princeton3
 1008   format(25x,f9.0,I7,f14.13,f8.6,f8.1,8x,a1)      ! nform=2 Parkes
 1009   format(9x,i5,f14.13,f6.2,f11.4,f10.6,2x,a2)     ! nform=3 Tempo ITOA
* nform = -1 -no format known
        nform= -1
* Timapr old format
        if(ST80(17:17).eq.'.')                            nform=5
* Tempo Princeton
        if((ST80(21:21).eq.'.').and.(ST80(31:31).eq.'.')) nform=0
* Tempo Princeton
        if((ST80(20:20).eq.'.').and.(ST80(30:30).eq.'.')) nform=1
* Tempo Princeton
        if((ST80(21:21).eq.'.').and.(ST80(30:30).eq.'.')) nform=4    ! <<
* Tempo Parkes
        if(ST80(42:42).eq.'.')                            nform=2
* Tempo Princeton ITOA
        if(ST80(15:15).eq.'.')                            nform=3
        
* unknown format - skip this line
        if(nform.lt.0) then
        i=i-1
        if(ST80(1:30).eq.'                              ') goto 19 ! (END)
        write(*,*) ' >> Unknown format in TOA line ! SKIPPED line:'
        write(*,'(a)') ST80(1:80)
        goto 18
        endif
* read TOA in corresponding format
        if (nform.eq.5)
     +  READ(ST80,1005,end=19) II,TMP,WEG(i),FR0D(i),dDM(i),ChSIT
        if (nform.eq.0)
     +  READ(ST80,1006,end=19) codobs,FR0D(i),II,TMP,WEG(i),dDM(i)
        if (nform.eq.1)
     +  READ(ST80,1007,end=19) codobs,FR0D(i),II,TMP,WEG(i),dDM(i)
        if (nform.eq.4)
     +  READ(ST80,1010,end=19) codobs,FR0D(i),II,TMP,WEG(i),dDM(i)    ! <<
        if (nform.eq.2) then
        READ(ST80,1008,end=19) FR0D(i),II,TMP,TMP1,WEG(i),codobs
        TMP=TMP+TMP1*p0/Sday
        dDM(i)=0.0
        endif
        if (nform.eq.3)
     +  READ(ST80,1009,end=19) II,TMP,WEG(i),FR0D(i),dDM(i),ChSIT

* clock correction applied if nform=0
        if((nform.eq.1).or.(nform.eq.2).or.(nform.eq.4)) nform=0
* clock correction NOT applied if nform=5
        if(nform.eq.3) nform=5
        
        if(II.lt.30000) II=II+39126 ! Convert 1966 days to MJD(Princeton)
        
            if(II.EQ.0) goto 19              ! exit from TOA file if blank
            RM(i) = TMP*Sday                 ! To UTC
* (Get Observatory code)
            IndSIT(I)=-1
            if (nform.ne.0) then             ! compare by 2-letter code
             do k=1,NSITES
             if(ChSIT.EQ.SiteC(k)) IndSIT(I)=k
             enddo!(k)
            else                             ! compare by 1-letter code
             do k=1,NSITES
             if(codobs.EQ.SiteL(k)) IndSIT(I)=k
             enddo!(k)
             ChSIT=codobs
            endif !(nform=0)
            if (IndSIT(I).EQ.-1) then
            write(*,*) ' >> ERROR: Unknown Observatory code ',
     &      ChSIT,', TOA N',i
            stop
            endif !(IndSIT(i)=-1)

        JDD(I) = 2400000.5D0 + dfloat(II)
        RM(i)  = RM(i) + dTOAs
        dNPHAS(I)=DNPHAS0
        
* (Apply observatory clock correction)
	
      
        if(nform.eq.0) then
        tmp1=dfloat(II)+RM(i)/Sday
        call UTCCOR(codobs,tmp1,dclk,clkdate,clockc,clsite,maxclk)
        RM(i)=RM(i)+dclk*1.0d-6
        endif
        
        Y(I)=float(II)+RM(I)/Sday  ! temporary used as TOA array
        if(weg(i).le.0.0) weg(i)=0.001
        if (weg00.ne.0.0d0) WEG(i)=weg00
        WEG(i)= dsqrt( (WEG(i) * facWEG)**2 + eqdWEG**2)
        
        if (WEG(i).LT.minWEG) Weg(i)=minWEG
	WGM=WEG(I)+WGM

* (check that values are in range of day)
        if (RM(i).GE.Sday) then
            RM(i)=RM(i)-Sday
	    JDD(i)=JDD(i)+1.d0
	endif
	if (RM(i).LT.0.0) then
            RM(i)=RM(i)+Sday
	    JDD(i)=JDD(i)-1.d0
	endif

*(skip data if frequency not in range FMIN-FMAX)
        if ((FR0D(i).LT.frqMIN).OR.(FR0D(i).GT.frqMAX)) then
        i=i-1
        goto 18
        endif!(skipTOA)
*(skip data if TOA Error bigger than EMAX)
        if (WEG(i).GT.maxWEG) then
        i=i-1
        goto 18
        endif!(skipTOA)
*(skip data if JD is not in range JDstart-Jdend)
        if ((Y(i).LT.JDstart).OR.(Y(i).GT.JDend)) then
        i=i-1
        goto 18
        endif!(skipTOA)
        
*(count jumped segments)
        if (JUMPED) then
        DNJUMP(i)=NJUMP+1
        else
        DNJUMP(i)=1
        endif!(JUMPED)

* (find first and last TOA in data array)
        TOAmin=dmin1(TOAmin,JDD(i))
        TOAmax=dmax1(TOAmax,JDD(i))
        goto 18

* LOOP through TOAs=================================================
 21     continue       ! tested
         i=i-1
 19	continue

* check if it was opened include-file and go back if yes
        if (toaun.eq.fuinc) then
          close(toaun)
          toaun=funit         ! go back to (2145.tpo)
          goto 18
        endif ! (toaun.eq.fuinc)
	CLOSE(unit=9)
        if(CLKCOR) CLOSE(124) ! close clock correction file
* ==================================================================
        
        NOPTS=I
        
        endif !(FAKET.OR.GENPER)

        open(331,file='jj')
c jjjjjjj
* re-calculate mean error, 'cause some toas skept
* check indexes in JUMP array dnjump(i) because
* numbering may be wrong due to omitted segments.
        
        wgm=weg(1)
c        Jcorr=0
        Jcorr=DNJUMP(1)-1
        NJUMP=NJUMP-Jcorr
        do i=2,nopts
        wgm=wgm+weg(i)
        DNJUMP(i)=DNJUMP(i)-JCORR
        if ((DNJUMP(i)-DNJUMP(i-1)).gt.1) then
        NJUMP=NJUMP-1
        JCORR=JCORR+1
        DNJUMP(i)=DNJUMP(i)-JCORR        
        endif
           write(331,*) i,dnjump(i)
        enddo
        close(331)
c jjjjjjj
        
* (delete indexes 1-75, 101-140  no parameter fit)
        
        if(FAKET.OR.GENPER.OR.(IE(52).eq.2)) then
        do i=1,75
        IE(i)=0
        enddo
        do i=101,140
        IE(i)=0
        enddo        
        IE(1)=1
        endif

        WRite(*,*) NOPTS,' TOAs read'

        
        WGM=WGM/DFLOAT(NOPTS)

* (Find sorting order of TOAs and put pointers into IndORD(i))

        if(EOP) call sortord(Y,IndORD,NOPTS)

* (Select timing model)
* (Damour & Deruelle inverese timing formula)
* (or Blandford & Teukolsky one)

       if (BINY) then
         k=IE(51)
         if((k.eq.0).or.(k.gt.6)) k=1
         if(k.eq.1) then
         write(*,*) 'DD timing model'
         mDADE=.TRUE.
         endif ! (k=1)
         if(k.eq.2) then
         write(*,*) 'BT timing model'
         mBT  =.TRUE.
         endif ! (k=2)
         if(k.eq.3) then
         write(*,*) 'BT+ timing model'
         mBTP =.TRUE.
         endif ! (k=3)
         if(k.eq.4) then
         write(*,*) 'BT++ timing model'
         mBTPP  =.TRUE.
         endif ! (k=4)
         if(k.eq.5) then
         write(*,*) 'DD+ timing model'
         mDADEP =.TRUE.
         endif ! (k=5)

*(make sure that correct parameters were set to fit for given timing model)

       if(mBT) then
        do k=1,norb
        n=(k-1)*13
        IE(18+n)=0
        IE(19+n)=0
        IE(24+n)=0
        enddo! k
       endif!mBT

       if(mBT.or.mBTP.or.mBTPP) then
        do k=1,norb
        n=(k-1)*13
        IE(24+n)=0
        enddo! k
       endif!mBT

       endif ! (BINY)

* (Earth's POLE motion (dX,dY) and dUT1 corrections)

        if(EOP) then
        write(*,*) 'Reading pole coor.dX,dY,dUT1'
        call geteop(Y,IndORD,NOPTS,DXPOLE,DYPOLE,DUT1,lack,fileop)
        endif ! EOP
*--------------------------------------------------------------
* Report used data

        if(EOP) then
        if(lack.eq.0)write(IO,'(1x,18HEOP from IERS data,15x,1h>,a)')
     &   fileop
        if(lack.eq.1) write(IO,'(/1x,
     &  48HEOP ignored, no data or there is no file EOP.DAT)')
        if (lack.ge.2)write(IO,'(1x,32HEOP from IERS. No data for toas:,
     &  i5,1H-,i5,1x,a)')lack,NOPTS,fileop
        write(IO,'(1x,34HJPL Ephemerides DE200/LE200      >,a)') fileph
        write(IO,'(1x,
     &  46HTime ephemerides Fairhead L.& Bretagnon P.1990/)')
        else
        write(IO,'(/1x,24HEOP (dUT1,dX,dY) ignored)')
        endif! EOP
        write(IO,'(80(1H-))')
*--------------------------------------------------------------

* (Select PSR parameters to be improved.)
*  -------------------------------------

* (delete indexes 51-59)
        do i=51,59
        IE(i)=0
        enddo
        
* (jumps are stored from 80 pointer)

        if(NJUMP.ne.0) then
        do I=80,79+NJUMP
        IE(I)=1
        ENDDO
        endif ! (NJUMP.ne.0)

* (Degree of system = number of adjustement parameters)

	IDEGEQ=0
	DO 217 I=1,IDEGS
 217	IDEGEQ=IDEGEQ+IE(I)
        if (IDEGEQ.EQ.1) NITER=1
        if (FAKET) NITER=50

* (Description of adjusted PSR parameters)

	K=1
	WR(ISM(IE,1))=WLET(1)
 216	ISS1=ISM(IE,k)
	IF(K.NE.1) then
		IF(ISS1.NE.ISS0) WR(ISS1)=WLET(k)
		endif
	IF(ISS1.EQ.IDEGEQ) GOTO 218
	k = k+1
        ISS0=ISS1
	if(K.LE.IDEGS) goto 216
 218	CONTINUE

* (Make vector of pointers)
	DO 219 I=1,IDEGS-1
	IF(IE(I).LE.IDEGEQ) IEIN=IE(I)
	IE(I+1)=IEIN+IE(I+1)
 219	IF(IE(I+1).EQ.IEIN) IE(I+1)=IDEGEQ+1

* (Start time of program)

	call gttime(ardt)
	ih0 =ardt(1)
	imm0=ardt(2)
	is0 =ardt(3)
	write(IO,'(2X,10HTime Beg.=,2(I2,1H:)I2)')
     &	IH0,IMM0,IS0

* (Print report)

	WRITE(IO,209) NOPTS
* (first/last TOAs)

 	CALL CALDJD(TOAmin,NN,IM,TMP,WW)
 	IYE=TMP
 	WRITE(IO,87) NN,IM,IYE,TOAmin

	CALL CALDJD(TOAmax,NN,IM,TMP,WW)
 	IYE=TMP
 	WRITE(IO,88) NN,IM,IYE,TOAmax
        
* (JD segment)
        
        if(JDSEGM) then
        WRITE(IO,'(2x,12HSelected MJD,f12.3,1H-,f12.3)') JDstart,JDend
        endif ! (JDSEGM)
        
* (report commands in TOA data file)

	WRITE(IO,'(2x,22HMedian TOA uncertainty,f10.3,1x,3Hmcs)') wgm
        if (NJUMP.ne.0) write(IO,
     & '(1x,i2,1x,38HTOA segments will be fitted for offset)') NJUMP
        if ((phORB1.ne.0.0d0).and.(phORB2.ne.1)) write(IO,
     &'(2x,28HOrbital phase ignored in fit,f5.2,1H-,f4.2)')phORB1,phORB2

* (way of fit)

        if ((IDEGEQ.NE.1).and.(.not.GENPER).and.(.not.faket)) then
        WRITE(IO,'(/2x,6HGlobal,I4,28H-th order fit of parameters,
     &   ,I5,20H degrees of freedom.)') IDEGEQ,NOPTS-IDEGEQ
                if (iweg) then
                write(IO,'(2x,12HWeighted fit)')
                else
                write(IO,'(2x,16HWeights=1 in fit)')
                endif
        else
        WRITE(IO,'(/2x,20HNo parameters fitted)')
        endif

         if (BINY) then
         if (mDADE)
     &   write(IO,'(2x,25HUsing timing formula <DD>)')
         if (mDADEP)
     &   write(IO,'(2x,26HUsing timing formula <DD+>)')
         if (mBT)
     &   write(IO,'(2x,25HUsing timing formula <BT>)')
         if (mBTP)
     &   write(IO,'(2x,26HUsing timing formula <BT+>)')
         if (mBTPP)
     &   write(IO,'(2x,27HUsing timing formula <BT++>)')
         endif !(BINY)

* (Beginning of iterations. Initializations.)

	KCNT=0
	R001=0.0D0
	R0=0.0D0
        if(NJUMP.NE.0) then
        do k=1,NJUMP
        R0N(k)=0.0d0
        enddo
        endif! (NJUMP.ne.0)
        D0=0.0D0
	CHI2=0.0

*********************************************************************
*********************************************************************
* do 62 ! loop for visual check of CHI2 surface near global minima
*       ! should be organised starting from this point.
        
                     R0=R001
	III = 1
        SIG0=1.0D+20
        
************************************************************************
c=======================================================================
*(Here we shall retun for next iteration of adjusted parameters)
*-----------------------------------------------------------------------
 1001   CONTINUE
					KCNT=KCNT+1
* (Clear variables)

	DO 9999 I=1,IDEGS
	DO 9998 J=1,IDEGS
	 ONE(I,J)=0.0D0
 9998	 SUM(I,J)=0.0D0
	 COEF(I)=0.0D0
	 ER(I)=0.0D0
	 RHS0(I)=0.0D0
 9999	 RHS(I)=0.0D0
	D0 = 0.0D0

* (Loop for TOAs set)
*  -----------------
 63     continue
        DO 60 II=1,NOPTS
	IPNT=II

* (Topocentric TOA in UTC scale)

        if(GENPER) then
            IPNT=1
            JD  = JDD00+dfloat(ii-1)*tstep
            JDDi= dint(JD)+0.5d0
            UTC = (JD-JDDi)*Sday
            if (UTC.LT.0.0d0) then
            JDDi =JDDi-1
            UTC=UTC + Sday
            endif! (UTC<0)
            JD=JDDi
            DXP=0.0d0
            DYP=0.0d0
            DUTC1=0.0d0
            XTOP =XTOPi(IndSIT(1))
            YTOP =YTOPi(IndSIT(1))
            ZTOP =ZTOPi(IndSIT(1))
        else
            JD   = JDD(IPNT)
	    UTC  = RM(IPNT)

* (Eearth Orientation Parameters)

        DXP  = DXPOLE(IPNT)*0.0001d0      ! Earth pole motion dX,dY - (rad.)
        DYP  = DYPOLE(IPNT)*0.0001d0      !
        DUTC1= DUT1(IPNT)                 ! DUTC1=UT1-UTC

        endif ! (GENPER)


* (Rectangular coordinates of PSR.)

	XP=DCOS(DPSR)*DCOS(APSR)
	YP=DCOS(DPSR)*DSIN(APSR)
	ZP=DSIN(DPSR)
        
        NJUMPN=DNJUMP(IPNT)-1 ! current jump No.

**         Time scales transformation.(1-3)
**	   -------------------------

        if (KCNT.EQ.1) then

* always OO 1-st observatory means ET - geocenter
	if (IndSIT(IPNT).eq.1) then
        EPHTM=.true.
	else
	EPHTM=.false.
        endif

* (Observatory coordinates)
        XTOP =XTOPi(IndSIT(IPNT))
        YTOP =YTOPi(IndSIT(IPNT))
        ZTOP =ZTOPi(IndSIT(IPNT))
        
* (1.from UTC to TDT/TAI)

	CALL TDTUTC(JD,TDT,UTC,filleap)

* (2.from TDT to TDB')

* calculate TB once and save difference TDB-TDT

*       CALL TTTB(JD,TDT,TDB)
	   if ( CDE405 ) then 
	   call tdbtrans(JD,TDT/86400.d0,TDB)
	   else
	   CALL TTTB(JD,TDT,TDB)
	   endif

        TDBSAV(IPNT) = TDB-UTC
        if (EPHTM) TDBSAV(IPNT) = 0.0d0
        endif ! (KCNT=1)

        JDB=JD
        TDB=UTC+TDBSAV(IPNT) ! TOA now is in TB scale

* (Earth's Coord. and Veloc. for time (JDB,TDB) of TOA from DE200/LE200)
	if ( CDE405 ) then
	CALL READEL405(JDB,TDB,IER,fileph)
	else
	CALL READEL(JDB,TDB,IER,fileph)
      endif

	IF(IER.EQ.0) GOTO 789
	ICW=0
        WRITE(*,*) 'Error in reading of Ephemeries ! =',IER
	STOP
 789	CONTINUE

        if (KCNT.EQ.1) then

* (Topocentric Roemer delay and EOP parameters)
        CALL DTTOP
     &(JD,UTC,DUTC1,XP,YP,ZP,XTOP,YTOP,ZTOP,TTS,SCTOP,DXP,DYP,SS,RTOPO)
        
        SCTOPS(IPNT)=SCTOP
 	TTSSAV(IPNT)=TTS

* (3.from TDB' to TDB diurnal term 1/c2*(V_earth,X_topo) in
*    relativistic time scale transformation ~2 mcs)
* (Norbert Wex.1996)

           TDB=TDB + (TABOUT(4,3)*RTOPO(1)+
     &       TABOUT(5,3)*RTOPO(2)+TABOUT(6,3)*RTOPO(3))*TAD

* (save the difference TDB-TDT)

           TDBSAV(IPNT) = TDBSAV(IPNT) +(TABOUT(4,3)*RTOPO(1)+
     &       TABOUT(5,3)*RTOPO(2)+TABOUT(6,3)*RTOPO(3))*TAD

        endif ! KCNT=1

* (Earth's coordinates XX,YY,ZZ-Rectang., AE,DE,RE-Spher.)

		XX=TABOUT(1,3)
		YY=TABOUT(2,3)
		ZZ=TABOUT(3,3)
		CALL SPHERC(XX,YY,ZZ,AE,DE,RE)

**	   Calculation of barycentric TOA.(1-7)
**	   ------------------------------------

* (1. Roemer term - 1)
* (   Earth's orbital motion )

	TOS=TA*(XX*XP+YY*YP+ZZ*ZP)
	TDB=TDB+TOS

* (2. Reomer's term - 2)
* (   Duirnial rotation of observer )

	TDB=TDB+TTSSAV(IPNT)

* (Doppler shift of frequency.)

	XX =TABOUT(4,3) ! Earth velocity in BSS
	YY =TABOUT(5,3)
	ZZ =TABOUT(6,3)
	SC =1.0D0-(XX*XP+YY*YP+ZZ*ZP)*TAD -SCTOPS(IPNT)

	FR =FR0D(IPNT)*0.01d0 * SC        ! frequency in [100MHz]

* (3. Dispersion delay. Go to 8 frequeency)

        ScDM=1.0d0/(FR*FR)/2.41D0
        TFR=(DM(1)+dDM(IPNT))*ScDM  ! Specific for Arecibo data
                                    ! format-excess DM.
        TMP1 = (JDB-JDBPOS)*Sday
c 14sec 
c       TFR=TFR+ ScDM * (                         ! DMdot
c    &   DM(2)*TMP1   +
c    &   Dm(3)*TMP1**2+
c    &   Dm(4)*TMP1**3+
c    &   Dm(5)*TMP1**4+
c    &   Dm(6)*TMP1**5+
c    &   Dm(7)*TMP1**6
c    &   )
c 14sec
        TFR=TFR+ScDM*(
     &  TMP1*(DM(2)+TMP1*(DM(3)+TMP1*(DM(4)+TMP1*(DM(5)+
     &  TMP1*(DM(6)+DM(7))))))
     &  )
        
        TDB=TDB-TFR

* (4. Shapiro delay for Sun & Jupiter.)

* (Vector c-o-m of Earth)

	XN=TABOUT(1,3)
	YN=TABOUT(2,3)
	ZN=TABOUT(3,3)

* (Vector c-o-m of SUN)

	XB=TABOUT(1,10)
	YB=TABOUT(2,10)
	ZB=TABOUT(3,10)
	XR=XN-XB
	YR=YN-YB
	ZR=ZN-ZB
	RNB=DSQRT(XR**2+YR**2+ZR**2)
	SCS=XR*XP+YR*YP+ZR*ZP
 	SCS=SCS+RNB
* (Shapiro for Sun)
	TSHAPS= SSUN*DLOG(SCS)
        
c	TSHAPS= SSUN*
c     &   DLOG(1.0+(XP*XB+YP*YB+ZP*YB)/DSQRT(XB**2+YB**2+ZB**2)/
c     & DSQRT(XP**2+YP**2+ZP**2))
        
* (Vector c-o-m of Jupiter)
        
	XB=TABOUT(1,5)
	YB=TABOUT(2,5)
	ZB=TABOUT(3,5)
	XR=XN-XB
	YR=YN-YB
	ZR=ZN-ZB
	RNB=DSQRT(XR**2+YR**2+ZR**2)
	SCJ=XR*XP+YR*YP+ZR*ZP
	SCJ=SCJ+RNB

* (Shapiro for Jupiter. Included, unless small)
        TSHAPJ= SJUP*DLOG(SCJ)
c oooooooo
c       U(IPNT,1) = - TSHAPJ * 1.0d3  ! in _kep.out ns
c       
* (Vector c-o-m of Saturn)
	XB=TABOUT(1,6)
	YB=TABOUT(2,6)
	ZB=TABOUT(3,6)
	XR=XN-XB
	YR=YN-YB
	ZR=ZN-ZB
	RNB=DSQRT(XR**2+YR**2+ZR**2)
	SCS=XR*XP+YR*YP+ZR*ZP
	SCS=SCS+RNB
* (Shapiro for Saturn. Max ~50 ns very small)
        TSHPSS= SSAT*DLOG(SCS)

* (5. Parallax of PSR.)

        VPR=(ZP*YN-YP*ZN)**2+(ZN*XP-XN*ZP)**2+(XN*YP-YN*XP)**2
        PR=-PRLX*VPR

        TDB = TDB+TSHAPS+TSHAPJ+PR
	JDB = JD

* (Save parameters from previous iteration.)

		APSR0 =APSR
		DPSR0 =DPSR
		PMRA0 =PMRA
		PMDEC0 =PMDEC
		R00   =R0
                if(NJUMP.NE.0) then
                do k=1,NJUMP
                R0N0(k)=R0N(k)
                enddo
                endif! (NJUMP.ne.0)
                do i=1,11
                f0(i)=f(i)
                enddo
                PRLX0=PRLX
                do i=1,11
                DM0(i)=DM(i)
                enddo
             if (BINY) then
               do k=1,norb
                RSG0(k) =RSG(k)
                X0(k)   =X(k)
                Om0(k)  =Om(k)
                Omdot0(k)=Omdot(k)
                Xdot0(k) =Xdot(k)
                Pbdot0(k)=Pbdot(k)
                Pbd20(k) =Pbd2(k) !<<
                Xd20(k)  =Xd2(k)  !<<
                Pbd30(k) =Pbd3(k) !<<                
                edot0(k) =edot(k)
                e0(k)   =e(k)
                Pb0(k)  =Pb(k)
                RG0(k)  =RG(k)
                Gc30(k) =Gc3(k)
                Si0(k)  =Si(k)
                dthet0(k)=dthet(k)
                enddo ! (k=1,norb)
             endif !(BINY)

* (6. PSR - proper motion)

	TPR=( PMDEC*( DCOS(DPSR)*DSIN(DE)-DSIN(DPSR)*DCOS(DE)*
     &	DCOS(AE-APSR))+
     &	PMRA*DCOS(DE)*DSIN(AE-APSR))
	TPR=TPR*TA*RE

* (7. Barycentric Time of Arrivial (JDB,TDB))

        TDB = TDB+TPR*(JDB-JDBPOS)*Sday
	JDB = JD

* (Time interval from epoch of period)

       DTB =(TDB)+(JDB-Pepoch)*Sday

c=======================================================================
* (*********************************)
* (This part is for binary PSRs only)

        DNP0=0.0d0
        if(GENPER) scl=0.0d0
*
       if (BINY) then
        if(GENPER) scl=0.0d0
        ShapT=0.0d0
*
        U(ipnt,1)=0.0d0
        U(ipnt,2)=0.0d0
        do k=1,norb

* (Time passed after perisatron passage)
        
        DTPER(k) =  (JDB-T0(k))*Sday+TDB

* (8. Keplerian motion of binary)

        ex(k)=(e(k)+edot(k)*DTPER(k))
        eSQ(k)=DSQRT(1.D0-ex(k)**2*(1.d0+dthet(k))**2)

* (Save parameters if we plan to use inverse timing formula and iterate it)

        DTBj=DTB
        TDBj=TDB
        DNPj=0.0d0
        nloops=0

 8000   DTPER(k) =  (JDB-T0(k))*Sday+TDBj
c ooooooooo????????????        
        RTP= DTPER(k) / (Pb(k) + 0.5d0*(Pbdot(k)+Pbd2(k)/2.0*DTPER(k)
     &  +Pbd3(k)/4.0*DTPER(k)**2)
     &  *DTPER(k) ) + RSG(k)
        RigKT(k)=RTP-RSG(k)                                   ! Need for LSQ
        RigKT(k)=RigKT(k)/Pb(k)

* (Eccentric anomaly EXC-Find it from Kepler's equation)

	NCIRC=DINT(RTP/PI2) ! No. of bin.periods ellapsed
        if(RTP.LT.0.0d0) NCIRC=NCIRC-1
	RTP=RTP-NCIRC*PI2
	EXC(k)=RTP

* (Iterate Kepler equation to get eccentric anomaly)
        
        
        if (ex(k).ne.0.0d0) then
 2378	TMP=EXC(k)
* (Newton's method)
	EXC(k)=EXC(k)-(EXC(k)-ex(k)*DSIN(EXC(k))-RTP)/
     &  (1.0d0-ex(k)*DCOS(EXC(k)))
	IF(DABS(EXC(k)-TMP).GT.1.D-14) GOTO 2378
        endif ! (ex=0)

        if (mDADE.OR.mBTP.OR.mBTPP.OR.mDADEP) then

* (Parameters using DD model)
        Ane(k)=2.0d0*datan(dsqrt((1.d0+ex(k))/(1.d0-ex(k)))
     &        * dtan(EXC(k)/2.0d0))
        if (Ane(k).LT.0.0d0) Ane(k)=Ane(k)+PI2
                Ane(k)=(Ane(k)+NCIRC*PI2)
                Om1=Om(k)+Omdot(k)*Ane(k)*Pb(k)
                Cw=dcos(Om1)
                Sw=dsin(Om1)
        else

* (Parameters using BT model)
                Cw=dcos(Om(k)+Omdot(k)*DTPER(k))
                Sw=dsin(Om(k)+Omdot(k)*DTPER(k))
        endif

        if (mBTPP.or.mDADEP) then
*  model BT++,DD+
        X1=X(k) + Xdot(k)*Ane(k)*Pb(k) + Xd2(k)/2.0*Ane(k)**2*Pb(k)**2
        else
        X1=X(k) + Xdot(k)*DTPER(k) + Xd2(k)/2.0*DTPER(k)**2
        endif

        AF(k)=X1*Sw
        BT(k)=eSQ(k)*X1*Cw

	Cu(k)=dcos(EXC(k))
	Su(k)=dsin(EXC(k))

        RDA=Cu(k)-ex(k)
        ABSI=AF(k)*RDA+BT(k)*Su(k)

* (Roemer plus Einstein time delay in binary system)

        if (mDADE) then
        RoemEp=-(AF(k)*RDA+(BT(k)+RG(k))*Su(k)) * F(1)
        else
        RoemEp=-F(1)*(AF(k)*RDA+(BT(k)+RG(k))*Su(k)
     & +   ABSI*(AF(k)*Su(k)-BT(k)*Cu(k))/
     &     (Pb(k)+0.5d0*(Pbdot(k)+Pbd2(k)/2.0*DTPER(k)
     & +    Pbd3(k)/4.0*DTPER(k)**2)
     & *    DTPER(k))/(1.D0-ex(k)*Cu(k)))
     & -   F(2)*DTBj*ABSI
        endif

* (Shapiro delay for companion)

        if (mBT) then
        ShapC=0.0d0
        Znn(k)=1.0d0
        else
        Znn(k) = 1.d0-ex(k)*Cu(k)-Si(k)*(Sw*RDA+eSQ(k)*Cw*Su(k))
        ShapC= 2.0d0*Gc3(k)*dlog(Znn(k))
        endif ! (mBT)

* (Phase due to classic/relativistic effects in binary system)

        DNP0k = RoemEp + ShapC*F(1)

* (Iterate inverse timing formula T=(t-DNP0k) in DD approach)

        if (mDADE) then                                  !(**)
                    if (dabs(DNP0k-DNPj).GT.1.d-14) then !(*)
                     DTBj=DTB  + DNP0k/F(1)
                     TDBj=TDB  + DNP0k/F(1)
                     DNPj=DNP0k
                     nloops=nloops+1
        if (nloops.LT.4000) then 
         DNPj=DNP0k
         goto 8000
        else
         write(*,*)
     &   'EXIT inverting timing formula(2000 loops) ',dabs(DNP0k-DNPj)
        endif !(nloops)
                    else                                !(*)
                     DTB =DTB  + DNP0k/F(1)
                     TDB =TDB  + DNP0k/F(1)
                     DNP0k=0.0d0
                    endif                               !(*)
       endif                                            !(**)

* (projected onto observer's l-o-s velocity of binary pulsar)

      if(GENPER) then
      xvel=X(k)/(Pb(k)+0.5d0*Pbdot(k)*DTPER(k))/esQ(k)
      TMP=(1.0d0-ex(k)*Cu(k))
      sFi=Su(k)*esQ(k)/TMP
      cFi=RDA/TMP
      scl=scl+xvel*(Cw*cFi-Sw*sFi+ex(k)*Cw)
      endif ! (GENPER)
      
       if(.not.GENPER) U(IPNT,1) = -(RoemEp/F(1))+U(ipnt,1)
       ShapT=ShapT+ShapC

       DNP0=DNP0+DNP0k
       
       enddo !(k=1,norb)
       endif !(BINY)
        
c=======================================================================
* (Phase and Residual of pulse.)
*  ----------------------------

* (Slow-down model of pulsar rotation - pulse phase)

        DNP=((((((((
     &    F(9)/362880.0d0   *DTB
     &   +F(8)/40320.0d0)  *DTB
     &   +F(7)/5040.0d0)  *DTB
     &   +F(6)/720.0d0)  *DTB
     &   +F(5)/120.0d0) *DTB
     &   +F(4)/24.0d0) *DTB
     &   +F(3)/6.0D0) *DTB
     &   +F(2)/2.0d0)*DTB
     &   +F(1))    *DTB   +R0
	DNP=DNP+DNP0

*(Calculate Apparent Period if requested)

        if(GENPER) then

*(Projected total velocity, Vobs=(SC-1)*SpDL,Vpsr=scl*SpDL)

        Vpsr=SCL*SpDL*0.001d0
        Vobs=(SC-1.0d0)*SpDL*0.001d0
        SCL=SCL+SC
        SCL=2.0d0-SCL
        Pbar=(0.5*P2*DTB+P1)*DTB+P0
        Pobs=Pbar/SCL

*(Vbar=(SCL-1)*SpDL)
        
        write(123,'(1x,i6,f10.3,f21.15,2f11.5$)')
     &   int(JD-2400000.5), UTC, Pobs,Vpsr,Vobs

*(Print orbital phase and Roemer delay if flag ORBSH=.true.)
        
        if (BINY.and.ORBSH) then
        do k=1,norb
                TMP=(JDB-T0(k))*Sday + TDB
                TMP=(TMP-DINT(TMP/Pb(k)/PI2)*Pb(k)*PI2)/Pb(k)/PI2
		if (TMP.LT.0.d0) TMP=TMP+1.d0
                write(123,'(f13.10,g17.6)') TMP,-(RoemEp/F(1))
        enddo !(k=1,norb)
        else
c           write(123,'(1x,c1)')
        endif !(BINY)
        goto 60
        endif !(GENPER)

*(add Jumps of data segments)

        if(NJUMPN.NE.0) then
         DNP=DNP-R0+R0N(NJUMPN)
        endif!(NJUMPN)

*(Add phase correction if needed)

        if (corrPHAS) DNP=DNP+dNPHAS(ipnt)

	DRN=DINT(DNP)
	DTI=DRN-DNP

	DRN=DRN+1.D0
	DT1=DRN-DNP
	DRN=DRN-2.D0
	DT2=DRN-DNP
	IF(DABS(DT1).LT.DABS(DTI)) DTI=DT1
	IF(DABS(DT2).LT.DABS(DTI)) DTI=DT2

*(Observed-Calculated, in seconds)-

	DTIS=-DTI/F(1)

* (Show the Residuals with Shapiro delay undeleted)

	if (BINY) U(IPNT,2) = (DTIS-ShapT)*1.d+06+U(IPNT,2) ! mcs

* (Construction of Syst.of Norm.Equations)
*  --------------------------------------

* (Assume weigths)
	if (WEG(IPNT).LT.1.d-6) then
		WEGHT=1.0
	else
		WEGHT=(WGM/WEG(IPNT))**2
	endif
	if (.NOT.iweg) WEGHT=1.0D0

	FRPDAT = F(1)+(F(2)+.5D0*F(3)*DTB)*DTB

c=======================================================================
* (Calculate partial derivatives and build matrix for LSQ)
*  ------------------------------------------------------
* (do not use in fit the points from phORB1 to phORB2)

        skipFIT=.false.

        if (biny) then
        do k=1,norb
        if (skipORB) then
        TMP=(DTPER(k)-DINT(DTPER(k)/Pb(k)/PI2)*Pb(k)*PI2)/Pb(k)/PI2
        if (TMP.LE.0.0d0) TMP=TMP+1.0d0
        if ((TMP.GE.phORB1).and.(TMP.LE.phORB2))
     &                    skipFIT=.true.
        endif! (skipORB)
        enddo! (norb)
        endif! (biny)

* (this procedure builds covar. matrix. Paramets - via CETBL7)

	if ((.not.skipFIT).and.(.not.FAKET).AND.(.not.GENPER))
     &    CALL LIN(NJUMPN,NORB)

        if(FAKET) RM(IPNT)=RM(IPNT) - DTIS

	DTIS       =DTIS*1.D+06
	Y(IPNT)    =SNGL(DTIS)
        TDBi(IPNT) =TDB

* (here is the end of loop over TOAs in set.)
* (In following LSQ solutions will be added)
* (to the previous values of parameters,)
* (and next iteration, if minima wasn't achieved, will be made again)

 60	CONTINUE ! loop for II (through whole TOAs set)
c=========================================================================
        
        if (GENPER) then        
        if (ORBSH) then
        St80=' _obsp.dat. Format:'//
     &   'MJD,UTC(s),Papp(s),Vpsr(km/s),Vobs(km/s)'
     &   //', OrbPhs, Roemer(s)'
        else
        St80=' _obsp.dat. Format:'//
     &   'MJD,UTC(s),Papp(s),Vpsr(km/s),Vobs(km/s)'
        endif ! (ORBSH)
        write(*,'(a80)') St80
        write(IO,'(a80)') St80
        close(123)
        stop
        endif ! (GENPER)

c=======================================================================
* (Calculate dispersion of residuals and chi2-value.)
*  --------------------------------------

* (delete constant term from residuals array Y(I))

	if (DELMEAN) then
                ymean=0.0
                DO I=1,NOPTS
                ymean= Y(I) + ymean
                ENDDO
                ymean=ymean/float(NOPTS)
                DO I=1,NOPTS
                Y(i) = Y(I)-ymean
                ENDDO
        endif


* (Calculate LSQ residuals SIG and CHI2 statistic)

        SIGu = 0.0D0
        SIG  = 0.0D0 ! weighted
        SWEGH= 0.0D0 !
	CHI2 = 0.0D0
        tmp=0.0d0
        DO I=1,NOPTS
		WEGHT=(WGM/WEG(I))**2
      	        if (.NOT.iweg) WEGHT=1.0D0 ! unweighted fit
                SWEGH=SWEGH+WEGHT
                SIGu=DBLE(Y(I))**2+SIGu
                SIG =(Y(I))**2*WEGHT + SIG 
                if(WEG(I).GT.1.d-06) CHI2=(Y(I)/WEG(I))**2+CHI2
        ENDDO
c
	SIGu=DSQRT(SIGu/DFLOAT(NOPTS-1))
	SIG =DSQRT(SIG /SWEGH)
        if(FAKET) goto 3003

* (Save right hand side of LSQ matrix which will be overriden in GAUSSD)

        DO 976 J=1,IDEGEQ
 976    RHS0(J)=RHS(J)

* (Solution of System Normal Equations by Gauss method.)
* (Vector of solutions is returned into COEF)

	if (NITER.NE.1) CALL GAUSSD(IDEGS,IDEGEQ,RHS,SUM,ONE,COEF)
        
* (Print PRECISE CORRELATION MATRIX)


        IF(KCNT.EQ.1.AND.iii.EQ.1.AND.NITER.NE.1) then
        if (MTRCR) then
        write(IO,*) ' Precise correlation matrix:'
	do i1=1,idegeq
	WRITE(IO,'(1x,I2,2h.),A7,1h:,26f16.12)') i1,WR(i1),
     &  (one(i1,j1)/dsqrt(dabs(one(i1,i1)*one(j1,j1))),j1=1,i1)
	enddo
	WRITE(IO,'(15x,26(i8,8x))') (i1,i1=1,idegeq)
        endif! (MTRCR)

* (Estimate accuracy of solutions and)
* (calculate correlations between parameters.)
        CALL ERR(IDEGS,IDEGEQ,NOPTS,COEF,ONE,ER,RHS0,D0,DISP)
        endif

* (Make sure we don't use elements we needn't)

        DO 3945 I=IDEGEQ+1,IDEGS
	ER(I)  = 0.0D0
 3945	COEF(I)= 0.0D0

* (Correction of parameters.)
*  -------------------------

 971    QM=1.0d0            ! set QM=1. usually, =0.9 for smooth iter.
        if (AKHILL) QM=0.9d0

* (Use steady corrections if need)

        DO 972 I=1,IDEGEQ
 972	COEF(I)=COEF(I)*QM

* (Correct values of parameters taken in previous iteration)

	  R0   = R00    + COEF(IE(1))*F(1)
        if (NJUMP.ne.0) then
        do k=1,NJUMP
        R0N(k)  = R0N0(k)+ COEF(IE(k+79))*F(1)
        enddo! (k)
        endif!(NJUMP.ne.0)

        do i=1,4
        F(i) = F0(i)  + COEF(IE(i+1))           
        enddo
        do i=5,9
        F(i) = F0(i)  + COEF(IE(55+i))
        enddo
        
	APSR  = APSR0  + COEF(IE(6))
	DPSR  = DPSR0  + COEF(IE(7))
        PMRA  = PMRA0  + COEF(IE(8))
        PMDEC = PMDEC0 + COEF(IE(9))
        PRLX  = PRLX0  + COEF(IE(10))
        
        DM(1) = DM0(1) + COEF(IE(11))
        do i=2,7
        DM(i) = DM0(i) + COEF(IE(68+i))
        enddo
        
       if (BINY) then
        do k=1,norb
        n=(k-1)*13
        X(k)    = X0(k)     + COEF(IE(12+n))
        Pb(k)   = Pb0(k)    + COEF(IE(13+n))
        e(k)    = e0(k)     + COEF(IE(14+n))
	RSG(k)  = RSG0(k)   + COEF(IE(15+n))
        Om(k)   = Om0(k)    + COEF(IE(16+n))
        RG(k)   = RG0(k)    + COEF(IE(17+n))
        Gc3(k)  = Gc30(k)   + COEF(IE(18+n))
        Si(k)   = Si0(k)    + COEF(IE(19+n))
        Omdot(k)= Omdot0(k) + COEF(IE(20+n))
        Xdot(k) = Xdot0(k)  + COEF(IE(21+n))
        Pbdot(k)= Pbdot0(k) + COEF(IE(22+n))
        Pbd2(k) = Pbd20(k)  + COEF(IE(101+n))  ! <<
        Xd2(k)  = Xd20(k)   + COEF(IE(102+n))  ! <<
        Pbd3(k) = Pbd30(k)  + COEF(IE(104+n))  ! <<        
        edot(k) = edot0(k)  + COEF(IE(23+n))
        dthet(k)= dthet0(k) + COEF(IE(24+n))
        enddo !(k=1,norb)
       endif !(BINY)

c=======================================================================
* (Output of data.)
* (Matrix of correlations.)

        IF(KCNT.EQ.1.AND.iii.EQ.1.AND.NITER.NE.1) then

        WRITE(IO,8988) 
 8988   FORMAT(/1X,80('-')
     &  /2X,' Matrix of correlations for fitted parameters:')
	DO 1002 IT=1,IDEGEQ
 1002   WRITE(IO,8989) IT,WR(IT),(-ONE(IT,KT),KT=1,IT)
        write(IO,'(14X,24(I2,2h.),2x))') (IT,IT=1,IDEGEQ)

* (Make correlations exceeding 0.95 labeled as '*')

        WRITE(IO,8990)

        DO 1003 IT=1,IDEGEQ
        do KT=1,IT
        if (DABS(ONE(IT,KT)).GE.0.950d0) then
             KS(KT:KT)='*'
             else
             KS(KT:KT)=':'
             endif
        enddo
        I=MOD(IT,10)
 1003   WRITE(IO,8991) IT,WR(IT),KS(1:IT)
        WRITE(IO,8992) (INT(MOD(IT,10)),IT=1,IDEGEQ)
        WRITE(IO,'(1X,80(1H-))/')
	write(*,'(/)')

* (Evaluate errors of adjustment parameters.)
* (Do it only once after first iteration)

* (Formal errors of parameters adjusted)

	ERCN =DABS(ER(IE(1)))
        ERP  =DABS(ER(IE(2)))/F(1)**2
        ERPT =DABS(ER(IE(3)))/F(1)**2
        ERPTT=DABS(ER(IE(4)))/F(1)**2
        ErF(1)=ER(IE(2))
        ErF(2)=ER(IE(3))
        ErF(3)=ER(IE(4))
        ErF(4)=ER(IE(5))
        ErF(5)=ER(IE(60))
        ErF(6)=ER(IE(61))
        ErF(7)=ER(IE(62))
        ErF(8)=ER(IE(63))
        ErF(9)=ER(IE(64))
        
        ERAL =ER(IE(6))*ArcRd/15.D0
        ERDL =ER(IE(7))*ArcRd
        ERALM=ER(IE(8))*ArcRd*Sday*365.25D0*1.D+03/DCOS(DPSR)
        ERDLM=ER(IE(9))*ArcRd*Sday*365.25D0*1.D+03
        ErPRLx=ER(IE(10))
        ERDM(1)=ER(IE(11))
        ERDM(2)=ER(IE(70))
        ERDM(3)=ER(IE(71))
        ERDM(4)=ER(IE(72))
        ERDM(5)=ER(IE(73))
        ERDM(6)=ER(IE(74))
        ERDM(7)=ER(IE(75))
        
        if (NJUMP.ne.0) then
        do k=1,NJUMP
        ErR0N(k)  = Er(IE(k+79))
        enddo! (k)
        endif!(NJUMP.ne.0)
        
       if (BINY) then
        do k=1,norb
        n=(k-1)*13
        ErX(k)  =ER(IE(12+n))
        ErPb(k) =ER(IE(13+n))
        Ere(k)  =ER(IE(14+n))
	ERSg(k) =ER(IE(15+n))
        ErOm(k) =ER(IE(16+n))
        ErRg(k) =ER(IE(17+n))
        ErGc3(k)=ER(IE(18+n))
        ErSi(k) =ER(IE(19+n))
        ErOmdot(k)  = ER(IE(20+n))*360.0/PI2*365.25*Sday   ! deg/yr
        ErXdot(k)   = ER(IE(21+n))*365.25*Sday*1.d+06      ! mcs/yr
        ErPbdot(k)  = ER(IE(22+n))*365.25*Sday*1.d+06*PI2  ! s/yr
        ErPbd2(k)   = ER(IE(101+n))                         ! <<
        ErXd2(k)    = ER(IE(102+n))                         ! <<
        ErPbd3(k)   = ER(IE(104+n))                         ! <<
        Eredot(k)   = ER(IE(23+n))                         ! 1/s
        Edthet(k)   = ER(IE(24+n))
        ErT0(k)     = ErSG(k)*Pb(k)/Sday
        enddo ! (k=1,norb)
       endif !(BINY)

       endif

* (Escape iterations if sigma oscillate near minima)     
* (This is point where iterative minimization will be stopped)

 3003   IF (SIG.GE.SIG0) GOTO 62
	SIG0=SIG

* (Next iteration step)

        III = III + 1

* (Display n.step of iteration and statistics)
        
        WRITE(IO,64) III-2,SIG,CHI2,NOPTS-IDEGEQ,SIGu
        WRITE(*,64) III-2,SIG,CHI2,NOPTS-IDEGEQ,SIGu

* (Escape iterations if the amount of loops exceeds maximal no.NITER)

        IF(III.GT.NITER) GOTO 62

* (Return back to make new iteration with new parameters)
*  -----------------------------------------------------
                GOTO 1001
 62		CONTINUE


* (------------------------------------------------------------------)
* (          Here is the end of LSQ minimization procedure           )
* (                 Below is output of data only                     )
* (------------------------------------------------------------------)

	if (KCNT.EQ.1) R001=R0

* (Use these data if you need to check CHI2 surface near global minima)
*        WRITE(IO,64) III-2,SIGu,CHI2,NOPTS-IDEGEQ

 3891   CONTINUE

c=======================================================================
* (Print results)
* --------------


* (Show run-time of program)

	call gttime(ardt)
	ih =ardt(1)
	imm=ardt(2)
	is =ardt(3)

        write(IO,'(/2X,10HTime End.=,2(I2,1H:)I2)')ih,imm,is
        imm=float(ih-ih0)*60.0 + float(imm-imm0)
        SEC= is-is0
        if (SEC.LT.0.0) then
                SEC=SEC+60.0
                imm=imm-1
                endif
       write(IO,'(2X,11HTotal Time=,I6,5H min.,F10.6,5H sec.)')imm,SEC

        if(FAKET) go to 3004

* (Reduce parameter's values to initial epoch)
* (and convert units).

        CALL RADEG(APSR,1,APS,WrsA)
        CALL RADEG(DPSR,0,DPS,WrsD)
        PMRA = PMRA*ArcRd*Sday*365.25D0*1.D+03
        PMDEC = PMDEC*ArcRd*Sday*365.25D0*1.D+03
        P0   =  1.D0/F(1)
	P1   = -F(2)/F(1)/F(1)
	P2   = -F(3)/F(1)/F(1)+F(2)**2*2.D0/F(1)/F(1)/F(1)
        if (F(3).EQ.0.0d0) P2=0.0d0
        Pepoch = Pepoch-2400000.5d0

        if (PRLX.NE.0.0) then
        DIST= 1.20962D-03/PRLX
        DISTMIN=1.20962D-03/(PRLX+ErPRLX)
        DISTMAX=1.20962D-03/(PRLX-ErPRLX)
        endif

        PRLX  =PRLX/1.20962d-03    *1.d+03
        ErPRLX=ErPRLX/1.20962d-03  *1.d+03

* (Print Astrometric and Spin parameters)

        WRITE(IO,74) R0/f(1),ERCN,WrsA,ERAL,WrsD,ERDL,
     &  PMRA,ERALM,PMDEC,ERDLM,P0,ERP,P1,ERPT,
     &  P2,ERPTT,
     &  Pepoch,DM(1),ErDM(1),Prlx,ErPrlx,DIST,DISTMIN,DISTMAX
        
* (spin and spin derivatives)
        if(F(1).ne.0.0d0) write(io,'(8x,6Hf0   :,1x,2(1pd24.16),1x,
     &  3H1/s)') f(1),erf(1)
        do i=2,9
        write(pp,'(i1)') i-1
        write(pw,'(i1)') i
        if(F(i).ne.0.0d0) write(io,'(8x,1Hf,a1,4Hdot:,1x,
     &  2(1pd24.12),1x,3H1/s,a1)') pp,f(i),erf(i),pw
        enddo ! i=2,9
* (DM derivatives )
        do i=2,7
        write(pp,'(i1)') i-1
        if(DM(i).ne.0.0d0)write(io,'(8x,3HDMd,a1,2H :,1x,2(1pd24.12),1x,
     &   9Hpc*cm3/yr,a1)') pp,dm(i)*Sday**(i-1)*365.25**(i-1),
     &   erdm(i)*Sday**(i-1)*365.25**(i-1),pp
        enddo ! i=2,7

*(Create parameter file _tim.inp if requested <ptimer -i inputfilename>)
************************************************************************
        if (MKINP) then
* (Read header of the parameter file)
        open(unit=10,file='_tim.inp',err=9007)
        goto 9008
 9007   write(*,'(2x,27H! _tim.inp already exist !!)')
        MKINP=.false.
        goto 9009
 9008   write(10,'(/a)') 'HEAD'
        write(10,*) ' '
        write(10,'(2a)') 'PSR            ',psrname
        pp='0'
        if (ErAl.ne.0.0) pp='1'
        write(10,'(a2,6x,a17,3x,a1,f20.8)') 'RA',WrsA,pp,ErAl
        pp='0'
        if (ErDl.ne.0.0) pp='1'
        write(10,'(a3,5x,a17,3x,a1,f20.8)') 'DEC',WrsD,pp,ErDl
        pp='0'
        if (ErAlm.ne.0.0) pp='1'
        write(10,'(a4,f21.4,3x,a1,f20.4)')  'PMRA',PmRA,pp,ErAlm
        pp='0'
        if (ErDLM.ne.0.0) pp='1'
        write(10,'(a5,f20.4,3x,a1,f20.4)')  'PMDEC',PmDEC,pp,ErDLM
        if (PRLX.GT.0.00001) then
        pp='0'
        if (ErPRLX.ne.0.0) pp='1'
        write(10,'(a2,f23.4,3x,a1,f20.4)')  'PX',PRLX,pp,ErPRLX
        endif

        pp='0'
        if (Erf(1).ne.0.0) pp='1'
        if(f(1).ne.0.0) write(10,'(a2,f24.16,2x,a1,f20.16)')
     &   'F0',f(1),pp,Erf(1)
        do i=2,9
        write(pw,'(i1)') i-1
        pp='0'
        if (Erf(i).ne.0.0) pp='1'
        if(f(i).ne.0.0) write(10,'(2a1,1pd24.12,2x,a1,1pd20.12)')
     &   'F',pw,f(i),pp,Erf(i)
        enddo ! i=2,9 
        write(10,'(6HPEPOCH,f20.3)') Pepoch
        write(10,'(5HSTART,f21.3)')  TOAmin-2400000.5d0
        write(10,'(6HFINISH,f20.3)') TOAmax-2400000.5d0
        pp='0'
        if (ErDM(1).ne.0.0) pp='1'
        write(10,'(2HDM,f24.6,2x,a,f20.6)') DM(1),pp,ErDM(1)
        do i=2,7
        write(pw,'(i1)') i-1
        pp='0'
        if (ErDM(i).ne.0.0) pp='1'
        if(DM(i).ne.0.0) write(10,'(a2,a1,1pd23.12,2x,a1,1pd20.12)')
     &   'DM',pw,DM(i)*Sday**(i-1)*365.25**(i-1),pp,
     &   ErDM(i)*Sday**(i-1)*365.25**(i-1)
        enddo ! i=2,7

        write(10,'(5HEPHEM,13x,5HDE200)')
        write(10,'(3HCLK,15x,9HUTC(NIST))')
        write(10,'(4HNITS,14x,1H3)')
        write(10,*) ' '
         if(BINY) then
         ww=' '
         if(mDADE) ww='DD'
         if(mBT)   ww='BT'
         if(mBTP)  ww='BT+'
         if(mBTPP) ww='BT++'
         if(mDADEP) ww='DT++'
         write(10,'(6HBINARY,12x,a)') ww
         endif ! (BINY)
        endif ! (MKINP)
************************************************************************
        
 9009   continue
       if (BINY) then
        do k=1,norb
	if (norb.gt.1) write(IO,'(/1x,i3,1x,5HOrbit)') k
        Pb(k)=Pb(k) * PI2
        ErPb(k)=ErPb(k) * PI2
        Pbdot(k)=Pbdot(k) * PI2
        Om(k)=Om(k)*360.0d0/PI2
                if (Om(k).LT.0.0d0) Om(k)=Om(k)+360.0d0
        ErOm(k)  = ErOm(k)*360.0/PI2
        Gc3(k)   = Gc3(k)*1.d+06                   ! mcs
	ErGc3(k) = ErGc3(k)*1.d+06
        RG(k)    = RG(k)*1.d+06                    ! mcs
	ErRG(k)  = ErRG(k)*1.d+06
        Omdot(k) = Omdot(k)*360.0d0/PI2*365.25*Sday    ! deg/yr
        Xdot(k)  = Xdot(k)*365.25*Sday*1.d+06          ! mcs/yr
        Pbdot(k) = Pbdot(k)*365.25*Sday*1.d+06         ! mcs/yr
        edot(k)  = edot(k)                             ! 1/s
        T0(k)    = T0(k)-RSG(k)*Pb(k)/PI2/Sday

* (Print binary parameters)

        WRITE(IO,75) Pb(k),ErPb(k),X(k),ErX(k),e(k),ErE(k),
     &  Om(k),ErOm(k),T0(k),ErT0(k),
     &  RG(k),ErRG(k),
     &  Omdot(k),ErOmdot(k),
     +  Omdot(k)/360.0*PI2/365.25/Sday,ErOmdot(k)/360.0*PI2/365.25/Sday,
     &  Pbdot(k),ErPbdot(k),
     +  Pbdot(k)/365.25/Sday*1.d-06,ErPbdot(k)/365.25/Sday*1.d-06,
     &  Xdot(k),ErXdot(k),
     +  Xdot(k)/365.25/Sday*1.d-06,ErXdot(k)/365.25/Sday*1.d-06,
     &  edot(k)*365.25*Sday,Eredot(k)*365.25*Sday,
     +  edot(k),Eredot(k),        
     &  Gc3(k),ErGc3(k),Si(k),
     &  nint(asin(Si(k))/PI2*360.0),ErSi(k),dthet(k),edthet(k)

        if (Pbd2(k).ne.0.0d0) write(io,'(1x,5HPbd2=,2(1pg24.12),
     &  1x,4H1/s2)')Pbd2(k),Erpbd2(k)
        if (Pbd3(k).ne.0.0d0) write(io,'(1x,5HPbd3=,2(1pg24.12),
     &  1x,4H1/s3)')Pbd3(k),Erpbd3(k)        
        if (Xd2(k).ne.0.0d0) write(io,'(1x,5H Xd2=,2(1pg24.12),
     &  1x,4H1/s2)') Xd2(k),ErXd2(k)

* (Add binary to parameter file <ptimer -i ...>)
************************************************************************
        if (MKINP) then
        pp='0'
        if (ErX(k).ne.0.0) pp='1'
        write(10,'(3HA1_,i1,f22.10,2x,a,f20.10)') k,X(k),pp,ErX(k)
        pp='0'
        if (ErE(k).ne.0.0) pp='1'
        write(10,'(2HE_,i1,f23.10,2x,a,f20.10)') k,E(k),pp,ErE(k)
        pp='0'
        if (ErT0(k).ne.0.0) pp='1'
        write(10,'(3HT0_,i1,f22.10,2x,a,f20.10)') k,T0(k),pp,ErT0(k)
        pp='0'
        if (ErPB(k).ne.0.0) pp='1'
        write(10,'(3HPB_,i1,f22.12,2x,a,f20.12)') k,PB(k)/Sday,pp,
     &   ErPB(k)/Sday
        pp='0'
        if (ErOM(k).ne.0.0) pp='1'
        write(10,'(3HOM_,i1,f22.10,2x,a,f20.10)') k,OM(k),pp,ErOM(k)
*(Xdot)
        if (Xdot(k).ne.0.0) then
        Xdot(k)=Xdot(k)/365.25/Sday*1.d+06
        ErXdot(k)=ErXdot(k)/365.25/Sday*1.d+06
        pp='0'
        if (ErXdot(k).ne.0.0) pp='1'
        write(10,'(5HXDOT_,i1,f20.10,2x,a,f20.10)') k,Xdot(k),
     &   pp,ErXdot(k)
        endif ! (Xdot)
*(Omdot)
        if (OMdot(k).ne.0.0) then
        pp='0'
        if (ErOMdot(k).ne.0.0) pp='1'
        write(10,'(6HOMDOT_,i1,f19.10,2x,a,f20.10)') k,OMdot(k),
     &   pp,ErOMdot(k)
        endif ! (OMdot)
*(Pbdot)        
        if (PBdot(k).ne.0.0) then
        PBdot(k)=PBdot(k)/365.25/Sday*1.d+06
        ErPBdot(k)=ErPBdot(k)/365.25/Sday*1.d+06
        pp='0'
        if (ErPBdot(k).ne.0.0) pp='1'
        write(10,'(6HPBDOT_,i1,f19.10,2x,a,f20.10)') k,PBdot(k),
     &   pp,ErPBdot(k)
        endif ! (PBdot)
*(Gamma)
        if (RG(k).ne.0.0) then
        RG(k)=RG(k)
        ErRG(k)=ErRG(k)
        pp='0'
        if (ErRG(k).ne.0.0) pp='1'
        write(10,'(6HGAMMA_,i1,f20.10,1x,a,f20.10)') k,RG(k),pp,ErRG(k)
        endif ! (Gamma)
*(r,s)
        if (Gc3(k).ne.0.0) then
        Gc3(k)=Gc3(k)/4.92549d0
        ErGc3(k)=ErGc3(k)/4.92549d0
        pp='0'
        if (ErGc3(k).ne.0.0) pp='1'
        write(10,'(3HM2_,i1,f20.10,4x,a,f20.10)') k,Gc3(k),pp,ErGc3(k)
        pp='0'
        if (ErSI(k).ne.0.0) pp='1'
        write(10,'(3HSI_,i1,f20.10,4x,a,f20.10)') k,SI(k),pp,ErSI(k)
        endif ! (M2)
        endif ! (MKINP)
************************************************************************
        enddo !(k=1,norb)       
        endif !(BINY)
        if (MKINP) close(10)
        
************************************************************************
* (Print TOA JUMPS)

        if(NJUMP.ne.0) then
        write(IO,*) ' '
        write(IO,*) '-- TOA offsets:'
        do k=1,NJUMP
        write(IO,'(1x,i2,1x,7HOffset=,f14.3,1x,1h(,f14.3,
     &  1h),1x,3Hmcs)') k,(R0N(k)-R0)*1.0d6/f(1),ErR0N(k)*1.d6
        enddo !(k)
        endif ! (NJUMP)
        
c=======================================================================
* (DEDUCE additional parameters)
*  ----------------------------

        AGE=0.0
        if(P1.NE.0.0) then
        AGE=(0.5*P0/P1/SDAY/365.25) ! char.age
        WRITE(IO,80) AGE
 80     FORMAT(/80('-')/
     &  9('.'),' DEDUCED Pulsar''s parameters:'/
     &  ' Char.age        tau=',1pg10.3,' yrs')
        if(P2.NE.0.0) write (IO,'(1x,19HBreak-down index n=,g12.3)')
     &  2.+P0*P2/P1**2
        endif

       if (BINY) then
        do k=1,norb
        mfunct = (Spdl*X(k))**3/(Pb(k)/PI2)**2 / 1.327d+20
        Write(IO,'(1x,19HMass function  fm= ,1pg19.7,2x,4HMsun)') mfunct

*(evaluate companion mass Mc for Mp=1.4Msun through 60 iterations)

        mc=0.0d0
	do i=1,60
	mc=mfunct**(1./3.)*(mc+1.4)**(2./3.)
	enddo ! i
        Write(IO,'(1x,19HCompanion mass Mc= ,1pg20.9,1x,4HMsun,
     &  24H/sin(i) (for Mp=1.4Msun))') Mc

        if (Si(k)*Gc3(k).ne.0.0) then
        Write(IO,'(1x,28Hfrom range and shape (r,s): )')

        ErOm(k)=Gc3(k)*1.d-06 * Spdl**3 / 1.327d+20
        tmp=ErOm(k)-ErGc3(k)*1.d-06 * Spdl**3 / 1.327d+20  ! min
        tmp1=ErOm(k)+ErGc3(k)*1.d-06 * Spdl**3 / 1.327d+20 ! max

        Write(IO,'(1x,19HCompanion     Mc = ,g20.9,1x,4HMsun,
     &  1H(,2g20.9,1H))')
     &  ErOm(k),tmp,tmp1
        Write(IO,'(1x,19HPulsar        Mp = ,g20.9,1x,4HMsun,
     &  1H(,2g20.9,1H))')
     &  dsqrt(dabs((Si(k)*ErOm(k))**3/mfunct-ErOm(k))),
     &  dsqrt(dabs((Si(k)*tmp)**3/mfunct-tmp)),
     &  dsqrt(dabs((Si(k)*tmp1)**3/mfunct-tmp1))
        ErOm(k)=-Erom(k)
        endif

*(parameters expected in General Relativity)
* -----------------------------------------

        Pb(k)=Pb(k)/PI2
        WRITE(IO,81)
 81     FORMAT(/9('.'),1x,
     &  'EXPECTED in General Relativity (for i=90,30 deg))')

        do i=1,2
        if (i.eq.2) Mc=Mc*2.0
        WRITE(IO,82) Mc
 82     FORMAT(/1x,'(Mp=1.40,Mc=',f5.2,')')

*(GR periastron advance)

        tmp=3.0d0/Pb(k)**(5.0/3.0)*(4.92549e-6*(1.4+Mc))
     &  **(2.0/3.0) / (1.0-e(k)**2)
        Write(IO,
     &  '(1x,20HPeriast.adv.  dW/dt=,f18.9,1x,
     &  7Hdeg/yr=,1pe9.2,1x,3H1/s)')
     &  tmp/PI2*360.0*Sday*365.25,tmp

*(GR Einstein delay Gamma)

        tmp=e(k)*Pb(k)**(1.0/3.0)*(4.92549e-6**(2.0/3.0))/
     &  (1.4+Mc)**(4.0/3.0)*Mc*(1.4+2.0*Mc)
        Write(IO,
     &  '(1x,20HEinstein del.Gamma =,f14.3,1x,3Hmcs)') tmp*1.0e+6

*(GR Orbital Period decay)

	tmp=-19.2*PI2*(4.92549e-6/Pb(k))**(5.0/3.0)*1.4*Mc/
     &  (1.4+Mc)**(1.0/3.0) /(1.0-e(k)**2)**(7.0/2.0)*
     &  (1.0+73.0/24.0*e(k)**2+37.0/96.0*e(k)**4)
	Write(IO,
     &  '(1x,20HOrb.per.decay dP/dt=,3x,1pe9.2,1x,
     &  1H=,0pf9.3,1x,6Hmcs/yr)')
     &  tmp,tmp*Sday*365.25*1.0e+6
        enddo ! (i=1,2)

*(GR Range of Shapiro delay)

	Write(IO,
     &  '(1x,20HRange of Shapiro r=,f9.3,7H/sin(i),1x,3Hmcs,
     &  /80(1H-))') Mc/2.0*4.92549

        enddo !(k=1,norb)
	endif !(BINY)

	CLOSE(unit= 9)
        
c=======================================================================
* (Write residuals and data for plot into files)
*  ------------------------------------------

	if(WRITRS) then
        write(IO,'(//2x,37HFile _res.out:(Date)vs(Residuals TOA))')

* (file _tim.out for UTC TOAs, BSS TOAs,residuals)
                OPEN(UNIT=4,FILE='_tim.out')
                PNM(1:7)=psrname(1:7)
                CALL OUT(NOPTS,JDD,TDBi,RM,Y,WEG,PNM,4)
		CLOSE(4)

* (file _res.out for (observed epoch,days) vs (residuals,mcs)
                open(4, file='_res.out')
                open(5, file='_kep.out')
		do i=1,nopts
*                write(4,'(f21.11,1x,g15.6,1x,g8.2)')
                write(4,'(f22.12,1x,g15.6,1x,g8.2)')         !AEA
     &          (JDD(i)-2400000.5d0)+TDBi(i)/Sday,Y(i),WEG(i)
                write(5,'(f21.11,1x,g15.6)')
     &          (JDD(i)-2400000.5d0)+TDBi(i)/Sday,U(i,1)*1.d6
		enddo
		close(4)
		close(5)

* (file _scl.res in input format for program to build time scale)

                open(4, file='_scl.res')
		do i=1,nopts
                write(4,'(f12.6,1x,es20.9,1x,es10.2)')         
     &          (JDD(i)-2400000.5d0-PEPOCH)+TDBi(i)/Sday,Y(i)/1000000.,
     &          WEG(i)/1000000.
          enddo
		close(4)


* (files _kep.orb, _res.orb, _res.shp for binaries)

        if (BINY) then
        write(IO,
     &  '(2x,47H_kep.orb:(Orb.ph)vs(Orb.Roemer[mcs])-predicted))')
        write(IO,'(2x,30H_res.orb:(Orb.ph)vs(Residuals))')
        write(IO,'(2x,34H_res.shp:(Orb.ph)vs(Shapiro [mcs]))')
        do k=1,norb
        Pb(k)=Pb(k)*PI2
        write(pp,'(i1)') k
                open(4, file='_res'//pp//'.orb')
                open(7, file='_res'//pp//'.shp')
                open(9, file='_kep'//pp//'.orb')
		do i=1,nopts
                TMP=(JDD(i)-T0(k))*Sday + TDBi(i)
                TMP=(TMP-DINT(TMP/Pb(k))*Pb(k))/Pb(k)
		if (TMP.LT.0.d0) TMP=TMP+1.0d0
                write(4,'(f15.9,g15.6,1x,g10.2)') TMP ,Y(i),WEG(i)
                write(9,'(f21.11,g17.6)') TMP,U(i,1)
                write(7,'(f15.9,g17.7,1x,g10.2)') TMP ,U(i,2),WEG(i)
		enddo
                close(9)
		close(7)
		close(4)
        enddo !(k=1,norb)
        endif !(BINY)
	endif


*
c=======================================================================
* (Write fake TOAs to file _fake.toa)

 3004   if(FAKET) then
        st80='Fake TOAS in file _fake.toa'
        write(IO,'(a80)') st80
        write(*,*) '  ',st80
        open(122,file='_fake.toa')
        do i=1,nopts
        JDint=dint(jdd(i)-2400000.5)
c        RM(i) = RM(i) + rndm(i*2+3)*WegRMS*1.d-06 !!! add random
         RM(i) = RM(i) + rndm(i)*WegRMS*1.d-06 !!! add random noise
c        RM(i) = RM(i) + Rand(i) * WegRMS*1.d-06 !!! add random noise

        TMP=RM(i)/Sday
           if (TMP.lt.0.0d0) then
           TMP=1+TMP
           JDint=JDint-1
           endif! RM(i)
           if (TMP.gt.1.0d0) then
           TMP=TMP-1.0d0
           JDint=JDint+1
           endif! RM(i)
        write(122,'(1x,a8,2x,i5,f14.13,f9.2,f12.5,f9.6,1x,a2)')
     & 'fakeTOAs',
     & JDint,TMP,WEG(i),FR0D(i),dDM(i),ChSIT
        enddo
        close(122)
        endif
c=======================================================================
        CLOSE(unit=IO)
	STOP ' '
*-----------------------------------------------------------------------
 5555   write(*,1) confile
        STOP ' '
*-----------------------------------------------------------------------
*  1     FORMAT(30X,'PTIMER v1.00' /80('-')/
*     &1x,'Pulsar Timing Program'/
*     &1X,'Positions for J2000.0;'/1x,'Planetary ephemeries DE200/'
*     &'LE200, DE405/LE405, EPM2015'
*     & /1x,'Time ephemeries: Fairhead L.& Bretagnon P.'/
*     &' Parameters read from ',a
*     & /80('-')/
*     &1X,'Timing model:Doroshenko O. & Kopeikin S.'/
*     &1X,'Binary timing model DD,DD+,BT, BT+'/
*     & 2X,' PTIMER <options> ... <file>'/
*     &19X,'<file>= pHHMMxDD with PSR name (f.eg. -p1937+21)'/
*     &19X,' (need 1937p21 TOA file and 1937p21.p par.file)'/
*     &19X,'<file>= file name (f.eg.1937.tpo) for free-format input'/
*     &9X,' <options>  - optional, are:'/
*     &19X,'-c -a priory values of PSR parameters will not be displayed'/
*     &19X,'-r -the mean offset of TOA residuals will not be removed'/
*     &19X,'-i -create input parameter file _tim.inp with fitted values'/
*     &19X,'-o -show orbital phase in _obsp.dat (in period gener.mode)'/
*     &19X,'-m -print precise correlation matrix'/
*     &19X,'-z -iterate with factor 0.9'/
*     &19X,'-h -show this comment'/
*     &19X,'-d<path> path to ptimer.cfg, changes default path'/
*     &19X,'-f fake TOAs generation mode'/
*     &10X,'Example: ptimer 1937.toa'/
*     &80('-'))


*-----------------------------------------------------------------------
  1     FORMAT(35X,'PTIMER v1.1'/80('-')/

     &1x,'Программа хронометрирования и уточнения параметров пульсаров'/
     &1X,'Эфемериды DE200/LE200, DE405/LE405, EPM2015'/
     &1X,'Преобразование шкал времени: Fairhead L.& Bretagnon P.'/
     &1X,'Модель хронометрирования : Doroshenko O. & Kopeikin S.'/
     &1X,'Модели двойного пульсара: DD,DD+,BT, BT+'/
     &' По умолчанию параметры считываются из: ',a
     & /80('-')/
     &2X,' ptimer <options> ... <file>'/
     &19X,'<file>=имя пульсара в формате pHHMMxDD (eg. -p1937+21)'/
     &19X,' (исп. файл МПИ 1937p21 и файл параметров 1937p21.p'/
     &19X,' фиксированной длины - формат TIMAPR1 )'/
     &19X,'<file>= имя файла (eg. 1937.tpo) - для свободного формата'/
     &19X,'(формат TIMAPR410/расширенный TEMPO2)'/
     &9X,' <options>  - опции:'/
     &19X,'-c -a не показывать начальные параметры'/
     &19X,'-r -не вычитать среднее значение ОУ МПИ'/
     &19X,'-i -вывод уточненных параметров в файл _tim.inp'/
*     &19X,'-o -вывод орб. фазы  в _obsp.dat(режим генерации периодов)'/
     &19X,'-m -вывод точной корреляционной матрицы'/
*     &19X,'-z -иетрация с фактором 0.9 (мягкая)'/
     &19X,'-h -показать help'/
     &19X,'-d<path> путь к файлу ptimer.cfg (вместо пути по умолчанию)'/
     &19X,'-f генерация расчетных МПИ из заданных параметров пульсара)'/ 
     &10X,'Примеры использования : ptimer 1937.toa; ptimer -p1913+16'/
     &80('-'))



 10     FORMAT(a5/2(F13.9/),1X,1PD23.16/2(1X,
     &  1PD17.10/),1X,5PD17.10,4(/1X,2PD12.5))
 11	FORMAT(//2(2PD15.6/),2X,5PD14.5)
 12	FORMAT(1X,/,1X,7PD19.8,1X,/,1X,5PD20.10)
 17	FORMAT(I5,f13.6,I5)
 20	FORMAT(/5X,I6/f11.7/)
c 21	FORMAT(/(2X,2F20.8))
 22	FORMAT(2X,I5,2F20.2)
 33	FORMAT(2X,'Input PSR - name [1919+21]                   >',$)
 64     FORMAT(I2,' sig=',
     &  F12.3,' us; Chi2=',f18.3,'(n=',i7,') Sig(unw)=',f12.3) ! PC
c64     FORMAT(I2,' Iterat., TOA Res.=',
c    &  F13.4,' mcs; Chi2=',f18.3,'(n=',i7,') Sig(unw)=',f13.4)
 97	FORMAT(1X,1PD23.15)
 87     FORMAT( 1x,':First TOA      :',2(I2,'.'),I4,'(JD=',f10.1,')')
 88     FORMAT( 1x,':Last TOA       :',2(I2,'.'),I4,'(JD=',f10.1,')'/)
 209    FORMAT(/1x,':Number of TOAs :',I6)
 74     FORMAT(/80('-')/
     &  9('.'),1x,'ADJUSTED Pulsar''s parameters:'/
     &  ' ** Spin parameters:'/
     &  ' Const       :',g21.8,'+/-',1pg8.1 's'/
     &  ' Alpha       :',A17,3x,1pg12.1,'s'/
     &  ' Delta       :',A17,3x,1pg12.1,'" '/
     &  ' M-Alpha     :',g12.5,8x,1pg12.1,' mas/yr'/
     &  ' M-Delta     :',g12.5,8x,1pg12.1,' mas/yr'/
     &  ' P           :',1Pg22.15,1pg10.1,' s'/
     &  ' dP/dt       :',1Pg22.15,1pg10.1/
     &  ' d2P/dt2     :',1Pg22.15,1pg10.1,' s-1'/
     &  ' MJD(P)      :',g22.12/
     &  ' DM          :',g22.9,1pg9.1,' pc cm-3'/
     &  ' Parallax    :',g22.5,1pg11.3,' mas'/
     &  ' =Distance   :',g22.6,' pc (',g12.3,'...',g12.3,' pc)')
 75     FORMAT(/
     &  ' ** Orbital parameters:'/
     &  ' * KEPLERIAN:',5('*')/
     &  ' Pb,     binary Period:',g25.17,2x,'+/-',1pg9.1,' s'/
     &  ' (a sini), proj.sem.ax:',g25.17,2x,1pg12.1,' s'/
     &  ' e,       eccentricity:',g25.15,2x,1pg12.1/
     &  ' OM,   long.periastron:',g25.15,2x,1pg12.1,' deg'/
     &  ' To,   time periastron:',g25.16,2x,1pg12.1,' JD'/
     &  ' * POST-KEPLERIAN:',5('*')/
     &  ' Gamma  Gr.redsh+Dopl2:',g23.9,4x,1pg12.1,' mcs'/
     &  ' dOM/dt periastron adv:',g23.9,4x,1pg12.1,' deg/yr =',
     +   g22.5,1pg9.1,' 1/s'/
     &  ' d(Pb)/dt     Pb-decay:',g23.9,4x,1pg12.1,' mcs/yr =',
     +   g22.5,1pg9.1,' s/s'/
     &  ' d(a sini)/dt  x-decay:',g23.9,4x,1pg12.1,' mcs/yr =',
     +   g22.5,1pg9.1,' s/s'/
     &  ' de/dt         e-decay:',g23.9,4x,1pg12.1,' 1/yr   =',
     +   g22.5,1pg9.1,' 1/s'/
     &  ' r,              range:',g23.9,4x,1pg12.1,' mcs'/
     &  ' s = sin(i)      shape:',g23.9,'(',i2,'deg)',1x,1pg8.1/
     &  ' dTHETA               :',g23.9,4x,1pg12.1,' ns')
 76     FORMAT(80('-')/
     &  9('.'),1x,'ASSUMED Pulsar parameters:'/
     &  ' ** SPIN PARAMETERS:'/
     &  ' Alpha(2000) :',A17,'s'/
     &  ' Delta(2000) :',A17,'"'/
     &  ' Alpha-Motion:',g22.7,' mas/yr'/
     &  ' Delta-Motion:',g22.7,' mas/yr'/
     &  ' P           :',1pg22.15,' s'/
     &  ' dP/dt       :',1pg22.15/
     &  ' d2P/dt2     :',1pg22.5,' s-2'/
     &  ' MJD(P)      :',g22.12/
     &  ' DM          :',g22.9,' pc cm-3'/
     &  ' Distance    :',g22.9,' pc')
 77       FORMAT(/
     &  ' ** Orbital parameters:'/
     &  ' * KEPLERIAN:',5('*')/
     &  ' Pb,     Binary Period:',g24.18,' s'/
     &  ' (a sini), Proj.sem.ax:',g24.12,' s'/
     &  ' e,       Eccentricity:',g24.9/
     &  ' OM,   Long.periastron:',g24.14,' deg'/
     &  ' To,   Time periastron:',g24.14,' JD'/
     &  ' * POST-KEPLERIAN:',5('*')/
     &  ' Gamma  Gr.redsh+Dopl2:',g23.9,' mcs'/
     &  ' OMdot Periastron adv.:',g23.9,' deg/yr'/
     &  ' dPb/dt,      Pb-decay:',g23.9,' mcs/yr'/
     &  ' r,              range:',g23.9,' mcs'/
     &  ' s = sin(i)           :',g23.9)

 8989   FORMAT(1X,I2,'.)',A7,'|',24(f6.3))
 8992   FORMAT(12x,250I1,/)
 8990   FORMAT(/1x,80('-')/2x,'Correlations exceeding 0.95=','''*''')
 8991   FORMAT(1x,I2,1x,A7,')',40A)
	stop
 9004	write(*,'(2x,23H! There is no file TOAs)')
	stop
 9005	write(*,'(2x,18H! There is no file,A50)') Filpar
        write(*,1) confile
	END

*========================================================================

        SUBROUTINE LIN(NJUMPN,NORB)

* Builds matrix SUM[i,j] of normal equation system
* COEF[i]*SUM[i,j]=RHS[i]
* from partial derivatives of adjusted parameters.

        implicit none
        INTEGER*4 IDEGS

        PARAMETER(IDEGS=224)

        INTEGER*2 NJUMPN
        REAL*8 DTB,FRP,RE,AE,DE,AP,DP,DTI,Weg,C,SUM,RHS,D0,TA,DTII,
     &     VPR,AF(3),BT(3),RigKT(3),e(3),Znn(3),X(3),Gc3(3),
     &     WW,eSQ(3),Su(3),Cu(3),ScDM,EXC(3),DTPER(3),Ane(3),PB(3),
     &     cosdp,sindp,cosde,sinde,sinape,cosape,Cue
	INTEGER*4 IE,I,J,IDEGEQ,NORB,N,K
        LOGICAL BINY,mDADE,mBTPP,mDADEP,mBT
	COMMON/CETBL6/C(IDEGS),SUM(IDEGS,IDEGS),RHS(IDEGS),D0,IE(IDEGS)
	COMMON/CETBL7a/BINY,mDADE,mBTPP,mDADEP,mBT
        COMMON/CETBL7b/IDEGEQ
        COMMON/CETBL7c/DTB,DTPER,FRP,RE,AE,DE,AP,DP,DTI,Weg,
     &     VPR,AF,BT,RigKT,Cu,Su,e,Znn,X,Gc3,esQ,ScDM,Ane,EXC,PB

* notations here:
*  EXC=u
*  esQ=DSQRT(1.d0-e**2)
*  Su = DSIN(EXC)      ! = sin(Ecc.anom)
*  Cu = DCOS(EXC)      ! = cos(Ecc.anom)
*  Pb = Pb/2pi

        TA   =  499.004783703D0
        DTII = -DTI
        cosdp=  DCOS(DP)
        sindp=  DSIN(DP)
        cosde=  DCOS(DE)
        sinde=  DSIN(DE)
        sinape= DSIN(AP-AE)
        cosape= DCOS(AP-AE)

* (Calculation of partial derivatives of fitted param.)

**-Epoch
        if(NJUMPN.EQ.0) then
            C(IE(1)) =  1.0D0
        else
            C(IE(1)) = 0.0d0
            do I=80,IDEGS
            C(IE(I)) = 0.0d0
            enddo !(I)
* (JUMP Pointers are from 80 )
            C(IE(79+NJUMPN))=1.0d0
        endif

**-Frequency
        C(IE(2)) = DTB/FRP
**-First deriv. of Freq.
        C(IE(3)) = C(IE(2))/2.D0*DTB
**-Second deriv. of Freq.
        C(IE(4)) = C(IE(3))/3.D0*DTB
**-Third deriv. of Freq.
        C(IE(5)) = C(IE(4))/4.D0*DTB
**-4th deriv. of Freq.
        C(IE(60)) = C(IE(5))/5.D0*DTB
**-5th deriv. of Freq.
        C(IE(61)) = C(IE(60))/6.D0*DTB
**-6th deriv. of Freq.
        C(IE(62)) = C(IE(61))/7.D0*DTB
**-7th deriv. of Freq.
        C(IE(63)) = C(IE(62))/8.D0*DTB
**-8th deriv. of Freq.
        C(IE(64)) = C(IE(63))/9.D0*DTB
        
**-Right accention
        C(IE(6)) = -RE*TA*COSDP*SINAPE*COSDE
**-Declination
	C(IE(7)) = -RE*(SINDP*COSAPE*COSDE-SINDE*COSDP)*TA
**-(Al)-proper motion
        C(IE(8)) = C(IE(6))/COSDP*DTB
**-(Dl)-proper motion
	C(IE(9)) = C(IE(7))*DTB
**-Parallax
        C(IE(10)) = -VPR
**-DM
        C(IE(11))= -ScDM
**-DMd
        C(IE(70))= -ScDM*DTB
**-DMd2 
        C(IE(71))=C(IE(70))*DTB
**-DMd3 
        C(IE(72))=C(IE(71))*DTB
**-DMd4 
        C(IE(73))=C(IE(72))*DTB
**-DMd5
        C(IE(74))=C(IE(73))*DTB
**-DMd6
        C(IE(75))=C(IE(74))*DTB

* (for binaries)

        if (BINY) then
        do k=1,norb
        n = (k-1)*13
        Cue=Cu(k)-e(k)
        WW=(AF(k)*Su(k)-BT(k)*Cu(k))/(1.0d0-e(k)*Cu(k))
**-X
        C(IE(12+n))=-(AF(k)*Cue+BT(k)*Su(k))/X(k)
**-Pb
        C(IE(13+n))=-WW*RigKT(k)
**-e
        C(IE(14+n))= (WW*Su(k)+AF(k)+BT(k)/eSQ(k)**2*e(k)*Su(k))
**-To(Sigma)
        C(IE(15+n))= WW
**-OM
        C(IE(16+n))=-(BT(k)/eSQ(k)*Cue-AF(k)*eSQ(k)*Su(k))
**-Gm
        C(IE(17+n))=-Su(k)
**-r
        C(IE(18+n))= 2.d0*dlog(Znn(k))
**-s
        C(IE(19+n))=-(AF(k)*Cue+BT(k)*Su(k))*Gc3(k)*2.d0/Znn(k)/X(k)

        if (mBT) then
**-OMdot in BT model
        C(IE(20+n))= C(IE(16+n))*DTPER(k)
        else
**-OMdot
        C(IE(20+n))=(Ane(k)-EXC(k)+e(k)*Su(k))*C(IE(16+n))*Pb(k)
        endif!(mDADE,mDADEP)

        if(mDADEP.or.mBTPP) then
**-Xdot in DD+,BT++ model
        C(IE(21+n))=(Ane(k)-EXC(k)+e(k)*Su(k))*C(IE(12+n))*Pb(k)
**-Xd2  in DD+,BT++ model
        C(IE(102+n))=C(IE(21+n)) *Pb(k)*Ane(k)/2.0
        else
**-Xdot
        C(IE(21+n))= C(IE(12+n))*DTPER(k)
**-Xd2
        C(IE(102+n))=C(IE(21+n))*DTPER(k)/2.0d0
        endif
**-Pbdot
        C(IE(22+n))= C(IE(13+n))*DTPER(k)/2.0d0
**-Pbdot-dot
        C(IE(101+n))= C(IE(22+n))*DTPER(k)/2.0 !<<
**-Pb-3dot
        C(IE(104+n))= C(IE(101+n))*DTPER(k)/2.0 !<<
**-edot
        C(IE(23+n))= C(IE(14+n))*DTPER(k)
**-dtheta
        C(IE(24+n))= BT(k)/eSQ(k)*e(k)**2*Su(k)

        enddo ! (k=1,norb)

        endif !B

	IF(IDEGEQ.LT.IDEGS) C(IDEGEQ+1)=0.0D0

* (Build Norm.Eq.Syst.-SUM(ij))

	DO 10 I=1,IDEGEQ
	DO 20 J=1,IDEGEQ
 20	SUM(I,J) =C(I)*C(J)*Weg + SUM(I,J)
 10	RHS(I)   =DTII*C(I)*Weg + RHS(I)
        D0 =(DTII)**2*Weg + D0
	RETURN
	END

*========================================================================
        integer*4 function ISM(IE,N)

* Summing of N elements of IE(N) vector

        dimension IE(1)
	integer*4 N
	ISM=0
	do 1 i=1,N
 1	ISM=IE(i)+ISM
	return
	end

*========================================================================
	SUBROUTINE TDTUTC(JDT,TDT,UTC,filleap)
        implicit none

*  Transrormation from UTC to TDT(TAI).
* modified by C.Lange-1998
      
      real*8 JDT,TDT,UTC
      integer*4 DTA,i,ileap(99),nleap
      integer*4 JDJ
      character filleap*80
      
      JDJ=JDT+UTC/86400.0D0 - 2400000.5D0
      open(17,file=filleap,status='OLD',err=200)
      do i=1,99
         read(17,*,END=100) ileap(i)
         if(ileap(i).gt.1000000) ileap(i)=ileap(i)-2400000
      enddo
 100  nleap=i-1
      close(17)
      do i=1,nleap
	IF(JDJ.GE.ileap(i)) DTA=10+i
      enddo

      TDT=UTC+32.184D0+DFLOAT(DTA)
	
      RETURN
 200  write(*,'(a)') ' >> leap second-file not found > '//filleap
      stop
      END
        
*========================================================================
        
        subroutine UTCCOR(codobs,MJD,cor,clkdate,clockc,clsite,maxclk)
*
* find clock correction to get UTC time scale
* The correctios read from file time.dat
* codobs - letter observatory code
* clkdate(MaxCLK) - MJDs of clock correction
* clksite(MaxCLK) - one-letter observatory code
* clockc(MaxCLK) - clock correction, us
* cor - output correction to UTC in us
        implicit none
        real*8 MJD,cor,cjd,cjd1,c2,c4,clkdate(1),clockc(1)
        character cod*1,codobs*1,clsite(1)*1
        integer*4 k,maxclk
        
        k=0
        cor =0.0
        cjd1=0.0
        c4  =0.0
        
 10     cjd=cjd1
        c2=c4
 50     continue
        k=k+1
        cjd1= clkdate(k)
        c4  = clockc(k)
        cod = clsite(k)

c        if (cod.eq.'a') write(*,*) cjd1,c4
        
        if (k.ge.maxclk) then
          write(*,'(f12.5,1x,a1,a)')
     +     MJD,codobs,' >> No UTC correction! d(UTC)=0'
           cor=0.0
           return
        endif
        
        if ((MJD.gt.cjd1).or.(cod.ne.codobs)) goto 10
        
*(assume linaer trend of clock)
        
        cor=(MJD-cjd)*(c4-c2)/(cjd1-cjd)+c2
        cor=-cor
        return
        end
      
*========================================================================
        
        subroutine READCLK(io,filutc,clkdate,clsite,clockc)
        
* read observatory clock correction file
* filutc - clock corrections (time.dat) file
* clkdate(MaxCLK) - MJDs of clock correction
* clksite(MaxCLK) - one-letter observatory code
* clockc(MaxCLK) - clock correction, us
* O.D.
        
        implicit none
        character filutc*80,fil1*80,dirclk*60,str50*50,clsite(1)*1,cod*1
        integer*4 clklen,i,filre,k,io
        real*8 clkdate(1),clockc(1),cjd1,c3,c4
        
        filre=124
        open(124,file=filutc,status='OLD')
        k=0
* find clock directory pointer in file name 'filutc'
* dirclk - directory of clock correction file
        
        clklen=1                    ! length of string before '/'
        do i=1,len(filutc)
        if((filutc(i:i).eq.'/').or.(filutc(i:i).eq.'\')) clklen=i
        enddo
        dirclk=filutc(1:clklen)     ! clock directory

* read data from clock correction file
        
 50     read(filre,'(a50)',end=30) str50
        
* if command 'INCLUDE' then open this file
        
        if(((str50(6:6).ne.'.').and.(str50(7:7).ne.'.')
     +   .and.(str50(8:8).ne.'.'))
     +   .and.(str50(1:7).ne.'INCLUDE')) goto 50
        
        if(str50(1:7).eq.'INCLUDE') then
        fil1=dirclk(1:clklen)//str50(9:50)
        filre=125
        close(filre)
        open(unit=filre,file=fil1,status='OLD',Err=60)
        goto 50
 60     continue
        write(*,'(a)') ' >> INCLUDE file not found:'//fil1
        filre=124                       ! go back to 'time.dat' file
        goto 50
        endif ! (='INCLUDE')
        
        if ((str50(6:6).eq.'.').or.(str50(7:7).eq.'.')
     +  .or.(str50(8:8).eq.'.')) then
        read(str50,'(f10.3,1x,f10.3,1x,f11.3,1x,a1)',end=50)
     +  cjd1,c3,c4,cod
        
        k = k+1
        if(c3.gt.800.0d0) c3=c3-818.8d0
        clkdate(k)=cjd1
        clockc(k)=c3-c4
        clsite(k)=cod
        goto 50
        else
        goto 50
        endif
* go back if end was in INCLUDE file
 30     if (filre.eq.125) then
        filre=124
        goto 50
        endif
        close(124)
        close(125)

        write(IO,'(i10,a)') k,' Clock corrections read'
        
        return
        end      
      
*========================================================================
      
	include 'timlib.f'
	include 'form1.f'
