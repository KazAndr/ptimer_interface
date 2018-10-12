c ver 1.01
*======================================================================
	SUBROUTINE CALDJD(JD,NUM,MON,YEA,WEEK )
* 
* using JD (INPUT),
* returns day/month/yaer/day of week
*
* O.D.

	implicit none

	INTEGER*4 NUM,MON,M
	CHARACTER WEEK*5
	REAL*8 JD,YEA,JD1,Z,F,A,AL,B,C,D,E,CC,REM,R,RMON

	JD1=JD+0.5D0
	Z = DINT(JD1+0.0001D0)
	F = JD1-Z
	A = Z
	IF (Z.LT.2299161.0D0) GOTO 10
	AL=DINT((Z-1867216.25D0)/36524.25D0)
	A=1.0D0+AL-DINT(AL/4.0D0)+Z
 10	B=A+1524.0D0
	C=DINT((B-122.1D0)/365.25D0)
	D=DINT(365.25D0*C)
	E=DINT((B-D)/30.6001D0)
	CC=B-D-DINT(30.6001D0*E)+F
	NUM=IDINT(CC)
	IF(E.LT.13.5D0) RMON=E-1.0D0
	IF(E.GT.13.5D0) RMON=E-13.0D0
	YEA=C-4716.0D0
	IF(RMON.LT.2.5D0) YEA=YEA+1.0D0
	MON=IDINT(RMON)
	JD1=JD1+1.0D0
	REM=JD1/7.0D0
	R=REM-DINT(REM)
	R=R*7.0D0
	M=IDINT(R)
	IF(M.EQ.0) WEEK='Sun.'
	IF(M.EQ.1) WEEK='Mon.'
	IF(M.EQ.2) WEEK='Tue.'
	IF(M.EQ.3) WEEK='Wed.'
	IF(M.EQ.4) WEEK='Thu.'
	IF(M.EQ.5) WEEK='Fri.'
	IF(M.EQ.6) WEEK='Sat.'
	RETURN
	END

*======================================================================
*
	SUBROUTINE JDCALD( DATA,JD )

* for given year/month/day DATA returns Julian date JD
*
* DATA :YYYY.MMDDDD (INPUT) 
* JD                (OUTPUT)
* O.D.

	implicit none
	REAL*8 DATA,JD,Y,RRM,RM,RD,A,B,C

	Y=DINT( DATA )
	RRM=(DATA-Y)*1.0D2
	RM=DINT(RRM)
	RD=((RRM-RM)*1.0D2)
	IF((RM-2.0D0)) 3,3,2
 3	Y=Y-1.0D0
	RM=RM+12.0D0
 2	A=DINT(Y/1.0D02)
	B=2.0D0-A+DINT(A/4.0D0)
	IF(DATA.LT.1582.1015D0) B=0.0D0
	A=DINT(30.6001D0*(RM+1.0D0))
	C=DINT(365.25D0*Y)
	JD=RD+C+A+B+1720994.5D0
	RETURN
	END

*======================================================================

	SUBROUTINE SPHERC ( X,Y,Z, AE,DE,RE )

* spherical coordinates (AE,DE,RE) (OUTPUT)
* from rectangular (X,Y,Z)         (INPUT)


	implicit none
*
	REAL*8 X,Y,Z,AE,DE,RE,PI2,PI
	PI=3.1415926535897932D0
	PI2=2.0D0*PI
	RE=DSQRT(X**2+Y**2+Z**2)

**	 Spherical coordinates
	AE=DATAN(Y/X)
	DE=DSQRT(X*X+Y*Y)
	DE=DATAN(Z/DE)
	IF(X.LT.0.0D0) AE=AE+PI
	IF((X.GT.0.0D0).AND.(Y.LT.0.0D0)) AE=AE+PI2
	IF((X.LT.0.0D0).AND.(Y.GT.0.0D0)) AE=AE+PI2
	IF(AE.GE.PI2) AE=AE-PI2
	RETURN
	END

*======================================================================

        SUBROUTINE RADIAN(RDEG,IND,RAD)

* converts angle RDEG
* 
*format of Input/Output data:
*
* ---INPUT-----  -OUTPUT-
* (IND)  (RDEG)  (RAD) 
*  0     DEG.'"   RAD 
*  1     H.mmss   RAD 
*  2     DEG.'"   DEG.XXX
*        h.mmss   HOUR.XXX 

	implicit none
	REAL*8 RDEG,RAD,PI,SEC
	INTEGER*4 IDEG,IMIN,IND

	PI=3.141592653589793D0

	IDEG = DINT(RDEG)
	IMIN = DINT((RDEG-IDEG)*(100.0+1d-12))
	SEC  = RDEG*1.0d4-IDEG*1.0d4-IMIN*1.0d2

        RAD=IDEG+IMIN/60.0d0+SEC/3600.0d0

	IF(IND.EQ.2) RETURN
	RAD=RAD*PI/12.0D0
	IF(IND.EQ.0) RAD=RAD/15.0D0
	RETURN
	END
 
*======================================================================

        SUBROUTINE RADEG(RAD,I,DEG,W)
*
* converts angle RAD and writes it in conventional manner to string W
*
* format of Input/Output data:
*
* ---INPUT-----  -OUTPUT------
* (I)    (RAD)   (DEG)      (W)
*  0      radians oo.''""   string oo '' "".""
*  1      =       hh.mmssss string hh mm ss.ds
*
* O.D.

	implicit none
        REAL*8 RAD,DEG,B,DG,DM,DS,DRM,DRS
        REAL*4 S
        integer*4 i
        character W*17,SS*1
        S=SIGN(1.0,RAD)
        B=DABS(RAD)*180.0D0/3.141592653589793d0
	IF(I.NE.0)B=B/15.0D0
	DG=DINT(B)
	DRM=(B-DG)*60.0D0
	DM=DINT(DRM)
	DRS=(DRM-DM)*60.0D0
	DS=DRS
        DEG=(DG+1.0D-02*DM+1.0D-04*DS)*S
        SS=' '
        if (S.LE.0.) SS='-'
        write(W,'(A1,i2,1H:,i2,1H:,f10.7)') SS,NINT(DG),NINT(DM),DS
	do i=2,9
	if (W(i:i).EQ.' ') W(i:i)='0'
	enddo
	RETURN
	END

*======================================================================

		SUBROUTINE READEL(JED,TSEC,IERR,fileph)

*
* Reading and interpolation of the JPL ephemerides DE200/LE200
* which are given as coefficents of Chebyshev's polinomas.
*
*
* For timing purposes the binary data file de200.bin is used
* with direct access of data records. Check the OPEN statement
* where you should specify correct record length (RECL=..) -
* it is compiler dependent. The record length should be equivalent
* to the length of (826 x real*8) numbers.
*
* The center of the reference frame BSS is choosen by setting ICENT=0
* Other values of ICENT may choosen according to values of 
* index (i), i=1..11 (see below)
*
* The choice of planets which positions/velocities to be calculated
* is made by setting of elements of array IREQ(i), i=1..13:
*
* IREQ(i)=0 - positions/velocities will not be calculated
* IREQ(i)=1 - positions will be calculated
* IREQ(i)=2 - velocities will be calculated
* 
* Index (i) refers to:
*
* 1-Mercury; 2-Venus; 3-Earth; 4-Mars; 5-Jupiter; 6-Saturn;
* 7-Uran; 8-Neptune; 9-Pluto; 10-Sun; 11-Moon; 12-Earth/Moon Barycenter;
* 13-Nutations
*
* The values of velocities/positions returned in array TABOUT(k,i)
* k=1,2,3 - X,Y,Z; k=4,5,6 - Xdot,Ydot,Zdot [AU and AU/DAY].
* i=1..12 - see index (i) above.
* The nutation coefficents will be returned in DNUT(i)
*
* (INPUT) JED, TSEC - time (TB) for which positions/velocities to be 
*                     calculated (JED - at midnight, TSEC - in seconds)
*
* The program transfer data through common blocks CETBL1..CETBL6
*
*
*

	implicit none
	INTEGER*4 LIST(11),LOC(3,12)
	REAL*8 JED,TSEC,DJE,BUF,PV,PVSUN,SS,T,EM,EMR,AUFAC,
     &	BF,TABOUT,NUT,DEJ
	REAL*4 EMRAT
	INTEGER*4 ICW,ICENT,IREQ,NZP,I,LST,IC,NR,NSK,LL,NM,J,ierr
	DIMENSION PV(6,13),PVSUN(6),SS(6),T(2)
	character fileph*80
	COMMON/CETBL1/EMRAT
	COMMON/CETBL2/ICW,ICENT,IREQ(13)
	COMMON/CETBL3/BUF(826)
	COMMON/CETBL4/TABOUT(6,12),NUT(4)
	COMMON/NZP/NZP
	
	DATA LOC/3,12,4,147,12,1,183,15,2,273,10,1,303,9,1,
     &	330,8,1,354,8,1,378,6,1,396,6,1,414,12,8,702,15,1,
     &	747,10,4/,SS/2439632.5D0,2458800.5D0,32.0D0,32.0D0,
     &	12.0D0,4.0D0/,BF/149597870.66D0/

* Note that above SS(1) and SS(2) are the first and the last
* Julian date of ephemerides file DE200.BIN

     	IERR=0
	EM=EMRAT
	IF(EM.EQ.0.0D0)EM=81.3007D0
	EMR=1.0D0/(1.0D0+EM)
	DJE=TSEC/86400.0D0
	IF(ICW.EQ.2) GOTO 1
	ICW=2
	BUF(1)=0.0D0
	BUF(2)=0.0D0
	NZP=0
 1	CONTINUE
	IF(ICENT.GE.0.AND.ICENT.LE.11) GOTO 2
	IERR=4
	RETURN
 2	CONTINUE
	DO 3 I=1,9
	IF(IREQ(I).GE.0.AND.IREQ(I).LE.2) GOTO 3
	IERR=3
	RETURN
 3	LIST(I)=IREQ(I)
	LIST(10)=IREQ(11)
	IF(ICENT.EQ.0.OR.ICENT.EQ.10) GOTO 5
	LST=0
	DO 4 I=1,10
 4	LST=MAX0(LST,LIST(I))
	IC=ICENT-ICENT/11
	LIST(IC)=MAX0(LST,IREQ(10))
 5	CONTINUE
	LIST(3)=MAX0(LIST(3),LIST(10))
	LIST(3)=MAX0(LIST(3),IREQ(12))
	LIST(10)=LIST(3)
	LIST(11)=IREQ(13)
	T(2)=SS(3)
	AUFAC=1.0D0/BF
	DEJ=JED+DJE
	IF(DEJ.LT.SS(1)) GOTO 99
	IF(DEJ.GT.SS(2)) GOTO 98
	NR=IDINT((DEJ-SS(1))/SS(4))+1


	NSK=NR-NZP
	IF(NSK) 18,19,18
  18	continue
	open(unit=12,file=fileph,access='DIRECT',
     &	status='OLD',form='UNFORMATTED',recl=6608,err=45)
*     &	blocksize=32768)
	goto 46
 45	write(*,*) ' ERROR: Not found',fileph
	write(*,*) ' Check config file ptimer.cfg'
	stop
	
 46	READ(12,REC=NR) BUF
	close(unit=12)
 19	CONTINUE

	NZP=NR
 6	IF(DEJ.GE.BUF(1).AND.DEJ.LE.BUF(2)) GOTO 7
	IERR=5
	RETURN
 7	T(1)=(DEJ-BUF(1))/T(2)
	LL=LOC(1,11)
	CALL XINTRP(BUF(LL),T,LOC(2,11),3,LOC(3,11),2,PVSUN)
	DO 9 I=1,10
	IF(LIST(I).EQ.0) GOTO 9
	LL=LOC(1,I)
	CALL XINTRP(BUF(LL),T,LOC(2,I),3,LOC(3,I),LIST(I),
     &	PV(1,I))
	NM=LIST(I)*3
	DO 8 J=1,NM
	IF(I.LE.9.AND.ICENT.NE.0)PV(J,I)=(PV(J,I)-PVSUN(J))
     &	*AUFAC
	IF(I.LE.9.AND.ICENT.EQ.0)PV(J,I)=PV(J,I)*AUFAC
	IF(I.EQ.10)PV(J,I)=PV(J,I)*AUFAC
 8	CONTINUE
 9	CONTINUE
	DO 10 J=1,6
 10	PVSUN(J)=PVSUN(J)*AUFAC
	IF(LIST(11).EQ.0) GOTO 11
	LL=LOC(1,12)
	CALL XINTRP(BUF(LL),T,LOC(2,12),2,LOC(3,12),
     &	LIST(11),NUT)
	GOTO 11
 99	CONTINUE
	IERR=1
	RETURN
 98	IERR=2
	RETURN
 11	IF(IERR.NE.0) RETURN
	DO 12 I=1,6
	PV(I,11)=PV(I,3)+EM*EMR*PV(I,10)
	PV(I,12)=PV(I,3)
	PV(I,3)=PV(I,3)-EMR*PV(I,10)
	PV(I,10)=0.0D0
 12	CONTINUE
	IF(ICENT.EQ.0) GOTO 14
	DO 13 I=1,12
	DO 13 J=1,3
	IF(IREQ(I).EQ.0) GOTO 13
	TABOUT(J,I)=PV(J,I)-PV(J,ICENT)
	IF(IREQ(I).EQ.1) GOTO 13
	TABOUT(J+3,I)=PV(J+3,I)-PV(J+3,ICENT)
 13	CONTINUE
	RETURN
 14	DO 15 I=1,12
	DO 15 J=1,3
	IF(IREQ(I).EQ.0) GOTO 15
	TABOUT(J,I)=PV(J,I)
	IF(I.EQ.10)TABOUT(J,10)=PVSUN(J)
	IF(IREQ(I).EQ.1) GOTO 15
	TABOUT(J+3,I)=PV(J+3,I)
	IF(I.EQ.10)TABOUT(J+3,10)=PVSUN(J+3)
 15	CONTINUE
	RETURN
	END
	
*======================================================================

	SUBROUTINE XINTRP(BUF,T,NCF,NCM,NA,FL,PVA)

c  This subroutine differentiates and interpolates a		
c  set of Chebyshev coefficients to give Position,		
c  Velocity,and Acceleration.					
c  Calling sequence parameters--				
c 								
c Input--BUF 	 ist location of array of DP chebyshev coef	
c ---------- 	 ficients of position.				
c 	     t	 T(1) is DP fractional time in interval co-	
c 		 vered by coefficients at which  interpola-	
c 		 tion is wanted(0.LE.T(1).LE.1). T(2) is DP	
c 		 length of whole.				
c 								
c 	   NCF	 no.of coefficients per component		
c 								
c 	   NCM	 no.of components per set of coefficients	
c 								
c 	    NA	 no.of sets of coefficients in full array	
c 		 (i.e.,no. of sub-intervals in full inter-	
c 		 val).						
c 								
c 	    FL	 INTEGER FLAG =1  for Positions only		
c 			      =2  for P,V  only			
c 			      =3  for P,V,A			
c 								
c Output--PVA  interpolated quantities requested.		
c 		 dimension expected is PVA(NCM,FL) , DP.	

	implicit none
	integer*4 NCF,NCM,NA
	REAL*8 BUF(NCF,NCM,1),T(2),PVA(NCM,1),PC(18),VC(18),
     &	AC(18),BMA,BMA2,TEMP,TC,TCL,TWOT,DINT,ARG
	EQUIVALENCE(TCL,PC(2))
	INTEGER*4 FL,L,NP,NV,NAC,M,I,J,JJ
	DATA PC(1)/1.0D0/,PC(2)/2.0D0/,VC(2)/1.0D0/,AC(3)/4.0D0/
	DINT(ARG)=ARG-DMOD(ARG,1.0D0)
*
*	   Calculate location of correct component arrays
*	   and scaled Chebyshev time in that interval
*
	TEMP=T(1)*DFLOAT(NA)
	L=IDINT(TEMP-DINT(T(1)))+1
	TC=2.0D0*((TEMP-DINT(TEMP))+DINT(T(1)))-1.0D0
*
*	   Check to see whether Chebyshev time has changed,
*	   and compute new polynomial values if it has.
*
	IF(TC.EQ.TCL) GOTO 1
	NP=2
	NV=3
	NAC=4
	TCL=TC
	TWOT=TC+TC
 1	CONTINUE
	IF(NP.GE.NCF) GOTO 2
	M=NP+1
	NP=NCF
	DO 3 I=M,NCF
 3	PC(I)=TWOT*PC(I-1)-PC(I-2)
 2	CONTINUE


*	   Interpolate to get position for each component


	DO 4 I=1,NCM
	PVA(I,1)=0.0D0
	DO 4 J=1,NCF
	JJ=NCF-J+1
	PVA(I,1)=PVA(I,1)+PC(JJ)*BUF(JJ,I,L)
 4	CONTINUE
	IF(FL.LE.1)RETURN

*	   Check velocity polynomial values,and generate
*	   new ones if required

	BMA=2.0D0*DFLOAT(NA)/T(2)
	VC(3)=TWOT+TWOT
	IF(NV.GE.NCF) GOTO 5
	M=NV+1
	NV=NCF
	DO 6 I=M,NCF
 6	VC(I)=TWOT*VC(I-1)+PC(I-1)+PC(I-1)-VC(I-2)

*	   Evaluate velocity for each component

 5	CONTINUE
	DO 7 I=1,NCM
	PVA(I,2)=0.0D0
	DO 11 J=2,NCF
	JJ=NCF-J+2
	PVA(I,2)=PVA(I,2)+VC(JJ)*BUF(JJ,I,L)
 11	CONTINUE
	PVA(I,2)=PVA(I,2)*BMA
 7	CONTINUE
	IF(FL.EQ.2) RETURN

*	   Check acceleration polynomial values , and
*	   re-do if necessary

	BMA2=BMA*BMA
	AC(4)=24.0D0*PC(2)
	IF(NAC.GE.NCF) GOTO 8
	M=NAC+1
	NAC=NCF
	DO 9 I=M,NCF
 9	AC(I)=TWOT*AC(I-1)+4.0D0*VC(I-1)-AC(I-2)
 8	CONTINUE

*	   Get acceleration for each component

	DO 10 I=1,NCM
	PVA(I,3)=0.0D0
	DO 12 J=3,NCF
	JJ=NCF-J+3
	PVA(I,3)=PVA(I,3)+AC(JJ)*BUF(JJ,I,L)
 12	CONTINUE
	PVA(I,3)=PVA(I,3)*BMA2
 10	CONTINUE
	RETURN
	END
	
*======================================================================
	SUBROUTINE TTTB(JD,TDT,TDB)

* Relativistic time scales transformation.
* Transformation from TDT to TDB is made to take into account
* the Einstein time delay and the second-order Doppler shift which
* observatory clock suffer in gravitational potential of
* the Solar system's bodies.
*
* for TDT->TDB transformation a semi-analytical formula of
* (Fairhead L.& Bretagnon P. (Astron.Aph.1990.v.229.pp.240-247)
* is used. The coefficents of the transformation merged in program
* were kindly provided by L.Fairhead, BDL, France.
*
* (INPUT)  UTC    - UTC time, seconds
*          JD     - Julian date corresponding to TDT value
*          TDT    - value of TAI to be transformed to TDB scale,seconds
* (OUTPUT) TDB    - TDB, seconds
*
* O.D.

	implicit none

	REAL*8 JD,TDT,TDB,PI,PI2
      REAL*8 AA(127),WW(127),FF(127),RAD,AP,T,A,W,F

        integer*4 IN,I
	DATA AA /
     & 1.43388D-01,1.3374D-02,1.6916D-02,1.9476D-02,1.22605D-01,
     & 3.70115D-01 ,4.32299D0,1.0084D-02,9.3570D-03,1.0143D-02,
     & 1.1081D-02 ,1.2152D-02,1.2869D-02,1.1774D-02,1.3276D-02,
     & 1.3794D-02 ,1.7925D-02,2.1774D-02,2.4976D-02,2.2997D-02,
     & 2.519D-02 , 2.7764D-02,2.9198D-02,3.3595D-02,3.2088D-02,
     & 3.4420D-02 ,5.9146D-02,5.4764D-02,7.7996D-02,2.1056D-01,
     & 2.6591D-01,2.69668D-01,1.706807D0,102.156724D0,9.37D-03,
     & 1.0099D-02, 9.858D-03, 9.963D-03,1.0825D-02,1.1622D-02,
     & 1.1847D-02,1.2420D-02,1.0453D-02,1.1707D-02,1.0396D-02,
     & 1.5078D-02,1.0962D-02,1.2462D-02,1.4318D-02,1.1942D-02,
     & 1.3671D-02,1.7806D-02,1.4223D-02,1.5974D-02,1.5949D-02,
     & 1.4751D-02,1.6155D-02,2.2567D-02,1.7673D-02,2.1691D-02,
     & 2.5196D-02,2.4816D-02,2.0322D-02,2.0937D-02,2.2509D-02,
     & 2.9862D-02,2.7567D-02,3.2423D-02,3.3529D-02,2.9247D-02,
     & 3.0215D-02,3.2438D-02,3.4867D-02,3.3477D-02,2.8244D-02,
     & 4.0480D-02,4.2560D-02,4.2732D-02,3.6955D-02,3.6507D-02,
     & 4.0759D-02,3.6564D-02,4.0184D-02,4.2411D-02,5.4139D-02,
     & 4.6551D-02,5.8844D-02,4.8373D-02,4.8042D-02,6.3814D-02,
     & 6.4397D-02,7.5019D-02,6.2617D-02,7.9645D-02,8.0164D-02,
     & 1.01868D-01,9.83580D-02,1.37927D-01,1.16120D-01,1.18971D-01,
     & 1.19979D-01, 1.5908D-01,1.43935D-01,2.03747D-01,2.30685D-01,
     & 1.73435D-01, 2.43085D-01,3.7551D-01,4.68597D-01,4.32392D-01,
     & 4.86306D-01,4.96817D-01,6.00309D-01,4.35206D-01,4.47061D-01,
     &  7.94185D-01,1.115322D0,1.193379D0,1.276839D0,1.554905D0,
     & 1.694205D0,2.256707D0,4.67674D0,4.770086D0,13.839792D0,
     & 22.417471D0,1656.674564D0/
	DATA WW /
     & 6283.075849991D0,-3.523118349D0,529.690965095D0,
     & 213.299095438D0,12566.151699983D0,0.0D0,6283.075849991D0,
     & 522.577418094D0,5486.777843175D0,4694.002954708D0,
     &-7.113547001D0,1059.381930189D0,6076.890301554D0,
     & 12036.460734888D0,6062.663207553D0,426.598190876D0,
     &-775.522611324D0,206.185548437D0,5760.498431898D0,
     &-796.298006816D0,5746.271337896D0,155.420399434D0,
     & 5856.477659115D0,5507.553238667D0,18849.227549974D0,
     &-398.149003408D0,26.2983198D0,1577.343542448D0,5223.693919802D0,
     &-3.523118349D0,529.690965095D0,213.299095438D0,12566.151699983D0,
     & 6283.075849991D0,149854.400134205D0,283.859318865D0,
     & 6309.374169791D0,149.563197135D0,553.569402842D0,
     & 5120.601145584D0,5643.178563677D0,4690.479836359D0,
     & 5863.591206116D0,-4705.732307544D0,951.718406251D0,
     & 19651.048481098D0,3.590428652D0,103.092774219D0,
     & 16730.463689596D0,8031.092263058D0,-536.804512095D0,
     & 73.297125859D0,17789.845619785D0,-2352.866153772D0,
     &-220.412642439D0,1349.867409659D0,10213.285546211D0,
     & 6133.512652857D0,6812.766815086D0,14143.495242431D0,
     & 1748.016413067D0,-1194.447010225D0,419.484643875D0,
     & 8429.241266467D0,10447.387839604D0,12139.553509107D0,
     & 6279.552731642D0,8827.390269875D0,9437.762934887D0,
     &-71430.695617928D0,7084.896781115D0,6076.890301554D0,
     & 522.577418094D0,6062.663207553D0,-6286.59896834D0,
     & 15720.838784878D0,161000.685737473D0,632.783739313D0,
     & 3154.687084896D0,801.820931124D0,12352.852604545D0,
     & 5088.628839767D0,-7.113547001D0,6275.962302991D0,
     & 17260.154654690D0,-9.803210680D-01,426.598190876D0,
     & 155.420399434D0,2146.165416475D0,5760.498431898D0,
     & 5746.271337896D0,2942.463423292D0,20.775395492D0,
     & 4694.002954708D0,206.185548437D0,-5573.142801634D0,
     & 2544.314419883D0,11790.629088659D0,1059.381930189D0,
     & 5486.777843175D0,38.133035638D0,10977.078804699D0,
     &-796.298006816D0,12036.460734888D0,5856.477659115D0,
     & 18849.227549974D0,-775.522611324D0,5507.553238667D0,
     & 6244.942814354D0,74.781598567D0,5884.926846583D0,
     & 6208.294251424D0,1577.343542448D0,-398.149003408D0,
     & 26.298319800D0,11506.769769794D0,3930.209696220D0,
     & 5223.693919802D0,7860.419392439D0,77713.771467920D0,
     &-3.523118349D0,213.299095438D0,6069.776754553D0,
     & 529.690965095D0,12566.151886066D0,5753.384970095D0,
     & 6283.075943033D0/
	DATA FF /
     & 1.131453581D0,1.502210314D0,4.510959344D0,1.642186981D0
     & ,2.438140634D0,4.712388980D0,2.642893748D0,7.493202620D-01
     & ,3.416081409D0,4.044013795D0,5.154724984D0,6.222874454D0
     & ,5.333425680D0,2.292832062D0,5.845801920D0,2.699831988D0
     & ,1.092065955D0,3.854787540D0,2.467913690D0,1.174411803D0
     & ,2.980330535D0,3.745318113D0,6.238118630D-01,5.980162321D0
     & ,4.162913471D0,5.980077351D0,1.083044735D0,4.53480017D0
     & ,4.670344204D0,6.262738348D0,5.836047367D0,3.400290479D0
     & ,4.205904248D0,4.249032005D0,6.738803950D-01,1.942176992D0
     & ,1.06181641D0,4.870690598D0,8.427150110D-01,4.863931876D0
     & ,5.489005403D0,4.734090399D0,1.91370455D0,2.654125618D0
     & ,5.717799605D0,3.96948077D0,2.196567739D0,1.737438797D0
     & ,3.016058075D0,2.053414715D0,5.971672571D0,3.475975097D0
     & ,2.104551349D0,6.145309371D0,4.00529827D0,4.308933301D0
     & ,1.331103168D0,3.307984806D0,3.186129845D0,5.952658009D0
     & ,2.901883301D0,1.087136918D0,3.735430632D0,6.523034140D-01
     & ,1.460726241D0,1.770181024D0,5.040846034D0,5.541473556D0
     & ,2.404714239D0,4.183178762D0,3.389610345D0,7.493174120D-01
     & ,5.210064075D0,4.144987272D0,5.069663519D0,2.546610123D0
     & ,1.270837679D0,5.720622217D0,5.071801441D0,6.248866009D0
     & ,3.981496998D0,3.324679049D0,3.565975565D0,2.869567043D0
     & ,3.411091093D0,9.215735391D-01,4.839650148D0,2.251573730D0
     & ,1.495846011D0,4.167901731D0,1.280308748D0,4.980931759D0
     & ,2.654394814D0,2.949233637D0,2.095377709D0,5.984503847D0
     & ,9.279388601D-02,1.135934669D0,8.735041230D-01,1.914547226D0
     & ,4.551585768D0,1.890075226D0,5.957517795D0,4.333987818D0
     & ,4.773852582D0,6.153743485D0,3.651837925D0,4.103476804D0
     & ,5.866398759D0,2.435898309D0,5.200071790D-01,5.696701824D0
     & ,2.678271909D0,4.349338347D0,3.615796498D0,2.322313077D0
     & ,1.422745069D0,3.649823730D0,5.988822341D0,5.198467090D0
     & ,5.025132748D0,5.543113262D0,4.021195093D0,4.444016030D-01
     & ,6.196904410D0,4.296977442D0,6.240054195D0/
	PI=3.1415926535897932D0
	PI2=2.0D0*PI
	RAD=206264.806D0
	AP=0.0D0

* calculate TDT->TDB

	T=(JD+TDT/86400.0d0-2451545.0d0)/(1000.0d0*365.25d0)
	IN=0
** Terms - T**3
 	DO 10 I=1,1
 	IN=IN+1
 	A=AA(IN)
 	W=WW(IN)
 	F=FF(IN)
 10	AP=T*T*T*A*DSIN(F+W*T)+AP
** - T**2
	DO 20 I=1,6
	IN=IN+1
	A=AA(IN)
	W=WW(IN)
	F=FF(IN)
 20	AP=T*T*A*DSIN(F+W*T)+AP
** - T**1
	DO 30 I=1,27
	IN=IN+1
	A=AA(IN)
	W=WW(IN)
	F=FF(IN)
 30	AP=T*A*DSIN(F+W*T)+AP
** - T**0
 	DO 40 I=1,93
 	IN=IN+1
 	A=AA(IN)
 	W=WW(IN)
	F=FF(IN)
 40	AP=A*DSIN(F+W*T)+AP
	
	TDB = TDT+AP*1.0D-06

	RETURN
	END

*======================================================================

	subroutine out(NOPTS,JDD,ZA,RM,Y,WEG,PNM,IO)
	implicit none

* Prints table of TOAs, residuals
* into logical device output IO
* O.D.

	REAL Y(1)
        REAL R22,R33,R15,R14
	REAL*8 JDD(1),RM(1),WEG(1),YE,R12,R13,ZA(1),JD,R,X,JDC,TT
        integer*4 NN,IM,IND,I11,NOPTS,IO,NW0,MON,II,NOM,KYA,NW,I12,I13
        CHARACTER WW*5,PNM*7,imnt*3,NWC*2

	IND=1
	open(12,file='_tim.res')

	NW0=0
	MON=0

	WRITE(IO,76)  PNM
	WRITE(IO,79)

	DO 100 II=1,NOPTS

	JD=JDD(II)
	R=RM(II)

	X=ZA(II)
	IF(R.GE.0.0D0) GOTO 102
	JD=JD-1.0D0
	R=R+86400.0D0
	X=X+86400.0D0
 102	CONTINUE
	JDC=JD+R/86400.0D0
	CALL CALDJD(JDC, NN , IM , YE , WW)
	NOM=IDINT(YE)
	KYA=1900
	NW =NOM-KYA
	if(NW.GE.100) NW=NW-100 ! y2k-problem solved here
	write(nwc,'(I2)') NW
	if(NWC(1:1).eq.' ') NWC(1:1)='0'

	IF(NW0.NE.NW) GOTO 444
	GOTO 445
 444	WRITE(IO,75)
 445    IF(IM.EQ.1) IMNT='Jan'
        IF(IM.EQ.2) IMNT='Feb'
        IF(IM.EQ.3) IMNT='Mar'
        IF(IM.EQ.4) IMNT='Apr'
        IF(IM.EQ.5) IMNT='Mai'
        IF(IM.EQ.6) IMNT='Jun'
        IF(IM.EQ.7) IMNT='Jul'
        IF(IM.EQ.8) IMNT='Aug'
        IF(IM.EQ.9) IMNT='Sep'
        IF(IM.EQ.10)IMNT='Oct'
        IF(IM.EQ.11)IMNT='Nov'
        IF(IM.EQ.12)IMNT='Dec'
        I11=DINT(JD-2400000.5D0+0.000000001d0)
	R12=R
	I12=IDINT(R12)
	R22=R12-DINT(R12)
	R13=X
	I13=IDINT(R13)
	R33=ABS(R13-DINT(R13))
        R14=Y(II)
        R15=WEG(II)
        WRITE(IO,77) II,NN,IMNT,NWC,I11,I12,R22,R15,I13,R33,R14
 77     FORMAT(1X,I4,')',1X,I2,'-',a3,'-',A2,I7,I7,F9.8,
     &  1X,f8.2,I7,F9.8,F11.2)
	MON=IM
	NW0=NW
		IF(IND.NE.1) GOTO 100
		TT=(JDD(II)-JDD(1)) ! sut
		TT=(ZA(II)-ZA(1))/86400.+TT
		write(12,'(1X,F12.6,2F12.3)') TT,Y(II),WEG(II)
 100	CONTINUE
	WRITE(IO,75)
		if(IND.EQ.1) close(12)
**
 75     FORMAT(2X,74('_'))
 76     FORMAT(2X,74('-'),/,22X
     &	,'pulsar  PSR ',A7)
 79     FORMAT(6X,'Date',3X,'JD-2400000.5',4X,'UTC[sec]',4X,'Err[mcs]',
     &  4X,'TB[sec]',6X,'Res[mcs ]')
	END


*======================================================================

	SUBROUTINE GAUSSD(NM,DEGP1,RHS,SUM,ONE,COEF)

*
*      Solution of system eq. using Gauss method.
*
* (INPUT DATA):
* -----------  NM - MAIN DIMENSION for SUM,ONE,RHS,COEF.
*              DEGR1-Degree of system
*              RHS(DEGP1)-Row of right-hand
*              SUM(DEGP1,DEGP1)-Matrix of Sist.N.Eq.
* (OUTPUT DATA):
* ------------ ONE(DEGP1,DEGP1)-Inverse matrix.
*              COEF(DEGP1)-Row of solutions.
*              SUM(DEGP1,DEGP1)*COEF(DEGP1)=RHS(DEGP1)
*              ONE(I,I)=1/P , where -weight of I-th term

	implicit none	
	INTEGER*4 NM,I,J,K,KPLUS1,L,IPLUS1,JJ
	
	REAL*8 RHS,SUM,ONE,COEF,TOT1
	REAL*8 ABS,DUMP,FACTOR,TOTAL,XX
	INTEGER*4 DEGP1,DEGREE
	DIMENSION RHS(NM),COEF(NM),SUM(NM,NM),ONE(NM,NM)
	DIMENSION TOT1(200) ! check if order maybe>200
 
	ABS(XX)=DABS(XX)
	DEGREE=DEGP1-1
 	
 	DO 3 I=1,DEGP1
	DO 3 J=1,DEGP1
 3	ONE(I,J)=0.0D0
	DO 4 I=1,DEGP1
 4	ONE(I,I)=1.0D0
 	
	DO 50 K=1,DEGREE
*	
*	  carry out the elimination process degree times

	KPLUS1=K+1
	L=K
	
*	   find terms of maximal value 

	DO 10 I=KPLUS1,DEGP1
	IF(ABS(SUM(I,K)).LE.ABS(SUM(L,K))) GOTO 10
	L=I
 10	CONTINUE

	IF(L.LE.K) GOTO 30

*	   if terms are already ordered then omit interchange

	DO 20 J=K,DEGP1
	DUMP=	SUM(L,J)
	SUM(L,J)=SUM(K,J)
	SUM(K,J)=DUMP
 20	CONTINUE

	DO 25 J=1,DEGP1
	DUMP=ONE(L,J)
	ONE(L,J)=ONE(K,J)
	ONE(K,J)=DUMP
 25	CONTINUE

	DUMP=RHS(K)
	RHS(K)=RHS(L)
	RHS(L)=DUMP

 30	DO 50 I=KPLUS1,DEGP1

*	   find factor which will eliminate term after
*	   subtraction.

	FACTOR=SUM(I,K)/SUM(K,K)
	SUM(I,K)=0.

	DO 40 J=KPLUS1,DEGP1

*	   compute other terms

 40	SUM(I,J)=SUM(I,J)-FACTOR*SUM(K,J)

	RHS(I)=RHS(I)-FACTOR*RHS(K)

	DO 45 J=1,DEGP1
	ONE(I,J)=ONE(I,J)-FACTOR*ONE(K,J)
 45	CONTINUE

 50	CONTINUE

*   	   compute solution of single eqation remaining

	COEF(DEGP1)=RHS(DEGP1)/SUM(DEGP1,DEGP1)

	DO 55 J=1,DEGP1
 55	ONE(DEGP1,J)=ONE(DEGP1,J)/SUM(DEGP1,DEGP1)

	I=DEGREE
 60	IPLUS1=I+1
	TOTAL=0.

*	   computes other solution by substitution

	DO 65 J=1,DEGP1
 65	TOT1(J)=0.

	DO 70 J=IPLUS1,DEGP1
	TOTAL=TOTAL+SUM(I,J)*COEF(J)
 70	CONTINUE

	DO 73 J=1,DEGP1
	DO 72 JJ=IPLUS1,DEGP1
 72	TOT1(J)=SUM(I,JJ)*ONE(JJ,J)+TOT1(J)
 73	CONTINUE

* 	   return solution in the vector 'coef'

	COEF(I)=(RHS(I)-TOTAL)/SUM(I,I)

	DO 75 J=1,DEGP1
 75	ONE(I,J)=(ONE(I,J)-TOT1(J))/SUM(I,I)

	I=I-1
	IF(I.GT.0) GOTO 60
	RETURN
	END

*======================================================================

*	Calculation of Matrix of correlations from
*	Covariance matrix.
*
* O.D.
	
	SUBROUTINE ERR(NM,N,NDAT,COEF,ONE,
     &		ER,RHS,R0,DISP)
        implicit none    
	INTEGER*4 NM,N,NDAT,I,J
	REAL*8 ONE(NM,NM),ER(1),RHS(1),COEF(1),R0
	REAL*8 DISP,D0,RR,AA,ON
	
	D0=-R0
	
	DO 10 I=1,N
 10	D0=RHS(I)*COEF(I)+D0
	
	RR=NDAT-N
	D0=D0/RR
	DISP=DSQRT(DABS(D0))

	DO 20 I=1,N
	AA=DABS(ONE(I,I)*D0)
 20	ER(I)=DSQRT(AA)

	DO 30 I=1,N
	DO 30 J=1,N
	ON=ONE(I,J)*D0
 30	ONE(I,J)=ON/ER(I)/ER(J)
	
	RETURN
	END


*======================================================================

	SUBROUTINE GEO(FI,LA,H,X,Y,Z)

*	
* Transformation from geodetical coordinates
* of an observatory (FI,LA,H) to topocentric
* rectangular coordinates (X,Y,Z) using
* standard Merit-83
*
* (INPUT)  FI,LA - geod.coord.,(rad)
* ------   H - elevation, (m.)
* (OUTPUT) X,Y,Z - topoc.coord.,(s) 
* -------
*
* O.D 1989

        implicit none
	REAL*8 FI,LA,H,X,Y,Z,A,B,C,CC,RN

	C=299792458.0D0
*(Merit-83)
	A=6378137.0D0
	B=6356752.3141D0
	CC=(DCOS(FI)*A)**2+(DSIN(FI)*B)**2
	RN=A**2/DSQRT(CC)
	X=(RN+H)*DCOS(FI)*DCOS(LA)
	Y=(RN+H)*DCOS(FI)*DSIN(LA)
	Z=((B/A)**2*RN+H)*DSIN(FI)
*(shift c-o-m according Merit-83)
c	X=X-80.0D0
c	Y=Y-101.0D0
c	Z=Z-119.0D0
	X=X/C
	Y=Y/C
	Z=Z/C
	RETURN
	END
**

*======================================================================

	SUBROUTINE XYZ2GEO(x,y,z,fi,la,h)
	
* Transforms rectangular coordinates (x,y,z) of observatory
* to geodetical coordinates (fi,la,h).
* For h < 1000m: error in h < 0.3cm and error in fi: < 0.003 arcsec
* For h < 3000m: error in h < 3cm   and error in fi: < 0.004 arcsec
*
* (input)
* x,y,z - rectangular coordinates of site, m
* (output)
* phi - latitude, rad
* la  - east longitude, rad
* h   - height (Sea level), m

* (Norbert Wex)

	implicit none
	real*8 a,b,q,x,y,z,fi,la,h,r
	data a,b,q/6378137.0d0,6356752.3141d0,0.003358431d0/

	r  = DSQRT(x*x+y*y+z*z)
	fi = DASIN(z/r)
	la = DATAN2(y,x)

	h=r-a*b/DSQRT((a*DSIN(fi))**2+(b*DCOS(fi))**2)
	fi=fi+q*(1.0d0-h/a)*DSIN(2.0*fi)+0.5d0*q**2*DSIN(4.0d0*fi)
	
	RETURN
	END

*===================================================================
        SUBROUTINE DTTOP
     &  (JD,UTC,dUT1,XP,YP,ZP,X,Y,Z,TTS,SC,DXP,DYP,GMAST,RTOPO)
*
*    Topocentric Reomer delay (TTS) and Doppler term (SC)
*    for a given source with coordinates (XP,YP,ZP) is
*    calculated taking into account:
*    1) diurnal rotation of an observer (X,Y,Z), 
*    2) precessional
*    3) nutational
*    4) Earth Orientation Parameters
* 
*(INPUT)  JD  - Julian date (days) at midnight
*------   UTC - time in UTC scale (s) 
*         XP,YP,ZP - rect.(Equatorial) BSS coord.of PSR (1)
*         X,Y,Z - rect.coord. of an observer (in light-sec). 
*         DXP,DYP - coordinates of Earth's pole (radians)
*(OUTPUT) TTS - Topoc.Roemer term (sec)   
* ------  GMAST - Grinvitch apparent sideral time
*         SC - Doppler term =(Vobs*Rpsr)/c (1)
*         RTOPO - topocentric coordinates of site in BRS
* O.D.

        implicit none
	REAL*8 TABOUT,DNUT,JD,UTC,dUT1,
     &	XP,YP,ZP,X,Y,Z,T0,T2,T3,DZITA,TETA,ZET,S,O(3),
     &	P(3,3),N(3,3),SY(3),AS,EPS,RJD,G,TTS,
     &  WR,VX,VY,SC,DXP,DYP,GMAST,RTOPO(3),
     &	sin dz,cos dz,sin z,cos z,sin th,cos th,
     &	sin ps,cos ps,sin ob,cos ob,sin ep,cos ep,
     &	cxp,sxp,cyp,syp,cs,ss
	INTEGER*4 J,L
	COMMON/CETBL4/TABOUT(6,12),DNUT(4)
**	
	AS=206264.806D0
	WR=7.29211514667D-05
	T0 = (JD-2451545.0D0)/36525.0D0
	
* (Lieske parameters)
	
	T2=T0*T0
	T3=T2*T0
	DZITA= 1.7998D-02*T3+3.0188D-01*T2+2306.2181D0*T0
	ZET  = 1.8203D-02*T3+1.09468D0 *T2+2306.2181D0*T0
	TETA =-4.1833D-02*T3-4.2665D-01*T2+2004.3109D0*T0
	
	DZITA = (-DZITA)/AS
	ZET   = (-ZET) /AS
	TETA  = (-TETA)/AS
	
*(Calculate inverse Precession matrix 1/P(Dzita,Zet,Teta)

	sin dz =DSIN(DZITA)
	cos dz =DCOS(DZITA)
	sin z  =DSIN(ZET)
	cos z  =DCOS(ZET)
	sin th =DSIN(TETA)
	cos th =DCOS(TETA)

	P(1,1) =-sin dz * sin z +cos dz * cos z * cos th
	P(2,1) = sin dz * cos z +cos dz * sin z * cos th
	P(3,1) = cos dz * sin th
	P(1,2) =-cos dz * sin z -sin dz * cos z * cos th
	P(2,2) = cos dz * cos z -sin dz * sin z * cos th
	P(3,2) =-sin dz * sin th
	P(1,3) =-cos z * sin th
	P(2,3) =-sin z * sin th
	P(3,3) = cos th

*(True obliquity for the date)

	EPS=1.813D-03*T3-0.00059D0*T2-46.815D0*T0+84381.448D0
	EPS=EPS/AS

*(Grinvitch apparent sideral Time s=GMAT for the time (JD,UTC))

	RJD= T0*36525.0D0+UTC/86400.0D0
	S  = 236.555367908D0*RJD+(0.093104D0
     &	   - 6.2D-06*(RJD/36525.0D0))*(RJD/36525.0D0)**2
	S  = (S+UTC+dUT1)/3600.0D0+
     &       (6.0D0+41.0D0/60.0D0+50.54841D0/3600.0D0)

*(Grinv.mean app.sid.time GMAST)
*(NB:local sideral time will be =s+Lambda)
	
        S  = S+DNUT(1)*AS*DCOS(EPS)/15.0d0/3600.0
	
	G=DINT(S/24.0D0)
	S=S-G*24.0D0
	IF(S.LT.0.0D0) S=S+24.0D0

	S=S/12.0D0*3.141592653589D0
        GMAST=S

*(Inverse matrix of Nutation N)

	sin ps =DSIN(DNUT(1))
	cos ps =DCOS(DNUT(1))
	sin ob =DSIN(EPS+DNUT(2))
	cos ob =DCOS(EPS+DNUT(2))
	sin ep =DSIN(EPS)
	cos ep =DCOS(EPS)

	N(1,1)= cos ps
	N(2,1)=-sin ps * cos ep
	N(3,1)=-sin ps * sin ep
	N(1,2)= sin ps * cos ob
	N(2,2)= cos ep * cos ob * cos ps + sin ep * sin ob
	N(3,2)= sin ep * cos ps * cos ob - cos ep * sin ob
	N(1,3)= sin ps * sin ob
	N(2,3)= cos ep * cos ps * sin ob - sin ep * cos ob
	N(3,3)= sin ep * cos ps * sin ob + cos ep * cos ob

* RTOPO is topocentric coordinates of observatory in BSS frame.
* (calculate R-vector in BSS as RTOPO(1,2,3)=(1/P)*(1/N)*SY(1,2,3)
* (i.e.product of inverse matrixes of precession P,)
* (nutation N, and vector of observer SY)

	cxp = dcos(DXP)
	sxp = dsin(DXP)
	cyp = dcos(DYP)
	syp = dsin(DYP)
        cs  = dcos(S)
        ss  = dsin(S)

	SY(1)= X*(cxp*cs-sxp*syp*ss) -Y*cyp*ss -Z*(sxp*cs-cxp*syp*ss)
	SY(2)= X*(cxp*ss+sxp*syp*cs) +Y*cyp*cs +Z*(-sxp*ss+cxp*syp*cs)
	SY(3)= X*sxp*cyp -Y*syp +Z*cxp*cyp

	do 30 J=1,3
	S=0.0D0
	do 25 L=1,3
   25	S=N(J,L)*SY(L)+S
   30	O(J)=S

	do 40 J=1,3
	S=0.0D0
	do 35 L=1,3
   35	S=P(J,L)*O(L)+S
	RTOPO(J)=S
   40	continue

*(Projection of velocity of observer onto the l-o-s of PSR)

c	VX=-WR* SY(2)
c	VY= WR* SY(1)
	
	vx=-wr*rtopo(2)
	vy= wr*rtopo(1)
	SC= VX*XP+VY*YP
	
*(Scalar product of RTOPO(1,2,3) and (XP,YP,ZP)
	
  	TTS=RTOPO(1)*XP+RTOPO(2)*YP+RTOPO(3)*ZP

	RETURN
	END

*======================================================================

      subroutine sortord(ar,ind,nopts)
      implicit none
c
c returns array of pointers IND(NOPTS) showing
c the sorting order of array AR(NOPTS).
c INPUT AR(NOPTS),NOPTS
c OUTPUT IND(NOPTS)
c Attention: AR(NOPTS) will be sorted too.

      real*4 ar(1),arless,arsav
      integer*4 nopts,j,i,numb,indsav
      integer*2 ind(1)
      do j=1,nopts
      ind(j)=j
      enddo

      do j=1,nopts
      arless=ar(j)
      numb=j

      do i=j,nopts
      arless= amin1(arless,ar(i))
      if (arless.eq.ar(i)) then
         numb=i
         endif
      enddo

      arsav=ar(numb)
      ar(numb)=ar(j)
      ar(j)=arsav

      indsav=ind(numb)
      ind(numb)=ind(j)
      ind(j)=indsav

      enddo

      return
      end

c====================================================================

        subroutine geteop(TOA,IndORD,Npoint,DXP,DYP,DUT1,lack,fileop)

* Interpolate EOP parameters for moments of TOAs.
* The values are read from file 'eop.dat', generated
* from original IERS data eopc04.XX
*
* INPUT:
* ------
* TOA(Npoint) - MJD.DD - toas array (must be sorted)
* IndORD(Npoint) - array of pointers to the original
*                  (unsorted) array of TOAs
* Npoint - number of TOAs
* fileop - full path to file 'eop.dat'
*
* OUTPUT:
* ------
* lack-flag:
* lack=0 - successful run
* lack=1 - no data or there is no file 'eop.dat'
* lack>1 - not enough data in 'eop.dat', starting from TOA(lack)
*
* dXp,dYp,dUT1 - Interpolated values of Earth Orientation parameters:
* DXp(Npoint)
* DYp(NPOINT)
* DUT1(Npoint)
*


        implicit none
        integer*4 datc0,Nread,Npoint,Lack,i
	integer*2 IndORD(1)
	real*4 N,X0,Y0,X1,X2,X3,Y1,Y2,Y3,u2,u3,u1,u0
	real*4 TOA(1),DXP(1),DYP(1),DUT1(1),d0,dat,datc,a,b,c
	character fileop*80

        do i=1,Npoint
        DXP(i)=0.0
        DYP(i)=0.0
        DUT1(i)=0.0
        enddo

	lack=0
	X2=0
	Y2=0
	X3=0
	Y3=0
	u2=0
	u3=0
	datc=0
	Nread=0
	
	open(10,file=fileop,status='OLD',err=40)
        go to 30
 40     lack=1
        return

*(Go through points)

 30     continue
	do i=1,Npoint
	dat = TOA(i)
	
	if (datc.GE.dat) goto 20

 10	X1=X2
	Y1=Y2
	u1=u2
	X2=X3
	Y2=Y3
	u2=u3
	
	d0=datc
	read(10,'(i6,2F9.6,f9.6)',end=300) datc0,X3,Y3,u3
        goto 33
 300	lack=i
	close(10)
        return
        
 33	datc=float(datc0)
	if (datc.GE.dat) goto 20
	goto 10
	
 20	continue

* Interpolate X0,Y0 on DAT

        n=dat-d0
* X0
	a=X2-X1
	b=X3-X2
	c=b-a
	X0=X2+n/2.*(a+b+n*c)
* Y0
	a=Y2-Y1
	b=Y3-Y2
	c=b-a
	Y0=Y2+n/2.*(a+b+n*c)
* u0
	a=u2-u1
	b=u3-u2
	c=b-a
	u0=u2+n/2.*(a+b+n*c)
	
        DXP(IndORD(i)) =X0
        DYP(IndORD(i)) =Y0
        DUT1(IndORD(i))=u0

        enddo ! i
        
 1000	close(10)
        return
	end

c====================================================================

	real*4 FUNCTION RNDM(I)
	integer*4 i
	L=1
	OUT=0.
	IF(i-1) 2,1,2
 1	U=.37843
 2	IF(L-1) 4,5,5
 4	K=1
	GOTO 6
 5	K=12
 6	DO 3 J=1,K
	F=37.*U
	U=F-AINT(F)
 3	OUT=OUT+U-.5
c	if(OUT.LT.0.0) OUT=-OUT
	RNDm=OUT
	RETURN
	END

*======================================================================

		SUBROUTINE READEL405(JED,TSEC,IERR,fileph)

*
* Reading and interpolation of the JPL ephemerides DE405/LE405
* which are given as coefficents of Chebyshev's polinomas.
*
*
* For timing purposes the binary data file de405.bin is used
* with direct access of data records. Check the OPEN statement
* where you should specify correct record length (RECL=..) -
* it is compiler dependent. The record length should be equivalent
* to the length of (1018 x real*8) numbers.
*				 
* The center of the reference frame BSS is choosen by setting ICENT=0
* Other values of ICENT may choosen according to values of 
* index (i), i=1..11 (see below)
*
* The choice of planets which positions/velocities to be calculated
* is made by setting of elements of array IREQ(i), i=1..13:
*
* IREQ(i)=0 - positions/velocities will not be calculated
* IREQ(i)=1 - positions will be calculated
* IREQ(i)=2 - velocities will be calculated
* 
* Index (i) refers to:
*
* 1-Mercury; 2-Venus; 3-Earth; 4-Mars; 5-Jupiter; 6-Saturn;
* 7-Uran; 8-Neptune; 9-Pluto; 10-Sun; 11-Moon; 12-Earth/Moon Barycenter;
* 13-Nutations
*
* The values of velocities/positions returned in array TABOUT(k,i)
* k=1,2,3 - X,Y,Z; k=4,5,6 - Xdot,Ydot,Zdot [AU and AU/DAY].
* i=1..12 - see index (i) above.
* The nutation coefficents will be returned in DNUT(i)
*
* (INPUT) JED, TSEC - time (TB) for which positions/velocities to be 
*                     calculated (JED - at midnight, TSEC - in seconds)
*
* The program transfer data through common blocks CETBL1..CETBL6
*
*


	implicit none
	INTEGER*4 LIST(11),LOC(3,12)
	REAL*8 JED,TSEC,DJE,BUF405,PV,PVSUN,SS,T,EM,EMR,AUFAC,
     &	BF,TABOUT,NUT,DEJ
	REAL*4 EMRAT
	INTEGER*4 ICW,ICENT,IREQ,NZP,I,LST,IC,NR,NSK,LL,NM,J,ierr
	DIMENSION PV(6,13),PVSUN(6),SS(6),T(2)
	character fileph*80
	COMMON/CETBL1/EMRAT
	COMMON/CETBL2/ICW,ICENT,IREQ(13)
c	COMMON/CETBL3/BUF(826)
	COMMON/CETBL4/TABOUT(6,12),NUT(4)
	COMMON/NZP/NZP
	COMMON/CETBL5/BUF405(1018)
	
	DATA LOC/3,14,4,171,10,2,231,13,2,309,11,1,342,8,1,
     &	366,7,1,387,6,1,405,6,1,423,6,1,441,13,8,753,11,2,
     &	819,10,4/,SS/2436912.5D0,2473424.5D0,32.0D0,32.0D0,
     &	12.0D0,4.0D0/,BF/149597870.691D0/

* Note that above SS(1) and SS(2) are the first and the last
* Julian date of ephemerides file DE405.BIN

     	IERR=0
	EM=EMRAT
	IF(EM.EQ.0.0D0)EM=81.30056D0
	EMR=1.0D0/(1.0D0+EM)
	DJE=TSEC/86400.0D0
	IF(ICW.EQ.2) GOTO 1
	ICW=2
	BUF405(1)=0.0D0
	BUF405(2)=0.0D0
	NZP=0
 1	CONTINUE
	IF(ICENT.GE.0.AND.ICENT.LE.11) GOTO 2
	IERR=4
	RETURN
 2	CONTINUE
	DO 3 I=1,9
	IF(IREQ(I).GE.0.AND.IREQ(I).LE.2) GOTO 3
	IERR=3
	RETURN
 3	LIST(I)=IREQ(I)
	LIST(10)=IREQ(11)
	IF(ICENT.EQ.0.OR.ICENT.EQ.10) GOTO 5
	LST=0
	DO 4 I=1,10
 4	LST=MAX0(LST,LIST(I))
	IC=ICENT-ICENT/11
	LIST(IC)=MAX0(LST,IREQ(10))
 5	CONTINUE
	LIST(3)=MAX0(LIST(3),LIST(10))
	LIST(3)=MAX0(LIST(3),IREQ(12))
	LIST(10)=LIST(3)
	LIST(11)=IREQ(13)
	T(2)=SS(3)
	AUFAC=1.0D0/BF
	DEJ=JED+DJE
	IF(DEJ.LT.SS(1)) GOTO 99
	IF(DEJ.GT.SS(2)) GOTO 98
	NR=IDINT((DEJ-SS(1))/SS(4))+1


	NSK=NR-NZP
	IF(NSK) 18,19,18
  18	continue
	open(unit=12,file=fileph,access='DIRECT',
     &	status='OLD',form='UNFORMATTED',recl=8144,err=45)
*     &	blocksize=32768)
	goto 46
 45	write(*,*) ' ERROR: Not found',fileph
	write(*,*) ' Check config file ptimer.cfg'
	stop
	
 46	READ(12,REC=NR) BUF405
	close(unit=12)
 19	CONTINUE

	NZP=NR
 6	IF(DEJ.GE.BUF405(1).AND.DEJ.LE.BUF405(2)) GOTO 7
	IERR=5
	RETURN
 7	T(1)=(DEJ-BUF405(1))/T(2)
	LL=LOC(1,11)
	CALL XINTRP(BUF405(LL),T,LOC(2,11),3,LOC(3,11),2,PVSUN)
	DO 9 I=1,10
	IF(LIST(I).EQ.0) GOTO 9
	LL=LOC(1,I)
	CALL XINTRP(BUF405(LL),T,LOC(2,I),3,LOC(3,I),LIST(I),
     &	PV(1,I))
	NM=LIST(I)*3
	DO 8 J=1,NM
	IF(I.LE.9.AND.ICENT.NE.0)PV(J,I)=(PV(J,I)-PVSUN(J))
     &	*AUFAC
	IF(I.LE.9.AND.ICENT.EQ.0)PV(J,I)=PV(J,I)*AUFAC
	IF(I.EQ.10)PV(J,I)=PV(J,I)*AUFAC
 8	CONTINUE
 9	CONTINUE
	DO 10 J=1,6
 10	PVSUN(J)=PVSUN(J)*AUFAC
	IF(LIST(11).EQ.0) GOTO 11
	LL=LOC(1,12)
	CALL XINTRP(BUF405(LL),T,LOC(2,12),2,LOC(3,12),
     &	LIST(11),NUT)
	GOTO 11
 99	CONTINUE
	IERR=1
	RETURN
 98	IERR=2
	RETURN
 11	IF(IERR.NE.0) RETURN
	DO 12 I=1,6
	PV(I,11)=PV(I,3)+EM*EMR*PV(I,10)
	PV(I,12)=PV(I,3)
	PV(I,3)=PV(I,3)-EMR*PV(I,10)
	PV(I,10)=0.0D0
 12	CONTINUE
	IF(ICENT.EQ.0) GOTO 14
	DO 13 I=1,12
	DO 13 J=1,3
	IF(IREQ(I).EQ.0) GOTO 13
	TABOUT(J,I)=PV(J,I)-PV(J,ICENT)
	IF(IREQ(I).EQ.1) GOTO 13
	TABOUT(J+3,I)=PV(J+3,I)-PV(J+3,ICENT)
 13	CONTINUE
	RETURN
 14	DO 15 I=1,12
	DO 15 J=1,3
	IF(IREQ(I).EQ.0) GOTO 15
	TABOUT(J,I)=PV(J,I)
	IF(I.EQ.10)TABOUT(J,10)=PVSUN(J)
	IF(IREQ(I).EQ.1) GOTO 15
	TABOUT(J+3,I)=PV(J+3,I)
	IF(I.EQ.10)TABOUT(J+3,10)=PVSUN(J+3)
 15	CONTINUE
	RETURN
	END
	
*======================================================================
      SUBROUTINE TDBTRANS(tdt,tdtsec,tdbsec)

c     Subroutine transforms TDT into TDB
c     enter: TDT in days, integer (tdt) and float (tdtsec) parts
c     result: tdbsec = tdtsec+ctatv (CTATV = TDB - TDT)

      implicit real*8(a-h,o-z)
      T=((tdt-2451545.d0)+tdtsec)/(1000.d0*365.25d0) 
      TT = T*T
      TTT = T*TT
      TTTT = T*TTT

c T**0
      t1= 1656.674564 * SIN(   6283.075849991 *T + 6.240054195 )+
     & 22.417471 * SIN(   5753.384884897 *T + 4.296977442 )  +
     & 13.839792 * SIN(  12566.151699983 *T + 6.196904410 )  +
     & 4.770086 * SIN(    529.690965095 *T + 0.444401603 )   +
     & 4.676740 * SIN(   6069.776754553 *T + 4.021195093 )   +
     & 2.256707 * SIN(    213.299095438 *T + 5.543113262 )   +
     & 1.694205 * SIN(     -3.523118349 *T + 5.025132748 )   +
     & 1.554905 * SIN(  77713.771467920 *T + 5.198467090 )   +
     & 1.276839 * SIN(   7860.419392439 *T + 5.988822341 )   +
     & 1.193379 * SIN(   5223.693919802 *T + 3.649823730 )   +
     & 1.115322 * SIN(   3930.209696220 *T + 1.422745069 )   +
     & 0.794185 * SIN(  11506.769769794 *T + 2.322313077 )   +
     & 0.447061 * SIN(     26.298319800 *T + 3.615796498 )   +
     & 0.435206 * SIN(   -398.149003408 *T + 4.349338347 )   +
     & 0.600309 * SIN(   1577.343542448 *T + 2.678271909 )   +
     & 0.496817 * SIN(   6208.294251424 *T + 5.696701824 )   +
     & 0.486306 * SIN(   5884.926846583 *T + 0.520007179 )   +
     & 0.432392 * SIN(     74.781598567 *T + 2.435898309 )   +
     & 0.468597 * SIN(   6244.942814354 *T + 5.866398759 )   +
     & 0.375510 * SIN(   5507.553238667 *T + 4.103476804 )   
      t2= 0.243085 * SIN(   -775.522611324 *T + 3.651837925 )   +
     & 0.173435 * SIN(  18849.227549974 *T + 6.153743485 )   +
     & 0.230685 * SIN(   5856.477659115 *T + 4.773852582 )   +
     & 0.203747 * SIN(  12036.460734888 *T + 4.333987818 )   +
     & 0.143935 * SIN(   -796.298006816 *T + 5.957517795 )   +
     & 0.159080 * SIN(  10977.078804699 *T + 1.890075226 )   +
     & 0.119979 * SIN(     38.133035638 *T + 4.551585768 )   +
     & 0.118971 * SIN(   5486.777843175 *T + 1.914547226 )   +
     & 0.116120 * SIN(   1059.381930189 *T + 0.873504123 )   +
     & 0.137927 * SIN(  11790.629088659 *T + 1.135934669 )   +
     & 0.098358 * SIN(   2544.314419883 *T + 0.092793886 )   +
     & 0.101868 * SIN(  -5573.142801634 *T + 5.984503847 )   +
     & 0.080164 * SIN(    206.185548437 *T + 2.095377709 )   +
     & 0.079645 * SIN(   4694.002954708 *T + 2.949233637 )   +
     & 0.062617 * SIN(     20.775395492 *T + 2.654394814 )   +
     & 0.075019 * SIN(   2942.463423292 *T + 4.980931759 )   +
     & 0.064397 * SIN(   5746.271337896 *T + 1.280308748 )   +
     & 0.063814 * SIN(   5760.498431898 *T + 4.167901731 )   +
     & 0.048042 * SIN(   2146.165416475 *T + 1.495846011 )   +
     & 0.048373 * SIN(    155.420399434 *T + 2.251573730 )   
      t3= 0.058844 * SIN(    426.598190876 *T + 4.839650148 )   +
     & 0.046551 * SIN(     -0.980321068 *T + 0.921573539 )   +
     & 0.054139 * SIN(  17260.154654690 *T + 3.411091093 )   +
     & 0.042411 * SIN(   6275.962302991 *T + 2.869567043 )   +
     & 0.040184 * SIN(     -7.113547001 *T + 3.565975565 )   +
     & 0.036564 * SIN(   5088.628839767 *T + 3.324679049 )   +
     & 0.040759 * SIN(  12352.852604545 *T + 3.981496998 )   +
     & 0.036507 * SIN(    801.820931124 *T + 6.248866009 )   +
     & 0.036955 * SIN(   3154.687084896 *T + 5.071801441 )   +
     & 0.042732 * SIN(    632.783739313 *T + 5.720622217 )   +
     & 0.042560 * SIN( 161000.685737473 *T + 1.270837679 )   +
     & 0.040480 * SIN(  15720.838784878 *T + 2.546610123 )   +
     & 0.028244 * SIN(  -6286.598968340 *T + 5.069663519 )   +
     & 0.033477 * SIN(   6062.663207553 *T + 4.144987272 )   +
     & 0.034867 * SIN(    522.577418094 *T + 5.210064075 )   +
     & 0.032438 * SIN(   6076.890301554 *T + 0.749317412 )   +
     & 0.030215 * SIN(   7084.896781115 *T + 3.389610345 )   +
     & 0.029247 * SIN( -71430.695617928 *T + 4.183178762 )   +
     & 0.033529 * SIN(   9437.762934887 *T + 2.404714239 )   +
     & 0.032423 * SIN(   8827.390269875 *T + 5.541473556 )   
      t4= 0.027567 * SIN(   6279.552731642 *T + 5.040846034 )   +
     & 0.029862 * SIN(  12139.553509107 *T + 1.770181024 )   +
     & 0.022509 * SIN(  10447.387839604 *T + 1.460726241 )   +
     & 0.020937 * SIN(   8429.241266467 *T + 0.652303414 )   +
     & 0.020322 * SIN(    419.484643875 *T + 3.735430632 )   +
     & 0.024816 * SIN(  -1194.447010225 *T + 1.087136918 )   +
     & 0.025196 * SIN(   1748.016413067 *T + 2.901883301 )   +
     & 0.021691 * SIN(  14143.495242431 *T + 5.952658009 )   +
     & 0.017673 * SIN(   6812.766815086 *T + 3.186129845 )   +
     & 0.022567 * SIN(   6133.512652857 *T + 3.307984806 )   +
     & 0.016155 * SIN(  10213.285546211 *T + 1.331103168 )   +
     & 0.014751 * SIN(   1349.867409659 *T + 4.308933301 )   +
     & 0.015949 * SIN(   -220.412642439 *T + 4.005298270 )   +
     & 0.015974 * SIN(  -2352.866153772 *T + 6.145309371 )   +
     & 0.014223 * SIN(  17789.845619785 *T + 2.104551349 )   +
     & 0.017806 * SIN(     73.297125859 *T + 3.475975097 )   +
     & 0.013671 * SIN(   -536.804512095 *T + 5.971672571 )   +
     & 0.011942 * SIN(   8031.092263058 *T + 2.053414715 )   +
     & 0.014318 * SIN(  16730.463689596 *T + 3.016058075 )   +
     & 0.012462 * SIN(    103.092774219 *T + 1.737438797 )   
      t5= 0.010962 * SIN(      3.590428652 *T + 2.196567739 )   +
     & 0.015078 * SIN(  19651.048481098 *T + 3.969480770 )   +
     & 0.010396 * SIN(    951.718406251 *T + 5.717799605 )   +
     & 0.011707 * SIN(  -4705.732307544 *T + 2.654125618 )   +
     & 0.010453 * SIN(   5863.591206116 *T + 1.913704550 )   +
     & 0.012420 * SIN(   4690.479836359 *T + 4.734090399 )   +
     & 0.011847 * SIN(   5643.178563677 *T + 5.489005403 )   +
     & 0.008610 * SIN(   3340.612426700 *T + 3.661698944 )   +
     & 0.011622 * SIN(   5120.601145584 *T + 4.863931876 )   +
     & 0.010825 * SIN(    553.569402842 *T + 0.842715011 )   +
     & 0.008666 * SIN(   -135.065080035 *T + 3.293406547 )   +
     & 0.009963 * SIN(    149.563197135 *T + 4.870690598 )   +
     & 0.009858 * SIN(   6309.374169791 *T + 1.061816410 )   +
     & 0.007959 * SIN(    316.391869657 *T + 2.465042647 )   +
     & 0.010099 * SIN(    283.859318865 *T + 1.942176992 )   +
     & 0.007147 * SIN(   -242.728603974 *T + 3.661486981 )   +
     & 0.007505 * SIN(   5230.807466803 *T + 4.920937029 )   +
     & 0.008323 * SIN(  11769.853693166 *T + 1.229392026 )   +
     & 0.007490 * SIN(  -6256.777530192 *T + 3.658444681 )   +
     & 0.009370 * SIN( 149854.400134205 *T + 0.673880395 )   
      t6= 0.007117 * SIN(     38.027672636 *T + 5.294249518 )   +
     & 0.007857 * SIN(  12168.002696575 *T + 0.525733528 )   +
     & 0.007019 * SIN(   6206.809778716 *T + 0.837688810 )   +
     & 0.006056 * SIN(    955.599741609 *T + 4.194535082 )   +
     & 0.008107 * SIN(  13367.972631107 *T + 3.793235253 )   +
     & 0.006731 * SIN(   5650.292110678 *T + 5.639906583 )   +
     & 0.007332 * SIN(     36.648562930 *T + 0.114858677 )   +
     & 0.006366 * SIN(   4164.311989613 *T + 2.262081818 )   +
     & 0.006858 * SIN(   5216.580372801 *T + 0.642063318 )   +
     & 0.006919 * SIN(   6681.224853400 *T + 6.018501522 )   +
     & 0.006826 * SIN(   7632.943259650 *T + 3.458654112 )   +
     & 0.005308 * SIN(  -1592.596013633 *T + 2.500382359 )   +
     & 0.005096 * SIN(  11371.704689758 *T + 2.547107806 )   +
     & 0.004841 * SIN(   5333.900241022 *T + 0.437078094 )   +
     & 0.005582 * SIN(   5966.683980335 *T + 2.246174308 )   +
     & 0.006304 * SIN(  11926.254413669 *T + 2.512929171 )   +
     & 0.006603 * SIN(  23581.258177318 *T + 5.393136889 )   +
     & 0.005123 * SIN(     -1.484472708 *T + 2.999641028 )   +
     & 0.004648 * SIN(   1589.072895284 *T + 1.275847090 )   +
     & 0.005119 * SIN(   6438.496249426 *T + 1.486539246 )   
      t7= 0.004521 * SIN(   4292.330832950 *T + 6.140635794 )   +
     & 0.005680 * SIN(  23013.539539587 *T + 4.557814849 )   +
     & 0.005488 * SIN(     -3.455808046 *T + 0.090675389 )   +
     & 0.004193 * SIN(   7234.794256242 *T + 4.869091389 )   +
     & 0.003742 * SIN(   7238.675591600 *T + 4.691976180 )   +
     & 0.004148 * SIN(   -110.206321219 *T + 3.016173439 )   +
     & 0.004553 * SIN(  11499.656222793 *T + 5.554998314 )   +
     & 0.004892 * SIN(   5436.993015240 *T + 1.475415597 )   +
     & 0.004044 * SIN(   4732.030627343 *T + 1.398784824 )   +
     & 0.004164 * SIN(  12491.370101415 *T + 5.650931916 )   +
     & 0.004349 * SIN(  11513.883316794 *T + 2.181745369 )   +
     & 0.003919 * SIN(  12528.018664345 *T + 5.823319737 )   +
     & 0.003129 * SIN(   6836.645252834 *T + 0.003844094 )   +
     & 0.004080 * SIN(  -7058.598461315 *T + 3.690360123 )   +
     & 0.003270 * SIN(     76.266071276 *T + 1.517189902 )   +
     & 0.002954 * SIN(   6283.143160294 *T + 4.447203799 )   +
     & 0.002872 * SIN(     28.449187468 *T + 1.158692983 )   +
     & 0.002881 * SIN(    735.876513532 *T + 0.349250250 )   +
     & 0.003279 * SIN(   5849.364112115 *T + 4.893384368 )   +
     & 0.003625 * SIN(   6209.778724132 *T + 1.473760578 )   
      t8= 0.003074 * SIN(    949.175608970 *T + 5.185878737 )   +
     & 0.002775 * SIN(   9917.696874510 *T + 1.030026325 )   +
     & 0.002646 * SIN(  10973.555686350 *T + 3.918259169 )   +
     & 0.002575 * SIN(  25132.303399966 *T + 6.109659023 )   +
     & 0.003500 * SIN(    263.083923373 *T + 1.892100742 )   +
     & 0.002740 * SIN(  18319.536584880 *T + 4.320519510 )   +
     & 0.002464 * SIN(    202.253395174 *T + 4.698203059 )   +
     & 0.002409 * SIN(      2.542797281 *T + 5.325009315 )   +
     & 0.003354 * SIN( -90955.551694697 *T + 1.942656623 )   +
     & 0.002296 * SIN(   6496.374945429 *T + 5.061810696 )   +
     & 0.003002 * SIN(   6172.869528772 *T + 2.797822767 )   +
     & 0.003202 * SIN(  27511.467873537 *T + 0.531673101 )   +
     & 0.002954 * SIN(  -6283.008539689 *T + 4.533471191 )   +
     & 0.002353 * SIN(    639.897286314 *T + 3.734548088 )   +
     & 0.002401 * SIN(  16200.772724501 *T + 2.605547070 )   +
     & 0.003053 * SIN( 233141.314403759 *T + 3.029030662 )   +
     & 0.003024 * SIN(  83286.914269554 *T + 2.355556099 )   +
     & 0.002863 * SIN(  17298.182327326 *T + 5.240963796 )   +
     & 0.002103 * SIN(  -7079.373856808 *T + 5.756641637 )   +
     & 0.002303 * SIN(  83996.847317911 *T + 2.013686814 )   
      t9= 0.002303 * SIN(  18073.704938650 *T + 1.089100410 )   +
     & 0.002381 * SIN(     63.735898303 *T + 0.759188178 )   +
     & 0.002493 * SIN(   6386.168624210 *T + 0.645026535 )   +
     & 0.002366 * SIN(      3.932153263 *T + 6.215885448 )   +
     & 0.002169 * SIN(  11015.106477335 *T + 4.845297676 )   +
     & 0.002397 * SIN(   6243.458341645 *T + 3.809290043 )   +
     & 0.002183 * SIN(   1162.474704408 *T + 6.179611691 )   +
     & 0.002353 * SIN(   6246.427287062 *T + 4.781719760 )   +
     & 0.002199 * SIN(   -245.831646229 *T + 5.956152284 )   +
     & 0.001729 * SIN(   3894.181829542 *T + 1.264976635 )   +
     & 0.001896 * SIN(  -3128.388765096 *T + 4.914231596 )   +
     & 0.002085 * SIN(     35.164090221 *T + 1.405158503 )   +
     & 0.002024 * SIN(  14712.317116458 *T + 2.752035928 )   +
     & 0.001737 * SIN(   6290.189396992 *T + 5.280820144 )   +
     & 0.002229 * SIN(    491.557929457 *T + 1.571007057 )   +
     & 0.001602 * SIN(  14314.168113050 *T + 4.203664806 )   +
     & 0.002186 * SIN(    454.909366527 *T + 1.402101526 )   +
     & 0.001897 * SIN(  22483.848574493 *T + 4.167932508 )   +
     & 0.001825 * SIN(  -3738.761430108 *T + 0.545828785 )   +
     & 0.001894 * SIN(   1052.268383188 *T + 5.817167450 )   
      t10= 0.001421 * SIN(     20.355319399 *T + 2.419886601 )   +
     & 0.001408 * SIN(  10984.192351700 *T + 2.732084787 )   +
     & 0.001847 * SIN(  10873.986030480 *T + 2.903477885 )   +
     & 0.001391 * SIN(  -8635.942003763 *T + 0.593891500 )   +
     & 0.001388 * SIN(     -7.046236698 *T + 1.166145902 )   +
     & 0.001810 * SIN( -88860.057071188 *T + 0.487355242 )   +
     & 0.001288 * SIN(  -1990.745017041 *T + 3.913022880 )   +
     & 0.001297 * SIN(  23543.230504682 *T + 3.063805171 )   +
     & 0.001335 * SIN(   -266.607041722 *T + 3.995764039 )   +
     & 0.001376 * SIN(  10969.965257698 *T + 5.152914309 )   +
     & 0.001745 * SIN( 244287.600007027 *T + 3.626395673 )   +
     & 0.001649 * SIN(  31441.677569757 *T + 1.952049260 )   +
     & 0.001416 * SIN(   9225.539273283 *T + 4.996408389 )   +
     & 0.001238 * SIN(   4804.209275927 *T + 5.503379738 )   +
     & 0.001472 * SIN(   4590.910180489 *T + 4.164913291 )   +
     & 0.001169 * SIN(   6040.347246017 *T + 5.841719038 )   +
     & 0.001039 * SIN(   5540.085789459 *T + 2.769753519 )   +
     & 0.001004 * SIN(   -170.672870619 *T + 0.755008103 )   +
     & 0.001284 * SIN(  10575.406682942 *T + 5.306538209 )   +
     & 0.001278 * SIN(     71.812653151 *T + 4.713486491 )   
      t11= 0.001321 * SIN(  18209.330263660 *T + 2.624866359 )   +
     & 0.001297 * SIN(  21228.392023546 *T + 0.382603541 )   +
     & 0.000954 * SIN(   6282.095528923 *T + 0.882213514 )   +
     & 0.001145 * SIN(   6058.731054289 *T + 1.169483931 )   +
     & 0.000979 * SIN(   5547.199336460 *T + 5.448375984 )   +
     & 0.000987 * SIN(  -6262.300454499 *T + 2.656486959 )   +
     & 0.001070 * SIN(-154717.609887482 *T + 1.827624012 )   +
     & 0.000991 * SIN(   4701.116501708 *T + 4.387001801 )   +
     & 0.001155 * SIN(    -14.227094002 *T + 3.042700750 )   +
     & 0.001176 * SIN(    277.034993741 *T + 3.335519004 )   +
     & 0.000890 * SIN(  13916.019109642 *T + 5.601498297 )   +
     & 0.000884 * SIN(  -1551.045222648 *T + 1.088831705 )   +
     & 0.000876 * SIN(   5017.508371365 *T + 3.969902609 )   +
     & 0.000806 * SIN(  15110.466119866 *T + 5.142876744 )   +
     & 0.000773 * SIN(  -4136.910433516 *T + 0.022067765 )   +
     & 0.001077 * SIN(    175.166059800 *T + 1.844913056 )   +
     & 0.000954 * SIN(  -6284.056171060 *T + 0.968480906 )   +
     & 0.000737 * SIN(   5326.786694021 *T + 4.923831588 )   +
     & 0.000845 * SIN(   -433.711737877 *T + 4.749245231 )   +
     & 0.000819 * SIN(   8662.240323563 *T + 5.991247817 )   
      t12= 0.000852 * SIN(    199.072001436 *T + 2.189604979 )   +
     & 0.000723 * SIN(  17256.631536341 *T + 6.068719637 )   +
     & 0.000940 * SIN(   6037.244203762 *T + 6.197428148 )   +
     & 0.000885 * SIN(  11712.955318231 *T + 3.280414875 )   +
     & 0.000706 * SIN(  12559.038152982 *T + 2.824848947 )   +
     & 0.000732 * SIN(   2379.164473572 *T + 2.501813417 )   +
     & 0.000764 * SIN(  -6127.655450557 *T + 2.236346329 )   +
     & 0.000908 * SIN(    131.541961686 *T + 2.521257490 )   +
     & 0.000907 * SIN(  35371.887265976 *T + 3.370195967 )   +
     & 0.000673 * SIN(   1066.495477190 *T + 3.876512374 )   +
     & 0.000814 * SIN(  17654.780539750 *T + 4.627122566 )   +
     & 0.000630 * SIN(     36.027866677 *T + 0.156368499 )   +
     & 0.000798 * SIN(    515.463871093 *T + 5.151962502 )   +
     & 0.000798 * SIN(    148.078724426 *T + 5.909225055 )   +
     & 0.000806 * SIN(    309.278322656 *T + 6.054064447 )   +
     & 0.000607 * SIN(    -39.617508346 *T + 2.839021623 )   +
     & 0.000601 * SIN(    412.371096874 *T + 3.984225404 )   +
     & 0.000646 * SIN(  11403.676995575 *T + 3.852959484 )   +
     & 0.000704 * SIN(  13521.751441591 *T + 2.300991267 )   +
     & 0.000603 * SIN( -65147.619767937 *T + 4.140083146 )   
      t13= 0.000609 * SIN(  10177.257679534 *T + 0.437122327 )   +
     & 0.000631 * SIN(   5767.611978898 *T + 4.026532329 )   +
     & 0.000576 * SIN(  11087.285125918 *T + 4.760293101 )   +
     & 0.000674 * SIN(  14945.316173554 *T + 6.270510511 )   +
     & 0.000726 * SIN(   5429.879468239 *T + 6.039606892 )   +
     & 0.000710 * SIN(  28766.924424484 *T + 5.672617711 )   +
     & 0.000647 * SIN(  11856.218651625 *T + 3.397132627 )   +
     & 0.000678 * SIN(  -5481.254918868 *T + 6.249666675 )   +
     & 0.000618 * SIN(  22003.914634870 *T + 2.466427018 )   +
     & 0.000738 * SIN(   6134.997125565 *T + 2.242668890 )   +
     & 0.000660 * SIN(    625.670192312 *T + 5.864091907 )   +
     & 0.000694 * SIN(   3496.032826134 *T + 2.668309141 )   +
     & 0.000531 * SIN(   6489.261398429 *T + 1.681888780 )   +
     & 0.000611 * SIN(-143571.324284214 *T + 2.424978312 )   +
     & 0.000575 * SIN(  12043.574281889 *T + 4.216492400 )   +
     & 0.000553 * SIN(  12416.588502848 *T + 4.772158039 )   +
     & 0.000689 * SIN(   4686.889407707 *T + 6.224271088 )   +
     & 0.000495 * SIN(   7342.457780181 *T + 3.817285811 )   +
     & 0.000567 * SIN(   3634.621024518 *T + 1.649264690 )   +
     & 0.000515 * SIN(  18635.928454536 *T + 3.945345892 )   
      t14= 0.000486 * SIN(   -323.505416657 *T + 4.061673868 )   +
     & 0.000662 * SIN(  25158.601719765 *T + 1.794058369 )   +
     & 0.000509 * SIN(    846.082834751 *T + 3.053874588 )   +
     & 0.000472 * SIN( -12569.674818332 *T + 5.112133338 )   +
     & 0.000461 * SIN(   6179.983075773 *T + 0.513669325 )   +
     & 0.000641 * SIN(  83467.156352816 *T + 3.210727723 )   +
     & 0.000520 * SIN(  10344.295065386 *T + 2.445597761 )   +
     & 0.000493 * SIN(  18422.629359098 *T + 1.676939306 )   +
     & 0.000478 * SIN(   1265.567478626 *T + 5.487314569 )   +
     & 0.000472 * SIN(    -18.159247265 *T + 1.999707589 )   +
     & 0.000559 * SIN(  11190.377900137 *T + 5.783236356 )   +
     & 0.000494 * SIN(   9623.688276691 *T + 3.022645053 )   +
     & 0.000463 * SIN(   5739.157790895 *T + 1.411223013 )   +
     & 0.000432 * SIN(  16858.482532933 *T + 1.179256434 )   +
     & 0.000574 * SIN(  72140.628666286 *T + 1.758191830 )   +
     & 0.000484 * SIN(  17267.268201691 *T + 3.290589143 )   +
     & 0.000550 * SIN(   4907.302050146 *T + 0.864024298 )   +
     & 0.000399 * SIN(     14.977853527 *T + 2.094441910 )   +
     & 0.000491 * SIN(    224.344795702 *T + 0.878372791 )   +
     & 0.000432 * SIN(  20426.571092422 *T + 6.003829241 )   
       t15= 0.000481 * SIN(   5749.452731634 *T + 4.309591964 )   +
     & 0.000480 * SIN(   5757.317038160 *T + 1.142348571 )   +
     & 0.000485 * SIN(   6702.560493867 *T + 0.210580917 )   +
     & 0.000426 * SIN(   6055.549660552 *T + 4.274476529 )   +
     & 0.000480 * SIN(   5959.570433334 *T + 5.031351030 )   +
     & 0.000466 * SIN(  12562.628581634 *T + 4.959581597 )   +
     & 0.000520 * SIN(  39302.096962196 *T + 4.788002889 )   +
     & 0.000458 * SIN(  12132.439962106 *T + 1.880103788 )   +
     & 0.000470 * SIN(  12029.347187887 *T + 1.405611197 )   +
     & 0.000416 * SIN(  -7477.522860216 *T + 1.082356330 )   +
     & 0.000449 * SIN(  11609.862544012 *T + 4.179989585 )   +
     & 0.000465 * SIN(  17253.041107690 *T + 0.353496295 )   +
     & 0.000362 * SIN(  -4535.059436924 *T + 1.583849576 )   +
     & 0.000383 * SIN(  21954.157609398 *T + 3.747376371 )   +
     & 0.000389 * SIN(     17.252277143 *T + 1.395753179 )   +
     & 0.000331 * SIN(  18052.929543158 *T + 0.566790582 )   +
     & 0.000430 * SIN(  13517.870106233 *T + 0.685827538 )   +
     & 0.000368 * SIN(  -5756.908003246 *T + 0.731374317 )   +
     & 0.000330 * SIN(  10557.594160824 *T + 3.710043680 )   +
     & 0.000332 * SIN(  20199.094959633 *T + 1.652901407 )   
      t16= 0.000384 * SIN(  11933.367960670 *T + 5.827781531 )   +
     & 0.000387 * SIN(  10454.501386605 *T + 2.541182564 )   +
     & 0.000325 * SIN(  15671.081759407 *T + 2.178850542 )   +
     & 0.000318 * SIN(    138.517496871 *T + 2.253253037 )   +
     & 0.000305 * SIN(   9388.005909415 *T + 0.578340206 )   +
     & 0.000352 * SIN(   5749.861766548 *T + 3.000297967 )   +
     & 0.000311 * SIN(   6915.859589305 *T + 1.693574249 )   +
     & 0.000297 * SIN(  24072.921469776 *T + 1.997249392 )   +
     & 0.000363 * SIN(   -640.877607382 *T + 5.071820966 )   +
     & 0.000323 * SIN(  12592.450019783 *T + 1.072262823 )   +
     & 0.000341 * SIN(  12146.667056108 *T + 4.700657997 )   +
     & 0.000290 * SIN(   9779.108676125 *T + 1.812320441 )   +
     & 0.000342 * SIN(   6132.028180148 *T + 4.322238614 )   +
     & 0.000329 * SIN(   6268.848755990 *T + 3.033827743 )   +
     & 0.000374 * SIN(  17996.031168222 *T + 3.388716544 )   +
     & 0.000285 * SIN(   -533.214083444 *T + 4.687313233 )   +
     & 0.000338 * SIN(   6065.844601290 *T + 0.877776108 )   +
     & 0.000276 * SIN(     24.298513841 *T + 0.770299429 )   +
     & 0.000336 * SIN(  -2388.894020449 *T + 5.353796034 )   +
     & 0.000290 * SIN(   3097.883822726 *T + 4.075291557 )   
      t17= 0.000318 * SIN(    709.933048357 *T + 5.941207518 )   +
     & 0.000271 * SIN(  13095.842665077 *T + 3.208912203 )   +
     & 0.000331 * SIN(   6073.708907816 *T + 4.007881169 )   +
     & 0.000292 * SIN(    742.990060533 *T + 2.714333592 )   +
     & 0.000362 * SIN(  29088.811415985 *T + 3.215977013 )   +
     & 0.000280 * SIN(  12359.966151546 *T + 0.710872502 )   +
     & 0.000267 * SIN(  10440.274292604 *T + 4.730108488 )   +
     & 0.000262 * SIN(    838.969287750 *T + 1.327720272 )   +
     & 0.000250 * SIN(  16496.361396202 *T + 0.898769761 )   +
     & 0.000325 * SIN(  20597.243963041 *T + 0.180044365 )   +
     & 0.000268 * SIN(   6148.010769956 *T + 5.152666276 )   +
     & 0.000284 * SIN(   5636.065016677 *T + 5.655385808 )   +
     & 0.000301 * SIN(   6080.822454817 *T + 2.135396205 )   +
     & 0.000294 * SIN(   -377.373607916 *T + 3.708784168 )   +
     & 0.000236 * SIN(   2118.763860378 *T + 1.733578756 )   +
     & 0.000234 * SIN(   5867.523359379 *T + 5.575209112 )   +
     & 0.000268 * SIN(-226858.238553767 *T + 0.069432392 )   +
     & 0.000265 * SIN( 167283.761587465 *T + 4.369302826 )   +
     & 0.000280 * SIN(  28237.233459389 *T + 5.304829118 )   +
     & 0.000292 * SIN(  12345.739057544 *T + 4.096094132 )   
      t18= 0.000223 * SIN(  19800.945956225 *T + 3.069327406 )   +
     & 0.000301 * SIN(  43232.306658416 *T + 6.205311188 )   +
     & 0.000264 * SIN(  18875.525869774 *T + 1.417263408 )   +
     & 0.000304 * SIN(  -1823.175188677 *T + 3.409035232 )   +
     & 0.000301 * SIN(    109.945688789 *T + 0.510922054 )   +
     & 0.000260 * SIN(    813.550283960 *T + 2.389438934 )   +
     & 0.000299 * SIN( 316428.228673312 *T + 5.384595078 )   +
     & 0.000211 * SIN(   5756.566278634 *T + 3.789392838 )   +
     & 0.000209 * SIN(   5750.203491159 *T + 1.661943545 )   +
     & 0.000240 * SIN(  12489.885628707 *T + 5.684549045 )   +
     & 0.000216 * SIN(   6303.851245484 *T + 3.862942261 )   +
     & 0.000203 * SIN(   1581.959348283 *T + 5.549853589 )   +
     & 0.000200 * SIN(   5642.198242609 *T + 1.016115785 )   +
     & 0.000197 * SIN(    -70.849445304 *T + 4.690702525 )   +
     & 0.000227 * SIN(   6287.008003254 *T + 2.911891613 )   +
     & 0.000197 * SIN(    533.623118358 *T + 1.048982898 )   +
     & 0.000205 * SIN(  -6279.485421340 *T + 1.829362730 )   +
     & 0.000209 * SIN( -10988.808157535 *T + 2.636140084 )   +
     & 0.000208 * SIN(   -227.526189440 *T + 4.127883842 )   +
     & 0.000191 * SIN(    415.552490612 *T + 4.401165650 )   
      t19= 0.000190 * SIN(  29296.615389579 *T + 4.175658539 )   +
     & 0.000264 * SIN(  66567.485864652 *T + 4.601102551 )   +
     & 0.000256 * SIN(  -3646.350377354 *T + 0.506364778 )   +
     & 0.000188 * SIN(  13119.721102825 *T + 2.032195842 )   +
     & 0.000185 * SIN(   -209.366942175 *T + 4.694756586 )   +
     & 0.000198 * SIN(  25934.124331089 *T + 3.832703118 )   +
     & 0.000195 * SIN(   4061.219215394 *T + 3.308463427 )   +
     & 0.000234 * SIN(   5113.487598583 *T + 1.716090661 )   +
     & 0.000188 * SIN(   1478.866574064 *T + 5.686865780 )   +
     & 0.000222 * SIN(  11823.161639450 *T + 1.942386641 )   +
     & 0.000181 * SIN(  10770.893256262 *T + 1.999482059 )   +
     & 0.000171 * SIN(   6546.159773364 *T + 1.182807992 )   +
     & 0.000206 * SIN(     70.328180442 *T + 5.934076062 )   +
     & 0.000169 * SIN(  20995.392966449 *T + 2.169080622 )   +
     & 0.000191 * SIN(  10660.686935042 *T + 5.405515999 )   +
     & 0.000228 * SIN(  33019.021112205 *T + 4.656985514 )   +
     & 0.000184 * SIN(  -4933.208440333 *T + 3.327476868 )   +
     & 0.000220 * SIN(   -135.625325010 *T + 1.765430262 )   +
     & 0.000166 * SIN(  23141.558382925 *T + 3.454132746 )   +
     & 0.000191 * SIN(   6144.558353121 *T + 5.020393445 )   
      t20= 0.000180 * SIN(   6084.003848555 *T + 0.602182191 )   +
     & 0.000163 * SIN(  17782.732072784 *T + 4.960593133 )   +
     & 0.000225 * SIN(  16460.333529525 *T + 2.596451817 )   +
     & 0.000222 * SIN(   5905.702242076 *T + 3.731990323 )   +
     & 0.000204 * SIN(    227.476132789 *T + 5.636192701 )   +
     & 0.000159 * SIN(  16737.577236597 *T + 3.600691544 )   +
     & 0.000200 * SIN(   6805.653268085 *T + 0.868220961 )   +
     & 0.000187 * SIN(  11919.140866668 *T + 2.629456641 )   +
     & 0.000161 * SIN(    127.471796607 *T + 2.862574720 )   +
     & 0.000205 * SIN(   6286.666278643 *T + 1.742882331 )   +
     & 0.000189 * SIN(    153.778810485 *T + 4.812372643 )   +
     & 0.000168 * SIN(  16723.350142595 *T + 0.027860588 )   +
     & 0.000149 * SIN(  11720.068865232 *T + 0.659721876 )   +
     & 0.000189 * SIN(   5237.921013804 *T + 5.245313000 )   +
     & 0.000143 * SIN(   6709.674040867 *T + 4.317625647 )   +
     & 0.000146 * SIN(   4487.817406270 *T + 4.815297007 )   +
     & 0.000144 * SIN(   -664.756045130 *T + 5.381366880 )   +
     & 0.000175 * SIN(   5127.714692584 *T + 4.728443327 )   +
     & 0.000162 * SIN(   6254.626662524 *T + 1.435132069 )   +
     & 0.000187 * SIN(  47162.516354635 *T + 1.354371923 )   
      t21= 0.000146 * SIN(  11080.171578918 *T + 3.369695406 )   +
     & 0.000180 * SIN(   -348.924420448 *T + 2.490902145 )   +
     & 0.000148 * SIN(    151.047669843 *T + 3.799109588 )   +
     & 0.000157 * SIN(   6197.248551160 *T + 1.284375887 )   +
     & 0.000167 * SIN(    146.594251718 *T + 0.759969109 )   +
     & 0.000133 * SIN(  -5331.357443741 *T + 5.409701889 )   +
     & 0.000154 * SIN(     95.979227218 *T + 3.366890614 )   +
     & 0.000148 * SIN(  -6418.140930027 *T + 3.384104996 )   +
     & 0.000128 * SIN(  -6525.804453965 *T + 3.803419985 )   +
     & 0.000130 * SIN(  11293.470674356 *T + 0.939039445 )   +
     & 0.000152 * SIN(  -5729.506447149 *T + 0.734117523 )   +
     & 0.000138 * SIN(    210.117701700 *T + 2.564216078 )   +
     & 0.000123 * SIN(   6066.595360816 *T + 4.517099537 )   +
     & 0.000140 * SIN(  18451.078546566 *T + 0.642049130 )   +
     & 0.000126 * SIN(  11300.584221356 *T + 3.485280663 )   +
     & 0.000119 * SIN(  10027.903195729 *T + 3.217431161 )   +
     & 0.000151 * SIN(   4274.518310832 *T + 4.404359108 )   +
     & 0.000117 * SIN(   6072.958148291 *T + 0.366324650 )   +
     & 0.000165 * SIN(  -7668.637425143 *T + 4.298212528 )   +
     & 0.000117 * SIN(  -6245.048177356 *T + 5.379518958 )   
      t22= 0.000130 * SIN(  -5888.449964932 *T + 4.527681115 )   +
     & 0.000121 * SIN(   -543.918059096 *T + 6.109429504 )   +
     & 0.000162 * SIN(   9683.594581116 *T + 5.720092446 )   +
     & 0.000141 * SIN(   6219.339951688 *T + 0.679068671 )   +
     & 0.000118 * SIN(  22743.409379516 *T + 4.881123092 )   +
     & 0.000129 * SIN(   1692.165669502 *T + 0.351407289 )   +
     & 0.000126 * SIN(   5657.405657679 *T + 5.146592349 )   +
     & 0.000114 * SIN(    728.762966531 *T + 0.520791814 )   +
     & 0.000120 * SIN(     52.596639600 *T + 0.948516300 )   +
     & 0.000115 * SIN(     65.220371012 *T + 3.504914846 )   +
     & 0.000126 * SIN(   5881.403728234 *T + 5.577502482 )   +
     & 0.000158 * SIN( 163096.180360983 *T + 2.957128968 )   +
     & 0.000134 * SIN(  12341.806904281 *T + 2.598576764 )   +
     & 0.000151 * SIN(  16627.370915377 *T + 3.985702050 )   +
     & 0.000109 * SIN(   1368.660252845 *T + 0.014730471 )   +
     & 0.000131 * SIN(   6211.263196841 *T + 0.085077024 )   +
     & 0.000146 * SIN(   5792.741760812 *T + 0.708426604 )   +
     & 0.000146 * SIN(    -77.750543984 *T + 3.121576600 )   +
     & 0.000107 * SIN(   5341.013788022 *T + 0.288231904 )   +
     & 0.000138 * SIN(   6281.591377283 *T + 2.797450317 )   
      t23= 0.000113 * SIN(  -6277.552925684 *T + 2.788904128 )   +
     & 0.000115 * SIN(   -525.758811831 *T + 5.895222200 )   +
     & 0.000138 * SIN(   6016.468808270 *T + 6.096188999 )   +
     & 0.000139 * SIN(  23539.707386333 *T + 2.028195445 )   +
     & 0.000146 * SIN(  -4176.041342449 *T + 4.660008502 )   +
     & 0.000107 * SIN(  16062.184526117 *T + 4.066520001 )   +
     & 0.000142 * SIN(  83783.548222473 *T + 2.936315115 )   +
     & 0.000128 * SIN(   9380.959672717 *T + 3.223844306 )   +
     & 0.000135 * SIN(   6205.325306007 *T + 1.638054048 )   +
     & 0.000101 * SIN(   2699.734819318 *T + 5.481603249 )   +
     & 0.000104 * SIN(   -568.821874027 *T + 2.205734493 )   +
     & 0.000103 * SIN(   6321.103522627 *T + 2.440421099 )   +
     & 0.000119 * SIN(   6321.208885629 *T + 2.547496264 )   +
     & 0.000138 * SIN(   1975.492545856 *T + 2.314608466 )   +
     & 0.000121 * SIN(    137.033024162 *T + 4.539108237 )   +
     & 0.000123 * SIN(  19402.796952817 *T + 4.538074405 )   +
     & 0.000119 * SIN(  22805.735565994 *T + 2.869040566 )   +
     & 0.000133 * SIN(  64471.991241142 *T + 6.056405489 )   +
     & 0.000129 * SIN(    -85.827298831 *T + 2.540635083 )   +
     & 0.000131 * SIN(  13613.804277336 *T + 4.005732868 )   
      t24= 0.000104 * SIN(   9814.604100291 *T + 1.959967212 )   +
     & 0.000112 * SIN(  16097.679950283 *T + 3.589026260 )   +
     & 0.000123 * SIN(   2107.034507542 *T + 1.728627253 )   +
     & 0.000121 * SIN(  36949.230808424 *T + 6.072332087 )   +
     & 0.000108 * SIN( -12539.853380183 *T + 3.716133846 )   +
     & 0.000113 * SIN(  -7875.671863624 *T + 2.725771122 )   +
     & 0.000109 * SIN(   4171.425536614 *T + 4.033338079 )   +
     & 0.000101 * SIN(   6247.911759770 *T + 3.441347021 )   +
     & 0.000113 * SIN(   7330.728427345 *T + 0.656372122 )   +
     & 0.000113 * SIN(  51092.726050855 *T + 2.791483066 )   +
     & 0.000106 * SIN(   5621.842923210 *T + 1.815323326 )   +
     & 0.000101 * SIN(    111.430161497 *T + 5.711033677 )   +
     & 0.000103 * SIN(    909.818733055 *T + 2.812745443 )   +
     & 0.000101 * SIN(   1790.642637886 *T + 1.965746028 )   +
c T**1
     & T* 102.156724 * SIN(   6283.075849991 *T + 4.249032005 ) +
     & T* 1.706807 * SIN(  12566.151699983 *T + 4.205904248 )   +
     & T* 0.269668 * SIN(    213.299095438 *T + 3.400290479 )   +
     & T* 0.265919 * SIN(    529.690965095 *T + 5.836047367 )   +
     & T* 0.210568 * SIN(     -3.523118349 *T + 6.262738348 )   +
     & T* 0.077996 * SIN(   5223.693919802 *T + 4.670344204 )   
      t25= T* 0.054764 * SIN(   1577.343542448 *T + 4.534800170 )   +
     & T* 0.059146 * SIN(     26.298319800 *T + 1.083044735 )   +
     & T* 0.034420 * SIN(   -398.149003408 *T + 5.980077351 )   +
     & T* 0.032088 * SIN(  18849.227549974 *T + 4.162913471 )   +
     & T* 0.033595 * SIN(   5507.553238667 *T + 5.980162321 )   +
     & T* 0.029198 * SIN(   5856.477659115 *T + 0.623811863 )   +
     & T* 0.027764 * SIN(    155.420399434 *T + 3.745318113 )   +
     & T* 0.025190 * SIN(   5746.271337896 *T + 2.980330535 )   +
     & T* 0.022997 * SIN(   -796.298006816 *T + 1.174411803 )   +
     & T* 0.024976 * SIN(   5760.498431898 *T + 2.467913690 )   +
     & T* 0.021774 * SIN(    206.185548437 *T + 3.854787540 )   +
     & T* 0.017925 * SIN(   -775.522611324 *T + 1.092065955 )   +
     & T* 0.013794 * SIN(    426.598190876 *T + 2.699831988 )   +
     & T* 0.013276 * SIN(   6062.663207553 *T + 5.845801920 )   +
     & T* 0.011774 * SIN(  12036.460734888 *T + 2.292832062 )   +
     & T* 0.012869 * SIN(   6076.890301554 *T + 5.333425680 )   +
     & T* 0.012152 * SIN(   1059.381930189 *T + 6.222874454 )   +
     & T* 0.011081 * SIN(     -7.113547001 *T + 5.154724984 )   +
     & T* 0.010143 * SIN(   4694.002954708 *T + 4.044013795 )   +
     & T* 0.009357 * SIN(   5486.777843175 *T + 3.416081409 )   
      t26= T* 0.010084 * SIN(    522.577418094 *T + 0.749320262 )   +
     & T* 0.008587 * SIN(  10977.078804699 *T + 2.777152598 )   +
     & T* 0.008628 * SIN(   6275.962302991 *T + 4.562060226 )   +
     & T* 0.008158 * SIN(   -220.412642439 *T + 5.806891533 )   +
     & T* 0.007746 * SIN(   2544.314419883 *T + 1.603197066 )   +
     & T* 0.007670 * SIN(   2146.165416475 *T + 3.000200440 )   +
     & T* 0.007098 * SIN(     74.781598567 *T + 0.443725817 )   +
     & T* 0.006180 * SIN(   -536.804512095 *T + 1.302642751 )   +
     & T* 0.005818 * SIN(   5088.628839767 *T + 4.827723531 )   +
     & T* 0.004945 * SIN(  -6286.598968340 *T + 0.268305170 )   +
     & T* 0.004774 * SIN(   1349.867409659 *T + 5.808636673 )   +
     & T* 0.004687 * SIN(   -242.728603974 *T + 5.154890570 )   +
     & T* 0.006089 * SIN(   1748.016413067 *T + 4.403765209 )   +
     & T* 0.005975 * SIN(  -1194.447010225 *T + 2.583472591 )   +
     & T* 0.004229 * SIN(    951.718406251 *T + 0.931172179 )   +
     & T* 0.005264 * SIN(    553.569402842 *T + 2.336107252 )   +
     & T* 0.003049 * SIN(   5643.178563677 *T + 1.362634430 )   +
     & T* 0.002974 * SIN(   6812.766815086 *T + 1.583012668 )   +
     & T* 0.003403 * SIN(  -2352.866153772 *T + 2.552189886 )   +
     & T* 0.003030 * SIN(    419.484643875 *T + 5.286473844 )   
      t27= T* 0.003210 * SIN(     -7.046236698 *T + 1.863796539 )   +
     & T* 0.003058 * SIN(   9437.762934887 *T + 4.226420633 )   +
     & T* 0.002589 * SIN(  12352.852604545 *T + 1.991935820 )   +
     & T* 0.002927 * SIN(   5216.580372801 *T + 2.319951253 )   +
     & T* 0.002425 * SIN(   5230.807466803 *T + 3.084752833 )   +
     & T* 0.002656 * SIN(   3154.687084896 *T + 2.487447866 )   +
     & T* 0.002445 * SIN(  10447.387839604 *T + 2.347139160 )   +
     & T* 0.002990 * SIN(   4690.479836359 *T + 6.235872050 )   +
     & T* 0.002890 * SIN(   5863.591206116 *T + 0.095197563 )   +
     & T* 0.002498 * SIN(   6438.496249426 *T + 2.994779800 )   +
     & T* 0.001889 * SIN(   8031.092263058 *T + 3.569003717 )   +
     & T* 0.002567 * SIN(    801.820931124 *T + 3.425611498 )   +
     & T* 0.001803 * SIN( -71430.695617928 *T + 2.192295512 )   +
     & T* 0.001782 * SIN(      3.932153263 *T + 5.180433689 )   +
     & T* 0.001694 * SIN(  -4705.732307544 *T + 4.641779174 )   +
     & T* 0.001704 * SIN(  -1592.596013633 *T + 3.997097652 )   +
     & T* 0.001735 * SIN(   5849.364112115 *T + 0.417558428 )   +
     & T* 0.001643 * SIN(   8429.241266467 *T + 2.180619584 )   +
     & T* 0.001680 * SIN(     38.133035638 *T + 4.164529426 )   +
     & T* 0.002045 * SIN(   7084.896781115 *T + 0.526323854 )   
      t28= T* 0.001458 * SIN(   4292.330832950 *T + 1.356098141 )   +
     & T* 0.001437 * SIN(     20.355319399 *T + 3.895439360 )   +
     & T* 0.001738 * SIN(   6279.552731642 *T + 0.087484036 )   +
     & T* 0.001367 * SIN(  14143.495242431 *T + 3.987576591 )   +
     & T* 0.001344 * SIN(   7234.794256242 *T + 0.090454338 )   +
     & T* 0.001438 * SIN(  11499.656222793 *T + 0.974387904 )   +
     & T* 0.001257 * SIN(   6836.645252834 *T + 1.509069366 )   +
     & T* 0.001358 * SIN(  11513.883316794 *T + 0.495572260 )   +
     & T* 0.001628 * SIN(   7632.943259650 *T + 4.968445721 )   +
     & T* 0.001169 * SIN(    103.092774219 *T + 2.838496795 )   +
     & T* 0.001162 * SIN(   4164.311989613 *T + 3.408387778 )   +
     & T* 0.001092 * SIN(   6069.776754553 *T + 3.617942651 )   +
     & T* 0.001008 * SIN(  17789.845619785 *T + 0.286350174 )   +
     & T* 0.001008 * SIN(    639.897286314 *T + 1.610762073 )   +
     & T* 0.000918 * SIN(  10213.285546211 *T + 5.532798067 )   +
     & T* 0.001011 * SIN(  -6256.777530192 *T + 0.661826484 )   +
     & T* 0.000753 * SIN(  16730.463689596 *T + 3.905030235 )   +
     & T* 0.000737 * SIN(  11926.254413669 *T + 4.641956361 )   +
     & T* 0.000694 * SIN(   3340.612426700 *T + 2.111120332 )   +
     & T* 0.000701 * SIN(   3894.181829542 *T + 2.760823491 )   
      t29= T* 0.000689 * SIN(   -135.065080035 *T + 4.768800780 )   +
     & T* 0.000700 * SIN(  13367.972631107 *T + 5.760439898 )   +
     & T* 0.000664 * SIN(   6040.347246017 *T + 1.051215840 )   +
     & T* 0.000654 * SIN(   5650.292110678 *T + 4.911332503 )   +
     & T* 0.000788 * SIN(   6681.224853400 *T + 4.699648011 )   +
     & T* 0.000628 * SIN(   5333.900241022 *T + 5.024608847 )   +
     & T* 0.000755 * SIN(   -110.206321219 *T + 4.370971253 )   +
     & T* 0.000628 * SIN(   6290.189396992 *T + 3.660478857 )   +
     & T* 0.000635 * SIN(  25132.303399966 *T + 4.121051532 )   +
     & T* 0.000534 * SIN(   5966.683980335 *T + 1.173284524 )   +
     & T* 0.000543 * SIN(   -433.711737877 *T + 0.345585464 )   +
     & T* 0.000517 * SIN(  -1990.745017041 *T + 5.414571768 )   +
     & T* 0.000504 * SIN(   5767.611978898 *T + 2.328281115 )   +
     & T* 0.000485 * SIN(   5753.384884897 *T + 1.685874771 )   +
     & T* 0.000463 * SIN(   7860.419392439 *T + 5.297703006 )   +
     & T* 0.000604 * SIN(    515.463871093 *T + 0.591998446 )   +
     & T* 0.000443 * SIN(  12168.002696575 *T + 4.830881244 )   +
     & T* 0.000570 * SIN(    199.072001436 *T + 3.899190272 )   +
     & T* 0.000465 * SIN(  10969.965257698 *T + 0.476681802 )   +
     & T* 0.000424 * SIN(  -7079.373856808 *T + 1.112242763 )   
      t30= T* 0.000427 * SIN(    735.876513532 *T + 1.994214480 )   +
     & T* 0.000478 * SIN(  -6127.655450557 *T + 3.778025483 )   +
     & T* 0.000414 * SIN(  10973.555686350 *T + 5.441088327 )   +
     & T* 0.000512 * SIN(   1589.072895284 *T + 0.107123853 )   +
     & T* 0.000378 * SIN(  10984.192351700 *T + 0.915087231 )   +
     & T* 0.000402 * SIN(  11371.704689758 *T + 4.107281715 )   +
     & T* 0.000453 * SIN(   9917.696874510 *T + 1.917490952 )   +
     & T* 0.000395 * SIN(    149.563197135 *T + 2.763124165 )   +
     & T* 0.000371 * SIN(   5739.157790895 *T + 3.112111866 )   +
     & T* 0.000350 * SIN(  11790.629088659 *T + 0.440639857 )   +
     & T* 0.000356 * SIN(   6133.512652857 *T + 5.444568842 )   +
     & T* 0.000344 * SIN(    412.371096874 *T + 5.676832684 )   +
     & T* 0.000383 * SIN(    955.599741609 *T + 5.559734846 )   +
     & T* 0.000333 * SIN(   6496.374945429 *T + 0.261537984 )   +
     & T* 0.000340 * SIN(   6055.549660552 *T + 5.975534987 )   +
     & T* 0.000334 * SIN(   1066.495477190 *T + 2.335063907 )   +
     & T* 0.000399 * SIN(  11506.769769794 *T + 5.321230910 )   +
     & T* 0.000314 * SIN(  18319.536584880 *T + 2.313312404 )   +
     & T* 0.000424 * SIN(   1052.268383188 *T + 1.211961766 )   +
     & T* 0.000307 * SIN(     63.735898303 *T + 3.169551388 )   
      t31= T* 0.000329 * SIN(     29.821438149 *T + 6.106912080 )   +
     & T* 0.000357 * SIN(   6309.374169791 *T + 4.223760346 )   +
     & T* 0.000312 * SIN(  -3738.761430108 *T + 2.180556645 )   +
     & T* 0.000301 * SIN(    309.278322656 *T + 1.499984572 )   +
     & T* 0.000268 * SIN(  12043.574281889 *T + 2.447520648 )   +
     & T* 0.000257 * SIN(  12491.370101415 *T + 3.662331761 )   +
     & T* 0.000290 * SIN(    625.670192312 *T + 1.272834584 )   +
     & T* 0.000256 * SIN(   5429.879468239 *T + 1.913426912 )   +
     & T* 0.000339 * SIN(   3496.032826134 *T + 4.165930011 )   +
     & T* 0.000283 * SIN(   3930.209696220 *T + 4.325565754 )   +
     & T* 0.000241 * SIN(  12528.018664345 *T + 3.832324536 )   +
     & T* 0.000304 * SIN(   4686.889407707 *T + 1.612348468 )   +
     & T* 0.000259 * SIN(  16200.772724501 *T + 3.470173146 )   +
     & T* 0.000238 * SIN(  12139.553509107 *T + 1.147977842 )   +
     & T* 0.000236 * SIN(   6172.869528772 *T + 3.776271728 )   +
     & T* 0.000296 * SIN(  -7058.598461315 *T + 0.460368852 )   +
     & T* 0.000306 * SIN(  10575.406682942 *T + 0.554749016 )   +
     & T* 0.000251 * SIN(  17298.182327326 *T + 0.834332510 )   +
     & T* 0.000290 * SIN(   4732.030627343 *T + 4.759564091 )   +
     & T* 0.000261 * SIN(   5884.926846583 *T + 0.298259862 )   
      t32= T* 0.000249 * SIN(   5547.199336460 *T + 3.749366406 )   +
     & T* 0.000213 * SIN(  11712.955318231 *T + 5.415666119 )   +
     & T* 0.000223 * SIN(   4701.116501708 *T + 2.703203558 )   +
     & T* 0.000268 * SIN(   -640.877607382 *T + 0.283670793 )   +
     & T* 0.000209 * SIN(   5636.065016677 *T + 1.238477199 )   +
     & T* 0.000193 * SIN(  10177.257679534 *T + 1.943251340 )   +
     & T* 0.000182 * SIN(   6283.143160294 *T + 2.456157599 )   +
     & T* 0.000184 * SIN(   -227.526189440 *T + 5.888038582 )   +
     & T* 0.000182 * SIN(  -6283.008539689 *T + 0.241332086 )   +
     & T* 0.000228 * SIN(  -6284.056171060 *T + 2.657323816 )   +
     & T* 0.000166 * SIN(   7238.675591600 *T + 5.930629110 )   +
     & T* 0.000167 * SIN(   3097.883822726 *T + 5.570955333 )   +
     & T* 0.000159 * SIN(   -323.505416657 *T + 5.786670700 )   +
     & T* 0.000154 * SIN(  -4136.910433516 *T + 1.517805532 )   +
     & T* 0.000176 * SIN(  12029.347187887 *T + 3.139266834 )   +
     & T* 0.000167 * SIN(  12132.439962106 *T + 3.556352289 )   +
     & T* 0.000153 * SIN(    202.253395174 *T + 1.463313961 )   +
     & T* 0.000157 * SIN(  17267.268201691 *T + 1.586837396 )   +
     & T* 0.000142 * SIN(  83996.847317911 *T + 0.022670115 )   +
     & T* 0.000152 * SIN(  17260.154654690 *T + 0.708528947 )   
      t33= T* 0.000144 * SIN(   6084.003848555 *T + 5.187075177 )   + 
     & T* 0.000135 * SIN(   5756.566278634 *T + 1.993229262 )   + 
     & T* 0.000134 * SIN(   5750.203491159 *T + 3.457197134 )   +
     & T* 0.000144 * SIN(   5326.786694021 *T + 6.066193291 )   +
     & T* 0.000160 * SIN(  11015.106477335 *T + 1.710431974 )   +
     & T* 0.000133 * SIN(   3634.621024518 *T + 2.836451652 )   +
     & T* 0.000134 * SIN(  18073.704938650 *T + 5.453106665 )   +
     & T* 0.000134 * SIN(   1162.474704408 *T + 5.326898811 )   +
     & T* 0.000128 * SIN(   5642.198242609 *T + 2.511652591 )   +
     & T* 0.000160 * SIN(    632.783739313 *T + 5.628785365 )   +
     & T* 0.000132 * SIN(  13916.019109642 *T + 0.819294053 )   +
     & T* 0.000122 * SIN(  14314.168113050 *T + 5.677408071 )   +
     & T* 0.000125 * SIN(  12359.966151546 *T + 5.251984735 )   +
     & T* 0.000121 * SIN(   5749.452731634 *T + 2.210924603 )   +
     & T* 0.000136 * SIN(   -245.831646229 *T + 1.646502367 )   +
     & T* 0.000120 * SIN(   5757.317038160 *T + 3.240883049 )   +
     & T* 0.000134 * SIN(  12146.667056108 *T + 3.059480037 )   +
     & T* 0.000137 * SIN(   6206.809778716 *T + 1.867105418 )   +
     & T* 0.000141 * SIN(  17253.041107690 *T + 2.069217456 )   +
     & T* 0.000129 * SIN(  -7477.522860216 *T + 2.781469314 )   
      t34= T* 0.000116 * SIN(   5540.085789459 *T + 4.281176991 )   +
     & T* 0.000116 * SIN(   9779.108676125 *T + 3.320925381 )   +
     & T* 0.000129 * SIN(   5237.921013804 *T + 3.497704076 )   +
     & T* 0.000113 * SIN(   5959.570433334 *T + 0.983210840 )   +
     & T* 0.000122 * SIN(   6282.095528923 *T + 2.674938860 )   +
     & T* 0.000140 * SIN(    -11.045700264 *T + 4.957936982 )   +
     & T* 0.000108 * SIN(  23543.230504682 *T + 1.390113589 )   +
     & T* 0.000106 * SIN( -12569.674818332 *T + 0.429631317 )   +
     & T* 0.000110 * SIN(   -266.607041722 *T + 5.501340197 )   +
     & T* 0.000115 * SIN(  12559.038152982 *T + 4.691456618 )   +
     & T* 0.000134 * SIN(  -2388.894020449 *T + 0.577313584 )   +
     & T* 0.000109 * SIN(  10440.274292604 *T + 6.218148717 )   +
     & T* 0.000102 * SIN(   -543.918059096 *T + 1.477842615 )   +
     & T* 0.000108 * SIN(  21228.392023546 *T + 2.237753948 )   +
     & T* 0.000101 * SIN(  -4535.059436924 *T + 3.100492232 )   +
     & T* 0.000103 * SIN(     76.266071276 *T + 5.594294322 )   +
     & T* 0.000104 * SIN(    949.175608970 *T + 5.674287810 )   +
     & T* 0.000101 * SIN(  13517.870106233 *T + 2.196632348 )   +
     & T* 0.000100 * SIN(  11933.367960670 *T + 4.056084160 )   +
c T**2
     & TT* 4.322990 * SIN(   6283.075849991 *T + 2.642893748 )   
      t35= TT* 0.406495 * SIN(      0.000000000 *T + 4.712388980 )   +
     & TT* 0.122605 * SIN(  12566.151699983 *T + 2.438140634 )   +
     & TT* 0.019476 * SIN(    213.299095438 *T + 1.642186981 )   +
     & TT* 0.016916 * SIN(    529.690965095 *T + 4.510959344 )   +
     & TT* 0.013374 * SIN(     -3.523118349 *T + 1.502210314 )   +
     & TT* 0.008042 * SIN(     26.298319800 *T + 0.478549024 )   +
     & TT* 0.007824 * SIN(    155.420399434 *T + 5.254710405 )   +
     & TT* 0.004894 * SIN(   5746.271337896 *T + 4.683210850 )   +
     & TT* 0.004875 * SIN(   5760.498431898 *T + 0.759507698 )   +
     & TT* 0.004416 * SIN(   5223.693919802 *T + 6.028853166 )   +
     & TT* 0.004088 * SIN(     -7.113547001 *T + 0.060926389 )   +
     & TT* 0.004433 * SIN(  77713.771467920 *T + 3.627734103 )   +
     & TT* 0.003277 * SIN(  18849.227549974 *T + 2.327912542 )   +
     & TT* 0.002703 * SIN(   6062.663207553 *T + 1.271941729 )   +
     & TT* 0.003435 * SIN(   -775.522611324 *T + 0.747446224 )   +
     & TT* 0.002618 * SIN(   6076.890301554 *T + 3.633715689 )   +
     & TT* 0.003146 * SIN(    206.185548437 *T + 5.647874613 )   +
     & TT* 0.002544 * SIN(   1577.343542448 *T + 6.232904270 )   +
     & TT* 0.002218 * SIN(   -220.412642439 *T + 1.309509946 )   +
     & TT* 0.002197 * SIN(   5856.477659115 *T + 2.407212349 )   
      t36= TT* 0.002897 * SIN(   5753.384884897 *T + 5.863842246 )   +
     & TT* 0.001766 * SIN(    426.598190876 *T + 0.754113147 )   +
     & TT* 0.001738 * SIN(   -796.298006816 *T + 2.714942671 )   +
     & TT* 0.001695 * SIN(    522.577418094 *T + 2.629369842 )   +
     & TT* 0.001584 * SIN(   5507.553238667 *T + 1.341138229 )   +
     & TT* 0.001503 * SIN(   -242.728603974 *T + 0.377699736 )   +
     & TT* 0.001552 * SIN(   -536.804512095 *T + 2.904684667 )   +
     & TT* 0.001370 * SIN(   -398.149003408 *T + 1.265599125 )   +
     & TT* 0.001889 * SIN(  -5573.142801634 *T + 4.413514859 )   +
     & TT* 0.001722 * SIN(   6069.776754553 *T + 2.445966339 )   +
     & TT* 0.001124 * SIN(   1059.381930189 *T + 5.041799657 )   +
     & TT* 0.001258 * SIN(    553.569402842 *T + 3.849557278 )   +
     & TT* 0.000831 * SIN(    951.718406251 *T + 2.471094709 )   +
     & TT* 0.000767 * SIN(   4694.002954708 *T + 5.363125422 )   +
     & TT* 0.000756 * SIN(   1349.867409659 *T + 1.046195744 )   +
     & TT* 0.000775 * SIN(    -11.045700264 *T + 0.245548001 )   +
     & TT* 0.000597 * SIN(   2146.165416475 *T + 4.543268798 )   +
     & TT* 0.000568 * SIN(   5216.580372801 *T + 4.178853144 )   +
     & TT* 0.000711 * SIN(   1748.016413067 *T + 5.934271972 )   +
     & TT* 0.000499 * SIN(  12036.460734888 *T + 0.624434410 )   
      t37= TT* 0.000671 * SIN(  -1194.447010225 *T + 4.136047594 )   +
     & TT* 0.000488 * SIN(   5849.364112115 *T + 2.209679987 )   +
     & TT* 0.000621 * SIN(   6438.496249426 *T + 4.518860804 )   +
     & TT* 0.000495 * SIN(  -6286.598968340 *T + 1.868201275 )   +
     & TT* 0.000456 * SIN(   5230.807466803 *T + 1.271231591 )   +
     & TT* 0.000451 * SIN(   5088.628839767 *T + 0.084060889 )   +
     & TT* 0.000435 * SIN(   5643.178563677 *T + 3.324456609 )   +
     & TT* 0.000387 * SIN(  10977.078804699 *T + 4.052488477 )   +
     & TT* 0.000547 * SIN( 161000.685737473 *T + 2.841633844 )   +
     & TT* 0.000522 * SIN(   3154.687084896 *T + 2.171979966 )   +
     & TT* 0.000375 * SIN(   5486.777843175 *T + 4.983027306 )   +
     & TT* 0.000421 * SIN(   5863.591206116 *T + 4.546432249 )   +
     & TT* 0.000439 * SIN(   7084.896781115 *T + 0.522967921 )   +
     & TT* 0.000309 * SIN(   2544.314419883 *T + 3.172606705 )   +
     & TT* 0.000347 * SIN(   4690.479836359 *T + 1.479586566 )   +
     & TT* 0.000317 * SIN(    801.820931124 *T + 3.553088096 )   +
     & TT* 0.000262 * SIN(    419.484643875 *T + 0.606635550 )   +
     & TT* 0.000248 * SIN(   6836.645252834 *T + 3.014082064 )   +
     & TT* 0.000245 * SIN(  -1592.596013633 *T + 5.519526220 )   +
     & TT* 0.000225 * SIN(   4292.330832950 *T + 2.877956536 )   
      t38= TT* 0.000214 * SIN(   7234.794256242 *T + 1.605227587 )   +
     & TT* 0.000205 * SIN(   5767.611978898 *T + 0.625804796 )   +
     & TT* 0.000180 * SIN(  10447.387839604 *T + 3.499954526 )   +
     & TT* 0.000229 * SIN(    199.072001436 *T + 5.632304604 )   +
     & TT* 0.000214 * SIN(    639.897286314 *T + 5.960227667 )   +
     & TT* 0.000175 * SIN(   -433.711737877 *T + 2.162417992 )   +
     & TT* 0.000209 * SIN(    515.463871093 *T + 2.322150893 )   +
     & TT* 0.000173 * SIN(   6040.347246017 *T + 2.556183691 )   +
     & TT* 0.000184 * SIN(   6309.374169791 *T + 4.732296790 )   +
     & TT* 0.000227 * SIN( 149854.400134205 *T + 5.385812217 )   +
     & TT* 0.000154 * SIN(   8031.092263058 *T + 5.120720920 )   +
     & TT* 0.000151 * SIN(   5739.157790895 *T + 4.815000443 )   +
     & TT* 0.000197 * SIN(   7632.943259650 *T + 0.222827271 )   +
     & TT* 0.000197 * SIN(     74.781598567 *T + 3.910456770 )   +
     & TT* 0.000138 * SIN(   6055.549660552 *T + 1.397484253 )   +
     & TT* 0.000149 * SIN(  -6127.655450557 *T + 5.333727496 )   +
     & TT* 0.000137 * SIN(   3894.181829542 *T + 4.281749907 )   +
     & TT* 0.000135 * SIN(   9437.762934887 *T + 5.979971885 )   +
     & TT* 0.000139 * SIN(  -2352.866153772 *T + 4.715630782 )   +
     & TT* 0.000142 * SIN(   6812.766815086 *T + 0.513330157 )   
      t39= TT* 0.000120 * SIN(  -4705.732307544 *T + 0.194160689 )   +
     & TT* 0.000131 * SIN( -71430.695617928 *T + 0.000379226 )   +
     & TT* 0.000124 * SIN(   6279.552731642 *T + 2.122264908 )   +
     & TT* 0.000108 * SIN(  -6256.777530192 *T + 0.883445696 )   +
c T**3
     & TTT* 0.143388 * SIN(   6283.075849991 *T + 1.131453581 )   +
     & TTT* 0.006671 * SIN(  12566.151699983 *T + 0.775148887 )   +
     & TTT* 0.001480 * SIN(    155.420399434 *T + 0.480016880 )   +
     & TTT* 0.000934 * SIN(    213.299095438 *T + 6.144453084 )   +
     & TTT* 0.000795 * SIN(    529.690965095 *T + 2.941595619 )   +
     & TTT* 0.000673 * SIN(   5746.271337896 *T + 0.120415406 )   +
     & TTT* 0.000672 * SIN(   5760.498431898 *T + 5.317009738 )   +
     & TTT* 0.000389 * SIN(   -220.412642439 *T + 3.090323467 )   +
     & TTT* 0.000373 * SIN(   6062.663207553 *T + 3.003551964 )   +
     & TTT* 0.000360 * SIN(   6076.890301554 *T + 1.918913041 )   +
     & TTT* 0.000316 * SIN(    -21.340641002 *T + 5.545798121 )   +
     & TTT* 0.000315 * SIN(   -242.728603974 *T + 1.884932563 )   +
     & TTT* 0.000278 * SIN(    206.185548437 *T + 1.266254859 )   +
     & TTT* 0.000238 * SIN(   -536.804512095 *T + 4.532664830 )   +
     & TTT* 0.000185 * SIN(    522.577418094 *T + 4.578313856 )   +
     & TTT* 0.000245 * SIN(  18849.227549974 *T + 0.587467082 )   
      t40= TTT* 0.000180 * SIN(    426.598190876 *T + 5.151178553 )   +
     & TTT* 0.000200 * SIN(    553.569402842 *T + 5.355983739 )   +
     & TTT* 0.000141 * SIN(   5223.693919802 *T + 1.336556009 )   +
     & TTT* 0.000104 * SIN(   5856.477659115 *T + 4.239842759 )   +
c T**4
     & TTTT* 0.003826 * SIN(   6283.075849991 *T + 5.705257275 )   +
     & TTTT* 0.000303 * SIN(  12566.151699983 *T + 5.407132842 )   +
     & TTTT* 0.000209 * SIN(    155.420399434 *T + 1.989815753 )   
      t41=0.00065*sin(6069.776754*T+4.021194)+
     &    0.00033*sin( 213.299095*T+5.543132)
     &   -0.00196*sin(6208.294251*T+5.696701)
     &   -0.00173*sin(  74.781599*T+2.435900)

      ctatv=t1+t2+t3+t4+t5+t6+t7+t8+t9+t10+t11+t12+t13+t14+t15
     &     +t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28+t29+t30
     &     +t31+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41
      ctatv=ctatv*1.d-6

      tdbsec=tdtsec*86400.d0+ctatv

      RETURN
      END 

