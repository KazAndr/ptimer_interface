c v 1.01
*
* reads new free-form format
*********************************************************************
*********************************************************************
*
* Read PSR parameters in free-form format
* from unit number funit untill option 'TOAS'
*    
       subroutine readpar(funit,psrname,found,norb,neph,f,
     +     p0,p1,p2,px,pepoch,aps,dps,pmra,pmdec,dm,
     +     a1,e,edot,t0,pb,om,omdot,gamma,pbdot,si,am,am2,
     +     a1dot,  
     +     bp,bpp,dr,dth,a0,b0,pbd2,pbd3,xd2,IE,
     *     NITER)
c
c funit    - log.unit number of parameter file
c
c output:
c psrname  - source name
c found   - .true. if record 'TOAS' os found in par.file
c norb    - number of components in system (single PSR norb=0)
c neph    - =0 for J2000 =1 for B1950
c f(1)=f0,f(2)=f1 - 1st deriv ,...
c p0,p1,p2- period and 1st and 2nd derivatives
c px      - distance to PSR, (pck)
c pepoch  - epoch of p0,...,f.., (MJD)
c aps     - R.A., (hh.mmssds)
c dps     - DEC,  (dd.mmssds)
c pmra    - p.m.R.A. (mas/yr)
c pmdec   - p.m.DEC  (mas/yr)
c dm      - DM, (pc*cm^-3)
c a1(i)   - proj.semi-major axis, (sec)
c e(i)    - eccentricity
c edot(i) - e-dot, (1/s)
c t0(i)   - time of periastron pass., (JD)
c pb(i)   - binary period, (sec)
c om(i)   - longitude of periastron, (deg)
c omdot(i)- rate of periastron advance, (deg/yr)
c a1dot(i)- rate of change of a1 (s/s)
c gamma(i)- grav.redshift and time delation, (s)
c pbdot(i)- orbital period change, (10^-12)
c si(i)   - sine of inclination angle
c am      - total system mass, (solar masses)
c am2(i)  - companion's mass, (solar masses)
c a1dot(i)- deriv. of x (10^-12)
c bp      - beta-prime
c bpp     - beta-prime-prime
c dr      - relativ.orbital deformation
c dth(i)  - relativ.orbital deformation (10^-6)
c a0,b0   - aberration parameter
c
c
c NITER  - number of iteration in fit, changed if command NITS found
c
c 
c bp,bpp,dr,a0,b0 - not yet used
c
       
      implicit none
      integer *4 IDEGS
      parameter(idegs=224)
      character chstr*80,valstr*80,psrname*9,chbinm*4
      integer*4 norb,neph,i,k,compMax,funit,ie(idegs),niter

      parameter (compMax=3)
      real*8 p0,p1,p2,px,pepoch,aps,dps,pmra,pmdec,dm(11),f(11)
      real*8 a1(compMax),e(compMax),t0(compMax),pb(compMax),om(compMax)
      real*8 omdot(compMax),gamma(compMax),pbdot(compMax),si(compMax),
     + a1dot(compMax),edot(compMax),am2(compMax),am,bp,bpp,dr,
     + dth(compMax),a0,b0, eps1,eps2,tasc,pi,pbd2(compMax),Xd2(compMax),
     + pbd3(compMax)
      logical found
      
      found=.false.
      PI=3.1415926535897932D0
         
c-- init variables
      chbinm='  '
      norb=0
      neph=0 ! default J2000
      p0=0.0
      p1=0.0
      p2=0.0
      eps1=0.0
      eps2=0.0
       do i=1,11
       f(i)=0.0d0
       dm(i)=0.0       
       enddo
      px=0.0d0
      pepoch=0.0
      aps=0.0
      dps=0.0
      pmra=0.0
      pmdec=0.0
      am=0.0
      bp=0.0
      bpp=0.0
      do i=1, compMax
         a1(i)=0.0
         e(i)=0.0
         t0(i)=0.0
         pb(i)=0.0
         om(i)=0.0
         omdot(i)=0.0
         a1dot(i)=0.0
         pbdot(i)=0.0
         gamma(i)=0.0
         am2(i)=0.0
         si(i)=0.0
      enddo
      IE(1)=1
 100  continue
 222  format(f24.12)
 223  format(f24.12,i1)
      
*(Read Pulsar parameters===========================================)
         do i=1,1000
         read(funit,'(a80)',end=800) chstr
* (convert string to upper case and format it to (f24.12,i1))
	 call upval(chstr,valstr)

          if(chstr(1:3).eq.'TOA') then ! exit when met 'TOAS' string
          found=.true.
          goto 800
         endif
         if(chstr(1:3).eq.'PSR')read(valstr,'(a9)')    psrname(1:9)! (PSRname)
         if(chstr(1:2).eq.'RA') read(valstr,223,err=3) aps,ie(6)   ! (RA)
         if(chstr(1:3).eq.'DEC')read(valstr,223,err=3) dps,ie(7)   ! (DEC)
         if(chstr(1:3).eq.'PMR')read(valstr,223,err=3) pmra,ie(8)  ! (PMRA)
         if(chstr(1:3).eq.'PMD')read(valstr,223,err=3) pmdec,ie(9) ! (PMDEC)
         if(chstr(1:2).eq.'PX') read(valstr,223,err=3) px,ie(10)   ! (PX)
         if(chstr(1:2).eq.'P0') read(valstr,223,err=3) p0,ie(2)    ! (P0/P)
         if(chstr(1:2).eq.'P ') read(valstr,223,err=3) p0,ie(2)    !  =
         if(chstr(1:2).eq.'P1') read(valstr,223,err=3) p1,ie(3)    ! (P1)
         if(chstr(1:2).eq.'P2') read(valstr,223,err=3) p2,ie(4)    ! (P2)
         if(chstr(1:2).eq.'F0') read(valstr,223,err=3) f(1),ie(2)  ! (F0/F)
         if(chstr(1:2).eq.'F ') read(valstr,223,err=3) f(1),ie(2)  !  =
         if(chstr(1:2).eq.'F1') read(valstr,223,err=3) f(2),ie(3)  ! (F1)
         if(chstr(1:2).eq.'F2') read(valstr,223,err=3) f(3),ie(4)  ! (F2)
         if(chstr(1:2).eq.'F3') read(valstr,223,err=3) f(4),ie(5)  ! (F3)
         if(chstr(1:2).eq.'F4') read(valstr,223,err=3) f(5),ie(60) ! (F4)
         if(chstr(1:2).eq.'F5') read(valstr,223,err=3) f(6),ie(61) ! (F5)
         if(chstr(1:2).eq.'F6') read(valstr,223,err=3) f(7),ie(62) ! (F6)
         if(chstr(1:2).eq.'F7') read(valstr,223,err=3) f(8),ie(63) ! (F7)
         if(chstr(1:2).eq.'F8') read(valstr,223,err=3) f(9),ie(64) ! (F8)
         if(chstr(1:2).eq.'F9') read(valstr,222,err=3) f(10)       ! (F9)
         if(chstr(1:3).eq.'F10')read(valstr,222,err=3) f(11)       ! (F10)
         if(chstr(1:3).eq.'PEP') then                              ! (PEPOCH)
         read(valstr,222,err=3) pepoch
         if(pepoch.gt.2400000.5d0) pepoch=pepoch-2400000.5d0
         endif
         if(chstr(1:3).eq.'DM ') read(valstr,223,err=3) dm(1),ie(11) ! (DM)
         if(chstr(1:3).eq.'DM1') read(valstr,223,err=3) dm(2),ie(70) ! (DMd)
         if(chstr(1:3).eq.'DM2') read(valstr,223,err=3) dm(3),ie(71) ! (DMd2)
         if(chstr(1:3).eq.'DM3') read(valstr,223,err=3) dm(4),ie(72) ! (DMd3)
         if(chstr(1:3).eq.'DM4') read(valstr,223,err=3) dm(5),ie(73) ! (DMd4)
         if(chstr(1:3).eq.'DM5') read(valstr,223,err=3) dm(6),ie(74) ! (DMd5)
         if(chstr(1:3).eq.'DM6') read(valstr,223,err=3) dm(7),ie(75) ! (DMd6)
         if(chstr(1:4).eq.'COOR') then                               ! (COORD)
          if(valstr(1:5).eq.'B1950') neph=1
         endif
         if(chstr(1:4).eq.'NITS') read(valstr,'(i9)',err=3) NITER    ! (NITS)
         
*(Binary model in CHBINM)
         if(chstr(1:4).eq.'BINA') read(valstr,'(a4)') chbinm         ! (BINARY)
*(a1,e,t0,pb,om) 1st companion
         if(chstr(2:2).ne.'_'.or.
     &     chstr(3:3).ne.'_'.or.chstr(4:4).ne.'_'.or.
     &     chstr(5:5).ne.'_'.or.chstr(6:6).ne.'_'.or.
     &     chstr(7:7).ne.'_') then
         if(chstr(1:3).eq.'A1 '.or.chstr(1:2).eq.'A ')
     &     read(valstr,223,err=3) a1(1),ie(12)                           ! (A1/A)
         if(chstr(1:3).eq.'E1 '.or.chstr(1:2).eq.'E ')
     &     read(valstr,223,err=3) e(1),ie(14)                            ! (E1/E)
         if(chstr(1:3).eq.'T0 ')  read(valstr,223,err=3) t0(1),ie(15)    ! (T0)
         if(chstr(1:3).eq.'PB ') read(valstr,223,err=3) pb(1),ie(13)     ! (PB)
         if(chstr(1:3).eq.'OM ') read(valstr,223,err=3) om(1),ie(16)     ! (OM)
         endif ! (chstr.ne.'_')
*(a1,e,t0,pb,om) 1st companion (alternative for above readings)
         if(chstr(1:4).eq.'A1_1') read(valstr,223,err=3) a1(1),ie(12)    ! (A1_1)
         if(chstr(1:3).eq.'E_1')  read(valstr,223,err=3) e(1), ie(14)    ! (E1_1)
         if(chstr(1:4).eq.'T0_1') read(valstr,223,err=3) t0(1),ie(15)    ! (T0_1)
         if(chstr(1:4).eq.'PB_1') read(valstr,223,err=3) pb(1),ie(13)    ! (PB_1)
         if(chstr(1:4).eq.'OM_1') read(valstr,223,err=3) om(1),ie(16)    ! (OM_1)
*(a1,e,t0,pb,om) 2nd companion
         if(chstr(1:4).eq.'A1_2') read(valstr,223,err=3) a1(2),ie(25)    ! (A1_2)
         if(chstr(1:3).eq.'E_2')  read(valstr,223,err=3) e(2), ie(27)    ! (E1_2)
         if(chstr(1:4).eq.'T0_2') read(valstr,223,err=3) t0(2),ie(28)    ! (T0_2)
         if(chstr(1:4).eq.'PB_2') read(valstr,223,err=3) pb(2),ie(26)    ! (PB_2)
         if(chstr(1:4).eq.'OM_2') read(valstr,223,err=3) om(2),ie(29)    ! (OM_2)
*(a1,e,t0,pb,om) 3rd companion
         if(chstr(1:4).eq.'A1_3') read(valstr,223,err=3) a1(3),ie(38)    ! (A1_3)
         if(chstr(1:3).eq.'E_3')  read(valstr,223,err=3) e(3), ie(40)    ! (E1_3)
         if(chstr(1:4).eq.'T0_3') read(valstr,223,err=3) t0(3),ie(41)    ! (T0_3)
         if(chstr(1:4).eq.'PB_3') read(valstr,223,err=3) pb(3),ie(39)    ! (PB_3)
         if(chstr(1:4).eq.'OM_3') read(valstr,223,err=3) om(3),ie(42)    ! (OM_3)
* (PN-parameters)
       if(chstr(1:6).eq.'OMDOT ')read(valstr,223,err=3) omdot(1),ie(20)   ! (OMDOT)
       if(chstr(1:7).eq.'OMDOT_1')read(valstr,223,err=3)omdot(1),ie(20)   ! (OMDOT_1)
       if(chstr(1:5).eq.'XDOT ') read(valstr,223,err=3) a1dot(1),ie(21)   ! (XDOT)
       if(chstr(1:6).eq.'XDOT_1')read(valstr,223,err=3) a1dot(1),ie(21)   ! (XDOT_1)
       if(chstr(1:6).eq.'PBDOT ')read(valstr,223,err=3) pbdot(1),ie(22)   ! (PBDOT)
       if(chstr(1:7).eq.'PBDOT_1')read(valstr,223,err=3)pbdot(1),ie(22)   ! (PBDOT_1)       
       if(chstr(1:5).eq.'EDOT ') read(valstr,223,err=3) edot(1), ie(23)   ! (EDOT)
       if(chstr(1:6).eq.'EDOT_1')read(valstr,223,err=3) edot(1), ie(23)   ! (EDOT_1)       
       if(chstr(1:7).eq.'PB2DOT ')read(valstr,223,err=3) pbd2(1),ie(101)  ! (PB2DOT)
       if(chstr(1:8).eq.'PB2DOT_1')read(valstr,223,err=3)pbd2(1),ie(101)  ! (PB2DOT_1)
       if(chstr(1:7).eq.'PB3DOT ')read(valstr,223,err=3) pbd3(1),ie(104)  ! (PB3DOT)
       if(chstr(1:8).eq.'PB3DOT_1')read(valstr,223,err=3)pbd3(1),ie(104)  ! (PB3DOT_1)       
       if(chstr(1:6).eq.'X2DOT ') read(valstr,223,err=3) Xd2(1),ie(102)   ! (X2DOT)
       if(chstr(1:7).eq.'X2DOT_1')read(valstr,223,err=3) Xd2(1),ie(102)   ! (X2DOT_1)
       if(chstr(1:5).eq.'DTHET')read(valstr,223,err=3) dth(1),  ie(24)    ! (DTHETA)
       if(chstr(1:4).eq.'GAMM') read(valstr,223,err=3) gamma(1),ie(17)    ! (GAMMA)
       if(chstr(1:7).eq.'GAMMA_1')read(valstr,223,err=3) gamma(1),ie(17)  ! (GAMMA_1)       
       if(chstr(1:2).eq.'M2')   read(valstr,223,err=3) am2(1),ie(18)      ! (M2)
       if(chstr(1:4).eq.'M2_1') read(valstr,223,err=3) am2(1),ie(18)      ! (M2_1)       
       if(chstr(1:2).eq.'SI')   read(valstr,223,err=3) si(1),ie(19)       ! (SI)
       if(chstr(1:4).eq.'SI_1') read(valstr,223,err=3) si(1),ie(19)       ! (SI_1)
C       if(chstr(1:4).eq.'PPNG') read(valstr,222,err=3) gamma?             ! (PPNGA)
       if(chstr(1:4).eq.'MTOT') read(valstr,222,err=3) am                 ! (MTOT)
       if(chstr(1:2).eq.'BP')   read(valstr,222,err=3) bp                 ! (BP)
       if(chstr(1:3).eq.'BPP')  read(valstr,222,err=3) bpp                ! (BPP)
       if(chstr(1:2).eq.'DR')   read(valstr,222,err=3) dr                 ! (DR)
       if(chstr(1:2).eq.'A0')   read(valstr,222,err=3) a0                 ! (A0)
       if(chstr(1:2).eq.'B0')   read(valstr,222,err=3) b0                 ! (B0)
*(Model of N.Wex-EPS1,EPS2,TASC)
         if(chstr(1:4).eq.'EPS1') read(valstr,223,err=3) eps1, ie(14)    ! (EPS1        E1_1)
         if(chstr(1:4).eq.'EPS2') read(valstr,223,err=3) eps2, ie(16)    ! (EPS2        OM_1)
         if(chstr(1:4).eq.'TASC') read(valstr,223,err=3) tasc, ie(15)    ! (TASC        T0  )
* (Another commands)
         if(chstr(1:4).eq.'GENP') IE(55)=1
         if(chstr(1:4).eq.'GENT') IE(55)=2
         if(IE(55).ge.1) then
          found=.true.
          goto 800               ! exit if command GENTOA or GENPER found 
          endif ! (IE(55)
         enddo ! (i=1,1000)
 800  continue
         
*(Count number of companions to PSR)
         if (a1(1).ne.0.0d0) norb=1
         if (a1(2).ne.0.0d0) norb=2
         if (a1(3).ne.0.0d0) norb=3
*(Convert f 2 p)
         if(p0.eq.0.0d0) then
         if (f(1).eq.0.0d0) f(1)=1.0d0
         p0=1.0d0/f(1)
         p1=-f(2)/f(1)/f(1)
         p2=-f(3)/f(1)/f(1)+2.0*f(2)*f(2)/f(1)/f(1)/f(1)
         if(f(3).eq.0.0d0) p2=0.0d0
         endif ! (p0)
*(Convert p 2 f)
         if(f(1).eq.0.0d0) then
         p1=p1*1.0d-15
         if (p0.eq.0.0d0) p0=1.0d0
         f(1)=1.0d0/p0
         if(f(2).eq.0.0) f(2)=-p1/p0/p0
         if(f(3).eq.0.0) f(3)=-p2/p0/p0+2.0*p1*p1/p0/p0/p0
         if(p2.eq.0.0d0) f(3)=0.0d0
         endif ! (p0)
         
* (Convert DM derivatives units from pc/cm3/yr**i, i=1..5 to pc/cm3/s**i)
        do i=2,7
        dm(i) = dm(i)/86400.0d0**(i-1)/365.25**(i-1)
        enddo ! i=2,7
*(Timing model)
         k=0
         if(chbinm(1:4). eq.'NONE') k=0
         if(chbinm(1:2). eq.'DD')   k=1
          if(chbinm(1:4).eq.'DDGR') k=1
          if(chbinm(1:3).eq.'DDT')  k=1
          if(chbinm(1:3).eq.'MSS')  k=1
          if(chbinm(1:4).eq.'ELL1') k=1
         if(chbinm(1:2). eq.'BT')   k=2
          if(chbinm(1:2).eq.'EH')   k=2
          if(chbinm(1:3).eq.'H88')  k=2
         if(chbinm(1:3). eq.'BT+')  k=3
         if(chbinm(1:4). eq.'BT++') k=4
          if(chbinm(1:4).eq.'BT1P') k=4
          if(chbinm(1:4).eq.'BT2P') k=4
         if(chbinm(1:3). eq.'DD+')  k=5
         IE(51)=k
*(Set parameters of pointer array)
         IE(52)=0                   ! Weighted LSQ
         IE(53)=0                   ! account EOP
         IE(54)=0                   ! write residuald 2 hdisk
         if (NITER.lt.1)  IE(52)=2  ! no fit if NITER=0
         if (NITER.eq.0)  NITER=2   ! fit jumps only
         if (NITER.eq.1)  NITER=2   ! fit jumps only
         
*(convert to necessary units)
         
* we don't use ELL1,ELL2,TASC just - convert them to om(1),e(1),t0(1)
* 
         if(eps2.ne.0.0d0) then
         om(1)=datan2(eps1,eps2)
         if (om(1).lt.0.0d0) om(1)=om(1)+pi*2.0d0
         e(1)=eps1/dsin(om(1))
         om(1)=om(1)*180.0d0/pi
         t0(1)=tasc+om(1)/360.0*pb(1)
         endif

         aps=aps*0.0001d0
         dps=dps*0.0001d0
         if (px.ne.0.0) then
            px=1.0d0/px * 1.0d+3 ! px to distance pc
            else
            px=1.0d9
         endif ! (px.ne.0)
         if(norb.ge.1) then
         do i=1,norb
          pb(i)   =pb(i)*86400.0d0
          a1dot(i)=a1dot(i)*1.0d-12            ! -> s/s
          edot(i) = edot(i)*1.0d-12            ! -> s/s
          if(t0(i).lt.2400000.5d0) t0(i)=t0(i)+2400000.5d0
         enddo
         endif ! (norb>=1)
      return
 900  continue
      found=.false.
      write(*,*) ' Parameter file - not found'
      return
 3    write(*,*) ' STOP: Error in line:',chstr
      STOP
      end

*********************************************************************

      subroutine upval(str,strval)

*  Converts string (str) to upper case
*  and copies its contents starting from
*  first blank into string (strval).
*  Symbols ':' are deleted in both strings
*
*  The format of strval is: (f24.12,i1)-value and flag
*  can be used for further formatted reading
      
      implicit none
      integer i,j,ileng,ipos1,ipos2,ipos3,ipos4
      character str*80,strval*80
      j=0

c delete ':' in string (str)
      do i=1,len(str)
      if (str(i:i).ne.':') then
       j=j+1
       strval(j:j)=str(i:i)
      endif
      enddo
      str=strval
c convert to upcase
      ileng=len(str)
      do i=1,ileng
      j=ichar(str(i:i))
      if(j.ge.97.and.j.le.122) str(i:i)=char(j-32)
      enddo
c copy string from 1st blank:-------------------
c DM         2.4567899123       1      0.000013
c   ipos1    ipos2       ipos3  ipos4
c-----------------------------------------------
      ipos1=1
      ipos2=1
      ipos3=1
      ipos4=1
      do i=1,ileng
      if(str(i:i).eq.' ') then
                 ipos1=i
                 goto 10
                 endif
      enddo
 10   do i=ipos1,ileng
      if(str(i:i).ne.' ') then
                 ipos2=i
                 goto 20
                 endif
      enddo
 20   do i=ipos2,ileng
      if(str(i:i).eq.' ') then
                 ipos3=i
                 goto 25
                 endif
      enddo
 25   do i=ipos3,ileng
      if(str(i:i).ne.' ') then
                 ipos4=i
                 goto 30
                 endif
      enddo
 30   continue
      strval='                                       '
      if (ipos2.le.ipos1) then
      strval='                                       '
      else
      strval(1:ipos3-ipos2+1)=str(ipos2:ipos3)
      if (ipos4.gt.1) strval(25:25)=str(ipos4:ipos4) ! copy flag 0/1
      endif
c check
c      write(*,'(4i3,a,a)') ipos1,ipos2,ipos3,ipos4,'>>',strval(1:26)
      return
      end
*********************************************************************

