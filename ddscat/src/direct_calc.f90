    SUBROUTINE DIRECT_CALC(IX,IY,IZ,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                           PYDDX,PZDDX,CXZC)
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_0,ONLY: IDIPINT
      USE DDCOMMON_10,ONLY: MYID
#ifdef openmp
      USE OMP_LIB   !! Art omp v2.0 supplies variable and function definitions
#endif
      IMPLICIT NONE

! arguments:

      INTEGER :: IPBC,IX,IY,IZ,NX,NY,NZ
      REAL(WP) :: AKD,AKD2,GAMMA,PYDDX,PZDDX
      REAL(WP) :: &
         AK(3),   &
         DX(3)
      COMPLEX (WP) ::                                             &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6)

! local variables

      CHARACTER :: CMSGNM*66
      INTEGER :: I,ICXZC,II,IMIN,J,JCXZC,                  &
                 JJ,JMIN,JPBC,JPY,JPY1,JPY2,JPZ,JPZ1,JPZ2, &
                 K,KCXZC,KMIN,M
      REAL(WP) :: BETA,CIMINUS,CIPLUS,COSKFR,COSKR,DAMP,    &
                  KF,KFR,KPERP,KPERP2,KR,KRF,               &
                  PHASF,PHASY,PHASYZ,PHI,PHIMAX,PI,         &
                  R,R02,R2,R3,RF,RJPY,RJPZ,RPERP,           &
                  SIMINUS,SIPLUS,SINKFR,SINKR,              &
                  T0,T1,T2,TERM0,TERM1,TERM2Y,TERM2Z,TERM3, &
                  X0,X2,X2Y2,XX,Y0,Y2,Y2Z2,YF,YHI,YLO,      &
                  Z0,Z2,ZF,ZHI,ZLO
      REAL(WP) :: X(3)
      COMPLEX(WP) :: CX2PIH,CXA0,CXA1,CXA2,CXEXPIKR,CXFAC, &
                     CXG0,CXG1,CXG2,CXI,CXIKR,CXPHAS,CXZERO
                     
      COMPLEX(WP) :: DCXSUM(6)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      SAVE CXZERO,CXI
      DATA CXI/(0._WP,1._WP)/,CXZERO/(0._WP,0._WP)/

!-----------------------------------------------------------------------
! subroutine DIRECT_CALC v8

! calculates electric field at (I,J,K) due to dipole at (1,1,1)
!                                      plus its replicas
! given:
!        IX,IY,IZ = 1, 1, 1 to do octant I>0, J>0, K>0 only
!                   1, 1,-1              I>0, J>0, all K
!                   1,-1, 1              I>0, K>0, all J
!                   1,-1,-1              I>0, all J, all K
!                  -1, 1, 1              all I, J>0, K>0
!                  -1, 1,-1              all I, J>0, all K
!                  -1,-1, 1              all I, all J, K>0
!                  -1,-1,-1              all I, all J, all K

!        NX,NY,NZ = size of first octant (I = 1 -> NX,
!                                         J = 1 -> NY,
!                                         K = 1 -> NZ)

!        IPBC     = 0 for isolated target
!                 = 1 to use periodic boundary conditions
!                     (y direction, z direction, or both y and z directions)
!                   N.B.: IPBC affects dimensioning of CXZC

!        GAMMA     = parameter controlling range of summations over
!                    replica dipoles.
!                    PHIMAX=1./gamma
!
!        PYDDX     = period of lattice in y direction/d
!        PZDDX     = period of lattice in z direction/d
!        DX(1-3)   = lattice spacing in x,y,z direction, in units of
!                    d = n**(-1/3).
!        AK(1-3)   = k(1-3)*d, where k = k vector in vacuo
!        AKD       = |k|d
!        CXZC      = array with dimensions
!                    (NX+1)*(NY+1)*(NZ+1)*6 if IPBC=0
!                    (2*NX)*(2*NY)*(2*NZ)*6 if IPBC=1
!
! and through module DDCOMMON_0:
!
!	 IDIPINT   = 0 for point dipole interaction method
!		   = 1 for filtered coupled dipole (FCD) interaction method 

! returns:
!        CXZC(I,J,K,M) = a_M to calculate E field at (I,J,K)
!                        produced by dipole P at (1,1,1) and replica dipoles
!                        if IPBC > 0
!                        -E_x = a_1*P_x + a_2*P_y + a_3*P_z
!                        -E_y = a_2*P_x + a_4*P_y + a_5*P_z
!                        -E_z = a_3*P_x + a_5*P_y + a_6*P_z 

! note that even for short-range interaction, for PBC calculations with
! periodicity in the y-direction, we need to extend sums
! from JPY=-1 to JPY=+1 to include "edge" interactions with next TUC
! and similarly if the target is periodic in the z-direction

! PYDDX=PYDDY=0. for nonperiodic case (IPBC=0)
! PYDDX>0 for PBC in y direction
! PZDDX>0 for PBC in z direction

! Adapted from original subroutine ESELF written by Jeremy Goodman
! history:
! 06.09.28 (BTD) DIRECT_CALC works for isolated target
! 06.09.30 (BTD) further modifications to DIRECT_CALC
!                change definition of BETA, so that given BETA 
!                determines range/wavelength
! 08.03.11 (BTD) v7.0.5
!                * eliminate BETA, introduce ALPHA=BETA**(0.25)
! 08.03.15 (BTD) * added code to estimate time to completion
!                  may need to change way this is done before running under MPI
! 08.03.15 (BTD) * added DDCOMMON_10 to communicate value of MYID.
! 08.04.20 (BTD) * changed notation: ALPHA -> GAMMA
! 08.06.05 (ASL) v7.0.6: 
!                * Arthur S. Lazanoff added OpenMP directives
! 08.07.02 (BTD) * Added call to CPU_TIME(T0) outside of PARALLEL region
! 12.04.16 (BTD) v7.2:
!                * add check to verify number of OpenMP threads               
! 12.07.06 (BTD) v7.2.3 edited comments
! 12.08.02 (IYW) added DIPINT to arg list and FCD Green func coefficients
! 12.08.10 (BTD) **** need to examine how these changes affect OpenMP ***
! 12.08.11 (BTD) v7.3 and eself_v7
!                * removed DIPINT from argument list
!                * added IDIPINT from module DDCOMMON_0
! 12.12.21 (BTD) * cosmetic changes, added comments
! 12.12.27 (BTD) * corrected errors in case IDPINT=1 (FILTRDD)
!                * added comments relating code to notation of
!                  Gay-Balmaz & Martin (2002) [Comp. Phys. Comm. 144, 111] 
!                  and Yurkin, Min & Hoekstra (2010) [PRE 82, 036703]
! 12.12.29 (BTD) * added phase correction for DIPINT=1
!                  (needed for PBC calculations)
! 13.01.05 (BTD) * OpenMP bug fix: add ICXZC to list of PRIVATE variables
!                  change
!                  PRIVATE(JCXZC,JPZM,KCXZC)               &
!                  to
!                  PRIVATE(ICXZC,JCXZC,JPZM,KCXZC)         &
! 13.01.07 (BTD) * corrected typo PRIVATE(ICXCZ -> PRIVATE(ICXZC
!                  (noted 13.01.07 by Zhenpeng Qin)
! 13.12.24 (BTD) v8 modifications to properly treat sums over replicas
!                   we now determine size of Fresnel zone for each
!                   dipole and its replicas
! 13.12.26 (BTD) * further changes to use new form of suppression factor
!                  exp[-beta*(phi-phiF)^4] with beta = 8
! 13.12.27 (BTD) * further modifications.  Now use exp(-beta*PHI^4)
!                  with beta chosen so that factor = 1e-7 when PHI=PHIMAX
! 13.12.29 (BTD) * set PHIMAX = 1/GAMMA so that PHIMAX can be
!                  modified from ddscat.par
! 19.07.10 (BTD) * added error-catcher suggested by Fangzhou Liu to
!                  ensure that IDIPINT is properly set in MPI
! end history
! Copyright (C) 1993,1994,1996,1997,1999,2000,2003,2004,2005,2006,2007,
!               2008,2011,2012,2013,2019 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------
! PI and KF are used if IDIPINT=1 (filtered coupled dipole)

      PI=4._WP*DATAN(1.D0)
      KF=PI

! maximum phase shift defining extent of integration over replica dipoles
! should do some experiments with different PHIMAX values
! BETA=16/PHIMAX**4: 
!    suppression factor = 1e-16 = 1.1e-7 for PHI=PHIMAX
!                       = 0.99840            PHI=0.1*PHIMAX = 5
!                       = 0.97472            PHI=0.2*PHIMAX = 10

! we will continue to use input parameter gamma
! and set PHIMAX = 1/gamma
! then gamma = 0.02 -> PHIMAX = 50
!              0.01 -> PHIMAX = 100

      PHIMAX=1._WP/GAMMA

      BETA=16._WP/PHIMAX**4

      IF(IPBC>0)THEN

! calculate various quantities that will be needed later when
! summing over replica dipoles
! JPBC=0 if not periodic
!      1 if periodic in y only
!      2 if periodic in z only
!      3 if periodic in both y and z

         IF(PYDDX<=0._WP)THEN
! periodic in z only
            JPBC=2
            KPERP2=AK(1)**2+AK(2)**2
         ELSE
            IF(PZDDX<=0._WP)THEN
! periodic in y only
               JPBC=1
               KPERP2=AK(1)**2+AK(3)**2
            ELSE
! periodic in both y and z
               JPBC=3
               KPERP2=AK(1)**2
            ENDIF
         ENDIF
! sanity check
         IF(KPERP2<=0._WP)CALL ERRMSG('FATAL','DIRECT_CALC',    &
                                      ' KPERP2 <= 0: unphysical')
         KPERP=SQRT(KPERP2)
         TERM0=AKD2*PHIMAX**2/KPERP2
         TERM1=2._WP*AKD*PHIMAX
         TERM2Y=AK(2)*PHIMAX/KPERP2
         TERM2Z=AK(3)*PHIMAX/KPERP2

      ELSE
! no periodicity
         JPBC=0
      ENDIF

!*** diagnostic
!      write(0,*)'direct_calc_v8 ckpt 1'
!***
! Compute 6 independent elements of 3x3 symmetric matrix A_jk, where
! A_jk*P_k = -electric field at location j due to dipole P_k at location

! A_jk = (a_1  a_2  a_3)
!        (a_2  a_4  a_5)
!        (a_3  a_5  a_6)_jk

! initialize CX(I,J,K,M) = a_M for electric field at (I,J,K)
!                          produced by a dipole at (1,1,1)
!                          and replica dipoles (if PYD or PYZ are
!                          nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! X0,Y0,Z0 = X(I,J,K) - X(1,1,1) = vector from dipole location (1,1,1)
!                                  to point where E is to be calculated
! IX = +1 -> IMIN = 1
! IX = -1 -> IMIN = 2-NX   (1-IMIN = NX-1)

      IMIN=1+(1-NX)*(1-IX)/2
      JMIN=1+(1-NY)*(1-IY)/2
      KMIN=1+(1-NZ)*(1-IZ)/2

!*** diagnostic
!      write(0,fmt='(a,3i5)')'direct_calc_v8 ckpt 2, imin,jmin,kmin=', &
!                            imin,jmin,kmin
!***

! Determine elapsed cpu time for these sums

      CALL CPU_TIME(T0)

#ifdef openmp
! fork a team of threads 
!$OMP PARALLEL               &
!$OMP&   PRIVATE(NTHREADS,TID)

      TID=OMP_GET_THREAD_NUM()

!*** diagnostic
!      write(0,*)'direct_calc_v8 ckpt 3: hello world from thread = ',TID
!***
! only master thread does this:

      IF(TID.EQ.0)THEN
         NTHREADS=OMP_GET_NUM_THREADS()
         WRITE(CMSGNM,FMT='(A,I3)')'number of OpenMP threads =',NTHREADS
         CALL WRIMSG('DIRECT_CALC',CMSGNM)
      ENDIF
!$OMP END PARALLEL
#endif

! 2013.01.04 (BTD) Zhenpen Qin reports error within following parallelization
! 2013.01.05 (BTD) ICXZC was inadvertently omitted from PRIVATE variables list
!                  add it: change
!                  PRIVATE(JCXZC,JPZM,KCXZC) to
!                  PRIVATE(ICXZC,JCXZC,JPZM,KCXZC)
#ifdef openmp
! 2012.04.27 (BTD) added R to private variables
!                  removed FLUSH directive (should not be needed)
! 2013.12.24 (BTD) v8: added a number of variables to private list
!$OMP PARALLEL                                           &
!$OMP&   PRIVATE(I,ICXZC,II)                             &
!$OMP&   PRIVATE(J,JCXZC,JJ,JPY,JPY1,JPY2,JPZ,JPZ1,JPZ2) &
!$OMP&   PRIVATE(K,KCXZC,M)                              &
!$OMP&   PRIVATE(CIMINUS,CIPLUS,COSKFR,COSKR,DAMP)       &
!$OMP&   PRIVATE(KFR,KR,KRF,PHASF,PHASY,PHASYZ,PHI)      &
!$OMP&   PRIVATE(R,R2,R3,RF,RJPY,RJPZ,RPERP)             &
!$OMP&   PRIVATE(SIMINUS,SIPLUS,SINKFR,SINKR,TERM3)      &
!$OMP&   PRIVATE(X,X0,X2,X2Y2,XX)                        &
!$OMP&   PRIVATE(Y0,Y2,Y2Z2,YF,YHI,YLO)                  &
!$OMP&   PRIVATE(Z0,Z2,ZF,ZHI,ZLO)                       &
!$OMP&   PRIVATE(CX2PIH,CXA0,CXA1,CXA2,CXEXPIKR,CXFAC)   &
!$OMP&   PRIVATE(CXG0,CXG1,CXG2,CXIKR,CXPHAS,DCXSUM)
!$OMP DO
#endif

      DO K=KMIN,NZ              !- loop over K

         Z0=REAL(K-1,KIND=WP)*DX(3)
         Z2=Z0**2
         IF(K>0)THEN
            KCXZC=K
         ELSE
            KCXZC=2*NZ+K
         ENDIF
         DO J=JMIN,NY              !-- loop over J
            Y0=REAL(J-1,KIND=WP)*DX(2)
            Y2=Y0**2
            Y2Z2=Y2+Z2
            IF(J>0)THEN
               JCXZC=J
            ELSE
               JCXZC=2*NY+J
            ENDIF
            DO I=IMIN,NX              !--- loop over I

! for first dipole, determine time to sum over replicas

               IF(I.EQ.IMIN.AND.J.EQ.JMIN.AND.          &
                  K.EQ.KMIN.AND.MYID==0)CALL CPU_TIME(T1)

               X0=REAL(I-1,KIND=WP)*DX(1)
               X2=X0**2

! calculate JPY1, JPY2, JPZ1, JPZ2
!           KRF = k*R_F, R_F=distance from (I,J,K) to center of Fresnel zone

               JPY1=0
               JPY2=0
               JPZ1=0
               JPZ2=0
               IF(JPBC==1)THEN
                  RPERP=SQRT(X2+Z2)
                  KRF=AKD2*RPERP/KPERP
                  YF=Y0-AK(2)*RPERP/KPERP
                  PHASF=AK(2)*YF
                  TERM3=SQRT(TERM0+TERM1*RPERP)/KPERP
                  YLO=(YF-TERM2Y-TERM3)/PYDDX
                  YHI=(YF-TERM2Y+TERM3)/PYDDX
                  JPY1=INT(YLO+0.999999_WP)
                  JPY2=INT(YHI)

!*** diagnostic
!       write(0,fmt='(a,3i5,a,1pe10.3,a,2i5)')                       &
!             'direct_calc_v8 ckpt 4: i,j,k=',i,j,k,' rperp=',rperp, &
!             ' jpy1,jpy2=',jpy1,jpy2
!***
               ELSEIF(JPBC==2)THEN
                  RPERP=SQRT(X2+Y2)
                  KRF=AKD2*RPERP/KPERP
                  ZF=Z0-AK(3)*RPERP/KPERP
                  PHASF=AK(3)*ZF
                  TERM3=SQRT(TERM0+TERM1*RPERP)/KPERP
                  ZLO=(ZF-TERM2Z-TERM3)/PZDDX
                  ZHI=(ZF-TERM2Z+TERM3)/PZDDX
                  JPZ1=INT(ZLO+0.999999_WP)
                  JPZ2=INT(ZHI)
               ELSEIF(JPBC==3)THEN
                  RPERP=ABS(X0)
                  KRF=AKD2*RPERP/KPERP
                  YF=Y0-AK(2)*RPERP/KPERP
                  ZF=Z0-AK(3)*RPERP/KPERP
                  PHASF=AK(2)*YF+AK(3)*ZF
                  TERM3=SQRT(TERM0+TERM1*RPERP)/KPERP
                  YLO=(YF-TERM2Y-TERM3)/PYDDX
                  YHI=(YF-TERM2Y+TERM3)/PYDDX
                  JPY1=INT(YLO+0.999999_WP)
                  JPY2=INT(YHI)
                  ZLO=(ZF-TERM2Z-TERM3)/PZDDX
                  ZHI=(ZF-TERM2Z+TERM3)/PZDDX
                  JPZ1=INT(ZLO+0.999999_WP)
                  JPZ2=INT(ZHI)
               ENDIF

               IF(I>0)THEN
                  ICXZC=I
               ELSE
                  ICXZC=2*NX+I
               ENDIF
               X(1)=X0
               DO M=1,6
                  DCXSUM(M)=CXZERO
               ENDDO

!*** diagnostic
!               write(0,*)'direct_calc_v8 ckpt 5'
!               write(0,*)'   x2=',x2
!***

! JPY=0, JPZ=0 gives E field from dipole at (1,1,1)
! general JPY,JPZ gives E field from dipole at
! (1,JPY*NPY+1,JPZ*NPZ+1)
! replica dipoles have same magnitude of polarization as dipole (1,1,1),
! but different phase.
! PHASYZ = phase of replica dipole - phase of dipole (1,1,1)

               DO JPY=JPY1,JPY2         !---- loop over JPY
                  RJPY=REAL(JPY,KIND=WP)*PYDDX*DX(2)
                  X(2)=Y0-RJPY
                  X2Y2=X2+X(2)**2
                  PHASY=AK(2)*RJPY

                  DO JPZ=JPZ1,JPZ2         !----- loop over JPZ

!*** diagnostic
!      write(0,*)'jpz=',jpz
!***

                     RJPZ=REAL(JPZ,KIND=WP)*PZDDX*DX(3)
                     X(3)=Z0-RJPZ
                     R2=X2Y2+X(3)**2

! skip the self-interaction (R=0) case (I=J=K=1 and JPY=JPZ=0)
! this also avoids possible problems if replica dipoles fall close
! to "ghost" dipoles in "padding" locations

                     IF(R2>1.E-6_WP)THEN

                        R=SQRT(R2)
                        R3=R*R2

! PHASYZ = phase at (1,JPY*NPY+1,JPZ*NPZ+1) - phase at (1,1,1)
! PHI = phase shift at (I,J,K) from replica of (1,1,1) relative to
!       phase shift from center of Fresnel zone

                        PHASYZ=PHASY+AK(3)*RJPZ
                        KR=AKD*R
                        CXIKR=CXI*KR
                        IF(JPBC==0)THEN
                           DAMP=0._WP
                        ELSE
                           PHI=KR-KRF+PHASYZ-PHASF
                           DAMP=BETA*PHI**4
                        ENDIF
!*** diagnostic
!                        write(0,fmt='(a,2i5,a,1pe10.3,a,1pe10.3)')    &
!                           'direct_calc_v8 ckpt 6: jpy,jpz=',jpy,jpz, &
!                           ' Phi=',PHI,' damp=',damp
!***
! include artificial factor exp[-DAMP] to assist convergence

			IF(IDIPINT==0)THEN
                           CXPHAS=EXP(CXIKR+CXI*PHASYZ-DAMP)/R3
                           CXFAC=(1._WP-CXIKR)/R2

!*** diagnostic
!                           write(0,*)'direct_calc_v8 ckpt 7'
!                           write(0,*)'  cxterm=',cxterm
!                           write(0,*)'  cxphas=',cxphas
!                           write(0,*)'   cxfac=',cxfac
!***
                           
! II=1      -> M=1   a_1 (xx)
!      JJ=2      2   a_2 (xy)
!         3      3   a_3 (xz)
!    2           4   a_4 (yy)
!         3      5   a_5 (yz)
!    3           6   a_6 (zz)

                           M=0
                           DO II=1,3              !------- loop over II
                              M=M+1
                              XX=X(II)**2
                              DCXSUM(M)=DCXSUM(M)-CXPHAS*(AKD2*(XX-R2)+ &
                                        CXFAC*(R2-3._WP*XX))
                              IF(II<3)THEN
                                 DO JJ=II+1,3        !------- loop over JJ
                                    M=M+1
                                    DCXSUM(M)=DCXSUM(M)-                       &
                                           CXPHAS*X(II)*X(JJ)*(AKD2-3._WP*CXFAC)
                                 ENDDO               !------- end loop over JJ
                              ENDIF
                           ENDDO                  !------ end loop over II

			ELSEIF(IDIPINT==1)THEN

! Expressions for the Green coefficients correspond to those in 
!   A library for computing the filtered and non-filtered 3D Greens
!   tensor associated with infinite homogeneous space and surfaces 
! by P. Gay-Balmaz and O. Martin (2002; Computer Physics Communications,
! 144, 111-120) apart from a factor of 4*pi.

                           CXEXPIKR=EXP(CXIKR)
                           CXPHAS=EXP(CXI*PHASYZ-DAMP)
                           KFR=KF*R
                           CALL CISI(KFR+KR,CIPLUS,SIPLUS)
                           CALL CISI(KFR-KR,CIMINUS,SIMINUS)
!*** diagnostic
!                           write(57,fmt='(1pe10.3,1p2e11.3)')kfr+kr,     &
!                                                             ciplus,siplus
!                           write(57,fmt='(1pe10.3,1p2e11.3)')kfr-kr,       &
!                                                             ciminus,siminus
!***
                           COSKR=COS(KR)
                           SINKR=SIN(KR)
                           COSKFR=COS(KFR)
                           SINKFR=SIN(KFR)

! CXA0 = alpha(R) of Gay-Balmaz & Martin 2002
! CXA1 = (d/dR) alpha   = alphaprime
! CXA2 = (d2/dR2) alpha = alphadoubleprime

                           CXA0=SINKR*(CIPLUS-CIMINUS)+ &
                                COSKR*(PI-SIPLUS-SIMINUS)
                           CXA1=AKD*SINKR*(-PI+SIPLUS+SIMINUS)+         &
                                AKD*COSKR*(CIPLUS-CIMINUS)-2._WP*SINKFR/R
                           CXA2=AKD2*(SINKR*(CIMINUS-CIPLUS)+     &
                                      COSKR*(SIPLUS+SIMINUS-PI))+ &
                                2._WP*(SINKFR-KFR*COSKFR)/R2

! Yurkin, Min & Hoekstra 2010 define G_ij such that G_ij*P_j = E field at i
! so that elements of G_ij are same as elements of our matrix CXZC
! notation:
! CXG0   = g_F(R)                     of YMH2010
! CXG1   = (d/dR) g_F   = g_Fprime    of YMH2010
! CXG2   = (d2/dR2) g_F = g_Fdblprime of YMH2010  
! CX2PIH = 2*pi*h(R)                  of YMH2010

! CXG0 = exp(ikR)/R - alpha/(pi*R)
! CXG1 = exp(ikR)*(-1+ikR)/R^2 - alphaprime/(pi*R) + alpha/(pi*R^2)
! CXG2 = exp(ikR)*(2-2ikR-(kR)^2)/R^3 - 
!       (2*alpha-2R*alphaprime+R^2*alphadoubleprime)/(pi*R^3)

                           CXG0=(CXEXPIKR-CXA0/PI)/R
                           CXG1=(CXEXPIKR*(CXIKR-1._WP)+(CXA0-CXA1*R)/PI)/R2
                           CXG2=(CXEXPIKR*(2._WP-2._WP*CXIKR-KR**2)- &
                                 (2._WP*(CXA0-R*CXA1)+R2*CXA2)/PI)/R3
                           CX2PIH=(SINKFR-KFR*COSKFR)/(PI*R3)
                           M=0
                           DO II=1,3   ! do II
                              M=M+1
! diagonal elements 
! II=1: M=1 xx
! II=2: M=4 yy
! II=3: M=6 zz
                              DCXSUM(M)=DCXSUM(M)+CXPHAS*(AKD2*CXG0+ &
                                        CXG1/R+(2._WP/3._WP)*CX2PIH+ &
                                        (CXG2/R2-CXG1/R3)*X(II)**2)
                              IF(II<3)THEN
                                 DO JJ=II+1,3  ! do JJ
                                    M=M+1
! off-diagonal elements
! II=1: M=2 xy
! II=1: M=3 xz
! II=2: M=5 yz
                                    DCXSUM(M)=DCXSUM(M)+CXPHAS*           &
                                              (CXG2/R2-CXG1/R3)*X(II)*X(JJ)
                                 ENDDO  ! enddo JJ
                              ENDIF
                           ENDDO  ! enddo II
                        ELSE
! 190710 error-catcher suggested by Fangzhou Liu:
                           WRITE(CMSGNM,FMT='(A,I4,A)')      &
                                'IDIPINT=',IDIPINT,' is invalid'
                           CALL ERRMSG('FATAL','DIRECT_CALC',CMSGNM)
			ENDIF  ! end elseif(IDIPINT==1)

                     ENDIF  ! endif(R2 > 1e-6)
                  ENDDO                     !----- end loop over JPZ
               ENDDO                     !---- end loop over JPY
!*** diagnostic
!               write(0,*)'direct_calc_v8 ckpt 8'
!               write(0,*)'   icxzc,jcxzc,kcxzc=',icxzc,jcxzc,kcxzc
!***
               DO M=1,6
                  CXZC(ICXZC,JCXZC,KCXZC,M)=DCXSUM(M)
!*** diagnostic
!                  write(0,fmt='(a,i1,a,2f10.4)')'   dcxsum(',m,')=',dcxsum(m)
!***
               ENDDO
               IF(I==IMIN.AND.J==JMIN.AND.K==KMIN.AND.MYID==0)THEN

! NB: this call to CPU_TIME occurs within the first thread
!     (I==IMIN, J==JMIN, K==KMIN)

                  CALL CPU_TIME(T2)
                  T2=REAL((NX-IMIN+1)*(NY-JMIN+1)*(NZ-KMIN+1),KIND=WP)*(T2-T1)
                  IF(T2<600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                         &
                        'Estimated total cputime required by DIRECT_CALC=', &
                        T2,' cpu-sec'
                  ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                         &
                        'Estimated total cputime required by DIRECT_CALC=', &
                        (T2/60._WP),' cpu-min'
                  ELSE
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                         &
                        'Estimated total cputime required by DIRECT_CALC=', &
                        (T2/3600._WP),' cpu-hr'
                  ENDIF
                  CALL WRIMSG('DIRECT_CALC',CMSGNM)
                  IF(PYDDX*PZDDX>0)THEN
                     WRITE(CMSGNM,FMT='(A)')          &
                          'cputime scales as 1/gamma^2'
                  ELSE
                     IF(PYDDX+PZDDX==0)THEN
                        WRITE(CMSGNM,FMT='(A)')               &
                              'cputime is independent of gamma'
                     ELSE
                        WRITE(CMSGNM,FMT='(A)')           &
                              '[cputime scales as 1/gamma]'
                     ENDIF
                  ENDIF
                  CALL WRIMSG('DIRECT_CALC',CMSGNM)
               ENDIF
            ENDDO                     !--- end loop over I
         ENDDO                     !-- end loop over J

! 2012.04.27 (BTD) it is not necessary to have a flush(cxzc) here
!                  flush is implied by both $omp END DO and 
!                  $omp END PARALLEL directives

      ENDDO                     !- end loop over K

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

!*** diagnostic
!      write(0,*)'direct_calc_v8 ckpt 9'
!***
      CALL CPU_TIME(T2)
!*** diagnostic
!      write(0,*)'direct_calc_v8 ckpt 10'
!***
      T2=T2-T0
      IF(MYID==0)THEN
         IF(T2<600.D0)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                   &
                  'Actual cputime to complete DIRECT_CALC=', &
                  T2,' cpu-sec'
         ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                   &
                  'Actual cputime to complete DIRECT_CALC=', &
                  (T2/60._WP),' cpu-min'
         ELSE
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                   &
                  'Actual cputime to complete DIRECT_CALC=', &
                  (T2/3600._WP),' cpu-hr'
         ENDIF
         CALL WRIMSG('DIRECT_CALC',CMSGNM)
      ENDIF

!*** diagnostic
!      write(0,*)'direct_calc_v8 ckpt 11'
!***

! If IPBC=1: set the elements with ICXZC=NX+1 or JCXZC=NY+1
!            or KCXZC=NZ+1 to zero

      IF(IPBC==1)THEN

         DO M=1,6

#ifdef openmp
!$OMP PARALLEL                    &
!$OMP&   PRIVATE(ICXZC,JCXZC,KCXZC)
!$OMP DO
#endif

            DO KCXZC=1,2*NZ
               DO JCXZC=1,2*NY
                  CXZC(NX+1,JCXZC,KCXZC,M)=CXZERO
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO KCXZC=1,2*NZ

               DO ICXZC=1,2*NX
                  CXZC(ICXZC,NY+1,KCXZC,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO JCXZC=1,2*NY

               DO ICXZC=1,2*NX
                  CXZC(ICXZC,JCXZC,NZ+1,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

         ENDDO
 
      ENDIF
      RETURN
 6700 FORMAT(F8.2,' cpu-sec')
 6701 FORMAT(F8.2,' cpu-min')
 6702 FORMAT(F8.2,' cpu-hr')
    END SUBROUTINE DIRECT_CALC
