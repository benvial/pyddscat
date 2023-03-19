    SUBROUTINE DIRECT_CALCB(IX,IY,IZ,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                           PYDDX,PZDDX,CXZG)
      USE DDPRECISION,ONLY: WP
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
         CXZG(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),3)

! local variables

      CHARACTER :: CMSGNM*66
      INTEGER :: I,ICXZG,II,IMIN,                                  &
                 J,JPBC,JCXZG,JJ,JMIN,JPY,JPY1,JPY2,JPZ,JPZ1,JPZ2, &
                 K,KCXZG,KMIN,M
      REAL(WP) :: BETA,DAMP,KPERP,KPERP2,KR,KRF,            &
                  PHASF,PHASY,PHASYZ,PHI,PHIMAX,            &
                  R,R2,RF,RJPY,RJPZ,RPERP,                  &
                  T0,T1,T2,TERM0,TERM1,TERM2Y,TERM2Z,TERM3, &
                  X0,X2,X2Y2,XX,Y0,Y2,Y2Z2,YF,YHI,YLO,      &
                  Z0,Z2,ZF,ZHI,ZLO
      REAL(WP) :: X(3)
      COMPLEX(WP) :: CXCOEFF,CXFAC,CXI,CXIKR,CXPHAS,CXZERO
      COMPLEX(WP) :: DCXSUM(6)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      SAVE CXZERO,CXI
      DATA CXI/(0._WP,1._WP)/,CXZERO/(0._WP,0._WP)/

!-----------------------------------------------------------------------
! subroutine DIRECT_CALCB_v5

! calculates magnetic field at (I,J,K) due to dipole at (1,1,1)
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

!        IPBC     = 0 if only doing first octant (IX=IY=IZ=1)
!                 = 1 otherwise
!                   N.B.: IPBC affects dimensioning of CXZG

!        GAMMA    = parameter controlling range of summations over
!                   replica dipoles.
!                   PHIMAX=1./gamma

!        NX,NY,NZ = size of first octant (I = 1 -> NX,
!                                         J = 1 -> NY,
!                                         K = 1 -> NZ)
!        PYDDX     = period of lattice in y direction/d
!        PZDDX     = period of lattice in z direction/d
!        DX(1-3)   = lattice spacing in x,y,z direction, in units of
!                    d = n**(-1/3).
!        AK(1-3)   = k(1-3)*d, where k = k vector in vacuo
!        AKD       = |k|d
!        AKD2      = |kd|^2
!        CXZG      = array with dimensions
!                    (NX+1)*(NY+1)*(NZ+1)*3 if IPBC=0
!                    (2*NX)*(2*NY)*(2*NZ)*3 if IPBC=1

! returns:
!        CXZG(I,J,K,M) = c_M to calculate magnetic field B at (I,J,K)
!                        contributed by dipole at (1,1,1) and
!                        replica dipoles (if IPBC=1)
!                        B_x =           c_1*P_y + c_2*P_z
!                        B_y = -c_1*P_x          + c_3*P_z
!                        B_z = -c_2*P_x -c_3*P_y

! note that even for short-range interaction, need to extend sums
! from JPY=-1 to JPY=+1 to include "edge" interactions with next TUC

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
! 12.07.07 (IW)  created DIRECT_CALCB from DIRECT_CALC
! 12.07.11 (BTD) v7.3 cosmetic changes to code written by Ian Wong
! 12.12.21 (BTD) bself_v3
!                * added comments
!                * corrected error when used for IPBC>0
! 13.12.26 (BTD) direct_calcb_v5 separated from bself
!                * new approach to summing over replica dipoles
!                  to ensure coverage of the Fresnel zone
! 13.12.27 (BTD) further modifications, paralleling direct_calc_v8.f90
! 13.12.29 (BTD) * set PHIMAX = 1/GAMMA so that PHIMAX can be
!                  modified from ddscat.par
! end history
! Copyright (C) 2012,2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!-----------------------------------------------------------------------

      PHIMAX=1./GAMMA

      BETA=16./PHIMAX**4

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
         IF(KPERP2<=0._WP)CALL ERRMSG('FATAL','DIRECT_CALCB',   &
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

! Compute 3 independent elements of 3x3 anti-symmetric matrix C_jk, where
! C_jk*P_k = magnetic field at location j due to dipole P_k at location

! C_jk = ( 0   c_1  c_2)
!        (-c_1  0   c_3)
!        (-c_2 -c_3  0 )_jk

! initialize CXG(I,J,K,M) = c_M for magnetic field at (I,J,K)
!                           produced by a dipole at (1,1,1)
!                           and replica dipoles (if PYD or PYZ are
!                           nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! X0,Y0,Z0 = X(I,J,K) - X(1,1,1) = vector from dipole location (1,1,1)
!                                  to point where B is to be calculated
! IX = +1 -> IMIN = 1
! IX = -1 -> IMIN = 2-NX   (1-IMIN = NX-1)
!
! similarly for JMIN and KMIN

      IMIN=1+(1-NX)*(1-IX)/2
      JMIN=1+(1-NY)*(1-IY)/2
      KMIN=1+(1-NZ)*(1-IZ)/2

! Determine elapsed cpu time for these sums

      CALL CPU_TIME(T0)

#ifdef openmp
! fork a team of threads 
!$OMP PARALLEL               &
!$OMP&   PRIVATE(NTHREADS,TID)

      TID=OMP_GET_THREAD_NUM()
!      WRITE(0,*)'direct_calcb_v5 ckpt 1: hello world from thread = ',TID

! only master thread does this:

      IF(TID.EQ.0)THEN
         NTHREADS=OMP_GET_NUM_THREADS()
         WRITE(CMSGNM,FMT='(A,I3)')'number of OpenMP threads =',NTHREADS
         CALL WRIMSG('DIRECT_CALCB',CMSGNM)
      ENDIF
!$OMP END PARALLEL
#endif

#ifdef openmp
! 2012.04.27 (BTD) added R to private variables
!                  removed FLUSH directive (should not be needed)
!$OMP PARALLEL                                        &
!$OMP&   PRIVATE(I,ICXZG,II)                          &
!$OMP&   PRIVATE(J,JCXZG,JPY,JPY1,JPY2,JPZ,JPZ1,JPZ2) &
!$OMP&   PRIVATE(K,KCXZG,M)                           &
!$OMP&   PRIVATE(DAMP,KR,KRF,PHASF,PHASY,PHASYZ,PHI)  &
!$OMP&   PRIVATE(R,R2,RF,RJPY,RJPZ,RPERP,TERM3)       &
!$OMP&   PRIVATE(X,X0,X2,X2Y2,XX)                     &
!$OMP&   PRIVATE(Y0,Y2,Y2Z2,YF,YHI,YLO)               &
!$OMP&   PRIVATE(Z0,Z2,ZF,ZHI,ZLO)                    &
!$OMP&   PRIVATE(CXCOEFF,CXFAC,CXIKR,CXPHAS,DCXSUM)
!$OMP DO
#endif

      DO K=KMIN,NZ              !- loop over K

         Z0=REAL(K-1,KIND=WP)*DX(3)
         Z2=Z0**2
         IF(K>0)THEN
            KCXZG=K
         ELSE
            KCXZG=2*NZ+K
         ENDIF
         DO J=JMIN,NY              !-- loop over J
            Y0=REAL(J-1,KIND=WP)*DX(2)
            Y2=Y0**2
            Y2Z2=Y2+Z2
            IF(J>0)THEN
               JCXZG=J
            ELSE
               JCXZG=2*NY+J
            ENDIF
            DO I=IMIN,NX              !--- loop over I

! for first dipole, determine time to sum over replicas

               IF(I.EQ.IMIN.AND.J.EQ.JMIN.AND.          &
                  K.EQ.KMIN.AND.MYID==0)CALL CPU_TIME(T1)

               X0=REAL(I-1,KIND=WP)*DX(1)
               X2=X0**2

! calculate JPY1, JPY2, JPZ1, JPZ2, 
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
                  ICXZG=I
               ELSE
                  ICXZG=2*NX+I
               ENDIF
               X(1)=X0
               DO M=1,3
                  DCXSUM(M)=CXZERO
               ENDDO

! JPY=0, JPZ=0 gives B field from dipole at (1,1,1)
! general JPY,JPZ gives B field from dipole at
! (1,JPY*NPY+1,JPZ*NPZ+1)
! replica dipoles have same magnitude of polarization as dipole (1,1,1),
! but different phase.
! PHASYZ = phase of replica dipole - phase of dipole (1,1,1)

               DO JPY=-JPY1,JPY2         !---- loop over JPY
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

                     IF(R2>1.E-6_WP)THEN

                        R=SQRT(R2)

! PHASYZ = phase at (1,JPY*NPY+1,JPZ*NPZ+1) - phase at (1,1,1)
! PHIF = phase shift at (I,J,K) from hypothetical dipole at center
!        of Fresnel zone
! PHI  = phase shift at (I,J,K) from replica of (1,1,1) relative to
!        phase shift from center of Fresnel zone
! include artificial factor exp[-DAMP] to assist convergence

                        PHASYZ=PHASY+AK(3)*RJPZ
                        KR=AKD*R
                        IF(JPBC==0)THEN
                           DAMP=0._WP
                        ELSE
                           PHI=KR-KRF+PHASYZ-PHASF
                           DAMP=BETA*PHI**4
                        ENDIF

!     k^2                    1
! B = --- * exp(ikr) * (1 - --- ) * r x p
!     r^2                   ikr 

! p = exp(i*phasyz) * p(1,1,1)

!     k^2                                    1
!   = --- * exp(ikr) * exp(i*phasyz) * (1 - --- ) * r x p(1,1,1)
!     r^2                                   ikr

! f(r) = (k/r)^2 * exp(ikr+i*phasyz) * [1 - 1/(ikr)]

! B = f(r) * r x p(1,1,1)     

! B_x =  c_1*p_y + c_2*p_z
! B_y = -c_1*p_x + c_3*p_z
! B_z = -c_2*p_x - c_3*p_y

! c_1 = -f * r_z
! c_2 =  f * r_y
! c_3 = -f * r_x

! include artificial factor exp[-DAMP] to assist convergence

                        CXIKR=CXI*KR
                        CXPHAS=EXP(CXIKR+CXI*PHASYZ-DAMP)
                        CXFAC=(1._WP-(1._WP/CXIKR))
			CXCOEFF=AKD2*CXFAC*CXPHAS/R2

                        DCXSUM(1)=DCXSUM(1)-CXCOEFF*X(3)
			DCXSUM(2)=DCXSUM(2)+CXCOEFF*X(2)
			DCXSUM(3)=DCXSUM(3)-CXCOEFF*X(1)
                      ENDIF   !----- endif (R2 > 1e-6)
                  ENDDO                     !----- end loop over JPZ
               ENDDO                     !---- end loop over JPY
               DO M=1,3
                  CXZG(ICXZG,JCXZG,KCXZG,M)=DCXSUM(M)
               ENDDO
               IF(I==IMIN.AND.J==JMIN.AND.K==KMIN.AND.MYID==0)THEN

! NB: this call to CPU_TIME occurs within the first thread
!     (I==IMIN, J==JMIN, K==KMIN)

                  CALL CPU_TIME(T2)
                  T2=REAL((NX-IMIN+1)*(NY-JMIN+1)*(NZ-KMIN+1),KIND=WP)*(T2-T1)
                  IF(T2<600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                          &
                        'Estimated total cputime required by DIRECT_CALCB=', &
                        T2,' cpu-sec'
                  ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                          &
                        'Estimated total cputime required by DIRECT_CALCB=', &
                        (T2/60._WP),' cpu-min'
                  ELSE
                     WRITE(CMSGNM,FMT='(A,F8.2,A)')                          &
                        'Estimated total cputime required by DIRECT_CALCB=', &
                        (T2/3600._WP),' cpu-hr'
                  ENDIF
                  CALL WRIMSG('DIRECT_CALCB',CMSGNM)
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
                  CALL WRIMSG('DIRECT_CALCB',CMSGNM)
               ENDIF
            ENDDO                     !--- end loop over I
         ENDDO                     !-- end loop over J

! 2012.04.27 (BTD) it is not necessary to have a flush(hxzc) here
!                  flush is implied by both $omp END DO and 
!                  $omp END PARALLEL directives

      ENDDO                     !- end loop over K

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

!*** diagnostic
!      write(0,*)'direct_calcb_v5 ckpt 2'
!***
      CALL CPU_TIME(T2)
!*** diagnostic
!      write(0,*)'direct_calcb_v5 ckpt 3'
!***
      T2=T2-T0
      IF(MYID==0)THEN
         IF(T2<600.D0)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                    &
                  'Actual cputime to complete DIRECT_CALCB=', &
                  T2,' cpu-sec'
         ELSEIF(T2>=600._WP.AND.T2<3600._WP)THEN
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                    &
                  'Actual cputime to complete DIRECT_CALCB=', &
                  (T2/60._WP),' cpu-min'
         ELSE
            WRITE(CMSGNM,FMT='(A,F8.2,A)')                    &
                  'Actual cputime to complete DIRECT_CALCB=', &
                  (T2/3600._WP),' cpu-hr'
         ENDIF
         CALL WRIMSG('DIRECT_CALCB',CMSGNM)
      ENDIF

!*** diagnostic
!      write(0,*)'direct_calcb_v5 ckpt 4'
!***

! If IPBC=1: set the elements with ICXZG=NX+1 or JCXZG=NY+1
!            or KCXZG=NZ+1 to zero

      IF(IPBC==1)THEN

         DO M=1,3

#ifdef openmp
!$OMP PARALLEL                    &
!$OMP&   PRIVATE(ICXZG,JCXZG,KCXZG)
!$OMP DO
#endif

            DO KCXZG=1,2*NZ
               DO JCXZG=1,2*NY
                  CXZG(NX+1,JCXZG,KCXZG,M)=CXZERO
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO KCXZG=1,2*NZ

               DO ICXZG=1,2*NX
                  CXZG(ICXZG,NY+1,KCXZG,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP DO
#endif

            DO JCXZG=1,2*NY

               DO ICXZG=1,2*NX
                  CXZG(ICXZG,JCXZG,NZ+1,M)=CXZERO
               ENDDO

            ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

         ENDDO   ! enddo M=1,3
 
      ENDIF   ! endif (IPBC==1)
      RETURN
 6700 FORMAT(F8.2,' cpu-sec')
 6701 FORMAT(F8.2,' cpu-min')
 6702 FORMAT(F8.2,' cpu-hr')
    END SUBROUTINE DIRECT_CALCB
