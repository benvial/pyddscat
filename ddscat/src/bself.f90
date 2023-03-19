    SUBROUTINE BSELF(CMETHD,CXZP,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK,AKD,DX, &
                     CXZG,CXZW,CXZB)
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_0,ONLY: AK2OLD_B,AK3OLD_B,NGRID,WOLD_B
      IMPLICIT NONE

!----------------------- bself v4 --------------------------------
! Arguments:

      CHARACTER(6) :: CMETHD
      INTEGER :: IPBC,NX,NY,NZ
      REAL(WP) :: AKD,GAMMA,PYD,PZD
      REAL(WP) :: AK(3),DX(3)
      COMPLEX(WP) ::                                                 &
         CXZB(NX,NY,NZ,3),                                           &
         CXZG(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),3), &
         CXZP(NX,NY,NZ,3),                                           &
         CXZW(2*NX,2*NY,2*NZ,*)

! NB: module DDCOMMON_0 must have previously set values of
!       AK2OLD_B,AK3OLD_B,WOLD_B
!    to be used by BSELF

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: I,IR,ISGN,J,JR,JX,JY,JZ,JSGN,K,KR,KSGN,M
      REAL(WP) :: AKD2,DTIME,FAC,PYDDX,PZDDX
      COMPLEX(WP) :: CXBX,CXBY,CXBZ,CXXY,CXXZ,CXYZ

! Local arrays:

      INTEGER :: ISYM(3)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      EXTERNAL CXFFTW,DIRECT_CALCB,EXTND,PAD,TIMEIT,TRIM
      INTRINSIC MIN,NINT,SIGN
 
!-----------------------------------------------------------------------
! Parameter GAMMA determines the range of the sums when periodic
! boundary conditions are employed.  Dipole-dipole interactions
! are screened by a factor exp[-(gamma*kr)^4]
! The effective
! range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
! range/lambda = 1/(2*pi*gamma) = 31.8
 
! The sums are actually continued out to
! r/d = 2*/(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
! [screening factor = exp(-16)=1.1e-7]

!-----------------------------------------------------------------------
! subroutine BSELF

!       Given the dipole moments, CXZP, at all points on
!       a rectangular grid, oscillating at frequency AKD,
!       compute the magnetic field amplitude, CXZB,
!       at each point produced by all the other dipoles except the one
!       at that point.
!       The relationship between the dipoles and the field
!       values at the grid points can be expressed as a convolution,
!       and the convolution is efficiently evaluated using 3D FFTs.
!
!       options for computation of 3-dimensional FFT:
!
!    if CMETHD='GPFAFT':
!          Use CXFFT3N interface to GPFA code of Temperton.
!          Good points:
!             -On CRAY the code is on average 30-40% faster in comparison
!              to TMPRTN (and 10-15 faster in comparison to BRENNR)
!             -On scalar machines the code is 2-8 faster in comparison
!              to BRENNR or TMPRTN
!             -Does not require additional storage
!          Limitations:
!              -Requires that NX,NY,NZ be of form (2**I)*(3**J)*(5**K);
!               subroutine EXTEND takes care of choosing suitable NX,NY,
!              -The choice of "lvr" variable in gpfa2f, gpfa3f, gpfa5f
!               depends on machine. WARNING: on C90 use lvr=128, on
!               all other CRAYs use lvr=64; wrong lvr will produce WRONG
!               RESULTS. On scalar machines optimal lvr depends on cache
!               length; sub-optimal choice degrades performance but still
!               produces correct results.
!
!   if CMETHD='FFTW21'
!      Use CXFFTW interface to FFTW (Fastest Fourier Transform in the
!               West) version 2.1.x from Frigo and Johnson.
!
!   if CMETHD='FFTMKL':
!      Use CXFFT3_MKL interface to Intel Math Kernel Library (MKL) FFT
!
! Input:
!
!    CXZP(I,J,K,L)   Lth cartesian component of the dipole
!                    moment at the grid point (I,J,K);
!                    the DIMENSIONed length of CXZP in
!                    the calling routine is CXZP(NX,NY,NZ,3)
!                    [or CXZP(NX*NY*NZ,3) or CXZP(3*NX*NY*NZ)]

!    NX,NY,NZ        Size of grid in x,y,z directions (INTEGER).

!    IPBC          = 0 for isolated target
!                    1 for periodic target

!    GAMMA         = coefficient used to assist convergence of sums
!                    over replica dipoles by suppressing long-range
!                    contributions with factor exp(-gamma*(kr)^2)
!                    typical value gamma = 0.005
!                    The effective
!                    range/d = 1/(gamma*k*d) = 400 if gamma=5e-3 and kd=0.5
!                    range/lambda = 1/(2*pi*gamma) = 31.8
!                    The sums are actually continued out to
!                    r/d = 2*/(gamma*kd) = 800 if gamma=5e-3 , kd=0.5
!                    [screening factor = exp(-16)=1.1e-7]

!       PYD          (Period of lattice in y direction)/DX(2)
!       PZD          (Period of lattice in z direction)/DX(3)

!       DX(1-3)      Lattice spacing in x,y,z directions, in units of
!                    n**(-1./3.) .  Note that with this normalization
!                       we have DX(1)*DX(2)*DX(3)=1.

!       AK(1-3)         k(1-3)*d, where k = k vector in vacuo, and
!                       d = effective lattice spacing = (dx*dy*dz)**(1/3)

!       AKD             = (omega/c)*d = k*d (dimensionless)

!       CXZG            (NX+1)*(NY+1)*(NZ+1)*3 array of Green
!                       function coefficients used
!                       internally by BSELF and
!                       *************NOT TO BE OVERWRITTEN***********
!                       between calls, because these coefficients are
!                       recomputed only if W has changed since the last
!                       call to BSELF.

!       CXZW            Complex, scratch-space vector of length:
!                       2*NX*2*NY*2*NZ*3
!                       See comment about FFT usage and CMETHD flag.
!                       Can be overwritten between calls to BSELF

! OUTPUT:

!       CXZB(I,J,K,L)   Lth component of dipole-generated magnetic field
!                       at grid point (I,J,K);
!                       the DECLARED length of ZB in the calling
!                       program is CXZB(NX,NY,NZ,3)
!                       [or CXZB(NX*NY*NZ,3) or CXZB(3*NX*NY*NZ)]
!
!=======================================================================
! subroutine BSELF
! history
! based on subroutine ESELF written originally by Jeremy Goodman
! BSELF created by Ian Wong, Princeton University, July 2012
! 12.07.07 (IW)  v1 written
! 12.07.11 (BTD) v2 created from v1
!                * corrected typos (HXZB -> CXZB)
!                * changed notation
!                  HXZC -> CXZG  (Green function)
!                  HXZW -> CXZW
!                  HXZB -> CXZB
!                * a few changes to comments
! 12.12.21 (BTD) v3
!                * added comments
!                * corrected error for periodic targets
! 12.12.22 (BTD) * revised comments, minor cleanup
! 13.01.03 (BTD) v4
!                * added AK2OLD_B,AK3OLD_B,WOLD_B from DDCOMMON_0
!                * modified to skip recomputation of Green-function
!                  coefficients on second call to BSELF
!                  NB: this requires that CXZG *not* be deallocated
!                  in subroutine NEARFIELD after first call to BSELF
! end history
! Copyright (C) 2013 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!=======================================================================        

! check if we can skip recomputation of Green-function coefficients

      IF(PYD.EQ.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD_B)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.EQ.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(3)-AK3OLD_B)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD_B-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD_B)<1.E-6_WP*AKD.AND. &
            ABS(AK(3)-AK3OLD_B)<1.E-6_WP*AKD)GOTO 70
      ENDIF

! Compute Green function coefficients 

! We have to compute the Green-function coefficients giving
! components of magnetic field strength at a each grid point R
! produced by unit-valued component of dipole moment at
! point Rprime, and then Fourier transform these components.

      WOLD_B=AKD
      AK2OLD_B=AK(2)
      AK3OLD_B=AK(3)
      NGRID=8*NX*NY*NZ
      AKD2=AKD*AKD

! We assume screening function exp(-(gamma*kr)^4) so
! range/d = 1/(gamma*kd) = 2000 if gamma=1e-3 and kd=0.5
! although the sums are actually continued out to
! r/d = 2/(gamma*kd) = 4000 if gamma=1e-3, kd=0.5
! [screening factor = exp(-16)=1.1e-7]

! PYD*DX(2) = periodicity in Y direction
! PZD*DX(3) = periodicity in Z direction

      IF(PYD>0._WP.OR.PZD>0._WP)THEN
         WRITE(CMSGNM,FMT='(A,2F8.2,A,1PE9.2)')'PBC with PYD, PZD=',PYD, &
                                               PZD,', GAMMA=',GAMMA
         CALL WRIMSG('BSELF',CMSGNM)
      ENDIF

      PYDDX=PYD*DX(2)
      PZDDX=PZD*DX(3)


! Compute 3 independent elements of 3x3 anti-symmetric matrix C_jk,where
! C_jk*P_k = magnetic field at location j due to dipole P_k at location

! C_jk = ( 0    c_1  c_2)
!        (-c_1  0    c_3)
!        (-c_2 -c_3   0 )_jk

      IF(IPBC==0)THEN

! initialize CXZG(I,J,K,M) = c_M for magnetic field at (I,J,K)
!                            produced by a dipole at (1,1,1)
!                            and replica dipoles (if PYD or PYZ are
!                            nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! Later obtain C_jk values for other octants by using symmetry

!*** diagnostic
!         write(0,*)'bself_v4 ckpt 1: call DIRECT_CALCB'
!***
        CALL DIRECT_CALCB(1,1,1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZG(1,1,1,1))
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 2'
!        write(0,*)'   returned from direct_calcb'
!        write(0,*)'   check for NaN...'
!        jr=0
!        do i=1,nx
!           do j=1,ny
!              do k=1,nz
!                 do m=1,3
!                    if(.not.(abs(hxzc(i,j,k,m))>=0.d0).or. &
!                        abs(hxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'hxzc=',hxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 3: cxy checked for NaN or overflow',jr, &
!                   ' instances found'
!***

! At this point, CXZG(I,J,K,1-3) contains the upper triangular part of the
! anti-symmetric 3 x 3 matrix giving the magnetic field at grid point (i,j,k)
! produced by a dipole at (1,1,1)

! Fill out CXZG to twice the size in each grid dimension to accomodate
! negative lags [periodicity in each dimension is assumed, so (e.g.)
! nx < i <= 2*nx is equivalent to -nx < i <= 0], exploiting symmetries,
! and then Fourier transform.

! If PYDDX=0 and PZDDX=0 , need only do direct calculation of C matrix
! for first octant, since remaining octants can be obtained by symmetry.
! After calculating C matrix, store only the first octant of the
! transform, since the others can be obtained by symmetry.

!-----------------------------------------------------------------------
! extend c_1 = c_xy(x,y,z) : c -> +c for x -> -x
!                                 +c     y -> -y
!                                 -c     z -> -z

        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=-1
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 4, about to call EXTND'
!***
        CALL EXTND(CXZG(1,1,1,1),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'bself_v4 ckpt 5, about to call cxfft3n'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!          write(0,*)'bself_v4 ckpt 6, returned from cxfft3n'
!***
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 7, about to call TRIM'
!***
        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZG(1,1,1,1))
!*** diagnostic
!        write(0,*)'returned from TRIM'
!***
!-----------------------------------------------------------------------
! extend c_2 = c_xz(x,y,z) : c -> +c for x -> -x
!                                 -c     y -> -y
!                                 +c     z -> -z

        ISYM(1)=1
        ISYM(2)=-1
        ISYM(3)=1
        CALL EXTND(CXZG(1,1,1,2),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZG(1,1,1,2))

!-----------------------------------------------------------------------
! extend c_3 = c_yz(x,y,z) : c -> -c for x -> -x
!                                 +c     y -> -y
!                                 +c     z -> -z
        ISYM(1)=-1
        ISYM(2)=1
        ISYM(3)=1
        CALL EXTND(CXZG(1,1,1,3),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZG(1,1,1,3))


      ELSEIF(IPBC==1)THEN

! This point is reached when PYDDX or PZDDX are nonzero.
! When PBC are used for general direction of incident wave,
! all octants of C matrix require direct calculation: symmetries valid
! for single target no longer apply because of position-dependent phases
! of replica dipoles.

! DIRECT_CALCB computes arrays (c_1)_jk , (c_2)_jk  (c_3)_jk
!
!               (   0   c_1  c_2 )
! where B(r_j)= ( -c_1   0   c_3 ) * P(r_k)
!               ( -c_2 -c_3   0  )
!
! CXZG(I,J,K,M) = c_M for (r_j - r_k)/d = (I-1)*xhat + (J-1)*yhat + (K-1)*zhat
!
! when IPBC=1, CXZG includes contribution to B from replica dipoles

!*** diagnostic
!         write(0,fmt='(a,a,4I4)')'bself_v4 ckpt 7.5',    &
!              ' call direct_calcb with NX,NY,NZ=',NX,NY,NZ
!***

        CALL DIRECT_CALCB(-1,-1,-1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZG(1,1,1,1))

!*** diagnostic
!        write(0,*)'bself_v4 ckpt 8'
!***
! The array CXZG(1-2*NX,1-2*NY,1-2*NZ,1-3) of C matrix coefficients
! now covers all octants.

! Fourier transform the C matrix CXZG:

        DO M=1,3
           IF(CMETHD=='GPFAFT')THEN
              CALL CXFFT3N(CXZG(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD=='FFTW21')THEN
              CALL CXFFTW(CXZG(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD.EQ.'FFTMKL')THEN
              CALL CXFFT3_MKL(CXZG(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ENDIF
        ENDDO

!*** diagnostic
!        write(0,*)'bself_v4 ckpt 9'
!***

! CXZG now contains the full Fourier transform of the C convolution
! and should not be overwritten between calls to BSELF

      ENDIF
!      CALL TIMEIT('BSELF (first call)',DTIME)
!      CALL TIMEIT('BSELF',DTIME)
!-----------------------------------------------------------------------

! End of computation of Green-function coefficients

70    CONTINUE

!*** diagnostic
!      write(0,*)'bself_v4 ckpt 10'
!****
! Fourier transform the polarizations

    DO M=1,3
       CALL PAD(CXZP(1,1,1,M),NX,NY,NZ,CXZW(1,1,1,M))
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 11: returned from PAD: ', &
!                  'check cxzw for NaN or overflow...'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 if(.not.(abs(cxzw(i,j,k,m))>=0.d0).or. &
!                     abs(cxzw(i,j,k,m))>=1.d100)then
!                    write(0,*)'i,j,k,m=',i,j,k,m,'cxzw=',cxzw(i,j,k,m)
!                    jr=jr+1
!                 endif
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 12: cxzw checked for NaN or overflow: ', &
!                  jr,' instances found'
!        if(jr>0)stop
!***

        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!        write(0,*)'bself_v4 ckpt 13: returned from CXFFT3N: ', &
!                  'check cxzw for NaN or overflow...'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 if(.not.(abs(cxzw(i,j,k,m))>=0.d0).or. &
!                     abs(cxzw(i,j,k,m))>=1.d100)then
!                    write(0,*)'i,j,k,m=',i,j,k,m,'cxzw=',cxzw(i,j,k,m)
!                    jr=jr+1
!                 endif
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 14: cxzw checked for NaN or overflow',jr, &
!                  ' instances found'
!        if(jr>0)stop
!***
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ENDIF
      ENDDO

!***********************************************************************

! Multiply by F.t. of Green function.

      IF(IPBC==0)THEN

!*** diagnostic
!         write(0,*)'bself_v4 ckpt 15, check hxzc(i,j,k,m=1,6)'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 do m=1,3
!                    if(.not.(abs(hxzc(i,j,k,m))>=0.d0).or. &
!                        abs(hxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'hxzc=',hxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 16: hxzc(i,j,k,m=1-3) checked for NaN ', &
!                  'or overflow: ',jr,' instances found'
!        if(jr>0)stop
!        write(0,*)'bself_v4 ckpt 17: check cxzw(i,j,k,m=1-3)'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 do m=1,3
!                    if(.not.(abs(cxzw(i,j,k,m))>=0.d0).or. &
!                        abs(cxzw(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'cxzw=',cxzw(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'bself_v4 ckpt 18: cxzw(i,j,k,m=1-3) checked for NaN or ', &
!                  'overflow: ',jr,' instances found'
!        if(jr>0)stop
!***
!***

! If IPBC=0, then only one octant of F.t. of Green function has been
!            stored, but can recover others using symmetry.

#ifdef openmp
!$OMP PARALLEL DO                                            &
!$OMP&   PRIVATE(K,J,I,KSGN,KR,JSGN,JR,ISGN,IR)              &
!$OMP&   PRIVATE(CXXY,CXXZ,CXYZ,CXBX,CXBY,CXBZ)
#endif

        DO K=1,2*NZ
           KSGN=NINT(SIGN(1._WP,NZ+1.5_WP-K))
           KR=MIN(K,2*NZ+2-K)
!*** diagnostic
!          write(0,*)'K,KSGN,KR=',K,KSGN,KR
!***
           DO J=1,2*NY
              JSGN=NINT(SIGN(1._WP,NY+1.5_WP-J))
              JR=MIN(J,2*NY+2-J)
!*** diagnostic
!            write(0,*)'   J,JSGN,JR=',J,JSGN,JR
!***
              DO I=1,2*NX
                 ISGN=NINT(SIGN(1._WP,NX+1.5_WP-I))
                 IR=MIN(I,2*NX+2-I)
!*** diagnostic
!              write(0,*)'       I,ISGN,IR=',I,ISGN,IR
!***
                 CXXY=CXZG(IR,JR,KR,1)*(ISGN*JSGN)
                 CXXZ=CXZG(IR,JR,KR,2)*(ISGN*KSGN)
                 CXYZ=CXZG(IR,JR,KR,3)*(JSGN*KSGN)
                 
!*** diagnostic
!              if(.not.(abs(cxxy)>=0.d0))write(0,*) &
!                'ir,jr,kr,hxzc(ir,jr,kr,1)=',ir,jr,kr,hxzc(ir,jr,kr,1)
!              if(.not.(abs(cxxz)>=0.d0))write(0,*) &
!                 'ir,jr,kr,hxzc(ir,jr,kr,2)=',ir,jr,kr,hxzc(ir,jr,kr,2)
!              if(.not.(abs(cxyz)>=0.d0))write(0,*) &
!                  'ir,jr,kr,hxzc(ir,jr,kr,3)=',ir,jr,kr,hxzc(ir,jr,kr,3)

!***
                 CXBX=CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                 CXBY=-CXXY*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,3)
                 CXBZ=-CXXZ*CXZW(I,J,K,1)-CXYZ*CXZW(I,J,K,2)
                      
!*** diagnostic
!              if(.not.(abs(cxbx+cxby+cxbz)>=0.d0))then
!                 write(0,*)'i,j,j,ir,jr,kr=',i,j,k,ir,jr,kr
!                 write(0,*)'hxzc(ir,jr,kr,1)=',hxzc(ir,jr,kr,1)
!                 write(0,*)'hxzc(ir,jr,kr,2)=',hxzc(ir,jr,kr,2)
!                 write(0,*)'hxzc(ir,jr,kr,3)=',hxzc(ir,jr,kr,3)
!                 write(0,*)'      cxbx=',cxbx
!                 write(0,*)'      cxby=',cxby
!                 write(0,*)'      cxbz=',cxbz
!                 write(0,*)'      cxxy=',cxxy
!                 write(0,*)'      cxxz=',cxxz
!                 write(0,*)'      cxyz=',cxyz
!                 write(0,*)'      cxzw(i,j,k,1)=',cxzw(i,j,k,1)
!                 write(0,*)'      cxzw(i,j,k,2)=',cxzw(i,j,k,2)
!                 write(0,*)'      cxzw(i,j,k,3)=',cxzw(i,j,k,3)
!                 stop
!              endif
!***

! now overwrite CXZW [which contained the Fourier transform of C matrix]
!                     with F.t. of B

                 CXZW(I,J,K,1)=CXBX
                 CXZW(I,J,K,2)=CXBY
                 CXZW(I,J,K,3)=CXBZ

! CXZW is now the Fourier transform of B

!*** diagnostic
!              if(.not.(abs(cxzw(i,j,k,1))>=0.d0))write(0,*) &
!                 'i,j,k,cxzw(i,j,k,1)=',i,j,k,cxzw(i,j,k,1)
!              if(.not.(abs(cxzw(i,j,k,2))>=0.d0))write(0,*) &
!                 'i,j,k,cxzw(i,j,k,2)=',i,j,k,cxzw(i,j,k,2)
!              if(.not.(abs(cxzw(i,j,k,3))>=0.d0))write(0,*) &
!                 'i,j,k,cxzw(i,j,k,3)=',i,j,k,cxzw(i,j,k,3)
!***
               ENDDO
            ENDDO
         ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ELSEIF(IPBC==1)THEN

! If IPBC=1, then the full F.t. of the Green function has been stored.

#ifdef openmp
!$OMP PARALLEL DO                                                  &
!$OMP&   PRIVATE(K,J,I,CXXY,CXXZ,CXYZ,CXBX,CXBY,CXBZ)
#endif

         DO K=1,2*NZ
            DO J=1,2*NY
               DO I=1,2*NX
                  CXXY=CXZG(I,J,K,1)
                  CXXZ=CXZG(I,J,K,2)
                  CXYZ=CXZG(I,J,K,3)
                  CXBX=CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                  CXBY=-CXXY*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,3)
                  CXBZ=-CXXZ*CXZW(I,J,K,1)-CXYZ*CXZW(I,J,K,2)

! now overwrite CXZW [which contained the F.t. of C matrix]
!                     with Fourier transform of B

                  CXZW(I,J,K,1)=CXBX
                  CXZW(I,J,K,2)=CXBY
                  CXZW(I,J,K,3)=CXBZ

! CXZW is now the Fourier transform of B

               ENDDO
            ENDDO
         ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ENDIF

! Inverse Fourier transform to obtain magnetic field:

      DO M=1,3
         IF(CMETHD=='GPFAFT')THEN
            CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ENDIF

!***********************************************************************

! Note: the Convex FFT routine already normalizes result.
!       For other FFT routines need to divide result by NGRID

         IF(CMETHD=='CONVEX')THEN
            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                     CXZB(I,J,K,M)=CXZW(I,J,K,M)
                  ENDDO
               ENDDO
            ENDDO
         ELSE
            FAC=1._WP/REAL(NGRID,KIND=WP)

#ifdef openmp
!$OMP PARALLEL DO     &
!$OMP&   PRIVATE(I,J,K)
#endif

            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                     CXZB(I,J,K,M)=FAC*CXZW(I,J,K,M)
                  ENDDO
               ENDDO
            ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

         ENDIF
      ENDDO
      RETURN
    END SUBROUTINE BSELF
