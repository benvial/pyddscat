    SUBROUTINE ESELF(CMETHD,CXZP,NX,NY,NZ,IPBC,GAMMA,PYD,PZD,AK,AKD, &
                     DX,CXZC,CXZW,CXZE)
      USE DDPRECISION,ONLY: WP
      USE DDCOMMON_0,ONLY: AK2OLD,AK3OLD,IDIPINT,NGRID,WOLD
      IMPLICIT NONE

!----------------------- eself v8 --------------------------------
! Arguments:

      CHARACTER(6) :: CMETHD
      INTEGER :: IPBC,NX,NY,NZ
      REAL(WP) :: AKD,GAMMA,PYD,PZD
      REAL(WP) :: AK(3),DX(3)
      COMPLEX(WP) ::                                                 &
         CXZC(NX+1+IPBC*(NX-1),NY+1+IPBC*(NY-1),NZ+1+IPBC*(NZ-1),6), &
         CXZE(NX,NY,NZ,3),                                           &
         CXZP(NX,NY,NZ,3),                                           &
         CXZW(2*NX,2*NY,2*NZ,*)

! NB: module DDCOMMON_0 must have previously set values of
!       AK2OLD,AK3OLD,WOLD
!    to be used by ESELF

! Local scalars:

      CHARACTER :: CMSGNM*70
      INTEGER :: I,IR,ISGN,J,JR,JX,JY,JZ,JSGN,K,KR,KSGN,M
      REAL(WP) :: AKD2,DTIME,FAC,PYDDX,PZDDX
      COMPLEX(WP) :: CXEX,CXEY,CXEZ,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ

! Local arrays:

      INTEGER :: ISYM(3)

#ifdef openmp
      INTEGER NTHREADS,TID
#endif

      EXTERNAL CXFFTW,DIRECT_CALC,EXTND,PAD,TIMEIT,TRIM
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
!
! 
!-----------------------------------------------------------------------
! subroutine ESELF

!       Given the dipole moments, CXZP, at all points on
!       a rectangular grid, oscillating at frequency AKD,
!       compute the electric field amplitude, CXZE,
!       at each point produced by all the other dipoles except the one
!       at that point.
!       The relationship between the dipoles and the field
!       values at the grid points can be expressed as a convolution,
!       and the convolution is efficiently evaluated using 3D FFTs.

!       options for computation of 3-dimensional FFT:

!    if CMETHD='GPFAFT':
!          Use CXFFT3N interface to GPFA code of Temperton.
!          Good points:
!             -On scalar machines the code is 2-8 faster in comparison
!              to BRENNR, and comparable in speed to FFTW 2.1
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
!   if CMETHD='FFTW21'
!      Use CXFFTW interface to FFTW (Fastest Fourier Transform in the
!               West) version 2.1.x from Frigo and Johnson.
!
!   if CMETHD='FFTMKL':
!      Use CXFFT3_MKL interface to Intel Math Kernel Library (MKL) FFT
!
! INPUT:

!       CXZP(I,J,K,L)   Lth cartesian component of the dipole
!                       moment at the grid point (I,J,K);
!                       the DIMENSIONed length of CXZP in
!                       the calling routine is CXZP(NX,NY,NZ,3)
!                       [or CXZP(NX*NY*NZ,3) or CXZP(3*NX*NY*NZ)]

!       NX,NY,NZ        Size of grid in x,y,z directions (INTEGER).

!       IPBC          = 0 for isolated target
!                     = 1 to use periodic boundary conditions
!                         (in either y direction, z direction, or both)
!       PYD             (Period of lattice in y direction)/DX(2)
!       PZD             (Period of lattice in z direction)/DX(3)

!       GAMMA         = coefficient used to assist convergence in
!                       lattice sums by suppressing long-range contributions
!                       with factor exp(-(gamma*k*r)^4)

!       DX(1-3)         Lattice spacing in x,y,z directions, in units of
!                       n**(-1./3.) .  Note that with this normalization
!                       we have DX(1)*DX(2)*DX(3)=1.

!       AK(1-3)         k(1-3)*d, where k = k vector in vacuo, and
!                       d = effective lattice spacing = (dx*dy*dz)**(1/3)

!       AKD           = (omega/c)*d = k*d (dimensionless)

!       CXZC            (NX+1)*(NY+1)*(NZ+1)*6 array of Green
!                       function coefficients used
!                       internally by ESELF and
!                       *************NOT TO BE OVERWRITTEN***********
!                       between calls, because these coefficients are
!                       recomputed only if W has changed since the last
!                       call to ESELF.

!       CXZW            Complex, scratch-space vector of length:
!                       2*NX*2*NY*2*NZ*3
!                       See comment about FFT usage and CMETHD flag.
!                       Can be overwritten between calls to ESELF
!
! OUTPUT:

!       CXZE(I,J,K,L)   Lth component of dipole-generated electric field
!                       at grid point (I,J,K);
!                       the DECLARED length of CXZE in the calling
!                       program is CXZE(NX,NY,NZ,3)
!                       [or CXZE(NX*NY*NZ,3) or CXZE(3*NX*NY*NZ)]

! Originally written by Jeremy Goodman, 
! Princeton Univ. Observatory, 90.09.22
! History:
! 90.11.29 (BTD): Modified to set untransformed ZC(1,1,1,1-6)=0.
! 90.11.29 (PJF): Modified to use FOURX and CXFFT99
! 90.12.05 (BTD): Modified for new ordering of elements of polarization
!                 and electric field vectors in calling program.  
!                 Modified ESELF and PAD to remove distinction between 
!                 NX,NY,NZ and dimensions of CXZE and CXZP, 
!                 since our new ordering always assumes this.
! 90.12.13 (BTD): Modified to include CONVEX option.
! 92.04.20 (BTD): removed ISYM from argument list of TRIM (was not used)
! 94.06.20 (PJF): modified to call CXFFT3N when CMETHD='NEWTMP'
! 96.10.18 (BTD): changed NEWTMP to GPFAFT
! 97.10.16 (BTD): added DX(1-3) to argument list to support use for
!                 noncubic rectangular lattices.
!                 Added DX to 3 lines computing X(1-3)
! 99.04.26 (BTD): changed notation: CXY -> CXZE
! 00.06.22 (BTD): modified to support option FFTWFJ, with calls to CXFFT
! 00.06.25 (BTD): modified to eliminate calls to FOURX (CMETHD=BRENNR)
!                 and CXFFT3 (CMETHD=TMPRTN), as these options are no
!                 longer supported
! 00.07.05 (BTD): further cleanup
! 03.07.13 (BTD): changed FFTWFJ to FFTW21 to allow future distinction
!                 between FFTW 2.1.x and FFTW 3.0.x
! 04.03.05 (BTD): modified to handle a periodic lattice of scatterers
!                 with periodicity NPY*D(2) in y direction and
!                 NPZ*D(3) in z direction
!                 Add parameter BETA to assist convergence
! 05.06.16 (BTD): Replaced integer NPY,NPZ by real PYD,PZD in
!                 argument list and in calculation
!                 corrected error in calculation of JPZM
! 05.07.08 (BTD): corrected error in calculation of JPZM
! 05.08.01 (BTD): changes to increase efficiency of summations
!                 required for PBC option (new variables AKD2,
!                 PYDDX,PZDDX,X0,Y0,Z0,X2,Y2,X2Y2,CXTERM)
! 05.08.03 (BTD): corrected typo.
! 05.08.04 (BTD): added AK(1-3) to argument list,
!                 added local variables PHASY and PHASYZ
!                 added phase shift exp(i*PHASYZ) for replica dipole
! 06.09.15 (BTD): added comments, reduced BETA to 1e-12 for improved
!                 accuracy (at expense of speed).
! 06.09.23 (BTD): corrected error: sign error in calculation of
!                 diplacements X(2) and X(3) for replica dipoles
! 06.09.23 (BTD): corrected error in computation of factor CXFAC
!                 appearing in summation for CXSUM
! 06.09.28 (BTD): eself v2.0 and DDSCAT 6.2.3:
!                 * put part of calculation of A matrix into subroutine
!                   DIRECT_CALC
!                 * modified to support PBC option
!                 * added IPBC to argument list to support PBC option
!                 * changed dimensioning of CXZC when IPBC=1
!                 * when IPBC=1, store full CXZC rather than only first
!                   octant
! 07.08.06 (BTD): Version 7.0.3
!                 eliminate option CONVEX -> CXC3DFFT
!                 No longer appear to be any sites where this would be
!                 useful.
! 07.10.04 (BTD): Recompute Green function coefficients only if
!                 frequency changes by more than 1 part in 1e6
!                 (previous condition of requiring frequency to be
!                 unchanged was evidently being confused by roundoff
!                 error, leading to unnecessary recalculation).
! 07.10.25 (BTD): Modified so that for IPBC>0, Green functions are
!                 recalculated when direction of incidence is
!                 changed.
! 08.03.11 (BTD): v7.0.5
!                 * added ALPHA to argument list
!                 * eliminated variable BETA
!                 * replaced exp(-beta*(k*r)^4) with exp(-(alpha*k*r)^4)
! 08.04.20 (BTD): * change notation: ALPHA -> GAMMA
! 08.05.12 (BTD): v7.0.6
!                 * added OpenMP directives as suggested by Art Lazanoff, eg.
!                      #ifdef eself_omp
!                      !$omp parallel do
!                      #endif
!                   which will be compiled when the preprocessor flag
!                      -fpp -Deself_omp
!                   is used.
! 08.06.05 (ASL,BTD): eself_v3
!                 * added call to new routine CXFFT3_MKL to use Intel MKL
!                   library for FFT.
!                   This will now be invoked with CMETHD='FFTMKL'
! 08.06.27 (BTD) : eself_v4
!                 * BTD learning the ropes...
! 08.08.07 (BTD): removed #ifdef debug / #endif statements inserted
!                 around !omp directives by Art Lazanoff for testing.
!                 omp directives appear to work properly, they can
!                 continued to be tested as a whole by compiling with or 
!                 without the -openmp flag, and it is desirable for
!                 export code to be able to compile without cpp insofar
!                 as possible
! 11.07.03 (PFJ,BTD): eself_v6
!                 * removed local initialization of AK2OLD,AK3OLD,WOLD
!                   so that these values are now passed from calling
!                   program through module DDCOMMON_0
!                 * NGRID is now passed through module DDCOMMON_0 so
!                   that NGRID may be provided to other routines
! 12.07.06 (BTD): edited comments
! 12.08.02 (IYW): v7.3 : eself_v7
!                 * added DIPINT to arg list
!		  * introduced possibility of calculating Green func 
!                   coefficients using "Filtered Coupled Dipole" (FCD) method
!                   (cf. ADDA)
!		  * added subroutine CISI for calculation of sine and cosine 
!                   integrals
! 12.08.11 (BTD): * removed DIPINT from arg list of ESELF
!                 * removed DIPINT from arg list of DIRECT_CALC
!                 * added new variable IDIPINT from module DDCOMMON_0
! 13.12.23 (BTD): v8 
!                 * removed DIRECT_CALC to direct_calc_v8.f90
! 14.10.30 (BTD): * removed PAD to pad.f90 so that it can be compiled
!                   with a different level of optimization
!                   need to do this because ifort -O2 was producing bad 
!                   code on mac powerbook ourania4 (ifort version 14.0.1
!                   under osx 10.9.5)
! end history
! Copyright (C) 1993,1994,1996,1997,1999,2000,2003,2004,2005,2006,2007,
!               2008,2011,2012,2013,2014 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

!--------------------------------------------

!*** diagnostic
!      write(0,*)'eself_v8 ckpt 0 with'
!      write(0,*)'       cmethd=',cmethd
!      write(0,*)'     nx,ny,nz=',nx,ny,nz
!      write(0,*)'      idipint=',idipint
!      write(0,*)'         ipbc=',ipbc
!      write(0,*)'        gamma=',gamma
!      write(0,*)'      pyd,pzd=',pyd,pzd
!      write(0,*)'      ak(1-3)=',ak
!      write(0,*)'          akd=',akd
!      write(0,*)'      dx(1-3)=',dx
!      write(0,*)'         wold=',wold
!      write(0,*)'       ak2old=',ak2old
!      write(0,*)'       ak3old=',ak3old
!      write(0,*)'  j jx jy jz k   cxzp(jx,jy,jz,k)'
!      do k=1,3
!         do jz=1,nz
!            do jy=1,ny
!               do jx=1,nx
!                  i=nz*ny*nx*(k-1)+ny*nx*(jz-1)+nx*(jy-1)+jx
!                  write(0,fmt='(i4,i3,i3,i3,i2,1p2e11.3)') &
!                               i,jx,jy,jz,k,cxzp(jx,jy,jz,k)
!               enddo
!            enddo
!         enddo
!      enddo
!***

! 14.10.16 (BTD) sanity check to find problem with SPHRN_PBC...
!*** diagnostic
!      do k=1,3
!         do jz=1,nz
!            do jy=1,nx
!               if(.not.(abs(cxzp(jx,jy,jz,k)).lt.1.e32))then
!                  write(0,*)'eself_v8 ckpt 1: *** problem ***'
!                  write(0,*)'  jx,jy,jz,k=',jx,jy,jz,k
!                  write(0,*)'  cxzp(jx,jy,jz,k)=',cxzp(jx,jy,jz,k)
!               endif
!            enddo
!         enddo
!      enddo
!***

! check if we can skip recomputation of Green-function coefficients

      IF(PYD.EQ.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.EQ.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.EQ.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(3)-AK3OLD)<1.E-6_WP*AKD)GOTO 70
      ELSEIF(PYD.NE.0._WP.AND.PZD.NE.0._WP)THEN
         IF(ABS(WOLD-AKD)<1.E-6_WP*AKD.AND.   &
            ABS(AK(2)-AK2OLD)<1.E-6_WP*AKD.AND. &
            ABS(AK(3)-AK3OLD)<1.E-6_WP*AKD)GOTO 70
      ENDIF

!*** diagnostic
!      write(0,*)'eself_v8 ckpt 2: recompute Green-function coefficients'
!***
! AKD.NE.WOLD :

! We have to recompute the Green-function coefficients giving
! components of field strength at a each grid point R
! produced by unit-valued component of dipole moment at
! point Rprime, and then Fourier transform these components.

      WOLD=AKD
      AK2OLD=AK(2)
      AK3OLD=AK(3)
      NGRID=8*NX*NY*NZ
      AKD2=AKD*AKD

! We assume screening function exp(-(gamma*kr)^4) so
! range/d = 1/(gamma*kd) = 400 if gamma=5e-3 and kd=0.5
! although the sums are actually continued out to
! r/d = 2/(gamma*kd) = 800 if gamma=5e-3, kd=0.5
! [screening factor = exp(-16)=1.1e-7]

! PYDDX=PYD*DX(2) = periodicity in Y direction
! PZDDX=PZD*DX(3) = periodicity in Z direction

      IF(PYD>0._WP.OR.PZD>0._WP)THEN
         WRITE(CMSGNM,FMT='(A,2F8.2,A,1PE9.2)')'PBC with PYD, PZD=',PYD, &
                                               PZD,', GAMMA=',GAMMA
         CALL WRIMSG('ESELF ',CMSGNM)
      ENDIF

      PYDDX=PYD*DX(2)
      PZDDX=PZD*DX(3)

!               Compute the actual coefficients:

! Compute 6 independent elements of 3x3 symmetric matrix A_jk,where
! A_jk*P_k = -electric field at location j due to dipole P at location k

! A_jk = (a_1  a_2  a_3)
!        (a_2  a_4  a_5)
!        (a_3  a_5  a_6)_jk

      IF(IPBC==0)THEN

! initialize CXZC(I,J,K,M) = a_M for electric field at (I,J,K)
!                            produced by a dipole at (1,1,1)
!                            and replica dipoles (if PYD or PYZ are
!                            nonzero).

! need to calculate this for all (I,J,K) for one octant:
! I running from 1 to NX, J from 1 to NY, K from 1 to NZ

! Later obtain A_jk values for other octants by using symmetry

!*** diagnostic
!         write(0,*)'eself_v8 ckpt 3: call DIRECT_CALC'
!***
        CALL DIRECT_CALC(1,1,1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZC(1,1,1,1))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 4'
!        write(0,*)'   returned from direct_calc'
!        write(0,*)'   check for NaN...'
!        jr=0
!        do i=1,nx
!           do j=1,ny
!              do k=1,nz
!                 do m=1,6
!                    if(.not.(abs(cxzc(i,j,k,m))>=0.d0).or. &
!                        abs(cxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'cxzc=',cxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'cxy checked for NaN or overflow',jr,' instances found'
!***

! At this point, CXZC(I,J,K,1-6) contains the upper triangular part of t
! symmetric 3 x 3 matrix giving the electric field at grid point (i,j,k)
! produced by a dipole at (1,1,1)

! Fill out CXZC to twice the size in each grid dimension to accomodate
! negative lags [periodicity in each dimension is assumed, so (e.g.)
! nx < i <= 2*nx is equivalent to -nx < i <= 0], exploiting symmetries,
! and then Fourier transform.

! If PYDDX=0 and PZDDX=0 , need only do direct calculation of A matrix
! for first octant, since remaining octants can be obtained by symmetry.
! After calculating A matrix, store only the first octant of the
! transform, since the others can be obtained by symmetry.

!-----------------------------------------------------------------------
! extend a_1 = a_xx(x,y,z) : a -> +a for x -> -x
!                                 +a     y -> -y
!                                 +a     z -> -z
        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=1
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 5, about to call EXTEND'
!***
        CALL EXTND(CXZC(1,1,1,1),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 6, returned from EXTEND'
!***
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'eself_v8 ckpt 7, cxzw(jx,jy,jz,m=1-6)'
!           do jx=1,2*NX
!              do jy=1,2*NY
!                 do jz=1,2*NZ
!                    write(0,fmt='(a,3i3,a,6f10.3)')                        &
!                       'eself_v8 ckpt 7',jx,jy,jz,' cxzw(jx,jy,jz,1-3)=',  &
!                       cxzw(jx,jy,jz,1),cxzw(jx,jy,jz,2),cxzw(jx,jy,jz,3)
!                 enddo
!              enddo
!           enddo
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!           write(0,*)'eself_v8 ckpt 8, cxzw(jx,jy,jz,m=1-6)'
!           do jx=1,2*NX
!              do jy=1,2*NY
!                 do jz=1,2*NZ
!                    write(0,fmt='(a,3i3,a,6f10.3)')                        &
!                       'eself_v8 ckpt 8',jx,jy,jz,' cxzw(jx,jy,jz,1-3)=',  &
!                       cxzw(jx,jy,jz,1),cxzw(jx,jy,jz,2),cxzw(jx,jy,jz,3)
!                 enddo
!              enddo
!           enddo                  ! ok thru here
!***
!**
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 9, about to call TRIM'
!***
        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,1))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 9.5, returned from TRIM'
!        do jz=1,nz
!           do jy=1,ny
!              do jx=1,nx
!                 write(0,fmt='(a,3i3,a,2f10.3)')                      &
!                    'eself_v8 ckpt 9.5',jx,jy,jz,' cxzc(jx,jy,jz,1)', &
!                    cxzc(jx,jy,jz,1)
!              enddo
!           enddo
!        enddo                 ! ok thru here
!***
!-----------------------------------------------------------------------
! extend a_2 = a_xy(x,y,z) : a -> -a for x -> -x
!                                 -a     y -> -y
!                                 +a     z -> -z

        ISYM(1)=-1
        ISYM(2)=-1
        ISYM(3)=1
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 10, about to call EXTND'
!***
        CALL EXTND(CXZC(1,1,1,2),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
!*** diagnostic
!           write(0,*)'eself_v8 ckpt 11, about to call cxfft3n'
!***
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!          write(0,*)'eself_v8 ckpt 12, returned from cxfft3n'
!***
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 13, about to call TRIM'
!***
        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,2))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 13.5, returned from TRIM'
!        do jz=1,nz
!           do jy=1,ny
!              do jx=1,nx
!                 write(0,fmt='(a,3i3,a,2f10.3)')                      &
!                    'eself_v8 ckpt 13.5',jx,jy,jz,' cxzc(jx,jy,jz,2)', &
!                    cxzc(jx,jy,jz,2)
!              enddo
!           enddo
!        enddo
!***

!-----------------------------------------------------------------------
! extend a_3 = a_xz(x,y,z) : a -> -a for x -> -x
!                                 +a     y -> -y
!                                 -a     z -> -z

        ISYM(1)=-1
        ISYM(2)=1
        ISYM(3)=-1
        CALL EXTND(CXZC(1,1,1,3),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,3))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 13.6, returned from TRIM'
!        do jz=1,nz
!           do jy=1,ny
!              do jx=1,nx
!                 write(0,fmt='(a,3i3,a,2f10.3)')                       &
!                    'eself_v8 ckpt 13.6',jx,jy,jz,' cxzc(jx,jy,jz,3)', &
!                    cxzc(jx,jy,jz,3)
!              enddo
!           enddo
!        enddo
!***

!-----------------------------------------------------------------------
! extend a_4 = a_yy(x,y,z) : a -> +a for x -> -x
!                                 +a     y -> -y
!                                 +a     z -> -z

        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=1
        CALL EXTND(CXZC(1,1,1,4),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,4))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 13.7, returned from TRIM'
!        do jz=1,nz
!           do jy=1,ny
!              do jx=1,nx
!                 write(0,fmt='(a,3i3,a,2f10.3)')                       &
!                    'eself_v8 ckpt 13.7',jx,jy,jz,' cxzc(jx,jy,jz,4)', &
!                    cxzc(jx,jy,jz,4)
!              enddo
!           enddo
!        enddo
!***
!-----------------------------------------------------------------------
! extend a_5 = a_yz(x,y,z) : a -> +a for x -> -x
!                                 -a     y -> -y
!                                 -a     z -> -z
        ISYM(1)=1
        ISYM(2)=-1
        ISYM(3)=-1
        CALL EXTND(CXZC(1,1,1,5),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,5))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 13.8, returned from TRIM'
!        do jz=1,nz
!           do jy=1,ny
!              do jx=1,nx
!                 write(0,fmt='(a,3i3,a,2f10.3)')                       &
!                    'eself_v8 ckpt 13.8',jx,jy,jz,' cxzc(jx,jy,jz,5)', &
!                    cxzc(jx,jy,jz,5)
!              enddo
!           enddo
!        enddo
!***
!-----------------------------------------------------------------------
! extend a_6 = a_zz(x,y,z) : a -> +a for x -> -x
!                                 +a     y -> -y
!                                 +a     z -> -z
        ISYM(1)=1
        ISYM(2)=1
        ISYM(3)=1
        CALL EXTND(CXZC(1,1,1,6),NX,NY,NZ,ISYM,CXZW(1,1,1,1))
        IF(CMETHD=='GPFAFT')THEN
           CALL CXFFT3N(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD=='FFTW21')THEN
           CALL CXFFTW(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ELSEIF(CMETHD.EQ.'FFTMKL')THEN
           CALL CXFFT3_MKL(CXZW(1,1,1,1),2*NX,2*NY,2*NZ,+1)
        ENDIF

        CALL TRIM(CXZW(1,1,1,1),NX,NY,NZ,CXZC(1,1,1,6))
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 13.9, returned from TRIM'
!        do jz=1,nz
!           do jy=1,ny
!              do jx=1,nx
!                 write(0,fmt='(a,3i3,a,2f10.3)')                       &
!                    'eself_v8 ckpt 13.9',jx,jy,jz,' cxzc(jx,jy,jz,6)', &
!                    cxzc(jx,jy,jz,6)
!              enddo
!           enddo
!        enddo           ! ok thru here
!***
      ELSEIF(IPBC==1)THEN

! This point is reached when PYDDX or PZDDX are nonzero.
! When PBC are used for general direction of incident wave,
! all octants of A matrix require direct calculation: symmetries valid
! for single target no longer apply because of position-dependent phases
! of replica dipoles.

!*** diagnostic
!         write(0,*)'eself_v8 ckpt 14'
!*** 
        CALL DIRECT_CALC(-1,-1,-1,NX,NY,NZ,IPBC,DX,AK,AKD,AKD2,GAMMA, &
                         PYDDX,PZDDX,CXZC(1,1,1,1))

!*** diagnostic
!        write(0,*)'eself_v8 ckpt 15'
!***
! The array CXZC(1-2*NX,1-2*NY,1-2*NZ,1-6) of A matrix coefficients
! now covers all octants.

! Fourier transform the A matrix CXZC:

        DO M=1,6
           IF(CMETHD=='GPFAFT')THEN
              CALL CXFFT3N(CXZC(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD=='FFTW21')THEN
              CALL CXFFTW(CXZC(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ELSEIF(CMETHD.EQ.'FFTMKL')THEN
              CALL CXFFT3_MKL(CXZC(1,1,1,M),2*NX,2*NY,2*NZ,+1)
           ENDIF
        ENDDO

!*** diagnostic
!        write(0,*)'eself_v8 ckpt 16'
!        do m=1,6
!           do jx=1,2*nx
!              do jy=1,2*ny
!                 do jz=1,2*nz
!                    if(.not.(abs(cxzc(jx,jy,jz,m)).lt.1.e30))then
!                       write(0,fmt='(a,4i4,a,1p2e10.3)')'jx,jy,jz,m=', &
!                          jx,jy,jz,m,' cxzc=',cxzc(jx,jy,jz,m)
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'eself_v8 ckpt 16.1'
!***

! CXZC now contains the full Fourier transform of the A convolution
! and should not be overwritten between calls to ESELF

      ENDIF
!      CALL TIMEIT('ESELF (first call)',DTIME)
!      CALL TIMEIT('ESELF',DTIME)
!-----------------------------------------------------------------------

! End of recomputation of Green-function coefficients

70    CONTINUE
!*** diagnostic
!      write(0,*)'eself_v8 ckpt 17'
!****
! Fourier transform the polarizations:

      DO M=1,3

!*** diagnostic
!         do jz=1,nz
!            do jy=1,ny
!               do jx=1,nx
!                  write(0,fmt='(a,4i3,a,2f10.3)')                &
!                     'eself_v8 ckpt 17, jx,jy,jz,m=',jx,jy,jz,m, &
!                     ' cxzp(jx,jy,jz,m)=',cxzp(jx,jy,jz,m)       !
!               enddo
!            enddo
!         enddo
!***
         CALL PAD(CXZP(1,1,1,M),NX,NY,NZ,CXZW(1,1,1,M))

!*** diagnostic
!         do jz=1,2*nz
!            do jy=1,2*ny
!               do jx=1,2*nx
!                  write(0,fmt='(a,4i3,a,2f10.3)')                &
!                     'eself_v8 ckpt 18, jx,jy,jz,m=',jx,jy,jz,m, &
!                     ' cxzw(jx,jy,jz,m)=',cxzw(jx,jy,jz,m)       !
!               enddo
!            enddo
!         enddo
!***
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 18.1: returned from PAD: ', &
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
!        write(0,*)'eself_v8 ckpt 18.2'
!        if(jr>0)then
!           write(0,*)'eself_v8 ckpt 19: cxzw checked for NaN or overflow: ', &
!              jr,' instances found'
!           stop
!        endif
!***

        IF(CMETHD=='GPFAFT')THEN

!*** diagnostic
!        write(0,*)'eself_v8 ckpt 19.9'
!        do k=1,2*nz
!           do j=1,2*ny
!              do i=1,2*nx
!                 write(0,fmt='(a,4i3,a,2f10.3)')            &
!                    'eself_v8 ckpt 19.9, i,j,k,m=',i,j,k,m, &
!                    ' cxzw(i,j,k,m)=',cxzw(i,j,k,m)         !  PROBLEM HERE
!              enddo
!           enddo
!        enddo
!***

           CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
!*** diagnostic
!        write(0,*)'eself_v8 ckpt 20'
!        do k=1,2*nz
!           do j=1,2*ny
!              do i=1,2*nx
!                 write(0,fmt='(a,4i3,a,2f10.3)')          &
!                    'eself_v8 ckpt 20, i,j,k,m=',i,j,k,m, &
!                    ' cxzw(i,j,k,m)=',cxzw(i,j,k,m)       ! PROBLEM
!              enddo
!           enddo
!        enddo
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
!        if(jr>0)then
!           write(0,*)'eself_v8 ckpt 21: cxzw checked for NaN or overflow', &
!                     jr,' instances found'
!           stop
!        endif
!***
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,+1)
         ENDIF
      ENDDO
!*** diagnostic
!      write(0,*)'eself_v8 ckpt 22, IPBC=',IPBC
!***
!***********************************************************************

! Multiply by F.t. of Green function.

      IF(IPBC==0)THEN

!*** diagnostic
!         write(0,*)'eself_v8 ckpt 23, check cxzc(i,j,k,m=1,6)'
!        jr=0
!        do i=1,2*nx
!           do j=1,2*ny
!              do k=1,2*nz
!                 do m=1,6
!                    if(.not.(abs(cxzc(i,j,k,m))>=0.d0).or. &
!                        abs(cxzc(i,j,k,m))>=1.d100)then
!                       write(0,*)'i,j,k,m=',i,j,k,m,'cxzc=',cxzc(i,j,k,m)
!                       jr=jr+1
!                    endif
!                 enddo
!              enddo
!           enddo
!        enddo
!        write(0,*)'eself_v8 ckpt 24: cxzc(i,j,k,m=1-6) checked for NaN ', &
!                  'or overflow: ',jr,' instances found'
!        if(jr>0)stop
!        write(0,*)'eself_v8 ckpt 25: check cxzw(i,j,k,m=1-3)'
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
!        write(0,*)'eself_v8 ckpt 26: cxzw(i,j,k,m=1-3) checked for NaN or ', &
!                  'overflow: ',jr,' instances found'
!        if(jr>0)stop
!***
!***

! If IPBC=0, then only one octant of F.t. of Green function has been
!            stored, but can recover others using symmetry.

#ifdef openmp
!$OMP PARALLEL DO                                            &
!$OMP&   PRIVATE(K,J,I,KSGN,KR,JSGN,JR,ISGN,IR)              &
!$OMP&   PRIVATE(CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

        DO K=1,2*NZ
           KSGN=NINT(SIGN(1._WP,NZ+1.5_WP-K))
           KR=MIN(K,2*NZ+2-K)
!*** diagnostic
!         write(0,*)'eself_v8 ckpt 27'
!***
           DO J=1,2*NY
              JSGN=NINT(SIGN(1._WP,NY+1.5_WP-J))
              JR=MIN(J,2*NY+2-J)
!*** diagnostic
!         write(0,*)'eself_v8 ckpt 28'
!***
              DO I=1,2*NX
                 ISGN=NINT(SIGN(1._WP,NX+1.5_WP-I))
                 IR=MIN(I,2*NX+2-I)
!*** diagnostic
!         write(0,*)'eself_v8 ckpt 29'
!***
                 CXXX=CXZC(IR,JR,KR,1)
                 CXXY=CXZC(IR,JR,KR,2)*(ISGN*JSGN)
                 CXXZ=CXZC(IR,JR,KR,3)*(ISGN*KSGN)
                 CXYY=CXZC(IR,JR,KR,4)
                 CXYZ=CXZC(IR,JR,KR,5)*(JSGN*KSGN)
                 CXZZ=CXZC(IR,JR,KR,6)
!*** diagnostic
!              if(.not.(abs(cxxx)>=0.d0))write(0,*) &
!                'ir,jr,kr,cxzc(ir,jr,kr,1)=',ir,jr,kr,cxzc(ir,jr,kr,1)
!              if(.not.(abs(cxxy)>=0.d0))write(0,*) &
!                'ir,jr,kr,cxzc(ir,jr,kr,2)=',ir,jr,kr,cxzc(ir,jr,kr,2)
!              if(.not.(abs(cxxz)>=0.d0))write(0,*) &
!                 'ir,jr,kr,cxzc(ir,jr,kr,3)=',ir,jr,kr,cxzc(ir,jr,kr,3)
!              if(.not.(abs(cxyy)>=0.d0))write(0,*) &
!                  'ir,jr,kr,cxzc(ir,jr,kr,4)=',ir,jr,kr,cxzc(ir,jr,kr,4)
!              if(.not.(abs(cxyz)>=0.d0))write(0,*) &
!                  'ir,jr,kr,cxzc(ir,jr,kr,5)=',ir,jr,kr,cxzc(ir,jr,kr,5)
!              if(.not.(abs(cxzz)>=0.d0))write(0,*) &
!                  'ir,jr,kr,cxzc(ir,jr,kr,6)=',ir,jr,kr,cxzc(ir,jr,kr,6)
!***
                 CXEX=CXXX*CXZW(I,J,K,1)+CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                 CXEY=CXXY*CXZW(I,J,K,1)+CXYY*CXZW(I,J,K,2)+CXYZ*CXZW(I,J,K,3)
                 CXEZ=CXXZ*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,2)+CXZZ*CXZW(I,J,K,3)
                      
!*** diagnostic
!              if(.not.(abs(cxex+cxey+cxez)>=0.d0))then
!                 write(0,*)'i,j,j,ir,jr,kr=',i,j,k,ir,jr,kr
!                 write(0,*)'cxzc(ir,jr,kr,1)=',cxzc(ir,jr,kr,1)
!                 write(0,*)'cxzc(ir,jr,kr,2)=',cxzc(ir,jr,kr,2)
!                 write(0,*)'cxzc(ir,jr,kr,3)=',cxzc(ir,jr,kr,3)
!                 write(0,*)'cxzc(ir,jr,kr,4)=',cxzc(ir,jr,kr,4)
!                 write(0,*)'cxzc(ir,jr,kr,5)=',cxzc(ir,jr,kr,5)
!                 write(0,*)'cxzc(ir,jr,kr,6)=',cxzc(ir,jr,kr,6)
!                 write(0,*)'      cxex=',cxex
!                 write(0,*)'      cxey=',cxey
!                 write(0,*)'      cxez=',cxez
!                 write(0,*)'      cxxx=',cxxx
!                 write(0,*)'      cxxy=',cxxy
!                 write(0,*)'      cxxz=',cxxz
!                 write(0,*)'      cxyy=',cxyy
!                 write(0,*)'      cxyz=',cxyz
!                 write(0,*)'      cxzz=',cxzz
!                 write(0,*)'      cxzw(i,j,k,1)=',cxzw(i,j,k,1)
!                 write(0,*)'      cxzw(i,j,k,2)=',cxzw(i,j,k,2)
!                 write(0,*)'      cxzw(i,j,k,3)=',cxzw(i,j,k,3)
!                 stop
!              endif
!***
                 CXZW(I,J,K,1)=CXEX
                 CXZW(I,J,K,2)=CXEY
                 CXZW(I,J,K,3)=CXEZ
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

!*** diagnostic
!         write(0,*)'eself_v8 ckpt 30'
!***

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ELSEIF(IPBC==1)THEN

!*** diagnostic
!         write(0,*)'eself_v8 ckpt 31'
!***

! If IPBC=1, then the full F.t. of the Green function has been stored.

#ifdef openmp
!$OMP PARALLEL DO                                                  &
!$OMP&   PRIVATE(K,J,I,CXXX,CXXY,CXXZ,CXYY,CXYZ,CXZZ,CXEX,CXEY,CXEZ)
#endif

         DO K=1,2*NZ
            DO J=1,2*NY
               DO I=1,2*NX
                  CXXX=CXZC(I,J,K,1)
                  CXXY=CXZC(I,J,K,2)
                  CXXZ=CXZC(I,J,K,3)
                  CXYY=CXZC(I,J,K,4)
                  CXYZ=CXZC(I,J,K,5)
                  CXZZ=CXZC(I,J,K,6)
                  CXEX=CXXX*CXZW(I,J,K,1)+CXXY*CXZW(I,J,K,2)+CXXZ*CXZW(I,J,K,3)
                  CXEY=CXXY*CXZW(I,J,K,1)+CXYY*CXZW(I,J,K,2)+CXYZ*CXZW(I,J,K,3)
                  CXEZ=CXXZ*CXZW(I,J,K,1)+CXYZ*CXZW(I,J,K,2)+CXZZ*CXZW(I,J,K,3)
                  CXZW(I,J,K,1)=CXEX
                  CXZW(I,J,K,2)=CXEY
                  CXZW(I,J,K,3)=CXEZ
               ENDDO
            ENDDO
         ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

      ENDIF

! Inverse Fourier transform to obtain electric field:

      DO M=1,3
         IF(CMETHD=='GPFAFT')THEN
            CALL CXFFT3N(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ELSEIF(CMETHD=='FFTW21')THEN
            CALL CXFFTW(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ELSEIF(CMETHD.EQ.'FFTMKL')THEN
            CALL CXFFT3_MKL(CXZW(1,1,1,M),2*NX,2*NY,2*NZ,-1)
         ENDIF

!*** diagnostic
!         write(0,*)'eself_v8 ckpt 32'
!***
!***********************************************************************

! Note: the Convex FFT routine already normalizes result.
!       For other FFT routines need to divide result by NGRID

         IF(CMETHD=='CONVEX')THEN
            DO K=1,NZ
               DO J=1,NY
                  DO I=1,NX
                     CXZE(I,J,K,M)=CXZW(I,J,K,M)
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
                     CXZE(I,J,K,M)=FAC*CXZW(I,J,K,M)
                  ENDDO
               ENDDO


            ENDDO

#ifdef openmp
!$OMP END PARALLEL DO
#endif

         ENDIF
      ENDDO

! 14.10.16 (BTD) sanity check to find problem with SPHRN_PBC...
!*** diagnostic
!      write(0,*)'eself_v8 ckpt 33'
!      do k=1,3
!         do jz=1,nz
!            do jy=1,ny
!               do jx=1,nx
!                  if(.not.(abs(cxze(jx,jy,jz,k)).lt.1.e32))then
!                     write(0,*)'eself_v8 ckpt 34: *** problem ***'
!                     write(0,*)'  jx,jy,jz,k=',jx,jy,jz,k
!                     write(0,*)'  cxze(jx,jy,jz,k)=',cxze(jx,jy,jz,k)
!                  endif
!               enddo
!            enddo
!         enddo
!      enddo
!      write(0,*)'eself_v8 ckpt 35'
!***

      RETURN
    END SUBROUTINE ESELF

!***********************************************************************

    SUBROUTINE EXTND(CXA,NX,NY,NZ,ISYM,CXB)
      USE DDPRECISION,ONLY: WP

#ifdef openmp
      USE OMP_LIB   !Art for OpenMP function declarations
#endif

      IMPLICIT NONE

!       Using symmetries, extend first octant of coefficient matrix to
!       other octants.

! Originally written by Jeremy Goodman
! 080620 (ASL) modified to use OpenMP
! 080627 (BTD) eself_v4: further mods related to OpenMP
!              reordered several nested do loops
!              introduced variables J2 and K2 to speed computations
! end history
! Copyright (C) 1993,2008 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: NX,NY,NZ
      COMPLEX(WP) :: CXA(NX+1,NY+1,NZ+1),CXB(2*NX,2*NY,2*NZ)
      INTEGER :: ISYM(3)

! Local variables:

      INTEGER :: I,J,J2,K,K2
      COMPLEX(WP) :: CXZERO

! SAVE statement:

      SAVE CXZERO

! DATA statement:

      DATA CXZERO/(0._WP,0._WP)/
!***********************************************************************

#ifdef openmp
!$OMP PARALLEL              & 
!$OMP&   PRIVATE(I,J,J2,K,K2)
!$OMP DO
#endif

      DO K=1,NZ

         DO J=1,NY
            DO I=1,NX
               CXB(I,J,K)=CXA(I,J,K)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP ENDDO
#endif

!-----------------------------------------------------------------------
! x -> -x

!btd 080627 moved I=NX+1 out of loop
! the SINGLE directive specifies that enclosed code is to be executed
! by only one thread in the team

#ifdef openmp
!$OMP SINGLE
#endif

      CXB(NX+1,1:NY,1:NZ)=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

      DO K=1,NZ

         DO J=1,NY
            DO I=NX+2,2*NX
               CXB(I,J,K)=CXA(2*NX+2-I,J,K)*ISYM(1)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
#endif

!-----------------------------------------------------------------------
! y -> -y

!btd 080627 moved J=NY+1 out of loop, switched order of loops I and J
! the SINGLE directive specifies that enclosed code is to be executed
! by only one thread in the team

#ifdef openmp
!$OMP SINGLE
#endif

      CXB(1:2*NX,NY+1,1:NZ)=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

      DO K=1,NZ

         DO J=NY+2,2*NY
            J2=2*NY+2-J
            DO I=1,2*NX
               CXB(I,J,K)=CXB(I,J2,K)*ISYM(2)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
#endif

!-----------------------------------------------------------------------
! z -> -z

!Art pulling the 3rd dimension to the outer most loop.
!Art we will do this expression in only the thread that has NZ+1
!btd 080627 reorder loops: J,I,K -> K,J,I
! the SINGLE directive specifies that enclosed code is to be executed
! by only one thread in the team

#ifdef openmp
!$OMP SINGLE
#endif

      CXB(1:2*NX,1:2*NY,NZ+1)=CXZERO

#ifdef openmp
!$OMP END SINGLE
!$OMP DO 
#endif

      DO K=NZ+2,2*NZ

         K2=2*NZ+2-K
         DO J=1,2*NY
            DO I=1,2*NX
               CXB(I,J,K)=CXB(I,J,K2)*ISYM(3)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
    END SUBROUTINE EXTND
!***********************************************************************

    SUBROUTINE TRIM(CXB,NX,NY,NZ,CXA)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Copy the first octant of CXB into CXA

! Originally written by Jeremy Goodman
! History:
! 92.04.20 (BTD) removed ISYM from argument list of TRIM:
! 06.09.28 (BTD) eself v2.0
! 08.06.27 (BTD) eself_v4
!                
! Copyright (C) 1993,2006 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:

      INTEGER :: NX,NY,NZ
      COMPLEX(WP) ::          &
         CXA(NX+1,NY+1,NZ+1), &
         CXB(2*NX,2*NY,2*NZ)

! Local variables:

      INTEGER :: I,J,K

!***********************************************************************

#ifdef openmp
!$OMP PARALLEL             &
!$OMP&    PRIVATE(I,J,K)   &
!$OMP&    SHARED(NX,NY,NZ) &
!$OMP&    SHARED(CXA,CXB)
!$OMP DO
#endif

      DO K=1,NZ+1

         DO J=1,NY+1
            DO I=1,NX+1
               CXA(I,J,K)=CXB(I,J,K)
            ENDDO
         ENDDO

      ENDDO

#ifdef openmp
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
    END SUBROUTINE TRIM

