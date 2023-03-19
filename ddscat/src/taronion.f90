    SUBROUTINE TARONION(A1,A2,ICORE,AOUT,RINROUT,DX,X0,CDESCR,IOSHP,        &
                        MXNAT,NAT,IXYZ,NCOMP_NEED,ICOMP,BETADF,PHIDF,THETADF)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!** Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: ICORE,IOSHP,MXNAT,NAT,NCOMP_NEED
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: AOUT,RINROUT
      REAL(WP) ::        &
         A1(3),          &
         A2(3),          &
         BETADF(MXNAT),  &
         DX(3),          &
         PHIDF(MXNAT),   &
         THETADF(MXNAT), &
         X0(3)

!** Local variables:

      INTEGER :: JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2
      REAL(WP) :: AX2,AY2,AZ2,PI,R,R2,RIN2,ROUT2,RSINTHETA,RYZ2,RZ2, &
                  X,XOFF,Y,YOFF,Z,ZOFF

!***********************************************************************
! Routine to construct onion-like sphere; 
!            inner core can be vacuum (if ICORE=0)
!            or isotropic material (if ICORE=1)
!         
!     (outer diameter)/d = AOUT
!     (inner radius/outer radius) = RINROUT
!     if ICORE=0: vacuum interior to inner radius
!                 shell: unixaxial material with ICOMP=1 for E parallel c
!                                                ICOMP=2 for E perp c
!                 NCOMP_NEED=2

!     if ICORE=1: core:  isotropic material with ICOMP=1
!                 shell: unixaxial material with ICOMP=2 for E parallel c
!                                                ICOMP=3 for E perp c
!                 NCOMP_NEED=3

!     target axes set to a1=x_TF
!                        a2=y_TF
! Input:
!        ICORE = 0 or 1 (inner void or inner material)
!        AOUT =(outer diameter)/d   (d=lattice spacing)
!        RINROUT=(inner radius)/(outer radius)
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT=dimensioning information (max number of atoms)
! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Fr
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Fr
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=(x-x0(1))/d,(y-x0(2))/d,(z-x0(3))/d
!                        for atoms of target
!        CDESCR=description of target (up to 67 characters)
!        NCOMP_NEED = number of dielectric functions needed
!                   = 2 for ICORE=0, 3 for ICORE=1
!        ICOMP(1-NAT,1-3)=composition identifier for dielectric
!                         components 1-3 at locations 1-NAT
!        BETADF(1-NAT)  = angle (radians) defining orientation of
!                         "Dielectric Frame" for locations 1-NAT
!        PHIDF(1-NAT)   = angle (radians) defining orientation of
!                         "Dielectric Frame" for locations 1-NAT
!        THETADF(1-NAT) = angle (radians) defining orientation of
!                         "Dielectric Frame" for locations 1-NAT
!        X0(1-3)=(location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).
!        The TF origin is set to be centroid of ellipsoid.

! B.T.Draine, Princeton Univ. Obs.
! History:
! 16.07.18 (BTD): written using tarell.f90 as starting point
! 20.04.11 (BTD): cosmetic change
! 20.07.16 (BTD): corrected error for case ICORE=1
! end history

! Copyright (C) 2016,2020
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Routine to construct pseudo-ellipsoidal target aray.
! With occupied array sites contained within ellipsoidal surface
! defined by (X/AX*d)**2+(Y/AY*d)**2+(Z/AZ*d)**2=0.25
! Ideal volume V=(pi/6)*AX*AY*AZ*d**3
! where d = effective lattice spacing

! Dipoles are located on lattice at sites
! (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers
!                                 XOFF,YOFF,ZOFF=constants

! For sphere: call with AX=AY=AZ
! For spheroid: call with AX=AY (or AY=AZ or AX=AZ)
! B.T.Draine, Princeton Univ. Obs., 88.08.12

! Criterion for choosing XOFF:
! If AX is close to an even integer, take XOFF=1/2
! If AX is close to an odd integer, take XOFF=0
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,fmt='(a)')'taronion ckpt 0'
!      write(0,fmt='(a,i10)')'   icore=',icore
!      write(0,fmt='(a,1pe10.3)')'   aout=',aout
!      write(0,fmt='(a,1pe10.3)')'   rinrout=',rinrout
!***
      PI=4.0_WP*ATAN(1.0_WP)
      JX=INT(AOUT/DX(1)+0.5_WP)
      IF(2*(JX/2)==JX)THEN
         XOFF=0.5_WP
      ELSE
         XOFF=0.0_WP
      ENDIF

      IF(ICORE==0)THEN
         NCOMP_NEED=2
      ELSE
         NCOMP_NEED=3
      ENDIF

      YOFF=XOFF
      ZOFF=XOFF

      LMX1=-INT(0.5_WP*AOUT/DX(1)+0.5_WP)
      LMX2=INT(0.5_WP*AOUT/DX(1)-0.5_WP)
      LMY1=-INT(0.5_WP*AOUT/DX(2)+0.5_WP)
      LMY2=INT(0.5_WP*AOUT/DX(2)-0.5_WP)
      LMZ1=-INT(0.5_WP*AOUT/DX(3)+0.5_WP)
      LMZ2=INT(0.5_WP*AOUT/DX(3)-0.5_WP)
      ROUT2=(0.5_WP*AOUT)**2
      RIN2=RINROUT**2*ROUT2
      NAT=0

! Specify target axes A1 and A2
! Convention: A1=(1,0,0) in target frame
!             A2=(0,1,0) in target frame

      DO JX=1,3
         A1(JX)=0.0_WP
         A2(JX)=0.0_WP
      ENDDO
      A1(1)=1.0_WP
      A2(2)=1.0_WP

! Determine list of occupied sites.

      DO JZ=LMZ1,LMZ2
         Z=(REAL(JZ,KIND=WP)+ZOFF)*DX(3)
         RZ2=Z*Z
         IF(RZ2<ROUT2)THEN
            DO JY=LMY1,LMY2
               Y=(REAL(JY,KIND=WP)+YOFF)*DX(2)
               RYZ2=RZ2+Y*Y
               IF(RYZ2<ROUT2)THEN
                  DO JX=LMX1,LMX2
                     X=(REAL(JX,KIND=WP)+XOFF)*DX(1)
                     R2=RYZ2+X*X

! Check whether site is in shell:

                     IF(R2<ROUT2)THEN
                        R=SQRT(R2)
                        IF(ICORE==0)THEN   ! vacuum core

                           IF(R2>RIN2)THEN

! Site is occupied:
                              NAT=NAT+1
                              IXYZ(NAT,1)=JX
                              IXYZ(NAT,2)=JY
                              IXYZ(NAT,3)=JZ

! Set composition:
! uniaxial material (e.g., graphite):

                              ICOMP(NAT,1)=1
                              ICOMP(NAT,2)=2
                              ICOMP(NAT,3)=2
                           ENDIF

                        ELSE   ! isotropic material 1 in core

                           NAT=NAT+1
                           IXYZ(NAT,1)=JX
                           IXYZ(NAT,2)=JY
                           IXYZ(NAT,3)=JZ

! Set composition: 
! core = isotropic material 1
! shell = uniaxial material (2 and 3)

                           IF(R2<=RIN2)THEN

! isotropic material 1:
                              ICOMP(NAT,1)=1
                              ICOMP(NAT,2)=1
                              ICOMP(NAT,3)=1

                           ELSE
! uniaxial material (e.g., graphite):

                              ICOMP(NAT,1)=2
                              ICOMP(NAT,2)=3
                              ICOMP(NAT,3)=3

                           ENDIF
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARONION',' NAT.GT.MXNAT ')
      ENDIF

! Set X0 so that TF origin is at centroid

      DO JY=1,3
         X0(JY)=0._WP
         DO JX=1,NAT
            X0(JY)=X0(JY)+REAL(IXYZ(JX,JY))
         ENDDO
         X0(JY)=-X0(JY)/REAL(NAT)
      ENDDO

! Determine rotation angles BETADF,PHIDF,THETADF to rotate the local
! Dielectric Frame to desired orientation in the Target Frame
! We want axis e1 of DF to be in radial direction:

      DO JX=1,NAT
         X=REAL(IXYZ(JX,1)+X0(1))*DX(1)
         Y=REAL(IXYZ(JX,2)+X0(2))*DX(2)
         Z=REAL(IXYZ(JX,3)+X0(3))*DX(3)
         RYZ2=Y**2+Z**2
         R2=X**2+RYZ2
         R=SQRT(R2)

         IF(R.GT.1.E-6_WP)THEN
            THETADF(JX)=ACOS(X/R)
            RSINTHETA=SQRT(RYZ2)
            PHIDF(JX)=ACOS(Y/RSINTHETA)      ! this gives 0 < phidf < pi

! if z < 0, then shift to other solution

            IF(Z<0.)PHIDF(JX)=2._WP*PI-PHIDF(JX)

         ELSE   ! special treatment for point at origin
            THETADF(JX)=0.
            PHIDF(JX)=0.
         ENDIF

! we don't care about rotations around axis e1, so set BETADF=0

         BETADF(JX)=0.
      ENDDO

!***********************************************************************
! Write target description into string CDESCR
      WRITE(CDESCR,FMT='(A,I7,A,F6.3,F6.3,F6.3,A)')' Onion,',NAT, &
            ' dipoles,',DX(1),DX(2),DX(3),'=x,y,z lattice spacing'
!***********************************************************************

      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)AOUT,RINROUT,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3), &
                                 THETADF(JX),PHIDF(JX),BETADF(JX)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARONION  onion-like grain; AOUT,RINROUT=',2F8.4,/,  &
         I10,' = NAT',/,                                             &
         3F10.6,' = A_1 vector',/,                                   &
         3F10.6,' = A_2 vector',/,                                   &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,             &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ',   &
               'for dipole 0 0 0',/,                                 &
         '     JA   IX   IY   IZ ICOMP  theta_DF   phi_DF   beta_DF')
9030  FORMAT(I7,3I5,3I2,3F10.7)
    END SUBROUTINE TARONION
