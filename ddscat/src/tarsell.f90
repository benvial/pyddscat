    SUBROUTINE TARSELL(A1,A2,AX,AY,AZ,RE,RN,DX,X0,CDESCR,IOSHP,MXNAT,&
                       NAT,IXYZ,ICOMP)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!** Arguments:

      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP, MXNAT, NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: AX, AY, AZ, RE, RN
      REAL(WP) :: &
         A1(3),   &
         A2(3),   &
         DX(3),   &
         X0(3)

!** Local variables:

      INTEGER :: JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2
      REAL(WP) :: AX2,AY2,AZ2,R,RYZ2,RZ2,X,XOFF,Y,YOFF,Z,ZOFF

!***********************************************************************
! Routine to construct "super ellipsoid" from "atoms"
! Input:
!        AX=(x-length)/d    (d=lattice spacing)
!        AY=(y-length)/d
!        AZ=(z-length)/d
!        RE = East-West exponent of [1], (e of [2]; epsilon2 of [3])
!        RN = North-South exponent of [1] (n of [2]; epsilon1 of [3])
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
!        ICOMP(1-NAT,1-3)=1 (composition identifier)
!        X0(1-3)=(location/d) in Target Frame corresponding to dipole with
!                IXYZ=(0,0,0).
!        TF origin is set to be centroid of ellipsoid.
!
! [1] Barr, A. H. (1981). Superquadrics and angle preserving transformation,
!              IEEE Computer Graphics and Applications, 1, 11-23.
! [2] Lin, W., Bi, L., & Dubovik, O. (2018).  Assessing superspheroids in
!              modeling the scattering matrices of dust aerosols, JGR
!              Atmospheres, 123, 13,917-13,943,
!              https://doi.org/10.1029/2018JD029464 
! [3] Jaklic A., Leonardis A., Solina F. (2000) Superquadrics and Their
!              Geometric Properties. In: Segmentation and Recovery of
!              Superquadrics. Computational Imaging and Vision, vol 20.
!              Springer, Dordrecht
!       (https://cse.buffalo.edu/~jryde/cse673/files/superquadrics.pdf)
!
! B.T.Draine, Princeton Univ. Obs.
! This routine contributed by M.J.Wolff
! History:
! 19.05.03 (mjw): created from TARELL
! end history

! Copyright (C) 1993,1996,1997,1998,2000,2007,2008,2019 
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Routine to construct pseudo-superellipsoidal target aray.
! With occupied array sites contained within ellipsoidal surface
! defined by equation 1 of [2]:
! ( (x/a)**(2/re) + (y/b)**(2/re) )**(re/rn)  +  (z/c)**(2/rn)  = 1
!  a = AX*d/2, b = AY*d/2, c = AZ*d/2
! Ideal volume, according to Wikipedia involves beta functions, so not
! worrying about this right now.
! where d = effective lattice spacing

! Dipoles are located on lattice at sites
! (x,y,z)=(I+XOFF,J+YOFF,Z+KOFF), I,J,K=integers
!                                 XOFF,YOFF,ZOFF=constants

! B.T.Draine, Princeton Univ. Obs.

! Criterion for choosing XOFF:
! If AX is close to an even integer, take XOFF=1/2
! If AX is close to an odd integer, take XOFF=0
!-----------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'entered tarsell,ax,ay,az=',ax,ay,az
!***
      JX=INT(AX/DX(1)+.5_WP)
      IF(2*(JX/2)==JX)THEN
         XOFF=0.5_WP
      ELSE
         XOFF=0._WP
      ENDIF

! Similar criterion for YOFF:

      JY=INT(AY/DX(2)+.5_WP)
      IF(2*(JY/2)==JY)THEN
         YOFF=0.5_WP
      ELSE
         YOFF=0._WP
      ENDIF

! Similar criterion for ZOFF:

      JZ=INT(AZ/DX(3)+.5_WP)
      IF(2*(JZ/2)==JZ)THEN
         ZOFF=0.5_WP
      ELSE
         ZOFF=0._WP
      ENDIF

      LMX1=-INT(.5_WP*AX/DX(1)+.5_WP)
      LMX2=INT(.5_WP*AX/DX(1)-.25_WP)
      LMY1=-INT(.5_WP*AY/DX(2)+.5_WP)
      LMY2=INT(.5_WP*AY/DX(2)-.25_WP)
      LMZ1=-INT(.5_WP*AZ/DX(3)+.5_WP)
      LMZ2=INT(.5_WP*AZ/DX(3)-.25_WP)
      AX2=AX*AX
      AY2=AY*AY
      AZ2=AZ*AZ
      NAT=0

! Specify target axes A1 and A2
! Convention: A1=(1,0,0) in target frame
!             A2=(0,0,1) in target frame

      DO JX=1,3
         A1(JX)=0._WP
         A2(JX)=0._WP
      ENDDO
      A1(1)=1._WP
      A2(3)=1._WP

!*** diagnostic
!      write(0,*)'tarsell, ckpt 2'
!***

! Determine list of occupied sites.

      DO JZ=LMZ1,LMZ2
         Z=(REAL(JZ,KIND=WP)+ZOFF)*DX(3)
         RZ2=(Z*Z/(AZ2/4._WP))**(1._WP / RN)
         IF(RZ2<1._WP)THEN
            DO JY=LMY1,LMY2
               Y=(REAL(JY,KIND=WP)+YOFF)*DX(2)
               RYZ2=(Y*Y/(AY2/4._WP))**(1._WP / RE)
               IF(RYZ2<1._WP)THEN
                  DO JX=LMX1,LMX2
                     X=(REAL(JX,KIND=WP)+XOFF)*DX(1)
                     R = (X*X/(AX2/4._WP))**(1._WP / RE)
                     R = (RYZ2 + R)**(RE/RN) + RZ2
                     IF(R<1._WP)THEN
! Site is occupied:
                        NAT=NAT+1
                        IXYZ(NAT,1)=JX
                        IXYZ(NAT,2)=JY
                        IXYZ(NAT,3)=JZ
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ENDIF
      ENDDO

!*** diagnostic
!      write(0,*)'tarsell, ckpt 3'
!***
      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARELL',' NAT.GT.MXNAT ')
      ENDIF

! Set composition

      DO JX=1,3
         DO JZ=1,NAT
            ICOMP(JZ,JX)=1
         ENDDO
      ENDDO

!*** diagnostic
!      write(0,*)'tarsell, ckpt 4'
!***
! Set X0 so that TF origin is at centroid

      DO JY=1,3
         X0(JY)=0._WP
         DO JX=1,NAT
            X0(JY)=X0(JY)+REAL(IXYZ(JX,JY))
         ENDDO
         X0(JY)=-X0(JY)/REAL(NAT)
      ENDDO
!*** diagnostic
!      write(0,*)'tarsell, ckpt 5'
!***
!***********************************************************************
! Write target description into string CDESCR
      WRITE(CDESCR,FMT='(A,I7,A,F6.3,F6.3,F6.3,A)')' S-Ellipse',NAT, &
          ' dipoles,',DX(1),DX(2),DX(3),'=x,y,z lattice spacing'
!***********************************************************************
!*** diagnostic
!      write(0,*)'tarsell,ckpt 6, ioshp=',ioshp
!***
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)AX,AY,AZ,RE,RN,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TARSELL  super ellipsoidal grain; AX,AY,AZ,E,N=',5F8.4,/,  &
         I10,' = NAT',/,                                           &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT(I7,3I5,3I2)
    END SUBROUTINE TARSELL
