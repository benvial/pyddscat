    SUBROUTINE TAROCT(A1,A2,AX,DX,X0,CDESCR,IOSHP,MXNAT,NAT,IXYZ,ICOMP)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

!** Arguments:
      CHARACTER :: CDESCR*67
      INTEGER :: IOSHP,MXNAT,NAT
      INTEGER*2 :: ICOMP(MXNAT,3)
      INTEGER :: IXYZ(MXNAT,3)
      REAL(WP) :: AX
      REAL(WP) :: A1(3),A2(3),DX(3),X0(3)
!** Local variables:
      INTEGER :: I,IMAX,IMIN,J,JMAX,JMIN,JX,K,KMAX,KMIN
      REAL(WP) :: FTHETA,FZ1,FZ2,FZ3,S,X,XOFF,Y,YMAX,YMIN,YOFF, &
                  Z,ZMAX,ZMAX0,ZMAX_0,ZOFF
!***********************************************************************

! Routine to construct regular octahedral target array
! Base of octrahedron is the xy-plane
! apex points on z-axis

! Input:
!        AX = length of one edge of octahedron (assume equilateral triangles
!             are used.
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        IOSHP=device number for "target.out" file
!             =-1 to suppress printing of "target.out"
!        MXNAT = dimensioning information
! Returns:
!        A1(1-3) = unit vector (1,0,0) (along one axis of tetrahedron)
!        A2(1-3) = unit vector (0,1,0)
!        CDESCR = string describing target (up to 67 chars.)
!        NAT = number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms in target
!        ICOMP(1-NAT,1-3)=1 (composition identifier)
!        X0(1-3) = location/d in TF of site with IXYZ=0 0 0
! Occupied array sites are those within octahedral surface
! Size:
!    S=length of each side of tetrahedron
!    Volume = S**3  *  sqrt(2)/3
! Orientation:
!    Center of mass is at origin.
!    Base is xy-plane.
!    S=length of one side (in lattice units)
!    Vertices are at
! Length in x direction = S
!           y           = S
!           z           = S * 2 / sqrt(2)
! Angle(AOB)=arccos(-1/3)=109.4712 deg.
! Occupied sites are assumed to be located at
! (X,Y,Z)=(I+XOFF,J+YOFF,K+ZOFF)
! where I,J,K are integers, and XOFF,YOFF,ZOFF are constants.
! Program sets XOFF,YOFF,ZOFF depending on choice of parameter S.

! Criterion for choosing XOFF and YOFF:
!    Base of octahedron is located in the xy plane
!    Choose XOFF and YOFF in order to have number of dipoles along
!    these edges as close as possible to S
!    e.g., if S=odd integer, then take {X,Y}OFF=0 to place dipoles
!          at integral locations
!          if S=even integer, then take {X,Y}OFF=1/2 to place dipoles
!          at half-integral locations
! Criterion for choosing ZOFF:
!    Choose ZOFF to have number of dipoles running through origin and along
!    along z-axis to have number of dipoles as close as possible to
!    S * 2/sqrt(2)....same odd,even scheme as above.
!    

! B.T.Draine, Princeton Univ. Obs.,
! History:
! 12.02.28 (MJW): written by Michael J. Wolff using tartet.f90
!                 as template
! end history

! Copyright (C) 1993,1995,1996,1997,1998,2000,2007,2008,2012
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! Current version of TAROCT is restricted to cubic lattices

      IF(DX(1)/=1._WP.OR.DX(2)/=1._WP)THEN
         CALL ERRMSG('FATAL','TAROCT',                          &
                     ' taroct does not support noncubic lattice')
      ENDIF
      S = AX

! Determine whether S is closest to even or odd integer.
!  (Temporarily let IMIN be integer which S is closest to)
      IMIN=INT(S+0.5_WP)
! If IMIN is even, then {X,Y}OFF=0.5
! If IMIN is odd, then {X,Y}OFF=0.


! Set XOFF (and IMIN,IMAX):
      XOFF=0._WP
      IF(IMIN-2*(IMIN/2)==0)XOFF=0.5_WP
      IMIN=-INT(0.5_WP*S+XOFF)
      IMAX=IMIN+INT(S+0.5_WP)-1

! Set YOFF (and JMIN,JMAX):
! assume same as X
      YOFF=XOFF
      JMIN=IMIN
      JMAX=IMAX

! Set ZOFF (and KMIN,KMAX):
! Determine whether ZMAX_0 = S * 2/sqrt(2) is closest to even or odd integer.
!  (Temporarily let KMIN be integer which ZMAX_0 is closest to)
      ZMAX_0 = S*2.0_WP/sqrt(2._WP)
      KMIN=INT(ZMAX_0+0.5_WP)

! If KMIN is even, then ZOFF=0.5
! If KMIN is odd, then ZOFF=0.
      ZOFF=0._WP
      IF(KMIN-2*(KMIN/2)==0)ZOFF=0.5_WP
      KMIN=-INT(0.5_WP*ZMAX_0+ZOFF)
      KMAX=KMIN+INT(ZMAX_0+0.5_WP)-1
      write(0,*) KMIN,KMAX,ZOFF

! Determine list of occupied sites.

      NAT=0
      DO I=IMIN,IMAX
         X=REAL(I,KIND=WP)+XOFF

         DO J=JMIN,JMAX
            Y=REAL(J,KIND=WP)+YOFF

! ZMAX=largest value of Z which can occur for this (X,Y)
! FTHETA = angle of between line connecting the origin to (X,Y) AND the axis
!                for the larger of the X,Y values
! FZ1 = radial distance of (X,Y)
! FZ2 = distance from origin to side of base (square) passing through (X,Y)
! FZ3 = distance of (X,Y) to base for segment defined in FZ2, to be used to
!           construct similar triangle relationship
            if (abs(X) >= abs(Y)) then
               FTHETA = atan2(Y,X)
            else
               FTHETA = atan2(X,Y)
            endif
            FZ1 = sqrt(X*X + Y*Y)
            FZ2 = abs( (S/2.0_WP) / cos(FTHETA) )
            FZ3 = abs(FZ2) - FZ1
            ZMAX= (ZMAX_0 / 2.0_WP) * (FZ3 / FZ2)

! alternate calculation that uses the equation of a line that
! connects the apex and the intersection of FZ2-base along face
!            ZMAX = -(ZMAX_0/2._WP/FZ2)*FZ1 + ZMAX_0/2._WP

            DO K=KMIN,KMAX
               Z=REAL(K,KIND=WP)+ZOFF
               IF(abs(Z) <= ZMAX) THEN
                  write(20,'(10(e12.5,1x))') x,y,z,zmax,fz1,fz2,fz3


! Site is occupied:
                  NAT=NAT+1
                  IF(NAT>MXNAT)THEN
                     CALL ERRMSG('FATAL','TAROCT',' NAT.GT.MXNAT ')
                  ENDIF
                  IXYZ(NAT,1)=I
                  IXYZ(NAT,2)=J
                  IXYZ(NAT,3)=K
               ENDIF
            ENDDO
         ENDDO
      ENDDO

! Initialize composition:
! Set X0 so that origin of TF is at centroid

      DO K=1,3
         X0(K)=0._WP
         DO I=1,NAT
            X0(K)=X0(K)+REAL(IXYZ(I,K))
            ICOMP(I,K)=1
         ENDDO
         X0(K)=-X0(K)/REAL(NAT)
      ENDDO

!*** Specify vectors A1 and A2 which will define target orientation
!    A1 is initially along x-axis
!    A2 is initially along y-axis

      A1(1)=1._WP
      A1(2)=0._WP
      A1(3)=0._WP
      A2(1)=0._WP
      A2(2)=1._WP
      A2(3)=0._WP

      WRITE(CDESCR,FMT='(A,I7,A)')' Octahedron of NAT=',NAT,' dipoles'
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)S,NAT,A1,A2,DX,X0
         DO JX=1,NAT
            WRITE(IOSHP,FMT=9030)JX,IXYZ(JX,1),IXYZ(JX,2),IXYZ(JX,3), &
                                 ICOMP(JX,1),ICOMP(JX,2),ICOMP(JX,3)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF
      RETURN
9020  FORMAT(' >TAROCT  octahedral grain: S=',F9.4,/,             &
         I10,' = NAT ',/,                                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')
9030  FORMAT (I7,3I4,3I2)
    END SUBROUTINE TAROCT
