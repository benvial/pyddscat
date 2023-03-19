    SUBROUTINE EVALE(CXE00,AKD,DX,X0,IXYZ0,MXNAT,MXN3,NAT,NAT0,NX,NY,NZ,CXE)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

!*** Arguments:

      INTEGER :: MXN3,MXNAT,NAT,NAT0,NX,NY,NZ
      INTEGER ::     &
         IXYZ0(NAT0,3)
      REAL(WP) :: &
         AKD(3),  &
         DX(3),   &
         X0(3)

! Note: CXE should be dimensioned to CXE(NAT,3) in this routine
!       so that first 3*NAT elements of CXE are employed.
!       XYZ0 should be dimensioned to

      COMPLEX(WP) :: & 
         CXE(NAT,3), &
         CXE00(3)

!***  Local variables:

      COMPLEX(WP) :: CXFAC,CXI
      REAL(WP) :: X,X1,X2
      INTEGER :: IA,IX,IY,IZ,M

!*** Intrinsic functions:

      INTRINSIC EXP,REAL

!***********************************************************************
! subroutine EVALE

! Given:   CXE00(1-3)=Incident E field at origin (complex) at t=0
!          AKD(1-3)=(kx,ky,kz)*d for incident wave (d=effective
!                    lattice spacing)
!          DX(1-3)=(dx/d,dy/d,dz/d) for lattice (dx,dy,dz=lattice
!                   spacings in x,y,z directions, d=(dx*dy*dz)**(1/3)
!          X0(1-3)=(x,y,z)location/(d*DX(1-3)) in TF of lattice site
!                  with IX=0,IY=0,IZ=0
!          IXYZ0(1-NAT0,3)=[x-x0(1)]/dx, [y-x0(2)]/dy, [z-x0(3)]/dz
!                  for each of NAT0 physical dipoles
!          MXNAT,MXN3=dimensioning information
!          NAT0=number of dipoles in physical target
!          NAT=number of locations at which to calculate CXE

! Returns: CXE(1-NAT,3)=incident E field at NAT locations at t=0

! B.T.Draine, Princeton Univ. Obs., 88.05.09
! History:
! 90.11.06 (BTD): Modified to pass array dimension.
! 90.11.29 (BTD): Modified to allow option for either
!                   physical locations only (NAT=NAT0), or
!                   extended dipole array (NAT>NAT0)
! 90.11.30 (BTD): Corrected error for case NAT>NAT0
! 90.11.30 (BTD): Corrected another error for case NAT>NAT0
! 90.12.03 (BTD): Change ordering of XYZ0 and CXE
! 90.12.05 (BTD): Corrected error in dimensioning of CXE
! 90.12.10 (BTD): Remove XYZ0, replace with IXYZ0
! 97.11.02 (BTD): Add DX to argument list to allow use with
!                 noncubic lattices.
! 07.06.20 (BTD): Add X0 to the argument list to specify location
!                 in TF corresponding to IX=0,IY=0,IZ=0
! 07.09.11 (BTD): Changed IXYZ0 from INTEGER*2 to INTEGER
! 08.03.14 (BTD): v7.05
!                 corrected dimensioning
!                 IXYZ0(MXNAT,3) -> IXYZO(NAT0,3)
! Copyright (C) 1993,1997,2007 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

      CXI=(0._WP,1._WP)

! Evaluate electric field vector at each dipole location.

! If NAT=NAT0, then evaluate E only at occupied sites.
! If NAT>NAT0, then evaluate E at all sites.

      IF(NAT==NAT0)THEN
         DO IA=1,NAT0
            X=0._WP
            DO M=1,3
               X=X+AKD(M)*DX(M)*(REAL(IXYZ0(IA,M),KIND=WP)+X0(M))
            ENDDO
            CXFAC=EXP(CXI*X)
            DO M=1,3
               CXE(IA,M)=CXE00(M)*CXFAC
            ENDDO
         ENDDO
      ELSE
         IA=0
         DO IZ=1,NZ
            X1=AKD(3)*DX(3)*(REAL(IZ,KIND=WP)+X0(3))
            DO IY=1,NY
               X2=X1+AKD(2)*DX(2)*(REAL(IY,KIND=WP)+X0(2))
               DO IX=1,NX
                  IA=IA+1
                  X=X2+AKD(1)*DX(1)*(REAL(IX,KIND=WP)+X0(1))
                  CXFAC=EXP(CXI*X)
                  DO M=1,3
                     CXE(IA,M)=CXE00(M)*CXFAC
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      RETURN
    END SUBROUTINE EVALE
