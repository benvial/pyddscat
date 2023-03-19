    SUBROUTINE PAD(CXA,NX,NY,NZ,CXB)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

! Pad the array CXA with zeros out to twice its length in each
! dimension, and put the result in CXB.

! Originally written by Jeremy Goodman, 90.09.22
! as part of Goodman & Draine 1991, Optics Letters 16, 1198
! history
! 90.09.22 (J.Goodman) first written 
! 90.11.29 (BTD) adapted to DDSCAT
! 14.10.30 (BTD) separated from eself_v8 to facilitate debugging
!                for some reason ifort -O2 under mac osx 10.9.5
!                on ourania4 is skipping the zero initialization loop
!                BEWARE: compilation with ifort -O2 will create
!                bad object code on some systems.
!                Observed for ifort version 14.0.1 on ourania4 running
!                osx 10.9.5
!
! Copyright (C) 1993,2014 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Arguments:
      INTEGER :: NX,NY,NZ
      COMPLEX(WP) ::       &
         CXA(NX,NY,NZ),    &
         CXB(2*NX,2*NY,2*NZ)
! Local variables:
      COMPLEX(WP) :: CXZERO
      INTEGER :: I,J,K
      SAVE CXZERO
      DATA CXZERO/(0._WP,0._WP)/
!***********************************************************************

#ifdef openmp
!$OMP PARALLEL        &
!$OMP&   PRIVATE(I,J,K)
!$OMP DO
#endif

! Caution: ifort -O2 may decide to skip this initialization loop...
!          [14.10.31: This is observed with ifort 14.0.1 
!                     under osx 10.9.5 on ourania4]
      DO K=1,2*NZ
         DO J=1,2*NY
            DO I=1,2*NX
               CXB(I,J,K)=CXZERO
            ENDDO
         ENDDO
      ENDDO

#ifdef openmp
!$OMP END DO
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
!$OMP END DO
!$OMP END PARALLEL
#endif

      RETURN
    END SUBROUTINE PAD
