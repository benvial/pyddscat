    SUBROUTINE EVALQ(CXADIA,CXAOFF,AK,NAT3,E02,CXE,CXP,CABS,CEXT,CPHA,MXN3, &
                     IMETHD)
      USE DDPRECISION,ONLY: WP
      IMPLICIT NONE

!*** Arguments:

      REAL(WP) :: CABS,CEXT,CPHA,E02
      INTEGER :: IMETHD,MXN3,NAT3
      COMPLEX(WP) ::   &
         CXADIA(MXN3), &
         CXAOFF(MXN3), &
         CXE(MXN3),    &
         CXP(MXN3)
      REAL(WP) :: &
         AK(3)

!*** Local variables:

      COMPLEX(WP) :: CXA,CXI,DCXA,RABS
      REAL(WP) :: AK1,AK2,AK3,PI
      INTEGER :: J1,J2,J3,NAT

!*** Intrinsic functions:

      INTRINSIC AIMAG,CONJG,REAL,SQRT

!*** SAVE statements:

      SAVE CXI

!*** Data statements:

      DATA CXI/(0._WP,1._WP)/

!***********************************************************************

! Given: 
!        MXN3=dimension parameter for CXADIAG,CXAOFF,CXE,CXP
!        CXADIA(J,1-3)=(a_11,a_22,a_33) for dipoles J=1-NAT,
!                  where symmetric 3x3 matrix a_ij is inverse of complex
!                  polarizability tensor alpha_ij for dipole J
!        CXAOFF(J,1-3)=(a_23,a_31,a_12) for dipoles J=1-NAT
!        AK(1-3) = (k_x,k_y,k_z)*d
!                  d = (d_x*d_y*d_z)**(1/3) = effective lattice spacing
!        NAT3 = 3*number of dipoles
!        E02 = |E_0|^2 , where E_0 = incident complex E field.
!        CXE(1-NAT3) = components of E field at each dipole, in order
!                      E_1x,E_2x,...,E_NATx,E_1y,E_2y,...,E_NATy,
!                      E_1z,E_2z,...,E_NATz
!        CXP(1-NAT3) = components of polarization vector at each dipole,
!                      in order
!                      P_1x,P_2x,...,P_NATx,P_1y,P_2y,...,P_NATy,
!                      P_1z,P_2z,...,P_NATz
!        IMETHD = 0 or 1
! Finds:
!        CEXT = extinction cross section /d**2
!  and, if IMETHD=1, also computes
!        CPHA = phase-lag cross section /d**2
!        CABS = absorption cross section /d**2

! Note: present definition of CPHA is 1/2 of Martin's definition

! B.T.Draine, Princeton Univ. Obs., 87/1/4

! History:
! 88.04.28 (BTD): modifications
! 90.11.02 (BTD): modified to allow use of vacuum sites (now pass E02
!                 from calling routine instead of evaluating it here)
! 90.12.13 (BTD): modified to use IMETHD flag, to allow "fast" calls
!                 in which only CEXT is computed.
! 97.12.26 (BTD): removed CXALPH from argument list; replaced with
!                 CXADIA and CXAOFF.
!                 CXADIA and CXAOFF are diagonal and off-diagonal
!                 elements of alpha^{-1} for each dipole.
!                 Modified to properly evaluate CABS
! 98.01.01 (BTD): Correct inconsistencies in assumed data ordering.
! 98.01.13 (BTD): Examine for possible error in evaluation of Qabs
! 98.04.27 (BTD): Minor polishing.
! 08.01.13 (BTD): cosmetic changes to f90 version
! End history

! Copyright (C) 1993,1997,1998,2008 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
! Initialization:
! Compute magnitude of kd

      CEXT=0._WP
      CABS=0._WP
      CPHA=0._WP
      PI=4._WP*ATAN(1._WP)
      AK2=0._WP
      DO J1=1,3
         AK2=AK2+AK(J1)*AK(J1)
      ENDDO
      AK1=SQRT(AK2)
      AK3=AK1*AK2

! Initialization complete.

      IF(IMETHD==0)THEN

!*** Compute CEXT:

         DO J1=1,NAT3
            CEXT=CEXT+AIMAG(CXP(J1))*REAL(CXE(J1))-REAL(CXP(J1))*AIMAG(CXE(J1))
         ENDDO
         CEXT=4._WP*PI*AK1*CEXT/E02
      ELSEIF(IMETHD==1)THEN

! Compute
! C_abs=(4*pi*k/|E_0|^2)*
!       sum_J { Im[P_J*conjg(a_J*P_J)] - (2/3)*k^3*|P_J|^2 }
!      =(4*pi*k/|E_0|^2)*
!       sum_J {-Im[conjg(P_J)*a_J*P_J] - Im[i*(2/3)*k^3*P_J*conjg(P_J)]}
!      =-(4*pi*k/|E_0|^2)*
!       Im{ sum_J [ conjg(P_J)*(a_J*P_J + i*(2/3)*k^3*P_J) ] }

         NAT=NAT3/3
         CXA=(0._WP,0._WP)
         RABS=CXI*AK3/1.5_WP
         DO J1=1,NAT
            J2=J1+NAT
            J3=J2+NAT
            DCXA=CONJG(CXP(J1))*((CXADIA(J1)+RABS)*CXP(J1)+              &
                                 CXAOFF(J2)*CXP(J3)+CXAOFF(J3)*CXP(J2))+ &
                 CONJG(CXP(J2))*((CXADIA(J2)+RABS)*CXP(J2)+              &
                                 CXAOFF(J3)*CXP(J1)+CXAOFF(J1)*CXP(J3))+ &
                 CONJG(CXP(J3))*((CXADIA(J3)+RABS)*CXP(J3)+              &
                                 CXAOFF(J1)*CXP(J2)+CXAOFF(J2)*CXP(J1))
            CXA=CXA+DCXA
         ENDDO
         CABS=-4._WP*PI*AK1*AIMAG(CXA)/E02
! Compute CEXT and CPHA
         CXA=(0._WP,0._WP)
         DO J1=1,NAT3
            CXA=CXA+CXP(J1)*CONJG(CXE(J1))
         ENDDO
         CEXT=4._WP*PI*AK1*AIMAG(CXA)/E02
         CPHA=2._WP*PI*AK1*REAL(CXA)/E02
      ENDIF
      RETURN

    END SUBROUTINE EVALQ
