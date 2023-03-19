    SUBROUTINE ROTATE(A1_TF,A2_TF,AK1,CXE01_LF,CXE02_LF,BETA,THETA,PHI,     &
         EN0_TF,CXE01_TF,CXE02_TF,MXSCA,MYID,NSCAT,ENSC_LF,EM1_LF,EM2_LF,   &
         AKS_TF,EM1_TF,EM2_TF,XLF_TF,YLF_TF,ZLF_TF)                         !
! -------------------------------- v5 ---------------------------------------
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: MXSCA,MYID,NSCAT
      REAL(WP) :: AK1,BETA,PHI,THETA
      REAL(WP) ::          &
         A1_TF(3),         &
         A2_TF(3),         &
         AKS_TF(3,MXSCA),  &
         EM1_LF(3,MXSCA),  &
         EM1_TF(3,MXSCA),  &
         EM2_LF(3,MXSCA),  &
         EM2_TF(3,MXSCA),  &
         EN0_TF(3),        &
         ENSC_LF(3,MXSCA), &
         XLF_TF(3),        &
         YLF_TF(3),        &
         ZLF_TF(3)         !
      COMPLEX(WP) ::  &
         CXE01_LF(3), &   
         CXE02_LF(3), &
         CXE01_TF(3), &
         CXE02_TF(3)  !

! Local variables

      CHARACTER :: CMSGNM*70
      INTEGER :: J
      REAL(WP) :: COSBET,COSPHI,COSTHE,DENOM,SINBET,SINPHI,SINTHE,SUM
      REAL(WP) ::  & 
         A3_TF(3), &
         RM(3,3),  &
         VEC(3)    !

!***********************************************************************
! Given:
!       A1_TF(1-3)=target axis 1 in TF=Target Frame
!       A2_TF(1-3)=target axis 2 in TF
!       AK1=magnitude of k vector = k*d
!       CXE01_LF(1-3)=incident polarization vector 1 in LF=Lab Frame
!       CXE02_LF(1-3)=incident polarization vector 2 in LF
!       BETA=angle (radians) through which target is to be rotated
!            around target axis A1
!       THETA        =angle (radians) between A1 and xLab
!       PHI=angle (radians) between xLab-A1 plane and xLab,yLab plane
!       ENSC_LF(1-3,1-NSCAT)=scattering directions in LF
!       EM1_LF(1-3,1-NSCAT)=scattering polarization vector 1 in LF
!       EM2_LF(1-3,1-NSCAT)=scattering polarization vector 2 in LF
!       MYID         =MPI process ID (=0 for master process)
! Returns:
!       EN0_TF(1-3)  =incident propagation vector in TF
!       CXE01_TF(1-3)=incident polariz. vector 1 in TF
!       CXE02_TF(1-3)=incident polariz. vector 2 in TF
!       AKS_TF(1-3,1-NSCAT)=scattering k vectors in TF
!       EM1_TF(1-3,1-NSCAT)=scattering pol. vector 1 in TF
!       EM2_TF(1-3,1-NSCAT)=scattering pol. vector 2 in TF

! Purpose:
!       In the Lab Frame, we hold the propagation direction (xLab-axis)
!       and incident polarization vectors (CXE01_LF and CXE02_LF) fixed,
!       and rotate target to an orientation specified by the three
!       angles BETA,THETA,PHI.
!       Computations in the main program DDSCAT are carried out
!       in the Target Frame, in which the target is held fixed and the
!       propagation vectors and polarization vectors are rotated.
!       This routine finds the incident propagation vector,
!       scattering propagation vectors, and associated polarization
!       vectors in the Target Frame.

! B. T. Draine, Princeton Univ. Obs., 89.11.20
! History:
! 91.04.22 (BTD) Eliminated vector A1R and removed superfluous line
!                computing A1R
! 96.11.03 (BTD) Corrected comments
! 03.06.06 (BTD) Corrected calculation of BETA0,PHI0 -- it appears that
!                these were not correct when a_1z < 0
!                                     and/or a_3x < 0
! 05.03.29 (BTD) Add a few lines of code to guard against possibility
!                that numerical roundoff will cause argument of
!                ACOS to have absolute value > 1.
!                Declare new variable TERM
! 15.03.13 (BTD) v2
!                * change notation
!                  A1     -> A1_TF
!                  A2     -> A2_TF
!                  AKSR   -> AKS_TF
!                  CXE01  -> CXE01_LF
!                  CXE02  -> CXE02_LF
!                  CXE01R -> CXE01_TF
!                  CXE02R -> CXE02_TF
!                  EM1    -> EM1_LF
!                  EM2    -> EM2_LF
!                  EM1R   -> EM1_TF
!                  EM2R   -> EM2_TF
!                  ENSC   -> ENSC_LF
! 15.03.14 (BTD) * add code to handle case of THETA0 close to 0 or pi
! 15.03.15 (BTD) v3 rewrite code determining phi0 and beta0
! 15.03.17 (BTD) v4 further changes in implementation of rotations
! 15.04.04 (BTD) v5 add XLF_TF,YLF_TF,ZLF_TF to argument list
! 15.04.05 (BTD) * major rewrite of ROTATE
! 15.04.06 (BTD) * add MYID to argument list to control output 
! end history

! Copyright (C) 1996,2003,2005,2015 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************

! find target axis a3:

      A3_TF(1)=A1_TF(2)*A2_TF(3)-A1_TF(3)*A2_TF(2)
      A3_TF(2)=A1_TF(3)*A2_TF(1)-A1_TF(1)*A2_TF(3)
      A3_TF(3)=A1_TF(1)*A2_TF(2)-A1_TF(2)*A2_TF(1)

! 2015.04.05 (BTD) obtain rotation matrix RM by direct algebra instead 
!                  of series of rotation operations

! RM(1,1) = x_LF dot x_TF
! RM(2,1) = y_LF dot x_TF
! RM(3,1) = z_LF dot x_TF

! RM(1,2) = x_LF dot y_TF
! RM(2,2) = y_LF dot y_TF
! RM(3,3) = z_LF dot y_TF

! RM(1,3) = x_LF dot z_TF
! RM(2,3) = y_LF dot z_TF
! RM(3,3) = z_LF dot x_TF

      COSTHE=COS(THETA)
      SINTHE=SIN(THETA)
      COSBET=COS(BETA)
      SINBET=SIN(BETA)
      COSPHI=COS(PHI)
      SINPHI=SIN(PHI)
      DO J=1,3
         RM(J,1)=A1_TF(J)*COSTHE        &
                -A2_TF(J)*SINTHE*COSBET &
                +A3_TF(J)*SINTHE*SINBET !

         RM(J,2)=A1_TF(J)*SINTHE*COSPHI                        &
                +A2_TF(J)*(COSTHE*COSBET*COSPHI-SINBET*SINPHI) &
                -A3_TF(J)*(COSTHE*SINBET*COSPHI+COSBET*SINPHI) !

         RM(J,3)=A1_TF(J)*SINTHE*SINPHI                        &
                +A2_TF(J)*(COSTHE*COSBET*SINPHI+SINBET*COSPHI) &
                -A3_TF(J)*(COSTHE*SINBET*SINPHI-COSBET*COSPHI) !
      ENDDO

!**** Now rotate all required vectors:
! x_LF,y_LF,z_LF in TF

      VEC(1)=1._WP
      VEC(2)=0._WP
      VEC(3)=0._WP

      CALL PROD3(RM,VEC,XLF_TF)
      VEC(1)=0._WP
      VEC(2)=1._WP
      CALL PROD3(RM,VEC,YLF_TF)
      VEC(2)=0._WP
      VEC(3)=1._WP
      CALL PROD3(RM,VEC,ZLF_TF)

! Incident propagation vector:

      DO J=1,3
         EN0_TF(J)=XLF_TF(J)
      ENDDO

! Incident polarization vectors:

      CALL PROD3C(RM,CXE01_LF,CXE01_TF)
      CALL PROD3C(RM,CXE02_LF,CXE02_TF)

! Scattering vectors:

      CALL PROD3V(RM,ENSC_LF,AKS_TF,NSCAT,MXSCA)

      DO J=1,NSCAT
         AKS_TF(1,J)=AK1*AKS_TF(1,J)
         AKS_TF(2,J)=AK1*AKS_TF(2,J)
         AKS_TF(3,J)=AK1*AKS_TF(3,J)
      ENDDO

! Scattering polarization vectors:

      CALL PROD3V(RM,EM1_LF,EM1_TF,NSCAT,MXSCA)
      CALL PROD3V(RM,EM2_LF,EM2_TF,NSCAT,MXSCA)

      IF(MYID==0)THEN
!*** diagnostic
!         write(0,*)'rotate_v5 ckpt 1: rotations complete:'
!***
         WRITE(CMSGNM,FMT='(A,3F12.8)')'      a_1 in TF:',A1_TF
         CALL WRIMSG('ROTATE',CMSGNM)
         WRITE(CMSGNM,FMT='(A,3F12.8)')'      a_2 in TF:',A2_TF
         CALL WRIMSG('ROTATE',CMSGNM)
         WRITE(CMSGNM,FMT='(A,3F12.8)')'      a_3 in TF:',A3_TF
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,3F12.8)')'     x_LF in TF:',XLF_TF
         CALL WRIMSG('ROTATE',CMSGNM)
         WRITE(CMSGNM,FMT='(A,3F12.8)')'     y_LF in TF:',YLF_TF
         CALL WRIMSG('ROTATE',CMSGNM)
         WRITE(CMSGNM,FMT='(A,3F12.8)')'     z_LF in TF:',ZLF_TF
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'n0=x_LF dot a_1=',            &
            (XLF_TF(1)*A1_TF(1)+XLF_TF(2)*A1_TF(2)+XLF_TF(3)*A1_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'n0=x_LF dot a_2=',            &
            (XLF_TF(1)*A2_TF(1)+XLF_TF(2)*A2_TF(2)+XLF_TF(3)*A2_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'n0=x_LF dot a_3=',            &
            (XLF_TF(1)*A3_TF(1)+XLF_TF(2)*A3_TF(2)+XLF_TF(3)*A3_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'   y_LF dot a_1=',            &
            (YLF_TF(1)*A1_TF(1)+YLF_TF(2)*A1_TF(2)+YLF_TF(3)*A1_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'   y_LF dot a_2=',            &
            (YLF_TF(1)*A2_TF(1)+YLF_TF(2)*A2_TF(2)+YLF_TF(3)*A2_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'   y_LF dot a_3=',            &
            (YLF_TF(1)*A3_TF(1)+YLF_TF(2)*A3_TF(2)+YLF_TF(3)*A3_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'   z_LF dot a_1=',            &
            (ZLF_TF(1)*A1_TF(1)+ZLF_TF(2)*A1_TF(2)+ZLF_TF(3)*A1_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'   z_LF dot a_2=',            &
            (ZLF_TF(1)*A2_TF(1)+ZLF_TF(2)*A2_TF(2)+ZLF_TF(3)*A2_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)

         WRITE(CMSGNM,FMT='(A,F12.8)')'   z_LF dot a_3=',            &
            (ZLF_TF(1)*A3_TF(1)+ZLF_TF(2)*A3_TF(2)+ZLF_TF(3)*A3_TF(3))
         CALL WRIMSG('ROTATE',CMSGNM)
      ENDIF
      RETURN
    END SUBROUTINE ROTATE
