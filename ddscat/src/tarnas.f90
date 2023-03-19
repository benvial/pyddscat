    SUBROUTINE TARNAS(A1,A2,DIAMX,PRINAX,DX,X0,BETADF,PHIDF,THETADF,   &
                      CFLSHP,CDESCR,IDVSHP,IOSHP,MXNAT,NAT,IXYZ,ICOMP, &
                      NCOMP_NEED,IANISO)
      USE DDPRECISION,ONLY : WP
      IMPLICIT NONE
!-------------------------------- tarnas v4 --------------------------------
!** Arguments:

      CHARACTER :: CDESCR*67,CFLSHP*80
      INTEGER :: IANISO,IDVSHP,IOSHP,MXNAT,NAT,NCOMP_NEED
      INTEGER*2 ::    &
         ICOMP(MXNAT,3)
      INTEGER ::     &
         IXYZ(MXNAT,3)
      REAL(WP) :: DIAMX,PRINAX
      REAL(WP) ::        &
         A1(3),          &
         A2(3),          &
         BETADF(MXNAT),  &
         DX(3),          &
         PHIDF(MXNAT),   &
         THETADF(MXNAT), &
         X0(3)

!** Local variables:

      INTEGER :: NSPHMX
      PARAMETER(NSPHMX=10000)
      CHARACTER :: CMSGNM*70
      LOGICAL :: OCC
      INTEGER :: JA,JS,JX,JY,JZ,LMX1,LMX2,LMY1,LMY2,LMZ1,LMZ2,NSPH
      INTEGER ::      &
         IC1(NSPHMX), &
         IC2(NSPHMX), &
         IC3(NSPHMX)
      REAL(WP) :: R2,SCALE,X,XMAX,XMIN,Y,YMAX,YMIN,Z,ZMAX,ZMIN
      REAL(WP) ::     &
         ALPHA(3),    &
         AS(NSPHMX),  &
         AS2(NSPHMX), &
         BE(NSPHMX),  &
         PH(NSPHMX),  &
         TH(NSPHMX),  &
         XS(NSPHMX),  &
         YS(NSPHMX),  &
         ZS(NSPHMX)

!***********************************************************************
! Routine to construct multisphere target of general anisotropic
! materials (target option 'NANSPH')

! Input:
!        DIAMX =max extent in X direction/d
!        PRINAX=0 to set A1=(1,0,0), A2=(0,1,0)
!              =1 to use principal axes for vectors A1 and A2
!                 as found here for dipole realization
!              =2 to use axes a1 and a2 read in from CFLSHP
!        DX(1-3)=(dx,dy,dz)/d where dx,dy,dz=lattice spacings in x,y,z
!                directions, and d=(dx*dy*dz)**(1/3)=effective lattice
!                spacing
!        CFLSHP= name of file containing locations and radii of spheres
!        IOSHP =device number for "target.out" file
!              =-1 to suppress printing of "target.out"
!        MXNAT =dimensioning information (max number of atoms)

! and, from input file CFLSHP:

!        N              = number of spheres
!        descriptive line which will be ignored (may be blank)
!        descriptive line which will be ignored (may be blank)
!        descriptive line which will be ignored (may be blank)
!        descriptive line which will be ignored (may be blank)
!        XS(1) YS(1) ZS(1) AS(1) IC1(1) IC2(1) IC3(1) TH(1) PH(1) BE(1)
!        XS(2) YS(2) ZS(2) AS(2) IC1(2) IC2(2) IC3(2) TH(2) PH(2) BE(2)
!        ...
!        XS(N) YS(N) ZS(N) AS(N) IC1(N) IC2(N) IC3(N) TH(N) PH(N) BE(N)

! where XS(J),YS(J),ZS(J) = x,y,z coordinates in target frame, in
!                                 arbitrary units, of center of sphere J
!        AS(J)            = radius of sphere J (same arbitrary units)
!        IC1(J),IC2(J),IC3(J) = dielectric function identifier
!                               for directions 1,2,3 in DF
!        TH(J),PH(J),BE(J) = angles theta_epsilon,phi_epsilon,beta_epsil
!                            specifying orientation of DF in TF

! units used for XS,YS,ZS, and AS are really arbitrary: actual size of
! overall target, in units of lattice spacing d, is controlled by
! parameter SHPAR(1) in ddscat.par

! Output:
!        A1(1-3)=(1,0,0)=unit vector defining target axis 1 in Target Frame
!        A2(1-3)=(0,1,0)=unit vector defining target axis 2 in Target Frame
!        X0(1-3)=location/d in TF of lattice site IXYZ=0 0 0
!                TF origin is taken to be located at centroid of N spheres.
!        NAT=number of atoms in target
!        IXYZ(1-NAT,1-3)=x/d,y/d,z/d for atoms of target
!        ICOMP(1-NAT,1-3)=composition identifier
!        IANISO=0 if all spheres are composed of isotropic materials
!               1 if some spheres are composed of anisotropic materials
!                 but with all optical axes parallel to x_TF,y_TF,z_TF
!                 [more specifically, BETADF=PHIDF=THETADF=0]
!               2 if one or more spheres is composed of anisotropic
!                 material with optical axes not aligned with x_TF,y_TF,z_TF
!        NCOMP_NEED=number of dielectric files needed
!        BETADF(1-NAT)=angle beta_epsilon (radians) specifying orientation
!                      of Dielectric Frame (DF) rel to Target Frame (TF)
!                      See UserGuide for definition of angles theta_epsilon
!                      phi_epsilon, beta_epsilon and discussion.
!        PHIDF(1-NAT)=angle phi_epsilon (radians) specifying orientation
!                     DF relative to TF
!        THETADF(1-NAT)=angle theta_epsilon (radians) specifying
!                       orientation of DF relative to TF
!        CDESCR=description of target (up to 67 characters)

! B.T.Draine, Princeton Univ. Obs., 2004.09.15

! History: 
! 04.09.15 (BTD) adapted from subroutine TARNSP 
! 04.09.16 (BTD) further work 
! 07.09.11 (BTD) changed IXYZ from INTEGER*2 to INTEGER 
! 08.08.08 (BTD) added X0 to argument list 
!                added code to set X0 so that TF origin is located at 
!                volume-weighted centroid of the N spheres.  
! 08.08.30 (BTD) modified format 9020 
! 08.09.17 (BTD) v2: fixed bugs 
!                target centroid and therefore X0 was not being 
!                computed correctly in v1: now corrected 
! 09.10.11 (BTD) v3, v7.1.0 
!                * changed number of unused descriptor lines
!                  from 1 to 4
!                  to provide compatibility with output of
!                  programs agglom.f and assign_comp.f
! 10.01.06 (BTD) * added sanity checks to detect and report input file
!                  incompatibilities
! 10.02.06 (BTD) * added rotation angles for dielectric tensor to 
!                  output file target.out
! 14.01.03 (BTD) * amended write statement to target.out to make consistent
!                  with UserGuide
!                  (error noted by Martin Ernesto Figuroa Delgadillo)
!                * changed NSPHMX from 4096 to 10000 to allow larger targets
! 14.10.15 (BTD) * modified to set alpha(1-3)=0 when called with
!                  PRINAX=0
! 14.10.26 (BTD) v4, v7.3.1
!                * added NCOMP_NEED and IANISO to argument list
!                * added code to determine NCOMP_NEED and IANISO
! 15.03.19 (BTD) * added code to support option PRINAX=2
!                  whereby a1,a2 are read in from target file
! end history

! Copyright (C) 2004,2007,2008,2009,2010,2014,2015
! B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 1'
!***

! Read parameters of spheres:

      WRITE(0,FMT='(A,A)')'>TARNAS_v4 open file=',CFLSHP
!--------------------------------------
      OPEN(UNIT=IDVSHP,FILE=CFLSHP)
      READ(IDVSHP,*,ERR=8100)NSPH
      IF(NSPH>NSPHMX)THEN
         WRITE(CMSGNM,FMT='(A,I6,A,I6,A)')'fatal error: NSPH=',NSPH, &
                                        ' > NSPHMX',NSPHMX
         CALL WRIMSG('TARNAS',CMSGNM)
         CALL ERRMSG('FATAL','TARNAS',' NSPH > NSPHMX')
         STOP
      ENDIF
      IF(NINT(PRINAX)<2)THEN
! skip 4 lines
         DO JA=1,4
            READ(IDVSHP,*)
         ENDDO
      ELSE
         READ(IDVSHP,*)
         READ(IDVSHP,*)A1
         READ(IDVSHP,*)A2
         READ(IDVSHP,*)
      ENDIF
      DO JA=1,NSPH
         READ(IDVSHP,*,ERR=8200,END=8300)XS(JA),YS(JA),ZS(JA),AS(JA), &
                       IC1(JA),IC2(JA),IC3(JA),TH(JA),PH(JA),BE(JA)
      ENDDO

!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 2, NSPH=',NSPH
!***

      CLOSE(IDVSHP)
      WRITE(0,FMT='(A,A)')'>TARNAS close file=',CFLSHP

! Determine number of compositions
!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 3, NSPH=',NSPH
!***
      NCOMP_NEED=0
!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 3.1, NSPH=',NSPH
!***
!      IANISO=0
!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 3.2, NSPH=',NSPH
!***
      IANISO=0
      DO JA=1,NSPH
!*** diagnostic
!         write(0,*)'tarnas_v4 ckpt 3.3: JA=',JA
!***
         IF(IC1(JA).GT.NCOMP_NEED)NCOMP_NEED=IC1(JA)
!*** diagnostic
!         write(0,*)'tarnas_v4 ckpt 3.4: NCOMP_NEED=',NCOMP_NEED
!***
         IF(IC2(JA).GT.NCOMP_NEED)NCOMP_NEED=IC2(JA)
         IF(IC3(JA).GT.NCOMP_NEED)NCOMP_NEED=IC3(JA)
         IF(IC1(JA).NE.IC2(JA))IANISO=1
         IF(IC2(JA).NE.IC3(JA))IANISO=1
      ENDDO

!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 4'
!      write(0,*)'  NCOMP_NEED=',NCOMP_NEED,' IANISO=',IANISO
!***
! Determine value of IANISO

      IF(IANISO==1)THEN
         DO JA=1,NSPH
            IF((TH(JA)/=0.).OR.(PH(JA)/=0.).OR.(TH(JA)/=0.))IANISO=2
         ENDDO
      ENDIF
!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 5: IANISO=',IANISO
!***
! Check for overlap and, if so, possible conflict in composition
! or (if anisotropic) orientation

      IF(NSPH>=2)THEN
         DO JA=2,NSPH
            DO JY=1,JA-1
               R2=(XS(JY)-XS(JA))**2+(YS(JY)-YS(JA))**2 + &
                  (ZS(JY)-ZS(JA))**2-(AS(JY)+AS(JA))**2

! exact R2 should be nonnegative, but allow for roundoff error

               IF(R2<-0.0001_WP*(AS(JY)+AS(JA)))THEN
                  IF(IC1(JY)/=IC1(JA).OR.IC2(JY)/=IC2(JA).OR. &
                     IC3(JY)/=IC3(JA))THEN
!*** diagnostic
!      write(0,*)'in tarnas, apparent overlap of spheres JY,JA=',jy,ja
!      write(0,*)'jy: x,y,z=',xs(jy),ys(jy),zs(jy)
!      write(0,*)'ja: x,y,z=',xs(ja),ys(ja),zs(ja)
!      write(0,*)'a(jy),a(ja)=',as(jy),as(ja)
!***
                     CALL ERRMSG('FATAL','TARNAS',                          &
                                 ' Sphere overlap but differing composition')
                  ENDIF
                  IF(IC1(JY)/=IC2(JY).OR.IC1(JY)/=IC3(JY))THEN
                     IF(BE(JY)/=BE(JA).OR.PH(JY)/=PH(JA).OR.TH(JY)/=TH(JA)) &
                        THEN
                        CALL ERRMSG('FATAL','TARNAS',                      &
                                ' Sphere overlap but differing orientation')
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF

! Dipoles are located at sites
! (x,y,z)=(I,J,K), I,J,K=integers

! Determine max extent in X direction:

      XMIN=XS(1)-AS(1)
      XMAX=XS(1)+AS(1)
      YMIN=YS(1)-AS(1)
      YMAX=YS(1)+AS(1)
      ZMIN=ZS(1)-AS(1)
      ZMAX=ZS(1)+AS(1)
      IF(NSPH>1)THEN
         DO JA=2,NSPH
            IF(XS(JA)-AS(JA)<XMIN)XMIN=XS(JA)-AS(JA)
            IF(XS(JA)+AS(JA)>XMAX)XMAX=XS(JA)+AS(JA)
            IF(YS(JA)-AS(JA)<YMIN)YMIN=YS(JA)-AS(JA)
            IF(YS(JA)+AS(JA)>YMAX)YMAX=YS(JA)+AS(JA)
            IF(ZS(JA)-AS(JA)<ZMIN)ZMIN=ZS(JA)-AS(JA)
            IF(ZS(JA)+AS(JA)>ZMAX)ZMAX=ZS(JA)+AS(JA)
         ENDDO
      ENDIF
      SCALE=DIAMX/(XMAX-XMIN)

! Now determine min,max values of I,J,K:

      LMX1=NINT(SCALE*XMIN/DX(1)-0.01_WP)
      LMX2=NINT(SCALE*XMAX/DX(1)+0.01_WP)
      LMY1=NINT(SCALE*YMIN/DX(2)-0.01_WP)
      LMY2=NINT(SCALE*YMAX/DX(2)+0.01_WP)
      LMZ1=NINT(SCALE*ZMIN/DX(3)-0.01_WP)
      LMZ2=NINT(SCALE*ZMAX/DX(3)+0.01_WP)

      DO JA=1,NSPH
         AS2(JA)=(SCALE*AS(JA))**2
         XS(JA)=SCALE*XS(JA)
         YS(JA)=SCALE*YS(JA)
         ZS(JA)=SCALE*ZS(JA)
      ENDDO

! Determine list of occupied sites

      NAT=0
      DO JZ=LMZ1,LMZ2
         Z=REAL(JZ,KIND=WP)*DX(3)
         DO JY=LMY1,LMY2
            Y=REAL(JY,KIND=WP)*DX(2)
            DO JX=LMX1,LMX2
               X=REAL(JX,KIND=WP)*DX(1)
               OCC=.FALSE.
               DO JA=1,NSPH
                  R2=(X-XS(JA))**2+(Y-YS(JA))**2+(Z-ZS(JA))**2
                  IF(R2<AS2(JA))THEN
                     OCC=.TRUE.
                     JS=JA
                  ENDIF
               ENDDO
               IF(OCC)THEN

! Site is occupied: assign location IXYZ, composition ICOMP, and
! material orientation BETADF,PHIDF,THETADF:

                  NAT=NAT+1
                  IXYZ(NAT,1)=JX
                  IXYZ(NAT,2)=JY
                  IXYZ(NAT,3)=JZ
                  ICOMP(NAT,1)=IC1(JS)
                  ICOMP(NAT,2)=IC2(JS)
                  ICOMP(NAT,3)=IC3(JS)
                  BETADF(NAT)=BE(JS)
                  PHIDF(NAT)=PH(JS)
                  THETADF(NAT)=TH(JS)
               ENDIF
            ENDDO
         ENDDO
      ENDDO

      IF(NAT>MXNAT)THEN
         CALL ERRMSG('FATAL','TARNAS',' NAT.GT.MXNAT ')
      ENDIF

! Specify target axes A1 and A2
! If PRINAX=0, then
!     A1=(1,0,0) in target frame
!     A2=(0,1,0) in target frame
!     set alpha(1-3)=0.
! If PRINAX=1., then
!     A1,A2 are principal axes of largest, second largest moment of
!     inertia
!     alpha(j)=I_j/(0.4*M*aeff^2) for principal axes j=1-3
! If PRINAX=2, then a1,a2 are already set from input file
      IF(NINT(PRINAX)<=0)THEN
         DO JX=1,3
            A1(JX)=0._WP
            A2(JX)=0._WP
            ALPHA(JX)=0._WP
         ENDDO
         A1(1)=1._WP
         A2(2)=1._WP
      ELSEIF(NINT(PRINAX)==1)THEN
!*** diagnostic
!         write(0,*)'tarnas_v4 ckpt 6'
!***
         CALL PRINAXIS(MXNAT,NAT,ICOMP,IXYZ,DX,A1,A2,ALPHA)
!*** diagnostic
!         write(0,*)'tarnas_v4 ckpt 7'
!         write(0,*)'  alpha(1-3)=',alpha
!***
      ELSEIF(NINT(PRINAX)>=2)THEN
!*** diagnostic
!         write(0,*)'tarnas_v4 ckpt 8: prinax=2'
!***
         DO JX=1,3
            ALPHA(JX)=0._WP
         ENDDO
      ENDIF

! Find volume-weighted centroid of the N spheres in IJK space
! Take this to be TF origin.
! Set X0 = -(IJK centroid)

      DO JX=1,3
         X0(JX)=0._WP
      ENDDO
      Z=0._WP

      DO JA=1,NSPH
         Z=Z+AS(JA)**3
      ENDDO

      DO JA=1,NSPH
         Y=AS(JA)**3/Z
         X0(1)=X0(1)-XS(JA)*Y
         X0(2)=X0(2)-YS(JA)*Y
         X0(3)=X0(3)-ZS(JA)*Y
      ENDDO

!*** diagnostic
!      write(0,*)'tarnas_v4 ckpt 8, x0=',x0
!***

!***********************************************************************
! Write target description into string CDESCR

      WRITE(CDESCR,FMT='(A,I7,A)')' Multisphere cluster containing',NAT, &
                                  ' dipoles'

!***********************************************************************

      CMSGNM=CDESCR
      CALL WRIMSG('TARNAS',CMSGNM)
      IF(IOSHP>=0)THEN
         OPEN(UNIT=IOSHP,FILE='target.out',STATUS='UNKNOWN')
         WRITE(IOSHP,FMT=9020)NSPH,DIAMX,NAT,ALPHA,A1,A2,DX,X0
         DO JA=1,NAT
            WRITE(IOSHP,FMT=9030)JA,IXYZ(JA,1),IXYZ(JA,2),IXYZ(JA,3), &
                                 ICOMP(JA,1),ICOMP(JA,2),ICOMP(JA,3), &
                                 THETADF(JA),PHIDF(JA),BETADF(JA)
! 2014.01.04 (BTD) corrected error to make consistent with UserGuide
!                  error noted by Martin Ernesto Figueroa Delgadillo 
!                  incorrect line follows:
!                                 BETADF(JA),PHIDF(JA),THETADF(JA)
         ENDDO
         CLOSE(UNIT=IOSHP)
      ENDIF

      RETURN
8100  WRITE(CMSGNM,FMT='(2A)')' Fatal error: unable to read number', &
                              ' N of spheres'
      CALL WRIMSG('TARNAS',CMSGNM)
      WRITE(CMSGNM,FMT='(2A)')' from file=',CFLSHP
      CALL WRIMSG('TARNAS',CMSGNM)
      STOP
8200  WRITE(CMSGNM,FMT='(2A,I7)')' Fatal error: error reading data for', &
                                 ' sphere',JA
      CALL WRIMSG('TARNAS',CMSGNM)
      STOP
8300  WRITE(CMSGNM,FMT='(A,I7,A)')' Fatal error: expected to read data for', &
                                  NSPH,' spheres'
      CALL WRIMSG('TARNAS',CMSGNM)
      WRITE(CMSGNM,FMT='(A,I7,A)')' but only found data for',JA-1,' spheres'
      CALL WRIMSG('TARNAS',CMSGNM)
      STOP

9020  FORMAT(' multisphere target composed of ',I3,                &
             ' spheres, DIAMX=',F8.4,/,                            &
         I10,3F8.4,' = NAT, alpha_1-3',/,                          &
         3F10.6,' = A_1 vector',/,                                 &
         3F10.6,' = A_2 vector',/,                                 &
         3F10.6,' = lattice spacings (d_x,d_y,d_z)/d',/,           &
         3F10.5,' = lattice offset x0(1-3) = (x_TF,y_TF,z_TF)/d ', &
               'for dipole 0 0 0',/,                               &
         '     JA  IX  IY  IZ ICOMP(x,y,z)')

9030  FORMAT(I10,3I10,3I2,3F10.6)
    END SUBROUTINE TARNAS
