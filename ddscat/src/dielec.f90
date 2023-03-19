    SUBROUTINE DIELEC(WAVE,IDVOUT,CFLEPS,CXEPS,MXCOMP,MXWAVT,NCOMP,E1A,E2A,WVA)
      USE DDPRECISION, ONLY : WP
      IMPLICIT NONE

! Arguments:

      INTEGER :: MXCOMP,MXWAVT
      CHARACTER(60) :: &
         CFLEPS(MXCOMP)
      COMPLEX(WP) :: &
         CXEPS(MXCOMP)
      INTEGER :: IDVOUT,NCOMP
      REAL(WP) :: WAVE
      REAL(WP) ::     &
         E1A(MXWAVT), &
         E2A(MXWAVT), &
         WVA(MXWAVT)

! Local Variables:

      LOGICAL :: INIT
      CHARACTER(70) :: CMSGNM
      CHARACTER(70) :: CDESCR
      COMPLEX(WP) :: CXI
      INTEGER :: I,ICEPS1,ICEPS2,ICIMN,ICOL,ICREN,ICWV, &
                 IEPS1,IEPS2,IIMN,IREN,IWV,J,JJ,NWAVT   !

      REAL(WP) :: DUM1,DUM2,E1,E2,TEMP
      REAL(WP) :: & 
         XIN(10)
      EXTERNAL ERRMSG
      SAVE CXI,ICEPS1,INIT,NWAVT
!**********************************************************************
! Given:
!       WAVE = wavelength (micron)
!       IDVOUT=output unit number
!       CFLEPS(1-NCOMP) = names of files containing dielectric data
!       MXCOMP,MXWAVT = dimensioning information
!       NCOMP = number of components
!       E1A, E2A, WVA = scratch arrays
! Returns:
!       CXEPS(1-NCOMP) = dielectric constant for materials
!                        1-NCOMP
! NOTE:
!       It is assumed that file(s) containing
!       the table(s) will have following format:
! line1 = description (CHARACTER*80)
! line2 = IWV, IREN, IIMN, IEPS1, IEPS2
!                IWV: if IWV > 0, wavelength is tabulated in col IWV
!                     if IWV < 0, energy (eV) is tabulated in col -IWV
!               IREN: if IREN = 0, Re(N) is not tabulated
!                     if IREN > 0: Re(N) is tabulated in col IREN
!                     if IREN < 0: Re(N-1) is tabulated in col -IREN
!               IIMN: if IIMN=0, Im(N) is not tabulated
!                     otherwise Im(N) is tabulated in col |IIMN|
!              IEPS1: if IEPS1 = 0: Re(EPS) is not tabulated
!                     if IEPS1 > 0: Re(EPS) is tabulated in col IEPS1
!                     if IEPS1 < 0: Re(EPS-1) is tabulaed in col -IEPS1
!              IEPS2: if IEPS2 = 0: Im(EPS) is not tabulated
!                     otherwise, Im(EPS) is tabulated in col |IEPS2|

! line3 = header line (will not be read)

! line4... = data, according to columns specified on line2
! wavelengths must be monotonic, either increasing or decreasing

! B.T.Draine, Princeton Univ. Obs.
! History:
! 90.12.04 (BTD): Rewritten to allow H2OICE and H2OLIQ options.
! 90.12.21 (BTD): Correct handling of 'TABLES' option when more than
!                 one wavelength is considered.
! 91.05.02 (BTD): Added IDVOUT to argument list for subr. INTERP
! 91.09.12 (BTD): Moved declaration of MXCOMP and MWAVT ahead of
!                 declaration of CFLEPS and CXEPS
! 93.12.06 (BTD): Added IEPS1,INIT, and NWAVT to SAVE statement
!                 (without this did not run properly on SGI Indigo)
! 96.12.16 (BTD): Corrected temperature to T=250K in output statement
!                 for H2OICE option
! 98.12.21 (BTD): changed dimension of CFLPAR from CHARACTER*40 to
!                 CHARACTER*60 to allow longer file names
!                 (also changed in reapar.f and DDSCAT.f)
! 04.10.14 (BTD): added check to ensure that user does not specify
!                 refractive index = 1 (this would cause division
!                 by zero elsewhere in code)
! 07.07.30 (PJF): Converted to f90
! 07.07.31 (BTD): Removed calls to REFWAT and REFICE
! 07.10.28 (BTD): Eliminated CDIEL -- reading from tables is standard
!                 Added output line just before table read as clue in
!                 case of failure during table read.
! 09.10.19 (BTD): ver7.0.8
!                 Added some more output to report begin/end reading file
! 19.02.05 (BTD): ver7.3.2
!                 modified to allow input of eps-1 or N-1
!                 by using IREN  = - column number that  Re(N-1)  is stored in
!                 or       IEPS1 = - column number that Re(eps-1) is stored in
! 19.08.31 (BTD): ver7.3.3
!                 corrected error in logic for case where input was eps-1
!                 (cases where epsilon or m or m-1 were read in were OK)
! end history
! Copyright (C) 1993,1996,1998,2004,2007,2009,2019
!               B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.
!***********************************************************************
      DATA CXI/(0._WP,1._WP)/,INIT/.TRUE./

!*** diagnostic
!      write(0,*)'dielec ckpt 0'
!      write(0,*)'  init=',init
!      write(0,*)'  ncomp=',ncomp
!      write(0,*)'  wave=',wave
!      write(0,*)'  idvout=',idvout
!      do j=1,ncomp
!         write(0,fmt='(a,i2,a,a)')'  cfleps(',j,')=',cfleps(j)
!      enddo
!***
      IF(INIT.OR.NCOMP>1)THEN
!*** diagnostic
!         write(0,*)'dielec ckpt 1'
!***
         DO J=1,NCOMP
            WRITE(CMSGNM,FMT='(A)')'about to read file ='
            CALL WRIMSG('DIELEC',CMSGNM)
            CALL WRIMSG('DIELEC',CFLEPS(J))
            OPEN(UNIT=3,FILE=CFLEPS(J),STATUS='OLD')

! Read header line:

            READ(3,9000)CDESCR                        ! line 1
            CALL WRIMSG('DIELEC',CDESCR)
!            WRITE(IDVOUT,9000)CDESCR

! Read line specifying columns for wavelength,Re(n),Im(n),Re(eps),Im(eps)
!            write(0,*)'dielec ckpt 2, j=',j

            READ(3,*)IWV,IREN,IIMN,IEPS1,IEPS2        ! line 2

!*** diagnostic
!            write(0,*)'dielec ckpt 3, iwv,iren,iimn,ieps1,ieps2=', &
!                      iwv,iren,iimn,ieps1,ieps2
!***
            ICWV=ABS(IWV)
            ICOL=ICWV

            ICREN=ABS(IREN)
            IF(ICREN>ICOL)ICOL=ICREN

            ICIMN=ABS(IIMN)
            IF(ICIMN>ICOL)ICOL=ICIMN

            ICEPS1=ABS(IEPS1)
            IF(ICEPS1>ICOL)ICOL=ICEPS1

            ICEPS2=ABS(IEPS2)
            IF(ICEPS2>ICOL)ICOL=ICEPS2

! Skip header line:

            READ(3,*)                                 ! line 3
            DO I=1,MXWAVT
               READ(3,*,END=600)(XIN(JJ),JJ=1,ICOL)   ! line 4,5,...
               IF(IWV>0)THEN
                  WVA(I)=XIN(ICWV)            ! XIN(ICWV) = wavelength
               ELSEIF(IWV<0)THEN
                  WVA(I)=1.23984/XIN(ICWV)   ! XIN(ICWV) = E(eV)
               ENDIF
               IF(IEPS1/=0)THEN
                  IF(IEPS1>0)THEN
                     E1A(I)=XIN(ICEPS1)
                  ELSE
                     E1A(I)=1.+XIN(ICEPS1)
                  ENDIF
                  E2A(I)=XIN(ICEPS2)
               ELSEIF(IIMN/=0)THEN
                  IF(IIMN>0)THEN
                     E1A(I)=XIN(ICREN)
                  ELSE
                     E1A(I)=1.+XIN(ICREN)
                  ENDIF
                  E2A(I)=XIN(ICIMN)
               ELSE
                  WRITE(CMSGNM,FMT='(A,I3,A,I3)')'fatal error: IREN=', &
                      IREN,' and IEPS1=',IEPS1                         !
                  CALL ERRMSG('FATAL','DIELEC',CMSGNM)
               ENDIF
               NWAVT=I
            ENDDO
            IEPS1=ABS(IEPS1)   ! IEPS1 > 0 if E1A,E2A=Re(epsilon),Im(epsilon)
                               ! IEPS1 = 0 if E1A,E2A=Re(index),Im(index)

! Check whether there is unread data remaining in file

            READ(3,*,END=600)(XIN(JJ),JJ=1,ICOL)

! If this point is reached, apparently unread data remains, so
! issue warning:

            CALL ERRMSG('WARNING','DIELEC', &
                        'parameter MXWAVT not large enough to read full dielec file')
600         CLOSE(3)

!*** diagnostic
!            write(0,*)'dielec ckpt 4, j=',j,' cfleps(j)=',cfleps(j)
!            write(0,*)'  nwavt=',nwavt 
!            do i=1,nwavt
!               write(0,fmt='(a,i2,a,1pe10.3)')'  wva(',i,')=',wva(i)
!            enddo
!***

            WRITE(CMSGNM,FMT='(A,A)')' completed reading file ='
            CALL WRIMSG('DIELEC',CMSGNM)
            CALL WRIMSG('DIELEC',CFLEPS(J))

!*** Have completed reading in table for this composition.
!    Now interpolate

!*** diagnostic
!            write(0,*)'dielec ckpt 5'
!***
            CALL INTERP(WAVE,E1,E2,WVA,E1A,E2A,IDVOUT,MXWAVT,NWAVT)
            IF(IEPS1.NE.0)THEN
               CXEPS(J)=E1+CXI*E2
            ELSE
               CXEPS(J)=(E1**2-E2**2)+2._WP*CXI*E1*E2
            ENDIF

! check that user has not specified refractive index = 1
! issue fatal warning if this occurs

            IF(CXEPS(J)==1._WP)CALL ERRMSG('FATAL','DIELEC', &
               'Refractive index = 1 is not allowed: should leave unoccupied')

         ENDDO ! enddo j=1,ncomp
         INIT=.FALSE.
      ELSE

!*** Perform this only if NCOMP=1 and previously initialized:

!*** diagnostic
!         write(0,*)'dielec ckpt 6'
!         write(0,*)'  nwavt=',nwavt
!         do j=1,nwavt
!            write(0,fmt='(a,i2,a,1pe10.3)')'  wva(',j,')=',wva(j)
!         enddo
!***

         CALL INTERP(WAVE,E1,E2,WVA,E1A,E2A,IDVOUT,MXWAVT,NWAVT)
         IF(ICEPS1>0)THEN
            CXEPS(1)=E1+CXI*E2
         ELSE
            CXEPS(1)=(E1**2-E2**2)+2._WP*CXI*E1*E2
         ENDIF
      ENDIF
      RETURN
9000  FORMAT(A70)
    END SUBROUTINE DIELEC
