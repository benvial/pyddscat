      PROGRAM CHECKER
      IMPLICIT NONE
      CHARACTER CHAR1
      LOGICAL NOPROBLEM
      INTEGER J,JPBC,JSKIP,JTRQ,NCOMP
      REAL QEXT,S1,S2,S3,S4,S5,S6,S7,TOL,V1,V2,V3,V4,V5,V6,V7

!=================================================================
! program checker v1
! purpose: 
!    read file w000r000k000.sca
!    compare selected entries to stored file w000r000k000.sca.sav
!    report discrepancies if fractional difference exceeds TOL
! intended to work for JPBC=0,1,2, or 3
!
! B.T. Draine, Princeton University, 2014.10.14
! history
! 14.10.14 (BTD) first written to use in "runexamples.btd" suite
! end history
!=================================================================
      TOL=5.E-4

      NOPROBLEM=.TRUE.
      OPEN(UNIT=3,FILE='w000r000k000.sca.sav')
      OPEN(UNIT=4,FILE='w000r000k000.sca')

      JSKIP=5
      DO J=1,JSKIP
         READ(3,*)   ! lines 1-5
         READ(4,*)   ! 
      ENDDO

      READ(3,*)JPBC  ! line 6
      READ(3,*)JTRQ  ! line 7
      READ(3,*)NCOMP ! line 8
      READ(4,*)
      READ(4,*)
      READ(4,*)

      IF(JPBC.EQ.0)THEN
         JSKIP=28+NCOMP
      ELSEIF(JPBC.EQ.1.OR.JPBC.EQ.2)THEN
         JSKIP=29+NCOMP
      ELSEIF(JPBC.EQ.3)THEN
         JSKIP=26+NCOMP
      ENDIF

! jpbc=0     : skip lines 9 - 36+NCOMP
! jpbc=1 or 2: skip lines 9 - 37+NCOMP
! jpbc=3     : skip lines 9 - 34+NCOMP

      DO J=1,JSKIP
         READ(3,*) !jpbc=0: 9 - 36+NCOMP ; JPBC=1: 9 - 37+NCOMP
         READ(4,*)
      ENDDO

      IF(JPBC.EQ.0)THEN
         READ(3,*)CHAR1,S1,S2,S3,S4,S5,S6,S7 ! line 37+NCOMP
         READ(4,*)CHAR1,V1,V2,V3,V4,V5,V6,V7
         QEXT=S1
         IF(.NOT.(ABS(V1-S1).LT.TOL*ABS(S1)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qext=',V1,' vs',S1,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
      
         IF(.NOT.(ABS(V2-S2).LT.TOL*ABS(S2)))THEN
             WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qabs=',V2,' vs',S2,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V3-S3).LT.TOL*ABS(S3)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qsca=',V3,' vs',S3,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V4-S4).LT.TOL))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: <cos>=',V4,' vs',S4,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V5-S5).LT.TOL))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: <cos^2>=',V5,' vs',S5,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V6-S6).LT.TOL*ABS(S6)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qbk=',V6,' vs',S6,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V7-S7).LT.TOL*ABS(S7)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qpha=',V7,' vs',S7,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF

         READ(3,*)CHAR1,S1,S2,S3,S4,S5,S6,S7 ! line 38+NCOMP
         READ(4,*)CHAR1,V1,V2,V3,V4,V5,V6,V7
         QEXT=0.5*(QEXT+S1)   ! pol-averaged Qext for future ref.
         IF(.NOT.(ABS(V1-S1).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qext=',V1,' vs',S1,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
      
         IF(.NOT.(ABS(V2-S2).LT.TOL*ABS(S2)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qabs=',V2,' vs',S2,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V3-S3).LT.TOL*ABS(S3)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qsca=',V3,' vs',S3,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V4-S4).LT.TOL))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: <cos>=',V4,' vs',S4,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V5-S5).LT.TOL))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: <cos^2>=',V5,' vs',S5,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V6-S6).LT.TOL*ABS(S6)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qbk=',V6,' vs',S6,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V7-S7).LT.TOL*ABS(S7)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qpha=',V7,' vs',S7,' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF

         JSKIP=3
         DO J=1,JSKIP
            READ(3,*) ! line 39+NCOMP - 41+NCOMP
            READ(4,*)
         ENDDO
!*** diagnostic
!         write(0,*)'checker ckpt 1, s1,s2,s3=',s1,s2,s3
!***
         READ(3,*)CHAR1,S1,S2,S3 ! line 42+NCOMP
!*** diagnostic
!         write(0,*)'checker ckpt 2, s1,s2,s3=',s1,s2,s3
!***
         READ(4,*)CHAR1,V1,V2,V3
         IF(.NOT.(ABS(V1-S1).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qsca*g(1)=',V1,' vs',S1,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V2-S2).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qsca*g(2)=',V2,' vs',S2,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V3-S3).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: Qsca*g(3)=',V3,' vs',S3,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         READ(3,*)CHAR1,S1,S2,S3 ! line 43+NCOMP
         READ(4,*)CHAR1,V1,V2,V3
         IF(.NOT.(ABS(V1-S1).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qsca*g(1)=',V1,' vs',S1,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V2-S2).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qsca*g(2)=',V2,' vs',S2,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V3-S3).LT.TOL*QEXT))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: Qsca*g(3)=',V3,' vs',S3,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         JSKIP=2
         DO J=1,JSKIP
            READ(3,*)
            READ(4,*)
         ENDDO

         IF(JTRQ.EQ.1)THEN
            READ(3,*)
            READ(4,*)
            READ(3,*)CHAR1,S1,S2,S3,S4,S5,S6 ! line 46+NCOMP
            READ(4,*)CHAR1,V1,V2,V3,V4,V5,V6
            IF(.NOT.(ABS(V1-S1).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=1: Qtrqab(1)=',V1,' vs',S1,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V2-S2).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=1: Qtrqab(2)=',V2,' vs',S2,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V3-S3).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=1: Qtrqab(3)=',V3,' vs',S3,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V4-S4).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=1: Qtrqsc(1)=',V4,' vs',S4,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V5-S5).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=1: Qtrqsc(2)=',V5,' vs',S5,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V6-S6).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=1: Qtrqsc(3)=',V6,' vs',S6,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            READ(3,*)CHAR1,S1,S2,S3,S4,S5,S6 ! line 47+NCOMP
            READ(4,*)CHAR1,V1,V2,V3,V4,V5,V6
            IF(.NOT.(ABS(V1-S1).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=2: Qtrqab(1)=',V1,' vs',S1,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V2-S2).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=2: Qtrqab(2)=',V2,' vs',S2,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V3-S3).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=2: Qtrqab(3)=',V3,' vs',S3,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V4-S4).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=2: Qtrqsc(1)=',V4,' vs',S4,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V5-S5).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=2: Qtrqsc(2)=',V5,' vs',S5,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V6-S6).LT.TOL*QEXT))THEN
               WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &              'problem for JO=2: Qtrqsc(3)=',V6,' vs',S6,
     &              ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
         ENDIF

         JSKIP=3
         DO J=1,JSKIP ! skip 44+NCOMP+4*JTRQ - 46+NCOMP+4*JRTQ
            READ(3,*)
            READ(4,*)
         ENDDO

! read scattering for theta=0

         READ(3,*)S1,S2,S3,S4,S5,S6,S7 ! line 47+NCOMP+4*JTRQ
         READ(4,*)V1,V2,V3,V4,V5,V6,V7
         IF(.NOT.(ABS(V3-S3).LT.TOL))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for theta,phi=',S1,S2,' Pol=',V3,' vs',S3,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V4-S4).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for theta,phi=',S1,S2,' S_11=',V4,' vs',S4,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V5-S5).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for theta,phi=',S1,S2,' S_12=',V5,' vs',S5,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V6-S6).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for theta,phi=',S1,S2,' S_21=',V6,' vs',S6,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V7-S7).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for theta,phi=',S1,S2,' S_22=',V7,' vs',S7,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF

      ELSEIF(JPBC.EQ.1.OR.JPBC.EQ.2)THEN

         READ(3,*)S1,S2,S3,S4,S5,S6,S7 ! line 38+NCOMP
         READ(4,*)V1,V2,V3,V4,V5,V6,V7
         IF(.NOT.(ABS(V3-S3).LT.TOL))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' Pol=',V3,' vs',S3,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V4-S4).LT.TOL*ABS(S4)))THEN
            WRITE(0,*)'problem with S_11 for alpha,zeta=',S1,S2
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_11=',V4,' vs',S4,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V5-S5).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_12=',V5,' vs',S5,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V6-S6).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_21=',V6,' vs',S6,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V7-S7).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_22=',V7,' vs',S7,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF

         JSKIP=10
         IF(JSKIP.GT.0)THEN
            DO J=1,JSKIP
               READ(3,*)
               READ(4,*)
            ENDDO
         ENDIF

         READ(3,*)S1,S2,S3,S4,S5,S6,S7 ! line 39+NCOMP+JSKIP
         READ(4,*)V1,V2,V3,V4,V5,V6,V7
         IF(.NOT.(ABS(V3-S3).LT.TOL))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' Pol=',V3,' vs',S3,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V4-S4).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_11=',V4,' vs',S4,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V5-S5).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_12=',V5,' vs',S5,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V6-S6).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_21=',V6,' vs',S6,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         IF(.NOT.(ABS(V7-S7).LT.TOL*ABS(S4)))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &         'problem for alpha,zeta=',S1,S2,' S_22=',V7,' vs',S7,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF

      ELSEIF(JPBC.EQ.3)THEN
         READ(3,*)CHAR1,S1 ! line 35+NCOMP
         READ(4,*)CHAR1,V1
         IF(.NOT.(ABS(V1-S1).LT.4.*TOL*ABS(S1)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=1: abs.coeff=',V1,' vs',S1,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         READ(3,*)CHAR1,S1 ! line 36+NCOMP
         READ(4,*)CHAR1,V1
         IF(.NOT.(ABS(V1-S1).LT.4.*TOL*ABS(S1)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for JO=2: abs.coeff=',V1,' vs',S1,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         READ(3,*)CHAR1,S1 ! line 37+NCOMP
         READ(4,*)CHAR1,V1
         IF(.NOT.(ABS(V1-S1).LT.4.*TOL*ABS(S1)))THEN
            WRITE(0,FMT='(A,1PE11.3,A,1PE11.3,A)')
     &         'problem for mean abs.coeff=',V1,' vs',S1,
     &         ' (previous calc.)'
            NOPROBLEM=.FALSE.
         ENDIF
         JSKIP=2
         DO J=1,JSKIP
            READ(3,*)
            READ(4,*)
         ENDDO
         DO J=1,2
            READ(3,*)S1,S2,S3,S4,S5,S6,S7 ! line 40+NCOMP
            READ(4,*)V1,V2,V3,V4,V5,V6,V7
            IF(.NOT.(ABS(V3-S3).LT.5.*TOL))THEN
            WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &            'problem for theta,phi=',S1,S2,' Pol=',V3,' vs',S3,
     &            ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V4-S4).LT.6.*TOL*ABS(S4)))THEN
               WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &            'problem for theta,phi=',S1,S2,' S_11=',V4,' vs',S4,
     &            ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V5-S5).LT.5.*TOL*ABS(S4)))THEN
               WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &            'problem for theta,phi=',S1,S2,' S_12=',V5,' vs',S5,
     &            ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V6-S6).LT.5.*TOL*ABS(S4)))THEN
               WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &            'problem for theta,phi=',S1,S2,' S_21=',V6,' vs',S6,
     &            ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
            IF(.NOT.(ABS(V7-S7).LT.6.*TOL*ABS(S4)))THEN
               WRITE(0,FMT='(A,2F7.2,A,1PE11.3,A,1PE11.3,A)')
     &            'problem for theta,phi=',S1,S2,' S_22=',V7,' vs',S7,
     &            ' (previous calc.)'
               NOPROBLEM=.FALSE.
            ENDIF
         ENDDO
      ENDIF

      IF(NOPROBLEM)THEN
         WRITE(0,*)'results in w000r000k000.sca look OK'
      ELSE
         WRITE(0,*)'       --- !!! WARNING !!! ---'
         WRITE(0,*)'results in w000r000k000.sca appear to be in error'
      ENDIF
      STOP
      END
