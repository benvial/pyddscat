    SUBROUTINE PRECOND
      RETURN
    END SUBROUTINE PRECOND

    SUBROUTINE ZBCG2(PRINT_RESID,L,N,X,NONZERO_X,RHS,MATVEC,PRECOND,TOLER, &
                     MXMATVEC,WORK,INFO)
      USE DDPRECISION,ONLY : WP
      USE DDCOMMON_9,ONLY: ITERMX,ITERN

! subroutine zbcg2 (print_resid,l,n,x,nonzero_x,rhs,matvec,precond,toler, &
!                   mxmatvec,work,info)

! Improved "vanilla" BiCGStab(2) iterative method

! Copyright (c) 2001 by M.A.Botchev (http://www.math.utwente.nl/~botchev/),
!                       University of Twente
! Permission to copy all or part of this work is granted,
! provided that the copies are not made or distributed
! for resale, and that the copyright notice and this
! notice are retained.

! This is the "vanilla" version of BiCGstab(\ell) as described
! in PhD thesis of D.R.Fokkema, Chapter 3 (also available as
! Preprint 976, Dept. of Mathematics, Utrecht University, URL
! http://www.math.uu.nl/publications/).  It includes two enhancements 
! to BiCGstab(\ell) proposed by G.Sleijpen and H.van der Vorst in
! 1) G.Sleijpen and H.van der Vorst "Maintaining convergence 
!    properties of BiCGstab methods in finite precision arithmetic",
!    Numerical Algorithms, 10, 1995, pp.203-223
! 2) G.Sleijpen and H.van der Vorst "Reliable updated residuals in
!    hybrid BiCG methods", Computing, 56, 1996, pp.141-163

! {{ This code based on original work of D.R.Fokkema:

! subroutine zbistbl v1.1 1998    
! Copyright (c) 1995-1998 by D.R. Fokkema.
! Permission to copy all or part of this work is granted, 
! provided that the copies are not made or distributed 
! for resale, and that the copyright notice and this 
! notice are retained.

! }}

! Your bug reports, comments, etc. are welcome: 
! m.a.botchev@math.utwente.nl

! ------------------------------
! Description of the parameters:
! ------------------------------

! print_resid (input) LOGICAL. If print_resid=.true. the number of 
!            matrix-vector multiplications done so far and residual norm will
!            be printed to the standard output each iteration

! l          (input) INTEGER the dimension \ell of BiCGstab(\ell)
!            in this simple version it is required that l <= 2
!            l=2 is often useful for systems with nonsymmetric matrices

! n          (input) INTEGER size of the linear system to solve 

! x          (input/output) COMPLEX*16 array dimension n
!            initial guess on input, solution on output

! rhs        (input) COMPLEX*16 array dimension n
!            the right-hand side (r.h.s.) vector

! matvec     (input) EXTERNAL name of matrix vector subroutine
!            to deliver y:=A*x by CALL matvec(n,x,y)

! nonzero_x  (input) LOGICAL tells
!            BiCGstab(\ell) if the initial guess x is zero or not. 
!            If nonzero_x is .FALSE., initial residual equals r.h.s. vector
!            and one MATVEC call is avoided

! toler      (input/output) DOUBLE PRECISION tolerance: the iterations are 
!            stopped as soon as || residual ||/|| initial residual|| <= toler,
!            the norm is Euclidean.  On output, if info>=0, the value of 
!            toler is set to the actually achieved residual reduction

! mxmatvec   (input/output) INTEGER.  On input: maximum number of matrix 
!            vector multiplications allowed to be done.  On output: 
!            if info>=0, mxmatvec is set to the actual number of matrix 
!            vector multiplications done

! work       (workspace) COMPLEX*16 array of dimension (n,2*l+5)

! info       (output) INTEGER.  info = 0 in case of succesful computations
!            and 
!            info = -m (<0) - means paramater number m has an illegal value
!            info = 1 - means no convergence achieved (stopping criterion
!            is not fulfilled)
!            info = 2 - means breakdown of the algorithm (taking a larger
!            value of parameter l usually helps)

! WARNING: If the iterations are ended normally (info=0 or info=1),
! the true residual norm is computed and returned as an output value 
! of the parameter toler.  The true residual norm can be slightly larger
! than the projected residual norm used by the algorithm to stop the
! iterations.  It may thus happen that on output info=0 but the value
! of toler is (slightly) larger than tolerance prescribed on input.
! ----------------------------------------------------------
! history:
! 08.03.14 (BTD) added declaration of dummy array IPAR(*)
!                changed
!                CALL MATVEC(X,WORK(1:N,R),N) -> 
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
! 08.05.12 (BTD) Following suggestion by Art Lazanoff, NASA Ames:
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
!                -> CALL MATVEC(X,WORK(1,R),IPAR) ->
!                RNRM0=DNORM2_BCG(N,WORK(1:N,R))
!                -> RNRM0=DNORM2_BCG(N,WORK(1,R))
!                RHO1=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,R+K-1))
!                -> RHO1=ZDOT_BCG(N,WORK(1,RR),WORK(1,R+K-1))
!                CALL MATVEC(WORK(1:N,U+K-1),WORK(1:N,U+K),IPAR)
!                -> CALL MATVEC(WORK(1,U+K-1),WORK(1,U+K),IPAR)
!                SIGMA=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,U+K))
!                -> SIGMA=ZDOT_BCG(N,WORK(1,RR),WORK(1,U+K))
!                CALL MATVEC(WORK(1:N,R+K-1),WORK(1:N,R+K),IPAR)
!                -> CALL MATVEC(WORK(1,R+K-1),WORK(1,R+K),IPAR)
!                MATRIX_Z(I,J)=
!                            CONJG(ZDOT_BCG(N,WORK(1:N,R+J-1),WORK(1:N,R+I-1)))
!                -> MATRIX_Z(I,J)=
!                            CONJG(ZDOT_BCG(N,WORK(1,R+J-1),WORK(1,R+I-1)))
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
!                -> CALL MATVEC(X,WORK(1,R),IPAR)
!                CALL MATVEC(X,WORK(1:N,R),IPAR)
!                -> CALL MATVEC(X,WORK(1,R),IPAR)
!                RNRM=DNORM2_BCG(N,WORK(1:N,R))
!                -> RNRM=DNORM2_BCG(N,WORK(1,R))
! 08.07.16 (BTD) Added DDCOMMON_9 to communicate ITERMX and ITERN
! 10.04.30 (BTD) Modified manner in which number of iterations is
!                limited (previous treatment had return without
!                error message when ITERN=ITERMX/4)
! end history
      IMPLICIT NONE

! Parameters:

      LOGICAL,INTENT(IN) :: PRINT_RESID,NONZERO_X
      INTEGER,INTENT(IN) :: L,N
      INTEGER,INTENT(INOUT) :: MXMATVEC
      INTEGER,INTENT(OUT) :: INFO
      COMPLEX(WP),INTENT(INOUT) :: X(N)
      COMPLEX(WP),INTENT(IN) :: RHS(N)
      REAL(WP),INTENT(INOUT) :: TOLER
      COMPLEX(WP),INTENT(OUT) :: WORK(N,3+2*(L+1))
      EXTERNAL MATVEC, PRECOND

! Local variables:

      CHARACTER :: CMSGNM*70
      COMPLEX(WP) :: MATRIX_Z(L+1,L+1),Y0(L+1),YL(L+1),ZY0(L+1),ZYL(L+1)
      LOGICAL :: RCMP,XPDT
      INTEGER :: I,J,K,NMATVEC
      integer :: jj
      COMPLEX(WP) :: ALPHA,BETA,OMEGA,RHO0,RHO1,SIGMA
      COMPLEX(WP) :: VARRHO,HATGAMMA
      REAL(WP) :: RNRM0,RNRM
      REAL(WP) :: MXNRMX,MXNRMR
      COMPLEX(WP) :: KAPPA0,KAPPAL

! Define dummy array IPAR

      INTEGER IPAR(13)

! Aliases for the parts of the work array:

      INTEGER :: RR,R,U,XP,BP

! Constants:

      REAL(WP),PARAMETER :: DELTA=1D-2
      COMPLEX(WP),PARAMETER :: ZZERO=(0D0,0D0),ZONE=(1D0,0D0)

! Functions:

      REAL(WP) :: DNORM2_BCG
      COMPLEX(WP) :: ZDOT_BCG

!---------------------------------------------------------------------------
!*** diagnostic
!      write(0,*)'zbcg2wp ckpt 1'
!      write(0,*)'   j  ------ x(j) --------'
!      do k=1,n
!         write(0,fmt='(i6,1p2e11.3)')k,x(k)
!      enddo
!***
      INFO=0

      IF(L<1 .OR. L>2)INFO=-2
      IF(N<1)INFO=-3
      IF(TOLER<=0E0_WP)INFO=-9
      IF(MXMATVEC<0)INFO=-10

      RR=1
      R=RR+1
      U=R+(L+1)
      XP=U+(L+1)
      BP=XP+1

      IF(INFO/=0)RETURN

!*** diagnostic
!      do j=1,n
!         if(.not.(abs(rhs(j)).lt.1.e32))then
!            write(0,*)'zbcg2wp ckpt 2 problem: rhs(j)=',rhs(j),' for j=',j
!         endif
!      enddo
!      do j=1,n
!         if(.not.(abs(x(j)).lt.1.e32))then
!            write(0,*)'zbcg2wp ckpt 3 problem: x(j)=',x(j),' for j=',j
!         endif
!      enddo
!***

! Initialize first residual

      IF(NONZERO_X)THEN
!*** diagnostic
!         write(0,*)'zbcg2wp ckpt 4'
!****
!        CALL MATVEC(X,WORK(1:N,R),N)
!         CALL MATVEC(X,WORK(1:N,R),IPAR)
         CALL MATVEC(X,WORK(1,R),IPAR)
!*** diagnostic
!         write(0,*)'zbcg2wp ckpt 5'
!****

         WORK(1:N,R)=RHS-WORK(1:N,R)
         NMATVEC=1
      ELSE
         WORK(1:N,R)=RHS
         NMATVEC=0
      ENDIF

!*** diagnostic
!      do j=1,N
!         if(.not.(abs(work(j,R)).lt.1.e32))then
!            write(0,*)'zbcg2wp ckpt 6 problem: work(j,R)=', &
!                      work(j,R),' for j=',j
!         endif
!      enddo
!***

!call precond (n,work(1:n,r))

! Initialize iteration loop

      WORK(1:N,RR)=WORK(1:N,R)
      WORK(1:N,BP)=WORK(1:N,R)
      WORK(1:N,XP)=X
      X=ZZERO

!      RNRM0=DNORM2_BCG(N,WORK(1:N,R))
      RNRM0=DNORM2_BCG(N,WORK(1,R))
      RNRM=RNRM0

      MXNRMX=RNRM0
      MXNRMR=RNRM0
      RCMP=.FALSE.
      XPDT=.FALSE.

      ALPHA=ZZERO
      OMEGA=ZONE
      SIGMA=ZONE
      RHO0=ZONE

!BTD 080716:
      ITERN=0
!-----------

!*** diagnostic
!      write(0,*)'zbcg2wp ckpt 7: itermx=',itermx
!***
!----------

! Iterate

! BTD 100430 change
!      DO WHILE (RNRM>TOLER*RNRM0 .AND. NMATVEC<MXMATVEC)
      DO WHILE (RNRM>TOLER*RNRM0 .AND. ITERN<=ITERMX)

!BTD 080716:
         ITERN=ITERN+1
!*** diagnostic
!         write(0,*)'zbcg2wp ckpt 8: itern=',itern,' itermx=',itermx
!***
         IF(ITERN>ITERMX)CALL ERRMSG('FATAL','ZBCG2WP',' ITERN>ITERMX')
!-----------
         
! =====================
! The BiCG part ---
! =====================

         RHO0=-OMEGA*RHO0
         DO K=1,L
!            RHO1=ZDOT_BCG(N,WORK(1:N,RR),WORK(1:N,R+K-1))
            RHO1=ZDOT_BCG(N,WORK(1,RR),WORK(1,R+K-1))
!*** diagnostic
!            write(0,fmt='(a,i3,a,1p2e12.4)')'ckpt 8.1, K=',k,' rho1=',rho1
!***
            IF(RHO0==ZZERO)THEN
               INFO=2
               TOLER=RNRM/RNRM0
               MXMATVEC=NMATVEC
               RETURN
            ENDIF
            BETA=ALPHA*(RHO1/RHO0)
            RHO0=RHO1
            DO J=0,K-1
               WORK(1:N,U+J)=WORK(1:N,R+J)-BETA*WORK(1:N,U+J)
            ENDDO
!*** diagnostic
!            do j=1,n
!               do jj=1,4
!                  if(.not.(abs(work(j,jj)).eq.0.))     &
!                    write(0,fmt='(a,i5,i2,a,1p2e12.4)') &
!                    'zbcg2wp ckpt 9, j,jj=',j,jj,      &
!                    ' work(j,jj)=',work(j,jj)
!               enddo
!           enddo
!***
           CALL MATVEC(WORK(1,U+K-1),WORK(1,U+K),IPAR)
!*** diagnostic
!           write(0,*)'zbcq2wp ckpt 10'
!            do j=1,n
!               do jj=1,4
!                  if(.not.(abs(work(j,jj)).eq.0.))     &
!                    write(0,fmt='(a,i5,i2,a,1p2e12.4)') &
!                    'zbcg2wp ckpt 10, j,jj=',j,jj,     &
!                    ' work(j,jj)=',work(j,jj)
!               enddo
!           enddo
!***
!***
!          call precond (n, work(1:n,u+k))
           NMATVEC=NMATVEC+1
!*** diagnostic
!           write(0,fmt='(a,i5,a,i5,a,i5)')'ckpt 10.4, n=',n,' rr=',rr, &
!              ' u+k=',(u+k)
!           do j=1,n
!              write(0,fmt='(a,i5,a,1p2e11.3,a,1p2e11.3,a,1p2e11.3)') &
!                 'j=',j,                                             &
!                 ' work(j,rr)=',work(j,rr),                          &
!                 ' work(j,u+k-1)=',work(j,u+k-1)                   
!                 ' work(j,u+k)=',work(j,u+k)
!           enddo
!***
           SIGMA=ZDOT_BCG(N,WORK(1,RR),WORK(1,U+K))
!*** diagnostic
!            write(0,fmt='(a,i3,a,1p2e12.4)')'ckpt 10.5, K=',k,' sigma=',sigma
!***

           IF(SIGMA==ZZERO)THEN
              INFO=2
              TOLER=RNRM/RNRM0
              MXMATVEC=NMATVEC
              RETURN
           ENDIF
           ALPHA=RHO1/SIGMA
!*** diagnostic
!            write(0,fmt='(a,i3,a,1p2e12.4)')'ckpt 10.6, K=',k,' alpha=',alpha
!***
           X(1:N)=ALPHA*WORK(1:N,U)+X(1:N)
!*** diagnostic
!           do j=1,n
!              if(.not.(abs(x(j)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 11 problem: x(j)=',x(j),' for j=',j
!              endif
!           enddo
!***
!*** diagnostic
!           write(0,fmt='(a,1p2e12.4)')'ckpt 11.1 alpha=',alpha
!           do j=1,n
!              do jj=1,4
!                 if(.not.(abs(work(j,jj)).eq.0.))       &
!                    write(0,fmt='(a,i5,i2,a,1p2e12.4)')  &
!                       'zbcg2wp ckpt 11.2, j,jj=',j,jj, &
!                       ' work(j,jj)=',work(j,jj)
!              enddo
!           enddo
!***

           DO J=0,K-1
              WORK(1:N,R+J)=-ALPHA*WORK(1:N,U+J+1)+WORK(1:N,R+J)
           ENDDO
!*** diagnostic
!           do j=1,n
!              do jj=1,4
!                 if(.not.(abs(work(j,jj)).eq.0.))       &
!                    write(0,fmt='(a,i5,i2,a,1p2e12.4)')  &
!                       'zbcg2wp ckpt 11.3, j,jj=',j,jj, &
!                       ' work(j,jj)=',work(j,jj)
!              enddo
!           enddo
!***
           CALL MATVEC(WORK(1,R+K-1),WORK(1,R+K),IPAR)
!*** diagnostic
!           do j=1,n
!              do jj=1,4
!                 if(.not.(abs(work(j,jj)).eq.0.))then
!                    write(0,fmt='(a,i5,i2,a,1p2e9.1)') &
!                    'zbcg2wp ckpt 11.6, j,jj=',j,jj,   &
!                    ' work(j,jj)=',work(j,jj)
!                 endif
!              enddo
!           enddo
!!!           stop
!***
!      call precond (n, work(1:n,r+k))
           NMATVEC=NMATVEC+1
           RNRM=DNORM2_BCG(N,WORK(1,R))
           MXNRMX=MAX(MXNRMX,RNRM)
           MXNRMR=MAX(MXNRMR,RNRM)
        ENDDO

! ==================================
! The convex polynomial part ---
! ==================================

!  --- Z = R'R

!*** diagnostic
!        write(0,*)'zbcg2wp ckpt 11.5, R=',R,' L=',L
!***
        DO I=1,L+1
           DO J=1,I

!*** diagnostic
!              if(.not.(abs(work(1,r+j-1)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 12 problem for i,j=',i,j,     &
!                           ' : work(1,r+j-1)=',                        &
!                           work(1,r+j-1),' for i,j=',i,j,' r+j-1=',r+j-1
!              endif
!              if(.not.(abs(work(1,r+i-1)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 13 problem for i,j=',i,j,     &
!                           ' : work(1,r+i-1)=',                        &
!                           work(1,r+i-1),' for i,j=',i,j,' r+i-1=',r+i-1
!              endif
!***
!*** diagnostic
!              write(0,*)'zbcg2wp ckpt 13.5, i,j=',i,j
!              do jj=1,n
!                 if(.not.(abs(work(jj,R+J-1)).lt.1.e10))then
!                    write(0,*)'work(j,k)=',work(jj,r+j-1),' for j,k=',jj,r+j-1
!                 endif
!                 if(.not.(abs(work(jj,R+I-1)).lt.1.e10))then
!                    write(0,*)'work(j,k)=',work(jj,r+i-1),' for j,k=',jj,r+i-1
!                 endif
!              enddo
!***

              MATRIX_Z(I,J)=CONJG(ZDOT_BCG(N,WORK(1,R+J-1),WORK(1,R+I-1)))

!*** diagnostic
!              if(.not.(abs(matrix_z(i,j)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 14 problem: matrix_z(i,j)=', &
!                           matrix_z(i,j),' for i,j=',i,j,' r=',r
!                 write(0,*)'   work(1,r+j-1)=',  &
!                           work(1,r+j-1),' for r+j-1=',r+j-1
!                 write(0,*)'   work(1,r+i-1)=',  &
!                           work(1,r+i-1),' for r+i-1=',r+i-1 
!              endif
!***
           ENDDO
        ENDDO

!  lower triangular part of Z is computed; compute the rest knowing that Z^H=Z
        DO J=2,L+1
           MATRIX_Z(1:J-1,J)=CONJG(MATRIX_Z(J,1:J-1))
        ENDDO

!  small vectors y0 and yl

        Y0(1)=-ZONE
        Y0(2)=(MATRIX_Z(2,1)/MATRIX_Z(2,2)) ! works only for l=2
        Y0(L+1)=ZZERO

!*** diagnostic
!        if(.not.(abs(y0(2)).lt.1.e32))then
!           write(0,*)'zbcg2wp ckpt 15 problem: y0(2)=',y0(2)
!        endif
!        if(.not.(abs(y0(l+1)).lt.1.e32))then
!           write(0,*)'zbcg2wp ckpt 16 problem: y0(l+1)=',y0(l+1), &
!                     ' for l+1=',(l+1)
!        endif
!***

        YL(1)=ZZERO
        YL(2)=(MATRIX_Z(2,3)/MATRIX_Z(2,2)) ! works only for l=2
        YL(L+1)=-ZONE

!*** diagnostic
!        do j=1,l+1
!           if(.not.(abs(yl(j)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 17 problem: yl(j)=',yl(j),' for j=',j
!           endif
!        enddo
!        do j=1,l+1
!           do jj=1,l+1
!              if(.not.(abs(matrix_z(j,jj)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 18 problem: matrix_z(j,jj)=', &
!                           matrix_z(j,jj),' for j,jj=',j,jj
!              endif
!           enddo
!        enddo
!***

!  --- Convex combination

! compute Z*y0 and Z*yl
        ZY0=ZZERO
        ZYL=ZZERO
        DO J=1,L+1
           ZY0=ZY0+MATRIX_Z(:,J)*Y0(J)
           ZYL=ZYL+MATRIX_Z(:,J)*YL(J)
        ENDDO

!*** diagnostic
!        do j=1,l+1
!           if(.not.(abs(zy0(j)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 19 problem: zy0(j)=',zy0(j),' for j=',j
!           endif
!           if(.not.(abs(zyl(j)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 20 problem: zyl(j)=',zyl(j),' for j=',j
!           endif
!        enddo
!***

        KAPPA0=SQRT(ABS(ZDOT_BCG(L+1,Y0,ZY0)))
        KAPPAL=SQRT(ABS(ZDOT_BCG(L+1,YL,ZYL)))

!*** diagnostic
!        if(.not.(abs(kappa0).lt.1.e32))then
!           write(0,*)'zbcg2wp ckpt 21 problem: kappa0=',kappa0
!        endif
!        if(.not.(abs(kappal).lt.1.e32))then
!           write(0,*)'zbcg2wp ckpt 22 problem: kappal=',kappal
!        endif
!***

        VARRHO=ZDOT_BCG(L+1,YL,ZY0)/(KAPPA0*KAPPAL)

        HATGAMMA=VARRHO/ABS(VARRHO)*MAX(ABS(VARRHO),7E-1_WP)

!*** diagnostic
!        if(.not.(abs(varrho).lt.1.e32))then
!           write(0,*)'zbcg2wp ckpt 23 problem: varrho=',varrho
!        endif
!        if(.not.(abs(hatgamma).lt.1.e32))then
!           write(0,*)'zbcg2wp ckpt 24 problem: hatgamma=',hatgamma
!        endif
!***

        Y0=Y0-(HATGAMMA*KAPPA0/KAPPAL)*YL

!*** diagnostic
!        do j=1,l+1
!           if(.not.(abs(y0(j)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 25 problem: y0(j)=',y0(j),' for j=',j
!           endif
!        enddo
!*** diagnostic

!  --- Update

        OMEGA=Y0(L+1)

        DO J=1,L

!*** diagnostic
!           do jj=1,n
!              if(.not.(abs(work(jj,u)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 26 problem: work(jj,u)=', &
!                           work(jj,u),' for jj=',jj
!              endif
!           enddo
!           if(.not.(abs(y0(j+1)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 27 problem: y0(j+1)=',y0(j+1), &
!                        ' for j+1=',j+1
!           endif
!           do jj=1,n
!              if(.not.(abs(work(jj,u+j)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 28 problem: work(jj,u+j)=', &
!                           work(jj,u+j),' for jj=',jj,' u+j=',(u+j)
!              endif
!           enddo
!***

           WORK(1:N,U)=WORK(1:N,U)- Y0(J+1)*WORK(1:N,U+J)

!*** diagnostic
!           do jj=1,n
!              if(.not.(abs(work(jj,u)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 29 problem: work(jj,u)=', &
!                           work(jj,u),' for jj=',jj
!              endif
!           enddo
!***

!*** diagnostic
!           if(.not.(abs(y0(j+1)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 30 problem: y0(j+1)=',y0(j+1), &
!                        ' for j+1=',(j+1)
!           endif
!           do jj=1,n
!              if(.not.(abs(work(jj,r+j-1)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 31 problem: work(jj,r+j-1)=', &
!                           work(j,r+j-1),' for jj=',jj
!              endif
!           enddo
!***

           X(1:N)=X(1:N)+Y0(J+1)*WORK(1:N,R+J-1)

!*** diagnostic
!           do jj=1,n
!              if(.not.(abs(x(jj)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 32 problem: x(jj)=',x(jj),' for j=',j
!              endif
!           enddo
!***

!*** diagnostic
!           do jj=1,n
!              if(.not.(abs(work(jj,r)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 33 problem: work(jj,r)=', &
!                           work(j,r),' for j=',j
!              endif
!           enddo
!***
           WORK(1:N,R)=WORK(1:N,R)- Y0(J+1)*WORK(1:N,R+J)

!*** diagnostic
!           do jj=1,n
!              if(.not.(abs(work(j,u)).lt.1.e32))then
!                 write(0,*)'zbcg2wp ckpt 34 problem: work(j,u)=', &
!                           work(j,u),' for j=',j
!              endif
!           enddo
!***
        ENDDO

!*** diagnostic
!        do j=1,n
!           if(.not.(abs(x(j)).lt.1.e32))then
!              write(0,*)'zbcg2wp ckpt 35 problem: x(j)=',x(j),' for j=',j
!           endif
!        enddo
!***
        
! y0 has changed; compute Z*y0 once more

        ZY0=ZZERO
        DO J=1,L+1
           ZY0=ZY0+MATRIX_Z(:,J)*Y0(J)
        ENDDO

        RNRM=SQRT(ABS(ZDOT_BCG(L+1,Y0,ZY0)))

! ================================
! The reliable update part ---
! ================================

        MXNRMX=MAX(MXNRMX,RNRM)
        MXNRMR=MAX(MXNRMR,RNRM)
        XPDT=(RNRM<DELTA*RNRM0 .AND. RNRM0<MXNRMX)
        RCMP=((RNRM<DELTA*MXNRMR .AND. RNRM0<MXNRMR).OR. XPDT)
        IF(RCMP)THEN
           CALL MATVEC(X,WORK(1,R),IPAR)
!   call precond (n, work(1:n,r))
           NMATVEC=NMATVEC+1
           WORK(1:N,R)=WORK(1:N,BP)- WORK(1:N,R)
           MXNRMR=RNRM
           IF(XPDT)THEN

              WORK(1:N,XP)=X(1:N)+WORK(1:N,XP)
              X=ZZERO
              WORK(1:N,BP)=WORK(1:N,R)

              MXNRMX=RNRM
           ENDIF
        ENDIF

!        IF(print_resid)PRINT *,nmatvec,' ',rnrm
        
        IF(PRINT_RESID)THEN
           WRITE(CMSGNM,FMT='(A,I8,A,1P,E10.3)') &
                 'IT=',ITERN,' f.err=',RNRM/RNRM0
! BTD 100430 why user NMATVEC?
!                 'IT=',NMATVEC/4,' f.err=',RNRM/RNRM0
           CALL WRIMSG('ZBCG2 ',CMSGNM)
         ENDIF

      ENDDO

!*** diagnostic
!      write(0,*)'zbcg2wp ckpt 36, itern=',itern
!***
! =========================
! End of iterations ---
! =========================

      X(1:N)=X(1:N)+WORK(1:N,XP)

      IF(RNRM>TOLER*RNRM0)INFO=1

! compute the true residual:

! ----------- One matvec can be saved by commenting out this: ---------
!      CALL MATVEC(X,WORK(1,R),N)
!----------------------------------------------------------------------

      CALL MATVEC(X,WORK(1,R),IPAR)
      WORK(1:N,R)=RHS(1:N)- WORK(1:N,R)
!      call precond (n,work(1:n,r))
      RNRM=DNORM2_BCG(N,WORK(1,R))
      NMATVEC=NMATVEC+1

      TOLER=RNRM/RNRM0
      MXMATVEC=NMATVEC

    END SUBROUTINE ZBCG2


    FUNCTION ZDOT_BCG(N,ZX,ZY)
      USE DDPRECISION,ONLY : WP

! complex inner product function

      IMPLICIT NONE
      INTEGER,INTENT (IN) :: N
      COMPLEX(WP),INTENT (IN) :: ZX(N),ZY(N)
      COMPLEX(WP) ZDOT_BCG

! 14.10.17 (BTD) this function appears to be sometimes mis-compiled by
!                g95 compiler on ourania4, giving bad results
!                experiment with replacing original f90 line with
!                simple do loop
!      ZDOT_BCG=SUM(CONJG(ZX)*ZY)

      integer :: j

      zdot_bcg=(0._wp,0._wp)

      do j=1,n
         zdot_bcg=zdot_bcg+conjg(zx(j))*zy(j)
      enddo

    END FUNCTION ZDOT_BCG


    FUNCTION DNORM2_BCG(N,ZX)
      USE DDPRECISION,ONLY : WP

! L2 norm function

      IMPLICIT NONE
      INTEGER,INTENT (IN) :: N
      COMPLEX(WP),INTENT (IN) :: ZX(N)
      COMPLEX(WP),EXTERNAL :: ZDOT_BCG
      REAL(WP) :: DNORM2_BCG
      DNORM2_BCG=SQRT(ABS(ZDOT_BCG(N,ZX,ZX)))

    END FUNCTION DNORM2_BCG
