    SUBROUTINE ERRMSG(CSTATS,CSUBRT,CMSGNM)
      IMPLICIT NONE
! Arguments:
      CHARACTER :: CMSGNM*(*), CSTATS*(*), CSUBRT*(*)

! if MPI is not installed, comment out the INCLUDE 'mpif.h' statement
#ifdef mpi
!      INCLUDE 'mpif.h'
#endif

#ifndef mpi
      INTEGER :: MPI_COMM_WORLD
#endif

! Local variables:
      INTEGER :: IOERR

!***********************************************************************
! Given:
!       CSTATS = 'WARNING' or 'FATAL'
!       CSUBRT = name of subroutine
!       CMESGN = message
! Prints a warning message in a standardized way,
! and STOPs if CSTATS='FATAL'

! History:
! 96.11.14 (PJF) Remove "getset" and hardwire "ioerr"
! 04.05.23 (BTD) cleanup
! 10.05.08 (BTD) modify output
! 19.07.10 (BTD) following suggestion from Fangzhou Liu:
!                add CALL MPI_ABORT
!                and CALL MPI_FINALIZE
!                so that master process will stop if slave process aborts 
!                need to also add INCLUDE 'mpif.h'
! 19.12.18 (BTD) need to comment out the INCLUDE 'mpif.h' when MPI is not
!                installed
! end history
! Copyright (C) 1993,1996,2004,2010,2019 B.T. Draine and P.J. Flatau
! This code is covered by the GNU General Public License.

!***********************************************************************
      DATA IOERR/6/

      IF(CSTATS=='FATAL')THEN
         WRITE(IOERR,FMT=9000)CSUBRT
         WRITE(IOERR,FMT=9010)CMSGNM
         WRITE(IOERR,FMT=9020)
! added by Fangzhou Liu to terminate all MPI processes
! otherwise when slave processes exit due to fatal errors,
! the master process will go on:
         CALL MPI_ABORT(MPI_COMM_WORLD,IOERR)
         CALL MPI_FINALIZE(IOERR)
         STOP
      ELSEIF(CSTATS=='WARNING')THEN
         WRITE(IOERR,FMT=9030)CSUBRT
         WRITE(IOERR,FMT=9010)CMSGNM
      ENDIF
      RETURN
9000  FORMAT(/' >>>>> FATAL ERROR IN PROCEDURE: ',A)
9010  FORMAT(' >>>>> ',A)
9020  FORMAT(' >>>>> EXECUTION ABORTED ')
9030  FORMAT(/' >>>>> WARNING IN PROCEDURE: ',A)
    END SUBROUTINE ERRMSG
