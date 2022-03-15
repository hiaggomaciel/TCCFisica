      PROGRAM circuito

      IMPLICIT none
      DOUBLE PRECISION:: A(3,3),W(3)
      DOUBLE PRECISION,ALLOCATABLE:: WORK(:)      
      INTEGER :: INFO,LWORK
      
      ! Determinando a matriz A
      A(1,1)=4.0D0
      A(1,2)=2.0D0
      A(1,3)=2.0D0
      A(2,2)=4.0D0
      A(2,3)=2.0D0
      A(3,3)=4.0D0
      
      ALLOCATE(WORK(1))
      LWORK=-1
      CALL DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
      LWORK=WORK(1)
      DEALLOCATE(WORK)
      ALLOCATE(WORK(LWORK))
      CALL DSYEV('V','U',3,A,3,W,WORK,LWORK,INFO)
      
      PRINT*,'lambda1=',W(1)
      PRINT*,'x1=(',A(1,1),',',A(2,1),',',A(3,1),')'
      PRINT*,
      PRINT*,'lambda2=',W(2)
      PRINT*,'x2=(',A(1,2),',',A(2,2),',',A(3,2),')'
      PRINT*,
      PRINT*,'lambda3=',W(3)
      PRINT*,'x3=(',A(1,3),',',A(2,3),',',A(3,3),')'
           

      END PROGRAM circuito
