PROGRAM BenzenoIF

IMPLICIT NONE
INTEGER::r,i,a,K
DOUBLE PRECISION::l, w, pi, rad ,x ,y ,z 
PRINT*, 'Insira a quantidade de Benzenos que voce deseja.'
READ*,k
k=k+1
l=1.42D0
z=0
a=4.0D0*k-2.0D0
OPEN (UNIT=1,FILE='BenzenoIF.xyz')
WRITE(1,*)a
WRITE(1,*)
  DO i=1,k
   DO r=1,4
	w=r*60.0D0+30.0D0
	pi=dacos(-1.0D0)
	rad= pi*w*(2.0D0/360.0D0)
	x=(i-1)*l*sqrt(3.0D0)+l*dcos(rad)
	y=l*dsin(rad)
     IF (i==k) THEN
      IF (r==1)THEN
	CYCLE
      IF (r==4)	THEN
        CYCLE
	WRITE(1,*) 'c',x,y,z
       END IF
      END IF
     END IF
	WRITE(1,*) 'c',x,y,z
    END DO
   END DO
 CLOSE(UNIT=1)
END PROGRAM
