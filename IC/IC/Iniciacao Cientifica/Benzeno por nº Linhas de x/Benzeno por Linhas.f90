PROGRAM BenzenoporLinhas
IMPLICIT NONE
INTEGER::i,k,r,a,f,n
DOUBLE PRECISION::w,pi,rad,x,y,z,l
PRINT*, 'Insira a quantidade de linhas que voce deseja'
READ*,k
l=1.42D0
z=0
a=k*2
OPEN (UNIT=1,FILE='BenzenoporLinhas.xyz')
WRITE(1,*)a
WRITE(1,*)
   DO i=1,k
    DO r=1,4
        w=r*60.0D0+30.0D0
	pi=dacos(-1.0D0)
	rad= pi*w*(2.0D0/360.0D0)

    IF (mod(i,2).eq.1) THEN
      IF (r==1) CYCLE
      IF (r==4) CYCLE
	x=(i-1)*l*sqrt(3.0D0)/(2.0D0)+l*dcos(rad)
	y=l*dsin(rad)
      WRITE (1,*) 'c',x,y,z
    ELSE IF (mod(i,2).eq.0) THEN
      IF (r==2) CYCLE
      IF (r==3) CYCLE
	x=(i-2)*l*sqrt(3.0D0)/(2.0D0)+l*dcos(rad)
	y=l*dsin(rad)
      WRITE (1,*) 'c',x,y,z
   END IF
  ! WRITE (1,*) 'c',x,y,z
   END DO
   END DO
   !WRITE (1,*) 'c',x,y,z
   CLOSE(UNIT=1)
END PROGRAM