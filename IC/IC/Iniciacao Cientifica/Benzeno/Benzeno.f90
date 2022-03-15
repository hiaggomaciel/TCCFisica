PROGRAM	Benzeno
IMPLICIT NONE
INTEGER::r, n, b, k, a
DOUBLE PRECISION:: pi, rad, y, x, w, z, acc,d
INTEGER,ALLOCATABLE:: aq(:,:)
d=0.0D0
acc=1.41D0
z=0.0D0
a=6
OPEN(UNIT=1,FILE='Benzeno.xyz')
OPEN(UNIT=2,FILE='Matriz')
WRITE(1,*) a
WRITE(1,*)
ALLOCATE(aq(a,a))
      DO r=0,5
	w=r*60.0D0-30.0D0
	pi=dacos(-1.0D0)
	rad=pi*w*(1.0D0/180.0D0)
	y=acc*dsin(rad)
	x=acc*dcos(rad)


       d=sqrt((x**2)+(y**2))-d
        IF (d==acc) THEN
        b=2
        ELSE IF (d==0) THEN
        b=1
        ELSE
        b=0
        END IF
     WRITE(2,*) b
     WRITE(1,*) 'c',x,y,z
      END DO
END PROGRAM Benzeno
