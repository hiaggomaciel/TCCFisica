PROGRAM	Naftaleno
IMPLICIT NONE
INTEGER::r,n
DOUBLE PRECISION:: pi, rad, y, x, w, z, l,a
l=1.41D0
z=0.0D0
OPEN(UNIT=1,FILE='Naftaleno.xyz')
WRITE(1,*) 10
WRITE(1,*)
DO r=0,5
	w=r*60.0D0-30.0D0
	pi=dacos(-1.0D0)
	rad= pi*w*(2.0D0/360.0D0)
	x=l*dcos(rad)
	y=l*dsin(rad)
WRITE(1,*) 'c',x,y,z
END DO
DO r=2,5
	w=r*60.0D0-30.0D0
	pi=dacos(-1.0D0)
	rad= pi*w*(2.0D0/360.0D0)
	x=l*dcos(rad)-l*sqrt(3.0D0)
	y=l*dsin(rad)
WRITE(1,*) 'c',x,y,z
END DO
 CLOSE(UNIT=1)
END PROGRAM
