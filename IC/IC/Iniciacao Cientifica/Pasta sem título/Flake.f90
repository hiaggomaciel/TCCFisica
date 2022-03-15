PROGRAM	Flake
IMPLICIT NONE
INTEGER::r, n, k, i, a, j, m, ac, b, iaq, jaq, INFO, LWORK, l, t
DOUBLE PRECISION,ALLOCATABLE:: AQ(:,:),x(:),y(:),z(:), WORK(:), W(:),Energy(:)
DOUBLE PRECISION:: pi, rad, acc, a1x, a1y, p, a2x, a2y, d, x1 ,y1, z1, gamma,e_site, Ef, Ef1, Ef2, GAP
CHARACTER*1,ALLOCATABLE:: atom(:)

      d=0
      acc=1.41D0
      z=0.0D0
      a1x=acc*1.5
      a1y=a1x*sqrt(3.0D0)
      a2x=acc*1.5
      a2y=-a2x*sqrt(3.0D0)
       PRINT*, 'Tamanho do floco'
       READ*, k
       OPEN(UNIT=1,FILE='Flake.xyz')
        a=0
            DO m=1,k
                a=a+((2*k)-m)*2
            END DO
                a=a-(2*k-1)
                a=a*6
       WRITE(1,*) a
       WRITE(1,*)
              DO i=0,k-1
                DO j=0,k-1
                  DO r=0,5
                     w=r*60.0D0
                     pi=dacos(-1.0D0)
                     rad=pi*p*(1.0D0/180.0D0)
                     x=acc*dcos(rad)
                     y=acc*dsin(rad)
                IF (j>i) CYCLE
       WRITE(1,*) 'c',x+(i+j)*a1x,y+(i-j)*a1y,z
                   END DO
                END DO
              END DO
       CLOSE(UNIT=1)
       
       ALLOCATE(x(a),y(a),z(a),atom(a))

                READ (1,*) a
                READ (1,*)
               
               DO i=1,a
                 READ(1,*) atom(i),x(i),y(i),z(i)
               END DO
             CLOSE (UNIT=1)

             OPEN (UNIT=2,FILE='matriz.mat')
             ALLOCATE(AQ(a,a),W(a))
	     
	      DO i=1,a
	        DO j=1,a
	          d=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
	        IF ((d.gt.1.4D0).and.(d.lt.1.45D0)) THEN
	           AQ(i,j)=gamma
	         ELSE IF (d.lt.1.0D0) THEN
	           AQ(i,j)=e_site
	         ELSE
	           AQ(i,j)=0.0D0
	         END IF
	         END DO 
                END DO
WRITE(2,*) AQ
       
       CLOSE(UNIT=2)
END PROGRAM