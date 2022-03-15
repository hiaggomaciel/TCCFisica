PROGRAM	flake
     IMPLICIT NONE
     INTEGER::r, n, k, i, a, j, m, ac, b, iaq, jaq, INFO, LWORK, l, t
     DOUBLE PRECISION,ALLOCATABLE:: AQ(:,:),x(:),y(:),z(:), WORK(:), W(:)
     DOUBLE PRECISION:: pi, rad, acc, a1x, a1y, p, a2x, a2y, d, x1 ,y1, z1, gamma,e_site, Ef, Ef1, Ef2, GAP
     DOUBLE PRECISION:: Energy
     CHARACTER*1,ALLOCATABLE:: atom(:)
     CHARACTER*100:: filename
       gamma=-3.0D0
       e_site=0.0D0
       acc= 1.41D0
       a1x= acc*1.5
       a1y= a1x*sqrt(3.0D0)
       a2x= acc*1.5
       a2y= -a2x*sqrt(3.0D0)
         PRINT*, 'Size of flake.'
         READ*, k
            IF(k.lt.1) PRINT*,'Enter with a number greater than zero!'
            IF(k.lt.1) STOP
              WRITE(filename,'(I2)')k
            IF(k.lt.10) filename='0'//TRIM(ADJUSTL(filename))
             filename='flake-'//TRIM(ADJUSTL(filename))//'.xyz'
            OPEN(UNIT=1,FILE=TRIM(filename))
	        a=0                                                                ! declaração do númedo de átomos
	        DO m=1,k                                                           
	           a=a+((2*k)-m)*2                                                 
	        END DO                                                             
	           a=a-(2*k-1)                                                     
	           a=a*6 
	        WRITE(1,*) a
                WRITE(1,*)
                  DO i=-(k-1),k-1                                                   !transposição vetorial dos átomos multiplicando de -(k-1)
	           DO j=-(k-1),k-1                                                  ! a (k-1) para que a estrutura se comporte como losango
	            DO r=0,5                                                        !benzeno central
	               p=r*60.0D0                                                   !
	               pi=dacos(-1.0D0)                                             !
		       rad=pi*p*(1.0D0/180.0D0)                                     !
	               x1=acc*dcos(rad)                                             !
	               y1=acc*dsin(rad)                                             !
	               z1=0.0D0                                                     !
	             IF ((i-j).le. -(k)) CYCLE                                      !para a criação do floco em si foi necessário que o programa
	             IF ((i-j).ge.  (k)) CYCLE                                      !retirasse o fim de cada losango de benzeno
	         WRITE(1,*) 'c',x1+(i+j)*a1x,y1+(i-j)*a1y,z1                        !escrevendo a estrutura por transposição
	            END DO
	           END DO
	          END DO
                 CLOSE(UNIT=1)
 
                OPEN(UNIT=1,FILE=TRIM(filename))
    
           ALLOCATE(x(a),y(a),z(a),atom(a))                                         !alocando os valores x, y e z em vetores
    
           READ (1,*) a
           READ (1,*) 
            DO i=1,a
              READ(1,*) atom(i),x(i),y(i),z(i)                                      !como x, y e z são números reias não daria pra calcular
            END DO                                                                  !a matriz, sendo obrigado a alocar os mesmo em vetor
                CLOSE (UNIT=1)

    
       WRITE(filename,'(I2)')k
         IF(k.lt.10) filename='0'//TRIM(ADJUSTL(filename))
          filename='matriz-'//TRIM(ADJUSTL(filename))//'.mat'    
            OPEN(UNIT=2,FILE=TRIM(filename))
             ALLOCATE(AQ(a,a),W(a))                                                 !após por os valores de x, y e z em vetores foi feita uma
	       DO i=1,a                                                             !situação que se pudesse calcular a distâcia entre os átomos
	        DO j=1,a                                                            !de carbono e assim escrever a matriz de tamanho axa
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

       ALLOCATE(WORK(1))
       LWORK=-1
       CALL DSYEV('N','L',a,AQ,a,W,WORK,LWORK,INFO)
       LWORK=WORK(1)
       DEALLOCATE(WORK)
       ALLOCATE(WORK(LWORK))
       CALL DSYEV('N','L',a,AQ,a,W,WORK,LWORK,INFO)
     
           DO l=1,a 
            WRITE(2,*) W(l)
           END DO
            
        CLOSE (UNIT=2)
        IF(k.lt.10) filename='0'//TRIM(ADJUSTL(filename))
         filename='info-'//TRIM(ADJUSTL(filename))//'.inf'
          OPEN(UNIT=3,FILE=TRIM(filename))
           !OPEN(UNIT=1,FILE=TRIM(filename))
             !READ(1,*) a
             !READ(1,*)
           OPEN(UNIT=2,FILE=TRIM(filename))
            
             DO t=1,a
              READ(2,*) Energy
               IF (t==(a/2)) THEN
                    t=Ef1
               ELSE IF (t==((a/2)+1)) THEN
                    t=Ef2
               END IF
              END DO
             Ef=((Ef1+Ef2)/2)
             GAP= Ef2-Ef1
           PRINT*, Ef, GAP 
        WRITE(3,*) Ef
        WRITE(3,*)
        WRITE(3,*) GAP
      
       WRITE(3,*) a/2
       CLOSE(UNIT=2)
       CLOSE(UNIT=3)

END PROGRAM