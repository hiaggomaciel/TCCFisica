PROGRAM flake
IMPLICIT NONE
INTEGER::r,k,i,a,j,m
DOUBLE PRECISION:: pi,rad,acc,a1x,a1y,p,a2x,a2y,x1,y1,z1
CHARACTER*100:: filename


! Parâmetros de rede
acc= 1.42D0
a1x= acc*1.5
a1y= a1x*sqrt(3.0D0)
a2x= acc*1.5
a2y= -a2x*sqrt(3.0D0)

! Tamanho do flake
PRINT*, 'Size of flake.'
READ*, k
IF(k.lt.1) PRINT*,'Enter with a number greater than zero!'
IF(k.lt.1) STOP

! Cálculo do número de átomos
a=0
DO m=1,k                                                           
  a=a+((2*k)-m)*2                                                 
END DO                                                             
a=a-(2*k-1)                                                     
a=a*6 

! Criando arquivo de entrada
WRITE(filename,'(I2)')k
IF(k.lt.10) filename='0'//TRIM(ADJUSTL(filename))
filename='flake-'//TRIM(ADJUSTL(filename))//'.xyz'
OPEN(UNIT=1,FILE=TRIM(filename))
WRITE(1,*) a
WRITE(1,*)
                
!Translação dos átomos da célula unitária
DO i=-(k-1),k-1 ! i*a1, com i no intervalo [-(k-1),(k-1)]
  DO j=-(k-1),k-1 ! j*a2, com j no intervalo [-(k-1),(k-1)]
    DO r=0,5
      p=r*60.0D0
      pi=dacos(-1.0D0)
      rad=pi*p*(1.0D0/180.0D0)
      x1=acc*dcos(rad)
      y1=acc*dsin(rad)
      z1=0.0D0
      IF ((i-j).le. -(k)) CYCLE ! Átomos fora do hexagono... 
      IF ((i-j).ge.  (k)) CYCLE ! ...são eliminados.
      WRITE(1,*) 'c',x1+(i+j)*a1x,y1+(i-j)*a1y,z1
    END DO
  END DO
END DO
WRITE(1,*)
CLOSE(UNIT=1)
 
END PROGRAM
