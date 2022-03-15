PROGRAM tb

IMPLICIT none
DOUBLE PRECISION:: hopping,epsil,pi,limite
DOUBLE PRECISION:: dis_tol,extreme,gap
DOUBLE PRECISION:: distancia,Ef,U,f,kt,distancia_ij,max_charge_dif,Efinal
DOUBLE PRECISION,ALLOCATABLE:: E_up(:),E_dw(:),R(:,:),RWORK(:),E_tot(:),charge_dif(:)
DOUBLE PRECISION,ALLOCATABLE:: charge_up(:),charge_dw(:),charge_up1(:),charge_dw1(:)
COMPLEX*16,ALLOCATABLE :: H_up(:,:),H_dw(:,:),WORK(:)
INTEGER :: INFO,LWORK,nat,i,j,LIWORK,LRWORK,l,nel,first,last,passo,color,q0option
INTEGER,ALLOCATABLE :: IWORK(:),M(:,:),n_vizinhos(:)
CHARACTER*1,ALLOCATABLE:: atom(:)
CHARACTER*100:: systemlabel
      
!Program to make the calculation of the electronic structure of systems based in graphene
!We take in consideration only the interaction between the first neighbors

!Initial parameters
hopping=3.0D0    ! Same as gamma
epsil=0.0D0      !-> Site energy parameter
U=3.00D0         !-> Hubbard parameter
pi=DACOS(-1.0D0)
dis_tol=0.20D0   ! Now we identify neighbors by a distance tolerance
kt=0.025D0

PRINT*,'************************************************************************'
PRINT*,'************************************************************************'
PRINT*,'\\Electronic structure of graphitic 0D systems by Tight Binding Method//'
PRINT*,'************************************************************************'
PRINT*,'************************************************************************'
!Reading the file of the coordinates of the unit cell of the system
PRINT*,'************************************************************************'
PRINT*,'\\\\\\\\\\\\\Reading input data for coordinates and charge//////////////'
PRINT*,'************************************************************************'
OPEN(UNIT=1,FILE='input.in')                           !!!!!!!!!!!!!!!!!!!!!!!!!!
READ(1,*) systemlabel                                  !!!!!!!!!!!!!!!!!!!!!!!!!!
READ(1,*) nel                                          !!!!!!!!!!!!!!!!!!!!!!!!!!
READ(1,*) q0option                                     !!!!!!!!!!!!!!!!!!!!!!!!!!
CLOSE(UNIT=1)                                          !!!!!!!!!!!!!!!!!!!!!!!!!!
                                                       !!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN (UNIT=1, FILE=TRIM(ADJUSTL(systemlabel))//'.xyz') !!!!!!!!!!!!!!!!!!!!!!!!!!
READ (1,*) nat                                         !!!!!!!!!!!!!!!!!!!!!!!!!!
READ (1,*)                                             !!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE (atom(nat),R(nat,3))                          !!!!!!!!!!!!!!!!!!!!!!!!!!
DO i=1,nat                                             !!!!!!!!!!!!!!!!!!!!!!!!!!
  READ(1,*) atom(i),R(i,:)                             !!!!!!!!!!!!!!!!!!!!!!!!!!
END DO                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!
CLOSE(UNIT=1)                                          !!!!!!!!!!!!!!!!!!!!!!!!!!
                                                       !!!!!!!!!!!!!!!!!!!!!!!!!!
ALLOCATE(charge_up(nat),charge_dw(nat))                !!!!!!!!!!!!!!!!!!!!!!!!!!
IF(q0option.eq.0)THEN                                  !!!!!!!!!!!!!!!!!!!!!!!!!!
  charge_up=0.5D0                                      !!!!!!!!!!!!!!!!!!!!!!!!!!
  charge_dw=0.5D0                                      !!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF(q0option.eq.1)THEN                             !!!!!!!!!!!!!!!!!!!!!!!!!!
  charge_up=1.0D0                                      !!!!!!!!!!!!!!!!!!!!!!!!!!
  charge_dw=0.0D0                                      !!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE IF(q0option.eq.-1)THEN                            !!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i=1,nat                                           !!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(MOD(i,2).eq.1)THEN                              !!!!!!!!!!!!!!!!!!!!!!!!!!
      charge_up(i)=1.0D0                               !!!!!!!!!!!!!!!!!!!!!!!!!!
      charge_dw(i)=0.0D0                               !!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSE IF(MOD(i,2).eq.0)THEN                         !!!!!!!!!!!!!!!!!!!!!!!!!!
      charge_up(i)=0.0D0                               !!!!!!!!!!!!!!!!!!!!!!!!!!
      charge_dw(i)=1.0D0                               !!!!!!!!!!!!!!!!!!!!!!!!!!
    END IF                                             !!!!!!!!!!!!!!!!!!!!!!!!!!
  END DO                                               !!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE                                                   !!!!!!!!!!!!!!!!!!!!!!!!!!
  OPEN (UNIT=4,FILE=TRIM(ADJUSTL(systemlabel))//'.q')  !!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i=1,nat                                           !!!!!!!!!!!!!!!!!!!!!!!!!!
    READ(4,*) charge_up(i),charge_dw(i)                !!!!!!!!!!!!!!!!!!!!!!!!!!
  END DO                                               !!!!!!!!!!!!!!!!!!!!!!!!!!
  CLOSE (UNIT=4)                                       !!!!!!!!!!!!!!!!!!!!!!!!!!
END IF                                                 !!!!!!!!!!!!!!!!!!!!!!!!!!

! Table of Neighbors
PRINT*,'************************************************************************'
PRINT*,'\\\\\\\\\Determining the neighborhood conditions of the lattice/////////'
PRINT*,'************************************************************************'
ALLOCATE (M(nat,12),n_vizinhos(nat))
n_vizinhos=0
DO i=1,nat
  DO j=1,nat
    distancia_ij=distancia(R(i,:),R(j,:))
    IF (ABS(distancia_ij-1.42).lt.dis_tol) THEN
      n_vizinhos(i)=n_vizinhos(i)+1
      M(i,n_vizinhos(i))=j
    END IF
  END DO
END DO
      
! Self-consistency cycle
limite=0.0000001D0
passo=0
DO 
  passo=passo+1
  PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
  PRINT*, 'integration step', passo
  !Compute the Hamiltonian to spin up
  ALLOCATE (H_up(nat,nat),E_up(nat))
  H_up=0.0D0 
  DO i=1,nat
    H_up(i,i)=epsil+U*charge_dw(i)
    DO j=1,n_vizinhos(i)
      l=M(i,j)
      H_up(i,l)=-hopping
    END DO
  END DO
  ALLOCATE (WORK(1),RWORK(1),IWORK(1))
  LWORK=-1
  CALL ZHEEVD( 'V', 'U', nat, H_up, nat, E_up, WORK, LWORK,RWORK,LRWORK, IWORK, LIWORK, INFO )
  LWORK=WORK(1)
  LRWORK=RWORK(1)
  LIWORK=IWORK(1)
  DEALLOCATE(WORK,IWORK,RWORK)
  ALLOCATE (WORK(LWORK),RWORK(LRWORK),IWORK(LIWORK))
  CALL ZHEEVD( 'V', 'U', nat, H_up, nat, E_up, WORK, LWORK,RWORK,LRWORK, IWORK, LIWORK, INFO )
  DEALLOCATE (WORK,RWORK,IWORK)
  !Compute the Hamiltonian to spin down
  ALLOCATE (H_dw(nat,nat),E_dw(nat))
  H_dw=0.0D0 
  DO i=1,nat
    H_dw(i,i)=epsil+U*charge_up(i)
    DO j=1,n_vizinhos(i)
      l=M(i,j)
      H_dw(i,l)=-hopping
    END DO
  END DO
  ALLOCATE (WORK(1),RWORK(1),IWORK(1))
  LWORK=-1
  CALL ZHEEVD( 'V', 'U', nat, H_dw, nat, E_dw, WORK, LWORK,RWORK,LRWORK, IWORK, LIWORK, INFO )
  LWORK=WORK(1)
  LRWORK=RWORK(1)
  LIWORK=IWORK(1)
  DEALLOCATE(WORK,IWORK,RWORK)
  ALLOCATE (WORK(LWORK),RWORK(LRWORK),IWORK(LIWORK))
  CALL ZHEEVD( 'V', 'U', nat, H_dw, nat, E_dw, WORK, LWORK,RWORK,LRWORK, IWORK, LIWORK, INFO )
  DEALLOCATE (WORK,RWORK,IWORK)
  !ordering the E_tot vector with the quick_sort subroutine
  ALLOCATE (E_tot(2*nat))
  E_tot(1:nat)=E_up            
  E_tot(nat+1:2*nat)=E_dw     
  first=1
  last=2*nat
  CALL quick_sort(E_tot,first,last,nat)
  
  !Energy Fermi and Energy Gap
  Ef=(E_tot(nel)+E_tot(nel+1))/2.0D0
  
  
  
  !Computing the charge in each site
  ALLOCATE(charge_up1(nat),charge_dw1(nat))
  charge_up1=0.0D0
  charge_dw1=0.0D0
  DO i=1,nat
    charge_up1(i)=0.0D0
    charge_dw1(i)=0.0D0
    DO j=1,nat
      f=1.0D0/(1.0D0+DEXP((E_up(j)-Ef)/kt))
      charge_up1(i)=charge_up1(i)+f*H_up(i,j)*CONJG(H_up(i,j))
      f=1.0D0/(1.0D0+DEXP((E_dw(j)-Ef)/kt))
      charge_dw1(i)=charge_dw1(i)+f*H_dw(i,j)*CONJG(H_dw(i,j))    
    END DO
  END DO
  WRITE (*,*) 'Total charge up=',SUM(charge_up1)
  WRITE (*,*) 'Total charge dw=',SUM(charge_dw1)
  
  !maximum difference between out and in charges 
  ALLOCATE (charge_dif(2*nat))
!   DO i=1,nat
!    charge_dif(i)=ABS(charge_up1(i)-charge_up(i))
!    charge_dif(nat+i)=ABS(charge_dw1(i)-charge_dw(i))     !!!!!!!!!!!!!!!!!! Corrigido
!    WRITE (12,*) charge_up1(i), charge_dw1(i) 
!   END DO
! Mas veja essa outra forma
  charge_dif(1:nat)=ABS(charge_up1-charge_up)
  charge_dif(nat+1:2*nat)=ABS(charge_dw1-charge_dw)
  max_charge_dif=MAXVAL(charge_dif(1:2*nat))
  PRINT*, 'Maximum dq=:', max_charge_dif
  charge_up=charge_up1
  charge_dw=charge_dw1
  IF (max_charge_dif.lt.limite) EXIT
  DEALLOCATE (H_up,E_up,H_dw,E_dw,E_tot,charge_up1,charge_dw1,charge_dif)
END DO
PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'  
PRINT*, 'Fermi Energy is: Ef=',Ef      
IF(MOD(nel,2).eq.0) gap=E_tot(nel+1)-E_tot(nel)
IF(MOD(nel,2).eq.1) gap=0.0D0
PRINT*, 'Gap is: Gap=',gap      
PRINT*, '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'    
PRINT*,'************************************************************************'
PRINT*,'\\\\\\\\\\\\\\\\\\\\\\\\Writing the Energy Bands////////////////////////'
PRINT*,'************************************************************************'
OPEN (UNIT=2,FILE=TRIM(ADJUSTL(systemlabel))//'.dat')
DO i=1,nat
  WRITE (2,*) 0,E_up(i)-Ef
  WRITE (2,*) 1,E_up(i)-Ef
  WRITE (2,*)  
END DO

DO i=1,nat
  WRITE (2,*) 0,E_dw(i)-Ef
  WRITE (2,*) 1,E_dw(i)-Ef
  WRITE (2,*)
END DO
CLOSE (UNIT=2)
      

! Energia total      
Efinal=0.0D0     
DO i=1,2*nat
  f=1.0D0/(1.0D0+DEXP((E_tot(i)-Ef)/kt))
  Efinal=Efinal+f*E_tot(i)
END DO
DO i=1,nat
  Efinal=Efinal-U*charge_up(i)*charge_dw(i)
END DO
PRINT*, 'Total Energy is: Etot=',Efinal      

      
OPEN (UNIT=4,FILE=TRIM(ADJUSTL(systemlabel))//'.mag')  !!!!!!!!!!!!!!!!!!!!!!!!!!
DO i=1,nat                                           !!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(4,*) charge_up(i),charge_dw(i)                !!!!!!!!!!!!!!!!!!!!!!!!!!
END DO                                               !!!!!!!!!!!!!!!!!!!!!!!!!!
WRITE(4,*)
CLOSE (UNIT=4)                                       !!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(UNIT=8,FILE=TRIM(ADJUSTL(systemlabel))//'.tcl')  !!!!!!!!!!!!!!!!!!!!!!!!!!
!write a tcl/tk file for rendering in VMD
        IF((ABS(minval(charge_up-charge_dw))).gt.(ABS(maxval(charge_up-charge_dw))))THEN
          extreme=ABS(minval(charge_up-charge_dw))
        ELSE
          extreme=ABS(maxval(charge_up-charge_dw))
        END IF
        IF(extreme.lt.1.0D-10) extreme=1.0
        DO i=1,nat
          !color go from 33 to 1056 (1024 color)
          !see http://www.ks.uiuc.edu/Research/vmd/current/ug/node81.html#ug:topic:coloring
          color=(extreme+charge_up(i)-charge_dw(i))/(extreme*2.0D0)*1023
          write(8,*) "graphics top color " , int(color)+33
          write(8,*) "graphics top sphere {",R(i,:),"} radius .5 resolution 80"
        END DO
        write(8,*) "graphics top color " , 8
        DO i=1,nat
          DO j=i+1,nat
            if(DSQRT(SUM((R(i,:)-R(j,:))**2)).gt. 1.65) cycle
            write(8,*) "graphics top cylinder {",R(i,:),"} {",R(j,:),&
            &"} radius .25 resolution 60 filled yes"
          end DO
        end DO
        WRITE(8,*)
        CLOSE(8)
        
        
PRINT*,'************************************************************************'
PRINT*,'************************************************************************'
PRINT*,'\\\\\\\\\\\\\\\\\\\\\\\\\Finalized calculation//////////////////////////'
PRINT*,'************************************************************************'
PRINT*,'************************************************************************'

END PROGRAM tb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!Function to compute the distance between two points of the lattice

FUNCTION distancia(R1,R2)
      
      IMPLICIT none
      DOUBLE PRECISION:: distancia,R1(3),R2(3)

      distancia=SQRT(SUM((R1-R2)**2))

END FUNCTION distancia


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1



RECURSIVE SUBROUTINE quick_sort(E_tot,first,last,nat)
IMPLICIT none
     INTEGER, INTENT (IN):: nat,first,last
     DOUBLE PRECISION, INTENT (IN OUT):: E_tot(2*nat)
     INTEGER:: i,j
     DOUBLE PRECISION:: x,t
     
     x=E_tot((first+last)/2)
     i=first
     j=last
     
     DO
       DO WHILE (E_tot(i).lt.x) 
       i=i+1
       END DO
       DO WHILE (x.lt.E_tot(j)) 
       j=j-1
       END DO
       IF (i.ge.j) EXIT
       t=E_tot(i)
       E_tot(i)=E_tot(j)
       E_tot(j)=t
       i=i+1
       j=j-1
     END DO
     IF (first.lt.(i-1)) CALL quick_sort(E_tot,first,i-1,nat)
     IF ((j+1).lt.last) CALL quick_sort(E_tot,j+1,last,nat)
     
END SUBROUTINE quick_sort

