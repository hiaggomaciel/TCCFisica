PROGRAM trigonometric

IMPLICIT NONE

DOUBLE PRECISION:: pi,x, rad, z, w
PRINT*, 'Insira o valor de um angulo qualquer para ter o seu seno e cosseno.'
READ*, x

pi=dacos(-1.0D0)

rad= pi*x*(2.0D0/360.0D0)

z= dsin(rad)

w= dcos(rad)


PRINT*, 'Seno =',z
PRINT*, 'Cosseno =', w

END PROGRAM
