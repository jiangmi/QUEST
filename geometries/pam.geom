#NDIM
2 
#PRIM
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 2.0
#SUPER
NCELL 0
0 NCELL
#ORB
s0  0.0d0  0.0d0  0.0d0 #0
s1  0.0d0  0.0d0  1.0d0 #1
#HAMILT
0 0   1.0 0.0 0.0   1.0  1.0  0.0
0 0   0.0 1.0 0.0   1.0  1.0  0.0
0 1   0.0 0.0 1.0   Vval Vval 0.0
0 1   1.0 0.0 1.0   Vpval Vpval 0.0
0 1   0.0 1.0 1.0   Vpval Vpval 0.0
0 1  -1.0 0.0 1.0   Vpval Vpval 0.0
0 1   0.0 -1.0 1.0  Vpval Vpval 0.0
0 0   0.0 0.0 0.0   0.0  0.0  0.0 
1 1   0.0 0.0 0.0   0.0  0.0  Uval
#SYMM
d  0.0d0 0.0d0 0.0d0 1.0d0 0.0d0 0.d0
d  0.0d0 0.0d0 0.0d0 0.0d0 1.0d0 0.d0
d  0.0d0 0.0d0 1.0d0 1.0d0 0.0d0 0.d0
d  0.0d0 0.0d0 1.0d0 0.0d0 1.0d0 0.d0
c4 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 1.d0
c4 0.0d0 0.0d0 1.0d0 0.0d0 0.0d0 1.d0
#bcond = boundary conditions: 0 is periodic 1 is antiperiodic
bcond  = 0,0
#END

