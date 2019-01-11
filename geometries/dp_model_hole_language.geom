#NDIM
2 
#PRIM
1.0 0.0
0.0 1.0
#SUPER
NCELL 0
0 NCELL
#ORB
s0  0.0d0  0.0d0  0.0d0 #0
s1  0.5d0  0.0d0  0.0d0 #1
s2  0.0d0  0.5d0  0.0d0 #2
#HAMILT
0 1   0.5 0.0 0.0    tpd  tpd 0.0
0 1  -0.5 0.0 0.0   -tpd -tpd 0.0                                      
0 2   0.0 0.5 0.0   -tpd -tpd 0.0
0 2   0.0 -0.5 0.0   tpd  tpd 0.0
1 2  -0.5 0.5 0.0   -tpp -tpp 0.0
1 2   0.5 0.5 0.0    tpp  tpp 0.0
2 1   0.5 0.5 0.0    tpp  tpp 0.0
2 1  -0.5 0.5 0.0   -tpp -tpp 0.0
0 0   0.0 0.0 0.0    ed  ed  Ud 
1 1   0.0 0.0 0.0    ep  ep  Up
2 2   0.0 0.0 0.0    ep  ep  Up
#bcond = boundary conditions: 0 is periodic 1 is antiperiodic
bcond  = 0,0
#END

