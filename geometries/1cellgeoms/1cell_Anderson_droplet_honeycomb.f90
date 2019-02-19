program geom
! Generate geometry file for Anderson droplet (single supercell as the whole lattice)
! the whole lattice is a honeycomb lattice as shown in Dirk Morr's paper
! then Anderson droplet put in the center
implicit none

! Variable declarations:

real :: s3 = 1.733  ! a little larger than sqrt(3)
real :: x, y, N=4.0, Nring=0.0, delta=0.0, U=4.0, b
integer :: i, j, k, a, flag
character*3 :: str  
character*18 :: s1, s2, s3
real*8, allocatable :: r(:)
real*8, allocatable :: ran(:)

open(unit=11,file='g_Anderson_droplet',status='replace', action='write')
write(11,"(A5)") "#NDIM"
write(11,'(A1)') "2"
write(11,'(A5)') "#PRIM"

! Ny should be a little larger than sqrt(3)/2*N
Ny = s3/2.0*N  
if (int(N)<10) then
  write(11,'(F3.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
  write(11,'(A3, A1, F3.1, A1, A3)') "0.0"," ",Ny," ","0.0"
else if (int(N)<100) then
  write(11,'(F4.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
  write(11,'(A3, A1, F4.1, A1, A3)') "0.0"," ",Ny," ","0.0"
else
  write(11,'(F5.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
  write(11,'(A3, A1, F5.1, A1, A3)') "0.0"," ",Ny," ","0.0"
endif

write(11,'(A3, A1, A3, A1, A3)') "0.0"," ","0.0"," ","2.0"
write(11,'(A6)') "#SUPER"
write(11,'(A3)') "1 0"
write(11,'(A3)') "0 1"
! orbitals: metal surface
write(11,'(A4)') "#ORB"
do j = 0, Ny-1
  do i = 0, N-1
    a = j*N+i
    if ((i*i+j*j)<=N*N) then
      if (a<10) then
        write(str,'(I1)') int(a)
        write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"
      else if (a<100) then
        write(str,'(I2)') int(a)
        write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"
      else
        write(str,'(I3)') int(a)
        write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"
      endif
    endif
  end do
end do

! orbitals: Anderson droplet
if (Nring==0) then
    a = N*N
    b = N/2
    if (a<10) then
      write(str,'(I1)') int(a)
      write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(b)," ",real(b)," ","1.0"
    else if (a<100) then
      write(str,'(I2)') int(a)
      write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(b)," ",real(b)," ","1.0"
        else
      write(str,'(I3)') int(a)
      write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(b)," ",real(b)," ","1.0"
    endif
endif

write(11,'(A30)') "#HAMILT            tup  tdn  U"
! hopping along x direction
do i = 1, N**2
  if (mod(i,int(N))/=0) then
    if(i<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i==10) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i==100) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
  else
    if(i-1<10 .and. i-int(N)<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-int(N)," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<10 .and. i-int(N)==10) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-int(N)," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<100 .and. i-int(N)<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-int(N)," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<100 .and. i-int(N)==100) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-int(N)," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-int(N)," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
  endif
end do

! hopping along y direction
do i = 1, N**2-N
    if(i-1<10 .and. i+int(N)-1<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<10 .and. i+int(N)-1<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<100 .and. i+int(N)-1<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<100 .and. i+int(N)-1==100) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
end do
do i = N**2-N+1, N**2
    if(i-1<10 .and. i-1-int(N*(N-1))<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1-int(N*(N-1))," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<10 .and. i-1-int(N*(N-1))<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1-int(N*(N-1))," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<100 .and. i-1-int(N*(N-1))<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1-int(N*(N-1))," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i-1<100 .and. i-1-int(N*(N-1))==100) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1-int(N*(N-1))," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1-int(N*(N-1))," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
end do

! local interaction
! consider addition potential: staggered for one (SDW) or two spins (CDW)
do j = 0, N-1
  do i = 0, N-1
    a = j*N+i
    if (mod(i+j,2)==0) then
      if (a<10) then
        write(11,'(I1, A1, I1, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",delta," ",delta," ",U
      else if (a<100) then
        write(11,'(I2, A1, I2, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",delta," ",delta," ",U
          else
        write(11,'(I3, A1, I3, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",delta," ",delta," ",U
      endif
    else
      if (a<10) then
        write(11,'(I1, A1, I1, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",-delta," ",-delta," ",U
      else if (a<100) then
        write(11,'(I2, A1, I2, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",-delta," ",-delta," ",U
          else
        write(11,'(I3, A1, I3, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",-delta," ",-delta," ",U
      endif
    endif
  end do
end do
write(11,'(A4)') "#END"
end program geom
