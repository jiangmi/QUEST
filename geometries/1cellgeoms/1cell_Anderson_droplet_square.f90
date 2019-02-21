program geom
! Generate geometry file for Anderson droplet on triangular lattice 
! (single supercell as the whole lattice)
! the whole lattice is a square lattice with additional diagonal hoppings so
! that triangular lattice, then Anderson droplet put in the center
! This program used PBC
implicit none

real :: x, y, N=6.0, U=4.0, p, V=1.2

! No. of atoms for each ring and their locations:
integer, dimension(5) :: Natom_ring = (/ 6,12,18,24,30 /)

! coordinates of droplet rings
real, dimension(6)  :: x1, y1
real, dimension(12) :: x2, y2
real, dimension(18) :: x3, y3
real, dimension(24) :: x4, y4
real, dimension(30) :: x5, y5
integer, allocatable :: idx_atom(:,:)  ! orbital index for droplet atoms and
                                       ! correponding metallic atoms

integer :: i, j, k, a, b, Nring=1, Nhop_droplet
character*3 :: str  
real*8, allocatable :: r(:)

write(str,'(I1)') Nring
open(unit=11,file='g_Anderson_droplet_Nring'//str,status='replace', action='write')
write(11,"(A5)") "#NDIM"
write(11,'(A1)') "2"
write(11,'(A5)') "#PRIM"

if (int(N)<10) then
  write(11,'(F3.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
  write(11,'(A3, A1, F3.1, A1, A3)') "0.0"," ",N," ","0.0"
else if (int(N)<100) then
  write(11,'(F4.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
  write(11,'(A3, A1, F4.1, A1, A3)') "0.0"," ",N," ","0.0"
else
  write(11,'(F5.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
  write(11,'(A3, A1, F5.1, A1, A3)') "0.0"," ",N," ","0.0"
endif

write(11,'(A3, A1, A3, A1, A3)') "0.0"," ","0.0"," ","2.0"
write(11,'(A6)') "#SUPER"
write(11,'(A3)') "1 0"
write(11,'(A3)') "0 1"

! orbitals: metal surface
write(11,'(A4)') "#ORB"
do j = 0, N-1
  do i = 0, N-1
    a = j*N+i
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
  end do
end do

! orbitals: Anderson droplet
Nhop_droplet = sum(Natom_ring(1:Nring))+1
allocate(idx_atom(Nhop_droplet, 2))
p = N/2   
x1 = (/ p-1, p,   p-1, p+1, p,   p+1 /)
y1 = (/ p-1, p-1, p,   p,   p+1, p+1 /)
x2 = (/ p-2, p-1, p,   p-2, p+1, p-2, p+2, p-1, p+2, p,   p+1, p+2 /)
y2 = (/ p-2, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+2 /)
x3 = (/ p-3, p-2, p-1, p,   p-3, p+1, p-3, p+2, p-3, p+3, p-2, p+3, p-1, p+3, p,   p+1, p+2, p+3 /)
y3 = (/ p-3, p-3, p-3, p-3, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+3, p+3 /)

! single atom in the center:
a = N*N
if (a<10) then
  write(str,'(I1)') int(a)
  write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",p," ",p," ","1.0"
else if (a<100) then
  write(str,'(I2)') int(a)
  write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",p," ",p," ","1.0"
    else
  write(str,'(I3)') int(a)
  write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",p," ",p," ","1.0"
endif

idx_atom(1,:) = (/ p*N+p, real(a) /)

! ring 1:
if (Nring>=1) then
  do i = 1,Natom_ring(1)
    a = N*N+i
    if (a<10) then
      write(str,'(I1)') int(a)
      write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x1(i)," ",y1(i)," ","1.0"
    else if (a<100) then
      write(str,'(I2)') int(a)
      write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x1(i)," ",y1(i)," ","1.0"
        else
      write(str,'(I3)') int(a)
      write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x1(i)," ",y1(i)," ","1.0"
    endif

    idx_atom(i+1,:) = (/ y1(i)*N+x1(i), real(a) /)
  enddo
endif

! ring 2:
if (Nring>=2) then
  do i = 1,Natom_ring(2)
    a = N*N+Natom_ring(1)+i
    if (a<10) then
      write(str,'(I1)') int(a)
      write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x2(i)," ",y2(i)," ","1.0"
    else if (a<100) then
      write(str,'(I2)') int(a)
      write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x2(i)," ",y2(i)," ","1.0"
        else
      write(str,'(I3)') int(a)
      write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x2(i)," ",y2(i)," ","1.0"
    endif

    idx_atom(i+7,:) = (/ y2(i)*N+x2(i), real(a) /)
  enddo
endif

! ring 3:
if (Nring>=3) then
  do i = 1,Natom_ring(3)
    a = N*N+Natom_ring(1)+Natom_ring(2)+i
    if (a<10) then
      write(str,'(I1)') int(a)
      write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x3(i)," ",y3(i)," ","1.0"
    else if (a<100) then
      write(str,'(I2)') int(a)
      write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x3(i)," ",y3(i)," ","1.0"
        else
      write(str,'(I3)') int(a)
      write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",x3(i)," ",y3(i)," ","1.0"
    endif

    idx_atom(i+19,:) = (/ y3(i)*N+x3(i), real(a) /)
  enddo
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
  ! PBC:
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
! PBC:
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

! hopping along (1,1) direction
do j = 0, N-2
  do i = 0, N-2
    a = j*N+i
    b = (j+1)*N+(i+1)
    if(a<10 .and. b<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
  enddo
enddo
! PBC:
do j = 0, N-2
   a = j*N+(N-1)
   b = a+1
    if(a<10 .and. b<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
end do
do i = 0, N-2
   a = (N-1)*N+i
   b = i+1
    if(a<10 .and. b<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
end do

! local hopping between Anderson droplet and metal
do i = 1, Nhop_droplet
    j = idx_atom(i,1)
    k = idx_atom(i,2)
    if(j<10 .and. k<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, F4.1, A2, F4.1, A2, A3)') j," ",k," ","0.0"," ","0.0"," ","1.0"," ",V," ",V," ","0.0"
    else if (j<10 .and. k<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, F4.1, A2, F4.1, A2, A3)') j," ",k," ","0.0"," ","0.0"," ","1.0"," ",V," ",V," ","0.0"
    else if (j<100 .and. k<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, F4.1, A2, F4.1, A2, A3)') j," ",k," ","0.0"," ","0.0"," ","1.0"," ",V," ",V," ","0.0"
    else if (j<100 .and. k<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, F4.1, A2, F4.1, A2, A3)') j," ",k," ","0.0"," ","0.0"," ","1.0"," ",V," ",V," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, F4.1, A2, F4.1, A2, A3)') j," ",k," ","0.0"," ","0.0"," ","1.0"," ",V," ",V," ","0.0"
    endif
end do

! local interaction for metal = 0
do j = 0, N-1
  do i = 0, N-1
    a = j*N+i
    if (mod(i+j,2)==0) then
      if (a<10) then
        write(11,'(I1, A1, I1, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0," ",0.0
      else if (a<100) then
        write(11,'(I2, A1, I2, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0," ",0.0
          else
        write(11,'(I3, A1, I3, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0," ",0.0
      endif
    else
      if (a<10) then
        write(11,'(I1, A1, I1, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",-0.0," ",0.0," ",0.0
      else if (a<100) then
        write(11,'(I2, A1, I2, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",-0.0," ",0.0," ",0.0
          else
        write(11,'(I3, A1, I3, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') a," ",a," ","0.0"," ","0.0" &
               ," ","0.0"," ",-0.0," ",0.0," ",0.0
      endif
    endif
  end do
end do

! local interaction for droplet = finite U
do i = 1, Nhop_droplet
    k = idx_atom(i,2)
      if (k<10) then
        write(11,'(I1, A1, I1, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') k," ",k," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0," ",U
      else if (k<100) then
        write(11,'(I2, A1, I2, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') k," ",k," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0," ",U
          else
        write(11,'(I3, A1, I3, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A2, F4.1)') k," ",k," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0," ",U
      endif
end do

! =======================================================
write(11,'(A5)') "#SYMM"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 0.0d0 1.0d0 0.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 0.0d0 0.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 0.0d0 1.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 0.0d0 -1.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 1.0d0 1.0d0 0.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 1.0d0 0.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 1.0d0 1.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F4.1,A25)') "d  ",p, " ", p, " 1.0d0 -1.0d0 1.0d0 0.0d0"

write(11,'(A4)') "#END"
end program geom
