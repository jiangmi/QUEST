program geom
! Generate geometry file for Anderson droplet on triangular lattice:
! (single supercell as the whole lattice)
! the whole unit cell is a triangular lattice, 
! then cut the three corner sites to form a hexogonal as Fig.2a in Morr's paper (quantum engineered Kondo lattices)
! This program used open boundary
implicit none

! Variable declarations:

real :: pa = 0.8660254  ! a little larger than sqrt(3)
real :: U=4.0, V=1.2, xmin, xmax, x, y, p, q

! No. of atoms for each ring and their locations:
integer, dimension(5) :: Natom_ring = (/ 6,12,18,24,30 /)

! coordinates of droplet rings
real :: x0, y0
real, dimension(6)  :: x1, y1
real, dimension(12) :: x2, y2
real, dimension(18) :: x3, y3
real, dimension(24) :: x4, y4
real, dimension(30) :: x5, y5
integer, allocatable :: idx_atom(:,:)  ! orbital index for droplet atoms and
                                       ! correponding metallic atoms

integer :: i, j, k, m, cnt=0, a, b, N=6, Nring=0, Nmetal, Nhop_droplet, Naj
character*3 :: str 
character*1 :: s1

! map from coordinate (x,y) into index of atom
integer, allocatable :: xy2idx(:,:)   
allocate(xy2idx(N, N))

write(s1,'(I1)') Nring
write(str,'(I3)') N
open(unit=11,file='g_Anderson_droplet_Nring'//s1//'_N'//adjustl(str),status='replace', action='write')
write(11,"(A5)") "#NDIM"
write(11,'(A1)') "2"
write(11,'(A5)') "#PRIM"

y = N*pa
if (N<10) then
  write(11,'(F3.1, A1, A3, A1, A3)') real(N),"","0.0"," ","0.0"
  write(11,'(F3.1, A1, F10.7, A1, A3)') 0.5*N," ",y," ","0.0"
else if (N/2<10 .and. N<100) then
  write(11,'(F4.1, A1, A3, A1, A3)') real(N),"","0.0"," ","0.0"
  write(11,'(F3.1, A1, F10.7, A1, A3)') 0.5*N," ",y," ","0.0"
else if (N/2<100 .and. N<100) then
write(11,'(F4.1, A1, A3, A1, A3)') real(N),"","0.0"," ","0.0"                                       
  write(11,'(F4.1, A1, F10.7, A1, A3)') 0.5*N," ",y," ","0.0"
endif

write(11,'(A3, A1, A3, A1, A3)') "0.0"," ","0.0"," ","2.0"
write(11,'(A6)') "#SUPER"
write(11,'(A3)') "1 0"
write(11,'(A3)') "0 1"

! orbitals: metal surface, (i,j) is coordinate of acquired atom
write(11,'(A4)') "#ORB"
do j = 0, N/3
  xmin = N/3.0-j*0.5
  xmax = N/1.5+j*0.5

  do i = 0, N/3+j
    ! record the map from (x,y) to atom/orbital index
    xy2idx(i,j) = cnt
    !write(*,*) i, " ", j, " ", xy2idx(i,j)  
  
    x = xmin+i
    if (cnt <10) then
      write(str,'(I1)') cnt
      write(11,'(A2, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",j*pa," ","0.0"
    else if (cnt <100) then
      write(str,'(I2)') cnt
      write(11,'(A3, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",j*pa," ","0.0"
    else
      write(str,'(I3)') cnt
      write(11,'(A4, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",j*pa," ","0.0"
    endif

    cnt = cnt+1
  end do
end do
do j = N/3+1, 2*N/3
  xmin =   N/6.0+(j-N/3)*0.5
  xmax = 5*N/6.0-(j-N/3)*0.5

  do i = 0, 2*N/3-(j-N/3-1)-1
    ! record the map from (x,y) to atom/orbital index
    xy2idx(i,j) = cnt
    !write(*,*) i, " ", j, " ", xy2idx(i,j)
    
    x = xmin+i

    if (cnt <10) then
      write(str,'(I1)') cnt
      write(11,'(A2, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",j*pa," ","0.0"
    else if (cnt <100) then
      write(str,'(I2)') cnt
      write(11,'(A3, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",j*pa," ","0.0"
    else
      write(str,'(I3)') cnt
      write(11,'(A4, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",j*pa," ","0.0"
    endif

    cnt = cnt+1
  end do
end do

! assert the total number of metallic sites should be N*N/3+N+1
Nmetal = N*N/3+N+1
if (cnt/=Nmetal) then
  write(*,*) 'There was an error of index of metallic sites!'
  stop
endif

! orbitals: Anderson droplet
Nhop_droplet = sum(Natom_ring(1:Nring))+1
allocate(idx_atom(Nhop_droplet, 2))
p = N/3

! following denotes the *-th atom at jth line
x0 = p
y0 = p
x1 = (/ p-1, p,   p-1, p+1, p-1, p   /)
y1 = (/ p-1, p-1, p,   p,   p+1, p+1 /)
x2 = (/ p-2, p-1, p,   p-2, p+1, p-2, p+2, p-2, p+1, p-2, p-1, p   /)
y2 = (/ p-2, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+2 /)
x3 = (/ p-3, p-2, p-1, p,   p-3, p+1, p-3, p+2, p-3, p+3, p-3, p+2, p-3, p+1, p-3, p-2, p-1, p   /)
y3 = (/ p-3, p-3, p-3, p-3, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+3, p+3 /)

! single atom in the center:
if (cnt<10) then
  write(str,'(I1)') cnt
  write(11,'(A2, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",N*0.5," ",p*pa," ","1.0"
else if (cnt<100) then
  write(str,'(I2)') cnt
  write(11,'(A3, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",N*0.5," ",p*pa," ","1.0"
else
  write(str,'(I3)') cnt
  write(11,'(A4, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",N*0.5," ",p*pa," ","1.0"
endif

idx_atom(1,:) = (/ xy2idx(N/3,N/3), cnt /)
write(*,*) xy2idx(N/3,N/3), cnt
cnt = cnt+1

! ring 1:
if (Nring>=1) then
  do i = 1,Natom_ring(1)
    if (y1(i)<=N/3) then
      xmin = N/3.0-y1(i)*0.5
    else
      xmin = N/6.0+(y1(i)-N/3)*0.5
    endif

    x = xmin + x1(i)
    y = y1(i)*pa

    if (cnt<10) then
      write(str,'(I1)') cnt
      write(11,'(A2, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    else if (cnt<100) then
      write(str,'(I2)') cnt
      write(11,'(A3, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    else
      write(str,'(I3)') cnt
      write(11,'(A4, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    endif

    idx_atom(i+1,:) = (/ xy2idx(x1(i),y1(i)), cnt /)
    write(*,*) xy2idx(x1(i),y1(i)), cnt
    cnt = cnt+1
  enddo
endif

! ring 2:
if (Nring>=2) then
  do i = 1,Natom_ring(2)
    if (y2(i)<=N/3) then
      xmin = N/3.0-y2(i)*0.5
    else
      xmin = N/6.0+(y2(i)-N/3)*0.5
    endif

    x = xmin + x2(i)
    y = y2(i)*pa

    if (cnt<10) then
      write(str,'(I1)') cnt
      write(11,'(A2, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    else if (cnt<100) then
      write(str,'(I2)') cnt
      write(11,'(A3, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    else
      write(str,'(I3)') cnt
      write(11,'(A4, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    endif

    idx_atom(i+7,:) = (/ xy2idx(x2(i),y2(i)), cnt /)
    write(*,*) xy2idx(x2(i),y2(i)), cnt
    cnt = cnt+1
  enddo
endif

! ring 3:
if (Nring>=3) then
  do i = 1,Natom_ring(3)
    if (y3(i)<=N/3) then
      xmin = N/3.0-y3(i)*0.5
    else
      xmin = N/6.0+(y3(i)-N/3)*0.5
    endif

    x = xmin + x3(i)
    y = y3(i)*pa

    if (cnt<10) then
      write(str,'(I1)') cnt
      write(11,'(A2, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    else if (cnt<100) then
      write(str,'(I2)') cnt
      write(11,'(A3, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    else
      write(str,'(I3)') cnt
      write(11,'(A4, A1, F4.1, A2, F10.7, A2, A3)') 's'//str," ",x," ",y," ","1.0"
    endif

    idx_atom(i+19,:) = (/ xy2idx(x3(i),y3(i)), cnt /)
    write(*,*) xy2idx(x3(i),y3(i)), cnt
    cnt = cnt+1
  enddo
endif

write(11,'(A30)') "#HAMILT            tup  tdn  U"
! hopping along x direction
do j = 0, 2*N/3
  ! get the No. of atoms for jth line
  if (j<=N/3) then
    Naj = N/3+j
  else
    Naj = 2*N/3-(j-N/3-1)-1
  endif

  do i = 0, Naj-1
    a = xy2idx(i,  j)
    b = xy2idx(i+1,j)

    if(a<10 .and. b<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
  enddo
end do

! hopping along 0.5*(1,sqrt(3)) direction
do j = 0, 2*N/3-1
  ! get the No. of atoms needed for hopping for jth line
  if (j<N/3) then
    Naj = N/3+j
  else
    Naj = 2*N/3-(j-N/3-1)-2
  endif

  do i = 0, Naj
    a = xy2idx(i,j)
    if (j<N/3) then
      b = xy2idx(i+1,j+1)
    else
      b = xy2idx(i,j+1)
    endif

    if(a<10 .and. b<10) then
      write(11,'(I1, A1, I1, A2, A3, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
  enddo
end do

! hopping along (-1/2,sqrt(3)/2) direction
do j = 0, 2*N/3-1
  ! get the No. of atoms needed for hopping for jth line
  if (j<N/3) then
    Naj = N/3+j
    k = 0
    m = Naj
  else
    Naj = 2*N/3-(j-N/3-1)-2
    k = 1
    m = Naj+1
  endif

  do i = k, m
    a = xy2idx(i,j)
    if (j<N/3) then
      b = xy2idx(i,j+1)
    else
      b = xy2idx(i-1,j+1)
    endif

    if(a<10 .and. b<10) then
      write(11,'(I1, A1, I1, A2, A4, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A4, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A4, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A4, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A4, A1, F10.7, A1, A3, A2, A3, A2, A3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
  enddo
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
do a = 0, Nmetal-1
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
x = N/2
y = p*pa
write(11,'(A3,F4.1, A2, F10.7,A24)') "c6 ",x, " ", y, " 0.0d0 0.0d0 0.0d0 1.0d0"
write(11,'(A3,F4.1, A2, F10.7,A24)') "c6 ",x, " ", y, " 1.0d0 0.0d0 0.0d0 1.0d0"
write(11,'(A3,F4.1, A2, F10.7,A24)') "d  ",x, " ", y, " 0.0d0 0.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F10.7,A24)') "d  ",x, " ", y, " 1.0d0 0.0d0 1.0d0 0.0d0"

write(11,'(A4)') "#END"
end program geom
