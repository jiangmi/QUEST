program geom
! Generate geometry file for Anderson droplet on triangular lattice:
! (single supercell as the whole lattice)
! the whole unit cell is a triangular lattice, 
! then cut the three corner sites to form a hexogonal as Fig.2a in Morr's paper (quantum engineered Kondo lattices)
! This program used open boundary
implicit none

! Variable declarations:

real :: pa = 0.8660254  ! sqrt(3)/2
real :: tt=1.0, xmin, xmax, x, y, p, q

! No. of atoms for each ring and their locations:
integer, dimension(9) :: Natom_ring = (/ 6,12,18,24,30,36,42,48,54 /)

! coordinates of droplet rings
character*1 :: AB = 'B'
integer :: dr = 4   ! control the distance between rings
integer :: N=30, Nring=1, pp
real :: x0, y0
real, dimension(6)  :: x1, y1
real, dimension(12) :: x2, y2
real, dimension(18) :: x3, y3
real, dimension(24) :: x4, y4
real, dimension(30) :: x5, y5
real, dimension(36) :: x6, y6
real, dimension(42) :: x7, y7
real, dimension(48) :: x8, y8
real, dimension(54) :: x9, y9
integer, allocatable :: idx_atom(:,:)  ! orbital index for droplet atoms and
                                       ! correponding metallic atoms

integer :: i, j, k, m, cnt=0, a, b, Nmetal, Nimp, Naj
character*3 :: str 
character*1 :: s1, s2

! map from coordinate (x,y) into index of atom
integer, allocatable :: xy2idx(:,:)   
allocate(xy2idx(0:N-1, 0:N-1))

pp = Nring*dr
Nimp = sum(Natom_ring(dr:pp:dr))+1
write(*,*) 'No. of droplet impurities = ', Nimp
allocate(idx_atom(1:Nimp, 2))
xy2idx = 0
idx_atom = 0

write(s1,'(I1)') Nring
write(s2,'(I1)') dr
write(str,'(I3)') N/3
open(unit=11,file='g_Anderson_droplet_'//AB//s2//'_Nr'//s1//'_L'//adjustl(str),status='replace', action='write')
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
write(*,*) N/3
do j = 0, N/3
  xmin = N/3.0-j*0.5
  xmax = N/1.5+j*0.5

  do i = 0, N/3+j
    ! record the map from (x,y) to atom/orbital index
    xy2idx(i,j) = cnt
    write(*,*) i, " ", j, " ", xy2idx(i,j)  
  
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
    write(*,*) i, " ", j, " ", xy2idx(i,j)
    
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
p = N/3

! following denotes the *-th atom at jth line
x0 = p
y0 = p
if (AB=='A') then
  x1 = (/ p-1, p,   p-1, p+1, p-1, p   /)
  y1 = (/ p-1, p-1, p,   p,   p+1, p+1 /)
  x2 = (/ p-2, p-1, p,   p-2, p+1, p-2, p+2, p-2, p+1, p-2, p-1, p   /)
  y2 = (/ p-2, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+2 /)
  x3 = (/ p-3, p-2, p-1, p,   p-3, p+1, p-3, p+2, p-3, p+3, p-3, p+2, p-3, p+1, p-3, p-2, p-1, p   /)
  y3 = (/ p-3, p-3, p-3, p-3, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+3, p+3 /)
  x4 = (/ p-4, p-3, p-2, p-1, p,   p-4, p+1, p-4, p+2, p-4, p+3, p-4, p+4, p-4, p+3, p-4, p+2, p-4, p+1, p-4, p-3, p-2, p-1, p   /)
  y4 = (/ p-4, p-4, p-4, p-4, p-4, p-3, p-3, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+4, p+4, p+4, p+4, p+4 /)
  x5 = (/ p-5, p-4, p-3, p-2, p-1, p,   p-5, p+1, p-5, p+2, p-5, p+3, p-5, p+4, p-5, p+5, p-5, p+4, p-5, p+3, p-5, p+2, p-5, p+1, p-5, p-4, p-3, p-2, p-1, p   /)
  y5 = (/ p-5, p-5, p-5, p-5, p-5, p-5, p-4, p-4, p-3, p-3, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+4, p+4, p+5, p+5, p+5, p+5, p+5, p+5 /)

  x6 = (/ p-6, p-5, p-4, p-3, p-2, p-1, p,   p-6, p+1, p-6, p+2, p-6, p+3, p-6, p+4, p-6, p+5, p-6, &
        p+6, p-6, p+5, p-6, p+4, p-6, p+3, p-6, p+2, p-6, p+1, p-6, p-5, p-4, p-3, p-2, p-1, p /)

  y6 = (/ p-6, p-6, p-6, p-6, p-6, p-6, p-6, p-5, p-5, p-4, p-4, p-3, p-3, p-2, p-2, p-1, p-1, p,   &
        p,   p+1, p+1, p+2, p+2, p+3, p+3, p+4, p+4, p+5, p+5, p+6, p+6, p+6, p+6, p+6, p+6, p+6 /)

  x7 = (/ p-7, p-6, p-5, p-4, p-3, p-2, p-1, p,   p-7, p+1, p-7, p+2, p-7, p+3, p-7, p+4, p-7, p+5, p-7, p+6, &
        p-7, p+7, p-7, p+6, p-7, p+5, p-7, p+4, p-7, p+3, p-7, p+2, p-7, p+1, p-7, p-6, p-5, p-4, p-3, p-2, p-1, p /)

  y7 = (/ p-7, p-7, p-7, p-7, p-7, p-7, p-7, p-7, p-6, p-6, p-5, p-5, p-4, p-4, p-3, p-3, p-2, p-2, p-1, p-1,   &
        p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+4, p+4, p+5, p+5, p+6, p+6, p+7, p+7, p+7, p+7, p+7, p+7, p+7, p+7 /)
    
  x8 = (/ p-8, p-7, p-6, p-5, p-4, p-3, p-2, p-1, p,   p-8, p+1, p-8, p+2, p-8, p+3, p-8, p+4, p-8, p+5, p-8, p+6, &
        p-8, p+7, p-8, p+8, p-8, p+7, p-8, p+6, p-8, p+5, p-8, p+4, p-8, p+3, p-8, p+2, p-8, p+1, &
        p-8, p-7, p-6, p-5, p-4, p-3, p-2, p-1, p /)
    
  y8 = (/ p-8, p-8, p-8, p-8, p-8, p-8, p-8, p-8, p-8, p-7, p-7, p-6, p-6, p-5, p-5, p-4, p-4, p-3, p-3, p-2, p-2,  &
        p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+4, p+4, p+5, p+5, p+6, p+6, p+7, p+7, &
        p+8, p+8, p+8, p+8, p+8, p+8, p+8, p+8, p+8 /)
    
  x9 = (/ p-9, p-8, p-7, p-6, p-5, p-4, p-3, p-2, p-1, p,   p-9, p+1, p-9, p+2, p-9, p+3, p-9, p+4, p-9, p+5, &
        p-9, p+6, p-9, p+7, p-9, p+8, p-9, p+9, p-9, p+8, p-9, p+7, p-9, p+6, p-9, p+5, p-9, p+4, p-9, p+3, &
        p-9, p+2, p-9, p+1, p-9, p-8, p-7, p-6, p-5, p-4, p-3, p-2, p-1, p /)
    
  y9 = (/ p-9, p-9, p-9, p-9, p-9, p-9, p-9, p-9, p-9, p-9, p-8, p-8, p-7, p-7, p-6, p-6, p-5, p-5, p-4, p-4, &
        p-3, p-3, p-2, p-2, p-1, p-1, p,   p,   p+1, p+1, p+2, p+2, p+3, p+3, p+4, p+4, p+5, p+5, p+6, p+6, &
        p+7, p+7, p+8, p+8, p+9, p+9, p+9, p+9, p+9, p+9, p+9, p+9, p+9, p+9 /)
elseif (AB=='B') then
  x1 = (/ p-1, p-2, p+1, p-2, p+1, p-1 /)
  y1 = (/ p-2, p-1, p-1, p+1, p+1, p+2 /)
  x2 = (/ p-2, p-3, p,   p-4, p+2, p-3, p+3, p-4, p+2, p-3, p,   p-2 /)
  y2 = (/ p-4, p-3, p-3, p-2, p-2, p,   p,   p+2, p+2, p+3, p+3, p+4 /)
  x3 = (/ p-3, p-4, p-1, p-5, p+1, p-6, p+3, p-5, p+4, p-5, p+4, p-6, p+3, p-5, p+1, p-4, p-1, p-3 /)
  y3 = (/ p-6, p-5, p-5, p-4, p-4, p-3, p-3, p-1, p-1, p+1, p+1, p+3, p+3, p+4, p+4, p+5, p+5, p+6 /)
  x4 = (/ p-4, p-5, p-2, p-6, p,   p-7, p+2, p-8, p+4, p-7, p+5, p-6, p+6, p-7, p+5, p-8, p+4, p-7, p+2, p-6, p,   p-5, p-2, p-4 /)
  y4 = (/ p-8, p-7, p-7, p-6, p-6, p-5, p-5, p-4, p-4, p-2, p-2, p,   p,   p+2, p+2, p+4, p+4, p+5, p+5, p+6, p+6, p+7, p+7, p+8 /)
  x5 = (/ p-5,  p-6, p-3, p-7, p-1, p-8, p+1, p-9, p+3, p-10, p+5, p-9, p+6, p-8, p+7, p-8, p+7, p-9, p+6, p-10, p+5, p-9, p+3, p-8, p+1, p-7, p-1, p-6, p-3, p-5  /)
  y5 = (/ p-10, p-9, p-9, p-8, p-8, p-7, p-7, p-6, p-6, p-5,  p-5, p-3, p-3, p-1, p-1, p+1, p+1, p+3, p+3, p+5,  p+5, p+6, p+6, p+7, p+7, p+8, p+8, p+9, p+9, p+10 /)
    
  x6 = (/ p-6,  p-7,  p-4,  p-8,  p-2,  p-9, p,   p-10, p+2,  p-11, p+4,  p-12, p+6, p-11, p+7, p-10,  p+8, &
        p-9,  p+9,  p-10, p+8,  p-11, p+7, p-12,p+6,  p-11, p+4,  p-10, p+2,  p-9, p,    p-8,  p-2,  p-7,  p-4, p-6 /)

  y6 = (/ p-12, p-11, p-11, p-10, p-10, p-9, p-9, p-8,  p-8,  p-7,  p-7,  p-6,  p-6, p-4,  p-4, p-2,  p-2,  &
        p,    p,    p+2,  p+2,  p+4,  p+4, p+6, p+6,  p+7,  p+7,  p+8,  p+8,  p+9, p+9,  p+10, p+10, p+11, p+11, p+12 /)
endif

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
write(*,*) 'ring 0:'
write(*,*) xy2idx(N/3,N/3), cnt
cnt = cnt+1

if (dr==1) then
    ! ring 1:
    write(*,*) 'ring 1:'
    if (Nring>=1) then
      do i = 1,Natom_ring(1)
        if (y1(i)<=N/3) then
          xmin = N/3.0-y1(i)*0.5
        else
          xmin = N/6.0+(y1(i)-N/3)*0.5
        endif

        x = xmin + x1(i)
        y = y1(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+1,:) = (/ xy2idx(x1(i),y1(i)), cnt /)
        write(*,*) xy2idx(x1(i),y1(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 2:
    write(*,*) 'ring 2:'
    if (Nring>=2) then
      do i = 1,Natom_ring(2)
        if (y2(i)<=N/3) then
          xmin = N/3.0-y2(i)*0.5
        else
          xmin = N/6.0+(y2(i)-N/3)*0.5
        endif

        x = xmin + x2(i)
        y = y2(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+7,:) = (/ xy2idx(x2(i),y2(i)), cnt /)
        write(*,*) xy2idx(x2(i),y2(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 3:
    write(*,*) 'ring 3:'
    if (Nring>=3) then
      do i = 1,Natom_ring(3)
        if (y3(i)<=N/3) then
          xmin = N/3.0-y3(i)*0.5
        else
          xmin = N/6.0+(y3(i)-N/3)*0.5
        endif

        x = xmin + x3(i)
        y = y3(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+19,:) = (/ xy2idx(x3(i),y3(i)), cnt /)
        write(*,*) xy2idx(x3(i),y3(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 4:
    write(*,*) 'ring 4:'
    if (Nring>=4) then
      do i = 1,Natom_ring(4)
        if (y4(i)<=N/3) then
          xmin = N/3.0-y4(i)*0.5
        else
          xmin = N/6.0+(y4(i)-N/3)*0.5
        endif

        x = xmin + x4(i)
        y = y4(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+37,:) = (/ xy2idx(x4(i),y4(i)), cnt /)
        write(*,*) xy2idx(x4(i),y4(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 5:
    write(*,*) 'ring 5:'
    if (Nring>=5) then
      do i = 1,Natom_ring(5)
        if (y5(i)<=N/3) then
          xmin = N/3.0-y5(i)*0.5
        else
          xmin = N/6.0+(y5(i)-N/3)*0.5
        endif

        x = xmin + x5(i)
        y = y5(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+61,:) = (/ xy2idx(x5(i),y5(i)), cnt /)
        write(*,*) xy2idx(x5(i),y5(i)), cnt
        cnt = cnt+1
      enddo
    endif
    ! ring 6:
    write(*,*) 'ring 6:'
    if (Nring>=6) then
      do i = 1,Natom_ring(6)
        if (y6(i)<=N/3) then
          xmin = N/3.0-y6(i)*0.5
        else
          xmin = N/6.0+(y6(i)-N/3)*0.5
        endif

        x = xmin + x6(i)
        y = y6(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+91,:) = (/ xy2idx(x6(i),y6(i)), cnt /)
        write(*,*) xy2idx(x6(i),y6(i)), cnt
        cnt = cnt+1
      enddo
    endif
    ! ring 7:
    write(*,*) 'ring 7:'
    if (Nring>=7) then
      do i = 1,Natom_ring(7)
        if (y7(i)<=N/3) then
          xmin = N/3.0-y7(i)*0.5
        else
          xmin = N/6.0+(y7(i)-N/3)*0.5
        endif

        x = xmin + x7(i)
        y = y7(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+127,:) = (/ xy2idx(x7(i),y7(i)), cnt /)
        write(*,*) xy2idx(x7(i),y7(i)), cnt
        cnt = cnt+1
      enddo
    endif
    ! ring 8:
    write(*,*) 'ring 8:'
    if (Nring>=8) then
      do i = 1,Natom_ring(8)
        if (y8(i)<=N/3) then
          xmin = N/3.0-y8(i)*0.5
        else
          xmin = N/6.0+(y8(i)-N/3)*0.5
        endif

        x = xmin + x8(i)
        y = y8(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+169,:) = (/ xy2idx(x8(i),y8(i)), cnt /)
        write(*,*) xy2idx(x8(i),y8(i)), cnt
        cnt = cnt+1
      enddo
    endif
    ! ring 9:
    write(*,*) 'ring 9:'
    if (Nring>=9) then
      do i = 1,Natom_ring(9)
        if (y9(i)<=N/3) then
          xmin = N/3.0-y9(i)*0.5
        else
          xmin = N/6.0+(y9(i)-N/3)*0.5
        endif

        x = xmin + x9(i)
        y = y9(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+217,:) = (/ xy2idx(x9(i),y9(i)), cnt /)
        write(*,*) xy2idx(x9(i),y9(i)), cnt
        cnt = cnt+1
      enddo
    endif
! ===========================
! dr=2
elseif (dr==2) then
    ! ring 1:
    write(*,*) 'ring 1:'
    if (Nring>=1) then
      do i = 1,Natom_ring(2)
        if (y2(i)<=N/3) then
          xmin = N/3.0-y2(i)*0.5
        else
          xmin = N/6.0+(y2(i)-N/3)*0.5
        endif

        x = xmin + x2(i)
        y = y2(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+1,:) = (/ xy2idx(x2(i),y2(i)), cnt /)
        write(*,*) xy2idx(x2(i),y2(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 2:
    write(*,*) 'ring 2:'
    if (Nring>=2) then
      do i = 1,Natom_ring(4)
        if (y4(i)<=N/3) then
          xmin = N/3.0-y4(i)*0.5
        else
          xmin = N/6.0+(y4(i)-N/3)*0.5
        endif

        x = xmin + x4(i)
        y = y4(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+13,:) = (/ xy2idx(x4(i),y4(i)), cnt /)
        write(*,*) xy2idx(x4(i),y4(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 3:
    write(*,*) 'ring 3:'
    if (Nring>=3) then
      do i = 1,Natom_ring(6)
        if (y6(i)<=N/3) then
          xmin = N/3.0-y6(i)*0.5
        else
          xmin = N/6.0+(y6(i)-N/3)*0.5
        endif

        x = xmin + x6(i)
        y = y6(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+37,:) = (/ xy2idx(x6(i),y6(i)), cnt /)
        write(*,*) xy2idx(x6(i),y6(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 4:
    write(*,*) 'ring 4:'
    if (Nring>=4) then
      do i = 1,Natom_ring(8)
        if (y8(i)<=N/3) then
          xmin = N/3.0-y8(i)*0.5
        else
          xmin = N/6.0+(y8(i)-N/3)*0.5
        endif

        x = xmin + x8(i)
        y = y8(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+73,:) = (/ xy2idx(x8(i),y8(i)), cnt /)
        write(*,*) xy2idx(x8(i),y8(i)), cnt
        cnt = cnt+1
      enddo
    endif
! ===========================
! dr=3
elseif (dr==3) then
    ! ring 1:
    write(*,*) 'ring 1:'
    if (Nring>=1) then
      do i = 1,Natom_ring(3)
        if (y3(i)<=N/3) then
          xmin = N/3.0-y3(i)*0.5
        else
          xmin = N/6.0+(y3(i)-N/3)*0.5
        endif

        x = xmin + x3(i)
        y = y3(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+1,:) = (/ xy2idx(x3(i),y3(i)), cnt /)
        write(*,*) xy2idx(x3(i),y3(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 2:
    write(*,*) 'ring 2:'
    if (Nring>=2) then
      do i = 1,Natom_ring(6)
        if (y6(i)<=N/3) then
          xmin = N/3.0-y6(i)*0.5
        else
          xmin = N/6.0+(y6(i)-N/3)*0.5
        endif

        x = xmin + x6(i)
        y = y6(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+19,:) = (/ xy2idx(x6(i),y6(i)), cnt /)
        write(*,*) xy2idx(x6(i),y6(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 3:
    write(*,*) 'ring 3:'
    if (Nring>=3) then
      do i = 1,Natom_ring(9)
        if (y9(i)<=N/3) then
          xmin = N/3.0-y9(i)*0.5
        else
          xmin = N/6.0+(y9(i)-N/3)*0.5
        endif

        x = xmin + x9(i)
        y = y9(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+55,:) = (/ xy2idx(x9(i),y9(i)), cnt /)
        write(*,*) xy2idx(x9(i),y9(i)), cnt
        cnt = cnt+1
      enddo
    endif
! ===========================
! dr=4
elseif (dr==4) then
    ! ring 1:
    write(*,*) 'ring 1:'
    if (Nring>=1) then
      do i = 1,Natom_ring(4)
        if (y4(i)<=N/3) then
          xmin = N/3.0-y4(i)*0.5
        else
          xmin = N/6.0+(y4(i)-N/3)*0.5
        endif

        x = xmin + x4(i)
        y = y4(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+1,:) = (/ xy2idx(x4(i),y4(i)), cnt /)
        write(*,*) xy2idx(x4(i),y4(i)), cnt
        cnt = cnt+1
      enddo
    endif

    ! ring 2:
    write(*,*) 'ring 2:'
    if (Nring>=2) then
      do i = 1,Natom_ring(8)
        if (y8(i)<=N/3) then
          xmin = N/3.0-y8(i)*0.5
        else
          xmin = N/6.0+(y8(i)-N/3)*0.5
        endif

        x = xmin + x8(i)
        y = y8(i)*pa
        call writefile(cnt, x, y)

        idx_atom(i+25,:) = (/ xy2idx(x8(i),y8(i)), cnt /)
        write(*,*) xy2idx(x8(i),y8(i)), cnt
        cnt = cnt+1
      enddo
    endif
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
      write(11,'(I1, A1, I1, A2, A3, A1, A3, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, A3, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, A3, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, A3, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","1.0"," ","0.0"&
               ," ","0.0"," ",tt," ",tt," ","0.0"
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
      write(11,'(I1, A1, I1, A2, A3, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A3, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A3, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A3, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A3, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
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
      write(11,'(I1, A1, I1, A2, A4, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<10 .and. b<100) then
      write(11,'(I1, A1, I2, A2, A4, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<100 .and. b<100) then
      write(11,'(I2, A1, I2, A2, A4, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else if (a<100 .and. b<1000) then
      write(11,'(I2, A1, I3, A2, A4, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    else
      write(11,'(I3, A1, I3, A2, A4, A1, F10.7, A1, A3, A2, F5.3, A2, F5.3, A2, A3)') a," ",b," ","-0.5"," ",pa &
               ," ","0.0"," ",tt," ",tt," ","0.0"
    endif
  enddo
end do

! local hopping between Anderson droplet and metal
do i = 1, Nimp
    j = idx_atom(i,1)
    k = idx_atom(i,2)
    if(j<10 .and. k<10) then
      write(11,'(I1, A1, I1, A27)') j," ",k," 0.0 0.0 1.0  Vval Vval 0.0"
    else if (j<10 .and. k<100) then
      write(11,'(I1, A1, I2, A27)') j," ",k," 0.0 0.0 1.0  Vval Vval 0.0"
    else if (j<100 .and. k<100) then
      write(11,'(I2, A1, I2, A27)') j," ",k," 0.0 0.0 1.0  Vval Vval 0.0"
    else if (j<100 .and. k<1000) then
      write(11,'(I2, A1, I3, A27)') j," ",k," 0.0 0.0 1.0  Vval Vval 0.0"
    else
      write(11,'(I3, A1, I3, A27)') j," ",k," 0.0 0.0 1.0  Vval Vval 0.0"
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
do i = 1, Nimp
    k = idx_atom(i,2)
      if (k<10) then
        write(11,'(I1, A1, I1, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A6)') k," ",k," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0,"  Uval"
      else if (k<100) then
        write(11,'(I2, A1, I2, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A6)') k," ",k," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0,"  Uval"
          else
        write(11,'(I3, A1, I3, A2, A3, A2, A3, A2, A3, A2, F4.1, A2, F4.1, A6)') k," ",k," ","0.0"," ","0.0" &
               ," ","0.0"," ",0.0," ",0.0,"  Uval"
      endif
end do

! symmetry:
write(11,'(A5)') "#SYMM"
x = N/2
y = p*pa
write(11,'(A3,F4.1, A2, F10.7,A24)') "c6 ",x, " ", y, " 0.0d0 0.0d0 0.0d0 1.0d0"
write(11,'(A3,F4.1, A2, F10.7,A24)') "c6 ",x, " ", y, " 1.0d0 0.0d0 0.0d0 1.0d0"
write(11,'(A3,F4.1, A2, F10.7,A24)') "d  ",x, " ", y, " 0.0d0 0.0d0 1.0d0 0.0d0"
write(11,'(A3,F4.1, A2, F10.7,A24)') "d  ",x, " ", y, " 1.0d0 0.0d0 1.0d0 0.0d0"

write(11,'(A4)') "#END"
end program geom

! =======================================================
subroutine writefile(cnt, x, y)
  integer :: cnt
  real    :: x, y
  character*3 :: str 

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
end subroutine writefile

