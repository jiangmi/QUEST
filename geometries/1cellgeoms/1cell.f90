program geom
! Generate geometry file for slanted interface (single supercell as the whole lattice)

implicit none

! Variable declarations:

real :: x, y, N=6.0
integer :: i, j, k, a, flag
character*3 :: str  
real*8, allocatable :: r(:)
real*8, allocatable :: ran(:)

open(unit=11,file='bilayerrandV',status='replace', action='write')
write(11,"(A5)") "#NDIM"
write(11,'(A1)') "1"
write(11,'(A5)') "#PRIM"
write(11,'(F4.1, A1, A3, A1, A3)') N,"","0.0"," ","0.0"
write(11,'(A3, A1, F4.1, A1, A3)') "0.0"," ",N," ","0.0"
write(11,'(A3, A1, A3, A1, A3)') "0.0"," ","0.0"," ","2.0"
write(11,'(A6)') "#SUPER"
write(11,'(A1)') "1"
! orbitals
write(11,'(A4)') "#ORB"
do i = 0, N-1
  do j = 0, N-1
    a = i*N+j
    if (a<10) then
      write(str,'(I1)') int(i*N+j)
      write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"
    else if (a<100) then
      write(str,'(I2)') int(i*N+j)
      write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"
        else
      write(str,'(I3)') int(i*N+j)
      write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"
    endif
  end do
end do
do i = 0, N-1
  do j = 0, N-1
    a = i*N+j+N**2
    if (a<10) then
      write(str,'(I1)') int(a)
      write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","1.0"
    else if (a<100) then
      write(str,'(I2)') int(a)
      write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","1.0"
        else
      write(str,'(I3)') int(a)
      write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","1.0"
    endif
  end do
end do

write(11,'(A30)') "#HAMILT            tup  tdn  U"
! hopping along y direction
do i = 1, 2*N**2-1
  x = mod(i, int(N))
  y = floor(i/N)
  ! hopping across the interface
  if (x>0) then
    if(i<10) then
      write(11,'(I1, A2, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i==10) then
      write(11,'(I1, A2, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i<100) then
      write(11,'(I2, A2, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else if (i==100) then
      write(11,'(I2, A2, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    else 
      write(11,'(I3, A2, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i," ","0.0"," ","1.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
    endif
 ! else 
 !   if(i<10) then
 !     write(11,'(I1, A2, I1, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",int(i-N)," ","0.0"," ","1.0"&
 !              ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
 !   else if (i==10) then
 !     write(11,'(I1, A2, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",int(i-N)," ","0.0"," ","1.0"&
 !              ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
 !   else if (i<100) then
 !     write(11,'(I2, A2, I2, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",int(i-N)," ","0.0"," ","1.0"&
 !              ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
 !   else if (i==100) then
 !     write(11,'(I2, A2, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",int(i-N)," ","0.0"," ","1.0"&
 !              ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
 !   else
 !     write(11,'(I3, A2, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",int(i-N)," ","0.0"," ","1.0"&
 !              ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
 !   endif
  endif
end do
! hopping along x direction
do i = 1, N**2-N
  write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
end do
do i = N**2+1, 2*N**2-N
  write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i+int(N)-1," ","1.0"," ","0.0"&
               ," ","0.0"," ","1.0"," ","1.0"," ","0.0"
end do
! hopping along z direction
! random nonzero t_perp
allocate(r(int(N*N)))
call random_seed()
call random_number(r)
do i = 1, N**2
  write(11,'(I3, A1, I3, A2, A3, A1, A3, A1, A3, A2, F6.3, A2, F6.3, A2, A3)') i-1," ",i+int(N**2)-1," ","0.0"," ","0.0"&
               ," ","1.0"," ",0.5*( 2*r(i)-1 )+1.7," ",0.5*( 2*r(i)-1 )+1.7," ","0.0"
end do

!do i = 0, N-1
!  do j = 0, N-1
!    flag = 0
!    do k = 1, size(r)
!      if (int(i*N+j)==r(k)) then
!        flag = 1
!      endif
!    enddo
!    if (flag == 1) then
!      write(11,'(I3, A2, I3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3)')
!      int(i*N+j)," ",int(i*N+j)," ","0.0"," ","0.0"&
!               ," ","0.0"," ","0.0"," ","0.0"," ","0.0"
!    else
!      write(11,'(I3, A2, I3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3)')
!      int(i*N+j)," ",int(i*N+j)," ","0.0"," ","0.0"&
!               ," ","0.0"," ","0.0"," ","0.0"," ","4.0"
!    endif
!  enddo
!enddo
! ================================================================
! local interaction
! Ordinary interface (along x direction):
!do i = 1, N**2/2
!  write(11,'(I3, A2, I3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1," ","0.0"," ","0.0"&
!               ," ","0.0"," ","0.0"," ","0.0"," ","0.0"
!enddo
do i = 1, 2*N**2
  write(11,'(I3, A2, I3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3)') i-1," ",i-1," ","0.0"," ","0.0"&
               ," ","0.0"," ","0.0"," ","0.0"," ","4.0"
enddo

! =======================================================
! random nonzero U sites in U=0 metal lattice
!allocate(ran(int(N*N*0.5)))
!allocate(r(int(N*N*0.5)))
!call random_seed()
!call random_number(ran)
!r=int(ran*N*N)+1
!do i = 0, N-1
!  do j = 0, N-1
!    flag = 0
!    do k = 1, size(r)
!      if (int(i*N+j)==r(k)) then
!        flag = 1
!      endif
!    enddo
!    if (flag == 1) then
!      write(11,'(I3, A2, I3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3)') int(i*N+j)," ",int(i*N+j)," ","0.0"," ","0.0"&
!               ," ","0.0"," ","0.0"," ","0.0"," ","0.0"
!    else
!      write(11,'(I3, A2, I3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3, A2, A3)') int(i*N+j)," ",int(i*N+j)," ","0.0"," ","0.0"&
!               ," ","0.0"," ","0.0"," ","0.0"," ","4.0"
!    endif
!  enddo
!enddo
!do i=1, size(r)
!  write(*,*) r(i)
!enddo
! =======================================================
write(11,'(A5)') "#SYMM"
write(11,'(A38)') "d  0.0d0 0.0d0 0.5d0 0.0d0 0.0d0 0.5d0"
! ======================================================
write(11,'(A6)') "#PHASE"
write(11,'(A34)') "1  # of cells to specify the phase"
do i = 0, N-1
  do j = 0, N-1
    if (mod(i+j,2)==0) then
      write(11,'(I3, A1, F4.1, A2, F4.1, A2, A3, A2, A4)') int(i*N+j)," ",real(i)," ",real(j)," ","0.0"," ","1.d0"
    else
      write(11,'(I3, A1, F4.1, A2, F4.1, A2, A3, A2, A5)') int(i*N+j)," ",real(i)," ",real(j)," ","0.0"," ","-1.d0"
    endif
  end do
end do
do i = 0, N-1
  do j = 0, N-1
    if (mod(i+j,2)==0) then
      write(11,'(I3, A1, F4.1, A2, F4.1, A2, A3, A2, A5)') int(i*N+j+N**2)," ",real(i)," ",real(j)," ","1.0"," ","-1.d0"
    else
      write(11,'(I3, A1, F4.1, A2, F4.1, A2, A3, A2, A4)') int(i*N+j+N**2)," ",real(i)," ",real(j)," ","1.0"," ","1.d0"
    endif
  end do
end do

write(11,'(A4)') "#END"
end program geom
