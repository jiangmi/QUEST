program geom
! Generate geometry file for slanted interface (single supercell as the whole lattice)

implicit none

! Variable declarations:

real :: x, y, N=6.0, delta=0.5, U=1.0
integer :: i, j, k, a, flag
character*3 :: str  
character*18 :: s1, s2, s3
real*8, allocatable :: r(:)
real*8, allocatable :: ran(:)

open(unit=11,file='SDWgeom',status='replace', action='write')
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

write(11,'(A3, A1, A3, A1, A3)') "0.0"," ","0.0"," ","1.0"
write(11,'(A6)') "#SUPER"
write(11,'(A3)') "1 0"
write(11,'(A3)') "0 1"
! orbitals
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

! =======================================================
write(11,'(A5)') "#SYMM"
write(s1,'(A18)') " 1.0d0 0.0d0 0.0d0"
write(s2,'(A18)') " 0.0d0 1.0d0 0.0d0"
write(s3,'(A18)') " 0.0d0 0.0d0 1.0d0"
do j = 0, N-1
  do i = 0, N-1
   ! if (i<10) then
      write(11,'(A3, F3.1, A2, F3.1, A2, A3, A18)') "d  ",real(i)," ",real(j)," ","0.0",s1
      write(11,'(A3, F3.1, A2, F3.1, A2, A3, A18)') "d  ",real(i)," ",real(j)," ","0.0",s2
   ! endif
  end do
end do
do j = 0, N-1
  do i = 0, N-1
   ! if (i<10) then
      write(11,'(A4, F3.1, A2, F3.1, A2, A3, A18)') "c4  ",real(i)," ",real(j)," ","0.0",s3
   ! endif
  end do
end do

! ======================================================
!write(11,'(A6)') "#PHASE"
!write(11,'(F4.1, A1, A3)') N,"","0.0"
!write(11,'(A3, A1, F4.1)') "0.0"," ",N
!do i = 0, N-1
!  do j = 0, N-1
!    a = i*N+j
!    if (mod(i+j,2)==0) then
!      if (a<10) then
!        write(str,'(I1)') int(a)
!        write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"," ","1.0"
!      else if (a<100) then
!        write(str,'(I2)') int(a)
!        write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"," ","1.0"
!          else
!        write(str,'(I3)') int(a)
!        write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3, A2, A3)') 's'//str," ",real(i)," ",real(j)," ","0.0"," ","1.0"
!      endif
!    else
!      if (a<10) then
!        write(str,'(I1)') int(a)
!        write(11,'(A2, A1, F4.1, A2, F4.1, A2, A3, A2, A4)') 's'//str," ",real(i)," ",real(j)," ","0.0"," ","-1.0"
!      else if (a<100) then
!        write(str,'(I2)') int(a)
!        write(11,'(A3, A1, F4.1, A2, F4.1, A2, A3, A2, A4)') 's'//str," ",real(i)," ",real(j)," ","0.0"," ","-1.0"
!          else
!        write(str,'(I3)') int(a)
!        write(11,'(A4, A1, F4.1, A2, F4.1, A2, A3, A2, A4)') 's'//str," ",real(i)," ",real(j)," ","0.0"," ","-1.0"
!      endif
!    endif
!  end do
!end do
write(11,'(A4)') "#END"
end program geom
