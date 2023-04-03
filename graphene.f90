implicit none
real*8 :: position(3,4),a1,a2,z1,z2,a0,dx_min,dx_max,dx,dy,diff,min_diff,min_diff_valx,min_diff_valy
real*8 :: best_x_size,best_y_size,use_grid_space
integer :: i1,i2,n1,n2,ii,i,j,min_diff_vali,min_diff_valj

a0=0.529177d0

! n1 is the number of units in the y direction
! n2 is the number of units in the x direction
n1=2
n2=2

ii=0
a1=2.46d0/a0
a2=4.26d0/a0

position(1,1)=0.d0
position(2,1)=0.d0
position(3,1)=0.d0
position(1,2)=0.d0
position(2,2)=2.84d0/a0
position(3,2)=0.d0
position(1,3)=1.23d0/a0
position(2,3)=2.13d0/a0
position(3,3)=0.d0
position(1,4)=1.23d0/a0
position(2,4)=0.71d0/a0
position(3,4)=0.d0

write(1,*)n1*n2*4,n1*n2*4

z1=-n1*a1/2.d0+a1/6.d0
z2=-n2*a2/2.d0+a2/6.d0

do i1=1,n1
  do i2=1,n2
    do i=1,4
      write(1,100)position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i),6,ii
!      write(1,100)position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i),6,ii
!      write(2,*)'C',position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i)
    end do
  end do
end do


100 format(3x,3f10.5,2i5)

! bx=nx*dx
! by=ny*dx
! sx=nu_x*a1
! sy=nu_y*a2
! 
! nu_x*a1=nx*dx
! nu_y*a2=ny*dx
!
!  nu_x * ny     a2
! ----------- = ----
!  nx * nu_y     a1
!
!  dx=nu_x*a1/nx
!  dy=nu_y*a2/ny
!   we want dx=dy . Find case where abs(dx-dy is the smallest) 

write(6,*)' Let us find a size and number of grid points that works well'
!write(6,*)'smallest size of grid points'
!read(5,*)dx_min
!write(6,*)'largest size of grid points'
!read(5,*)dx_max
dx_min=0.4
dx_max=0.55

min_diff=9.d9
min_diff_valx=0.0
min_diff_valy=0.0
min_diff_vali=0
min_diff_valj=0
do i=1,999999
  dx=n1*a1/i
  if((dx>=dx_min).AND.(dx<=dx_max)) then
    do j=1,999999
      dy=n2*a2/j
	  diff=abs(dx-dy)
	  if(diff<min_diff) then
		min_diff_valx=dx
        min_diff_valy=dy
		min_diff_vali=i
        min_diff_valj=j
		min_diff=diff
	  end if
	end do
  end if
end do

use_grid_space=(min_diff_valx+min_diff_valy)/2.d0

write(6,*)' Best match'
write(6,*)' dx=',min_diff_valx
write(6,*)' dy=',min_diff_valy
write(6,*)' with a difference = ',min_diff
write(6,*)' use grid_space=',use_grid_space
write(6,*)' nx = ',min_diff_vali
write(6,*)' ny = ',min_diff_valj

!z1=n1*a1
!write(6,*)z1
!write(6,*)'# of grid points?'
!read(5,*)ii
!write(6,*)'dx',z1/ii
!z2=n2*a2
!write(6,*)z2
!write(6,*)'# of grid points?'
!read(5,*)ii
!write(6,*)'dx',z2/ii

best_x_size=min_diff_vali*use_grid_space
best_y_size=min_diff_valj*use_grid_space
write(6,*) "***adjusting lattice***"
write(6,*) "old lattice height=",a1
a1=best_x_size/n1
write(6,*) "new lattice height=",a1
write(6,*) "old lattice width=",a2
a2=best_y_size/n2
write(6,*) "new lattice width=",a2

position(1,1)=0.d0
position(2,1)=0.d0
position(3,1)=0.d0
position(1,2)=0.d0
position(2,2)=2.84d0/a0
position(3,2)=0.d0
position(1,3)=1.23d0/a0
position(2,3)=2.13d0/a0
position(3,3)=0.d0
position(1,4)=1.23d0/a0
position(2,4)=0.71d0/a0
position(3,4)=0.d0

write(2,*)n1*n2*4,n1*n2*4
open(199,file=trim(adjustl("atom.xyz")),form='formatted',status='replace',action='write')
write(199,*)n1*n2*4
write(199,*) ""
open(299,file=trim(adjustl("atom_periodic.xyz")),form='formatted',status='replace',action='write')
write(299,*)n1*n2*4*5
write(299,*) ""

z1=-n1*a1/2.d0+a1/6.d0
z2=-n2*a2/2.d0+a2/6.d0

do i1=1,n1
  do i2=1,n2
    do i=1,4
      write(2,100)position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i),6,ii
	  write(199,'(a,3F14.6)') "C ",position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i)
	  write(299,'(a,3F14.6)') "C ",position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i)
	  write(299,'(a,3F14.6)') "C ",position(3,i),z1+(i1-1)*a1+position(1,i)+best_x_size,z2+(i2-1)*a2+position(2,i)
	  write(299,'(a,3F14.6)') "C ",position(3,i),z1+(i1-1)*a1+position(1,i)-best_x_size,z2+(i2-1)*a2+position(2,i)
	  write(299,'(a,3F14.6)') "C ",position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i)+best_y_size
	  write(299,'(a,3F14.6)') "C ",position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i)-best_y_size
!      write(1,100)position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i),6,ii
!      write(2,*)'C',position(3,i),z1+(i1-1)*a1+position(1,i),z2+(i2-1)*a2+position(2,i)
    end do
  end do
end do

end
