implicit none
integer :: i,j,length
real*8 :: A,pi = 3.14159265359,l,AH,sss
real*8 :: position(8,3),center_x,center_y,center_z,basis2

! Jason
! Take notice of how this script looks identical to fcc_sheet.f90. It's virtually
! identical to that program, only this generates Si structures. Follow the comments
! carefully.



sss = .0 !change to lattice parameter
A = (5.43102511 - sss)/0.52917721067121 ! lattice constant for Si in atomic units
! 5.4310... is the lattice constant for Si in Angstroms. The denominator is the
! Angstrom --> atomic units conversion factor.

 ! lattice constant in atomic units
! the numerator is the lattice constance in angstroms
! the denominator is the Angstroms --> atomic units conversion factor
! The "lattice constant" is the average distance between adjacent atoms in any structure
! it is unique for every material





! These are the basis atoms. Their location determines the overall structure of the
! Si. Atoms may only occupy positions determined by linear combinations of these
! "linearly independent" vectors. The structure of Si is "fcc diamond" -- front-centered
! cubic diamond. All "fcc-d" structures have "basis" atoms at the locations defined below.
! They're annoying, more annoying than ordinary fcc because fcc-d has TWO basis sets.


! build the primitive lattice ----------------------------------

! the first basis (0,0,0)
 position(1,1) = 0
 position(1,2) = 0
 position(1,3) = 0
 position(2,1) = A/2
 position(2,2) = A/2
 position(2,3) = 0
 position(3,1) = A/2
 position(3,2) = 0
 position(3,3) = A/2
 position(4,1) = 0
 position(4,2) = A/2
 position(4,3) = A/2

! the second basis (1/4,1/4,1/4)
	basis2 = A/4
 position(5,1) = basis2
 position(5,2) = basis2
 position(5,3) = basis2
 position(6,1) = basis2+A/2
 position(6,2) = basis2+A/2
 position(6,3) = basis2
 position(7,1) = basis2+A/2
 position(7,2) = basis2
 position(7,3) = basis2+A/2
 position(8,1) = basis2
 position(8,2) = basis2+A/2
 position(8,3) = basis2+A/2
	



!--------------------------------------------------------------------



!This part of the program centers the rod along the cartesian axes.
! determines the length of the rod. The units here are undefined.
length = 1
l = (7./4.+(length-1))*A - A/4.
center_x= l/2
center_y = 7*A/8 !3/4 * A
center_z = 7*A/8 !3/4 * A


open(1,file='atom.inp',action ='write')
open(2,file = 'atom.xyz',action='write')
 !name of file coordinates are written to. RENAME THIS.



! Don't mess with any of this...
!------------------------------------------------------------------------------------

do i = 0,length
 do j = 1,5
  write(1,*)position(j,1)+i*A-center_x,position(j,2)-center_y,position(j,3)-center_z, "14    0"
! these can be used to generate n x n structures. before using, they should be edited
!  if(j<6) then
   write(1,*)position(j,1)+i*A-center_x,position(j,2)-center_y+A,position(j,3)-center_z, "14   0"
   write(1,*)position(j,1)+i*A-center_x,position(j,2)-center_y,position(j,3)-center_z+A, "14    0"
   write(1,*)position(j,1)+i*A-center_x,position(j,2)-center_y+A,position(j,3)-center_z+A, "14   0"
    
  write(2,*)"Si   ",position(j,1)+i*A-center_x,position(j,2)-center_y,position(j,3)-center_z
! these can be used to generate n x n structures. before using, they should be edited
!  if(j<6) then
   write(2,*)"Si   ",position(j,1)+i*A-center_x,position(j,2)-center_y+A,position(j,3)-center_z
   write(2,*)"Si   ",position(j,1)+i*A-center_x,position(j,2)-center_y,position(j,3)-center_z+A
   write(2,*)"Si   ",position(j,1)+i*A-center_x,position(j,2)-center_y+A,position(j,3)-center_z+A
  ! write(1,*)position(j,1)+i*A-center_x,position(j,2)-center_y+A,position(j,3)-center_z
  ! write(1,*)position(j,1)+i*A-center_x,position(j,2)-center_y+AH,position(j,3)-center_z+AH
 ! endif
1 format(3x,3f10.5,i5)
enddo
end do
close(1)
!-----------------------------------------------------------------------------------------

! write out length of rod
write(*,*) "Length of rod in..." 
write(*,*) "...Angstroms", l*0.52917721067121
write(*,*) "...Bohr", l
write(*,*) "...Nanometers", l*0.52917721067121*0.1
write(*,*) "Length of side in..."
write(*,*) "...Bohr", (7./4.)*A 
!write(*,*) "Use this grid spacing", A/16
   !denominator is grid step number in y and z








!use this if nxn structure is used
write(*,*) "Use this grid spacing", (2.*A)/45.
   !denominator is grid step number in y and z

!This will tell you what grid spacing you should use DEPENDING on your
! grid point number in the y and z directions, which should be the same
! if you're making a sheet from a rod.
! Below is the case if Ny=Nz=32 









!write(*,*) 2.*A

!write(*,*) (7.*A/8. + (1./8.)*A)*2.

end
