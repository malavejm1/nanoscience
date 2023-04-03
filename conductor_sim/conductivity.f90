MODULE conducting_surface
 USE global_stuff
 implicit none

 ! the main variable
  real*8,allocatable  :: conductivity(:)

contains

subroutine RS_assign_conductivity
 implicit none
  integer :: i,j,k,l,geo
  real*8  :: R,R1,R2,R3,D,C
  real*8  :: mask_width
  real*8  :: ii,jj,kk,f1,f2,f3
  real*8  :: shiftx,shifty,shiftz
  real*8  :: sx,sy,sz
  real*8  :: dR1,dR2,dR3,decR1,decR2,decR3



  C=0.d0 ! conductivity
  R=0.d0  !
  mask_width = 0.d0
  allocate(conductivity(N_L_points_ME))
  conductivity=0.d0
  R1 = 0.d0
  R2 = 0.d0
  R3 = 0.d0

  write(6,*) "----- CONDUCTING SURFACE PARAMETERS -----"  

  open(199,file="conductivity.inp",form='formatted')
  read(199,*) geo
        if(geo==1) write(6,*) "geometry: Sphere"
        if(geo==2) write(6,*) "geometry: Sheet"
  read(199,*) R1,R2,R3
        if(geo==1) write(6,*) "R1 only: Radius =", R1
 	if(geo==2) write(6,*) "R1,R2,R3: sidex,sidey,sidez =",R1,R2,R3       
  read(199,*) mask_width
        if(geo==2) write(6,*) "Mask width (only for sheet) =", mask_width
  read(199,*) C
	write(6,*) "Conducitivty (between 0 and 1) =", C
  read(199,*)shiftx,shifty,shiftz
	write(6,*) "Surface center with respect to grid origin: (",shiftx,shifty,shiftz,")"
        if(shiftx.ne.0.or.shifty.ne.0.or.shiftz.ne.0) then
         write(6,*) "--- CAUTION ----"
         write(6,*) "Translation parameters of conducting surface are nonzero."
	 write(6,*) "Imagery slices of displaced surface may improperly represent conductivity density."
	 write(6,*) "Before shifting, visualization of conducting surface with vanishing translation values is advised."  
	endif
  close(199)  





	! CONDUCTING SPHERE
if(geo==1) then
   
   R = R1
   sx = 0.; sy = 0.; sz = 0.
   open(199,file="conductivity.out",form='formatted',status='replace',action='write')
   do i=1,N_L_points_ME

     sx = grid_point_ME(1,i) - shiftx
     sy = grid_point_ME(2,i) - shifty
     sz = grid_point_ME(3,i) - shiftz

     D = sqrt(sx**2 + sy**2 + sz**2)
!     D=sqrt(grid_point_ME(1,i)**2+grid_point_ME(2,i)**2+grid_point_ME(3,i)**2)
        if(D<R) then
          conductivity(i)=C
        else
      conductivity(i)=C*exp(-((D-R)/0.5d0)**2) ! smooth out edge (not discontinuous as we traverse interfaces)
      end if
        if(conductivity(i)>1.d-5) write(199,'(3es10.3e1,1es8.1e1)') grid_point_ME(1,i),grid_point_ME(2,i),grid_point_ME(3,i),conductivity(i)
     end do
     close(199)
     call save_array_text_ME(conductivity,"cond",0,1,1) ! save

endif


if(geo==2) then

open(199,file="conductivity.out",form='formatted',status='replace',action='write')

   mask_width = nint(mask_width)
   


! Ensuring steep density decay off rectangular conductors requires a mask function frequency commensurate
! with the conductor's spatial dimensions. The formula determining frequency compatbility is:

!	(side/2)/(2*(non-zero integer)) = width in which density decays steeply off conductor edges 

! The code-readable statement rounding user input to the nearest integer is justified by this compatiblity
! equation. Ask Justin about its derivation.

    dR1 = .5*R1/(2.*mask_width)
    dR2 = .5*R2/(2.*mask_width)
    dR3 = .5*R3/(2.*mask_width)

    decR1 = dR1+.5*R1
    decR2 = dR2+.5*R2
    decR3 = dR3+.5*R3

    f1 = 0; f2 =0; f3 =0
    sx = 0; sy =0; sz =0

    do i=1,N_L_ME(1),1
     ii = grid_x_ME(i)
     do j=1,N_L_ME(2),1
      jj = grid_y_ME(j)
      do k=1,N_L_ME(3),1
       kk = grid_z_ME(k)

       l = Lattice_inv_ME(i,j,k)


       sx = ii - shiftx
       if(abs(sx).le.(R1/2)) then
        f1 = 1.
       elseif(abs(sx).gt.(R1/2.).and.abs(sx).le.decR1) then
        f1 = 1.- sin((pi/2.)*(abs(sx-.5*R1)/abs(decR1-.5*R1)))**2.
       else
        f1 = 0.
       endif

       sy = jj - shifty
       if(abs(sy).le.(R2/2)) then
        f2 = 1.
       elseif(abs(sy).gt.(R2/2).and.abs(sy).le.decR2) then
        f2 = 1. - sin((pi/2.)*(abs(sy-.5*R2)/abs(decR2-.5*R2)))**2.
       else
        f2=0.
       endif

       sz = kk - shiftz
       if(abs(sz).le.(R3/2)) then
        f3 = 1.
       elseif(abs(sz).gt.(R3/2.).and.abs(sz).le.decR3) then
        f3 = 1. - sin((pi/2.)*(abs(sz-.5*R3)/abs(decR3-.5*R3)))**2.
       else
        f3 = 0.
       endif
       conductivity(l) = C*f1*f2*f3

        if(conductivity(l)>1.d-5) write(199,'(3es10.3e1,1es8.1e1)') grid_point_ME(1,l),grid_point_ME(2,l),grid_point_ME(3,l),conductivity(l)

     end do
     end do
     end do
     close(199)
     call save_array_text_ME(conductivity,"cond",0,1,1) ! save

endif

end subroutine
end
