MODULE ReimannSilberstein
!
! Maxwell module
!
  Use global_stuff
  Use FD_Operators
  USE conducting_surface

  implicit none
  double precision               :: dx,dt,time_RS
  integer                        :: RS_iter_per_TS=20 ! # of FDTD steps per Taylor time step
  integer                        :: RS_HelmhMODE=1
  integer                        :: laser_shape
  ! to use custom laser
  integer                        :: use_custom_laser=0 
  integer                        :: custom_laser_type=1
  real*8                         :: custom_laser_parm1,custom_laser_parm2,custom_laser_parm3
  real*8                         :: custom_laser_parm4,custom_laser_parm5,custom_laser_parm6
  
  !grid step, same as SE part
  real*8                         :: scale_current=1.d0

  real*8,parameter               :: convert_rF2E=sqrt(2.d0/eps0)
  real*8,parameter               :: convert_E2rF=sqrt(eps0/2.d0)
  real*8,parameter               :: convert_iF2B=sqrt(2*mu0)
  real*8,parameter               :: convert_B2iF=1.d0/convert_iF2B  

!  include 'fftw3.f'
!  integer,parameter      :: dirnF_FFT=FFTW_FORWARD
!  integer,parameter      :: dirnB_FFT=FFTW_BACKWARD    
!  integer*8              :: RS_planF,RS_planB
  
! laser parameters
  real*8,parameter              :: min_omega=1.d-4 ! if omega is less than this, it is single cycle
  double precision              :: E_0,omega,tt,laser_wavelength
  double precision              :: laser_pulse_shift1,laser_pulse_width1,ramp_time
  !
  integer                     :: RS_freq_save_wall
  real*8                      :: Wid_s,C_s
  integer,parameter           :: RS_N_taylor=4
  complex*16                  :: tcRS(0:RS_N_taylor),Fextern(3) ! homogeneous external EM field
  complex*16,allocatable                 ::  Fxyz(:,:),Finduced(:,:),Fdiv(:,:),Fkeep(:,:)
  complex*16,allocatable                 ::  Fnew(:,:),HemF(:,:),curAdd2F(:,:)
  complex*16,allocatable                 ::  deldotF(:)
  real*8,allocatable                     ::  divE(:),V_divE(:),Efield(:,:),Efield_init(:,:)
  real*8,allocatable                     ::  cur(:,:),cur_nuc(:,:),AvecF(:,:),AoldF(:,:),tmp_real3_N(:,:),tmp_realN_3(:,:),tmp2_realN_3(:,:)
  real*8,allocatable                     ::  conductivity_core(:),E_core_max(:)
  real*8                                 ::  F_energy,F_norm,source,source_prev
  real*8                                 :: ramp,Fene_init,Fene_final,F_eneinduced
  complex*16,allocatable                 :: Fcap(:),Fcap_ax(:),Fcap_ay(:),Fcap_az(:),Fcap_hard(:)
  complex*16,allocatable                 :: RS_cap_x(:),RS_cap_y(:),RS_cap_z(:)
  integer                                :: RS_cap_thickness_ngp(3)
  real*8,allocatable                     :: Ef1dx(:) ! 1d E field from laser
  real*8,allocatable                     :: divA(:),density_out_ME(:)


  real*8,parameter :: RS_capP(5)=(/0.00127967d0,0.000486973d0,0.00978732d0,0.000277563d0,10.d0/)
  complex*16,allocatable                  :: Fc3rs(:,:,:),Fc3k(:,:,:)
  real*8                 :: faceidx(6),xpidx(6)
  logical                :: usecurlcur=.TRUE.

  integer,parameter      :: wall_left_file=1301301
  integer,parameter      :: wall_right_file=1301302
  integer,parameter      :: wall_right_fmulti=1301303
  integer,parameter      :: wall_right_multi_dist=5
  
  !
  logical                :: set_initial_fstate=.FALSE.
  real*8                 :: set_initial_f_e0=0.d0
  real*8                 :: set_initial_f_wid=0.d0
  real*8                 :: set_initial_f_shift=0.d0
 
  
  logical                :: use_core_conductivity=.FALSE.
  logical		 :: conducting_shape, mix_conductor_schrodinger
  
  CONTAINS

  include 'custom_laser.f90'
  
subroutine init_RS
implicit none
  integer                     :: i1,i2,i3,num,i,k,j
  real*8                      :: cap_depth_x,cap_depth_y,cap_depth_z,init_F_max,init_F_maxpos,FEt0
  
  if(RS_iter_per_TS/=0.d0) then
    RS_time_step=time_step/dfloat(RS_iter_per_TS)
  else
    RS_time_step=0.d0
  end if
  ! allocate RS stuff
  allocate(Fxyz(N_L_Points_ME,3),Finduced(N_L_Points_ME,3),Fdiv(N_L_Points_ME,3),Fnew(N_L_Points_ME,3),HemF(N_L_Points_ME,3),cur(3,N_L_Points_ME),curAdd2F(N_L_Points_ME,3),Fkeep(N_L_Points_ME,3))
  allocate(cur_nuc(3,N_L_Points_ME))
  allocate(AvecF(3,N_L_Points_ME))
  allocate(Fcap(N_L_Points_ME),Fcap_ax(N_L_Points_ME),Fcap_ay(N_L_Points_ME),Fcap_az(N_L_Points_ME),Fcap_hard(N_L_Points_ME))
  allocate(RS_cap_x(N_L_ME(1)),RS_cap_y(N_L_ME(2)),RS_cap_z(N_L_ME(3)))
  allocate(Fc3rs(N_L_ME(1),N_L_ME(2),N_L_ME(3)),Fc3k(N_L_ME(1),N_L_ME(2),N_L_ME(3)))

  allocate(Ef1dx(N_L_ME(1)))
  allocate(divA(N_L_Points_ME))
  allocate(density_out_ME(N_L_Points_ME))
  allocate(tmp_real3_N(3,N_L_Points_ME),tmp_realN_3(N_L_Points_ME,3),tmp2_realN_3(N_L_Points_ME,3))
  allocate(wf_ME(-N_d:N_L_ME(1)+N_d,-N_d:N_L_ME(2)+N_d,-N_d:N_L_ME(3)+N_d))
  
  allocate(deldotF(N_L_Points_ME),divE(N_L_Points_ME),V_divE(N_L_Points_ME),Efield(3,N_L_Points_ME),Efield_init(3,N_L_Points_ME))
  
  if(time_order_mode==12) allocate(AoldF(3,N_L_Points_ME))
  
  call sleep(1)
  
  Fxyz=0.d0
  Finduced=0.d0
  Fdiv=0.d0
  Fnew=0.d0
  HemF=0.d0
  cur=0.d0
  curAdd2F=0.d0
  Fkeep=0.d0
  cur_nuc=0.d0
  AvecF=0.d0
  divA=0.d0
  Fcap=0.d0
  Fcap_ax=0.d0
  Fcap_ay=0.d0
  Fcap_az=0.d0
  RS_cap_x=0.d0
  RS_cap_y=0.d0
  RS_cap_z=0.d0
  Fc3rs=0.d0
  Fc3k=0.d0
  Ef1dx=0.d0
  tmp_real3_N=0.d0
  tmp_realN_3=0.d0
  
  ! init RS cap
  cap_depth_x=RS_cap_thickness_ngp(1)*grid_step_ME(1)
  if(N_L_ME(1)>2*RS_cap_thickness_ngp(1).AND.(RS_cap_thickness_ngp(1)>0)) then
    call polycap(RS_cap_x,grid_x_ME,N_L_ME(1),cap_depth_x,grid_step_ME(1),RS_capP(1),RS_capP(2),RS_capP(3),RS_capP(4),RS_capP(5))
  end if
  cap_depth_y=RS_cap_thickness_ngp(2)*grid_step_ME(2)
  if(N_L_ME(2)>2*RS_cap_thickness_ngp(2).AND.(RS_cap_thickness_ngp(2)>0)) then
    call polycap(RS_cap_y,grid_y_ME,N_L_ME(2),cap_depth_y,grid_step_ME(2),RS_capP(1),RS_capP(2),RS_capP(3),RS_capP(4),RS_capP(5))
  end if
  
  cap_depth_z=RS_cap_thickness_ngp(3)*grid_step_ME(3)
  if(N_L_ME(3)>2*RS_cap_thickness_ngp(3).AND.(RS_cap_thickness_ngp(3)>0)) then
    call polycap(RS_cap_z,grid_z_ME,N_L_ME(3),cap_depth_z,grid_step_ME(3),RS_capP(1),RS_capP(2),RS_capP(3),RS_capP(4),RS_capP(5))
  end if
  
  write(6,*) "cap depth xyz ",cap_depth_x,cap_depth_y,cap_depth_z
  
  open(386,file="cap_F.dat",form='formatted',status='replace',action='write')
  write(386,*) "#cap_p1 ",RS_capP(1)
  write(386,*) "#cap_p2 ",RS_capP(2)
  write(386,*) "#cap_p3 ",RS_capP(3)
  write(386,*) "#cap_p4 ",RS_capP(4)
  write(386,*) "#extent  ",RS_capP(5)
  do j=1,N_L_Points_ME
    Fcap_ax(j)=RS_cap_x(lattice_ME(1,j))
    Fcap_ay(j)=RS_cap_y(lattice_ME(2,j))
    Fcap_az(j)=RS_cap_z(lattice_ME(3,j))
    if(abs(imag(Fcap_ax(j)+Fcap_ay(j)+Fcap_az(j)))>1.d-6) write(386,'(3f15.3,3es16.6e3)') grid_point_ME(1,j),grid_point_ME(2,j),grid_point_ME(3,j),imag(Fcap_ax(j)+Fcap_ay(j)+Fcap_az(j)),Fcap_hard(j)
  end do
  close(386)  
  
  do i=1,RS_N_taylor
    tcRS(i)=(-zi*RS_time_step)**i
    do k=1,i
      tcRS(i)=tcRS(i)/k
    end do
  end do
  
  open(wall_left_file,file="wall_avg_left.txt",form='formatted',status='replace',action='write')
  write(wall_left_file,'(a)') "# time      EfAVGx           EfAVGy           EfAVGz           Ef1d_z          Enet_z"
  write(wall_left_file,'(a,f12.5)') "#Xpos=",grid_x_ME(RS_cap_thickness_ngp(1)+1)
  open(wall_right_file,file="wall_avg_right.txt",form='formatted',status='replace',action='write')
  write(wall_right_file,'(a)') "# time      EfAVGx           EfAVGy           EfAVGz           Ef1d_z          Enet_z"
  if(N_L_ME(1)-RS_cap_thickness_ngp(1)-1>0) write(wall_left_file,'(a,f12.5)') "#Xpos=",grid_x_ME(N_L_ME(1)-RS_cap_thickness_ngp(1)-1)
  open(wall_right_fmulti,file="wall_avg_multiR.txt",form='formatted',status='replace',action='write')  
  write(wall_right_fmulti,'(a)') "# time      EfAVGx1           EfAVGx2           EfAVGx3           Ef1d_x1          Ef1d_x2          Ef1d_x3"
  i=N_L_ME(1)-RS_cap_thickness_ngp(1)-1

  write(wall_right_fmulti,'(a,3f12.5)') "#Xpos=",grid_x_ME(i-2*wall_right_multi_dist),grid_x_ME(i-wall_right_multi_dist),grid_x_ME(i) 
  
  if(use_core_conductivity==.TRUE.) call init_core_conductivity

  if(conducting_shape==.TRUE.) call RS_assign_conductivity


 
end subroutine init_RS


function Ez_atX(time,x)
implicit none
  real*8   :: time,x
  real*8   :: tf,pf,KWAVE,xc,ef,T0c,WID,tcorr,rfac
  real*8   :: Ez_atX
  
  T0c=laser_pulse_shift1!/c_speed
  WID=laser_pulse_width1
  tcorr=time-x/c_speed
  
  rfac=1.d0
  if(tcorr<ramp_time) rfac=tcorr/ramp_time
  if(tcorr<0.d0) rfac=0.d0
  if(omega<min_omega) then
    tf=1.d0
  else
    tf=sin(omega*tcorr)
  end if
  pf=exp(-((tcorr - T0c)/WID)**2)
  !if(laser_shape==2) pf=sin(pi*(tcorr-T0c)/T0c+pi)
  Ez_atX=E_0*tf*pf*rfac

end function Ez_atX

subroutine gen_1dE_addA(time)
implicit none
  integer :: i,j
  real*8  :: time
  
  ! use_custom_laser==0 default laser profile
  ! use_custom_laser==1 custom laser profile
  ! use_custom_laser==2 3d laser field
  
  if(use_custom_laser==2) then
      ! set tmp_real3_N
      call Exyz_custom(time)
      AvecF=AvecF-tmp_real3_N*time_step
  else
      do i=1,N_L_ME(1)
        if(use_custom_laser==0) then ! default laser profile
          Ef1dx(i)=Ez_atX(time,grid_x_ME(i))
    	else if(use_custom_laser==1) then ! custom profile in 1D
    	  Ef1dx(i)=Ez_custom(time,grid_x_ME(i))
    	end if
      end do
      
      do i=1,N_L_Points_ME
        AvecF(3,i)=AvecF(3,i)-Ef1dx(lattice_ME(1,i))*time_step
      end do
  end if
end subroutine gen_1dE_addA

subroutine polycap(capval,axis,n,capdist,gridspacing,pc1,pc2,pc3,pc4,extent)
implicit none
  integer  :: n,i
  complex*16   :: capval(n),captmp1,captmp2
  real*8   :: axis(n)
  real*8   :: capdist,gridspacing
  real*8   :: pc1,pc2,pc3,pc4,extent
  real*8   :: axmax,axmin
  !
  real*8   :: x,xl,xr,yr,yl,dr1,c1,c2
  real*8   :: cap_parameter=0.3d0

  axmin=axis(1)
  axmax=axis(n)
  c1=axmin+capdist
  c2=axmax-capdist
  write(*,'(a,f15.5,a,f15.5)') "generate Cap from ",c2," to ",axmax
  write(*,'(a,f15.5,a,f15.5)') "generate Cap from ",c1," to ",axmin
  write(*,'(a,4f20.10)') "polynomial coefficients ",pc1,pc2,pc3,pc4

  do i=1,n
    x=axis(i)
    captmp1=0.d0
    captmp2=0.d0
    xl=abs(x-c1)*extent/capdist
    yl=pc1*xl + pc2*xl**2 + pc3*xl**3 + pc4*xl**4
    xr=abs(x-c2)*extent/capdist
    yr=pc1*xr + pc2*xr**2 + pc3*xr**3 + pc4*xr**4
    if(x<=c1) captmp1=-zi*h2m*(2*3.14159265359)**2*yl
    if(x>=c2) captmp2=-zi*h2m*(2*3.14159265359)**2*yr
    capval(i)=capval(i)+captmp1+captmp2
  end do
	
end subroutine polycap

subroutine RS_timestep
implicit none
  integer     :: n,i,j,NL,i2,i3,step
  logical  :: dopropF
  real*8   :: mixc,R1

  
  NL=N_L_Points_ME
  
  ! split operator type ABCs
  Finduced(:,1)=Finduced(:,1)*exp(-zi*(Fcap_ay(:)+Fcap_az(:))*time_step/2.d0)
  Finduced(:,2)=Finduced(:,2)*exp(-zi*(Fcap_ax(:)+Fcap_az(:))*time_step/2.d0)
  Finduced(:,3)=Finduced(:,3)*exp(-zi*(Fcap_ax(:)+Fcap_ay(:))*time_step/2.d0)
  do j=1,3
    Finduced(:,j)=Finduced(:,j)*exp(-zi*Fcap_hard(:)*time_step/2.d0)
  end do
  
 ! conducting sheet

  if(conducting_shape==.TRUE.) then
   if(mix_conductor_schrodinger==.false.) cur = 0.d0
    do i=1,N_L_points_ME
         if(RS_HelmhMODE==1) then
          cur(1,i)=cur(1,i)+conductivity(i)*(real(Finduced(i,1)*convert_rF2E))
          cur(2,i)=cur(2,i)+conductivity(i)*(real(Finduced(i,2)*convert_rF2E))
          cur(3,i)=cur(3,i)+conductivity(i)*((real(Finduced(i,3)*convert_rF2E)+Ef1dx(Lattice_ME(1,i))))

         else
         write(*,*) "conductivity restricted to calculations including RS_HelmhMODE==1"
		stop
         endif
    enddo
   endif
  
  ! core conductivity
  if(use_core_conductivity) call gen_core_current
  
  ! add current
   do i=1,N_L_Points_ME
     do j=1,3
!         --------- BUG FIX  :: old was curAdd2F(i,j)=-cur(j,i)*RS_time_step/eps0
       curAdd2F(i,j)=-(cur(j,i)+cur_nuc(j,i))*RS_time_step/sqrt(2*eps0) ! use different time step
!       curAdd2F(i,j)=-cur(j,i)*RS_time_step/sqrt(2*eps0) ! use different time step
     end do
   end do

   call dofftvector(curAdd2F,planForward)
   call dofftvector(Finduced,planForward)

  ! do Helmholtz decomposition in k-space
  if(RS_HelmhMODE==2) then
    call Helmholtz_decomp_FC(curAdd2F,Fdiv)
    curAdd2F=Fdiv
  end if   
  
  do step=1,RS_iter_per_TS
! apply time propagator for 
    Finduced=Finduced+curAdd2F
    Fnew=Finduced
    do n=1,RS_N_taylor
      call do_Hem_F(Fnew)
      do j=1,3
	Finduced(:,j)=Finduced(:,j)+tcRS(n)*(HemF(:,j))
      end do
      Fnew=HemF
    end do
  end do
  
  ! do Helmholtz decomposition in k-space
  if(RS_HelmhMODE==1) then
    call Helmholtz_decomp_FC(Finduced,Fdiv,deldotF)
    call dofftvector(Fdiv,planBackward)
    if((TD_MODE==4).or.(conducting_shape==.true.)) then
! need divE potential
! deldotF is essentially div(E)=rho/eps_0
! for potential need  F'{e*F{div(E)}/|k**2|}
! solve by solving Poisson equation. Already in K space so just multiply
      do i=1,N_L_Points_ME
        R1=ElectronCharge/kvec_xyz_lenSq_ME(i)
	if(kvec_xyz_lenSq_ME(i)<1d-10) R1=0.d0
	deldotF(i)=deldotF(i)*R1
      end do
      call dofftscalar(deldotF,planBackward)
      V_divE=real(deldotF)*convert_rF2E
    end if
  else if(RS_HelmhMODE==3) then
    call Helmholtz_decomp_FC(Finduced,Fdiv)
    Finduced=Finduced-Fdiv
  end if

  call dofftvector(Finduced,planBackward)
  
! split operator type ABCs
  Finduced(:,1)=Finduced(:,1)*exp(-zi*(Fcap_ay(:)+Fcap_az(:))*time_step/2.d0)
  Finduced(:,2)=Finduced(:,2)*exp(-zi*(Fcap_ax(:)+Fcap_az(:))*time_step/2.d0)
  Finduced(:,3)=Finduced(:,3)*exp(-zi*(Fcap_ax(:)+Fcap_ay(:))*time_step/2.d0)
  do j=1,3
    Finduced(:,j)=Finduced(:,j)*exp(-zi*Fcap_hard(:)*time_step/2.d0)
  end do

! energy density 
  F_eneinduced=0.d0
  do j=1,3
    F_eneinduced=F_eneinduced+sum(real(conjg(Finduced(:,j))*Finduced(:,j)))
  end do
  F_eneinduced=F_eneinduced*grid_volume_ME

! Vector potential without divergence part of E
  do i=1,N_L_Points_ME
    do j=1,3
      if(RS_HelmhMODE==1) then
	AvecF(j,i)=AvecF(j,i)-real(Fxyz(i,j)+Finduced(i,j)-Fdiv(i,j))*time_step*convert_rF2E
      else if((RS_HelmhMODE==2).OR.(RS_HelmhMODE==3)) then
        AvecF(j,i)=AvecF(j,i)-real(Fxyz(i,j)+Finduced(i,j))*time_step*convert_rF2E
      end if
    end do
  end do
  
  call divergence(AvecF,divA) ! calculate the divergence of A Just in case...

end subroutine RS_timestep

subroutine do_Hem_F(Fuse)
implicit none
  complex*16 :: Fuse(N_L_Points_ME,3)
  complex*16 :: sf
  integer    :: i
    
  sf=1.d0
  HemF=(0.d0,0.d0)

  do i=1,N_L_Points_ME
! HemF1= d2 F3 - d3 F2
    HemF(i,1)=HemF(i,1)+zi*kvec_xyz_ME(2,i)*Fuse(i,3)-zi*kvec_xyz_ME(3,i)*Fuse(i,2)	
! HemF2= d3 F1 - d1 F3
    HemF(i,2)=HemF(i,2)+zi*kvec_xyz_ME(3,i)*Fuse(i,1)-zi*kvec_xyz_ME(1,i)*Fuse(i,3)	
! HemF3= d1 f2 - d2 f1
    HemF(i,3)=HemF(i,3)+zi*kvec_xyz_ME(1,i)*Fuse(i,2)-zi*kvec_xyz_ME(2,i)*Fuse(i,1)	
  end do
  HemF=HemF*c_speed
end subroutine do_Hem_F


subroutine Helmholtz_decomp_FC(Vin,Vout,a3d) ! works on fourier transformed vectors
implicit none
  integer  :: i,j
  real*8   :: khat,kmag
  complex*16,intent(in)              :: Vin(N_L_Points_ME,3)
  complex*16,intent(out)             :: Vout(N_L_Points_ME,3)
  complex*16,optional,intent(out)    :: a3d(N_L_Points_ME) ! divergence of vector field
  
  Vout=0.d0
  a3d=0.d0
! F'div=k`(k`*F') where k`=k/|k|
! 1st do divergence
  do j=1,3
    do i=1,N_L_Points_ME
      kmag=sqrt(kvec_xyz_ME(1,i)**2+kvec_xyz_ME(2,i)**2+kvec_xyz_ME(3,i)**2)
      if(kmag/=0.d0) then
	khat=kvec_xyz_ME(j,i)/kmag
      else
	khat=0.d0
      end if
      a3d(i)=a3d(i)+Vin(i,j)*khat
    end do
  end do
  ! do gradient
  do j=1,3
    do i=1,N_L_Points_ME
      kmag=sqrt(kvec_xyz_ME(1,i)**2+kvec_xyz_ME(2,i)**2+kvec_xyz_ME(3,i)**2)
      if(kmag/=0.d0) then
	khat=kvec_xyz_ME(j,i)/kmag
      else
	khat=0.d0
      end if
      Vout(i,j)=a3d(i)*khat
    end do
  end do  
  
end subroutine Helmholtz_decomp_FC

subroutine RS_save_wall(iter)
 ! saves average value of E field at wall (just before cap begins)
implicit none 
  integer :: iter
  integer :: i1,i2,i3,i,j,npts
  real*8  :: EfSUM(3)
  real*8  :: zcomp_multi(3),zcomp_1dmulti(3)
  
  ! left
  EfSUM=0.d0
  npts=0
  i1=RS_cap_thickness_ngp(1)+1
  do i2=RS_cap_thickness_ngp(2)+1,N_L_ME(2)-RS_cap_thickness_ngp(2)
    do i3=RS_cap_thickness_ngp(3)+1,N_L_ME(3)-RS_cap_thickness_ngp(3)
      EfSUM(:)=EfSUM(:)+real(Finduced(lattice_inv_ME(i1,i2,i3),:))*convert_rF2E
      npts=npts+1
    end do
  end do
  EfSUM=EfSUM/dfloat(npts)
  write(wall_left_file,'(f10.3,5es16.7e3)') iter*time_step,EfSUM(1:3),Ef1dx(i1),EfSUM(3)+Ef1dx(i1)
  
  ! right !!!!! (also sets stuff below)
  EfSUM=0.d0
  i1=N_L_ME(1)-RS_cap_thickness_ngp(1)-1
  do i2=RS_cap_thickness_ngp(2)+1,N_L_ME(2)-RS_cap_thickness_ngp(2)
    do i3=RS_cap_thickness_ngp(3)+1,N_L_ME(3)-RS_cap_thickness_ngp(3)
      EfSUM(:)=EfSUM(:)+real(Finduced(lattice_inv_ME(i1,i2,i3),:))*convert_rF2E
    end do
  end do
  EfSUM=EfSUM/dfloat(npts)
  write(wall_right_file,'(f10.3,5es16.7e3)') iter*time_step,EfSUM(1:3),Ef1dx(i1),EfSUM(3)+Ef1dx(i1)
  ! right multiple points !!!!!! SET ABOVE BY EfSUM and i1
  zcomp_1dmulti(3)=Ef1dx(i1)
  zcomp_multi(3)=EfSUM(3) ! closest to the wall
  EfSUM=0.d0
  i1=i1-wall_right_multi_dist ! points closer to center
  do i2=RS_cap_thickness_ngp(2)+1,N_L_ME(2)-RS_cap_thickness_ngp(2)
    do i3=RS_cap_thickness_ngp(3)+1,N_L_ME(3)-RS_cap_thickness_ngp(3)
      EfSUM(3)=EfSUM(3)+real(Finduced(lattice_inv_ME(i1,i2,i3),3))*convert_rF2E
    end do
  end do
  EfSUM=EfSUM/dfloat(npts) 
  zcomp_multi(2)=EfSUM(3)
  zcomp_1dmulti(2)=Ef1dx(i1)
  EfSUM=0.d0
  i1=i1-wall_right_multi_dist ! wall_right_multi_dist points in
  do i2=RS_cap_thickness_ngp(2)+1,N_L_ME(2)-RS_cap_thickness_ngp(2)
    do i3=RS_cap_thickness_ngp(3)+1,N_L_ME(3)-RS_cap_thickness_ngp(3)
      EfSUM(3)=EfSUM(3)+real(Finduced(lattice_inv_ME(i1,i2,i3),3))*convert_rF2E
    end do
  end do
  EfSUM=EfSUM/dfloat(npts) 
  zcomp_multi(1)=EfSUM(3)
  zcomp_1dmulti(1)=Ef1dx(i1)
  
  write(wall_right_fmulti,'(f10.3,6es16.7e3)') iter*time_step,zcomp_multi(1:3),zcomp_1dmulti(1:3)
  
end subroutine RS_save_wall

subroutine init_core_conductivity
implicit none
integer             :: N_core,i,j,k1,k2,k3
real*8              :: dx,dy,dz,r1
real*8,allocatable  :: core_pos(:,:),core_rad(:),core_cond(:),core_max(:)

  allocate (conductivity_core(N_L_Points_ME),E_core_max(N_L_Points_ME))
  write(6,*) "initialize conductivity of for core electrons."
  open(1,file='cond.inp')
  read(1,*)N_core
  allocate(core_pos(3,N_core),core_rad(N_core),core_cond(N_core),core_max(N_core))
  do i=1,N_core
      read(1,*)(core_pos(j,i),j=1,3),core_rad(i),core_cond(i),core_max(i)
  end do
  close(1)
  conductivity_core=0.d0
  E_core_max=0.d0
  do i=1,N_core
    do k1=1,N_L_ME(1)
	  dx=grid_x_ME(k1)-core_pos(1,i)
	  if(abs(dx)>core_rad(i)) CYCLE
	  do k2=1,N_L_ME(2)
	    dy=grid_y_ME(k2)-core_pos(2,i)
	    if(abs(dy)>core_rad(i)) CYCLE
	    do k3=1,N_L_ME(3)
		  dz=grid_z_ME(k3)-core_pos(3,i)
		  if(abs(dz)>core_rad(i)) CYCLE
		  r1=sqrt(dx**2+dy**2+dz**2)
		  if(r1<core_rad(i)) then
		    j=Lattice_inv_ME(k1,k2,k3)
		    if(conductivity_core(j)>1.d-10) then
			  write(6,*) "core conductivity has already been assigned for point"
			  write(6,*) grid_x_ME(k1),grid_y_ME(k2),grid_z_ME(k3)
			  write(6,*) "core conductivity spheres should not overlap"
			  stop
			end if
			conductivity_core(j)=core_cond(i)
			E_core_max(j)=core_max(i)
		  end if
		end do
	  end do
	end do
  end do
  
  call save_real_array_bov_ME(conductivity_core,"core",0,0.d0)
  
end subroutine init_core_conductivity

subroutine gen_core_current 
! use E and conductivity_core at each point to generate current 
! from core electrons
implicit none
! NOTE current is added to cur_nuc
integer :: i,j

  tmp_realN_3=real(Finduced)*convert_rF2E  
   do i=1,N_L_Points_ME
     tmp_realN_3(i,3)=tmp_realN_3(i,3)+Ef1dx(Lattice_ME(1,i))
     do j=1,3
    	if(abs(tmp_realN_3(i,j))>E_core_max(i)) then ! just set to +- E_max
    	  tmp_realN_3(i,j)=E_core_max(i)*tmp_realN_3(i,j)/abs(tmp_realN_3(i,j))
    	end if
		cur_nuc(j,i)=cur_nuc(j,i)+conductivity_core(i)*tmp_realN_3(i,j)
     end do
   end do
end subroutine gen_core_current

END MODULE ReimannSilberstein
