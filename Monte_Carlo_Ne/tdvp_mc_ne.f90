module tdvp_vmc
  use linalg  
  implicit none      
! parameters
  integer, parameter           :: ntherm=10
  integer                      :: N_particle,n_par,N_basis,n_var,n_params,Nt
  double precision,parameter   :: h2m=0.5d0,fd_h=0.001d0,pi=3.1415926d0
  double precision             :: step,lambda,omega,dt
! 
  double precision,allocatable :: r(:,:),hmat(:,:),omat(:,:),eval(:),evec(:,:),params(:),c_basis(:), &
&   mmat(:,:),hvec(:),ekin(:), epot(:), photons(:), vmat(:,:),dh(:,:)

  contains
  
  subroutine init
  implicit none
  integer        :: i,k
  N_particle=10
  n_par=0

  open(1,file='mc.dat')
  read(1,*)N_basis
  n_params=n_basis*n_par
  n_var=n_basis*n_par+n_basis
  allocate(params(n_params),photons(n_basis))
  do i=1,n_basis
!    read(1,*)(params((i-1)*n_par+k),k=1,n_par),photons(i)
    read(1,*)photons(i)
  end do
  write(6,*)params
  read(1,*)step,dt,Nt
  read(1,*)omega
  read(1,*)lambda
  lambda=sqrt(2.d0)*lambda
 
  close(1)
    
  allocate(r(3,n_particle),hmat(n_basis,n_basis),omat(n_basis,n_basis),eval(n_basis), &
&   evec(n_basis,n_basis),c_basis(n_basis),mmat(n_var,n_var),hvec(n_var), &
&   ekin(n_basis), epot(n_basis), vmat(n_basis,n_basis),dh(n_basis,n_basis))


      c_basis=1.0d0
      c_basis(1)=1.d0


  end subroutine init

  subroutine mcwalk(naccept,par1,par2,wfold)
    
    implicit none
    integer                    :: naccept,i,k,itherm
    double precision           :: wfnew, wfold, rn, par1(n_params), par2(n_params), rp(3,n_particle)

    naccept=0    
    do itherm=1,ntherm      
!   attempt to move 
    do i=1,N_particle
      do k=1,3
        call random_number(rn)
        rp(k,i)=r(k,i)+step*(rn-0.5d0)
        end do
      end do
      wfnew=overlap(rp,par1,par2)
!     metropolis:
      call random_number(rn)
      if(wfnew/wfold.ge.rn) then
!     accepted:
        naccept=naccept+1
        r=rp
        wfold=wfnew 
      endif
    end do

  end subroutine mcwalk  

subroutine potential_energy(r,e_p)
! calculate the potential energy
implicit none
  real*8             :: r(3,N_particle)
  integer            :: i,j
  real*8             :: e_p,rr
  real*8,parameter   :: Z=10.d0

  e_p=0.d0
  do i=1,N_particle
    rr=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
    e_p=e_p-Z/rr
  end do
  do i=1,N_particle
    do j=i+1,N_particle
      rr=sqrt((r(1,i)-r(1,j))**2+(r(2,i)-r(2,j))**2+(r(3,i)-r(3,j))**2)
      e_p=e_p+1.d0/rr
    end do
  end do
  
end subroutine potential_energy

  subroutine elocal(par1,par2,energy,ov)
    
    implicit none
    double precision   :: par1(n_params),par2(n_params),energy,ov
    integer            :: i,j,k,n,n1,n2
    double precision   :: wf1(n_basis), wf2(n_basis), wfm(n_basis), wfp(n_basis)
    double precision   :: rp(3,N_particle), rm(3,N_particle)
    double precision   :: vp, rr, ke, pe, cc, d, dd

    call wave_fun(r,par1,wf1)
    call wave_fun(r,par2,wf2)

!   kinetic energy
    rp=r
    rm=r
    ekin=0.d0
    do i=1,N_particle
      do k=1,3
        rm(k,i)=r(k,i)-fd_h
        rp(k,i)=r(k,i)+fd_h
        call wave_fun(rp,par2,wfp)
        call wave_fun(rm,par2,wfm)
        ekin=ekin-h2m*(wfp+wfm-2.d0*wf2)/fd_h**2
        rp(k,i)=r(k,i)
        rm(k,i)=r(k,i)
      end do
    end do

    call potential_energy(r,vp)

    d=0.d0
    do i=1,n_particle
!      d=d+lambda*r(3,i)
      d=d+lambda*r(1,i)+lambda*r(2,i)
    end do
    vp=vp+0.5d0*d**2
    
    vmat=0.d0
    do i=1,n_basis
      do j=1,n_basis
        n1=photons(i)
        n2=photons(j)
        if(n1==n2) then
          vmat(i,j)=vp+omega*n1
        endif
        if(n1==n2+1) then
          vmat(i,j)=sqrt(omega/2.d0)*d*sqrt(1.d0*(n1))
        endif
        if(n1==n2-1) then
          vmat(i,j)=sqrt(omega/2.d0)*d*sqrt(1.d0*(n1+1))
        endif
      end do
    end do
 
    epot=matmul(vmat,wf2)

    ke=0.d0
    pe=0.d0
    ov=0.d0
    do i=1,n_basis
      do j=1,n_basis
        dd=0.d0
        if(photons(i)==photons(j)) dd=1.d0
        cc=c_basis(i)*c_basis(j)
        ov=ov+wf1(i)*wf2(j)*cc*dd
        ke=ke+wf1(i)*ekin(j)*cc*dd
        pe=pe+wf1(i)*vmat(i,j)*wf2(j)*cc
      end do
    end do
    energy=ke+pe

  end subroutine elocal

  function overlap(p,par1,par2)

    implicit none
    double precision :: p(3,N_particle),par1(n_params),par2(n_params),overlap
    double precision :: wf1(n_basis),wf2(n_basis),dd
    integer          :: i,j

    call wave_fun(p,par1,wf1)
    call wave_fun(p,par2,wf2)

    overlap=0.d0
    do i=1,n_basis
      do j=1,n_basis
        dd=0.d0
        if(photons(i)==photons(j)) dd=1.d0
        overlap=overlap+wf1(i)*wf2(j)*c_basis(i)*c_basis(j)*dd
      end do      
    end do
          
  end function overlap 

    
  subroutine wave_fun(p,par,wf)
    implicit none
    double precision :: p(3,n_particle), par(n_params), wf(n_basis)

    call wf_Ne(p,par,wf)

  end subroutine wave_fun



  subroutine wf_Ne(p,par,wf)
!     wave function
      implicit none

      double precision :: p(3,n_particle), par(n_params), wf(n_basis)
      real*8             :: fact(0:4),rr(10),phi_up(5,5),phi_down(5,5),Cbs(10,5),xi(10),d_up,d_down
      integer            :: i,j,nb(10),lb(10),mb(10),spin(10)
      real*8             :: a,b,rij,r2x,r2y,r2z,jf,wave_fun,w,d
      integer            :: k,ij,Nbs,f(0:4),kk




      
    d=0.d0
    do i=1,n_particle
!      d=d+lambda*p(3,i)
      d=d+lambda*r(1,i)+lambda*r(2,i)
    end do
    wf=0.d0

    do kk=1,n_basis


Nbs=10

f(0)=1
f(1)=1
f(2)=2
f(3)=6
f(4)=24

cbs=0.d0
cbs(1,1)=1.
cbs(2,1)=.6281849
cbs(3,1)=-.00350573
cbs(2,2)=-.54502867
cbs(3,2)=.58266109
cbs(4,2)=1.
cbs(5,3)=.43026168
cbs(6,3)=1.
cbs(7,4)=.43026168
cbs(8,4)=1.
cbs(9,5)=.43026168
cbs(10,5)=1.



xi(1)=10.69407156
xi(2)=8.41060152
xi(3)=4.36838035
xi(4)=2.82036562
xi(5)=5.45912448
xi(6)=2.41353954
xi(7)=5.45912448
xi(8)=2.41353954
xi(9)=5.45912448
xi(10)=2.41353954

nb(1)=1
lb(1)=0
mb(1)=0
nb(2)=1
lb(2)=0
mb(2)=0
nb(3)=2
lb(3)=0
mb(3)=0
nb(4)=2
lb(4)=0
mb(4)=0
nb(5)=2
lb(5)=1
mb(5)=-1
nb(6)=2
lb(6)=1
mb(6)=-1
nb(7)=2
lb(7)=1
mb(7)=0
nb(8)=2
lb(8)=1
mb(8)=0
nb(9)=2
lb(9)=1
mb(9)=1
nb(10)=2
lb(10)=1
mb(10)=1

!     wave function
      
      wave_fun=0.d0
      jf=1.d0
      b=2.454995
      spin(1)=1
      spin(2)=1
      spin(3)=1
      spin(4)=1
      spin(5)=1
      spin(6)=-1
      spin(7)=-1
      spin(8)=-1
      spin(9)=-1
      spin(10)=-1
      ij=0      
      do i=1,N_particle-1      
       do j=i+1,N_particle
         ij=ij+1
         r2x=(p(1,i)-p(1,j))**2
         r2y=(p(2,i)-p(2,j))**2
         r2z=(p(3,i)-p(3,j))**2
         rij=sqrt(r2x+r2y+r2z)
         if(spin(i).eq.spin(j)) then
           a=0.25d0
         else
           a=0.5d0
         endif
         jf=jf*exp(a*rij/(1.d0+b*rij))
       end do
     end do


do i=1,10
  rr(i)=sqrt(p(1,i)**2+p(2,i)**2+p(3,i)**2)
end do

do k=1,10
  do i=1,5
    w=0.d0
    do j=1,Nbs
      w=w+Cbs(j,i)*sqrt((2*xi(j))**(2*nb(j)+1)/f(2*nb(j)))*rr(k)**(nb(j)-1)*exp(-xi(j)*rr(k))*ylm(p(:,k),lb(j),mb(j))
    end do
    if(k.le.5) then
      phi_up(k,i)=w
    else
      phi_down(k-5,i)=w
    endif
  end do
end do
call det(phi_up,5,d_up)
call det(phi_down,5,d_down)
wf(kk)=jf*d_up*d_down*d**photons(kk)

    end do
    

  end subroutine wf_Ne


function ylm(rv,l,m)
  real*8  :: rv(3),ylm,rr
  integer :: l,m
  rr=sqrt(rv(1)**2+rv(2)**2+rv(3)**2)
  if(l.eq.0.and.m.eq.0)  ylm=1.d0
  if(l.eq.1.and.m.eq.-1) ylm=rv(1)/rr
  if(l.eq.1.and.m.eq.0)  ylm=rv(3)/rr
  if(l.eq.1.and.m.eq.1)  ylm=rv(2)/rr
end function ylm


  
   
  subroutine mcstep(par1,par2,nmoves)
    implicit none
    double precision   :: par1(n_params),par2(n_params)
    integer            :: naccept,imoves,itherm,iaccept,i,j,k,kk,ki,kj,kki,kkj,nmoves
    double precision   :: sume, sume2, wfold, e, o, energy, e2, rn, dx, dd, d1,d2
    double precision   :: p1(n_params), p2(n_params)
    double precision   :: wf1(n_basis), wf2(n_basis), wf(n_basis)
    double precision   :: d_wf1(n_basis*n_par), d_wf2(n_basis*n_par) 

    naccept=0
    sume=0.d0
    sume2=0.d0
    wfold=1.d-50
    dx=0.001d0

    mmat=0.d0
    hvec=0.d0
    hmat=0.d0
    dh=0.d0
    do i=1,N_particle
      do k=1,3
        call random_number(rn)
        r(k,i)=step*(rn-0.5d0)
      end do
    end do

!   MC moves  

    do imoves=1,nmoves
      call mcwalk(iaccept,par1,par2,wfold)
      naccept=naccept+iaccept
!     calculate the energy for the new position
      call elocal(par1,par2,e,o)
      energy=e/o
      sume=sume+energy
      sume2=sume2+energy**2      
!     calculate the M matrix 
      call wave_fun(r,par1,wf1)
      call wave_fun(r,par2,wf2)
!      call d_wave_fun(r,par1,d_wf1)
!      call d_wave_fun(r,par2,d_wf2)
      k=0
      do j=1,n_basis
        do i=1,n_par
          k=k+1
          p1=par1
          p1(k)=par1(k)+dx
          call wave_fun(r,p1,wf)
          d_wf1(k)=(wf(j)-wf1(j))/dx
          p2=par2
          p2(k)=par2(k)+dx
          call wave_fun(r,p2,wf)
          d_wf2(k)=(wf(j)-wf2(j))/dx
        end do
      end do
      
      do i=1,n_basis
        do j=1,n_basis
          dd=0.d0
          if(photons(i)==photons(j)) dd=1.d0
          mmat(i,j)=mmat(i,j)+wf1(i)*wf2(j)/o*dd
          hvec(i)=hvec(i)+wf1(i)*c_basis(j)*(ekin(j)*dd+vmat(i,j)*wf2(j))/o
          hmat(i,j)=hmat(i,j)+wf1(i)*(ekin(j)*dd+vmat(i,j)*wf2(j))/o
        end do
      end do
      do i=1,n_basis
        do j=1,n_basis
          dd=0.d0
          if(photons(i)==photons(j)) dd=1.d0
          do k=1,n_par
            kk=(j-1)*n_par+k
            mmat(i,n_basis+kk)=mmat(i,n_basis+kk)+wf1(i)*d_wf2(kk)/o*c_basis(j)*dd
            mmat(n_basis+kk,i)=mmat(i,n_basis+kk)
          end do
        end do
      end do
      do i=1,n_basis
        do j=1,n_basis
          dd=0.d0
          if(photons(i)==photons(j)) dd=1.d0
          do ki=1,n_par
            kki=(i-1)*n_par+ki
            dh(i,j)=dh(i,j)+d_wf1(kki)*(ekin(j)*dd+vmat(i,j)*wf2(j))/o
            hvec(n_basis+kki)=hvec(n_basis+kki)+d_wf1(kki)*c_basis(i)*c_basis(j)*(ekin(j)*dd+vmat(i,j)*wf2(j))/o
            do kj=1,n_par
              kkj=(j-1)*n_par+kj
              mmat(n_basis+kki,n_basis+kkj)=mmat(n_basis+kki,n_basis+kkj)+d_wf1(kki)*d_wf2(kkj)/o*c_basis(i)*c_basis(j)*dd
            end do
          end do
        end do
      end do
      
      end do
      e=sume/nmoves
      e2=sume2/nmoves      
      write(6,*)'>>>>>>>>>>>>>  imoves:'      
      write(6,*)'energy',e
      write(6,*)'error ',sqrt((e2-e*e)/imoves)      
      write(6,*)'accept',(100.d0*naccept)/(imoves*ntherm)      
      mmat=mmat/nmoves
      hmat=hmat/nmoves
      hvec=hvec/nmoves
      dh=dh/nmoves
      write(16,*)'mmat'
      do i=1,n_var
        write(16,'(6d16.8)')(mmat(j,i),j=1,n_var)
      end do
      write(16,*)'hvec',hvec     

!      stop
      
    end subroutine mcstep
    

end module tdvp_vmc
  

PROGRAM MC
  USE tdvp_vmc  
  use linalg  
  implicit none
 
  integer                                  :: naccept,idm,k,i,imoves,itherm,iaccept,j,n1,n2
  real*8                                   :: sume,sume2,wfnew,wfold,z,e,energy,e2,rn
  real*8, allocatable                      :: update(:),mi(:,:),fot(:)


  call init
  allocate(update(n_var),mi(n_var,n_var),fot(0:10))

      naccept=0
      sume=0.d0
      sume2=0.d0
      wfnew=0.d0
      wfold=1.d-50
      do i=1,N_particle
        do k=1,3
          call random_number(rn)
          r(k,i)=step*(rn-0.5d0)
        end do
      end do
!

      
      do k=1,Nt

      call mcstep(params,params,250000)
      do i=1,n_basis
        do j=1,n_basis
          omat(i,j)=mmat(i,j)
        end do
      end do
 
      write(98,*)hmat
      write(98,*)omat
      

      
      call diag1(hmat,omat,n_basis,eval,evec)
      write(6,*)'energy diag',eval(1)
      c_basis(:)=evec(:,1)     
      call inv(mmat,n_var,mi)

      update=matmul(mi,hvec)
!      do i=1,n_basis
!        c_basis(i)=c_basis(i)-update(i)*dt
!      end do
      

      do i=1,n_basis*n_par
        params(i)=params(i)-update(n_basis+i)*dt
      end do
!      write(6,*)'basis',c_basis
      write(6,*)'params'
      write(6,*)params
      do i=1,n_basis
        do j=1,n_basis
          omat(i,j)=mmat(i,j)
        end do
      end do

      fot=0.d0
      do i=1,n_basis
        do j=1,n_basis
          n1=photons(i)
          n2=photons(j)
          if(n1==n2) fot(n1)=fot(n1)+omat(i,j)*evec(i,1)*evec(j,1)
        end do
      end do
      
      write(6,*)'photons',(fot(i),i=0,2)

      end do

!     final check
      
      call mcstep(params,params,1000000)

      call diag1(hmat,omat,n_basis,eval,evec)
      write(6,*)'energy',eval(1)
      do i=1,n_basis
        do j=1,n_basis
          omat(i,j)=mmat(i,j)
        end do
      end do

      fot=0.d0
      do i=1,n_basis
        do j=1,n_basis
          n1=photons(i)
          n2=photons(j)
          if(n1==n2) fot(n1)=fot(n1)+omat(i,j)*evec(i,1)*evec(j,1)
        end do
      end do
      
      write(6,*)'photons',(fot(i),i=0,2)
      
      end


