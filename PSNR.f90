! peak signal-to-noise calculation that uses hard coding
# to extract results regarding the coherence of spectra
# calculated using FFTs


implicit none
! define variables
integer :: l=514,i,j
double precision :: t1,t2,t3,i1,i2,i3,m1,m2,m3,t4,i4,m4
real*8,dimension(15) :: p1,p2,p3,p4
real*8 :: av1,av2,av3,avp1,avp2,avp3,qqq,avp4,av4
real*8 :: s1,s2,s3,sigma1,sigma2,sigma3,s4,sigma4
!0_space,2_space,4_space

p1= (/17.6296280855614,0.795741187328543,0.189836783089217,0.02974771000362832, &
	0.171710867434058,2.44246575543044,5.516945558113,   &
	1.29,1.00,1.27,0.76,0.793,1.02,1.0448,1.0244/)
p2= (/25.7943926406712,10.5849292202837,14.3045616137613,8.67239842404742, &
        10.0,29.0,23.,   &
	25.37,13.511,12.7,7.33,7.48,6.65,4.94,5.25/)
p3= (/26.4629106278310,9.10255957011416,11.7343194408946,4.58723342632496, &
        10.0,22.0,29.0,   &
	16.13,9.66,9.85,10.88,4.94,5.04,5.57,5.909/)
p4= (/22.5,1.607,3.18,3.34,10.9,12.0,13.9,16.89,18.15,19.35,2.5,7.96,15.31,6.80,6.65/)


write(*,*) p1(15),p2(15),p3(15)
!AVERAGE OF THE INTENSITIES
avp1=SUM(p1)/15
avp2=SUM(p2)/15
avp3=SUM(p3)/15
avp4=SUM(p4)/15

m1=0.d0
m2=0.d0
m3=0.d0
m4=0.d0
open(1,file="0_space")
open(2,file="2_space")
open(3,file="4_space")
open(4,file="2_space_0.01")

write(6,*) "calculate variance"
do i = 1,l
 read(1,*)t1,i1
 read(2,*)t2,i2
 read(3,*)t3,i3
 read(4,*)t4,i4
 m1=m1+i1 
 m2=m2+i2
 m3=m3+i3
 m4=m4+i4
enddo
write(6,*) 'ok'
close(1);close(2);close(3);close(4)
av1 = (m1/l)!**0.5d0
av2 = (m2/l)!**0.5d0
av3 = (m3/l)!**0.5d0
av4 = (m4/l)
write(6,*) 'ok2'

s1=0.d0
s2=0.d0
s3=0.d0
s4=0.d0



open(1,file="0_space")
open(2,file="2_space")
open(3,file="4_space")
open(4,file="2_space_0.01")

i1=0.d0
i2=0.d0
i3=0.d0
i4=0.d0
write(6,*) 'ok3'
do j = 1,l
 read(1,*)t1,i1
 read(2,*)t2,i2
 read(3,*)t3,i3
 read(4,*)t4,i4
 s1 = s1 + (i1-av1)**2.d0
 s2 = s2 + (i2-av2)**2.d0
 s3 = s3 + (i3-av3)**2.d0
 s4 = s4 + (i4-av4)**2.d0
enddo
write(6,*) 'ok4'
sigma1=s1/(l-1)
sigma2=s2/(l-1)
sigma3=s3/(l-1)
sigma4=s4/(l-1)

write(*,*) "Variance calculated"

qqq = avp1/sigma1
!PSNR calculation
write(6,*) "PSNR"
write(6,*) "uncoupled:", (avp1/sigma1)/qqq
write(6,*) "2-space:",((avp2)/(sigma2))/qqq
write(6,*) "4-space:",((avp3)/(sigma3))/qqq
write(6,*) "2-space_0.01:",((avp4)/(sigma4))/qqq
write(6,*) "   "
write(6,*) "-------------------------------"
write(6,*) "AHI     AOI     sigma"
write(6,*) avp1,"     ",av1,"     ",sigma1
write(6,*) avp2,"     ",av2,"     ",sigma2
write(6,*) avp3,"     ",av3,"     ",sigma3

