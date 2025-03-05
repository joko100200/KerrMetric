program main
implicit none

integer :: error,T,i,ttemp,cut,temptemp
real*8,dimension(20) :: connection
real*8,dimension(2) :: r,theta,phi,vr,vtheta,vphi
real*8 :: ar, atheta, aphi, vtemp, a, E, delta, phiterm
real*8 :: dt,r0,theta0,phi0,v01,v02,v03,M,boxsize,J, tc
namelist /inputs/ dt,r0,theta0,phi0,v01,v02,v03,M,J
namelist /simspecs/ boxsize,T,cut

open(10,file='propagator/parameter.inp',status='old',&
        action = 'read', position = 'rewind', iostat = error)
if(error .ne. 0) stop
read(10,inputs)
read(10,simspecs)
close(10)

connection(:) = 0.d0 !Initializing connection array to 0

!------Variable Explinations------!
!T Total integration time
!cut This integer cuts the number of entries written to rthetaphi.out from T to int(T/cut)
!boxsize is the kill barrier for the if r greater than, integrator stops 
!M is mass of black hole making event horizon at 2*M (meters). The mass of the sun is 1480 
!ar,atheta,aphi accelerations in coord directions
!v01,v02,v03 initial velocity in r theta and phi respecivly
!ttemp is a temp variable for the cutting procedure
!connection array stores connections look at subroutine at bottom for further info
!a is the angular momentum per unit mass
!E and delta are length scales
!J is angular momentum
!tc is timecorrection since we are integrating over coordinate time and not proper time

a = J/M

!Initial conditions into arrays
r(1) = r0
theta(1) = theta0
phi(1) = phi0

!Radius check (shouldn't be in the blackhole)
if(r(1)<=2.d0*M) then
write(*,*) "Initial radius inside black hole"
stop
endif

!Initial conditions into arrays
vr(1) = v01
vtheta(1) = v02
vphi(1) = v03

!Velocity Checking should be less than c which is 1
E = r(1)*r(1)+a*a*cos(theta(1))*cos(theta(1))
delta = r(1)*r(1)-2.d0*M*r(1)+a*a
phiterm = (sin(theta(1))*sin(theta(1)))*(r(1)*r(1)+a*a+((2.d0*M*r(1)*a*a)/E)*sin(theta(1))*sin(theta(1))) !g33
vtemp = sqrt(E/delta*vr(1)*vr(1)+E*vtheta(1)*vtheta(1)+phiterm*vphi(1)*vphi(1))
write(*,*) "Initial velocity is (units of c)"
write(*,*) vtemp
if(vtemp >= 1.d0) then
write(*,*) "initial velocity faster than light lmaoo"
stop
endif



open(11,file='rthetaphi.out',status='unknown',action='write',position='rewind')
write(11,*) r(1),theta(1),phi(1) !Writing initial coords to file

!------------Integration loop------------!
do i = 1, T, 1
call connections(connection,r(1),theta(1),M,a) !Calling connections at point in space

!---Time Correction---!
!Since this integrator is propegating using coordinate time it needs an additional 
tc = 2.d0*(connection(2)*vr(1)+connection(3)*vtheta(1)+connection(9)*vr(1)*vtheta(1)+connection(20)*vtheta(1)*vphi(1)) 

!---Radial Integration---!
ar = -connection(1)-2.d0*(connection(4)*vphi(1)+connection(6)*vr(1)*vtheta(1))-connection(5)*vr(1)*vr(1)&
-connection(7)*vtheta(1)*vtheta(1)-connection(11)*vphi(1)*vphi(1)+tc*vr(1)
vr(2) = ar*dt+vr(1)
r(2) = vr(2)*dt+r(1)

do temptemp = 1, 20
write(*,*) temptemp,connection(temptemp)
enddo
!Inside black hole checker
if(r(2) <= 2.d0*M) then
write(*,*) "Fell into black hole"
write(*,*) "Integration terminated"
close(11)
stop
endif

!Out of box checker
if(r(2)>boxsize) then
write(*,*) r(2)
ttemp = int(i/cut)
close(11)
goto 50
endif

!---Theta Integration---!
atheta = -connection(12)*vphi(1)*vphi(1)-connection(13)-2.d0*(connection(16)*vphi(1)+connection(18)*vtheta(1)*vr(1))&
-connection(17)*vr(1)*vr(1)-connection(19)*vtheta(1)*vtheta(1)+tc*vtheta(1)
vtheta(2)=atheta*dt+vtheta(1)
theta(2)=vtheta(2)*dt+theta(1)

!---Phi Integration---!
aphi = -2.d0*(connection(8)*vtheta(1)*vphi(1)+connection(10)*vr(1)*vphi(1)+connection(14)*vr(1)+connection(15)*vtheta(1))&
+tc*vphi(1)
vphi(2)=aphi*dt+vphi(1)
phi(2)=vphi(2)*dt+phi(1)


!writing new positions to file
if(mod(i,cut) .eq. 0) then
write(11,*) r(2),theta(2),phi(2)
endif

!updating variables for new step
vr(1)=vr(2)
r(1)=r(2)
vtheta(1)=vtheta(2)
theta(1)=theta(2)
vphi(1)=vphi(2)
phi(1)=phi(2)
connection(:) = 0.d0
enddo
close(11) !Closing output file


!---rthetaphi->xyz---!
ttemp = int(T/cut)
50 open(11,file='rthetaphi.out',status='old',action='read',iostat=error)
if(error .ne. 0) stop
open(12,file='xyz.out',status='unknown',action='write')
do i=1,ttemp,1
read(11,*)r(1),theta(1),phi(1)
write(12,*)sqrt(r(1)*r(1)+a*a)*sin(theta(1))*cos(phi(1)),sqrt(r(1)*r(1)+a*a)*sin(theta(1))*sin(phi(1)),r(1)*cos(theta(1))
enddo
close(12)
close(11)

end program main


!Connection subroutine
subroutine connections(connection,r,theta,M,a)
real*8, intent(in) :: r, theta, M, a
real*8, intent(out), dimension(20) :: connection
real*8 :: rs, invr, E, delta, invE, invdelta, costheta, sintheta

costheta = cos(theta)
sintheta = sin(theta)
E = r*r + a*a*costheta
invE = 1.d0/E
rs = 2.d0*M
delta = r*r - rs*r + a*a
invdelta = 1.d0/delta

connection(1) = rs*delta*(r*r-a*a*costheta*costheta)*0.5d0*invE**3 
connection(2) = rs*(r*r+a*a)*(r*r-a*a*costheta*costheta)*0.5d0*(invE**2)*invdelta
connection(3) = -rs*a*a*r*sintheta*costheta*invE*invE
connection(4) = delta*rs*a*sintheta*sintheta*(r*r-a*a*costheta*costheta)*0.5d0*invE**3
connection(5) = 2.d0*r*a*a*sintheta*sintheta-rs*(r*r-a*a*costheta)*0.5d0*invE*invdelta
connection(6) = -a*a*sintheta*costheta*invE
connection(7) = -r*delta*invE
connection(8) = (costheta/sintheta)*invE**2*(E*E+rs*a*a*r*sintheta*sintheta) 
connection(9) = rs*a*sintheta*sintheta*(a*a*costheta*costheta*(a*a-r*r)-r*r*(a*a+3.d0*r*r))*0.5d0*invE*invE*invdelta
connection(10) = (2.d0*r*E*E+rs*((a**4)*sintheta*sintheta*costheta*costheta-r*r*(E+r*r+a*a)))*0.5d0*invE*invE*invdelta
connection(11) = delta*sintheta*sintheta*0.5d0*(invE**3)*(-2.d0*r*E*E+rs*a*a*sintheta*sintheta*(r*r-a*a*costheta*costheta))
connection(12) = -sintheta*costheta*(invE**3)*((r*r+a*a)**2-a*a*delta*sintheta*sintheta)*E+(r*r+a*a)*rs*a*a*r*sintheta*sintheta
connection(13) = rs*a*a*r*sintheta*costheta*(invE**3) 
connection(14) = rs*a*(r*r-a*a*costheta*costheta)*0.5d0*invE*invE*invdelta
connection(15) = -rs*a*r*(costheta/sintheta)*invE*invE 
connection(16) = rs*a*r*(r*r+a*a)*sintheta*costheta*(invE**3)
connection(17) = a*a*sintheta*costheta*invE*invdelta
connection(18) = r*invE
connection(19) = -a*a*sintheta*costheta*invE
connection(20) = rs*a*a*a*r*sintheta*sintheta*sintheta*costheta*invE*invE


!t=0 r=1 theta=2 phi=3
end subroutine connections
