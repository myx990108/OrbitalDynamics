! *********************************************************************
! orbel_xv2el.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Given the cartesian position and velocity of an orbit, compute
!     the osculating orbital elements in the two-body problem.
!
! Remarks:
!     If the inclination INC is less than TINY, we
!     arbitrarily choose the longitude of the ascending node LGNODE
!     to be 0.0 (so the ascending node is then along the X axis).  If 
!     the  eccentricity E is less than SQRT(TINY), we arbitrarily
!     choose the argument of perihelion to be 0.
! *********************************************************************
! Input:
!     gmsum  ==>  G*(M1+M2)
!     x,y,z,vx,vy,vz  ==> Cartesian position and velovity
! Output:
!     a  ==>  semimajor axis
!     e  ==>  orbital eccentricity
!     inc  ==>  orbital inclination
!     capom  ==> longitude of ascending node (rad)
!     omega  ==> argument of pericenter (rad)
!     capm   ==> orbital mean anomaly (rad)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine orbel_xv2el(gmsum,x,y,z,vx,vy,vz,a,e,inc,capom,omega,capm)
    implicit none 
    ! Input/Output param
    real*8  gmsum,x,y,z,vx,vy,vz 
    real*8  a,e,inc,capom,omega,capm 
    ! Local param
    integer ialpha
    real*8  TINY,PI,TWOPI
    real*8  hx,hy,hz,h2,h,r,v2,v,vdotr,energy
    real*8  fac,u,face,cape,capf,tmpf,cw,sw,w
    parameter(TINY = 4.d-15)
    parameter(PI = 3.1415926535897932384d0)
    parameter(TWOPI = 2.d0 * PI)
    ! ... Execution ...
    ! Compute the angular momentum H, and thereby the inclination inc.
    hx = y*vz - z*vy 
    hy = z*vx - x*vz 
    hz = x*vy - y*vx 
    h2 = hx*hx + hy*hy + hz*hz 
    h = sqrt(h2)
    if(hz.gt.h) then 
        hz = h 
        hx = 0.d0 
        hy = 0.d0
    endif
    inc = acos(hz/h)
    ! Compute longitude of ascending node capom and the argument of latitude u
    fac = sqrt(hx*hx + hy*hy)/h 
    if((fac.lt.TINY).or.(inc.eq.0.d0)) then
        capom = 0.d0
        u = atan2(y,x)
        if(abs(inc-PI).lt.(10.d0*TINY)) u = -u 
    else
        capom = atan2(hx,-hy)
        u = atan2(z/sin(inc),x*cos(capom)+y*sin(capom))
    endif
    if(capom.lt.0.d0) capom = capom + TWOPI 
    if(u.lt.0.d0) u = u + TWOPI 
    ! Compute the radius R and velocity squared V2, and the dot
    ! product RDOTV , the energy pre unit mass ECERGY
    r = sqrt(x*x + y*y + z*z)
    v2 = vx*vx + vy*vy + vz*vz 
    v = sqrt(v2)
    vdotr = x*vx + y*vy + z*vz 
    energy = 0.5d0 * v2 - gmsum / r 
    ! Determine type of conic section and label it via IALPHA
    if(abs(energy*r/gmsum).lt.sqrt(TINY)) then 
        ialpha = 0
    else 
        if(energy.lt.0.d0) ialpha = -1
        if(energy.gt.0.d0) ialpha = +1
    endif
    ! Depending on the conic type, determine the remaining elements
    ! Ellipse orbits
    if(ialpha.eq.-1) then 
        a = -0.5*gmsum/energy 
        fac = 1.d0 - h2/(gmsum*a)
        if(fac.gt.TINY) then 
            e = sqrt(fac)
            face = (a-r)/(a*e) ! cos(E)
            if(face.gt.1.d0) then 
                cape = 0.d0    ! eccentric anomaly
            else 
                if(face.gt.-1.d0) then 
                    cape = acos(face)
                else 
                    cape = PI 
                endif
            endif
            if(vdotr.lt.0.d0) cape = TWOPI - cape 
            cw = (cos(cape)-e)/(1.d0-e*cos(cape)) ! x/r
            sw = sqrt(1.d0-e*e)*sin(cape)/(1.d0-e*cos(cape)) ! y/r
            w = atan2(sw,cw)
            if(w.lt.0.d0) w = w + TWOPI
        else 
            e = 0.d0 
            w = u 
            cape = u
        endif
        capm = cape - e * sin(cape)
        omega = u - w 
        if(omega.lt.0.d0) omega = omega + TWOPI 
        omega = omega - int(omega/TWOPI)*TWOPI
    endif
    ! Hyperbola orbits
    if(ialpha.eq.+1) then 
        a = +0.5d0*gmsum/energy
        fac = h2/(gmsum*a)
        if(fac.gt.TINY) then 
            e = sqrt(1.d0+fac)
            tmpf = (a+r)/(a*e)
            if(tmpf.lt.1.d0) tmpf = 1.d0
            capf = log(tmpf + sqrt(tmpf*tmpf - 1.d0))
            if(vdotr.lt.0.d0) capf = -capf 
            cw = (e - cosh(capf))/(e*cosh(capf) - 1.d0)
            sw = sqrt(e*e - 1.d0)*sinh(capf)/(e*cosh(capf) - 1.d0)
            w = atan2(sw,cw)
            if(w.lt.0.d0) w = w + TWOPI 
        else 
            e = 1.d0 
            tmpf = 0.5d0*h2/gmsum 
            w = acos(2.d0*tmpf/r - 1.d0)
            if(vdotr.lt.0.d0) w = TWOPI - w 
            tmpf = (a+r)/(a*e)
            capf = log(tmpf + sqrt(tmpf*tmpf -1.d0))
        endif
        capm = e*sinh(capf) - capf 
        omega = u - w 
        if(omega.lt.0.d0) omega = omega + TWOPI 
        omega = omega - int(omega/TWOPI)*TWOPI 
    endif
    ! Parabolic orbits : ( NOTE - in this case we use "a" to mean pericentric distance : P/2)
    if(ialpha.eq.0) then 
        a = 0.5d0*h2/gmsum 
        e = 1.d0 
        w = acos(2.d0*a/r-1.d0)
        if(vdotr.lt.0.d0) w = TWOPI - w 
        tmpf = tan(0.5d0 * w)
        capm = tmpf * (1.d0 + tmpf*tmpf/3.d0) !! according to orbel_zget(q)
        omega = u - w 
        if(omega .lt. 0.d0) omega = omega + 2.d0*PI
	    omega = omega - int(omega/(2.d0*PI))*2.d0*PI 	 
    endif
    return
end subroutine

! *********************************************************************
! orbel_el2xv.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Given the osculating orbital elements of an orbit, compute
!     the cartesian position and velocity in the two-body problem.
!
! Remarks:
!     If the inclination INC is less than TINY, we
!     arbitrarily choose the longitude of the ascending node LGNODE
!     to be 0.0 (so the ascending node is then along the X axis).  If 
!     the  eccentricity E is less than SQRT(TINY), we arbitrarily
!     choose the argument of perihelion to be 0.
! *********************************************************************
! Input:
!     gmsum  ==>  G*(M1+M2)
!     a  ==>  semimajor axis
!     e  ==>  orbital eccentricity
!     inc  ==>  orbital inclination
!     capom  ==> longitude of ascending node (rad)
!     omega  ==> argument of pericenter (rad)
!     capm   ==> orbital mean anomaly (rad)
! Output:
!     x,y,z,vx,vy,vz ==> cartesian positions and velocities
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine orbel_el2xv(gmsum,a,e,inc,capom,omega,capm,x,y,z,vx,vy,vz)
    implicit none 
    ! Input/Output param
    real*8  gmsum,a,e,inc,capom,omega,capm 
    real*8  x,y,z,vx,vy,vz 
    ! Local param
    real*8  orbel_ehybrid,orbel_fhybrid,orbel_zget
    integer ialpha
    real*8  TINY 
    real*8  em1,sp,cp,so,co,si,ci
    real*8  cape,scap,ccap,sqe,sqgma
    real*8  capf,shcap,chcap,zpara
    real*8  d11,d12,d13,d21,d22,d23
    real*8  xfac1,xfac2,vfac1,vfac2,ri
    parameter(TINY = 4.d-15)
    ! ... Execution ...
    ! Check whether eccentricity e < 0
    if(e.lt.0.d0) then 
        write(6,"(a)") "ORBEL_EL2XV : ERROE FOR ECC < 0, SETTING ECC = 0!"
        e = 0.d0 
    endif
    ! Check for the type of conic section
    em1 = e - 1.d0
    if(abs(em1).lt.TINY) ialpha = 0
    if(em1.gt.TINY) ialpha = +1
    if(em1.lt.-TINY) ialpha = -1
    ! Generate rotation matrices
    sp = sin(omega)
    cp = cos(omega)
    so = sin(capom)
    co = cos(capom)
    si = sin(inc)
    ci = cos(inc)
    d11 = cp*co - sp*so*ci
	d12 = cp*so + sp*co*ci
	d13 = sp*si
	d21 = -sp*co - cp*so*ci
	d22 = -sp*so + cp*co*ci
	d23 = cp*si
    ! Get the other quantities depending on the type of orbital type (ialpha)
    ! Elliptic orbits
    if(ialpha.eq.-1) then 
        cape = orbel_ehybrid(e,capm) ! Mean anomaly is normalized to [0,2pi]
        scap = sin(cape)
        ccap = cos(cape)
        sqe = sqrt(1.d0 - e*e)
        sqgma = sqrt(gmsum*a)
        xfac1 = a*(ccap - e)
        xfac2 = a*sqe*scap 
        ri = 1.d0 / (a*(1.d0 - e*ccap))
        vfac1 = -ri * sqgma * scap 
        vfac2 =  ri * sqgma * sqe * ccap 
    endif
    ! Hyperbolic orbits
    if(ialpha.eq.+1) then 
        capf = orbel_fhybrid(e,capm)
        shcap = sinh(capf)
        chcap = cosh(capf)
        sqe = sqrt(e*e - 1.d0)
        sqgma = sqrt(gmsum*a)
        xfac1 = a*(e - chcap)
        xfac2 = a*sqe*shcap 
        ri = 1.d0 / (a*(e*chcap - 1.d0))
        vfac1 = -ri * sqgma * shcap
	    vfac2 =  ri * sqgma * sqe * chcap
    endif
    ! Parabolic orbits
    if(ialpha.eq.0) then 
        zpara = orbel_zget(capm)
        sqgma = sqrt(2.d0*gmsum*a)
        xfac1 = a*(1.d0 - zpara*zpara)
        xfac2 = 2.d0*a*zpara 
        ri = 1.d0 / (a*(1.d0 + zpara*zpara))
        vfac1 = -ri * sqgma * zpara
        vfac2 =  ri * sqgma
    endif

    x =  d11*xfac1 + d21*xfac2
	y =  d12*xfac1 + d22*xfac2
	z =  d13*xfac1 + d23*xfac2
	vx = d11*vfac1 + d21*vfac2
	vy = d12*vfac1 + d22*vfac2
	vz = d13*vfac1 + d23*vfac2
    return
end subroutine

! *********************************************************************
! orbel_esolmd.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for elliptic orbits with 0 < e < 1 using 
!     some sort of quartic convergence from Wisdom.
! 
! Kepler function is E - e*sin(E) = M = n*(t-\tau)
! This function is only good for low values of e (0 < e < 0.18) since 
! it only iterates once.
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (0 < e < 0.18)
!     capm  ==>  ellipse mean anomaly
! Output:
!     orbel_esolmd  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_esolmd(e,capm)
    implicit none 
    ! Input/Output param
    real*8  e,capm 
    real*8  orbel_esolmd
    ! Local param
    real*8  x,sm,cm,sx,cx,es,ec
    real*8  f,fp,fpp,fppp,dx
    ! ... Execution ...
    ! Function to solve Kepler's eqn for E (here called
    ! x) for given e and M. returns value of x.
    sm = sin(capm)
    cm = cos(capm)
    ! Begin with a guess value accurate to order eccen**3.d0
    x = capm + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))
    ! Iteration once with Newton's method
    sx = sin(x)
    cx = cos(x)
    es = e*sx 
    ec = e*cx 
    f = x - es - capm 
    fp = 1.d0 - ec 
    fpp = es 
    fppp = ec 
    dx = -f/fp 
    dx = -f/(fp + dx*fpp/2.d0)
    dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
    orbel_esolmd = x + dx 
    return
end function

! *********************************************************************
! orbel_eget.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for elliptic orbits with 0 < e < 1 using 
!     quartic convergence from Danby.
! 
! Kepler function is E - e*sin(E) = M = n*(t-\tau)
! This function is only good for low values of e (0.15 < e < 0.8) since 
! it iterates twice.
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (0.15 < e < 0.8)
!     capm  ==>  ellipse mean anomaly
! Output:
!     orbel_eget  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_eget(e,capm)
    implicit none 
    ! Input/Output param
    real*8  e,capm 
    real*8  orbel_eget 
    ! Local param
    real*8  x,sm,cm,sx,cx,es,ec
    real*8  f,fp,fpp,fppp,dx
    ! ... Execution ...
    ! Function to solve Kepler's eqn for E (here called
    ! x) for given e and M. returns value of x.
    sm = sin(capm)
    cm = cos(capm)
    ! Begin with a guess value accurate to order eccen**3.d0
    x = capm + e*sm*( 1.d0 + e*( cm + e*( 1.d0 -1.5d0*sm*sm)))
    ! Iterate once
    sx = sin(x)
    cx = cos(x)
    es = e*sx
    ec = e*cx 
    f = x - es - capm
    fp = 1.d0 - ec 
    fpp = es 
    fppp = ec 
    dx = -f/fp 
    dx = -f/(fp + dx*fpp/2.d0)
    dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
    orbel_eget = x + dx 
    ! Iterate twice
    x = orbel_eget
    sx = sin(x)
    cx = cos(x)
    es = e*sx
    ec = e*cx 
    f = x - es - capm
    fp = 1.d0 - ec 
    fpp = es 
    fppp = ec 
    dx = -f/fp 
    dx = -f/(fp + dx*fpp/2.d0)
    dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
    orbel_eget = x + dx 
    return
end function

! *********************************************************************
! orbel_ehie.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for elliptic orbits with 0 < e < 1 using 
!     Danby's quartic convergence method for 3 iterations. Eq is 
!     f(x) = x - e*sin(x+M), where E = x + M
! 
! Kepler function is E - e*sin(E) = M = n*(t-\tau)
! This function is only good for low values of e (e > 0.8) since 
! it iterates twice.
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (e > 0.8)
!     capm  ==>  ellipse mean anomaly
! Output:
!     orbel_ehie  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_ehie(e,capm)
    implicit none
    ! Input/Output param
    real*8  e,capm 
    real*8  orbel_ehie
    ! Local param
    integer iflg,nper,niter,NMAX 
    real*8  m,dx,x,sa,ca,esa,eca,f,fp
    real*8  PI,TWOPI 
    parameter(NMAX = 3)
    parameter(PI = 3.1415926535897932384d0)
    parameter(TWOPI = 2.d0 * PI)
    ! ... Execution ...
    ! In this section, bring M into the range (0,TWOPI) and if
    ! the result is greater than PI, solve for (TWOPI - M).
    iflg = 0
    nper = capm/TWOPI
    m = capm - nper*TWOPI 
    if(m.lt.0.d0) m = m + TWOPI 
    if(m.gt.PI) then 
        iflg = 1
        m = TWOPI - m 
    endif
    ! Make a first guess that works well for e near 1
    x = (6.d0*m)**(1.d0/3.d0) - m
    niter = 0
    ! Iteration Loop
    do niter = 1,NMAX 
        sa = sin(x+m)
        ca = cos(x+m)
        esa = e*sa 
        eca = e*ca 
        f = x - esa 
        fp = 1.d0 - eca 
        dx = -f/fp
        dx = -f/(fp + dx*esa/2.d0)
        dx = -f/(fp + dx*esa/2.d0 + dx*dx*eca/6.d0)
        x = x + dx
    enddo
    orbel_ehie = x + m 
    if(iflg.eq.1) then 
        orbel_ehie = TWOPI - orbel_ehie 
        m = TWOPI - m 
    endif
    capm = m
    return
end function

! *********************************************************************
! orbel_ehybrid.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for elliptic orbits with 0 < e < 1 using 
!     Danby's quartic convergence method.
! 
! Kepler function is E - e*sin(E) = M = n*(t-\tau)
! This function calls orbel_esolmd when 0 < e < 0.18, calls orbel_eget
! when 0.18 < e < 0.8, and calls orbel_ehie when e > 0.8
!
! The mean anomaly is normalized to [0,TWOPI]
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (0 < e < 1)
!     capm  ==>  ellipse mean anomaly
! Output:
!     orbel_ehie  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_ehybrid(e,capm)
    implicit none 
    ! Input/Output param
    real*8  e,capm
    real*8  orbel_ehybrid 
    ! Local param
    real*8  orbel_esolmd,orbel_eget,orbel_ehie 
    real*8  PI,TWOPI 
    parameter(PI = 3.1415926535897932384d0)
    parameter(TWOPI = 2.d0 * PI)
    ! ... Execution ...
    ! Normalization
    capm = capm - int(capm/TWOPI) * TWOPI
    if(capm.lt.0.d0) capm = capm + TWOPI 
    ! Check the value of eccectricity
    if(e.lt.0.18d0) then 
        orbel_ehybrid = orbel_esolmd(e,capm)
    elseif(e.lt.0.8d0) then 
        orbel_ehybrid = orbel_eget(e,capm)
    else 
        orbel_ehybrid = orbel_ehie(e,capm)
    endif
end function

! *********************************************************************
! orbel_flon.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for hyperbolic orbits with e > 1 using 
!     power series for capn in terms of F and Newton's method
!
! Kepler function is  e*sinh(F) - F = n*(t - \tau)
! This function is only good for low values of capn (capn < 0.636*e - 0.6)
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (e > 1)
!     capn  ==>  hyperbola mean anomaly
! Output:
!     orbel_flon  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_flon(e,capn)
    implicit none
    ! Input/Output param
    real*8  e,capn 
    real*8  orbel_flon 
    ! Local param
    integer iflg,i,IMAX
    real*8  TINY
    real*8  a0,a1,a3,a5,a7,a9,a11
    real*8  b1,b3,b5,b7,b9,b11
    real*8  a,b,sq,biga,bigb
    real*8  x,x2,f,fp,dx,diff
    parameter(TINY = 4.d-15)
    parameter(IMAX = 10)
    parameter(a11 = 156.d0,a9 = 17160.d0,a7 = 1235520.d0)
	parameter(a5 = 51891840.d0,a3 = 1037836800.d0)
	parameter(b11 = 11.d0*a11,b9 = 9.d0*a9,b7 = 7.d0*a7)
	parameter(b5 = 5.d0*a5, b3 = 3.d0*a3)
    ! ... Execution ...
    ! Function to solve "Kepler's eqn" for F (here called x) for given e and CAPN.
    ! Only good for smallish CAPN 
    iflg = 0
    if(capn.lt.0.d0) then 
        iflg = 1
        capn = -capn
    endif
    ! Initial guess value
    a1 = 6227020800.d0 * (1.d0 - 1.d0/e)
    a0 = -6227020800.d0*capn/e
    b1 = a1 

    a = 6.d0*(e-1.d0)/e
    b = -6.d0*capn/e
    sq = sqrt(0.25*b*b +a*a*a/27.d0)
    biga = (-0.5*b + sq)**0.3333333333333333d0
    bigb = -(+0.5*b + sq)**0.3333333333333333d0
    x = biga + bigb 
    orbel_flon = x
    ! When capn is TINY (or zero), no need to go further than cubic 
    ! even for e = 1.
    if(capn.lt.TINY) then 
        if(iflg.eq.1) then 
            capn = -capn 
            orbel_flon = -orbel_flon
        endif
        return
    endif
    ! Iterations for Newton's method
    do i = 1,IMAX
        x2 = x * x 
        f = a0 +x*(a1+x2*(a3+x2*(a5+x2*(a7+x2*(a9+x2*(a11+x2))))))
        fp = b1 +x2*(b3+x2*(b5+x2*(b7+x2*(b9+x2*(b11 + 13.d0*x2)))))
        dx = -f/fp
        orbel_flon = x + dx 
        ! If converged, then return
        if(abs(dx).le.TINY) then 
            if(iflg.eq.1) then 
                capn = -capn 
                orbel_flon = -orbel_flon
            endif
            return 
        endif
        ! Continue to iterate
        x = orbel_flon
    enddo
    ! Abnormal return here - we've gone through the loop IMAX 
    ! times without convergence
    write(6,"(a)") "ORBEL_FLON : RETURNING WITHOUT COMPLETE CONVERGENCE"
	if(iflg.eq.1) then
	    orbel_flon = -orbel_flon
	    capn = -capn
	endif
	return
end function

! *********************************************************************
! orbel_fget.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for hyperbolic orbits with e > 1 using
!     hybrid approach
!
! Kepler function is  e*sinh(F) - F = n*(t - \tau)
! This function is good for high values of capn (capn > 0.636*e - 0.6)
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (e > 1)
!     capn  ==>  hyperbola mean anomaly
! Output:
!     orbel_fget  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_fget(e,capn)
    implicit none
    ! Input/Output param
    real*8  e,capn 
    real*8  orbel_fget 
    ! Local param
    integer i,IMAX 
    real*8  TINY
    real*8  tmp,x,shx,chx
    real*8  esh,ech,f,fp,fpp,fppp,dx
    parameter(IMAX = 10)
    parameter(TINY = 4.d-15)
    ! ... Execution ...
    ! Function to solve "Kepler's eqn" for F (here called x) for given e and CAPN.
    ! Begin with a guess proposed by Danby
    if(capn.lt.0.d0) then 
        tmp = -2.d0*capn/e + 1.8d0
        x = -log(tmp)
    else 
        tmp = +2.d0*capn/e + 1.8d0 
        x = log(tmp)
    endif
    orbel_fget = x 
    ! Iterations for Newton's method
    do i = 1,IMAX 
        shx = sinh(x)
        chx = cosh(x)
        esh = e*shx 
        ech = e*chx 
        f = esh - x - capn 
        fp = ech - 1.d0 
        fpp = esh 
        fppp = ech 
        dx = -f/fp 
        dx = -f/(fp + dx*fpp/2.d0)
        dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
        orbel_fget = x + dx 
        ! If converged, then return
        if(abs(dx).le.TINY) return 
        ! Continue to iterate
        x = orbel_fget
    enddo
    write(6,"(a)") "ORBEL_FGET : RETURNING WITHOUT COMPLETE CONVERGENCE!"
    return
end function

! *********************************************************************
! orbel_fhybrid.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for hyperbolic orbits with e > 1 using
!     hybrid approach
!
! Kepler function is  e*sinh(F) - F = n*(t - \tau)
! This function calls orbel_flon when capn < 0.636*e - 0.6, and calls
! orbel_fget when capn > 0.636*e - 0.6
! *********************************************************************
! Input:
!     e  ==>  orbital eccentricity (e > 1)
!     capn  ==>  hyperbola mean anomaly
! Output:
!     orbel_fhybrid  ==> eccentricity anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_fhybrid(e,capn)
    implicit none 
    ! Input/Output param
    real*8  e,capn
    real*8  orbel_fhybrid 
    ! Local param
    real*8  orbel_fget,orbel_flon 
    real*8  abn
    ! ... Execution ...
    ! Check the value of capn
    if(abs(capn).lt.0.636d0*e - 0.6d0) then 
        orbel_fhybrid = orbel_flon(e,capn)
    else 
        orbel_fhybrid = orbel_fget(e,capn)
    endif
    return
end function

! *********************************************************************
! orbel_zget.f90 (function)
! *********************************************************************
! Purpose:
!     Solve Kepler function for parabolic orbits with e = 1
!
! Kepler function is tan(f/2) + 1/3*tan^3(f/2) = M 
! where M = sqrt(mu/2*q**3)*(t-\tau), the solution is E = tan(f/2)
! *********************************************************************
! Input:
!     q  ==>  parabola mean anomaly
! Output:
!     orbel_zget ==> eccentric anomaly
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
function orbel_zget(q)
    implicit none
    ! Input/Output param
    real*8  q,orbel_zget 
    ! Local param
    integer iflg 
    real*8  x,tmp
    ! ... Execution ...
    ! Check the sign of q
    iflg = 0
    if(q.lt.0.d0) then 
        q = -q 
        iflg = 1
    endif
    ! Analytical solution
    if(q.lt.1.d-3) then 
        orbel_zget = q*(1.d0 - (q*q/3.d0)*(1.d0 - q*q))
    else 
        x = 0.5d0 * (3.d0*q + sqrt(9.d0*q*q + 4.d0))
        tmp = x ** (1.d0/3.d0)
        orbel_zget = tmp - 1.d0 / tmp 
    endif
    ! Recover the sign of q 
    if(iflg.eq.1) then 
        orbel_zget = -orbel_zget
        q = -q 
    endif 
    return
end function