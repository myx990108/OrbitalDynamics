! *********************************************************************
! drift_one.f90 (subroutine)
! *********************************************************************
! Purpose:
!     This subroutine does the danby-type drift for one particle, using 
!     appropriate vbles and redoing a drift if the accuracy is too poor 
!     (as flagged by the integer iflg).
! *********************************************************************
! Input:
!     mu  ==> mass
!     x,y,z,vx,vy,vz  ==> position and velocity of the particle
!     dt  ==> time step
! Output:
!     iflg  ==>  (==0 if succeeded)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_one(mu,x,y,z,vx,vy,vz,dt,iflg)
    implicit none 
    ! Input/Output param
    real*8  mu
    real*8  x,y,z,vx,vy,vz,dt
    integer iflg 
    ! Local param
    real*8  dttmp 
    integer i 
    ! ... Execution ...
    call drift_danby(mu,x,y,z,vx,vy,vz,dt,iflg)
    if(iflg.ne.0) then 
        do i = 1,10
            dttmp = dt/10.d0 
            call drift_danby(mu,x,y,z,vx,vy,vz,dttmp,iflg)
        enddo
    endif
    return
end subroutine

! *********************************************************************
! drift_danby.f90 (subroutine)
! *********************************************************************
! Purpose:
!     This subroutine does the Danby and decides which variables to use
! *********************************************************************
! Input:
!     mu  ==> mass of the body
!     x0,y0,z0  ==> initial position in jacobi coords
!     vx0,vy0,vz0  ==> initial velocity in jacobi coords
!     dt0  ==> time step
! Output:
!     iflg  ==> integer flag (==0 if converged)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_danby(mu,x0,y0,z0,vx0,vy0,vz0,dt0,iflg)
    implicit none 
    ! Input/Output param
    real*8  mu,x0,y0,z0,vx0,vy0,vz0,dt0
    integer iflg
    ! Local param
    real*8  PI,TWOPI,DANBYB
    parameter(PI = 3.1415926535897932384d0)
    parameter(TWOPI = 2.d0*PI)
    parameter(DANBYB = 1.d-13)
    real*8  dt
    real*8  r0,v0s,u,alpha
    real*8  a,asq,en,ec,es,esq,dm,fp
    real*8  f,g,fdot,gdot,c1,c2,c3
    real*8  x,y,z,vx,vy,vz 
    real*8  xkep,fchk,s,c
    ! ... Execution ...
    ! Set dt = dt0 to be sure timestep is not altered while solving
    ! for new coords.
    dt = dt0 
    iflg = 0
    r0 = sqrt(x0*x0 + y0*y0 + z0*z0)
    v0s = vx0*vx0 + vy0*vy0 + vz0*vz0 
    u = x0*vx0 + y0*vy0 + z0*vz0 
    alpha = 2.d0*mu/r0 - v0s 
    ! Check for the binding energy
    if(alpha.gt.0.d0) then 
        a = mu/alpha 
        asq = a*a 
        en = sqrt(mu/(a*asq))
        ec = 1.d0 - r0/a
        es = u/(en*asq)
        esq = ec*ec + es*es 
        dm = en*dt - int(en*dt/TWOPI)*TWOPI 
        dt = dm/en
        if((dm*dm.gt.0.16d0).or.(esq.gt.0.36d0)) then 
            call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
            if(iflg.eq.0) then 
                f = 1.d0 - (mu/r0)*c2
                g = dt - mu*c3 
                fdot = -(mu/(fp*r0))*c1
                gdot = 1. - (mu/fp)*c2

                x = x0*f + vx0*g
                y = y0*f + vy0*g
                z = z0*f + vz0*g
                vx = x0*fdot + vx0*gdot 
                vy = y0*fdot + vy0*gdot 
                vz = z0*fdot + vz0*gdot 

                x0 = x 
                y0 = y 
                z0 = z 
                vx0 = vx 
                vy0 = vy
                vz0 = vz 
            endif
            return
        endif

        if(esq*dm*dm.lt.0.0016d0) then 
            call drift_kepmd(dm,es,ec,xkep,s,c)
            fchk = (xkep - ec*s + es*(1.d0-c) - dm)
            if(fchk*fchk.gt.DANBYB*DANBYB) then 
                iflg = 1
                return
            endif
            fp = 1.d0 - ec*c + es*s 
            f = (a/r0) * (c-1.d0) + 1.d0 
            g = dt + (s-xkep)/en 
            fdot = - (a/(r0*fp))*en*s
            gdot = (c-1.d0)/fp + 1.d0 

            x = x0*f + vx0*g
            y = y0*f + vy0*g
            z = z0*f + vz0*g
            vx = x0*fdot + vx0*gdot 
            vy = y0*fdot + vy0*gdot 
            vz = z0*fdot + vz0*gdot 

            x0 = x 
            y0 = y 
            z0 = z 
            vx0 = vx 
            vy0 = vy
            vz0 = vz 
            iflg = 0
        endif
    endif

    call drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
    if(iflg.eq.0) then 
        f = 1.d0 - (mu/r0)*c2
        g = dt - mu*c3 
        fdot = -(mu/(fp*r0))*c1
        gdot = 1. - (mu/fp)*c2

        x = x0*f + vx0*g
        y = y0*f + vy0*g
        z = z0*f + vz0*g
        vx = x0*fdot + vx0*gdot 
        vy = y0*fdot + vy0*gdot 
        vz = z0*fdot + vz0*gdot 

        x0 = x 
        y0 = y 
        z0 = z 
        vx0 = vx 
        vy0 = vy
        vz0 = vz 
    endif
    return
end subroutine

! *********************************************************************
! drift_kepu.f90 (subroutine)
! *********************************************************************
! Purpose:
!     subroutine for solving kepler's equation using universal variables.
! *********************************************************************
! Input:
!     dt  ==> time step
!     r0  ==> Distance between the center body and particles
!     mu  ==> Reduces mass of system
!     alpha  ==> twice the binding energy
!     u  ==> Vel. dot radial vector
! Output:
!     fp  ==> f' from p170 of Danby.1988
!     c1,c2,c3  ==> c's from p171-172 of Danby.1988
!     iflg  ==> convergence flag (==0 if succeeded)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu(dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
    implicit none 
    ! Input/Output param
    integer iflg
    real*8  dt,r0,mu,alpha,u 
    real*8  fp,c1,c2,c3
    ! Local param
    real*8  s,st,fo,fn
    ! ... Execution ...
    ! Get the guess initial value for s
    call drift_kepu_guess(dt,r0,mu,alpha,u,s)
    ! Store initial guess for possible use later in
    ! laguerre's method, in case newton's method fails.
    st = s 
    call drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
    if(iflg.ne.0) then 
        call drift_kepu_fchk(dt,r0,mu,alpha,u,st,fo)
        call drift_kepu_fchk(dt,r0,mu,alpha,u,s,fn)
        if(abs(fo).lt.abs(fn)) then 
            s = st 
        endif
        call drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
    endif
    return
end subroutine

! *********************************************************************
! drift_kepmd.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Solve Kepler function in difference form for an ellipse, given 
!     small dm and small eccentricity. See drift_dan.f90 for the criteria.
! 
! The Kepler function in difference form is:
!     x - e*cos(E0)*sin(x) + e*sin(E0)*(1-cos(x)) - dm = 0
! *********************************************************************
! Input:
!     dm  ==>  increment in mean anomaly M
!     es,ec  ==>  e*sin(E0), e*cos(E0)
! Output:
!     x  ==> solution to Kepler function in difference form
!     s,c  ==> sin(x), cos(x)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepmd(dm,es,ec,x,s,c)
    implicit none 
    ! Input/Output param
    real*8  dm,es,ec,x,s,c 
    ! Local param
    real*8  dx,fac1,fac2,q,y 
    real*8  f,fp,fpp,fppp
    ! ... Execution ...
    ! Compute initial guess for root x
    fac1 = 1.d0 / (1.d0 - ec)
    q = fac1 * dm 
    fac2 = es*es*fac1 - ec/3.d0
    x = q*(1.d0 - 0.5d0*fac1*q*(es - q*fac2))
    ! Sin and Cos value of x 
    s = sin(x)
    c = cos(x)
    ! Compute better value for the root using quartic Newton's method
    f = x - ec*s + es*(1.d0-c) - dm
    fp = 1.d0 - ec*c + es*s
    fpp = ec*s + es*c 
    fppp = ec*c - es*s
    dx = -f/fp 
    dx = -f/(fp + dx*fpp/2.d0)
    dx = -f/(fp + dx*fpp/2.d0 + dx*dx*fppp/6.d0)
    x = x + dx 
    ! Sin and Cos value of x 
    s = sin(x)
    c = cos(x)
    return
end subroutine

! *********************************************************************
! drift_kepu_stumpff.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Subroutine for the calculation of stumpff functions.
! 
! The stumpff functions are as follows:
!     See Danby.1988 p.172  equations 6.9.15
! *********************************************************************
! Input:
!     x  ==> argument
! Output:
!     c0,c1,c2,c3  ==> c's from p171-172 (Danby.1988)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu_stumpff(x,c0,c1,c2,c3)
    implicit none 
    ! Input/Output param
    real*8  x 
    real*8  c0,c1,c2,c3 
    ! Local param
    integer n,i 
    real*8  xm
    ! ... Execution ...
    n = 0
    xm = 0.1
    do while(abs(x).ge.xm)
        n = n + 1
        x = x / 4.d0
    enddo
    ! Find c2 and c3 using series
    c2 = (1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x/182.d0)/132.d0)/90.d0)/56.d0)/30.d0)/12.d0)/2.d0
    c3 = (1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x*(1.d0-x/210.d0)/156.d0)/110.d0)/72.d0)/42.d0)/20.d0)/6.d0
    ! Find c1 and c0
    c1 = 1.d0 - x*c3 
    c0 = 1.d0 - x*c2 
    if(n.ne.0) then 
        do i = n,1,-1
            c3 = (c2 + c0*c3)/4.d0
            c2 = c1*c1/2.d0
            c1 = c0*c1 
            c0 = 2*c0*c0 - 1.d0
            x = x * 4.d0
        enddo
    endif
    return 
end subroutine

! *********************************************************************
! drift_kepu_p3solve.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Returns the real root of cubic often found in solving kepler
!     problem in universal variables.
! *********************************************************************
! Input:
!     dt  ==> time step
!     r0  ==> Distance between the center body and particles
!     mu  ==> Reduces mass of system
!     alpha  ==> Twice the binding energy
!     u  ==> Vel. dot radial vector
! Output:
!     s  ==> solution of cubic eqn for the universal variable
!     iflg  ==> success flag ( = 0 if O.K.) (integer)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
    implicit none 
    ! Input/Output param
    integer iflg
    real*8  dt,r0,mu,alpha,u 
    real*8  s
    ! Local param
    real*8  denom,a2,a1,a0
    real*8  q,r,sq2,sq,p1,p2
    ! ... Execution ...
    ! The coefficient a2, a1, a0 of the polynomial f2(s) from Danby.1988 P177.6.9.38
    denom = (mu - alpha*r0)/6.d0
    a2 = 0.5d0*u/denom 
    a1 = r0/denom
    a0 = -dt/denom 
    ! The process of solution
    q = (a1 - a2*a2/3.d0)/3.d0
    r = (a1*a2 - 3.d0*a0)/6.d0 - a2*a2*a2/27.d0
    sq2 = q*q*q + r*r 
    if(sq2.ge.0.d0) then 
        sq = sqrt(sq2)
        if((r+sq).le.0.d0) then 
            p1 =  -(-(r+sq))**(1.d0/3.d0)
        else 
            p1 = (r+sq)**(1.d0/3.d0)
        endif
        if((r-sq).le.0.d0) then 
            p2 = -(-(r-sq))**(1.d0/3.d0)
        else 
            p2 = (r-sq)**(1.d0/3.d0)
        endif
        iflg = 0
        s = p1 + p2 - a2/3.d0 
    else 
        iflg = 1
        s = 0
    endif
    return
end subroutine

! *********************************************************************
! drift_kepu_guess.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Initial guess for solving kepler's function using universal variables
!     See Danby.1988 P176
! 
! Select for elliptic case and hyperbolic case by alpha (binding energy)
! *********************************************************************
! Input:
!     dt  ==> time step
!     r0  ==> Distance between the center body and particles
!     mu  ==> Reduces mass of system
!     alpha  ==> Twice the binding energy
!     u  ==> Vel. dot radial vector
! Output:
!     s  ==> Initial guess for the value of universal variable
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu_guess(dt,r0,mu,alpha,u,s)
    implicit none 
    ! Input/Output param
    real*8  dt,r0,mu,alpha,u
    real*8  s 
    ! Local param
    integer iflg
    real*8  a,en,ec,es,e,y,sy,cy
    real*8  sigma,x
    ! ... Execution ...
    ! Find initial guess for elliptic orbits
    if(alpha.gt.0.d0) then 
        if(dt/r0 .le. 0.4d0)  then
            s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0)
	        return
        else    
            a = mu/alpha
            en = sqrt(mu/(a*a*a))
            ec = 1.d0 - r0/a  ! e * cos(E)
            es = u/(en*a*a)   ! e * sin(E)
            e = sqrt(ec*ec + es*es)
            y = en*dt - es
            sy = sin(y)
            cy = cos(y)
            sigma = dsign(1.d0,(es*cy + ec*sy))
            x = y + sigma*0.85d0*e
            s = x/sqrt(alpha)
        endif
    else 
    ! Find initial guess for hyperbolic motion
        call drift_kepu_p3solve(dt,r0,mu,alpha,u,s,iflg)
        if(iflg.ne.0) s = dt/r0
    endif
    return
end subroutine

! *********************************************************************
! drift_kepu_new.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Subroutine for solving kepler's equation in universal variables.
!     using NEWTON'S METHOD
! 
! Algorithm: Danby.1988 Quartic Convergence
! *********************************************************************
! Input:
!     s  ==> initial value of universal variable
!     dt  ==> time step
!     r0  ==> Distance between the center body and particles
!     mu  ==> Reduces mass of system
!     alpha  ==> Twice the binding energy
!     u  ==> Vel. dot radial vector
! Output:
!     s  ==> final value of universal variable
!     fp  ==> f' from p170 of Danby.1988
!     c1,c2,c3  ==> c's from p171-172 of Danby.1988
!     iflg  ==> convergence flag (=0 if succeeded)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu_new(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
    implicit none 
    ! Input/Output param
    real*8  s,dt,r0,mu,alpha,u 
    real*8  fp,c1,c2,c3 
    integer iflg 
    ! Local param
    integer nc
    real*8  x,c0
    real*8  f,fpp,fppp,ds,fdt
    real*8  DANBYB 
    parameter(DANBYB = 1.0d-13)
    ! ... Execution ...
    do nc = 0,6
        x = alpha * s * s 
        call drift_kepu_stumpff(x,c0,c1,c2,c3)
        c1 = c1*s 
        c2 = c2*s*s 
        c3 = c3*s*s*s 
        f = r0*c1 + u*c2 + mu*c3 - dt 
        fp = r0*c0 + u*c1 + mu*c2 
        fpp = (-r0*alpha + mu)*c1 + u*c0 
        fppp = (-r0*alpha + mu)*c0 - u*alpha*c1 
        ds = -f/fp 
        ds = -f/(fp + ds*fpp/2.d0)
        ds = -f/(fp + ds*fpp/2.d0 + ds*ds*fppp/6.d0)
        s = s + ds 
        fdt = f/dt
    ! Quartic convergence
        if(fdt*fdt.lt.DANBYB*DANBYB) then 
            iflg = 0
            return
        endif
    enddo
    iflg = 1
    return
end subroutine

! *********************************************************************
! drift_kepu_lag.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Subroutine for solving kepler's equation in universal variables.
!     using LAGUERRE'S METHOD
! 
! Algorithm: Danby.1988 Quartic Convergence
! *********************************************************************
! Input:
!     s  ==> initial value of universal variable
!     dt  ==> time step
!     r0  ==> Distance between the center body and particles
!     mu  ==> Reduces mass of system
!     alpha  ==> Twice the binding energy
!     u  ==> Vel. dot radial vector
! Output:
!     s  ==> final value of universal variable
!     fp  ==> f' from p170 of Danby.1988
!     c1,c2,c3  ==> c's from p171-172 of Danby.1988
!     iflg  ==> convergence flag (=0 if succeeded)
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu_lag(s,dt,r0,mu,alpha,u,fp,c1,c2,c3,iflg)
    implicit none 
    ! Input/Output param
    real*8  s,dt,r0,mu,alpha,u 
    real*8  fp,c1,c2,c3
    integer iflg 
    ! Local param
    integer NLAG2 
    parameter(NLAG2 = 400)
    integer ncmax,nc 
    real*8  ln
    real*8  x,c0,f,fpp,ds,fdt
    real*8  DANBYB 
    parameter(DANBYB = 1.0d-13)
    ! ... Execution ...
    ! To get close approach needed to take lots of iterations if alpha < 0
    if(alpha.lt.0.d0) then 
        ncmax = NLAG2 
    else
        ncmax = NLAG2 
    endif
    ! Start laguere's method
    ln = 5.d0
    do nc = 0,ncmax 
        x = alpha * s * s 
        call drift_kepu_stumpff(x,c0,c1,c2,c3)
        c1 = c1*s 
        c2 = c2*s*s 
        c3 = c3*s*s*s 
        f = r0*c1 + u*c2 + mu*c3 - dt 
        fp = r0*c0 + u*c1 + mu*c2 
        fpp = (-40.d0*alpha + mu)*c1 + u*c0
        ds = - ln*f/(fp + dsign(1.d0,fp)*sqrt(abs((ln - 1.0)*(ln - 1.0)*fp*fp - (ln - 1.0)*ln*f*fpp)))
        s = s + ds 
        fdt = f/dt 
    ! Quartic convergence
        if(fdt*fdt.lt.DANBYB*DANBYB) then 
            iflg = 0
            return 
        endif
    enddo
    iflg = 2
    return
end subroutine

! *********************************************************************
! drift_kepu_fchk.f90 (subroutine)
! *********************************************************************
! Purpose:
!     Returns the value of the function f of which we are trying to find
!     the root in universal variables.
! 
! The f function is as follows:
!     f(s) = r0*s*c1(alpha*s2) + r0*r0_dot*s2*c2(alpha*s2) + \mu*s3*c3(alpha*s2) - dt
! *********************************************************************
! Input:
!     dt  ==> time step
!     r0  ==> Distance between the center body and particles
!     mu  ==> Reduces mass of system
!     alpha  ==> Twice the binding energy
!     u  ==> Vel. dot radial vector
!     s  ==> Approximate root of function f
! Output:
!     f  ==> function value
! *********************************************************************
! Date written:
!     2023.10.31 by Miao Yuxuan
!     Revision from M. Duncan 
subroutine drift_kepu_fchk(dt,r0,mu,alpha,u,s,f)
    implicit none 
    ! Input/Output param
    real*8  dt,r0,mu,alpha,u,s
    real*8  f 
    ! Local param
    real*8  x,c0,c1,c2,c3
    ! ... Execution ...
    x = alpha * s * s 
    call drift_kepu_stumpff(x,c0,c1,c2,c3)
    c1 = c1 * s 
    c2 = c2 * s * s 
    c3 = c3 * s * s * s 
    f = r0 * c1 + u * c2 + mu * c3 - dt 
    return
end subroutine