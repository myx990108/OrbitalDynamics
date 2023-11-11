program main 
    implicit none 
    real*8  x,y,z,vx,vy,vz
    real*8  a,e,inc,capom,omega,capm0,capm1,capm2,en
    real*8  dt,t0
    integer iflg
    real*8  PI,TWOPI
    parameter(PI = 3.1415926535897932384d0)
    parameter(TWOPI = 2.d0*PI)
    a = 1.43d0 
    e = 0.999999999999999d0
    inc = 0.12d0
    capom = 0.32d0 
    omega = 0.11d0
    capm0 = 0.23d0 
    t0 = 0.d0
    write(6,"(a,f8.3)") "Initial orbit element at t0 = ", t0
    write(6,"(5(f12.8),f20.16)") a,e,inc,capom,omega,capm0

    en = sqrt(1.d0/a/a/a)
    dt = 0.1d0
    call orbel_el2xv(1.d0,a,e,inc,capom,omega,capm0,x,y,z,vx,vy,vz)
    call drift_one(1.d0,x,y,z,vx,vy,vz,dt,iflg)
    call orbel_xv2el(1.d0,x,y,z,vx,vy,vz,a,e,inc,capom,omega,capm1)
    write(6,"(a)") "Orbit element after dt computed by danby is: "
    write(6,"(5(f12.8),f20.16)") a,e,inc,capom,omega,capm1

    capm2 = capm0 + en*dt 
    capm2 = capm2 - int(capm2/TWOPI)*TWOPI
    if(capm2.lt.0.d0) capm2 = capm2 + TWOPI
    write(6,"(a)") "Orbit element after dt theoretically is: "
    write(6,"(5(f12.8),f20.16)") a,e,inc,capom,omega,capm2

    write(6,"(a,1p1e12.5)") "The difference is: ",abs(capm2 - capm1)
end program
