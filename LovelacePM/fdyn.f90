!solve aircraft positions as they vary in time
subroutine tstep_solve(nstep, dt, rho, U0, g, Sref, cref, bref, perturbation, onboard, inertia, m, &
coeffs, derivs, external_history, alpha_history, beta_history, euler_history)
    !for pertutbation, here we use the notation employed by Drela in Flight Viechle Aerodynamics
    !equations are also directly obtained from there
    integer, intent(IN) :: nstep
    real(8), intent(IN) :: dt, rho, U0, g, Sref, cref, bref, perturbation(1:12), onboard(1:3), &
    inertia(1:3, 1:3), m, coeffs(1:6), derivs(1:5, 1:6)
    !notation for coefs: (CX, CY, CZ, Cl, Cm, Cn)
    !notation for derivs: (CX, CY, CZ, Cl, Cm, Cn), first row: a derivatives; second, b derivatives; third, p; 4, q; 5, r
    real(8), intent(OUT) :: external_history(1:nstep, 1:3), alpha_history(1:nstep), beta_history(1:nstep), &
    euler_history(1:nstep, 1:3)

    integer :: i, j

    real(8) :: Uinf, xe, ye, ze, phi, theta, psi, u, v, w, p, q, r, &
    forces(1:3), moments(1:3), invI(1:3, 1:3), det_l, a, b, mtemp(1:3)

    xe=perturbation(1)
    ye=perturbation(2)
    ze=perturbation(3)
    phi=perturbation(4)
    theta=perturbation(5)
    psi=perturbation(6)
    u=perturbation(7)+U0
    v=perturbation(8)
    w=perturbation(9)
    p=perturbation(10)*2*U0/bref
    q=perturbation(11)*2*U0/cref
    r=perturbation(12)*2*U0/bref !denormalizing angular velocities

    invI(1:3, 1)=inertia(1:3, 1)
    invI(1:3, 2)=inertia(1:3, 2)
    invI(1:3, 3)=inertia(1:3, 3)
    call invert3(invI, det_l)

    do i=1, nstep
        Uinf=sqrt(u**2+v**2+w**2)
        a=atan(w/u)
        b=atan(v/sqrt(u**2+w**2))

        do j=1, 3
            forces(j)=(coeffs(j)+derivs(1, j)*a+derivs(2, j)*b+derivs(3, j)*p*bref/(2*Uinf)+&
            derivs(4, j)*q*cref/(2*Uinf)+derivs(5, j)*r*bref/(2*Uinf))*Sref*Uinf**2*rho/2
        end do
        moments(1)=(coeffs(4)+derivs(1, 4)*a+derivs(2, 4)*b+derivs(3, 4)*p*bref/(2*Uinf)+&
        derivs(4, 4)*q*cref/(2*Uinf)+derivs(5, 4)*r*bref/(2*Uinf))*Sref*bref*Uinf**2*rho/2
        moments(2)=(coeffs(5)+derivs(1, 5)*a+derivs(2, 5)*b+derivs(3, 5)*p*bref/(2*Uinf)+&
        derivs(4, 5)*q*cref/(2*Uinf)+derivs(5, 5)*r*bref/(2*Uinf))*Sref*cref*Uinf**2*rho/2
        moments(3)=(coeffs(6)+derivs(1, 6)*a+derivs(2, 6)*b+derivs(3, 6)*p*bref/(2*Uinf)+&
        derivs(4, 6)*q*cref/(2*Uinf)+derivs(5, 6)*r*bref/(2*Uinf))*Sref*bref*Uinf**2*rho/2

        xe=xe+dt*(cos(theta)*cos(psi)*u+(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))*v+&
        (cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))*w)
        ye=ye+dt*(cos(theta)*sin(psi)*u+(sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi))*v+&
        (cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi))*w)
        ze=ze+dt*(-sin(theta)*u+sin(phi)*cos(theta)*v+cos(phi)*cos(theta)*w)

        u=u+dt*(forces(1)/m-g*sin(theta)-q*w+r*v)
        v=v+dt*(forces(2)/m+g*sin(phi)*cos(theta)-r*u+p*w)
        w=w+dt*(forces(3)/m+g*cos(phi)*cos(theta)-p*v+q*u)

        mtemp(1)=moments(1)-(inertia(3, 3)-inertia(2, 2))*q*r-inertia(2, 3)*(q**2-r**2)-&
        inertia(1, 3)*p*q+inertia(1, 2)*p*r-onboard(3)*q+onboard(2)*r
        mtemp(2)=moments(2)-(inertia(1, 1)-inertia(3, 3))*r*p-inertia(1, 3)*(r**2-p**2)-&
        inertia(1, 2)*q*r+inertia(2, 3)*q*p-onboard(1)*r+onboard(3)*p
        mtemp(3)=moments(3)-(inertia(3, 3)-inertia(1, 1))*p*q-inertia(1, 2)*(p**2-q**2)-&
        inertia(2, 3)*r*p+inertia(1, 3)*r*q-onboard(2)*p+onboard(1)*q

        mtemp=matmul(invI, mtemp)

        p=p+dt*mtemp(1)
        q=q+dt*mtemp(2)
        r=r+dt*mtemp(3)

        phi=phi+dt*(p+q*sin(phi)*tan(theta)+r*cos(phi)*tan(theta))
        theta=theta+dt*(q*cos(phi)-r*sin(phi))
        psi=psi+dt*(q*sin(phi)/cos(theta)+r*cos(phi)/cos(theta))

        alpha_history(i)=a
        beta_history(i)=b
        euler_history(i, 1)=phi
        euler_history(i, 2)=theta
        euler_history(i, 3)=psi
        external_history(i, 1)=xe
        external_history(i, 2)=ye
        external_history(i, 3)=ze
    end do
end subroutine tstep_solve

subroutine invert3(a, det_l) !simple matrix inverter
    double precision, intent(inout) :: a (1:3, 1:3)
    real(8), intent(out) :: det_l
    real(8) :: b(1:3, 1:3)

    det_l = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) &
        -a(1,2)*(a(2,1)*a(3,3)-a(2,3)*a(3,1)) &
        +a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    
    b=a

    a(1,1) =  b(2,2)*b(3,3) - b(2,3)*b(3,2)
    a(2,1) =  b(2,3)*b(3,1) - b(2,1)*b(3,3)
    a(3,1) =  b(2,1)*b(3,2) - b(2,2)*b(3,1)

    a(1,2) =  b(1,3)*b(3,2) - b(1,2)*b(3,3)
    a(2,2) =  b(1,1)*b(3,3) - b(1,3)*b(3,1)
    a(3,2) =  b(1,2)*b(3,1) - b(1,1)*b(3,2)

    a(1,3) =  b(1,2)*b(2,3) - b(1,3)*b(2,2)
    a(2,3) =  b(1,3)*b(2,1) - b(1,1)*b(2,3)
    a(3,3) =  b(1,1)*b(2,2) - b(1,2)*b(2,1)

    a=a/det_l
end