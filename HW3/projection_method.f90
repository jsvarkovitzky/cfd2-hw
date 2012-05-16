! ============================================================================
!  Program:     /Users/mandli/Documents/School/grad/ns_course
!  File:        projection_method
!  Created:     2010-05-02
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-05-02 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

module grid_module

    implicit none
    
    ! Grid parameters
    integer :: N_x,N_y
    integer, parameter :: mbc = 1
    double precision, parameter :: gamma = 0.66d0
    double precision :: L_x,L_y,dx,dzeta
    
    ! Grid arrays
    ! i-1,j  |  i,j  |  i+1,j           u <-- x_edge, matches p indices
    ! ---v---|---v---|---v---       
    !        |       |
    !  i-1,j |  i,j  |  i+1,j       
    !    p   u   p   u   p    <--- x_center, y_center, zeta_center, F_center
    !        |       |             
    ! i-1,j-1| i,j-1 | i+1,j-1 
    ! ---v---|---v---|---v--- <--- y_edge, zeta_edge, F_edge
    !        |       |
    double precision, allocatable :: x_edge(:),x_center(:)
    double precision, allocatable :: zeta_edge(:),zeta_center(:)
    double precision, allocatable :: y_edge(:),y_center(:)
    double precision, allocatable :: F_edge(:),F_center(:)
    
contains

    ! Setup and allocate all grid quantities
    subroutine setup_grid()

        implicit none
        
        ! Local
        integer :: i
        
        ! Grid spacing
        dx = (L_x - 0.d0) / (N_x)
        dzeta = (L_y - 0.d0) / (N_y)
        
        ! Grid array allocation
        allocate(x_edge(0:N_x+1),x_center(0:N_x+1))
        allocate(zeta_edge(0:N_y+1),zeta_center(0:N_y+1))
        allocate(y_edge(0:N_y+1),y_center(0:N_y+1))
        allocate(F_edge(0:N_y+1),F_center(0:N_y+1))

        ! Calculate grid arrays
        forall (i=1-mbc:N_x+mbc) 
            x_edge(i) = 0.d0 + i * dx
            x_center(i) = 0.d0 + (i-0.5d0) * dx
        end forall
        forall (i=1-mbc:N_y+mbc)
            zeta_edge(i) = 0.d0 + i * dzeta
            zeta_center(i) = 0.d0 + (i-0.5d0) * dzeta
            F_edge(i) = F(zeta_edge(i))
            F_center(i) = F(zeta_center(i))
        end forall
        ! Really don't need these in this implementation
        y_edge = L_y * (1.d0 - tanh(gamma * (L_y - zeta_edge)) / tanh(gamma*L_y))
        y_center = L_y * (1.d0 - tanh(gamma * (L_y - zeta_center)) / tanh(gamma*L_y))
        
    end subroutine setup_grid
    
    ! Grid stretching function
    double precision pure function F(zeta)
        implicit none
        double precision, intent(in) :: zeta
        F =  tanh(gamma*L_y) / (gamma*L_y) * cosh(gamma*(L_y - zeta))**2
    end function F

    ! Output grid
    subroutine output_grid(frame,t,u,v,p)

        implicit none
        
        ! Arguments
        integer, intent(in) :: frame
        double precision, intent(in) :: t
        double precision, dimension(0:N_x+1,0:N_y+1) :: u,v,p
        
        ! Locals
        integer :: i,j

    
! Open output file and write file header:

        if (t==0) then
        open(unit=70,file='_output/UVP.dat',access='sequential',status='unknown')
!        write(70,*) ' VARIABLES= "x", "y", "u", "v", "p"'
!        write(70,100) t,N_X,N_Y
        endif

 100 FORMAT('ZONE T="t = ',e26.16,'"',' F=POINT, I=',I5,' J=', I5)

        if (t>0) then
!        write(70,100) t,N_X,N_Y
        endif
        
        ! Write out data
        do j=1,N_y
        do i=1,N_x

        ! Reset to zero if exponent is too large
        if (abs(u(i,j)) < 1d-99) u(i,j) = 0.d0
        if (abs(v(i,j)) < 1d-99) v(i,j) = 0.d0
        if (abs(p(i,j)) < 1d-99) p(i,j) = 0.d0

        write(70,"(5e26.16)") x_edge(i),y_edge(j),u(i,j),v(i,j),p(i,j)

        enddo
        enddo

    end subroutine output_grid

end module grid_module

! ============================================================================
!  Program:     /Users/mandli/Documents/School/grad/ns_course
!  File:        main
!  Created:     2010-05-02
!  Author:      Kyle Mandli
! ============================================================================
!      Copyright (C) 2010-05-02 Kyle Mandli <mandli@amath.washington.edu>
!
!  Distributed under the terms of the Berkeley Software Distribution (BSD) 
!  license
!                     http://www.opensource.org/licenses/
! ============================================================================

program main
    
    use grid_module
    
    implicit none
    
    ! ========================================================================
    ! Solver parameters
    integer, parameter :: MAX_ITERATIONS = 1000000
    double precision, parameter :: TOLERANCE = 1.d-6, CFL = 0.8d0
    logical, parameter :: write_star = .false.
    integer :: n_steps

    ! ========================================================================
    ! Physics
    double precision :: U_inf = 1.d0
    double precision :: rho,nu,Re !rho = 1.d0, nu=1.d-3

    ! ========================================================================
    ! Velocity and pressures
    double precision, allocatable :: u(:,:),v(:,:),p(:,:),u_star(:,:),v_star(:,:)
    double precision, allocatable :: u_old(:,:), v_old(:,:)
    
    ! ========================================================================
    ! Locals
    character*20 :: arg
    integer :: i,j,n,m,frame,i_R,j_R
    double precision :: R,t,dt,a
    double precision, allocatable :: Flux_ux(:,:),Flux_uy(:,:),Flux_vy(:,:)
    double precision :: uu_x,uv_y,uv_x,vv_y,u_xx,u_yy,v_xx,v_yy
    double precision, allocatable :: Q(:,:),b(:),cp(:),cm(:)
    ! ========================================================================
    
    ! Get command line arguments
!    if (iargc() /= 6) then
!        print *,"Wrong number of command line arguments, expected 6."
!        print *,"  blasius Nx Ny Lx Ly Nsteps Re"
!        print *,"    Nx - Number of grid points in x"
!        print *,"    Ny - Number of grid points in y"
!        print *,"    Lx - Length of domain in x"
!        print *,"    Ly - Length of domain in y"
!        print *,"    Nsteps - Number of steps in between output"
!        print *,"    Re - Reynolds number of flow"
!        stop
!    else
!        call getarg(1,arg)
!        read (arg,'(I10)') N_x
!        call getarg(2,arg)
!        read (arg,'(I10)') N_y
!        call getarg(3,arg)
!        read (arg,'(e16.8)') L_x
!        call getarg(4,arg)
!        read (arg,'(e16.8)') L_y
!        call getarg(5,arg)
!        read (arg,'(I10)') n_steps
!        call getarg(6,arg)
!        read (arg,'(e16.8)') Re
!    endif


!!!!!!!!!!!!!!!!!!
!!  Parameters: !!
!!!!!!!!!!!!!!!!!!

    N_x=300  !Number of grid points in x-direction
    N_y=60   !Number of grid points in y-direction
    L_x=10.0 !Length of box in x-direction
    L_y=4.0  !Length of box in y-direction
    n_steps=100000 !Interval that u,v and p are printed to UVP.dat
    Re=100.0   !Reynolds number

    print *,"Running blasius with following parameters: "
    print "(a,i3,a,i3,a)"," (N_x,N_y) = (",N_x,",",N_y,")"
    print "(a,e16.8,a,e16.8,a)"," (L_x,L_y) = (",L_x,",",L_y,")"
    print "(a,i4,a)"," Output every ",n_steps," steps"
    print "(a,e16.8)"," Reynold's number = ",Re
    
    ! ========================================================================
    ! Setup grid and arrays
    call setup_grid()
    allocate(Flux_ux(1:N_x+1,0:N_y+1))
    allocate(Flux_uy(0:N_x,0:N_y))
    allocate(Flux_vy(0:N_x+1,0:N_y))
    allocate(Q(1:N_x,1:N_y))
    allocate(b(1:N_y),cp(1:N_y),cm(1:N_y))
    
    ! Calculate matrix coefficients for Poisson solve
    a = 1.d0 / dx**2
    forall (j=1:N_y)
        b(j) = - (2.d0 / dx**2 + F_center(j) / dzeta**2 * (F_edge(j) + F_edge(j-1))) 
        cp(j) = (F_center(j) * F_edge(j)) / dzeta**2
        cm(j) = (F_center(j) * F_edge(j-1)) / dzeta**2
    end forall
        
    ! Velocity and pressure arrays
    allocate(u(1-mbc:N_x+mbc,1-mbc:N_y+mbc),u_star(1-mbc:N_x+mbc,1-mbc:N_y+mbc))
    allocate(v(1-mbc:N_x+mbc,1-mbc:N_y+mbc),v_star(1-mbc:N_x+mbc,1-mbc:N_y+mbc))
    allocate(p(1-mbc:N_x+mbc,1-mbc:N_y+mbc))
    allocate(u_old(1:N_x,1:N_y),v_old(1:N_x,1:N_y))
    
    ! Inital conditions
    u = U_inf
    v = 0.d0
    p = 0.d0
    dt = CFL * dx / (Re * U_inf)
    t = 0.d0
    frame = 0

    nu = 1.d-3
    rho = 1.d0


    ! Output inital condition
    call output_grid(frame,t,u,v,p)
    print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",0," t=",t
    
    ! Open up file to store residual information in
    open(unit=13, file='residual.dat', status="unknown", action="write")
    
    ! ========================================================================
    ! Main algorithm loop
    do n=1,MAX_ITERATIONS
        ! Store old step for convergence test
        u_old = u(1:N_x,1:N_y)
        v_old = v(1:N_x,1:N_y)

        ! ====================================================================
        ! Apply BC  V = x, U = o, P = +
        call bc(u,v,U_inf)
        
        ! ====================================================================
        ! Step 1: Update velocity to intermediate step
        ! Calculate fluxes at each boundary
        !
        !  U-Fluxes
        !                                   |   u   |
        !         FUY_i,j                i,j| i,j+1 |i+1,j
        ! -------|---o---|-------    -------v-------v-------
        !        |       |            i-1,j |  i,j  |  i+1,j
        !    FUX_i,j    FUX_i+1,j       u   |   u   |   u
        !        |       |                  |       |
        ! -------|---o---|-------    -------v-------v-------
        !         FUY_i,j-1            i,j-1| i,j-1 | i+1,j-1
        !                                   |   u   |  
        !  
        !  V-Fluxes
        !                                   |   v   |
        !         FVY_i,j            i-1,j+1| i,j+1 |i,j+1
        ! -------|---o---|-------    -------u---o---u-------
        !        |       |            i-1,j |  i,j  |  i+1,j
        !    FUY_i-1,j  FUY_i,j         v   o   v   o   v
        !        |       |                  |       |
        ! -------|---o---|-------    -------u---o---u-------
        !         FVY_i,j-1            i-1,j| i,j-1 |i,j
        !                                   |   v   |    
        forall (i=1:N_x+1,j=0:N_y+1) Flux_ux(i,j) = (u(i-1,j)+u(i,j))**2 / 4.d0
        forall (i=0:N_x,j=0:N_y) Flux_uy(i,j) = (u(i,j)+u(i,j+1)) * (v(i+1,j)+v(i,j)) / 4.d0
        forall (i=0:N_x+1,j=0:N_y) Flux_vy(i,j) = (v(i,j+1) + v(i,j))**2 / 4.d0
    
        do j=1,N_y
            do i=1,N_x

!********************************************
!***  ADD RHS TO THE FOLLOWING FORMULAS:  ***
!********************************************

!***  Source of formulas are CDS  derivatives from Lecture 5 pg 8 notes ***
                ! Advective terms
                uu_x = 1.d0/dx*(Flux_ux(i+1,j)-Flux_ux(i,j))
                uv_y = F_center(j)/dzeta*(Flux_uy(i,j)-Flux_uy(i,j-1))
                
                uv_x = 1.d0/dx*(Flux_uy(i,j)-Flux_uy(i-1,j))
                vv_y = F_edge(j)/dzeta*(Flux_vy(i,j)-Flux_vy(i,j-1))
                
                ! Diffusive terms

                !!!CHECK THAT FLUXes USED CORRECTLY AND IF F() SHOULD BE USED INSTEAD!!!

                u_xx = 1.d0/dx**2*(u(i+1,j)-2*u(i,j)+u(i-1,j))
                u_yy = F_center(j)/dzeta**2*(F_edge(j)*(u(i,j+1)-u(i,j))-F_edge(j-1)*(u(i,j)-u(i,j-1)))
                v_xx = 1.d0/dx**2*(v(i+1,j)-2*v(i,j)+v(i-1,j))
                v_yy = F_edge(j)/dzeta**2*(F_center(j+1)*(v(i,j+1)-v(i,j))-F_center(j)*(v(i,j)-v(i,j-1)))
                
                ! Update to u* and v* value
                u_star(i,j) = u(i,j) + dt*(-(uu_x+uv_y)+nu*(u_xx+u_yy))
                v_star(i,j) = v(i,j) + dt*(-(vv_y+uv_x)+nu*(v_xx+v_yy))
            enddo
        enddo
        
        ! Debug, save out u_star,v_star,p
        if (write_star) then
            frame = frame + 1
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
            call output_grid(frame,t,u_star,v_star,p)
        endif
        
        ! ====================================================================
        ! Step 2: Solve projection poisson problem
        call bc(u_star,v_star,U_inf)
        forall(i=1:N_x,j=1:N_y)
            Q(i,j) = 1.d0/dt * ((u_star(i,j)-u_star(i-1,j)) / dx + (v_star(i,j)-v_star(i,j-1)) / dzeta * F_center(j))
        end forall
        ! Solve poisson problem
        call solve_poisson(p,Q,a,b,cm,cp)
        
        ! ====================================================================
        ! Step 3: Update velocities to end time
        !**********************************************
        !***  ADD RHS TO VELOCITY UPDATE FORLUMAS:  ***
        !**********************************************
        forall (i=1:N_x,j=1:N_y)
            u(i,j) = u_star(i,j) - dt/(rho*dx)*(p(i+1,j)-p(i,j))
            v(i,j) = v_star(i,j) - dt/(rho*dzeta)*(p(i,j+1)-p(i,j))*F_edge(j)
        end forall
        
        ! ====================================================================
        ! Step 4: Check convergence
        R = 0.d0
        do j=1,N_y
            do i=1,N_x
                if (R < abs(u(i,j)-u_old(i,j)) .or. R < abs(v(i,j)-v_old(i,j))) then
                    R = max(R,abs(u(i,j)-u_old(i,j)),abs(v(i,j)-v_old(i,j)))
                    i_R = i
                    j_R = j
                endif
            enddo
        enddo
        
        ! Finish up loop
        print "(a,i4,a,i3,a,i3,a,e16.8)","Loop ",n,": (",i_R,",",j_R,") - R = ",R   
        write (13,"(i4,i4,i4,e16.8)") n,i_R,j_R,R
        ! Write out u,v,p every n_steps 
        if (mod(n,n_steps) == 0) then
            frame = frame + 1
            call output_grid(frame,t,u,v,p)
            print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
        endif
        ! Check tolerance
        if (R < TOLERANCE) then
            print *, "Convergence reached, R = ",R
    call output_grid(frame,t,u,v,p)
    print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
            exit
        endif
        ! We did not reach our tolerance, iterate again
        t = t + dt
    enddo
    if (R > TOLERANCE) then
        print "(a,e16.8)","Convergence was never reached, R = ", R
    endif
!    call output_grid(frame,t,u,v,p)
!    print "(a,i3,a,i4,a,e16.8)","Writing frame ",frame," during step n=",n," t=",t
    close(13)
end program main


! ============================================================================
! Boundary Conditions                                                        
! ============================================================================

subroutine bc(u,v,U_inf)

    use grid_module

    implicit none
    
    ! Input/Output arguments
    double precision, intent(in) :: U_inf
    double precision, intent(inout) :: u(1-mbc:N_x+mbc,1-mbc:N_y+mbc)
    double precision, intent(inout) :: v(1-mbc:N_x+mbc,1-mbc:N_y+mbc)

    ! Locals
    integer :: i,j
    
    ! Lower Boundary         Upper Boundary     
    !  i,1  i,1    |          
    !   x    +     x i,1     |   o   | 
    !   |          |         |i,N_y+1|
    ! -------o--------       x   +   x 
    !  |    i,0    |         |       |
    !  x     +     x        -|---o---|-
    !  |    i,0   i,0        | i,N_y |
    !  |           |         x   +   x   
    !                               i,N_y
    ! Wall                  Free shear
    ! u(i,0) = -u(i,1)      u(i,N_y+1) = u(i,N_y)
    ! v(i,0) = 0.d0         v(i,N_y+1) = v(i,N_y)
    ! p(i,0) = p(i,1)       p(i,N_y+1) = 0.d0
    forall(i=N_x/2:N_x+1)
        u(i,0) = -u(i,1)
        v(i,0) = 0.d0     
        u(i,N_y+1) = u(i,N_y)
        v(i,N_y+1) = v(i,N_y)
    end forall

    !Adding first half of domain with free shear BCs
    forall(i=0:N_x/2)
        u(i,0) = u(i,1)
        v(i,0) = v(i,1)     
        u(i,N_y+1) = u(i,N_y)
        v(i,N_y+1) = v(i,N_y)
    end forall

    ! Left Boundaries
    !   x 0,j  |      x 1,j             
    !  0,j     |                        u(0,j) = U_inf
    !   +      o 0,j  + 1,j   o 1,j     v(0,j) = 0.d0
    !          |                        p(0,j) = p(1,j)  (P_x = 0)
    ! Right Boundaries (outflow)
    !     x N_x,j   |         x N_x+1,j             u(N_x+1,j) = 2 u(N_x,j) - u(N_x-1,j)
    !               |                               v(N_x+1,j) = 2 v(N_x,j) - v(N_x-1,j)
    !     + N_x,j   o N_x,j   + N_x+1,j  o N_x+1,j  p(N_x+1,j) = p(N_x,j)
    !               |
    forall(j=1:N_y+1)
        u(0,j) = U_inf
        v(0,j) = 0.d0
!         u(N_x+1,j) = 2.d0*u(N_x,j) - u(N_x-1,j)
        u(N_x+1,j) = u(N_x,j)
        v(N_x+1,j) = v(N_x,j)
!         v(N_x+1,j) = 2.d0*v(N_x,j) - v(N_x-1,j)
    end forall

end subroutine bc


! ============================================================================
! Solve the poisson problem
!  laplace(P)_ij = Q_ij
! ============================================================================

subroutine solve_poisson(P,Q,a,b,cm,cp)

    use grid_module

    implicit none
    
    ! Input
    double precision, intent(in) :: Q(1:N_x,1:N_y),a
    double precision, intent(in) :: b(1:N_y),cm(1:N_y),cp(1:N_y)
    double precision, intent(inout) :: P(1-mbc:N_x+mbc,1-mbc:N_y+mbc)

    ! Solver parameters
    logical, parameter :: verbose = .false.
    integer, parameter :: MAX_ITERATIONS = 100
    double precision, parameter :: TOLERANCE = 10.d-6
    double precision, parameter :: w = 1.6d0 ! 1 (GS) < w < 2

    ! Local variables
    integer :: i,j,n
    double precision :: R,P_old
    
    do n=1,MAX_ITERATIONS
        R = 0.d0
        ! Boundary conditions
        forall (j=0:N_y+1)
            P(0,j) = P(1,j)        ! Left
            P(N_x+1,j) = P(N_x,j)  ! Right
        end forall
        forall (i=0:N_x+1)
            P(i,0) = P(i,1)        ! Bottom wall
            P(i,N_y+1) = 0.d0      ! Free stream
        end forall
        do j=1,N_y
            do i=1,N_x
                ! Update p
                P_old = P(i,j)
                P(i,j) = (Q(i,j) - (a * (P(i+1,j) + P(i-1,j)) + cp(j)*P(i,j+1) + cm(j)*P(i,j-1))) / b(j)
                ! relaxation
                P(i,j) = w * P(i,j) + (1-w) * P_old
                ! Calculate Residual
                R = max(R,abs(P(i,j)-P_old))
            enddo
        enddo
        
        ! Print out convergence and exit
        if (verbose) then
            print "(a,i5,a,e16.8,a)", "(",n,",",R,")"
        endif
        if (R < TOLERANCE) then
            exit
        endif
    enddo
    ! Check to see if we converged and quit if we did not
    if (n > MAX_ITERATIONS) then
        print *,"Poisson solver did not converge, exiting!"
        print "(a,e16.8)","R = ",R
!        stop
    else
        print "(a,i4,a,e16.8)", "Solver converged in ",n," steps: R = ",R
        ! Boundary conditions
        forall (j=0:N_y+1)
            P(0,j) = P(1,j)        ! Left
            P(N_x+1,j) = P(N_x,j)  ! Right
        end forall
        forall (i=0:N_x+1)
            P(i,0) = P(i,1)        ! Bottom wall
            P(i,N_y+1) = 0.d0      ! Free stream
        end forall
    endif

end subroutine solve_poisson

