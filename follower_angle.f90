program test
  use constants
  implicit none
  integer :: i,j,k,iter, counter, store, per
  integer, parameter :: N=64;
  integer, parameter :: alpha= 1;
  real(kind=8) :: delt, domain, x0(alpha), v
  real(kind=8) :: psiprev(alpha,N), psinew(alpha,N)
  real(kind=8) :: ten(alpha,N), F(alpha,N)
  real(kind=8) :: inter(alpha,2,N), ufiber(2,N)
  real(kind=8) :: X(alpha,N), Y(alpha,N), Frep(alpha,2,N)
  real(kind=8) :: fx(alpha,N), fy(alpha,N)
  character*40 filename_velocity
  character*40 filename_position

  ! Initialize few things and set time step
  inter = 0.d0;  ten = 0.d0; ufiber = 0.d0;
  delt  = 4.d-7; store = 50;
  ! Initial condition for the problem
  call initial(N,alpha,domain,x0,psiprev)
  ! Forward euler step to begin with
  do k = 1,alpha
    call tension(N,psiprev(k,:),ufiber,Frep(k,:,:),ten(k,:))
    call feul(N,psiprev(k,:),ten(k,:),ufiber,F(k,:),psinew(k,:))
  end do

  do counter = 1,10

    write(filename_position,'(A9,I3.3,A4)')'position_',int(counter),'.txt'
    open(2,file=filename_position)

    write(filename_velocity,'(A9,I3.3,A4)')'velocity_',int(counter),'.txt'
    open(3,file=filename_velocity)

  do iter = 1,25000

  ! Compute replusive forces
    per = 1;
    call position (N,alpha,psinew,x0,X,Y)
    call repelforce(N,alpha,per,domain,psinew,X,Y,Frep)
  ! Compute hydrodynamic interaction
  ! if (iter > 500 .or. counter > 1) then
  !   per = 1;
  !   call wallinterper(N,alpha,per,domain,x0,psinew,ten,inter)
  ! end if

    !$OMP PARALLEL PRIVATE(ufiber)
    !$OMP DO
    do k = 1,alpha
      ufiber = inter(k,:,:);
      ! Solve for tension
      call tension(N,psinew(k,:),ufiber,Frep(k,:,:),ten(k,:))
      ! March forward in time
      call beul(N,delt,ten(k,:),ufiber,Frep(k,:,:), &
              & F(k,:),psiprev(k,:),psinew(k,:))
    end do
    !$OMP END DO
    !$OMP END PARALLEL

! Write data to file
    if (mod(iter,store)==0) then
      print*, counter, iter
      do j=1,N
        write(2,*)  j, (real(psinew(k,j)), k = 1,alpha)
        write(3,*)  j, (real(inter(k,1,j)), real(inter(k,2,j)), k = 1,alpha)
     end do
   end if

! Check for code stability
   v = maxval(abs(psinew))
   if (isnan(v)) then
      print*, "NAN eroor"
      stop
   end if
 end do
 close(2)
end do

end program test
!*********************************************************************
! INTIAIL CONDITION ON THE ANGLE
!*********************************************************************
subroutine initial(N,alpha,dom,x0,phi)
  use constants
  implicit none
  integer, intent(in) :: N,alpha
  real(kind=8), intent(out) :: x0(alpha), phi(alpha,N), dom
  real(kind=8) :: h, dis,s
  integer :: i,j,k, idum

! Define spacing and random seed
  h = 1.d0/(N-1.d0);
  ! Create the initian angle of the filament
  do k = 1,alpha
    do i = 1,N
      s = (i-1.d0)*h
      if (mod(k,2) == 0) then
        phi(k,i) = pi/2.d0+ 2*k*h*s;
      else
        phi(k,i) = pi/2.d0+ 2*k*h*s;
      end if
    end do
  end do

!! Create the x base point for the filaments
  dis = 0.2d0;
  do k = 1,alpha
    x0(k) = (k-1.d0)*dis
  end do
  x0  = x0 + dis/2.d0
  dom = alpha*dis;

end subroutine initial
!*********************************************************************
! Position in physical space (x and y coordinates)
!*********************************************************************
subroutine position(N,alpha,phi,x0,X,Y)
  use constants
  implicit none
  integer, intent(in) :: N,alpha
  real(kind=8), intent(in) :: x0(alpha), phi(alpha,N) ! Initial angle at the base
  real(kind=8), intent(out) :: X(alpha,N), Y(alpha,N)
  real(kind=8) :: h, out, tx(N), ty(N)
  integer :: i,k

  h = 1.d0/(N-1.d0);
  do k = 1,alpha
    do i = 1,N
  ! X co-ordinate of the filament
      tx(i) = cos(phi(k,i)); ! tangent vector in the x dir
      call trapz(i,out,h,tx(1:i))
      X(k,i)= out + x0(k); ! x-cord

  ! Y co-ordinate of the filament
      ty(i) = sin(phi(k,i)); ! tangent vector in the y dir
      call trapz(i,out,h,ty(1:i))
      Y(k,i)= out; ! y-cord
    end do
  end do

end subroutine position
!***********************************************************************
! Solve for the tension in the problem
!*************************************************************************
subroutine tension(N,psi,inter,Frep,ten)
  use constants
  implicit none
  integer, intent(in) :: N
  real(kind=8), intent(in) :: psi(N), inter(2,N), Frep(2,N)
  real(kind=8), intent(out):: ten(N)
  real(kind=8) :: Ut(N), Un(N), Uts(N)
  real(kind=8) :: ps1(N), ps2(N), ps3(N)
  real(kind=8) :: T1(N-2), T2(N-2), T3(N-2), b(N-2), Ttemp(N-2)
  real(kind=8) :: h, h2, cpll, cperp, beta
  real(kind=8) :: Ft(N), Fn(N), Fts(N)
  integer :: i, indx(N)
  !real(kind=8) :: T1(N-2), T2(N-2), T3(N-2), b(N-2), Ttemp(N-2)

  ! Define friction coefficient and spacing
  h = 1.d0/(N-1.d0); h2 = h**2;
  cpll = -1.d0/(2.d0*c); cperp = 1.d0/(2.d0-c);
  beta = cperp/cpll;
  ! Define the interaction velocities
  Ut(:) = inter(1,:); call gder(N,1,Ut,Uts)
  Un(:) = inter(2,:);
  ! The repulsive force and its derivatives
  Ft(:) = Frep(1,:); call gder(N,1,Ft,Fts)
  Fn(:) = Frep(2,:);
  ! Take derivatives of the tangent angle
  call gder(N,1,psi,ps1); call gder(N,2,psi,ps2); call gder(N,3,psi,ps3);
  ! Form the matrix for tension
  b = 0.d0;
  ! Form the interior of the matrix
  do i = 2,N-1
    T1(i-1) = 1.d0/h2;
    T2(i-1) =-2.d0/h2-ps1(i)**2/beta;
    T3(i-1) = 1.d0/h2;

    b(i-1)  = (ps1(i)**4)/beta - 3.d0*ps2(i)**2 - (3.d0+1.d0/beta)*ps1(i)*ps3(i);
    b(i-1)  = b(i-1) + cpll*(ps1(i)*Un(i)-Uts(i));
    b(i-1)  = b(i-1) - (Fts(i)*beta-Fn(i)*ps1(i))/beta;
  end do
  ! Boundary condition at s=0; dT/ds = (sigma-3*ps2*ps1)
  T2(1) = T2(1) + (4.d0/3.d0)/h2;
  T3(1) = T3(1) - (1.d0/3.d0)/h2;
  b(1)  = b(1)  + (2.d0*h/3.d0)*(sigma-3.d0*ps1(1)*ps2(1))/h2

  ! Solve using thomas algorithm decomposition
  call tridag(T1,T2,T3,b,Ttemp,N-2)
  ! Substitte
  ten(2:N-1) = Ttemp;
  ten(1)     = (4.d0*ten(2)-ten(3))/(3.d0) - (2.d0*h/3.d0)*(sigma-3.d0*ps1(1)*ps2(1))
  ten(N)     = 0.d0;

end subroutine tension
!***********************************************************************
! Forward Euler marching for the first time step
!*************************************************************************
subroutine feul(N,psi,ten,inter,Fprev,psinew)
  use constants
  implicit none
  integer, intent(in) :: N
  real(kind=8), intent(in) :: psi(N), ten(N), inter(2,N)
  real(kind=8), intent(out):: Fprev(N),psinew(N)
  real(kind=8) :: Ut(N), Un(N), Uns(N)
  real(kind=8) :: ps1(N), ps2(N), ps3(N), Ts(N)
  real(kind=8) :: cpll, cperp, beta, factor
  real(kind=8) :: h,h2,h3,h4
  real(kind=8) :: Amat(N,N), b(N), d
  real(kind=8) :: comp(N,9), AL(N,9)
  integer :: i, j, indx(N),k

  ! Define friction coefficient and spacing
  h = 1.d0/(N-1.d0); h2 = h**2; h3 = h**3; h4 = h**4;
  cpll = -1.d0/(2.d0*c); cperp = 1.d0/(2.d0-c);
  beta = cperp/cpll;
  ! Define the interaction velocities
  Ut(:) = inter(1,:); Un(:) = inter(2,:); call gder(N,1,Un,Uns)
  factor = (2.d0-c)*dteu; Fprev = 0.d0;

  ! Take derivatives of the tangent angle and tension
  call gder(N,1,psi,ps1); call gder(N,2,psi,ps2)
  call gder(N,1,ten,Ts);  call gder(N,3,psi,ps3)

  ! First form the interior of the matrix and the RHS
  b = 0.d0; Amat = 0.d0;
  do i = 3,N-2
    b(i) = (Ts(i)*ps1(i)+3.d0*ps2(i)*ps1(i)**2)*(1.d0+beta) - beta*sigma*ps1(i)
    b(i) = b(i) + cperp*(Ut(i)*ps1(i)+Uns(i));
    Fprev(i) = b(i);
    b(i) = b(i)*factor + psi(i);


    Amat(i,i-2) = 1.d0*factor/h4;
    Amat(i,i-1) =-4.d0*factor/h4 -      ten(i)*factor/h2;
    Amat(i,i)   = 6.d0*factor/h4 + 2.d0*ten(i)*factor/h2 + 1.d0;
    Amat(i,i+1) =-4.d0*factor/h4 -      ten(i)*factor/h2;
    Amat(i,i+2) = 1.d0*factor/h4;
  end do
  ! Apply the boundary conditions
  ! EBC3
  b(1)       = -ps1(1)**3;
  Amat(1,1)  = (-3.d0*ten(1)/(2.d0*h)-(-2.5d0/h3));
  Amat(1,2)  = ( 4.d0*ten(1)/(2.d0*h)-( 9.0d0/h3));
  Amat(1,3)  = (-1.d0*ten(1)/(2.d0*h)-(-12.d0/h3));
  Amat(1,4)  = (0.d0                 -( 7.0d0/h3));
  Amat(1,5)  = (0.d0                 -(-1.5d0/h3));
  ! EBC5
  Amat(2,1)  = 1.d0; b(2)  = pi/2.d0;
  ! EBC6
  Amat(N-1,N)= 3.d0; Amat(N-1,N-1)= -4.d0; Amat(N-1,N-2) = 1.d0
  b(N-1) = 0.d0;
  ! EBC4
  Amat(N,N)   = 2.d0; Amat(N,N-1) =-5.d0; Amat(N,N-2) = 4.d0; Amat(N,N-3) =-1.d0;
  b(N) = 0.d0;

  ! Solve the band diagonal system
    ! call ludcmp(Amat,N,N,indx,d)
    ! call lubksb(Amat,N,N,indx,b)

    ! Solve the band diagonal system
       do i=1,N
          do j= 1,9
             k= i-5+j
             if(k.ge.1.and.k.le.N) comp(i,j)=Amat(i,k)
          enddo
       enddo
      call bandec(comp,N,4,4,N,9,AL,9,indx,d)
      call banbks(comp,N,4,4,N,9,AL,9,indx,b)

    ! Substitute solution as the new angle
  !psitemp = b;

 ! Assign the solution to new tangent angle
    psinew = b;

end subroutine feul
!***********************************************************************
! Forward Euler marching for the first time step
!*************************************************************************
subroutine beul(N,delt,ten,inter,Frep,Fprev,psiprev,psinew)
  use constants
  integer, intent(in) :: N
  real(kind=8), intent(in) :: delt, ten(N), inter(2,N), Frep(2,N)
  real(kind=8), intent(inout):: psiprev(N), psinew(N), Fprev(N)
  real(kind=8) :: Fnew(N), psP1(N)
  real(kind=8) :: Ut(N), Un(N), Uns(N)
  real(kind=8) :: ps1(N), ps2(N), ps3(N), Ts(N)
  real(kind=8) :: cpll, cperp, beta
  real(kind=8) :: Ft(N), Fn(N), Fns(N)
  real(kind=8) :: h,h2,h3,h4
  real(kind=8) :: Amat(N,N), b(N), d, factor
  real(kind=8) :: comp(N,9), AL(N,9)
  integer :: i, indx(N)

  ! Define friction coefficient and spacing
  h = 1.d0/(N-1.d0); h2 = h**2; h3 = h**3; h4 = h**4;
  cpll = -1.d0/(2.d0*c); cperp = 1.d0/(2.d0-c);
  beta = cpll/cperp;
  ! Define the interaction velocities
  Ut(:) = inter(1,:);
  Un(:) = inter(2,:); call gder(N,1,Un,Uns)
  ! The repulsive force and its derivatives
  Ft(:) = Frep(1,:);
  Fn(:) = Frep(2,:); call gder(N,1,Fn,Fns)
  ! Take derivatives of the tangent angle and tension
  call gder(N,1,psinew,ps1);   call gder(N,2,psinew,ps2);
  call gder(N,1,psiprev,psP1); call gder(N,1,ten,Ts);
  factor = 2.d0*(2.d0-c)*delt; Fnew = 0.d0;

  ! First form the interior of the matrix and the RHS
  b = 0.d0; Amat = 0.d0;
  do i = 3,N-2
    Fnew(i) = (Ts(i)*ps1(i)+3.d0*ps2(i)*ps1(i)**2)*(1.d0+beta) - beta*sigma*ps1(i)
    Fnew(i) = Fnew(i) + cperp*(Ut(i)*ps1(i)+Uns(i));
    Fnew(i) = Fnew(i) + (Fns(i)+beta*Ft(i)*ps1(i));

    b(i)    = 2.d0*Fnew(i)-Fprev(i);
    b(i)    = b(i)*factor + (4.d0*psinew(i)-psiprev(i));

    Fprev(i) = Fnew(i);


    Amat(i,i-2) = 1.d0*factor/h4;
    Amat(i,i-1) =-4.d0*factor/h4 -      ten(i)*factor/h2;
    Amat(i,i)   = 6.d0*factor/h4 + 2.d0*ten(i)*factor/h2 + 3.d0;
    Amat(i,i+1) =-4.d0*factor/h4 -      ten(i)*factor/h2;
    Amat(i,i+2) = 1.d0*factor/h4;
  end do
  ! Apply the boundary conditions
  ! EBC3
  b(1)       = -(2.d0*ps1(1)**3-psP1(1)**3);
  Amat(1,1)  = (-3.d0*ten(1)/(2.d0*h)-(-2.5d0/h3));
  Amat(1,2)  = ( 4.d0*ten(1)/(2.d0*h)-( 9.0d0/h3));
  Amat(1,3)  = (-1.d0*ten(1)/(2.d0*h)-(-12.d0/h3));
  Amat(1,4)  = (0.d0                 -( 7.0d0/h3));
  Amat(1,5)  = (0.d0                 -(-1.5d0/h3));
  ! EBC5
  !Amat(2,1)  = -3.d0/(2.d0*h); Amat(2,3) = 4.d0/(2.d0*h); Amat(2,3) =-1.d0/(2.d0*h);
  Amat(2,1)  = Amat(2,1) - kp; b(2)  = -kp*pi/2.d0;
  ! EBC6
  Amat(N-1,N)= 3.d0; Amat(N-1,N-1)= -4.d0; Amat(N-1,N-2) = 1.d0
  b(N-1) = 0.d0;
  ! EBC4
  Amat(N,N)   = 2.d0; Amat(N,N-1) =-5.d0; Amat(N,N-2) = 4.d0; Amat(N,N-3) =-1.d0;
  b(N) = 0.d0;


  ! Solve the band diagonal system
     do i=1,N
        do j= 1,9
           k= i-5+j
           if(k.ge.1.and.k.le.N) comp(i,j)=Amat(i,k)
        enddo
     enddo
    call bandec(comp,N,4,4,N,9,AL,9,indx,d)
    call banbks(comp,N,4,4,N,9,AL,9,indx,b)

  ! Substitute solution as the new angle
  psiprev = psinew; psinew = b;

end subroutine beul
!******************************************************************************
! Force acting on the filament
!*****************************************************************************
subroutine force(N,alpha,psi,ten,fx,fy)
  use constants;
  implicit none
  integer, intent(in) :: N,alpha
  real(kind=8), intent(in) :: psi(alpha,N), ten(alpha,N)
  real(kind=8), intent(out):: fx(alpha,N), fy(alpha,N)
  real(kind=8) :: Ts(N), ps1(N), ps2(N), ps3(N), fT(N),fN(N)
  integer :: j,k

  do k = 1,alpha
     call gder(N,1,ten(k,:) ,Ts )
     call gder(N,1,psi(k,:) ,ps1)
     call gder(N,2,psi(k,:) ,ps2)
     call gder(N,3,psi(k,:) ,ps3)
    do j = 1,N   ! At the jth point
      ! Obtain the tangential and normal force this is actually force per unit length
      fT(j) =  3.d0*ps1(j)*ps2(j) + Ts(j) !+ sigma*ps1(i)**2-gamma;
      fN(j) = -ps3(j)+ps1(j)**3 + ten(k,j)*ps1(j)!-sigma*ps2(i);

      ! Obtain the x and y force
      fx(k,j) = fT(j)*cos(psi(k,j)) - fN(j)*sin(psi(k,j));
      fy(k,j) = fT(j)*sin(psi(k,j)) + fN(j)*cos(psi(k,j));
    end do
  end do

end subroutine force
!*********************************************************************
! Steric interaction between filaments
!*********************************************************************
subroutine repelforce(N,alpha,per,dom,psi,X,Y,Frepel)
  use constants; use omp_lib
  implicit none;
  integer, intent(in) :: N,alpha,per
  real(kind=8), intent(in) :: psi(alpha,N), X(alpha,N), Y(alpha,N),dom
  real(kind=8), intent(out):: Frepel(alpha,2,N)
  real(kind=8) :: h, s(N), cutoff, epsi
  real(kind=8) :: freg(N), fs0(N), s0
  real(kind=8) :: dis, nhat(2), Vrep, amp
  real(kind=8) :: Fx(alpha,N), Fy(alpha,N)
  real(kind=8) :: FtR(alpha,N), FnR(alpha,N)
  real(kind=8) :: Xper(alpha,N), Yper(alpha,N),xsrc,ysrc
  real(kind=8) :: ftxk(N), ftyk(N), ftxl(N), ftyl(N), dbase
  real(kind=8) :: vsrc(2,N), veval(2,N), delv(2), vdot, gdis
  integer :: i,j,k,l,p,Q,ind

  ! Dissipation coefficient
  gdis = 100.d0;

  ! Initialize
  Fx = 0.d0; Fy = 0.d0; freg = 0.d0
  FtR= 0.d0; FnR= 0.d0;

  dbase = abs(X(1,1)-X(2,1));

  ! Define arc length and other quantities
  h = 1.d0/(N-1.d0);
  do i = 1,N
    s(i) = h*(i-1.d0);
  end do
  ! Define parameters for repulsion
  ! p = decay rate; cutoff = radius of influence; amp = amplitude
  ! epsi = the regularization length of the delta function approximation
  p = 12; cutoff  = 2*h; amp = 0.05d0; epsi = 8*h;

  ! Compute repulsive forces
  do Q = -per,per
    Xper = X + dom*Q; Yper = Y;
    !!$OMP PARALLEL PRIVATE(xsrc,ysrc,dis,nhat,Vrep,epsi,s0,freg) shared(Fx,Fy)
    !!$OMP DO REDUCTION(+:Fx,Fy)
        do k = 1,alpha ! We will obtain the force on the k filament
          do l = k+1,alpha ! Will consider repulsion from l filament
            do i = 1,N ! This is the point on which the repulsion will act
              do j = 1,N
            if (abs(X(k,1)-Xper(l,1)) < 1.5d0*dbase) then
              xsrc = Xper(l,j); ysrc = Yper(l,j);
              ! distance between points
              dis = sqrt((xsrc-X(k,i))**2 + (ysrc-Y(k,i))**2);
            if (dis < cutoff) then
              nhat(1) = (xsrc-X(k,i))/dis;
              nhat(2) = (ysrc-Y(k,i))/dis;

              ! Compute repulsive force and adaptively change epsi
              Vrep  = amp*(cutoff/dis)**p;
              if (Vrep > 1e3) then
                epsi = 6*h;
              else
                epsi = 6*h;
              end if

              ! Compute the regularization vector and add repulsion
              s0 = s(i); call regdelta(N,s,freg,epsi,s0);
              Fx(k,:)= Fx(k,:) - freg*Vrep*nhat(1);
              Fy(k,:)= Fy(k,:) - freg*Vrep*nhat(2);


              s0 = s(j); call regdelta(N,s,freg,epsi,s0);
              Fx(l,:) = Fx(l,:) + freg*Vrep*nhat(1);
              Fy(l,:) = Fy(l,:) + freg*Vrep*nhat(2);
            end if
          end if ! end of checking between base distance
        end do ! end of j loop
      end do ! end of i loop
    end do ! end of periodic box loop
   end do ! end of l loop
   !!$OMP END DO
   !!$OMP END PARALLEL
  end do ! end of k loop


  ! Now project the x and y forces on tangent and normal
  do k = 1,alpha
    Frepel(k,1,:)= Fx(k,:)*cos(psi(k,:)) + Fy(k,:)*sin(psi(k,:));
    Frepel(k,2,:)=-Fx(k,:)*sin(psi(k,:)) + Fy(k,:)*cos(psi(k,:));
 end do

end subroutine repelforce
!*****************************************************************************************
! Create regularized delta function array around a point
!*****************************************************************************************
subroutine regdelta(N,s,freg,epsi,s0)
  use constants;
  implicit none
  integer, intent(in) :: N
  real(kind=8), intent(in) :: s(N),epsi,s0
  real(kind=8), intent(out):: freg(N)
  real(kind=8) :: zeta, eps, h
  integer :: i,j,k

 ! Initialize regularization array
 freg = 0.d0

  ! Create the arc-length array
  h = 1.d0/(N-1.d0);

  ! This is the coordinate of s0
  j = floor(s0/h)+1;
  k = floor(epsi/h)-1;

  ! This is delta function
  do i = j-k,j+k
    if (i > 0 .and. i < N+1) then
      zeta = abs((s0-s(i)))/epsi
      freg(i) = 0.5*(1.d0+cos(pi*zeta));
    end if
  end do

end subroutine regdelta
!*********************************************************************
! Calculate first or second derivative of any given vector
!*********************************************************************
subroutine gder(N,order,v,der)
  implicit none
  integer, intent(in) :: N, order
  real(kind=8), intent(in) :: v(N)
  real(kind=8), intent(out) :: der(N)
  real(kind=8) :: h, h2, h3, h4
  integer :: i,j

  h = 1.d0/(N-1.d0); h2 = h**2.d0; h3 = h**3; h4 = h**4;

  if (order == 1) then  ! First derivative
    do i = 1,N
      if (i == 1) then
        der(1) = (-v(3)+4.d0*v(2)-3.d0*v(1))/(2.d0*h);
      elseif (i == N) then
        der(i) = (v(i-2)-4.d0*v(i-1)+3.d0*v(i))/(2.d0*h);
      else
        der(i) = (v(i+1) - v(i-1))/(2.d0*h);
      end if
    end do
  end if

  if (order == 2) then  ! Second derivative
    do i = 1,N
      if (i == 1) then
        der(i) = (-v(i+3)+4.d0*v(i+2)-5.d0*v(i+1) + 2.d0*v(i))/h2;
      elseif (i == N) then
        der(i) = (-v(i-3)+4.d0*v(i-2)-5.d0*v(i-1) + 2.d0*v(i))/h2;
      else
        der(i) = (v(i+1) -2.d0*v(i) + v(i-1))/h2;
      end if
    end do
  end if

  if (order == 3) then  ! Third derivative
    do j = 1,N
    if (j.eq.1) then
          der(j)=(-1.5d0*v(j+4)+7.d0*v(j+3)-12.d0*v(j+2)+9.d0*v(j+1)-2.5d0*v(j))/h3
    elseif (j.eq.2) then
          der(j)=(-0.5d0*v(j+3)+3.d0*v(j+2)-6.d0*v(j+1)+5.d0*v(j)-1.5d0*v(j-1))/h3
    elseif (j.eq.N-1) then
          der(j)=(0.5d0*v(j-3)-3.d0*v(j-2)+6.d0*v(j-1)-5.d0*v(j)+1.5d0*v(j+1))/h3
    elseif (j.eq.N) then
          der(j)=(1.5d0*v(j-4)-7.d0*v(j-3)+12.d0*v(j-2)-9.d0*v(j-1)+2.5d0*v(j))/h3
    else
          der(j)=(v(j+2)-2.d0*v(j+1)+2.d0*v(j-1)-v(j-2))/(2*h3)
    endif
   end do
 end if

if (order == 4) then ! Fourth derivative
 do j = 1,N
   if (j.eq.1) then
         der(j)=(-2.d0*v(j+5)+11.d0*v(j+4)-24.d0*v(j+3)+26.d0*v(j+2)-14.d0*v(j+1)+3.d0*v(j))/h4
   elseif (j.eq.2) then
         der(j)=(-1.d0*v(j+4)+6.d0*v(j+3)-14.d0*v(j+2)+16.d0*v(j+1)-9.d0*v(j)+2.d0*v(j-1))/h4
   elseif (j.eq.N-1) then
         der(j)=(-1.d0*v(j-4)+6.d0*v(j-3)-14.d0*v(j-2)+16.d0*v(j-1)-9.d0*v(j)+2.d0*v(j+1))/h4
   elseif (j.eq.N) then
         der(j)=(-2.d0*v(j-5)+11.d0*v(j-4)-24.d0*v(j-3)+26.d0*v(j-2)-14.d0*v(j-1)+3.d0*v(j))/h4
   else
         der(j)=(v(j+2)-4.d0*v(j+1)+6.d0*v(j)-4.d0*v(j-1)+v(j-2))/h4
   endif
 end do
end if

end subroutine gder
! !------------TRAPEZOIDAL INTEGRAL ROUTINE-------------------------
!****************************************************************************************
! Trapzoidal integration
!****************************************************************************************
subroutine trapz(N,out,h,v)
  use constants
  implicit none
  integer, intent(in) :: N
  real(kind=8), intent(in) :: h,v(N)
  real(kind=8), intent(out) :: out
  integer :: i

  out = 0.d0
  if (N == 1) then
    out = 0.d0
  elseif (N==2) then
    out = h*(v(1)+v(2))/2.d0
  else
    do i = 1,N
      out = out + h*v(i);
    end do
    out = out - h/24.d0*(15.d0*(v(1)+v(N))-4.d0*(v(2)+v(N-1))+(v(3)+v(N-2)));
  end if

end subroutine trapz
!********************************************************************************************
! Subroutine for LU decomposition from NR
! The subroutine takes the following inputs
! a = matrix NXN
! n = size of the matrix a
! np = physical dimension that are used for memory allocation, 500 for 's code style
! indx = this is an integer array that records the permutation order of the LU subst
! d = it takes the value of either +- 1 depending on the odd-even permutation
!*********************************************************************************************
SUBROUTINE ludcmp(a,n,np,indx,d)
implicit none
INTEGER, intent(in) :: n,np
integer, intent(out) :: indx(n)
integer :: NMAX
real(kind=8), intent(inout) :: a(np,np)
real(kind=8), intent(out) :: d
real(kind=8) :: TINY
PARAMETER (NMAX=500,TINY=1.d-20)
INTEGER i,imax,j,k
real(kind=8) aamax,dum,sum,vv(NMAX)
d=1.d0
do i=1,n
   aamax=0.d0
   do j=1,n
      if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
   enddo
   !if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
   vv(i)=1.d0/aamax
enddo
do j=1,n
   do i=1,j-1
      sum=a(i,j)
      do k=1,i-1
         sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
   enddo
   aamax=0.d0
   do i=j,n
      sum=a(i,j)
      do k=1,j-1
         sum=sum-a(i,k)*a(k,j)
      enddo
      a(i,j)=sum
      dum=vv(i)*abs(sum)
      if (dum.ge.aamax) then
         imax=i
         aamax=dum
      endif
   enddo
   if (j.ne.imax)then
      do k=1,n
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
      enddo
      d=-d
      vv(imax)=vv(j)
   endif
   indx(j)=imax
   if(a(j,j).eq.0.d0) a(j,j)=TINY
   if(j.ne.n) then
      dum=1.d0/a(j,j)
      do i=j+1,n
         a(i,j)=a(i,j)*dum
      enddo
   endif
enddo
end SUBROUTINE ludcmp
!*************************************************************************************
! SUBROUTINE for LU backsubstitution
! This routine has to be use in conjunction with ludcmp
! a = the input is the matrix from output of lubksb
! n = size of the matrix a
! np = physical dimension that are used for memory allocation, 500 for 's code style
! indx = this is an integer array that records the permutation order of the LU subst
! b = This is input/output of the subroutine. b comes as the input of the RHS of the eqn Ax = b.
!     Then it gets modified and the solution x of the equation is stored in b
!*************************************************************************************
SUBROUTINE lubksb(a,n,np,indx,b)
INTEGER, intent(in) :: n,np,indx(n)
real(kind=8), intent(in) :: a(np,np)
real(kind=8), intent(inout) ::  b(n)
INTEGER i,ii,j,ll
real(kind=8) :: sum
ii=0
do i=1,n
   ll=indx(i)
   sum=b(ll)
   b(ll)=b(i)
   if (ii.ne.0)then
      do j=ii,i-1
         sum=sum-a(i,j)*b(j)
      enddo
   elseif (sum.ne.(0.d0)) then
      ii=i
   endif
   b(i)=sum
enddo
do i=n,1,-1
   sum=b(i)
   do j=i+1,n
      sum=sum-a(i,j)*b(j)
   enddo
   b(i)=sum/a(i,i)
enddo
end SUBROUTINE lubksb
!***********************************************************************
! Bandec
! LU decomposition of banded matrix, from Numerical recipes, O(N)
!***********************************************************************
!tab
subroutine bandec(a,n,m1,m2,np,mp,al,mpl,indx,d)
  implicit none
  integer, intent(in) :: m1,m2,mp,mpl,n,np
  real(kind=8), intent(inout) :: a(np,mp)
  real(kind=8), intent(out) :: al(np,mpl),d
  integer, intent(out) :: indx(n)
  real(kind=8) :: TINY
  parameter(TINY=1.d-20)
  integer i,j,k,l,mm
  real(kind=8) :: dum

  mm=m1+m2+1
  ! if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
  l=m1
  do i=1,m1
  	do j=m1+2-i,mm
  		a(i,j-l)=a(i,j)
  	enddo
  	l=l-1
  	do j=mm-l,mm
  		a(i,j)=0
  	end do
  enddo

  d=1.d0
  l=m1
  do k=1,n
  	dum=a(k,1)
  	i=k
  	if(l.lt.n) l=l+1
  	do j=k+1,l
  		if(abs(a(j,1)).gt.abs(dum)) then
  			dum=a(j,1)
  			i=j
  		endif
  	enddo
  	indx(k)=i
  	if(dum.eq.0.d0) a(k,1)=TINY
  	if(i.ne.k) then
  		d=-d
  		do j=1,mm
  			dum=a(k,j)
  			a(k,j)=a(i,j)
  			a(i,j)=dum
  		enddo
  	endif
  	do i=k+1,l
  		dum=a(i,1)/a(k,1)
  		al(k,i-k)=dum
  		do j=2,mm
  			a(i,j-1)=a(i,j)-dum*a(k,j)
  		enddo
  		a(i,mm)=0
  	enddo
  enddo

end subroutine bandec


!***********************************************************************
!c Banbks
!c 		Back substitution of banded matrix, from Numerical recipes, O(N)
!***********************************************************************

subroutine banbks(a,n,m1,m2,np,mp,al,mpl,indx,b)
implicit none

integer, intent(in) :: m1,m2,mp,mpl,n,np,indx(n)
real(kind=8), intent(in) :: a(np,mp),al(np,mpl)
real(kind=8), intent(out) :: b(n)
integer i,k,l,mm
real(kind=8) :: dum

mm=m1+m2+1
! if(mm.gt.mp.or.m1.gt.mpl.or.n.gt.np) pause 'bad args in bandec'
l=m1
do k=1,n
	i=indx(k)
	if(i.ne.k) then
		dum=b(k)
		b(k)=b(i)
		b(i)=dum
	endif
	if(l.lt.n) l=l+1
	do i=k+1,l
		b(i)=b(i)-al(k,i-k)*b(k)
	enddo
enddo
l=1
do i=n,1,-1
	dum=b(i)
	do k=2,l
		dum=dum-a(i,k)*b(k+i-1)
	enddo
	b(i)=dum/a(i,1)
	if(l.lt.mm) l=l+1
enddo

!return
end
! ***********************************************************************
! * GASDEV returns double precision figures from a Gaussian Distribution
! * with mean 0 and variance 1. 'IDUM' is set to a negative number to
! * initialize and then not changed. Reset IDUM every time, with a new
! * negative number (iteration number X -1, for instance) to avoid
! * repetition of saem sequence
! ***********************************************************************
    real(kind=8) function gasdev(idum)
		integer, intent(inout):: idum
		integer :: iset
		real(kind=8) :: fac,gset,rsq,v1,v2
    real(kind=8), external :: ran1

		!common/rndm/idum
		SAVE iset,gset
		DATA iset/0/
		if (idum.lt.0) iset=0
		if (iset.eq.0) then
  1	 	v1=2.d0*ran1(idum)-1.d0
			v2=2.d0*ran1(idum)-1.d0
			rsq=v1**2.d0+v2**2.d0
			if(rsq>=1.d0 .or. rsq==0.d0) go to 1
			fac=sqrt(-2.d0*log(rsq)/rsq)
			gset=v1*fac
			gasdev=v2*fac
			iset=1
		else
			gasdev=gset
			iset=0
		endif
		return
		END
    !--------------------------------------------------------------------------------------

    real(kind=8) function ran1(idum)
		!common/rndm/idum
		INTEGER, intent(inout):: idum
		integer, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32,NDIV=1+(IM-1)/NTAB
    real(kind=8), parameter ::  AM=1.d0/IM,EPS=1.2d-7,RNMX=1.d0-EPS
		INTEGER :: j,k,iv(NTAB),iy
		SAVE iv,iy
		DATA iv /NTAB*0/, iy /0/
		if (idum<=0 .or. iy==0) then
			idum=max(-idum,1)
			do j=NTAB+8,1,-1
				k=idum/IQ
				idum=IA*(idum-k*IQ)-IR*k
				if (idum<0) idum=idum+IM
				if (j<=NTAB) iv(j)=idum
			enddo
			iy=iv(1)
		endif
		k=idum/IQ
		idum=IA*(idum-k*IQ)-IR*k
		if(idum<0) idum=idum+IM
		j=1+iy/NDIV
		iy=iv(j)
		iv(j)=idum
		ran1=min(AM*iy,RNMX)
		return
		END
!***********************************************************************
!		Tridiagonal linear system solver, from Numerical recipes, O(N)
! 		Arguments: a,b,c are diagonals, r is RHS, u is solution
!***********************************************************************
subroutine tridag(a,b,c,r,u,n)
  implicit none
  integer, intent(in) :: N
  real(kind=8), intent(in) :: a(N), b(N), c(N),r(N)
  real(kind=8), intent(out):: u(N)
  real(kind=8) :: gam(N), bet
  integer :: nmax,j

  bet=b(1)
  u(1)=r(1)/bet
  do j=2,n
     gam(j)=c(j-1)/bet
     bet=b(j)-a(j)*gam(j)
     !if(bet.eq.0) pause 'tridag failed'
     u(j)=(r(j)-a(j)*u(j-1))/bet
  enddo
  do j=n-1,1,-1
     u(j)=u(j)-gam(j+1)*u(j+1)
  enddo


end subroutine tridag
! !************************************************************************************
! ! Prepare the nonlinear equation vector/residue vector
! !************************************************************************************
! subroutine nonlineq(N,delt,psinew,psiold,psiold2,ten,ufiber,nonlin)
!   use constants
!   implicit none
!   integer, intent(in) :: N
!   real(kind=8), intent(in) :: psinew(N), psiold(N), psiold2(N)
!   real(kind=8), intent(in) :: delt, ten(N), ufiber(2,N)
!   real(kind=8), intent(out):: nonlin(2*N)
!   real(kind=8) :: h, h2, beta, cpll, cperp
!   real(kind=8) :: Ts(N),  Tss(N)
!   real(kind=8) :: ps1(N), ps2(N), ps3(N), ps4(N), psidot(N)
!   real(kind=8) :: Ut(N), Un(N), Uts(N), Uns(N)
!   real(kind=8) :: E1(N), E2(N)
!   integer :: i
!
!   ! Define relevant geometric quantities and required derivatives
!   h = 1.d0/(N-1.d0); h2 = h**2;
!   cpll = -1.d0/(2.d0*c); cperp = 1.d0/(2.d0-c);
!   beta = cpll/cperp;
!   ! Calculate derivatives of the angle
!   call gder(N,1,psinew,ps1); call gder(N,2,psinew,ps2)
!   call gder(N,3,psinew,ps3); call gder(N,4,psinew,ps4)
!   ! Calculate derivatives of the tension
!   call gder(N,1,ten,Ts);     call gder(N,2,ten,Tss)
!   ! Define components of disturbance flow and motor related stuff
!   Ut = ufiber(1,:);    Un = ufiber(2,:);
!   call gder(N,1,Ut,Uts); call gder(N,1,Un,Uns)
!   psidot = (3.d0*psinew-4.d0*psiold+psiold2)/(2.d0*delt)
!   ! Now create nonliear equations as eq = 0
!
!   ! FIRST WE CREATE THE ARRAY E1
!   E1(1)= Ts(1) + ps2(1)*ps1(1); ! EBC1
!   do i = 2,N-1
!     E1(i) = Tss(i) + (1.d0+beta)*ps3(i)*ps1(i) + ps2(i)**2 - &
!           & beta*ten(i)*(ps1(i))**2 - cpll*(ps1(i)*Un(i)-Uts(i));
!   end do
!   E1(N)= Ts(N) + ps2(N)*ps1(N) + Vel*cpll ! EBC2
!
!   ! NEXT WE CREATE E2
!   E2(1)= -ps3(1) + ten(1)*ps1(1); ! EBC3
!   E2(2)= psinew(1)-pi/2.d0; ! EBC5
!   do i = 3,N-2
!     E2(i) =-ps4(i) + ps2(i)*(ps1(i)**2)/beta + ten(i)*ps2(i) + &
!           & (1.d0+1.d0/beta)*Ts(i)*ps1(i) - cperp*psidot(i)  + &
!           & cperp*(Ut(i)*ps1(i)+Uns(i));
!   end do
!   E2(N-1)= psinew(N)-pi/2.d0 ! EBC6
!   E2(N)  =-ps3(N) + ten(N)*ps1(N); ! EBC4
!
!   ! Now assemble all the equations into one big array
!   nonlin(1:N) = E1; nonlin(N+1:2*N) = E2;
!
! end subroutine nonlineq
! !*************************************************************************************
! ! Jacobian for Newton solve
! !*************************************************************************************
! subroutine jacobian(N,delt,psi,ten,ufiber,JAC)
!   use constants
!   implicit none
!   integer, intent(in) :: N
!   real(kind=8), intent(in) :: delt, psi(N), ten(N), ufiber(2,N)
!   real(kind=8), intent(out):: JAC(2*N,2*N)
!   real(kind=8) :: h, h2, h3, h4, beta, cpll, cperp
!   real(kind=8) :: Ts(N),Ut(N), Un(N)
!   real(kind=8) :: ps1(N), ps2(N), ps3(N)
!   real(kind=8) :: de1_dT(N,N), de1_dpsi(N,N)
!   real(kind=8) :: de2_dT(N,N), de2_dpsi(N,N)
!   integer :: i
!
!   ! Define relevant geometric quantities and required derivatives
!   h = 1.d0/(N-1.d0); h2 = h**2; h3 = h**3; h4 = h**4;
!   cpll = -1.d0/(2.d0*c); cperp = 1.d0/(2.d0-c);
!   beta = cpll/cperp;
!
!   ! Take derivatives
!   call gder(N,1,psi,ps1); call gder(N,2,psi,ps2)
!   call gder(N,3,psi,ps3); call gder(N,1,ten,Ts)
!
!   ! Define disturbance velocities
!   Ut = ufiber(1,:); Un = ufiber(2,:);
!
!   ! THE ORDER OF EQUATIONS:
!   ! E1 = [EBC1 E1 EBC2]
!   ! E2 = [EBC3 EBC5 E2 EBC6 EBC4]
!
!   ! LIST OF BCS:
!   ! EBC1:= Ts + Pss*Ps = 0    at s =  0
!   ! EBC2:= TS + Pss*Ps =-v    at s = 1
!   ! EBC3:=  -Psss + T*Ps  = 0 at s = 0
!   ! EBC4:=  -Psss + T*Ps  = 0 at s = 1
!   ! EBC5:=  P(1) = pi/2
!   ! EBC6:=  P(N) = pi/2
!
!   ! ORDER OF VARIABLE ARRAGNEMENT
!   ! X = [T(1),T(2) ... T(N); psi(1),psi(2) ... psi(N)]
!
!   ! START CONSTRUCTING JACOBIANS
!
!   !---------------- A: Jacobian components related to equation 1----------------
!   ! de1/dT
!   de1_dT = 0.d0;
!   de1_dT(1,1) = -3.d0/(2.d0*h); de1_dT(1,2) = 4.d0/(2.d0*h); de1_dT(1,3) = -1.d0/(2.d0*h);   ! Comes from the BC: EBC1
!   do i = 2,N-1
!     de1_dT(i,i-1)=  1.d0/h2;
!     de1_dT(i,i)  = -2.d0/h2 - beta*(ps1(i))**2; !beta = zeta_pll/zeta_perp
!     de1_dT(i,i+1)=  1.d0/h2;
!   end do
!   de1_dT(N,N) = 3.d0/(2.d0*h); de1_dT(N,N-1) =-4.d0/(2.d0*h);  de1_dT(N,N-2) = 1.d0/(2.d0*h) ! Comes from the BC: EBC2
!
!   ! de1/dpsi
!   de1_dpsi = 0.d0;
!   ! First apply the boundary condition EBC1
!   de1_dpsi(1,1) = (-3.d0*ps2(1)/(2.d0*h)+2.d0*ps1(1)/h2);
!   de1_dpsi(1,2) = ( 4.d0*ps2(1)/(2.d0*h)-5.d0*ps1(1)/h2);
!   de1_dpsi(1,3) = (-1.d0*ps2(1)/(2.d0*h)+4.d0*ps1(1)/h2);
!   de1_dpsi(1,4) = (0.d0                 -1.d0*ps1(1)/h2);
!   ! Firt add all the terms that do not involve the 3rd derivative of phi
!   do i = 2,N-1
!     de1_dpsi(i,i-1)=  2.d0*ps2(i)/h2 + 2.d0*beta*ten(i)*ps1(i)/(2.d0*h) + cpll*Un(i)/(2.d0*h);
!     de1_dpsi(i,i)  = -4.d0*ps2(i)/h2;
!     de1_dpsi(i,i+1)=  2.d0*ps2(i)/h2 - 2.d0*beta*ten(i)*ps1(i)/(2.d0*h) - cpll*Un(i)/(2.d0*h);
!
!     if ((i > 2) .and. (i < N-1)) then
!       de1_dpsi(i,i-2) = de1_dpsi(i,i-2) + (1.d0+beta)*ps1(i)*(-1.d0/(2*h3));
!       de1_dpsi(i,i-1) = de1_dpsi(i,i-1) + (1.d0+beta)*(-ps3(i)/(2.d0*h)+2.d0*ps1(i)/(2*h3));
!       de1_dpsi(i,i+1) = de1_dpsi(i,i+1) + (1.d0+beta)*( ps3(i)/(2.d0*h)-2.d0*ps1(i)/(2*h3));
!       de1_dpsi(i,i+2) = de1_dpsi(i,i+2) + (1.d0+beta)*ps1(i)*( 1.d0/(2*h3));
!     elseif (i == 2) then
!       de1_dpsi(i,i-1) = de1_dpsi(i,i-1) + (1.d0+beta)*(-ps3(i)/(2.d0*h)-1.5d0*ps1(i)/(h3));
!       de1_dpsi(i,i)   = de1_dpsi(i,i)   + (1.d0+beta)*ps1(i)*(  5.d0/(h3));
!       de1_dpsi(i,i+1) = de1_dpsi(i,i+1) + (1.d0+beta)*( ps3(i)/(2.d0*h)- 6.d0*ps1(i)/(h3));
!       de1_dpsi(i,i+2) = de1_dpsi(i,i+2) + (1.d0+beta)*ps1(i)*(  3.d0/(h3));
!       de1_dpsi(i,i+3) = de1_dpsi(i,i+3) + (1.d0+beta)*ps1(i)*(-0.5d0/(h3));
!     elseif (i==N-1) then
!       de1_dpsi(i,i+1) = de1_dpsi(i,i-1) + (1.d0+beta)*( ps3(i)/(2.d0*h)+1.5d0*ps1(i)/(h3));
!       de1_dpsi(i,i)   = de1_dpsi(i,i)   + (1.d0+beta)*ps1(i)*( -5.d0/(h3));
!       de1_dpsi(i,i-1) = de1_dpsi(i,i-1) + (1.d0+beta)*(-ps3(i)/(2.d0*h)+6.d0*ps1(i)/(h3));
!       de1_dpsi(i,i-2) = de1_dpsi(i,i-2) + (1.d0+beta)*ps1(i)*( -3.d0/(h3));
!       de1_dpsi(i,i-3) = de1_dpsi(i,i-3) + (1.d0+beta)*ps1(i)*( 0.5d0/(h3));
!     end if
!   end do
!   ! Finally apply the boundary condition EBC2
!   de1_dpsi(N,N)   = ( 3.d0*ps2(N)/(2.d0*h)+2.d0*ps1(N)/h2);
!   de1_dpsi(N,N-1) = (-4.d0*ps2(N)/(2.d0*h)-5.d0*ps1(N)/h2);
!   de1_dpsi(N,N-2) = ( 1.d0*ps2(N)/(2.d0*h)+4.d0*ps1(N)/h2);
!   de1_dpsi(N,N-3) = (0.d0                 -1.d0*ps1(N)/h2);
!
!   !---------------- B: Jacobian components related to equation 2----------------
!   ! de2/dT
!   de2_dT = 0.d0;
!   de2_dT(1,1) = ps1(1);  !From BC: EBC3
!   ! No contribution from EBC5
!   do i = 3,N-2
!     de2_dT(i,i-1)= -(1.d0+1.d0/beta)*ps1(i)/(2.d0*h);
!     de2_dT(i,i)  =  ps2(i);
!     de2_dT(i,i+1)=  (1.d0+1.d0/beta)*ps1(i)/(2.d0*h);
!   end do
!   de2_dT(N,N) = ps1(N); ! From BC: EBC4
!
!   ! de2/dpsi
!   de2_dpsi  = 0.d0;
!   ! Contribution from EBC3
!   de2_dpsi(1,1)  = (-3.d0*ten(1)/(2.d0*h)-(-2.5d0/h3));
!   de2_dpsi(1,2)  = ( 4.d0*ten(1)/(2.d0*h)-( 9.0d0/h3));
!   de2_dpsi(1,3)  = (-1.d0*ten(1)/(2.d0*h)-(-12.d0/h3));
!   de2_dpsi(1,4)  = (0.d0                 -( 7.0d0/h3));
!   de2_dpsi(1,5)  = (0.d0                 -(-1.5d0/h3));
!   ! Contribution from EBC5
!   de2_dpsi(2,1)  = 1.d0;
!   ! Structure for the actual equations
!   do i = 3,N-2
!     ! first write down contribution from -phi_{ssss}
!     de2_dpsi(i,i-2) =-1.d0/h4;
!     de2_dpsi(i,i-1) = 4.d0/h4;
!     de2_dpsi(i,i)   =-6.d0/h4;
!     de2_dpsi(i,i+1) = 4.d0/h4;
!     de2_dpsi(i,i+2) =-1.d0/h4;
!
!     ! Contribution from: T*pss + (1+1/beta)*ps*Ts-cp*phi_t
!     de2_dpsi(i,i-1)=  de2_dpsi(i,i-1) + ten(i)/h2 - (1.d0+1.d0/beta)*Ts(i)/(2.d0*h) - cperp*Ut(i)/(2.d0*h)
!     de2_dpsi(i,i)  =  de2_dpsi(i,i)   + -2.d0*ten(i)/h2 - 3.d0*cperp/(2.d0*delt);
!     de2_dpsi(i,i+1)=  de2_dpsi(i,i+1) + ten(i)/h2 + (1.d0+1.d0/beta)*Ts(i)/(2.d0*h) + cperp*Ut(i)/(2.d0*h)
!
!     ! Contribution from: pss*ps1^2/beta
!     de2_dpsi(i,i-1)= de2_dpsi(i,i-1) + (-2.d0*ps2(i)*ps1(i)/(2.d0*h) + ps1(i)**2/h2)/beta;
!     de2_dpsi(i,i)  = de2_dpsi(i,i)   + (-2.d0*ps1(i)**2/h2)/beta;
!     de2_dpsi(i,i+1)= de2_dpsi(i,i+1) + ( 2.d0*ps2(i)*ps1(i)/(2.d0*h) + ps1(i)**2/h2)/beta;
!   end do
!   ! Contribution from EBC6
!   de2_dpsi(N-1,N) = 1.d0;
!   ! Contribution from EBC4
!   de2_dpsi(N,N)  = ( 3.d0*ten(N)/(2.d0*h)-( 2.5d0/h3));
!   de2_dpsi(N,N-1)= (-4.d0*ten(N)/(2.d0*h)-(-9.0d0/h3));
!   de2_dpsi(N,N-2)= ( 1.d0*ten(N)/(2.d0*h)-( 12.d0/h3));
!   de2_dpsi(N,N-3)= (0.d0                 -(-7.0d0/h3));
!   de2_dpsi(N,N-4)= (0.d0                 -( 1.5d0/h3));
! !----------------- ASSEMBLE THE WHOLE JACOBIAN MATRIX-------------------------
!   JAC = 0.d0;
!   ! Fill up in blocks
!
!   ! First N rows to E1
!   JAC(1:N,1:N)      = de1_dT;
!   JAC(1:N,N+1:2*N)  = de1_dpsi;
!
!   ! Second N rows to E2
!   JAC(N+1:2*N,1:N)      = de2_dT;
!   JAC(N+1:2*N,N+1:2*N)  = de2_dpsi;
!   !----------------------------------------------------------------------------
!
!
! end subroutine jacobian
