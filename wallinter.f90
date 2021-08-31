!*****************************************************************************
! Subroutine to calculate interaction due to wall or free space and self
!*****************************************************************************
subroutine wallinterper(N,alpha,per,domain,x0,psi,ten,HI)
  use constants; use omp_lib
  implicit none
  integer, intent(in) :: N, alpha,per
  real(kind=8), intent(in) :: domain,x0(alpha), psi(alpha,N)
  real(kind=8), intent(in) :: ten(alpha,N)
  real(kind=8), intent(out):: HI(alpha,2,N)
  real(kind=8) :: h, X(alpha,N), Y(alpha,N), theta(N)
  real(kind=8) :: fx(alpha,N), fy(alpha,N)
  real(kind=8) :: inter(alpha,2,N), xeval, yeval, xsrc, ysrc
  real(kind=8) :: R(2), Rhat(2), Rmag, dot, Ux, Uy
  real(kind=8) :: adx(N), ady(N), dbx(N), dby(N)
  real(kind=8) :: Xtar(2), Xsource(2), f0(2), U(2), dis
  real(kind=8) :: ftR(alpha,N), fnR(alpha,N), FxR(alpha,N), FyR(alpha,N)
  integer :: i,j,k,l,flag

 !print*, "hi"
! Define spacing
  h = 1.d0/(N-1.d0); inter = 0.d0;
! First calculate the positions and set domain
call position(N,alpha,psi,x0,X,Y)
! Now calculate the forces in the x and y direction
call force(N,alpha,psi,ten,fx,fy)

 ! Now calculate the interactions
 inter = 0.d0;
 !$OMP PARALLEL PRIVATE(xeval,yeval,xsrc,ysrc,Xtar,Xsource,f0,U)
 !$OMP DO
 do k = 1,alpha ! where we want disturbance (target)
   do l = 1,alpha ! SOURCE of disturbance
    !if (k .ne. l) then
     do i = 1,N
       xeval = X(k,i); yeval = Y(k,i) ! where you want vel (target)
       do j = 1,N
         xsrc = X(l,j); ysrc = Y(l,j) ! where the point forces are(source)

          ! Correct for periodicity of the domain
          if (per == 1) then
            xeval = xeval - domain*floor(xeval/domain)
            xsrc  = xsrc  - domain*floor(xsrc/domain)
          end if

          ! Create data in form of greens subroutine
          Xtar(1)    = xeval; Xtar(2)    = yeval;
          Xsource(1) = xsrc;  Xsource(2) = ysrc;

          ! Now calculate disturbance velocity
          f0(1) = 8.d0*pi*h*fx(l,j);
          f0(2) = 8.d0*pi*h*fy(l,j);
          call greens(per,domain,Xtar,Xsource,f0,U)

          ! Add the velocity to interaction
          inter(k,1,i) = inter(k,1,i) + U(1);
          inter(k,2,i) = inter(k,2,i) + U(2);
        end do
      end do
    !end if
   end do ! loop of l ends (sources)
 end do   ! loop of k ends (targets)
 !$OMP END DO
 !$OMP END PARALLEL


  !! Now convert to tangent and normal direction
  !$OMP PARALLEL
  !$OMP DO
  do k = 1,alpha
    do i =1,N
      HI(k,1,i)= inter(k,1,i)*cos(psi(k,i)) + inter(k,2,i)*sin(psi(k,i));
      HI(k,2,i)=-inter(k,1,i)*sin(psi(k,i)) + inter(k,2,i)*cos(psi(k,i));
    end do
  end do
  !$OMP END DO
  !$OMP END PARALLEL


end subroutine wallinterper
!******************************************************************************
! Subroutine greens
!******************************************************************************
subroutine greens(per,domain,Xtar,Xsource,f0,U)
  use constants
  implicit none
  integer, intent(in) :: per
  real(kind=8), intent(in) :: domain,Xtar(2), Xsource(2)
  real(kind=8), intent(out):: f0(2), U(2)
  real(kind=8) :: xt,yt,xsrc,ysrc,yimage
  real(kind=8) :: fi(2),g(2),h
  real(kind=8) :: x0p,xip,r0(2),ri(2),R0m,Rim
  real(kind=8) :: dot1,dot2,dot3,e(2)
  real(kind=8) :: us0(2),usi(2),upd(2),usd(2)
  integer :: i,j,p,Q

  ! Define forces for image
  fi=-f0;
  g = f0; g(1) = -f0(1);
  ! Define image points
  xt    = Xtar(1);      yt  = Xtar(2);
  xsrc  = Xsource(1);   ysrc= Xsource(2);
  yimage= -ysrc
  ! Initialize velocity
  U = 0.d0;
  ! Set range of periodicity
  Q  = 0
  if (per == 1) then
    Q = 5;
  end if

do p = -Q,Q,1
  x0p = xsrc + p*domain; xip = x0p;
  ! Define position vectors
  r0(1)= xt-x0p; r0(2) = yt-ysrc;  ! Between actual points
  ri(1)= xt-xip; ri(2) = yt-yimage; ! Target and image points
  ! Define magnitude of position vectors
  R0m  = sqrt(r0(1)**2 + r0(2)**2);
  Rim  = sqrt(ri(1)**2 + ri(2)**2);
  ! Define force and position dot products
  dot1 = f0(1)*r0(1) + f0(2)*r0(2);
  dot2 = fi(1)*ri(1) + fi(2)*ri(2);
  dot3 = g(1)*ri(1)  + g(2)*ri(2);
  h = abs(yimage); e(1) = 0.d0; e(2) = 1.d0;

  do i = 1,2
    ! Stokeslet contribution
    us0(i) =  f0(i)/R0m + dot1*r0(i)/(R0m**3);
    ! Velocities from image
    usi(i) =  fi(i)/Rim + dot2*ri(i)/(Rim**3);

    upd(i) = -(h**2)*(g(i)/(Rim**3) -3*(dot3*ri(i))/(Rim**5));

    usd(i) = f0(2)*ri(i)/(Rim**3) + ri(2)*g(i)/(Rim**3) - dot3*e(i)/(Rim**3) &
           & - 3.d0*(ri(2)*dot3*ri(i))/(Rim**5);

  end do

  us0 = us0/(8.d0*pi);
  usi = usi/(8.d0*pi);
  upd = upd/(4.d0*pi);
  usd = usd*2.d0*h/(8.d0*pi);

  ! Check for self interaction
  if (p == 0 .and.  R0m < 1.d-3) then
    us0 = 0.d0;
  end if

  ! This happens at the anchor point of the filaments
  if (p .eq. 0 .and. Rim < 1.d-3) then
      usi = 0.d0; upd = 0.d0; usd = 0.d0;
  end if

  U = U + us0 + (usi+upd+usd);
end do

end subroutine greens
