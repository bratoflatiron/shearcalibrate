!*********************************************************************************************
!A module containing all the constant values
!*********************************************************************************************
module constants
implicit none
real(kind=8), parameter :: pi   = 4.d0*atan(1.d0);
real(kind=8), parameter :: dteu = 1.d-7 ! A very small time step chosen for explicit euler
real(kind=8), parameter :: c    = log(exp(1.d0)*1d-4)!log(exp(1.d0)*(8.d-3/lpol)**2) ! slenderness parameter
real(kind=8), parameter :: sigma= 0.0d0;
real(kind=8), parameter :: gamma= 1.0d2;
end module constants
