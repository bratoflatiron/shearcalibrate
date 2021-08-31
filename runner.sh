rm a.out
gfortran -c constants.f90
gfortran -finit-real=zero follower_angle.f90 -O3
#./a.out


# rm *.so
# gfortran -c constants.f90
# f2py -c -m wrinkle --opt='-O3' constants.f90 newton_angle2.f90
# mv wrinkle.cpython-38-darwin.so wrinkle.so
