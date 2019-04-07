module parameters
  implicit none
  save
  integer, parameter :: sp = kind(1.e0)
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: qp = selected_real_kind(30) !
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
  real(dp), parameter :: mu = 4.d-7 * pi
  real(dp), parameter :: epsilon = 8.85d-12
  real(dp), parameter :: eps = 1.d-7
  real(dp), parameter :: Iw = 1.d0
  real(dp), parameter :: dsx = 1.d0
  real(dp), parameter :: dsy = 1.d0
  real(dp), parameter :: dsz = 1.d0
  real(dp), parameter :: mx = 1.d0
  real(dp), parameter :: my = 1.d0
  real(dp), parameter :: mz = 1.d0
end module parameters
