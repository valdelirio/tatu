module utils
use parameters
contains
subroutine sanitizedata(n, h0, z, esp, camadT, camad, h, prof)
  implicit none
  integer :: n
  real(dp), intent(in) :: h0, z, esp(:)
  integer, intent(out) :: camadT, camad
  real(dp), dimension(:), allocatable, intent(out) :: h, prof

  integer :: i, j, k

  allocate( h(0:n), prof(-1:n) )
  if (size(esp) == n) then
    h(0)=0.d0
    h(1:n)=esp
  else
    h(0)=0.d0
    h(1:n-1)=esp
    h(n)=1.d300
  end if
  ! create depths array that suits any pathological situation
  prof(-1)=-1.d300
  prof(0)=0.d0
  if (n > 1) then
    prof(1) = h(1)
    if (n > 2) then
      do k=2,n-1
        prof(k) = prof(k-1) + h(k)
      end do
    end if
  end if
  prof(n)=1.d300

  ! find the layer where the receiver is
  camad = 0
  if (z < 0.d0) then
    camad=0
  else if (z >= prof(n-1)) then
    camad=n
  else
    do i=n-1,1,-1
      if (z >= prof(i-1)) then
        camad=i
        exit
      end if
    end do
  end if

  ! find the layer where the transmitter is
  camadT = 0
  if (h0 < 0.d0) then
    camadT = 0
  else if (h0 >= prof(n-1)) then
    camadT = n
  else
    do j=n-1,1,-1
      if (h0 >= prof(j-1)) then
        camadT = j
        exit
      end if
    end do
  end if
end subroutine sanitizedata

function int2str(num)
  implicit none
  integer, intent(in) :: num
  character(len=16) :: int2str
  write(int2str, '(i16)') num
end function int2str

function real2str(real, format)
  implicit none
  real(sp), intent(in) :: real
  character(len=*), intent(in), optional :: format
  character(len=16) :: real2str
  if (present(format)) then
    write(real2str, format) real
  else
    write(real2str, '(f6.2)') real
  end if
end function real2str

function real2strPerc(real)
  implicit none
  real(sp), intent(in) :: real
  character(len=3) :: real2strPerc
  write(real2strPerc, '(i3)') nint(real)
end function real2strPerc

subroutine commonarraysED(n, npt, r, krJ0J1, zeta, eta0, h0, h, prof, camadT, eta, u, uh, AdmInt, ImpInt, &
                          RTEdw, RTEup, RTMdw, RTMup, FEdw, FEup, AMdw, AMup, AMdwz, AMupz)
  implicit none
  integer, intent(in) :: n, npt, camadT
  real(dp), intent(in) :: r, krJ0J1(npt), h0, h(0:n), prof(-1:n)
  complex(dp), intent(in) :: eta(n), zeta, eta0
  complex(dp), dimension(npt), intent(out) :: AMdw, AMup, FEdw, FEup
  complex(dp), optional, dimension(npt), intent(out) :: AMdwz, AMupz
  complex(dp), dimension(npt,0:n), intent(out) :: u, AdmInt, ImpInt, uh, RTEdw, RTMdw, RTEup, RTMup

  integer :: i
  real(dp) :: kr(npt)
  complex(dp) :: den(npt), wvnb2(0:n), tgh(npt,0:n), AdmApdw(npt,1:n), ImpApdw(npt,1:n), AdmApup(npt,0:n-1), ImpApup(npt,0:n-1)

  kr = krJ0J1 / r
  wvnb2(0) = -zeta * eta0
  u(:,0) = sqrt(kr * kr - wvnb2(0))
  AdmInt(:,0) = u(:,0) / zeta
  ImpInt(:,0) = u(:,0) / eta0
  uh(:,0) = u(:,0) * h(0)
  tgh(:,0) = (1.0 - exp(-2.0 * uh(:,0))) / (1.0 + exp(-2.0 * uh(:,0)))
  do i = 1, n
    wvnb2(i) = -zeta * eta(i)
    u(:,i) = sqrt(kr * kr - wvnb2(i))
    AdmInt(:,i) = u(:,i) / zeta
    ImpInt(:,i) = u(:,i) / eta(i)
    uh(:,i) = u(:,i) * h(i)
    tgh(:,i) = (1.0 - exp(-2.0 * uh(:,i))) / (1.0 + exp(-2.0 * uh(:,i)))
  end do

  AdmApdw(:,n) = AdmInt(:,n)
  ImpApdw(:,n) = ImpInt(:,n)
  RTEdw(:,n) = (0.0,0.0)
  RTMdw(:,n) = (0.0,0.0)
  do i = n-1, 1, -1
    AdmApdw(:,i) = AdmInt(:,i) * (AdmApdw(:,i+1) + AdmInt(:,i) * &
                  tgh(:,i)) / (AdmInt(:,i) + AdmApdw(:,i+1) * tgh(:,i))
    ImpApdw(:,i) = ImpInt(:,i) * (ImpApdw(:,i+1) + ImpInt(:,i) * &
                  tgh(:,i)) / (ImpInt(:,i) + ImpApdw(:,i+1) * tgh(:,i))
    RTEdw(:,i) = (AdmInt(:,i) - AdmApdw(:,i+1)) / (AdmInt(:,i) + AdmApdw(:,i+1))
    RTMdw(:,i) = (ImpInt(:,i) - ImpApdw(:,i+1)) / (ImpInt(:,i) + ImpApdw(:,i+1))
  end do
  RTEdw(:,0) = (AdmInt(:,0) - AdmApdw(:,1)) / (AdmInt(:,0) + AdmApdw(:,1))
  RTMdw(:,0) = (ImpInt(:,0) - ImpApdw(:,1)) / (ImpInt(:,0) + ImpApdw(:,1))

  AdmApup(:,0) = AdmInt(:,0)
  ImpApup(:,0) = ImpInt(:,0)
  RTEup(:,0) = (0.0,0.0)
  RTMup(:,0) = (0.0,0.0)
  do i = 1, n-1
    AdmApup(:,i) = AdmInt(:,i) * (AdmApup(:,i-1) + AdmInt(:,i) * &
                tgh(:,i)) / (AdmInt(:,i) + AdmApup(:,i-1) * tgh(:,i))
    ImpApup(:,i) = ImpInt(:,i) * (ImpApup(:,i-1) + ImpInt(:,i) * &
                tgh(:,i)) / (ImpInt(:,i) + ImpApup(:,i-1) * tgh(:,i))
    RTEup(:,i) = (AdmInt(:,i) - AdmApup(:,i-1)) / (AdmInt(:,i) + AdmApup(:,i-1))
    RTMup(:,i) = (ImpInt(:,i) - ImpApup(:,i-1)) / (ImpInt(:,i) + ImpApup(:,i-1))
  end do
  RTEup(:,n) = (AdmInt(:,n) - AdmApup(:,n-1)) / (AdmInt(:,n) + AdmApup(:,n-1))
  RTMup(:,n) = (ImpInt(:,n) - ImpApup(:,n-1)) / (ImpInt(:,n) + ImpApup(:,n-1))

  den = 1.0 - RTEup(:,camadT) * RTEdw(:,camadT) * exp(-2.0 * uh(:,camadT))
  FEdw = (exp(-u(:,camadT) * (prof(camadT)-h0)) + RTEup(:,camadT) * &
        exp(u(:,camadT) * (prof(camadT-1)-h(camadT)-h0))) / den
  FEup = (exp(u(:,camadT) * (prof(camadT-1)-h0)) + RTEdw(:,camadT) * &
        exp(-u(:,camadT) * (prof(camadT)+h(camadT)-h0))) / den

  den = 1.0 - RTMup(:,camadT) * RTMdw(:,camadT) * exp(-2.0 * uh(:,camadT))
  AMdw = (exp(-u(:,camadT) * (prof(camadT) - h0)) - RTMup(:,camadT) * &
        exp(u(:,camadT) * (prof(camadT-1) - h(camadT)-h0))) / den
  AMup = (exp(u(:,camadT) * (prof(camadT-1)-h0)) - RTMdw(:,camadT) * &
        exp(-u(:,camadT) * (prof(camadT)+h(camadT)-h0))) / den

  if (present(AMdwz) .and. present(AMupz)) then
    AMdwz = (exp(-u(:,camadT) * (prof(camadT)-h0)) + RTMup(:,camadT) * &
          exp(u(:,camadT) * (prof(camadT-1) - h(camadT) - h0))) / den
    AMupz = (exp(u(:,camadT) * (prof(camadT-1)-h0)) + RTMdw(:,camadT) * &
          exp(-u(:,camadT) * (prof(camadT)+h(camadT)-h0))) / den
  end if
end subroutine commonarraysED

subroutine commonarraysMD(n, npt, hordist, krJ0J1, zeta, eta0, h0, h, prof, camadT, eta, u, uh, AdmInt, ImpInt, &
                          RTEdw, RTEup, RTMdw, RTMup, FEdw, FEup, AMdw, AMup, FEdwz, FEupz)
  implicit none
  integer, intent(in) :: n, npt, camadT
  real(dp), intent(in) :: hordist, krJ0J1(npt), h0, h(0:n), prof(-1:n)
  complex(dp), intent(in) :: eta(n), zeta, eta0
  complex(dp), dimension(npt), intent(out) :: AMdw, AMup, FEdw, FEup
  complex(dp), optional, dimension(npt), intent(out) :: FEdwz, FEupz
  complex(dp), dimension(npt,0:n), intent(out) :: u, AdmInt, ImpInt, uh, RTEdw, RTMdw, RTEup, RTMup

  integer :: i
  real(dp) :: r, kr(npt)
  complex(dp) :: den(npt), wvnb2(0:n), tgh(npt,0:n), AdmApdw(npt,1:n), ImpApdw(npt,1:n), AdmApup(npt,0:n-1), ImpApup(npt,0:n-1)

  if (hordist < eps) then
    r = 1.d-2
  else
    r = hordist
  end if
  kr = krJ0J1 / r
  wvnb2(0) = -zeta * eta0
  u(:,0) = sqrt(kr * kr - wvnb2(0))
  AdmInt(:,0) = u(:,0) / zeta
  ImpInt(:,0) = u(:,0) / eta0
  uh(:,0) = u(:,0) * h(0)
  tgh(:,0) = (1.0 - exp(-2.0 * uh(:,0))) / (1.0 + exp(-2.0 * uh(:,0)))
  do i = 1, n
    wvnb2(i) = -zeta * eta(i)
    u(:,i) = sqrt(kr * kr - wvnb2(i))
    AdmInt(:,i) = u(:,i) / zeta
    ImpInt(:,i) = u(:,i) / eta(i)
    uh(:,i) = u(:,i) * h(i)
    tgh(:,i) = (1.0 - exp(-2.0 * uh(:,i))) / (1.0 + exp(-2.0 * uh(:,i)))
  end do

  AdmApdw(:,n) = AdmInt(:,n)
  ImpApdw(:,n) = ImpInt(:,n)
  RTEdw(:,n) = (0.0,0.0)
  RTMdw(:,n) = (0.0,0.0)
  do i = n-1, 1, -1
    AdmApdw(:,i) = AdmInt(:,i) * (AdmApdw(:,i+1) + AdmInt(:,i) * &
                  tgh(:,i)) / (AdmInt(:,i) + AdmApdw(:,i+1) * tgh(:,i))
    ImpApdw(:,i) = ImpInt(:,i) * (ImpApdw(:,i+1) + ImpInt(:,i) * &
                  tgh(:,i)) / (ImpInt(:,i) + ImpApdw(:,i+1) * tgh(:,i))
    RTEdw(:,i) = (AdmInt(:,i) - AdmApdw(:,i+1)) / (AdmInt(:,i) + AdmApdw(:,i+1))
    RTMdw(:,i) = (ImpInt(:,i) - ImpApdw(:,i+1)) / (ImpInt(:,i) + ImpApdw(:,i+1))
  end do
  RTEdw(:,0) = (AdmInt(:,0) - AdmApdw(:,1)) / (AdmInt(:,0) + AdmApdw(:,1))
  RTMdw(:,0) = (ImpInt(:,0) - ImpApdw(:,1)) / (ImpInt(:,0) + ImpApdw(:,1))

  AdmApup(:,0) = AdmInt(:,0)
  ImpApup(:,0) = ImpInt(:,0)
  RTEup(:,0) = (0.0,0.0)
  RTMup(:,0) = (0.0,0.0)
  do i = 1, n-1
    AdmApup(:,i) = AdmInt(:,i) * (AdmApup(:,i-1) + AdmInt(:,i) * &
                tgh(:,i)) / (AdmInt(:,i) + AdmApup(:,i-1) * tgh(:,i))
    ImpApup(:,i) = ImpInt(:,i) * (ImpApup(:,i-1) + ImpInt(:,i) * &
                tgh(:,i)) / (ImpInt(:,i) + ImpApup(:,i-1) * tgh(:,i))
    RTEup(:,i) = (AdmInt(:,i) - AdmApup(:,i-1)) / (AdmInt(:,i) + AdmApup(:,i-1))
    RTMup(:,i) = (ImpInt(:,i) - ImpApup(:,i-1)) / (ImpInt(:,i) + ImpApup(:,i-1))
  end do
  RTEup(:,n) = (AdmInt(:,n) - AdmApup(:,n-1)) / (AdmInt(:,n) + AdmApup(:,n-1))
  RTMup(:,n) = (ImpInt(:,n) - ImpApup(:,n-1)) / (ImpInt(:,n) + ImpApup(:,n-1))

  den = 1.0 - RTMup(:,camadT) * RTMdw(:,camadT) * exp(-2.0 * uh(:,camadT))
  AMdw = (exp(-u(:,camadT) * (prof(camadT) - h0)) + RTMup(:,camadT) * &
        exp(u(:,camadT) * (prof(camadT-1) - h(camadT)-h0))) / den

  AMup = (exp(u(:,camadT) * (prof(camadT-1)-h0)) + RTMdw(:,camadT) * &
      exp(-u(:,camadT) * (prof(camadT)+h(camadT)-h0))) / den

  den = 1.0 - RTEup(:,camadT) * RTEdw(:,camadT) * exp(-2.0 * uh(:,camadT))
  FEdw = (exp(-u(:,camadT) * (prof(camadT)-h0)) - RTEup(:,camadT) * &
        exp(u(:,camadT) * (prof(camadT-1)-h(camadT)-h0))) / den

  FEup = (exp(u(:,camadT) * (prof(camadT-1)-h0)) - RTEdw(:,camadT) * &
      exp(-u(:,camadT) * (prof(camadT)+h(camadT)-h0))) / den

  if (present(FEdwz) .and. present(FEupz)) then
    FEdwz = (exp(-u(:,camadT) * (prof(camadT)-h0)) + RTEup(:,camadT) * &
          exp(u(:,camadT) * (prof(camadT-1) - h(camadT) - h0))) / den

    FEupz = (exp(u(:,camadT) * (prof(camadT-1)-h0)) + RTEdw(:,camadT) * &
          exp(-u(:,camadT) * (prof(camadT) + h(camadT) - h0))) / den
  end if
end subroutine commonarraysMD

end module utils
