module hedx
use parameters
contains

subroutine hedx_optimized(Tx, Ty, h0, n, camadR, camadT, npt, krJ0J1, wJ0, wJ1, h, prof, eta, zeta, eta0, cx, cy, z, &
                            u, uh, AdmInt, ImpInt, RTEdw, RTEup, RTMdw, RTMup, FEdw, FEup, AMdw, AMup, &
                            Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
  implicit none
  integer,intent(in) :: n, camadT, camadR, npt
  real(dp), intent(in) :: Tx, Ty, h0, h(0:n), prof(-1:n), krJ0J1(npt), wJ0(npt), wJ1(npt), cx, cy, z
  complex(dp), intent(in) :: eta(1:n), zeta, eta0, AMdw(npt), AMup(npt), FEdw(npt), FEup(npt)
  complex(dp), dimension(npt,0:n), intent(in) :: u, uh, AdmInt, ImpInt, RTEdw, RTMdw, RTEup, RTMup
  complex(dp), intent(out) :: Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p

  integer :: j
  real(dp) :: x, y, r, kr(npt), x2_r2, y2_r2, xy_r2, twox2_r2m1, twoy2_r2m1, twopir
  complex(dp), dimension(npt) :: Ktmdz, Ktmdz_J0, Ktmdz_J1, Ktm, Ktm_J0, Ktm_J1
  complex(dp), dimension(npt) :: Kte, Kte_J0, Kte_J1, Ktedz, Ktedz_J0, Ktedz_J1
  complex(dp) :: kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1, kernelEzJ1
  complex(dp) :: kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1
  complex(dp), dimension(:,:), allocatable :: TMdw, TEdw, TMup, TEup

  if ( dabs(cx - Tx) < eps ) then
    x = 0.0
  else
    x = cx - Tx
  end if
  if ( dabs(cy - Ty) < eps ) then
    y = 0.0
  else
    y = cy - Ty
  end if
  r = dsqrt( x ** 2 + y ** 2 )

  kr = krJ0J1 / r

  ! To workaround the warning: ... may be used uninitialized in this function
  allocate(TMdw(1,1), TEdw(1,1))
  allocate(TMup(1,1), TEup(1,1))

  if (camadR > camadT) then
    deallocate(TMdw, TEdw)
    allocate(TMdw(npt,camadT:camadR),TEdw(npt,camadT:camadR))
    do j=camadT,camadR
      if (j == camadT) then
        TMdw(:,j)= - Iw * dsx / 2.0
        TEdw(:,j)= - Iw * dsx / (2.0 * AdmInt(:,camadT))
      elseif (j == (camadT + 1) .and. j == n) then
        TMdw(:,j)=TMdw(:,j-1)*(exp(-u(:,camadT)*(prof(camadT)-h0)) - &
                    RTMup(:,camadT)*AMup(:)*exp(-uh(:,camadT)) + &
                    RTMdw(:,camadT)*AMdw(:))
        TEdw(:,j)=TEdw(:,j-1)*(exp(-u(:,camadT)*(prof(camadT)-h0)) + &
                    RTEup(:,camadT)*FEup(:)*exp(-uh(:,camadT)) + &
                    RTEdw(:,camadT)*FEdw(:))
      elseif (j == (camadT + 1) .and. j /= n) then
        TMdw(:,j)=TMdw(:,j-1)*(exp(-u(:,camadT)*(prof(camadT)-h0)) - &
                    RTMup(:,camadT)*AMup(:)*exp(-uh(:,camadT)) + &
                    RTMdw(:,camadT)*AMdw(:)) / (1.0 + RTMdw(:,j)*exp(-2.0*uh(:,j)))
        TEdw(:,j)=TEdw(:,j-1)*(exp(-u(:,camadT)*(prof(camadT)-h0)) + &
                    RTEup(:,camadT)*FEup(:)*exp(-uh(:,camadT)) + &
                    RTEdw(:,camadT)*FEdw(:)) / (1.0 + RTEdw(:,j)*exp(-2.0*uh(:,j)))
      elseif (j /= n) then
        TMdw(:,j)=TMdw(:,j-1)*(1.0+RTMdw(:,j-1))*exp(-uh(:,j-1))/ &
            (1.0 + RTMdw(:,j)*exp(-2.0*uh(:,j)))
        TEdw(:,j)=TEdw(:,j-1)*(1.0+RTEdw(:,j-1))*exp(-uh(:,j-1))/ &
            (1.0 + RTEdw(:,j)*exp(-2.0*uh(:,j)))
      elseif (j==n) then
        TMdw(:,j)=TMdw(:,j-1)*(1.0+RTMdw(:,j-1))*exp(-uh(:,j-1))
        TEdw(:,j)=TEdw(:,j-1)*(1.0+RTEdw(:,j-1))*exp(-uh(:,j-1))
      end if
    end do
  elseif (camadR < camadT) then
    deallocate(TMup, TEup)
    allocate(TMup(npt,camadR:camadT),TEup(npt,camadR:camadT))
    do j=camadT,camadR,-1
      if (j == camadT) then
        TMup(:,j)= Iw * dsx / 2.0
        TEup(:,j)= - Iw * dsx / (2.0 * AdmInt(:,camadT))
      elseif (j == (camadT - 1) .and. j == 0) then
        TMup(:,j)=TMup(:,j+1)*(exp(-u(:,camadT)*h0) + &
                    RTMup(:,camadT)*AMup(:) - RTMdw(:,camadT)*AMdw(:)* &
                    exp(-uh(:,camadT)))
        TEup(:,j)=TEup(:,j+1)*(exp(-u(:,camadT)*h0) + &
                    RTEup(:,camadT)*FEup(:) + RTEdw(:,camadT)*FEdw(:)* &
                    exp(-uh(:,camadT)))
      elseif (j == (camadT - 1) .and. j /= 0) then
        TMup(:,j)=TMup(:,j+1)*(exp(u(:,camadT)*(prof(camadT-1)-h0)) + &
                    RTMup(:,camadT)*AMup(:) - RTMdw(:,camadT)*AMdw(:)* &
                    exp(-uh(:,camadT))) / (1.0+RTMup(:,j)*exp(-2.0*uh(:,j)))
        TEup(:,j)=TEup(:,j+1)*(exp(u(:,camadT)*(prof(camadT-1)-h0)) + &
                    RTEup(:,camadT)*FEup(:) + RTEdw(:,camadT)*FEdw(:)* &
                    exp(-uh(:,camadT))) / (1.0+RTEup(:,j)*exp(-2.0*uh(:,j)))
      elseif (j /= 0) then
        TMup(:,j)=TMup(:,j+1)*(1.0+RTMup(:,j+1))*exp(-uh(:,j+1)) / &
                    (1.0 + RTMup(:,j)*exp(-2.0*uh(:,j)))
        TEup(:,j)=TEup(:,j+1)*(1.0+RTEup(:,j+1))*exp(-uh(:,j+1)) / &
                    (1.0 + RTEup(:,j)*exp(-2.0*uh(:,j)))
      elseif (j == 0) then
        TMup(:,j)=TMup(:,1)*(1.0+RTMup(:,1))*exp(-uh(:,1))
        TEup(:,j)=TEup(:,1)*(1.0+RTEup(:,1))*exp(-uh(:,1))
      end if
    end do
  else
    deallocate(TMdw, TEdw, TMup, TEup)
    allocate(TMdw(npt,camadT:camadR), TEdw(npt,camadT:camadR))
    allocate(TMup(npt,camadR:camadT), TEup(npt,camadR:camadT))
    TMdw(:,camadR)= - Iw * dsx / 2.0
    TEdw(:,camadR)= - Iw * dsx / (2.0 * AdmInt(:,camadT))
    TMup(:,camadR)= Iw * dsx / 2.0
    TEup(:,camadR)= - Iw * dsx / (2.0 * AdmInt(:,camadT))
  end if

  x2_r2 = x * x / r / r
  y2_r2 = y * y / r / r
  xy_r2 = x * y / r / r
  twox2_r2m1 = 2.0 * x2_r2 - 1.0
  twoy2_r2m1 = 2.0 * y2_r2 - 1.0
  twopir = 2.0 * pi * r
  if (camadR == 0 .and. camadT /= 0) then
    Kte = TEup(:,0) * exp(u(:,0)*z)
    Kte_J0 = Kte * wJ0(:)
    Kte_J1 = Kte * wJ1(:)
    Ktedz = AdmInt(:,0) * Kte
    Ktedz_J0 = Ktedz * wJ0(:)
    Ktedz_J1 = Ktedz * wJ1(:)
    Ktm = TMup(:,0) * exp(u(:,0)*z)
    Ktm_J0 = Ktm * wJ0(:)
    Ktm_J1 = Ktm * wJ1(:)
    Ktmdz = ImpInt(:,0) * Ktm
    Ktmdz_J0 = Ktmdz * wJ0(:)
    Ktmdz_J1 = Ktmdz * wJ1(:)

    kernelExJ1 = (twox2_r2m1 * sum(Ktmdz_J1) - twoy2_r2m1 * sum(Kte_J1)) / r
    kernelExJ0 = sum((x2_r2 * Ktmdz_J0 - y2_r2 * Kte_J0) * kr)
    Ex_p = (kernelExJ1 - kernelExJ0) / twopir

    kernelEyJ1 = 2.0 * xy_r2 * sum(Ktmdz_J1 + Kte_J1) / r
    kernelEyJ0 = xy_r2 * sum((Ktmdz_J0 + Kte_J0) * kr)
    Ey_p = (kernelEyJ1 - kernelEyJ0) / twopir

    kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta0
    Ez_p = - kernelEzJ1 / twopir

    kernelHxJ1 = 2.0 * xy_r2 * sum(Ktm_J1 + Ktedz_J1) / r
    kernelHxJ0 = xy_r2 * sum((Ktm_J0 + Ktedz_J0) * kr)
    Hx_p = (kernelHxJ1 - kernelHxJ0) / twopir

    kernelHyJ1 = (-twox2_r2m1 * sum(Ktm_J1) + twoy2_r2m1 * sum(Ktedz_J1)) / r
    kernelHyJ0 = sum((x2_r2 * Ktm_J0 - y2_r2 * Ktedz_J0) * kr)
    Hy_p = (kernelHyJ1 + kernelHyJ0) / twopir

    kernelHzJ1 = y * sum(Kte_J1 * kr * kr) / r / zeta
    Hz_p = -kernelHzJ1 / twopir
  elseif (camadR < camadT) then !camada k
    Kte = TEup(:,camadR)*(exp(u(:,camadR)*(z-prof(camadR))) + &
          RTEup(:,camadR)*exp(-u(:,camadR)*(z-prof(camadR-1)+h(camadR))))
    Kte_J0 = Kte * wJ0(:)
    Kte_J1 = Kte * wJ1(:)
    Ktedz = AdmInt(:,camadR)*TEup(:,camadR)*(exp(u(:,camadR)*(z-prof(camadR))) - &
          RTEup(:,camadR)*exp(-u(:,camadR)*(z-prof(camadR-1)+h(camadR))))
    Ktedz_J0 = Ktedz * wJ0(:)
    Ktedz_J1 = Ktedz * wJ1(:)
    Ktm = TMup(:,camadR)*(exp(u(:,camadR)*(z-prof(camadR))) + &
          RTMup(:,camadR)*exp(-u(:,camadR)*(z-prof(camadR-1)+h(camadR))))
    Ktm_J0 = Ktm * wJ0(:)
    Ktm_J1 = Ktm * wJ1(:)
    Ktmdz = ImpInt(:,camadR)*TMup(:,camadR)*(exp(u(:,camadR)*(z-prof(camadR))) - &
            RTMup(:,camadR)*exp(-u(:,camadR)*(z-prof(camadR-1)+h(camadR))))
    Ktmdz_J0 = Ktmdz * wJ0(:)
    ktmdz_J1 = Ktmdz * wJ1(:)

    kernelExJ1 = (twox2_r2m1 * sum(Ktmdz_J1) - twoy2_r2m1 * sum(Kte_J1)) / r
    kernelExJ0 = sum((x2_r2 * Ktmdz_J0 - y2_r2 * Kte_J0) * kr)
    Ex_p = (kernelExJ1 - kernelExJ0) / twopir

    kernelEyJ1 = 2.0 * xy_r2 * sum(Ktmdz_J1 + Kte_J1) / r
    kernelEyJ0 = xy_r2 * sum((Ktmdz_J0 + Kte_J0) * kr)
    Ey_p = (kernelEyJ1 - kernelEyJ0) / twopir

    kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta(camadR)
    Ez_p = -kernelEzJ1 / twopir

    kernelHxJ1 = 2.0 * xy_r2 * sum(Ktm_J1 + Ktedz_J1) / r
    kernelHxJ0 = xy_r2 * sum((Ktm_J0 + Ktedz_J0) * kr)
    Hx_p = (kernelHxJ1 - kernelHxJ0) / twopir

    kernelHyJ1 = (-twox2_r2m1 * sum(Ktm_J1) + twoy2_r2m1 * sum(Ktedz_J1)) / r
    kernelHyJ0 = sum((x2_r2 * Ktm_J0 - y2_r2 * Ktedz_J0) * kr)
    Hy_p = (kernelHyJ1 + kernelHyJ0) / twopir

    kernelHzJ1 = y * sum(Kte_J1 * kr * kr) / r / zeta
    Hz_p = -kernelHzJ1 / twopir
  elseif (camadR == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
    Kte = TEup(:,camadR)*(exp(u(:,camadR)*(z-h0)) + &
          RTEup(:,camadR)*FEup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
          RTEdw(:,camadR)*FEdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Kte_J0 = Kte * wJ0(:)
    Kte_J1 = Kte * wJ1(:)
    Ktedz = AdmInt(:,camadR)*TEup(:,camadR)*(exp(u(:,camadR)*(z-h0)) - &
            RTEup(:,camadR)*FEup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
            RTEdw(:,camadR)*FEdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Ktedz_J0 = Ktedz * wJ0(:)
    Ktedz_J1 = Ktedz * wJ1(:)
    Ktm = TMup(:,camadR)*(exp(u(:,camadR)*(z-h0)) + &
          RTMup(:,camadR)*AMup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) - &
          RTMdw(:,camadR)*AMdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Ktm_J0 = Ktm * wJ0(:)
    Ktm_J1 = Ktm * wJ1(:)
    Ktmdz = ImpInt(:,camadR)*TMup(:,camadR)*(exp(u(:,camadR)*(z-h0)) - &
            RTMup(:,camadR)*AMup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) - &
            RTMdw(:,camadR)*AMdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Ktmdz_J0 = Ktmdz * wJ0(:)
    Ktmdz_J1 = Ktmdz * wJ1(:)

    kernelExJ1 = (twox2_r2m1 * sum(Ktmdz_J1) - twoy2_r2m1 * sum(Kte_J1)) / r
    kernelExJ0 = sum((x2_r2 * Ktmdz_J0 - y2_r2 * Kte_J0) * kr)
    Ex_p = (kernelExJ1 - kernelExJ0) / twopir

    kernelEyJ1 = 2.0 * xy_r2 * sum(Ktmdz_J1 + Kte_J1) / r
    kernelEyJ0 = xy_r2 * sum((Ktmdz_J0 + Kte_J0) * kr)
    Ey_p = (kernelEyJ1 - kernelEyJ0) / twopir

    if (camadR /= 0) then
      kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta(camadR)
    else
      kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta0
    end if
    Ez_p = - kernelEzJ1 / twopir

    kernelHxJ1 = 2.0 * xy_r2 * sum(Ktm_J1 + Ktedz_J1) / r
    kernelHxJ0 = xy_r2 * sum((Ktm_J0 + Ktedz_J0) * kr)
    Hx_p = (kernelHxJ1 - kernelHxJ0) / twopir

    kernelHyJ1 = (-twox2_r2m1 * sum(Ktm_J1) + twoy2_r2m1 * sum(Ktedz_J1)) / r
    kernelHyJ0 = sum((x2_r2 * Ktm_J0 - y2_r2 * Ktedz_J0) * kr)
    Hy_p = (kernelHyJ1 + kernelHyJ0) / twopir

    kernelHzJ1 = y * sum(Kte_J1 * kr * kr) / r / zeta
    Hz_p = -kernelHzJ1 / twopir
  elseif (camadR == camadT .and. z > h0) then  !na mesma camada do transmissor mas abaixo dele
    Kte = TEdw(:,camadR)*(exp(-u(:,camadR)*(z-h0)) + &
          RTEup(:,camadR)*FEup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
          RTEdw(:,camadR)*FEdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Kte_J0 = Kte * wJ0(:)
    Kte_J1 = Kte * wJ1(:)
    Ktedz = AdmInt(:,camadR)*TEdw(:,camadR)*(exp(-u(:,camadR)*(z-h0)) + &
            RTEup(:,camadR)*FEup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) - &
            RTEdw(:,camadR)*FEdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Ktedz_J0 = Ktedz * wJ0(:)
    Ktedz_J1 = Ktedz * wJ1(:)
    Ktm = TMdw(:,camadR)*(exp(-u(:,camadR)*(z-h0)) - &
          RTMup(:,camadR)*AMup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
          RTMdw(:,camadR)*AMdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Ktm_J0 = Ktm * wJ0(:)
    Ktm_J1 = Ktm * wJ1(:)
    Ktmdz = ImpInt(:,camadR)*TMdw(:,camadR)*(exp(-u(:,camadR)*(z-h0)) - &
            RTMup(:,camadR)*AMup(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) - &
            RTMdw(:,camadR)*AMdw(:)*exp(u(:,camadR)*(z-prof(camadR))))
    Ktmdz_J0 = Ktmdz * wJ0(:)
    Ktmdz_J1 = Ktmdz * wJ1(:)

    kernelExJ1 = -(twox2_r2m1 * sum(Ktmdz_J1) + twoy2_r2m1 * sum(Kte_J1)) / r
    kernelExJ0 = sum((x2_r2 * Ktmdz_J0 + y2_r2 * Kte_J0) * kr)
    Ex_p = (kernelExJ1 + kernelExJ0) / twopir

    kernelEyJ1 = 2.0 * xy_r2 * sum(-Ktmdz_J1 + Kte_J1) / r
    kernelEyJ0 = xy_r2 * sum((-Ktmdz_J0 + Kte_J0) * kr)
    Ey_p = (kernelEyJ1 - kernelEyJ0) / twopir

    if (camadR /= 0) then
      kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta(camadR)
    else
      kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta0
    end if
    Ez_p = - kernelEzJ1 / twopir

    kernelHxJ1 = 2.0 * xy_r2 * sum(Ktm_J1 - Ktedz_J1) / r
    kernelHxJ0 = xy_r2 * sum((Ktm_J0 - Ktedz_J0) * kr)
    Hx_p = (kernelHxJ1 - kernelHxJ0) / twopir

    kernelHyJ1 = -(twox2_r2m1 * sum(Ktm_J1) + twoy2_r2m1 * sum(Ktedz_J1)) / r
    kernelHyJ0 = sum((x2_r2 * Ktm_J0 + y2_r2 * Ktedz_J0) * kr)
    Hy_p = (kernelHyJ1 + kernelHyJ0) / twopir

    kernelHzJ1 = y * sum(Kte_J1 * kr * kr) / r / zeta
    Hz_p = -kernelHzJ1 / twopir
  elseif (camadR > camadT .and. camadR /= n) then !camada j
    Kte = TEdw(:,camadR)*(exp(-u(:,camadR)*(z-prof(camadR-1))) + &
          RTEdw(:,camadR)*exp(u(:,camadR)*(z-prof(camadR)-h(camadR))))
    Kte_J0 = Kte * wJ0(:)
    Kte_J1 = Kte * wJ1(:)
    Ktedz = AdmInt(:,camadR)*TEdw(:,camadR)*(exp(-u(:,camadR)*(z-prof(camadR-1))) - &
          RTEdw(:,camadR)*exp(u(:,camadR)*(z-prof(camadR)-h(camadR))))
    Ktedz_J0 = Ktedz * wJ0(:)
    Ktedz_J1 = Ktedz * wJ1(:)
    Ktm = TMdw(:,camadR)*(exp(-u(:,camadR)*(z-prof(camadR-1))) + &
          RTMdw(:,camadR)*exp(u(:,camadR)*(z-prof(camadR)-h(camadR))))
    Ktm_J0 = Ktm * wJ0(:)
    Ktm_J1 = Ktm * wJ1(:)
    Ktmdz = ImpInt(:,camadR)*TMdw(:,camadR)*(exp(-u(:,camadR)*(z-prof(camadR-1))) - &
            RTMdw(:,camadR)*exp(u(:,camadR)*(z-prof(camadR)-h(camadR))))
    Ktmdz_J0 = Ktmdz * wJ0(:)
    Ktmdz_J1 = Ktmdz * wJ1(:)

    kernelExJ1 = -(twox2_r2m1 * sum(Ktmdz_J1) + twoy2_r2m1 * sum(Kte_J1)) / r
    kernelExJ0 = sum((x2_r2 * Ktmdz_J0 + y2_r2 * Kte_J0) * kr)
    Ex_p = (kernelExJ1 + kernelExJ0) / twopir

    kernelEyJ1 = 2.0 * xy_r2 * sum(-Ktmdz_J1 + Kte_J1) / r
    kernelEyJ0 = xy_r2 * sum((-Ktmdz_J0 + Kte_J0) * kr)
    Ey_p = (kernelEyJ1 - kernelEyJ0) / twopir

    kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta(camadR)
    Ez_p = -kernelEzJ1 / twopir

    kernelHxJ1 = 2.0 * xy_r2 * sum(Ktm_J1 - Ktedz_J1) / r
    kernelHxJ0 = xy_r2 * sum((Ktm_J0 - Ktedz_J0) * kr)
    Hx_p = (kernelHxJ1 - kernelHxJ0) / twopir

    kernelHyJ1 = -(twox2_r2m1 * sum(Ktm_J1) + twoy2_r2m1 * sum(Ktedz_J1)) / r
    kernelHyJ0 = sum((x2_r2 * Ktm_J0 + y2_r2 * Ktedz_J0) * kr)
    Hy_p = (kernelHyJ1 + kernelHyJ0) / twopir

    kernelHzJ1 = y * sum(Kte_J1 * kr * kr) / r / zeta
    Hz_p = -kernelHzJ1 / twopir
  else  !camada n
    Kte = TEdw(:,n)*exp(-u(:,n)*(z-prof(n-1)))
    Kte_J0 = Kte * wJ0(:)
    Kte_J1 = Kte * wJ1(:)
    Ktedz_J0 = AdmInt(:,n) * Kte * wJ0(:)
    Ktedz_J1 = AdmInt(:,n) * Kte * wJ1(:)
    Ktm = TMdw(:,n)*exp(-u(:,n)*(z-prof(n-1)))
    Ktm_J0 = Ktm * wJ0(:)
    Ktm_J1 = Ktm * wJ1(:)
    Ktmdz = ImpInt(:,n) * Ktm
    Ktmdz_J0 = Ktmdz * wJ0(:)
    Ktmdz_J1 = Ktmdz * wJ1(:)

    kernelExJ1 = -(twox2_r2m1 * sum(Ktmdz_J1) + twoy2_r2m1 * sum(Kte_J1)) / r
    kernelExJ0 = sum((x2_r2 * Ktmdz_J0 + y2_r2 * Kte_J0) * kr)
    Ex_p = (kernelExJ1 + kernelExJ0) / twopir

    kernelEyJ1 = 2.0 * xy_r2 * sum(-Ktmdz_J1 + Kte_J1) / r
    kernelEyJ0 = xy_r2 * sum((-Ktmdz_J0 + Kte_J0) * kr)
    Ey_p = (kernelEyJ1 - kernelEyJ0) / twopir

    kernelEzJ1 = x * sum(Ktm_J1 * kr * kr) / r / eta(n)
    Ez_p = - kernelEzJ1 / twopir

    kernelHxJ1 = 2.0 * xy_r2 * sum(Ktm_J1 - Ktedz_J1) / r
    kernelHxJ0 = xy_r2 * sum((Ktm_J0 - Ktedz_J0) * kr)
    Hx_p = (kernelHxJ1 - kernelHxJ0) / twopir

    kernelHyJ1 = -(twox2_r2m1 * sum(Ktm_J1) + twoy2_r2m1 * sum(Ktedz_J1)) / r
    kernelHyJ0 = sum((x2_r2 * Ktm_J0 + y2_r2 * Ktedz_J0) * kr)
    Hy_p = (kernelHyJ1 + kernelHyJ0) / twopir

    kernelHzJ1 = y * sum(Kte_J1 * kr * kr) / r / zeta
    Hz_p = -kernelHzJ1 / twopir
  end if
end subroutine hedx_optimized
end module hedx
