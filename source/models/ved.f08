module ved
   use parameters
contains

   subroutine ved_optimized(Tx, Ty, h0, n, camadR, camadT, npt, krJ0J1, wJ0, wJ1, h, prof, eta, eta0, cx, cy, z, &
      u, uh, ImpInt, RTMdw, RTMup, AMdwz, AMupz, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
      implicit none
      integer,intent(in) :: n, camadT, camadR, npt
      real(dp), intent(in) :: Tx, Ty, h0, h(0:n), prof(-1:n), krJ0J1(npt), wJ0(npt), wJ1(npt), cx, cy, z
      complex(dp), intent(in) :: eta(1:n), eta0, AMdwz(npt), AMupz(npt)
      complex(dp), dimension(npt,0:n), intent(in) :: u, uh, ImpInt, RTMdw, RTMup
      complex(dp), intent(out) :: Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p

      integer :: j
      real(dp) :: x, y, r, kr(npt), twopir

      complex(dp), dimension(npt) :: fac, KtmzJ0, KtmzJ1, KtmdzzJ1
      complex(dp) :: kernelEJ1, kernelEJ0, kernelHJ1
      complex(dp), dimension(:,:), allocatable :: TMdwz, TMupz

      Hz_p = (0.0,0.0)

      if (dabs(cx - Tx) < eps .and. dabs(cy - Ty) < eps) then
         x = 0.0    !assigning null values to very small coordinates
         y = 0.0    !assigning null values to very small coordinates
         r = 1.d-2  !value to avoid division by zero in preliminary field determination steps
      elseif (dabs(cx - Tx) < eps) then
         x = 0.0  !assigning null value to very small abscissa
         y = cy - Ty
         r = dabs(y)
      elseif (dabs(cy - Ty) < eps) then
         x = cx - Tx
         y = 0.0  !assigning null value to very small ordinate
         r = dabs(x)
      else
         x = cx - Tx
         y = cy - Ty
         r = dsqrt( x ** 2 + y ** 2 )
      end if
      kr = krJ0J1 / r

      ! To workaround the warning: ... may be used uninitialized in this function
      allocate(TMdwz(1,1), TMupz(1,1))
      if (camadR > camadT) then
         deallocate(TMdwz)
         allocate(TMdwz(npt,camadT:camadR))
         do j = camadT, camadR
            if (j == camadT) then
               TMdwz(:,j) = Iw * dsz / (2.0 * u(:,camadT))
            elseif (j == (camadT + 1) .and. j == n) then
               TMdwz(:,j) = TMdwz(:,j-1) * (exp(-u(:,camadT) * (prof(camadT)-h0)) + &
                  RTMup(:,camadT) * AMupz(:) * exp(-uh(:,camadT)) + RTMdw(:,camadT) * AMdwz(:))
            elseif (j == (camadT + 1) .and. j /= n) then
               TMdwz(:,j) = TMdwz(:,j-1) * (exp(-u(:,camadT) * (prof(camadT)-h0)) + &
                  RTMup(:,camadT) * AMupz(:) * exp(-uh(:,camadT)) + &
                  RTMdw(:,camadT) * AMdwz(:)) / (1.0 + RTMdw(:,j) * exp(-2.0 * uh(:,j)))
            elseif (j /= n) then
               TMdwz(:,j) = TMdwz(:,j-1)*(1.0+RTMdw(:,j-1))*exp(-uh(:,j-1))/ &
                  (1.0 + RTMdw(:,j)*exp(-2.0*uh(:,j)))
            elseif (j==n) then
               TMdwz(:,j) = TMdwz(:,j-1)*(1.0+RTMdw(:,j-1))*exp(-uh(:,j-1))
            end if
         end do
      elseif (camadR < camadT) then
         deallocate(TMupz)
         allocate(TMupz(npt,camadR:camadT))
         do j=camadT,camadR,-1
            if (j == camadT) then
               TMupz(:,j) = Iw * dsz / (2.0 * u(:,camadT))
            elseif (j == (camadT - 1) .and. j == 0) then
               TMupz(:,j) = TMupz(:,j+1)*(exp(-u(:,camadT)*h0) + &
                  RTMup(:,camadT) * AMupz(:) + RTMdw(:,camadT) * AMdwz(:) * exp(-uh(:,camadT)))
            elseif (j == (camadT - 1) .and. j /= 0) then
               TMupz(:,j) = TMupz(:,j+1)*(exp(u(:,camadT)*(prof(camadT-1)-h0)) + &
                  RTMup(:,camadT)*AMupz(:) + RTMdw(:,camadT)*AMdwz(:)* &
                  exp(-uh(:,camadT))) / (1.0+RTMup(:,j)*exp(-2.0*uh(:,j)))
            elseif (j /= 0) then
               TMupz(:,j) = TMupz(:,j+1)*(1.0+RTMup(:,j+1))*exp(-uh(:,j+1)) / &
                  (1.0 + RTMup(:,j)*exp(-2.0*uh(:,j)))
            elseif (j == 0) then
               TMupz(:,j) = TMupz(:,1)*(1.0+RTMup(:,1))*exp(-uh(:,1))
            end if
         end do
      else
         deallocate(TMdwz, TMupz)
         allocate(TMdwz(npt,camadT:camadR), TMupz(npt,camadR:camadT))
         TMdwz(:,camadR)= Iw * dsz / (2.0*u(:,camadT))
         TMupz(:,camadR)= Iw * dsz / (2.0*u(:,camadT))
      end if

      twopir = 2.0 * pi * r
      if (camadR == 0 .and. camadT /= 0) then
         fac = TMupz(:,0) * exp(u(:,0)*z)
         KtmzJ0 = fac * wJ0
         KtmzJ1 = fac * wJ1
         KtmdzzJ1 = ImpInt(:,0) * fac * wJ1

         kernelEJ1 = sum(KtmdzzJ1 * kr * kr) / r
         Ex_p = -x * kernelEJ1 / twopir

         Ey_p = -y * kernelEJ1 / twopir

         kernelEJ0 = sum(KtmzJ0 * kr * kr * kr)
         Ez_p = kernelEJ0 / twopir / eta0

         kernelHJ1 = sum(KtmzJ1 * kr * kr) / r
         Hx_p = -y * kernelHJ1 / twopir

         Hy_p = x * kernelHJ1 / twopir
      elseif (camadR < camadT) then !camada k
         fac = TMupz(:,camadR) * (exp(u(:,camadR) * (z-prof(camadR))) + &
            RTMup(:,camadR) * exp(-u(:,camadR) * (z-prof(camadR-1) + h(camadR))))
         KtmzJ0 = fac * wJ0(:)
         KtmzJ1 = fac * wJ1(:)
         KtmdzzJ1 = (ImpInt(:,camadR) * TMupz(:,camadR) * (exp(u(:,camadR) * (z-prof(camadR))) - &
            RTMup(:,camadR) * exp(-u(:,camadR) * (z-prof(camadR-1)+h(camadR))))) * wJ1(:)

         kernelEJ1 = sum(KtmdzzJ1 * kr * kr) / r
         Ex_p = -x * kernelEJ1 / twopir

         Ey_p = -y * kernelEJ1 / twopir

         kernelEJ0 = sum(KtmzJ0 * kr * kr * kr)
         Ez_p = kernelEJ0 / twopir / eta(camadR)

         kernelHJ1 = sum(KtmzJ1 * kr * kr) / r
         Hx_p = -y * kernelHJ1 / twopir

         Hy_p = x * kernelHJ1 / twopir
      elseif (camadR == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
         fac = TMupz(:,camadR) * (exp(u(:,camadR) * (z-h0)) + &
            RTMup(:,camadR) * AMupz(:) * exp(-u(:,camadR) * (z-prof(camadR-1))) + &
            RTMdw(:,camadR) * AMdwz(:) * exp(u(:,camadR) * (z-prof(camadR))))
         KtmzJ0 = fac * wJ0(:)
         KtmzJ1 = fac * wJ1(:)
         KtmdzzJ1 = (ImpInt(:,camadR) * TMupz(:,camadR) * (exp(u(:,camadR) * (z-h0)) - &
            RTMup(:,camadR) * AMupz(:) * exp(-u(:,camadR) * (z-prof(camadR-1))) + &
            RTMdw(:,camadR) * AMdwz(:) * exp(u(:,camadR) * (z-prof(camadR))))) * wJ1(:)

         kernelEJ1 = sum(KtmdzzJ1 * kr * kr) / r
         Ex_p = -x * kernelEJ1 / twopir

         Ey_p = -y * kernelEJ1 / twopir

         kernelEJ0 = sum(KtmzJ0 * kr * kr * kr)
         if (camadR /= 0) then
            Ez_p = kernelEJ0 / twopir / eta(camadR)
         else
            Ez_p = kernelEJ0 / twopir / eta0
         end if

         kernelHJ1 = sum(KtmzJ1 * kr * kr) / r
         Hx_p = -y * kernelHJ1 / twopir

         Hy_p = x * kernelHJ1 / twopir
      elseif (camadR == camadT .and. z > h0) then !na mesma camada do transmissor mas abaixo dele
         fac = TMdwz(:,camadR) * (exp(-u(:,camadR) * (z-h0)) + &
            RTMup(:,camadR) * AMupz(:) * exp(-u(:,camadR) * (z-prof(camadR-1))) + &
            RTMdw(:,camadR) * AMdwz(:) * exp(u(:,camadR) * (z-prof(camadR))))
         KtmzJ0 = fac * wJ0(:)
         KtmzJ1 = fac * wJ1(:)
         KtmdzzJ1 = (ImpInt(:,camadR) * TMdwz(:,camadR) * (exp(-u(:,camadR) * (z-h0)) + &
            RTMup(:,camadR) * AMupz(:) * exp(-u(:,camadR) * (z-prof(camadR-1))) - &
            RTMdw(:,camadR) * AMdwz(:) * exp(u(:,camadR) * (z-prof(camadR))))) * wJ1(:)

         kernelEJ1 = sum(KtmdzzJ1 * kr * kr) / r
         Ex_p = x * kernelEJ1 / twopir

         Ey_p = y * kernelEJ1 / twopir

         kernelEJ0 = sum(KtmzJ0 * kr * kr * kr)
         if (camadR /= 0) then
            Ez_p = kernelEJ0 / twopir / eta(camadR)
         else
            Ez_p = kernelEJ0 / twopir / eta0
         end if

         kernelHJ1 = sum(KtmzJ1 * kr * kr) / r
         Hx_p = -y * kernelHJ1 / twopir

         Hy_p = x * kernelHJ1 / twopir
      elseif (camadR > camadT .and. camadR /= n) then !camada j
         fac = TMdwz(:,camadR) * (exp(-u(:,camadR) * (z-prof(camadR-1))) + &
            RTMdw(:,camadR) * exp(u(:,camadR) * (z-prof(camadR)-h(camadR))))
         KtmzJ0 = fac * wJ0(:)
         KtmzJ1 = fac * wJ1(:)
         KtmdzzJ1=(ImpInt(:,camadR) * TMdwz(:,camadR) * (exp(-u(:,camadR) * (z-prof(camadR-1))) - &
            RTMdw(:,camadR) * exp(u(:,camadR) * (z-prof(camadR)-h(camadR)))))*wJ1(:)

         kernelEJ1 = sum(KtmdzzJ1 * kr * kr) / r
         Ex_p = x * kernelEJ1 / twopir

         Ey_p = y * kernelEJ1 / twopir

         kernelEJ0 = sum(KtmzJ0 * kr * kr * kr)
         Ez_p = kernelEJ0 / twopir / eta(camadR)

         kernelHJ1 = sum(KtmzJ1 * kr * kr) / r
         Hx_p = -y * kernelHJ1 / twopir

         Hy_p = x * kernelHJ1 / twopir
      else !camada n
         fac = TMdwz(:,n) * exp(-u(:,n) * (z-prof(n-1)))
         KtmzJ0 = fac * wJ0(:)
         KtmzJ1 = fac * wJ1(:)
         KtmdzzJ1 = ImpInt(:,n) * fac * wJ1(:)

         kernelEJ1 = sum(KtmdzzJ1 * kr * kr) / r
         Ex_p = x * kernelEJ1 / twopir

         Ey_p = y * kernelEJ1 / twopir

         kernelEJ0 = sum(KtmzJ0 * kr * kr * kr)
         Ez_p = kernelEJ0 / twopir / eta(camadR)

         kernelHJ1 = sum(KtmzJ1 * kr * kr) / r
         Hx_p = -y * kernelHJ1 / twopir

         Hy_p = x * kernelHJ1 / twopir
      end if
   end subroutine ved_optimized

end module ved
