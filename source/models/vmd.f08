module vmd
   use parameters
contains

   subroutine vmd_optimized(Tx, Ty, h0, n, camadR, camadT, npt, krJ0J1, wJ0, wJ1, h, prof, zeta, cx, cy, z, &
      u, uh, AdmInt, RTEdw, RTEup, FEdwz, FEupz, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p)
      implicit none
      integer, intent(in) :: n, camadT, camadR, npt
      real(dp), intent(in) :: Tx, Ty, h0, h(0:n), prof(-1:n), krJ0J1(npt), wJ0(npt), wJ1(npt), cx, cy, z
      complex(dp), intent(in) :: zeta, FEdwz(npt), FEupz(npt)
      complex(dp), dimension(npt,0:n), intent(in) :: u, uh, AdmInt, RTEdw, RTEup
      complex(dp), intent(out) :: Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p

      integer :: j
      real(dp) :: x, y, r, kr(npt), twopir

      complex(dp), dimension(npt) :: fac, KtezJ0, KtezJ1, KtedzzJ1
      complex(dp) :: kernelEx, kernelEy, kernelHx, kernelHy, kernelHz
      complex(dp), dimension(:,:), allocatable :: TEdwz, TEupz

      Ez_p = (0.0,0.0)

      if (dabs(cx - Tx) < eps .and. dabs(cy - Ty) < eps) then
         ! x = dsign(1.d-2,cx)
         x = 0.0    !assigning null values to very small coordinates
         ! y = dsign(1.d-2,cy)
         y = 0.0    !assigning null values to very small coordinates
         r = 1.d-2  !value to avoid division by zero in preliminary field determination steps
      elseif (dabs(cx - Tx) < eps) then
         x = dsign(1.d-2,cx)
         y = cy - Ty
         r = dabs(y)
      elseif (dabs(cy - Ty) < eps) then
         x = cx - Tx
         y = dsign(1.d-2,cy)
         r = dabs(x)
      else
         x = cx - Tx
         y = cy - Ty
         r = dsqrt( x ** 2 + y ** 2 )
      end if
      kr = krJ0J1 / r

      ! To workaround the warning: ... may be used uninitialized in this function
      allocate(TEdwz(1,1), TEupz(1,1))
      if (camadR > camadT) then
         deallocate(TEdwz)
         allocate(TEdwz(npt,camadT:camadR))
         do j = camadT, camadR
            if (j == camadT) then
               TEdwz(:,j) = zeta * mz / ( 2.0 * u(:,j) )
            elseif (j == (camadT + 1) .and. j == n) then
               TEdwz(:,j) = TEdwz(:,j-1)*(exp(-u(:,camadT)*(prof(camadT)-h0)) + &
                  RTEup(:,camadT)*FEupz*exp(-uh(:,camadT))+RTEdw(:,camadT)*FEdwz)
            elseif (j == (camadT + 1) .and. j /= n) then
               TEdwz(:,j) = TEdwz(:,j-1)*(exp(-u(:,camadT)*(prof(camadT)-h0)) + &
                  RTEup(:,camadT)*FEupz(:)*exp(-uh(:,camadT)) + &
                  RTEdw(:,camadT)*FEdwz(:))/(1.0+RTEdw(:,j)*exp(-2.0*uh(:,j)))
            elseif (j /= n) then
               TEdwz(:,j) = TEdwz(:,j-1)*(1.0+RTEdw(:,j-1))*exp(-uh(:,j-1)) / &
                  (1.0+RTEdw(:,j)*exp(-2.0*uh(:,j)))
            elseif (j==n) then
               TEdwz(:,j) = TEdwz(:,j-1)*(1.0+RTEdw(:,j-1))*exp(-uh(:,j-1))
            end if
         end do
      elseif (camadR < camadT) then
         deallocate(TEupz)
         allocate(TEupz(npt,camadR:camadT))
         do j=camadT,camadR,-1
            if (j == camadT) then
               TEupz(:,j) = zeta * mz / (2.0 * u(:,j))
            elseif (j == (camadT - 1) .and. j == 0) then
               TEupz(:,j) = TEupz(:,j+1)*(exp(-u(:,camadT)*h0) + &
                  RTEup(:,camadT)*FEupz(:)+RTEdw(:,camadT)*FEdwz(:)*exp(-uh(:,camadT)))
            elseif (j == (camadT - 1) .and. j /= 0) then
               TEupz(:,j) = TEupz(:,j+1)*(exp(u(:,camadT)*(prof(camadT-1)-h0)) + &
                  RTEup(:,camadT)*FEupz(:)+RTEdw(:,camadT)*FEdwz(:) * &
                  exp(-uh(:,camadT)))/(1.0+RTEup(:,j)*exp(-2.0*uh(:,j)))
            elseif (j /= 0) then
               TEupz(:,j) = TEupz(:,j+1)*(1.0+RTEup(:,j+1))*exp(-uh(:,j+1)) / &
                  (1.0+RTEup(:,j)*exp(-2.0*uh(:,j)))
            elseif (j == 0) then
               TEupz(:,j) = TEupz(:,1) * (1.0+RTEup(:,1))*exp(-uh(:,1))
            end if
         end do
      else
         deallocate(TEdwz, TEupz)
         allocate(TEdwz(npt,camadT:camadR), TEupz(npt,camadR:camadT))
         TEdwz(:,camadR) = zeta * mz / (2.0 * u(:,camadT))
         TEupz(:,camadR) = TEdwz(:,camadR)
      end if

      twopir = 2.0 * pi * r
      if (camadR == 0 .and. camadT /= 0) then
         fac = TEupz(:,0)*exp(u(:,0)*z)
         KtezJ0 = fac * wJ0
         KtezJ1 = fac * wJ1
         KtedzzJ1 = AdmInt(:,0) * KtezJ1

         kernelEx = y * sum(KtezJ1 * kr * kr) / twopir
         Ex_p = kernelEx / r

         kernelEy = - x * sum(KtezJ1 * kr * kr) / twopir
         Ey_p = kernelEy / r

         kernelHx = - x * sum(KtedzzJ1 * kr * kr) / twopir
         Hx_p = kernelHx / r

         kernelHy = - y * sum(KtedzzJ1 * kr * kr) / twopir
         Hy_p = kernelHy / r

         kernelHz = sum(KtezJ0 * kr * kr * kr) / 2.0 / pi / zeta
         Hz_p = kernelHz / r
      elseif (camadR < camadT) then !camada k
         fac = TEupz(:,camadR)*(exp(u(:,camadR)*(z-prof(camadR))) + &
            RTEup(:,camadR)*exp(-u(:,camadR)*(z-prof(camadR-1)+h(camadR))))
         KtezJ0 = fac * wJ0
         KtezJ1 = fac * wJ1
         KtedzzJ1 = (AdmInt(:,camadR) * TEupz(:,camadR)*(exp(u(:,camadR)*(z-prof(camadR))) - &
            RTEup(:,camadR)*exp(-u(:,camadR)*(z-prof(camadR-1)+h(camadR))))) * wJ1

         kernelEx = y * sum(KtezJ1 * kr * kr) / twopir
         Ex_p = kernelEx / r

         kernelEy = - x * sum(KtezJ1 * kr * kr) / twopir
         Ey_p = kernelEy / r

         kernelHx = - x * sum(KtedzzJ1 * kr * kr) / twopir
         Hx_p = kernelHx / r

         kernelHy = - y * sum(KtedzzJ1 * kr * kr) / twopir
         Hy_p = kernelHy / r

         kernelHz = sum(KtezJ0 * kr * kr * kr) / 2.0 / pi / zeta
         Hz_p = kernelHz / r
      elseif (camadR == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
         fac = TEupz(:,camadR)*(exp(u(:,camadR)*(z-h0)) + &
            RTEup(:,camadR)*FEupz(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
            RTEdw(:,camadR)*FEdwz(:)*exp(u(:,camadR)*(z-prof(camadR))))
         KtezJ0 = fac * wJ0
         KtezJ1 = fac * wJ1
         KtedzzJ1 = (AdmInt(:,camadR)*TEupz(:,camadR)*(exp(u(:,camadR)*(z-h0)) - &
            RTEup(:,camadR)*FEupz(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
            RTEdw(:,camadR)*FEdwz(:)*exp(u(:,camadR)*(z-prof(camadR))))) * wJ1

         kernelEx = y * sum(KtezJ1 * kr * kr) / twopir
         Ex_p = kernelEx / r

         kernelEy = - x * sum(KtezJ1 * kr * kr) / twopir
         Ey_p = kernelEy / r

         kernelHx = - x * sum(KtedzzJ1 * kr * kr) / twopir
         Hx_p = kernelHx / r

         kernelHy = - y * sum(KtedzzJ1 * kr * kr) / twopir
         Hy_p = kernelHy / r

         kernelHz = sum(KtezJ0 * kr * kr * kr) / 2.0 / pi / zeta
         Hz_p = kernelHz / r
      elseif (camadR == camadT .and. z > h0) then !na mesma camada do transmissor mas abaixo dele
         fac = TEdwz(:,camadR)*(exp(-u(:,camadR)*(z-h0)) + &
            RTEup(:,camadR)*FEupz(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) + &
            RTEdw(:,camadR)*FEdwz(:)*exp(u(:,camadR)*(z-prof(camadR))))
         KtezJ0 = fac * wJ0
         KtezJ1 = fac * wJ1
         KtedzzJ1 = (-AdmInt(:,camadR)*TEdwz(:,camadR)*(exp(-u(:,camadR)*(z-h0)) + &
            RTEup(:,camadR)*FEupz(:)*exp(-u(:,camadR)*(z-prof(camadR-1))) - &
            RTEdw(:,camadR)*FEdwz(:)*exp(u(:,camadR)*(z-prof(camadR))))) * wJ1

         kernelEx = y * sum(KtezJ1 * kr * kr) / twopir
         Ex_p = kernelEx / r

         kernelEy = - x * sum(KtezJ1 * kr * kr) / twopir
         Ey_p = kernelEy / r

         kernelHx = - x * sum(KtedzzJ1 * kr * kr) / twopir
         Hx_p = kernelHx / r

         kernelHy = - y * sum(KtedzzJ1 * kr * kr) / twopir
         Hy_p = kernelHy / r

         kernelHz = sum(KtezJ0 * kr * kr * kr) / 2.0 / pi / zeta
         Hz_p = kernelHz / r
      elseif (camadR > camadT .and. camadR /= n) then !camada j
         fac = TEdwz(:,camadR)*(exp(-u(:,camadR)*(z-prof(camadR-1))) + &
            RTEdw(:,camadR)*exp(u(:,camadR)*(z-prof(camadR)-h(camadR))))
         KtezJ0 = fac * wJ0
         KtezJ1 = fac * wJ1
         KtedzzJ1 = (-AdmInt(:,camadR)*TEdwz(:,camadR)*(exp(-u(:,camadR)*(z-prof(camadR-1))) - &
            RTEdw(:,camadR)*exp(u(:,camadR)*(z-prof(camadR)-h(camadR))))) * wJ1

         kernelEx = y * sum(KtezJ1 * kr * kr) / twopir
         Ex_p = kernelEx / r

         kernelEy = - x * sum(KtezJ1 * kr * kr) / twopir
         Ey_p = kernelEy / r

         kernelHx = - x * sum(KtedzzJ1 * kr * kr) / twopir
         Hx_p = kernelHx / r

         kernelHy = - y * sum(KtedzzJ1 * kr * kr) / twopir
         Hy_p = kernelHy / r

         kernelHz = sum(KtezJ0 * kr * kr * kr) / 2.0 / pi / zeta
         Hz_p = kernelHz / r
      else !camada n
         fac = TEdwz(:,n)*exp(-u(:,n)*(z-prof(n-1)))
         KtezJ0 = fac * wJ0
         KtezJ1 = fac * wJ1
         KtedzzJ1 = -AdmInt(:,n) * fac * wJ1

         kernelEx = y * sum(KtezJ1 * kr * kr) / twopir
         Ex_p = kernelEx / r

         kernelEy = - x * sum(KtezJ1 * kr * kr) / twopir
         Ey_p = kernelEy / r

         kernelHx = - x * sum(KtedzzJ1 * kr * kr) / twopir
         Hx_p = kernelHx / r

         kernelHy = - y * sum(KtedzzJ1 * kr * kr) / twopir
         Hy_p = kernelHy / r

         kernelHz = sum(KtezJ0 * kr * kr * kr) / 2.0 / pi / zeta
         Hz_p = kernelHz / r
      end if
   end subroutine vmd_optimized

end module vmd
