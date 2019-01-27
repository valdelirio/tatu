module vmd
use parameters
use escolha_do_filtro
    contains
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine vmd_xyz_loops( Tx, Ty, h0, n, esp, condut, neta, zeta, cx, cy, z, Ex_p, Ey_p, Hx_p, Hy_p, Hz_p )
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: Tx, Ty, h0, esp(:), condut(1:n), cx, cy, z
  complex(dp), intent(in) :: zeta, neta
  complex(dp), intent(out) :: Ex_p, Ey_p, Hx_p, Hy_p, Hz_p

  integer :: i, j, k, camad, camadT, filtro, idtfcd_cJ0, ident_fJ0, nJ0, idtfcd_cJ1, ident_fJ1, nJ1
  real(dp) :: x, y, r
  real(dp), dimension(:), allocatable :: h, krJ0, krJ1, w_J0, w_J1, prof

! Para uso de loops:
  complex(dp), dimension(:), allocatable :: wvnb2, FEdwJ0, FEdwJ1, FEupJ0, FEupJ1
  complex(dp), dimension(:,:), allocatable :: uJ0, uJ1, AdmIntJ0, AdmIntJ1, uhJ0, uhJ1, tghJ0, tghJ1
  complex(dp), dimension(:,:), allocatable :: AdmApdwJ0, AdmApdwJ1, AdmApupJ0, AdmApupJ1
  complex(dp), dimension(:,:), allocatable :: RTEdwJ0, RTEdwJ1, RTEupJ0, RTEupJ1
  complex(dp), dimension(:,:), allocatable :: TEdwJ0, TEdwJ1, TEupJ0, TEupJ1

  complex(dp), dimension(:), allocatable :: Kte_J0, Kte_J1, Ktedz_J1
  complex(dp), dimension(:), allocatable :: kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1
  complex(dp), dimension(:), allocatable :: kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1
  real(dp), parameter :: eps = 1.d-7

  if ( dabs(cx - Tx) < eps ) then
    x = dsign( 1.d-1, Tx )
  else
    x = cx - Tx
  end if

  y = cy - Ty
  r = dsqrt( x ** 2 + y ** 2 )

  allocate( h(0 : n), prof(-1 : n) )
  if ( size(esp) == n ) then
    h(0) = 0.d0
    h(1 : n) = esp
  else
    h(0) = 0.d0
    h(1 : n - 1) = esp
    h(n) = 1.d300
  end if
! criando um novo vetor de profundidades que se adeque à qualquer situação patológica
  prof(-1) = -1.d300
  prof(0) = 0.d0
  if ( n > 1 ) then
        prof(1) = h(1)
        if ( n > 2 ) then
            do k = 2, n - 1
                prof(k) = prof(k - 1) + h(k)
            end do
        end if
  end if
  prof(n) = 1.d300
camad = 0
!para descobrir em que camada está a observação
  if ( z <= 0.d0 ) then
    camad = 0
  else if ( z > prof(n - 1) ) then
    camad = n
  else
        do i = n - 1, 1, -1
            if ( z > prof(i - 1) ) then
                camad = i
                exit
            end if
        end do
  end if
camadT = 0
!para descobrir em que camada está o transmissor
  if ( h0 <= 0.d0 ) then
        camadT = 0
  else if ( h0 > prof(n - 1) ) then
    camadT = n
  else
        do j = n - 1, 1, -1
            if ( h0 > prof(j - 1) ) then
                camadT = j
                exit
            end if
        end do
  end if

!!  write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!  read(*,*)idtfcd_cJ0
  filtro = 0  !esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!  call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
!!  write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!  read(*,*)idtfcd_cJ1
!!  call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)

  idtfcd_cJ0 = 5 !3
  ident_fJ0 = 0
  nJ0 = 201  !241
  idtfcd_cJ1 = 5 !3
  ident_fJ1 = 1
  nJ1 = 201  !241

  allocate( KrJ0(nJ0), KrJ1(nJ1), w_J0(nJ0), w_J1(nJ1) )

  call constfiltro( filtro, idtfcd_cJ0, ident_fJ0, nJ0, r, KrJ0, w_J0 )
  call constfiltro( filtro, idtfcd_cJ1, ident_fJ1, nJ1, r, KrJ1, w_J1 )

  allocate( wvnb2(0 : n), uJ0(nJ0, 0 : n), uJ1(nJ1, 0 : n), AdmIntJ0(nJ0, 0 : n), AdmIntJ1(nJ1, 0 : n) )
  allocate( uhJ0(nJ0, 0 : n), uhJ1(nJ1, 0 : n), tghJ0(nJ0, 0 : n), tghJ1(nJ1, 0 : n) )
! work around the warning: ... may be used uninitialized in this function
  allocate( TEupJ0(1,1), TEupJ1(1,1) )
  allocate( TEdwJ0(1,1), TEdwJ1(1,1) )
!
  do i = 0, n
        if ( i == 0 ) then
            wvnb2(i) = -zeta * neta
            uJ0(:,i) = sqrt( krJ0 * krJ0 - wvnb2(i) )
            uJ1(:,i) = sqrt( krJ1 * krJ1 - wvnb2(i) )
            AdmIntJ0(:,i) = uJ0(:,i) / zeta
            AdmIntJ1(:,i) = uJ1(:,i) / zeta
            uhJ0(:,i) = uJ0(:,i) * h(i)
            uhJ1(:,i) = uJ1(:,i) * h(i)
            tghJ0(:,i) = ( 1.d0 - exp( -2.d0 * uhJ0(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhJ0(:,i) ) )
            tghJ1(:,i) = ( 1.d0 - exp( -2.d0 * uhJ1(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhJ1(:,i) ) )
    else
            wvnb2(i) = -zeta * condut(i)
            uJ0(:,i) = sqrt( krJ0 * krJ0 - wvnb2(i) )
            uJ1(:,i) = sqrt( krJ1 * krJ1 - wvnb2(i) )
            AdmIntJ0(:,i) = uJ0(:,i) / zeta
            AdmIntJ1(:,i) = uJ1(:,i) / zeta
            uhJ0(:,i) = uJ0(:,i) * h(i)
            uhJ1(:,i) = uJ1(:,i) * h(i)
            tghJ0(:,i) = ( 1.d0 - exp( -2.d0 * uhJ0(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhJ0(:,i) ) )
            tghJ1(:,i) = ( 1.d0 - exp( -2.d0 * uhJ1(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhJ1(:,i) ) )
    end if
  end do

  allocate( AdmApdwJ0(nJ0, 1 : n), AdmApdwJ1(nJ1, 1 : n) )
  allocate( RTEdwJ0(nJ0, 0 : n), RTEdwJ1(nJ1, 0 : n) )

  do i = n, 1, -1
        if ( i == n ) then
            AdmApdwJ0(:,i) = AdmIntJ0(:,i)
            AdmApdwJ1(:,i) = AdmIntJ1(:,i)
            RTEdwJ0(:,i) = (0.d0, 0.d0)
            RTEdwJ1(:,i) = (0.d0, 0.d0)
    else
            AdmApdwJ0(:,i) = AdmIntJ0(:,i) * ( AdmApdwJ0(:,i + 1) + AdmIntJ0(:,i) * &
              tghJ0(:,i) ) / ( AdmIntJ0(:,i) + AdmApdwJ0(:,i + 1) * tghJ0(:,i) )
            AdmApdwJ1(:,i) = AdmIntJ1(:,i) * ( AdmApdwJ1(:,i + 1) + AdmIntJ1(:,i) * &
              tghJ1(:,i) ) / ( AdmIntJ1(:,i) + AdmApdwJ1(:,i + 1) * tghJ1(:,i) )
            RTEdwJ0(:,i) = ( AdmIntJ0(:,i) - AdmApdwJ0(:,i + 1) ) / ( AdmIntJ0(:,i) + AdmApdwJ0(:,i + 1) )
            RTEdwJ1(:,i) = ( AdmIntJ1(:,i) - AdmApdwJ1(:,i + 1) ) / ( AdmIntJ1(:,i) + AdmApdwJ1(:,i + 1) )
    end if
  end do
    RTEdwJ0(:,0) = ( AdmIntJ0(:,0) - AdmApdwJ0(:,1) ) / ( AdmIntJ0(:,0) + AdmApdwJ0(:,1) )
    RTEdwJ1(:,0) = ( AdmIntJ1(:,0) - AdmApdwJ1(:,1) ) / ( AdmIntJ1(:,0) + AdmApdwJ1(:,1) )

  allocate( AdmApupJ0(nJ0,0 : n - 1), AdmApupJ1(nJ1, 0 : n - 1) )
  allocate( RTEupJ0(nJ0, 0 : n), RTEupJ1(nJ1,0 : n) )

  do i = 0, n - 1
        if ( i == 0 ) then
            AdmApupJ0(:,i) = AdmIntJ0(:,i)
            AdmApupJ1(:,i) = AdmIntJ1(:,i)
            RTEupJ0(:,i) = (0.d0,0.d0)
            RTEupJ1(:,i) = (0.d0,0.d0)
    else
            AdmApupJ0(:,i) = AdmIntJ0(:,i) * ( AdmApupJ0(:,i - 1) + AdmIntJ0(:,i) * &
              tghJ0(:,i) ) / ( AdmIntJ0(:,i) + AdmApupJ0(:,i - 1) * tghJ0(:,i) )
            AdmApupJ1(:,i) = AdmIntJ1(:,i) * ( AdmApupJ1(:,i - 1) + AdmIntJ1(:,i) * &
              tghJ1(:,i) ) / ( AdmIntJ1(:,i) + AdmApupJ1(:,i - 1) * tghJ1(:,i) )
            RTEupJ0(:,i) = ( AdmIntJ0(:,i) - AdmApupJ0(:,i - 1) ) / ( AdmIntJ0(:,i) + AdmApupJ0(:,i - 1) )
            RTEupJ1(:,i) = ( AdmIntJ1(:,i) - AdmApupJ1(:,i - 1) ) / ( AdmIntJ1(:,i) + AdmApupJ1(:,i - 1) )
    end if
  end do
    RTEupJ0(:,n) = ( AdmIntJ0(:,n) - AdmApupJ0(:,n - 1) ) / ( AdmIntJ0(:,n) + AdmApupJ0(:,n - 1) )
    RTEupJ1(:,n) = ( AdmIntJ1(:,n) - AdmApupJ1(:,n - 1) ) / ( AdmIntJ1(:,n) + AdmApupJ1(:,n - 1) )

  allocate( FEdwJ0(nJ0), FEdwJ1(nJ1), FEupJ0(nJ0), FEupJ1(nJ1) )

  FEdwJ0 = ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0) ) + RTEupJ0(:,camadT) * &
            exp( uJ0(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ0(:,camadT) * RTEdwJ0(:,camadT) * exp( -2.d0 * uhJ0(:,camadT) ) )
  FEdwJ1 = ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0) ) + RTEupJ1(:,camadT) * &
            exp( uJ1(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ1(:,camadT) * RTEdwJ1(:,camadT) * exp( -2.d0 * uhJ1(:,camadT) ) )

  FEupJ0 = ( exp( uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTEdwJ0(:,camadT) * &
            exp( -uJ0(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ0(:,camadT) * RTEdwJ0(:,camadT) * exp( -2.d0 * uhJ0(:,camadT) ) )
  FEupJ1 = ( exp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTEdwJ1(:,camadT) * &
            exp( -uJ1(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ1(:,camadT) * RTEdwJ1(:,camadT) * exp( -2.d0 * uhJ1(:,camadT) ) )

  if ( camad > camadT ) then
        deallocate(TEdwJ0, TEdwJ1)
        allocate( TEdwJ0(nJ0,camadT : camad), TEdwJ1(nJ1,camadT : camad) )
        do j = camadT, camad
            if ( j == camadT ) then

                TEdwJ0(:,j)= zeta * mz / ( 2.d0 * uJ0(:,j) )
                TEdwJ1(:,j)= zeta * mz / ( 2.d0 * uJ1(:,j) )

      else if ( j == (camadT + 1) .and. j == n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + &
                        RTEupJ0(:,camadT) * FEupJ0(:) * exp( -uhJ0(:,camadT) ) + RTEdwJ0(:,camadT) * FEdwJ0(:) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0) ) + &
                        RTEupJ1(:,camadT) * FEupJ1(:) * exp( -uhJ1(:,camadT) ) + RTEdwJ1(:,camadT) * FEdwJ1(:) )

      else if ( j == (camadT + 1) .and. j /= n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + &
                        RTEupJ0(:,camadT) * FEupJ0(:) * exp( -uhJ0(:,camadT) ) + &
                        RTEdwJ0(:,camadT) * FEdwJ0(:) ) / ( 1.d0 + RTEdwJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) + &
                        RTEupJ1(:,camadT) * FEupJ1(:) * exp( -uhJ1(:,camadT) ) + &
                        RTEdwJ1(:,camadT) * FEdwJ1(:) ) / ( 1.d0 + RTEdwJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )

      else if ( j /= n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( 1.d0 + RTEdwJ0(:,j - 1) ) * exp( -uhJ0(:,j - 1) ) / &
                        ( 1.d0 + RTEdwJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( 1.d0 + RTEdwJ1(:,j - 1) ) * exp( -uhJ1(:,j - 1) ) / &
                        ( 1.d0 + RTEdwJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )

      else if ( j == n ) then

                TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( 1.d0 + RTEdwJ0(:,j - 1) ) * exp( -uhJ0(:,j - 1) )
                TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( 1.d0 + RTEdwJ1(:,j - 1) ) * exp( -uhJ1(:,j - 1) )

      end if
    end do
  else if ( camad < camadT ) then
        deallocate(TEupJ0, TEupJ1)
        allocate( TEupJ0(nJ0,camad : camadT), TEupJ1(nJ1,camad : camadT) )
        do j = camadT, camad, -1
            if ( j == camadT ) then

                TEupJ0(:,j) = zeta * mz / ( 2.d0 * uJ0(:,j) )
                TEupJ1(:,j) = zeta * mz / ( 2.d0 * uJ1(:,j) )

      else if ( j == (camadT - 1) .and. j == 0 ) then

                TEupJ0(:,j) = TEupJ0(:,j + 1) * ( exp( -uJ0(:,camadT) * h0 ) + &
                    RTEupJ0(:,camadT) * FEupJ0(:) + RTEdwJ0(:,camadT) * FEdwJ0(:) * exp( -uhJ0(:,camadT) ) )
                TEupJ1(:,j) = TEupJ1(:,j + 1) * ( exp( -uJ1(:,camadT) * h0 ) + &
                    RTEupJ1(:,camadT) * FEupJ1(:) + RTEdwJ1(:,camadT) * FEdwJ1(:) * exp( -uhJ1(:,camadT) ) )

      else if ( j == (camadT - 1) .and. j /= 0 ) then

                TEupJ0(:,j) = TEupJ0(:,j + 1) * ( exp(uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                    RTEupJ0(:,camadT) * FEupJ0(:) + RTEdwJ0(:,camadT) * FEdwJ0(:) * &
                    exp( -uhJ0(:,camadT) ) ) / ( 1.d0 + RTEupJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
                TEupJ1(:,j) = TEupJ1(:,j + 1) * ( exp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                    RTEupJ1(:,camadT) * FEupJ1(:) + RTEdwJ1(:,camadT) * FEdwJ1(:) * &
                    exp( -uhJ1(:,camadT) ) ) / ( 1.d0 + RTEupJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )

      else if ( j /= 0 ) then

                TEupJ0(:,j) = TEupJ0(:,j + 1) * ( 1.d0 + RTEupJ0(:,j + 1) ) * exp( -uhJ0(:,j + 1) ) / &
                    ( 1.d0 + RTEupJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
                TEupJ1(:,j) = TEupJ1(:,j + 1) * ( 1.d0 + RTEupJ1(:,j + 1) ) * exp( -uhJ1(:,j + 1) ) / &
                    ( 1.d0 + RTEupJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )

      else if ( j == 0 ) then

                TEupJ0(:,j) = TEupJ0(:,1) * ( 1.d0 + RTEupJ0(:,1) ) * exp( -uhJ0(:,1) )
                TEupJ1(:,j) = TEupJ1(:,1) * ( 1.d0 + RTEupJ1(:,1) ) * exp( -uhJ1(:,1) )

      end if
    end do
  else
    deallocate(TEdwJ0, TEdwJ1, TEupJ0, TEupJ1)
    allocate( TEdwJ0(nJ0,camadT : camad), TEdwJ1(nJ1,camadT : camad) )
    allocate( TEupJ0(nJ0,camad : camadT), TEupJ1(nJ1,camad : camadT) )

    TEdwJ0(:,camad) = zeta * mz / ( 2.d0 * uJ0(:,camadT) )
    TEdwJ1(:,camad) = zeta * mz / ( 2.d0 * uJ1(:,camadT) )

    TEupJ0(:,camad) = TEdwJ0(:,camad)
    TEupJ1(:,camad) = TEdwJ1(:,camad)

  end if

  allocate( Kte_J0(nJ0), Kte_J1(nJ1), Ktedz_J1(nJ1) )
  allocate( kernelExJ0(nJ0), kernelExJ1(nJ1), kernelEyJ0(nJ0), kernelEyJ1(nJ1) )
  allocate( kernelHxJ0(nJ0), kernelHxJ1(nJ1), kernelHyJ0(nJ0), kernelHyJ1(nJ1), kernelHzJ1(nJ1) )
  if ( camad == 0 .and. camadT /= 0 ) then

    Kte_J0 = TEupJ0(:,0) * exp( uJ0(:,0) * z ) * w_J0(:)
    Kte_J1 = TEupJ1(:,0) * exp( uJ1(:,0) * z ) * w_J1(:)
    Ktedz_J1 = AdmIntJ1(:,0) * Kte_J1

    kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
    Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

    else if ( camad < camadT ) then !camada k

    Kte_J0 = ( TEupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - prof(camad) ) ) + &
                RTEupJ0(:,camad) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - prof(camad) ) ) + &
                RTEupJ1(:,camad) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)
    Ktedz_J1 = ( AdmIntJ1(:,camad) * TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - prof(camad) ) ) - &
                RTEupJ1(:,camad) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelEyJ1 = - x * Kte_J1 * krJ1 ** 2 / ( 2.d0 * pi * r )
    Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
    Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

  else if ( camad == camadT .and. z <= h0 ) then  !na mesma camada do transmissor mas acima dele

    Kte_J0 = ( TEupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - h0 ) ) + &
                RTEupJ0(:,camad) * FEupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * FEdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEupJ1(:,camad) * ( exp( uJ1 (:,camad) * ( z - h0 ) ) + &
                RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
    Ktedz_J1 = ( AdmIntJ1(:,camad) * TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - h0 ) ) - &
                RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
    Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

  else if ( camad == camadT .and. z > h0 ) then !na mesma camada do transmissor mas abaixo dele

    Kte_J0 = ( TEdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - h0 ) ) + &
                RTEupJ0(:,camad) * FEupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * FEdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - h0 ) ) + &
                RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
    Ktedz_J1 = ( - AdmIntJ1(:,camad) * TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - h0 ) ) + &
                RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
    Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

  else if ( camad > camadT .and. camad /= n ) then !camada j

    Kte_J0 = ( TEdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * exp( uJ0(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * exp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)
    Ktedz_J1 = ( - AdmIntJ1(:,camad) * TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwJ1(:,camad) * exp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
    Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

  else  !camada n

    Kte_J0 = ( TEdwJ0(:,n) * exp( -uJ0(:,n) * ( z - prof(n - 1) ) ) ) * w_J0(:)
    Kte_J1 = ( TEdwJ1(:,n) * exp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)
    Ktedz_J1=( - AdmIntJ1(:,n) * TEdwJ1(:,n) * exp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)

    kernelExJ1 = y * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ex_p = sum( kernelExJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelEyJ1 = - x * Kte_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Ey_p = sum( kernelEyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHxJ1 = - x * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hx_p = sum( kernelHxJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHyJ1 = - y * Ktedz_J1 * krJ1 * krJ1 / ( 2.d0 * pi * r )
    Hy_p = sum( kernelHyJ1 ) / r     !este último r é decorrente do uso dos filtros

    kernelHzJ1 = Kte_J0 * krJ0 * krJ0 * krJ0 / ( 2.d0 * pi * zeta )
    Hz_p = sum( kernelHzJ1 ) / r     !este último r é decorrente do uso dos filtros

  end if

  deallocate( h, KrJ0, KrJ1, w_J0, w_J1 )
  deallocate( wvnb2, uJ0, uJ1, AdmIntJ0, AdmIntJ1, uhJ0, uhJ1, tghJ0, tghJ1 )
  deallocate( AdmApdwJ0, AdmApdwJ1, RTEdwJ0, RTEdwJ1 )
  deallocate( AdmApupJ0, AdmApupJ1, RTEupJ0, RTEupJ1 )
  deallocate( FEdwJ0, FEdwJ1, FEupJ0, FEupJ1 )
  deallocate( Kte_J0, Kte_J1, Ktedz_J1 )
  deallocate( kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1 )
  deallocate( kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1 )
end subroutine vmd_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine vmd_xkyz_loops( Tx, ky, h0, n, esp, condut, neta, zeta, cx, z, Ex_ky, Ey_ky, Hx_ky, Hy_ky, Hz_ky )
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: Tx, ky, h0, esp(:), condut(1 : n), cx, z
  complex(dp), intent(in) :: zeta, neta
  complex(dp), intent(out) :: Ex_ky, Ey_ky, Hx_ky, Hy_ky, Hz_ky

  integer :: i, j, k, camad, camadT, autor, filtro, npts, nptc, funs, func
  real(dp) :: x
  real(dp), dimension(:), allocatable :: h, kxsen, kxcos, kr2sen, kr2cos, w_sen, w_cos, prof

  complex(dp), dimension(:), allocatable :: wvnb2, FEdwSen, FEdwCos, FEupSen, FEupCos
  complex(dp), dimension(:,:), allocatable :: uSen, uCos, AdmIntSen, AdmIntCos
  complex(dp), dimension(:,:), allocatable :: uhSen, uhCos, tghSen, tghCos
  complex(dp), dimension(:,:), allocatable :: AdmApdwSen, AdmApdwCos
  complex(dp), dimension(:,:), allocatable :: RTEdwSen, RTEdwCos
  complex(dp), dimension(:,:), allocatable :: AdmApupSen, AdmApupCos
  complex(dp), dimension(:,:), allocatable :: RTEupSen, RTEupCos
  complex(dp), dimension(:,:), allocatable :: TEdwSen, TEdwCos, TEupSen, TEupCos
  complex(dp), dimension(:), allocatable :: Kte_Sen, Kte_Cos, Ktedz_Sen, Ktedz_Cos
  complex(dp), dimension(:), allocatable :: kernelEx, kernelEy, kernelHx, kernelHy, kernelHz
  real(dp), parameter :: eps = 1.d-7

  if ( dabs(cx - Tx) < eps ) then
    x = dsign( 1.d-1, Tx )
  else
    x = cx - Tx
  end if

  allocate( h(0 : n), prof(-1 : n) )
  if ( size(esp) == n ) then
    h(0) = 0.d0
    h(1 : n) = esp
  else
    h(0) = 0.d0
    h(1 : n - 1) = esp
    h(n) = 1.d300
  end if
! criando um novo vetor de profundidades que se adeque à qualquer situação patológica
  prof(-1) = -1.d300
  prof(0) = 0.d0
  if ( n > 1 ) then
        prof(1) = h(1)
        if ( n > 2 ) then
            do k = 2, n - 1
                prof(k) = prof(k-1) + h(k)
            end do
        end if
  end if
  prof(n) = 1.d300
  camad = 0
!para descobrir em que camada está a observação
  if ( z < 0.d0 ) then
        camad = 0
  else if ( z >= prof(n - 1) ) then
    camad = n
  else
        do i = n - 1, 1, -1
            if ( z >= prof(i - 1) ) then
                camad = i
                exit
            end if
        end do
  end if
  camadT = 0
!para descobrir em que camada está o transmissor
  if ( h0 < 0.d0 ) then
        camadT = 0
  else if ( h0 >= prof(n - 1) ) then
    camadT = n
  else
        do j = n - 1, 1, -1
            if ( h0 >= prof(j - 1) ) then
                camadT = j
                exit
            end if
        end do
  end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! work around the warning: ... may be used uninitialized in this function
allocate( TEupSen(1,1), TEupCos(1,1) )
allocate( TEdwSen(1,1), TEdwCos(1,1) )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  filtro = 1    !Designa o tipo de filtro usado na subrotina de pesos e abscisas de vários filtros.
          !O algarismo 0 é usado para J0 e J1, enquanto 1 é para seno e cosseno.
  autor = 4   !Designa o criador do filtro. No caso de seno ou cosseno os que possuo são os do Frayzer (1) e do Kerry Key (4)
  funs = 2    !Designa se é filtro seno (2) ou filtro cosseno (3)
  func = 3    !Designa se é filtro seno (2) ou filtro cosseno (3)
  npts = 81   !Designa o número de pontos usado no filtro seno.
  nptc = 81   !Designa o número de pontos usado no filtro cosseno.

  allocate( kxsen(npts), kxcos(nptc), w_sen(npts), w_cos(nptc) )

  call constfiltro( filtro, autor, funs, npts, x, Kxsen, w_sen )
  call constfiltro( filtro, autor, func, nptc, x, Kxcos, w_cos )

  allocate( wvnb2(0 : n), uSen(npts,0 : n), uCos(nptc,0 : n), AdmIntSen(npts,0 : n), AdmIntCos(nptc,0 : n) )
  allocate( uhSen(npts,0 : n), uhCos(nptc,0 : n), tghSen(npts,0 : n), tghCos(nptc,0 : n), kr2sen(npts), kr2cos(nptc) )

  kr2sen = kxsen * kxsen + ky * ky
  kr2cos = kxcos * kxcos + ky * ky
  do i = 0, n
    if ( i == 0 ) then
            wvnb2(i) = -zeta * neta
            uSen(:,i) = sqrt( kr2sen - wvnb2(i) )
            uCos(:,i) = sqrt( kr2cos - wvnb2(i) )
            AdmIntSen(:,i) = uSen(:,i) / zeta
            AdmIntCos(:,i) = uCos(:,i) / zeta
            uhSen(:,i) = uSen(:,i) * h(i)
            uhCos(:,i) = uCos(:,i) * h(i)
            tghSen(:,i) = ( 1.d0 - exp( -2.d0 * uhSen(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhSen(:,i) ) )
            tghCos(:,i) = ( 1.d0 - exp( -2.d0 * uhCos(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhCos(:,i) ) )
    else
            wvnb2(i) = -zeta * condut(i)
            uSen(:,i) = sqrt( kr2Sen - wvnb2(i) )
            uCos(:,i) = sqrt( kr2Cos - wvnb2(i) )
            AdmIntSen(:,i) = uSen(:,i) / zeta
            AdmIntCos(:,i) = uCos(:,i) / zeta
            uhSen(:,i) = uSen(:,i) * h(i)
            uhCos(:,i) = uCos(:,i) * h(i)
            tghSen(:,i) = ( 1.d0 - exp( -2.d0 * uhSen(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhSen(:,i) ) )
            tghCos(:,i) = ( 1.d0 - exp( -2.d0 * uhCos(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhCos(:,i) ) )
    end if
  end do

  allocate( AdmApdwSen(npts,1 : n), AdmApdwCos(nptc,1 : n) )
  allocate( RTEdwSen(npts,0 : n), RTEdwCos(nptc,0 : n) )

  do i = n, 1, -1
    if ( i == n ) then
            AdmApdwSen(:,i) = AdmIntSen(:,i)
            AdmApdwCos(:,i) = AdmIntCos(:,i)
            RTEdwSen(:,i) = (0.d0,0.d0)
            RTEdwCos(:,i) = (0.d0,0.d0)
    else
            AdmApdwSen(:,i) = AdmIntSen(:,i) * ( AdmApdwSen(:,i + 1) + AdmIntSen(:,i) * &
                tghSen(:,i) ) / ( AdmIntSen(:,i) + AdmApdwSen(:,i + 1) * tghSen(:,i) )
            AdmApdwCos(:,i) = AdmIntCos(:,i) * ( AdmApdwCos(:,i + 1) + AdmIntCos(:,i) * &
              tghCos(:,i) ) / ( AdmIntCos(:,i) + AdmApdwCos(:,i + 1) * tghCos(:,i) )
            RTEdwSen(:,i) = ( AdmIntSen(:,i) - AdmApdwSen(:,i + 1) ) / ( AdmIntSen(:,i) + AdmApdwSen(:,i + 1) )
            RTEdwCos(:,i) = ( AdmIntCos(:,i) - AdmApdwCos(:,i + 1) ) / ( AdmIntCos(:,i) + AdmApdwCos(:,i + 1) )
    end if
  end do
  RTEdwSen(:,0) = ( AdmIntSen(:,0) - AdmApdwSen(:,1) ) / ( AdmIntSen(:,0) + AdmApdwSen(:,1) )
    RTEdwCos(:,0) = ( AdmIntCos(:,0) - AdmApdwCos(:,1) ) / ( AdmIntCos(:,0) + AdmApdwCos(:,1) )

  allocate( AdmApupSen(npts,0 : n - 1), AdmApupCos(nptc,0 : n - 1) )
  allocate( RTEupSen(npts,0 : n), RTEupCos(nptc,0 : n) )

  do i = 0, n - 1
    if ( i == 0 ) then
            AdmApupSen(:,i) = AdmIntSen(:,i)
            AdmApupCos(:,i) = AdmIntCos(:,i)
            RTEupSen(:,i) = (0.d0,0.d0)
            RTEupCos(:,i) = (0.d0,0.d0)
    else
            AdmApupSen(:,i) = AdmIntSen(:,i) * ( AdmApupSen(:,i - 1) + AdmIntSen(:,i) * &
                tghSen(:,i) ) / ( AdmIntSen(:,i) + AdmApupSen(:,i - 1) * tghSen(:,i) )
            AdmApupCos(:,i) = AdmIntCos(:,i) * ( AdmApupCos(:,i - 1) + AdmIntCos(:,i) * &
              tghCos(:,i) ) / ( AdmIntCos(:,i) + AdmApupCos(:,i - 1) * tghCos(:,i) )
            RTEupSen(:,i) = ( AdmIntSen(:,i) - AdmApupSen(:,i - 1) ) / ( AdmIntSen(:,i) + AdmApupSen(:,i - 1) )
            RTEupCos(:,i) = ( AdmIntCos(:,i) - AdmApupCos(:,i - 1) ) / ( AdmIntCos(:,i) + AdmApupCos(:,i - 1) )
    end if
  end do
  RTEupSen(:,n) = ( AdmIntSen(:,n) - AdmApupSen(:,n - 1) ) / ( AdmIntSen(:,n) + AdmApupSen(:,n - 1) )
    RTEupCos(:,n) = ( AdmIntCos(:,n) - AdmApupCos(:,n - 1) ) / ( AdmIntCos(:,n) + AdmApupCos(:,n - 1) )

    allocate( FEdwSen(npts), FEdwCos(nptc), FEupSen(npts), FEupCos(nptc) )

  FEdwSen = ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) + RTEupSen(:,camadT) * &
                exp( uSen(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
        ( 1.d0 - RTEupSen(:,camadT) * RTEdwSen(:,camadT) * exp( -2.d0 * uhSen(:,camadT) ) )
  FEdwCos = ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) + RTEupCos(:,camadT) * &
        exp( uCos(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
        ( 1.d0 - RTEupCos(:,camadT) * RTEdwCos(:,camadT) * exp( -2.d0 * uhCos(:,camadT) ) )

  FEupSen = ( exp( uSen(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTEdwSen(:,camadT) * &
        exp( -uSen(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
        ( 1.d0 - RTEupSen(:,camadT) * RTEdwSen(:,camadT) * exp( -2.d0 * uhSen(:,camadT) ) )
  FEupCos = ( exp( uCos(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTEdwCos(:,camadT) * &
        exp( -uCos(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
        ( 1.d0 - RTEupCos(:,camadT) * RTEdwCos(:,camadT) * exp( -2.d0 * uhCos(:,camadT) ) )

  if ( camad > camadT ) then
        deallocate(TEdwSen,TEdwCos)
        allocate( TEdwSen(npts,camadT : camad), TEdwCos(nptc,camadT : camad) )
        do j = camadT, camad
            if ( j == camadT ) then

                TEdwSen(:,j) = zeta * mz / ( 2.d0 * uSen(:,camadT) )
                TEdwCos(:,j) = zeta * mz / ( 2.d0 * uCos(:,camadT) )

      else if ( j == (camadT + 1) .and. j == n ) then

                TEdwSen(:,j) = TEdwSen(:,j - 1) * ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) + &
                    RTEupSen(:,camadT) * FEupSen(:) * exp( -uhSen(:,camadT) ) + RTEdwSen(:,camadT) * FEdwSen(:) )
                TEdwCos(:,j) = TEdwCos(:,j - 1) * ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) + &
                    RTEupCos(:,camadT) * FEupCos(:) * exp( -uhCos(:,camadT) ) + RTEdwCos(:,camadT) * FEdwCos(:) )

      else if ( j == (camadT + 1) .and. j /= n ) then

                TEdwSen(:,j) = TEdwSen(:,j - 1) * ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) + &
                    RTEupSen(:,camadT) * FEupSen(:) * exp( -uhSen(:,camadT) ) + &
                    RTEdwSen(:,camadT) * FEdwSen(:) ) / ( 1.d0 + RTEdwSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
                TEdwCos(:,j) = TEdwCos(:,j - 1) * ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) + &
                    RTEupCos(:,camadT) * FEupCos(:) * exp( -uhCos(:,camadT) ) + &
                    RTEdwCos(:,camadT) * FEdwCos(:) ) / ( 1.d0 + RTEdwCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )

      else if ( j /= n ) then

            TEdwSen(:,j) = TEdwSen(:,j - 1) * ( 1.d0 + RTEdwSen(:,j - 1) ) * exp( -uhSen(:,j - 1) )/ &
                    ( 1.d0 + RTEdwSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
            TEdwCos(:,j) = TEdwCos(:,j - 1) * ( 1.d0 + RTEdwCos(:,j - 1) ) * exp( -uhCos(:,j - 1) ) / &
                    ( 1.d0 + RTEdwCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )

      else if ( j == n ) then

          TEdwSen(:,j) = TEdwSen(:,j - 1) * ( 1.d0 + RTEdwSen(:,j - 1) ) * exp( -uhSen(:,j - 1) )
          TEdwCos(:,j) = TEdwCos(:,j - 1) * ( 1.d0 + RTEdwCos(:,j - 1) ) * exp( -uhCos(:,j - 1) )

      end if
    end do
  else if ( camad < camadT ) then
        deallocate( TEupSen, TEupCos)
        allocate( TEupSen(npts,camad : camadT), TEupCos(nptc,camad : camadT) )
        do j = camadT, camad, -1
            if ( j == camadT ) then

                TEupSen(:,j) = zeta * mz / ( 2.d0 * uSen(:,camadT) )
                TEupCos(:,j) = zeta * mz / ( 2.d0 * uCos(:,camadT) )

      else if ( j == (camadT - 1) .and. j == 0 ) then

                TEupSen(:,j) = TEupSen(:,j + 1) * ( exp( -uSen(:,camadT) * h0 ) + &
                    RTEupSen(:,camadT) * FEupSen(:) + RTEdwSen(:,camadT) * FEdwSen(:) * exp( -uhSen(:,camadT) ) )
                TEupCos(:,j) = TEupCos(:,j + 1) * ( exp( -uCos(:,camadT) * h0 ) + &
                    RTEupCos(:,camadT) * FEupCos(:) + RTEdwCos(:,camadT) * FEdwCos(:) * exp( -uhCos(:,camadT) ) )

      else if ( j == (camadT - 1) .and. j /= 0 ) then

                TEupSen(:,j) = TEupSen(:,j + 1) * ( exp( uSen(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                    RTEupSen(:,camadT) * FEupSen(:) + RTEdwSen(:,camadT) * FEdwSen(:) * &
                    exp( -uhSen(:,camadT) ) ) / ( 1.d0 + RTEupSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
                TEupCos(:,j) = TEupCos(:,j + 1) * ( exp( uCos(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                    RTEupCos(:,camadT) * FEupCos(:) + RTEdwCos(:,camadT) * FEdwCos(:) * &
                    exp( -uhCos(:,camadT) ) ) / ( 1.d0 + RTEupCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )

      else if ( j /= 0 ) then

                TEupSen(:,j) = TEupSen(:,j + 1) * ( 1.d0 + RTEupSen(:,j + 1) ) * exp( -uhSen(:,j + 1) ) / &
                    ( 1.d0 + RTEupSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
                TEupCos(:,j) = TEupCos(:,j + 1) * ( 1.d0 + RTEupCos(:,j + 1) ) * exp( -uhCos(:,j + 1) ) / &
                    ( 1.d0 + RTEupCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )
      else if ( j == 0 ) then

                TEupSen(:,j) = TEupSen(:,1) * ( 1.d0 + RTEupSen(:,1) ) * exp( -uhSen(:,1) )
                TEupCos(:,j) = TEupCos(:,1) * ( 1.d0 + RTEupCos(:,1) ) * exp( -uhCos(:,1) )

      end if
    end do
  else
        deallocate( TEupSen, TEupCos, TEdwSen, TEdwCos)
        allocate( TEdwSen(npts,camadT : camad), TEdwCos(nptc,camadT : camad) )
        allocate( TEupSen(npts,camad : camadT), TEupCos(nptc,camad : camadT) )

    TEdwSen(:,camad) = zeta * mz / ( 2.d0 * uSen(:,camadT) )
    TEdwCos(:,camad) = zeta * mz / ( 2.d0 * uCos(:,camadT) )

    TEupSen(:,camad) = TEdwSen(:,camad)
    TEupCos(:,camad) = TEdwCos(:,camad)

  end if

  allocate( Kte_Sen(npts), Kte_Cos(nptc), Ktedz_Sen(npts), Ktedz_Cos(nptc) )
  allocate( kernelEx(nptc), kernelEy(npts) )
  allocate( kernelHx(npts), kernelHy(nptc), kernelHz(nptc) )

  if ( camad == 0 .and. camadT /= 0 ) then

    Kte_Sen = TEupSen(:,0) * exp( uSen(:,0) * z )
    Kte_Cos = TEupCos(:,0) * exp( uCos(:,0) * z )

    Ktedz_Sen = AdmIntSen(:,0) * Kte_Sen
    Ktedz_Cos = AdmIntCos(:,0) * Kte_Cos

    kernelEx = cmplx(0.d0,-ky,kind=dp) * Kte_Cos * w_cos
    Ex_ky = sum( kernelEx ) / ( pi * dabs(x) )

    kernelEy = - kxsen * Kte_Sen * w_sen
    Ey_ky = sum( kernelEy ) / ( pi * x )

    kernelHx = - kxsen *  Ktedz_Sen * w_sen
    Hx_ky = sum( kernelHx ) / ( pi * x )

    kernelHy = cmplx(0.d0,ky,kind=dp) * Ktedz_Cos * w_cos
    Hy_ky = sum( kernelHy ) / ( pi * dabs(x) )

    kernelHz = ( kxcos * kxcos + ky * ky ) * Kte_Cos * w_cos / zeta
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )

  else if ( camad < camadT ) then !layer k

    Kte_Sen = TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - prof(camad) ) ) + &
                RTEupSen(:,camad) * exp( -uSen(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    Kte_Cos = TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - prof(camad) ) ) + &
                RTEupCos(:,camad) * exp( -uCos(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )

    Ktedz_Sen = AdmIntSen(:,camad) * TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - prof(camad) ) ) - &
                RTEupSen(:,camad) * exp( -uSen(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    Ktedz_Cos = AdmIntCos(:,camad) * TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - prof(camad) ) ) - &
                RTEupCos(:,camad) * exp( -uCos(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )

    kernelEx = cmplx(0.d0,-ky,kind=dp) * Kte_Cos * w_cos
    Ex_ky = sum( kernelEx ) / ( pi * dabs(x) )

    kernelEy = - kxsen * Kte_Sen * w_sen
    Ey_ky = sum( kernelEy ) / ( pi * x )

    kernelHx = - kxsen *  Ktedz_Sen * w_sen
    Hx_ky = sum( kernelHx ) / ( pi * x )

    kernelHy = cmplx(0.d0,ky,kind=dp) * Ktedz_Cos * w_cos
    Hy_ky = sum( kernelHy ) / ( pi * dabs(x) )

    kernelHz = ( kxcos * kxcos + ky * ky ) * Kte_Cos * w_cos / zeta
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )

  else if ( camad == camadT .and. z <= h0 ) then  !in the same layer, but receiver above transmitter

    Kte_Sen = TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - h0 ) ) + &
                RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Kte_Cos = TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - h0 ) ) + &
                RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )
        Ktedz_Sen = AdmIntSen(:,camad) * TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - h0 ) ) - &
                RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktedz_Cos = AdmIntCos(:,camad) * TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - h0 ) ) - &
                RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    kernelEx = cmplx(0.d0,-ky,kind=dp) * Kte_Cos * w_cos
    Ex_ky = sum( kernelEx ) / ( pi * dabs(x) )

    kernelEy = - kxsen * Kte_Sen * w_sen
    Ey_ky = sum( kernelEy ) / ( pi * x )

    kernelHx = - kxsen *  Ktedz_Sen * w_sen
    Hx_ky = sum( kernelHx ) / ( pi * x )

    kernelHy = cmplx(0.d0,ky,kind=dp) * Ktedz_Cos * w_cos
    Hy_ky = sum( kernelHy ) / ( pi * dabs(x) )

    kernelHz = ( kxcos * kxcos + ky * ky ) * Kte_Cos * w_cos / zeta
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )

  else if ( camad == camadT .and. z > h0 ) then !in the same layer, but receiver below transmitter

    Kte_Sen = TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - h0 ) ) + &
      RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
      RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Kte_Cos = TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - h0 ) ) + &
      RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
      RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    Ktedz_Sen = - AdmIntSen(:,camad) * TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - h0 ) ) + &
      RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
      RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktedz_Cos = - AdmIntCos(:,camad) * TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - h0 ) ) + &
      RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
      RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    kernelEx = cmplx(0.d0,-ky,kind=dp) * Kte_Cos * w_cos
    Ex_ky = sum( kernelEx ) / ( pi * dabs(x) )

    kernelEy = - kxsen * Kte_Sen * w_sen
    Ey_ky = sum( kernelEy ) / ( pi * x )

    kernelHx = - kxsen *  Ktedz_Sen * w_sen
    Hx_ky = sum( kernelHx ) / ( pi * x )

    kernelHy = cmplx(0.d0,ky,kind=dp) * Ktedz_Cos * w_cos
    Hy_ky = sum( kernelHy ) / ( pi * dabs(x) )

    kernelHz = ( kxcos * kxcos + ky * ky ) * Kte_Cos * w_cos / zeta
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )

  else if ( camad > camadT .and. camad /= n ) then !layer j

    Kte_Sen = TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
      RTEdwSen(:,camad) * exp( uSen(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Kte_Cos = TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
      RTEdwCos(:,camad) * exp( uCos(:,camad) * ( z - prof(camad) - h(camad) ) ) )

    Ktedz_Sen = - AdmIntSen(:,camad) * TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
      RTEdwSen(:,camad) * exp( uSen(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Ktedz_Cos = - AdmIntCos(:,camad) * TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
      RTEdwCos(:,camad) * exp( uCos(:,camad) * ( z - prof(camad) - h(camad) ) ) )

    kernelEx = cmplx(0.d0,-ky,kind=dp) * Kte_Cos * w_cos
    Ex_ky = sum( kernelEx ) / ( pi * dabs(x) )

    kernelEy = - kxsen * Kte_Sen * w_sen
    Ey_ky = sum( kernelEy ) / ( pi * x )

    kernelHx = - kxsen *  Ktedz_Sen * w_sen
    Hx_ky = sum( kernelHx ) / ( pi * x )

    kernelHy = cmplx(0.d0,ky,kind=dp) * Ktedz_Cos * w_cos
    Hy_ky = sum( kernelHy ) / ( pi * dabs(x) )

    kernelHz = ( kxcos * kxcos + ky * ky ) * Kte_Cos * w_cos / zeta
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )

  else  !layer n

    Kte_Sen = TEdwSen(:,n) * exp( -uSen(:,n) * ( z - prof(n - 1) ) )
    Kte_Cos = TEdwCos(:,n) * exp( -uCos(:,n) * ( z - prof(n - 1) ) )

    Ktedz_Sen = - AdmIntSen(:,n) * Kte_Sen
    Ktedz_Cos = - AdmIntCos(:,n) * Kte_Cos

    kernelEx = cmplx(0.d0,-ky,kind=dp) * Kte_Cos * w_cos
    Ex_ky = sum( kernelEx ) / ( pi * dabs(x) )

    kernelEy = - kxsen * Kte_Sen * w_sen
    Ey_ky = sum( kernelEy ) / ( pi * x )

    kernelHx = - kxsen *  Ktedz_Sen * w_sen
    Hx_ky = sum( kernelHx ) / ( pi * x )

    kernelHy = cmplx(0.d0,ky,kind=dp) * Ktedz_Cos * w_cos
    Hy_ky = sum( kernelHy ) / ( pi * dabs(x) )

    kernelHz = ( kxcos * kxcos + ky * ky ) * Kte_Cos * w_cos / zeta
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )

  end if

  deallocate( h, KxSen, KxCos, w_Sen, w_Cos )
  deallocate( wvnb2, uSen, uCos, AdmIntSen, AdmIntCos, uhSen, uhCos, tghSen, tghCos )
  deallocate( AdmApdwSen, AdmApdwCos, RTEdwSen, RTEdwCos )
  deallocate( AdmApupSen, AdmApupCos, RTEupSen, RTEupCos )
  deallocate( FEdwSen, FEdwCos, FEupSen, FEupCos )
  deallocate( Kte_Sen, Kte_Cos, Ktedz_Sen, Ktedz_Cos )
  deallocate( kernelEx, kernelEy, kernelHx, kernelHy, kernelHz )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine vmd_xkyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end module vmd
