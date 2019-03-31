module hmdy
use parameters
use utils
use select_filter
    contains
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine hmdy_xyz_loops( Tx, Ty, h0, n, esp, condut, neta, zeta, cx, cy, z, Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p )
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: Tx, Ty, h0, esp(:), condut(1 : n), cx, cy, z
  complex(dp), intent(in) :: zeta, neta
  complex(dp), intent(out) :: Ex_p, Ey_p, Ez_p, Hx_p, Hy_p, Hz_p

  integer :: i, j, camad, camadT, filtro, idtfcd_cJ0, ident_fJ0, nJ0, idtfcd_cJ1, ident_fJ1, nJ1
  real(dp) :: x, y, r
  real(dp), dimension(:), allocatable :: h, krJ0, krJ1, w_J0, w_J1, prof

! Para uso de loops:
  complex(dp), dimension(:), allocatable :: wvnb2, AMdwJ0, AMdwJ1, AMupJ0, AMupJ1, FEdwJ0, FEdwJ1, FEupJ0, FEupJ1
  complex(dp), dimension(:,:), allocatable :: uJ0, uJ1, AdmIntJ0, AdmIntJ1, ImpIntJ0, ImpIntJ1, uhJ0, uhJ1, tghJ0, tghJ1
  complex(dp), dimension(:,:), allocatable :: AdmApdwJ0, AdmApdwJ1, ImpApdwJ0, ImpApdwJ1, AdmApupJ0, AdmApupJ1, ImpApupJ0, ImpApupJ1
  complex(dp), dimension(:,:), allocatable :: RTEdwJ0, RTEdwJ1, RTMdwJ0, RTMdwJ1, RTEupJ0, RTEupJ1, RTMupJ0, RTMupJ1
  complex(dp), dimension(:,:), allocatable :: TMdwJ0, TMdwJ1, TEdwJ0, TEdwJ1, TMupJ0, TMupJ1, TEupJ0, TEupJ1

  complex(dp), dimension(:), allocatable :: Ktmdz_J0, Ktmdz_J1, Ktm_J0, Ktm_J1, Kte_J0, Kte_J1, Ktedz_J0, Ktedz_J1
  complex(dp), dimension(:), allocatable :: kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1, kernelEzJ1
  complex(dp), dimension(:), allocatable :: kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1

  if ( dabs(cx - Tx) < eps .and. dabs(Tx) > eps ) then
    x = dsign( 1.d-1, Tx )
  elseif ( dabs(cx - Tx) < eps .and. dabs(Tx) < eps ) then
    x = 1.d-1
  else
    x = cx - Tx
  end if
  if ( dabs(cy - Ty) < eps .and. dabs(Ty) > eps ) then
    y = dsign(1.d-1,Ty)
!   else if ( dabs(cy - Ty) < eps .and. dabs(Ty) < eps ) then
!     y = 1.d-2
  else
    y = cy - Ty
  end if
  r = dsqrt( x ** 2 + y ** 2 )

  call sanitizedata(n, h0, z, esp, camadT, camad, h, prof)

!!  write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3), Key(4) ou Werthmüller (5)'
!!  read(*,*)idtfcd_cJ0
  filtro = 0  !esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
!!  call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
!!  write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3), Key(4) ou Werthmüller (5)'
!!  read(*,*)idtfcd_cJ1
!!  call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)

  idtfcd_cJ0 = 3 !5 !
  ident_fJ0 = 0
  nJ0 = 241 !201  !
  idtfcd_cJ1 = 3 !5 !
  ident_fJ1 = 1
  nJ1 = 241 !201  !

  allocate( KrJ0(nJ0), KrJ1(nJ1), w_J0(nJ0), w_J1(nJ1) )

  call constfiltro( filtro, idtfcd_cJ0, ident_fJ0, nJ0, r, KrJ0, w_J0 )
  call constfiltro( filtro, idtfcd_cJ1, ident_fJ1, nJ1, r, KrJ1, w_J1 )

  allocate( wvnb2(0 : n), uJ0(nJ0,0 : n), uJ1(nJ1,0 : n), AdmIntJ0(nJ0,0 : n), AdmIntJ1(nJ1,0 : n) )
  allocate( ImpIntJ0(nJ0,0 : n), ImpIntJ1(nJ1,0 : n) )
  allocate( uhJ0(nJ0,0 : n), uhJ1(nJ1,0 : n), tghJ0(nJ0,0 : n), tghJ1(nJ1,0 : n) )
! To workaround the warning: ... may be used uninitialized in this function
  allocate( TMdwJ0(1,1), TMdwJ1(1,1), TEdwJ0(1,1), TEdwJ1(1,1) )
  allocate( TMupJ0(1,1), TMupJ1(1,1), TEupJ0(1,1), TEupJ1(1,1) )
!
  wvnb2(0) = -zeta * neta
  uJ0(:,0) = sqrt( krJ0 * krJ0 - wvnb2(0) )
  uJ1(:,0) = sqrt( krJ1 * krJ1 - wvnb2(0) )
  AdmIntJ0(:,0) = uJ0(:,0) / zeta
  AdmIntJ1(:,0) = uJ1(:,0) / zeta
  ImpIntJ0(:,0) = uJ0(:,0) / neta
  ImpIntJ1(:,0) = uJ1(:,0) / neta
  uhJ0(:,0) = uJ0(:,0) * h(0)
  uhJ1(:,0) = uJ1(:,0) * h(0)
  tghJ0(:,0) = ( 1.d0 - exp( -2.d0 * uhJ0(:,0) ) ) / ( 1.d0 + exp( -2.d0 * uhJ0(:,0) ) )
  tghJ1(:,0) = ( 1.d0 - exp( -2.d0 * uhJ1(:,0) ) ) / ( 1.d0 + exp( -2.d0 * uhJ1(:,0) ) )
  do i = 1, n
    wvnb2(i) = -zeta * condut(i)
    uJ0(:,i) = sqrt( krJ0 * krJ0 - wvnb2(i) )
    uJ1(:,i) = sqrt( krJ1 * krJ1 - wvnb2(i) )
    AdmIntJ0(:,i) = uJ0(:,i) / zeta
    AdmIntJ1(:,i) = uJ1(:,i) / zeta
    ImpIntJ0(:,i) = uJ0(:,i) / condut(i)
    ImpIntJ1(:,i) = uJ1(:,i) / condut(i)
    uhJ0(:,i) = uJ0(:,i) * h(i)
    uhJ1(:,i) = uJ1(:,i) * h(i)
    tghJ0(:,i) = ( 1.d0 - exp( -2.d0 * uhJ0(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhJ0(:,i) ) )
    tghJ1(:,i) = ( 1.d0 - exp( -2.d0 * uhJ1(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhJ1(:,i) ) )
  end do

  allocate( AdmApdwJ0(nJ0,1 : n), AdmApdwJ1(nJ1,1 : n), ImpApdwJ0(nJ0,1 : n), ImpApdwJ1(nJ1,1 : n) )
  allocate( RTEdwJ0(nJ0,0:n), RTEdwJ1(nJ1,0 : n), RTMdwJ0(nJ0,0 : n), RTMdwJ1(nJ1,0 : n) )

  AdmApdwJ0(:,n) = AdmIntJ0(:,n)
  AdmApdwJ1(:,n) = AdmIntJ1(:,n)
  ImpApdwJ0(:,n) = ImpIntJ0(:,n)
  ImpApdwJ1(:,n) = ImpIntJ1(:,n)
  RTEdwJ0(:,n) = (0.d0,0.d0)
  RTEdwJ1(:,n) = (0.d0,0.d0)
  RTMdwJ0(:,n) = (0.d0,0.d0)
  RTMdwJ1(:,n) = (0.d0,0.d0)
  do i = n-1, 1, -1
    AdmApdwJ0(:,i) = AdmIntJ0(:,i) * ( AdmApdwJ0(:,i + 1) + AdmIntJ0(:,i) * &
                     tghJ0(:,i) ) / ( AdmIntJ0(:,i) + AdmApdwJ0(:,i + 1) * tghJ0(:,i) )
    AdmApdwJ1(:,i) = AdmIntJ1(:,i) * ( AdmApdwJ1(:,i + 1) + AdmIntJ1(:,i) * &
                     tghJ1(:,i) ) / ( AdmIntJ1(:,i) + AdmApdwJ1(:,i + 1) * tghJ1(:,i) )
    ImpApdwJ0(:,i) = ImpIntJ0(:,i) * ( ImpApdwJ0(:,i + 1) + ImpIntJ0(:,i) * &
                     tghJ0(:,i) ) / ( ImpIntJ0(:,i) + ImpApdwJ0(:,i + 1) * tghJ0(:,i) )
    ImpApdwJ1(:,i) = ImpIntJ1(:,i) * ( ImpApdwJ1(:,i + 1) + ImpIntJ1(:,i) * &
                     tghJ1(:,i) ) / ( ImpIntJ1(:,i) + ImpApdwJ1(:,i + 1) * tghJ1(:,i) )
    RTEdwJ0(:,i) = ( AdmIntJ0(:,i) - AdmApdwJ0(:,i + 1) ) / ( AdmIntJ0(:,i) + AdmApdwJ0(:,i + 1) )
    RTEdwJ1(:,i) = ( AdmIntJ1(:,i) - AdmApdwJ1(:,i + 1) ) / ( AdmIntJ1(:,i) + AdmApdwJ1(:,i + 1) )
    RTMdwJ0(:,i) = ( ImpIntJ0(:,i) - ImpApdwJ0(:,i + 1) ) / ( ImpIntJ0(:,i) + ImpApdwJ0(:,i + 1) )
    RTMdwJ1(:,i) = ( ImpIntJ1(:,i) - ImpApdwJ1(:,i + 1) ) / ( ImpIntJ1(:,i) + ImpApdwJ1(:,i + 1) )
  end do
  RTEdwJ0(:,0) = ( AdmIntJ0(:,0) - AdmApdwJ0(:,1) ) / ( AdmIntJ0(:,0) + AdmApdwJ0(:,1) )
  RTEdwJ1(:,0) = ( AdmIntJ1(:,0) - AdmApdwJ1(:,1) ) / ( AdmIntJ1(:,0) + AdmApdwJ1(:,1) )
  RTMdwJ0(:,0) = ( ImpIntJ0(:,0) - ImpApdwJ0(:,1) ) / ( ImpIntJ0(:,0) + ImpApdwJ0(:,1) )
  RTMdwJ1(:,0) = ( ImpIntJ1(:,0) - ImpApdwJ1(:,1) ) / ( ImpIntJ1(:,0) + ImpApdwJ1(:,1) )

  allocate( AdmApupJ0(nJ0,0 : n - 1), AdmApupJ1(nJ1,0 : n - 1), ImpApupJ0(nJ0,0 : n - 1), ImpApupJ1(nJ1,0 : n - 1) )
  allocate( RTEupJ0(nJ0,0 : n), RTEupJ1(nJ1,0 : n), RTMupJ0(nJ0,0 : n), RTMupJ1(nJ1,0 : n) )

  AdmApupJ0(:,0) = AdmIntJ0(:,0)
  AdmApupJ1(:,0) = AdmIntJ1(:,0)
  ImpApupJ0(:,0) = ImpIntJ0(:,0)
  ImpApupJ1(:,0) = ImpIntJ1(:,0)
  RTEupJ0(:,0) = (0.d0,0.d0)
  RTEupJ1(:,0) = (0.d0,0.d0)
  RTMupJ0(:,0) = (0.d0,0.d0)
  RTMupJ1(:,0) = (0.d0,0.d0)
  do i = 1, n - 1
    AdmApupJ0(:,i) = AdmIntJ0(:,i) * ( AdmApupJ0(:,i - 1) + AdmIntJ0(:,i) * &
                    tghJ0(:,i) ) / ( AdmIntJ0(:,i) + AdmApupJ0(:,i - 1) * tghJ0(:,i) )
    AdmApupJ1(:,i) = AdmIntJ1(:,i) * ( AdmApupJ1(:,i - 1) + AdmIntJ1(:,i) * &
                    tghJ1(:,i) ) / ( AdmIntJ1(:,i) + AdmApupJ1(:,i - 1) * tghJ1(:,i) )
    ImpApupJ0(:,i) = ImpIntJ0(:,i) * ( ImpApupJ0(:,i - 1) + ImpIntJ0(:,i) * &
                    tghJ0(:,i) ) / ( ImpIntJ0(:,i) + ImpApupJ0(:,i - 1) * tghJ0(:,i) )
    ImpApupJ1(:,i) = ImpIntJ1(:,i) * ( ImpApupJ1(:,i - 1) + ImpIntJ1(:,i) * &
                    tghJ1(:,i) ) / ( ImpIntJ1(:,i) + ImpApupJ1(:,i - 1) * tghJ1(:,i) )
    RTEupJ0(:,i) = ( AdmIntJ0(:,i) - AdmApupJ0(:,i - 1) ) / ( AdmIntJ0(:,i) + AdmApupJ0(:,i - 1) )
    RTEupJ1(:,i) = ( AdmIntJ1(:,i) - AdmApupJ1(:,i - 1) ) / ( AdmIntJ1(:,i) + AdmApupJ1(:,i - 1) )
    RTMupJ0(:,i) = ( ImpIntJ0(:,i) - ImpApupJ0(:,i - 1) ) / ( ImpIntJ0(:,i) + ImpApupJ0(:,i - 1) )
    RTMupJ1(:,i) = ( ImpIntJ1(:,i) - ImpApupJ1(:,i - 1) ) / ( ImpIntJ1(:,i) + ImpApupJ1(:,i - 1) )
  end do
  RTEupJ0(:,n) = ( AdmIntJ0(:,n) - AdmApupJ0(:,n - 1) ) / ( AdmIntJ0(:,n) + AdmApupJ0(:,n - 1) )
  RTEupJ1(:,n) = ( AdmIntJ1(:,n) - AdmApupJ1(:,n - 1) ) / ( AdmIntJ1(:,n) + AdmApupJ1(:,n - 1) )
  RTMupJ0(:,n) = ( ImpIntJ0(:,n) - ImpApupJ0(:,n - 1) ) / ( ImpIntJ0(:,n) + ImpApupJ0(:,n - 1) )
  RTMupJ1(:,n) = ( ImpIntJ1(:,n) - ImpApupJ1(:,n - 1) ) / ( ImpIntJ1(:,n) + ImpApupJ1(:,n - 1) )

  allocate( AMdwJ0(nJ0), AMdwJ1(nJ1), AMupJ0(nJ0), AMupJ1(nJ1), FEdwJ0(nJ0), FEdwJ1(nJ1), FEupJ0(nJ0), FEupJ1(nJ1) )

  AMdwJ0 = ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + RTMupJ0(:,camadT) * &
            exp( uJ0(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupJ0(:,camadT) * RTMdwJ0(:,camadT) * exp( -2.d0 * uhJ0(:,camadT) ) )
  AMdwJ1 = ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) + RTMupJ1(:,camadT) * &
            exp( uJ1(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupJ1(:,camadT) * RTMdwJ1(:,camadT) * exp( -2.d0 * uhJ1(:,camadT) ) )

  AMupJ0 = ( exp( uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTMdwJ0(:,camadT) * &
            exp( -uJ0(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupJ0(:,camadT) * RTMdwJ0(:,camadT) * exp( -2.d0 * uhJ0(:,camadT) ) )
  AMupJ1 = ( exp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTMdwJ1(:,camadT) * &
            exp( -uJ1(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupJ1(:,camadT) * RTMdwJ1(:,camadT) * exp( -2.d0 * uhJ1(:,camadT) ) )

  FEdwJ0 = ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) - RTEupJ0(:,camadT) * &
            exp( uJ0(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ0(:,camadT) * RTEdwJ0(:,camadT) * exp( -2.d0 * uhJ0(:,camadT) ) )
  FEdwJ1 = ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) - RTEupJ1(:,camadT) * &
            exp( uJ1(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ1(:,camadT) * RTEdwJ1(:,camadT) * exp( -2.d0 * uhJ1(:,camadT) ) )

  FEupJ0 = ( exp( uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) - RTEdwJ0(:,camadT) * &
            exp( -uJ0(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ0(:,camadT) * RTEdwJ0(:,camadT) * exp( -2.d0 * uhJ0(:,camadT) ) )
  FEupJ1 = ( exp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) - RTEdwJ1(:,camadT) * &
            exp( -uJ1(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupJ1(:,camadT) * RTEdwJ1(:,camadT) * exp( -2.d0 * uhJ1(:,camadT) ) )

  if ( camad > camadT ) then
    deallocate( TMdwJ0, TMdwJ1, TEdwJ0, TEdwJ1)
    allocate( TMdwJ0(nJ0,camadT : camad), TMdwJ1(nJ1,camadT : camad), TEdwJ0(nJ0,camadT : camad), TEdwJ1(nJ1,camadT : camad) )
    do j = camadT, camad
      if ( j == camadT ) then
        TMdwJ0(:,j) = zeta * my / ( 2.d0 * ImpIntJ0(:,camadT) )
        TMdwJ1(:,j) = zeta * my / ( 2.d0 * ImpIntJ1(:,camadT) )
        TEdwJ0(:,j) = - zeta * my / 2.d0
        TEdwJ1(:,j) = - zeta * my / 2.d0
      else if ( j == (camadT + 1) .and. j == n ) then
        TMdwJ0(:,j) = TMdwJ0(:,j - 1) * ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupJ0(:,camadT) * AMupJ0(:) * exp( -uhJ0(:,camadT) ) + RTMdwJ0(:,camadT) * AMdwJ0(:) )
        TMdwJ1(:,j) = TMdwJ1(:,j - 1) * ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupJ1(:,camadT) * AMupJ1(:) * exp( -uhJ1(:,camadT) ) + RTMdwJ1(:,camadT) * AMdwJ1(:) )
        TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) - &
                      RTEupJ0(:,camadT) * FEupJ0(:) * exp( -uhJ0(:,camadT) ) + RTEdwJ0(:,camadT) * FEdwJ0(:) )
        TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) - &
                      RTEupJ1(:,camadT) * FEupJ1(:) * exp( -uhJ1(:,camadT) ) + RTEdwJ1(:,camadT) * FEdwJ1(:) )
      else if ( j == (camadT + 1) .and. j /= n ) then
        TMdwJ0(:,j) = TMdwJ0(:,j - 1) * ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupJ0(:,camadT) * AMupJ0(:) * exp( -uhJ0(:,camadT) ) + &
                      RTMdwJ0(:,camadT) * AMdwJ0(:) ) / ( 1.d0 + RTMdwJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TMdwJ1(:,j) = TMdwJ1(:,j - 1) * ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupJ1(:,camadT) * AMupJ1(:) * exp( -uhJ1(:,camadT) ) + &
                      RTMdwJ1(:,camadT) * AMdwJ1(:) ) / ( 1.d0 + RTMdwJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
        TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( exp( -uJ0(:,camadT) * ( prof(camadT) - h0 ) ) - &
                      RTEupJ0(:,camadT) * FEupJ0(:) * exp( -uhJ0(:,camadT) ) + &
                      RTEdwJ0(:,camadT) * FEdwJ0(:) ) / ( 1.d0 + RTEdwJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( exp( -uJ1(:,camadT) * ( prof(camadT) - h0 ) ) - &
                      RTEupJ1(:,camadT) * FEupJ1(:) * exp( -uhJ1(:,camadT) ) + &
                      RTEdwJ1(:,camadT) * FEdwJ1(:) ) / ( 1.d0 + RTEdwJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
      else if ( j /= n ) then
        TMdwJ0(:,j) = TMdwJ0(:,j - 1) * ( 1.d0 + RTMdwJ0(:,j - 1) ) * exp( -uhJ0(:,j - 1) ) / &
                      ( 1.d0 + RTMdwJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TMdwJ1(:,j) = TMdwJ1(:,j - 1) * ( 1.d0 + RTMdwJ1(:,j - 1) ) * exp( -uhJ1(:,j - 1) ) / &
                      ( 1.d0 + RTMdwJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
        TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( 1.d0 + RTEdwJ0(:,j - 1) ) * exp( -uhJ0(:,j - 1) ) / &
                      ( 1.d0 + RTEdwJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( 1.d0 + RTEdwJ1(:,j - 1) ) * exp( -uhJ1(:,j - 1) ) / &
                      ( 1.d0 + RTEdwJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
      else if ( j == n ) then
        TMdwJ0(:,j) = TMdwJ0(:,j - 1) * ( 1.d0 + RTMdwJ0(:,j - 1) ) * exp( -uhJ0(:,j - 1) )
        TMdwJ1(:,j) = TMdwJ1(:,j - 1) * ( 1.d0 + RTMdwJ1(:,j - 1) ) * exp( -uhJ1(:,j - 1) )
        TEdwJ0(:,j) = TEdwJ0(:,j - 1) * ( 1.d0 + RTEdwJ0(:,j - 1) ) * exp( -uhJ0(:,j - 1) )
        TEdwJ1(:,j) = TEdwJ1(:,j - 1) * ( 1.d0 + RTEdwJ1(:,j - 1) ) * exp( -uhJ1(:,j - 1) )
      end if
    end do
  else if ( camad < camadT ) then
    deallocate( TMupJ0, TMupJ1, TEupJ0, TEupJ1)
    allocate( TMupJ0(nJ0,camad : camadT), TMupJ1(nJ1,camad : camadT), TEupJ0(nJ0,camad : camadT), TEupJ1(nJ1,camad : camadT) )
    do j = camadT, camad, -1
      if ( j == camadT ) then
        TMupJ0(:,j) = zeta * my / ( 2.d0 * ImpIntJ0(:,camadT) )
        TMupJ1(:,j) = zeta * my / ( 2.d0 * ImpIntJ1(:,camadT) )
        TEupJ0(:,j) = zeta * my / 2.d0
        TEupJ1(:,j) = zeta * my / 2.d0
      else if ( j == (camadT - 1) .and. j == 0 ) then
        TMupJ0(:,j) = TMupJ0(:,j + 1) * ( exp( -uJ0(:,camadT) * h0 ) + &
                      RTMupJ0(:,camadT) * AMupJ0(:) + RTMdwJ0(:,camadT) * AMdwJ0(:) * exp( -uhJ0(:,camadT) ) )
        TMupJ1(:,j) = TMupJ1(:,j + 1) * ( exp( -uJ1(:,camadT) * h0 ) + &
                      RTMupJ1(:,camadT) * AMupJ1(:) + RTMdwJ1(:,camadT) * AMdwJ1(:) * exp( -uhJ1(:,camadT) ) )
        TEupJ0(:,j) = TEupJ0(:,j + 1) * ( exp( -uJ0(:,camadT) * h0 ) + &
                      RTEupJ0(:,camadT) * FEupJ0(:) - RTEdwJ0(:,camadT) * FEdwJ0(:) * exp( -uhJ0(:,camadT) ) )
        TEupJ1(:,j) = TEupJ1(:,j + 1) * ( exp( -uJ1(:,camadT) * h0 ) + &
                      RTEupJ1(:,camadT) * FEupJ1(:) - RTEdwJ1(:,camadT) * FEdwJ1(:) * exp( -uhJ1(:,camadT) ) )
      else if ( j == (camadT - 1) .and. j /= 0 ) then
        TMupJ0(:,j) = TMupJ0(:,j + 1) * ( exp( uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                      RTMupJ0(:,camadT) * AMupJ0(:) + RTMdwJ0(:,camadT) * AMdwJ0(:)* &
                      exp( -uhJ0(:,camadT) ) ) / ( 1.d0 + RTMupJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TMupJ1(:,j) = TMupJ1(:,j + 1) * ( exp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                      RTMupJ1(:,camadT) * AMupJ1(:) + RTMdwJ1(:,camadT) * AMdwJ1(:) * &
                      exp( -uhJ1(:,camadT) ) ) / ( 1.d0 + RTMupJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
        TEupJ0(:,j) = TEupJ0(:,j + 1) * ( exp( uJ0(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                      RTEupJ0(:,camadT) * FEupJ0(:) - RTEdwJ0(:,camadT) * FEdwJ0(:) * &
                      exp( -uhJ0(:,camadT) ) ) / ( 1.d0 + RTEupJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TEupJ1(:,j) = TEupJ1(:,j + 1) * ( exp( uJ1(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
                      RTEupJ1(:,camadT) * FEupJ1(:) - RTEdwJ1(:,camadT) * FEdwJ1(:) * &
                      exp( -uhJ1(:,camadT) ) ) / ( 1.d0 + RTEupJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
      else if ( j /= 0 ) then
        TMupJ0(:,j) = TMupJ0(:,j + 1) * ( 1.d0 + RTMupJ0(:,j + 1) ) * exp( -uhJ0(:,j + 1) ) / &
                      ( 1.d0 + RTMupJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TMupJ1(:,j) = TMupJ1(:,j + 1) * ( 1.d0 + RTMupJ1(:,j + 1) ) * exp( -uhJ1(:,j + 1) ) / &
                      ( 1.d0 + RTMupJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
        TEupJ0(:,j) = TEupJ0(:,j + 1) * ( 1.d0 + RTEupJ0(:,j + 1) ) * exp( -uhJ0(:,j + 1) ) / &
                      ( 1.d0 + RTEupJ0(:,j) * exp( -2.d0 * uhJ0(:,j) ) )
        TEupJ1(:,j) = TEupJ1(:,j + 1) * ( 1.d0 + RTEupJ1(:,j + 1) ) * exp( -uhJ1(:,j + 1) ) / &
                      ( 1.d0 + RTEupJ1(:,j) * exp( -2.d0 * uhJ1(:,j) ) )
      else if ( j == 0 ) then
        TMupJ0(:,j) = TMupJ0(:,1) * ( 1.d0 + RTMupJ0(:,1) ) * exp( -uhJ0(:,1) )
        TMupJ1(:,j) = TMupJ1(:,1) * ( 1.d0 + RTMupJ1(:,1) ) * exp( -uhJ1(:,1) )
        TEupJ0(:,j) = TEupJ0(:,1) * ( 1.d0 + RTEupJ0(:,1) ) * exp( -uhJ0(:,1) )
        TEupJ1(:,j) = TEupJ1(:,1) * ( 1.d0 + RTEupJ1(:,1) ) * exp( -uhJ1(:,1) )
      end if
    end do
  else
    deallocate( TMdwJ0, TMdwJ1, TEdwJ0, TEdwJ1)
    allocate( TMdwJ0(nJ0,camadT : camad), TMdwJ1(nJ1,camadT : camad), TEdwJ0(nJ0,camadT : camad), TEdwJ1(nJ1,camadT : camad) )
    allocate( TMupJ0(nJ0,camad : camadT), TMupJ1(nJ1,camad : camadT), TEupJ0(nJ0,camad : camadT), TEupJ1(nJ1,camad : camadT) )

      TMdwJ0(:,camad) = zeta * my / ( 2.d0 * ImpIntJ0(:,camadT) )
      TMdwJ1(:,camad) = zeta * my / ( 2.d0 * ImpIntJ1(:,camadT) )

      TEdwJ0(:,camad) = - zeta * my / 2.d0
      TEdwJ1(:,camad) = - zeta * my / 2.d0

      TMupJ0(:,camad) = TMdwJ0(:,camad)
      TMupJ1(:,camad) = TMdwJ1(:,camad)

      TEupJ0(:,camad) = - TEdwJ0(:,camad)
      TEupJ1(:,camad) = - TEdwJ1(:,camad)
  end if

  allocate( Ktmdz_J0(nJ0), Ktmdz_J1(nJ1), Ktm_J0(nJ0), Ktm_J1(nJ1), Kte_J0(nJ0), Kte_J1(nJ1), Ktedz_J0(nJ0), Ktedz_J1(nJ1) )
  allocate( kernelExJ0(nJ0), kernelExJ1(nJ1), kernelEyJ0(nJ0), kernelEyJ1(nJ1), kernelEzJ1(nJ1))
  allocate( kernelHxJ0(nJ0), kernelHxJ1(nJ1), kernelHyJ0(nJ0), kernelHyJ1(nJ1), kernelHzJ1(nJ1))
  if ( camad == 0 .and. camadT /= 0 ) then
    Ktmdz_J0 = ( ImpIntJ0(:,0) * TMupJ0(:,0) * exp( uJ0(:,0) * z ) ) * w_J0(:)
    Ktmdz_J1 = ( ImpIntJ1(:,0) * TMupJ1(:,0) * exp( uJ1(:,0) * z ) ) * w_J1(:)
    Kte_J0 = ( TEupJ0(:,0) * exp( uJ0(:,0) * z ) ) * w_J0(:)
    Kte_J1 = ( TEupJ1(:,0) * exp( uJ1(:,0) * z ) ) * w_J1(:)

    Ktm_J0 = ( TMupJ0(:,0) * exp( uJ0(:,0) * z ) ) * w_J0(:)
    Ktm_J1 = ( TMupJ1(:,0) * exp( uJ1(:,0) * z ) ) * w_J1(:)
    Ktedz_J0 = ( AdmIntJ0(:,0) * TEupJ0(:,0) * exp( uJ0(:,0) * z ) ) * w_J0(:)
    Ktedz_J1 = ( AdmIntJ1(:,0) * TEupJ1(:,0) * exp( uJ1(:,0) * z ) ) * w_J1(:)

    kernelExJ1 = - y * x * ( Ktmdz_J1 - Kte_J1 ) / ( r ** 3 )
    kernelExJ0 = - y * x * ( Ktmdz_J0 - Kte_J0 ) * krJ0 / ( 2.d0 * r ** 2 )
    Ey_p = ( sum( kernelExJ1 ) - sum( kernelExJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelEyJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktmdz_J1 + ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Kte_J1
    kernelEyJ0 = ( x * x * Ktmdz_J0 + y * y * Kte_J0 ) * krJ0 / r ** 2
    Ex_p = ( sum( kernelEyJ1 ) + sum( kernelEyJ0 ) ) / ( 2.d0 * pi * r )  !este último r é decorrente do uso dos filtros

    kernelEzJ1 = ( x / ( r * neta ) ) * Ktm_J1 * krJ1 * krJ1
    Ez_p = sum( kernelEzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros

    kernelHxJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktm_J1 + ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Ktedz_J1
    kernelHxJ0 = ( x * x * Ktm_J0 + y * y * Ktedz_J0 ) * krJ0 / r ** 2
    Hy_p = -( sum( kernelHxJ1 ) + sum( kernelHxJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelHyJ1 = y * x * ( Ktm_J1 - Ktedz_J1 ) / r ** 3
    kernelHyJ0 = y * x * ( Ktm_J0 - Ktedz_J0 ) * krJ0 / ( 2.d0 * r * r )
    Hx_p = ( - sum( kernelHyJ1 ) + sum( kernelHyJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelHzJ1 = y/(r*zeta) * Kte_J1 * krJ1 * krJ1
    Hz_p = - sum( kernelHzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros
  else if ( camad < camadT ) then !camada k
    ktmdz_J0 = ( ImpIntJ0(:,camad) * TMupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - prof(camad) ) ) - &
              RTMupJ0(:,camad) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J0(:)
    ktmdz_J1 = ( ImpIntJ1(:,camad) * TMupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - prof(camad) ) ) - &
              RTMupJ1(:,camad) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)
    Kte_J0 = ( TEupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - prof(camad) ) ) + &
              RTEupJ0(:,camad) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - prof(camad) ) ) + &
              RTEupJ1(:,camad) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)

    ktm_J0 = ( TMupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - prof(camad) ) ) + &
              RTMupJ0(:,camad) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J0(:)
    ktm_J1 = ( TMupJ1(:,camad) * ( exp( uJ1(:,camad)  * ( z - prof(camad) ) ) + &
              RTMupJ1(:,camad) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)
    Ktedz_J0 = ( AdmIntJ0(:,camad) * TEupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - prof(camad) ) ) - &
              RTEupJ0(:,camad) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J0(:)
    Ktedz_J1 = ( AdmIntJ1(:,camad) * TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - prof(camad) ) ) - &
              RTEupJ1(:,camad) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = -y * x * ( Ktmdz_J1 - Kte_J1 ) / ( r * r * r )
    kernelExJ0 = y * x * ( Ktmdz_J0 - Kte_J0 ) * krJ0 / ( 2.d0 * r * r )
    Ey_p = ( sum( kernelExJ1 ) + sum( kernelExJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelEyJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktmdz_J1 + ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Kte_J1
    kernelEyJ0 = ( x * x * Ktmdz_J0 + y * y * Kte_J0 ) * krJ0 / r ** 2
    Ex_p = ( sum( kernelEyJ1 ) + sum( kernelEyJ0 ) ) / ( 2.d0 * pi * r )  !este último r é decorrente do uso dos filtros

    kernelEzJ1 = ( x / ( r * condut( camad ) ) ) * Ktm_J1 * krJ1 * krJ1
    Ez_p = sum( kernelEzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros

    kernelHxJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktm_J1 + ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Ktedz_J1
    kernelHxJ0 = ( x * x * Ktm_J0 + y * y * Ktedz_J0 ) * krJ0 / r ** 2
    Hy_p = -( sum( kernelHxJ1 ) + sum( kernelHxJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelHyJ1 = y * x * ( Ktm_J1 - Ktedz_J1 ) / r ** 3
    kernelHyJ0 = y * x * ( Ktm_J0 - Ktedz_J0 ) * krJ0 / ( 2.d0 * r * r )
    Hx_p = ( -sum( kernelHyJ1 ) + sum( kernelHyJ0 ) ) / ( pi * r )  !este último r é decorrente do uso dos filtros

    kernelHzJ1 = y / ( r * zeta ) * Kte_J1 * krJ1 * krJ1
    Hz_p = - sum( kernelHzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros
  else if ( camad == camadT .and. z <= h0 ) then  !na mesma camada do transmissor mas acima dele
    Ktmdz_J0 = ( ImpIntJ0(:,camad) * TMupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - h0 ) ) - &
              RTMupJ0(:,camad) * AMupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwJ0(:,camad) * AMdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Ktmdz_J1 = ( ImpIntJ1(:,camad) * TMupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - h0 ) ) - &
              RTMupJ1(:,camad) * AMupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwJ1(:,camad) * AMdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
    Kte_J0 = ( TEupJ0(:,camad) * ( exp(uJ0(:,camad) * ( z - h0 ) ) + &
              RTEupJ0(:,camad) * FEupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwJ0(:,camad) * FEdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - h0 ) ) + &
              RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    Ktm_J0 = ( TMupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - h0 ) ) + &
              RTMupJ0(:,camad) * AMupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwJ0(:,camad) * AMdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Ktm_J1 = ( TMupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - h0 ) ) + &
              RTMupJ1(:,camad) * AMupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwJ1(:,camad) * AMdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)
    Ktedz_J0 = ( AdmIntJ0(:,camad) * TEupJ0(:,camad) * ( exp( uJ0(:,camad) * ( z - h0 ) ) - &
              RTEupJ0(:,camad) * FEupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwJ0(:,camad) * FEdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Ktedz_J1 = ( AdmIntJ1(:,camad) * TEupJ1(:,camad) * ( exp( uJ1(:,camad) * ( z - h0 ) ) - &
              RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = -y * x * ( Ktmdz_J1 - Kte_J1 ) / (r * r * r)
    kernelExJ0 = y * x * ( Ktmdz_J0 - Kte_J0 ) * krJ0 / ( 2.d0 * r * r )
    Ey_p = ( sum( kernelExJ1 ) + sum( kernelExJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelEyJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktmdz_J1 + ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Kte_J1
    kernelEyJ0 = ( x * x * Ktmdz_J0 + y * y * Kte_J0 ) * krJ0 / r ** 2
    Ex_p = ( sum( kernelEyJ1 ) + sum( kernelEyJ0 ) ) / ( 2.d0 * pi * r )  !este último r é decorrente do uso dos filtros

    if ( camad /= 0 ) then
      kernelEzJ1 = ( x / ( r * condut(camad) ) ) * Ktm_J1 * krJ1 * krJ1
    else
      kernelEzJ1 = ( x / ( r * neta ) ) * Ktm_J1 * krJ1 * krJ1
    end if
    Ez_p = sum( kernelEzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros

    kernelHxJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktm_J1 + ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Ktedz_J1
    kernelHxJ0 = ( x * x * Ktm_J0 + y * y * Ktedz_J0 ) * krJ0 / r ** 2
    Hy_p = -( sum( kernelHxJ1 ) + sum( kernelHxJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelHyJ1 = y * x * ( Ktm_J1 - Ktedz_J1 ) / r ** 3
    kernelHyJ0 = y * x * ( Ktm_J0 - Ktedz_J0 ) * krJ0 / ( 2.d0 * r * r )
    Hx_p = ( -sum( kernelHyJ1 ) + sum( kernelHyJ0 ) ) / ( pi * r )  !este último r é decorrente do uso dos filtros

    kernelHzJ1 = y / ( r * zeta ) * Kte_J1 * krJ1 * krJ1
    Hz_p = - sum( kernelHzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros
  else if ( camad == camadT .and. z > h0 ) then !na mesma camada do transmissor mas abaixo dele
    Ktmdz_J0 = ( ImpIntJ0(:,camad) * TMdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - h0 ) ) + &
               RTMupJ0(:,camad) * AMupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) - &
               RTMdwJ0(:,camad) * AMdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Ktmdz_J1 = ( ImpIntJ1(:,camad) * TMdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - h0 ) ) + &
               RTMupJ1(:,camad) * AMupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTMdwJ1(:,camad) * AMdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    Kte_J0 = ( TEdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - h0 ) ) - &
              RTEupJ0(:,camad) * FEupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTEdwJ0(:,camad) * FEdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - h0 ) ) - &
              RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    Ktm_J0 = ( TMdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - h0 ) ) + &
              RTMupJ0(:,camad) * AMupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwJ0(:,camad) * AMdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Ktm_J1 = ( TMdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - h0 ) ) + &
              RTMupJ1(:,camad) * AMupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwJ1(:,camad) * AMdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    Ktedz_J0 = ( AdmIntJ0(:,camad) * TEdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - h0 ) ) - &
              RTEupJ0(:,camad) * FEupJ0(:) * exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwJ0(:,camad) * FEdwJ0(:) * exp( uJ0(:,camad) * ( z - prof(camad) ) ) ) ) * w_J0(:)
    Ktedz_J1 = ( AdmIntJ1(:,camad) * TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - h0 ) ) - &
              RTEupJ1(:,camad) * FEupJ1(:) * exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwJ1(:,camad) * FEdwJ1(:) * exp( uJ1(:,camad) * ( z - prof(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = y * x * ( Ktmdz_J1 + Kte_J1 ) / ( r * r * r )
    kernelExJ0 = y * x * ( Ktmdz_J0 + Kte_J0 ) * krJ0 / ( 2.d0 * r * r )
    Ey_p = ( sum( kernelExJ1 ) - sum( kernelExJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelEyJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktmdz_J1 - ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Kte_J1
    kernelEyJ0 = ( x * x * Ktmdz_J0 - y * y * Kte_J0 ) * krJ0 / r ** 2
    Ex_p = -( sum( kernelEyJ1 ) + sum( kernelEyJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    if ( camad /= 0 ) then
      kernelEzJ1 = ( x / ( r * condut(camad) ) ) * Ktm_J1 * krJ1 * krJ1
    else
      kernelEzJ1 = ( x / ( r * neta ) ) * Ktm_J1 * krJ1 * krJ1
    end if
    Ez_p = sum( kernelEzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros

    kernelHxJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktm_J1 - ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Ktedz_J1
    kernelHxJ0 = ( x * x * Ktm_J0 - y * y * Ktedz_J0 ) * krJ0 / r ** 2
    Hy_p = -( sum( kernelHxJ1 ) + sum( kernelHxJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelHyJ1 = y * x * ( Ktm_J1 + Ktedz_J1 ) / r ** 3
    kernelHyJ0 = y * x * ( Ktm_J0 + Ktedz_J0 ) * krJ0 / ( 2.d0 * r * r )
    Hx_p = ( -sum( kernelHyJ1 ) + sum( kernelHyJ0 ) ) / ( pi * r )  !este último r é decorrente do uso dos filtros

    kernelHzJ1 = y / ( r * zeta ) * Kte_J1 * krJ1 * krJ1
    Hz_p = - sum( kernelHzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros
  else if ( camad > camadT .and. camad /= n ) then !camada j
    Ktmdz_J0 = ( ImpIntJ0(:,camad) * TMdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTMdwJ0(:,camad) * exp( uJ0(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J0(:)
    Ktmdz_J1 = ( ImpIntJ1(:,camad) * TMdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTMdwJ1(:,camad) * exp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)

    Kte_J0 = ( TEdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ0(:,camad) * exp( uJ0(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J0(:)
    Kte_J1 = ( TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwJ1(:,camad) * exp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)

    Ktm_J0 = ( TMdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTMdwJ0(:,camad) * exp( uJ0(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J0(:)
    Ktm_J1 = ( TMdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTMdwJ1(:,camad) * exp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)

    Ktedz_J0 = ( AdmIntJ0(:,camad) * TEdwJ0(:,camad) * ( exp( -uJ0(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwJ0(:,camad) * exp( uJ0(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J0(:)
    Ktedz_J1 = ( AdmIntJ1(:,camad) * TEdwJ1(:,camad) * ( exp( -uJ1(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwJ1(:,camad) * exp( uJ1(:,camad) * ( z - prof(camad) - h(camad) ) ) ) ) * w_J1(:)

    kernelExJ1 = y * x * ( Ktmdz_J1 + Kte_J1 ) / ( r * r * r )
    kernelExJ0 = y * x * ( Ktmdz_J0 + Kte_J0 ) * krJ0 / ( 2.d0 * r * r )
    Ey_p = ( sum( kernelExJ1 ) - sum( kernelExJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelEyJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktmdz_J1 - ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Kte_J1
    kernelEyJ0 = ( x * x * Ktmdz_J0 - y * y * Kte_J0 ) * krJ0 / r ** 2
    Ex_p = -( sum( kernelEyJ1 ) + sum( kernelEyJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelEzJ1 = ( x / ( r * condut(camad) ) ) * Ktm_J1 * krJ1 * krJ1
    Ez_p = sum( kernelEzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros

    kernelHxJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktm_J1 - ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Ktedz_J1
    kernelHxJ0 = ( x * x * Ktm_J0 - y * y * Ktedz_J0 ) * krJ0 / r ** 2
    Hy_p = -( sum( kernelHxJ1 ) + sum( kernelHxJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelHyJ1 = y * x * ( Ktm_J1 + Ktedz_J1 ) / r ** 3
    kernelHyJ0 = y * x * ( Ktm_J0 + Ktedz_J0 ) * krJ0 / ( 2.d0 * r * r )
    Hx_p = ( -sum( kernelHyJ1 ) + sum( kernelHyJ0 ) ) / ( pi * r )  !este último r é decorrente do uso dos filtros

    kernelHzJ1 = y / ( r * zeta ) * Kte_J1 * krJ1 * krJ1
    Hz_p = - sum( kernelHzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros
  else  !camada n
    Ktmdz_J0 = ( ImpIntJ0(:,n) * TMdwJ0(:,n) * exp( -uJ0(:,n) * ( z - prof(n - 1) ) ) ) * w_J0(:)
    Ktmdz_J1 = ( ImpIntJ1(:,n) * TMdwJ1(:,n) * exp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)

    Kte_J0 = ( TEdwJ0(:,n) * exp( -uJ0(:,n) * ( z - prof(n - 1) ) ) ) * w_J0(:)
    Kte_J1 = ( TEdwJ1(:,n) * exp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)

    Ktm_J0 = ( TMdwJ0(:,n) * exp( -uJ0(:,n) * ( z - prof(n - 1) ) ) ) * w_J0(:)
    Ktm_J1 = ( TMdwJ1(:,n) * exp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)

    Ktedz_J0 = ( AdmIntJ0(:,n) * TEdwJ0(:,n) * exp( -uJ0(:,n) * ( z - prof(n - 1) ) ) ) * w_J0(:)
    Ktedz_J1 = ( AdmIntJ1(:,n) * TEdwJ1(:,n) * exp( -uJ1(:,n) * ( z - prof(n - 1) ) ) ) * w_J1(:)

    kernelExJ1 = y * x * ( Ktmdz_J1 + Kte_J1 ) / ( r * r * r )
    kernelExJ0 = y * x * ( Ktmdz_J0 + Kte_J0 ) * krJ0 / ( 2.d0 * r * r )
    Ey_p = ( sum( kernelExJ1 ) - sum( kernelExJ0 ) ) / ( pi * r ) !este último r é decorrente do uso dos filtros

    kernelEyJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktmdz_J1 - ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Kte_J1
    kernelEyJ0 = ( x * x * Ktmdz_J0 - y * y * Kte_J0 ) * krJ0 / r ** 2
    Ex_p = -( sum( kernelEyJ1 ) + sum( kernelEyJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelEzJ1 = ( x / ( r * condut(n) ) ) * Ktm_J1 * krJ1 * krJ1
    Ez_p = sum( kernelEzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros

    kernelHxJ1 = ( 1.d0 / r - 2.d0 * x ** 2 / r ** 3 ) * Ktm_J1 - ( 1.d0 / r - 2.d0 * y ** 2 / r ** 3 ) * Ktedz_J1
    kernelHxJ0 = ( x * x * Ktm_J0 - y * y * Ktedz_J0 ) * krJ0 / r ** 2
    Hy_p = -( sum( kernelHxJ1 ) + sum( kernelHxJ0 ) ) / ( 2.d0 * pi * r ) !este último r é decorrente do uso dos filtros

    kernelHyJ1 = y * x * ( Ktm_J1 + Ktedz_J1 ) / r ** 3
    kernelHyJ0 = y * x * ( Ktm_J0 + Ktedz_J0 ) * krJ0 / ( 2.d0 * r * r )
    Hx_p = ( -sum( kernelHyJ1 ) + sum( kernelHyJ0 ) ) / ( pi * r )  !este último r é decorrente do uso dos filtros

    kernelHzJ1 = y / ( r * zeta ) * Kte_J1 * krJ1 * krJ1
    Hz_p = - sum( kernelHzJ1 ) / ( 2.d0 * pi * r )        !este último r é decorrente do uso dos filtros
  end if

  deallocate( h, KrJ0, KrJ1, w_J0, w_J1 )
  deallocate( wvnb2, uJ0, uJ1, AdmIntJ0, AdmIntJ1, ImpIntJ0, ImpIntJ1, uhJ0, uhJ1, tghJ0, tghJ1 )
  deallocate( AdmApdwJ0, AdmApdwJ1, ImpApdwJ0, ImpApdwJ1, RTEdwJ0, RTEdwJ1, RTMdwJ0, RTMdwJ1 )
  deallocate( AdmApupJ0, AdmApupJ1, ImpApupJ0, ImpApupJ1, RTEupJ0, RTEupJ1, RTMupJ0, RTMupJ1 )
  deallocate( AMdwJ0, AMdwJ1, AMupJ0, AMupJ1, FEdwJ0, FEdwJ1, FEupJ0, FEupJ1 )
  deallocate( Ktmdz_J0, Ktmdz_J1, Ktm_J0, Ktm_J1, Kte_J0, Kte_J1, Ktedz_J0, Ktedz_J1 )
  deallocate( kernelExJ0, kernelExJ1, kernelEyJ0, kernelEyJ1, kernelEzJ1 )
  deallocate( kernelHxJ0, kernelHxJ1, kernelHyJ0, kernelHyJ1, kernelHzJ1 )
end subroutine hmdy_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine hmdy_xkyz_loops( Tx, ky, h0, n, esp, condut, neta, zeta, cx, z, Ex_ky, Ey_ky, Ez_ky, Hx_ky, Hy_ky, Hz_ky )
  implicit none
  integer, intent(in) :: n
  real(dp), intent(in) :: Tx, ky, h0, esp(:), condut(1 : n), cx, z
  complex(dp), intent(in) :: zeta, neta
  complex(dp), intent(out) :: Ex_ky, Ey_ky, Ez_ky, Hx_ky, Hy_ky, Hz_ky

  integer :: i, j, camad, camadT, autor, filtro, npts, nptc, funs, func
  real(dp) :: x
  real(dp), dimension(:), allocatable :: h, kxsen, kxcos, kr2sen, kr2cos, w_sen, w_cos, prof

  complex(dp), dimension(:), allocatable :: wvnb2, AMdwSen, AMdwCos, AMupSen, AMupCos, FEdwSen, FEdwCos, FEupSen, FEupCos
  complex(dp), dimension(:,:), allocatable :: uSen, uCos, AdmIntSen, AdmIntCos, ImpIntSen, ImpIntCos
  complex(dp), dimension(:,:), allocatable :: uhSen, uhCos, tghSen, tghCos
  complex(dp), dimension(:,:), allocatable :: AdmApdwSen, AdmApdwCos, ImpApdwSen, ImpApdwCos
  complex(dp), dimension(:,:), allocatable :: RTEdwSen, RTEdwCos, RTMdwSen, RTMdwCos
  complex(dp), dimension(:,:), allocatable :: AdmApupSen, AdmApupCos, ImpApupSen, ImpApupCos
  complex(dp), dimension(:,:), allocatable :: RTEupSen, RTEupCos, RTMupSen, RTMupCos
  complex(dp), dimension(:,:), allocatable :: TMdwSen, TMdwCos, TEdwSen, TEdwCos, TMupSen, TMupCos, TEupSen, TEupCos
  complex(dp), dimension(:), allocatable :: Ktmdz_Sen, Ktmdz_Cos, Ktm_Sen, Ktm_Cos, Kte_Sen, Kte_Cos, Ktedz_Sen, Ktedz_Cos
  complex(dp), dimension(:), allocatable :: kernelEx, kernelEy, kernelEz, kernelHx, kernelHy, kernelHz

  if ( dabs(cx - Tx) < eps .and. dabs(Tx) > eps ) then
    x = dsign( 1.d-1, Tx )
  elseif ( dabs(cx - Tx) < eps .and. dabs(Tx) < eps ) then
    x = 1.d-1
  else
    x = cx - Tx
  end if

  call sanitizedata(n, h0, z, esp, camadT, camad, h, prof)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! to workaround the warning: ... may be used uninitialized in this function
  allocate( TMupSen(1,1), TMupCos(1,1) )
  allocate( TEupSen(1,1), TEupCos(1,1) )
  allocate( TMdwSen(1,1), TMdwCos(1,1) )
  allocate( TEdwSen(1,1), TEdwCos(1,1) )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  filtro = 1  !Designa o tipo de filtro usado na subrotina de pesos e abscisas de vários filtros.
        !O algarismo 0 é usado para J0 e J1, enquanto 1 é para seno e cosseno.
  autor = 4   !Designa o criador do filtro. No caso de seno ou cosseno os que possuo são os do Frayzer (1) e do Kerry Key (4)
  funs = 2    !Designa se é filtro seno (2) ou filtro cosseno (3)
  func = 3    !Designa se é filtro seno (2) ou filtro cosseno (3)
  npts = 241  !Designa o número de pontos usado no filtro seno.
  nptc = 241    !Designa o número de pontos usado no filtro cosseno.

  allocate( kxsen(npts), kxcos(nptc), w_sen(npts), w_cos(nptc) )

  call constfiltro( filtro, autor, funs, npts, x, Kxsen, w_sen )
  call constfiltro( filtro, autor, func, nptc, x, Kxcos, w_cos )

  allocate( wvnb2(0 : n), uSen(npts,0 : n), uCos(nptc,0 : n), AdmIntSen(npts,0 : n), AdmIntCos(nptc,0 : n) )
  allocate( ImpIntCos(nptc,0 : n), ImpIntSen(npts,0 : n) )
  allocate( uhSen(npts,0 : n), uhCos(nptc,0 : n), tghSen(npts,0 : n), tghCos(nptc,0 : n), kr2sen(npts), kr2cos(nptc) )

  kr2sen = kxsen * kxsen + ky * ky
  kr2cos = kxcos * kxcos + ky * ky

  wvnb2(0) = -zeta * neta
  uSen(:,0) = sqrt( kr2sen - wvnb2(0) )
  uCos(:,0) = sqrt( kr2cos - wvnb2(0) )
  AdmIntSen(:,0) = uSen(:,0) / zeta
  AdmIntCos(:,0) = uCos(:,0) / zeta
  ImpIntSen(:,0) = uSen(:,0) / neta
  ImpIntCos(:,0) = uCos(:,0) / neta
  uhSen(:,0) = uSen(:,0) * h(0)
  uhCos(:,0) = uCos(:,0) * h(0)
  tghSen(:,0) = ( 1.d0 - exp( -2.d0 * uhSen(:,0) ) ) / ( 1.d0 + exp( -2.d0 * uhSen(:,0) ) )
  tghCos(:,0) = ( 1.d0 - exp( -2.d0 * uhCos(:,0) ) ) / ( 1.d0 + exp( -2.d0 * uhCos(:,0) ) )
  do i = 1, n
    wvnb2(i) = -zeta * condut(i)
    uSen(:,i) = sqrt( kr2Sen - wvnb2(i) )
    uCos(:,i) = sqrt( kr2Cos - wvnb2(i) )
    AdmIntSen(:,i) = uSen(:,i) / zeta
    AdmIntCos(:,i) = uCos(:,i) / zeta
    ImpIntSen(:,i) = uSen(:,i) / condut(i)
    ImpIntCos(:,i) = uCos(:,i) / condut(i)
    uhSen(:,i) = uSen(:,i) * h(i)
    uhCos(:,i) = uCos(:,i) * h(i)
    tghSen(:,i) = ( 1.d0 - exp( -2.d0 * uhSen(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhSen(:,i) ) )
    tghCos(:,i) = ( 1.d0 - exp( -2.d0 * uhCos(:,i) ) ) / ( 1.d0 + exp( -2.d0 * uhCos(:,i) ) )
  end do

  allocate( AdmApdwSen(npts,1 : n), AdmApdwCos(nptc,1 : n), ImpApdwSen(npts,1 : n), ImpApdwCos(nptc,1 : n) )
  allocate( RTEdwSen(npts,0 : n), RTEdwCos(nptc,0 : n), RTMdwSen(npts,0 : n), RTMdwCos(nptc,0 : n) )

  AdmApdwSen(:,n) = AdmIntSen(:,n)
  AdmApdwCos(:,n) = AdmIntCos(:,n)
  ImpApdwSen(:,n) = ImpIntSen(:,n)
  ImpApdwCos(:,n) = ImpIntCos(:,n)
  RTEdwSen(:,n) = (0.d0,0.d0)
  RTEdwCos(:,n) = (0.d0,0.d0)
  RTMdwSen(:,n) = (0.d0,0.d0)
  RTMdwCos(:,n) = (0.d0,0.d0)
  do i = n-1, 1, -1
    AdmApdwSen(:,i) = AdmIntSen(:,i) * ( AdmApdwSen(:,i + 1) + AdmIntSen(:,i) * &
                     tghSen(:,i) ) / ( AdmIntSen(:,i) + AdmApdwSen(:,i + 1) * tghSen(:,i) )
    AdmApdwCos(:,i) = AdmIntCos(:,i) * ( AdmApdwCos(:,i + 1) + AdmIntCos(:,i) * &
                     tghCos(:,i) ) / ( AdmIntCos(:,i) + AdmApdwCos(:,i + 1) * tghCos(:,i) )
    ImpApdwSen(:,i) = ImpIntSen(:,i) * ( ImpApdwSen(:,i + 1) + ImpIntSen(:,i) * &
                     tghSen(:,i) ) / ( ImpIntSen(:,i) + ImpApdwSen(:,i + 1) * tghSen(:,i) )
    ImpApdwCos(:,i) = ImpIntCos(:,i) * ( ImpApdwCos(:,i + 1) + ImpIntCos(:,i) * &
                     tghCos(:,i) ) / ( ImpIntCos(:,i) + ImpApdwCos(:,i + 1) * tghCos(:,i) )
    RTEdwSen(:,i)=(AdmIntSen(:,i)-AdmApdwSen(:,i+1))/(AdmIntSen(:,i)+AdmApdwSen(:,i+1))
    RTEdwCos(:,i)=(AdmIntCos(:,i)-AdmApdwCos(:,i+1))/(AdmIntCos(:,i)+AdmApdwCos(:,i+1))
    RTMdwSen(:,i)=(ImpIntSen(:,i)-ImpApdwSen(:,i+1))/(ImpIntSen(:,i)+ImpApdwSen(:,i+1))
    RTMdwCos(:,i)=(ImpIntCos(:,i)-ImpApdwCos(:,i+1))/(ImpIntCos(:,i)+ImpApdwCos(:,i+1))
  end do
  RTEdwSen(:,0) = ( AdmIntSen(:,0) - AdmApdwSen(:,1) ) / ( AdmIntSen(:,0) + AdmApdwSen(:,1) )
  RTEdwCos(:,0) = ( AdmIntCos(:,0) - AdmApdwCos(:,1) ) / ( AdmIntCos(:,0) + AdmApdwCos(:,1) )
  RTMdwSen(:,0) = ( ImpIntSen(:,0) - ImpApdwSen(:,1) ) / ( ImpIntSen(:,0) + ImpApdwSen(:,1) )
  RTMdwCos(:,0) = ( ImpIntCos(:,0) - ImpApdwCos(:,1) ) / ( ImpIntCos(:,0) + ImpApdwCos(:,1) )

  allocate( AdmApupSen(npts,0 : n - 1), AdmApupCos(nptc,0 : n - 1), ImpApupSen(npts,0 : n - 1), ImpApupCos(nptc,0 : n - 1) )
  allocate( RTEupSen(npts,0 : n), RTEupCos(nptc,0 : n), RTMupSen(npts,0 : n), RTMupCos(nptc,0 : n) )

  AdmApupSen(:,0) = AdmIntSen(:,0)
  AdmApupCos(:,0) = AdmIntCos(:,0)
  ImpApupSen(:,0) = ImpIntSen(:,0)
  ImpApupCos(:,0) = ImpIntCos(:,0)
  RTEupSen(:,0) = (0.d0,0.d0)
  RTEupCos(:,0) = (0.d0,0.d0)
  RTMupSen(:,0) = (0.d0,0.d0)
  RTMupCos(:,0) = (0.d0,0.d0)
  do i = 1, n - 1
    AdmApupSen(:,i) = AdmIntSen(:,i) * ( AdmApupSen(:,i - 1) + AdmIntSen(:,i) * &
                      tghSen(:,i) ) / ( AdmIntSen(:,i) + AdmApupSen(:,i - 1) * tghSen(:,i) )
    AdmApupCos(:,i) = AdmIntCos(:,i) * ( AdmApupCos(:,i - 1) + AdmIntCos(:,i) * &
                      tghCos(:,i) ) / ( AdmIntCos(:,i) + AdmApupCos(:,i - 1) * tghCos(:,i) )
    ImpApupSen(:,i) = ImpIntSen(:,i) * ( ImpApupSen(:,i - 1) + ImpIntSen(:,i) * &
                      tghSen(:,i) ) / ( ImpIntSen(:,i) + ImpApupSen(:,i - 1) * tghSen(:,i) )
    ImpApupCos(:,i) = ImpIntCos(:,i) * ( ImpApupCos(:,i - 1) + ImpIntCos(:,i) * &
                      tghCos(:,i) ) / ( ImpIntCos(:,i) + ImpApupCos(:,i - 1) * tghCos(:,i) )
    RTEupSen(:,i) = ( AdmIntSen(:,i) - AdmApupSen(:,i - 1) ) / ( AdmIntSen(:,i) + AdmApupSen(:,i - 1) )
    RTEupCos(:,i) = ( AdmIntCos(:,i) - AdmApupCos(:,i - 1) ) / ( AdmIntCos(:,i) + AdmApupCos(:,i - 1) )
    RTMupSen(:,i) = ( ImpIntSen(:,i) - ImpApupSen(:,i - 1) ) / ( ImpIntSen(:,i) + ImpApupSen(:,i - 1) )
    RTMupCos(:,i) = ( ImpIntCos(:,i) - ImpApupCos(:,i - 1) ) / ( ImpIntCos(:,i) + ImpApupCos(:,i - 1) )
  end do
  RTEupSen(:,n)=(AdmIntSen(:,n)-AdmApupSen(:,n-1))/(AdmIntSen(:,n)+AdmApupSen(:,n-1))
  RTEupCos(:,n)=(AdmIntCos(:,n)-AdmApupCos(:,n-1))/(AdmIntCos(:,n)+AdmApupCos(:,n-1))
  RTMupSen(:,n)=(ImpIntSen(:,n)-ImpApupSen(:,n-1))/(ImpIntSen(:,n)+ImpApupSen(:,n-1))
  RTMupCos(:,n)=(ImpIntCos(:,n)-ImpApupCos(:,n-1))/(ImpIntCos(:,n)+ImpApupCos(:,n-1))

  allocate( AMdwSen(npts), AMdwCos(nptc), AMupSen(npts), AMupCos(nptc) )
  allocate( FEdwSen(npts), FEdwCos(nptc), FEupSen(npts), FEupCos(nptc) )

  AMdwSen = ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) + RTMupSen(:,camadT) * &
            exp( uSen(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupSen(:,camadT) * RTMdwSen(:,camadT) * exp( -2.d0 * uhSen(:,camadT) ) )
  AMdwCos = ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) + RTMupCos(:,camadT) * &
            exp( uCos(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupCos(:,camadT) * RTMdwCos(:,camadT) * exp( -2.d0 * uhCos(:,camadT) ) )

  AMupSen = ( exp( uSen(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTMdwSen(:,camadT) * &
            exp( -uSen(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupSen(:,camadT) * RTMdwSen(:,camadT) * exp( -2.d0 * uhSen(:,camadT) ) )
  AMupCos = ( exp( uCos(:,camadT) * ( prof(camadT - 1) - h0 ) ) + RTMdwCos(:,camadT) * &
            exp( -uCos(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTMupCos(:,camadT) * RTMdwCos(:,camadT) * exp( -2.d0 * uhCos(:,camadT) ) )

  FEdwSen = ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) - RTEupSen(:,camadT) * &
            exp( uSen(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupSen(:,camadT) * RTEdwSen(:,camadT) * exp( -2.d0 * uhSen(:,camadT) ) )
  FEdwCos = ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) - RTEupCos(:,camadT) * &
            exp( uCos(:,camadT) * ( prof(camadT - 1) - h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupCos(:,camadT) * RTEdwCos(:,camadT) * exp( -2.d0 * uhCos(:,camadT) ) )

  FEupSen = ( exp( uSen(:,camadT) * ( prof(camadT - 1) - h0 ) ) - RTEdwSen(:,camadT) * &
            exp( -uSen(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupSen(:,camadT) * RTEdwSen(:,camadT) * exp( -2.d0 * uhSen(:,camadT) ) )
  FEupCos = ( exp( uCos(:,camadT) * ( prof(camadT - 1) - h0 ) ) - RTEdwCos(:,camadT) * &
            exp( -uCos(:,camadT) * ( prof(camadT) + h(camadT) - h0 ) ) ) / &
            ( 1.d0 - RTEupCos(:,camadT) * RTEdwCos(:,camadT) * exp( -2.d0 * uhCos(:,camadT) ) )

  if ( camad > camadT ) then
    deallocate( TMdwSen, TMdwCos, TEdwSen, TEdwCos )
    allocate( TMdwSen(npts,camadT : camad), TMdwCos(nptc,camadT : camad) )
    allocate( TEdwSen(npts,camadT : camad), TEdwCos(nptc,camadT : camad) )
    do j = camadT, camad
      if ( j == camadT ) then
        TMdwSen(:,j) = zeta * my / ( 2.d0 * ImpIntSen(:,camadT) )
        TMdwCos(:,j) = zeta * my / ( 2.d0 * ImpIntCos(:,camadT) )
        TEdwSen(:,j) = - zeta * my / 2.d0
        TEdwCos(:,j) = - zeta * my / 2.d0
      else if ( j == (camadT + 1) .and. j == n ) then
        TMdwSen(:,j) = TMdwSen(:,j - 1) * ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupSen(:,camadT) * AMupSen(:) * exp( -uhSen(:,camadT) ) + RTMdwSen(:,camadT) * AMdwSen(:) )
        TMdwCos(:,j) = TMdwCos(:,j - 1) * ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupCos(:,camadT) * AMupCos(:) * exp( -uhCos(:,camadT) ) + RTMdwCos(:,camadT) * AMdwCos(:) )
        TEdwSen(:,j) = TEdwSen(:,j - 1) * ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) - &
                      RTEupSen(:,camadT) * FEupSen(:) * exp( -uhSen(:,camadT) ) + RTEdwSen(:,camadT) * FEdwSen(:) )
        TEdwCos(:,j) = TEdwCos(:,j - 1) * ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) - &
                      RTEupCos(:,camadT) * FEupCos(:) * exp( -uhCos(:,camadT) ) + RTEdwCos(:,camadT) * FEdwCos(:) )
      else if ( j == (camadT + 1) .and. j /= n ) then
        TMdwSen(:,j) = TMdwSen(:,j - 1) * ( exp( -uSen(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupSen(:,camadT) * AMupSen(:) * exp( -uhSen(:,camadT) ) + &
                      RTMdwSen(:,camadT) * AMdwSen(:) ) / ( 1.d0 + RTMdwSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TMdwCos(:,j) = TMdwCos(:,j - 1) * ( exp( -uCos(:,camadT) * ( prof(camadT) - h0 ) ) + &
                      RTMupCos(:,camadT) * AMupCos(:) * exp( -uhCos(:,camadT) ) + &
                      RTMdwCos(:,camadT) * AMdwCos(:) ) / ( 1.d0 + RTMdwCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )

        TEdwSen(:,j)=TEdwSen(:,j-1)*(exp(-uSen(:,camadT)*(prof(camadT)-h0)) - &
                      RTEupSen(:,camadT)*FEupSen(:)*exp(-uhSen(:,camadT)) + &
                      RTEdwSen(:,camadT)*FEdwSen(:)) / (1.d0 + RTEdwSen(:,j)*exp(-2.d0*uhSen(:,j)))
        TEdwCos(:,j)=TEdwCos(:,j-1)*(exp(-uCos(:,camadT)*(prof(camadT)-h0)) - &
                      RTEupCos(:,camadT)*FEupCos(:)*exp(-uhCos(:,camadT)) + &
                      RTEdwCos(:,camadT)*FEdwCos(:)) / (1.d0 + RTEdwCos(:,j)*exp(-2.d0*uhCos(:,j)))
      else if ( j /= n ) then
        TMdwSen(:,j) = TMdwSen(:,j - 1) * ( 1.d0 + RTMdwSen(:,j - 1) ) * exp( -uhSen(:,j - 1) ) / &
                      ( 1.d0 + RTMdwSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TMdwCos(:,j) = TMdwCos(:,j - 1) * ( 1.d0 + RTMdwCos(:,j - 1) ) * exp( -uhCos(:,j - 1) ) / &
                      ( 1.d0 + RTMdwCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )
        TEdwSen(:,j) = TEdwSen(:,j - 1) * ( 1.d0 + RTEdwSen(:,j - 1) ) * exp( -uhSen(:,j - 1) ) / &
                      ( 1.d0 + RTEdwSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TEdwCos(:,j) = TEdwCos(:,j - 1) * ( 1.d0 + RTEdwCos(:,j - 1) ) * exp( -uhCos(:,j - 1) ) / &
                      ( 1.d0 + RTEdwCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )
      else if ( j == n ) then
        TMdwSen(:,j) = TMdwSen(:,j - 1) * ( 1.d0 + RTMdwSen(:,j - 1) ) * exp( -uhSen(:,j - 1) )
        TMdwCos(:,j) = TMdwCos(:,j - 1) * ( 1.d0 + RTMdwCos(:,j - 1) ) * exp( -uhCos(:,j - 1) )
        TEdwSen(:,j) = TEdwSen(:,j - 1) * ( 1.d0 + RTEdwSen(:,j - 1) ) * exp( -uhSen(:,j - 1) )
        TEdwCos(:,j) = TEdwCos(:,j - 1) * ( 1.d0 + RTEdwCos(:,j - 1) ) * exp( -uhCos(:,j - 1) )
      end if
    end do
  else if ( camad < camadT ) then
    deallocate( TMupSen, TMupCos,TEupSen, TEupCos )
    allocate( TMupSen(npts,camad : camadT), TMupCos(nptc,camad : camadT) )
    allocate( TEupSen(npts,camad : camadT), TEupCos(nptc,camad : camadT) )
    do j = camadT, camad, -1
      if ( j == camadT ) then
        TMupSen(:,j) = zeta * my / ( 2.d0 * ImpIntSen(:,camadT) )
        TMupCos(:,j) = zeta * my / ( 2.d0 * ImpIntCos(:,camadT) )
        TEupSen(:,j) = zeta * my / 2.d0
        TEupCos(:,j) = zeta * my / 2.d0
      else if ( j == (camadT - 1) .and. j == 0 ) then
        TMupSen(:,j) = TMupSen(:,j + 1) * ( exp( -uSen(:,camadT) * h0 ) + &
          RTMupSen(:,camadT) * AMupSen(:) + RTMdwSen(:,camadT) * AMdwSen(:) * exp( -uhSen(:,camadT) ) )
          TMupCos(:,j) = TMupCos(:,j + 1) * ( exp( -uCos(:,camadT) * h0 ) + &
          RTMupCos(:,camadT) * AMupCos(:) + RTMdwCos(:,camadT) * AMdwCos(:) * exp( -uhCos(:,camadT) ) )
        TEupSen(:,j) = TEupSen(:,j + 1) * ( exp( -uSen(:,camadT) * h0 ) + &
          RTEupSen(:,camadT) * FEupSen(:) - RTEdwSen(:,camadT) * FEdwSen(:) * exp( -uhSen(:,camadT) ) )
        TEupCos(:,j) = TEupCos(:,j + 1) * ( exp( -uCos(:,camadT) * h0 ) + &
          RTEupCos(:,camadT) * FEupCos(:) - RTEdwCos(:,camadT) * FEdwCos(:) * exp( -uhCos(:,camadT) ) )
      else if ( j == (camadT - 1) .and. j /= 0 ) then
        TMupSen(:,j) = TMupSen(:,j + 1) * ( exp( uSen(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
          RTMupSen(:,camadT) * AMupSen(:) + RTMdwSen(:,camadT) * AMdwSen(:) * &
          exp( -uhSen(:,camadT) ) ) / ( 1.d0 + RTMupSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TMupCos(:,j) = TMupCos(:,j + 1) * ( exp( uCos(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
          RTMupCos(:,camadT) * AMupCos(:) + RTMdwCos(:,camadT) * AMdwCos(:) * &
          exp( -uhCos(:,camadT) ) ) / ( 1.d0 + RTMupCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )

        TEupSen(:,j) = TEupSen(:,j + 1) * ( exp( uSen(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
          RTEupSen(:,camadT) * FEupSen(:) - RTEdwSen(:,camadT) * FEdwSen(:) * &
          exp( -uhSen(:,camadT) ) ) / ( 1.d0 + RTEupSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TEupCos(:,j) = TEupCos(:,j + 1) * ( exp( uCos(:,camadT) * ( prof(camadT - 1) - h0 ) ) + &
          RTEupCos(:,camadT) * FEupCos(:) - RTEdwCos(:,camadT) * FEdwCos(:) * &
          exp( -uhCos(:,camadT) ) ) / ( 1.d0 + RTEupCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )
      else if ( j /= 0 ) then
        TMupSen(:,j) = TMupSen(:,j + 1) * ( 1.d0 + RTMupSen(:,j + 1) ) * exp( -uhSen(:,j + 1) ) / &
                      ( 1.d0 + RTMupSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TMupCos(:,j) = TMupCos(:,j + 1) * ( 1.d0 + RTMupCos(:,j + 1) ) * exp( -uhCos(:,j + 1) ) / &
                      ( 1.d0 + RTMupCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )
        TEupSen(:,j) = TEupSen(:,j + 1) * ( 1.d0 + RTEupSen(:,j + 1) ) * exp( -uhSen(:,j + 1) ) / &
                      ( 1.d0 + RTEupSen(:,j) * exp( -2.d0 * uhSen(:,j) ) )
        TEupCos(:,j) = TEupCos(:,j + 1) * ( 1.d0 + RTEupCos(:,j + 1) ) * exp( -uhCos(:,j + 1) ) / &
                      ( 1.d0 + RTEupCos(:,j) * exp( -2.d0 * uhCos(:,j) ) )
      else if ( j == 0 ) then
        TMupSen(:,j) = TMupSen(:,1) * ( 1.d0 + RTMupSen(:,1) ) * exp( -uhSen(:,1) )
        TMupCos(:,j) = TMupCos(:,1) * ( 1.d0 + RTMupCos(:,1) ) * exp( -uhCos(:,1) )
        TEupSen(:,j) = TEupSen(:,1) * ( 1.d0 + RTEupSen(:,1) ) * exp( -uhSen(:,1) )
        TEupCos(:,j) = TEupCos(:,1) * ( 1.d0 + RTEupCos(:,1) ) * exp( -uhCos(:,1) )
      end if
    end do
  else
    deallocate( TMupSen, TMupCos, TEupSen, TEupCos, TMdwSen, TMdwCos, TEdwSen, TEdwCos )
    allocate( TMdwSen(npts,camadT : camad), TMdwCos(nptc,camadT : camad) )
    allocate( TEdwSen(npts,camadT : camad), TEdwCos(nptc,camadT : camad) )
    allocate( TMupSen(npts,camad : camadT), TMupCos(nptc,camad : camadT) )
    allocate( TEupSen(npts,camad : camadT), TEupCos(nptc,camad : camadT) )

    TMdwSen(:,camad) = zeta * my / (2.d0 * ImpIntSen(:,camadT))
    TMdwCos(:,camad) = zeta * my / (2.d0 * ImpIntCos(:,camadT))

    TEdwSen(:,camad) = - zeta * my / 2.d0
    TEdwCos(:,camad) = - zeta * my / 2.d0

    TMupSen(:,camad) = TMdwSen(:,camad)
    TMupCos(:,camad) = TMdwCos(:,camad)

    TEupSen(:,camad) = - TEdwSen(:,camad)
    TEupCos(:,camad) = - TEdwCos(:,camad)
  end if

  allocate(Ktmdz_Sen(npts),Ktmdz_Cos(nptc),Ktm_Sen(npts),Ktm_Cos(nptc),Kte_Sen(npts),Kte_Cos(nptc),Ktedz_Sen(npts),Ktedz_Cos(nptc))
  allocate(kernelEx(nptc),kernelEy(npts),kernelEz(npts))
  allocate(kernelHx(npts),kernelHy(nptc),kernelHz(nptc))

  if ( camad == 0 .and. camadT /= 0 ) then
    Ktmdz_Sen = ImpIntSen(:,0) * TMupSen(:,0) * exp( uSen(:,0) * z )
    Ktmdz_Cos = ImpIntCos(:,0) * TMupCos(:,0) * exp( uCos(:,0) * z )
    Kte_Sen = TEupSen(:,0) * exp( uSen(:,0) * z )
    Kte_Cos = TEupCos(:,0) * exp( uCos(:,0) * z )

    Ktm_Sen = TMupSen(:,0) * exp( uSen(:,0) * z )
    Ktm_Cos = TMupCos(:,0) * exp( uCos(:,0) * z )
    Ktedz_Sen = AdmIntSen(:,0) * TEupSen(:,0) * exp( uSen(:,0) * z )
    Ktedz_Cos = AdmIntCos(:,0) * TEupCos(:,0) * exp( uCos(:,0) * z )

    kernelEx = ( kxsen * ky * ( Ktmdz_Sen - Kte_Sen ) / kr2sen ) * w_sen
    Ey_ky = (0.d0,1.d0) * sum( kernelEx ) / ( pi * x )

    kernelEy = ( ( kxcos * kxcos * Ktmdz_Cos + ky * ky * Kte_Cos ) / kr2cos ) * w_cos
    Ex_ky = sum( kernelEy ) / ( pi * dabs(x) )

    kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen  * w_sen
    Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * neta )

    kernelHx = ( ( kxcos * kxcos * Ktm_Cos + ky * ky * Ktedz_Cos ) / kr2cos ) * w_cos
    Hy_ky = -sum( kernelHx ) / ( pi * dabs(x) )

    kernelHy = ( kxsen * ky * ( Ktm_Sen - Ktedz_Sen ) / kr2sen ) * w_sen
    Hx_ky = (0.d0,1.d0) * sum( kernelHy ) / ( pi * x )

    kernelHz = ( cmplx(0.d0,ky,kind=dp) / zeta * Kte_Cos ) * w_cos
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )
  else if ( camad < camadT ) then !camada k
    Ktmdz_Sen = ImpIntSen(:,camad) * TMupSen(:,camad) * ( exp( uSen(:,camad) * ( z - prof(camad) ) ) - &
                    RTMupSen(:,camad) * exp( -uSen(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    ktmdz_Cos = ImpIntCos(:,camad) * TMupCos(:,camad) * ( exp( uCos(:,camad) * ( z - prof(camad) ) ) - &
                    RTMupCos(:,camad) * exp( -uCos(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    Kte_Sen = TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - prof(camad) ) ) + &
                    RTEupSen(:,camad) * exp( -uSen(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    Kte_Cos = TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - prof(camad) ) ) + &
                    RTEupCos(:,camad) * exp( -uCos(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )

    ktm_Sen = TMupSen(:,camad) * ( exp( uSen(:,camad) * ( z - prof(camad) ) ) + &
                    RTMupSen(:,camad) * exp( -uSen(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    ktm_Cos = TMupCos(:,camad) * ( exp( uCos(:,camad) * ( z - prof(camad) ) ) + &
                    RTMupCos(:,camad) * exp( -uCos(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    Ktedz_Sen = AdmIntSen(:,camad) * TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - prof(camad) ) ) - &
                    RTEupSen(:,camad) * exp( -uSen(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )
    Ktedz_Cos = AdmIntCos(:,camad) * TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - prof(camad) ) ) - &
                    RTEupCos(:,camad) * exp( -uCos(:,camad) * ( z - prof(camad - 1) + h(camad) ) ) )

    kernelEx = ( kxsen * ky * ( Ktmdz_Sen - Kte_Sen ) / kr2sen ) * w_sen
    Ey_ky = (0.d0,1.d0) * sum( kernelEx ) / ( pi * x )

    kernelEy = ( ( kxcos * kxcos * Ktmdz_Cos + ky * ky * Kte_Cos ) / kr2cos ) * w_cos
    Ex_ky = sum( kernelEy ) / ( pi * dabs(x) )

    kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen  * w_sen
    Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * condut(camad) )

    kernelHx = ( ( kxcos * kxcos * Ktm_Cos + ky * ky * Ktedz_Cos ) / kr2cos ) * w_cos
    Hy_ky = -sum( kernelHx ) / ( pi * dabs(x) )

    kernelHy = ( kxsen * ky * ( Ktm_Sen - Ktedz_Sen ) / kr2sen ) * w_sen
    Hx_ky = (0.d0,1.d0) * sum( kernelHy ) / ( pi * x )

    kernelHz = ( cmplx(0.d0,ky,kind=dp) / zeta * Kte_Cos ) * w_cos
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )
  else if ( camad == camadT .and. z <= h0 ) then  !na mesma camada do transmissor mas acima dele
    Ktmdz_Sen = ImpIntSen(:,camad) * TMupSen(:,camad) * ( exp( uSen(:,camad) * ( z - h0 ) ) - &
                    RTMupSen(:,camad) * AMupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                    RTMdwSen(:,camad) * AMdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktmdz_Cos = ImpIntCos(:,camad) * TMupCos(:,camad) * ( exp( uCos(:,camad) * ( z - h0 ) ) - &
                    RTMupCos(:,camad) * AMupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                    RTMdwCos(:,camad) * AMdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )
    Kte_Sen = TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - h0 ) ) + &
                    RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
                    RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Kte_Cos = TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - h0 ) ) + &
                    RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
                    RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    Ktm_Sen = TMupSen(:,camad) * ( exp( uSen(:,camad) * ( z - h0 ) ) + &
              RTMupSen(:,camad) * AMupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwSen(:,camad) * AMdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktm_Cos = TMupCos(:,camad) * ( exp( uCos(:,camad) * ( z - h0 ) ) + &
              RTMupCos(:,camad) * AMupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
              RTMdwCos(:,camad) * AMdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )
    Ktedz_Sen = AdmIntSen(:,camad) * TEupSen(:,camad) * ( exp( uSen(:,camad) * ( z - h0 ) ) - &
              RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktedz_Cos = AdmIntCos(:,camad) * TEupCos(:,camad) * ( exp( uCos(:,camad) * ( z - h0 ) ) - &
              RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
              RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    kernelEx = ( kxsen * ky * ( Ktmdz_Sen - Kte_Sen ) / kr2sen ) * w_sen
    Ey_ky = (0.d0,1.d0) * sum( kernelEx ) / ( pi * x )

    kernelEy = ( ( kxcos * kxcos * Ktmdz_Cos + ky * ky * Kte_Cos ) / kr2cos ) * w_cos
    Ex_ky = sum( kernelEy ) / ( pi * dabs(x) )

    if ( camad /= 0 ) then
      kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen  * w_sen
      Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * condut(camad) )
    else
      kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen  * w_sen
      Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * neta )
    end if

    kernelHx = ( ( kxcos * kxcos * Ktm_Cos + ky * ky * Ktedz_Cos ) / kr2cos ) * w_cos
    Hy_ky = -sum( kernelHx ) / ( pi * dabs(x) )

    kernelHy = ( kxsen * ky * ( Ktm_Sen - Ktedz_Sen ) / kr2sen ) * w_sen
    Hx_ky = (0.d0,1.d0) * sum( kernelHy ) / ( pi * x )

    kernelHz = ( cmplx(0.d0,ky,kind=dp) / zeta * Kte_Cos ) * w_cos
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )
  else if ( camad == camadT .and. z > h0 ) then !na mesma camada do transmissor mas abaixo dele
    Ktmdz_Sen = ImpIntSen(:,camad) * TMdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - h0 ) ) + &
                RTMupSen(:,camad) * AMupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTMdwSen(:,camad) * AMdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktmdz_Cos = ImpIntCos(:,camad) * TMdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - h0 ) ) + &
                RTMupCos(:,camad) * AMupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTMdwCos(:,camad) * AMdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )
    Kte_Sen = TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - h0 ) ) - &
                RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Kte_Cos = TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - h0 ) ) - &
                RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    Ktm_Sen = TMdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - h0 ) ) + &
                RTMupSen(:,camad) * AMupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTMdwSen(:,camad) * AMdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktm_Cos = TMdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - h0 ) ) + &
                RTMupCos(:,camad) * AMupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                RTMdwCos(:,camad) * AMdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )
    Ktedz_Sen = AdmIntSen(:,camad) * TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - h0 ) ) - &
                RTEupSen(:,camad) * FEupSen(:) * exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwSen(:,camad) * FEdwSen(:) * exp( uSen(:,camad) * ( z - prof(camad) ) ) )
    Ktedz_Cos = AdmIntCos(:,camad) * TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - h0 ) ) - &
                RTEupCos(:,camad) * FEupCos(:) * exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
                RTEdwCos(:,camad) * FEdwCos(:) * exp( uCos(:,camad) * ( z - prof(camad) ) ) )

    kernelEx = ( -kxsen * ky * ( Ktmdz_Sen + Kte_Sen ) / kr2sen ) * w_sen
    Ey_ky = (0.d0,1.d0) * sum( kernelEx ) / ( pi * x )

    kernelEy = ( ( -kxcos * kxcos * Ktmdz_Cos + ky * ky * Kte_Cos ) / kr2cos ) * w_cos
    Ex_ky = sum( kernelEy ) / ( pi * dabs(x) )

    if ( camad /= 0 ) then
            kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen * w_sen
            Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * condut(camad) )
    else
            kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen * w_sen
            Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * neta )
    end if

    kernelHx =  ( ( -kxcos * kxcos * Ktm_Cos + ky * ky * Ktedz_Cos ) / kr2cos ) * w_cos
    Hy_ky = sum( kernelHx ) / ( pi * dabs(x) )

    kernelHy = ( kxsen * ky * ( Ktm_Sen + Ktedz_Sen ) / kr2sen ) * w_sen
    Hx_ky = (0.d0,1.d0) * sum( kernelHy ) / ( pi * x )

    kernelHz = ( cmplx(0.d0,ky,kind=dp) / zeta * Kte_Cos ) * w_cos
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )
  else if ( camad > camadT .and. camad /= n ) then !camada j
    Ktmdz_Sen = ImpIntSen(:,camad) * TMdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
                    RTMdwSen(:,camad) * exp( uSen(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Ktmdz_Cos = ImpIntCos(:,camad) * TMdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
                    RTMdwCos(:,camad) * exp( uCos(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Kte_Sen = TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                    RTEdwSen(:,camad) * exp( uSen(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Kte_Cos = TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                    RTEdwCos(:,camad) * exp( uCos(:,camad) * ( z - prof(camad) - h(camad) ) ) )

    Ktm_Sen = TMdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) + &
                    RTMdwSen(:,camad) * exp( uSen(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Ktm_Cos = TMdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) + &
                    RTMdwCos(:,camad) * exp( uCos(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Ktedz_Sen = AdmIntSen(:,camad) * TEdwSen(:,camad) * ( exp( -uSen(:,camad) * ( z - prof(camad - 1) ) ) - &
                    RTEdwSen(:,camad) * exp( uSen(:,camad) * ( z - prof(camad) - h(camad) ) ) )
    Ktedz_Cos = AdmIntCos(:,camad) * TEdwCos(:,camad) * ( exp( -uCos(:,camad) * ( z - prof(camad - 1) ) ) - &
                    RTEdwCos(:,camad) * exp( uCos(:,camad) * ( z - prof(camad) - h(camad) ) ) )

    kernelEx = ( -kxsen * ky * ( Ktmdz_Sen + Kte_Sen ) / kr2sen ) * w_sen
    Ey_ky = (0.d0,1.d0) * sum( kernelEx ) / ( pi * x )

    kernelEy = ( ( -kxcos * kxcos * Ktmdz_Cos + ky * ky * Kte_Cos ) / kr2cos ) * w_cos
    Ex_ky = sum( kernelEy ) / ( pi * dabs(x) )

    kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen  * w_sen
    Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * condut(camad) )

    kernelHx =  ( ( -kxcos * kxcos * Ktm_Cos + ky * ky * Ktedz_Cos ) / kr2cos ) * w_cos
    Hy_ky = sum( kernelHx ) / ( pi * dabs(x) )

    kernelHy = ( kxsen * ky * ( Ktm_Sen + Ktedz_Sen ) / kr2sen ) * w_sen
    Hx_ky = (0.d0,1.d0) * sum( kernelHy ) / ( pi * x )

    kernelHz = ( cmplx(0.d0,ky,kind=dp) / zeta * Kte_Cos ) * w_cos
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )
  else  !camada n
    Ktmdz_Sen = ImpIntSen(:,n) * TMdwSen(:,n) * exp( -uSen(:,n) * ( z - prof(n - 1) ) )
    Ktmdz_Cos = ImpIntCos(:,n) * TMdwCos(:,n) * exp( -uCos(:,n) * ( z - prof(n - 1) ) )
    Kte_Sen = TEdwSen(:,n) * exp( -uSen(:,n) * ( z - prof(n - 1) ) )
    Kte_Cos = TEdwCos(:,n) * exp( -uCos(:,n) * ( z - prof(n - 1) ) )

    Ktm_Sen = TMdwSen(:,n) * exp( -uSen(:,n) * ( z - prof(n - 1) ) )
    Ktm_Cos = TMdwCos(:,n) * exp( -uCos(:,n) * ( z - prof(n - 1) ) )
    Ktedz_Sen = AdmIntSen(:,n) * TEdwSen(:,n) * exp( -uSen(:,n) * ( z - prof(n - 1) ) )
    Ktedz_Cos = AdmIntCos(:,n) * TEdwCos(:,n) * exp( -uCos(:,n) * ( z - prof(n - 1) ) )

    kernelEx = ( -kxsen * ky * ( Ktmdz_Sen + Kte_Sen ) / kr2sen ) * w_sen
    Ey_ky = (0.d0,1.d0) * sum( kernelEx ) / ( pi * x )

    kernelEy = ( ( -kxcos * kxcos * Ktmdz_Cos + ky * ky * Kte_Cos ) / kr2cos ) * w_cos
    Ex_ky = sum( kernelEy ) / ( pi * dabs(x) )

    kernelEz = cmplx(0.d0,-kxsen,kind=dp) * Ktm_Sen  * w_sen
    Ez_ky = (0.d0,1.d0) * sum( kernelEz ) / ( pi * x * condut(camad) )

    kernelHx =  ( ( -kxcos * kxcos * Ktm_Cos + ky * ky * Ktedz_Cos ) / kr2cos ) * w_cos
    Hy_ky = sum( kernelHx ) / ( pi * dabs(x) )

    kernelHy = ( kxsen * ky * ( Ktm_Sen + Ktedz_Sen ) / kr2sen ) * w_sen
    Hx_ky = (0.d0,1.d0) * sum( kernelHy ) / ( pi * x )

    kernelHz = ( cmplx(0.d0,ky,kind=dp) / zeta * Kte_Cos ) * w_cos
    Hz_ky = sum( kernelHz ) / ( pi * dabs(x) )
  end if

  deallocate( h, KxSen, KxCos, w_Sen, w_Cos )
  deallocate( wvnb2, uSen, uCos, AdmIntSen, AdmIntCos, ImpIntSen, ImpIntCos, uhSen, uhCos, tghSen, tghCos )
  deallocate( AdmApdwSen, AdmApdwCos, ImpApdwSen, ImpApdwCos, RTEdwSen, RTEdwCos, RTMdwSen, RTMdwCos )
  deallocate( AdmApupSen, AdmApupCos, ImpApupSen, ImpApupCos, RTEupSen, RTEupCos, RTMupSen, RTMupCos )
  deallocate( AMdwSen, AMdwCos, AMupSen, AMupCos, FEdwSen, FEdwCos, FEupSen, FEupCos )
  deallocate( Ktmdz_Sen, Ktmdz_Cos, Ktm_Sen, Ktm_Cos, Kte_Sen, Kte_Cos, Ktedz_Sen, Ktedz_Cos )
  deallocate( kernelEx, kernelEy, kernelEz, kernelHx, kernelHy, kernelHz )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
end subroutine hmdy_xkyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end module hmdy
