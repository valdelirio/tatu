module ved
use parameters
use escolha_do_filtro
contains
  subroutine ved_xyz_loops(Tx,Ty,h0,n,esp,condut,neta,zeta,cx,cy,z,Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p)
  implicit none
  integer,intent(in)::n
  real(dp),intent(in)::Tx,Ty,h0,esp(:),condut(1:n),cx,cy,z !,neta
  complex(dp),intent(in)::zeta,neta
  complex(dp),intent(out)::Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p

  integer::i,j,k,camad,camadT,filtro,idtfcd_cJ0,ident_fJ0,nJ0,idtfcd_cJ1,ident_fJ1,nJ1
  real(dp)::x,y,r
  real(dp),dimension(:),allocatable::h,krJ0,krJ1,w_J0,w_J1,prof

! Para uso de loops:
  complex(dp),dimension(:),allocatable::wvnb2,AMdwJ0,AMdwJ1,AMupJ0,AMupJ1
  complex(dp),dimension(:,:),allocatable::uJ0,uJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1
  complex(dp),dimension(:,:),allocatable::ImpApdwJ0,ImpApdwJ1,ImpApupJ0,ImpApupJ1
  complex(dp),dimension(:,:),allocatable::RTMdwJ0,RTMdwJ1,RTMupJ0,RTMupJ1
  complex(dp),dimension(:,:),allocatable::TMdwJ0,TMdwJ1,TMupJ0,TMupJ1

  complex(dp),dimension(:),allocatable::Ktmdz_J1,Ktm_J0,Ktm_J1
  complex(dp),dimension(:),allocatable::kernelExJ1,kernelEyJ1,kernelEzJ0
  complex(dp),dimension(:),allocatable::kernelHxJ1,kernelHyJ1
  real(dp), parameter :: eps = 1.d-7
  Hz_p = (0.d0,0.d0)

  if ( dabs(cx-Tx) < eps ) then
    x=1.d-1
  else
    x=cx-Tx
  end if
  y=cy-Ty
  r = dsqrt(x**2 + y**2)

  allocate(h(0:n),prof(-1:n))
  if (size(esp)==n) then
    h(0)=0.d0
    h(1:n)=esp
  else
    h(0)=0.d0
    h(1:n-1)=esp
    h(n)=1.d300
  end if
! criando um novo vetor de profundidades que se adeque à qualquer situação patológica
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

!para descobrir em que camada está a observação
  camad = 0
  if (z <= 0.d0) then
    camad=0
  else if (z > prof(n-1)) then
    camad=n
  else
  do i=n-1,1,-1
  if (z > prof(i-1)) then
    camad=i
    exit
  end if
  end do
  end if

!para descobrir em que camada está o transmissor
  camadT = 0
  if (h0 <= 0.d0) then
    camadT = 0
  else if (h0 > prof(n-1)) then
    camadT = n
  else
  do j=n-1,1,-1
  if (h0 > prof(j-1)) then
    camadT = j
    exit
  end if
  end do
  end if
!
! To workaround the warning: ... may be used uninitialized in this function
allocate( TMdwJ0(1,1), TMdwJ1(1,1), TMupJ0(1,1), TMupJ1(1,1) )
!
!!  write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
!!  read(*,*)idtfcd_cJ0
  filtro=0  !esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
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

  allocate(KrJ0(nJ0),KrJ1(nJ1),w_J0(nJ0),w_J1(nJ1))

  call constfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0,r,KrJ0,w_J0)
  call constfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1,r,KrJ1,w_J1)

  allocate(wvnb2(0:n),uJ0(nJ0,0:n),uJ1(nJ1,0:n),ImpIntJ0(nJ0,0:n),ImpIntJ1(nJ1,0:n))
  allocate(uhJ0(nJ0,0:n),uhJ1(nJ1,0:n),tghJ0(nJ0,0:n),tghJ1(nJ1,0:n))
!
  do i=0,n
    if (i==0) then
    wvnb2(i)=-zeta*neta
    uJ0(:,i)=sqrt(krJ0*krJ0-wvnb2(i))
    uJ1(:,i)=sqrt(krJ1*krJ1-wvnb2(i))
    ImpIntJ0(:,i)=uJ0(:,i)/neta
    ImpIntJ1(:,i)=uJ1(:,i)/neta
    uhJ0(:,i)=uJ0(:,i)*h(i)
    uhJ1(:,i)=uJ1(:,i)*h(i)
    tghJ0(:,i)=(1.d0-exp(-2.d0*uhJ0(:,i)))/(1.d0+exp(-2.d0*uhJ0(:,i)))
    tghJ1(:,i)=(1.d0-exp(-2.d0*uhJ1(:,i)))/(1.d0+exp(-2.d0*uhJ1(:,i)))
    else
    wvnb2(i) = -zeta*condut(i)
    uJ0(:,i)=sqrt(krJ0*krJ0-wvnb2(i))
    uJ1(:,i)=sqrt(krJ1*krJ1-wvnb2(i))
    ImpIntJ0(:,i)=uJ0(:,i)/condut(i)
    ImpIntJ1(:,i)=uJ1(:,i)/condut(i)
    uhJ0(:,i)=uJ0(:,i)*h(i)
    uhJ1(:,i)=uJ1(:,i)*h(i)
    tghJ0(:,i)=(1.d0-exp(-2.d0*uhJ0(:,i)))/(1.d0+exp(-2.d0*uhJ0(:,i)))
    tghJ1(:,i)=(1.d0-exp(-2.d0*uhJ1(:,i)))/(1.d0+exp(-2.d0*uhJ1(:,i)))
    end if
  end do

  allocate(ImpApdwJ0(nJ0,1:n),ImpApdwJ1(nJ1,1:n),RTMdwJ0(nJ0,0:n),RTMdwJ1(nJ1,0:n))

  do i=n,1,-1
    if (i==n) then
    ImpApdwJ0(:,i)=ImpIntJ0(:,i)
    ImpApdwJ1(:,i)=ImpIntJ1(:,i)
    RTMdwJ0(:,i)=(0.d0,0.d0)
    RTMdwJ1(:,i)=(0.d0,0.d0)
    else
    ImpApdwJ0(:,i)=ImpIntJ0(:,i)*(ImpApdwJ0(:,i+1)+ImpIntJ0(:,i)* &
      tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1)*tghJ0(:,i))
    ImpApdwJ1(:,i)=ImpIntJ1(:,i)*(ImpApdwJ1(:,i+1)+ImpIntJ1(:,i)* &
      tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1)*tghJ1(:,i))
    RTMdwJ0(:,i)=(ImpIntJ0(:,i)-ImpApdwJ0(:,i+1))/(ImpIntJ0(:,i)+ImpApdwJ0(:,i+1))
    RTMdwJ1(:,i)=(ImpIntJ1(:,i)-ImpApdwJ1(:,i+1))/(ImpIntJ1(:,i)+ImpApdwJ1(:,i+1))
    end if
  end do
    RTMdwJ0(:,0)=(ImpIntJ0(:,0)-ImpApdwJ0(:,1))/(ImpIntJ0(:,0)+ImpApdwJ0(:,1))
    RTMdwJ1(:,0)=(ImpIntJ1(:,0)-ImpApdwJ1(:,1))/(ImpIntJ1(:,0)+ImpApdwJ1(:,1))

  allocate(ImpApupJ0(nJ0,0:n-1),ImpApupJ1(nJ1,0:n-1),RTMupJ0(nJ0,0:n),RTMupJ1(nJ1,0:n))

  do i=0,n-1
    if (i==0) then
    ImpApupJ0(:,i)=ImpIntJ0(:,i)
    ImpApupJ1(:,i)=ImpIntJ1(:,i)
    RTMupJ0(:,i)=(0.d0,0.d0)
    RTMupJ1(:,i)=(0.d0,0.d0)
    else
    ImpApupJ0(:,i)=ImpIntJ0(:,i)*(ImpApupJ0(:,i-1)+ImpIntJ0(:,i)* &
      tghJ0(:,i))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1)*tghJ0(:,i))
    ImpApupJ1(:,i)=ImpIntJ1(:,i)*(ImpApupJ1(:,i-1)+ImpIntJ1(:,i)* &
      tghJ1(:,i))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1)*tghJ1(:,i))
    RTMupJ0(:,i)=(ImpIntJ0(:,i)-ImpApupJ0(:,i-1))/(ImpIntJ0(:,i)+ImpApupJ0(:,i-1))
    RTMupJ1(:,i)=(ImpIntJ1(:,i)-ImpApupJ1(:,i-1))/(ImpIntJ1(:,i)+ImpApupJ1(:,i-1))
    end if
  end do
    RTMupJ0(:,n)=(ImpIntJ0(:,n)-ImpApupJ0(:,n-1))/(ImpIntJ0(:,n)+ImpApupJ0(:,n-1))
    RTMupJ1(:,n)=(ImpIntJ1(:,n)-ImpApupJ1(:,n-1))/(ImpIntJ1(:,n)+ImpApupJ1(:,n-1))

  allocate(AMdwJ0(nJ0),AMdwJ1(nJ1),AMupJ0(nJ0),AMupJ1(nJ1))

  AMdwJ0=(exp(-uJ0(:,camadT)*(prof(camadT)-h0)) + RTMupJ0(:,camadT)* &
        exp(uJ0(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
        (1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*exp(-2.d0*uhJ0(:,camadT)))
  AMdwJ1=(exp(-uJ1(:,camadT)*(prof(camadT)-h0)) + RTMupJ1(:,camadT)* &
        exp(uJ1(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
        (1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*exp(-2.d0*uhJ1(:,camadT)))
  AMupJ0=(exp(uJ0(:,camadT)*(prof(camadT-1)-h0))+RTMdwJ0(:,camadT)* &
        exp(-uJ0(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
        (1.d0-RTMupJ0(:,camadT)*RTMdwJ0(:,camadT)*exp(-2.d0*uhJ0(:,camadT)))
  AMupJ1=(exp(uJ1(:,camadT)*(prof(camadT-1)-h0))+RTMdwJ1(:,camadT)* &
        exp(-uJ1(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
        (1.d0-RTMupJ1(:,camadT)*RTMdwJ1(:,camadT)*exp(-2.d0*uhJ1(:,camadT)))

  if (camad > camadT) then
    deallocate( TMdwJ0, TMdwJ1 )
    allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad))
    do j=camadT,camad
      if (j == camadT) then
      TMdwJ0(:,j)= Iw * dsz / (2.d0*uJ0(:,camadT))
      TMdwJ1(:,j)= Iw * dsz / (2.d0*uJ1(:,camadT))
      elseif (j == (camadT + 1) .and. j == n) then
      TMdwJ0(:,j)=TMdwJ0(:,j-1)*(exp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
        RTMupJ0(:,camadT)*AMupJ0(:)*exp(-uhJ0(:,camadT)) + &
        RTMdwJ0(:,camadT)*AMdwJ0(:))
      TMdwJ1(:,j)=TMdwJ1(:,j-1)*(exp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
        RTMupJ1(:,camadT)*AMupJ1(:)*exp(-uhJ1(:,camadT)) + &
        RTMdwJ1(:,camadT)*AMdwJ1(:))
      elseif (j == (camadT + 1) .and. j /= n) then
      TMdwJ0(:,j)=TMdwJ0(:,j-1)*(exp(-uJ0(:,camadT)*(prof(camadT)-h0)) + &
        RTMupJ0(:,camadT)*AMupJ0(:)*exp(-uhJ0(:,camadT)) + &
        RTMdwJ0(:,camadT)*AMdwJ0(:)) / (1.d0 + RTMdwJ0(:,j)*exp(-2.d0*uhJ0(:,j)))
      TMdwJ1(:,j)=TMdwJ1(:,j-1)*(exp(-uJ1(:,camadT)*(prof(camadT)-h0)) + &
        RTMupJ1(:,camadT)*AMupJ1(:)*exp(-uhJ1(:,camadT)) + &
        RTMdwJ1(:,camadT)*AMdwJ1(:)) / (1.d0 + RTMdwJ1(:,j)*exp(-2.d0*uhJ1(:,j)))
      elseif (j /= n) then
      TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*exp(-uhJ0(:,j-1))/ &
            (1.d0 + RTMdwJ0(:,j)*exp(-2.d0*uhJ0(:,j)))
      TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*exp(-uhJ1(:,j-1))/ &
            (1.d0 + RTMdwJ1(:,j)*exp(-2.d0*uhJ1(:,j)))
      elseif (j==n) then
      TMdwJ0(:,j)=TMdwJ0(:,j-1)*(1.d0+RTMdwJ0(:,j-1))*exp(-uhJ0(:,j-1))
      TMdwJ1(:,j)=TMdwJ1(:,j-1)*(1.d0+RTMdwJ1(:,j-1))*exp(-uhJ1(:,j-1))
      end if
    end do
  elseif (camad < camadT) then
    deallocate(TMupJ0, TMupJ1)
    allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT))
    do j=camadT,camad,-1
      if (j == camadT) then
      TMupJ0(:,j)= Iw * dsz / (2.d0*uJ0(:,camadT))
      TMupJ1(:,j)= Iw * dsz / (2.d0*uJ1(:,camadT))
      elseif (j == (camadT - 1) .and. j == 0) then
      TMupJ0(:,j)=TMupJ0(:,j+1)*(exp(-uJ0(:,camadT)*h0) + &
        RTMupJ0(:,camadT)*AMupJ0(:) + RTMdwJ0(:,camadT)*AMdwJ0(:)* &
        exp(-uhJ0(:,camadT)))
      TMupJ1(:,j)=TMupJ1(:,j+1)*(exp(-uJ1(:,camadT)*h0) + &
        RTMupJ1(:,camadT)*AMupJ1(:) + RTMdwJ1(:,camadT)*AMdwJ1(:)* &
        exp(-uhJ1(:,camadT)))
      elseif (j == (camadT - 1) .and. j /= 0) then
      TMupJ0(:,j)=TMupJ0(:,j+1)*(exp(uJ0(:,camadT)*(prof(camadT-1)-h0)) + &
        RTMupJ0(:,camadT)*AMupJ0(:) + RTMdwJ0(:,camadT)*AMdwJ0(:)* &
        exp(-uhJ0(:,camadT))) / (1.d0+RTMupJ0(:,j)*exp(-2.d0*uhJ0(:,j)))
      TMupJ1(:,j)=TMupJ1(:,j+1)*(exp(uJ1(:,camadT)*(prof(camadT-1)-h0)) + &
        RTMupJ1(:,camadT)*AMupJ1(:) + RTMdwJ1(:,camadT)*AMdwJ1(:)* &
        exp(-uhJ1(:,camadT))) / (1.d0+RTMupJ1(:,j)*exp(-2.d0*uhJ1(:,j)))
      elseif (j /= 0) then
      TMupJ0(:,j)=TMupJ0(:,j+1)*(1.d0+RTMupJ0(:,j+1))*exp(-uhJ0(:,j+1)) / &
        (1.d0 + RTMupJ0(:,j)*exp(-2.d0*uhJ0(:,j)))
      TMupJ1(:,j)=TMupJ1(:,j+1)*(1.d0+RTMupJ1(:,j+1))*exp(-uhJ1(:,j+1)) / &
        (1.d0 + RTMupJ1(:,j)*exp(-2.d0*uhJ1(:,j)))
      elseif (j == 0) then
      TMupJ0(:,j)=TMupJ0(:,1)*(1.d0+RTMupJ0(:,1))*exp(-uhJ0(:,1))
      TMupJ1(:,j)=TMupJ1(:,1)*(1.d0+RTMupJ1(:,1))*exp(-uhJ1(:,1))
      end if
    end do
  else
    deallocate( TMdwJ0, TMdwJ1, TMupJ0,TMupJ1 )
    allocate(TMdwJ0(nJ0,camadT:camad),TMdwJ1(nJ1,camadT:camad))
    allocate(TMupJ0(nJ0,camad:camadT),TMupJ1(nJ1,camad:camadT))
      TMdwJ0(:,camad)= Iw * dsz / (2.d0*uJ0(:,camadT))
      TMdwJ1(:,camad)= Iw * dsz / (2.d0*uJ1(:,camadT))
      TMupJ0(:,camad)= Iw * dsz / (2.d0*uJ0(:,camadT))
      TMupJ1(:,camad)= Iw * dsz / (2.d0*uJ1(:,camadT))
  end if

  allocate(Ktmdz_J1(nJ1),Ktm_J0(nJ0),Ktm_J1(nJ1))
  allocate(kernelExJ1(nJ1),kernelEyJ1(nJ1),kernelEzJ0(nJ0))
  allocate(kernelHxJ1(nJ1),kernelHyJ1(nJ1))
  if (camad == 0 .and. camadT /= 0) then
    Ktmdz_J1=(ImpIntJ1(:,0)*TMupJ1(:,0)*exp(uJ1(:,0)*z))*w_J1(:)
    Ktm_J0=(TMupJ0(:,0)*exp(uJ0(:,0)*z))*w_J0(:)
    Ktm_J1=(TMupJ1(:,0)*exp(uJ1(:,0)*z))*w_J1(:)

    kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
    Ex_p = - sum(kernelExJ1) / (2.d0*pi*r)    !este último r é decorrente do uso dos filtros

    kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
    Ey_p = - sum(kernelEyJ1) / (2.d0*pi*r)    !este último r é decorrente do uso dos filtros

    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*neta) !este último r é decorrente do uso dos filtros

    kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
    Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
    Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

  elseif (camad < camadT) then !camada k
    ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(exp(uJ1(:,camad)*(z-prof(camad))) - &
      RTMupJ1(:,camad)*exp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)
    ktm_J0=(TMupJ0(:,camad)*(exp(uJ0(:,camad)*(z-prof(camad))) + &
      RTMupJ0(:,camad)*exp(-uJ0(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J0(:)
    ktm_J1=(TMupJ1(:,camad)*(exp(uJ1(:,camad)*(z-prof(camad))) + &
      RTMupJ1(:,camad)*exp(-uJ1(:,camad)*(z-prof(camad-1)+h(camad)))))*w_J1(:)

    kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
    Ex_p = - sum(kernelExJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
    Ey_p = - sum(kernelEyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))  !este último r é decorrente do uso dos filtros

    kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
    Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
    Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

  elseif (camad == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
    Ktmdz_J1=(ImpIntJ1(:,camad)*TMupJ1(:,camad)*(exp(uJ1(:,camad)*(z-h0)) - &
      RTMupJ1(:,camad)*AMupJ1(:)*exp(-uJ1(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ1(:,camad)*AMdwJ1(:)*exp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
    Ktm_J0=(TMupJ0(:,camad)*(exp(uJ0(:,camad)*(z-h0)) + &
      RTMupJ0(:,camad)*AMupJ0(:)*exp(-uJ0(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ0(:,camad)*AMdwJ0(:)*exp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
    Ktm_J1=(TMupJ1(:,camad)*(exp(uJ1(:,camad)*(z-h0)) + &
      RTMupJ1(:,camad)*AMupJ1(:)*exp(-uJ1(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ1(:,camad)*AMdwJ1(:)*exp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

    kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
    Ex_p = - sum(kernelExJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
    Ey_p = - sum(kernelEyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    if (camad /= 0) then
    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))  !este último r é decorrente do uso dos filtros
    else
    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*neta) !este último r é decorrente do uso dos filtros
    end if

    kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
    Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
    Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

  elseif (camad == camadT .and. z > h0) then  !na mesma camada do transmissor mas abaixo dele
    Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(exp(-uJ1(:,camad)*(z-h0)) + &
      RTMupJ1(:,camad)*AMupJ1(:)*exp(-uJ1(:,camad)*(z-prof(camad-1))) - &
      RTMdwJ1(:,camad)*AMdwJ1(:)*exp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)
    Ktm_J0=(TMdwJ0(:,camad)*(exp(-uJ0(:,camad)*(z-h0)) + &
      RTMupJ0(:,camad)*AMupJ0(:)*exp(-uJ0(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ0(:,camad)*AMdwJ0(:)*exp(uJ0(:,camad)*(z-prof(camad)))))*w_J0(:)
    Ktm_J1=(TMdwJ1(:,camad)*(exp(-uJ1(:,camad)*(z-h0)) + &
      RTMupJ1(:,camad)*AMupJ1(:)*exp(-uJ1(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ1(:,camad)*AMdwJ1(:)*exp(uJ1(:,camad)*(z-prof(camad)))))*w_J1(:)

    kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
    Ex_p = sum(kernelExJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
    Ey_p = sum(kernelEyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    if (camad /= 0) then
    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))  !este último r é decorrente do uso dos filtros
    else
    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*neta) !este último r é decorrente do uso dos filtros
    end if

    kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
    Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
    Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

  elseif (camad > camadT .and. camad /= n) then !camada j
    Ktmdz_J1=(ImpIntJ1(:,camad)*TMdwJ1(:,camad)*(exp(-uJ1(:,camad)*(z-prof(camad-1))) - &
      RTMdwJ1(:,camad)*exp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)
    Ktm_J0=(TMdwJ0(:,camad)*(exp(-uJ0(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ0(:,camad)*exp(uJ0(:,camad)*(z-prof(camad)-h(camad)))))*w_J0(:)
    Ktm_J1=(TMdwJ1(:,camad)*(exp(-uJ1(:,camad)*(z-prof(camad-1))) + &
      RTMdwJ1(:,camad)*exp(uJ1(:,camad)*(z-prof(camad)-h(camad)))))*w_J1(:)

    kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
    Ex_p = sum(kernelExJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
    Ey_p = sum(kernelEyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))  !este último r é decorrente do uso dos filtros

    kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
    Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
    Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

  else  !camada n
    Ktmdz_J1=(ImpIntJ1(:,n)*TMdwJ1(:,n)*exp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)
    Ktm_J0=(TMdwJ0(:,n)*exp(-uJ0(:,n)*(z-prof(n-1))))*w_J0(:)
    Ktm_J1=(TMdwJ1(:,n)*exp(-uJ1(:,n)*(z-prof(n-1))))*w_J1(:)

    kernelExJ1 = x/r * Ktmdz_J1 * krJ1 * krJ1
    Ex_p = sum(kernelExJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEyJ1 = y/r * Ktmdz_J1 * krJ1 * krJ1
    Ey_p = sum(kernelEyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelEzJ0 = Ktm_J0 * krJ0 * krJ0 * krJ0
    Ez_p = sum(kernelEzJ0) / (2.d0*pi*r*condut(camad))  !este último r é decorrente do uso dos filtros

    kernelHxJ1 = y/r * Ktm_J1 * krJ1 * krJ1
    Hx_p = - sum(kernelHxJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

    kernelHyJ1 = x/r * Ktm_J1 * krJ1 * krJ1
    Hy_p = sum(kernelHyJ1) / (2.d0*pi*r)  !este último r é decorrente do uso dos filtros

  end if

  deallocate(h,KrJ0,KrJ1,w_J0,w_J1)
  deallocate(wvnb2,uJ0,uJ1,ImpIntJ0,ImpIntJ1,uhJ0,uhJ1,tghJ0,tghJ1)
  deallocate(ImpApdwJ0,ImpApdwJ1,RTMdwJ0,RTMdwJ1)
  deallocate(ImpApupJ0,ImpApupJ1,RTMupJ0,RTMupJ1)
  deallocate(AMdwJ0,AMdwJ1,AMupJ0,AMupJ1)
  deallocate(Ktmdz_J1,Ktm_J0,Ktm_J1)
  deallocate(kernelExJ1,kernelEyJ1,kernelEzJ0)
  deallocate(kernelHxJ1,kernelHyJ1)
  end subroutine ved_xyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
  subroutine ved_xkyz_loops(Tx,ky,h0,n,esp,condut,neta,zeta,cx,z,Ex_ky,Ey_ky,Ez_ky,Hx_ky,Hy_ky,Hz_ky)
  implicit none
  integer,intent(in)::n
  real(dp),intent(in)::Tx,ky,h0,esp(:),condut(1:n),cx,z    !,neta
  complex(dp),intent(in)::zeta,neta
  complex(dp),intent(out)::Ex_ky,Ey_ky,Ez_ky,Hx_ky,Hy_ky,Hz_ky

  integer::i,j,k,camad,camadT,autor,filtro,npts,nptc,funs,func
  real(dp)::x
  real(dp),dimension(:),allocatable::h,kxsen,kxcos,kr2sen,kr2cos,w_sen,w_cos,prof

  complex(dp),dimension(:),allocatable::wvnb2,AMdwSen,AMdwCos,AMupSen,AMupCos
  complex(dp),dimension(:,:),allocatable::uSen,uCos,AdmIntSen,AdmIntCos,ImpIntSen,ImpIntCos
  complex(dp),dimension(:,:),allocatable::uhSen,uhCos,tghSen,tghCos
  complex(dp),dimension(:,:),allocatable::ImpApdwSen,ImpApdwCos
  complex(dp),dimension(:,:),allocatable::RTMdwSen,RTMdwCos
  complex(dp),dimension(:,:),allocatable::ImpApupSen,ImpApupCos
  complex(dp),dimension(:,:),allocatable::RTMupSen,RTMupCos
  complex(dp),dimension(:,:),allocatable::TMdwSen,TMdwCos,TMupSen,TMupCos
  complex(dp),dimension(:),allocatable::Ktmdz_Sen,Ktmdz_Cos,Ktm_Sen,Ktm_Cos,Kte_Sen,Kte_Cos,Ktedz_Sen,Ktedz_Cos
  complex(dp),dimension(:),allocatable::kernelEx,kernelEy,kernelEz,kernelHx,kernelHy
  real(dp), parameter :: eps = 1.d-7
  Hz_ky = (0.d0,0.d0)
  if ( dabs(cx-Tx) < eps ) then
    x=1.d-1
  else
    x=cx-Tx
  end if

  allocate(h(0:n),prof(-1:n))
  if (size(esp)==n) then
    h(0)=0.d0
    h(1:n)=esp
  else
    h(0)=0.d0
    h(1:n-1)=esp
    h(n)=1.d300
  end if
! criando um novo vetor de profundidades que se adeque à qualquer situação patológica
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

!para descobrir em que camada está a observação
  camad = 0
  if (z <= 0.d0) then
    camad=0
  else if (z > prof(n-1)) then
    camad=n
  else
  do i=n-1,1,-1
  if (z > prof(i-1)) then
    camad=i
    exit
  end if
  end do
  end if

!para descobrir em que camada está o transmissor
  camadT = 0
  if (h0 <= 0.d0) then
    camadT = 0
  else if (h0 > prof(n-1)) then
    camadT = n
  else
  do j=n-1,1,-1
  if (h0 > prof(j-1)) then
    camadT = j
    exit
  end if
  end do
  end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! to workaround the warning: ... may be used uninitialized in this function
allocate( TMdwSen(1,1), TMdwCos(1,1) )
allocate( TMupSen(1,1), TMupCos(1,1) )
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  filtro=1  !Designa o tipo de filtro usado na subrotina de pesos e abscisas de vários filtros.
        !O algarismo 0 é usado para J0 e J1, enquanto 1 é para seno e cosseno.
  autor = 4   !Designa o criador do filtro. No caso de seno ou cosseno os que possuo são os do Frayzer (1) e do Kerry Key (4)
  funs = 2    !Designa se é filtro seno (2) ou filtro cosseno (3)
  func = 3    !Designa se é filtro seno (2) ou filtro cosseno (3)
  npts = 241    !Designa o número de pontos usado no filtro seno.
  nptc = 241    !Designa o número de pontos usado no filtro cosseno.

  allocate(kxsen(npts),kxcos(nptc),w_sen(npts),w_cos(nptc))

  call constfiltro(filtro,autor,funs,npts,x,Kxsen,w_sen)
  call constfiltro(filtro,autor,func,nptc,x,Kxcos,w_cos)

  allocate(wvnb2(0:n),uSen(npts,0:n),uCos(nptc,0:n),AdmIntSen(npts,0:n),AdmIntCos(nptc,0:n),ImpIntSen(npts,0:n),ImpIntCos(nptc,0:n))
  allocate(uhSen(npts,0:n),uhCos(nptc,0:n),tghSen(npts,0:n),tghCos(nptc,0:n),kr2sen(npts),kr2cos(nptc))

  kr2sen = kxsen*kxsen+ky*ky
  kr2cos = kxcos*kxcos+ky*ky
  do i=0,n
    if (i==0) then
    wvnb2(i)=-zeta*neta
    uSen(:,i)=sqrt(kr2sen-wvnb2(i))
    uCos(:,i)=sqrt(kr2cos-wvnb2(i))
    ImpIntSen(:,i)=uSen(:,i)/neta
    ImpIntCos(:,i)=uCos(:,i)/neta
    uhSen(:,i)=uSen(:,i)*h(i)
    uhCos(:,i)=uCos(:,i)*h(i)
    tghSen(:,i)=(1.d0-exp(-2.d0*uhSen(:,i)))/(1.d0+exp(-2.d0*uhSen(:,i)))
    tghCos(:,i)=(1.d0-exp(-2.d0*uhCos(:,i)))/(1.d0+exp(-2.d0*uhCos(:,i)))
    else
    wvnb2(i) = -zeta*condut(i)
    uSen(:,i)=sqrt(kr2Sen-wvnb2(i))
    uCos(:,i)=sqrt(kr2Cos-wvnb2(i))
    ImpIntSen(:,i)=uSen(:,i)/condut(i)
    ImpIntCos(:,i)=uCos(:,i)/condut(i)
    uhSen(:,i)=uSen(:,i)*h(i)
    uhCos(:,i)=uCos(:,i)*h(i)
    tghSen(:,i)=(1.d0-exp(-2.d0*uhSen(:,i)))/(1.d0+exp(-2.d0*uhSen(:,i)))
    tghCos(:,i)=(1.d0-exp(-2.d0*uhCos(:,i)))/(1.d0+exp(-2.d0*uhCos(:,i)))
    end if
  end do

  allocate(ImpApdwSen(npts,1:n),ImpApdwCos(nptc,1:n))
  allocate(RTMdwSen(npts,0:n),RTMdwCos(nptc,0:n))

  do i=n,1,-1
    if (i==n) then
    ImpApdwSen(:,i)=ImpIntSen(:,i)
    ImpApdwCos(:,i)=ImpIntCos(:,i)
    RTMdwSen(:,i)=(0.d0,0.d0)
    RTMdwCos(:,i)=(0.d0,0.d0)
    else
    ImpApdwSen(:,i)=ImpIntSen(:,i)*(ImpApdwSen(:,i+1)+ImpIntSen(:,i)* &
      tghSen(:,i))/(ImpIntSen(:,i)+ImpApdwSen(:,i+1)*tghSen(:,i))
    ImpApdwCos(:,i)=ImpIntCos(:,i)*(ImpApdwCos(:,i+1)+ImpIntCos(:,i)* &
      tghCos(:,i))/(ImpIntCos(:,i)+ImpApdwCos(:,i+1)*tghCos(:,i))
    RTMdwSen(:,i)=(ImpIntSen(:,i)-ImpApdwSen(:,i+1))/(ImpIntSen(:,i)+ImpApdwSen(:,i+1))
    RTMdwCos(:,i)=(ImpIntCos(:,i)-ImpApdwCos(:,i+1))/(ImpIntCos(:,i)+ImpApdwCos(:,i+1))
    end if
  end do
    RTMdwSen(:,0)=(ImpIntSen(:,0)-ImpApdwSen(:,1))/(ImpIntSen(:,0)+ImpApdwSen(:,1))
    RTMdwCos(:,0)=(ImpIntCos(:,0)-ImpApdwCos(:,1))/(ImpIntCos(:,0)+ImpApdwCos(:,1))

  allocate(ImpApupSen(npts,0:n-1),ImpApupCos(nptc,0:n-1))
  allocate(RTMupSen(npts,0:n),RTMupCos(nptc,0:n))

  do i=0,n-1
    if (i==0) then
    ImpApupSen(:,i)=ImpIntSen(:,i)
    ImpApupCos(:,i)=ImpIntCos(:,i)
    RTMupSen(:,i)=(0.d0,0.d0)
    RTMupCos(:,i)=(0.d0,0.d0)
    else
    ImpApupSen(:,i)=ImpIntSen(:,i)*(ImpApupSen(:,i-1)+ImpIntSen(:,i)* &
      tghSen(:,i))/(ImpIntSen(:,i)+ImpApupSen(:,i-1)*tghSen(:,i))
    ImpApupCos(:,i)=ImpIntCos(:,i)*(ImpApupCos(:,i-1)+ImpIntCos(:,i)* &
      tghCos(:,i))/(ImpIntCos(:,i)+ImpApupCos(:,i-1)*tghCos(:,i))
    RTMupSen(:,i)=(ImpIntSen(:,i)-ImpApupSen(:,i-1))/(ImpIntSen(:,i)+ImpApupSen(:,i-1))
    RTMupCos(:,i)=(ImpIntCos(:,i)-ImpApupCos(:,i-1))/(ImpIntCos(:,i)+ImpApupCos(:,i-1))
    end if
  end do
    RTMupSen(:,n)=(ImpIntSen(:,n)-ImpApupSen(:,n-1))/(ImpIntSen(:,n)+ImpApupSen(:,n-1))
    RTMupCos(:,n)=(ImpIntCos(:,n)-ImpApupCos(:,n-1))/(ImpIntCos(:,n)+ImpApupCos(:,n-1))

  allocate(AMdwSen(npts),AMdwCos(nptc),AMupSen(npts),AMupCos(nptc))

  AMdwSen=(exp(-uSen(:,camadT)*(prof(camadT)-h0)) + RTMupSen(:,camadT)* &
        exp(uSen(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
        (1.d0-RTMupSen(:,camadT)*RTMdwSen(:,camadT)*exp(-2.d0*uhSen(:,camadT)))
  AMdwCos=(exp(-uCos(:,camadT)*(prof(camadT)-h0)) + RTMupCos(:,camadT)* &
        exp(uCos(:,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
        (1.d0-RTMupCos(:,camadT)*RTMdwCos(:,camadT)*exp(-2.d0*uhCos(:,camadT)))
  AMupSen=(exp(uSen(:,camadT)*(prof(camadT-1)-h0)) + RTMdwSen(:,camadT)* &
        exp(-uSen(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
        (1.d0-RTMupSen(:,camadT)*RTMdwSen(:,camadT)*exp(-2.d0*uhSen(:,camadT)))
  AMupCos=(exp(uCos(:,camadT)*(prof(camadT-1)-h0)) + RTMdwCos(:,camadT)* &
        exp(-uCos(:,camadT)*(prof(camadT)+h(camadT)-h0))) / &
        (1.d0-RTMupCos(:,camadT)*RTMdwCos(:,camadT)*exp(-2.d0*uhCos(:,camadT)))

  if (camad > camadT) then
    deallocate( TMdwSen, TMdwCos )
    allocate(TMdwSen(npts,camadT:camad),TMdwCos(nptc,camadT:camad))
    do j=camadT,camad
      if (j == camadT) then
      TMdwSen(:,j)= Iw * dsz / (2.d0*uSen(:,camadT))
      TMdwCos(:,j)= Iw * dsz / (2.d0*uCos(:,camadT))
      elseif (j == (camadT + 1) .and. j == n) then
      TMdwSen(:,j)=TMdwSen(:,j-1)*(exp(-uSen(:,camadT)*(prof(camadT)-h0)) + &
        RTMupSen(:,camadT)*AMupSen(:)*exp(-uhSen(:,camadT)) + &
        RTMdwSen(:,camadT)*AMdwSen(:))
      TMdwCos(:,j)=TMdwCos(:,j-1)*(exp(-uCos(:,camadT)*(prof(camadT)-h0)) + &
        RTMupCos(:,camadT)*AMupCos(:)*exp(-uhCos(:,camadT)) + &
        RTMdwCos(:,camadT)*AMdwCos(:))
      elseif (j == (camadT + 1) .and. j /= n) then
      TMdwSen(:,j)=TMdwSen(:,j-1)*(exp(-uSen(:,camadT)*(prof(camadT)-h0)) + &
        RTMupSen(:,camadT)*AMupSen(:)*exp(-uhSen(:,camadT)) + &
        RTMdwSen(:,camadT)*AMdwSen(:)) / (1.d0 + RTMdwSen(:,j)*exp(-2.d0*uhSen(:,j)))
      TMdwCos(:,j)=TMdwCos(:,j-1)*(exp(-uCos(:,camadT)*(prof(camadT)-h0)) + &
        RTMupCos(:,camadT)*AMupCos(:)*exp(-uhCos(:,camadT)) + &
        RTMdwCos(:,camadT)*AMdwCos(:)) / (1.d0 + RTMdwCos(:,j)*exp(-2.d0*uhCos(:,j)))
      elseif (j /= n) then
      TMdwSen(:,j)=TMdwSen(:,j-1)*(1.d0+RTMdwSen(:,j-1))*exp(-uhSen(:,j-1))/ &
            (1.d0 + RTMdwSen(:,j)*exp(-2.d0*uhSen(:,j)))
      TMdwCos(:,j)=TMdwCos(:,j-1)*(1.d0+RTMdwCos(:,j-1))*exp(-uhCos(:,j-1))/ &
            (1.d0 + RTMdwCos(:,j)*exp(-2.d0*uhCos(:,j)))
      elseif (j==n) then
      TMdwSen(:,j)=TMdwSen(:,j-1)*(1.d0+RTMdwSen(:,j-1))*exp(-uhSen(:,j-1))
      TMdwCos(:,j)=TMdwCos(:,j-1)*(1.d0+RTMdwCos(:,j-1))*exp(-uhCos(:,j-1))
      end if
    end do
  elseif (camad < camadT) then
    deallocate( TMupSen, TMupCos )
    allocate(TMupSen(npts,camad:camadT),TMupCos(nptc,camad:camadT))
    do j=camadT,camad,-1
      if (j == camadT) then
      TMupSen(:,j)= Iw * dsz / (2.d0*uSen(:,camadT))
      TMupCos(:,j)= Iw * dsz / (2.d0*uCos(:,camadT))
      elseif (j == (camadT - 1) .and. j == 0) then
      TMupSen(:,j)=TMupSen(:,j+1)*(exp(-uSen(:,camadT)*h0) + &
        RTMupSen(:,camadT)*AMupSen(:) + RTMdwSen(:,camadT)*AMdwSen(:)* &
        exp(-uhSen(:,camadT)))
      TMupCos(:,j)=TMupCos(:,j+1)*(exp(-uCos(:,camadT)*h0) + &
        RTMupCos(:,camadT)*AMupCos(:) + RTMdwCos(:,camadT)*AMdwCos(:)* &
        exp(-uhCos(:,camadT)))
      elseif (j == (camadT - 1) .and. j /= 0) then
      TMupSen(:,j)=TMupSen(:,j+1)*(exp(uSen(:,camadT)*(prof(camadT-1)-h0)) + &
        RTMupSen(:,camadT)*AMupSen(:) + RTMdwSen(:,camadT)*AMdwSen(:)* &
        exp(-uhSen(:,camadT))) / (1.d0+RTMupSen(:,j)*exp(-2.d0*uhSen(:,j)))
      TMupCos(:,j)=TMupCos(:,j+1)*(exp(uCos(:,camadT)*(prof(camadT-1)-h0)) + &
        RTMupCos(:,camadT)*AMupCos(:) + RTMdwCos(:,camadT)*AMdwCos(:)* &
        exp(-uhCos(:,camadT))) / (1.d0+RTMupCos(:,j)*exp(-2.d0*uhCos(:,j)))
      elseif (j /= 0) then
      TMupSen(:,j)=TMupSen(:,j+1)*(1.d0+RTMupSen(:,j+1))*exp(-uhSen(:,j+1)) / &
        (1.d0 + RTMupSen(:,j)*exp(-2.d0*uhSen(:,j)))
      TMupCos(:,j)=TMupCos(:,j+1)*(1.d0+RTMupCos(:,j+1))*exp(-uhCos(:,j+1)) / &
        (1.d0 + RTMupCos(:,j)*exp(-2.d0*uhCos(:,j)))
      elseif (j == 0) then
      TMupSen(:,j)=TMupSen(:,1)*(1.d0+RTMupSen(:,1))*exp(-uhSen(:,1))
      TMupCos(:,j)=TMupCos(:,1)*(1.d0+RTMupCos(:,1))*exp(-uhCos(:,1))
      end if
    end do
  else
    deallocate( TMdwSen, TMdwCos, TMupSen,TMupCos )
    allocate(TMdwSen(npts,camadT:camad),TMdwCos(nptc,camadT:camad))
    allocate(TMupSen(npts,camad:camadT),TMupCos(nptc,camad:camadT))
      TMdwSen(:,camad)= Iw * dsz / (2.d0*uSen(:,camadT))
      TMdwCos(:,camad)= Iw * dsz / (2.d0*uCos(:,camadT))
      TMupSen(:,camad)= Iw * dsz / (2.d0*uSen(:,camadT))
      TMupCos(:,camad)= Iw * dsz / (2.d0*uCos(:,camadT))
  end if

  allocate(Ktmdz_Sen(npts),Ktmdz_Cos(nptc),Ktm_Sen(npts),Ktm_Cos(nptc),Kte_Sen(npts),Kte_Cos(nptc),Ktedz_Sen(npts),Ktedz_Cos(nptc))
  allocate(kernelEx(npts),kernelEy(nptc),kernelEz(nptc))
  allocate(kernelHx(nptc),kernelHy(npts))
  if (camad == 0 .and. camadT /= 0) then
    Ktmdz_Sen=ImpIntSen(:,0)*TMupSen(:,0)*exp(uSen(:,0)*z)
    Ktmdz_Cos=ImpIntCos(:,0)*TMupCos(:,0)*exp(uCos(:,0)*z)
    Ktm_Sen=TMupSen(:,0)*exp(uSen(:,0)*z)
    Ktm_Cos=TMupCos(:,0)*exp(uCos(:,0)*z)

    kernelEx = (0.d0,1.d0) * kxsen * Ktmdz_Sen * w_sen
    Ex_ky = (0.d0,1.d0) * sum(kernelEx) / (pi*x)

    kernelEy = (0.d0,1.d0) * ky * Ktmdz_Cos * w_cos
    Ey_ky = sum(kernelEy) / (pi*dabs(x))

    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*neta)

    kernelHx =  (0.d0,1.d0) * ky * Ktm_Cos * w_cos
    Hx_ky = sum(kernelHx) / (pi*dabs(x))

    kernelHy = (0.d0,-1.d0) * kxsen * Ktm_Sen * w_sen
    Hy_ky = (0.d0,1.d0) * sum(kernelHy) / (pi*x)

  elseif (camad < camadT) then !camada k
    ktmdz_Sen=ImpIntSen(:,camad)*TMupSen(:,camad)*(exp(uSen(:,camad)*(z-prof(camad))) - &
      RTMupSen(:,camad)*exp(-uSen(:,camad)*(z-prof(camad-1)+h(camad))))
    ktmdz_Cos=ImpIntCos(:,camad)*TMupCos(:,camad)*(exp(uCos(:,camad)*(z-prof(camad))) - &
      RTMupCos(:,camad)*exp(-uCos(:,camad)*(z-prof(camad-1)+h(camad))))
    ktm_Sen=TMupSen(:,camad)*(exp(uSen(:,camad)*(z-prof(camad))) + &
      RTMupSen(:,camad)*exp(-uSen(:,camad)*(z-prof(camad-1)+h(camad))))
    ktm_Cos=TMupCos(:,camad)*(exp(uCos(:,camad)*(z-prof(camad))) + &
      RTMupCos(:,camad)*exp(-uCos(:,camad)*(z-prof(camad-1)+h(camad))))

    kernelEx = (0.d0,1.d0) * kxsen * Ktmdz_Sen * w_sen
    Ex_ky = (0.d0,1.d0) * sum(kernelEx) / (pi*x)

    kernelEy = (0.d0,1.d0) * ky * Ktmdz_Cos * w_cos
    Ey_ky = sum(kernelEy) / (pi*dabs(x))

    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*condut(camad))

    kernelHx =  (0.d0,1.d0) * ky * Ktm_Cos * w_cos
    Hx_ky = sum(kernelHx) / (pi*dabs(x))

    kernelHy = (0.d0,-1.d0) * kxsen * Ktm_Sen * w_sen
    Hy_ky = (0.d0,1.d0) * sum(kernelHy) / (pi*x)

  elseif (camad == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
    Ktmdz_Sen=ImpIntSen(:,camad)*TMupSen(:,camad)*(exp(uSen(:,camad)*(z-h0)) - &
      RTMupSen(:,camad)*AMupSen(:)*exp(-uSen(:,camad)*(z-prof(camad-1))) + &
      RTMdwSen(:,camad)*AMdwSen(:)*exp(uSen(:,camad)*(z-prof(camad))))
    Ktmdz_Cos=ImpIntCos(:,camad)*TMupCos(:,camad)*(exp(uCos(:,camad)*(z-h0)) - &
      RTMupCos(:,camad)*AMupCos(:)*exp(-uCos(:,camad)*(z-prof(camad-1))) + &
      RTMdwCos(:,camad)*AMdwCos(:)*exp(uCos(:,camad)*(z-prof(camad))))
    Ktm_Sen=TMupSen(:,camad)*(exp(uSen(:,camad)*(z-h0)) + &
      RTMupSen(:,camad)*AMupSen(:)*exp(-uSen(:,camad)*(z-prof(camad-1))) + &
      RTMdwSen(:,camad)*AMdwSen(:)*exp(uSen(:,camad)*(z-prof(camad))))
    Ktm_Cos=TMupCos(:,camad)*(exp(uCos(:,camad)*(z-h0)) + &
      RTMupCos(:,camad)*AMupCos(:)*exp(-uCos(:,camad)*(z-prof(camad-1))) + &
      RTMdwCos(:,camad)*AMdwCos(:)*exp(uCos(:,camad)*(z-prof(camad))))

    kernelEx = (0.d0,1.d0) * kxsen * Ktmdz_Sen * w_sen
    Ex_ky = (0.d0,1.d0) * sum(kernelEx) / (pi*x)

    kernelEy = (0.d0,1.d0) * ky * Ktmdz_Cos * w_cos
    Ey_ky = sum(kernelEy) / (pi*dabs(x))

    if (camad /= 0) then
    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*condut(camad))
    else
    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*neta)
    end if

    kernelHx =  (0.d0,1.d0) * ky * Ktm_Cos * w_cos
    Hx_ky = sum(kernelHx) / (pi*dabs(x))

    kernelHy = (0.d0,-1.d0) * kxsen * Ktm_Sen * w_sen
    Hy_ky = (0.d0,1.d0) * sum(kernelHy) / (pi*x)

  elseif (camad == camadT .and. z > h0) then  !na mesma camada do transmissor mas abaixo dele
    Ktmdz_Sen=ImpIntSen(:,camad)*TMdwSen(:,camad)*(exp(-uSen(:,camad)*(z-h0)) + &
      RTMupSen(:,camad)*AMupSen(:)*exp(-uSen(:,camad)*(z-prof(camad-1))) - &
      RTMdwSen(:,camad)*AMdwSen(:)*exp(uSen(:,camad)*(z-prof(camad))))
    Ktmdz_Cos=ImpIntCos(:,camad)*TMdwCos(:,camad)*(exp(-uCos(:,camad)*(z-h0)) + &
      RTMupCos(:,camad)*AMupCos(:)*exp(-uCos(:,camad)*(z-prof(camad-1))) - &
      RTMdwCos(:,camad)*AMdwCos(:)*exp(uCos(:,camad)*(z-prof(camad))))
    Ktm_Sen=TMdwSen(:,camad)*(exp(-uSen(:,camad)*(z-h0)) + &
      RTMupSen(:,camad)*AMupSen(:)*exp(-uSen(:,camad)*(z-prof(camad-1))) + &
      RTMdwSen(:,camad)*AMdwSen(:)*exp(uSen(:,camad)*(z-prof(camad))))
    Ktm_Cos=TMdwCos(:,camad)*(exp(-uCos(:,camad)*(z-h0)) + &
      RTMupCos(:,camad)*AMupCos(:)*exp(-uCos(:,camad)*(z-prof(camad-1))) + &
      RTMdwCos(:,camad)*AMdwCos(:)*exp(uCos(:,camad)*(z-prof(camad))))

    kernelEx = (0.d0,1.d0) * kxsen * Ktmdz_Sen * w_sen
    Ex_ky = (0.d0,-1.d0) * sum(kernelEx) / (pi*x)

    kernelEy = (0.d0,-1.d0) * ky * Ktmdz_Cos * w_cos
    Ey_ky = sum(kernelEy) / (pi*dabs(x))

    if (camad /= 0) then
    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*condut(camad))
    else
    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*neta)
    end if

    kernelHx =  (0.d0,1.d0) * ky * Ktm_Cos * w_cos
    Hx_ky = sum(kernelHx) / (pi*dabs(x))

    kernelHy = (0.d0,-1.d0) * kxsen * Ktm_Sen * w_sen
    Hy_ky = (0.d0,1.d0) * sum(kernelHy) / (pi*x)

  elseif (camad > camadT .and. camad /= n) then !camada j
    Ktmdz_Sen=ImpIntSen(:,camad)*TMdwSen(:,camad)*(exp(-uSen(:,camad)*(z-prof(camad-1))) - &
      RTMdwSen(:,camad)*exp(uSen(:,camad)*(z-prof(camad)-h(camad))))
    Ktmdz_Cos=ImpIntCos(:,camad)*TMdwCos(:,camad)*(exp(-uCos(:,camad)*(z-prof(camad-1))) - &
      RTMdwCos(:,camad)*exp(uCos(:,camad)*(z-prof(camad)-h(camad))))
    Ktm_Sen=TMdwSen(:,camad)*(exp(-uSen(:,camad)*(z-prof(camad-1))) + &
      RTMdwSen(:,camad)*exp(uSen(:,camad)*(z-prof(camad)-h(camad))))
    Ktm_Cos=TMdwCos(:,camad)*(exp(-uCos(:,camad)*(z-prof(camad-1))) + &
      RTMdwCos(:,camad)*exp(uCos(:,camad)*(z-prof(camad)-h(camad))))

    kernelEx = (0.d0,1.d0) * kxsen * Ktmdz_Sen * w_sen
    Ex_ky = (0.d0,-1.d0) * sum(kernelEx) / (pi*x)

    kernelEy = (0.d0,-1.d0) * ky * Ktmdz_Cos * w_cos
    Ey_ky = sum(kernelEy) / (pi*dabs(x))

    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*condut(camad))

    kernelHx =  (0.d0,1.d0) * ky * Ktm_Cos * w_cos
    Hx_ky = sum(kernelHx) / (pi*dabs(x))

    kernelHy = (0.d0,-1.d0) * kxsen * Ktm_Sen * w_sen
    Hy_ky = (0.d0,1.d0) * sum(kernelHy) / (pi*x)

  else  !camada n
    Ktmdz_Sen=ImpIntSen(:,n)*TMdwSen(:,n)*exp(-uSen(:,n)*(z-prof(n-1)))
    Ktmdz_Cos=ImpIntCos(:,n)*TMdwCos(:,n)*exp(-uCos(:,n)*(z-prof(n-1)))
    Ktm_Sen=TMdwSen(:,n)*exp(-uSen(:,n)*(z-prof(n-1)))
    Ktm_Cos=TMdwCos(:,n)*exp(-uCos(:,n)*(z-prof(n-1)))

    kernelEx = (0.d0,1.d0) * kxsen * Ktmdz_Sen * w_sen
    Ex_ky = (0.d0,-1.d0) * sum(kernelEx) / (pi*x)

    kernelEy = (0.d0,-1.d0) * ky * Ktmdz_Cos * w_cos
    Ey_ky = sum(kernelEy) / (pi*dabs(x))

    kernelEz = kr2cos * Ktm_Cos * w_cos
    Ez_ky = sum(kernelEz) / (pi*dabs(x)*condut(camad))

    kernelHx =  (0.d0,1.d0) * ky * Ktm_Cos * w_cos
    Hx_ky = sum(kernelHx) / (pi*dabs(x))

    kernelHy = (0.d0,-1.d0) * kxsen * Ktm_Sen * w_sen
    Hy_ky = (0.d0,1.d0) * sum(kernelHy) / (pi*x)

  end if

  deallocate(h,KxSen,KxCos,w_Sen,w_Cos)
  deallocate(wvnb2,uSen,uCos,ImpIntSen,ImpIntCos,uhSen,uhCos,tghSen,tghCos)
  deallocate(ImpApdwSen,ImpApdwCos,RTMdwSen,RTMdwCos)
  deallocate(ImpApupSen,ImpApupCos,RTMupSen,RTMupCos)
  deallocate(AMdwSen,AMdwCos,AMupSen,AMupCos)
  deallocate(Ktmdz_Sen,Ktmdz_Cos,Ktm_Sen,Ktm_Cos)
  deallocate(kernelEx,kernelEy,kernelEz,kernelHx,kernelHy)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  end subroutine ved_xkyz_loops
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
subroutine ved_xkyz(Tx,ky,h0,n,esp,condut,neta,zeta,cx,z,Ex_ky,Ey_ky,Ez_ky,Hx_ky,Hy_ky,Hz_ky)
implicit none
integer,intent(in)::n
real(dp),intent(in)::Tx,ky,h0,esp(:),condut(1:n),cx,z    !,neta
complex(dp),intent(in)::zeta,neta
complex(dp),intent(out)::Ex_ky,Ey_ky,Ez_ky,Hx_ky,Hy_ky,Hz_ky

integer::i,j,k,camad,camadT,autor,filtro,npts,nptc,funs,func
real(dp)::x,kx
real(dp),dimension(:),allocatable::h,kxsen,kxcos,w_sen,w_cos,prof

real(dp), parameter :: eps = 1.d-7

complex(dp) :: kerEx, kerEy, kerEz, kerHx, kerHy
Hz_ky = (0.d0,0.d0)

if ( dabs(cx - Tx) < eps ) then
  x=1.d-1
else
  x=cx-Tx
end if

allocate(h(0:n),prof(-1:n))
if (size(esp)==n) then
  h(0)=0.d0
  h(1:n)=esp
else
  h(0)=0.d0
  h(1:n-1)=esp
  h(n)=1.d300
end if
! criando um novo vetor de profundidades que se adeque à qualquer situação patológica
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

!para descobrir em que camada está a observação
if (z <= 0.d0) then
  camad=0
else if (z > prof(n-1)) then
  camad=n
else
do i=n-1,1,-1
if (z > prof(i-1)) then
  camad=i
  exit
end if
end do
end if

!para descobrir em que camada está o transmissor
if (h0 <= 0.d0) then
  camadT = 0
else if (h0 > prof(n-1)) then
  camadT = n
else
do j=n-1,1,-1
if (h0 > prof(j-1)) then
  camadT = j
  exit
end if
end do
end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
filtro=1  !Designa o tipo de filtro usado na subrotina de pesos e abscisas de vários filtros.
      !O algarismo 0 é usado para J0 e J1, enquanto 1 é para seno e cosseno.
autor = 1   !Designa o criador do filtro. No caso de seno ou cosseno os que possuo são os do Frayzer (1) e do Kerry Key (4)
funs = 2    !Designa se é filtro seno (2) ou filtro cosseno (3)
func = 3    !Designa se é filtro seno (2) ou filtro cosseno (3)
npts = 30   !Designa o número de pontos usado no filtro seno.
nptc = 30   !Designa o número de pontos usado no filtro cosseno.

allocate(kxsen(npts),kxcos(nptc),w_sen(npts),w_cos(nptc))

call constfiltro(filtro,autor,funs,npts,x,Kxsen,w_sen)
call constfiltro(filtro,autor,func,nptc,x,Kxcos,w_cos)
kerEx = 0.d0; kerEy = 0.d0; kerEz = 0.d0;
kerHx = 0.d0; kerHy = 0.d0;
if (camad == 0 .and. camadT /= 0) then

  kerEx = (0.d0,0.d0)
  kerHy = (0.d0,0.d0)
    do i=1,npts
      kx=kxsen(i)
      kerEx = kerEx + KEx0(kx,ky)*w_sen(i)
      kerHy = kerHy + KHy0(kx,ky)*w_sen(i)
    end do
  kerEy = (0.d0,0.d0)
  kerEz = (0.d0,0.d0)
  kerHx = (0.d0,0.d0)
    do j=1,nptc
      kx=kxcos(j)
      kerEy = kerEy + KEy0(kx,ky)*w_cos(j)
      kerEz = kerEz + KEz0(kx,ky)*w_cos(j)
      kerHx = kerHx + KHx0(kx,ky)*w_cos(j)
    end do

elseif (camad < camadT) then !camada k

  kerEx = (0.d0,0.d0)
  kerHy = (0.d0,0.d0)
    do i=1,npts
      kx=kxsen(i)
      kerEx = kerEx + KExk(kx,ky)*w_sen(i)
      kerHy = kerHy + KHyk(kx,ky)*w_sen(i)
    end do
  kerEy = (0.d0,0.d0)
  kerEz = (0.d0,0.d0)
  kerHx = (0.d0,0.d0)
    do j=1,nptc
      kx=kxcos(j)
      kerEy = kerEy + KEyk(kx,ky)*w_cos(j)
      kerEz = kerEz + KEzk(kx,ky)*w_cos(j)
      kerHx = kerHx + KHxk(kx,ky)*w_cos(j)
    end do

elseif (camad == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele

  kerEx = (0.d0,0.d0)
  kerHy = (0.d0,0.d0)
    do i=1,npts
      kx=kxsen(i)
      kerEx = kerEx + KExl_up(kx,ky)*w_sen(i)
      kerHy = kerHy + KHyl_up(kx,ky)*w_sen(i)
    end do
  kerEy = (0.d0,0.d0)
  kerEz = (0.d0,0.d0)
  kerHx = (0.d0,0.d0)
    do j=1,nptc
      kx=kxcos(j)
      kerEy = kerEy + KEyl_up(kx,ky)*w_cos(j)
      kerEz = kerEz + KEzl_up(kx,ky)*w_cos(j)
      kerHx = kerHx + KHxl_up(kx,ky)*w_cos(j)
    end do

elseif (camad == camadT .and. z > h0) then  !na mesma camada do transmissor mas abaixo dele

  kerEx = (0.d0,0.d0)
  kerHy = (0.d0,0.d0)
    do i=1,npts
      kx=kxsen(i)
      kerEx = kerEx + KExl_dw(kx,ky)*w_sen(i)
      kerHy = kerHy + KHyl_dw(kx,ky)*w_sen(i)
    end do
  kerEy = (0.d0,0.d0)
  kerEz = (0.d0,0.d0)
  kerHx = (0.d0,0.d0)
    do j=1,nptc
      kx=kxcos(j)
      kerEy = kerEy + KEyl_dw(kx,ky)*w_cos(j)
      kerEz = kerEz + KEzl_dw(kx,ky)*w_cos(j)
      kerHx = kerHx + KHxl_dw(kx,ky)*w_cos(j)
    end do

elseif (camad > camadT .and. camad /= n) then !camada j

  kerEx = (0.d0,0.d0)
  kerHy = (0.d0,0.d0)
    do i=1,npts
      kx=kxsen(i)
      kerEx = kerEx + KExj(kx,ky)*w_sen(i)
      kerHy = kerHy + KHyj(kx,ky)*w_sen(i)
    end do
  kerEy = (0.d0,0.d0)
  kerEz = (0.d0,0.d0)
  kerHx = (0.d0,0.d0)
    do j=1,nptc
      kx=kxcos(j)
      kerEy = kerEy + KEyj(kx,ky)*w_cos(j)
      kerEz = kerEz + KEzj(kx,ky)*w_cos(j)
      kerHx = kerHx + KHxj(kx,ky)*w_cos(j)
    end do

elseif (camad == n .and. camadT /= n) then  !camada n

  kerEx = (0.d0,0.d0)
  kerHy = (0.d0,0.d0)
    do i=1,npts
      kx=kxsen(i)
      kerEx = kerEx + KExn(kx,ky)*w_sen(i)
      kerHy = kerHy + KHyn(kx,ky)*w_sen(i)
    end do
  kerEy = (0.d0,0.d0)
  kerEz = (0.d0,0.d0)
  kerHx = (0.d0,0.d0)
    do j=1,nptc
      kx=kxcos(j)
      kerEy = kerEy + KEyn(kx,ky)*w_cos(j)
      kerEz = kerEz + KEzn(kx,ky)*w_cos(j)
      kerHx = kerHx + KHxn(kx,ky)*w_cos(j)
    end do

end if

  Ex_ky = (0.d0,1.d0) * kerEx / (pi*x)
  Ey_ky = kerEy / (pi*dabs(x))
  Ez_ky = kerEz / (pi*dabs(x))
  Hx_ky = kerHx / (pi*dabs(x))
  Hy_ky = (0.d0,1.d0) * kerHy / (pi*x)

contains
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function wvnb_2(cam)
implicit none
integer,intent(in)::cam
complex(dp)::wvnb_2

if (cam == 0) then
wvnb_2 = -zeta*neta !consideramos aqui a contribuição apenas da parte imaginária da admitividade, já que sigma é zero.
else
wvnb_2 = -zeta*condut(cam)
end if

return
end function wvnb_2
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function u(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::u
! Não mude a natureza. Use essa fórmula, mesmo no ar em EI. Mas, talvez seja bom modificar em EF.
u=sqrt(kr**2.d0-wvnb_2(cam))

return
end function u
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function ImpInt(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::ImpInt

if (cam /= 0) then
ImpInt=u(kr,cam)/condut(cam)
else
ImpInt=u(kr,cam)/neta
end if

return
end function ImpInt
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function uh(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::uh
! nao esquecer de construir h(0)=0
uh=u(kr,cam)*h(cam)

return
end function uh
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function tgh(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::tgh

tgh = (1.d0 - exp(-2.d0 * uh(kr,cam))) / (1.d0 + exp(-2.d0 * uh(kr,cam)))

return
end function tgh
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Imp_Ap_dw(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
integer::j
complex(dp)::Imp_Ap_dw

Imp_Ap_dw = ImpInt(kr,n)
do j=n-1,cam,-1
Imp_Ap_dw = ImpInt(kr,j) * (Imp_Ap_dw + ImpInt(kr,j) * &
    tgh(kr,j)) / (ImpInt(kr,j) + Imp_Ap_dw * tgh(kr,j))
end do

return
end function Imp_Ap_dw
!   recursive function Imp_Ap_dw(kr,cam) result(ImpAp_dw)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::ImpAp_dw
!
!   if (cam == n) then
!   ImpAp_dw = ImpInt(kr,n)
!   else
! ! do j=ncam-1,1,-1
!   ImpAp_dw = ImpInt(kr,cam) * (Imp_Ap_dw(kr,cam+1) + ImpInt(kr,cam) * &
!       tgh(kr,cam)) / (ImpInt(kr,cam) + Imp_Ap_dw(kr,cam+1) * tgh(kr,cam))
!   end if
!
!   return
!   end function Imp_Ap_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Imp_Ap_up(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
integer::j
complex(dp)::Imp_Ap_up

Imp_Ap_up = ImpInt(kr,0)
do j=1,cam
Imp_Ap_up = ImpInt(kr,j) * (Imp_Ap_up + ImpInt(kr,j) * &
    tgh(kr,j)) / (ImpInt(kr,j) + Imp_Ap_up * tgh(kr,j))
end do

return
end function Imp_Ap_up
!   recursive function Imp_Ap_up(kr,cam) result(ImpAp_up)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::ImpAp_up
!
!   if (cam == 0) then
!   ImpAp_up = ImpInt(kr,0)
!   else
! ! do j=1,ncam
!   ImpAp_up = ImpInt(kr,cam) * (Imp_Ap_up(kr,cam-1) + ImpInt(kr,cam) * &
!       tgh(kr,cam)) / (ImpInt(kr,cam) + Imp_Ap_up(kr,cam-1) * tgh(kr,cam))
!   end if
!
!   return
!   end function Imp_Ap_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function RTM_dw(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::RTM_dw

if (cam==n) then
RTM_dw = (0.d0,0.d0)
else
RTM_dw = (ImpInt(kr,cam) - Imp_Ap_dw(kr,cam+1)) / (ImpInt(kr,cam) + Imp_Ap_dw(kr,cam+1))
end if

return
end function RTM_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function RTM_up(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::RTM_up

if (cam == 0) then
RTM_up = (0.d0,0.d0)
else
RTM_up = (ImpInt(kr,cam) - Imp_Ap_up(kr,cam-1)) / (ImpInt(kr,cam) + Imp_Ap_up(kr,cam-1))
end if

return
end function RTM_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function AM_dw(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::AM_dw

AM_dw = (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*exp(u(kr,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
    (1.d0 - RTM_up(kr,camadT)*RTM_dw(kr,camadT)*exp(-2.d0*uh(kr,camadT)))

return
end function AM_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function AM_up(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::AM_up

AM_up = (exp(u(kr,camadT)*(prof(camadT-1)-h0)) + RTM_dw(kr,camadT)*exp(-u(kr,camadT)*(prof(camadT)+h(camadT)-h0))) / &
    (1.d0 - RTM_up(kr,camadT)*RTM_dw(kr,camadT)*exp(-2.d0*uh(kr,camadT)))

return
end function AM_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function TM_dw(Kr,cam)
implicit none
real(dp),intent(in)::Kr
integer,intent(in)::cam
integer::j
complex(dp)::TM_dw

TM_dw = Iw * dsz / (2.d0*u(kr,camadT))
do j=camadT+1,cam
  if (j == (camadT + 1) .and. j == n) then
  TM_dw = TM_dw * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
      RTM_dw(kr,camadT)*AM_dw(kr))
  elseif (j == (camadT + 1) .and. j /= n) then
  TM_dw = TM_dw * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
      RTM_dw(kr,camadT)*AM_dw(kr)) / (1.d0 + RTM_dw(kr,camadT+1)*exp(-2.d0*uh(kr,camadT+1)))
  elseif (j /= n) then
  TM_dw = TM_dw * (1.d0 + RTM_dw(kr,j-1)) * exp(-uh(kr,j-1)) / (1.d0 + RTM_dw(kr,j) * exp(-2.d0*uh(kr,j)))
  else
  TM_dw = TM_dw * (1.d0 + RTM_dw(kr,n-1)) * exp(-uh(kr,n-1))
  end if
end do

return
end function
!   recursive function TM_dw(Kr,cam) result(TMdw)
!   implicit none
!   real(dp),intent(in)::Kr
!   integer,intent(in)::cam
!   complex(dp)::TMdw
!
!   if (cam == camadT) then
!   TMdw = Iw * dsz / (2.d0*u(kr,camadT))
!   elseif (cam == (camadT + 1) .and. cam == n) then
!   TMdw = TM_dw(kr,cam-1) * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
!       RTM_dw(kr,camadT)*AM_dw(kr))
!   elseif (cam == (camadT + 1) .and. cam /= n) then
!   TMdw = TM_dw(kr,cam-1) * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
!       RTM_dw(kr,camadT)*AM_dw(kr)) / (1.d0 + RTM_dw(kr,camadT+1)*exp(-2.d0*uh(kr,camadT+1)))
!   elseif (cam /= n) then
!   TMdw = TM_dw(kr,cam-1) * (1.d0 + RTM_dw(kr,cam-1)) * exp(-uh(kr,cam-1)) / (1.d0 + RTM_dw(kr,cam) * exp(-2.d0*uh(kr,cam)))
!   else
!   TMdw = TM_dw(kr,n-1) * (1.d0 + RTM_dw(kr,n-1)) * exp(-uh(kr,n-1))
!   end if
!
!   return
!   end function TM_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function TM_up(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
integer::j
complex(dp)::TM_up

TM_up = Iw * dsz / (2.d0*u(kr,camadT))
do j=camadT-1,cam,-1
  if (j == (camadT - 1) .and. j == 0) then
  TM_up = TM_up * (exp(-u(kr,camadT)*h0) + RTM_up(kr,camadT)*AM_up(kr) + RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT)))
  elseif (j == (camadT - 1) .and. j /= 0) then
  TM_up = TM_up * (exp(u(kr,camadT)*(prof(camadT-1)-h0)) + RTM_up(kr,camadT)*AM_up(kr) + &
      RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT))) / (1.d0 + RTM_up(kr,camadT-1)*exp(-2.d0*uh(kr,camadT-1)))
  elseif (j /= 0) then
  TM_up = TM_up * (1.d0 + RTM_up(kr,j+1)) * exp(-uh(kr,j+1)) / (1.d0 + RTM_up(kr,j) * exp(-2.d0*uh(kr,j)))
  else
  TM_up = TM_up * (1.d0 + RTM_up(kr,1)) * exp(-uh(kr,1))
  end if
end do

return
end function TM_up
!   recursive function TM_up(kr,cam) result(TMup)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::TMup
!
!   if (cam == camadT) then
!   TMup = Iw * dsz / 2.d0
!   elseif (cam == (camadT - 1) .and. cam == 0) then
!   TMup = TM_up(kr,cam+1) * (exp(-u(kr,camadT)*h0) + RTM_up(kr,camadT)*AM_up(kr) + RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT)))
!   elseif (cam == (camadT - 1) .and. cam /= 0) then
!   TMup = TM_up(kr,cam+1) * (exp(u(kr,camadT)*(prof(camadT-1)-h0)) + RTM_up(kr,camadT)*AM_up(kr) + &
!       RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT))) / (1.d0 + RTM_up(kr,camadT-1)*exp(-2.d0*uh(kr,camadT-1)))
!   elseif (cam /= 0) then
!   TMup = TM_up(kr,cam+1) * (1.d0 + RTM_up(kr,cam+1)) * exp(-uh(kr,cam+1)) / (1.d0 + RTM_up(kr,cam) * exp(-2.d0*uh(kr,cam)))
!   else
!   TMup = TM_up(kr,1) * (1.d0 + RTM_up(kr,1)) * exp(-uh(kr,1))
!   end if
!
!   return
!   end function TM_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktmdz_0(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktmdz_0

Ktmdz_0 = ImpInt(kr,0) * TM_up(kr,0) * exp(u(kr,0) * z)

return
end function Ktmdz_0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktmdz_k(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::Ktmdz_k

ktmdz_k = ImpInt(kr,cam) * TM_up(kr,cam) * (exp(u(kr,cam) * (z-prof(cam))) - &
       RTM_up(kr,cam) * exp(-u(kr,cam) * (z-prof(cam-1)+h(cam))))

return
end function Ktmdz_k
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktmdz_l_up(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktmdz_l_up

Ktmdz_l_up = ImpInt(kr,camadT) * TM_up(kr,camadT) * (exp(u(kr,camadT) * (z-h0)) - &
      RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) + &
      RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))

return
end function Ktmdz_l_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktmdz_l_dw(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktmdz_l_dw

Ktmdz_l_dw = ImpInt(kr,camadT) * TM_dw(kr,camadT) * (exp(-u(kr,camadT) * (z-h0)) + &
  RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) - &
  RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))

return
end function Ktmdz_l_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktmdz_j(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::Ktmdz_j

Ktmdz_j = ImpInt(kr,cam) * TM_dw(kr,cam) * (exp(-u(kr,cam) * (z-prof(cam-1))) - &
    RTM_dw(kr,cam) * exp(u(kr,cam) * (z-prof(cam)-h(cam))))

return
end function Ktmdz_j
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktmdz_n(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktmdz_n

Ktmdz_n = ImpInt(kr,n) * TM_dw(kr,n) * exp(-u(kr,n) * (z-prof(n-1)))

return
end function Ktmdz_n
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktm_0(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktm_0

Ktm_0 = TM_up(kr,0) * exp(u(kr,0) * z)

return
end function Ktm_0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktm_k(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::Ktm_k

ktm_k = TM_up(kr,cam) * (exp(u(kr,cam) * (z-prof(cam))) + &
       RTM_up(kr,cam) * exp(-u(kr,cam) * (z-prof(cam-1)+h(cam))))

return
end function Ktm_k
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktm_l_up(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktm_l_up

Ktm_l_up = TM_up(kr,camadT) * (exp(u(kr,camadT) * (z-h0)) + &
      RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) + &
      RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))

return
end function Ktm_l_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktm_l_dw(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktm_l_dw

Ktm_l_dw = TM_dw(kr,camadT) * (exp(-u(kr,camadT) * (z-h0)) + &
  RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) + &
  RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))

return
end function Ktm_l_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktm_j(kr,cam)
implicit none
real(dp),intent(in)::kr
integer,intent(in)::cam
complex(dp)::Ktm_j

Ktm_j = TM_dw(kr,cam) * (exp(-u(kr,cam) * (z-prof(cam-1))) + &
    RTM_dw(kr,cam) * exp(u(kr,cam) * (z-prof(cam)-h(cam))))

return
end function Ktm_j
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function Ktm_n(kr)
implicit none
real(dp),intent(in)::kr
complex(dp)::Ktm_n

Ktm_n = TM_dw(kr,n) * exp(-u(kr,n) * (z-prof(n-1)))

return
end function Ktm_n
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEx0(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEx0

kr=dsqrt(kx*kx+ky*ky)
KEx0 = cmplx(0.d0,kx,kind=dp) * Ktmdz_0(kr)

return
end function KEx0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KExk(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KExk

kr=dsqrt(kx*kx+ky*ky)
KExk = cmplx(0.d0,kx,kind=dp) * Ktmdz_k(kr,camad)

return
end function KExk
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KExl_up(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KExl_up

kr=dsqrt(kx*kx+ky*ky)
KExl_up = cmplx(0.d0,kx,kind=dp) * Ktmdz_l_up(kr)

return
end function KExl_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KExl_dw(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KExl_dw

kr=dsqrt(kx*kx+ky*ky)
KExl_dw = cmplx(0.d0,-kx,kind=dp) * Ktmdz_l_dw(kr)

return
end function KExl_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KExj(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KExj

kr=dsqrt(kx*kx+ky*ky)
KExj = cmplx(0.d0,-kx,kind=dp) * Ktmdz_j(kr,camad)

return
end function KExj
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KExn(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KExn

kr=dsqrt(kx*kx+ky*ky)
KExn = cmplx(0.d0,-kx,kind=dp) * Ktmdz_n(kr)

return
end function KExn
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEy0(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEy0

kr=dsqrt(kx*kx+ky*ky)
KEy0 = cmplx(0.d0,ky,kind=dp) * Ktmdz_0(kr)

return
end function KEy0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEyk(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEyk

kr=dsqrt(kx*kx+ky*ky)
KEyk = cmplx(0.d0,ky,kind=dp) * Ktmdz_k(kr,camad)

return
end function KEyk
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEyl_up(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEyl_up

kr=dsqrt(kx*kx+ky*ky)
KEyl_up = cmplx(0.d0,ky,kind=dp) * Ktmdz_l_up(kr)

return
end function KEyl_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEyl_dw(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEyl_dw

kr=dsqrt(kx*kx+ky*ky)
KEyl_dw = cmplx(0.d0,-ky,kind=dp) * Ktmdz_l_dw(kr)

return
end function KEyl_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEyj(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEyj

kr=dsqrt(kx*kx+ky*ky)
KEyj = cmplx(0.d0,-ky,kind=dp) * Ktmdz_j(kr,camad)

return
end function KEyj
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEyn(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEyn

kr=dsqrt(kx*kx+ky*ky)
KEyn = cmplx(0.d0,-ky,kind=dp) * Ktmdz_n(kr)

return
end function KEyn
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEz0(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEz0

kr=dsqrt(kx*kx+ky*ky)
KEz0 = Ktm_0(kr) * kr * kr / neta

return
end function KEz0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEzk(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEzk

kr=dsqrt(kx*kx+ky*ky)
KEzk = Ktm_k(kr,camad) * kr * kr / condut(camad)

return
end function KEzk
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEzl_up(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEzl_up

kr=dsqrt(kx*kx+ky*ky)
if (camad /= 0) then
KEzl_up = Ktm_l_up(kr) * kr * kr / condut(camad)
else
KEzl_up = Ktm_l_up(kr) * kr * kr / neta
end if

return
end function KEzl_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEzl_dw(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEzl_dw

kr=dsqrt(kx*kx+ky*ky)
if (camad /= 0) then
KEzl_dw = Ktm_l_dw(kr) * kr * kr / condut(camad)
else
KEzl_dw = Ktm_l_dw(kr) * kr * kr / neta
end if

return
end function KEzl_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEzj(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEzj

kr=dsqrt(kx*kx+ky*ky)
KEzj = Ktm_j(kr,camad) * kr * kr / condut(camad)

return
end function KEzj
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KEzn(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KEzn

kr=dsqrt(kx*kx+ky*ky)
KEzn = Ktm_n(kr) * kr * kr / condut(n)

return
end function KEzn
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHx0(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHx0

kr=dsqrt(kx*kx+ky*ky)
KHx0 = cmplx(0.d0,ky,kind=dp) * Ktm_0(kr)

return
end function KHx0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHxk(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHxk

kr=dsqrt(kx*kx+ky*ky)
KHxk = cmplx(0.d0,ky,kind=dp) * Ktm_k(kr,camad)

return
end function KHxk
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHxl_up(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHxl_up

kr=dsqrt(kx*kx+ky*ky)
KHxl_up = cmplx(0.d0,ky,kind=dp) * Ktm_l_up(kr)

return
end function KHxl_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHxl_dw(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHxl_dw

kr=dsqrt(kx*kx+ky*ky)
KHxl_dw = cmplx(0.d0,ky,kind=dp) * Ktm_l_dw(kr)

return
end function KHxl_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHxj(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHxj

kr=dsqrt(kx*kx+ky*ky)
KHxj = cmplx(0.d0,ky,kind=dp) * Ktm_j(kr,camad)

return
end function KHxj
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHxn(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHxn

kr=dsqrt(kx*kx+ky*ky)
KHxn = cmplx(0.d0,ky,kind=dp) * Ktm_n(kr)

return
end function KHxn
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHy0(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHy0

kr=dsqrt(kx*kx+ky*ky)
KHy0 = cmplx(0.d0,-kx,kind=dp) * Ktm_0(kr)

return
end function KHy0
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHyk(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHyk

kr=dsqrt(kx*kx+ky*ky)
KHyk = cmplx(0.d0,-kx,kind=dp) * Ktm_k(kr,camad)

return
end function KHyk
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHyl_up(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHyl_up

kr=dsqrt(kx*kx+ky*ky)
KHyl_up = cmplx(0.d0,-kx,kind=dp) * Ktm_l_up(kr)

return
end function KHyl_up
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHyl_dw(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHyl_dw

kr=dsqrt(kx*kx+ky*ky)
KHyl_dw = cmplx(0.d0,-kx,kind=dp) * Ktm_l_dw(kr)

return
end function KHyl_dw
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHyj(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHyj

kr=dsqrt(kx*kx+ky*ky)
KHyj = cmplx(0.d0,-kx,kind=dp) * Ktm_j(kr,camad)

return
end function KHyj
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
function KHyn(kx,ky)
implicit none
real(dp),intent(in)::kx,ky
real(dp)::kr
complex(dp)::KHyn

kr=dsqrt(kx*kx+ky*ky)
KHyn = cmplx(0.d0,-kx,kind=dp) * Ktm_n(kr)

return
end function KHyn
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
end subroutine ved_xkyz
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!   subroutine ved_xyz(filAnd,Tx,Ty,h0,n,esp,condut,neta,zeta,cx,cy,z,Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p)
!   implicit none
!   integer,intent(in)::filAnd,n
!   real(dp),intent(in)::Tx,Ty,h0,esp(:),condut(1:n),cx,cy,z !,neta
!   complex(dp),intent(in)::zeta,neta
!   complex(dp),intent(out)::Ex_p,Ey_p,Ez_p,Hx_p,Hy_p,Hz_p
!
!!   real(dp),parameter::pi=3.141592653589793238462643383279502884197d0
! real(dp), parameter :: dsz = 1.d0 !momento
! real(dp), parameter :: Iw = 1.d0 !corrente eletrica
!   integer::i,j,k,camad,camadT,filtro,idtfcd_cJ0,ident_fJ0,nJ0,idtfcd_cJ1,ident_fJ1,nJ1
!   real(dp)::x,y,r,k_r
!   real(dp),dimension(:),allocatable::h,krJ0,krJ1,w_J0,w_J1,prof
!
!   complex(dp)::kerEx_J1,kerEy_J1,kerEz_J0
!   complex(dp)::kerHx_J1,kerHy_J1
!   real(dp)::TOL1
!   integer::NF
!   complex(dp)::DZHANK
!   real(dp), parameter :: eps = 1.d-7
!   Hz_p = (0.d0,0.d0)
!
!   if ( dabs(cx-Tx) < eps ) then
!     x=1.d-2
!   else
!     x=cx-Tx
!   end if
!   y=cy-Ty
!   r = dsqrt(x**2 + y**2)
!
!   allocate(h(0:n),prof(-1:n))
!   if (size(esp)==n) then
!     h(0)=0.d0
!     h(1:n)=esp
!   else
!     h(0)=0.d0
!     h(1:n-1)=esp
!     h(n)=1.d300
!   end if
! ! criando um novo vetor de profundidades que se adeque à qualquer situação patológica
!   prof(-1)=-1.d300
!   prof(0)=0.d0
!   if (n > 1) then
!    prof(1) = h(1)
!    if (n > 2) then
!     do k=2,n-1
!      prof(k) = prof(k-1) + h(k)
!     end do
!    end if
!   end if
!   prof(n)=1.d300
!
! !para descobrir em que camada está a observação
!   if (z <= 0.d0) then
!     camad=0
!   else if (z > prof(n-1)) then
!     camad=n
!   else
!   do i=n-1,1,-1
!   if (z > prof(i-1)) then
!     camad=i
!     exit
!   end if
!   end do
!   end if
!
! !para descobrir em que camada está o transmissor
!   if (h0 <= 0.d0) then
!     camadT = 0
!   else if (h0 > prof(n-1)) then
!     camadT = n
!   else
!   do j=n-1,1,-1
!   if (h0 > prof(j-1)) then
!     camadT = j
!     exit
!   end if
!   end do
!   end if
!
! if (filAnd == 1) then !1 caso se queira usar os filtros de Anderson
!
!   TOL1 = 0.d0
!   if (camad == 0 .and. camadT /= 0) then
!
!       kerEx_J1=DZHANK(1,r,KEx0_J1,TOL1,NF,1)
!       kerEy_J1=DZHANK(1,r,KEy0_J1,TOL1,NF,1)
!       kerEz_J0=DZHANK(0,r,KEz0_J0,TOL1,NF,1)
!     Ex_p = -kerEx_J1 / (2.d0*pi)
!     Ey_p = -kerEy_J1 / (2.d0*pi)
!     Ez_p = kerEz_J0 / (2.d0*pi)
!
!       kerHx_J1=DZHANK(1,r,KHx0_J1,TOL1,NF,1)
!       kerHy_J1=DZHANK(1,r,KHy0_J1,TOL1,NF,1)
!     Hx_p = -kerHx_J1 / (2.d0*pi)
!     Hy_p = kerHy_J1 / (2.d0*pi)
!
!   elseif (camad < camadT) then !camada k
!
!       kerEx_J1=DZHANK(1,r,KExk_J1,TOL1,NF,1)
!       kerEy_J1=DZHANK(1,r,KEyk_J1,TOL1,NF,1)
!       kerEz_J0=DZHANK(0,r,KEzk_J0,TOL1,NF,1)
!   Ex_p = -kerEx_J1 / (2.d0*pi)
!   Ey_p = -kerEy_J1 / (2.d0*pi)
!   Ez_p = kerEz_J0 / (2.d0*pi)
!
!       kerHx_J1=DZHANK(1,r,KHxk_J1,TOL1,NF,1)
!       kerHy_J1=DZHANK(1,r,KHyk_J1,TOL1,NF,1)
!   Hx_p = -kerHx_J1 / (2.d0*pi)
!   Hy_p = kerHy_J1 / (2.d0*pi)
!
!   elseif (camad == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
!
!       kerEx_J1=DZHANK(1,r,KExl_up_J1,TOL1,NF,1)
!       kerEy_J1=DZHANK(1,r,KEyl_up_J1,TOL1,NF,1)
!       kerEz_J0=DZHANK(0,r,KEzl_up_J0,TOL1,NF,1)
!   Ex_p = -kerEx_J1 / (2.d0*pi)
!   Ey_p = -kerEy_J1 / (2.d0*pi)
!   Ez_p = kerEz_J0 / (2.d0*pi)
!
!       kerHx_J1=DZHANK(1,r,KHxl_up_J1,TOL1,NF,1)
!       kerHy_J1=DZHANK(1,r,KHyl_up_J1,TOL1,NF,1)
!   Hx_p = -kerHx_J1 / (2.d0*pi)
!   Hy_p = kerHy_J1 / (2.d0*pi)
!
!   elseif (camad == camadT .and. z > h0) then  !na mesma camada do transmissor mas abaixo dele
!
!       kerEx_J1=DZHANK(1,r,KExl_dw_J1,TOL1,NF,1)
!       kerEy_J1=DZHANK(1,r,KEyl_dw_J1,TOL1,NF,1)
!       kerEz_J0=DZHANK(0,r,KEzl_dw_J0,TOL1,NF,1)
!   Ex_p = kerEx_J1 / (2.d0*pi)
!   Ey_p = kerEy_J1 / (2.d0*pi)
!   Ez_p = kerEz_J0 / (2.d0*pi)
!
!       kerHx_J1=DZHANK(1,r,KHxl_dw_J1,TOL1,NF,1)
!       kerHy_J1=DZHANK(1,r,KHyl_dw_J1,TOL1,NF,1)
!   Hx_p = -kerHx_J1 / (2.d0*pi)
!   Hy_p = kerHy_J1 / (2.d0*pi)
!
!   elseif (camad > camadT .and. camad /= n) then !camada j
!
!       kerEx_J1=DZHANK(1,r,KExj_J1,TOL1,NF,1)
!       kerEy_J1=DZHANK(1,r,KEyj_J1,TOL1,NF,1)
!       kerEz_J0=DZHANK(0,r,KEzj_J0,TOL1,NF,1)
!   Ex_p = kerEx_J1 / (2.d0*pi)
!   Ey_p = kerEy_J1 / (2.d0*pi)
!   Ez_p = kerEz_J0 / (2.d0*pi)
!
!       kerHx_J1=DZHANK(1,r,KHxj_J1,TOL1,NF,1)
!       kerHy_J1=DZHANK(1,r,KHyj_J1,TOL1,NF,1)
!   Hx_p = -kerHx_J1 / (2.d0*pi)
!   Hy_p = kerHy_J1 / (2.d0*pi)
!
!   elseif (camad == n .and. camadT /= n) then  !camada n
!
!       kerEx_J1=DZHANK(1,r,KExn_J1,TOL1,NF,1)
!       kerEy_J1=DZHANK(1,r,KEyn_J1,TOL1,NF,1)
!       kerEz_J0=DZHANK(0,r,KEzn_J0,TOL1,NF,1)
!   Ex_p = kerEx_J1 / (2.d0*pi)
!   Ey_p = kerEy_J1 / (2.d0*pi)
!   Ez_p = kerEz_J0 / (2.d0*pi)
!
!       kerHx_J1=DZHANK(1,r,KHxn_J1,TOL1,NF,1)
!       kerHy_J1=DZHANK(1,r,KHyn_J1,TOL1,NF,1)
!   Hx_p = -kerHx_J1 / (2.d0*pi)
!   Hy_p = kerHy_J1 / (2.d0*pi)
!
!   end if
!
! else
! !!  write(*,*)'Entre com o criador dos filtros J0: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
! !!  read(*,*)idtfcd_cJ0
!   filtro=0  !esta variável direciona o uso de filtros J0 e J1 em vez de seno e cosseno
! !!  call identfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0)
! !!  write(*,*)'Entre com o criador dos filtros J1: Rijo(0), Frayzer(1), Guptasarma(2), Kong(3) ou Key(4)'
! !!  read(*,*)idtfcd_cJ1
! !!  call identfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1)
!
!   idtfcd_cJ0 = 3
!   ident_fJ0 = 0
!   nJ0 = 241
!   idtfcd_cJ1 = 3
!   ident_fJ1 = 1
!   nJ1 = 241
!
!   allocate(KrJ0(nJ0),KrJ1(nJ1),w_J0(nJ0),w_J1(nJ1))
!
!   call constfiltro(filtro,idtfcd_cJ0,ident_fJ0,nJ0,r,KrJ0,w_J0)
!   call constfiltro(filtro,idtfcd_cJ1,ident_fJ1,nJ1,r,KrJ1,w_J1)
!
!   if (camad == 0 .and. camadT /= 0) then
!
!     kerEx_J1 = (0.d0,0.d0)
!     kerEy_J1 = (0.d0,0.d0)
!     kerHx_J1 = (0.d0,0.d0)
!     kerHy_J1 = (0.d0,0.d0)
!       do i=1,nJ1
!         k_r=krJ1(i)
!         kerEx_J1 = kerEx_J1 + KEx0_J1(k_r)*w_J1(i)
!         kerEy_J1 = kerEy_J1 + KEy0_J1(k_r)*w_J1(i)
!         kerHx_J1 = kerHx_J1 + KHx0_J1(k_r)*w_J1(i)
!         kerHy_J1 = kerHy_J1 + KHy0_J1(k_r)*w_J1(i)
!       end do
!     Ex_p = -kerEx_J1 / (2.d0*pi*r)
!     Ey_p = -kerEy_J1 / (2.d0*pi*r)
!     Hx_p = -kerHx_J1 / (2.d0*pi*r)
!     Hy_p = kerHy_J1 / (2.d0*pi*r)
!
!     kerEz_J0 = (0.d0,0.d0)
!       do j=i,nJ0
!         k_r=krJ0(j)
!         kerEz_J0 = kerEz_J0 + KEz0_J0(k_r)*w_J0(j)
!       end do
!     Ez_p = kerEz_J0 / (2.d0*pi*r)
!
!   elseif (camad < camadT) then !camada k
!
!     kerEx_J1 = (0.d0,0.d0)
!     kerEy_J1 = (0.d0,0.d0)
!     kerHx_J1 = (0.d0,0.d0)
!     kerHy_J1 = (0.d0,0.d0)
!       do i=1,nJ1
!         k_r=krJ1(i)
!         kerEx_J1 = kerEx_J1 + KExk_J1(k_r)*w_J1(i)
!         kerEy_J1 = kerEy_J1 + KEyk_J1(k_r)*w_J1(i)
!         kerHx_J1 = kerHx_J1 + KHxk_J1(k_r)*w_J1(i)
!         kerHy_J1 = kerHy_J1 + KHyk_J1(k_r)*w_J1(i)
!       end do
!     Ex_p = -kerEx_J1 / (2.d0*pi*r)
!     Ey_p = -kerEy_J1 / (2.d0*pi*r)
!     Hx_p = -kerHx_J1 / (2.d0*pi*r)
!     Hy_p = kerHy_J1 / (2.d0*pi*r)
!
!     kerEz_J0 = (0.d0,0.d0)
!       do j=1,nJ0
!         k_r=krJ0(j)
!         kerEz_J0 = kerEz_J0 + KEzk_J0(k_r)*w_J0(j)
!       end do
!     Ez_p = kerEz_J0 / (2.d0*pi*r)
!
!   elseif (camad == camadT .and. z <= h0) then !na mesma camada do transmissor mas acima dele
!
!     kerEx_J1 = (0.d0,0.d0)
!     kerEy_J1 = (0.d0,0.d0)
!     kerHx_J1 = (0.d0,0.d0)
!     kerHy_J1 = (0.d0,0.d0)
!       do i=1,nJ1
!         k_r=krJ1(i)
!         kerEx_J1 = kerEx_J1 + KExl_up_J1(k_r)*w_J1(i)
!         kerEy_J1 = kerEy_J1 + KEyl_up_J1(k_r)*w_J1(i)
!         kerHx_J1 = kerHx_J1 + KHxl_up_J1(k_r)*w_J1(i)
!         kerHy_J1 = kerHy_J1 + KHyl_up_J1(k_r)*w_J1(i)
!       end do
!     Ex_p = -kerEx_J1 / (2.d0*pi*r)
!     Ey_p = -kerEy_J1 / (2.d0*pi*r)
!     Hx_p = -kerHx_J1 / (2.d0*pi*r)
!     Hy_p = kerHy_J1 / (2.d0*pi*r)
!
!     kerEz_J0 = (0.d0,0.d0)
!       do j=1,nJ0
!         k_r=krJ0(j)
!         kerEz_J0 = kerEz_J0 + KEzl_up_J0(k_r)*w_J0(j)
!       end do
!     Ez_p = kerEz_J0 / (2.d0*pi*r)
!
!   elseif (camad == camadT .and. z > h0) then  !na mesma camada do transmissor mas abaixo dele
!
!     kerEx_J1 = (0.d0,0.d0)
!     kerEy_J1 = (0.d0,0.d0)
!     kerHx_J1 = (0.d0,0.d0)
!     kerHy_J1 = (0.d0,0.d0)
!       do i=1,nJ1
!         k_r=krJ1(i)
!         kerEx_J1 = kerEx_J1 + KExl_dw_J1(k_r)*w_J1(i)
!         kerEy_J1 = kerEy_J1 + KEyl_dw_J1(k_r)*w_J1(i)
!         kerHx_J1 = kerHx_J1 + KHxl_dw_J1(k_r)*w_J1(i)
!         kerHy_J1 = kerHy_J1 + KHyl_dw_J1(k_r)*w_J1(i)
!       end do
!     Ex_p = kerEx_J1 / (2.d0*pi*r)
!     Ey_p = kerEy_J1 / (2.d0*pi*r)
!     Hx_p = -kerHx_J1 / (2.d0*pi*r)
!     Hy_p = kerHy_J1 / (2.d0*pi*r)
!
!     kerEz_J0 = (0.d0,0.d0)
!       do j=1,nJ0
!         k_r=krJ0(j)
!         kerEz_J0 = kerEz_J0 + KEzl_dw_J0(k_r)*w_J0(j)
!       end do
!     Ez_p = kerEz_J0 / (2.d0*pi*r)
!
!   elseif (camad > camadT .and. camad /= n) then !camada j
!
!     kerEx_J1 = (0.d0,0.d0)
!     kerEy_J1 = (0.d0,0.d0)
!     kerHx_J1 = (0.d0,0.d0)
!     kerHy_J1 = (0.d0,0.d0)
!       do i=1,nJ1
!         k_r=krJ1(i)
!         kerEx_J1 = kerEx_J1 + KExj_J1(k_r)*w_J1(i)
!         kerEy_J1 = kerEy_J1 + KEyj_J1(k_r)*w_J1(i)
!         kerHx_J1 = kerHx_J1 + KHxj_J1(k_r)*w_J1(i)
!         kerHy_J1 = kerHy_J1 + KHyj_J1(k_r)*w_J1(i)
!       end do
!     Ex_p = kerEx_J1 / (2.d0*pi*r)
!     Ey_p = kerEy_J1 / (2.d0*pi*r)
!     Hx_p = -kerHx_J1 / (2.d0*pi*r)
!     Hy_p = kerHy_J1 / (2.d0*pi*r)
!
!     kerEz_J0 = (0.d0,0.d0)
!       do j=1,nJ0
!         k_r=krJ0(j)
!         kerEz_J0 = kerEz_J0 + KEzj_J0(k_r)*w_J0(j)
!       end do
!     Ez_p = kerEz_J0 / (2.d0*pi*r)
!
!   elseif (camad == n .and. camadT /= n) then  !camada n
!
!     kerEx_J1 = (0.d0,0.d0)
!     kerEy_J1 = (0.d0,0.d0)
!     kerHx_J1 = (0.d0,0.d0)
!     kerHy_J1 = (0.d0,0.d0)
!       do i=1,nJ1
!         k_r=krJ1(i)
!         kerEx_J1 = kerEx_J1 + KExn_J1(k_r)*w_J1(i)
!         kerEy_J1 = kerEy_J1 + KEyn_J1(k_r)*w_J1(i)
!         kerHx_J1 = kerHx_J1 + KHxn_J1(k_r)*w_J1(i)
!         kerHy_J1 = kerHy_J1 + KHyn_J1(k_r)*w_J1(i)
!       end do
!     Ex_p = kerEx_J1 / (2.d0*pi*r)
!     Ey_p = kerEy_J1 / (2.d0*pi*r)
!     Hx_p = -kerHx_J1 / (2.d0*pi*r)
!     Hy_p = kerHy_J1 / (2.d0*pi*r)
!
!     kerEz_J0 = (0.d0,0.d0)
!       do j=1,nJ0
!         k_r=krJ0(j)
!         kerEz_J0 = kerEz_J0 + KEzn_J0(k_r)*w_J0(j)
!       end do
!     Ez_p = kerEz_J0 / (2.d0*pi*r)
!
!   end if
! end if
!
! contains
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function wvnb_2(cam)
!   implicit none
!   integer,intent(in)::cam
!   complex(dp)::wvnb_2
!
!   if (cam == 0) then
!   wvnb_2 = -zeta*neta !consideramos aqui a contribuição apenas da parte imaginária da admitividade, já que sigma é zero.
!   else
!   wvnb_2 = -zeta*condut(cam)
!   end if
!
!   return
!   end function wvnb_2
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function u(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::u
! ! Não mude a natureza. Use essa fórmula, mesmo no ar em EI. Mas, talvez seja bom modificar em EF.
!   u=sqrt(kr**2.d0-wvnb_2(cam))
!
!   return
!   end function u
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function ImpInt(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::ImpInt
!
!   if (cam /= 0) then
!   ImpInt=u(kr,cam)/condut(cam)
!   else
!   ImpInt=u(kr,cam)/neta
!   end if
!
!   return
!   end function ImpInt
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function uh(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::uh
! ! nao esquecer de construir h(0)=0
!   uh=u(kr,cam)*h(cam)
!
!   return
!   end function uh
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function tgh(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::tgh
!
!   tgh = (1.d0 - exp(-2.d0 * uh(kr,cam))) / (1.d0 + exp(-2.d0 * uh(kr,cam)))
!
!   return
!   end function tgh
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Imp_Ap_dw(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   integer::j
!   complex(dp)::Imp_Ap_dw
!
!   Imp_Ap_dw = ImpInt(kr,n)
!   do j=n-1,cam,-1
!   Imp_Ap_dw = ImpInt(kr,j) * (Imp_Ap_dw + ImpInt(kr,j) * &
!       tgh(kr,j)) / (ImpInt(kr,j) + Imp_Ap_dw * tgh(kr,j))
!   end do
!
!   return
!   end function Imp_Ap_dw
! !   recursive function Imp_Ap_dw(kr,cam) result(ImpAp_dw)
! !   implicit none
! !   real(dp),intent(in)::kr
! !   integer,intent(in)::cam
! !   complex(dp)::ImpAp_dw
! !
! !   if (cam == n) then
! !   ImpAp_dw = ImpInt(kr,n)
! !   else
! ! ! do j=ncam-1,1,-1
! !   ImpAp_dw = ImpInt(kr,cam) * (Imp_Ap_dw(kr,cam+1) + ImpInt(kr,cam) * &
! !       tgh(kr,cam)) / (ImpInt(kr,cam) + Imp_Ap_dw(kr,cam+1) * tgh(kr,cam))
! !   end if
! !
! !   return
! !   end function Imp_Ap_dw
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Imp_Ap_up(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   integer::j
!   complex(dp)::Imp_Ap_up
!
!   Imp_Ap_up = ImpInt(kr,0)
!   do j=1,cam
!   Imp_Ap_up = ImpInt(kr,j) * (Imp_Ap_up + ImpInt(kr,j) * &
!       tgh(kr,j)) / (ImpInt(kr,j) + Imp_Ap_up * tgh(kr,j))
!   end do
!
!   return
!   end function Imp_Ap_up
! !   recursive function Imp_Ap_up(kr,cam) result(ImpAp_up)
! !   implicit none
! !   real(dp),intent(in)::kr
! !   integer,intent(in)::cam
! !   complex(dp)::ImpAp_up
! !
! !   if (cam == 0) then
! !   ImpAp_up = ImpInt(kr,0)
! !   else
! ! ! do j=1,ncam
! !   ImpAp_up = ImpInt(kr,cam) * (Imp_Ap_up(kr,cam-1) + ImpInt(kr,cam) * &
! !       tgh(kr,cam)) / (ImpInt(kr,cam) + Imp_Ap_up(kr,cam-1) * tgh(kr,cam))
! !   end if
! !
! !   return
! !   end function Imp_Ap_up
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function RTM_dw(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::RTM_dw
!
!   if (cam==n) then
!   RTM_dw = (0.d0,0.d0)
!   else
!   RTM_dw = (ImpInt(kr,cam) - Imp_Ap_dw(kr,cam+1)) / (ImpInt(kr,cam) + Imp_Ap_dw(kr,cam+1))
!   end if
!
!   return
!   end function RTM_dw
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function RTM_up(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::RTM_up
!
!   if (cam == 0) then
!   RTM_up = (0.d0,0.d0)
!   else
!   RTM_up = (ImpInt(kr,cam) - Imp_Ap_up(kr,cam-1)) / (ImpInt(kr,cam) + Imp_Ap_up(kr,cam-1))
!   end if
!
!   return
!   end function RTM_up
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function AM_dw(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::AM_dw
!
!   AM_dw = (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*exp(u(kr,camadT)*(prof(camadT-1)-h(camadT)-h0))) / &
!       (1.d0 - RTM_up(kr,camadT)*RTM_dw(kr,camadT)*exp(-2.d0*uh(kr,camadT)))
!
!   return
!   end function AM_dw
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function AM_up(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::AM_up
!
!   AM_up = (exp(u(kr,camadT)*(prof(camadT-1)-h0)) + RTM_dw(kr,camadT)*exp(-u(kr,camadT)*(prof(camadT)+h(camadT)-h0))) / &
!       (1.d0 - RTM_up(kr,camadT)*RTM_dw(kr,camadT)*exp(-2.d0*uh(kr,camadT)))
!
!   return
!   end function AM_up
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function TM_dw(Kr,cam)
!   implicit none
!   real(dp),intent(in)::Kr
!   integer,intent(in)::cam
!   integer::j
!   complex(dp)::TM_dw
!
!   TM_dw = Iw * dsz / (2.d0*u(kr,camadT))
!   do j=camadT+1,cam
!     if (j == (camadT + 1) .and. j == n) then
!     TM_dw = TM_dw * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
!         RTM_dw(kr,camadT)*AM_dw(kr))
!     elseif (j == (camadT + 1) .and. j /= n) then
!     TM_dw = TM_dw * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
!         RTM_dw(kr,camadT)*AM_dw(kr)) / (1.d0 + RTM_dw(kr,camadT+1)*exp(-2.d0*uh(kr,camadT+1)))
!     elseif (j /= n) then
!     TM_dw = TM_dw * (1.d0 + RTM_dw(kr,j-1)) * exp(-uh(kr,j-1)) / (1.d0 + RTM_dw(kr,j) * exp(-2.d0*uh(kr,j)))
!     else
!     TM_dw = TM_dw * (1.d0 + RTM_dw(kr,n-1)) * exp(-uh(kr,n-1))
!     end if
!   end do
!
!   return
!   end function
! !   recursive function TM_dw(Kr,cam) result(TMdw)
! !   implicit none
! !   real(dp),intent(in)::Kr
! !   integer,intent(in)::cam
! !   complex(dp)::TMdw
! !
! !   if (cam == camadT) then
! !   TMdw = Iw * dsz / (2.d0*u(kr,camadT))
! !   elseif (cam == (camadT + 1) .and. cam == n) then
! !   TMdw = TM_dw(kr,cam-1) * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
! !       RTM_dw(kr,camadT)*AM_dw(kr))
! !   elseif (cam == (camadT + 1) .and. cam /= n) then
! !   TMdw = TM_dw(kr,cam-1) * (exp(-u(kr,camadT)*(prof(camadT)-h0)) + RTM_up(kr,camadT)*AM_up(kr)*exp(-uh(kr,camadT)) + &
! !       RTM_dw(kr,camadT)*AM_dw(kr)) / (1.d0 + RTM_dw(kr,camadT+1)*exp(-2.d0*uh(kr,camadT+1)))
! !   elseif (cam /= n) then
! !   TMdw = TM_dw(kr,cam-1) * (1.d0 + RTM_dw(kr,cam-1)) * exp(-uh(kr,cam-1)) / (1.d0 + RTM_dw(kr,cam) * exp(-2.d0*uh(kr,cam)))
! !   else
! !   TMdw = TM_dw(kr,n-1) * (1.d0 + RTM_dw(kr,n-1)) * exp(-uh(kr,n-1))
! !   end if
! !
! !   return
! !   end function TM_dw
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function TM_up(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   integer::j
!   complex(dp)::TM_up
!
!   TM_up = Iw * dsz / (2.d0*u(kr,camadT))
!   do j=camadT-1,cam,-1
!     if (j == (camadT - 1) .and. j == 0) then
!     TM_up = TM_up * (exp(-u(kr,camadT)*h0) + RTM_up(kr,camadT)*AM_up(kr) + RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT)))
!     elseif (j == (camadT - 1) .and. j /= 0) then
!     TM_up = TM_up * (exp(u(kr,camadT)*(prof(camadT-1)-h0)) + RTM_up(kr,camadT)*AM_up(kr) + &
!         RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT))) / (1.d0 + RTM_up(kr,camadT-1)*exp(-2.d0*uh(kr,camadT-1)))
!     elseif (j /= 0) then
!     TM_up = TM_up * (1.d0 + RTM_up(kr,j+1)) * exp(-uh(kr,j+1)) / (1.d0 + RTM_up(kr,j) * exp(-2.d0*uh(kr,j)))
!     else
!     TM_up = TM_up * (1.d0 + RTM_up(kr,1)) * exp(-uh(kr,1))
!     end if
!   end do
!
!   return
!   end function TM_up
! !   recursive function TM_up(kr,cam) result(TMup)
! !   implicit none
! !   real(dp),intent(in)::kr
! !   integer,intent(in)::cam
! !   complex(dp)::TMup
! !
! !   if (cam == camadT) then
! !   TMup = Iw * dsz / 2.d0
! !   elseif (cam == (camadT - 1) .and. cam == 0) then
! !   TMup = TM_up(kr,cam+1) * (exp(-u(kr,camadT)*h0) + RTM_up(kr,camadT)*AM_up(kr) + RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT)))
! !   elseif (cam == (camadT - 1) .and. cam /= 0) then
! !   TMup = TM_up(kr,cam+1) * (exp(u(kr,camadT)*(prof(camadT-1)-h0)) + RTM_up(kr,camadT)*AM_up(kr) + &
! !       RTM_dw(kr,camadT)*AM_dw(kr)*exp(-uh(kr,camadT))) / (1.d0 + RTM_up(kr,camadT-1)*exp(-2.d0*uh(kr,camadT-1)))
! !   elseif (cam /= 0) then
! !   TMup = TM_up(kr,cam+1) * (1.d0 + RTM_up(kr,cam+1)) * exp(-uh(kr,cam+1)) / (1.d0 + RTM_up(kr,cam) * exp(-2.d0*uh(kr,cam)))
! !   else
! !   TMup = TM_up(kr,1) * (1.d0 + RTM_up(kr,1)) * exp(-uh(kr,1))
! !   end if
! !
! !   return
! !   end function TM_up
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktmdz_0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktmdz_0
!
!   Ktmdz_0 = ImpInt(kr,0) * TM_up(kr,0) * exp(u(kr,0) * z)
!
!   return
!   end function Ktmdz_0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktmdz_k(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::Ktmdz_k
!
!   ktmdz_k = ImpInt(kr,cam) * TM_up(kr,cam) * (exp(u(kr,cam) * (z-prof(cam))) - &
!          RTM_up(kr,cam) * exp(-u(kr,cam) * (z-prof(cam-1)+h(cam))))
!
!   return
!   end function Ktmdz_k
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktmdz_l_up(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktmdz_l_up
!
!   Ktmdz_l_up = ImpInt(kr,camadT) * TM_up(kr,camadT) * (exp(u(kr,camadT) * (z-h0)) - &
!         RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) + &
!         RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))
!
!   return
!   end function Ktmdz_l_up
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktmdz_l_dw(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktmdz_l_dw
!
!   Ktmdz_l_dw = ImpInt(kr,camadT) * TM_dw(kr,camadT) * (exp(-u(kr,camadT) * (z-h0)) + &
!     RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) - &
!     RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))
!
!   return
!   end function Ktmdz_l_dw
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktmdz_j(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::Ktmdz_j
!
!   Ktmdz_j = ImpInt(kr,cam) * TM_dw(kr,cam) * (exp(-u(kr,cam) * (z-prof(cam-1))) - &
!       RTM_dw(kr,cam) * exp(u(kr,cam) * (z-prof(cam)-h(cam))))
!
!   return
!   end function Ktmdz_j
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktmdz_n(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktmdz_n
!
!   Ktmdz_n = ImpInt(kr,n) * TM_dw(kr,n) * exp(-u(kr,n) * (z-prof(n-1)))
!
!   return
!   end function Ktmdz_n
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktm_0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktm_0
!
!   Ktm_0 = TM_up(kr,0) * exp(u(kr,0) * z)
!
!   return
!   end function Ktm_0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktm_k(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::Ktm_k
!
!   ktm_k = TM_up(kr,cam) * (exp(u(kr,cam) * (z-prof(cam))) + &
!          RTM_up(kr,cam) * exp(-u(kr,cam) * (z-prof(cam-1)+h(cam))))
!
!   return
!   end function Ktm_k
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktm_l_up(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktm_l_up
!
!   Ktm_l_up = TM_up(kr,camadT) * (exp(u(kr,camadT) * (z-h0)) + &
!         RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) + &
!         RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))
!
!   return
!   end function Ktm_l_up
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktm_l_dw(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktm_l_dw
!
!   Ktm_l_dw = TM_dw(kr,camadT) * (exp(-u(kr,camadT) * (z-h0)) + &
!     RTM_up(kr,camadT) * AM_up(kr) * exp(-u(kr,camadT) * (z-prof(camadT-1))) + &
!     RTM_dw(kr,camadT) * AM_dw(kr) * exp(u(kr,camadT) * (z-prof(camadT))))
!
!   return
!   end function Ktm_l_dw
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktm_j(kr,cam)
!   implicit none
!   real(dp),intent(in)::kr
!   integer,intent(in)::cam
!   complex(dp)::Ktm_j
!
!   Ktm_j = TM_dw(kr,cam) * (exp(-u(kr,cam) * (z-prof(cam-1))) + &
!       RTM_dw(kr,cam) * exp(u(kr,cam) * (z-prof(cam)-h(cam))))
!
!   return
!   end function Ktm_j
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function Ktm_n(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::Ktm_n
!
!   Ktm_n = TM_dw(kr,n) * exp(-u(kr,n) * (z-prof(n-1)))
!
!   return
!   end function Ktm_n
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEx0_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEx0_J1
!
!   KEx0_J1 = x/r * Ktmdz_0(kr) * kr * kr
!
!   return
!   end function KEx0_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KExk_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KExk_J1
!
!   KExk_J1 = x/r * Ktmdz_k(kr,camad) * kr * kr
!
!   return
!   end function KExk_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KExl_up_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KExl_up_J1
!
!   KExl_up_J1 = x/r * Ktmdz_l_up(kr) * kr * kr
!
!   return
!   end function KExl_up_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KExl_dw_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KExl_dw_J1
!
!   KExl_dw_J1 = x/r * Ktmdz_l_dw(kr) * kr * kr
!
!   return
!   end function KExl_dw_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KExj_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KExj_J1
!
!   KExj_J1 = x/r * Ktmdz_j(kr,camad) * kr * kr
!
!   return
!   end function KExj_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KExn_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KExn_J1
!
!   KExn_J1 = x/r * Ktmdz_n(kr) * kr * kr
!
!   return
!   end function KExn_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEy0_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEy0_J1
!
!   KEy0_J1 = y/r * Ktmdz_0(kr) * kr * kr
!
!   return
!   end function KEy0_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEyk_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEyk_J1
!
!   KEyk_J1 = y/r * Ktmdz_k(kr,camad) * kr * kr
!
!   return
!   end function KEyk_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEyl_up_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEyl_up_J1
!
!   KEyl_up_J1 = y/r * Ktmdz_l_up(kr) * kr * kr
!
!   return
!   end function KEyl_up_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEyl_dw_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEyl_dw_J1
!
!   KEyl_dw_J1 = y/r * Ktmdz_l_dw(kr) * kr * kr
!
!   return
!   end function KEyl_dw_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEyj_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEyj_J1
!
!   KEyj_J1 = y/r * Ktmdz_j(kr,camad) * kr * kr
!
!   return
!   end function KEyj_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEyn_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEyn_J1
!
!   KEyn_J1 = y/r * Ktmdz_n(kr) * kr * kr
!
!   return
!   end function KEyn_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEz0_J0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEz0_J0
!
!   KEz0_J0 = Ktm_0(kr) * kr * kr * kr / neta
!
!   return
!   end function KEz0_J0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEzk_J0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEzk_J0
!
!   KEzk_J0 = Ktm_k(kr,camad) * kr * kr * kr / condut(camad)
!
!   return
!   end function KEzk_J0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEzl_up_J0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEzl_up_J0
!
!   if (camad /= 0) then
!   KEzl_up_J0 = Ktm_l_up(kr) * kr * kr * kr / condut(camad)
!   else
!   KEzl_up_J0 = Ktm_l_up(kr) * kr * kr * kr / neta
!   end if
!
!   return
!   end function KEzl_up_J0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEzl_dw_J0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEzl_dw_J0
!
!   if (camad /= 0) then
!   KEzl_dw_J0 = Ktm_l_dw(kr) * kr * kr * kr / condut(camad)
!   else
!   KEzl_dw_J0 = Ktm_l_dw(kr) * kr * kr * kr / neta
!   end if
!
!   return
!   end function KEzl_dw_J0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEzj_J0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEzj_J0
!
!   KEzj_J0 = Ktm_j(kr,camad) * kr * kr * kr / condut(camad)
!
!   return
!   end function KEzj_J0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KEzn_J0(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KEzn_J0
!
!   KEzn_J0 = Ktm_n(kr) * kr * kr * kr / condut(n)
!
!   return
!   end function KEzn_J0
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHx0_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHx0_J1
!
!   KHx0_J1 = y/r * Ktm_0(kr) * kr * kr
!
!   return
!   end function KHx0_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHxk_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHxk_J1
!
!   KHxk_J1 = y/r * Ktm_k(kr,camad) * kr * kr
!
!   return
!   end function KHxk_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHxl_up_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHxl_up_J1
!
!   KHxl_up_J1 = y/r * Ktm_l_up(kr) * kr * kr
!
!   return
!   end function KHxl_up_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHxl_dw_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHxl_dw_J1
!
!   KHxl_dw_J1 = y/r * Ktm_l_dw(kr) * kr * kr
!
!   return
!   end function KHxl_dw_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHxj_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHxj_J1
!
!   KHxj_J1 = y/r * Ktm_j(kr,camad) * kr * kr
!
!   return
!   end function KHxj_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHxn_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHxn_J1
!
!   KHxn_J1 = y/r * Ktm_n(kr) * kr * kr
!
!   return
!   end function KHxn_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHy0_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHy0_J1
!
!   KHy0_J1 = x/r * Ktm_0(kr) * kr * kr
!
!   return
!   end function KHy0_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHyk_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHyk_J1
!
!   KHyk_J1 = x/r * Ktm_k(kr,camad) * kr * kr
!
!   return
!   end function KHyk_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHyl_up_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHyl_up_J1
!
!   KHyl_up_J1 = x/r * Ktm_l_up(kr) * kr * kr
!
!   return
!   end function KHyl_up_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHyl_dw_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHyl_dw_J1
!
!   KHyl_dw_J1 = x/r * Ktm_l_dw(kr) * kr * kr
!
!   return
!   end function KHyl_dw_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHyj_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHyj_J1
!
!   KHyj_J1 = x/r * Ktm_j(kr,camad) * kr * kr
!
!   return
!   end function KHyj_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   function KHyn_J1(kr)
!   implicit none
!   real(dp),intent(in)::kr
!   complex(dp)::KHyn_J1
!
!   KHyn_J1 = x/r * Ktm_n(kr) * kr * kr
!
!   return
!   end function KHyn_J1
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! !-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!   end subroutine ved_xyz
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end module ved
