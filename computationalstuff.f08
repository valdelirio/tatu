module computational_stuff
  !     Modulo com tarefas computacionais bastantes corriqueiras.
  !     Autor: Valdelirio
  !     Ultima atualizacao: 17/10/2016
  implicit none
  save
  ! integer, parameter :: spi=selected_int_kind(3)
  ! integer, parameter :: dpi=selected_int_kind(9)
  integer, parameter :: sp = kind(1.e0)
  integer, parameter :: dp = kind(1.d0)
  integer, parameter :: qp = selected_real_kind(30) !
  real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
  real(dp), parameter :: mu = 4.d-7 * pi
  real(dp), parameter :: epsilon = 8.85d-12
  real(dp), parameter :: Iw = 1.d0
  real(dp), parameter :: dsx = 1.d0
  real(dp), parameter :: dsy = 1.d0
  real(dp), parameter :: dsz = 1.d0
  real(dp), parameter :: mx = 1.d0
  real(dp), parameter :: my = 1.d0
  real(dp), parameter :: mz = 1.d0
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
contains
!=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  !     Altera, se necessario, a ordem das colunas informadas em ja, em cada linha, para ordem
  !     crescente e, reordena a matriz Kg (no formato CSR) conforme a modificacao de ja
  subroutine OrdCres_ja( n, vetin, vetord, M )
    !   Autor: Valdelirio
    implicit none
    integer, intent(in) :: n
    integer, dimension(n), intent(in) :: vetin
    integer, dimension(n), intent(out) :: vetord
    real(dp), dimension(n), intent(inout) :: M

    integer :: i, j, mini, maxi, minvalatual, maxvalatual, ind(n), old(n)

    ind = (/(i,i=1,n)/)
    old = ind
    vetord = vetin

    do i = 2, n
       do j = i, n
          if ( vetord(i-1) < vetord(j) ) THEN
             mini = i - 1
             maxi = j
             minvalatual = vetord(mini)
             maxvalatual = vetord(maxi)
          else
             mini = j
             maxi = i - 1
             minvalatual = vetord(mini)
             maxvalatual = vetord(maxi)
          end if
          vetord(i-1) = minvalatual
          vetord(j) = maxvalatual
          ind(j) = old(maxi)
          ind(i-1) = old(mini)
          old = ind
       end do
    end do
    !   print*,ind
    !   print*,vetord
    !   read(*,*)
    M = M(ind)

  end subroutine OrdCres_ja
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  !     Cria um vetor com ordem crescente de valores a partir de um vetor qualquer, mas sem eliminar os possíveis valores repetidos
  subroutine OrdCresVal( n, vetin, vetord )
    !   Autor: Valdelirio
    implicit none
    integer, intent(in) :: n
    real(dp), intent(in) :: vetin(n)
    real(dp), intent(out) :: vetord(n)

    integer :: i, j
    real(dp) :: minvalatual

    vetord = vetin

    do i = 2, n
       do j = i, n
          minvalatual = dmin1( vetord(i - 1), vetord(j) )
          vetord(j) = dmax1( vetord(i - 1), vetord(j) )
          vetord(i - 1) = minvalatual
       end do
    end do

  end subroutine OrdCresVal
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  !   Redimensiona as coordenadas de um vetor no caso da existência de valores repetidos
  subroutine redim_coord( ncoord_tmp, coord_tmp, ncoord_def, indret )
    !   Autor: Valdelirio
    implicit none
    integer, intent(in) :: ncoord_tmp
    real(dp), intent(in) :: coord_tmp(ncoord_tmp)
    integer, intent(out) :: ncoord_def, indret(ncoord_tmp)

    real(dp), parameter :: eps = 1.d-7
    integer :: i, j, indrep( ncoord_tmp, ncoord_tmp )
    integer :: vetmax( ncoord_tmp, ncoord_tmp ), indmax( ncoord_tmp )

    indrep = 0
    vetmax = 0

    do i = 1, ncoord_tmp
       do j = i, ncoord_tmp
          if (i == j) THEN
             indrep(i,j) = 0
          elseif ( abs(coord_tmp(i) - coord_tmp(j)) < eps ) THEN
             indrep(i,j) = 1
          else
             indrep(i,j) = 0
          end if
       end do
       vetmax(i,:) = MAXLOC( indrep(i,:) )
    end do

    j = 0
    do i = 1, ncoord_tmp
       indmax(i) = vetmax(i,1)
       if ( indmax(i) /= 1 ) THEN
          indret(i) = indmax(i)
          j = j + 1
       else
          indret(i) = 0
       end if
       !        Há um erro horrorendo no uso de vetmax, mas que dá respostas corretas no primeiro elemento de cada linha, daí o código
       !        faz a coisa certa. Devo indubitavelmente verificar que forma devo endereçar as respostas de maxloc, já que no manual diz que é um
       !        array com o índice do maior elemento numa dada dimensão, no entanto vetmax não apresenta tais índices repetidos em seus elementos.
       !print*,vetmax(i,1)
       !pause
    end do

    ncoord_def = ncoord_tmp - j

  end subroutine redim_coord
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  !   Mediante um vetor com possíveis valores repetidos determina-se o vetor com apenas valores distintos. É dependente da subrotina
  !   redim_coord que informa os índices de localização dos valores repetidos
  subroutine renov_coord( ncoord_tmp, coord_tmp, ncoord_def, indret, coord_def )
    !   Autor: Valdelirio
    implicit none
    integer, intent(in) :: ncoord_tmp, ncoord_def
    integer, intent(in) :: indret( ncoord_tmp )
    real(dp), intent(in) :: coord_tmp( ncoord_tmp )
    real(dp), intent(out) :: coord_def( ncoord_def )

    integer(4) :: i, j

    j = 0
    do i = 1, ncoord_tmp
       if ( indret(i) == 0 ) THEN
          j = j + 1
          coord_def(j) = coord_tmp(i)
       else
          j = j
       end if
    end do

  end subroutine renov_coord
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine cont_nos_dist( nnos_tmp, nos_tmp, nnos_def )
    !   contabiliza o numero de inteiros distintos dentro de um vetor
    implicit none
    integer, intent(in) :: nnos_tmp, nos_tmp(:)
    integer, intent(out) :: nnos_def

    integer :: i, j
    real(dp), parameter :: epsilon = 5.d-1

    j = 1
    do i = 2, nnos_tmp
      if ( ABS( nos_tmp(i) - nos_tmp(i - 1) ) > epsilon ) THEN
        j = j + 1
      end if
    end do

    nnos_def = j

  end subroutine cont_nos_dist
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine nos_dist( nnos_tmp, nos_tmp, nnos_def, nos_def )
    !   retorna um vetor com apenas inteiros distintos de um vetor de inteiros
    implicit none
    integer, intent(in) :: nnos_tmp, nnos_def, nos_tmp(:)
    integer, intent(out) :: nos_def(nnos_def)

    integer :: i, j, k
    real(dp), parameter :: epsilon = 5.d-1


    j = 2
    nos_def(1) = nos_tmp(1)
    do k = 2, nnos_def
       loop_b:  do i = j, nnos_tmp
          if ( ABS( nos_tmp(i) - nos_tmp(i - 1) ) > epsilon ) THEN
             nos_def(k) = nos_tmp(i)
             j = i + 1
             EXIT loop_b
          end if
       end do loop_b
    end do

  end subroutine nos_dist
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine qualcamada(z,n,h,camada)
    !   Autor: Valdelirio
    implicit none
    integer,intent(in) :: n
    real(dp), intent(in) :: z, h(n-1)
    integer, intent(out) :: camada
    integer :: i
    real(dp) :: prof(0:n-1)

    prof(0) = 0.d0      !So para efeito de facilitacao de calculo
    if ( n > 1 )  THEN
       prof(1) = h(1)
       if ( n > 2 ) THEN
          do i = 2,n-1
             prof(i) = prof(i-1) + h(i)
          end do
       end if
    end if

    if (z < 0.d0) THEN
       camada = 0
    elseif ( z >= prof(n-1) ) THEN
       camada = n
    else
       do i = n-1,1,-1
          if ( z >= prof(i-1) ) THEN
             camada = i
             EXIT
          end if
       end do
    end if

  end subroutine qualcamada
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine info_rcptrs_com_offset( os1, osn, nT, nR, T_x, R_x, nTa, nTx_a, Tx_a, iTxa )
    !   Autor: Valdelirio
    implicit none
    real(dp), intent(in) :: os1, osn
    integer, intent(in) :: nT, nR
    real(dp), intent(in) :: T_x(nT), R_x(nR)

    integer, intent(out) :: nTa
    integer, dimension(:), allocatable, intent(out) :: nTx_a, iTxa
    real(dp), dimension(:), allocatable, intent(out) :: Tx_a

    integer :: i, j, k, contnTa

    ALLOCATE( nTx_a(nT) )
    k = 0                                                               !contador total de transmissores adjuntos (receptores)
    do i = 1,nT
       contnTa = 0                                              !contador dos transmissores adjuntos (receptores) de cada transmissor
       do j = 1,nR
          if ( dabs(R_x(j) - T_x(i)) >= os1 .AND. dabs(R_x(j) - T_x(i)) <= osn ) THEN
             contnTa = contnTa + 1
             k = k + 1
          end if
       end do
       nTx_a(i) = contnTa                               !registro do numero de transmissores adjuntos (receptores) para o transmissor i
    end do
    nTa = sum(nTx_a)                                    !numero total de transmissores adjuntos (receptores) para todos os transmissores
    ALLOCATE( Tx_a(nTa), iTxa(nTa) )    !vetor para registrar os indices e quais receptores serao de dados para cada transmissor

    k = 1
    do i = 1,nT
       do j = 1,nR
          if ( dabs(R_x(j) - T_x(i)) >= os1 .AND. dabs(R_x(j) - T_x(i)) <= osn ) THEN
             iTxa(k) = j
             Tx_a(k) = R_x(j)
             k = k + 1
          end if
       end do
    end do

  end subroutine info_rcptrs_com_offset
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine indcs_receptoresnomar( nf, nT, nR, nTa, nTx_a, iTxa, indrec )
    !   Autor: Valdelirio
    implicit none
    integer, intent(in) :: nf, nT, nR, nTa
    integer, intent(in) :: nTx_a(nT), iTxa(nTa)
    integer, dimension(:), allocatable, intent(out) :: indrec

    integer :: i, j, k, l, m

    ALLOCATE(indrec(nTa*nf))
    k = 1
    do i = 1,nf
       l = 0
       m = 1
       do j = 1,nT
          indrec(k:k+nTx_a(j)-1) = l + iTxa(m:m+nTx_a(j)-1) + (i - 1) * nR * nT
          m = m + nTx_a(j)
          l = l + nR
          k = k + nTx_a(j)
       end do
    end do

  end subroutine indcs_receptoresnomar
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine indcs_camposseparados( nf, nT, nR, nv_des, ind_mar, indvtcs )
    implicit none
    integer, intent(in) :: nf, nT, nR, nv_des
    integer, dimension(:), allocatable, intent(out) :: ind_mar, indvtcs

    integer :: i, j, k, l, m, n, p

    ALLOCATE( ind_mar(nf*nT*nR), indvtcs(nf*nT*nv_des) )

    k = 1
    l = 1
    m = 1
    p = nR + 1
    do i = 1,nf
       do j = 1, nT
          ind_mar(k:k+nR-1) = (/(n,n=l,l+nR-1)/)
          indvtcs(m:m+nv_des-1) = (/(n,n=p,p+nv_des-1)/)
          k = k + nR
          l = l + nR + nv_des
          m = m + nv_des
          p = p + nR + nv_des
       end do
    end do

  end subroutine indcs_camposseparados
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  ! subroutine suavizafase(idtfc,z,phi)
  !   !   Baseado no paper: Abbas, Kattoush. A new recurrent approach for phase unwrapping.
  !   !   International Journal of Applied Science and Engineering.
  !   ! implicit none
  !   integer(4), intent(in) :: idtfc
  !   complex(dp), dimension(:), intent(in) :: z
  !   real(dp), dimension(size(z)), intent(out) :: phi
  !
  !   real(dp), parameter :: pi = 3.141592653589793238462643383279502884197d0
  !   integer(4) :: n, i
  !
  !   n = SIZE(z)
  !   phi = atan2(dimag(z),real(z))
  !   do i = 2,n
  !      if ( dabs(phi(i) - phi(i-1)) > pi .AND. phi(i) < phi(i-1) ) THEN
  !         phi(i:n) = phi(i:n) + 2 * pi
  !      elseif ( dabs(phi(i) - phi(i-1)) > pi .AND. phi(i) > phi(i-1) ) THEN
  !         phi(i:n) = phi(i:n) - 2 * pi
  !      end if
  !   end do
  !   if ( idtfc /= 0 ) THEN
  !      phi = 18.d1 * phi / pi
  !   end if
  !
  ! end subroutine suavizafase
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  subroutine cria_unid_arq( prenumunid, sufnumunid, prenome, sufnome, numarq, filename )
    !   Autor: Valdelirio
    !   Entradas:
    !   prenumunid: prefixo com o numero da unidade que se deseja criar o arquivo
    !   sufnumunid: sufixo com o numero da unidade que se deseja criar o arquivo
    !   prenome: prefixo do nome do arquivo que se deseja criar
    !   sufbnome: sufixo do nome do arquivo que se deseja criar (geralmente .extensao)
    !   Saidas:
    !   numarq: numero da unidade do arquivo. Inteiro com a concatenacao do prenumunid e sufnumunid
    !   filename: nome da unidade do arquivo. String associando o prefnome com numero sufnumunid e a extensao dada em sufnome
    implicit none
    character(len=18), intent(in) :: prenome, sufnome
    integer, intent(in) :: prenumunid, sufnumunid
    integer, intent(out) :: numarq
    character(len=28), intent(out) :: filename

    character(len=18) :: prenarq, sufnarq
    character(len=18) :: strarq

    prenarq = int2str(prenumunid)
    sufnarq = int2str(sufnumunid)
    strarq = ADJUSTL(prenarq//TRIM(ADJUSTL(sufnarq)))
    numarq = str2int(strarq)

    filename = TRIM(prenome)//TRIM(ADJUSTL(sufnarq))//sufnome

  end subroutine cria_unid_arq
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  function int2str(num)
    !   Autor: Valdelirio
    implicit none
    integer, intent(in) :: num
    character(len=18) :: int2str

    write(int2str,'(I18)')num

    return
  end function int2str
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
  function str2int(str)
    !   Autor: Valdelirio
    implicit none
    character(len=18), intent(in) :: str
    integer :: str2int

    READ(str,'(I18)')str2int

    return
  end function str2int
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-132
  function real2str(num)
    !   Autor: Valdelirio
    implicit none
    real(dp), intent(in) :: num
    character(len=32) :: real2str

    write(real2str,'(F32.20)')num

    return
  end function real2str
  !=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=132
end MODULE computational_stuff
