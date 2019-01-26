module escolha_do_filtro
  use computational_stuff
  contains
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine identfiltro( tipofiltro, criador, idtfcd_f, np_f )
  implicit none
  integer, intent(in) :: tipofiltro, criador
  integer, intent(out) :: idtfcd_f, np_f

  integer :: idtfcdf, n_J0, n_J1, n_sen, n_cos

  write(*, *) 'Que filtro será usado, 0 (J0), 1 (J1), 2 (seno) ou 3 (cosseno)?'
  read(*, *) idtfcdf
  idtfcd_f = idtfcdf
! tipofiltro é uma variável que identifica se os filtros a serem usados serão para transformadas de Hankel (0) ou de Fourier (1)
  if ( tipofiltro == 0 ) then !Para uso de filtros J0, J1 ou J2
    if ( criador == 0 ) then !Para o uso dos filtros de Rijo
      if ( idtfcd_f == 0 ) then
        np_f = 61
      elseif ( idtfcd_f == 1 ) then
        np_f = 47
      end if
    elseif ( criador == 1 ) then !Para o uso dos filtros de Frayzer
      if ( idtfcd_f == 0 ) then
        np_f = 37
      elseif ( idtfcd_f == 1 ) then
        np_f = 27
      end if
    elseif ( criador == 2 ) then !Para o uso dos filtros de Guptasarma
      if ( idtfcd_f == 0 ) then
        write(*, *) 'Quantos pontos serão considerados, 61 ou 120 ?'
        read(*, *) n_J0
        np_f = n_J0
      elseif ( idtfcd_f == 1 ) then
        write(*, *) 'Quantos pontos serão considerados, 47 ou 140 ?'
        read(*, *) n_J1
        np_f = n_J1
      end if
    elseif ( criador == 3 ) then !Para o uso dos filtros de Kong
      if ( idtfcd_f == 0 ) then
        write(*, *) 'Quantos pontos serão considerados, 61 ou 241 ?'
        read(*, *) n_J0
        np_f = n_J0
      elseif ( idtfcd_f == 1 ) then
        write(*, *) 'Quantos pontos serão considerados, 61 ou 241 ?'
        read(*, *) n_J1
        np_f = n_J1
      end if
    elseif ( criador == 4 ) then !Para o uso dos filtros de Key
      if ( idtfcd_f == 0 ) then
        write(*, *) 'Quantos pontos serão considerados, 101, 201 ou 401 ?'
        read(*, *) n_J0
        np_f = n_J0
      elseif ( idtfcd_f == 1 ) then
        write(*, *) 'Quantos pontos serão considerados, 101, 201 ou 401 ?'
        read(*, *) n_J1
        np_f = n_J1
      end if
    elseif ( criador == 5 ) then !Para o uso dos filtros de Werthmuller
      if ( idtfcd_f == 0 ) then
        write(*, *) 'Quantos pontos serão considerados 201 ou 2001 ?'
        read(*, *) n_J0
        np_f = n_J0
      elseif ( idtfcd_f == 1 ) then
        write(*, *) 'Quantos pontos serão considerados 201 ou 2001 ?'
        read(*, *) n_J1
        np_f = n_J1
      end if
    end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  elseif ( tipofiltro == 1 ) then !Para uso de filtros seno ou cosseno
    if ( criador == 1 ) then !Para o uso dos filtros seno e cosseno do Frayzer
      if ( idtfcd_f == 2 ) then
        write(*, *) 'Quantos pontos serão considerados, 19, 30 ou 40 ?'
        read(*, *) n_sen
        np_f = n_sen
      elseif ( idtfcd_f == 3 ) then
        write(*, *) 'Quantos pontos serão considerados, 19, 30 ou 40 ?'
        read(*, *) n_cos
        np_f=n_cos
      end if
    elseif ( criador == 4 ) then !Para o uso dos filtros seno e cosseno do Kerry Key
      if ( idtfcd_f == 2 ) then
        write(*, *) 'Quantos pontos serão considerados, 81, 101, 201, 241 ou 601 ?'
        read(*, *) n_sen
        np_f = n_sen
      elseif ( idtfcd_f == 3 ) then
        write(*, *) 'Quantos pontos serão considerados, 81, 101, 201, 241 ou 601 ?'
        read(*, *) n_cos
        np_f = n_cos
      end if
    elseif ( criador == 5 ) then !Para o uso dos filtros seno e cosseno do Werthmuller
      if ( idtfcd_f == 2 ) then
        np_f = 201
      elseif ( idtfcd_f == 3 ) then
        np_f = 201
      end if
    end if
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  else
    print*,'especificação incorreta do tipo de filtro'
    stop
  end if

  end subroutine
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
  subroutine constfiltro( tipofiltro, criador, idtfcd_f, np_f, r, kr, w )
  use filtros_J0_J1_J2_sen_cos

  implicit none
  real(dp), intent(in) :: r
  integer, intent(in) :: tipofiltro, criador, idtfcd_f, np_f
  real(dp), intent(out) :: Kr( np_f ), w( np_f )
  real(dp) :: abs_f( np_f )

  real(dp), dimension(:), allocatable :: lr_abs_J061, lr_pes_J061, lr_abs_J147, lr_pes_J147
  real(dp), dimension(:), allocatable :: fl_abs_J037, fl_pes_J037, fl_abs_J127, fl_pes_J127, fl_abs_J219, fl_pes_J219
  real(dp), dimension(:), allocatable :: gs_abs_J061, gs_pes_J061, gs_abs_J0120, gs_pes_J0120
  real(dp), dimension(:), allocatable :: gs_abs_J147, gs_pes_J147, gs_abs_J1140, gs_pes_J1140
  real(dp), dimension(:), allocatable :: k_abs_J0161, k_pes_J061, k_pes_J161, k_abs_J01241, k_pes_J0241, k_pes_J1241
  real(dp), dimension(:), allocatable :: kk_abs_J01101, kk_pes_J0101, kk_pes_J1101, kk_abs_J01201, kk_pes_J0201, kk_pes_J1201
  real(dp), dimension(:), allocatable :: kk_abs_J01401, kk_pes_J0401, kk_pes_J1401
  real(dp), dimension(:), allocatable :: wm_abs_J01201, wm_pes_J0201, wm_pes_J1201, wm_abs_J012001, wm_pes_J02001, wm_pes_J12001

  real(dp), dimension(:), allocatable :: fl_abs_s19, fl_pes_sen19, fl_abs_s30, fl_pes_sen30, fl_abs_s40, fl_pes_sen40
  real(dp), dimension(:), allocatable :: fl_abs_c19, fl_pes_cos19, fl_abs_c30, fl_pes_cos30, fl_abs_c40, fl_pes_cos40
  real(dp), dimension(:), allocatable :: kk_abs_81, kk_pes_sen81, kk_pes_cos81, kk_abs_241, kk_pes_sen241, kk_pes_cos241
  real(dp), dimension(:), allocatable :: kk_abs_601, kk_pes_sen601, kk_pes_cos601, kk_abs_101, kk_pes_sen101, kk_pes_cos101
  real(dp), dimension(:), allocatable :: kk_abs_201, kk_pes_sen201, kk_pes_cos201
  real(dp), dimension(:), allocatable :: wm_abs_201, wm_pes_sen201, wm_pes_cos201

  real(dp), dimension(:), allocatable :: fl_abs_J070, fl_pes_J070, fl_abs_J170, fl_pes_J170

! tipofiltro é uma variável que identifica se os filtros a serem usados serão para transformadas de Hankel (0) ou de Fourier (1)
if ( tipofiltro == 0 ) then !Para uso de filtros J0, J1 ou J2
  if ( criador == 0 ) then !Para o uso dos filtros de Rijo
    allocate( lr_abs_J061(61), lr_pes_J061(61), lr_abs_J147(47), lr_pes_J147(47) )
    call J0J1Rijo( lr_abs_J061, lr_pes_J061, lr_abs_J147, lr_pes_J147 )
    if ( idtfcd_f == 0 ) then
      abs_f = lr_abs_J061
      w = lr_pes_J061
    elseif ( idtfcd_f == 1 ) then
      abs_f = lr_abs_J147
      w = lr_pes_J147
    end if
    deallocate( lr_abs_J061, lr_pes_J061, lr_abs_J147, lr_pes_J147 )
    kr = dexp(abs_f) / dabs(r)
  elseif ( criador == 1 ) then !Para o uso dos filtros de Frayzer
    allocate( fl_abs_J037(37), fl_pes_J037(37), fl_abs_J127(27), fl_pes_J127(27), fl_abs_J219(19), fl_pes_J219(19) )
    call J0J1J2Frayzer( fl_abs_J037, fl_pes_J037, fl_abs_J127, fl_pes_J127, fl_abs_J219, fl_pes_J219 , &
                        fl_abs_J070, fl_pes_J070, fl_abs_J170, fl_pes_J170 )
    if ( idtfcd_f == 0 ) then
      abs_f = fl_abs_J037
      w = fl_pes_J037
    elseif ( idtfcd_f == 1 ) then
      abs_f = fl_abs_J127
      w = fl_pes_J127
    end if
    deallocate( fl_abs_J037, fl_pes_J037, fl_abs_J127, fl_pes_J127, fl_abs_J219, fl_pes_J219 )
    kr = dexp(abs_f) / dabs(r)
  elseif ( criador == 2 ) then !Para o uso dos filtros de Guptasarma
    allocate( gs_abs_J061(61), gs_pes_J061(61), gs_abs_J0120(120), gs_pes_J0120(120), gs_abs_J147(47), gs_pes_J147(47))
    allocate( gs_abs_J1140(140), gs_pes_J1140(140) )
    call J0J1GS( gs_abs_J061, gs_pes_J061, gs_abs_J0120, gs_pes_J0120, gs_abs_J147, gs_pes_J147, gs_abs_J1140, gs_pes_J1140 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 61 ) then
        abs_f = gs_abs_J061
        w = gs_pes_J061
      elseif ( np_f == 120 ) then
        abs_f = gs_abs_J0120
        w = gs_pes_J0120
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 47 ) then
        abs_f = gs_abs_J147
        w = gs_pes_J147
      elseif ( np_f == 140 ) then
        abs_f = gs_abs_J1140
        w = gs_pes_J1140
      end if
    end if
    deallocate( gs_abs_J061, gs_pes_J061, gs_abs_J0120, gs_pes_J0120, gs_abs_J147, gs_pes_J147, gs_abs_J1140, gs_pes_J1140 )
    kr = ( 10.d0 ** (abs_f) ) / dabs(r)
  elseif ( criador == 3 ) then !Para o uso dos filtros de Kong
    allocate( k_abs_J0161( 61 ), k_pes_J061(61), k_pes_J161(61), k_abs_J01241(241), k_pes_J0241(241), k_pes_J1241(241) )
    call J0J1Kong( k_abs_J0161, k_pes_J061, k_pes_J161, k_abs_J01241, k_pes_J0241, k_pes_J1241 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 61 ) then
        abs_f = k_abs_J0161
        w = k_pes_J061
      elseif ( np_f == 241 ) then
        abs_f = k_abs_J01241
        w = k_pes_J0241
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 61 ) then
        abs_f = k_abs_J0161
        w = k_pes_J161
      elseif ( np_f == 241 ) then
        abs_f = k_abs_J01241
        w = k_pes_J1241
      end if
    end if
    deallocate( k_abs_J0161, k_pes_J061, k_pes_J161, k_abs_J01241, k_pes_J0241, k_pes_J1241 )
    kr = abs_f / dabs(r)
  elseif ( criador == 4 ) then !Para o uso dos filtros de Key
    allocate( kk_abs_J01101(101), kk_pes_J0101(101), kk_pes_J1101(101), kk_abs_J01201(201), kk_pes_J0201(201) )
    allocate( kk_pes_J1201(201), kk_abs_J01401(401), kk_pes_J0401(401), kk_pes_J1401(401) )
    call J0J1Key( kk_abs_J01101, kk_pes_J0101, kk_pes_J1101, kk_abs_J01201, kk_pes_J0201, kk_pes_J1201, &
                  kk_abs_J01401, kk_pes_J0401, kk_pes_J1401 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 101 ) then
        abs_f = kk_abs_J01101
        w = kk_pes_J0101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_J01201
        w = kk_pes_J0201
      elseif ( np_f == 401 ) then
        abs_f = kk_abs_J01401
        w = kk_pes_J0401
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 101 ) then
        abs_f = kk_abs_J01101
        w = kk_pes_J1101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_J01201
        w = kk_pes_J1201
      elseif ( np_f == 401 ) then
        abs_f = kk_abs_J01401
        w = kk_pes_J1401
      end if
    end if
    deallocate( kk_abs_J01101, kk_pes_J0101, kk_pes_J1101, kk_abs_J01201, kk_pes_J0201 )
    deallocate( kk_pes_J1201, kk_abs_J01401, kk_pes_J0401, kk_pes_J1401 )
    kr = abs_f / dabs(r)
  elseif ( criador == 5 ) then !Para o uso dos filtros de Werthmuller
    allocate( wm_abs_J01201(201), wm_pes_J0201(201), wm_pes_J1201(201))
    allocate( wm_abs_J012001(2001), wm_pes_J02001(2001), wm_pes_J12001(2001) )
    call J0J1Wer( wm_abs_J01201, wm_pes_J0201, wm_pes_J1201, wm_abs_J012001, wm_pes_J02001, wm_pes_J12001 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 201 ) then
        abs_f = wm_abs_J01201
        w = wm_pes_J0201
      elseif ( np_f == 2001 ) then
        abs_f = wm_abs_J012001
        w = wm_pes_J02001
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 201 ) then
        abs_f = wm_abs_J01201
        w = wm_pes_J1201
      elseif ( np_f == 2001 ) then
        abs_f = wm_abs_J012001
        w = wm_pes_J12001
      end if
    end if
    deallocate( wm_abs_J01201, wm_pes_J0201, wm_pes_J1201, wm_abs_J012001, wm_pes_J02001, wm_pes_J12001 )
    kr = abs_f / dabs(r)

  end if
elseif ( tipofiltro == 1 ) then !Para uso de filtros seno ou cosseno
  if ( criador == 1 ) then !Para o uso dos filtros seno e cosseno do Frayzer
    allocate( fl_abs_s19(19), fl_pes_sen19(19), fl_abs_s30(30), fl_pes_sen30(30), fl_abs_s40(40), fl_pes_sen40(40) )
    allocate( fl_abs_c19(19), fl_pes_cos19(19), fl_abs_c30(30), fl_pes_cos30(30), fl_abs_c40(40), fl_pes_cos40(40) )
    call SenCosFrayzer( fl_abs_s19, fl_pes_sen19, fl_abs_s30, fl_pes_sen30, fl_abs_s40, fl_pes_sen40, &
                        fl_abs_c19, fl_pes_cos19, fl_abs_c30, fl_pes_cos30, fl_abs_c40, fl_pes_cos40 )
    if ( idtfcd_f == 2 ) then
      if ( np_f == 19 ) then
        abs_f = fl_abs_s19
        w = fl_pes_sen19
      elseif ( np_f == 30 ) then
        abs_f = fl_abs_s30
        w = fl_pes_sen30
      elseif ( np_f == 40 ) then
        abs_f = fl_abs_s40
        w = fl_pes_sen40
      end if
    elseif ( idtfcd_f == 3 ) then
      if ( np_f == 19 ) then
        abs_f = fl_abs_c19
        w = fl_pes_cos19
      elseif ( np_f == 30 ) then
        abs_f = fl_abs_c30
        w = fl_pes_cos30
      elseif ( np_f == 40 ) then
        abs_f = fl_abs_c40
        w = fl_pes_cos40
      end if
    end if
    deallocate( fl_abs_s19, fl_pes_sen19, fl_abs_s30, fl_pes_sen30, fl_abs_s40, fl_pes_sen40 )
    deallocate( fl_abs_c19, fl_pes_cos19, fl_abs_c30, fl_pes_cos30, fl_abs_c40, fl_pes_cos40 )
    kr = dexp(abs_f) / dabs(r)
  elseif ( criador == 4 ) then !Para o uso dos filtros seno e cosseno do Kerry Key
    allocate( kk_abs_81(81), kk_pes_sen81(81), kk_pes_cos81(81))
    allocate( kk_abs_241(241), kk_pes_sen241(241), kk_pes_cos241(241))
    allocate( kk_abs_601(601), kk_pes_sen601(601), kk_pes_cos601(601))
    allocate( kk_abs_101(101), kk_pes_sen101(101), kk_pes_cos101(101))
    allocate( kk_abs_201(201), kk_pes_sen201(201), kk_pes_cos201(201))
    call SenCosKey( kk_abs_81, kk_pes_sen81, kk_pes_cos81, kk_abs_241, kk_pes_sen241, kk_pes_cos241, &
                    kk_abs_601, kk_pes_sen601, kk_pes_cos601, kk_abs_101, kk_pes_sen101, kk_pes_cos101, &
                    kk_abs_201, kk_pes_sen201, kk_pes_cos201 )
    if ( idtfcd_f == 2 ) then
      if ( np_f == 81 ) then
        abs_f = kk_abs_81
        w = kk_pes_sen81
      elseif ( np_f == 101 ) then
        abs_f = kk_abs_101
        w = kk_pes_sen101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_201
        w = kk_pes_sen201
      elseif ( np_f == 241 ) then
        abs_f = kk_abs_241
        w = kk_pes_sen241
      elseif ( np_f == 601 ) then
        abs_f = kk_abs_601
        w = kk_pes_sen601
      end if
    elseif ( idtfcd_f == 3 ) then
      if ( np_f == 81 ) then
        abs_f = kk_abs_81
        w = kk_pes_cos81
      elseif ( np_f == 101 ) then
        abs_f = kk_abs_101
        w = kk_pes_cos101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_201
        w = kk_pes_cos201
      elseif ( np_f == 241 ) then
        abs_f = kk_abs_241
        w = kk_pes_cos241
      elseif ( np_f == 601 ) then
        abs_f = kk_abs_601
        w = kk_pes_cos601
      end if
    end if
    deallocate( kk_abs_81, kk_pes_sen81, kk_pes_cos81, kk_abs_241, kk_pes_sen241, kk_pes_cos241 )
    deallocate( kk_abs_601, kk_pes_sen601, kk_pes_cos601, kk_abs_101, kk_pes_sen101, kk_pes_cos101 )
    deallocate( kk_abs_201, kk_pes_sen201, kk_pes_cos201 )
    kr = abs_f / dabs(r)
  elseif ( criador == 5 ) then !Para o uso dos filtros seno e cosseno do Werthmuller
    allocate( wm_abs_201(201), wm_pes_sen201(201), wm_pes_cos201(201))
    call SenCosWer( wm_abs_201, wm_pes_sen201, wm_pes_cos201 )
    if ( idtfcd_f == 2 ) then
      abs_f = wm_abs_201
      w = wm_pes_sen201
    elseif ( idtfcd_f == 3 ) then
      abs_f = wm_abs_201
      w = wm_pes_cos201
    end if
    deallocate( wm_abs_201, wm_pes_sen201, wm_pes_cos201 )
    kr = abs_f / dabs(r)
  end if
else
  print*,'especificação incorreta do tipo de filtro'
  stop
end if

  end subroutine
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
subroutine constfiltro2( uso, tipofiltro, criador, idtfcd_f, np_f, r, kr, w )
  use filtros_J0_J1_J2_sen_cos

  implicit none
  real(dp), intent(in) :: r
  integer, intent(in) :: uso, tipofiltro, criador, idtfcd_f, np_f
  real(dp), intent(out) :: Kr( np_f ), w( np_f )
  real(dp) :: abs_f( np_f )

  real(dp), dimension(:), allocatable :: lr_abs_J061, lr_pes_J061, lr_abs_J147, lr_pes_J147
  real(dp), dimension(:), allocatable :: fl_abs_J037, fl_pes_J037, fl_abs_J127, fl_pes_J127, fl_abs_J219, fl_pes_J219
  real(dp), dimension(:), allocatable :: gs_abs_J061, gs_pes_J061, gs_abs_J0120, gs_pes_J0120
  real(dp), dimension(:), allocatable :: gs_abs_J147, gs_pes_J147, gs_abs_J1140, gs_pes_J1140
  real(dp), dimension(:), allocatable :: k_abs_J0161, k_pes_J061, k_pes_J161, k_abs_J01241, k_pes_J0241, k_pes_J1241
  real(dp), dimension(:), allocatable :: kk_abs_J01101, kk_pes_J0101, kk_pes_J1101, kk_abs_J01201, kk_pes_J0201, kk_pes_J1201
  real(dp), dimension(:), allocatable :: kk_abs_J01401, kk_pes_J0401, kk_pes_J1401
  real(dp), dimension(:), allocatable :: wm_abs_J01201, wm_pes_J0201, wm_pes_J1201, wm_abs_J012001, wm_pes_J02001, wm_pes_J12001

  real(dp), dimension(:), allocatable :: fl_abs_s19, fl_pes_sen19, fl_abs_s30, fl_pes_sen30, fl_abs_s40, fl_pes_sen40
  real(dp), dimension(:), allocatable :: fl_abs_c19, fl_pes_cos19, fl_abs_c30, fl_pes_cos30, fl_abs_c40, fl_pes_cos40
  real(dp), dimension(:), allocatable :: kk_abs_81, kk_pes_sen81, kk_pes_cos81, kk_abs_241, kk_pes_sen241, kk_pes_cos241
  real(dp), dimension(:), allocatable :: kk_abs_601, kk_pes_sen601, kk_pes_cos601, kk_abs_101, kk_pes_sen101, kk_pes_cos101
  real(dp), dimension(:), allocatable :: kk_abs_201, kk_pes_sen201, kk_pes_cos201
  real(dp), dimension(:), allocatable :: wm_abs_201, wm_pes_sen201, wm_pes_cos201

  real(dp), dimension(:), allocatable :: fl_abs_J070, fl_pes_J070, fl_abs_J170, fl_pes_J170

! tipofiltro é uma variável que identifica se os filtros a serem usados serão para transformadas de Hankel (0) ou de Fourier (1)
if ( tipofiltro == 0 ) then !Para uso de filtros J0, J1 ou J2
  if ( criador == 0 ) then !Para o uso dos filtros de Rijo
    allocate( lr_abs_J061(61), lr_pes_J061(61), lr_abs_J147(47), lr_pes_J147(47) )
    call J0J1Rijo( lr_abs_J061, lr_pes_J061, lr_abs_J147, lr_pes_J147 )
    if ( idtfcd_f == 0 ) then
      abs_f = lr_abs_J061
      w = lr_pes_J061
    elseif ( idtfcd_f == 1 ) then
      abs_f = lr_abs_J147
      w = lr_pes_J147
    end if
    deallocate( lr_abs_J061, lr_pes_J061, lr_abs_J147, lr_pes_J147 )
    kr = dexp(abs_f) / dabs(r)
  elseif ( criador == 1 ) then !Para o uso dos filtros de Frayzer
    allocate( fl_abs_J037(37), fl_pes_J037(37), fl_abs_J127(27), fl_pes_J127(27), fl_abs_J219(19), fl_pes_J219(19) )
    allocate( fl_abs_J070(70), fl_pes_J070(70), fl_abs_J170(70), fl_pes_J170(70) )
    call J0J1J2Frayzer( fl_abs_J037, fl_pes_J037, fl_abs_J127, fl_pes_J127, fl_abs_J219, fl_pes_J219, &
                        fl_abs_J070, fl_pes_J070, fl_abs_J170, fl_pes_J170 )
    if ( idtfcd_f == 0 .and. uso == 1 ) then
      abs_f = fl_abs_J037
      w = fl_pes_J037
    elseif ( idtfcd_f == 1 .and. uso == 1 ) then
      abs_f = fl_abs_J127
      w = fl_pes_J127
    elseif ( idtfcd_f == 0 .and. uso == 0 ) then
      abs_f = fl_abs_J070
      w = fl_pes_J070
    elseif ( idtfcd_f == 1 .and. uso == 0 ) then
      abs_f = fl_abs_J170
      w = fl_pes_J170
    end if
    deallocate( fl_abs_J037, fl_pes_J037, fl_abs_J127, fl_pes_J127, fl_abs_J219, fl_pes_J219 )
    deallocate( fl_abs_J070, fl_pes_J070, fl_abs_J170, fl_pes_J170 )
    kr = dexp(abs_f) / dabs(r)
  elseif ( criador == 2 ) then !Para o uso dos filtros de Guptasarma
    allocate( gs_abs_J061(61), gs_pes_J061(61), gs_abs_J0120(120), gs_pes_J0120(120), gs_abs_J147(47), gs_pes_J147(47) )
    allocate( gs_abs_J1140(140), gs_pes_J1140(140) )
    call J0J1GS( gs_abs_J061, gs_pes_J061, gs_abs_J0120, gs_pes_J0120, gs_abs_J147, gs_pes_J147, gs_abs_J1140, gs_pes_J1140 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 61 ) then
        abs_f = gs_abs_J061
        w = gs_pes_J061
      elseif ( np_f == 120 ) then
        abs_f = gs_abs_J0120
        w = gs_pes_J0120
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 47 ) then
        abs_f = gs_abs_J147
        w = gs_pes_J147
      elseif ( np_f == 140 ) then
        abs_f = gs_abs_J1140
        w = gs_pes_J1140
      end if
    end if
    deallocate( gs_abs_J061, gs_pes_J061, gs_abs_J0120, gs_pes_J0120, gs_abs_J147, gs_pes_J147, gs_abs_J1140, gs_pes_J1140 )
    kr = ( 10.d0 ** (abs_f) ) / dabs(r)
  elseif ( criador == 3 ) then !Para o uso dos filtros de Kong
    allocate( k_abs_J0161(61), k_pes_J061(61), k_pes_J161(61), k_abs_J01241(241), k_pes_J0241(241), k_pes_J1241(241) )
    call J0J1Kong( k_abs_J0161, k_pes_J061, k_pes_J161, k_abs_J01241, k_pes_J0241, k_pes_J1241 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 61 ) then
        abs_f = k_abs_J0161
        w = k_pes_J061
      elseif ( np_f == 241 ) then
        abs_f = k_abs_J01241
        w = k_pes_J0241
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 61 ) then
        abs_f = k_abs_J0161
        w = k_pes_J161
      elseif ( np_f == 241 ) then
        abs_f = k_abs_J01241
        w = k_pes_J1241
      end if
    end if
    deallocate( k_abs_J0161, k_pes_J061, k_pes_J161, k_abs_J01241, k_pes_J0241, k_pes_J1241 )
    kr = abs_f / dabs(r)
  elseif ( criador == 4 ) then !Para o uso dos filtros de Key
    allocate( kk_abs_J01101(101), kk_pes_J0101(101), kk_pes_J1101(101), kk_abs_J01201(201), kk_pes_J0201(201) )
    allocate( kk_pes_J1201(201), kk_abs_J01401(401), kk_pes_J0401(401), kk_pes_J1401(401) )
    call J0J1Key( kk_abs_J01101, kk_pes_J0101, kk_pes_J1101, kk_abs_J01201, kk_pes_J0201, kk_pes_J1201, &
                  kk_abs_J01401, kk_pes_J0401, kk_pes_J1401 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 101 ) then
        abs_f = kk_abs_J01101
        w = kk_pes_J0101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_J01201
        w = kk_pes_J0201
      elseif ( np_f == 401 ) then
        abs_f = kk_abs_J01401
        w = kk_pes_J0401
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 101 ) then
        abs_f = kk_abs_J01101
        w = kk_pes_J1101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_J01201
        w = kk_pes_J1201
      elseif ( np_f == 401 ) then
        abs_f = kk_abs_J01401
        w = kk_pes_J1401
      end if
    end if
    deallocate( kk_abs_J01101, kk_pes_J0101, kk_pes_J1101, kk_abs_J01201, kk_pes_J0201 )
    deallocate( kk_pes_J1201, kk_abs_J01401, kk_pes_J0401, kk_pes_J1401 )
    kr = abs_f / dabs(r)
  elseif ( criador == 5 ) then !Para o uso dos filtros de Key
    allocate( wm_abs_J01201(201), wm_pes_J0201(201), wm_pes_J1201(201) )
    allocate( wm_abs_J012001(2001), wm_pes_J02001(2001), wm_pes_J12001(2001) )
    call J0J1Wer( wm_abs_J01201, wm_pes_J0201, wm_pes_J1201, wm_abs_J012001, wm_pes_J02001, wm_pes_J12001 )
    if ( idtfcd_f == 0 ) then
      if ( np_f == 201 ) then
        abs_f = wm_abs_J01201
        w = wm_pes_J0201
      elseif ( np_f == 2001 ) then
        abs_f = wm_abs_J012001
        w = wm_pes_J02001
      end if
    elseif ( idtfcd_f == 1 ) then
      if ( np_f == 201 ) then
        abs_f = wm_abs_J01201
        w = wm_pes_J1201
      elseif ( np_f == 2001 ) then
        abs_f = wm_abs_J012001
        w = wm_pes_J12001
      end if
    end if
    deallocate( wm_abs_J01201, wm_pes_J0201, wm_pes_J1201, wm_abs_J012001, wm_pes_J02001, wm_pes_J12001 )
    kr = abs_f / dabs(r)
  end if
elseif ( tipofiltro == 1 ) then !Para uso de filtros seno ou cosseno
  if ( criador == 1 ) then !Para o uso dos filtros seno e cosseno do Frayzer
    allocate( fl_abs_s19(19), fl_pes_sen19(19), fl_abs_s30(30), fl_pes_sen30(30), fl_abs_s40(40), fl_pes_sen40(40) )
    allocate( fl_abs_c19(19), fl_pes_cos19(19), fl_abs_c30(30), fl_pes_cos30(30), fl_abs_c40(40), fl_pes_cos40(40) )
    call SenCosFrayzer( fl_abs_s19, fl_pes_sen19, fl_abs_s30, fl_pes_sen30, fl_abs_s40, fl_pes_sen40, &
                        fl_abs_c19, fl_pes_cos19, fl_abs_c30, fl_pes_cos30, fl_abs_c40, fl_pes_cos40 )
    if ( idtfcd_f == 2 ) then
      if ( np_f == 19 ) then
        abs_f = fl_abs_s19
        w = fl_pes_sen19
      elseif ( np_f == 30 ) then
        abs_f = fl_abs_s30
        w = fl_pes_sen30
      elseif ( np_f == 40 ) then
        abs_f = fl_abs_s40
        w = fl_pes_sen40
      end if
    elseif ( idtfcd_f == 3 ) then
      if ( np_f == 19 ) then
        abs_f = fl_abs_c19
        w = fl_pes_cos19
      elseif ( np_f == 30 ) then
        abs_f = fl_abs_c30
        w = fl_pes_cos30
      elseif ( np_f == 40 ) then
        abs_f = fl_abs_c40
        w = fl_pes_cos40
      end if
    end if
    deallocate( fl_abs_s19, fl_pes_sen19, fl_abs_s30, fl_pes_sen30, fl_abs_s40, fl_pes_sen40 )
    deallocate( fl_abs_c19, fl_pes_cos19, fl_abs_c30, fl_pes_cos30, fl_abs_c40, fl_pes_cos40 )
    kr = dexp(abs_f) / dabs(r)
  elseif ( criador == 4 ) then !Para o uso dos filtros seno e cosseno do Kerry Key
    allocate( kk_abs_81(81), kk_pes_sen81(81), kk_pes_cos81(81))
    allocate( kk_abs_241(241), kk_pes_sen241(241), kk_pes_cos241(241))
    allocate( kk_abs_601(601), kk_pes_sen601(601), kk_pes_cos601(601))
    allocate( kk_abs_101(101), kk_pes_sen101(101), kk_pes_cos101(101))
    allocate( kk_abs_201(201), kk_pes_sen201(201), kk_pes_cos201(201))
    call SenCosKey( kk_abs_81, kk_pes_sen81, kk_pes_cos81, kk_abs_241, kk_pes_sen241, kk_pes_cos241, &
                    kk_abs_601, kk_pes_sen601, kk_pes_cos601, kk_abs_101, kk_pes_sen101, kk_pes_cos101, &
                    kk_abs_201, kk_pes_sen201, kk_pes_cos201 )
    if ( idtfcd_f == 2 ) then
      if ( np_f == 81 ) then
        abs_f = kk_abs_81
        w = kk_pes_sen81
      elseif ( np_f == 101 ) then
        abs_f = kk_abs_101
        w = kk_pes_sen101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_201
        w = kk_pes_sen201
      elseif ( np_f == 241 ) then
        abs_f = kk_abs_241
        w = kk_pes_sen241
      elseif ( np_f == 601 ) then
        abs_f = kk_abs_601
        w = kk_pes_sen601
      end if
    elseif ( idtfcd_f == 3 ) then
      if ( np_f == 81 ) then
        abs_f = kk_abs_81
        w = kk_pes_cos81
      elseif ( np_f == 101 ) then
        abs_f = kk_abs_101
        w = kk_pes_cos101
      elseif ( np_f == 201 ) then
        abs_f = kk_abs_201
        w = kk_pes_cos201
      elseif ( np_f == 241 ) then
        abs_f = kk_abs_241
        w = kk_pes_cos241
      elseif ( np_f == 601 ) then
        abs_f = kk_abs_601
        w = kk_pes_cos601
      end if
    end if
    deallocate( kk_abs_81, kk_pes_sen81, kk_pes_cos81, kk_abs_241, kk_pes_sen241, kk_pes_cos241 )
    deallocate( kk_abs_601, kk_pes_sen601, kk_pes_cos601, kk_abs_101, kk_pes_sen101, kk_pes_cos101 )
    deallocate( kk_abs_201, kk_pes_sen201, kk_pes_cos201 )
    kr = abs_f / dabs(r)
  elseif ( criador == 5 ) then !Para o uso dos filtros seno e cosseno do Werthmuller
    allocate( wm_abs_201(201), wm_pes_sen201(201), wm_pes_cos201(201))
    call SenCosWer( wm_abs_201, wm_pes_sen201, wm_pes_cos201 )
    if ( idtfcd_f == 2 ) then
        abs_f = wm_abs_201
        w = wm_pes_sen201
    elseif ( idtfcd_f == 3 ) then
      abs_f = wm_abs_201
      w = wm_pes_cos201
    end if
    deallocate( wm_abs_201, wm_pes_sen201, wm_pes_cos201 )
    kr = abs_f / dabs(r)
  end if
else
  print*,'especificação incorreta do tipo de filtro'
  stop
end if

  end subroutine constfiltro2

end module escolha_do_filtro
