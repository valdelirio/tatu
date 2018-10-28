program InvDipolos
use DEHx
use DEHy
use DEV
use escolha_do_filtro
	implicit none
	real(8), parameter :: pi = 3.1415926535897932384626433832795d0
	real(8), parameter :: mu = 4.d0*pi*1d-7
	real(8), parameter :: eps = 8.85d0*1.d-12
	real(8), parameter :: Iw = 1.d0
	real(8) :: freq, xini, xfin, yini, yfin, zini, zfin, px, py, pz, Tx, Ty, Tz
	real(8) :: omega, t1, t2, dsx, dsy, dsz, ky !,neta0
	real(8), dimension(:), allocatable :: sigmas, h, resist, Rx, Ry, Rz
	complex*16 :: zeta, neta0
	integer(4) :: ncam, nmed, i, j, k, filAnd, autor, filtro, npts, nptc, funs, func
	real(8), dimension(:), allocatable :: kysen, kycos, w_sen, w_cos
	complex*16 :: kerEx, kerEy, kerEz, kerHx, kerHy, kerHz
	complex*16 :: Exky, Eyky, Ezky, Hxky, Hyky, Hzky, Ex, Ey, Ez, Hx, Hy, Hz

	call cpu_time( t1 )
	open(unit = 100, file = 'param.in', status = 'old', action = 'read')

	read(100,*) ncam, freq
	
	allocate( resist(ncam), sigmas(ncam), h(ncam) )

	do i = 1, ncam
	read(100,*) resist(i)
	sigmas(i) = 1.d0 / resist(i)
	end do

	if ( ncam > 1 ) then
		do j = 1, ncam - 1
            read(100,*) h(j)
		end do
	endif

	read(100,*) xini, xfin, px
	read(100,*) yini, yfin, py
	read(100,*) zini, zfin, pz
	read(100,*) Tx, Ty, Tz, dsy	!dsz	!dsx	!
!	Para facilitação nos cálculos de dois semi-espaços
	h(ncam) = 1.d300
	close(100)

	dsx = 1.d0
	dsy = 1.d0
	dsz = 1.d0

	omega = 2.d0 * pi * freq
!	Curiosamente se as respostas do campo elétrico no ar não forem corretas, use neta0 = 1.d-7.
	neta0 = (0.d0,1.d0) * omega * eps	!dcmplx(1.d-7,0.d0)		!
	zeta = (0.d0, 1.d0) * omega * mu
!print*,'i*omega*epsilon =',(0.d0,1.d0)*omega*eps

	if ( xfin > xini .and. (yfin <= yini .or. zfin <= zini) ) then
		nmed = floor( (xfin - xini) / px + 1 )
		allocate ( Rx(nmed), Ry(nmed), Rz(nmed) )
		Rx = xini + (/(i , i = 0, nmed - 1)/) * px
		Ry = yini
		Rz = zini
	else if ( yfin > yini .and. (xfin <= xini .or. zfin <= zini) ) then
		nmed = floor( (yfin - yini) / py + 1 )
		allocate ( Rx(nmed), Ry(nmed), Rz(nmed) )
		Rx = xini
		Ry = yini + (/(j , j = 0, nmed - 1)/) * py
		Rz = zini
	else if ( zfin > zini .and. (xfin <= xini .or. yfin <= yini) ) then
		nmed = floor( (zfin - zini) / pz + 1 )
		allocate ( Rx(nmed), Ry(nmed), Rz(nmed) )
		Rx = xini
		Ry = yini
		Rz = zini + (/(k , k = 0, nmed - 1)/) * pz
    else
        write(*,*)'As posições de medida estão incompatíveis'
        stop
	end if
	10 format(6(1PG24.15E3))
	20 format(3(1PG24.15E3))
	open(unit = 1, file = 'obs.out', status = 'replace', action = 'write')
	open(unit = 2, file = 'Ep.out',  status = 'replace', action = 'write')
	open(unit = 3, file = 'Hp.out',  status = 'replace', action = 'write')

! 	filAnd = 0

do i = 1, nmed
!	call dehx_xyz(filAnd, Tx, Ty, Iw, dsx, Tz, ncam, h, sigmas, neta0, zeta, Rx(i), Ry(i), Rz(i), Ex, Ey, Ez, Hx, Hy, Hz)
	call dehx_xyz_loops(Tx, Ty, Iw, dsx, Tz, ncam, h, sigmas, neta0, zeta, Rx(i), Ry(i), Rz(i), Ex, Ey, Ez, Hx, Hy, Hz)
!	call dehy_xyz(filAnd, Tx, Ty, Iw, dsy, Tz, ncam, h, sigmas, neta0, zeta, Rx(i), Ry(i), Rz(i), Ex, Ey, Ez, Hx, Hy, Hz)
! 	call dehy_xyz_loops(Tx, Ty, Iw, dsy, Tz, ncam, h, sigmas, neta0, zeta, Rx(i), Ry(i), Rz(i), Ex, Ey, Ez, Hx, Hy, Hz)
!	call dev_xyz(filAnd, Tx, Ty, Iw, dsz, Tz, ncam, h, sigmas, neta0, zeta, Rx(i), Ry(i), Rz(i), Ex, Ey, Ez, Hx, Hy, Hz)
!	call dev_xyz_loops(Tx, Ty, Iw, dsz, Tz, ncam, h, sigmas, neta0, zeta, Rx(i), Ry(i), Rz(i), Ex, Ey, Ez, Hx, Hy, Hz)
   	write(1,20) Rx(i), Ry(i), Rz(i)
   	write(2,10) dreal(Ex), dimag(Ex), dreal(Ey), dimag(Ey), dreal(Ez), dimag(Ez)
   	write(3,10) dreal(Hx), dimag(Hx), dreal(Hy), dimag(Hy), dreal(Hz), dimag(Hz)
end do

!	Utilização da transformada inversa dupla de Fourier
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! 	filtro = 1	!Designa o tipo de filtro usado na subrotina de pesos e abscisas de vários filtros.
!!		O algarismo 0 é usado para J0 e J1, enquanto 1 é para seno e cosseno.
! 	autor = 1		!Designa o criador do filtro. No caso de seno ou cosseno os que possuo são os do Frayzer (1) e do Kerry Key (4)
! 	funs = 2		!Designa se é filtro seno (2) ou filtro cosseno (3)
! 	func = 3		!Designa se é filtro seno (2) ou filtro cosseno (3)
! 	npts = 19		!Designa o número de pontos usado no filtro seno.
! 	nptc = 19		!Designa o número de pontos usado no filtro cosseno.!
! 
! 	allocate(kysen( npts ), kycos( nptc ), w_sen( npts ), w_cos( nptc ))
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! do k = 1, nmed
! 	call constfiltro(filtro, autor, funs, npts, Ry(k) - Ty, Kysen, w_sen)
! 	call constfiltro(filtro, autor, func, nptc, Ry(k) - Ty, Kycos, w_cos)
!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!	Dipolo Elétrico Horizontal em x
! 	kerEy = (0.d0, 0.d0)
! 	kerHx = (0.d0, 0.d0)
! 	kerHz = (0.d0, 0.d0)
! 	do i = 1, npts
! 		ky = kysen(i)
!!			call dehx_xkyz(Tx, ky, Iw, dsx, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		call dehx_xkyz_loops(Tx, ky, Iw, dsx, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		kerEy = kerEy + Eyky * w_sen(i)
! 		kerHx = kerHx + Hxky * w_sen(i)
! 		kerHz = kerHz + Hzky * w_sen(i)
! 	end do
! 	Ey = (0.d0,1.d0) * kerEy / (pi * (Ry(k) - Ty))
! 	Hx = (0.d0,1.d0) * kerHx / (pi * (Ry(k) - Ty))
! 	Hz = (0.d0,1.d0) * kerHz / (pi * (Ry(k) - Ty))
! 
! 	kerEx = (0.d0,0.d0)
! 	kerEz = (0.d0,0.d0)
! 	kerHy = (0.d0,0.d0)
! 	do j = 1, nptc
! 		ky = kycos(j)
!!			call dehx_xkyz(Tx, ky, Iw, dsx, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		call dehx_xkyz_loops(Tx, ky, Iw, dsx, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		kerEx = kerEx + Exky * w_cos(j)
! 		kerEz = kerEz + Ezky * w_cos(j)
! 		kerHy = kerHy + Hyky * w_cos(j)
! 	end do
! 	Ex = kerEx / (pi * dabs(Ry(k) - Ty))
! 	Ez = kerEz / (pi * dabs(Ry(k) - Ty))
! 	Hy = kerHy / (pi * dabs(Ry(k) - Ty))
!!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!	Dipolo Elétrico Horizontal em y
! 	kerEx = (0.d0,0.d0)
! 	kerEz = (0.d0,0.d0)
! 	kerHy = (0.d0,0.d0)
! 	do i = 1, npts
! 		ky = kysen(i)
!!			call dehy_xkyz(Tx, ky, Iw, dsy, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		call dehy_xkyz_loops(Tx, ky, Iw, dsy, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		kerEx = kerEx + Exky * w_sen(i)
! 		kerEz = kerEz + Ezky * w_sen(i)
! 		kerHy = kerHy + Hyky * w_sen(i)
! 	end do
! 	Ex = (0.d0,1.d0) * kerEx / (pi * (Ry(k) - Ty))
! 	Ez = (0.d0,1.d0) * kerEz / (pi * (Ry(k) - Ty))
! 	Hy = (0.d0,1.d0) * kerHy / (pi * (Ry(k) - Ty))
! 
! 	kerEy = (0.d0,0.d0)
! 	kerHx = (0.d0,0.d0)
! 	kerHz = (0.d0,0.d0)
! 	do j = 1, nptc
! 		ky = kycos(j)
!!			call dehy_xkyz(Tx, ky, Iw, dsy, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		call dehy_xkyz_loops(Tx, ky, Iw, dsy, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
! 		kerEy = kerEy + Eyky * w_cos(j)
! 		kerHx = kerHx + Hxky * w_cos(j)
! 		kerHz = kerHz + Hzky * w_cos(j)
! 	end do
! 	Ey = kerEy / (pi * dabs(Ry(k) - Ty))
! 	Hx = kerHx / (pi * dabs(Ry(k) - Ty))
! 	Hz = kerHz / (pi * dabs(Ry(k) - Ty))
!!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
!!	Dipolo Elétrico Vertical
!!		kerEy = (0.d0,0.d0)
!!		kerHx = (0.d0,0.d0)
!!		do i = 1, npts
!!			ky = kysen(i)
!!			call dev_xkyz(Tx, ky, Iw, dsz, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
!!			call dev_xkyz_loops(Tx, ky, Iw, dsz, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
!!			kerEy = kerEy + Eyky * w_sen(i)
!!			kerHx = kerHx + Hxky * w_sen(i)
!!		end do
!!		Ey = (0.d0,1.d0) * kerEy / (pi * (Ry(k) - Ty))
!!		Hx = (0.d0,1.d0) * kerHx / (pi * (Ry(k) - Ty))
!!
!!		kerEx = (0.d0,0.d0)
!!		kerEz = (0.d0,0.d0)
!!		kerHy = (0.d0,0.d0)
!!		do j = 1, nptc
!!			ky = kycos(j)
!!			call dev_xkyz(Tx, ky, Iw, dsz, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
!!			call dev_xkyz_loops(Tx, ky, Iw, dsz, Tz, ncam, h, sigmas, neta0, zeta, Rx(k), Rz(k), Exky, Eyky, Ezky, Hxky, Hyky, Hzky)
!!			kerEx = kerEx + Exky * w_cos(j)
!!			kerEz = kerEz + Ezky * w_cos(j)
!!			kerHy = kerHy + Hyky * w_cos(j)
!!		end do
!!		Ex = kerEx / (pi * dabs(Ry(k) - Ty))
!!		Ez = kerEz / (pi * dabs(Ry(k) - Ty))
!!		Hy = kerHy / (pi * dabs(Ry(k) - Ty))
!!-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
! 	write(1,20) Rx(k), Ry(k), Rz(k)
! 	write(2,10) dreal(Ex), dimag(Ex), dreal(Ey), dimag(Ey), dreal(Ez), dimag(Ez)
! 	write(3,10) dreal(Hx), dimag(Hx), dreal(Hy), dimag(Hy), dreal(Hz), dimag(Hz)
! end do

	call cpu_time(t2)

	write(*,*)'tempo de processamento',(t2-t1),'segundos'

end program InvDipolos
