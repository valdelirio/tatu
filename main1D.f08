program main1Dmod
use computational_stuff
use DMHx
use DMHy
use DMV
use DEHx
use DEHy
use DEV
use csv_file
implicit none
character(len=20) :: sourcetype, sourceT, iniposT, finposT, orienT, direcT, stepT
character(len=20) :: iniposR, finposR, orienR, direcR, stepR
character(len=20) :: inif, finf, numbf, numlay, hlay, reslay
logical :: advance
! character(len=10) :: namesrc, namefrq, namercv
! character(len=10) :: namerEx, namerEy, namerEz, namerHx, namerHy, namerHz
! character(len=10) :: nameiEx, nameiEy, nameiEz, nameiHx, nameiHy, nameiHz
character(len=123) :: info
integer :: ncam, nT, nR, nf, i, j, k, numfile, p !, q
real(dp) :: t1, t2
real(dp) :: Tx1, Ty1, Tz1, pT, Tfin, Rx1, Ry1, Rz1, pR, Rfin, fini, ffin
real(dp) :: pf, Tx, Ty, Tz, f, Rx, Ry, Rz, w !, myout(15)
real(dp), dimension(:), allocatable :: resist, sigmas, h, freq, mydirecT, mydirecR
real(dp), dimension(:,:), allocatable :: tmt, rcv, myout
complex(dp) :: eta0, zeta, Exp, Eyp, Ezp, Hxp, Hyp, Hzp

! write(*,*) sp,dp, selected_real_kind(20) !
! read(*,*)

call cpu_time(t1)
!reading input file
open( unit = 100, file = 'param.in', status = 'old', action = 'read' )
read(100,*)sourcetype, sourceT

read(100,*)iniposT, Tx1, Ty1, Tz1
read(100,*)orienT, direcT
read(100,*)stepT, pT
read(100,*)finposT, Tfin

read(100,*)iniposR, Rx1, Ry1, Rz1
read(100,*)orienR, direcR
read(100,*)stepR, pR
read(100,*)finposR, Rfin

read(100,*)inif, fini
read(100,*)numbf, nf
read(100,*)finf, ffin

read(100,*)numlay, ncam
allocate(resist(ncam), sigmas(ncam), h(ncam-1))
read(100,*)hlay, (h(i),i=1,ncam-1)
read(100,*)reslay, (resist(i),i=1,ncam)
close(100)
sigmas = 1.d0 / resist
!constructing array of transmitters:
select case (direcT)
case ('x')
  nT = floor((Tfin-Tx1)/pT) + 1
  allocate(tmt(nT,3), mydirecT(nT))
  tmt(:,1) = Tx1 + (/(i,i=0,nT-1)/) * pT
  tmt(:,2) = Ty1
  tmt(:,3) = Tz1
  mydirecT = tmt(:,1)
case ('y')
  nT = floor((Tfin-Ty1)/pT) + 1
  allocate(tmt(nT,3), mydirecT(nT))
  tmt(:,1) = Tx1
  tmt(:,2) = Ty1 + (/(i,i=0,nT-1)/) * pT
  tmt(:,3) = Tz1
  mydirecT = tmt(:,2)
case ('z')
  nT = floor((Tfin-Tz1)/pT) + 1
allocate(tmt(nT,3), mydirecT(nT))
  tmt(:,1) = Tx1
  tmt(:,1) = Ty1
  tmt(:,3) = Tz1 + (/(i,i=0,nT-1)/) * pT
  mydirecT = tmt(:,3)
case default
  stop'Error associated to transmitters! In direction, step, first or final coordinate!'
end select

!constructing array of receivers:
select case (direcR)
case ('x')
  nR = floor((Rfin-Rx1)/pR) + 1
  allocate(rcv(nR,3), mydirecR(nR))
  rcv(:,1) = Rx1 + (/(i,i=0,nR-1)/) * pR
  rcv(:,2) = Ry1
  rcv(:,3) = Rz1
  mydirecR = rcv(:,1)
case ('y')
  nR = floor((Rfin-Ry1)/pR) + 1
allocate(rcv(nR,3), mydirecR(nR))
  rcv(:,1) = Rx1
  rcv(:,2) = Ry1 + (/(i,i=0,nR-1)/) * pR
  rcv(:,3) = Rz1
  mydirecR = rcv(:,2)
case ('z')
  nR = floor((Rfin-Rz1)/pR) + 1
  allocate(rcv(nR,3), mydirecR(nR))
  rcv(:,1) = Rx1
  rcv(:,2) = Ry1
  rcv(:,3) = Rz1 + (/(i,i=0,nR-1)/) * pR
  mydirecR = rcv(:,3)
case default
  stop'Error associated to receivers! In direction, step, first or final coordinate!'
end select

!constructing array of frequencies:
allocate(freq(nf))
pf = 10**( dlog10(ffin/fini) / (nf-1) )
freq = fini * pf ** (/(i,i=0,nf-1)/)
! allocates dimension of output file
allocate(myout(nT*nf*nR,15))
!
! information about output file:
info = 'source, frequency, receiver, ExReal, ExImag, EyReal, EyImag, EzReal, EzImag, HxReal, HxImag, HyReal, HyImag, HzReal, HzImag'
!advance = .true. faz com que o último argumento escrito no arquivo não tenha virgula depois dele.
!Se advance = .false. então depois do argumento aparece uma virgula
advance = .true.
! numfile = 1000
open(newunit=numfile,file='output.csv',status='replace',action='write')
! written in first line of output file:
! call csv_write_char( numfile, info, advance )
write(numfile,*)info
! the others lines are will insert in each iteration

!selects source and determines eletromagnetic fields in receivers positions
select case (sourceT)
  case ('dehx')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,nf
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dehx_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
                        real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          p = p + 1
          ! myout = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
          !           real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          ! call csv_write_dble_1d(numfile, myout, advance)
        end do
      end do
    end do
  case ('dehy')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,nf
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dehy_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
                        real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          p = p + 1
        end do
      end do
    end do
  case ('dev')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,nf
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dev_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
                        real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          p = p + 1
        end do
      end do
    end do
  case ('dmhx')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,nf
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dmhx_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
                        real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          p = p + 1
        end do
      end do
    end do
  case ('dmhy')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,nf
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dmhy_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
                        real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          p = p + 1
        end do
      end do
    end do
  case ('dmv')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,nf
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dmv_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),imag(Exp),real(Eyp),imag(Eyp),real(Ezp),imag(Ezp), &
                        real(Hxp),imag(Hxp),real(Hyp),imag(Hyp),real(Hzp),imag(Hzp)/)
          p = p + 1
        end do
      end do
    end do
  case default
    stop'Source was entered incorretly'
  end select

call csv_write_dble_2d(numfile, transpose(myout))

call cpu_time(t2)
write(*,*)'Tempo processamento',t2-t1,'segundos'

end program main1Dmod
