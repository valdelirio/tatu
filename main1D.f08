program main1Dmod
use cli
use json_io
use computational_stuff
use DMHx
use DMHy
use DMV
use DEHx
use DEHy
use DEV
use csv_file
implicit none
integer :: ncam, nT, nR, i, j, k, numfile, p
real(dp) :: t1, t2, pf, Tx, Ty, Tz, f, Rx, Ry, Rz, w
real(dp), dimension(:), allocatable :: sigmas, h, freq, mydirecT, mydirecR
real(dp), dimension(:,:), allocatable :: tmt, rcv, myout
complex(dp) :: eta0, zeta, Exp, Eyp, Ezp, Hxp, Hyp, Hzp
logical :: advance

character(len=:), allocatable :: input_file, output_file, info
character(len=20), dimension(15) :: labels
type(json_input) :: in

call cpu_time(t1)

input_file = cli_get_option_value('-i')
in = json_io_read_input(input_file)

! getting input parameters
ncam = in%layers%number
h = in%layers%thickness

sigmas = 1.d0 / in%layers%resistivity

!constructing array of transmitters:
select case (in%transmitter%direction)
  case ('x')
    nT = floor((in%transmitter%final - in%transmitter%initial%x) / in%transmitter%step) + 1
    allocate(tmt(nT,3), mydirecT(nT))
    tmt(:,1) = in%transmitter%initial%x + (/(i,i=0,nT-1)/) * in%transmitter%step
    tmt(:,2) = in%transmitter%initial%y
    tmt(:,3) = in%transmitter%initial%z
    mydirecT = tmt(:,1)
  case ('y')
    nT = floor((in%transmitter%final - in%transmitter%initial%y) / in%transmitter%step) + 1
    allocate(tmt(nT,3), mydirecT(nT))
    tmt(:,1) = in%transmitter%initial%x
    tmt(:,2) = in%transmitter%initial%y + (/(i,i=0,nT-1)/) * in%transmitter%step
    tmt(:,3) = in%transmitter%initial%z
    mydirecT = tmt(:,2)
  case ('z')
    nT = floor((in%transmitter%final - in%transmitter%initial%z) / in%transmitter%step) + 1
  allocate(tmt(nT,3), mydirecT(nT))
    tmt(:,1) = in%transmitter%initial%x
    tmt(:,1) = in%transmitter%initial%y
    tmt(:,3) = in%transmitter%initial%z + (/(i,i=0,nT-1)/) * in%transmitter%step
    mydirecT = tmt(:,3)
  case default
    stop 'Error associated to transmitters! In direction, step, first or final coordinate!'
end select

!constructing array of receivers:
select case (in%receiver%direction)
  case ('x')
    nR = floor((in%receiver%final - in%receiver%initial%x) / in%receiver%step) + 1
    allocate(rcv(nR,3), mydirecR(nR))
    rcv(:,1) = in%receiver%initial%x + (/(i,i=0,nR-1)/) * in%receiver%step
    rcv(:,2) = in%receiver%initial%y
    rcv(:,3) = in%receiver%initial%z
    mydirecR = rcv(:,1)
  case ('y')
    nR = floor((in%receiver%final - in%receiver%initial%y) / in%receiver%step) + 1
  allocate(rcv(nR,3), mydirecR(nR))
    rcv(:,1) = in%receiver%initial%x
    rcv(:,2) = in%receiver%initial%y + (/(i,i=0,nR-1)/) * in%receiver%step
    rcv(:,3) = in%receiver%initial%z
    mydirecR = rcv(:,2)
  case ('z')
    nR = floor((in%receiver%final - in%receiver%initial%z) / in%receiver%step) + 1
    allocate(rcv(nR,3), mydirecR(nR))
    rcv(:,1) = in%receiver%initial%x
    rcv(:,2) = in%receiver%initial%y
    rcv(:,3) = in%receiver%initial%z + (/(i,i=0,nR-1)/) * in%receiver%step
    mydirecR = rcv(:,3)
  case default
    stop 'Error associated to receivers! In direction, step, first or final coordinate!'
end select

!constructing array of frequencies:
allocate(freq(in%frequency%samples))
pf = 10**( dlog10(in%frequency%final / in%frequency%initial) / (in%frequency%samples - 1) )
freq = in%frequency%initial * pf ** (/(i, i=0, in%frequency%samples - 1)/)
! allocates dimension of output file
allocate(myout(nT*in%frequency%samples*nR,15))
!
! information about output file:
info = 'transmitter, frequency, receiver, '&
     //'ExReal, ExImag, EyReal, EyImag, EzReal, EzImag, '&
     //'HxReal, HxImag, HyReal, HyImag, HzReal, HzImag'
labels(1) = 'transmitter'
labels(2) = 'frequency'
labels(3) = 'receiver'
labels(4) = 'ExReal'
labels(5) = 'ExImag'
labels(6) = 'EyReal'
labels(7) = 'EyImag'
labels(8) = 'EzReal'
labels(9) = 'EzImag'
labels(10) = 'HxReal'
labels(11) = 'HxImag'
labels(12) = 'HyReal'
labels(13) = 'HyImag'
labels(14) = 'HzReal'
labels(15) = 'HzImag'

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
select case (in%transmitter%model)
  case ('dehx')
    p = 1
    do i = 1,nT
      Tx = tmt(i,1)
      Ty = tmt(i,2)
      Tz = tmt(i,3)
      do j = 1,in%frequency%samples
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dehx_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                        real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
          p = p + 1
          ! myout = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
          !           real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
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
      do j = 1,in%frequency%samples
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dehy_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                        real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
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
      do j = 1,in%frequency%samples
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dev_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                        real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
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
      do j = 1,in%frequency%samples
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dmhx_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                        real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
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
      do j = 1,in%frequency%samples
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dmhy_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                        real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
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
      do j = 1,in%frequency%samples
        f = freq(j)
        w = 2 * pi * f
        eta0 = cmplx(0,1,kind=dp) * w * epsilon  !cmplx(1.d-12,0.d0,kind=dp) !cmplx(1.d-7,0.d0,kind=dp) !
        zeta = cmplx(0,1,kind=dp) * w * mu
        do k = 1,nR
          Rx = rcv(k,1)
          Ry = rcv(k,2)
          Rz = rcv(k,3)
          call dmv_xyz_loops(Tx,Ty,Tz,ncam,h,sigmas,eta0,zeta,Rx,Ry,Rz,Exp,Eyp,Hxp,Hyp,Hzp)
          myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                        real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
          p = p + 1
        end do
      end do
    end do
  case default
    stop 'Source was entered incorretly'
  end select
! the file output is stored in csv format
! The first line gives information about the columns.
! Each line registers the field components in each receiver at a certain frequency,
! for a given transmitter position
call csv_write_dble_2d(numfile, transpose(myout))

output_file = cli_get_option_value('-o')
call json_io_write_output(output_file, in, labels, myout, tmt, freq, rcv)

call cpu_time(t2)
write(*,*) t2-t1

end program main1Dmod
