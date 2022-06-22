program tatu
   use clifor
   use tatu_io
   use parameters
   use filterscommonbase
   use utils
   use hmdx
   use hmdy
   use vmd
   use hedx
   use hedy
   use ved
   implicit none
   logical :: progress
   integer :: ncam, nT, nR, i, j, k, p, npt, camadT, camadR
   real(sp) :: start, finish
   real(dp) :: pf, Tx, Ty, Tz, f, Rx, Ry, Rz, w, r
   real(dp), dimension(:), allocatable :: sigmas, esp, freq, mydirecT, mydirecR
   real(dp), dimension(:), allocatable :: krJ0J1, wJ0, wJ1, h, prof
   real(dp), dimension(:,:), allocatable :: tmt, rcv, myout
   complex(dp) :: eta0, zeta, Exp, Eyp, Ezp, Hxp, Hyp, Hzp
   complex(dp), dimension(:), allocatable :: eta, AMdw, AMup, FEdw, FEup
   complex(dp), dimension(:), allocatable :: AMdwz, AMupz, FEdwz, FEupz
   complex(dp), dimension(:,:), allocatable :: u, AdmInt, ImpInt, uh, RTEdw, RTMdw, RTEup, RTMup

   character(len=:), allocatable :: input_file, output_file, output_type
   character(len=20), dimension(15) :: labels
   type(json_input) :: in

   character(len=15) :: Ttext
   character(len=12) :: Ftext
   character(len=11) :: Rtext
   character(len=16) :: Ttot, Ftot, Rtot
   integer :: ilenTtot, ilenFtot, ilenRtot
   character(len=:), allocatable :: iterT, iterF, iterR, clenTtot, clenFtot, clenRtot
   character(len=128) :: info
   real(sp) :: rPerc, incPerc
   character(len=3) :: cPerc

   call cpu_time(start)

   Ttext = '%: Transmitter '
   Ftext = '; Frequency '
   Rtext = '; Receptor '

   call clifor_set_program_info( &
      name='tatu', &
      version='1.0.0', &
      pretty_name='Tatu', &
      description='Geophysics Electromagnetic Modeling in 1D Layered Media' &
      )

   call clifor_create_option('v', 'version', 'Show version and exit')
   call clifor_create_option('h', 'help', 'Show this help message')
   call clifor_create_option('i', 'input-file', 'File to read the input data', &
      required=.true., need_value=.true., value_name='FILEPATH')
   call clifor_create_option('o', 'output-file', 'File to write the output data', &
      required=.true., need_value=.true., value_name='FILEPATH')
   call clifor_create_option('t', 'output-type', 'Output file type: json (default), ssv or all', &
      need_value=.true., value_name='FILETYPE')
   call clifor_create_option('p', 'progress', 'Show progress')

   call clifor_read_command_line

   if (clifor_flag_was_provided('help')) call clifor_show_program_help
   if (clifor_flag_was_provided('version')) call clifor_show_program_version

   call clifor_ensure_required_options

   allocate(character(1) :: input_file)
   input_file = clifor_get_value_from_option('input-file')

   allocate(character(1) :: output_file)
   output_file = clifor_get_value_from_option('output-file')

   allocate(character(1) :: output_type)
   output_type = clifor_get_value_from_option('output-type')

   progress = clifor_flag_was_provided('progress')

   call clifor_finalizer

   in = tatu_io_read_input(input_file)

! getting input parameters
   ncam = in%layers%number
   if ( ncam == 1 ) call clifor_write_warning('The thickness was ignored!')
   esp = in%layers%thickness

   allocate(sigmas(ncam), eta(ncam))
   sigmas = 1.d0 / in%layers%resistivity

!constructing array of transmitters:
   select case (in%transmitter%direction)
    case ('x')
      if ( dabs(in%transmitter%step) < eps ) then
         nT = 1
      else
         nT = floor((in%transmitter%final - in%transmitter%initial%x) / in%transmitter%step) + 1
      end if
      allocate(tmt(nT,3), mydirecT(nT))
      tmt(:,1) = in%transmitter%initial%x + (/(i,i=0,nT-1)/) * in%transmitter%step
      tmt(:,2) = in%transmitter%initial%y
      tmt(:,3) = in%transmitter%initial%z
      mydirecT = tmt(:,1)
    case ('y')
      if ( dabs(in%transmitter%step) < eps ) then
         nT = 1
      else
         nT = floor((in%transmitter%final - in%transmitter%initial%y) / in%transmitter%step) + 1
      end if
      allocate(tmt(nT,3), mydirecT(nT))
      tmt(:,1) = in%transmitter%initial%x
      tmt(:,2) = in%transmitter%initial%y + (/(i,i=0,nT-1)/) * in%transmitter%step
      tmt(:,3) = in%transmitter%initial%z
      mydirecT = tmt(:,2)
    case ('z')
      if ( dabs(in%transmitter%step) < eps ) then
         nT = 1
      else
         nT = floor((in%transmitter%final - in%transmitter%initial%z) / in%transmitter%step) + 1
      end if
      allocate(tmt(nT,3), mydirecT(nT))
      tmt(:,1) = in%transmitter%initial%x
      tmt(:,1) = in%transmitter%initial%y
      tmt(:,3) = in%transmitter%initial%z + (/(i,i=0,nT-1)/) * in%transmitter%step
      mydirecT = tmt(:,3)
    case default
      call clifor_write_error('Error associated to transmitters! In direction, step, first or final coordinate!')
      stop 1
   end select
   Ttot = int2str(nT)
   ilenTtot = len(trim(adjustl(Ttot)))
   allocate(character(len=1) :: clenTtot)
   clenTtot = trim(adjustl(int2str(ilenTtot)))
   allocate(character(len=2*ilenTtot+1) :: iterT)
!constructing array of receivers:
   select case (in%receiver%direction)
    case ('x')
      if ( dabs(in%receiver%step) < eps ) then
         nR = 1
      else
         nR = floor((in%receiver%final - in%receiver%initial%x) / in%receiver%step) + 1
      end if
      allocate(rcv(nR,3), mydirecR(nR))
      rcv(:,1) = in%receiver%initial%x + (/(i,i=0,nR-1)/) * in%receiver%step
      rcv(:,2) = in%receiver%initial%y
      rcv(:,3) = in%receiver%initial%z
      mydirecR = rcv(:,1)
    case ('y')
      if ( dabs(in%receiver%step) < eps ) then
         nR = 1
      else
         nR = floor((in%receiver%final - in%receiver%initial%y) / in%receiver%step) + 1
      end if
      allocate(rcv(nR,3), mydirecR(nR))
      rcv(:,1) = in%receiver%initial%x
      rcv(:,2) = in%receiver%initial%y + (/(i,i=0,nR-1)/) * in%receiver%step
      rcv(:,3) = in%receiver%initial%z
      mydirecR = rcv(:,2)
    case ('z')
      if ( dabs(in%receiver%step) < eps ) then
         nR = 1
      else
         nR = floor((in%receiver%final - in%receiver%initial%z) / in%receiver%step) + 1
      end if
      allocate(rcv(nR,3), mydirecR(nR))
      rcv(:,1) = in%receiver%initial%x
      rcv(:,2) = in%receiver%initial%y
      rcv(:,3) = in%receiver%initial%z + (/(i,i=0,nR-1)/) * in%receiver%step
      mydirecR = rcv(:,3)
    case default
      call clifor_write_error('Error associated to receivers! In direction, step, first or final coordinate!')
      stop 1
   end select
   Rtot = int2str(nR)
   ilenRtot = len(trim(adjustl(Rtot)))
   allocate(character(len=1) :: clenRtot)
   clenRtot = trim(adjustl(int2str(ilenRtot)))
   allocate(character(len=2*ilenRtot+1) :: iterR)
!constructing array of frequencies:
   allocate(freq(in%frequency%samples))
   if (in%frequency%samples == 1 ) then
      freq = in%frequency%initial
      call clifor_write_warning('The final frequency was ignored!')
   else
      pf = 10**( dlog10(in%frequency%final / in%frequency%initial) / (in%frequency%samples - 1) )
      freq = in%frequency%initial * pf ** (/(i, i=0, in%frequency%samples - 1)/)
   end if
   Ftot = int2str(in%frequency%samples)
   ilenFtot = len(trim(adjustl(Ftot)))
   allocate(character(len=1) :: clenFtot)
   clenFtot = trim(adjustl(int2str(ilenFtot)))
   allocate(character(len=2*ilenFtot+1) :: iterF)
! allocates dimension of output file
   allocate(myout(nT*in%frequency%samples*nR,15))
!
! information about output file:
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
! Werthmüller's filters:
   npt = 201
   call J0J1Wer(npt, krJ0J1, wJ0, wJ1)
! Alocating common arrays to any dipole:
   allocate(FEdw(npt), FEup(npt), AMdw(npt), AMup(npt))
   allocate(u(npt,0:ncam), uh(npt,0:ncam), AdmInt(npt,0:ncam), ImpInt(npt,0:ncam))
   allocate(RTEdw(npt,0:ncam), RTEup(npt,0:ncam), RTMdw(npt,0:ncam), RTMup(npt,0:ncam))

!selects source and determines eletromagnetic fields in receivers positions
   rPerc = 0
   incPerc = 100.0/nT/in%frequency%samples/nR
   select case (in%transmitter%model)
    case ('hedx')
      p = 1
      do i = 1,nT
         write(iterT, '(a'//clenTtot//',a)') trim(adjustl(int2str(i))), '/'//trim(adjustl(Ttot))
         Tx = tmt(i,1)
         Ty = tmt(i,2)
         Tz = tmt(i,3)
         do j = 1,in%frequency%samples
            write(iterF, '(a'//clenFtot//',a)') trim(adjustl(int2str(j))), '/'//trim(adjustl(Ftot))
            f = freq(j)
            w = 2 * pi * f
            ! eta0 = cmplx(1.d-12,w*epsilon,kind=dp)
            eta0 = cmplx(1.d-12,0.d0,kind=dp) !here, zero imaginary part facilitates convergence at high frequencies
            zeta = cmplx(0,1,kind=dp) * w * mu
            eta = cmplx(sigmas,w*epsilon,kind=dp)
            do k = 1,nR
               write(iterR, '(a'//clenRtot//',a)') trim(adjustl(int2str(k))), '/'//trim(adjustl(Rtot))
               rPerc = rPerc + incPerc
               write(cPerc, '(a3)') trim(adjustl(real2strPerc(rPerc)))
               info = cPerc//Ttext//trim(iterT)//Ftext//trim(iterF)//Rtext//trim(iterR)//'.'
               Rx = rcv(k,1)
               Ry = rcv(k,2)
               Rz = rcv(k,3)
               ! Calculating fields of hedx in most optimized manner:
               call sanitizedata(ncam,Tz,Rz,esp,camadT,camadR,h,prof)
               r = sqrt((Rx - Tx)**2 + (Ry - Ty)**2)
               call commonarraysED(ncam,npt,r,krJ0J1,zeta,eta0,Tz,h,prof,camadT,eta, &
                  u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw,RTMup,FEdw,FEup,AMdw,AMup)
               call hedx_optimized(Tx,Ty,Tz,ncam,camadR,camadT,npt,krJ0J1,wJ0,wJ1,h, &
                  prof,eta,zeta,eta0,Rx,Ry,Rz,u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw, &
                  RTMup,FEdw,FEup,AMdw,AMup,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)

               myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                  real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
               if (progress) call clifor_write_progress(info)
               p = p + 1
            end do
         end do
      end do
    case ('hedy')
      p = 1
      do i = 1,nT
         write(iterT, '(a'//clenTtot//',a)') trim(adjustl(int2str(i))), '/'//trim(adjustl(Ttot))
         Tx = tmt(i,1)
         Ty = tmt(i,2)
         Tz = tmt(i,3)
         do j = 1,in%frequency%samples
            write(iterF, '(a'//clenFtot//',a)') trim(adjustl(int2str(j))), '/'//trim(adjustl(Ftot))
            f = freq(j)
            w = 2 * pi * f
            ! eta0 = cmplx(1.d-12,w*epsilon,kind=dp)
            eta0 = cmplx(1.d-12,0.d0,kind=dp) !here, zero imaginary part facilitates convergence at high frequencies
            zeta = cmplx(0,1,kind=dp) * w * mu
            eta = cmplx(sigmas,w*epsilon,kind=dp)
            do k = 1,nR
               write(iterR, '(a'//clenRtot//',a)') trim(adjustl(int2str(k))), '/'//trim(adjustl(Rtot))
               rPerc = rPerc + incPerc
               write(cPerc, '(a3)') trim(adjustl(real2strPerc(rPerc)))
               info = cPerc//Ttext//trim(iterT)//Ftext//trim(iterF)//Rtext//trim(iterR)//'.'
               Rx = rcv(k,1)
               Ry = rcv(k,2)
               Rz = rcv(k,3)
               ! Calculating fields of hedy in most optimized manner:
               call sanitizedata(ncam,Tz,Rz,esp,camadT,camadR,h,prof)
               r = sqrt((Rx - Tx)**2 + (Ry - Ty)**2)
               call commonarraysED(ncam,npt,r,krJ0J1,zeta,eta0,Tz,h,prof,camadT,eta, &
                  u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw,RTMup,FEdw,FEup,AMdw,AMup)
               call hedy_optimized(Tx,Ty,Tz,ncam,camadR,camadT,npt,krJ0J1,wJ0,wJ1,h, &
                  prof,eta,zeta,eta0,Rx,Ry,Rz,u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw, &
                  RTMup,FEdw,FEup,AMdw,AMup,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)

               myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                  real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
               if (progress) call clifor_write_progress(info)
               p = p + 1
            end do
         end do
      end do
    case ('ved')
      allocate(AMdwz(npt),AMupz(npt))
      p = 1
      do i = 1,nT
         write(iterT, '(a'//clenTtot//',a)') trim(adjustl(int2str(i))), '/'//trim(adjustl(Ttot))
         Tx = tmt(i,1)
         Ty = tmt(i,2)
         Tz = tmt(i,3)
         do j = 1,in%frequency%samples
            write(iterF, '(a'//clenFtot//',a)') trim(adjustl(int2str(j))), '/'//trim(adjustl(Ftot))
            f = freq(j)
            w = 2 * pi * f
            ! eta0 = cmplx(1.d-12,w*epsilon,kind=dp)
            eta0 = cmplx(1.d-12,0.d0,kind=dp) !here, zero imaginary part facilitates convergence at high frequencies
            zeta = cmplx(0,1,kind=dp) * w * mu
            eta = cmplx(sigmas,w*epsilon,kind=dp)
            do k = 1,nR
               write(iterR, '(a'//clenRtot//',a)') trim(adjustl(int2str(k))), '/'//trim(adjustl(Rtot))
               rPerc = rPerc + incPerc
               write(cPerc, '(a3)') trim(adjustl(real2strPerc(rPerc)))
               info = cPerc//Ttext//trim(iterT)//Ftext//trim(iterF)//Rtext//trim(iterR)//'.'
               Rx = rcv(k,1)
               Ry = rcv(k,2)
               Rz = rcv(k,3)
               ! Calculating fields of ved in most optimized manner:
               call sanitizedata(ncam,Tz,Rz,esp,camadT,camadR,h,prof)
               r = sqrt((Rx - Tx)**2 + (Ry - Ty)**2)
               call commonarraysED(ncam,npt,r,krJ0J1,zeta,eta0,Tz,h,prof,camadT,eta, &
                  u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw,RTMup,FEdw,FEup,AMdw,AMup,AMdwz,AMupz)
               call ved_optimized(Tx,Ty,Tz,ncam,camadR,camadT,npt,krJ0J1,wJ0,wJ1,h, &
                  prof,eta,eta0,Rx,Ry,Rz,u,uh,ImpInt,RTMdw,RTMup,AMdwz,AMupz, &
                  Exp,Eyp,Ezp,Hxp,Hyp,Hzp)

               myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                  real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
               if (progress) call clifor_write_progress(info)
               p = p + 1
            end do
         end do
      end do
    case ('hmdx')
      p = 1
      do i = 1,nT
         write(iterT, '(a'//clenTtot//',a)') trim(adjustl(int2str(i))), '/'//trim(adjustl(Ttot))
         Tx = tmt(i,1)
         Ty = tmt(i,2)
         Tz = tmt(i,3)
         do j = 1,in%frequency%samples
            write(iterF, '(a'//clenFtot//',a)') trim(adjustl(int2str(j))), '/'//trim(adjustl(Ftot))
            f = freq(j)
            w = 2 * pi * f
            ! eta0 = cmplx(1.d-12,w*epsilon,kind=dp)
            eta0 = cmplx(1.d-12,0.d0,kind=dp) !here, zero imaginary part facilitates convergence at high frequencies
            zeta = cmplx(0,1,kind=dp) * w * mu
            eta = cmplx(sigmas,w*epsilon,kind=dp)
            do k = 1,nR
               write(iterR, '(a'//clenRtot//',a)') trim(adjustl(int2str(k))), '/'//trim(adjustl(Rtot))
               rPerc = rPerc + incPerc
               write(cPerc, '(a3)') trim(adjustl(real2strPerc(rPerc)))
               info = cPerc//Ttext//trim(iterT)//Ftext//trim(iterF)//Rtext//trim(iterR)//'.'
               Rx = rcv(k,1)
               Ry = rcv(k,2)
               Rz = rcv(k,3)
               ! Calculating fields of hmdx in most optimized manner:
               call sanitizedata(ncam,Tz,Rz,esp,camadT,camadR,h,prof)
               r = sqrt((Rx - Tx)**2 + (Ry - Ty)**2)
               call commonarraysMD(ncam,npt,r,krJ0J1,zeta,eta0,Tz,h,prof,camadT,eta, &
                  u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw,RTMup,FEdw,FEup,AMdw,AMup)
               call hmdx_optimized(Tx,Ty,Tz,ncam,camadR,camadT,npt,krJ0J1,wJ0,wJ1,h, &
                  prof,eta,zeta,eta0,Rx,Ry,Rz,u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw, &
                  RTMup,FEdw,FEup,AMdw,AMup,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)

               myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                  real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
               if (progress) call clifor_write_progress(info)
               p = p + 1
            end do
         end do
      end do
    case ('hmdy')
      p = 1
      do i = 1,nT
         write(iterT, '(a'//clenTtot//',a)') trim(adjustl(int2str(i))), '/'//trim(adjustl(Ttot))
         Tx = tmt(i,1)
         Ty = tmt(i,2)
         Tz = tmt(i,3)
         do j = 1,in%frequency%samples
            write(iterF, '(a'//clenFtot//',a)') trim(adjustl(int2str(j))), '/'//trim(adjustl(Ftot))
            f = freq(j)
            w = 2 * pi * f
            ! eta0 = cmplx(1.d-12,w*epsilon,kind=dp)
            eta0 = cmplx(1.d-12,0.d0,kind=dp) !here, zero imaginary part facilitates convergence at high frequencies
            zeta = cmplx(0,1,kind=dp) * w * mu
            eta = cmplx(sigmas,w*epsilon,kind=dp)
            do k = 1,nR
               write(iterR, '(a'//clenRtot//',a)') trim(adjustl(int2str(k))), '/'//trim(adjustl(Rtot))
               rPerc = rPerc + incPerc
               write(cPerc, '(a3)') trim(adjustl(real2strPerc(rPerc)))
               info = cPerc//Ttext//trim(iterT)//Ftext//trim(iterF)//Rtext//trim(iterR)//'.'
               Rx = rcv(k,1)
               Ry = rcv(k,2)
               Rz = rcv(k,3)
               ! Calculating fields of hmdy in most optimized manner:
               call sanitizedata(ncam,Tz,Rz,esp,camadT,camadR,h,prof)
               r = sqrt((Rx - Tx)**2 + (Ry - Ty)**2)
               call commonarraysMD(ncam,npt,r,krJ0J1,zeta,eta0,Tz,h,prof,camadT,eta, &
                  u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw,RTMup,FEdw,FEup,AMdw,AMup)
               call hmdy_optimized(Tx,Ty,Tz,ncam,camadR,camadT,npt,krJ0J1,wJ0,wJ1,h, &
                  prof,eta,zeta,eta0,Rx,Ry,Rz,u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw, &
                  RTMup,FEdw,FEup,AMdw,AMup,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)

               myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                  real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp),aimag(Hzp)/)
               if (progress) call clifor_write_progress(info)
               p = p + 1
            end do
         end do
      end do
    case ('vmd')
      allocate(FEdwz(npt), FEupz(npt))
      p = 1
      do i = 1,nT
         write(iterT, '(a'//clenTtot//',a)') trim(adjustl(int2str(i))), '/'//trim(adjustl(Ttot))
         Tx = tmt(i,1)
         Ty = tmt(i,2)
         Tz = tmt(i,3)
         do j = 1,in%frequency%samples
            write(iterF, '(a'//clenFtot//',a)') trim(adjustl(int2str(j))), '/'//trim(adjustl(Ftot))
            f = freq(j)
            w = 2 * pi * f
            ! eta0 = cmplx(1.d-12,w*epsilon,kind=dp)
            eta0 = cmplx(1.d-12,0.d0,kind=dp) !here, zero imaginary part facilitates convergence at high frequencies
            zeta = cmplx(0,1,kind=dp) * w * mu
            eta = cmplx(sigmas,w*epsilon,kind=dp)
            do k = 1,nR
               write(iterR, '(a'//clenRtot//',a)') trim(adjustl(int2str(k))), '/'//trim(adjustl(Rtot))
               rPerc = rPerc + incPerc
               write(cPerc, '(a3)') trim(adjustl(real2strPerc(rPerc)))
               info = cPerc//Ttext//trim(iterT)//Ftext//trim(iterF)//Rtext//trim(iterR)//'.'
               Rx = rcv(k,1)
               Ry = rcv(k,2)
               Rz = rcv(k,3)
               ! Calculating fields of vmd in most optimized manner:
               call sanitizedata(ncam,Tz,Rz,esp,camadT,camadR,h,prof)
               r = sqrt((Rx - Tx)**2 + (Ry - Ty)**2)
               call commonarraysMD(ncam,npt,r,krJ0J1,zeta,eta0,Tz,h,prof,camadT,eta, &
                  u,uh,AdmInt,ImpInt,RTEdw,RTEup,RTMdw,RTMup,FEdw,FEup,AMdw,AMup,FEdwz,FEupz)
               call vmd_optimized(Tx,Ty,Tz,ncam,camadR,camadT,npt,krJ0J1,wJ0,wJ1,h, &
                  prof,zeta,Rx,Ry,Rz,u,uh,AdmInt,RTEdw,RTEup,FEdwz,FEupz,Exp,Eyp,Ezp,Hxp,Hyp,Hzp)

               myout(p,:) = (/mydirecT(i),f,mydirecR(k),real(Exp),aimag(Exp),real(Eyp),aimag(Eyp),real(Ezp),aimag(Ezp), &
                  real(Hxp),aimag(Hxp),real(Hyp),aimag(Hyp),real(Hzp), aimag(Hzp)/)
               if (progress) call clifor_write_progress(info)
               p = p + 1
            end do
         end do
      end do
   end select

   if ( output_type /= '' ) then
      select case (output_type)
       case ('json')
         call tatu_io_write_output_json(output_file, in, labels, myout, tmt, freq, rcv)
       case ('ssv')
         call tatu_io_write_output_ssv(output_file, myout)
       case ('all')
         call tatu_io_write_output_json(output_file, in, labels, myout, tmt, freq, rcv)
         call tatu_io_write_output_ssv(output_file, myout)
       case default
         call clifor_write_error('Invalid type. Use "json" or "ssv" only')
         stop
      end select
   else
      call tatu_io_write_output_json(output_file, in, labels, myout, tmt, freq, rcv)
   end if

   call cpu_time(finish)
   write(*,*) achar(13)//achar(11)//achar(0), '[ TIME ]', finish-start

end program tatu
