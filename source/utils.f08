module utils
use parameters
contains
subroutine sanitizedata(n, h0, z, esp, camadT, camad, h, prof)
  implicit none
  integer :: n
  real(dp), intent(in) :: h0, z, esp(:)
  integer, intent(out) :: camadT, camad
  real(dp), dimension(:), allocatable, intent(out) :: h, prof

  integer :: i, j, k

  allocate( h(0:n), prof(-1:n) )
  if (size(esp) == n) then
    h(0)=0.d0
    h(1:n)=esp
  else
    h(0)=0.d0
    h(1:n-1)=esp
    h(n)=1.d300
  end if
  ! create depths array that suits any pathological situation
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

  ! find the layer where the receiver is
  camad = 0
  if (z < 0.d0) then
    camad=0
  else if (z >= prof(n-1)) then
    camad=n
  else
    do i=n-1,1,-1
      if (z >= prof(i-1)) then
        camad=i
        exit
      end if
    end do
  end if

  ! find the layer where the transmitter is
  camadT = 0
  if (h0 < 0.d0) then
    camadT = 0
  else if (h0 >= prof(n-1)) then
    camadT = n
  else
    do j=n-1,1,-1
      if (h0 >= prof(j-1)) then
        camadT = j
        exit
      end if
    end do
  end if
end subroutine sanitizedata

function int2str(num)
  implicit none
  integer, intent(in) :: num
  character(len=16) :: int2str
  write(int2str, '(i16)') num
end function int2str

function real2str(real, format)
  implicit none
  real(sp), intent(in) :: real
  character(len=*), intent(in), optional :: format
  character(len=16) :: real2str
  if (present(format)) then
    write(real2str, format) real
  else
    write(real2str, '(f6.2)') real
  end if
end function real2str

function real2strPerc(real)
  implicit none
  real(sp), intent(in) :: real
  character(len=3) :: real2strPerc
  integer :: percentage
  if (real > 99.1 .and. real < 100.9) then
    percentage = 100
  else
    percentage = int(real)
  end if
  write(real2strPerc, '(i3)') percentage
end function real2strPerc

end module utils
