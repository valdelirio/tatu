module test_types
  use fruit
  implicit none
contains

  subroutine test_types_point
    use types, only: point
    type(point) :: initial = point(1d0, 2d0, 3d0)
    call assert_equals(1d0, initial%x)
    call assert_equals(2d0, initial%y)
    call assert_equals(3d0, initial%z)
  end subroutine test_types_point

  subroutine test_types_transmitter
    use types, only: point, transmitter
    type(point) :: initial
    type(transmitter) :: src
    initial = point(-3d0, -2d0, -1d0)
    src = transmitter('dehx', 'x', initial, 5d-1, 4d0)
    call assert_equals('dehx', src%model)
    call assert_equals('x', src%direction)
    call assert_equals(-3d0, src%initial%x)
    call assert_equals(-2d0, src%initial%y)
    call assert_equals(-1d0, src%initial%z)
    call assert_equals(5d-1, src%step)
    call assert_equals(4d0, src%final)
  end subroutine test_types_transmitter

  subroutine test_types_receiver
    use types, only: point, receiver
    type(point) :: initial
    type(receiver) :: rec
    initial = point(-3d0, -2d0, -1d0)
    rec = receiver('x', initial, 5d-1, 4d0)
    call assert_equals('x', rec%direction)
    call assert_equals(-3d0, rec%initial%x)
    call assert_equals(-2d0, rec%initial%y)
    call assert_equals(-1d0, rec%initial%z)
    call assert_equals(5d-1, rec%step)
    call assert_equals(4d0, rec%final)
  end subroutine test_types_receiver

  ! subroutine test_types_input
  !   use types
  !   type(point) :: initial
  !   type(transmitter) :: src
  !   type(input) :: in
  !   initial = point(8d-1, 6d-1, 4d-1)
  !   src = transmitter('dehx', 'x', initial, 3d-1, 9d0)
  !   in = input(src)
  !   call assert_equals('dehx', in%transmitter%model)
  !   call assert_equals('x', in%transmitter%direction)
  !   call assert_equals(8d-1, in%transmitter%initial%x)
  !   call assert_equals(6d-1, in%transmitter%initial%y)
  !   call assert_equals(4d-1, in%transmitter%initial%z)
  !   call assert_equals(3d-1, in%transmitter%step)
  !   call assert_equals(9d0, in%transmitter%final)
  ! end subroutine test_types_input

end module test_types
