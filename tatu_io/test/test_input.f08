module test_input
  use fruit
  implicit none

contains

  subroutine test_input_transmitter_model
    use input, only: tatu_io_input_check_transmitter_model
    integer :: i
    logical :: isvalid
    character(len=4), dimension(6) :: models

    models = (/'hedx', 'hedy', ' ved', 'hmdx', 'hmdy', ' vmd'/)
    do i=1,6
      isvalid = tatu_io_input_check_transmitter_model(adjustl(models(i)))
      call assert_equals(.true., isvalid)
    end do

    isvalid = tatu_io_input_check_transmitter_model('nonexistent')
    call assert_equals(.false., isvalid)
  end subroutine test_input_transmitter_model


  subroutine test_input_direction
    use input, only: tatu_io_input_check_direction
    integer :: i
    logical :: isvalid
    character(len=1), dimension(3) :: directions

    directions = (/'x', 'y', 'z'/)
    do i=1,3
      isvalid = tatu_io_input_check_direction(directions(i))
      call assert_equals(.true., isvalid)
    end do

    isvalid = tatu_io_input_check_direction('nonexistent')
    call assert_equals(.false., isvalid)
  end subroutine test_input_direction


  subroutine test_input_step
    use input, only: tatu_io_input_check_step
    use types, only: real_dp
    logical :: isvalid

    isvalid = tatu_io_input_check_step(1.d0, 2.d0, 5.d0)
    call assert_equals(.true., isvalid)

    isvalid = tatu_io_input_check_step(1.d0, -2.d0, 5.d0)
    call assert_equals(.false., isvalid)

    isvalid = tatu_io_input_check_step(5.d0, 2.d0, 1.d0)
    call assert_equals(.false., isvalid)

    isvalid = tatu_io_input_check_step(5.d0, -2.d0, 1.d0)
    call assert_equals(.true., isvalid)

    isvalid = tatu_io_input_check_step(4.d0, 2.d0, 4.d0)
    call assert_equals(.true., isvalid)

    isvalid = tatu_io_input_check_step(4.d0, -2.d0, 4.d0)
    call assert_equals(.true., isvalid)

    isvalid = tatu_io_input_check_step(8.d0, 0.d0, 8.d0)
    call assert_equals(.true., isvalid)
  end subroutine test_input_step


  subroutine test_input_get_type
    use input, only: tatu_io_input_get_type
    character(len=:), allocatable :: type

    type = tatu_io_input_get_type('Égua')
    call assert_equals('string', type)

    type = tatu_io_input_get_type(3)
    call assert_equals('integer', type)

    type = tatu_io_input_get_type(3.d0)
    call assert_equals('realdp', type)

  end subroutine test_input_get_type
end module
