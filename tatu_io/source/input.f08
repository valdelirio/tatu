module input
  use, intrinsic :: iso_fortran_env, only: error_unit
  use json_module
  use types

  interface tatu_io_get
    module procedure tatu_io_get_string
    module procedure tatu_io_get_integer
    module procedure tatu_io_get_real
    module procedure tatu_io_get_real_array
    module procedure tatu_io_get_real_point
  end interface tatu_io_get

  interface tatu_io_input_get_type
    module procedure tatu_io_input_get_type_string
    module procedure tatu_io_input_get_type_integer
    module procedure tatu_io_input_get_type_realdp
  end interface tatu_io_input_get_type

contains

  ! IDEA Maybe create a generic interface for error
  subroutine error(message, filename)
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: filename
    if (present(filename)) then
      write(error_unit, '(a)') '[ ERROR ] '//message//': '//filename
    else
      write(error_unit, '(a)') '[ ERROR ] '//message
    end if
    stop 1
  end subroutine

  subroutine warning(message)
    character(len=*), intent(in) :: message
    write(error_unit, '(a)') '[ WARNING ] '//message
  end subroutine


  subroutine path_not_found(path)
    character(len=*), intent(in) :: path
    call error('"'//path//'": expected a numeric value')
  end subroutine path_not_found


  ! TODO Split and separeted procedures to each input file session
  subroutine tatu_io_get_input(input_file, in)
    character(len=*), intent(in) :: input_file
    type(json_input), intent(out) :: in
    type(json_file) :: json
    real(real_dp) :: initial_value_at_direction
    real(real_dp), parameter :: eps = 1.d-12
    logical :: valid

    call json%initialize()
    call json%load_file(filename = input_file)
    if (json%failed()) call error('no such file', input_file)

    call tatu_io_get(json, 'transmitter.model', in%transmitter%model)
    valid = tatu_io_input_check_transmitter_model(in%transmitter%model)
    if ( .not. valid ) then
      call error('Invalid transmitter model: '//trim(in%transmitter%model)// &
                 '. Available models: hedx, hedy, ved, hmdx, hmdy, vmd')
    end if

    call tatu_io_get(json, 'transmitter.direction', in%transmitter%direction)
    valid = tatu_io_input_check_direction(in%transmitter%direction)
    if ( .not. valid ) then
      call error('Invalid transmitter direction: '//in%transmitter%direction// &
                 '. Available directions: x, y, z')
    end if

    call tatu_io_get(json, 'transmitter.initial', in%transmitter%initial)
    call tatu_io_get(json, 'transmitter.step', in%transmitter%step)
    call tatu_io_get(json, 'transmitter.final', in%transmitter%final)
    select case (in%transmitter%direction)
      case ('x')
        initial_value_at_direction = in%transmitter%initial%x
      case ('y')
        initial_value_at_direction = in%transmitter%initial%y
      case ('z')
        initial_value_at_direction = in%transmitter%initial%z
    end select
    valid = tatu_io_input_check_step(initial_value_at_direction, in%transmitter%step, in%transmitter%final)
    if ( .not. valid ) then
      if ( in%transmitter%step < 0) then
        call error('You must provide transmitter step > 0 when initial < final')
      elseif ( in%transmitter%step > 0) then
        call error('You must provide transmitter step < 0 when initial > final')
      end if
    elseif (abs(in%transmitter%final - initial_value_at_direction) < eps) then
      call warning('Transmitter step will be ignored, because initial = final.')
    end if

    call tatu_io_get(json, 'receiver.direction', in%receiver%direction)
    valid = tatu_io_input_check_direction(in%receiver%direction)
    if ( .not. valid ) then
      call error('Invalid receiver direction: '//in%receiver%direction// &
                 '. Available directions: x, y, z')
    end if

    call tatu_io_get(json, 'receiver.initial', in%receiver%initial)
    call tatu_io_get(json, 'receiver.step', in%receiver%step)
    call tatu_io_get(json, 'receiver.final', in%receiver%final)
    select case (in%receiver%direction)
      case ('x')
        initial_value_at_direction = in%receiver%initial%x
      case ('y')
        initial_value_at_direction = in%receiver%initial%y
      case ('z')
        initial_value_at_direction = in%receiver%initial%z
    end select
    valid = tatu_io_input_check_step(initial_value_at_direction, in%receiver%step, in%receiver%final)
    if ( .not. valid ) then
      if ( in%receiver%step < 0) then
        call error('You must provide receiver step > 0 when initial < final')
      elseif ( in%receiver%step > 0) then
        call error('You must provide receiver step < 0 when initial > final')
      end if
    elseif (abs(in%receiver%final - initial_value_at_direction) < eps) then
      call warning('Receiver step will be ignored, because initial = final.')
    end if

    call tatu_io_get(json, 'frequency.initial', in%frequency%initial)
    call tatu_io_get(json, 'frequency.samples', in%frequency%samples)
    call tatu_io_get(json, 'frequency.final', in%frequency%final)

    call tatu_io_get(json, 'layers.number', in%layers%number)
    call tatu_io_get(json, 'layers.resistivity', in%layers%resistivity)
    call tatu_io_get(json, 'layers.thickness', in%layers%thickness)

    if (size(in%layers%resistivity) /= in%layers%number) then
      call error('The resistivity array size must be equal to layers number.')
    end if

    if (int(in%layers%number) > 1 .and. size(in%layers%thickness) /= int(in%layers%number) - 1) then
      call error('The thickness array size must be equal to layers number minus 1.')
    end if

    call json%destroy()
    if (json%failed()) call error('error on destroy')
  end subroutine tatu_io_get_input


  subroutine tatu_io_get_string(json, path, tatu_io_string)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: path
    character(len=*), intent(out) :: tatu_io_string
    character(len=:), allocatable :: temp
    logical :: found
    call json%get(path, temp, found)
    tatu_io_string = temp
    if (.not. found) call path_not_found(path)
  end subroutine tatu_io_get_string


  subroutine tatu_io_get_integer(json, path, tatu_io_integer)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: path
    integer, intent(out) :: tatu_io_integer
    logical :: found
    call json%get(path, tatu_io_integer, found)
    if (.not. found) call path_not_found(path)
  end subroutine tatu_io_get_integer


  subroutine tatu_io_get_real(json, path, tatu_io_real)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: path
    real(real_dp), intent(out) :: tatu_io_real
    logical :: found
    call json%get(path, tatu_io_real, found)
    if (.not. found) call path_not_found(path)
  end subroutine tatu_io_get_real


  subroutine tatu_io_get_real_array(json, path, tatu_io_real_array)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: path
    real(real_dp), dimension(:), allocatable, intent(out) :: tatu_io_real_array
    logical :: found
    call json%get(path, tatu_io_real_array, found)
    if (.not. found) call path_not_found(path)
  end subroutine tatu_io_get_real_array


  subroutine tatu_io_get_real_point(json, path, tatu_io_real_point)
    type(json_file), intent(inout) :: json
    character(len=*), intent(in) :: path
    type(point), intent(out) :: tatu_io_real_point
    real(real_dp), dimension(:), allocatable :: array
    logical :: found
    call json%get(path, array, found)
    if (.not. found) call path_not_found(path)
    tatu_io_real_point = point(array(1), array(2), array(3))
  end subroutine tatu_io_get_real_point


  function tatu_io_input_check_transmitter_model(model) result(isvalid)
    character(len=*), intent(in) :: model
    logical :: isvalid
    select case (model)
    case ('hedx', 'hedy', 'ved', 'hmdx', 'hmdy', 'vmd')
        isvalid = .true.
      case default
        isvalid = .false.
    end select
  end function tatu_io_input_check_transmitter_model


  function tatu_io_input_check_direction(direction) result(isvalid)
    character(len=*), intent(in) :: direction
    logical :: isvalid
    select case (direction)
    case ('x', 'y', 'z')
        isvalid = .true.
      case default
        isvalid = .false.
    end select
  end function tatu_io_input_check_direction


  function tatu_io_input_check_step(initial, step, final) result(isvalid)
    real(real_dp) :: initial, step, final
    logical :: isvalid
    if ( (initial < final .and. step > 0) .or. (initial > final .and. step < 0) ) then
      isvalid = .true.
    elseif ( initial < final .and. step < 0 ) then
      isvalid = .false.
    elseif ( initial > final .and. step > 0 ) then
      isvalid = .false.
    else
      isvalid = .true.
    end if
  end function tatu_io_input_check_step


  function tatu_io_input_get_type_string(value) result(type)
    character(len=*), intent(in) :: value
    character(len=128) :: unused
    character(len=6) :: type
    type = 'string'
    unused = value
  end function tatu_io_input_get_type_string


  function tatu_io_input_get_type_integer(value) result(type)
    integer, intent(in) :: value
    integer :: unused
    character(len=7) :: type
    type = 'integer'
    unused = value
  end function tatu_io_input_get_type_integer


  function tatu_io_input_get_type_realdp(value) result(type)
    real(real_dp), intent(in) :: value
    real(real_dp) :: unused
    character(len=6) :: type
    type = 'realdp'
    unused = value
  end function tatu_io_input_get_type_realdp
end module input
