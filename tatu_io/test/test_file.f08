module test_file
  use fruit
  implicit none
contains

  subroutine test_file_already_exist
    use file, only: file_already_exist
    logical :: already_exist
    already_exist = file_already_exist('./data/input.json')
    call assert_equals(.true., already_exist)
    already_exist = file_already_exist('inexistent file.not')
    call assert_equals(.false., already_exist)
  end subroutine test_file_already_exist

  subroutine test_file_already_open
    use file, only: file_already_open
    logical :: already_open
    integer :: unit
    character(len=:), allocatable :: filepath
    filepath = './data/input.json'
    open(newunit=unit, file=filepath)
    already_open = file_already_open(filepath)
    call assert_equals(.true., already_open)
    close(unit)
    already_open = file_already_open(filepath)
    call assert_equals(.false., already_open)
  end subroutine test_file_already_open

  subroutine test_file_is_readable
    use file, only: file_is_readable
    logical :: is_readable
    is_readable = file_is_readable('./data/input.json')
    call assert_equals(.true., is_readable)
    is_readable = file_is_readable('inexistent file.not')
    call assert_equals(.false., is_readable)
  end subroutine test_file_is_readable

  subroutine test_file_is_writable
    use file, only: file_is_writable
    logical :: is_writable
    is_writable = file_is_writable('./data/input.json')
    call assert_equals(.true., is_writable)
    is_writable = file_is_writable('inexistent file.not')
    call assert_equals(.false., is_writable)
  end subroutine test_file_is_writable

end module test_file
