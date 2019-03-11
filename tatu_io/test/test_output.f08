module test_output
  use fruit

contains

  subroutine test_output_write
    use output, only: output_write
    ! character(len=:), allocatable :: argument
    ! argument = cli_get_argument(0)
    ! call assert_equals('./driver.bin', argument)
  end subroutine test_output_write

end module test_output
