module test_utils
  use fruit

  implicit none

contains

  subroutine utils_int2str
    use utils, only: int2str
    call set_unit_name('int2str')
    call set_case_name('negative number')
    call assert_equals('-100', int2str(-100))
    call set_case_name('zero')
    call assert_equals('0', int2str(0))
    call set_case_name('positive number')
    call assert_equals('999100', int2str(999100))
  end subroutine utils_int2str

end module test_utils
