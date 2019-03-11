module file

contains

  function file_already_exist(filepath) result(already_exist)
    character(len=*), intent(in) :: filepath
    logical :: already_exist
    inquire(file=filepath, exist=already_exist)
  end function file_already_exist

  function file_already_open(filepath) result(already_open)
    character(len=*), intent(in) :: filepath
    logical :: already_open
    inquire(file=filepath, opened=already_open)
  end function file_already_open

  function file_is_readable(filepath) result(is_readable)
    character(len=*), intent(in) :: filepath
    character(len=7) :: readable
    logical :: is_readable
    inquire(file=filepath, read=readable)
    is_readable = merge(.true., .false., readable == 'YES')
  end function file_is_readable

  function file_is_writable(filepath) result(is_writable)
    character(len=*), intent(in) :: filepath
    character(len=7) :: writable
    logical :: is_writable
    inquire(file=filepath, write=writable)
    is_writable = merge(.true., .false., writable == 'YES')
  end function file_is_writable

end module file
