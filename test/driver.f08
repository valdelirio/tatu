program fruit_driver
  use fruit
  use test_utils

  call init_fruit(1)
  call fruit_show_dots


  call reset_unit_name
  call utils_int2str



  call fruit_summary
  call fruit_finalize
end program fruit_driver
