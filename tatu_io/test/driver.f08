program fruit_driver
  use fruit
  use test_file
  use test_types
  use test_input
  use test_output
  call init_fruit(1)


  call test_file_already_exist
  call test_file_already_open
  call test_file_is_readable
  call test_file_is_writable

  call test_types_point
  call test_types_transmitter
  call test_types_receiver

  call test_input_transmitter_model
  call test_input_direction
  call test_input_step
  call test_input_get_type


  call fruit_summary
  call fruit_finalize
end program fruit_driver
