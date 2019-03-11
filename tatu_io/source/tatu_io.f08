module tatu_io
  use json_module
  use types
  use input
  use output

contains

  function tatu_io_read_input(input_file)
    character(len=*), intent(in) :: input_file
    type(json_input) :: tatu_io_read_input
    call tatu_io_get_input(input_file, tatu_io_read_input)
  end function tatu_io_read_input

  subroutine tatu_io_write_output(output_file, in, labels, values, u_transmitter, u_frequency, u_receiver)
    character(len=*), intent(in) :: output_file
    type(json_input), intent(in) :: in
    character(len=*), dimension(:), intent(in) :: labels
    real(real_dp), dimension(:,:), intent(in) :: values
    real(real_dp), dimension(:,:), intent(in) :: u_transmitter
    real(real_dp), dimension(:), intent(in) :: u_frequency
    real(real_dp), dimension(:,:), intent(in) :: u_receiver
    call output_write(output_file, in, labels, values, u_transmitter, u_frequency, u_receiver)
  end subroutine tatu_io_write_output

  subroutine tatu_io_write_output_ssv(output_file, labels, output_data)
    use, intrinsic :: iso_fortran_env, only: error_unit
    character(len=*), intent(in) :: output_file
    character(len=*), dimension(:), intent(in) :: labels
    real(real_dp), dimension(:,:), intent(in) :: output_data
    integer :: i, unit_file, write_status

    10 format( 15(G24.15E3) )
    open(newunit = unit_file, file = output_file//'.ssv', status = 'replace', action = 'write', iostat = write_status)
    if (write_status /= 0) write(error_unit,'(a)') 'On save: '//output_file
    write(unit_file,'(15a)') (/ (labels(i), i=1, size(labels)) /)
    write(unit_file,10)transpose(output_data)
    close(unit_file)
  end subroutine tatu_io_write_output_ssv

end module tatu_io
