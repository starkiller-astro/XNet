Program preprocess_direct
  Use, Intrinsic :: iso_fortran_env, Only: output_unit
  Use xnet_preprocess, Only: net_preprocess
  Implicit None

  Character(80) :: data_desc
  Character(256) :: data_dir

  If ( command_argument_count() /= 2 ) Then
    Write(output_unit,*) 'usage: preprocess_direct DATA_DIR DESCRIPTION'
    Stop 1
  EndIf
  Call get_command_argument(1,data_dir)
  Call get_command_argument(2,data_desc)
  Call net_preprocess(output_unit,trim(data_dir),data_desc)

End Program preprocess_direct
