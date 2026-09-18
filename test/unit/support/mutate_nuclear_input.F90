Program mutate_nuclear_input
  Implicit None

  Character(32) :: mode
  Character(512) :: data_dir

  If ( command_argument_count() /= 2 ) Then
    Write(*,*) 'usage: mutate_nuclear_input MODE DATA_DIR'
    Stop 2
  EndIf
  Call get_command_argument(1,mode)
  Call get_command_argument(2,data_dir)

  Select Case (trim(mode))
  Case ('nets4-count','nets4-order','nets4-order-tail')
    Call mutate_nets4(trim(data_dir),trim(mode))
  Case ('match-header')
    Call empty_match_data(trim(data_dir))
  Case ('match-count-1')
    Call mutate_match_data(trim(data_dir),1)
  Case ('match-count-2')
    Call mutate_match_data(trim(data_dir),2)
  Case ('match-count-3')
    Call mutate_match_data(trim(data_dir),3)
  Case ('match-count-4')
    Call mutate_match_data(trim(data_dir),4)
  Case Default
    Write(*,*) 'unsupported nuclear-input mutation: ',trim(mode)
    Stop 2
  End Select

Contains

  Subroutine mutate_nets4(directory,mutation)
    Implicit None

    Character(*), Intent(in) :: directory, mutation

    Character(5), Allocatable :: names(:)
    Character(5) :: saved_name
    Integer :: lun, ny_file

    Open(newunit=lun,file=trim(directory)//'/nets4',form='unformatted',status='old',action='read')
    Read(lun) ny_file
    If ( mutation == 'nets4-count' ) Then
      Close(lun)
      Open(newunit=lun,file=trim(directory)//'/nets4',form='unformatted',status='replace',action='write')
      Write(lun) ny_file + 1
      Close(lun)
      Return
    EndIf

    Allocate (names(ny_file))
    Read(lun) names
    Close(lun)
    If ( ny_file < 2 ) Then
      Write(*,*) 'nets4 order mutation requires at least two nuclei'
      Stop 2
    EndIf
    If ( mutation == 'nets4-order' ) Then
      saved_name = names(1)
      names(1) = names(2)
      names(2) = saved_name
    Else
      names(ny_file) = 'ne20'
    EndIf
    Open(newunit=lun,file=trim(directory)//'/nets4',form='unformatted',status='replace',action='write')
    Write(lun) ny_file
    Write(lun) names
    Close(lun)
    Deallocate (names)

    Return
  End Subroutine mutate_nets4

  Subroutine empty_match_data(directory)
    Implicit None

    Character(*), Intent(in) :: directory

    Integer :: lun

    Open(newunit=lun,file=trim(directory)//'/match_data',form='unformatted',status='replace',action='write')
    Close(lun)

    Return
  End Subroutine empty_match_data

  Subroutine mutate_match_data(directory,group)
    Implicit None

    Character(*), Intent(in) :: directory
    Integer, Intent(in) :: group

    Integer :: lun, mflx_file, nr(4)

    Open(newunit=lun,file=trim(directory)//'/match_data',form='unformatted',status='old',action='read')
    Read(lun) mflx_file, nr
    Close(lun)
    nr(group) = nr(group) + 1
    Open(newunit=lun,file=trim(directory)//'/match_data',form='unformatted',status='replace',action='write')
    Write(lun) mflx_file, nr
    Close(lun)

    Return
  End Subroutine mutate_match_data

End Program mutate_nuclear_input
