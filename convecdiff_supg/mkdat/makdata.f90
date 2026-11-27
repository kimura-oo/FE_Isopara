include '../../tools/supg_fields.f90'

program makdata
    use supg_fields, only: prepare_fields
    implicit none
    call prepare_fields('right','zero','open_right')
end program
