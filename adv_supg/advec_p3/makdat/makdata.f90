include '../../../tools/supg_fields.f90'

program makdata
    use supg_fields, only: prepare_fields
    implicit none
    call prepare_fields('down','box','all')
end program
