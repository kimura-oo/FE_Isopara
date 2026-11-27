include '../src/mesh_elevation.f90'
program mesh_p1_to_p2
    use mesh_elevation, only: convert_mesh
    implicit none
    call convert_mesh(2)
end program
