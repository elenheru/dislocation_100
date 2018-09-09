
!MS$DECLARE

module array_parameters_mod
    save
    integer, parameter  ::  x_layers= 18,y_layers= 18,z_layers= 4
    integer, parameter  ::  &
    atoms_max_array = (2*x_layers+1)*(2*y_layers+1)*(2*z_layers+1)*3
endmodule array_parameters_mod

module comp_parameters_mod
    use array_parameters_mod
    save
    integer, parameter  ::  beauty_denominator=72 !number of shots during iterative wall shift refinement
    integer             ::  vestaWOcounter=0
    integer             ::  xyz__WOcounter=0
    integer             ::  matrixWOcounter=0
    integer             ::  burgers_poisson_WOcounter=0
    integer, parameter  ::  total_relax_passages    = 16!10000 !
    integer, parameter  ::  anisotropy_passages     = 8!10000 !
    !integer, parameter  ::  x_layers= 11,y_layers= 11,z_layers= 4
    !!@layers are 1+2*layers. see spawn_bcc_rectangular_100 for details
    integer, parameter  ::  x_edges = 2,y_edges = 2,z_edges = 1!@from 4 edges of cell how much of a0 are borders
    integer             ::  atoms__in_total
    integer             ::  first_relaxable,last__relaxable
    integer             ::  first_z__period,last__z__period
    integer             ::  first______wall,last_______wall
    real(8), parameter  ::  pi = 3.14159265d0
endmodule comp_parameters_mod

module phys_parameters_mod
    !use comp_parameters_mod
    save
    real(8), parameter  ::  a0 = 2.8600d0
    real(8), parameter  ::  poisson=293d-3
    real(8), parameter  ::  burgers=-a0*50d-2*2d0!dislocation is doubled
    real(8), parameter  ::  core_sign= 1d0
    !core is compressed if sign=1 or decomressed if sign=-1
endmodule phys_parameters_mod

module      positions_mod   !angstrems
    use array_parameters_mod
    save
    real(8), dimension(1:3,atoms_max_array) ::  R_perf,R_curr
endmodule   positions_mod

module      periodic_conditions_mod
    use array_parameters_mod
    save
    integer             ::  periodic_and_body(1:2,atoms_max_array/4)
    integer             ::  max_periodic_pair
endmodule   periodic_conditions_mod

module      interaction_mod
    use array_parameters_mod
    save
!    real(8), parameter ::   cutoff = 5.6_8 !sphere where the neighbor is
    real(8), parameter ::   cutoff = 3.5d0 !sphere where the neighbor is
    integer, parameter ::   neibors_max = 150
    integer verlet_list(0:neibors_max,atoms_max_array)
    real(8) distan_list(  neibors_max,atoms_max_array)
endmodule   interaction_mod

module      cgm_storage_mod
    use array_parameters_mod
    save
    real(8), dimension(1:3,atoms_max_array) ::  F_curr,D_prev,D_curr!,F_prev
    real(8), parameter  ::  system_gradient_small_enough=1d-9
    real(8), parameter  ::  cgm_denominator             =1d0/b"10"
    integer, parameter  ::  cgm_direction_switch        =19
endmodule   cgm_storage_mod

module      matrixes_mod
    use array_parameters_mod
    save
    integer,dimension(-x_layers:x_layers,-y_layers:y_layers) ::&
    matrix_z00
    integer,dimension(-x_layers:x_layers-1,-y_layers:y_layers-1) ::&
    matrix_z05
    !!atoms are in plane with z=const z=0 or z=0.5*a0
    real(8),dimension(-2*x_layers:2*x_layers,-2*y_layers:2*y_layers) ::&
    matrix_united_ux,matrix_united_uy,matrix_united_uz
    !matrix_united_ux(2*i  ,2*j  )=ux(matrix_z00(i,j))
    !matrix_united_ux(2*i+1,2*j+1)=ux(matrix_z05(i,j))
endmodule   matrixes_mod

module      anisotropy_mod
    use array_parameters_mod
    save
    integer ::  max_layer_list_index
    real(8) ::  burgers_list(atoms_max_array/4)
    real(8) ::  poisson_list(atoms_max_array/4)
    real(8) ::  phi_____list(atoms_max_array/4)
    integer ::  layer_list(  atoms_max_array/4)
endmodule   anisotropy_mod
