
addpath compute_lattice.m
addpath display_lattice.m

CL = compute_lattice( ...
    6, ...          % lattice size in x
    6, ...          % lattice size in y
    0.2000, ...     % bilayer relative twisting angle (in rad)
    2.36, ...       % bilayer total seperation distance between layers (in units of the bond length)
    2.7, ...        % bond energy
    0.39);          % cutoff energy

DL = display_lattice( ...
    -0.2, ...       % Energy range min                                                       
    0.5, ...        % Energy range max
    1.0e-6, ...     % Imaginary component of the energy (small for fine detail)
    200 );          % Energy plot resolution

CL.build_H();

DL.show_Hamiltonian_SL(CL);     % Visualize the Hamiltonian for the single layer and bilayer cases
DL.show_Hamiltonian_BL(CL); 
DL.trace_GreensF_SL(CL);        % Plot Im tr G for both cases
DL.trace_GreensF_BL(CL);
DL.show_lattice3D(CL);          % Plot the lattice in 3D (disable for large lattices)

% TO DO: 
% - Generate particle_pos_L1, particle_pos_L2 (in CL) procedurally, no need
% to be stored in memory.
% - Implement curve fitting to reduce the # of times G must be inverted.
% - Arbitrary lattice center
