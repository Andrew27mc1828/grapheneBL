classdef compute_lattice < handle

    properties
        L_sizeX;
        L_sizeY;
        H_dim;

        twist_angle_rad;        
        layer_seperation;  
        default_E;
        cutoff_E;

        Hamiltonian_SL;
        Hamiltonian_BL;

        index_array_L1;
        particle_pos_L1;
        particle_pos_L2;

        L1_center;
        L2_center;        

        rotation_matrix;
        rotation_matrix_I;
    end

    methods (Access = public)
        function obj = compute_lattice( ...
                lattice_dimX, lattice_dimY, twist_angle, distance_between_layers, ...
                bond_energy, cutoff_energy)
            
            obj.L_sizeX = lattice_dimX;
            obj.L_sizeY = lattice_dimY;

            basis_num = obj.L_sizeX * obj.L_sizeY;
            obj.H_dim = 4 * basis_num;

            obj.twist_angle_rad = twist_angle;
            obj.layer_seperation = distance_between_layers;
            obj.default_E = bond_energy;
            obj.cutoff_E = cutoff_energy;

            fprintf('Single Layer Hamiltonian dimension: %i \n', obj.H_dim);

        end

        function build_H(obj)

            obj.index_bases();      % Make an index array for tracking the bases
            obj.build_HSL();        % Build a Hamiltonian for a single layer
            obj.build_HBL();        % Build a Hamiltonian for a bilayer
        
        end
    end

    methods (Access = private)
        
        function index_bases( obj )
            
            % Index everything
            obj.index_array_L1 = zeros(obj.L_sizeX, obj.L_sizeY);

            counter = 1;
            for i = 1:obj.L_sizeX
                for j = 1:obj.L_sizeY
                    obj.index_array_L1(i,j) = counter;
                    counter = counter + 4;
                end
            end

            % Keep track of particle positions
            obj.particle_pos_L1 = zeros( obj.H_dim, 3);
            obj.particle_pos_L2 = zeros( obj.H_dim, 3);
            
            baseA = [0.0; 0.0; 0.0];
            baseB = [1.0; 0.0; 0.0];
            baseC = [1.5; 0.5*sqrt(3); 0.0];
            baseD = [2.5; 0.5*sqrt(3); 0.0];
            
            aX = [3.0; 0.0; 0.0];
            aY = [0.0; sqrt(3); 0.0];
            
            obj.L1_center = [0.0; 0.0; -0.5 * obj.layer_seperation];
            obj.L2_center = [0.0; 0.0; 0.5 * obj.layer_seperation];
            
            obj.L1_center = obj.L1_center - 0.5*( obj.L_sizeX * aX ) - 0.5*( obj.L_sizeY * aY );
            obj.L2_center = obj.L2_center - 0.5*( obj.L_sizeX * aX ) - 0.5*( obj.L_sizeY * aY );
            
            obj.rotation_matrix = [ cos(obj.twist_angle_rad),   -sin(obj.twist_angle_rad),  0.0;
                                    sin(obj.twist_angle_rad),   cos(obj.twist_angle_rad),   0.0;
                                    0.0,                        0.0,                        1.0 ];
            
            obj.rotation_matrix_I = inv(obj.rotation_matrix);

            counter = 1;
            for i = 1:obj.L_sizeX
                for j = 1:obj.L_sizeY
                    
                    center = ((i - 1) * aX) + ((j - 1) * aY) + obj.L1_center;
                    obj.particle_pos_L1( counter    , : ) = center + baseA;
                    obj.particle_pos_L1( counter + 1, : ) = center + baseB;
                    obj.particle_pos_L1( counter + 2, : ) = center + baseC;
                    obj.particle_pos_L1( counter + 3, : ) = center + baseD;
            
                    center = ((i - 1) * aX) + ((j - 1) * aY) + obj.L2_center;
                    obj.particle_pos_L2( counter    , : ) = obj.rotation_matrix * (center + baseA);
                    obj.particle_pos_L2( counter + 1, : ) = obj.rotation_matrix * (center + baseB);
                    obj.particle_pos_L2( counter + 2, : ) = obj.rotation_matrix * (center + baseC);
                    obj.particle_pos_L2( counter + 3, : ) = obj.rotation_matrix * (center + baseD);
            
                    counter = counter + 4;
                end
            end
        end

        function build_HSL( obj )

            obj.Hamiltonian_SL = zeros( obj.H_dim, obj.H_dim );

            % Construct the Single Layer Hamiltonian
            for i = 1:obj.L_sizeX
                for j = 1:obj.L_sizeY
                    
                    particle1 = obj.index_array_L1(i, j);
                    
                    % Connect particles within (i,j) basis
                    obj.Hamiltonian_SL( particle1,     particle1 + 1 ) = obj.default_E;
                    obj.Hamiltonian_SL( particle1 + 1, particle1 + 2 ) = obj.default_E;
                    obj.Hamiltonian_SL( particle1 + 2, particle1 + 3 ) = obj.default_E;
            
                    % Plus h.c.
                    obj.Hamiltonian_SL( particle1 + 1, particle1 )     = obj.default_E;
                    obj.Hamiltonian_SL( particle1 + 2, particle1 + 1 ) = obj.default_E;
                    obj.Hamiltonian_SL( particle1 + 3, particle1 + 2 ) = obj.default_E;
            
                    % Connect to other bases
                    if i ~= 1
                        particle2 = obj.index_array_L1(i - 1, j);            
                        obj.Hamiltonian_SL( particle1,     particle2 + 3 ) = obj.default_E;
                        obj.Hamiltonian_SL( particle2 + 3, particle1     ) = obj.default_E; % Plus h.c.
                
                        if j ~= 1
                            particle2 = obj.index_array_L1(i - 1, j - 1);
                            obj.Hamiltonian_SL( particle1    , particle2 + 3 ) = obj.default_E;
                            obj.Hamiltonian_SL( particle2 + 3, particle1     ) = obj.default_E; % Plus h.c.
                        end
                    end
            
                    if (i ~= obj.L_sizeX) && (j ~= obj.L_sizeY)
                        particle2 = obj.index_array_L1(i + 1, j + 1);
                        obj.Hamiltonian_SL( particle1 + 3, particle2     ) = obj.default_E;
                        obj.Hamiltonian_SL( particle2    , particle1 + 3 ) = obj.default_E; % Plus h.c.
                    end
            
                    if j ~= 1
                        particle2 = obj.index_array_L1(i, j - 1);
                        obj.Hamiltonian_SL( particle1 + 1, particle2 + 2 ) = obj.default_E;
                        obj.Hamiltonian_SL( particle2 + 2, particle1 + 1 ) = obj.default_E; % Plus h.c.
                    end
            
                end
            end

            % Done making SL Hamiltonian
        end
    
        function build_HBL( obj )
            
            obj.Hamiltonian_BL = zeros( 2*obj.H_dim, 2*obj.H_dim );

            % Direct sum H_SL
            for i = 1:obj.H_dim
                for j = 1:obj.H_dim
            
                    obj.Hamiltonian_BL(i,j) = obj.Hamiltonian_SL(i,j);
                    obj.Hamiltonian_BL(i + obj.H_dim, j + obj.H_dim) = obj.Hamiltonian_SL(i,j);
                end
            end
               
            % Compute inter-layer interactions
            for i = 1:obj.L_sizeX
                for j = 1:obj.L_sizeY
                    
                    particle1_index = obj.index_array_L1(i, j);
                    
                    obj.connect_particle_to_layer2( particle1_index );
                    obj.connect_particle_to_layer2( particle1_index + 1 );
                    obj.connect_particle_to_layer2( particle1_index + 2 );
                    obj.connect_particle_to_layer2( particle1_index + 3 );
                end
            end

            % Done building BL Hamiltonian
        end

        function connect_particle_to_layer2(obj, p1_index)
            
            particle1_pos_bA = obj.particle_pos_L1( p1_index, : );
                    
            particle2_pos_bA_approx = obj.rotation_matrix_I * transpose( particle1_pos_bA );
            particle2_c0 = particle2_pos_bA_approx - obj.L2_center;
            
            i2 = round( particle2_c0(1) / 3.0 ) + 1;
            j2 = round( particle2_c0(2) / sqrt(3) ) + 1;
            
            obj.connect_base(i2, j2, p1_index);
            obj.connect_base(i2 + 1, j2 + 1, p1_index);
            obj.connect_base(i2 + 1, j2 - 1, p1_index);
            obj.connect_base(i2 - 1, j2 + 1, p1_index);
            obj.connect_base(i2 - 1, j2 - 1, p1_index);
            
        end

        function connect_base( obj, i_L2, j_L2, p1_index )

            if i_L2 < 1
                i_L2 = 1;
            end
            
            if j_L2 < 1
                j_L2 = 1;
            end
            
            if i_L2 > obj.L_sizeX 
                i_L2 = obj.L_sizeX;
            end
            
            if j_L2 > obj.L_sizeY 
                j_L2 = obj.L_sizeY;
            end

            p2_indexA = obj.index_array_L1(i_L2, j_L2);
            
            obj.connect_index( p1_index, p2_indexA );
            obj.connect_index( p1_index, p2_indexA + 1 );
            obj.connect_index( p1_index, p2_indexA + 2 );
            obj.connect_index( p1_index, p2_indexA + 3 );
            
        end

        function connect_index(obj, p1_index, p2_index ) 
            
            E_interaction = obj.Energy( ...
                obj.particle_pos_L1( p1_index, : ), ...
                obj.particle_pos_L2( p2_index, : ) ); 
    
            if E_interaction > obj.cutoff_E
                obj.Hamiltonian_BL(p1_index, p2_index + obj.H_dim) = E_interaction;
                obj.Hamiltonian_BL(p2_index + obj.H_dim, p1_index) = E_interaction;
            end
        end

        function E = Energy( obj, r1, r2 )
            E = obj.default_E * exp( 1.2 * ( 1.0 - norm(r1 - r2) ) );
        end
    
    end
end