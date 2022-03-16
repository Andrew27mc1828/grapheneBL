classdef display_lattice < handle
    
    properties
        E0;
        E1;
        res;

        energy_res;
    end

    methods (Access = public)
        
        function obj = display_lattice( E_min, E_max, energy_resolution, plot_resolution )
            obj.E0 = E_min;
            obj.E1 = E_max;
            obj.energy_res = energy_resolution;
            obj.res = plot_resolution;
        end

        function show_Hamiltonian_SL(obj, compL)
            
            figure(1);
            spy(compL.Hamiltonian_SL,'b');
            grid on;
            title('Hamiltonian Plot Single Layer');
        end

        function show_Hamiltonian_BL(obj, compL)
            
            figure(2);
            spy(compL.Hamiltonian_BL,'b');
            grid on;
            title('Hamiltonian Plot Bilayer');
        end

        function trace_GreensF_SL(obj, compL)

            tic; % Time this precess
            [x3, y3] = obj.trace_GreensF(compL.Hamiltonian_SL, compL.H_dim);
            
            fprintf('elapsed time for SL solution is: %.2f seconds. \n', toc);
            
            figure(3);
            plot(x3, y3);
            ylim( [-0.05, 0.2] );
            grid on;
            title('Signle Layer Energy Plot');
            xlabel('Energy (eV)');
            ylabel('Im[tr(G)]');
        end

        function trace_GreensF_BL(obj, compL)

            tic; % Time this precess
            [x4, y4] = obj.trace_GreensF(compL.Hamiltonian_BL, 2*compL.H_dim);
            
            fprintf('elapsed time for BL solution is: %.2f seconds. \n', toc);
            
            figure(4);
            plot(x4, y4);
            ylim( [-0.05, 0.2] );
            grid on;
            title('Bilayer Energy Plot');
            xlabel('Energy (eV)');
            ylabel('Im[tr(G)]');
        end

        function show_lattice3D(obj, compL)

            aX = [3.0; 0.0; 0.0];
            aY = [0.0; sqrt(3); 0.0];
            
            figure(5);
            
            plot3( ...
                compL.particle_pos_L1(:,1), ...
                compL.particle_pos_L1(:,2), ...
                compL.particle_pos_L1(:,3), ...
                '.', 'Color', 'black', 'MarkerSize', 12 );
            hold on;
            
            plot3( ...
                compL.particle_pos_L2(:,1), ...
                compL.particle_pos_L2(:,2), ...
                compL.particle_pos_L2(:,3), ...
                '.', 'Color', 'black', 'MarkerSize', 12 );
            hold on;
            
            zlim( [-0.8 * compL.layer_seperation, 0.8 * compL.layer_seperation] );
            grid on;
            
            % L1 bounding box
            bbox_corner1 = ( compL.L_sizeX * aX ) + ( compL.L_sizeY * aY );
            bbox_corner2 = ( compL.L_sizeX * aX ) - ( compL.L_sizeY * aY );
            
            bbox_x = 0.5 * [ bbox_corner1(1), bbox_corner2(1), -bbox_corner1(1), -bbox_corner2(1), bbox_corner1(1) ];
            bbox_y = 0.5 * [ bbox_corner1(2), bbox_corner2(2), -bbox_corner1(2), -bbox_corner2(2), bbox_corner1(2) ];
            bbox_z = -0.5 * [ compL.layer_seperation, compL.layer_seperation, compL.layer_seperation, compL.layer_seperation, compL.layer_seperation ];
            
            plot3(bbox_x, bbox_y, bbox_z, 'Color', '#77AC30');
            hold on;
            
            % L2 bounding box
            bbox_corner1 = compL.rotation_matrix * (( compL.L_sizeX * aX ) + ( compL.L_sizeY * aY ));
            bbox_corner2 = compL.rotation_matrix * (( compL.L_sizeX * aX ) - ( compL.L_sizeY * aY ));
            
            bbox_x = 0.5 * [ bbox_corner1(1), bbox_corner2(1), -bbox_corner1(1), -bbox_corner2(1), bbox_corner1(1) ];
            bbox_y = 0.5 * [ bbox_corner1(2), bbox_corner2(2), -bbox_corner1(2), -bbox_corner2(2), bbox_corner1(2) ];
            bbox_z = 0.5 * [ compL.layer_seperation, compL.layer_seperation, compL.layer_seperation, compL.layer_seperation, compL.layer_seperation ];
            
            plot3(bbox_x, bbox_y, bbox_z, 'Color', '#77AC30');
            hold on;
            
            for m = 1:compL.H_dim
                for n = 1:compL.H_dim
                    
                    if compL.Hamiltonian_SL(m,n) ~= 0.0
            
                        % Plot connections in L1
                        edge_x = [ compL.particle_pos_L1( m, 1 ), compL.particle_pos_L1( n, 1 ) ];
                        edge_y = [ compL.particle_pos_L1( m, 2 ), compL.particle_pos_L1( n, 2 ) ];
                        edge_z = [ compL.particle_pos_L1( m, 3 ), compL.particle_pos_L1( n, 3 ) ];
            
                        plot3( edge_x, edge_y, edge_z, 'Color', '#0072BD' );
                        hold on;
            
                        % Plot connections in L2
                        edge_x = [ compL.particle_pos_L2( m, 1 ), compL.particle_pos_L2( n, 1 ) ];
                        edge_y = [ compL.particle_pos_L2( m, 2 ), compL.particle_pos_L2( n, 2 ) ];
                        edge_z = [ compL.particle_pos_L2( m, 3 ), compL.particle_pos_L2( n, 3 ) ];
            
                        plot3( edge_x, edge_y, edge_z, 'Color', '#0072BD' );
                        hold on;
            
                    end
            
                    % Plot inter-layer connections
                    if compL.Hamiltonian_BL(m, n + compL.H_dim) ~= 0.0
                        
                        edge_x = [ compL.particle_pos_L1( m, 1 ), compL.particle_pos_L2( n, 1 ) ];
                        edge_y = [ compL.particle_pos_L1( m, 2 ), compL.particle_pos_L2( n, 2 ) ];
                        edge_z = [ compL.particle_pos_L1( m, 3 ), compL.particle_pos_L2( n, 3 ) ];
            
                        pl_edge = plot3( edge_x, edge_y, edge_z, 'Color', '#0072BD' );
                        hold on;
                    end
                end
            end
            
            % Done showing lattice 3D
        end

    end

    methods (Access = private)
         
        function [x,y] = trace_GreensF(obj, H, dim)

            energy_scale = linspace(obj.E0, obj.E1, obj.res);
            dos = energy_scale;
            
            for n = 1:length(energy_scale)
                
                E_real = energy_scale(n);
                lambda = complex(E_real, obj.energy_res);
                G = inv( H - (lambda * eye( dim, dim )) );
                dos(n) = imag(trace(G));
            end

            x = energy_scale;
            y = dos;
        end
    end
end