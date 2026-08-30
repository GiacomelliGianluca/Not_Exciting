function  [A, B, C, D] = system_dynamics(sys_id)
    % Discrete-time system simulation (handles input as a sequence)
    % Inputs:
    %   x0 - Initial state (n×1)
    %   u_seq - Control input sequence (1×T)
    % Outputs:
    %   x_seq - State sequence (n×T)
    %   y_seq - Output sequence (p×T)
    % n: system order
    % System matrices
      switch sys_id
        case 1
            % system 1
            A = [0.7326, -0.0861;
                 0.1722,  0.9909];
            B = [0.0609;
                 0.0064];
            C = [0, 1.4142];

            D = 0;

        case 2
            % system 2
            A = [0.85, 0.30;
                 -0.20, 0.75];
            B = [1.0;
                 0.2];
            C = [0.0, 1.0];

            D = 0;

        case 3
            % system 3
            A = [0.9000, -0.4000;
                 0.6000,  0.7000];
            B = [0.0000;
                 1.0000];
            C = [1.0000, 0.0000];
            D = 0;

            
        case 4
            % system 4
            A = [  2.20632892, -2.29641076,  1.49610313, -0.5184;
                    1           0           0           0;
                    0           1           0           0;
                    0           0           1           0 ];
            B = [1; 0; 0; 0];
            C = [ 0, 1, -1.17557050, 1 ];
            D = 0;


            
      
            
       case 5
            % system 5
            A = [  1.01192665, -0.31215078,  0.33767358, -0.50765625;
                     1           0           0           0;
                     0           1           0           0;
                     0           0           1           0 ];
            B = [1; 0; 0; 0];
            C = [ 0, 1, -0.21900825, 0.49 ];

            D = 0;


            
        otherwise
            error('invalid sys_id');
     end


end
