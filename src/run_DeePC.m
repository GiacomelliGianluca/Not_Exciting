function [u_seq, y_seq,x_seq,y_pred_seq, is_infeasible, u_ini_0, y_ini_0] = run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R, x0, T, sys, flag_SPC, randseed) 
rng(randseed)
    % Runs the DeePC algorithm in a receding horizon manner.
    % Inputs:
    %   U_p, Y_p, U_f, Y_f - Hankel matrices
    %   desired_trajectory - Full reference trajectory
    %   T_ini - Initial condition length
    %   N - Prediction horizon
    % -  Q, R - Performance weight matrices
    % - T           - Total simulation steps
    % - flag_SPC: flag to use SPC instead 
    % Output:
    %   u_seq - Optimized control sequence over time
    
    % Initialize input and output sequences
    u_seq = zeros(T, 1);  % Stores the control inputs applied
    nx=size(x0,1);

    y_seq = zeros(T, 1);  % Stores the measured system outputs
    x_seq = zeros(nx, T+1);  % Stores the system states
    y_pred_seq = zeros(T, N);     


   % 1. Initialize past input and output data
   u_ini=2*rand(1, T_ini)' -1;
   y_ini = zeros(T_ini, 1);  


    % Compute initial y_ini from system
    x_tmp = x0;                 
    for i = 1:T_ini
        [y_ini(i), x_tmp] = system_dynamics_global_no_noise(u_ini(i), x_tmp,sys);
    end
    x_seq(:,1) = x_tmp; 
    u_ini = u_ini(end - (T_ini) + 1 : end);
    y_ini = y_ini(end - (T_ini) + 1 : end);


    u_ini_0 = u_ini;
    y_ini_0 = y_ini;
 

    if flag_SPC
        P_SPC = Y_f*pinv([Y_p;U_p; U_f]);
    else
        P_SPC = [];
    end

    % 2. Control loop (rolling horizon optimization)
    for t = 1:T
        if t > 1
            u_ini = [u_ini(2:end); u_seq(t-1)];
            y_ini = [y_ini(2:end); y_seq(t-1)];
        end
        
        r = desired_trajectory(t:t+N-1);
        
        [~, y_pred, u_opt, is_infeasible] = solve_DeePC(U_p, Y_p, U_f, Y_f, u_ini, y_ini, r, Q, R, flag_SPC, P_SPC);
        
        if is_infeasible
            warning('DeePC infeasibile at t = %d. Break.\n', t);
            u_seq = u_seq(1:t-1);
            y_seq = y_seq(1:t-1);
            x_seq = x_seq(:, 1:t);
            break;
        end
        
        u_seq(t) = u_opt(1); 
        y_pred_seq(t,:) = y_pred';
        
        [y_seq(t), x_next] = system_dynamics_global_no_noise(u_seq(t), x_seq(:,t), sys);
        x_seq(:,t+1) = x_next; 
    end
    y_seq = y_seq';


    % Return the computed control input sequence
    disp('DeePC optimization completed.');
end






