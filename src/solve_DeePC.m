

function [g_opt, y_pred,u_opt, inf] = solve_DeePC(U_p, Y_p, U_f, Y_f, u_ini, y_ini, r, Q, R, flag_SPC, P_SPC)



    % minus_Pi: matrix for subspace regularization
    % lambda_Pi: hyperparameters for subspace regularization


    % Sizes
    n_uf = size(U_f, 1); % Number of future inputs
    n_yf = size(Y_f, 1); % Number of future outputs
    n_g = size(U_p, 2); % Number of selector elements
    ny_p = size(Y_p, 1); % Past-output stacked length


   inf = 0; % infeasibility flag



try
 cvx_begin quiet

        % Declare the optimization variables
        variable u(n_uf,1)
        variable y(n_yf,1)
        variable g(n_g)

        % Cost function initialization
        J = 0;
        
        for ii=1:length(r)
            
            e = y(ii) - r(ii);
            J = J + e'*Q*e + u(ii)'*R*u(ii);

        end

        minimize(J)
        subject to
            
        if flag_SPC ~=1
            % Model constraints
            U_p*g==u_ini;
            U_f*g==u;
            Y_p*g == y_ini;
            Y_f*g==y;
        end

        if flag_SPC == 1
           y == P_SPC*[y_ini; u_ini; u];
        end

            u>=-2;
            u<=2;
    cvx_end
   
catch
    inf = 1;
    % keyboard;
end

    % Extract the optimal control action
    g_opt = g;
    u_opt = u;
    y_pred = y;   
    % fprintf('g = [');
    % fprintf(' %.6e', g_opt);
    % fprintf(' ]\n');




end
