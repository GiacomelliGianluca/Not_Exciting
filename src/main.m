clc
clear
close all

%% Parameters
T_D   = 200;     % Dataset dimension
Ts    = 0.01;    % Sampling time
fs    = 1/Ts;    % Sampling frequency
% Hankel parameters
T_ini = 8;              % Past horizon
N     = 15;             % Future horizon
L     = T_ini + N;
T     = 10*T_D;         % Extended data collection to remove transient (only the last T_D samples will be used)
T_c   = 100;           % Control period
% DeePC Weights
Q = 100;
R = 1;
flag_SPC = 1; % 1: use SPC instead
% Target reference signal frequency [Hz]
f_target = 2;     
MC = 2; % Number of MC simulations
data_dir = fullfile(pwd, 'data_pred_15/SPC_test');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end

%% Monte Carlo loop
for sysnum = 1:5
    if ismember(sysnum, [1, 2, 3])
        dim = 2;
    else
        dim = 4;
    end
    x0 = zeros(dim,1);
    for si = 0:6 
        fname_big = sprintf('s%02d_r%02d_case%02d.mat', sysnum, f_target, si);
        fname_big = fullfile(data_dir, fname_big);
        
        % ======= (only used by multisine) =======
        k    = 0:(T-1);     % Time index 
        Nper = T_D;        
        f0   = fs / Nper;   % f0 = fs/N
        rng(1);
        
        switch si
            case 0 % WN
                rng(1);
                u_seq = randn(1, T);
            case 1 % IB
                Jmax = 51; 
                u_seq = zeros(size(k));
                for j = 1:Jmax
                    phi_j = -pi * j*(j-1) / Jmax;  
                    A_j   = 1;                     
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j );
                end
            case 2 % IBN
                Jmax = 37;               
                u_seq = zeros(size(k));
                for j = 1:Jmax
                    phi_j = -pi * j*(j-1) / Jmax;  
                    A_j   = 1;                      
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j );
                end
            case 3 % OB
                j_min = 30;
                j_max = 80;
                Jnum  = j_max - j_min + 1;
                u_seq = zeros(size(k));
                idx = 0;
                for j = j_min:j_max
                    idx = idx + 1;
                    phi_j = -pi * idx*(idx-1) / Jnum;  
                    A_j   = 1;                         
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j );
                end
            case 4 % OBN 
                j_min = 40;
                j_max = 76;
                Jnum  = j_max - j_min + 1;

                u_seq = zeros(size(k));
                idx = 0;
                for j = j_min:j_max
                    idx = idx + 1;
                    phi_j = -pi * idx*(idx-1) / Jnum;  
                    A_j   = 1;                         
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j );
                end
             case 5 % IBW
                Jmax = 61;                
                u_seq = zeros(size(k));
                for j = 1:Jmax
                    phi_j = -pi * j*(j-1) / Jmax;  
                    A_j   = 1;                     
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j );
                end
            case 6 % OBW
                j_min = 35;
                j_max = 95;
                Jnum  = j_max - j_min + 1;
                u_seq = zeros(size(k));
                idx = 0;
                for j = j_min:j_max
                    idx = idx + 1;
                    phi_j = -pi * idx*(idx-1) / Jnum;  
                    A_j   = 1;                         
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j );
                end
            otherwise
                error('Unknown si = %d', si);
        end
        
        [y_seq, ~, n] = system_dynamics_no_noise(u_seq, x0, sysnum);
        u_seq_cut = u_seq(end - T_D + 1 : end);
        y_seq_cut = y_seq(end - T_D + 1 : end);
        desired_trajectory = sin( f_target * 2*pi / fs * (0:(T-1)) );
        [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_seq_cut, y_seq_cut, T_ini, N, dim);
        H_u_ext = hankel(u_seq_cut(1:L+n), u_seq_cut(L+n:length(u_seq_cut)));
        sv_H_u_ext = svd(H_u_ext);
        sv_H_w = svd([U_p; U_f; Y_p; Y_f]);
 
        run_prefix = sprintf('s%02d_r%02d_case%02d_run', sysnum, f_target, si);
        last_run   = 0;
        big        = struct();
        
        if isfile(fname_big)
            big = load(fname_big);
            fn  = fieldnames(big);
            run_nums = [];
            for kk = 1:numel(fn)
                if startsWith(fn{kk}, run_prefix)
                    run_nums(end+1) = str2double(extractAfter(fn{kk}, run_prefix)); 
                end
            end
            if ~isempty(run_nums)
                last_run = max(run_nums);
            end
        end
        
        if last_run >= MC
            fprintf(' -> [sys=%d si=%d] Already present %d/%d MC runs. Skipping.\n', sysnum, si, last_run, MC);
            continue
        elseif last_run > 0
            fprintf(' -> [sys=%d si=%d] Found %d runs. Adding from %d to %d.\n', sysnum, si, last_run, last_run + 1, MC);
        else
            fprintf(' -> [sys=%d si=%d] Starting new runs from 1 to %d.\n', sysnum, si, MC);
        end
        
        runs_to_do  = (last_run + 1):MC;
        run_results = cell(MC, 1);
        
        parfor i = runs_to_do
            % ---------- DeePC ----------
            tic;
            [u_seq_opt, y_seq_opt, x_seq, y_pred_seq, is_infeasible] = run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R, x0, T_c, sysnum, flag_SPC, i);
            elapsed_time = toc;
            fprintf('[sys=%d si=%d] Run %d | execution time: %.6f s\n', sysnum, si, i, elapsed_time);
            
            if ~is_infeasible 
                error = desired_trajectory(1:T_c) - y_seq_opt; 
                rms_error = sqrt(mean(error.^2)); 
            else
                error = NaN; 
                rms_error = NaN; 
            end
            
            run_results{i} = struct( ...
                'rms_error',     rms_error, ...
                'error',         error(:), ...
                'y_seq_opt',     y_seq_opt(:), ...
                'u_seq_opt',     u_seq_opt(:), ...
                'elapsed_time',  elapsed_time, ...
                'sv_H_u_ext',    sv_H_u_ext(:), ...
                'sv_H_w',        sv_H_w(:), ...
                'is_infeasible', is_infeasible ...
            );
        end
        
        for i = runs_to_do
            varname = sprintf('s%02d_r%02d_case%02d_run%02d', sysnum, f_target, si, i);
            big.(varname) = run_results{i};
        end
        
        big.params = struct( ...
            'T_D',      T_D, ...
            'Ts',       Ts, ...
            'fs',       fs, ...
            'T_ini',    T_ini, ...
            'N',        N, ...
            'L',        L, ...
            'T',        T, ...
            'T_c',      T_c, ...
            'Q',        Q, ...
            'R',        R, ...
            'flag_SPC', flag_SPC, ...
            'f_target', f_target, ...
            'MC',       MC, ...
            'sysnum',   sysnum, ...
            'si',       si, ...
            'dim',      dim, ...
            'x0',       x0 ...
        );
        
        if ~isempty(runs_to_do)
            save(fname_big, '-struct', 'big');
            fprintf('Saved: %s\n', fname_big);
        end
    end
end