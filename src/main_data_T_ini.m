clc
clear
close all

%% User Parameters (Grouped for clean saving)
user_params = struct();
user_params.Ts    = 0.01;                  % Sampling time
user_params.fs    = 1/user_params.Ts;      % Sampling frequency
user_params.N     = 4;                     % Future horizon

% Predefined ranges
user_params.T_ini_vec = [4, 6, 8, 10];     % [4, 6, 8, 10, 12] 

% --- ORIGINAL BAND PARAMETERS ---
user_params.band_width_pct = 0.3;          % Band width (30% of Nyquist)
user_params.band_step_pct  = 0.1;          % Band shift (10% of Nyquist)

% CHOOSE THE BAND TO TEST: 1, 2, 3... 
% (White Noise is now automatically run as a benchmark)
user_params.chosen_band_idx = 3; 
user_params.T_c = 100;

% DeePC Weights
user_params.Q = 100;
user_params.R = 1;   
user_params.flag_SPC = 0;                  % 1: use SPC instead

% Target frequency parameters
user_params.w_target_vec_pi = 0.04;
user_params.f_target_vec = user_params.w_target_vec_pi * (user_params.fs / 2);
user_params.MC = 2;                       % Number of MC simulations

% Variable extraction for code readability
fs = user_params.fs;
N = user_params.N;
T_c = user_params.T_c;
Q = user_params.Q;
R = user_params.R;
flag_SPC = user_params.flag_SPC;
MC = user_params.MC;
chosen_band = user_params.chosen_band_idx;

data_dir = fullfile(pwd, 'data_pred_4/Rho/paper');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

% --- START PARALLEL POOL ---
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end

%% Simulation Loops 
for sysnum = 1:1
    if ismember(sysnum, [1, 2, 3])
        dim = 2;
    else
        dim = 4;
    end
    x0 = zeros(dim,1);

    % =========================================================================
    % PAST HORIZON LOOP (T_ini)
    % =========================================================================
    for idx_ini = 1:length(user_params.T_ini_vec)
        T_ini = user_params.T_ini_vec(idx_ini);
        L = T_ini + N;
        
        % T_D vector dynamically calculated based on current T_ini
        T_D_vec = [2*(L+2), 4*(L+2), 6*(L+2), 8*(L+2), 10*(L+2)]; % Systems 1, 2, 3
        % T_D_vec = [2*(L+4), 2*(2*(L+4)), 4*(2*(L+4)), 6*(2*(L+4)), 8*(2*(L+4)), 10*(2*(L+4))]; % Systems 4, 5
        
        fprintf('\n******************************************************\n');
        fprintf(' STARTING TEST WITH T_ini = %d (L = %d)\n', T_ini, L);
        fprintf('******************************************************\n');

        % --- EXPERIMENT SETUP STRUCTURING ---
        test_setups = struct('type', {}, 'T_D', {});
        
        % 1. White Noise Setup (executed ONCE per T_ini with max T_D for robustness)
        test_setups(1).type = 'WN';
        test_setups(1).T_D  = max(T_D_vec);
        
        % 2. Multisine Setup for all requested T_D values
        for td_idx = 1:length(T_D_vec)
            test_setups(end+1).type = 'Multisine';
            test_setups(end).T_D = T_D_vec(td_idx);
        end

        % =========================================================================
        % SETUP LOOP (White Noise + Multisine across various T_D)
        % =========================================================================
        for test_idx = 1:length(test_setups)
            
            current_type = test_setups(test_idx).type;
            T_D = test_setups(test_idx).T_D;
            T   = 10 * max(T_D_vec); 
            
            k = 0:(T-1);
            
            rng(1)
            % Signal Generation
            if strcmp(current_type, 'WN')
                fprintf('\n======================================================\n');
                fprintf('[sys=%d | T_ini=%d] WHITE NOISE Test (Baseline Comparison) | T_D_used = %d\n', sysnum, T_ini, T_D);
                rng(1)
                u_seq    = randn(1, T); % Persistently exciting broad-spectrum signal
                j_min    = 0;
                j_max    = 0;
                w_center = NaN;         % Uniform spectrum across all frequencies
                band_label = 'WN';
                
                fprintf('Signal: White Noise\n');
                fprintf('======================================================\n');
                
            else
                % --- HARMONIC CALCULATION ACCORDING TO PERCENTAGE LOGIC ---
                j_max_tot = floor(T_D/2);                                
                Jnum      = max(1, round(user_params.band_width_pct * j_max_tot));   
                Jstep     = max(1, round(user_params.band_step_pct  * j_max_tot));   
                j_min_vec = 1 : Jstep : (j_max_tot - Jnum + 1);          
                n_bands   = numel(j_min_vec); 
                
                if chosen_band > n_bands
                    fprintf('\n [!] WARNING: With T_D=%d band %d does not exist (max bands: %d). Skipping.\n', T_D, chosen_band, n_bands);
                    continue; 
                end
                
                fprintf('\n======================================================\n');
                fprintf('[sys=%d | T_ini=%d] MULTISINE Test | T_D = %d | Selected Band: %d\n', sysnum, T_ini, T_D, chosen_band);
                
                j_min = j_min_vec(chosen_band);
                j_max = j_min + Jnum - 1;
                Nper = T_D; 
                f0   = fs / Nper;
                Jn    = j_max - j_min + 1;
                u_seq = zeros(size(k));
                idx   = 0;
                
                for j = j_min:j_max
                    idx   = idx + 1;
                    phi_j = -pi * idx*(idx-1) / Jn;   
                    A_j   = 1;
                    u_seq = u_seq + A_j * sin( 2*pi*j*(f0/fs)*k + phi_j ); 
                end
                w_center = 2*pi*((j_min+j_max)/2) / Nper;
                
                f_min_real = j_min * f0;
                f_max_real = j_max * f0;
                band_label = sprintf('bandIdx%02d', chosen_band);
                
                fprintf('Harmonics: [%d - %d] -> Real Freq: [%.2f, %.2f] Hz\n', j_min, j_max, f_min_real, f_max_real);
                fprintf('======================================================\n');
            end

            % Dynamics and Hankel
            [y_seq, ~, n] = system_dynamics_no_noise(u_seq, x0, sysnum);
            
            u_seq_cut = u_seq(end - T_D + 1 : end);
            y_seq_cut = y_seq(end - T_D + 1 : end);
            [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_seq_cut, y_seq_cut, T_ini, N, dim);
            
            H_u_ext = hankel(u_seq_cut(1:L+n), u_seq_cut(L+n:length(u_seq_cut)));
            sv_H_u_ext = svd(H_u_ext);
            
            H_w = [U_p; U_f; Y_p; Y_f];
            sv_H_w = svd(H_w);
            
            is_rank_deficient = (rank(H_w) < dim + L);

            % =================================================================
            % TARGET FREQUENCY LOOP
            % =================================================================
            for idx_f = 1:length(user_params.f_target_vec)
                
                f_target = user_params.f_target_vec(idx_f);
                ftgt_str = strrep(sprintf('%.2f', f_target), '.', 'p');
                
                % Dynamic filename creation: includes test type (WN or selected band)
                fname_big = sprintf('sys%02d_Tini%02d_TD%04d_ftgt%sHz_%s.mat', ...
                                    sysnum, T_ini, T_D, ftgt_str, band_label);
                fname_big = fullfile(data_dir, fname_big);
                
                % Rank Deficiency Handling
                if is_rank_deficient
                    fprintf('  [sys=%d T_ini=%d TD=%d] Rank deficiency! Skipping Monte Carlo.\n', sysnum, T_ini, T_D);
                    big_rd = struct();
                    if isfile(fname_big)
                        big_rd = load(fname_big);
                    end
                    big_rd.is_rank_deficient = true;
                    big_rd.sv_H_w = sv_H_w(:);
                    big_rd.T_D_used = T_D;
                    big_rd.T_ini_used = T_ini;
                    big_rd.user_params = user_params;
                    save(fname_big, '-struct', 'big_rd');
                    continue;
                end
                
                desired_trajectory = sin( f_target * 2*pi / fs * (0:(T-1)) );

                % ---- Resume / File Saving Handling ----
                run_prefix = sprintf('run_'); 
                last_run   = 0;
                big        = struct();
                
                big.user_params = user_params;
                if isfile(fname_big)
                    big = load(fname_big);
                    fn  = fieldnames(big);
                    run_nums = [];
                    for kk = 1:numel(fn)
                        if startsWith(fn{kk}, run_prefix)
                            run_nums(end+1) = str2double(extractAfter(fn{kk}, run_prefix)); %#ok<AGROW>
                        end
                    end
                    if ~isempty(run_nums)
                        last_run = max(run_nums);
                    end
                end
                if last_run >= MC
                    fprintf('  [sys=%d T_ini=%d TD=%d] Completed (%d/%d runs). Skipping.\n', sysnum, T_ini, T_D, last_run, MC);
                    continue
                elseif last_run > 0
                    fprintf('  [sys=%d T_ini=%d TD=%d] Resuming MC from %d to %d\n', sysnum, T_ini, T_D, last_run + 1, MC);
                else
                    fprintf('  [sys=%d T_ini=%d TD=%d] Starting Monte Carlo (run 1 to %d)...\n', sysnum, T_ini, T_D, MC);
                end

                % =============================================================
                % MONTE CARLO LOOP (PARALLELIZED)
                % =============================================================
                runs_to_do = (last_run+1):MC;
                run_results = cell(MC, 1);
                parfor i = runs_to_do
                    tic;
                    [u_seq_opt, y_seq_opt, x_seq, y_pred_seq, is_infeasible, y_ini_0, u_ini_0] = ...
                        run_DeePC(U_p, Y_p, U_f, Y_f, desired_trajectory, T_ini, N, Q, R, x0, T_c, sysnum, flag_SPC, i);
                    elapsed_time = toc;
                    
                    if ~is_infeasible
                        error = desired_trajectory(1:T_c) - y_seq_opt;
                        rms_error = sqrt(mean(error.^2));
                    else
                        error = NaN;
                        rms_error = NaN;
                    end
                    run_results{i} = struct( ...
                        'rms_error',         rms_error, ...
                        'error',             error(:), ...
                        'y_seq_opt',         y_seq_opt(:), ...
                        'u_seq_opt',         u_seq_opt(:), ...
                        'elapsed_time',      elapsed_time, ...
                        'sv_H_u_ext',        sv_H_u_ext(:), ...
                        'sv_H_w',            sv_H_w(:), ...
                        'is_infeasible',     is_infeasible, ...
                        'j_min',             j_min, ...
                        'j_max',             j_max, ...
                        'w_center',          w_center, ...
                        'f_target',          f_target, ...
                        'T_D_used',          T_D, ... 
                        'T_ini_used',        T_ini, ...
                        'y_ini_0',           y_ini_0(:), ...
                        'u_ini_0',           u_ini_0(:), ...
                        'is_rank_deficient', is_rank_deficient ...
                    );
                end 

                % ---- Dump Results and Save ----
                for i = runs_to_do
                    varname = sprintf('%s%02d', run_prefix, i);
                    big.(varname) = run_results{i};
                end
                
                if ~isempty(runs_to_do)
                    big.is_rank_deficient = false;
                    save(fname_big, '-struct', 'big');
                    fprintf('  [sys=%d T_ini=%d TD=%d] Saved: %s\n', sysnum, T_ini, T_D, fname_big);
                end
            end % End f_target loop
        end % End setup test loop
    end % End T_ini loop
end % End sysnum loop