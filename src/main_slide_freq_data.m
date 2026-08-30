clc
clear
close all

%% Parameters
Ts    = 0.01;    % Sampling time
fs    = 1/Ts;    % Sampling frequency
% Hankel parameters
T_ini = 8;       % Past horizon
N     = 4;       % Future horizon
L     = T_ini + N;
% T_D_vec = [2*(L+2), 2*(2*(L+2)), 4*(2*(L+2)), 6*(2*(L+2)), 8*(2*(L+2)), 10*(2*(L+2))]; % Systems 1,2,3
T_D_vec = [2*(L+4), 2*(2*(L+4)), 4*(2*(L+4)), 6*(2*(L+4)), 8*(2*(L+4)), 10*(2*(L+4))]; % Systems 4,5
T_D_max = max(T_D_vec);              % Maximum size to guarantee coherent excitation
T   = 10*T_D_max; % Extended data collection to remove transient (calculated on max)
T_c = 100;
% DeePC Weights
Q = 100;
R = 1;   
flag_SPC = 0; % 1: use SPC instead
% --- TARGET FREQUENCY PARAMETERS ---
w_target_vec_pi = 0.04;
% Conversion to physical frequencies [Hz] 
f_target_vec = w_target_vec_pi * (fs / 2);

%% Sliding-band multisine excitation parameters
band_width_pct = 0.3;    % Band width
band_step_pct  = 0.1;    % Band shift
MC = 2;                 % Number of MC simulations
data_dir = fullfile(pwd, 'data_pred_4\DeePC_test_numb_data\Diff_sgn\Sys_45_sliding_band_30_10');
if ~exist(data_dir, 'dir')
    mkdir(data_dir);
end

%% 1. DETERMINATION OF MAXIMUM NUMBER OF BANDS
max_n_bands = 0;
for TD = T_D_vec
    j_max_tot = floor(TD / 2);
    J_num  = max(1, round(band_width_pct * j_max_tot));
    J_step = max(1, round(band_step_pct  * j_max_tot));
    n_b = length(1 : J_step : (j_max_tot - J_num + 1));
    if n_b > max_n_bands
        max_n_bands = n_b;
    end
end

%% 2. DYNAMIC CALCULATION OF FREQUENCY INTERVALS [f_min, f_max]
table_data = repmat({"—"}, max_n_bands, length(T_D_vec));
col_names  = cell(1, length(T_D_vec));
for i_TD = 1:length(T_D_vec)
    TD = T_D_vec(i_TD);
    col_names{i_TD} = sprintf('TD_%d', TD);
    
    f0 = fs / TD;                         % Spectral resolution
    j_max_tot = floor(TD / 2);
    J_num  = max(1, round(band_width_pct * j_max_tot));
    J_step = max(1, round(band_step_pct  * j_max_tot));
    
    j_start_vec = 1 : J_step : (j_max_tot - J_num + 1);
    
    for b = 1:length(j_start_vec)
        j_start = j_start_vec(b);
        f_min = j_start * f0;
        f_max = (j_start + J_num - 1) * f0;
        
        % Cell formatting with endpoints in Hz
        table_data{b, i_TD} = sprintf('[%.2f, %.2f]', f_min, f_max);
    end
end

%% 3. CREATION AND DISPLAY OF THE MATLAB TABLE
Bands = strcat("Band ", string(pad(string((1:max_n_bands)'), 2, 'left', '0')));
T_summary = cell2table(table_data, 'VariableNames', col_names);
T_summary = [table(Bands) T_summary];

% Display on Command Window
disp('========================================================================================');
disp('                  FREQUENCY INTERVALS TABLE [f_min, f_max] in Hz                        ');
disp('========================================================================================');
disp(T_summary);

%% 4. AUTOMATIC EXPORT (EXCEL / PDF)
% Export directly to Excel
writetable(T_summary, 'Calculated_Bands_Table.xlsx');
disp('Table successfully saved to "Calculated_Bands_Table.xlsx"');

% --- START PARALLEL POOL ---
poolobj = gcp('nocreate');
if isempty(poolobj)
    parpool('local');
end

%% Simulation Loops 
for sysnum = 5:5
    if ismember(sysnum, [1, 2, 3])
        dim = 2;
    else
        dim = 4;
    end
    x0 = zeros(dim,1);

    % =========================================================================
    % LOOP OVER DATA DIMENSION FOR THE SPECIFIC EXPERIMENT (T_D)
    % =========================================================================
    for td_idx = 1:length(T_D_vec)
        
        T_D = T_D_vec(td_idx);
        T   = 10 * T_D; % Extension to extinguish transient for current T_D
        
        % Band parameters calculated on current T_D (ensures periodicity)
        j_max_tot = floor(T_D/2);                                % Nyquist for T_D
        Jnum      = max(1, round(band_width_pct * j_max_tot));   % Number of harmonics
        Jstep     = max(1, round(band_step_pct  * j_max_tot));   % Sliding step
        j_min_vec = 1 : Jstep : (j_max_tot - Jnum + 1);          
        n_bands   = numel(j_min_vec);                            
        
        % Include band 0 (White Noise) + multisine bands (1:n_bands)
        bands_to_run = [0, 1:n_bands];
        
        fprintf('\n======================================================\n');
        fprintf('[sys=%d] Test with T_D = %d (%d multisine bands + White Noise band00)\n', sysnum, T_D, n_bands);
        fprintf('======================================================\n');

        for band_idx = 1:n_bands
            k = 0:(T-1);
            if band_idx == 0
                % =============================================================
                % 1A. WHITE NOISE GENERATION (BAND 00)
                % =============================================================
                u_seq    = randn(1, T); % Persistently exciting broad-spectrum signal
                j_min    = 0;
                j_max    = 0;
                w_center = NaN;         % Uniform spectrum across all frequencies
            else
                % =============================================================
                % 1B. PERIODIC MULTISINE GENERATION DEDICATED TO T_D
                % =============================================================
                j_min = j_min_vec(band_idx);
                j_max = j_min + Jnum - 1;
                Nper = T_D;             % Fundamental period equal exactly to T_D
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
            end

            % Response simulation to reach steady state
            [y_seq, ~, n] = system_dynamics_no_noise(u_seq, x0, sysnum);
            
            % Extract last T_D samples
            u_seq_cut = u_seq(end - T_D + 1 : end);
            y_seq_cut = y_seq(end - T_D + 1 : end);

            % Hankel Matrix Construction
            [U_p, U_f, Y_p, Y_f] = construct_Hankel(u_seq_cut, y_seq_cut, T_ini, N, dim);
            
            H_u_ext = hankel(u_seq_cut(1:L+n), u_seq_cut(L+n:length(u_seq_cut)));
            sv_H_u_ext = svd(H_u_ext);
            
            H_w = [U_p; U_f; Y_p; Y_f];
            sv_H_w = svd(H_w);
            
            % --- RIGOROUS RANK DEFICIENCY CHECK ---
            is_rank_deficient = (rank(H_w) < dim + L);

            % =================================================================
            % 2. TARGET FREQUENCY LOOP
            % =================================================================
            for idx_f = 1:length(f_target_vec)
                
                f_target = f_target_vec(idx_f);
                w_target_pi = w_target_vec_pi(idx_f);
                
                w_str = strrep(sprintf('%g', w_target_pi), '.', 'p');
                
                % Saving with %02d: band_idx = 0 produces "band00.mat"
                fname_big = sprintf('s%02d_TD%04d_w%s_band%02d.mat', sysnum, T_D, w_str, band_idx);
                fname_big = fullfile(data_dir, fname_big);
                
                % If rank deficiency occurs, save outcome and skip Monte Carlo optimization
                if is_rank_deficient
                    fprintf('  [sys=%d TD=%d band=%02d] Rank deficiency detected! Skipping Monte Carlo.\n', sysnum, T_D, band_idx);
                    big_rd = struct();
                    if isfile(fname_big)
                        big_rd = load(fname_big);
                    end
                    big_rd.is_rank_deficient = true;
                    big_rd.sv_H_w = sv_H_w(:);
                    big_rd.T_D_used = T_D;
                    save(fname_big, '-struct', 'big_rd');
                    continue;
                end
                
                desired_trajectory = sin( f_target * 2*pi / fs * (0:(T-1)) );

                % ---- Resume / File Saving Handling ----
                run_prefix = sprintf('s%02d_TD%04d_w%s_band%02d_run', sysnum, T_D, w_str, band_idx);
                last_run   = 0;
                big        = struct();
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

                % Print Monte Carlo simulation status (Resuming / Completed / Starting)
                if last_run >= MC
                    fprintf('  [sys=%d TD=%d band=%02d] Simulations already completed (%d/%d runs). Skipping.\n', ...
                        sysnum, T_D, band_idx, last_run, MC);
                    continue
                elseif last_run > 0
                    fprintf('  [sys=%d TD=%d band=%02d] Resuming Monte Carlo simulation from run %d to %d (completed %d/%d)\n', ...
                        sysnum, T_D, band_idx, last_run + 1, MC, last_run, MC);
                else
                    fprintf('  [sys=%d TD=%d band=%02d] Starting Monte Carlo simulation (run 1 to %d)...\n', ...
                        sysnum, T_D, band_idx, MC);
                end

                % =============================================================
                % 3. MONTE CARLO LOOP (PARALLELIZED)
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
                        'w_target_pi',       w_target_pi, ...
                        'T_D_used',          T_D, ... 
                        'y_ini_0',           y_ini_0(:), ...
                        'u_ini_0',           u_ini_0(:), ...
                        'is_rank_deficient', is_rank_deficient ...
                    );
                end 

                % ---- Dump Results and Save ----
                for i = runs_to_do
                    varname = sprintf('s%02d_TD%04d_w%s_band%02d_run%02d', sysnum, T_D, w_str, band_idx, i);
                    big.(varname) = run_results{i};
                end
                
                if ~isempty(runs_to_do)
                    big.is_rank_deficient = false;
                    save(fname_big, '-struct', 'big');
                    fprintf('  [sys=%d TD=%d band=%02d] Saved: %s\n', sysnum, T_D, band_idx, fname_big);
                end
            end % End f_target loop
        end % End band loop
    end % End T_D loop
end % End sysnum loop