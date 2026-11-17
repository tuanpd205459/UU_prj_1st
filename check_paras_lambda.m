clc; clear; close all;
tic

%% 1. SYSTEM & SENSOR PARAMETERS (FIXED ACROSS SWEEP)
z1 = 350e-6;                    % Fixed distance z1 (m)
dz = 50e-6;                    % Fixed plane spacing z2 - z1 (m)
z2 = z1 + dz;

% Define sweep configurations: [lambda_0 (m), refractive_index]
configs = [
    405e-9, 1.000;              % Case 2: 405 nm in air
    405e-9, 1.515               % Case 3: 405 nm in immersion medium (n = 1.515)
];
short_labels = {'405nm (n=1)', '267nm (n=1.515)'};

area = 1000e-6;                  % FOV (m)
dx = 200e-9;                    % Pixel pitch (m)
N = round(area/dx);
dy = dx;

% Geometric NA depends only on FOV and z (independent of lambda)
NA_geo_z1_val = 1 / sqrt(1 + (2 * z1 / area)^2);
NA_geo_z2_val = 1 / sqrt(1 + (2 * z2 / area)^2);

% Define spatial coordinates
x = (-N/2 : N/2 - 1) * dx;
y = (-N/2 : N/2 - 1) * dy;
x_um = x * 1e6;
y_um = y * 1e6;
[~, mid_y] = min(abs(y));

% ROI settings for cross-section evaluation at y = 0
x_roi_um = 3.0;                 % Zoomed half-width (+/- 3 um) around center
idx_roi  = abs(x_um) <= x_roi_um;
x_cut    = x_um(idx_roi);

%% 2. CREATE OBJECT (COMPUTED ONCE OUTSIDE LOOP)
object_width_nm    = 200;       % Width of each slit (nm)
object_length_nm   = 5000;      % Length of each slit (nm)
center_distance_nm = 400;       % Center-to-center distance (nm) -> 200 nm gap

object_width    = object_width_nm * 1e-9;
object_length   = object_length_nm * 1e-9;
center_distance = center_distance_nm * 1e-9;

[X, Y] = meshgrid(x, y);
object = ones(N, N);
mask1 = abs(X + center_distance/2) <= object_width/2 & abs(Y) <= object_length/2;
mask2 = abs(X - center_distance/2) <= object_width/2 & abs(Y) <= object_length/2;
object(mask1 | mask2) = 0;

u_object = object;              % Incident amplitude A0 = 1
am = abs(u_object);
profile_true = am(mid_y, idx_roi);
U_obj_freq = FT2Dc(u_object);   % Pre-compute object FFT to save time

fprintf('================ FIXED SYSTEM & OBJECT PARAMETERS ================\n');
fprintf('Grid size (N x N)   : %d x %d pixels | FOV = %.0f um\n', N, N, area*1e6);
fprintf('Pixel pitch (dx)    : %.1f nm | z1 = %.0f um, z2 = %.0f um (dz = %.0f um)\n', dx*1e9, z1*1e6, z2*1e6, dz*1e6);
fprintf('NA_geo (z1 / z2)    : %.4f / %.4f\n', NA_geo_z1_val, NA_geo_z2_val);
fprintf('Object              : Two slits (Width = %.0f nm, Length = %.0f nm, Gap = %.0f nm)\n', ...
        object_width_nm, object_length_nm, center_distance_nm - object_width_nm);
fprintf('==================================================================\n');

% Create output directory for saving figures
output_dir = 'Lambda_Sweep_Results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end


%% 3. PRE-ALLOCATE STORAGE FOR SWEEP METRICS
num_cases = size(configs, 1);

lambda_eff_nm_arr = zeros(1, num_cases);
NA_Nyq_arr        = zeros(1, num_cases);
NA_geo_z1_arr     = zeros(1, num_cases);
NA_geo_z2_arr     = zeros(1, num_cases);
NA_eff_z1_arr     = zeros(1, num_cases);
NA_eff_z2_arr     = zeros(1, num_cases);
R_z1_nm_arr       = zeros(1, num_cases);
R_z2_nm_arr       = zeros(1, num_cases);
z_crit_um_arr     = zeros(1, num_cases);
z_max_um_arr      = zeros(1, num_cases);
rmse_arr          = zeros(1, num_cases);
regime_list       = strings(1, num_cases);

% Store cross-section profiles
x_cut_cell     = cell(1, num_cases);
prof_true_cell = cell(1, num_cases);
prof_rec_cell  = cell(1, num_cases);

Iterations = 100;

%% 4. PARAMETER SWEEP LOOP OVER WAVELENGTH & REFRACTIVE INDEX
fprintf('\nStarting sweep over %d configurations (Iterations per case: %d)...\n\n', num_cases, Iterations);

for k = 1:num_cases
    lambda_0          = configs(k, 1);
    diffractive_index = configs(k, 2);
    lambda            = lambda_0 / diffractive_index; % Effective wavelength in medium
    
    % Wavelength-dependent limits
    z_max  = N * dx^2 / lambda;
    NA_Nyq = lambda / (2 * dx);
    
    % Critical distance for R = dx
    if NA_Nyq < 1
        z_crit = (area / 2) * sqrt(1 - NA_Nyq^2) / NA_Nyq;
    else
        z_crit = 0;
    end
    
    % Calculate NA and Theoretical Resolution
    NA_geo_z1 = NA_geo_z1_val;
    NA_geo_z2 = NA_geo_z2_val;
    
    NA_eff_z1 = min([1, NA_geo_z1, NA_Nyq]);
    NA_eff_z2 = min([1, NA_geo_z2, NA_Nyq]);
    
    R_z1 = lambda / (2 * NA_eff_z1);
    R_z2 = lambda / (2 * NA_eff_z2);
    
    % Store metrics
    lambda_eff_nm_arr(k) = lambda * 1e9;
    NA_Nyq_arr(k)        = NA_Nyq;
    NA_geo_z1_arr(k)     = NA_geo_z1;
    NA_geo_z2_arr(k)     = NA_geo_z2;
    NA_eff_z1_arr(k)     = NA_eff_z1;
    NA_eff_z2_arr(k)     = NA_eff_z2;
    R_z1_nm_arr(k)       = R_z1 * 1e9;
    R_z2_nm_arr(k)       = R_z2 * 1e9;
    z_crit_um_arr(k)     = z_crit * 1e6;
    z_max_um_arr(k)      = z_max * 1e6;
    
    % Determine limiting regime
    if NA_Nyq < min(NA_geo_z1, NA_geo_z2)
        regime_list(k) = "Sampling-limited";
    elseif max(z1, z2) > z_crit
        regime_list(k) = "Aperture-limited";
    else
        regime_list(k) = "Diffraction-limited";
    end
    
    % Hologram Simulation
    prop_z1 = Propagator_BLASM(N, lambda, dx, z1);
    prop_z2 = Propagator_BLASM(N, lambda, dx, z2);
    
    amp1_measured = abs(IFT2Dc(U_obj_freq .* prop_z1));
    amp2_measured = abs(IFT2Dc(U_obj_freq .* prop_z2));
    
    % Plane-to-Plane Gerchberg-Saxton Phase Retrieval
    prop_dz       = Propagator_BLASM(N, lambda, dx, dz);
    prop_minus_dz = Propagator_BLASM(N, lambda, dx, -dz);
    prop_minus_z1 = Propagator_BLASM(N, lambda, dx, -z1);
    
    u1 = amp1_measured;
    for iter = 1:Iterations
        u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);
        u2     = amp2_measured .* exp(1i * angle(u2_sim));
        
        u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz);
        u1     = amp1_measured .* exp(1i * angle(u1_sim));
    end
    
    u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);
    amp_rec   = abs(u_obj_rec);
    
    % Extract cross-section at y = 0 and store into cells
    prof_rec = amp_rec(mid_y, idx_roi);
    x_cut_cell{k}     = x_cut;
    prof_true_cell{k} = profile_true;
    prof_rec_cell{k}  = prof_rec;
    
    % Calculate RMSE on ROI cross-section
    rmse_arr(k) = sqrt(mean((prof_rec - profile_true).^2));
    
    % Print log line-by-line
    fprintf('--- [%d/%d] %s ---\n', k, num_cases, short_labels{k});
    fprintf('  Lambda_eff : %.1f nm (NA_Nyq = %.4f)\n', lambda*1e9, NA_Nyq);
    fprintf('  z_crit     : %.1f um | z_max = %.2f mm\n', z_crit*1e6, z_max*1e3);
    fprintf('  At z1      : NA_eff = %.4f | R1 = %.1f nm\n', NA_eff_z1, R_z1*1e9);
    fprintf('  At z2      : NA_eff = %.4f | R2 = %.1f nm\n', NA_eff_z2, R_z2*1e9);
    fprintf('  RMSE       : %.4f\n', rmse_arr(k));
    fprintf('  Regime     : %s\n\n', regime_list(k));
    
    % Save individual evaluation figure
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 450]);
    
    subplot(1, 2, 1);
    imagesc(x_um, y_um, amp_rec);
    axis image;
    set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 14);
    colormap gray; colorbar; hold on;
    plot([min(x_cut), max(x_cut)], [0, 0], 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlim([-15 15]); ylim([-15 15]);
    title(sprintf('Reconstructed Amplitude (%s)', short_labels{k}));
    xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel('$y\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    
    subplot(1, 2, 2);
    stairs(x_cut, profile_true, 'Color', '#0072BD', 'LineWidth', 1.8); hold on;
    plot(x_cut, prof_rec, 'Color', '#D95319', 'LineWidth', 1.8, 'LineStyle', '--');
    grid on; grid minor;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
    xlim([-x_roi_um, x_roi_um]);
    ylim([0, 1.2]);
    xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel('Amplitude (a.u.)', 'Interpreter', 'latex');
    legend('Ground Truth', sprintf('Rec (R1=%.1fnm, RMSE=%.3f)', R_z1*1e9, rmse_arr(k)), 'Location', 'southeast');
    title(sprintf('Cross-section y = 0 (%s | %s)', short_labels{k}, regime_list(k)));
    
    fname = sprintf('Rec_Case%d_%dnm_n%.3f.png', k, round(lambda_0*1e9), diffractive_index);
    saveas(fig, fullfile(output_dir, fname));
    close(fig);
end

elapsed_time = toc;
fprintf('Completed in %.2f seconds.\n', elapsed_time);

%% 5. SUMMARY PLOTS ACROSS ALL CONFIGURATIONS
figure('Color', 'w', 'Name', 'Sweep Summary vs Wavelength', 'Position', [100 100 1350 480]);

% Panel 1: Overlay of cross-sections at y = 0
subplot(1, 2, 1);
stairs(x_cut_cell{1}, prof_true_cell{1}, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground Truth'); hold on;
colors = lines(num_cases);
for k = 1:num_cases
    plot(x_cut_cell{k}, prof_rec_cell{k}, 'LineWidth', 1.8, 'Color', colors(k, :), ...
         'DisplayName', sprintf('%s (R1=%.0f nm)', short_labels{k}, R_z1_nm_arr(k)));
end
grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlim([-2.5, 2.5]);
ylim([0, 1.2]);
xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel('Amplitude (a.u.)', 'Interpreter', 'latex');
title('Cross-Section Amplitude at y = 0');
legend('Location', 'southeast', 'FontSize', 12);

% Panel 2: Theoretical Resolution R_z1 & R_z2
subplot(1, 2, 2);
x_idx = 1:num_cases;
plot(x_idx, R_z1_nm_arr, '-o', 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', 'R at z_1'); hold on;
plot(x_idx, R_z2_nm_arr, '-s', 'LineWidth', 2, 'MarkerSize', 7, 'DisplayName', 'R at z_2');
yline(200, 'r--', 'Object Gap (200 nm)', 'LineWidth', 1.2, 'DisplayName', 'Gap (200 nm)');
ylabel('Theoretical Resolution $R$ (nm)', 'Interpreter', 'latex');
ylim([0, max(R_z2_nm_arr)*1.2]);
grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14, ...
         'XTick', x_idx, 'XTickLabel', short_labels);
xlim([0.7, num_cases + 0.3]);
xlabel('Configuration', 'Interpreter', 'latex');
title('Theoretical Resolution (z_1 & z_2)');
legend('Location', 'northeast', 'FontSize', 12);


%% =========================================================================
%                                FUNCTIONS
% =========================================================================
function p = Propagator_BLASM(N, lambda, dx, z)
    df = 1 / (N * dx);
    f = (-N/2 : N/2 - 1) * df;
    [FX, FY] = meshgrid(f, f);

    alpha = lambda * FX;
    beta  = lambda * FY;

    limit = 1 / sqrt((2 * abs(z) / (N * dx))^2 + 1);
    mask = (alpha.^2 + beta.^2) <= limit^2;

    p = zeros(N, N);
    p(mask) = exp(1i * 2 * pi * (z / lambda) * ...
              sqrt(1 - alpha(mask).^2 - beta(mask).^2));
end

function out = FT2Dc(in)
    out = fftshift(fft2(ifftshift(in)));
end

function out = IFT2Dc(in)
    out = fftshift(ifft2(ifftshift(in)));
end