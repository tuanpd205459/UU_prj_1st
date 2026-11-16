clc; clear; close all;
tic

%% 1. FIXED SYSTEM & OBJECT PARAMETERS
lambda_0 = 405e-9;          % Wavelength (m)
diffractive_index = 1;
lambda = lambda_0 / diffractive_index;

z1 = 100e-6;               % Plane 1 (m)
z2 = 150e-6;               % Plane 2 (m)
dz = z2 - z1;

dx = 100e-9;                % Pixel pitch (m)
dy = dx;
Iterations = 100;           % GS Iterations

% Object parameters (Two dark slits)
object_width_nm    = 200;   % Width of each slit (nm)
object_length_nm   = 5000;  % Length of each slit (nm)
center_distance_nm = 200 + 200;

object_width    = object_width_nm    * 1e-9;
object_length   = object_length_nm   * 1e-9;
center_distance = center_distance_nm * 1e-9;

NA_Nyq = lambda / (2 * dx);

%% 2. SWEEP PARAMETERS & STORAGE INITIALIZATION
FOV_list_um = [100, 200, 300, 500];
num_fov = length(FOV_list_um);

NA_geo_z1_arr  = zeros(1, num_fov);
NA_geo_z2_arr  = zeros(1, num_fov);
NA_eff_z1_arr  = zeros(1, num_fov);
NA_eff_z2_arr  = zeros(1, num_fov);
R_z1_nm_arr    = zeros(1, num_fov);
R_z2_nm_arr    = zeros(1, num_fov);
z_crit_um_arr  = zeros(1, num_fov);
visibility_arr = zeros(1, num_fov);
rmse_arr       = zeros(1, num_fov);
regime_list    = strings(1, num_fov);

% Tạo thư mục lưu ảnh kết quả
output_dir = 'FOV_Sweep_Results';
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Lưu trữ các profile mặt cắt (dùng cell vì N thay đổi theo FOV)
x_roi_um = 3.0; % Vùng cắt ngang +/- 3 um quanh tâm
x_cut_cell     = cell(1, num_fov);
prof_true_cell = cell(1, num_fov);
prof_rec_cell  = cell(1, num_fov);

fprintf('================ FIXED SYSTEM & OBJECT PARAMETERS ================\n');
fprintf('Wavelength (lambda) : %.1f nm (n = %.2f) | NA_Nyq = %.5f\n', lambda*1e9, diffractive_index, NA_Nyq);
fprintf('Pixel pitch (dx)    : %.1f nm | z1 = %.0f um, z2 = %.0f um (dz = %.0f um)\n', dx*1e9, z1*1e6, z2*1e6, dz*1e6);
fprintf('Object              : Two slits (Width = %.0f nm, Gap = %.0f nm)\n', object_width_nm, center_distance_nm - object_width_nm);
fprintf('==================================================================\n');

%% 3. PARAMETER SWEEP LOOP OVER FOV
fprintf('\nStarting sweep over %d FOVs (Iterations per FOV: %d)...\n', num_fov, Iterations);

for k = 1:num_fov
    area = FOV_list_um(k) * 1e-6;
    N = round(area / dx);
    
    % Critical distance for R = dx
    if NA_Nyq < 1
        z_crit = (area / 2) * sqrt(1 - NA_Nyq^2) / NA_Nyq;
    else
        z_crit = 0;
    end
    
    % Calculate NA and Theoretical Resolution
    NA_geo_z1 = 1 / sqrt(1 + (2 * z1 / area)^2);
    NA_geo_z2 = 1 / sqrt(1 + (2 * z2 / area)^2);
    
    NA_eff_z1 = min([diffractive_index, NA_geo_z1, NA_Nyq]);
    NA_eff_z2 = min([diffractive_index, NA_geo_z2, NA_Nyq]);
    
    R_z1 = lambda / (2 * NA_eff_z1);
    R_z2 = lambda / (2 * NA_eff_z2);
    
    % Store metrics
    NA_geo_z1_arr(k) = NA_geo_z1;
    NA_geo_z2_arr(k) = NA_geo_z2;
    NA_eff_z1_arr(k) = NA_eff_z1;
    NA_eff_z2_arr(k) = NA_eff_z2;
    R_z1_nm_arr(k)   = R_z1 * 1e9;
    R_z2_nm_arr(k)   = R_z2 * 1e9;
    z_crit_um_arr(k) = z_crit * 1e6;
    
    % Determine limiting regime
    if NA_Nyq < min(NA_geo_z1, NA_geo_z2)
        regime_list(k) = "Sampling-limited";
    elseif max(z1, z2) > z_crit
        regime_list(k) = "Aperture-limited";
    else
        regime_list(k) = "Diffraction-limited";
    end
    
    fprintf('[%d/%d] FOV = %3.0f um (N=%4d) | z_crit = %6.1f um | NA_z1 = %.4f (R1 = %6.1f nm) | NA_z2 = %.4f (R2 = %6.1f nm) | %s\n', ...
            k, num_fov, FOV_list_um(k), N, z_crit*1e6, NA_eff_z1, R_z1*1e9, NA_eff_z2, R_z2*1e9, regime_list(k));
        
    % Build Coordinates & Object for current FOV
    x = (-N/2 : N/2 - 1) * dx;
    y = (-N/2 : N/2 - 1) * dy;
    x_um = x * 1e6;
    y_um = y * 1e6;
    [X, Y] = meshgrid(x, y);
    
    x1_obj = -center_distance / 2;
    x2_obj =  center_distance / 2;
    
    u_object = ones(N, N);
    mask1 = abs(X - x1_obj) <= object_width/2 & abs(Y) <= object_length/2;
    mask2 = abs(X - x2_obj) <= object_width/2 & abs(Y) <= object_length/2;
    u_object(mask1 | mask2) = 0;
    
    am = abs(u_object);
    [~, mid_y] = min(abs(y));
    
    % Hologram Simulation
    U_obj_freq = FT2Dc(u_object);
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
    
    % Extract cross-section at y = 0
    idx_roi      = abs(x_um) <= x_roi_um;
    x_cut        = x_um(idx_roi);
    profile_true = am(mid_y, idx_roi);
    prof_rec     = amp_rec(mid_y, idx_roi);
    
    x_cut_cell{k}     = x_cut;
    prof_true_cell{k} = profile_true;
    prof_rec_cell{k}  = prof_rec;
    
    % Calculate contrast/visibility between the two dark slits & RMSE
    [~, idx_gap]   = min(abs(x_cut - 0));
    [~, idx_slit1] = min(abs(x_cut - x1_obj*1e6));
    [~, idx_slit2] = min(abs(x_cut - x2_obj*1e6));
    
    I_peak_gap = prof_rec(idx_gap);
    I_trough   = mean([prof_rec(idx_slit1), prof_rec(idx_slit2)]);
    visibility_arr(k) = max(0, (I_peak_gap - I_trough) / (I_peak_gap + I_trough));
    rmse_arr(k)       = sqrt(mean((profile_true - prof_rec).^2));
    
    % Save individual evaluation figure (invisible mode for fast execution)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 450]);
    
    subplot(1, 2, 1);
    imagesc(x_um, y_um, amp_rec);
    axis image;
    set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 11);
    colormap gray; colorbar; hold on;
    plot([min(x_cut), max(x_cut)], [0, 0], 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlim([-15 15]); ylim([-15 15]);
    title(sprintf('Reconstructed Amplitude (FOV = %.0f \\mum)', FOV_list_um(k)), 'FontName', 'Times New Roman');
    xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel('$y\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    
    subplot(1, 2, 2);
    stairs(x_cut, profile_true, 'Color', '#0072BD', 'LineWidth', 1.8); hold on;
    plot(x_cut, prof_rec, 'Color', '#D95319', 'LineWidth', 1.8, 'LineStyle', '--');
    grid on; grid minor;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
    xlim([-x_roi_um, x_roi_um]);
    ylim([0, 1.25]);
    xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel('Amplitude (a.u.)', 'Interpreter', 'latex');
    legend('Ground Truth', sprintf('Rec (R_{z1} = %.1f nm)', R_z1*1e9), 'Location', 'southeast');
    title(sprintf('Cross-section at y = 0 (NA = %.3f, V = %.2f)', NA_eff_z1, visibility_arr(k)), 'FontName', 'Times New Roman');
    
    % Export individual figure at 600 DPI
    exportgraphics(fig, fullfile(output_dir, sprintf('Eval_FOV_%dum.png', FOV_list_um(k))), 'Resolution', 600);
    close(fig);
    
    % Free large matrices before next FOV
    clear X Y u_object am U_obj_freq prop_z1 prop_z2 prop_dz prop_minus_dz prop_minus_z1 u1 u2 u1_sim u2_sim u_obj_rec amp_rec;
end

elapsed_time = toc;
fprintf('Sweep completed in %.2f seconds! Individual figures saved in "%s/".\n', elapsed_time, output_dir);

%% 4. COMBINED OVERLAY & VISIBILITY SUMMARY FIGURE
fig_summary = figure('Name', 'FOV Sweep Summary', 'Color', 'w', 'Position', [100 100 1200 480]);

% Panel 1: Overlay all cross-sections on a single plot
subplot(1, 2, 1);
stairs(x_cut_cell{1}, prof_true_cell{1}, 'Color', '#0072BD', 'LineWidth', 2.0); hold on;
colors = {'#D95319', '#EDB120', '#7E2F8E', '#77AC30'};
styles = {'--', '-.', ':', '-'};
leg_str = cell(1, num_fov + 1);
leg_str{1} = 'Ground Truth';

for k = 1:num_fov
    plot(x_cut_cell{k}, prof_rec_cell{k}, 'Color', colors{k}, 'LineStyle', styles{k}, 'LineWidth', 1.8);
    leg_str{k+1} = sprintf('FOV = %d \\mum (R_1=%.0f nm, V=%.2f)', FOV_list_um(k), R_z1_nm_arr(k), visibility_arr(k));
end
grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlim([-x_roi_um, x_roi_um]);
ylim([0, 1.25]);
xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel('Amplitude (a.u.)', 'Interpreter', 'latex');
title('Overlay Cross-sections at $y = 0$', 'Interpreter', 'latex');
legend(leg_str, 'Location', 'southeast', 'FontSize', 9.5);

% Panel 2: Contrast/Visibility & Theoretical Resolution vs FOV
subplot(1, 2, 2);
yyaxis left
plot(FOV_list_um, visibility_arr, '-o', 'Color', '#0072BD', 'LineWidth', 2.0, 'MarkerSize', 7, 'MarkerFaceColor', '#0072BD');
ylabel('Slit Contrast / Visibility ($V$)', 'Interpreter', 'latex');
ylim([0, 1.05]);

yyaxis right
plot(FOV_list_um, R_z1_nm_arr, '--s', 'Color', '#D95319', 'LineWidth', 1.8, 'MarkerSize', 7);
ylabel('Theoretical Resolution $R_{z1}$ (nm)', 'Interpreter', 'latex');

grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('Field of View, FOV ($\mu\mathrm{m}$)', 'Interpreter', 'latex');
title('Resolution & Slit Visibility vs. FOV', 'Interpreter', 'latex');
legend('Measured Visibility (V)', 'Theoretical R_{z1} (nm)', 'Location', 'east');

exportgraphics(fig_summary, fullfile(output_dir, 'Summary_FOV_Sweep.png'), 'Resolution', 600);

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