clc; clear; close all;
tic

%% 1. SYSTEM & SENSOR PARAMETERS (FIXED ACROSS SWEEP)
lambda = 405e-9;                % Wavelength (m)
diffractive_index = 1;
lambda = lambda / diffractive_index;

dz = 500e-6;                    % Fixed plane spacing z2 - z1 (m)
z1_list_um = [100, 250, 500, 1000, 1500, 2000]; % Sweep values for z1 (um)
z1_list = z1_list_um * 1e-6;    % Convert to meters

area = 1000e-6;                  % FOV (m)
dx = 200e-9;                    % Pixel pitch (m)
N = round(area/dx);
dy = dx;

z_max = N * dx^2 / lambda;      % Fresnel/ASM aliasing transition distance
NA_Nyq = lambda / (2 * dx);

if NA_Nyq < 1
    z_crit = (area / 2) * sqrt(1 - NA_Nyq^2) / NA_Nyq;
else
    z_crit = 0;                 % Evanescent regime (dx <= lambda/2)
end

% Define spatial coordinates
x = (-N/2 : N/2 - 1) * dx;
y = (-N/2 : N/2 - 1) * dy;
x_um = x * 1e6;
y_um = y * 1e6;
[~, mid_y] = min(abs(y));

% ROI settings for cross-section evaluation at y = 0
x_roi_um = 5;                   % Zoomed half-width (um) around center
idx_roi  = abs(x_um) <= x_roi_um;
x_cut    = x_um(idx_roi);

%% 2. DISPLAY FIXED SIMULATION PARAMETERS
fprintf('=====================================================\n');
fprintf('           FIXED SIMULATION PARAMETERS               \n');
fprintf('=====================================================\n');
fprintf('Wavelength       : %.1f nm\n', lambda*1e9);
fprintf('Refractive index : %.2f\n', diffractive_index);
fprintf('Grid size        : %d x %d pixels\n', N, N);
fprintf('Pixel pitch      : %.4f um (%.1f nm)\n', dx*1e6, dx*1e9);
fprintf('FOV              : %.2f um\n', area*1e6);
fprintf('Plane spacing dz : %.2f um\n', dz*1e6);
fprintf('z1 sweep range   : [%s] um\n', num2str(z1_list_um));
fprintf('z_max (Aliasing) : %.2f mm (%.1f um)\n', z_max*1e3, z_max*1e6);
fprintf('z_crit (R = dx)  : %.2f um\n', z_crit*1e6);
fprintf('NA_Nyq           : %.5f\n', NA_Nyq);
fprintf('=====================================================\n');

%% 3. CREATE OBJECT (COMPUTED ONCE OUTSIDE LOOP)
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

fprintf('\n================ OBJECT PARAMETERS ==================\n');
fprintf('Shape            : Elongated double rectangle\n');
fprintf('Width            : %.1f nm (%.2f pixels)\n', object_width_nm, object_width/dx);
fprintf('Length           : %.1f nm (%.2f pixels)\n', object_length_nm, object_length/dx);
fprintf('Center distance  : %.1f nm\n', center_distance_nm);
fprintf('Edge-to-edge gap : %.1f nm\n', center_distance_nm - object_width_nm);
fprintf('=====================================================\n');

% Create output directory for saving figures & logs
out_dir = 'Sweep_Results_z';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

%% 4. PRE-ALLOCATE STORAGE FOR SWEEP METRICS
num_z = length(z1_list);
z2_list_um     = z1_list_um + dz*1e6;
NA_geo_z1_arr  = zeros(1, num_z);
NA_geo_z2_arr  = zeros(1, num_z);
NA_eff_z1_arr  = zeros(1, num_z);
NA_eff_z2_arr  = zeros(1, num_z);
R_z1_nm_arr    = zeros(1, num_z);
R_z2_nm_arr    = zeros(1, num_z);
visibility_arr = zeros(1, num_z);
regime_list    = strings(1, num_z);
profiles_rec   = zeros(num_z, length(x_cut));

% Indices to measure dip (center x=0) and slit minima (x = +/- 200 nm)
[~, idx_gap]   = min(abs(x_cut - 0));
[~, idx_slit1] = min(abs(x_cut + 0.2));
[~, idx_slit2] = min(abs(x_cut - 0.2));

Iterations = 100;

%% 5. PARAMETER SWEEP LOOP OVER z1
fprintf('\nStarting sweep over %d distances (Iterations per z: %d)...\n', num_z, Iterations);

for k = 1:num_z
    z1 = z1_list(k);
    z2 = z1 + dz;
    
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
    
    % Determine limiting regime
    if NA_Nyq < min(NA_geo_z1, NA_geo_z2)
        regime_list(k) = "Sampling-limited";
    elseif max(z1, z2) > z_crit
        regime_list(k) = "Aperture-limited";
    else
        regime_list(k) = "Diffraction-limited";
    end
    
    fprintf('[%d/%d] z1 = %4.0f um | z2 = %4.0f um | NA_z1 = %.4f (R1 = %6.1f nm) | NA_z2 = %.4f (R2 = %6.1f nm) | %s\n', ...
            k, num_z, z1*1e6, z2*1e6, NA_eff_z1, R_z1*1e9, NA_eff_z2, R_z2*1e9, regime_list(k));
        
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
    
    % Extract cross-section at y = 0
    prof_rec = amp_rec(mid_y, idx_roi);
    profiles_rec(k, :) = prof_rec;
    
    % Calculate contrast/visibility between the two dark slits
    I_peak_gap = prof_rec(idx_gap);
    I_trough   = mean([prof_rec(idx_slit1), prof_rec(idx_slit2)]);
    visibility_arr(k) = max(0, (I_peak_gap - I_trough) / (I_peak_gap + I_trough));
    
    % Save individual evaluation figure (invisible mode for fast execution)
    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1100 450]);
    
    subplot(1, 2, 1);
    imagesc(x_um, y_um, amp_rec);
    axis image;
    set(gca, 'YDir', 'normal', 'FontName', 'Times New Roman', 'FontSize', 11);
    colormap gray; colorbar; hold on;
    plot([min(x_cut), max(x_cut)], [0, 0], 'r-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    xlim([-15 15]); ylim([-15 15]);
    title(sprintf('Reconstructed Amplitude (z_1 = %.0f \\mum)', z1*1e6));
    xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel('$y\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    
    subplot(1, 2, 2);
    stairs(x_cut, profile_true, 'Color', '#0072BD', 'LineWidth', 1.8); hold on;
    plot(x_cut, prof_rec, 'Color', '#D95319', 'LineWidth', 1.8, 'LineStyle', '--');
    grid on; grid minor;
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
    xlim([-x_roi_um, x_roi_um]);
    ylim([0, 1.2]);
    xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
    ylabel('Amplitude (a.u.)', 'Interpreter', 'latex');
    legend('Ground Truth', sprintf('Rec (R_{theo} = %.1f nm)', R_z1*1e9), 'Location', 'southeast');
    title(sprintf('Cross-section at y = 0 (NA = %.3f, V = %.2f)', NA_eff_z1, visibility_arr(k)));
    
    saveas(fig, fullfile(out_dir, sprintf('Rec_z1_%04dum.png', round(z1*1e6))));
    close(fig);
end

elapsed_time = toc;


%% 7. SUMMARY PLOTS ACROSS ALL z1 VALUES
figure('Color', 'w', 'Name', 'Sweep Summary vs z1', 'Position', [150 150 1200 450]);

% Panel 1: Overlay of cross-sections at y = 0
subplot(1, 2, 1);
stairs(x_cut, profile_true, 'k-', 'LineWidth', 2, 'DisplayName', 'Ground Truth'); hold on;
colors = lines(num_z);
for k = 1:num_z
    plot(x_cut, profiles_rec(k, :), 'LineWidth', 1.5, 'Color', colors(k, :), ...
         'DisplayName', sprintf('z_1 = %d \\mum (R=%.0f nm)', z1_list_um(k), R_z1_nm_arr(k)));
end
grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlim([-2.5, 2.5]);
ylim([0, 1.2]);
xlabel('$x\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
ylabel('Amplitude (a.u.)', 'Interpreter', 'latex');
title('Cross-section Reconstructed Amplitude vs. Propagation Distance z_1');
legend('Location', 'southeast', 'FontSize', 9);

% Panel 2: Theoretical Resolution & Measured Slit Contrast vs z1
subplot(1, 2, 2);
yyaxis left;
plot(z1_list_um, R_z1_nm_arr, '-o', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Theoretical Resolution $R$ (nm)', 'Interpreter', 'latex');
yline(200, 'r--', 'Object Gap (200 nm)', 'LineWidth', 1.2);

yyaxis right;
plot(z1_list_um, visibility_arr, '-s', 'LineWidth', 2, 'MarkerSize', 6);
ylabel('Measured Slit Contrast (Visibility)', 'Interpreter', 'latex');

grid on; grid minor;
set(gca, 'FontName', 'Times New Roman', 'FontSize', 11);
xlabel('Propagation Distance $z_1\ (\mu\mathrm{m})$', 'Interpreter', 'latex');
title('Resolution Limit & Contrast vs. Distance z_1');

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