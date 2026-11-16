clear; clc; close all;

%% 1. SYSTEM PARAMETERS SETUP
z1 = 350e-6;                    % Fixed distance z1 (m)
dz = 50e-6;                     % Fixed plane spacing z2 - z1 (m)
z2 = z1 + dz;

% Define sweep configurations: [lambda_0 (m), refractive_index]
configs = [
    405e-9, 1.000;              % Case 1: 405 nm in air
    405e-9, 1.515               % Case 2: 405 nm in immersion medium (n = 1.515)
];
short_labels = {'405nm (n=1)', '267nm (n=1.515)'};

area = 500e-6;                  % FOV (m)
dx = 100e-9;                    % Pixel pitch (m)
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
x_roi_um = 10.0;                 % Zoomed half-width (+/- 3 um) around center
idx_roi  = abs(x_um) <= x_roi_um;
x_cut    = x_um(idx_roi);

%% 2. CREATE DOUBLE-SLIT OBJECT
object_width_nm    = 200;       % Width of each slit (nm)
object_length_nm   = 5000;      % Length of each slit (nm)
center_distance_nm = 400;       % Center-to-center distance (nm) -> 200 nm gap

% Convert dimensions to meters
w_m = object_width_nm * 1e-9;
l_m = object_length_nm * 1e-9;
d_m = center_distance_nm * 1e-9;

% Create 2D spatial grid and object transmission mask (1 = transparent, 0 = opaque)
[X, Y] = meshgrid(x, y);
slit1 = (abs(X + d_m/2) <= w_m/2) & (abs(Y) <= l_m/2);
slit2 = (abs(X - d_m/2) <= w_m/2) & (abs(Y) <= l_m/2);
object_mask = double(slit1 | slit2);

% Print theoretical resolution limits for default configurations to Command Window
fprintf('=== THEORETICAL RESOLUTION LIMITS AT DEFAULT CONFIGURATION ===\n');
fprintf('Geometric NA at z1 (350 um): %.4f\n', NA_geo_z1_val);
fprintf('Geometric NA at z2 (400 um): %.4f\n', NA_geo_z2_val);
fprintf('Nyquist Limit (2*dx):        %.1f nm\n', 2 * dx * 1e9);
for i = 1:size(configs, 1)
    lam_eff = configs(i,1) / configs(i,2);
    res_z1  = lam_eff / (2 * NA_geo_z1_val);
    fprintf('Case %d [%s]: Diffraction limit at z1 = %.1f nm | System limit = %.1f nm\n', ...
        i, short_labels{i}, res_z1*1e9, max(res_z1, 2*dx)*1e9);
end

%% 3. CALCULATE & PLOT PARAMETER EFFECTS ON THEORETICAL RESOLUTION
% Baseline reference parameters
lambda0_ref = configs(1,1);     % 405 nm
n_ref       = configs(1,2);     % 1.000 (Air)
z_ref       = z1;               % 350 um
fov_ref     = area;             % 500 um
dx_ref      = dx;               % 100 nm

% Function for geometric diffraction-limited resolution: delta_diff = (lambda0 / n) / (2 * NA_geo)
calc_res_diff = @(lam, n, z, fov) (lam ./ n) ./ (2 .* (1 ./ sqrt(1 + (2 .* z ./ fov).^2)));

% Parameter sweep ranges
dx_range     = linspace(20e-9, 500e-9, 200);    % Pixel size: 20 nm -> 500 nm
fov_range    = linspace(100e-6, 2000e-6, 200);  % FOV: 100 um -> 2000 um
lambda_range = linspace(200e-9, 800e-9, 200);   % Wavelength: 200 nm -> 800 nm
z_range = linspace(5e-6, 1500e-6, 200); % Bắt đầu từ 5 um thay vì 50 um
% Initialize figure window
figure('Name', 'Parametric Resolution Analysis', 'Color', 'w', 'Position', [80, 80, 1150, 780]);

% --- Subplot 1: Effect of Pixel Size (dx) ---
subplot(2, 2, 1);
res_diff_const = calc_res_diff(lambda0_ref, n_ref, z_ref, fov_ref) * ones(size(dx_range));
res_nyquist    = 2 * dx_range;
res_total_1    = max(res_diff_const, res_nyquist);

plot(dx_range*1e9, res_diff_const*1e9, '--b', 'LineWidth', 1.5); hold on;
plot(dx_range*1e9, res_nyquist*1e9, '--r', 'LineWidth', 1.5);
plot(dx_range*1e9, res_total_1*1e9, 'k', 'LineWidth', 2.2);
xline(dx_ref*1e9, ':g', 'LineWidth', 1.8);
grid on;
title('1. Effect of Pixel Size (\Delta x)');
xlabel('Pixel size \Delta x (nm)');
ylabel('Theoretical resolution \delta (nm)');
legend('Geometric NA diffraction limit', 'Nyquist limit (2\Delta x)', 'System resolution', 'Current (100nm)', 'Location', 'northwest');

% --- Subplot 2: Effect of Field of View (FOV) ---
subplot(2, 2, 2);
res_diff_fov  = calc_res_diff(lambda0_ref, n_ref, z_ref, fov_range);
res_nyq_const = 2 * dx_ref * ones(size(fov_range));
res_total_2   = max(res_diff_fov, res_nyq_const);

plot(fov_range*1e6, res_diff_fov*1e9, '--b', 'LineWidth', 1.5); hold on;
plot(fov_range*1e6, res_nyq_const*1e9, '--r', 'LineWidth', 1.5);
plot(fov_range*1e6, res_total_2*1e9, 'k', 'LineWidth', 2.2);
xline(fov_ref*1e6, ':g', 'LineWidth', 1.8);
grid on;
title('2. Effect of Field of View (FOV)');
xlabel('FOV (\mum)');
ylabel('Theoretical resolution \delta (nm)');
legend('Geometric NA diffraction limit', 'Nyquist limit (2\Delta x)', 'System resolution', 'Current (500\mum)', 'Location', 'northeast');

% --- Subplot 3: Effect of Propagation Distance (z) ---
% --- Subplot 3: Effect of Propagation Distance (z) ---
subplot(2, 2, 3);
res_diff_z  = calc_res_diff(lambda0_ref, n_ref, z_range, fov_ref);
res_nyq_z   = 2 * dx_ref * ones(size(z_range));
res_total_3 = max(res_diff_z, res_nyq_z);

% 1. Vẽ đường tổng hợp nét liền dày bên dưới
plot(z_range*1e6, res_total_3*1e9, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 4.0); hold on;

% 2. Vẽ 2 đường giới hạn nét đứt mảnh đè lên trên
plot(z_range*1e6, res_diff_z*1e9, '--b', 'LineWidth', 1.6);
plot(z_range*1e6, res_nyq_z*1e9, '--r', 'LineWidth', 1.6);

set(gca, 'XScale', 'log');
xlim([5, 1500]);
xticks([5, 10, 36.2, 100, 350, 1000]);
xticklabels({'5', '10', '36.2', '100', '350', '1000'});

grid on;
title('3. Effect of Propagation Distance z (Log Scale)');
xlabel('Propagation distance z (\mum)');
ylabel('Theoretical resolution \delta (nm)');
legend('System resolution', 'Geometric NA diffraction limit', 'Nyquist limit (2\Delta x)', 'Location', 'northwest');

% --- Subplot 4: Effect of Wavelength (\lambda_0) & Refractive Index (n) ---
subplot(2, 2, 4);
res_diff_air  = calc_res_diff(lambda_range, 1.000, z_ref, fov_ref);
res_diff_oil  = calc_res_diff(lambda_range, 1.515, z_ref, fov_ref);
res_total_air = max(res_diff_air, 2 * dx_ref);
res_total_oil = max(res_diff_oil, 2 * dx_ref);

plot(lambda_range*1e9, res_total_air*1e9, 'b', 'LineWidth', 2.2); hold on;
plot(lambda_range*1e9, res_total_oil*1e9, 'c', 'LineWidth', 2.2);
plot(lambda_range*1e9, 2*dx_ref*ones(size(lambda_range))*1e9, '--r', 'LineWidth', 1.5);
xline(lambda0_ref*1e9, ':g', 'LineWidth', 1.8);
grid on;
title('4. Effect of Wavelength \lambda_0 & Medium n');
xlabel('Vacuum wavelength \lambda_0 (nm)');
ylabel('Theoretical resolution \delta (nm)');
legend('Resolution in Air (n = 1.000)', 'Glycerin (n = 1.515)', 'Nyquist limit (2\Delta x)', 'Current (405nm)', 'Location', 'northwest');

sgtitle('Theoretical Resolution Analysis (Abbe / Nyquist Criteria)', 'FontSize', 14, 'FontWeight', 'bold');

%% 4. VERIFY INITIALIZED OBJECT MASK IN ROI
figure('Name', 'Object Mask ROI', 'Color', 'w', 'Position', [150, 150, 500, 420]);
imagesc(x_cut, y_um(idx_roi), object_mask(idx_roi, idx_roi));
colormap(gray); axis image; colorbar;
title('Double-Slit Object (Width = 200nm, Gap = 200nm)');
xlabel('x (\mum)');
ylabel('y (\mum)');
