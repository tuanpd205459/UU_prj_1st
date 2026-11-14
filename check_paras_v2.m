clc; clear; close all;

%% 1. SYSTEM PARAMETERS
lambda = 633e-9;       % Wavelength (m)
z1 = 1500e-6;          % Plane 1 (m)
z2 = 2000e-6;          % Plane 2 (m)
dz = z2 - z1;

%% 2. SENSOR / SAMPLING PARAMETERS
area = 500e-6;        % FOV (m)
dx = 100e-9;           % Pixel pitch (m)
N = round(area/dx);
dy = dx;

z_max = N * dx^2 / lambda; % Fresnel/ASM aliasing transition distance

% Critical distance z_crit to achieve Nyquist resolution limit (R = dx)
NA_Nyq = lambda / (2 * dx);
if NA_Nyq < 1
    z_crit = (area / 2) * sqrt(1 - NA_Nyq^2) / NA_Nyq;
    % Tuong duong: z_crit = (N * dx^2 / lambda) * sqrt(1 - (lambda/(2*dx))^2);
else
    z_crit = 0; % Evanescent regime (dx <= lambda/2)
end

% Define coordinates
x = (-N/2 : N/2 - 1) * dx;
y = (-N/2 : N/2 - 1) * dy;

x_um = x * 1e6;
y_um = y * 1e6;

fprintf('Wavelength  : %.1f nm\n', lambda*1e9);
fprintf('Grid size   : %d x %d pixels\n', N, N);
fprintf('Pixel pitch : %.4f um\n', dx*1e6);
fprintf('FOV         : %.2f um\n', area*1e6);
fprintf('z1          : %.2f um\n', z1*1e6);
fprintf('z2          : %.2f um\n', z2*1e6);
fprintf('dz          : %.2f um\n', dz*1e6);
fprintf('z_max       : %.2f mm\n', z_max*1e3);
fprintf('z_crit      : %.2f um (Max z for R = dx)\n', z_crit*1e6);

%% 3. DIFFRACTION LIMIT CALCULATION
NA1 = sin(atan(area/(2*z1)));
NA2 = sin(atan(area/(2*z2)));

% Effective NA (bounded by both sensor aperture and Nyquist sampling)
NA1_eff = min(NA1, NA_Nyq);
NA2_eff = min(NA2, NA_Nyq);

d_Abbe_1 = lambda / (2*NA1_eff);
d_Abbe_2 = lambda / (2*NA2_eff);

fprintf('NA_Nyq                 : %.5f\n', NA_Nyq);
fprintf('NA at z1               : %.5f (Effective: %.5f)\n', NA1, NA1_eff);
fprintf('NA at z2               : %.5f (Effective: %.5f)\n', NA2, NA2_eff);
fprintf('Resolution at z1       : %.2f nm\n', d_Abbe_1*1e9);
fprintf('Resolution at z2       : %.2f nm\n', d_Abbe_2*1e9);

if z1 > z_crit || z2 > z_crit
    fprintf('[WARNING] z1 or z2 > z_crit: Resolution is aperture-limited (> %.1f nm).\n', dx*1e9);
else
    fprintf('[OK] Both planes <= z_crit: Nyquist limit (R = %.1f nm) is achievable.\n', dx*1e9);
end


%% 4. CREATE OBJECT (Two Point Objects)
object_shape = 'rectangle';      % 'circle' or 'square'

object_width_nm  = 500;       % Width of each slit (nm)
object_length_nm = 5000;      % Length of each slit (nm)

center_distance_nm = 500 + 500;

object_width  = object_width_nm  * 1e-9;
object_length = object_length_nm * 1e-9;
center_distance = center_distance_nm * 1e-9;

[X, Y] = meshgrid(x, y);

% Define object centers
x1_obj = -center_distance/2;
x2_obj =  center_distance/2;

y1_obj = 0;
y2_obj = 0;

object = ones(N, N);

% Rectangular elongated slits
mask1 = abs(X - x1_obj) <= object_width/2 & ...
        abs(Y - y1_obj) <= object_length/2;

mask2 = abs(X - x2_obj) <= object_width/2 & ...
        abs(Y - y2_obj) <= object_length/2;

object(mask1 | mask2) = 0;

fprintf('\n========== OBJECT PARAMETERS ==========\n');
fprintf('Shape              : Elongated rectangle\n');
fprintf('Width              : %.1f nm\n', object_width_nm);
fprintf('Length             : %.1f nm\n', object_length_nm);
fprintf('Center distance    : %.1f nm\n', center_distance_nm);
fprintf('Edge-to-edge gap   : %.1f nm\n', ...
        center_distance_nm - object_width_nm);
fprintf('Width (pixels)     : %.2f pixels\n', object_width/dx);
fprintf('Length (pixels)    : %.2f pixels\n', object_length/dx);
fprintf('=======================================\n');

%% 5. TRANSMISSION FUNCTION

I0 = 1;
A0 = sqrt(I0);

t = object;
U0 = A0;

u_object = U0 * t;

am = abs(u_object);


%% 6. DISPLAY OBJECT

figure('Name', 'Object Details', 'Color', 'w');

imagesc(x_um, y_um, am);

axis image;
colormap gray;
colorbar;

title('Object Amplitude Distribution');
xlabel('x (\mum)');
ylabel('y (\mum)');


%% 7. ORIGINAL OBJECT CROSS-SECTION

[~, mid_y] = min(abs(y));

profile_obj_am = am(mid_y, :);

figure('Name', 'Ground Truth Cross-section', 'Color', 'w');

stairs(x_um, profile_obj_am, 'LineWidth', 2);

grid on;
grid minor;

xlabel('x (\mum)');
ylabel('Amplitude (a.u.)');

title('Ground Truth Amplitude Cross-section (y = 0)');


%% 8. HOLOGRAM SIMULATION

prop_z1 = Propagator_BLASM(N, lambda, dx, z1);
prop_z2 = Propagator_BLASM(N, lambda, dx, z2);

% Propagation from object plane to z1 and z2

field_z1_true = IFT2Dc(FT2Dc(u_object) .* prop_z1);

field_z2_true = IFT2Dc(FT2Dc(u_object) .* prop_z2);

% Calculate holograms

hologram1 = abs(field_z1_true).^2;
hologram2 = abs(field_z2_true).^2;

% Measured amplitudes

amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);


%% 9. VISUALIZATION (HOLOGRAMS)

figure('Name', 'Measured Holograms', ...
       'Color', 'w', ...
       'Position', [100 100 1000 400]);

subplot(1,2,1);

imagesc(x_um, y_um, hologram1);

axis image;
title(sprintf('Hologram 1 (z_1 = %d \\mum)', round(z1*1e6)));

xlabel('x (\mum)');
ylabel('y (\mum)');

colorbar;
colormap gray;


subplot(1,2,2);

imagesc(x_um, y_um, hologram2);

axis image;
title(sprintf('Hologram 2 (z_2 = %d \\mum)', round(z2*1e6)));

xlabel('x (\mum)');
ylabel('y (\mum)');

colorbar;
colormap gray;


%% 10. HOLOGRAM CROSS-SECTION

profile_holo1 = hologram1(mid_y, :);
profile_holo2 = hologram2(mid_y, :);

figure('Name', 'Hologram Cross-Sections', ...
       'Color', 'w', ...
       'Position', [150 150 1000 400]);

subplot(1,2,1);

plot(x_um, profile_holo1, 'b', 'LineWidth', 1.2);

grid on;
grid minor;

xlabel('x (\mum)');
ylabel('Intensity (a.u.)');

title(sprintf('Hologram 1 Cross-section (z_1 = %.1f \\mum)', z1*1e6));


subplot(1,2,2);

plot(x_um, profile_holo2, 'r', 'LineWidth', 1.2);

grid on;
grid minor;

xlabel('x (\mum)');
ylabel('Intensity (a.u.)');

title(sprintf('Hologram 2 Cross-section (z_2 = %.1f \\mum)', z2*1e6));


%% 11. PHASE RETRIEVAL
% Plane-to-Plane Gerchberg-Saxton Algorithm

Iterations = 100;

prop_dz       = Propagator_BLASM(N, lambda, dx, dz);
prop_minus_dz = Propagator_BLASM(N, lambda, dx, -dz);
prop_minus_z1 = Propagator_BLASM(N, lambda, dx, -z1);

% Initial phase

u1 = amp1_measured;

fprintf('\nRunning Phase Retrieval loop (%d iterations)...\n', Iterations);

figure('Name', 'Reconstruction Progress', ...
       'Position', [200 200 1000 450], ...
       'Color', 'w');


for iter = 1:Iterations

    % Propagation from z1 to z2

    u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);

    u2 = amp2_measured .* exp(1i * angle(u2_sim));


    % Back propagation from z2 to z1

    u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz);

    u1 = amp1_measured .* exp(1i * angle(u1_sim));


    % Display reconstruction progress

    if mod(iter, 10) == 0 || iter == 1 || iter == Iterations

        u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);

        subplot(1,2,1);

        imagesc(x_um, y_um, abs(u_obj_rec));

        axis image;
        colormap gray;

        title(sprintf('Reconstructed Amplitude (Iter: %d)', iter));

        colorbar;

        xlabel('x (\mum)');
        ylabel('y (\mum)');


        subplot(1,2,2);

        imagesc(x_um, y_um, angle(u_obj_rec));

        axis image;
        colormap gray;

        title(sprintf('Reconstructed Phase (Iter: %d)', iter));

        colorbar;

        xlabel('x (\mum)');
        ylabel('y (\mum)');

        drawnow;

    end

end

fprintf('Reconstruction complete!\n');


%% 12. INTERACTIVE CROSS-SECTION EVALUATION

% Reconstructed amplitude

amp_rec = abs(u_obj_rec);


% Display reconstructed amplitude

fig = figure('Color', 'w', ...
             'Name', 'Select Cross-section', ...
             'Position', [300 300 800 600]);

ax = axes(fig);

imagesc(ax, x_um, y_um, amp_rec);

axis(ax, 'image');
colormap(ax, gray);
colorbar(ax);

title(ax, 'Reconstructed Amplitude - Select Cross-section');

xlabel(ax, 'x (\mum)');
ylabel(ax, 'y (\mum)');


% Select cross-section using mouse

h = drawline(ax, 'Color', 'r', 'LineWidth', 1.5);

wait(h);


if isvalid(h)

    % Get selected line coordinates

    pos = h.Position;

    x_line = linspace(pos(1,1), pos(2,1), 500);
    y_line = linspace(pos(1,2), pos(2,2), 500);


    % Ground Truth
    % Nearest-neighbor interpolation to preserve pixel values

    profile_true = interp2(x_um, y_um, am, ...
                           x_line, y_line, 'nearest');


    % Reconstruction
    % Linear interpolation for continuous profile

    profile_rec = interp2(x_um, y_um, amp_rec, ...
                          x_line, y_line, 'linear');


    % Distance along selected cross-section

    distance = hypot(x_line - x_line(1), ...
                     y_line - y_line(1));


    % Plot comparison

    figure('Color', 'w', ...
           'Name', 'Cross-section Comparison');

    % Ground Truth: Staircase plot

    stairs(distance, profile_true, ...
           'Color', '#0072BD', ...
           'LineWidth', 1.5);

    hold on;


    % Reconstruction: Continuous plot

    plot(distance, profile_rec, ...
         'Color', '#D95319', ...
         'LineWidth', 1.5, ...
         'LineStyle', '--');


    grid on;
    grid minor;

    xlabel('Distance (\mum)');
    ylabel('Amplitude (a.u.)');

    legend('Ground Truth', 'Reconstruction');

    title('Cross-section Comparison');

end


%% =========================================================================
%                               FUNCTIONS
% =========================================================================

function p = Propagator_BLASM(N, lambda, dx, z)

    % Spatial frequency grid

    df = 1 / (N * dx);

    f = (-N/2 : N/2 - 1) * df;

    [FX, FY] = meshgrid(f, f);


    alpha = lambda * FX;
    beta  = lambda * FY;


    % Band-Limited Angular Spectrum Method

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
