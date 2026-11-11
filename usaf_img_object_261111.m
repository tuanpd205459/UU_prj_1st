clc; clear; close all;

%% 1. SETTINGS & PARAMETERS
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
           'defaultTextFontName', 'Times New Roman', ...
           'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultTextInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex');

lambda = 633e-9;        % Wavelength (m)
z1 = 1500e-6;           % Distance to Plane 1 (m)
z2 = 1600e-6;           % Distance to Plane 2 (m)
dz = z2 - z1;

% Read USAF image
img = imread('USAF_220micro_104nm.png');

% Convert to grayscale
if ndims(img) == 3
    img = rgb2gray(img);
end

% ĐẢM BẢO ẢNH VUÔNG (Cắt ảnh về dạng N x N)
sz = min(size(img));
img = img(1:sz, 1:sz);

% Normalize image
object = double(img);
object = object / max(object(:));

% Image dimensions
[Ny, Nx] = size(object);
N = Nx;

% Pixel size
dx = 104e-9;
dy = dx;
area = dx * Nx; % Field of View

% Binary object
threshold = 0.5;
object_mask = object > threshold;

% Amplitude
amplitude = double(object_mask);

% Phase
phase_object = 1;  % rad
ph = double(object_mask) * phase_object;

% Complex transmission
T = amplitude .* exp(1i*ph);

%% 3. COORDINATE GRID & INFO
x = ((0:Nx-1)-(Nx-1)/2)*dx;
y = ((0:Ny-1)-(Ny-1)/2)*dy;
[X,Y] = meshgrid(x,y);

% KHAI BÁO x_um và y_um (bị thiếu ở code cũ)
x_um = x * 1e6;
y_um = y * 1e6;

FOV_x = Nx*dx;
FOV_y = Ny*dy;

fprintf('\n========== SIMULATION PARAMETERS ==========\n');
fprintf('Object mode    : USAF Target Image\n');
fprintf('Grid size      : %d x %d pixels\n', Nx, Ny);
fprintf('Pixel size     : %.2f nm\n', dx*1e9);
fprintf('FOV            : %.3f x %.3f um\n', FOV_x*1e6, FOV_y*1e6);
fprintf('Wavelength     : %.1f nm\n', lambda*1e9);
fprintf('Plane z1       : %.3f um\n', z1*1e6);
fprintf('Plane z2       : %.3f um\n', z2*1e6);
fprintf('Distance dz    : %.3f um\n', dz*1e6);


%% 4. HOLOGRAM SIMULATION
prop_z1 = Propagator(N, lambda, area, z1);
prop_z2 = Propagator(N, lambda, area, z2);

field_z1_true = IFT2Dc(FT2Dc(T) .* prop_z1);
field_z2_true = IFT2Dc(FT2Dc(T) .* prop_z2);

% Cường độ (hologram) và biên độ đo được
hologram1 = abs(field_z1_true).^2; 
hologram2 = abs(field_z2_true).^2;

amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

%% 5. VISUALIZATION (OBJECT & HOLOGRAMS)
figure('Name', 'Object and Measured Holograms', 'Color', 'w', 'Position', [100 100 1200 800]);

% Object amplitude
subplot(2,2,1);
imagesc(x_um, y_um, amplitude); axis image;
title('Object Amplitude (USAF Target)'); xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)');
colorbar; clim([0 1]); colormap gray;

% Object phase
subplot(2,2,2);
imagesc(x_um, y_um, ph); axis image;
title('Object Phase'); xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)');
colorbar; colormap gray;

% Hologram 1
subplot(2,2,3); 
imagesc(x_um, y_um, hologram1); axis image;
title(sprintf('Hologram 1 ($z_1 = %d$ $\\mu$m)', round(z1*1e6))); 
xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)'); colorbar; colormap gray;

% Hologram 2
subplot(2,2,4); 
imagesc(x_um, y_um, hologram2); axis image;
title(sprintf('Hologram 2 ($z_2 = %d$ $\\mu$m)', round(z2*1e6))); 
xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)'); colorbar; colormap gray;

%% 6. PHASE RETRIEVAL (Plane-to-Plane)
Iterations = 300;
prop_dz = Propagator(N, lambda, area, dz);
prop_minus_dz = Propagator(N, lambda, area, -dz);
prop_minus_z1 = Propagator(N, lambda, area, -z1);

u1 = amp1_measured .* exp(1i * zeros(N, N));

fprintf('Đang chạy vòng lặp phục hồi pha (Phase Retrieval) %d iterations...\n', Iterations);
figure('Name', 'Reconstruction Progress', 'Position', [200, 200, 1000, 450]);

for iter = 1:Iterations
    % Lan truyền z1 -> z2
    u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);
    u2 = amp2_measured .* exp(1i * angle(u2_sim)); 
    
    % Lan truyền ngược z2 -> z1
    u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz); 
    u1 = amp1_measured .* exp(1i * angle(u1_sim));
    
    % Hiển thị tiến trình
    if mod(iter, 10) == 0 || iter == 1 || iter == Iterations
        u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);
        
        subplot(1,2,1); imagesc(x_um, y_um, abs(u_obj_rec)); axis image; colormap gray;
        title(sprintf('Rec. Amplitude (Iter: %d)', iter)); colorbar;
        xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)');
        
        subplot(1,2,2); imagesc(x_um, y_um, angle(u_obj_rec)); axis image; colormap gray;
        title(sprintf('Rec. Phase (Iter: %d)', iter)); colorbar;
        xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)');
        drawnow;
    end
end
fprintf('Hoàn thành phục hồi!\n');

%% 7. CROSS-SECTION EVALUATION
amp_rec = abs(u_obj_rec);

fig = figure('Color','w', 'Name','Select Cross-section');
ax = axes(fig);
imagesc(ax, x_um, y_um, amp_rec); axis(ax,'image'); colormap(ax,gray); colorbar(ax);
xlabel(ax,'$x$ ($\mu$m)'); ylabel(ax,'$y$ ($\mu$m)');
title(ax,'Draw a line between TWO points');

h = drawline(ax,'Color','r','LineWidth',1.5);
wait(h);

pos = h.Position;
x_click = pos(:,1)';
y_click = pos(:,2)';

num_points = 500;
x_line = linspace(x_click(1), x_click(2), num_points);
y_line = linspace(y_click(1), y_click(2), num_points);

[X_grid, Y_grid] = meshgrid(x_um, y_um);
profile_rec = interp2(X_grid, Y_grid, amp_rec, x_line, y_line, 'linear');
profile_true = interp2(X_grid, Y_grid, amplitude, x_line, y_line, 'linear');
distance = hypot(x_line - x_line(1), y_line - y_line(1));

figure('Color','w');
plot(distance, profile_true, 'Color', '#0072BD', 'LineWidth', 1.5); hold on;
plot(distance, profile_rec, 'Color', '#D95319', 'LineWidth', 1.5);
grid on;
xlabel('Distance ($\mu$m)'); ylabel('Amplitude');
legend('Ground Truth', 'Reconstruction');
title('Amplitude Cross-section');

%% =========================================================================
%                               FUNCTIONS
% =========================================================================
function p = Propagator(N, lambda, area, z)
    [X, Y] = meshgrid(1:N, 1:N);
    alpha = lambda * (Y - N/2 - 1) / area;
    beta  = lambda * (X - N/2 - 1) / area;
    mask = (alpha.^2 + beta.^2) <= 1;
    
    p = zeros(N, N);
    p(mask) = exp(2 * pi * 1i * z * sqrt(1 - alpha(mask).^2 - beta(mask).^2) / lambda);
end

function out = FT2Dc(in)
    [Nx, Ny] = size(in);
    vx = (-1).^(1:Nx).'; 
    vy = (-1).^(1:Ny);   
    f1 = vx * vy;        
    out = f1 .* fft2(f1 .* in);
end

function out = IFT2Dc(in)
    [Nx, Ny] = size(in);
    vx = (-1).^(1:Nx).'; 
    vy = (-1).^(1:Ny);   
    f1 = vx * vy;        
    out = f1 .* ifft2(f1 .* in);
end