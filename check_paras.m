clc; clear; close all;

%% 1. SYSTEM PARAMETERS
lambda = 633e-9;       % Wavelength (m)
z1 = 1500e-6;          % Plane 1 (m)
z2 = 2000e-6;          % Plane 2 (m)
dz = z2 - z1;

%% 2. SENSOR / SAMPLING PARAMETERS
area = 1000e-6;        % FOV (m)
dx = 500e-9;           % Pixel pitch (m)
N = round(area/dx);
dy = dx;
z_max = N * dx^2 / lambda;

% Define Coordinates (chuẩn hóa chính xác cho lưới chẵn N)
x = (-N/2 : N/2 - 1) * dx;
y = (-N/2 : N/2 - 1) * dy;
x_um = x * 1e6;
y_um = y * 1e6;

fprintf('====================================\n');
fprintf('       SYSTEM PARAMETERS\n');
fprintf('====================================\n');
fprintf('Wavelength  : %.1f nm\n', lambda*1e9);
fprintf('Grid size   : %d x %d pixels\n', N, N);
fprintf('Pixel pitch : %.4f um\n', dx*1e6);
fprintf('FOV         : %.2f um\n', area*1e6);
fprintf('z1          : %.2f um\n', z1*1e6);
fprintf('z2          : %.2f um\n', z2*1e6);
fprintf('dz          : %.2f um\n', dz*1e6);
fprintf('z_max       : %.2f mm\n', z_max*1e3);

%% 3. CREATE OBJECT (Two Point Objects)
object_size_nm = 1000;       % Diameter / width of each object (nm)
center_distance_nm = 2000;  % Center-to-center distance (nm)
object_shape = 'square';    % 'circle' or 'square'

object_size = object_size_nm * 1e-9;
center_distance = center_distance_nm * 1e-9;

[X, Y] = meshgrid(x, y);

% Define object centers (symmetric around origin along X-axis)
x1 = -center_distance/2;
x2 =  center_distance/2;
y1 = 0;
y2 = 0;

object = ones(N, N);
switch lower(object_shape)
    case 'circle'
        mask1 = (X - x1).^2 + (Y - y1).^2 <= (object_size/2)^2;
        mask2 = (X - x2).^2 + (Y - y2).^2 <= (object_size/2)^2;
    case 'square'
        mask1 = abs(X - x1) <= object_size/2 & abs(Y - y1) <= object_size/2;
        mask2 = abs(X - x2) <= object_size/2 & abs(Y - y2) <= object_size/2;
    otherwise
        error('Unknown object shape.');
end
object(mask1 | mask2) = 0;

fprintf('\n========== OBJECT PARAMETERS ==========\n');
fprintf('Shape              : %s\n', object_shape);
fprintf('Object size        : %.1f nm\n', object_size_nm);
fprintf('Center distance    : %.1f nm\n', center_distance_nm);
fprintf('Edge-to-edge gap   : %.1f nm\n', center_distance_nm - object_size_nm);
fprintf('Object size (pixel): %.2f pixels\n', object_size/dx);
fprintf('Center distance(px): %.2f pixels\n', center_distance/dx);
fprintf('=======================================\n');

%% 4. TRANSMISSION FUNCTION
I0 = 1;
A0 = sqrt(I0);
t = object;
U0 = A0;
u_object = U0 * t;     % Trường phức ngay sau vật thể (biên độ nền = A0)
am = abs(u_object);

%% 5. DISPLAY OBJECT
% Mở rộng vùng hiển thị [-3, 3] um để nhìn trọn vẹn 2 vật thể tại +/- 1.2 um
% view_lim = 3; 
figure('Name', 'Object Details', 'Color', 'w');
imagesc(x_um, y_um, am);
axis image; colormap gray; colorbar;
% xlim([-view_lim, view_lim]);
% ylim([-view_lim, view_lim]);
title('Object Amplitude Distribution');
xlabel('x (\mum)'); ylabel('y (\mum)');

%% 6. ORIGINAL OBJECT CROSS-SECTION
[~, mid_y] = min(abs(y));
profile_obj_am = am(mid_y, :); 

figure('Name', 'Ground Truth Cross-section', 'Color', 'w');
stairs(x_um, am(mid_y,:), 'LineWidth', 2);
grid on; grid minor;
xlabel('x (\mum)'); ylabel('Amplitude (a.u.)');
title('Object Amplitude Cross-section (y = 0)');
% xlim([-view_lim, view_lim]);
% ylim([-0.5, A0 + 1]);

%% 7. HOLOGRAM SIMULATION
prop_z1 = Propagator_BLASM(N, lambda, dx, z1);
prop_z2 = Propagator_BLASM(N, lambda, dx, z2);

% SỬA LỖI: Truyền trường u_object (đã nhân A0) thay vì chỉ truyền t
field_z1_true = IFT2Dc(FT2Dc(u_object) .* prop_z1);
field_z2_true = IFT2Dc(FT2Dc(u_object) .* prop_z2);

hologram1 = abs(field_z1_true).^2; 
hologram2 = abs(field_z2_true).^2;

amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

%% 8. VISUALIZATION (HOLOGRAMS)
figure('Name', 'Measured Holograms', 'Color', 'w', 'Position', [100 100 1000 400]);
subplot(1,2,1); 
imagesc(x_um, y_um, hologram1); axis image;
title(sprintf('Hologram 1 ($z_1 = %d$ $\\mu$m)', round(z1*1e6)), 'Interpreter', 'latex'); 
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex'); 
colorbar; colormap gray;

subplot(1,2,2); 
imagesc(x_um, y_um, hologram2); axis image;
title(sprintf('Hologram 2 ($z_2 = %d$ $\\mu$m)', round(z2*1e6)), 'Interpreter', 'latex'); 
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex'); 
colorbar; colormap gray;

%% 9. HOLOGRAM CROSS-SECTION
profile_holo1 = hologram1(mid_y, :);
profile_holo2 = hologram2(mid_y, :);

figure('Name', 'Hologram Cross-Sections', 'Color', 'w', 'Position', [150 150 1000 400]);
subplot(1,2,1);
plot(x_um, profile_holo1, 'b', 'LineWidth', 1.2); grid on;
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('Intensity (a.u.)');
title(sprintf('Hologram 1 Cross-section ($z_1 = %.1f$ $\\mu$m)', z1*1e6), 'Interpreter', 'latex');

subplot(1,2,2);
plot(x_um, profile_holo2, 'r', 'LineWidth', 1.2); grid on;
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('Intensity (a.u.)');
title(sprintf('Hologram 2 Cross-section ($z_2 = %.1f$ $\\mu$m)', z2*1e6), 'Interpreter', 'latex');

%% 10. PHASE RETRIEVAL (Plane-to-Plane Gerchberg-Saxton)
Iterations = 100;
prop_dz       = Propagator_BLASM(N, lambda, dx,  dz);
prop_minus_dz = Propagator_BLASM(N, lambda, dx, -dz);
prop_minus_z1 = Propagator_BLASM(N, lambda, dx, -z1);

% Khởi tạo pha ban đầu bằng 0
u1 = amp1_measured; 

fprintf('\nRunning Phase Retrieval loop (%d iterations)...\n', Iterations);
figure('Name', 'Reconstruction Progress', 'Position', [200, 200, 1000, 450], 'Color', 'w');

for iter = 1:Iterations
    % Truyền từ z1 -> z2
    u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);
    u2 = amp2_measured .* exp(1i * angle(u2_sim)); 
    
    % Truyền ngược từ z2 -> z1
    u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz); 
    u1 = amp1_measured .* exp(1i * angle(u1_sim));
    
    if mod(iter, 10) == 0 || iter == 1 || iter == Iterations
        u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);
        
        subplot(1,2,1); 
        imagesc(x_um, y_um, abs(u_obj_rec)); axis image; colormap gray;
%         xlim([-view_lim, view_lim]); ylim([-view_lim, view_lim]);
        title(sprintf('Rec. Amplitude (Iter: %d)', iter)); colorbar;
        xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex');
        
        subplot(1,2,2); 
        imagesc(x_um, y_um, angle(u_obj_rec)); axis image; colormap gray;
%         xlim([-view_lim, view_lim]); ylim([-view_lim, view_lim]);
        title(sprintf('Rec. Phase (Iter: %d)', iter)); colorbar;
        xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex');
        drawnow;
    end
end
fprintf('Reconstruction complete!\n');

%% 11. CROSS-SECTION EVALUATION (Tự động tại y = 0 & Tùy chọn vẽ tay)
amp_rec = abs(u_obj_rec);
profile_rec_auto = amp_rec(mid_y, :);

% Đánh giá tự động qua tâm hai vật thể (tránh lỗi treo khi chưa vẽ tay)
figure('Name', 'Automatic Cross-section Evaluation', 'Color', 'w');
plot(x_um, profile_obj_am, 'Color', '#0072BD', 'LineWidth', 1.8); hold on;
plot(x_um, profile_rec_auto, 'Color', '#D95319', 'LineWidth', 1.5, 'LineStyle', '--');
grid on; 
% xlim([-view_lim, view_lim]);
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('Amplitude (a.u.)');
legend('Ground Truth', 'Reconstruction');
title('Amplitude Cross-section Evaluation (y = 0)');

% Nếu vẫn muốn vẽ đường cắt tương tác bằng chuột, mở comment khối dưới đây:

fig = figure('Color','w', 'Name','Select Cross-section', 'Position', [300, 300, 800, 600]);
ax = axes(fig);
imagesc(ax, x_um, y_um, amp_rec); axis(ax,'image'); colormap(ax,gray); colorbar(ax);
% xlim(ax, [-view_lim, view_lim]); ylim(ax, [-view_lim, view_lim]);
title(ax,'Reconstructed amplitude (Cross-section)');
h = drawline(ax,'Color','r','LineWidth',1.5);
wait(h);
if isvalid(h)
    pos = h.Position;
    x_line = linspace(pos(1,1), pos(2,1), 500);
    y_line = linspace(pos(1,2), pos(2,2), 500);
    profile_rec = interp2(x_um, y_um, amp_rec, x_line, y_line, 'linear');
    profile_true = interp2(x_um, y_um, am, x_line, y_line, 'linear'); 
    distance = hypot(x_line - x_line(1), y_line - y_line(1));
    figure('Color','w');
    plot(distance, profile_true, 'Color', '#0072BD', 'LineWidth', 1.5); hold on;
    plot(distance, profile_rec, 'Color', '#D95319', 'LineWidth', 1.5, 'LineStyle', '--');
    grid on; xlabel('Distance ($\mu$m)', 'Interpreter', 'latex'); ylabel('Amplitude (a.u.)');
    legend('Ground Truth', 'Reconstruction');
end


%% =========================================================================
%                               FUNCTIONS
% =========================================================================
function p = Propagator_BLASM(N, lambda, dx, z)
    % Lưới tần số không gian chuẩn hóa khớp với fftshift
    df = 1 / (N * dx);
    f = (-N/2 : N/2 - 1) * df;
    [FX, FY] = meshgrid(f, f);
    
    alpha = lambda * FX;
    beta  = lambda * FY;
    
    % Band-Limited Angular Spectrum Method (BL-ASM)
    limit = 1 / sqrt((2 * abs(z) / (N * dx))^2 + 1);
    mask = (alpha.^2 + beta.^2) <= limit^2;
    
    p = zeros(N, N);
    p(mask) = exp(1i * 2 * pi * (z / lambda) * sqrt(1 - alpha(mask).^2 - beta(mask).^2));
end

function out = FT2Dc(in)
    out = fftshift(fft2(ifftshift(in)));
end

function out = IFT2Dc(in)
    out = fftshift(ifft2(ifftshift(in)));
end