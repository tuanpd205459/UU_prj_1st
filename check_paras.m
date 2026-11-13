clc; clear; close all;

%% 1. SETTINGS & PARAMETERS
set(groot, 'defaultAxesFontName', 'Times New Roman', ...
           'defaultTextFontName', 'Times New Roman', ...
           'defaultAxesTickLabelInterpreter', 'latex', ...
           'defaultTextInterpreter', 'latex', ...
           'defaultLegendInterpreter', 'latex');

% Thông số hệ thống
lambda = 633e-9;         % Bước sóng: 633 nm
z1 = 1500e-6;              % SỬA LỖI: Thu hẹp z1 về 10 um để phù hợp với FOV nano
z2 = 2000e-6;              % SỬA LỖI: Thu hẹp z2 về 15 um
dz = z2 - z1;            

% Thông số cảm biến
dx = 200e-9;             % Kích thước pixel: 100 nm (0.1 um)
dy = dx;
% N = 1024;                 % Số lượng pixel (N x N)
% area = N * dx;           % FOV (m)
area = 110e-6;
N = round(area/dx);
%% 2. TẠO VẬT THỂ (Custom Resolution Target)
% Tùy chỉnh kích thước vạch (micromet)
bar_width_nm = 200;      % Độ rộng mỗi vạch (200 nm = 2 px)
bar_gap_nm = 1000;        % Khoảng cách giữa các vạch (400 nm = 4 px)
bar_length_um = 20;      % SỬA LỖI: Giảm chiều dài xuống 20 um (200 px) để lọt vào FOV 51.2 um
num_bars = 3;            % Số lượng vạch (chuẩn 3-bar target)

% Chuyển đổi sang số lượng pixel
bar_width_px = round((bar_width_nm * 1e-9) / dx);
bar_gap_px = round((bar_gap_nm * 1e-9) / dx);
bar_length_px = round((bar_length_um * 1e-6) / dy);

% Khởi tạo nền tối (giá trị 0)
object = zeros(N, N);

% Tính toán vị trí bắt đầu để căn giữa ma trận
total_width_px = num_bars * bar_width_px + (num_bars - 1) * bar_gap_px;
start_x = round(N/2) - floor(total_width_px / 2);
start_y = round(N/2) - floor(bar_length_px / 2);

% Vẽ các vạch (giá trị 1)
for i = 1:num_bars
    current_x = start_x + (i-1) * (bar_width_px + bar_gap_px);
    % Tránh vượt quá giới hạn ma trận
    if current_x + bar_width_px - 1 <= N && start_y + bar_length_px - 1 <= N
        object(start_y : start_y + bar_length_px - 1, current_x : current_x + bar_width_px - 1) = 1;
    end
end

% Hàm truyền đạt phức (Transmission function)
am = exp(-1.6 * object);      
ph0 = 3;                      
ph = ph0 * object;            
t = am .* exp(1i * ph);       

% Tọa độ hiển thị (micromet)
x_um = ((1:N) - N/2 - 1) * dx * 1e6; % SỬA LỖI: 1e6 để ra micromet
y_um = ((1:N) - N/2 - 1) * dy * 1e6; % SỬA LỖI: 1e6 để ra micromet

% --- Hiển thị Biên độ và Pha ban đầu ---
figure('Name', 'Original Object (Ground Truth)', 'Position', [100, 550, 800, 350], 'Color', 'w');
subplot(1,2,1); imagesc(x_um, y_um, am); colormap gray; axis image;
title('Original Amplitude'); xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)'); colorbar; 
subplot(1,2,2); imagesc(x_um, y_um, ph); colormap gray; axis image;
title('Original Phase'); xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)'); colorbar; 

%% 3. ORIGINAL OBJECT CROSS-SECTION
mid_y = round(N/2); 
profile_obj_am = am(mid_y, :);
profile_obj_ph = ph(mid_y, :);

figure('Name', 'Ground Truth Cross-Section', 'Color', 'w', 'Position', [150 500, 1000, 400]);
subplot(1,2,1);
plot(x_um, profile_obj_am, 'k', 'LineWidth', 1.5);
grid on; grid minor; 
xlabel('$x$ ($\mu$m)'); ylabel('Amplitude');
title('Amplitude Cross-section at center ($y=0$)');
xlim([-total_width_px*dx*1e6 * 2, total_width_px*dx*1e6 * 2]); % Zoom vừa phải

subplot(1,2,2);
plot(x_um, profile_obj_ph, 'b', 'LineWidth', 1.5);
grid on; grid minor;
xlabel('$x$ ($\mu$m)'); ylabel('Phase (rad)');
title('Phase Cross-section at center ($y=0$)');
xlim([-total_width_px*dx*1e6 * 2, total_width_px*dx*1e6 * 2]);

%% 4. DISPLAY SYSTEM PARAMETERS
fprintf('\n====================================\n');
fprintf('        SIMULATION PARAMETERS\n');
fprintf('====================================\n');
fprintf('Wavelength       : %.1f nm\n', lambda*1e9);
fprintf('Plane z1         : %.1f um\n', z1*1e6);
fprintf('Plane z2         : %.1f um\n', z2*1e6);
fprintf('Distance dz      : %.1f um\n', dz*1e6);
fprintf('\n[Target Generation]\n');
fprintf('Bar Width        : %.1f nm (%d px)\n', bar_width_nm, bar_width_px);
fprintf('Bar Gap          : %.1f nm (%d px)\n', bar_gap_nm, bar_gap_px);
fprintf('Bar Length       : %.1f um (%d px)\n', bar_length_um, bar_length_px);
fprintf('\n[Sampling & FOV]\n');
fprintf('Pixel size       : %.3f um (%.1f nm)\n', dx*1e6, dx*1e9);
fprintf('Grid size        : %d x %d pixels\n', N, N);
fprintf('Total FOV        : %.1f x %.1f um^2\n', area*1e6, area*1e6);
fprintf('====================================\n\n');

%% 5. HOLOGRAM SIMULATION
prop_z1 = Propagator_BLASM(N, lambda, area, dx, z1);
prop_z2 = Propagator_BLASM(N, lambda, area, dx, z2);

field_z1_true = IFT2Dc(FT2Dc(t) .* prop_z1);
field_z2_true = IFT2Dc(FT2Dc(t) .* prop_z2);

hologram1 = abs(field_z1_true).^2; 
hologram2 = abs(field_z2_true).^2;

amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

%% 6. VISUALIZATION (HOLOGRAMS)
figure('Name', 'Measured Holograms', 'Color', 'w', 'Position', [100 100 1000 400]);
subplot(1,2,1); 
imagesc(x_um, y_um, hologram1); axis image;
title(sprintf('Hologram 1 ($z_1 = %d$ $\\mu$m)', round(z1*1e6))); 
xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)'); colorbar; colormap gray;
subplot(1,2,2); 
imagesc(x_um, y_um, hologram2); axis image;
title(sprintf('Hologram 2 ($z_2 = %d$ $\\mu$m)', round(z2*1e6))); 
xlabel('$x$ ($\mu$m)'); ylabel('$y$ ($\mu$m)'); colorbar; colormap gray;

%% 7. HOLOGRAM CROSS-SECTION
profile_holo1 = hologram1(mid_y, :);
profile_holo2 = hologram2(mid_y, :);

figure('Name', 'Hologram Cross-Sections', 'Color', 'w', 'Position', [150 150 1000 400]);
subplot(1,2,1);
plot(x_um, profile_holo1, 'b', 'LineWidth', 1.2); grid on;
xlabel('$x$ ($\mu$m)'); ylabel('Intensity (a.u.)');
title(sprintf('Hologram 1 Cross-section ($z_1 = %.1f$ $\\mu$m)', z1*1e6));
subplot(1,2,2);
plot(x_um, profile_holo2, 'r', 'LineWidth', 1.2); grid on;
xlabel('$x$ ($\mu$m)'); ylabel('Intensity (a.u.)');
title(sprintf('Hologram 2 Cross-section ($z_2 = %.1f$ $\\mu$m)', z2*1e6));

%% 8. PHASE RETRIEVAL (Plane-to-Plane)
Iterations = 300;
prop_dz = Propagator_BLASM(N, lambda, area, dx, dz);
prop_minus_dz = Propagator_BLASM(N, lambda, area, dx, -dz);
prop_minus_z1 = Propagator_BLASM(N, lambda, area, dx, -z1);

u1 = amp1_measured .* exp(1i * zeros(N, N));

fprintf('Đang chạy vòng lặp phục hồi pha (Phase Retrieval) %d iterations...\n', Iterations);
figure('Name', 'Reconstruction Progress', 'Position', [200, 200, 1000, 450], 'Color', 'w');

for iter = 1:Iterations
    u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);
    u2 = amp2_measured .* exp(1i * angle(u2_sim)); 
    
    u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz); 
    u1 = amp1_measured .* exp(1i * angle(u1_sim));
    
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

%% 9. CROSS-SECTION EVALUATION
amp_rec = abs(u_obj_rec);
fig = figure('Color','w', 'Name','Select Cross-section', 'Position', [300, 300, 800, 600]);
ax = axes(fig);
imagesc(ax, x_um, y_um, amp_rec); axis(ax,'image'); colormap(ax,gray); colorbar(ax);
xlabel(ax,'$x$ ($\mu$m)'); ylabel(ax,'$y$ ($\mu$m)');
title(ax,'Draw a line between TWO points to evaluate cross-section');
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
profile_true = interp2(X_grid, Y_grid, am, x_line, y_line, 'linear'); 
distance = hypot(x_line - x_line(1), y_line - y_line(1));

figure('Color','w');
plot(distance, profile_true, 'Color', '#0072BD', 'LineWidth', 1.5); hold on;
plot(distance, profile_rec, 'Color', '#D95319', 'LineWidth', 1.5, 'LineStyle', '--');
grid on;
xlabel('Distance ($\mu$m)'); ylabel('Amplitude (a.u.)');
legend('Ground Truth', 'Reconstruction');
title('Amplitude Cross-section Evaluation');

%% =========================================================================
%                               FUNCTIONS
% =========================================================================
% SỬA LỖI: Thay thế Propagator thường bằng Band-Limited ASM để chống nhiễu Aliasing
function p = Propagator_BLASM(N, lambda, area, dx, z)
    [X, Y] = meshgrid(1:N, 1:N);
    alpha = lambda * (Y - N/2 - 1) / area;
    beta  = lambda * (X - N/2 - 1) / area;
    
    % Tính giới hạn tần số để cắt nhiễu aliasing (Band-Limited ASM)
    limit = 1 / sqrt((2 * abs(z) / (N * dx))^2 + 1);
    mask = (alpha.^2 + beta.^2) <= limit^2;
    
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