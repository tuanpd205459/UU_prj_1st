%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SO SÁNH INLINE DIGITAL HOLOGRAPHY: SINGLE-WAVELENGTH LASER VS. SLD
% OBJECT: DOUBLE SLITS (2 THANH SONG SONG) - CHUẨN ĐƠN VỊ SI (m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. SYSTEM PARAMETERS (Đơn vị chuẩn: mét - m)

% Bước sóng trung tâm & thông số phổ SLD
lambda_c   = 405e-9;       % 405 nm (Bước sóng của Laser & tâm phổ SLD)
bandwidth  = 40e-9;         % 4 nm FWHM (Băng thông SLD)
num_lambda = 100;

% Lấy mẫu bước sóng (Mở rộng dải lấy mẫu phổ +- 1.5*FWHM để bao phủ đủ đuôi Gaussian)
span    = 3 * bandwidth;
lambda  = linspace(lambda_c - span/2, lambda_c + span/2, num_lambda);
S       = exp(-4 * log(2) * ((lambda - lambda_c) / bandwidth).^2);
weights = S / sum(S);

% Cảm biến (Sensor)
N    = 512;
dx   = 2e-6;               % Kích thước pixel: 2 um
area = dx * N;             % Kích thước vùng cảm biến: 1.024 mm
x_axis_um = ((1:N) - N/2 - 1) * dx * 1e6; % Trục tọa độ X (đơn vị: um)

% Khoảng cách truyền sóng
z1 = 1.5e-3;               % 1.5 mm
z2 = 2.0e-3;               % 2.0 mm
dz = z2 - z1;              % 0.5 mm

%% 2. OBJECT: 2 THANH CÁCH NHAU (DOUBLE SLITS)

slit_width  = 10e-6;       % Độ rộng mỗi thanh: 10 um
slit_length = 0.5e-3;      % Chiều dài mỗi thanh: 0.5 mm
slit_gap    = 30e-6;       % Khoảng cách mép trong giữa 2 thanh: 30 um

slit_w_px   = max(1, round(slit_width / dx));
slit_l_px   = max(1, round(slit_length / dx));
gap_px      = round(slit_gap / dx);

object = zeros(N, N);
cx = round(N/2);
cy = round(N/2);

% Tọa độ theo trục Y (chiều dài)
y1 = cy - floor(slit_l_px/2);
y2 = y1 + slit_l_px - 1;

% Thanh bên trái (Left bar)
x1_left = cx - floor(gap_px/2) - slit_w_px;
x2_left = x1_left + slit_w_px - 1;

% Thanh bên phải (Right bar)
x1_right = cx + ceil(gap_px/2);
x2_right = x1_right + slit_w_px - 1;

object(y1:y2, x1_left:x2_left)   = 1;
object(y1:y2, x1_right:x2_right) = 1;

% Hàm truyền qua phức của vật
am = 1 - 0.2 * object;     % Nền = 1.0, Vật = 0.8
ph = 1.0 * object;         % Nền = 0 rad, Vật = 1.0 rad
t  = am .* exp(1i * ph);  

figure('Name', 'Double Slits Object', 'Position', [50, 500, 800, 350]);
subplot(1, 2, 1); 
imagesc(x_axis_um, x_axis_um, object); axis image; colormap gray; colorbar;
title('Object Mask (2 Thanh Song Song)'); xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(1, 2, 2); 
imagesc(x_axis_um, x_axis_um, am); axis image; colormap gray; colorbar;
title('Object Amplitude'); xlabel('x (\mum)'); ylabel('y (\mum)');

%% 3A. TÍNH HOLOGRAM CHO SINGLE-WAVELENGTH LASER

disp('Calculating Single-Wavelength Laser holograms...');

prop_laser_z1 = Propagator(N, lambda_c, area, z1);
prop_laser_z2 = Propagator(N, lambda_c, area, z2);

field_laser_z1 = IFT2Dc(FT2Dc(t) .* prop_laser_z1);
field_laser_z2 = IFT2Dc(FT2Dc(t) .* prop_laser_z2);

holo_laser_1 = abs(field_laser_z1).^2;
holo_laser_2 = abs(field_laser_z2).^2;

%% 3B. TÍNH HOLOGRAM CHO SLD (SPECTRAL SUPERPOSITION)

disp('Calculating SLD holograms...');

I1_total = zeros(N, N);
I2_total = zeros(N, N);

for j = 1:num_lambda
    lambda_j = lambda(j);

    prop_z1 = Propagator(N, lambda_j, area, z1);
    prop_z2 = Propagator(N, lambda_j, area, z2);

    field_z1 = IFT2Dc(FT2Dc(t) .* prop_z1);
    field_z2 = IFT2Dc(FT2Dc(t) .* prop_z2);

    I1_total = I1_total + weights(j) * (abs(field_z1).^2);
    I2_total = I2_total + weights(j) * (abs(field_z2).^2);
end

holo_sld_1 = I1_total;
holo_sld_2 = I2_total;

%% 4. SO SÁNH ẢNH HOLOGRAM 2D (LASER VS. SLD)

figure('Name', 'Hologram Comparison: Laser vs SLD', 'Position', [100, 100, 900, 750]);

% Hàng 1: Tại mặt phẳng z1 = 1.5 mm
subplot(2, 2, 1);
imagesc(x_axis_um, x_axis_um, holo_laser_1); axis image; colormap gray; colorbar;
title(sprintf('Laser Hologram (z1 = %.1f mm)', z1*1e3));
xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2, 2, 2);
imagesc(x_axis_um, x_axis_um, holo_sld_1); axis image; colormap gray; colorbar;
title(sprintf('SLD Hologram (z1 = %.1f mm)', z1*1e3));
xlabel('x (\mum)'); ylabel('y (\mum)');

% Hàng 2: Tại mặt phẳng z2 = 2.0 mm
subplot(2, 2, 3);
imagesc(x_axis_um, x_axis_um, holo_laser_2); axis image; colormap gray; colorbar;
title(sprintf('Laser Hologram (z2 = %.1f mm)', z2*1e3));
xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2, 2, 4);
imagesc(x_axis_um, x_axis_um, holo_sld_2); axis image; colormap gray; colorbar;
title(sprintf('SLD Hologram (z2 = %.1f mm)', z2*1e3));
xlabel('x (\mum)'); ylabel('y (\mum)');

%% 5. SO SÁNH MẶT CẮT NGANG HOLOGRAM (1D CROSS-SECTION TẠI Y = 0)

profile_laser_1 = holo_laser_1(cy, :);
profile_sld_1   = holo_sld_1(cy, :);

profile_laser_2 = holo_laser_2(cy, :);
profile_sld_2   = holo_sld_2(cy, :);

figure('Name', 'Hologram Cross-Section Comparison', 'Position', [150, 150, 1000, 600]);

% Mặt cắt ngang tại z1
subplot(2, 1, 1);
plot(x_axis_um, profile_laser_1, 'b-', 'LineWidth', 1.2); hold on;
plot(x_axis_um, profile_sld_1, 'r--', 'LineWidth', 1.5);
title(sprintf('Mặt cắt ngang Hologram tại z1 = %.1f mm (y = 0)', z1*1e3));
xlabel('Vị trí ngang x (\mum)');
ylabel('Cường độ sáng (Intensity)');
legend('Single-Wavelength Laser (405 nm)', 'SLD (FWHM = 4 nm)', 'Location', 'best');
xlim([-350, 350]); grid on;

% Mặt cắt ngang tại z2
subplot(2, 1, 2);
plot(x_axis_um, profile_laser_2, 'b-', 'LineWidth', 1.2); hold on;
plot(x_axis_um, profile_sld_2, 'r--', 'LineWidth', 1.5);
title(sprintf('Mặt cắt ngang Hologram tại z2 = %.1f mm (y = 0)', z2*1e3));
xlabel('Vị trí ngang x (\mum)');
ylabel('Cường độ sáng (Intensity)');
legend('Single-Wavelength Laser (405 nm)', 'SLD (FWHM = 4 nm)', 'Location', 'best');
xlim([-350, 350]); grid on;

%% 6. MULTI-PLANE PHASE RETRIEVAL (SO SÁNH TÁI TẠO LASER VS. SLD)

Iterations = 300; % Chạy 300 vòng lặp để so sánh nhanh kết quả tái tạo

rec_laser = RunPhaseRetrieval(sqrt(holo_laser_1), sqrt(holo_laser_2), N, lambda_c, area, z1, dz, Iterations);
rec_sld   = RunPhaseRetrieval(sqrt(holo_sld_1),   sqrt(holo_sld_2),   N, lambda_c, area, z1, dz, Iterations);

figure('Name', 'Phase Retrieval Comparison: Laser vs SLD', 'Position', [200, 100, 950, 700]);

subplot(2, 2, 1);
imagesc(x_axis_um, x_axis_um, abs(rec_laser)); axis image; colormap gray; colorbar;
title('Laser - Reconstructed Amplitude'); xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2, 2, 2);
imagesc(x_axis_um, x_axis_um, angle(rec_laser)); axis image; colormap gray; colorbar;
title('Laser - Reconstructed Phase (rad)'); xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2, 2, 3);
imagesc(x_axis_um, x_axis_um, abs(rec_sld)); axis image; colormap gray; colorbar;
title('SLD - Reconstructed Amplitude'); xlabel('x (\mum)'); ylabel('y (\mum)');

subplot(2, 2, 4);
imagesc(x_axis_um, x_axis_um, angle(rec_sld)); axis image; colormap gray; colorbar;
title('SLD - Reconstructed Phase (rad)'); xlabel('x (\mum)'); ylabel('y (\mum)');

%% 7. HELPER FUNCTIONS

function u_obj_rec = RunPhaseRetrieval(amp1, amp2, N, lambda_c, area, z1, dz, Iterations)
    prop_dz       = Propagator(N, lambda_c, area, dz);
    prop_minus_dz = Propagator(N, lambda_c, area, -dz);
    prop_minus_z1 = Propagator(N, lambda_c, area, -z1);

    u1 = amp1;
    for iter = 1:Iterations
        u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);
        u2     = amp2 .* exp(1i * angle(u2_sim));
        u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz);
        u1     = amp1 .* exp(1i * angle(u1_sim));
    end
    u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);
end

function p = Propagator(N, lambda, area, z)
    [X, Y] = meshgrid(1:N, 1:N);
    alpha = lambda * (Y - N/2 - 1) / area;
    beta  = lambda * (X - N/2 - 1) / area;

    mask = (alpha.^2 + beta.^2) <= 1;
    p = zeros(N, N);
    p(mask) = exp(2 * pi * 1i * z .* sqrt(1 - alpha(mask).^2 - beta(mask).^2) / lambda);
end

function out = FT2Dc(in)
    [Nx, Ny] = size(in);
    f1 = ((-1).^(1:Nx).') * ((-1).^(1:Ny));
    out = f1 .* fft2(f1 .* in);
end

function out = IFT2Dc(in)
    [Nx, Ny] = size(in);
    f1 = ((-1).^(1:Nx).') * ((-1).^(1:Ny));
    out = f1 .* ifft2(f1 .* in);
end