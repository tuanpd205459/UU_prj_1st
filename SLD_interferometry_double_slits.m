%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLD INLINE DIGITAL HOLOGRAPHY 
% SPECTRAL SUPERPOSITION + MULTI-PLANE PHASE RETRIEVAL
% OBJECT: DOUBLE SLITS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. SYSTEM PARAMETERS (Đơn vị chuẩn: mét - m)

% Bước sóng trung tâm và băng thông 3 dB (FWHM)
lambda_c  = 405e-9;        % 405 nm
bandwidth = 4e-9;          % 4 nm FWHM

num_lambda = 100;

% Lấy mẫu bước sóng
lambda = linspace(lambda_c - bandwidth/2, ...
                  lambda_c + bandwidth/2, num_lambda);

% Phân bố phổ Gaussian & chuẩn hóa trọng số
S       = exp(-4 * log(2) * ((lambda - lambda_c) / bandwidth).^2);
weights = S / sum(S);

% Cảm biến (Sensor)
N    = 512;
dx   = 2e-6;             
area = dx * N;            

% Khoảng cách truyền sóng
z1 = 1.5e-3;               % 1.5 mm
z2 = 2.0e-3;               % 2.0 mm
dz = z2 - z1;            

%% 2. OBJECT (DOUBLE SLITS - 2 THANH SONG SONG)

slit_width  = 10e-6;       % Độ rộng mỗi thanh: 10 um
slit_length = 0.5e-3;      % Chiều dài mỗi thanh: 0.5 mm
slit_gap    = 300e-6;       % Khoảng cách mép trong giữa 2 thanh: 30 um

slit_width_px  = max(1, round(slit_width / dx));
slit_length_px = max(1, round(slit_length / dx));
gap_px         = round(slit_gap / dx);

object = zeros(N, N);
cx = round(N/2);
cy = round(N/2);

% Tọa độ theo trục Y (chiều dài thanh)
y1 = cy - floor(slit_length_px/2);
y2 = y1 + slit_length_px - 1;

% Tọa độ thanh bên trái
x1_left = cx - floor(gap_px/2) - slit_width_px;
x2_left = x1_left + slit_width_px - 1;

% Tọa độ thanh bên phải
x1_right = cx + ceil(gap_px/2);
x2_right = x1_right + slit_width_px - 1;

% Gán giá trị = 1 cho cả 2 thanh
object(y1:y2, x1_left:x2_left)   = 1;
object(y1:y2, x1_right:x2_right) = 1;

% Hàm truyền qua
am = 1 - 0.2 * object;    % Nền = 1.0, Vật = 0.8
ph = 1.0 * object;        % Nền = 0 rad, Vật = 1.0 rad
t  = am .* exp(1i * ph);  

figure('Name', 'Double Narrow Slits');
subplot(1, 2, 1); imagesc(object); axis image; colormap gray; colorbar; title('Double Slits');
subplot(1, 2, 2); imagesc(am);     axis image; colormap gray; colorbar; title('Object Amplitude');

%% 3. SPECTRAL SUPERPOSITION

I1_total = zeros(N, N);
I2_total = zeros(N, N);

disp('Calculating SLD holograms...');

for j = 1:num_lambda
    lambda_j = lambda(j);

    prop_z1 = Propagator(N, lambda_j, area, z1);
    prop_z2 = Propagator(N, lambda_j, area, z2);

    field_z1 = IFT2Dc(FT2Dc(t) .* prop_z1);
    field_z2 = IFT2Dc(FT2Dc(t) .* prop_z2);

    I1 = abs(field_z1).^2;
    I2 = abs(field_z2).^2;

    I1_total = I1_total + weights(j) * I1;
    I2_total = I2_total + weights(j) * I2;
end

hologram1 = I1_total;
hologram2 = I2_total;

amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

figure('Name', 'SLD Holograms');
subplot(1, 2, 1); imagesc(hologram1); axis image; colormap gray; colorbar; title('SLD Hologram (z1 = 1.5 mm)');
subplot(1, 2, 2); imagesc(hologram2); axis image; colormap gray; colorbar; title('SLD Hologram (z2 = 2.0 mm)');

%% 4. MULTI-PLANE PHASE RETRIEVAL

% Dùng bước sóng trung tâm lambda_c làm mô hình tái tạo xấp xỉ
Iterations = 300;

prop_dz       = Propagator(N, lambda_c, area, dz);
prop_minus_dz = Propagator(N, lambda_c, area, -dz);
prop_minus_z1 = Propagator(N, lambda_c, area, -z1);

u1 = amp1_measured;

figure('Name', 'SLD Reconstruction');
disp('Starting phase retrieval...');

for iter = 1:Iterations
    % A. Truyền xuôi z1 -> z2
    u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);

    % B. Ép biên độ đo được tại z2 (giữ lại pha ước lượng)
    u2 = amp2_measured .* exp(1i * angle(u2_sim));

    % C. Truyền ngược z2 -> z1
    u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz);

    % D. Ép biên độ đo được tại z1 (giữ lại pha ước lượng)
    u1 = amp1_measured .* exp(1i * angle(u1_sim));

    % Hiển thị quá trình tái tạo
    if mod(iter, 10) == 0 || iter == 1 || iter == Iterations
        u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);

        subplot(1, 2, 1);
        imagesc(abs(u_obj_rec)); axis image; colormap gray; colorbar;
        title(sprintf('Amplitude - Iter %d', iter));

        subplot(1, 2, 2);
        imagesc(angle(u_obj_rec)); axis image; colormap gray; colorbar;
        title(sprintf('Phase - Iter %d', iter));

        drawnow;
    end
end

disp('Reconstruction completed!');

%% 5. SPECTRAL DISTRIBUTION

% Mở rộng dải lấy mẫu phổ lên +- 2.5 * FWHM (tổng cộng 5 * bandwidth)
span   = 1 * bandwidth; 
lambda = linspace(lambda_c - span/2, lambda_c + span/2, num_lambda);

% Tính lại phổ S và trọng số
S       = exp(-4 * log(2) * ((lambda - lambda_c) / bandwidth).^2);
weights = S / sum(S);

% Vẽ đồ thị phổ với vùng quan sát rộng
figure('Name', 'SLD Spectrum');
plot(lambda * 1e9, S, 'o-', 'LineWidth', 1.5);
xlabel('Wavelength (nm)');
ylabel('Relative Spectral Intensity');
title('Gaussian SLD Spectrum');
xlim([(lambda_c - span/2)*1e9, (lambda_c + span/2)*1e9]);
% ylim([-0.05, 1.05]);
grid on;

%% 6. HELPER FUNCTIONS

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