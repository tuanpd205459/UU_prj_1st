%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLD INLINE DIGITAL HOLOGRAPHY - ADJUSTABLE COHERENCE LENGTH (Lc)
% Nhập trực tiếp Lc (um) để tự động tính bề rộng phổ và mô phỏng
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. SYSTEM PARAMETERS (mm)

lambda_c = 0.405e-3;       % Bước sóng trung tâm: 405 nm (0.405e-3 mm)

% === THAY ĐỔI CHIỀU DÀI KẾT HỢP TẠI ĐÂY ===
Lc_target_um = 30;        % Chiều dài kết hợp mong muốn Lc (đơn vị: um)
                           % Gợi ý thử nghiệm:
                           %   Lc = 1.5 um -> SLD siêu băng rộng (vân tắt cực nhanh)
                           %   Lc = 5.0 um -> SLD trung bình
                           %   Lc = 25  um -> Gần giống Laser (vân lan rộng)
% ==========================================

Lc = Lc_target_um * 1e-3;  % Đổi sang đơn vị mm

% Tính ngược bề rộng phổ FWHM (bandwidth) từ Lc cho phổ Gaussian
bandwidth = (2 * log(2) / pi) * (lambda_c^2 / Lc);

% Tự động điều chỉnh số điểm lấy mẫu phổ để tránh aliasing khi Lc nhỏ (phổ rộng)
num_lambda = max(60, ceil(bandwidth / 0.5e-6)); 

fprintf('Chieu dai ket hop Lc  = %.2f um\n', Lc_target_um);
fprintf('Be rong pho FWHM      = %.2f nm\n', bandwidth * 1e6);
fprintf('So diem lay mau pho   = %d\n', num_lambda);

% Lấy mẫu phổ trong khoảng +- 1.5 * FWHM (phủ >99.9% năng lượng phổ Gaussian)
span   = 3 * bandwidth;
lambda = linspace(lambda_c - span/2, lambda_c + span/2, num_lambda);

% Kiểm tra an toàn không để bước sóng <= 0 nếu Lc đặt quá nhỏ
if any(lambda <= 0)
    error('Lc qua nho khien dai pho vuot qua buoc song trung tam! Hay tang Lc_target_um.');
end

% Phân bố phổ Gaussian & chuẩn hóa trọng số
S       = exp(-4 * log(2) * ((lambda - lambda_c) / bandwidth).^2);
weights = S / sum(S);

% Cảm biến (Sensor)
N    = 800;
dx   = 1e-3;               % Kích thước pixel: 1 um
area = dx * N;             % Vùng mô phỏng: 0.8 x 0.8 mm

% Khoảng cách truyền sóng
z1 = 2.0;                  % 2.0 mm

%% 2. OBJECT (IN-LINE HOLOGRAPHY SETUP)

slit_width  = 10e-3;       % Độ rộng vật mẫu: 10 um
slit_length = 0.4;         % Chiều dài vật mẫu: 0.4 mm

slit_width_px  = max(1, round(slit_width / dx));
slit_length_px = max(1, round(slit_length / dx));

object = zeros(N, N);
cx = round(N/2);
cy = round(N/2);

x1 = cx - floor(slit_width_px/2);
x2 = x1 + slit_width_px - 1;
y1 = cy - floor(slit_length_px/2);
y2 = y1 + slit_length_px - 1;

object(y1:y2, x1:x2) = 1;

% Vật hấp thụ một phần trên nền sáng (In-line Hologram)
t = 1 - 0.8 * object;

%% 3. MÔ PHỎNG HOLOGRAM: SLD VS. LASER ĐƠN SẮC

% A. Laser đơn sắc (Lc = vô cùng)
prop_mono = Propagator(N, lambda_c, area, z1);
I_mono    = abs(IFT2Dc(FT2Dc(t) .* prop_mono)).^2;

% B. Nguồn SLD với chiều dài kết hợp Lc_target_um
I_sld = zeros(N, N);
for j = 1:num_lambda
    prop_j  = Propagator(N, lambda(j), area, z1);
    field_j = IFT2Dc(FT2Dc(t) .* prop_j);
    I_sld   = I_sld + weights(j) * (abs(field_j).^2);
end

%% 4. HIỂN THỊ & ĐỐI CHIẾU

x_axis = ((1:N) - cx) * dx * 1e3; % Trục tọa độ ngang (um)
row_center = cy;

figure('Name', sprintf('SLD Hologram (Lc = %.1f um)', Lc_target_um), 'Position', [100, 100, 1150, 680]);

% Ảnh 2D Laser
subplot(2, 2, 1);
imagesc(x_axis, x_axis, I_mono);
axis image; colormap gray; colorbar;
caxis([0.4 1.6]);
title('1. Laser đơn sắc (L_c = \infty)');
xlabel('x (\mum)'); ylabel('y (\mum)');

% Ảnh 2D SLD
subplot(2, 2, 2);
imagesc(x_axis, x_axis, I_sld);
axis image; colormap gray; colorbar;
caxis([0.4 1.6]);
title(sprintf('2. Nguồn SLD (L_c = %.1f \\mum, \\Delta\\lambda = %.1f nm)', Lc_target_um, bandwidth*1e6));
xlabel('x (\mum)'); ylabel('y (\mum)');

% Lát cắt ngang qua tâm
subplot(2, 1, 2);
plot(x_axis, I_mono(row_center, :), 'b--', 'LineWidth', 1.1); hold on;
plot(x_axis, I_sld(row_center, :),  'r-',  'LineWidth', 1.6);
yline(1.0, 'k:', 'LineWidth', 1.0);
grid on;
xlim([-350 350]);
xlabel('Vị trí ngang x (\mum)');
ylabel('Cường độ chuẩn hóa');
title(sprintf('3. Lát cắt ngang qua tâm với L_c = %.1f \\mum', Lc_target_um));
legend('Laser đơn sắc (L_c = \infty)', sprintf('Nguồn SLD (L_c = %.1f \\mum)', Lc_target_um), 'Nền tham chiếu (I_R = 1)', 'Location', 'northeast');

%% 5. HELPER FUNCTIONS (BAND-LIMITED ASM)

function p = Propagator(N, lambda, area, z)
    [X, Y] = meshgrid(1:N, 1:N);
    fx = (X - N/2 - 1) / area;
    fy = (Y - N/2 - 1) / area;
    
    alpha = lambda * fx;
    beta  = lambda * fy;
    
    f_limit = 1 / (lambda * sqrt(1 + (2 * z / area)^2));
    mask = (alpha.^2 + beta.^2 <= 1) & (abs(fx) <= f_limit) & (abs(fy) <= f_limit);
    
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