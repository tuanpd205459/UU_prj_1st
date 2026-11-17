%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLD INLINE DIGITAL HOLOGRAPHY - COHERENCE VERIFICATION (OPTIMIZED)
% So sánh trực tiếp SLD (Low-Coherence) và Laser (Monochromatic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. SYSTEM PARAMETERS (mm)

lambda_c   = 0.405e-3;     % Bước sóng trung tâm: 405 nm
bandwidth  = 0.030e-3;     % Bề rộng phổ FWHM: 30 nm
num_lambda = 60;           % Số điểm lấy mẫu phổ

% Lấy mẫu phổ trong khoảng +- 1.5 * FWHM (phủ >99.9% năng lượng phổ Gaussian)
span   = 3 * bandwidth;
lambda = linspace(lambda_c - span/2, lambda_c + span/2, num_lambda);

% Phân bố phổ Gaussian & chuẩn hóa trọng số
S       = exp(-4 * log(2) * ((lambda - lambda_c) / bandwidth).^2);
weights = S / sum(S);

% Tính chiều dài kết hợp lý thuyết Lc (mm)
Lc = (2 * log(2) / pi) * (lambda_c^2 / bandwidth);
fprintf('Chieu dai ket hop ly thuyet Lc = %.2f um\n', Lc * 1e3);

% Cảm biến (Sensor) - Tăng độ phân giải lấy mẫu để tránh nhiễu aliasing ở rìa
N    = 800;
dx   = 1e-3;               % Kích thước pixel: 1 um
area = dx * N;             % Kích thước vùng mô phỏng: 0.8 x 0.8 mm

% Khoảng cách truyền sóng
z1 = 2.0;                  % 2.0 mm

%% 2. OBJECT (IN-LINE HOLOGRAPHY SETUP)

slit_width  = 10e-3;       % Độ rộng thanh vật mẫu: 10 um
slit_length = 0.4;         % Chiều dài thanh vật mẫu: 0.4 mm

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

% Chọn loại vật mẫu in-line (nền xung quanh = 1 đóng vai trò sóng tham chiếu):
% Cách 1: Vật hấp thụ một phần (bán trong suốt, biên độ truyền qua = 0.2)
t = 1 - 0.8 * object;

% Cách 2: Vật thuần pha (bỏ dấu % ở dòng dưới nếu muốn thử nghiệm vật pha)
% t = exp(1i * 1.5 * object);

%% 3. MÔ PHỎNG HOLOGRAM: SLD VS. LASER ĐƠN SẮC

% A. Trường hợp chuẩn đối chiếu: Laser đơn sắc tại lambda_c
prop_mono = Propagator(N, lambda_c, area, z1);
I_mono    = abs(IFT2Dc(FT2Dc(t) .* prop_mono)).^2;

% B. Trường hợp nguồn kết hợp thấp (SLD): Cộng dồn cường độ không kết hợp theo phổ
I_sld = zeros(N, N);
for j = 1:num_lambda
    prop_j  = Propagator(N, lambda(j), area, z1);
    field_j = IFT2Dc(FT2Dc(t) .* prop_j);
    I_sld   = I_sld + weights(j) * (abs(field_j).^2);
end

% Kiểm tra bảo toàn năng lượng
fprintf('Ty so bao toan nang luong (SLD / Object) = %.5f\n', sum(I_sld(:)) / sum(abs(t(:)).^2));

%% 4. HIỂN THỊ & ĐỐI CHIẾU KIỂM TRA (VISUAL VERIFICATION)

x_axis = ((1:N) - cx) * dx * 1e3; % Trục tọa độ ngang (đơn vị: um)
row_center = cy;

figure('Name', 'Kiem tra mo phong In-line Hologram (SLD vs Laser)', 'Position', [100, 100, 1150, 680]);

% Ảnh 2D của Laser
subplot(2, 2, 1);
imagesc(x_axis, x_axis, I_mono);
axis image; colormap gray; colorbar;
caxis([0.4 1.6]); % Cố định thang màu để so sánh độ tương phản vân công bằng
title('1. In-line Hologram - Laser đơn sắc (\Delta\lambda = 0)');
xlabel('x (\mum)'); ylabel('y (\mum)');

% Ảnh 2D của SLD
subplot(2, 2, 2);
imagesc(x_axis, x_axis, I_sld);
axis image; colormap gray; colorbar;
caxis([0.4 1.6]);
title(sprintf('2. In-line Hologram - Nguồn SLD (\\Delta\\lambda = %.0f nm)', bandwidth*1e6));
xlabel('x (\mum)'); ylabel('y (\mum)');

% Lát cắt ngang qua tâm vật thể
subplot(2, 1, 2);
plot(x_axis, I_mono(row_center, :), 'b--', 'LineWidth', 1.1); hold on;
plot(x_axis, I_sld(row_center, :),  'r-',  'LineWidth', 1.6);
yline(1.0, 'k:', 'LineWidth', 1.0); % Đường mức nền tham chiếu I = 1
grid on;
xlim([-350 350]);
xlabel('Vị trí ngang x (\mum)');
ylabel('Cường độ chuẩn hóa');
title('3. Lát cắt ngang qua tâm vật: Vân SLD (đỏ) tắt hẳn về mức nền 1.0 ở rìa xa');
legend('Laser đơn sắc (Độ kết hợp cao)', 'Nguồn SLD (Độ kết hợp thấp)', 'Cường độ nền tham chiếu (I_R = 1)', 'Location', 'northeast');

%% 5. HELPER FUNCTIONS (BAND-LIMITED ASM)

function p = Propagator(N, lambda, area, z)
    [X, Y] = meshgrid(1:N, 1:N);
    fx = (X - N/2 - 1) / area;
    fy = (Y - N/2 - 1) / area;
    
    alpha = lambda * fx;
    beta  = lambda * fy;
    
    % Giới hạn băng tần (Band-Limited ASM) tránh nhiễu aliasing ở biên
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