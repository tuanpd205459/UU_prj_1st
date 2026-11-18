%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE-WAVELENGTH LASER INLINE DIGITAL HOLOGRAPHY (CHUẨN ĐƠN VỊ SI: m)
% MULTI-PLANE PHASE RETRIEVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. SYSTEM PARAMETERS (Đơn vị chuẩn: mét - m)

% Bước sóng đơn sắc của Laser
lambda = 405e-9;           % 405 nm

% Cảm biến (Sensor)
N    = 512;
dx   = 2e-6;             
area = dx * N;            

% Khoảng cách truyền sóng
z1 = 1.5e-3;               % 1.5 mm
z2 = 2.0e-3;               % 2.0 mm
dz = z2 - z1;            

%% 2. OBJECT (SINGLE NARROW SLIT)

slit_width  = 10e-6;       % Độ rộng khe: 10 um
slit_length = 0.5e-3;      % Chiều dài khe: 0.5 mm

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

% Hàm truyền qua
am = 1 - 0.2 * object;    % Nền = 1.0, Vật = 0.8
ph = 1.0 * object;        % Nền = 0 rad, Vật = 1.0 rad
t  = am .* exp(1i * ph);  

figure('Name', 'Single Narrow Slit');
subplot(1, 2, 1); imagesc(object); axis image; colormap gray; colorbar; title('Single Slit');
subplot(1, 2, 2); imagesc(am);     axis image; colormap gray; colorbar; title('Object Amplitude');

%% 3. HOLOGRAM RECORDING (MONOCHROMATIC LASER)

disp('Calculating Laser holograms...');

prop_z1 = Propagator(N, lambda, area, z1);
prop_z2 = Propagator(N, lambda, area, z2);

field_z1 = IFT2Dc(FT2Dc(t) .* prop_z1);
field_z2 = IFT2Dc(FT2Dc(t) .* prop_z2);

hologram1 = abs(field_z1).^2;
hologram2 = abs(field_z2).^2;

amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

figure('Name', 'Laser Holograms');
subplot(1, 2, 1); imagesc(hologram1); axis image; colormap gray; colorbar; title('Laser Hologram (z1 = 1.5 mm)');
subplot(1, 2, 2); imagesc(hologram2); axis image; colormap gray; colorbar; title('Laser Hologram (z2 = 2.0 mm)');

%% 4. MULTI-PLANE PHASE RETRIEVAL

Iterations = 1000;

prop_dz       = Propagator(N, lambda, area, dz);
prop_minus_dz = Propagator(N, lambda, area, -dz);
prop_minus_z1 = Propagator(N, lambda, area, -z1);

u1 = amp1_measured;

figure('Name', 'Laser Reconstruction');
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

%% 5. HELPER FUNCTIONS

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