clc; clear; close all;
%% 1. SYSTEM PARAMETERS
N = 1024;              % (Sửa) Cố định N = 1024 giống code 2
lambda = 500e-9;       % (Sửa) Bước sóng 500nm giống code 2
area = 2000e-6;        % FOV 0.002m
dx = area / N;         % Tính lại dx theo N
dy = dx;

z1 = 4000e-6;          % Plane 1 (m) tương đương z = 0.004 của code 2
z2 = 4500e-6;          % Plane 2 (m)
dz = z2 - z1;
z_max = N * dx^2 / lambda;

% Define Coordinates
x = linspace(-area/2, area/2, N);
y = linspace(-area/2, area/2, N);
x_um = x * 1e6;
y_um = y * 1e6;

fprintf('====================================\n');
fprintf('       SYSTEM PARAMETERS (UPDATED)\n');
fprintf('====================================\n');
fprintf('Wavelength  : %.1f nm\n', lambda*1e9);
fprintf('Grid size   : %d x %d pixels\n', N, N);
fprintf('Pixel pitch : %.4f um\n', dx*1e6);
fprintf('FOV         : %.2f um\n', area*1e6);
fprintf('z1          : %.2f um\n', z1*1e6);

%% 3. CREATE OBJECT (Two Point Objects)
% Uniform background
object = ones(N,N);
% Two opaque points (Sửa bỏ hàm round vì N=1024 là số chẵn)
object(N/2+1, N/2+1) = 0;
object(N/2+7, N/2+1) = 0;

%% 4. TRANSMISSION FUNCTION
% (Sửa) Bỏ dấu chuyển vị ', dùng trực tiếp object giống Code 2
t = object; 
am = abs(t); 

%% 5. DISPLAY OBJECT
figure('Name', 'Object Details', 'Color', 'w', 'Position', [100 200, 500, 450]);
imagesc(x_um, y_um, object); axis image; colormap gray; colorbar;
title('Binary Object / Transmission Function'); 
xlabel('x (\mum)'); ylabel('y (\mum)');

%% 6. ORIGINAL OBJECT CROSS-SECTION
mid_y = N/2; 
profile_obj_am = am(:, mid_y); % Lấy mặt cắt theo trục y (cột giữa) vì 2 điểm xếp dọc
figure('Name', 'Ground Truth Cross-Section', 'Color', 'w', 'Position', [150 500, 600, 400]);
plot(y_um, profile_obj_am, 'k', 'LineWidth', 1.5);
grid on; grid minor; 
xlabel('$y$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('Amplitude');
title('Object Amplitude Cross-section at center ($x=0$)', 'Interpreter', 'latex');
xlim([-50, 50]); 

%% 7. HOLOGRAM SIMULATION
% (Sửa) Dùng hàm Propagator tiêu chuẩn (ASM) không có lọc băng thông giống Code 2
prop_z1 = Propagator(N, lambda, area, z1);
prop_z2 = Propagator(N, lambda, area, z2);

field_z1_true = IFT2Dc(FT2Dc(t) .* prop_z1);
field_z2_true = IFT2Dc(FT2Dc(t) .* prop_z2);

hologram1 = abs(field_z1_true).^2;
hologram2 = abs(field_z2_true).^2;
amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);
figure;
imshow(flipud(rot90(hologram1)), []);
title('simulated hologram / a.u.')
xlabel({'x / px'})
ylabel({'y / px'})
axis on
set(gca,'YDir','normal')
colormap('gray')
colorbar;
%% 8. VISUALIZATION (HOLOGRAMS)
figure('Name', 'Measured Holograms', 'Color', 'w', 'Position', [100 100 1000 400]);
subplot(1,2,1); 
imagesc(x_um, y_um, hologram1); axis image;
title(sprintf('Hologram 1 ($z_1 = %d$ $\\mu$m)', round(z1*1e6)), 'Interpreter', 'latex'); 
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex'); colorbar; colormap gray;
subplot(1,2,2); 
imagesc(x_um, y_um, hologram2); axis image;
title(sprintf('Hologram 2 ($z_2 = %d$ $\\mu$m)', round(z2*1e6)), 'Interpreter', 'latex'); 
xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex'); colorbar; colormap gray;

%% 10. PHASE RETRIEVAL (Plane-to-Plane)
Iterations = 300;
prop_dz = Propagator(N, lambda, area, dz);
prop_minus_dz = Propagator(N, lambda, area, -dz);
prop_minus_z1 = Propagator(N, lambda, area, -z1);
u1 = amp1_measured .* exp(1i * zeros(N, N));

fprintf('\nRunning Phase Retrieval loop (%d iterations)...\n', Iterations);
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
        xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex');
        
        subplot(1,2,2); imagesc(x_um, y_um, angle(u_obj_rec)); axis image; colormap gray;
        title(sprintf('Rec. Phase (Iter: %d)', iter)); colorbar;
        xlabel('$x$ ($\mu$m)', 'Interpreter', 'latex'); ylabel('$y$ ($\mu$m)', 'Interpreter', 'latex');
        drawnow;
    end
end
fprintf('Reconstruction complete!\n');

%% =========================================================================
%                               FUNCTIONS
% =========================================================================
% (Sửa) Hàm Propagator tiêu chuẩn (Standard ASM) giống lõi của Code 2
function [p] = Propagator(N, lambda, area, z)
p = zeros(N,N);
for ii = 1:N; 
    for jj = 1:N
        alpha = lambda*(ii - N/2 -1)/area;
        beta = lambda*(jj - N/2 -1)/area;
        if ((alpha^2 + beta^2) <= 1)
        p(ii,jj) = exp(2*pi*i*z*sqrt(1 - alpha^2 - beta^2)/lambda);
        end; % if
    end
end;
end


function [out] = IFT2Dc(in)

[Nx Ny] = size(in);

f1 = zeros(Nx,Ny);

for ii = 1:Nx
    for jj = 1:Ny
        f1(ii, jj) = exp(-i*pi*(ii + jj));
    end
end

FT = ifft2(f1.*in);

out = f1.*FT;
end


function [out] = FT2Dc(in)

[Nx Ny] = size(in);

f1 = zeros(Nx,Ny);

for ii = 1:Nx
    for jj = 1:Ny
         f1(ii,jj) = exp(i*pi*(ii + jj));
    end
end

FT = fft2(f1.*in);

out = f1.*FT;
end