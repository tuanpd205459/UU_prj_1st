
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. SIMULATION PARAMETERS

lambda = 633e-9;        % Wavelength (m)
z1 = 1500e-6;            % Propagation distance to Plane 1 (m)
z2 = 2000e-6;            % Propagation distance to Plane 2 (m)
dz = z2 - z1;           % Distance between two planes

object_mode = 'object'; % 'object' or '2dot'

%% 2. GENERATE OBJECT

switch lower(object_mode)

    case 'object'
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

    case '2dot'
        % Grid
        Nx = 512;
        Ny = 512;
        N = Nx;           % KHAI BÁO THÊM N
        
        dx = 20e-9;
        dy = dx;
        area = dx * Nx;   % KHAI BÁO THÊM area

        % Geometry
        diameter = 200e-9;
        radius = diameter/2;
        center_distance = 400e-9;

        % Coordinate grid
        x = ((0:Nx-1)-(Nx-1)/2)*dx;
        y = ((0:Ny-1)-(Ny-1)/2)*dy;

        [X,Y] = meshgrid(x,y);

        % Two dot positions
        x1 = -center_distance/2;
        x2 =  center_distance/2;

        % Circular masks
        circle1 = (X-x1).^2 + Y.^2 <= radius^2;
        circle2 = (X-x2).^2 + Y.^2 <= radius^2;

        object_mask = circle1 | circle2;

        % Optical properties
        A_object = 0.8;
        phase_object = pi/4;

        % Background transmission = 1
        T = ones(Ny,Nx);

        % Object transmission
        T(object_mask) = A_object*exp(1i*phase_object);

        % Ground truth
        amplitude = abs(T);
        ph = angle(T);

    otherwise
        error('Choose object_mode: object or 2dot.');
end

%% 3. COORDINATE GRID & INFO

x = ((0:Nx-1)-(Nx-1)/2)*dx;
y = ((0:Ny-1)-(Ny-1)/2)*dy;
[X,Y] = meshgrid(x,y);

FOV_x = Nx*dx;
FOV_y = Ny*dy;

fprintf('\n========== SIMULATION PARAMETERS ==========\n');
fprintf('Object mode    : %s\n',object_mode);
fprintf('Grid size      : %d x %d pixels\n',Nx,Ny);
fprintf('Pixel size     : %.2f nm\n',dx*1e9);
fprintf('FOV            : %.3f x %.3f um\n', FOV_x*1e6,FOV_y*1e6);
fprintf('Wavelength     : %.1f nm\n',lambda*1e9);
fprintf('Plane z1       : %.3f um\n',z1*1e6);
fprintf('Plane z2       : %.3f um\n',z2*1e6);
fprintf('Distance dz    : %.3f um\n',dz*1e6);

%% 4. MÔ PHỎNG 2 HOLOGRAM TẠI Z1 VÀ Z2
% Tạo bộ nhân tố truyền sóng (Propagator)
prop_z1 = Propagator(N, lambda, area, z1);
prop_z2 = Propagator(N, lambda, area, z2);

% Lan truyền sóng đến z1 và z2 (Sửa 't' thành 'T')
field_z1_true = IFT2Dc(FT2Dc(T) .* prop_z1);
field_z2_true = IFT2Dc(FT2Dc(T) .* prop_z2);

% Holograms (Cường độ sáng thu được từ cảm biến)
hologram1 = abs(field_z1_true).^2; 
hologram2 = abs(field_z2_true).^2;

% Biên độ đo được (CCD chỉ ghi nhận được cường độ)
amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

%% 5. VISUALIZATION (VẬT THỂ & HOLOGRAM)
x_um = x*1e6;
y_um = y*1e6;

figure('Name','Object and Measured Holograms','Color','w', 'Position',[100 100 1200 800]);

% Object amplitude
subplot(2,2,1);
imagesc(x_um,y_um,amplitude); axis image;
title('Object Amplitude (Ground Truth)'); xlabel('x (\mum)'); ylabel('y (\mum)');
colorbar; clim([0 1]); colormap gray;

% Object phase
subplot(2,2,2);
imagesc(x_um,y_um,ph); axis image;
title('Object Phase (Ground Truth)'); xlabel('x (\mum)'); ylabel('y (\mum)');
colorbar; colormap gray;

% Hologram 1
subplot(2,2,3); 
imagesc(x_um,y_um, hologram1); axis image;
title(sprintf('Hologram 1 (z_1 = %d \\mum)', z1*1e6)); 
xlabel('x (\mum)'); ylabel('y (\mum)'); colorbar; colormap gray;

% Hologram 2
subplot(2,2,4); 
imagesc(x_um,y_um, hologram2); axis image;
title(sprintf('Hologram 2 (z_2 = %d \\mum)', z2*1e6)); 
xlabel('x (\mum)'); ylabel('y (\mum)'); colorbar; colormap gray;


%% 6. PHỤC HỒI ẢNH (PLANE-TO-PLANE PHASE RETRIEVAL)
Iterations = 500;              

% Các nhân tố truyền sóng cần thiết cho vòng lặp
prop_dz       = Propagator(N, lambda, area, dz);    % Truyền xuôi: z1 -> z2
prop_minus_dz = Propagator(N, lambda, area, -dz);   % Truyền ngược: z2 -> z1
prop_minus_z1 = Propagator(N, lambda, area, -z1);   % Truyền ngược: z1 -> z=0 (mặt vật thể)

% Khởi tạo trường phức tại mặt z1 (Dùng biên độ đo được, pha ban đầu = 0)
u1 = amp1_measured .* exp(1i * zeros(N, N));

fprintf('\nĐang chạy vòng lặp phục hồi pha (Phase Retrieval) %d iterations...\n', Iterations);
figure('Name', 'Reconstruction Progress', 'Position', [200, 200, 1000, 450]);

for iter = 1:Iterations
    % [A] Lan truyền thuận: z1 -> z2
    u2_sim = IFT2Dc(FT2Dc(u1) .* prop_dz);
    
    % [B] Ràng buộc cường độ tại z2
    u2 = amp2_measured .* exp(1i * angle(u2_sim)); 
    
    % [C] Lan truyền ngược: z2 -> z1
    u1_sim = IFT2Dc(FT2Dc(u2) .* prop_minus_dz); 
    
    % [D] Ràng buộc cường độ tại z1
    u1 = amp1_measured .* exp(1i * angle(u1_sim));
    
    % Hiển thị tiến trình (mỗi 10 vòng lặp)
    if mod(iter, 10) == 0 || iter == 1 || iter == Iterations
        % Truyền trường u1 về lại z=0 để xem ảnh phục hồi
        u_obj_rec = IFT2Dc(FT2Dc(u1) .* prop_minus_z1);
        
        subplot(1,2,1); imagesc(x_um, y_um, abs(u_obj_rec)); axis image; colormap gray;
        title(sprintf('Rec. Amplitude (Iter: %d)', iter)); colorbar;
        xlabel('x (\mum)'); ylabel('y (\mum)');
        
        subplot(1,2,2); imagesc(x_um, y_um, angle(u_obj_rec)); axis image; colormap gray;
        title(sprintf('Rec. Phase (Iter: %d)', iter)); colorbar;
        xlabel('x (\mum)'); ylabel('y (\mum)');
        
        drawnow;
    end
end
fprintf('Hoàn thành phục hồi!\n');
%%
amp_rec = abs(u_obj_rec);

% Create figure
fig = figure('Color','w', ...
             'Name','Select Cross-section', ...
             'WindowStyle','normal');

ax = axes(fig);

imagesc(ax, x_um, y_um, amp_rec);

axis(ax,'image');
colormap(ax,gray);
colorbar(ax);

xlabel(ax,'x (\mum)');
ylabel(ax,'y (\mum)');
title(ax,'Draw a line between TWO points');

% Enable figure toolbar
fig.ToolBar = 'figure';

% Disable other mouse modes
zoom(fig,'off');
pan(fig,'off');
rotate3d(fig,'off');

% Select two endpoints by drawing a line
disp('Draw a line connecting two points on the image.');

h = drawline(ax,'Color','r','LineWidth',1.5);

% Wait until line is created
wait(h);

% Get two endpoint coordinates
pos = h.Position;

x_click = pos(:,1)';
y_click = pos(:,2)';

% =========================================
% EXTRACT PROFILE
% =========================================

num_points = 500;

x_line = linspace(x_click(1),x_click(2),num_points);
y_line = linspace(y_click(1),y_click(2),num_points);

[X_grid,Y_grid] = meshgrid(x_um,y_um);

profile_rec = interp2(X_grid,Y_grid,amp_rec, ...
                      x_line,y_line,'linear');

profile_true = interp2(X_grid,Y_grid,amplitude, ...
                       x_line,y_line,'linear');

distance = hypot(x_line-x_line(1), ...
                 y_line-y_line(1));

% =========================================
% PLOT PROFILE
% =========================================

figure('Color','w');

plot(distance,profile_true,'g','LineWidth',1.5);
hold on;

plot(distance,profile_rec,'r','LineWidth',1.5);

grid on;
xlabel('Distance (\mum)');
ylabel('Amplitude');

legend('Ground Truth','Reconstruction');
title('Amplitude Cross-section');

%% =========================================================================
%                             HÀM PHỤ TRỢ (FUNCTIONS)
% =========================================================================

function p = Propagator(N, lambda, area, z)
    % Hàm lan truyền phổ góc (Angular Spectrum Propagator)
    [X, Y] = meshgrid(1:N, 1:N);
    alpha = lambda * (Y - N/2 - 1) / area;
    beta  = lambda * (X - N/2 - 1) / area;
    
    % Mặt nạ lọc sóng suy biến (Evanescent waves)
    mask = (alpha.^2 + beta.^2) <= 1;
    
    p = zeros(N, N);
    p(mask) = exp(2 * pi * 1i * z * sqrt(1 - alpha(mask).^2 - beta(mask).^2) / lambda);
end

function out = FT2Dc(in)
    % Biến đổi Fourier 2D chuẩn hóa ở giữa tâm (Dùng nhân chập thay vì fftshift)
    [Nx, Ny] = size(in);
    vx = (-1).^(1:Nx).'; 
    vy = (-1).^(1:Ny);   
    f1 = vx * vy;        
    out = f1 .* fft2(f1 .* in);
end

function out = IFT2Dc(in)
    % Biến đổi Fourier ngược 2D chuẩn hóa ở giữa tâm
    [Nx, Ny] = size(in);
    vx = (-1).^(1:Nx).'; 
    vy = (-1).^(1:Ny);   
    f1 = vx * vy;        
    out = f1 .* ifft2(f1 .* in);
end