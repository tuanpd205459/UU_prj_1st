% 
% % ============================================================
% % USAF 1951: DISPLAY + COORDINATE SELECTION + PHYSICAL CROP
% % Pixel size = 104 nm/pixel
% % ============================================================
% 
% clear; clc; close all;
% 
% %% 1. PIXEL SIZE
% pixel_size_nm = 104;
% pixel_size_um = pixel_size_nm / 1000;
% 
% %% 2. READ IMAGE
% img = imread('USAF_98area_104nm.png');
% 
% [height, width, ~] = size(img);
% 
% fprintf('Original image: %d x %d pixels\n', width, height);
% fprintf('Pixel size: %.1f nm\n', pixel_size_nm);
% fprintf('Physical FOV: %.3f x %.3f um\n', ...
%     width*pixel_size_um, height*pixel_size_um);
% 
% %% 3. DISPLAY ORIGINAL IMAGE WITH PHYSICAL COORDINATES
% 
% x_thuc = (0:width-1) * pixel_size_um;
% y_thuc = (0:height-1) * pixel_size_um;
% 
% figure('Name','Original USAF Image');
% 
% imshow(img, 'XData', x_thuc, 'YData', y_thuc);
% 
% axis on;
% axis image;
% xlabel('X (\mum)');
% ylabel('Y (\mum)');
% title('Original USAF 1951 ');
% 
%%
% ============================================================
% USAF 1951: DISPLAY + COORDINATE SELECTION + PHYSICAL CROP
% Original pixel size = 104 nm/pixel
% Crop N x N pixels WITHOUT resizing
% ============================================================

clear; clc; close all;

%% 1. PIXEL SIZE
pixel_size_nm = 104;
pixel_size_um = pixel_size_nm / 1000;

%% 2. READ IMAGE
img = imread('USAF-1951_104nm_tilt_C.png');

[height, width, ~] = size(img);

fprintf('Original image: %d x %d pixels\n', width, height);
fprintf('Pixel size: %.1f nm\n', pixel_size_nm);
fprintf('Physical FOV: %.3f x %.3f um\n', ...
    width*pixel_size_um, height*pixel_size_um);

%% 3. DISPLAY ORIGINAL IMAGE WITH PHYSICAL COORDINATES

x_thuc = (0:width-1) * pixel_size_um;
y_thuc = (0:height-1) * pixel_size_um;

figure('Name','Original USAF Image');

imshow(img, 'XData', x_thuc, 'YData', y_thuc);
axis on;
axis image;
xlabel('X (\mum)');
ylabel('Y (\mum)');
title('Original USAF 1951');


N = 2000;  % 2000 x 2000 pixels

%% 5.

figure('Name','Select Crop Region');

imshow(img);
axis image;
hold on;

x_center = round(width/2);
y_center = round(height/2);

h = drawrectangle('Position', ...
    [x_center-N/2, y_center-N/2, N, N], ...
    'Color','r', ...
    'LineWidth',2, ...
    'FixedAspectRatio',true, ...
    'AspectRatio',1);

disp('Move the crop window and double-click to confirm.');

wait(h);

pos = round(h.Position);

x_start = max(1, min(pos(1), width-N+1));
y_start = max(1, min(pos(2), height-N+1));

x_end = x_start + N - 1;
y_end = y_start + N - 1;

%% 6. CROP WITHOUT RESIZING

img_crop = img(y_start:y_end, x_start:x_end, :);

%% 9. DISPLAY CROPPED IMAGE

figure('Name','Cropped USAF Image');

imshow(img_crop);
axis image;
title(sprintf('Cropped Image: %d x %d pixels', N, N));

%% 10. PHYSICAL DIMENSIONS

crop_width_um = N * pixel_size_um;
crop_height_um = N * pixel_size_um;

fprintf('\n========== CROP INFORMATION ==========\n');
fprintf('Crop size: %d x %d pixels\n', N, N);
fprintf('Pixel size: %.1f nm/pixel\n', pixel_size_nm);
fprintf('Physical size: %.3f x %.3f um\n', ...
    crop_width_um, crop_height_um);

%% 11. 

imwrite(img_crop, 'USAF_crop.png');

disp('Cropped image saved as USAF_crop.png');