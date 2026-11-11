
clear; clc; close all;

%% 1. PARAMETERS
dx = 100e-9;          % Pixel size: 100 nm
bar_width = 200e-9;   % Bar width: 200 nm
gap = 500e-9;         % Gap: 500 nm

bar_length = 1e-6;    % Bar length: 1 um
N = 256;              % Image size: 256 x 256 pixels

%% 2. CREATE OBJECT
object = zeros(N,N);

bar_px = round(bar_width/dx);
gap_px = round(gap/dx);
length_px = round(bar_length/dx);

total_width = 3*bar_px + 2*gap_px;

% Center the pattern
start_col = floor((N-total_width)/2) + 1;
start_row = floor((N-length_px)/2) + 1;

for k = 0:2
    c1 = start_col + k*(bar_px+gap_px);
    c2 = c1 + bar_px - 1;

    object(start_row:start_row+length_px-1,c1:c2) = 1;
end

%% 3. DISPLAY
x = (0:N-1)*dx*1e6;
y = (0:N-1)*dx*1e6;

figure;
imagesc(x,y,object);
colormap(gray);
axis image;
set(gca,'YDir','reverse');

xlabel('X (\mum)');
ylabel('Y (\mum)');
title('Three-bar target');
colorbar;

%% 4. PHYSICAL PARAMETERS
FOV_x = N*dx;
FOV_y = N*dx;

actual_bar_width = bar_px*dx;
actual_gap = gap_px*dx;
actual_length = length_px*dx;

pattern_width = 3*actual_bar_width + 2*actual_gap;
pattern_height = actual_length;

center_distance = actual_bar_width + actual_gap;

%% 5. DISPLAY PARAMETERS
fprintf('\n====================================\n');
fprintf('       OBJECT PARAMETERS\n');
fprintf('====================================\n');

fprintf('\n[Sampling]\n');
fprintf('Pixel size       : %.2f nm\n',dx*1e9);
fprintf('Image size       : %d x %d pixels\n',N,N);

fprintf('\n[Field of View]\n');
fprintf('FOV X            : %.3f um\n',FOV_x*1e6);
fprintf('FOV Y            : %.3f um\n',FOV_y*1e6);
fprintf('Total FOV        : %.3f x %.3f um^2\n',...
    FOV_x*1e6,FOV_y*1e6);

fprintf('\n[Individual Bar]\n');
fprintf('Number of bars   : 3\n');
fprintf('Bar width        : %.2f nm\n',actual_bar_width*1e9);
fprintf('Bar length       : %.2f um\n',actual_length*1e6);
fprintf('Bar pixels       : %d x %d pixels\n',...
    bar_px,length_px);

fprintf('\n[Spacing]\n');
fprintf('Edge-to-edge gap : %.2f nm\n',actual_gap*1e9);
fprintf('Center distance  : %.2f nm\n',center_distance*1e9);

fprintf('\n[Total Pattern]\n');
fprintf('Pattern width    : %.2f um\n',pattern_width*1e6);
fprintf('Pattern height   : %.2f um\n',pattern_height*1e6);
fprintf('Pattern size     : %.2f x %.2f um^2\n',...
    pattern_width*1e6,pattern_height*1e6);

fprintf('\n====================================\n');
%% SAVE OBJECT
save('three_bar_object.mat', ...
    'object', 'dx', 'N', ...
    'bar_width', 'gap', 'bar_length');

fprintf('\nObject saved: three_bar_object.mat\n');

