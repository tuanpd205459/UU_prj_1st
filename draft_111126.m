
clear; clc; close all;

%% 1. PARAMETERS
filename = 'USAF_220micro_104nm.png';
dx = 104e-9;                 % Pixel size: 104 nm

%% 2. READ IMAGE
img = imread(filename);

if size(img,3) == 3
    img = rgb2gray(img);
end

amplitude = double(img);

% Normalize amplitude
amplitude = amplitude / max(amplitude(:));

[Ny,Nx] = size(amplitude);

%% 3. PHYSICAL COORDINATES
x_um = (0:Nx-1)*dx*1e6;
y_um = (0:Ny-1)*dx*1e6;

FOV_x = Nx*dx*1e6;
FOV_y = Ny*dx*1e6;

%% 4. CROSS-SECTION EVALUATION

fig = figure('Color','w', ...
    'Name','Select Cross-section');

ax = axes(fig);

imagesc(ax,x_um,y_um,amplitude);
axis(ax,'image');
set(ax,'YDir','reverse');

colormap(ax,gray);
colorbar(ax);

xlabel(ax,'$x$ ($\mu$m)','Interpreter','latex');
ylabel(ax,'$y$ ($\mu$m)','Interpreter','latex');

title(ax,'Draw a line between TWO points');

% Draw adjustable line
h = drawline(ax,'Color','r','LineWidth',1.5);
wait(h);

% Get selected coordinates
pos = h.Position;

x_click = pos(:,1)';
y_click = pos(:,2)';

%% 5. EXTRACT AMPLITUDE PROFILE

num_points = 500;

x_line = linspace(x_click(1),x_click(2),num_points);
y_line = linspace(y_click(1),y_click(2),num_points);

[X_grid,Y_grid] = meshgrid(x_um,y_um);

profile_true = interp2( ...
    X_grid,Y_grid,amplitude, ...
    x_line,y_line,'linear');

% Physical distance along selected line
distance = hypot( ...
    x_line-x_line(1), ...
    y_line-y_line(1));

%% 6. PLOT AMPLITUDE CROSS-SECTION

figure('Color','w','Name','Amplitude Cross-section');

plot(distance,profile_true, ...
    'Color','#0072BD','LineWidth',1.5);

grid on;

xlabel('Distance ($\mu$m)','Interpreter','latex');
ylabel('Amplitude','Interpreter','latex');

title('USAF Amplitude Cross-section');

xlim([0 distance(end)]);

%% 7. DISPLAY PARAMETERS

fprintf('\n========== IMAGE PARAMETERS ==========\n');

fprintf('Image size     : %d x %d pixels\n',Nx,Ny);
fprintf('Pixel size     : %.2f nm\n',dx*1e9);

fprintf('FOV X          : %.3f um\n',FOV_x);
fprintf('FOV Y          : %.3f um\n',FOV_y);

fprintf('\n========== CROSS-SECTION ==========\n');

fprintf('Start point    : (%.3f, %.3f) um\n',...
    x_click(1),y_click(1));

fprintf('End point      : (%.3f, %.3f) um\n',...
    x_click(2),y_click(2));

fprintf('Profile length : %.3f um\n',distance(end));
fprintf('Samples        : %d\n',num_points);

%% 6. BAR WIDTH MEASUREMENT

% Threshold: midpoint between minimum and maximum amplitude
threshold = (min(profile_true) + max(profile_true))/2;

% Binary profile
binary_profile = profile_true > threshold;

% Find bar regions
edges = diff([0, binary_profile, 0]);

bar_starts = find(edges == 1);
bar_ends   = find(edges == -1) - 1;

% Calculate width of each bar
num_bars = length(bar_starts);
bar_widths = zeros(1,num_bars);

for k = 1:num_bars
    bar_widths(k) = distance(bar_ends(k)) - ...
                    distance(bar_starts(k));
end

% Calculate gaps
gaps = zeros(1,max(0,num_bars-1));

for k = 1:num_bars-1
    gaps(k) = distance(bar_starts(k+1)) - ...
              distance(bar_ends(k));
end

%% DISPLAY MEASUREMENT
fprintf('\n========== BAR WIDTH MEASUREMENT ==========\n');
fprintf('Threshold = %.4f\n',threshold);

for k = 1:num_bars
    fprintf('Bar %d width = %.2f nm\n', ...
        k,bar_widths(k)*1000);
end

for k = 1:length(gaps)
    fprintf('Gap %d = %.2f nm\n', ...
        k,gaps(k)*1000);
end

%% PLOT THRESHOLD AND BAR WIDTH
figure('Color','w');

plot(distance,profile_true,'b','LineWidth',1.5);
hold on;

yline(threshold,'r--','Threshold');

for k = 1:num_bars
    xline(distance(bar_starts(k)),'g--');
    xline(distance(bar_ends(k)),'g--');
end

grid on;
xlabel('Distance (\mum)');
ylabel('Amplitude');
title('Bar Width Measurement');
