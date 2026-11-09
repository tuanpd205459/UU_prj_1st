
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORWARD SIMULATION OF NANO-OBJECT IN-LINE HOLOGRAPHY
%
% Object plane -> ASM propagation -> Computational CCD
%
% 1. USAF 1951 or Two Nano Dots
% 2. No cropping or resizing
% 3. USAF pixel size = 104 nm
% 4. Standard Angular Spectrum Method
% 5. Ideal computational CCD (no pixel integration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. SIMULATION PARAMETERS

lambda = 633e-9;        % Wavelength (m)
z = 0.005;              % Propagation distance (m)

object_mode = 'object'; % 'object' or '2dot'

%% 2. GENERATE OBJECT

switch lower(object_mode)

    case 'object'

        % Read USAF image
        img = imread('USAF_760micro_104nm.png');

        % Convert to grayscale
        if ndims(img) == 3
            img = rgb2gray(img);
        end

        % Normalize image
        object = double(img);
        object = object / max(object(:));

        % Image dimensions
        [Ny, Nx] = size(object);

        % Pixel size
        dx = 104e-9;
        dy = dx;

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

        dx = 20e-9;
        dy = dx;

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

%% 3. COORDINATE GRID

x = ((0:Nx-1)-(Nx-1)/2)*dx;
y = ((0:Ny-1)-(Ny-1)/2)*dy;

[X,Y] = meshgrid(x,y);

FOV_x = Nx*dx;
FOV_y = Ny*dy;

fprintf('\n========== SIMULATION PARAMETERS ==========\n');

fprintf('Object mode    : %s\n',object_mode);
fprintf('Grid size      : %d x %d pixels\n',Nx,Ny);
fprintf('Pixel size     : %.2f nm\n',dx*1e9);
fprintf('FOV            : %.3f x %.3f um\n', ...
    FOV_x*1e6,FOV_y*1e6);
fprintf('Wavelength     : %.1f nm\n',lambda*1e9);
fprintf('Propagation    : %.3f um\n',z*1e6);

%% 4. OBJECT FIELD

% Incident plane wave
U_inc = ones(Ny,Nx);

% Field after object
U0 = U_inc .* T;

%% 5. ANGULAR SPECTRUM METHOD

% Spatial frequency coordinates
fx = (-floor(Nx/2):ceil(Nx/2)-1)/FOV_x;
fy = (-floor(Ny/2):ceil(Ny/2)-1)/FOV_y;

[FX,FY] = meshgrid(fx,fy);

% Longitudinal wavevector
KZ = sqrt(complex(1/lambda^2 - FX.^2 - FY.^2));

% ASM transfer function
H = exp(1i*2*pi*z.*KZ);

% Forward propagation
A0 = fftshift(fft2(U0));

Az = A0 .* H;

Uz = ifft2(ifftshift(Az));

%% 6. COMPUTATIONAL CCD

CCD_amplitude = abs(Uz);
CCD_phase = angle(Uz);

CCD_intensity = abs(Uz).^2;
CCD_intensity_norm = CCD_intensity/max(CCD_intensity(:));

fprintf('\n========== CCD OUTPUT ==========\n');

fprintf('Amplitude min : %.6f\n',min(CCD_amplitude(:)));
fprintf('Amplitude max : %.6f\n',max(CCD_amplitude(:)));

fprintf('Intensity min : %.6f\n',min(CCD_intensity(:)));
fprintf('Intensity max : %.6f\n',max(CCD_intensity(:)));

%% 7. VISUALIZATION

x_um = x*1e6;
y_um = y*1e6;

figure('Name','Forward Simulation','Color','w', ...
    'Position',[100 100 1200 700]);

% Object mask
subplot(2,3,1);
imagesc(x_um,y_um,object_mask);
axis image;
title('Object Mask');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% Object amplitude
subplot(2,3,2);
imagesc(x_um,y_um,amplitude);
axis image;
title('Object Amplitude');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;
caxis([0 1]);

% Object phase
subplot(2,3,3);
imagesc(x_um,y_um,ph);
axis image;
title('Object Phase');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% CCD amplitude
subplot(2,3,4);
imagesc(x_um,y_um,CCD_amplitude);
axis image;
title('CCD Amplitude');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% CCD phase
subplot(2,3,5);
imagesc(x_um,y_um,CCD_phase);
axis image;
title('CCD Phase');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

% CCD intensity
subplot(2,3,6);
imagesc(x_um,y_um,CCD_intensity_norm);
axis image;
title('CCD Intensity (Normalized)');
xlabel('x (\mum)');
ylabel('y (\mum)');
colorbar;

colormap gray;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
