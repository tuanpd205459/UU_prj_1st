
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SLD INLINE DIGITAL HOLOGRAPHY
% SPECTRAL SUPERPOSITION + MULTI-PLANE PHASE RETRIEVAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% 1. SYSTEM PARAMETERS (mm)

% Central wavelength and spectral bandwidth
lambda_c = 0.405e-3;       % 633 nm
bandwidth = 0.004e-3;      % 20 nm FWHM

num_lambda = 50;

% Wavelength sampling
lambda = linspace(lambda_c-bandwidth/2, ...
                  lambda_c+bandwidth/2, num_lambda);

% Gaussian spectral distribution
S = exp(-4*log(2)*((lambda-lambda_c)/bandwidth).^2);

% Normalize spectral weights
weights = S/sum(S);

% Sensor
N = 500;
dx = 2e-3;                 % 2 um
area = dx*N;

% Propagation distances
z1 = 1.0;
z2 = 1.5;
dz = z2-z1;

%% 2. OBJECT (GROUND TRUTH)

% try
%     object0 = imread('a_object.jpg');
%     object = double(object0(:,:,1));
%     object = rot90(object,-1);
% 
%     object = imresize(object,[N N]);
% 
% catch
%     warning('Using simulated object.');
%     object = phantom('Modified Shepp-Logan',N);
% end
% 
% % Normalize object
% object = double(object);
% object = (object-min(object(:))) / ...
%          (max(object(:))-min(object(:)));
% 
% % Complex transmission function
% am = exp(-1.6*object);
% 
% ph0 = 3;
% ph = ph0*object;
% 
% t = am.*exp(1i*ph);
% 
% figure('Name','Ground Truth');
% 
% subplot(1,2,1);
% imagesc(am);
% axis image;
% colormap gray;
% colorbar;
% title('Object Amplitude');
% 
% subplot(1,2,2);
% imagesc(ph);
% axis image;
% colormap gray;
% colorbar;
% title('Object Phase');

%% 2. OBJECT (SINGLE NARROW SLIT)

% Slit parameters
slit_width = 10e-3;       % Slit width: 10 um (mm)
slit_length = 0.5;        % Slit length: 0.5 mm

% Convert physical dimensions to pixels
slit_width_px = max(1,round(slit_width/dx));
slit_length_px = max(1,round(slit_length/dx));

% Initialize object
object = zeros(N,N);

% Center position
cx = round(N/2);
cy = round(N/2);

% Slit boundaries
x1 = cx-floor(slit_width_px/2);
x2 = x1+slit_width_px-1;

y1 = cy-floor(slit_length_px/2);
y2 = y1+slit_length_px-1;

% Create vertical slit
object(y1:y2,x1:x2) = 1;

% Transmission function
am = object;
ph = zeros(N,N);

t = am.*exp(1i*ph);

% Display ground truth
figure('Name','Single Narrow Slit');

subplot(1,2,1);
imagesc(object);
axis image;
colormap gray;
colorbar;
title('Single Slit');

subplot(1,2,2);
imagesc(am);
axis image;
colormap gray;
colorbar;
title('Object Amplitude');


%% 3. SPECTRAL SUPERPOSITION

I1_total = zeros(N,N);
I2_total = zeros(N,N);

disp('Calculating SLD holograms...');

for j = 1:num_lambda

    % Current wavelength
    lambda_j = lambda(j);

    % Propagators
    prop_z1 = Propagator(N,lambda_j,area,z1);
    prop_z2 = Propagator(N,lambda_j,area,z2);

    % Propagation
    field_z1 = IFT2Dc(FT2Dc(t).*prop_z1);
    field_z2 = IFT2Dc(FT2Dc(t).*prop_z2);

    % Intensity holograms
    I1 = abs(field_z1).^2;
    I2 = abs(field_z2).^2;

    % Spectral intensity accumulation
    I1_total = I1_total + weights(j)*I1;
    I2_total = I2_total + weights(j)*I2;

end

% Final SLD holograms
hologram1 = I1_total;
hologram2 = I2_total;

% Measured amplitudes
amp1_measured = sqrt(hologram1);
amp2_measured = sqrt(hologram2);

figure('Name','SLD Holograms');

subplot(1,2,1);
imagesc(hologram1);
axis image;
colormap gray;
colorbar;
title('SLD Hologram (z1)');

subplot(1,2,2);
imagesc(hologram2);
axis image;
colormap gray;
colorbar;
title('SLD Hologram (z2)');

%% 4. MULTI-PLANE PHASE RETRIEVAL

% IMPORTANT:
% This reconstruction uses the central wavelength as an
% approximate single-wavelength reconstruction model.

Iterations = 50;

prop_dz = Propagator(N,lambda_c,area,dz);
prop_minus_dz = Propagator(N,lambda_c,area,-dz);
prop_minus_z1 = Propagator(N,lambda_c,area,-z1);

% Initial complex field
u1 = amp1_measured;

figure('Name','SLD Reconstruction');

disp('Starting phase retrieval...');

for iter = 1:Iterations

    % A. Forward propagation z1 -> z2
    u2_sim = IFT2Dc(FT2Dc(u1).*prop_dz);

    % B. Amplitude constraint at z2
    u2 = amp2_measured.*exp(1i*angle(u2_sim));

    % C. Backward propagation z2 -> z1
    u1_sim = IFT2Dc(FT2Dc(u2).*prop_minus_dz);

    % D. Amplitude constraint at z1
    u1 = amp1_measured.*exp(1i*angle(u1_sim));

    % Display reconstruction
    if mod(iter,10)==0 || iter==1 || iter==Iterations

        u_obj_rec = IFT2Dc(FT2Dc(u1).*prop_minus_z1);

        subplot(1,2,1);
        imagesc(abs(u_obj_rec));
        axis image;
        colormap gray;
        colorbar;
        title(sprintf('Amplitude - Iter %d',iter));

        subplot(1,2,2);
        imagesc(angle(u_obj_rec));
        axis image;
        colormap gray;
        colorbar;
        title(sprintf('Phase - Iter %d',iter));

        drawnow;
    end

end

disp('Reconstruction completed!');

%% 5. SPECTRAL DISTRIBUTION

figure('Name','SLD Spectrum');

plot(lambda*1e6,S,'o-','LineWidth',1.5);
xlabel('Wavelength (um)');
ylabel('Relative Spectral Intensity');
title('Gaussian SLD Spectrum');
grid on;

%% 6. HELPER FUNCTIONS

function p = Propagator(N,lambda,area,z)

    [X,Y] = meshgrid(1:N,1:N);

    alpha = lambda*(Y-N/2-1)/area;
    beta  = lambda*(X-N/2-1)/area;

    mask = (alpha.^2+beta.^2)<=1;

    p = zeros(N,N);

    p(mask) = exp(2*pi*1i*z.* ...
        sqrt(1-alpha(mask).^2-beta(mask).^2)/lambda);

end

function out = FT2Dc(in)

    [Nx,Ny] = size(in);

    vx = (-1).^(1:Nx).';
    vy = (-1).^(1:Ny);

    f1 = vx*vy;

    out = f1.*fft2(f1.*in);

end

function out = IFT2Dc(in)

    [Nx,Ny] = size(in);

    vx = (-1).^(1:Nx).';
    vy = (-1).^(1:Ny);

    f1 = vx*vy;

    out = f1.*ifft2(f1.*in);

end
