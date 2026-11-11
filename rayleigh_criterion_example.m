
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAYLEIGH CRITERION USING ANGULAR SPECTRUM METHOD (ASM)
%
% Object plane -> ASM -> Thin lens -> ASM -> Image plane
%
% Two incoherent point sources
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc; close all;

%% 1. SYSTEM PARAMETERS

lambda = 633e-9;       % Wavelength (m)
NA = 0.05;             % Numerical aperture
n = 1.0;               % Refractive index

f = 1e-3;              % Focal length (m)

% Rayleigh criterion
d_R = 0.61*lambda/NA;

fprintf('Rayleigh resolution = %.3f um\n',d_R*1e6);

%% 2. NUMERICAL GRID

N = 1024;
dx = 0.5e-6;

L = N*dx;

x = (-N/2:N/2-1)*dx;
[X,Y] = meshgrid(x,x);

k = 2*pi*n/lambda;

% Frequency coordinates
df = 1/L;

fx = (-N/2:N/2-1)*df;
[FX,FY] = meshgrid(fx,fx);

%% 3. ANGULAR SPECTRUM TRANSFER FUNCTION

kz = k*sqrt(complex( ...
    1-(lambda/n)^2*(FX.^2+FY.^2)));

H = exp(1i*kz*(2*f));

%% 4. THIN LENS AND CIRCULAR APERTURE

% Aperture radius from NA = n*sin(theta)
theta = asin(NA/n);
R = f*tan(theta);

aperture = double(X.^2+Y.^2 <= R^2);

% Thin lens phase
lens = exp(-1i*k/(2*f)*(X.^2+Y.^2));

% Combined lens and aperture
Lens = lens.*aperture;

fprintf('Aperture diameter = %.2f um\n',2*R*1e6);

%% 5. DEFINE TWO POINT SOURCES

% Test different source separations
separations = [0.5 1.0 1.22 1.5 2.0]*d_R;

% Store intensity profiles
I_profiles = zeros(length(separations),N);

%% 6. SIMULATION LOOP

for s = 1:length(separations)

    d = separations(s);

    % Point source positions
    x1 = -d/2;
    x2 =  d/2;

    % Initialize object fields
    U01 = zeros(N,N);
    U02 = zeros(N,N);

    % Source 1
    [~,ix1] = min(abs(x-x1));
    U01(N/2+1,ix1) = 1;

    % Source 2
    [~,ix2] = min(abs(x-x2));
    U02(N/2+1,ix2) = 1;

    %% Source 1 propagation

    U1 = asm_propagation(U01,H);
    U1 = U1.*Lens;
    U1 = asm_propagation(U1,H);

    %% Source 2 propagation

    U2 = asm_propagation(U02,H);
    U2 = U2.*Lens;
    U2 = asm_propagation(U2,H);

    %% Incoherent intensity addition

    I1 = abs(U1).^2;
    I2 = abs(U2).^2;

    I_total = I1+I2;

    % Extract central cross-section
    profile = I_total(N/2+1,:);

    % Normalize
    profile = profile/max(profile);

    I_profiles(s,:) = profile;

end

%% 7. PLOT INTENSITY PROFILES

figure('Color','w');

for s = 1:length(separations)

    subplot(3,2,s);

    plot(x*1e6,I_profiles(s,:), ...
        'b','LineWidth',1.5);

    hold on;

    title(sprintf('d = %.2f d_R', ...
        separations(s)/d_R));

    xlabel('Position (um)');
    ylabel('Normalized intensity');

    grid on;
    xlim([-20 20]);
    ylim([0 1.1]);

end

sgtitle('Rayleigh Criterion using ASM');

%% 8. OVERLAY COMPARISON

figure('Color','w');
hold on;

for s = 1:length(separations)

    plot(x*1e6,I_profiles(s,:), ...
        'LineWidth',1.5, ...
        'DisplayName',sprintf('%.2f d_R', ...
        separations(s)/d_R));

end

xlabel('Position (um)');
ylabel('Normalized intensity');

title('Comparison of Two-Point Resolution');
legend('Location','best');
grid on;
xlim([-20 20]);
%% Check khoang cach giua 2 vat
lambda = 633e-9;     % Wavelength (m)
b = 106e-6;          % Sensor width (m)
z = 1e-3;            % Object-sensor distance (m)
n = 1.0;             % Refractive index

D = 200e-9;          % Object diameter (m)
d = 400e-9;          % Center-to-center distance (m)

% CALCULATE NA
theta = atan(b/(2*z));
NA = n*sin(theta);

% RAYLEIGH CRITERION
d_R = 0.61*lambda/NA;

% OBJECT DISTANCE
gap = d-D;

% RESULTS
fprintf('NA                  = %.4f\n', NA);
fprintf('Rayleigh distance   = %.2f nm\n', d_R*1e9);
fprintf('Object diameter     = %.2f nm\n', D*1e9);
fprintf('Center distance     = %.2f nm\n', d*1e9);
fprintf('Edge-to-edge gap    = %.2f nm\n', gap*1e9);

% CHECK
if d >= d_R
    fprintf('\nPASS: Meets Rayleigh criterion\n');
else
    fprintf('\nFAIL: Below Rayleigh criterion\n');
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION: ANGULAR SPECTRUM PROPAGATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Uz = asm_propagation(U0,H)

    % Centered Fourier transform
    A = fftshift(fft2(ifftshift(U0)));

    % Angular spectrum propagation
    B = A.*H;

    % Inverse Fourier transform
    Uz = fftshift(ifft2(ifftshift(B)));

end
