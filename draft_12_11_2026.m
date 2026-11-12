
clc; clear; close all;

%% 1. INPUT PARAMETERS
WAVELENGTH = 633;       % Wavelength (nm)
N_PIXELS = 1024;        % Sensor size: 1024 x 1024 pixels

% Pixel sizes to investigate (um)
PIXEL_SIZES = [0.1, 0.25, 0.5, 1.0, 1.67, 2.0, 3.45, 5.0, 10.0];

%% 2. UNIT CONVERSION
wl = WAVELENGTH * 1e-9;       % Wavelength (m)
dx = PIXEL_SIZES * 1e-6;      % Pixel pitch (m)

%% 3. NYQUIST NA
na_nyquist = wl ./ (2 .* dx);

%% 4. EFFECTIVE NA
na_effective = min(na_nyquist, 1.0);

%% 5. OPTIMAL RECORDING DISTANCE
z_opt = (N_PIXELS .* dx.^2) ./ wl;

%% 6. ABBE RESOLUTION LIMIT
resolution = wl ./ (2 .* na_effective);

%% 7. PHYSICAL STATUS CLASSIFICATION
status = strings(size(PIXEL_SIZES));

for k = 1:length(PIXEL_SIZES)

    if na_nyquist(k) > 1.0
        status(k) = "Sóng evanescent / ngoài miền propagating";

    elseif abs(na_nyquist(k) - 1.0) < 1e-12
        status(k) = "Ngưỡng NA = 1";

    else
        status(k) = "Sóng propagating";
    end

end

%% 8. CREATE RESULTS TABLE
Results = table( ...
    PIXEL_SIZES(:), ...
    round(na_nyquist(:), 3), ...
    round(na_effective(:), 3), ...
    round(z_opt(:) * 1e3, 3), ...
    round(resolution(:) * 1e6, 3), ...
    status(:), ...
    'VariableNames', { ...
    'Pixel_Size_um', ...
    'NA_Nyquist', ...
    'Effective_NA', ...
    'Z_opt_mm', ...
    'Resolution_um', ...
    'Physical_Status'});

%% 9. DISPLAY RESULTS
fprintf('\n============================================\n');
fprintf('       HOLOGRAPHY SYSTEM PARAMETERS\n');
fprintf('============================================\n');

fprintf('Wavelength       : %d nm\n', WAVELENGTH);
fprintf('Sensor size      : %d x %d pixels\n', N_PIXELS, N_PIXELS);

fprintf('============================================\n\n');

disp(Results);
