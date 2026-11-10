clear; clc; close all;

%% PARAMETERS
lambda = 633e-9;       % Wavelength (m)
b = 220e-6;            % FOV (m)
dx = 104e-9;             % Sensor pixel size (m)

z = linspace(100e-9, 10e-3, 1000);  % Object-sensor distance (m)

%% 1. STANDING WAVE PERIODICITY
theta = atan(b ./ (2*z));

p = lambda ./ (2*sin(theta));

%% 2. NYQUIST-SHANNON CONDITION
% Pixel size must be <= half of periodicity

p_limit = 2*dx;

sampling_condition = dx <= p/2;

%% 3. FIND MINIMUM DISTANCE
idx = find(sampling_condition, 1, 'first');

if isempty(idx)
    z_min = NaN;
    fprintf('Sampling condition is NOT satisfied.\n');
else
    z_min = z(idx);
    fprintf('Minimum distance z_min = %.4f mm\n', z_min*1e3);
end

%% 4. PLOT
figure;

plot(z*1e3, p*1e6, 'b', 'LineWidth', 2);
hold on;

yline(p_limit*1e6, '--r', 'Nyquist limit: 2dx');

xlabel('Object-Sensor Distance z (mm)');
ylabel('Standing Wave Periodicity p (\mum)');

title('Nyquist-Shannon Sampling Condition');
legend('Standing wave periodicity', 'Nyquist limit');
grid on;

%% 5. DISPLAY RESULT
fprintf('Wavelength       = %.1f nm\n', lambda*1e9);
fprintf('FOV              = %.1f um\n', b*1e6);
fprintf('Pixel size       = %.2f um\n', dx*1e6);
fprintf('Required p       >= %.2f um\n', 2*dx*1e6);




%% Calculate z_min
lambda = 633e-9;  % 633 nm
b = 220e-6;       % 760 um
dx = 104e-9;        % 2 um

z_min = (b/2)*sqrt((4*dx/lambda)^2 - 1);

fprintf('z_min = %.3f mm\n', z_min*1e3);


