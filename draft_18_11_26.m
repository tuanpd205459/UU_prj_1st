
clc; clear; close all;

%% 1. SLD PARAMETERS
lambda0 = 633e-9;       % Central wavelength (m)
FWHM = 20e-9;           % Spectral bandwidth (m)

Nlambda = 501;
lambda = linspace(lambda0-3*FWHM, ...
                  lambda0+3*FWHM, Nlambda);

% Gaussian spectral power density
S = exp(-4*log(2)*((lambda-lambda0)/FWHM).^2);
S = S / sum(S);

%% 2. TEMPORAL COHERENCE FUNCTION
OPD = linspace(-50e-6,50e-6,1001);

gamma = zeros(size(OPD));

for k = 1:length(OPD)
    gamma(k) = sum(S .* exp(1i*2*pi*OPD(k)./lambda));
end

gamma = gamma / sum(S);

% Approximate 1/e coherence length for Gaussian spectrum
Lc = 2*sqrt(log(2))/pi * lambda0^2/FWHM;

%% 3. INTERFERENCE SIMULATION
Nx = 2000;
x = linspace(-1e-3,1e-3,Nx);

theta = 0.5*pi/180;
Aobj = 0.8;
Aref = 1;

OPD0 = [0, 20e-6];

I_SLD = zeros(2,Nx);

for q = 1:2

    opd_map = OPD0(q) + x*sin(theta);

    for k = 1:Nlambda

        phase = 2*pi*opd_map/lambda(k);

        Uref = Aref;
        Uobj = Aobj .* exp(1i*phase);

        % Incoherent spectral intensity integration
        I_SLD(q,:) = I_SLD(q,:) + ...
            S(k)*abs(Uref + Uobj).^2;
    end
end

% Monochromatic laser comparison
I_Laser = zeros(2,Nx);

for q = 1:2
    opd_map = OPD0(q) + x*sin(theta);

    I_Laser(q,:) = abs(Aref + ...
        Aobj*exp(1i*2*pi*opd_map/lambda0)).^2;
end

%% 4. PLOT RESULTS
figure('Color','w','Position',[100 100 1200 800]);

subplot(2,2,1);
plot(lambda*1e9,S,'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Normalized power');
title('SLD Spectrum');
grid on;

subplot(2,2,2);
plot(OPD*1e6,abs(gamma),'LineWidth',2);
hold on;
xline(Lc*1e6,'--r','1/e coherence length');
xlabel('OPD (\mum)');
ylabel('|\gamma(OPD)|');
title('Temporal Coherence Function');
grid on;

subplot(2,2,3);
plot(x*1e3,I_Laser(1,:),'--','LineWidth',1.5);
hold on;
plot(x*1e3,I_SLD(1,:),'LineWidth',1.5);
xlabel('Position (mm)');
ylabel('Intensity');
title('OPD = 0 \mum');
legend('Monochromatic laser','SLD');
grid on;

subplot(2,2,4);
plot(x*1e3,I_Laser(2,:),'--','LineWidth',1.5);
hold on;
plot(x*1e3,I_SLD(2,:),'LineWidth',1.5);
xlabel('Position (mm)');
ylabel('Intensity');
title('OPD = 20 \mum');
legend('Monochromatic laser','SLD');
grid on;

sgtitle('Low Temporal Coherence Validation');

fprintf('Approximate 1/e coherence length: %.2f um\n',Lc*1e6);
fprintf('Coherence at OPD = 20 um: %.4f\n', ...
    abs(sum(S.*exp(1i*2*pi*20e-6./lambda))/sum(S)));
% 