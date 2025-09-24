% RADARSAT-1 Data Processing: RDA Imaging with SRC2
% Accurate implementation of RDA with secondary range compression
% Author: Hanyang Xu

clc; clear all; close all;

%% Extract raw data
% run specify_parameters.m
% run extract_data.m
clear all;

%% Parameter setup
% Radar parameters
PLOT = 0;
C = 2.9979e8;                                                               % Speed of light
f0 = 5.3e9;                                                                 % Center frequency
lambda = C / f0;                                                            % Wavelength
% Platform parameters
V = 7062;                                                                   % Platform velocity
% Range parameters
Kr = -0.72135e12;                                                           % Chirp rate
Tr = 41.75e-6;                                                              % Pulse duration
Fr = 32.317e6;                                                              % Range sampling rate
dtau = 1 / Fr;                                                              % Range sampling interval
% Azimuth parameters
Ka = 1733;                                                                  % Azimuth chirp rate
Fa = 1256.98;                                                               % Azimuth sampling rate
fc = -6900;                                                                 % Doppler centroid frequency
deta = 1 / Fa;                                                              % Azimuth sampling interval
% Other parameters
t0 = 6.5956e-3;                                                             % Data acquisition time
R0 = t0 * C / 2;                                                            % Nearest slant range

%% Load raw data
load CDdata1.mat
data = double(data);

% load CDdata2.mat
% data2 = double(data);
% data = [data1; data2];
[Na, Nr] = size(data);

%% Doppler centroid estimation
% fc = -6900;                                                              % Given by the CD
PRF = Fa;
DeltaR = C / (2 * Fr);
chirp_rg_BW = Kr * Tr;
La = 2 * DeltaR;
beam_width_az = lambda / La * 180 / pi;
% Za = 800;
% ta = (Na + Za - 1) * deta;
% fc = DopplerCenterEstimation_Eny(data, PRF)               % Good: -7050

%% Determine reference range and reference time
Rref = (t0 + Nr / 2 * dtau) * C / 2;                                                  % Reference range selected at the center of the imaging area

%% Zero-padding at the end of the data
% Za = 800;                                                                   % Azimuth zero-padding number
Zr = ceil(Tr / dtau);                                                         % Range zero-padding number
data = cat(2, data, zeros(Na, Zr));                                            % Range zero-padding
% data = cat(1, zeros(Za, Nr + Zr), data);                                         % Azimuth zero-padding
% Na = Na + Za;
Nr = Nr + Zr;

% figure, imagesc(abs(data)); axis image; set(gcf, 'Color', 'w');
% title('Time Domain: Amplitude of the Signal After Zero-Padding');
% xlabel('Range (Samples)'); ylabel('Azimuth (Samples)');

%% Time and frequency axis setup
tau = t0 + (0:Nr-1) * dtau;                                                     % Range time axis
ftau = ((0:Nr-1) - Nr / 2) / Nr * Fr;                                               % Range frequency axis
eta = ((0:Na-1) - Na / 2) * deta;                                                 % Azimuth time axis
feta = fc + ((0:Na-1) - Na / 2) / Na * Fa;                                            % Azimuth frequency axis

%% Align to the Doppler centroid

S0 = data .* exp(1i * 2 * pi * fc * (eta' * ones(1, Nr)));                              % Shift to the Doppler centroid frequency
clear data

%% Range Compression (RC)
Tref = t0 + Nr / 2 * dtau;                                                        % Reference time in range
ht_rc = (abs((ones(Na, 1) * tau) - Tref) < Tr / 2) .* exp(1i * pi * Kr * ((ones(Na, 1) * tau) - Tref).^2); % Matched filter (Method 2)
Hf_rc = conj(fftshift(fft(fftshift(ht_rc, 2), [], 2));                      % Range FFT of the matched filter
clear ht_rc
S0_rc = fftshift(fft(fftshift(S0, 2), [], 2);                               % Range FFT of the signal
clear S0
S1 = ifftshift(ifft(ifftshift(S0_rc .* Hf_rc, 2), [], 2);                     % Range IFFT after matched filtering
clear Hf_rc S0_rc
if PLOT
figure, imagesc(abs(S1)); axis image; set(gcf, 'Color', 'w');
title('Time Domain: Signal Amplitude After Range Compression');
xlabel('Range (Samples)'); ylabel('Azimuth (Samples)');
end

%% Secondary Range Compression (SRC)
Sff = fftshift(fft2(fftshift(S1)));                                         % 2D FFT of the signal
WindowR = ones(Na, 1) * kaiser(Nr, 2.5)';                                       % Range window
WindowA = kaiser(Na, 2.5) * ones(1, Nr);                                        % Azimuth window
Sff = Sff .* WindowR .* WindowA;                                                % Apply windows
clear WindowR WindowA
D = sqrt(1 - lambda^2 * feta.^2 / (4 * V^2));                                       % Range migration factor
Ksrc = (2 * V^2 * f0^3 * (D' * ones(1, Nr)).^3) ./ (C * R0 * (feta' * ones(1, Nr)).^2);       % SRC filter chirp rate
Hsrc = exp(-1i * pi * ((ones(Na, 1) * ftau).^2) ./ Ksrc);                            % SRC filter
clear Ksrc
S1 = Sff .* Hsrc;                                                             % Filtering
clear Sff Hsrc
Srd = ifftshift(ifft(ifftshift(S1, 2), [], 2);                              % Range IFFT after SRC correction, transforming to the range-Doppler domain
clear S1
if PLOT
figure, imagesc(abs(Srd)); axis image; set(gcf, 'Color', 'w');
title('Range-Doppler Domain: Signal Amplitude After SRC');
xlabel('Range (Samples)'); ylabel('Azimuth (Samples)');
end

%% Range Cell Migration Correction (RCMC)
deltaR = dtau * C / 2;                                                          % Range sampling interval
deltaRCM = Rref * (1 ./ D - 1);                                                   % Range migration amount
deltaRCM = deltaRCM / deltaR;                                                 % Quantization of range migration amount
IntNum = ceil(deltaRCM);                                                    % Round up the range migration amount
DecNum = IntNum - deltaRCM;                                                   % Fractional part of the range migration amount quantized to (1/16) multiples
Srcmc = zeros(Na, Nr);
h = waitbar(0, 'Sinc Interpolation');
for ii = 1:Na
    SincVal = sinc(DecNum(ii) - 4:1:DecNum(ii) + 3)';                           % Generate interpolation kernel
    for jj = 1:Nr
        kk = jj + IntNum(ii);                                                 % Position of range migration in the original matrix
        if (5 <= kk && kk <= Nr - 3)                                               % Prevent index overflow
            Srcmc(ii, jj) = Srd(ii, kk - 4:kk + 3) * SincVal;                       % Interpolation
        end
    end
    waitbar(ii / Na);
end
close(h);
% clear deltaRCM IntNum DecNum Srd
if PLOT
figure, imagesc(abs(Srcmc)); axis image; set(gcf, 'Color', 'w');
title('Range-Doppler Domain: Signal Amplitude After RCMC');
xlabel('Range (Samples)'); ylabel('Azimuth (Samples)');
end

%% Azimuth Compression (AC)
% Filtering first
% load Srcmc.mat
% S7 = xuhanyang(eta, V, lambda, Fa, Nr, Na, Srcmc);
Haz = exp(1i * 4 * pi * R0 * (D' * ones(1, Nr)) * f0 / C);                                 % Azimuth matched filter
clear D
S2 = Srcmc .* Haz;                                                            % Matched filtering
clear Haz 
Stt = ifftshift(ifft(ifftshift(S2, 1), [], 1);                              % Azimuth IFFT
clear S2

%% Adjust the image
% Img = (abs(Stt));                                                     % Flip the image vertically (flipud)
% Ntau = round((Tr / dtau / 2 - 1) / 2);                                              % Calculate the length of the range guard zone
% 
% Img = Img(1:Na, Ntau + 1:Ntau + Nr - Zr);                                       % Crop the effective region of the image
% Img = Img / max(max(Img));
% Img = 20 * log10(Img + eps);
% Img(Img < -50) = -50;
% Img = uint8((Img + 50) / 50 * 255);
% 
% figure, imagesc(Img)); axis image; set(gcf, 'Color', 'w');
% xlabel('Range');
% ylabel('Azimuth');
% imwrite(Img, 'SAR_RDA_2.bmp', 'bmp');

Img = (abs(Stt));                                                     % Flip the image vertically (flipud)
Ntau = round((Tr / dtau / 2 - 1) / 2);                                              % Calculate the length of the range guard zone
Img = Img(1:Na, Ntau + 1:Ntau + Nr - Zr);                                       % Crop the effective region of the image
Img = Img / max(max(Img));
Img = 20 * log10(Img);
Img(Img < -50) = -50;
% Img = uint8((Img + 50) / 50 * 255);

% Extract part of the scene
% Img2 = Img(97:680, 666:980);
figure, imagesc(Img, [-50 0]); axis image; set(gcf, 'Color', 'w');
title('RDA Imaging Result')
colormap(turbo);
colorbar
% figure, imagesc(Img(688:944, :)); axis image; set(gcf, 'Color', 'w');
% imwrite(Img, 'SAR_CSA.bmp', 'bmp'); % imwrite(A, filename, fmt), A is the image data, filename is the target image name, fmt is the format of the image to be generated
