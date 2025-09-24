% RADARSAT-1 Data Processing: CSA Imaging
% Author: Hanyang Xu

% clc; clear all; close all;
% extract = 0;

%% Extract raw data
% if extract
% run specify_parameters.m
% run extract_data.m
% end
clear all;
PLOT = 0;

%% Parameter setup
% Radar parameters
C = 2.9979e8;                                                               % Speed of light
f0 = 5.3e9;                                                                 % Center frequency
lambda = C / f0;

% Platform parameters
V = 7062;                                                                   % Effective platform velocity

% Range parameters
Kr = -0.72135e12;                                                           % Chirp rate
Tr = 41.75e-6;                                                              % Pulse duration
Fr = 32.317e6;                                                              % Range sampling rate
dtau = 1 / Fr;                                                              % Range sampling interval

% Azimuth parameters
Ka = 1733;                                                                  % Azimuth chirp rate
Fa = 1256.98;                                                               % Azimuth sampling rate
fc = -6900;                                                                 % Doppler centroid frequency, squint angle
deta = 1 / Fa;                                                              % Azimuth sampling interval

% Other parameters
t0 = 6.5956e-3;                                                             % Data acquisition time, start sampling time
t0 = t0 + (900 - 1) * dtau;
R0 = t0 * C / 2;                                                            % Nearest slant range

%% Explanation
% Reference target: Selected at the center of the swath
% Reference range: Nearest slant range of the reference target
% Reference frequency: Center frequency of the reference target

%% Load raw data
% Change the current directory to newfolder
% oldFolder = cd('.\scene01');
load CDdata1.mat
% data1 = data;
% load CDdata2.mat
% data2 = data;
% data = [data1; data2];
% cd(oldFolder);
data = double(data);
[Na, Nr] = size(data); % Azimuth * Range

%% Reference data setup
Rref = (t0 + Nr / 2 * dtau) * C / 2;                                        % Reference range, taken as the middle value of the sampled range
Vref = V;                                                                   % Reference velocity
fref = fc;                                                                  % Reference frequency, Doppler centroid frequency

%% Zero-padding at the end of the data
% Za = 800;                                                                 % Azimuth zero-padding number
Zr = ceil(Tr / dtau);                                                       % Range zero-padding number
data = cat(2, data, zeros(Na, Zr));                                         % Range zero-padding
% data = cat(1, zeros(Za, Nr + Zr), data);                                  % Azimuth zero-padding
% Na = Na + Za;
Nr = Nr + Zr;
if PLOT
    figure, imagesc(abs(data)); axis image; set(gcf, 'Color', 'w');
    title('Time Domain: Amplitude of the Signal After Zero-Padding');
    xlabel('Range (Samples)'); ylabel('Azimuth (Samples)');
end

%% Time and frequency axis setup
tau = t0 + (0:Nr-1) * dtau;                                                 % Range time axis
ftau = ((0:Nr-1) - Nr / 2) / Nr * Fr;                                       % Range frequency axis, equivalent to dividing the sampling frequency into N parts
eta = ((0:Na-1) - Na / 2) * deta;                                           % Azimuth time axis, half division
feta = fc + ((0:Na-1) - Na / 2) / Na * Fa;                                  % Azimuth frequency axis

%% Intermediate variable setup
D_feta_V = sqrt(1 - C^2 * feta.^2 / (4 * V^2 * f0^2));                      % D(feta, V)
D_feta_Vref = sqrt(1 - C^2 * feta.^2 / (4 * Vref^2 * f0^2));                % D(feta, Vref)
D_fref_V = sqrt(1 - C^2 * fref^2 / (4 * V^2 * f0^2));                       % D(fref, V)
D_fref_Vref = sqrt(1 - C^2 * fref^2 / (4 * Vref^2 * f0^2));                 % D(fref, Vref)

Ksrc = (2 * V^2 * f0^3 * D_feta_V.^3) ./ (C * R0 * feta.^2);                % SRC filter chirp rate
Km = Kr ./ (1 - Kr ./ Ksrc);                                                % Modified range chirp rate
clear Ksrc
RCMbulk = (1 ./ D_feta_Vref - 1 / D_fref_Vref) * Rref;                      % Bulk RCM
alpha = D_fref_Vref ./ D_feta_Vref - 1;                                     % Frequency scaling parameter

%% Align to the Doppler centroid
S0 = data .* exp(-1i * 2 * pi * fc * (eta' * ones(1, Nr)));                 % Shift to the Doppler centroid frequency
% clear data

%% Transform to the range-Doppler domain, implement scaling phase multiplication
Srd = fftshift(fft(fftshift(S0, 1), [], 1);                                 % Azimuth FFT of the original signal
clear S0
tt = 2 / C * (R0 / D_fref_V + RCMbulk) - 2 * Rref ./ (C * D_feta_Vref);     % P205 (7.26) (7.27)
Ssc = exp(1i * pi * Km .* alpha .* tt.^2);                                  % Scaling equation P207 (7.30)
clear tt
S1 = Srd .* (Ssc' * ones(1, Nr));                                           % Scaling phase multiplication
clear Srd Ssc

%% Transform to the 2D frequency domain, implement RC, SRC, and bulk RCMC
S2 = fftshift(fft(fftshift(S1, 2), [], 2);                                  % Transform the signal to the 2D frequency domain
clear S1
WindowR = ones(Na, 1) * kaiser(Nr, 2.4)';                                   % Range window, higher kaiser value means sharper window
WindowA = kaiser(Na, 2.4) * ones(1, Nr);                                    % Azimuth window
S2 = S2 .* WindowR .* WindowA;                                              % Apply windows
clear WindowR WindowA
Hm = exp(1i * pi ./ ((Km' * ones(1, Nr)) .* (1 + alpha' * ones(1, Nr))) .* (ones(Na, 1) * ftau).^2); % Combined range compression and bulk RCMC filter
clear alpha
Hrcm = exp(1i * 4 * pi / C * (RCMbulk' * ones(1, Nr)) .* (ones(Na, 1) * ftau));             % SRC filter
clear RCMbulk
S3 = S2 .* Hm .* Hrcm;                                                      % Phase multiplication, implement filtering
clear S2 Hm Hrcm

%% Transform to the range-Doppler domain, implement azimuth compression and additional phase correction
S4 = ifftshift(ifft(ifftshift(S3, 2), [], 2);                               % Range IFFT of the signal
clear S3
if PLOT
    figure, imagesc(abs(S4)); axis image; set(gcf, 'Color', 'w');
    title('Range-Doppler Domain: Signal Amplitude After RC, SRC, and Bulk RCMC');
    xlabel('Range (Samples)'); ylabel('Azimuth (Samples)');
end
Hac = exp(-1i * 4 * pi * R0 * f0 * D_feta_V / C);                           % Azimuth compression filter
Hpc = exp(-1i * (4 * pi * Km / C^2) .* (1 - D_feta_Vref / D_fref_Vref) .* (R0 ./ D_feta_V - Rref ./ D_feta_V).^2); % Phase correction filter
S5 = S4 .* (Hac' * ones(1, Nr)) .* (Hpc' * ones(1, Nr));                    % Phase multiplication, implement filtering
clear S4 Hac Hpc Km

%% Transform to the image domain
Stt = ifftshift(ifft(ifftshift(S5, 1), [], 1);                              % Azimuth IFFT of the signal
clear S5

%% Adjust the image
Img = (abs(Stt));                                                           % Flip the image vertically (flipud)
Ntau = round((Tr / dtau / 2 - 1) / 2);                                      % Calculate the length of the range guard zone
Img = Img(1:Na, Ntau + 1:Ntau + Nr - Zr);                                   % Crop the effective region of the image
Img = Img / max(max(Img));
Img = 20 * log10(Img);
Img(Img < -50) = -50;
Img2 = Img(1:797, :);
figure, imagesc(Img2, [-50 0]);
% figure, imagesc(tau * C / 2, eta * V, Img2, [-50 0]); axis image; set(gcf, 'Color', 'w');
colormap('turbo');
colorbar
xlabel("Range (m)", "FontName", "Times New Roman", "FontSize", 13);
ylabel("Azimuth (m)", "FontName", "Times New Roman", "FontSize", 13);
% Img = uint8((Img + 50) / 50 * 255);

% load('radarsat-1_data.mat');
% figure, imagesc(Img, [-50 0]); axis image; set(gcf, 'Color', 'w');
% title('RDA Imaging Result')
%
% load('radarsat-1_csdata.mat');
% figure, imagesc(Img, [-50 0]); axis image; set(gcf, 'Color', 'w');
% title('CS Imaging Result')

% figure, imagesc(Img(688:944, :)); axis image; set(gcf, 'Color', 'w');
% imwrite(Img, 'SAR_CSA.bmp', 'bmp');
