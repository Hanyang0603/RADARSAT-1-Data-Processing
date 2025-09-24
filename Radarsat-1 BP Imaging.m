%% Processing Radarsat-1 Data Using BP Algorithm %%
% Author: Hanyang Xu
% The last modified date: 2022/09/10

clear all;
close all;
clc;

select = 0;                                                                   % 1 for non-uniform resolution, 0 for uniform resolution
PLO = 0;                                                                      % Whether to plot

%% Extract Raw Data
extract = 0;                                                                  % Whether to extract new raw data
if extract
    run specify_parameters.m
    run extract_data.m
end

first_r = 760;                                                                % Position of the first range sample
load CDdata1.mat                                                              % Load raw data
data = double(data.');                                                        % Convert data format and transpose
[Nr, Na] = size(data);

%% Input Parameters
C = 2.9979e8;                                                                 % Speed of light
f0 = 5.3e9;                                                                   % Center frequency
lambda = C / f0;                                                              % Wavelength
clear f0
V = 7062;                                                                     % Platform velocity
H = 798505;                                                                   % Platform altitude

% Range Parameters
Kr = -0.72135e12;                                                             % Chirp rate
Tr = 41.75e-6;                                                                % Pulse duration
Fr = 32.317e6;                                                                % Range sampling rate
dtau = 1 / Fr;                                                                % Range sampling interval
t0 = 6.5956e-3;                                                               % Data acquisition time, start sampling time
t0 = t0 + (first_r - 1) * dtau;                                               % Update t0, actual start time
W = (Nr - 1) * C * dtau / 2;                                                  % Swath width corresponding to the data
Rnear = t0 * C / 2;                                                           % Nearest slant range calculated from the actual start sampling time
Rfar = Rnear + W;
Rref = (Rnear + Rfar) / 2;                                                    % Reference slant range calculated from the actual start sampling time
theta_near = acosd(H / Rnear);                                                % Near-range incidence angle
theta_far = acosd(H / Rfar);                                                  % Far-range incidence angle
theta_mid = acosd(H / Rref);                                                  % Mid-range incidence angle, approximately equal to the average of the previous two

% Azimuth Parameters
Fa = 1256.98;                                                                 % Azimuth sampling rate, PRF
fdc = -7045.1;                                                                % Doppler centroid frequency, estimated Doppler center frequency
deta = 1 / Fa;                                                                % Azimuth sampling interval, PRT
fd_min = -7557.91;                                                            % Minimum Doppler frequency
fd_max = -6532.3;                                                             % Maximum Doppler frequency
theta_sq_min = -asind(fd_max * lambda / 2 / V);                               % Minimum squint angle determined by the maximum Doppler frequency
theta_sq_max = -asind(fd_min * lambda / 2 / V);                               % Maximum squint angle determined by the minimum Doppler frequency
theta_sq_c = (theta_sq_min + theta_sq_max) / 2;                               % Center squint angle
beta = theta_sq_max - theta_sq_min;                                           % Azimuth beamwidth (3dB)
clear first_r theta_sq_max theta_sq_min fd_min fd_max

%% Time and Frequency Axis Setup
tau = t0 + (0:Nr-1) * dtau;                                                   % Range time axis
clear t0
eta = ((0:Na-1) - Na/2) * deta;                                               % Azimuth time axis
feta = fdc + ((0:Na-1) - Na/2) / Na * Fa;                                     % Azimuth frequency axis

%% Projection Plane Coordinate Axis Setup
sarposi_x = eta .* V * cosd(theta_sq_c);                                      % SAR position in the projection plane coordinate system (x-coordinate)
sarposi_y = -eta .* V * sind(theta_sq_c);                                     % SAR position in the projection plane coordinate system (y-coordinate)
Xmin = min(sarposi_x) + Rref * tand(beta / 2);
Xmax = max(sarposi_x) - Rref * tand(beta / 2);
Lsar = 2 * Rref * tand(beta / 2);                                             % Synthetic aperture length

%% Doppler Centroid Estimation
%{
% Since the spectrum amplitude is symmetric about the peak, the Doppler centroid can be simply determined by energy balancing;
% That is, find the frequency point that divides the spectrum energy equally;
% Energy balancing can be achieved by convolving the received average power spectrum with an energy balancing filter.
Srnm_Fr = fftshift(fft(fftshift(data, 2), [], 2);                             % Azimuth FFT
S_aps = mean(abs(Srnm_Fr).^2, 1);                                             % Azimuth average power spectrum
% Filter = [ones(1, Na/2), -1 .* ones(1, Na/2)];                               % Energy balancing filter
% S_aps1 = conv(S_aps, Filter);
% S_aps1 = S_aps1(1, Na/2+1:Na/2+Na);
% S_aps1 = cconv(S_aps, Filter, 2*Na-1);

% According to the Doppler model, find the true Doppler frequency spectrum using fractional PRF
Srnm_Fr2 = fftshift(Srnm_Fr, 2);
S_aps2 = mean(abs(Srnm_Fr2).^2, 1);
S_aps3 = 10 .* log10(S_aps2 ./ max(S_aps2));
M_amb = -6;
f = linspace(M_amb * Fa - Fa/2, M_amb * Fa + Fa/2, Na);

figure;
subplot(3,1,1)
plot((abs(Srnm_Fr(200, :))).^2);
grid on
xlabel("Azimuth frequency (cells)", "FontName", "Times New Roman", "FontSize", 11);
ylabel("Power", "FontName", "Times New Roman", "FontSize", 11);
title("Spectrum of one line of received signal", "FontName", "Times New Roman", "FontSize", 11);

subplot(3,1,2)
plot(S_aps);
grid on
xlabel("Azimuth frequency (cells)", "FontName", "Times New Roman", "FontSize", 11);
ylabel("Power", "FontName", "Times New Roman", "FontSize", 11);
title("Observed average spectrum", "FontName", "Times New Roman", "FontSize", 11);

subplot(3,1,3)
plot(f, S_aps3);
axis([M_amb * Fa - Fa/2, M_amb * Fa + Fa/2, -8, 0]);
hold on
plot([M_amb * Fa - Fa/2, M_amb * Fa + Fa/2], [-6, -6], 'r');
grid on
xlabel("Azimuth frequency (Hz)", "FontName", "Times New Roman", "FontSize", 11);
ylabel("Power", "FontName", "Times New Roman", "FontSize", 11);
title("Unaliased average spectrum", "FontName", "Times New Roman", "FontSize", 11);

S_aps3 = 10 .* log10(S_aps ./ max(S_aps));
f = linspace(M_amb * Fa - Fa/2, M_amb * Fa + Fa/2, Na);
plot(S_aps3);
grid on
xlabel("Azimuth frequency (cells)", "FontName", "Times New Roman", "FontSize", 11);
ylabel("Amplitude", "FontName", "Times New Roman", "FontSize", 11);
title("Output of power balancing filter", "FontName", "Times New Roman", "FontSize", 11);
%}

%% Range Spectrum Windowing and Pulse Compression
Srnm_Fr_1 = fftshift(fft(fftshift(data, 1), [], 1);                           % Range FFT of raw data
WindowR = kaiser(Nr, 2.4) * ones(1, Na);                                      % Range window, kaiser window
Srnm_Fr_2 = Srnm_Fr_1 .* WindowR;
data = fftshift(ifft(fftshift(Srnm_Fr_2, 1), [], 1);                           % Range IFFT

clear WindowR Srnm_Fr_2 Srnm_Fr_1

% Frequency domain filtering method 1: Time-reversed conjugate of the replica pulse, zero-padded DFT
Srnm = [data; zeros(ceil(Tr / dtau), Na)];                                    % Zero-padding
Srnm_Fr = fftshift(fft(fftshift(Srnm, 1), [], 1);                             % Range FFT
[Nr1, ~] = size(Srnm_Fr);                                                     % Nr1 = Nr + ceil(Tr / dtau)
tt0 = linspace(-Tr/2, Tr/2 - Tr/ceil(Tr * Fr), ceil(Tr * Fr));                % Time axis of the original pulse in the time domain
temp_x = exp(1i * pi * Kr * tt0.^2);                                          % Sampled original pulse
temp_x = [conj(temp_x), zeros(1, Nr1 - ceil(Tr * Fr))];                       % Zero-padding to the same length
tempF = ifftshift(fft(fftshift(temp_x)));                                     % ifftshift shifts the last element to the zero-frequency position
pcCoe = repmat(tempF.', [1, Na]);
Srnm_Fr = Srnm_Fr .* pcCoe;
Srnm_Tr = fftshift(ifft(fftshift(Srnm_Fr, 1), [], 1), 1);                     % Range IFFT, compressed signal in the time domain
Srnm_Tr = [Srnm_Tr(end - Nr/2 + 1:end, :); Srnm_Tr(1:Nr/2, :)];               % Compressed signal in the time domain

if PLO
    plotImageGray(Srnm_Tr, 'After range compression');
end
clear Srnm_Fr pcCoe tempF temp tt0 Nr1 Srnm Kr Tr data

if select                                                                     % If non-uniform
    data = data .* exp(1i * 2 * pi * fdc * (eta .* ones(Nr, 1)));             % Shift to Doppler centroid
end

%% Subregion and Resolution Setup
if select
    S_x = Xmax - Xmin;                                                        % Azimuth swath width
    M = 8;                                                                    % Number of subregions, coarse division
    step = S_x / M;                                                           % Azimuth width of one subregion under coarse division
    clear S
    for i = 1:M-1                                                             % First coarse division, then fine division for the last part
        preset2(i, :) = [Xmin + (i-1) * step, Xmin + i * step];               % Azimuth coordinate range of each subregion
    end
    temp_x = max(preset2(:, 2));                                              % Update the maximum boundary after coarse division
    S_x = Xmax - temp_x;                                                      % Remaining swath width
    M = 8;                                                                    % Fine division number
    step = S_x / M;                                                           % Azimuth width of one fine subregion
    for i = 1:M
        preset2 = [preset2; temp_x + (i-1) * step, temp_x + i * step];
    end
    % preset2 = [xmin, xmax, xcenter]
    preset2 = [preset2, (preset2(:, 1) + preset2(:, 2)) ./ 2];                % Minimum, maximum, and center values of the azimuth coordinates of a subregion
    M = size(preset2, 1);                                                     % Total number of subregions after coarse and fine division
    clear step S
    
    %% Parameters Required for Each Subregion Resolution
    % Assume the platform flies from top to bottom
    Lsar = [preset2(:, 2) + r1, preset2(:, 1) + r2];                         % Actual platform position, from the start of the first beam to the end of the last beam
    % SAR swath does not exceed the azimuth sampling position
    Lsar(:, 1) = max(Lsar(:, 1), min(eta * V));                               % SAR start point
    Lsar(:, 2) = min(Lsar(:, 2), max(eta * V));                               % SAR end point
    
    % Calculate the highest achievable resolution for the current subregion
    % Use the center of the subregion as the reference point to calculate the resolution accumulation angle
    theta1 = asin((Lsar(:, 1) - preset2(:, 3)) ./ Rref);                      % Smaller angle
    theta2 = asin((Lsar(:, 2) - preset2(:, 3)) ./ Rref);                      % Larger angle
    Lreso = lambda ./ 4 ./ sin(abs((theta2 - theta1)) ./ 2);                   % Highest resolution limit
    Lreso(1:6) = 20;
    Lreso(3:4) = 12;
    reso = Lreso;                                                             % Actual resolution requirement
    % clear Lreso
    delta = 2 .* asin(lambda ./ 4 ./ reso);                                   % rad, calculate accumulation angle
    
    Lsar_new = [Lsar(:, 1), preset2(:, 3) + Rref * sin(theta1 + delta)];      % Original SAR start point remains unchanged, only the end point changes
    clear delta theta1 theta2 reso
    
    Numposi = round(Lsar_new ./ V ./ deta) + Na / 2 + 1;                      % Corresponding azimuth sample index range
    
    %%%%%%%%%%%%%%%%%%%%%% Update Global SAR Index %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do not consider the case where the start point of the next subregion exceeds the start point of the previous subregion, default to step distribution
    [temp_x, b] = sort(Numposi(:, 1));                                        % temp is the azimuth sample start index of each subregion, b is the corresponding subregion index
    td = Numposi(b(1), 1):1:Numposi(b(1), 2);                                % Azimuth transmission time index, first use subregion 1's sample index as the initial value
    for i = 2:M                                                              % Traverse each subregion
        if td(end) < temp_x(i)                                                % If the start point of the next subregion is greater than the current SAR end point, i.e., no overlap
            td = [td, Numposi(b(i), 1):1:Numposi(b(i), 2)];
        else                                                                  % SAR overlap
            td = [td, td(end) + 1:1:max(td(end), Numposi(b(i), 2))];
        end
    end
    clear b temp
    sarposi = eta(td) * V;                                                    % Actual azimuth position axis
    % td = (min(Numposi) - 1 - Na / 2) * PRT:PRT:(max(Numposi) - 1 - Na / 2) * PRT; % Azimuth transmission (sampling) time axis
    tdd = td2tdd2(td, Numposi, M);                                            % Transmission matrix
    clear M
    
    %% Subaperture Filtering
    % First, do not change PRT, only select the corresponding illumination points and illumination range
    Srnm_rd = fftshift(fft(fftshift(Srnm_Tr, 2), [], 2);                      % Transform to range-Doppler domain
    
    S7 = [];
    df1 = min(feta);
    df2 = max(feta);
    df3 = max(feta) - min(feta);
    for i = 1:length(tdd)                                                     % Traverse within the new transmission axis range
        temp2 = find(tdd(:, i) == 1);                                         % Find the index of the required illumination region
        WindowA = zeros(Nr, Na);
        theta1 = asind((sarposi(i) - preset2(temp2, 1)) ./ Rref);             % Minimum boundary angle with respect to the velocity direction
        fd1 = 2 * V / lambda .* cosd(theta1 + 90);                            % More negative
        theta2 = asind((sarposi(i) - preset2(temp2, 2)) ./ Rref);             % Maximum boundary angle with respect to the velocity direction
        fd2 = 2 * V / lambda .* cosd(theta2 + 90);                            % Larger
        rank1 = floor((max(fd1, df1) + abs(df1)) ./ df3 .* Na) + 1;           % Position of the smaller Doppler frequency
        rank2 = floor((min(fd2, df2) + abs(df1)) ./ df3 .* Na);               % Position of the larger Doppler frequency
        for k = 1:length(temp2)
            WindowA(:, rank1(k):rank2(k)) = 1;                                % Rectangular window
        end
        S6 = Srnm_rd .* WindowA;                                              % Filter according to the current transmission position
        Stt = ifftshift(ifft(ifftshift(S6, 2), [], 2);                        % Azimuth IFFT
        Stt = Stt .* exp(-1i * 2 * pi * fdc * (eta .* ones(Nr, 1)));
        S7 = [S7, Stt(:, td(i))];
    end
    
    clear Srnm_rd td WindowA df1 df2 df3
end

%% Backprojection
tic
if select == 1                                                                % Non-uniform case
    Na = length(tdd);
    Naa = 2000;
    Nup1 = 10;                                                                % Range frequency domain zero-padding factor
    Nr_up = Nr * Nup1;
    Srnm_Fr2 = fft(S7, [], 1);                                                % Range FFT
    Srnm_up = [Srnm_Fr2(1:Nr/2,:);zeros(Nr_up - Nr,Na);Srnm_Fr2(Nr/2+1:end,:)];
    Srnm_Tr_up = ifft(Srnm_up, [], 1);
    dTr = dtau / Nup1;                                                        % Frequency domain interpolation = time domain interpolation
    clear Nup1
    
    ran = tau * C / 2;
    clear tau eta
    Srnm_TaTr = zeros(Nr, Naa);
    td3 = linspace(Xmin, Xmax, Naa);
    for i = 1:Na
        temp_x = (sarposi(i) - td3); % Difference between each position in the scene and the current transmission point
        iran = sqrt(repmat(ran.', [1, Naa]).^2 + repmat(temp_x.^2, [Nr, 1])); % Relative range
        t_ij = 2 * iran / C; % Corresponding time delay for each position in the scene, matrix Nr * Na * Nup2, approximately covering the range scope
        t_ij = round((t_ij - (2 * Rnear / C)) / dTr); % Calculate the index in the current signal column based on the time delay of each position
        it_ij = (t_ij > 0 & t_ij <= Nr_up); % Output 1 if within the range scope, otherwise 0
        t_ij = t_ij .* it_ij + Nr_up * (1 - it_ij); % Replace zeros in the matrix with Nr_up, all unsuitable elements are normalized to the last row
        clear it_ij
        % The matrix elements are now only t_ij and Nr_up
        % During projection, only consider projections within the imaging scope, i.e., the scope of Na
        sig_rdta = Srnm_Tr_up(:, i); % The i-th column of the time-domain signal after interpolation
        sig_rdta(Nr_up) = 0; % Set the last row to zero, column vector
        clear temp
        illu = find(tdd(:, i) == 1); % Find the region illuminated by the current azimuth
        td2 = zeros(1, Naa);
        for k = 1:length(illu) % Traverse all regions that need to be illuminated
            td2 = td2 + (preset2(illu(k), 2) >= td3 & td3 >= (preset2(illu(k), 1))); % Find the imaging scope illuminated by the current azimuth
        end
        clear illu
        td2 = ones(Nr, 1) .* td2; % Expand td2 to Nr * Na * Nup matrix
        Srnm_TaTr = Srnm_TaTr + td2 .* sig_rdta(t_ij) .* exp(1i * 4 * pi / lambda .* iran); % Indexing, exp(1i * 4 * pi / lambda .* iran) performs phase compensation, making each phase consistent
        % If sig_rdta(t_ij) indexes a position outside the imaging scope, set it to zero
    end
    
else                                                                        % Uniform resolution
    Nup1 = 10;                                                              % Range interpolation factor
    Nr_up = Nr * Nup1;
    Srnm_Fr2 = fft(Srnm_Tr, [], 1);
    Srnm_up = [Srnm_Fr2(1:Nr/2, :); zeros(Nr_up - Nr, Na); Srnm_Fr2(Nr/2+1:end, :)]; % Middle interpolation, frequency domain interpolation = time domain interpolation, equivalent to shortening the sampling period
    Srnm_Tr_up = ifft(Srnm_up, [], 1);                                        % Range IFFT
    dTr = dtau / Nup1;                                                        % Update sampling interval
    clear Nup1 Srnm_Fr2 Srnm_up dtau
    
    ran = tau * C / 2;                                                          % Slant range scope
    clear tau eta
    Srnm_TaTr = zeros(Nr, Na);                                               % Imaging matrix
    td3 = linspace(Xmin, Xmax, Na);                                             % Imaging scene scope azimuth axis
    for i = 1:Na
        temp_x = (sarposi_x(i) - td3);                                        % In the projection plane coordinate system, the x difference between each resolution unit in the scene and the current transmission point
        temp_y = (sarposi_y(i) - ran);                                          % In the projection plane coordinate system, the y difference between each resolution unit in the scene and the current transmission point
        td2 = repmat((temp_x <= Rref * tand(beta / 2)) & (temp_x >= -Rref * tand(beta / 2)), [Nr, 1]); % Original td2
        iran = sqrt(repmat((temp_y.^2).', [1, Na]) + repmat(temp_x.^2, [Nr, 1])); % Current slant range to each resolution unit in the scene, matrix [Nr, Na]
        t_ij = 2 * iran / C;                                                    % Corresponding time delay for each position in the scene, matrix Nr * Na * Nup2, approximately covering the range scope
        t_ij = round((t_ij - (2 * Rnear / C)) / dTr);                             % Calculate the index in the current signal column based on the time delay of each position
        it_ij = (t_ij > 0 & t_ij <= Nr_up);                                         % Output 1 if within the range scope, otherwise 0
        t_ij = t_ij .* it_ij + Nr_up * (1 - it_ij);                                   % Replace zeros in the matrix with Nr_up, all unsuitable elements are normalized to the last row
        clear it_ij
        % The matrix elements are now only t_ij and Nr_up
        % During projection, only consider projections within the imaging scope, i.e., the scope of Na
        sig_rdta = Srnm_Tr_up(:, i);                                           % The i-th column of the time-domain signal after interpolation
        sig_rdta(Nr_up) = 0;                                                  % Set the last row to zero, column vector
        Srnm_TaTr = Srnm_TaTr + td2 .* sig_rdta(t_ij) .* exp(1i * 4 * pi / lambda .* iran); % Indexing, exp(1i * 4 * pi / lambda .* iran) performs phase compensation, making each phase consistent
    end
    
    Img = abs(Srnm_TaTr).';
    Img = flipud(Img);
    Img = Img / max(max(Img));
    Img = 20 * log10(Img);
    Img(Img > -0) = 0;
    Img(Img < -50) = -50;
    td3 = td3 - min(td3);
    % td3 = fliplr(td3);
    
    if PLO
        figure;
        imagesc(ran, td3, (Img), [-50, 0]); %
        colormap('turbo');
        colorbar
        xlabel("Range (m)", "FontName", "Times New Roman", "FontSize", 13);
        ylabel("Azimuth (m)", "FontName", "Times New Roman", "FontSize", 13);
        axis xy;
        % title("Spectrum of one line of received signal", "FontName", "Times New Roman", "FontSize", 11);
    end
    
end
toc

Img = abs(Srnm_TaTr).';
Img2 = Img(1:500, :);
Img4 = Img(1001:1501, :);
Img6 = [Img2; Img3; Img4];
Img2 = Img2 / max(max(Img6));
Img2 = 20 * log10(Img2);
% Img2(Img2 < -50) = -50;

Img3 = Img(501:1000, :);
Img3 = Img3 / max(max(Img3));
Img3 = 20 * log10(Img3);
% Img3(Img3 < -50) = -50;

Img4 = Img4 / max(max(Img6));
Img4 = 20 * log10(Img4);
% Img4(Img4 < -50) = -50;

Img5 = Img(1502:end, :);
Img5 = Img5 / max(max(Img));
Img5 = 20 * log10(Img5);
% Img5(Img5 < -50) = -50;

Img = [Img2; Img3; Img4; Img5];
Img(Img > -5) = 0;
Img(Img < -40) = -40;
figure;
imagesc((Img), [-40, 0]);
ylabel('Azimuth (m)');
xlabel('Range (m)'); axis xy;
