%% Delay and Sum Beamformer with a uniform linear array capturing plane wave sound sources
% By Ferdinando Olivieri (f.olivieri@ieee.org)
% https://github.com/F-Olivieri/delay-and-sum-tutorial

% Copyright (c) 2017 Ferdinando Olivieri
%
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% Should you use this code, please send an acknowledgement email to the author.

% Clear Matlab's workspace
close all; clearvars; clc;

%% Customizable parameters

% Uniform Linear Microphone array
ElementSpacing = 0.06; % the spacing between transducers [m]
NumberOfSensors = 5; % number of transducers
N_fft = 2048; % length of the impulse responses

c0 = 343; % Speed of sound (meters/second)

% Sound source
fs = 44100; % Sampling frequency (Hz)
source_direction_deg = 60;%deg % Direction of the plane wave
duration_signal = N_fft + 512; % larger than N_fft
source_signal = randn(duration_signal, 1); % It could be a wavefile too, just make sure you modify the code

% Delay-and-sum Beamformer options
steering_direction_deg = 60; % Steering direction of the beamformer (may be the same as source_direction_deg)

% Plotting options
export_plots_flag = 0; % save plots as specified in graphic_format
% If export_plots_flag = 1 then set also the output format
graphic_format = '-depsc';
PW_frequency_plotting = 1000;%Hz (Monocromatic plane wave, for plotting only)
side = 5;  % Monocromatic plane wave: One-side of the square area in meters
res = 200; % Monocromatic plane wave: Spatial resolution for one-side

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% START PROCESSING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqvect = ((0:N_fft + 1))*fs/N_fft;
lengthfreqvect = length(freqvect);
%% Coordinates of the mic array
[x_c, y_c, ElementNumbering] = get_linear_array_coordinates(NumberOfSensors, ElementSpacing);
plot_linear_array(x_c, y_c, export_plots_flag); % Plot a ULA

%% Plane wave

[Pt, sx, sy] = planewave_soundfield(source_direction_deg, 0, PW_frequency_plotting, c0, side, res, 1);
Pt = real(Pt);

[H_PW, h_IR_PW] = calculate_transfer_function_plane_wave(freqvect, c0, source_direction_deg, N_fft, x_c);

mic_sigs = fftfilt(h_IR_PW, source_signal);
%% plane wave from 0deg and array

%Plot the microphone array to verify the correctness of the geometry
fig1 = figure;
plotPlaneWaveAndMicArray(sx, sy, Pt, x_c,y_c, ElementNumbering, source_direction_deg);
if export_plots_flag, print(fig1,graphic_format,['PW_' num2str(source_direction_deg) 'deg_MicArrayAndPW']); end

fig2 = plotimpulseresponsemics(h_IR_PW);
if export_plots_flag, print(fig2,graphic_format,['PW_' num2str(source_direction_deg) 'deg_MicArrayIR']); end

% %% Summing the mic signals
% h_sum_0deg = sum(h_IR_PW, 2);
% fig5 = figure; plot(h_sum_0deg);
% % xlim([1020 1028]);
% grid on;
% xlabel('Samples'); ylabel('Amplitude');
% if export_plots_flag, print(fig5,graphic_format,'PW_0deg_sumMicArray'); end

%% Delay and Sum beamforming

% Filter Design
W_PW = zeros(lengthfreqvect, NumberOfSensors); % H is a frequency dependent M X N matrix
for m = 1:NumberOfSensors % for each sensor
    for freq_idx = 1:lengthfreqvect % for each angular frequency
        angularfreq = 2*pi*freqvect(freq_idx);
        k_scalar = angularfreq/c0;
        k_vect = k_scalar*sin(deg2rad(steering_direction_deg));%
        W_PW(freq_idx,m) = H_PW(freq_idx, m)*exp(-1i*dot(k_vect, x_c(m)));
    end
end
W_PW = W_PW/NumberOfSensors; % Normalizing the output
w_IR_DS = ifft(W_PW, N_fft, 'symmetric'); % Inverse Fourier transform

beamformer_output = fftfilt(w_IR_DS, mic_sigs);

idx_sample_inf = 1; %N_fft/2 - 20
idx_sample_sup = N_fft; %N_fft/2 + 20
sample_interval = [idx_sample_inf,idx_sample_sup];

[X,Y] = meshgrid(1:NumberOfSensors, 1:N_fft);

fig7 = figure;
plot3(Y, X, circshift(w_IR_DS, N_fft/2), 'b');
xlim(sample_interval)
grid on; xlabel('Samples'); ylabel('Microphone number');
zlabel('Impulse Response Amplitude');
title('Response D&S Beamformer: signal shifting (modeling delay applied)');
view([-20 25]); %view([0 50]);
if export_plots_flag, print(fig7, graphic_format, 'PW_30deg_DS'); end

% fig8 = figure;
% plot(sum(w_IR_DS, 2), 'b');
% grid on; xlabel('Samples'); ylabel('Amplitude');
% xlim([1014 1034]);
% %zlabel('Impulse Response Amplitude');
% title('Response D&S Beamformer: signal summing');
% %view([0 50]);
% if export_plots_flag, print(fig8,graphic_format,'PW_30deg_DS_sum'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x_c, y_c, ElementNumbering] = get_linear_array_coordinates(NumberOfSensors, ElementSpacing)
ElementNumbering = (NumberOfSensors - 1)/2:-1:-(NumberOfSensors - 1)/2;
x_c = ElementSpacing*ElementNumbering; % Microphone positions
y_c = zeros(1, NumberOfSensors);
end

function plot_linear_array(x_c, y_c, export_plots_flag)

figs = figure;
scatter(x_c, y_c,'k'); xlabel('x, [m]'); ylabel('y, [m]'); axis equal; grid on;
ylim([-0.1 0.8]); %xlim([-0.2 0.2]);
yL = get(gca, 'YLim'); yline = line([0 0], yL, 'Color','b'); line2arrow(yline);
xL = get(gca, 'XLim'); xline = line(xL, [0 0], 'Color','b'); line2arrow(xline);
if export_plots_flag, print(figs,'-depsc','PW_0deg_MicArrayCoordinates.eps'); end

end

function [Pt, sx, sy] = planewave_soundfield(AzimuthAngle, ElevationAngle, f, c0, side, res, MeshgridFlag)
if nargin < 7, MeshgridFlag = 0; end

sy = linspace(-side/2, side/2, res + 1);
sx = linspace(-side/2, side/2, res + 1);

%%% Parameters - Frequency %%%
omega = 2*pi*f;       % Angular frequency
k = omega/c0;           % Wave number

%%% Virtual Source Locations
thet_q = ElevationAngle;    % Coelevation angle
phi_q = AzimuthAngle;     % Azimuth angle

% Convert the angle in radians
thet_q = deg2rad(thet_q);
phi_q = deg2rad(phi_q);
% Target fields induced by the Target source

if MeshgridFlag, [sx, sy] = meshgrid(sx, sy); end %%% Meshgrid of the simulated plane for the target fields %%%
sz = zeros(size(sx));      % Meshgrids of the simulated plane

phi_q = phi_q + pi;     % Azimuth angle

% kx=cos(phi_q).*sin(thet_q); %Original
% ky=sin(phi_q).*sin(thet_q);
% kz=cos(thet_q);

kx = sin(phi_q).*cos(thet_q);
ky = cos(phi_q).*cos(thet_q);
kz = sin(thet_q);
Pt = exp(-1i*k*(kx*sx + ky*sy + kz*sz));  % E. Williams Book page 23
end


function plotPlaneWaveAndMicArray(sx, sy, Pt, x_c,y_c, ElementNumbering, theta)

num_sensors = length(x_c);
hold all; pcolor(sx, sy, Pt); shading interp; asse=[-2 2]; caxis(asse); axis square;
colorbar;
%image(x,y,factor*(Z + shift)); colormap jet; shading interp;
hold all; scatter(x_c,y_c,'k'); xlabel('x, [m]'); ylabel('y, [m]'); axis equal; grid on;
for sensor = 1:num_sensors
    hold all; text(x_c(sensor), y_c(sensor) - 0.02,['$m_{' num2str(ElementNumbering(sensor)) '}$']);
end

ylim([-0.1 0.5]); %xlim([-0.2 0.2]);
yL = get(gca,'YLim'); yline = line([0 0],yL,'Color','b'); line2arrow(yline);
xL = get(gca,'XLim'); xline = line(xL,[0 0],'Color','b'); line2arrow(xline);

hold all;
x1 = sind(theta)/5; x2 = 0;
y1 = cosd(theta)/5; y2 = 0;
plot_arrow(x2, y2, x1, y1,'linewidth',.1,'color',[0 0 0],'facecolor',[0 0 0]);
hold all;
%ArcBetweenTwoPoints([x2;y2;0], [x1;y1;0], [0;0;0]);

text('Interpreter','latex',...
    'String','$$\bf{\hat{k}}$$',...
    'Position',[x1 (y1 - 0.03)],...
    'FontSize',16)

end

function [H_PW, h_IR_PW] = calculate_transfer_function_plane_wave(freqvect, c0, directionPW, N_fft, x_c)

num_freq_bins = length(freqvect);
NumberOfSensors = length(x_c);
% TRANSFER FUNCTION H
H_PW = zeros(num_freq_bins, NumberOfSensors); % H is a frequency dependent M X N matrix
for m = 1:NumberOfSensors % for each sensor
    for freq_idx = 1:num_freq_bins % for each angular frequency
        angularfreq = 2*pi*freqvect(freq_idx);
        k_scalar = angularfreq/c0;
        k_vect = -k_scalar*sin(deg2rad(directionPW));
        %exp_modellingdelay = exp(-1i*angularfreq*((N_fft/2 - 1)/fs));
        H_PW(freq_idx,m) = exp(1i*dot(k_vect,x_c(m))); %*exp_modellingdelay
    end
end
h_IR_PW = ifft(H_PW, N_fft, 'symmetric'); % Impulse response matrix

end

function [figs, h] = plotimpulseresponsemics(h)
%% Impulse responses at zero degree

[N_fft, NumberOfSensors] = size(h);

% idx_sample_inf = N_fft/2 - 20;
% idx_sample_sup = N_fft/2 + 20;
idx_sample_inf = 1;
idx_sample_sup = N_fft;
sample_interval = idx_sample_inf:idx_sample_sup;

[X,Y] = meshgrid(1:NumberOfSensors, sample_interval);

figs = figure;
plot3(Y, X, circshift(h( sample_interval, :), N_fft/2), 'b');
grid on; xlabel('Samples'); ylabel('Microphone number');
zlabel('Impulse Response Amplitude');
title('Impulse Responses (modeling delay applied)');
% view([0 50]);
view([-20 25]); %view([0 50]);

end


%% CODE ARCHIVE

% %% Direction 30deg
%
% directionPW = 30;
% [Pt, sx, sy] = planewave_soundfield(directionPW, 0, PW_frequency, c0,side, res, 1);
% Pt = real(Pt);
%
% fig3 = figure;
% plotPlaneWaveAndMicArray(sx, sy, Pt, x_c,y_c, ElementNumbering, k, directionPW);
% if export_plots_flag, print(fig3,graphic_format,['PW_' num2str(directionPW) 'deg_MicArrayAndPW']); end
%
% [fig4, h_30deg, H_30deg] = plotimpulseresponsemics(NumberOfSensors, x_c, directionPW, c0, N_fft, fs, freqvect);
% if export_plots_flag, print(fig4,graphic_format,['PW_' num2str(directionPW) 'deg_MicArrayIR']); end

% h_sum_30deg = sum(h_30deg,2);
% fig6 = figure; plot(h_sum_30deg); xlim([988 1060]); grid on;
% xlabel('Samples'); ylabel('Amplitude');
% if export_plots_flag, print(fig6,graphic_format,'PW_30deg_sumMicArray'); end

