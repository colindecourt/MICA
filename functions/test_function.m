clear all; 
close all;
clc;


%% Load a signal


signal = load('ecg_normal_1.mat');
data = signal.ecg;
Fs = signal.Fs; % Sampling frequency

%% Test functions

Smwi = pam_filter(data,Fs);
[TRESH1, TRESH2] = pam_tresholding(Smwi, data);
[ P_wave, P_wave_abs, Q_peak, Q_peak_abs, R_peak, R_peak_abs, S_peak, S_peak_abs, T_wave, T_wave_abs ] = pqrst_peak( Smwi, TRESH1, data);

figure 
plot(data/max(abs(data)));
hold on
plot(TRESH1/(max(abs(TRESH1))));
hold off