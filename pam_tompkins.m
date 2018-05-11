clear all; 
close all;
clc;

%% Load a signal
[file,path] = uigetfile('*.mat', 'rt');
signal = load(fullfile(path, file));
data = signal.ecg; % Your ecg data
Fs = signal.Fs; % Sampling frequency
N = size(data,2); % Data length

%% Amplifier 


%% Band-pass filter 

% Low-pass filter 

% High-pass filter 

%% Differenciation 

%% Intensification of local etrema

%% Moving Window Integration

%% Threesholding 

%% Detection of maxima 



