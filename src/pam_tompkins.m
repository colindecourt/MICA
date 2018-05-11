clear all; 
close all;
clc;

%% Load a signal


signal = load('ecg_normal_1.mat');
data = signal.ecg;
Fs = signal.Fs; % Sampling frequency



%% Band-pass filter 

% Low-pass filter 

b1=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a1=[1 -2 1];

Y1 = filter(b1, a1, data); 

% High-pass filter 

b2=zeros(1, 33);
b2(1)=-1;
b2(17)=32;
b2(18)=-32;
b2(33)=1;
a2=[1 -1];

Y2=filter(b2,a2,Y1); % combination of the two filters resulting in the band-pass filter 


% subplot(2,2,1);
% plot(data(1:500), 'blue');
% 
% subplot(2,2,2);
% plot(Y1(1:500), 'red');
% 
% subplot(2,2,3);
% plot(Y2(1:500), 'green');



%% Differenciation 

b3=Fs/8.*[1 2 0 -2 1];
a3=[1]; % the delay of 2 samples is for now ignored 

Y3=filter(b3, a3, Y2);

% subplot(2,2,4);
% plot(Y3(1:500));


%% Intensification of local etrema

Ssq = abs(Y3).^2;

% figure(2) 
% plot(Ssq(1:500));

%% Moving Window Integration

N=30; % average QRS window length


b4=1/N.*ones(1,N);
a4=[1]

Smwi=filter(b4, a4, Ssq);

plot(Smwi(1:500));

%% Threesholding 




%% Detection of maxima 



