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

Y2=filter(b2,a2,Y1); % signal filtered 

figure(1)
hold on
plot(Y2(1:500))



%% Differenciation 

b3=Fs/8.*[1 2 0 -2 1];
a3=[1]; % the delay of 2 samples is for now ignored 

Y3=filter(b3, a3, Y2);

%% Intensification of local etrema

Ssq = abs(Y3).^2;

%% Moving Window Integration

N=30; % average QRS window length

b4=1/N.*ones(1,N);
a4=[1]
Smwi=filter(b4, a4, Ssq);


plot((1/10^8)*Smwi(1:500));


%% Threesholding

PEAK=max(Smwi(1));

for i=1:200:length(data)-200
    [R ,x] = max(Smwi(i:i+200)); % Search the max during a period of 200 samples
    m = mean(Smwi(i:i+200)); %Search the mean to don't detect noise peak (as P wave for example)
    for k=i:i+200
     A(k) = 0.875*R+0.125*PEAK; 
     B(k) = x;
     M(k) = 0.875*m+0.125*PEAK;
    end
end


TRESH1 = M+0.25*(A-M);
TRESH2 = 0.5*TRESH1;


% plot(A(50780:51999));
% plot(M(50780:51999));
hold off;
%% Detection of maxima 



