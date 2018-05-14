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


% subplot(2,2,1);
% plot(data(1:500), 'blue');
% 
% subplot(2,2,2);
% plot(Y1(1:500), 'red');
% 
% subplot(2,2,3);
figure(1)
plot(Y2(1:1000))



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

figure(2)
hold on
plot(Smwi(1:500));

%% Threesholding

% M�thode 1 de ne marche pas car seuil adaptafif seulement jusqu'� trouver
% le premier max %

% SPKI=zeros(1,length(data));
% NPKI=zeros(1,length(data));
% SPKI(1)=max(Smwi(1:2*Fs));
% NPKI(1)=mean(Smwi(1:2*Fs));
% PEAK=max(Smwi(1));
% for j=1:length(data)
%     if(Smwi(j)>PEAK)
%         PEAK=Smwi(j);
%     end
%    SPKI(j+1)=0.125*PEAK+0.875*SPKI(j);
%    NPKI(j+1)=0.125*PEAK+0.875*NPKI(j);
% end
% 
% TRESH1 = NPKI+0.25*(SPKI-NPKI);
% 
% TRESH2 = 0.5*TRESH1;

% M�thode 2 %

for i=1:200
    [R , x] = max(Smwi(1:200));
    m = mean(Smwi(1:200));
end

plot(x,R,'*');
plot(m,'x');
hold off;
%% Detection of maxima 



