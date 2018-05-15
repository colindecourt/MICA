clear all; 
close all;
clc;

%% Load a signal


signal = load('ecg_normal_1.mat');
data = signal.ecg;
Fs = signal.Fs; % Sampling frequency

%figure(4)
%plot(data(1:500));


%% Band-pass filter 

% Low-pass filter 

b1=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a1=[1 -2 1];

Y1 = filter(b1, a1, data); 
% fvtool(b1, a1); %enables to see the delay introduced by the first filter 


% High-pass filter 

b2=zeros(1, 33);
b2(1)=-1;
b2(17)=32;
b2(18)=-32;
b2(33)=1;
a2=[1 -1];

Y2=filter(b2,a2,Y1); % signal filtered 
% fvtool(b2, a2); %enables to see the delay introduced by the second filter



% plot(Y2(1:500));
figure(1)
hold on
plot(Y2(1:500)/max(abs(Y2(1:500)))); % normalized signal 


%% Differenciation 

b3=(Fs/8).*[1 2 0 -2 -1];
a3=[1]; % the delay of 2 samples is for now ignored 

Y3=filter(b3, a3, Y2);

%% Intensification of local etrema

Ssq = abs(Y3).*abs(Y3);

%% Moving Window Integration

N=30; % average QRS window length

b4=ones(1,N);
%a4=[1]
Smwi=1/N*conv(b4,Ssq);


plot(Smwi(1:1000));

%plot(Smwi(17:517)/max(abs(Smwi(17:517)))); % normalized signal


%% Threesholding


PEAK=max(Smwi(1));

for i=1:400:length(data)-400
    R = max(Smwi(i:i+400)); % Search the max during a period of 200 samples
    m = mean(Smwi(i:i+400)); %Search the mean to don't detect noise peak (as P wave for example)
    for k=i:i+400
     A(k) = 0.875*R+0.125*PEAK; 
     
     M(k) = 0.875*m+0.125*PEAK;
    end
end


TRESH1 = M+0.25*(A-M);
TRESH2 = 0.5*TRESH1;

plot(TRESH1(1:1000));
plot(TRESH2(1:1000));
%plot(TRESH1(1:500)/max(abs(TRESH1(1:500))), 'yellow'); % normalized tresh
%plot(TRESH2(1:500)/max(abs(TRESH2(1:500))), 'green'); % normalize tresh
hold off;


delay=38 % delay in samples 



%% Detection of maxima 

X=[0];
for i=1:length(TRESH1)
    if (0.85*TRESH1(i) < Smwi(i) && Smwi(i)< 1.15*TRESH1(i))
        X=[X,i];
        
    end
end





