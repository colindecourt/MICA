clear all; 
close all;
clc;

%% Load a signal


signal = load('ecg_noiseBL.mat');
data = signal.ecg;
Fs = signal.Fs; % Sampling frequency

figure
plot(data);


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


%plot(Y2(1:500)/max(abs(Y2(1:500)))); % normalized signal 


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

% plot(TRESH1(1:1000));
% plot(TRESH2(1:1000));
%plot(TRESH1(1:500)/max(abs(TRESH1(1:500))), 'yellow'); % normalized tresh
%plot(TRESH2(1:500)/max(abs(TRESH2(1:500))), 'green'); % normalize tresh


delay=38 % delay in samples 



%% Detection of maxima 

%Detection of the QRS location
T=zeros(size(Smwi));
for i=1:length(TRESH1)
    if(Smwi(i)>TRESH1(i))
        T(i)=1;
    else
        T(i)=0;
    end
end

index=[0];
T1=diff(T);
for k=1:length(T1)
    if(abs(T1(k))==1)
        index=[index,k];
    end
end
index=index(2:end)-delay; %index of QRS location

%R detection 
R_peak=[0];
R_peak_abs=[0];

for i=1:2:length(index)
    temp=max(data(index(i):index(i+1)));
    R_peak=[R_peak,temp];
end
R_peak=R_peak(2:end);

for j=1:length(R_peak)
    temp1=find(data==R_peak(j));
    R_peak_abs=[R_peak_abs,temp1];
end
R_peak_abs = R_peak_abs(2:end);
hold on
plot(R_peak_abs,R_peak,'*');
    





