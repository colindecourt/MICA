function [ Smwi] = pam_filter(data, Fs)
%PAM_FILTER Summary of this function goes here
%   Detailed explanation goes here
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

%% Differenciation

b3=(Fs/8).*[1 2 0 -2 -1];
a3=[1]; % the delay of 2 samples is for now ignored

Y3=filter(b3, a3, Y2);

%% Intensification of local etrema

Ssq = abs(Y3).*abs(Y3);

%% Moving Window Integration

N=30; % average QRS window length

b4=ones(1,N);
Smwi=1/N*conv(b4,Ssq);



end

