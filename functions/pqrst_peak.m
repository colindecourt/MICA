function [ P_wave, P_wave_abs, Q_peak, Q_peak_abs, R_peak, R_peak_abs, S_peak, S_peak_abs, T_wave, T_wave_abs ] = pqrst_peak( Smwi, TRESH1, data)
%PQRST_PEAK Summary of this function goes here
%   Detailed explanation goes here
%% Detection of maxima
delay =38; 
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
    if (index(i)<0)
        index(i)=abs(index(i));
    end
    
    temp=max(data(index(i):index(i+1)));
    R_peak=[R_peak,temp];
    
    
end


R_peak=R_peak(2:end);

for j=1:length(R_peak)
    temp1=find(data==R_peak(j));
    R_peak_abs=[R_peak_abs,temp1];
end
R_peak_abs = R_peak_abs(2:end);

%Q detection
Q_peak_abs=[0];
c=R_peak_abs(1);
for i=2:length(R_peak_abs)
    while(data(c)>data(c-1))
        c=c-1;     
    end
    Q_peak_abs=[Q_peak_abs,c];
    c=R_peak_abs(i);
end
Q_peak_abs=Q_peak_abs(2:end);

Q_peak=data(Q_peak_abs);

%S detection
S_peak_abs=[0];

d=R_peak_abs(1);
for i=2:length(R_peak_abs)
    while(data(d)>data(d+1))
        d=d+1;     
    end
    S_peak_abs=[S_peak_abs,d];
    d=R_peak_abs(i);
end
S_peak_abs=S_peak_abs(2:end);

S_peak=data(S_peak_abs);

%% P and T wave detection 
b5=[1 0 0 0 0 0 -1];
a5=[1];

G1=filter(b5,a5,data); %differenciator 

b6=[1 0 0 0 0 0 0 0 -1];
a6=[1 -1];

G2=filter(b6,a6,G1);


delay2 = 6; %the delay is seen in the filtered signal each time the decreasing slope after the R peak crosses the level 0

index_RR_interval=R_peak_abs+6;
RR_interval=diff(index_RR_interval); %length of each RR interval

% T waves detection (we omit the first T wave)
T_wave_f=0;
for k=1:length(index_RR_interval)-1
    T_wave_f=[T_wave_f, max(((G2((index_RR_interval(k)+floor(0.1*RR_interval(k))):(index_RR_interval(k)+floor(0.65*RR_interval(k)))))))];
end

T_wave_f=T_wave_f(2:end);
T_wave_abs=0;
for j=1:length(T_wave_f)
    temp1=find(G2==T_wave_f(j));
    T_wave_abs=[T_wave_abs,temp1];
end

T_wave_abs=T_wave_abs(2:end);


% Detection of P waves
P_wave_f=0;
for k=1:length(index_RR_interval)-1
    P_wave_f=[P_wave_f, max(((G2((index_RR_interval(k)+floor(0.71*RR_interval(k))):(index_RR_interval(k)+floor(0.90*RR_interval(k)))))))];
end

P_wave_f=P_wave_f(2:end);
P_wave_abs=0;
for j=1:length(P_wave_f)
    temp1=find(G2==P_wave_f(j));
    P_wave_abs=[P_wave_abs,temp1];
end

P_wave_abs=P_wave_abs(2:end);

P_wave = data(P_wave_abs);
T_wave = data(T_wave_abs);


end

