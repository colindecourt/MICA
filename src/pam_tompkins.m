clear all;
close all;
clc;

%% Load a signal


signal = load('ecg_AF.mat');
data = signal.ecg;
Fs = signal.Fs; % Sampling frequency

%figure
%plot(data);


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
figure
plot(Smwi(1:1000));

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
hold on 
plot(TRESH1(1:1000));
plot(TRESH2(1:1000));
%plot(TRESH1(1:500)/max(abs(TRESH1(1:500))), 'yellow'); % normalized tresh
%plot(TRESH2(1:500)/max(abs(TRESH2(1:500))), 'green'); % normalize tresh


delay=38; % delay in samples



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
%R_peak_abs = R_peak_abs(2:end);
%hold on
%plot(R_peak_abs,R_peak,'o');


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

%plot(Q_peak_abs,Q_peak,'o');

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

%plot(S_peak_abs,S_peak,'o');

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

%plot(P_wave_abs,P_wave,'o');
%plot(T_wave_abs,T_wave,'o');

%% Pathologies detection

%Average BPM
RR_mean=mean(RR_interval);
BPM=60./((RR_mean)/200)

% Tachycardia/Bradycardia 

time_total=length(data)/200;
nb_R_peak=length(RR_interval)+1;
incr=floor(10*nb_R_peak/time_total); %determines the step we need to chose in order to have a mean every 10sec

moyenne=0;
moyenne2=0;
for k=1:incr:length(RR_interval)-incr-1
    RR_interval_mean1=mean(RR_interval(k:k+incr));
    moyenne=[moyenne,RR_interval_mean1];
    for i=k:k+incr
        temp=RR_interval_mean1;
        moyenne2(i)=temp; %the vector will be used after
    end
end



moyenne=moyenne(2:end); %Represents the average RR_interval during each period


% BPM for each 
BPM_moy=60./(moyenne/200);

% number of periods during which the patient presents bradycardia or tachycardia 
nb_b1=0;
ind_b1=0;
nb_b2=0;
ind_b2=0;
nb_t1=0;
ind_t1=0;
nb_t2=0;
ind_t2=0;
for i=1:length(BPM_moy)
    if BPM_moy(i)<50
        nb_b2=nb_b2+1;
        ind_b2=[ind_b2,i];
    elseif (BPM_moy(i)>50 && BPM_moy(i)<60) 
        nb_b1=nb_b1+1;
        ind_b1=[ind_b1,i];
    elseif (BPM_moy(i)>100 && BPM_moy(i)<110) 
        nb_t1=nb_t1+1;
        ind_t1=[ind_t1,i];
    elseif BPM_moy(i)>110
        nb_t2=nb_t2+1;
        ind_t2=[ind_t2,i];
    end 
end



% bradycardia rate for the whole signal
brad_rate=(nb_b1+nb_b2)*100/length(moyenne);

% tachycardia rate for the whole signal
tach_rate=(nb_t1+nb_t2)*100/length(moyenne);

% Warning display depending on the pathology (tachycardia od bradycardia)
% detected

if nb_b2/length(moyenne)>0.05
    disp('WARNING : Very low cardiac rhythm at some moments')
end

if (nb_b1+nb_b2)/length(moyenne)>0.5
    disp('General low cardiac rhythm')
end

if nb_t2/length(moyenne)>0.05
    disp('WARNING : Very high cardiac rhythm at some moments')
end

if (nb_t1+nb_t2)/length(moyenne)>0.5
    disp('General high cardiac rhythm')
end 


%Ectopic heart beat 
Ect_diff=abs(diff(RR_interval));
moyenne2=moyenne2(2:end);
l1=length(Ect_diff);
l2=length(moyenne2);
if l2<l1
    for i=1:l1-l2
    moyenne2=[moyenne2, moyenne2(l2)];
    end
end


nb_ect=0;
for i=1:l1
    if Ect_diff(i)>0.5*moyenne2(i) 
        nb_ect=nb_ect+1;
    end
end

%Warning display 
if nb_ect>5
    disp('Ectobeat heartbeat')
end


%Atrial fibrilation 

figure(4) 
h=histogram(RR_interval)

% 
% 
% for i=1:length(RR_interval)-1
%     AF_0(i)=(RR_interval(i)-RR_mean)*(RR_interval(i)-RR_mean);
% end
% 
% for i=1:length(RR_interval)-2
%     AF_1(i)=(RR_interval(i+1)-RR_mean)*(RR_interval(i)-RR_mean);
% end
% 
% R_0=sum(AF_0)/((length(RR_interval)-1)*RR_mean^2);
% R_1=sum(AF_1)/((length(RR_interval)-1)*RR_mean^2);
% 
% 
% 
