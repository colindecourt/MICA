function [nb_b1, nb_b2, nb_t1, nb_t2, ind_b1, ind_b2, ind_t1, ind_t2, brad_rate, tach_rate]=brad_tach(data, RR_interval)
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
end 