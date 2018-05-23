function [nb_ect]= ectopic(data, RR_interval)

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

end