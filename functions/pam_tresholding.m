function [ TRESH1, TRESH2 ] = pam_tresholding( Smwi, data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
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

end

