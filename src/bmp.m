function [BPM] = bpm(RR_interval)
RR_mean=mean(RR_interval);
BPM=60./((RR_mean)/200)
end 