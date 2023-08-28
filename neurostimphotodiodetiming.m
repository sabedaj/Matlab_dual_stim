%find timing of photodiode neurostim
FS=30000;%sample rate
num_channels = length(1); % ADC input info from header file
fileinfo = dir('analogin.dat');
if ~isempty(fileinfo)
num_samples = fileinfo.bytes/(num_channels * 2); % uint16 = 2 bytes
fid = fopen('analogin.dat', 'r');
v = fread(fid, [num_channels, num_samples], 'uint16');
end
timingSignal=0:1/FS:length(v)/FS-1/FS;
maxVal=max(v);
sampletimemax=find(v==maxVal);
signaldigitised=zeros(length(v),1);
signaldigitised(v==maxVal)=1;
zero_period_length = 10000;
zero_periods = [];
current_zero_count = 0;
for i = 1:length(signaldigitised)
if signaldigitised(i) == 0
current_zero_count = current_zero_count + 1;
if current_zero_count == zero_period_length
zero_periods = [zero_periods, i - zero_period_length + 1];
end
else
current_zero_count = 0;
end
end
onsetstimulus=length(signaldigitised)-zero_periods(2:end);%calculates stimulus onset sample
figure; plot(signaldigitised); hold on
scatter(onsetstimulus,ones(length(onsetstimulus),1),'r')