function plotTrialData(dName,chn,trialID,repeatNum)
%dName is the file name e.g. amplifier or xxxx.mu or xxxx_dn 
    %amplifier is raw, xxxx.mu if filtered
%chn is channel from original order e.g. 1==32
%trialID is the ID of the condition tested - can use loadtrialinfo
%repeatNum is the repeat of the condition
FS=30000;
nChn=128;
shortbytes=2;
TimeBeforeTrig=0.05;%time in ms
TimeAfterTrig=0.1;%time in ms
filepath=pwd;
vFID = fopen([filepath filesep dName '.dat'],'r');
TP = loadTrialParams;
trialinfo=loadTrialInfo;
trialinfo=cell2mat(trialinfo(2:end,[1 2 4:end]));
numsim_elect=length(trialinfo(:,2))/length(unique(trialinfo(:,1)));
tID = find(cell2mat(TP(numsim_elect:numsim_elect:end,2)) == trialID);
trig=loadTrig(0);
offset=trig(tID(repeatNum))*nChn*shortbytes-TimeBeforeTrig*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
fseek(vFID,offset,'bof');
v = fread(vFID,[nChn, (TimeBeforeTrig+TimeAfterTrig)*FS],'short') ./ 10;
figure
plot(-TimeBeforeTrig:1/FS:TimeAfterTrig-1/FS,v(chn,:))
xlabel('Time(s)')
ylabel('Amplitude (\muV)')
fclose('all');
end