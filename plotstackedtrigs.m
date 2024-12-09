function plotstackedtrigs(chsp, tID, trigs, DataType)
closenessvar=120;
nChn=128;
FS=30000;
timebefore=0.09;%in seconds 0.09
timeafter=0.09;%in seconds0.09
trig=loadTrig;
TrialParams=loadTrialParams;
numelect=find(diff(cell2mat(TrialParams(:,1)))~=0,1,'first');
TrialParamstID = find(cell2mat(TrialParams(1:numelect:end,2)) == tID); %identifies trial row matching trial ID
savespktimes=[];
trigtID = trig(TrialParamstID);
trigtID(trigtID==-500)=[];
trigtID=trigtID(trigs);
filepath = pwd;
[filepathm,name,ext] = fileparts(filepath);
name = name(1:end-14);

if strcmp(DataType,'amplifier')
    fileID=fopen('amplifier.dat','r');
elseif strcmp(DataType,'dn')
    fileID=fopen('amplifier_dn_sab.dat','r');
elseif strcmp(DataType,'mu')
    fileID=fopen([name '.mu_sab.dat'],'r');
elseif strcmp(DataType,'DT')
    fileID=fopen([name '_DT.mu.dat'],'r');
end
try
    ftell(fileID)
catch
    return
end
figure
    hold on
shortbytes=2;
loop=0;
for indT=1:length(trigs)
    
    offset=trigtID(indT)*nChn*shortbytes-timebefore*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
    fseek(fileID,offset,'bof');
    %fseek(fileID2,offset,'bof');
    ftell(fileID)
    if strcmp(DataType,'amplifier') || strcmp(DataType,'dn')
        vblankmu = fread(fileID,[nChn, (timebefore*FS+timeafter*FS)],'int16') .* 0.195;
    else 
        vblankmu = fread(fileID,[nChn, (timebefore*FS+timeafter*FS)],'short')./10; %plots one second from trigger and 250ms brefore
    end

    plot(-timebefore*1000:1000/FS:timeafter*1000-1000/FS,vblankmu(chsp,:)+closenessvar*loop,'k')
    loop=loop+1;
end
xlabel('Time (ms)')
plot([60 60],[0 500],'b','LineWidth',2)
text(65,250,'500\muV','Color','b')
ylim([-100 4800])
xlim([-90 90])
set(gca, 'YTick', [])
set(gca,'TickDir','out');
xline(0,'r');

end