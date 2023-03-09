function plotRawdata(chsp,tID,DataType,plotspktime)
%DataType dn mu or amplifier
nChn=128;
FS=30000;
timebefore=0.09;%in seconds
timeafter=0.09;%in seconds
trig=loadTrig;
TrialParams=loadTrialParams;
numelect=find(diff(cell2mat(TrialParams(:,1)))~=0,1,'first');
TrialParamstID = find(cell2mat(TrialParams(1:numelect:end,2)) == tID); %identifies trial row matching trial ID
savespktimes=[];
trigtID = trig(TrialParamstID);
trigtID(trigtID==-500)=[];

filepath = pwd;
[filepathm,name,ext] = fileparts(filepath);
name = name(1:end-14);
if plotspktime==1
    sp=load([name '.sp.mat'],'sp');
    sp=sp.sp{chsp};
    thresh=load([name '.sp.mat'],'thresh');
     thresh=thresh.thresh{chsp};
end
if strcmp(DataType,'amplifier')
    fileID=fopen('amplifier.dat','r');
elseif strcmp(DataType,'dn')
    fileID=fopen('amplifier_dn_sab.dat','r');
elseif strcmp(DataType,'mu')
    fileID=fopen([name '.mu_sab.dat'],'r');
end

try
    ftell(fileID)
catch
    return
end
shortbytes=2;
for indT=1:length(trigtID)
    offset=trigtID(indT)*nChn*shortbytes-timebefore*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
    fseek(fileID,offset,'bof');
    %fseek(fileID2,offset,'bof');
    ftell(fileID)
    if strcmp(DataType,'amplifier') || strcmp(DataType,'dn')
        vblankmu = fread(fileID,[nChn, (timebefore*FS+timeafter*FS)],'int16') .* 0.195;
    elseif strcmp(DataType,'mu')
        vblankmu = fread(fileID,[nChn, (timebefore*FS+timeafter*FS)],'short')./10; %plots one second from trigger and 250ms brefore
    end
    
    
    %
    %if any(vblankmu(chsp,:)<-50)
    %         figure
    %         plot (-0.02*1000:1000/FS:0.02*1000-1000/FS,vblankmu2(chsp,:))
    %         title(['Channel ' num2str(chsp)])
    %         xlabel('Time (ms)')
    %         ylabel('Voltage (uV)')

    figure
    hold on
    plot (-timebefore*1000:1000/FS:timeafter*1000-1000/FS,vblankmu(chsp,:))
    title(['Channel ' num2str(chsp)])
    xlabel('Time (ms)')
    ylabel('Voltage (uV)')
    if plotspktime==1
        yline(thresh);
        offsetspk=trigtID(indT)*1000./FS;
        spktimes=sp(((sp(:,1)>(offsetspk-timebefore*1000))&(sp(:,1)<(offsetspk+timeafter*1000))),:);
        for numspk=1:size(spktimes,1)
        scatter(spktimes(numspk,1)-offsetspk,min(spktimes(numspk,2:end)),'r');
        end
        
        savespktimes=[savespktimes; (spktimes(:,1)-offsetspk)];
    end

    %end
end
fclose(fileID);
end