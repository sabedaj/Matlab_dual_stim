function plotRawdata(chsp,tID)

nChn=128;
FS=30000;
trig=loadTrig;
TrialParams=loadTrialParams;
TrialParamstID = find(cell2mat(TrialParams(1:end,2)) == tID); %identifies trial row matching trial ID
trigtID = trig(TrialParamstID)./(FS/1000);
trigtID(trigtID==-500/(FS/1000))=[];

filepath = pwd;
[filepathm,name,ext] = fileparts(filepath);
name = name(1:end-14);

fileID=fopen([name '.mu_sab.dat'],'r');
fileID2=fopen(['amplifier_dn_sab.dat'],'r');
try
    ftell(fileID)
catch
    return
end
shortbytes=2;
for indT=1:40
    offset=trig(TrialParamstID(indT))*nChn*shortbytes-0.020*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
    fseek(fileID,offset,'bof');
     fseek(fileID2,offset,'bof');
    ftell(fileID)
    vblankmu = fread(fileID,[nChn, (0.020*FS+0.02*FS)],'short')./10; %plots one second from trigger and 250ms brefore
    vblankmu2 = fread(fileID2,[nChn, (0.020*FS+0.02*FS)],'int16') .* 0.195;
    %if any(vblankmu(chsp,:)<-50)
        figure
        plot (-0.02*1000:1000/FS:0.02*1000-1000/FS,vblankmu2(chsp,:))
        title(['Channel ' num2str(chsp)])
        xlabel('Time (ms)')
        ylabel('Voltage (uV)')
        
            figure
        plot (-0.02*1000:1000/FS:0.02*1000-1000/FS,vblankmu(chsp,:))
        title(['Channel ' num2str(chsp)])
        xlabel('Time (ms)')
        ylabel('Voltage (uV)')
    %end
end
fclose(fileID);
end