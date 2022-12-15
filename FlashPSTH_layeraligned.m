%%Layeraligned psth for each animal

%setup data structure
positionchn=fliplr(string(0:50:1600));
for dist_chn=1:33
Flashresp_depthsep.(strcat('D', positionchn(dist_chn)))=[];
end

fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir('info.rhs');
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=1;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
BIN = [-200 200]; MAX = 400;
FS = 30000; 
d = Depth(E_Mapnumber); 

%% Run for each respective animal
trig=loadTrig(0);%cleanTrig_flash; if no trig file saved
trig(trig==-500)=[];
%trig=trig(700:900); % to limit the number of trigs
theseTrig = trig./30;
nT=length(trig);
loadNREP;

%hist data
for chn=1:nChn
    xdata = [];
    ydata = [];
    chn1 = d(chn);
    sp = loadSpikes(chn1);
    for tr = 1:nT
        theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
        for i = 1:length(theseSp)
            xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
            ydata = [ydata, tr*(MAX/nT)];
        end
    end

    Z = hist(xdata,0:400); %#ok<HIST>
    if chn<17
        shankplot=1;
        chn_shank=chn;
    elseif chn<33
        shankplot=2;
        chn_shank=chn-16;
    elseif chn<49
        shankplot=3;
        chn_shank=chn-32;
    else
        shankplot=4;
        chn_shank=chn-48;
    end
    shanklayclass=ElectLayerClass(1+((shankplot-1)*16):(shankplot*16));
    WMborder=find(shanklayclass==4,1,'last');
    if isempty(WMborder)
        GIborder=find(shanklayclass==3,1,'last');
        if isempty(GIborder)
            dist=chn_shank+(24-SGborder);%SG electrode at e24
            SGborder=find(shanklayclass==2,1,'last');
            Flashresp_depthsep.(strcat('D', positionchn(dist_chn)))=[Flashresp_depthsep.(strcat('D', positionchn(dist_chn))); Z];%pos 24 is gran pos 25 is supra
        else
            dist_chn=chn_shank+(21-GIborder);%GI electrode at e21
            Flashresp_depthsep.(strcat('D', positionchn(dist_chn)))=[Flashresp_depthsep.(strcat('D', positionchn(dist_chn))); Z];%pos 21 is infra, pos 22 is gran
        end
    else
        dist_chn=chn_shank+(7-WMborder);%WM electrode at e7
        Flashresp_depthsep.(strcat('D', positionchn(dist_chn)))=[Flashresp_depthsep.(strcat('D', positionchn(dist_chn))); Z];%pos 7 is WM, pos 8 is infra
    end
end
beep