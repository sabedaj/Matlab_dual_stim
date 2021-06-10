function [IDstruct, baslinespikestruct] = sortTrials_SM(startpointms,mstoanalyse,trig,printspiking,starttrial,trialjump,endtrial,varargin)
%Sort into trial IDs and plot example spikes 
%OUTPUT - the structure containing spiking information for each trial/repeat where
%each cell relates to one trial ID. The array in this cell is in the format
%of (electrodes,trials/repeats) where the total number of spikes per 
%repeat is calculated for the given time

%INPUT - starting point after trigegr to begin counting spikes
%(startpointms) in ms, end point after trigger to stop counting spikes
%(mstoanalyse) in ms, trigger information (trig)

%Elecrode properties
filepath = pwd;
fourShank_cutoff = datetime('03-Aug-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
if (datetime(fileinfo.date) < fourShank_cutoff)
    nChn=32;
    E_Mapnumber=0;
else
    E_Mapnumber=loadMapNum;
    if E_Mapnumber>0
        nChn=64;
    else
        nChn=32;
    end
end
E_MAP = Depth(E_Mapnumber);
FS=30000;

loadThreshold;
TrialParams = loadTrialParams;
loadNREP;%number repeats
spike =0;
[filepathm,name,ext] = fileparts(filepath);
name = name(1:end-14);
load([name '.sp.mat'])
maxtid=max(cell2mat(TrialParams(:,2)));
nospI=[];
baslinespikestruct=[];
IDstruct=[];
Spike_trialstruct=[];

% C5 = intersect(spikedetailstrig(:,1),spikedetailstrig1(:,1)); %middle electrodes
% test3= intersect(C5,C3);

% remove values of common spike times across the array
C1 = intersect(sp{1}(:,1),sp{2}(:,1));%checks later electrodes
C2 = intersect(sp{24}(:,1),sp{6}(:,1)); %early electrodes
C3 = intersect(sp{21}(:,1),sp{5}(:,1)); %middle electrodes
test1= intersect(C1,C2);
test2= intersect(C2,C3);
check=intersect(test1,test2);
for chncount=1:nChn
    [C,r,c]=(intersect(sp{chncount}(:,1),check));
    sp{chncount}(r,:)=[];
end
dispstat('','init');
dispstat(sprintf('Working through trial: '),'keepthis','n');
% Sort into trial IDs and plot spikes for trial 1
for tID=starttrial:trialjump:endtrial
    dispstat(sprintf('%d',tID));
    if isempty(TrialParams(cell2mat(TrialParams(1:end,2)) == tID))
        fprintf('No trial ID: %d\n',tID)
    else
        TrialParamstID = find(cell2mat(TrialParams(1:end,2)) == tID); %identifies trial row matching trial ID
        TrialParamstID(1:2:end)=[];
        TrialParamstID=TrialParamstID./2;
        trigtID = trig(TrialParamstID)./(FS/1000);
        trigtID(trigtID==-500/(FS/1000))=[];
        nTrig = length(trigtID);
        for indT=1:nTrig
            for chsp=1:1:nChn
                v = sp{chsp};
                spikedetailstrig=v(((v(:,1)>(trigtID(indT)+startpointms))&(v(:,1)<(trigtID(indT)+mstoanalyse))),:); %v(((v(:,1)>0)&(v(:,1)<20000)),:);
                %spikedetailstrig=v(((v(:,1)>0)&(v(:,1)<20000)),:);
                timems=mstoanalyse-startpointms;
                avgtimebs=10;
                baslinespiketrig=v(((v(:,1)>(trigtID(indT)-5-avgtimebs*timems))&(v(:,1)<(trigtID(indT)-5))),:); %v(((v(:,1)>0)&(v(:,1)<20000)),:);

                %spikedetailstrig=v(((v(:,1)>(trigtID(indT)+startpointms)/1000)&(v(:,1)<(trigtID(indT)+mstoanalyse)/1000)),:);
                if ~isempty(spikedetailstrig)
                    if any(spikedetailstrig(:,2:end)>(500))%used to check if any are above the specified threshold
                        %Use this to plot the .mu data for any trial
                        fileID=fopen([name '.mu_sab.dat'],'r');
                        shortbytes=2;
                        offset=trig(TrialParamstID(indT))*nChn*shortbytes-0.065*FS*shortbytes*nChn;%offset from beginning of file to trigger-250ms
                        ftell(fileID)
                        fseek(fileID,offset,'bof');
                        ftell(fileID)
                        vblankmu = fread(fileID,[nChn, (0.008*FS+0.065*FS)],'short')./10; %plots one second from trigger and 250ms brefore
                        figure 
                        plot (-0.065*1000:1000/FS:0.008*1000-1000/FS,vblankmu(chsp,:))
                        title(['Channel ' num2str(chsp)])
                        xlabel('Time (ms)')
                        ylabel('Voltage (uV)')
                        fclose(fileID);
                        
                        fileID=fopen('amplifier.dat','r');
                        
                        shortbytes=2;
                        offset=trig(TrialParamstID(indT))*nChn*shortbytes-0.5*FS*shortbytes*nChn;%offset from beginning of file to trigger
                        ftell(fileID)
                        fseek(fileID,offset,'bof');
                        ftell(fileID)
                        vblankmu = fread(fileID,[nChn, (0.5*FS+0.01*FS)],'int16') .* 0.195;
                        figure 
                        plot (-0.5*1000:1000/FS:0.01*1000-1000/FS,vblankmu(chsp,:))
                        title(['Channel ' num2str(chsp)])
                        xlabel('Time (ms)')
                        ylabel('Voltage (uV)')
                        fclose(fileID);
                        
                       
                        fileID=fopen('amplifier_dn_sab.dat','r');
                        shortbytes=2;
                        offset=trig(TrialParamstID(indT))*nChn*shortbytes-0.25*FS*shortbytes*nChn;%offset from beginning of file to trigger
                        ftell(fileID)
                        fseek(fileID,offset,'bof');
                        ftell(fileID)
                        vblankmu = fread(fileID,[nChn, (1*FS+0.25*FS)],'int16') .* 0.195;
                        figure 
                        plot (-0.25*1000:1000/FS:1*1000-1000/FS,vblankmu(chsp,:))
                        title(['Channel ' num2str(chsp)])
                        xlabel('Time (ms)')
                        ylabel('Voltage (uV)')
                        fclose(fileID);
                        
                        fileID=fopen('amplifier.dat','r');
                        shortbytes=2;
                        offset=trig(TrialParamstID(indT))*nChn*shortbytes-0.25*FS*shortbytes*nChn;%offset from beginning of file to trigger
                        ftell(fileID)
                        fseek(fileID,offset,'bof');
                        ftell(fileID)
                        vblankmu = fread(fileID,[nChn, (1*FS+0.25*FS)],'int16') .* 0.195;
                        figure 
                        plot (-0.25*1000:1000/FS:1*1000-1000/FS,vblankmu(chsp,:))
                        title(['Channel ' num2str(chsp)])
                        xlabel('Time (ms)')
                        ylabel('Voltage (uV)')
                        fclose(fileID);
                        
                    end
                    spike=size(spikedetailstrig,1);
                    
                    if printspiking>0
                        if (size(varargin,2)==0) || (cell2mat(varargin(1))==0)
                            if ((tID>15) && (tID<37)) && (((indT>=10) && (indT<=13))) %%||( (tID>15) && (tID<19)) %
                                figure(find(E_MAP==chsp))
                                hold on
                                title(['Channel ' num2str(find(E_MAP==chsp))])
                                ylabel('Spike amplitude (uV)')
                                plot(1*1000/FS:1000/FS:49*1000/FS,spikedetailstrig(:,2:50),'k')
                                xlabel('Time (ms)')
                                
                            end
                        else
                            if ((tID>0) && (tID<=maxtid)) && (indT<=1)% prints all spikes after trigger
                                figure(find(E_MAP==chsp))
                                hold on
                                plot(1*1000/FS:1000/FS:49*1000/FS,spikedetailstrig(:,2:50),'k')
                                title(['Channel ' num2str(find(E_MAP==chsp))])
                                ylabel('Spike amplitude (uV)')
                                xlabel('Time (ms)')
                            end
                        end
                    end
                else
                    spike=0;
                end
                if ~isempty(baslinespiketrig)
                    spikebaseline=size(baslinespiketrig,1)/avgtimebs;
                else
                    spikebaseline=0;
                end
                mchsp=E_MAP(chsp);
                cnum=['Chn_' num2str(mchsp)];
                IDnum=['ID_' num2str(tID)];
                tnum=['Trial_' num2str(indT)];
                if isfield(Spike_trialstruct,cnum)&& isfield(Spike_trialstruct.(cnum),(IDnum)) && isfield(Spike_trialstruct.(cnum).(IDnum),(indT))
                        Spike_trialstruct.(cnum).(IDnum).(tnum)=[Spike_trialstruct.cnum.IDnum.tnum; spikedetailstrig];
                        baslinespike_trialstruct.(cnum).(IDnum).(tnum)=[ baslinespike_trialstruct.cnum.IDnum.tnum;  baslinespiketrig];
                else
                    Spike_trialstruct.(cnum).(IDnum).(tnum)=spikedetailstrig;
                    baslinespike_trialstruct.(cnum).(IDnum).(tnum)=baslinespiketrig;
                end

                nospI(chsp,indT)=spike;
                basespI(chsp,indT)=spikebaseline;
                spike=0;
            end
        end
        IDstruct=StructureIDgeneration(nospI, tID, IDstruct);
        baslinespikestruct=StructureIDgeneration(basespI, tID, baslinespikestruct);
        basespI=[];
        nospI=[];
    end
end

save('Spikes_trialsorted.mat','Spike_trialstruct','baslinespike_trialstruct')
end

