function DepthChangeingSpiking_SM(avgnospT, starttrial, trialjump, varargin)

%   Plots the response curve according to trial change over time. 
%   Analyses spiking according to trial ID and Time. Creates a 3
%   dimensional plot with three axes. Intpus:
% startpointseconds - How long after the trigger do you want skip spike analysis(ms)? 
% secondstoanalyse - How long after the trigger do you want to analyse spikes for(ms)? 
%TimeStep -Timestep of spike analysis between TimeBegin and TimeEnd(ms)? 
%Eletrode of interest number
%Start trial, jump trial(usually = cond)
%End trial
%%
trig = loadTrig(0);
theseTrig = trig;
trialinfo=TrialInfo_sab;
trialinfo(1,:)=[];
cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition
TrialParams = loadTrialParams;
maxtid=max(cell2mat(TrialParams(:,2)));
loadStimChn;
filepath = pwd;
fourShank_cutoff = datetime('04-Aug-2020 00:00:00');
fileinfo = dir([filepath,'\info.rhs']);
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
if nargin ==4
    endtrial=cell2mat(varargin(1));
else
    endtrial=maxtid;
end
    x=zeros(nChn,floor((endtrial-starttrial)/trialjump));
    y=zeros(nChn,floor((endtrial-starttrial)/trialjump));
    z=zeros(nChn,floor((endtrial-starttrial)/trialjump));
    for trial=1:ceil((endtrial-starttrial)/trialjump)
        z(:,trial)=avgnospT(:,starttrial+(trial-1)*trialjump);
        if cond~=starttrial
            y(:,trial)=cell2mat(trialinfo(trialjump*trial*2+floor(starttrial/cond)*cond,18));
        else
            y(:,trial)=cell2mat(trialinfo(trialjump*trial*2,18));
        end
        if y(:,trial)==-1
            y(:,trial)=0;
        end
        x(:,trial)=1:nChn;
    end
  
%     figure
%     hold on
%     surf(x,y,z)
%     xlabel(x_namelabel)
%     ylabel('Average spiking')
%     zlabel('Time (ms)')
%     title(title_namelabel)
%     xticks(1:size(xticks_names(1,:),2))
%     xticklabels(xticks_names(1,:))
%     xlim([1 size(xticks_names(1,:),2)])
%     
    figure
    hold on
    surf(y,x,z)
    xlabel('Amplitude (uA)')
    zlabel('Average spiking')
    ylabel('Channel number')
    if (cond~=starttrial) && (starttrial~=1)
        title_namelabel=['Channel changes in spiking. Stimchn: ' num2str(stimChn(1)) ' ' num2str(stimChn(2)) ' @ ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))*100/(cell2mat(trialinfo((starttrial+trialjump)*2-1,18))+cell2mat(trialinfo((starttrial+trialjump)*2,18)))) '/' num2str(cell2mat(trialinfo((starttrial+trialjump)*2,18))*100/((cell2mat(trialinfo((starttrial+trialjump)*2-1,18)))+cell2mat(trialinfo((starttrial+trialjump)*2,18))))];
    else
        title_namelabel=['Channel changes in spiking. Stimchn: ' num2str(cell2mat(trialinfo((starttrial+trialjump)*2-1,2)))];
    end
    title(title_namelabel)
    colorbar
    ylim([1 nChn])
    caxis([0, max(z,[],'all')]);
    
%     for loopnumcount=1:maxtid/cond
%         hold on
%         xline((cond*loopnumcount)+1);
%     end
 
  
    
    %%
end