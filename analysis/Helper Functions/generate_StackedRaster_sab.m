function generate_StackedRaster_sab(chn)
dbstop if error
%% Load in data from current directory
loadVarAmp;
plotHist=0;
plot_individualhist=0;
plot_raster=1;
loadNREP;
nT=n_REP_true;
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
d = Depth(E_Mapnumber); chn = d(chn);
sp = loadSpikes(chn);
sp = denoiseSpikes(sp,chn);
TP = loadTrialParams;
trialinfo=TrialInfo_sab;
trialinfo(1,:)=[];

cond= find(diff(cell2mat(trialinfo(:,18))),1,'first')/2; %condition 
totalnumbercond=length(trialinfo(:,2))/(2*cond); %gives the total number of condition changes i.e. 100/0 75/25 but with same electrode are classed as one condition change
loadNORECORDELECT;
loadStimChn;
loadVarAmp;
loadSingleDual;
checkmaxsize=0;
check=0;
for count=1:length(CHN)
    for recordcount=1:length(NORECORDELECT)
        clear relatedtrials
        desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)); %finds trials with desired initial electrode
        %single stim
        if (desiredchanneltrial(end)==length(trialinfo)) %need to remove the last trial in case there is a stim electrode pair resulting in 1 as the last position
            desiredchanneltrial(end)=[];
        end
        desiredchanneltrial_plus=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(CHN(count)+NORECORDELECT(recordcount)+1)),1);% finds triials with matching recording electrode spacing
        if (VarAmp==1)||(singanddual==1)
            modified_desiredchanneltrial_plus= desiredchanneltrial_plus;
            for condition=1:totalnumbercond
                if any(modified_desiredchanneltrial_plus<(cond*condition*2)) && (check>0)
                    check=check+1;
                    relatedtrials(:,check)=modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)));
                elseif any(modified_desiredchanneltrial_plus<(cond*condition*2))
                    relatedtrials=modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)));
                    check=1;
                end
                modified_desiredchanneltrial_plus((modified_desiredchanneltrial_plus<(cond*condition*2)))=[];
                
                if  isempty(modified_desiredchanneltrial_plus)
                    break;
                end
            end
            if exist('relatedtrials', 'var')
                relatedtrials( :, all(~relatedtrials,1) ) = [];
                desiredchannel__singleamp=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                desiredchannel__singleamp(desiredchannel__singleamp<checkmaxsize)=[];
                desiredchanneltrial=find(cell2mat(trialinfo(:,2))==CHN(count)+NORECORDELECT(recordcount)+1); %finds trials with desired initial electrode
                if desiredchanneltrial(end)==length(trialinfo)
                    desiredchanneltrial(end)=[];
                end
                desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                if VarAmp == 1
                    desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<max(relatedtrials,[],'all'));
                end
                if length(desiredchannel__singleampmath)~=size(relatedtrials,2)
                    desiredchannel__singleampmath=desiredchanneltrial((cell2mat(trialinfo(desiredchanneltrial+1,2))==(0)),1);% finds triials with matching recording electrode spacing
                    desiredchannel__singleampmath=desiredchannel__singleampmath(desiredchannel__singleampmath<(cond*2*length(relatedtrials)));
                end
                relatedtrials_all=[desiredchannel__singleampmath';relatedtrials; desiredchannel__singleamp'];
                relatedtrials_all=(relatedtrials_all+1)./2;
            end
        else
            relatedtrials_all=(desiredchanneltrial_plus+1)./2;
        end
    
    if VarAmp~=0
        checkmaxsize=max(relatedtrials,[],'all');
    end
    sum_xhist=0;
    figure; hold on; axM = gca;
    %%%%%%%%%%%%%raster input
    for counter_trials=1:size(relatedtrials_all,2)
        for counter_trials_row=1:size(relatedtrials_all,1)
            ID=relatedtrials_all(counter_trials_row,counter_trials);
            tID = find((cell2mat(TP(:,2)) == ID));
            tID(1:2:end)=[];
            tID=tID./2;
            trig = loadTrig(0);

%             trig(550/2:588/2)=0;
%             trig(7244/2:7246/2)=0;
%             b = [-500];
%             k = 550/2; %row position, can be 0,1,2 or 3 in this case
%             trig = [trig(1,1:k) b trig(1,k+1:end)];
%             k = 7244/2; %row position, can be 0,1,2 or 3 in this case
%             trig = [trig(1,1:k) b trig(1,k+1:end)];
%             if tID(end)>length(trig)
%                 tID(end)=[];
%             end
            
            trig = trig(tID);
            trig(trig==-500)=[];
            theseTrig = trig./30;
            %lfp = loadLFP(chn);
            %% Set up the raster data structure
            BIN = [-200 200];
            SPACING = 200; SMOOTHING = 2; MAX = 400;
            xdata = [];
            ydata = [];
            % Grab Data for the Inset Figure
            FS = 30000; file = pwd;
            mDIR = dir([file '\*mu_sab.dat']);
            mNAME = mDIR.name;
            mFID = fopen([file '\' mNAME],'r');
            X = BIN(1):1/30:BIN(2) - 1/30;
            mu = zeros(nT,(FS/1e3)*diff(BIN));
            for tr = 1:length(trig)
                theseSp = (sp(sp > theseTrig(tr)+BIN(1) & sp < theseTrig(tr)+BIN(2)) - theseTrig(tr));
                for i = 1:length(theseSp)
                    xdata = [xdata, (theseSp(i) + abs(BIN(1)))]; %#ok<*AGROW>
                    ydata = [ydata, tr*(MAX/nT)];
                end
                % Also, extract MUA waveforms
%                 OFFSET = cast(nChn*2*(trig(tr)+(BIN(1)*FS/1e3)),'int64');
%                 fseek(mFID,OFFSET,'bof');
%                 v = fread(mFID,[nChn, (FS/1e3)*diff(BIN)],'short') ./ 10;
%                 mu(tr,:) = v(chn,:);
                % Also, calculate phase for each trial
                %     phase = generatePhaseVector(lfp,[4,5,6,10,20,80]);
                %     for b = 1:6
                %         thisPhase(tr,b) = phase{b}(cast(theseTrig(tr)/30 + offset(b),'int64'));
                %     end
            end
            %% Plot the raster

            yScale = MAX./nT;
            if plot_raster==1
                for i = 1:size(xdata,2)
                    line([xdata(i) xdata(i)],[ydata(i)-yScale/2.3+(counter_trials_row-1)*MAX+(counter_trials-1)*MAX*size(relatedtrials_all,1) ydata(i)+yScale/2.3+(counter_trials_row-1)*MAX+(counter_trials-1)*MAX*size(relatedtrials_all,1)],'Color','k','LineWidth',7);
                end
            end
            line([200 200],[0 MAX*size(relatedtrials_all,1)*size(relatedtrials_all,2)],'Color','b');
            yline(MAX*size(relatedtrials_all,1)*counter_trials,'LineWidth',3);
            yline((counter_trials_row)*MAX+(counter_trials-1)*MAX*size(relatedtrials_all,1),'LineWidth',1);
            % Add an inset to show it's not artefact
%             axI = axes('Position',[.175 .7 .2 .2]); hold on;
%             xlim(axI,[-10 20]); %xlabel('Time (ms)');
%             for tr = 1:nT
%                 plot(axI,X,mu(tr,:)+(tr-1)*SPACING,'Color','k');
%                 %line(axI,X,thresh*ones(1,length(X))+(tr-1)*SPACING,'Color','r');
%             end
%             ylim(axI,[3500 5000]); yticks(''); xticklabels('');
%             line(axI,[0 0],[0 nT*SPACING],'Color','b');
            % Add the convolved spikerate
            Z = hist(xdata,0:400); %#ok<HIST>
            maxval=max(Z(150:250));
            if plotHist==1
                if plot_individualhist==1
                    window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
                    rate = (1000/nT)*conv(Z,window);
                    rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
                    %rate=Z;
                    plot(axM,rate+MAX*(counter_trials_row-1)+(counter_trials-1)*MAX*size(relatedtrials_all,1),'b','LineWidth',2);
                else
                    sum_xhist=Z+sum_xhist;
                end
            end
            txt = strcat('Max: ', string(maxval), ' Sp/s');
            text(161,MAX*(counter_trials_row-1)+size(relatedtrials_all,1)*MAX*(counter_trials-1)+150,txt,'FontSize',10,'Color','k')
            
        end

            if cell2mat(trialinfo(ID*2,18))==-1
                txt = 'No Stim';
            else
                txt = strcat(string(trialinfo(ID*2,18)), 'uA');
            end
            text(151,MAX*size(relatedtrials_all,1)*counter_trials-150,txt,'FontSize',14)
        
        if (plot_individualhist==0) && (plotHist==1)
            window = normpdf(-3*SMOOTHING:3*SMOOTHING,0,SMOOTHING);
            rate = (1000/nT)*conv(sum_xhist,window);
            rate = rate(3*SMOOTHING+1:end-3*SMOOTHING);
            plot(axM,rate+MAX*size(relatedtrials_all,1)*(counter_trials-1),'b','LineWidth',2);
            sum_xhist=0;
        end

            %text(151,MAX*size(relatedtrials_all,1)*counter_trials-150,txt,'FontSize',14)
    end
    ylim(axM,[0 MAX*size(relatedtrials_all,1)*size(relatedtrials_all,2)]); xlim(axM,[150 300]);
    % Tidy Up
    
    %ylabel(axM,'Firing Rate (Sp/s)');
    xticks(axM,150:50:300);
    xticklabels(axM,-50:50:100);
    %beautifyPlot(30,axI);
    beautifyPlot(13,axM);
    set(gca,'yticklabel',{[]})
    title(['Electrodes ', (num2str(CHN(count))),' and ', (num2str(CHN(count)+NORECORDELECT(recordcount)+1))])
    if VarAmp==1
        title(['Channel ', num2str(find(d==chn)), ' w/ StimChn ', num2str(cell2mat(trialinfo(ID*2-3,2))), ' ', num2str(cell2mat(trialinfo(ID*2-2,2)))])
    elseif singanddual==1
        title(['Channel ', num2str(find(d==chn)), ' w/ StimChn ', num2str(cell2mat(relatedtrials_all(2,1)*2-1,2)), ' Dualstim ', num2str(cell2mat(trialinfo(relatedtrials_all(2,1)*2,2)))])
    else
        title(['Channel ', num2str(find(d==chn)), ' w/ StimChn ', num2str(cell2mat(trialinfo(ID*2,2))), ' ', num2str(cell2mat(trialinfo(ID*2+1,2)))])
    end
    xlabel(axM,'Time (ms)');
    end
end





end