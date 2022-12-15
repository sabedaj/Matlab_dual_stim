%% save array layer classification
ratN='Rat_006';
superL=[];% electrodes in superficial layers
inputL=[15:16 15+32:32+16 15+48:16+48 15+16:16+16];%electrodes in input layers
deepL=[1:14 (1+32):(14+32) (1+48):(14+48) (1+16):(14+16)];% electrodes in deep layers
WM=[];

ElectLayerClass=zeros(64,1);
if any(superL~=0)
    ElectLayerClass(superL)=1;
end
if any(inputL~=0)
    ElectLayerClass(inputL)=2;
end
if any(deepL~=0)
    ElectLayerClass(deepL)=3;
end
if any(WM~=0)
    ElectLayerClass(WM)=4;
end
layers_depthdefinition=[368 522 1240];
Depthestimate.s1=zeros(16,1);
Depthestimate.s2=zeros(16,1);
Depthestimate.s3=zeros(16,1);
Depthestimate.s4=zeros(16,1);
shankorder=[1 4 2 3];
for shankiterate=0:3
    checkshank=['s' num2str(shankorder(shankiterate+1))];
    if any(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==4)
        firste1=find(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==3,1,'first');
        Depthestimate.(checkshank)(1:firste1-1)=layers_depthdefinition(end)+50*(firste1-1):-50:layers_depthdefinition(end)+50;
        Depthestimate.(checkshank)(firste1:end)=layers_depthdefinition(end):-50:layers_depthdefinition(end)-50*(16-firste1);
    else
        firste1=find(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==2,1,'first');
        Depthestimate.(checkshank)(1:firste1-1)=layers_depthdefinition(2)+50*(firste1-1):-50:layers_depthdefinition(2)+50;
        Depthestimate.(checkshank)(firste1:end)=layers_depthdefinition(2):-50:layers_depthdefinition(2)-50*(16-firste1);
    end
end
% firste2=1077;
% firste3=1024;
% firste4=1039;
% Depthestimate.s1=firste1:-50:firste1-50*16;
% Depthestimate.s2=firste2:-50:firste2-50*16;
% Depthestimate.s3=firste3:-50:firste3-50*16;
% Depthestimate.s4=firste4:-50:firste4-50*16;
% layers_depthdefinition=[382 555 1240];
% layers_depthdefinition_start=[0 190+353+149+165 525+700+190+353+149+165];

save(['E:\DATA\' ratN filesep 'ElectLayerClass.mat'],'ElectLayerClass')

%%
%r23
% superL=[15:16 (15+32):(16+32) (15+48):(16+48) (16+16)];% electrodes in superficial layers
% inputL=[11:14 (11+32):(14+32) (12+48):(14+48) (13+16):(15+16)];%electrodes in input layers
% deepL=[1:10 (1+32):(10+32) (1+48):(11+48) (1+16):(12+16)];% electrodes in deep layers
% WM=[];

%r22
%%LFP
% superL=[15:16 (15+32):(16+32) (16+48) (16+16)];% electrodes in superficial layers
% inputL=[11:14 (11+32):(14+32) (13+48):(15+48) (13+16):(15+16)];%electrodes in input layers
% deepL=[1:10 (1+32):(10+32) (1+48):(12+48) (1+16):(12+16)];% electrodes in deep layers
% WM=[];
%%hist only
% superL=[12:16 (12+32):(16+32) (12+48):(16+48) (13+16):(16+16)];% electrodes in superficial layers
% inputL=[8:11 (9+32):(11+32) (9+48):(11+48) (10+16):(12+16)];%electrodes in input layers
% deepL=[1:7 (1+32):(8+32) (1+48):(8+48) (1+16):(9+16)];% electrodes in deep layers
% WM=[];

%R21
% superL=[];% electrodes in superficial layers
% inputL=[];%electrodes in input layers
% deepL=[2:16 (2+32):(16+32) (3+48):(16+48) (3+16):(16+16)];% electrodes in deep layers
% WM=[1 (1+32) (1+48):(2+48) (1+16):(2+16)];

%R13
% superL=[];% electrodes in superficial layers
% inputL=[];%electrodes in input layers
% deepL=[8:16 (8+32):(16+32) (8+48):(16+48) (8+16):(16+16)];% electrodes in deep layers
% WM=[1:7 (1+32):(7+32) (1+48):(7+48) (1+16):(7+16)];

%R12
%LFP



%R9
% superL=[];% electrodes in superficial layers
% inputL=[];%electrodes in input layers
% deepL=[6:16 (6+32):(16+32) (6+48):(16+48) (6+16):(16+16)];% electrodes in deep layers
% WM=[1:5 (1+32):(5+32) (1+48):(5+48) (1+16):(5+16)];

%R8
% superL=[];% electrodes in superficial layers
% inputL=[16 32 48 64];%electrodes in input layers
% deepL=[2:15 18:31 34:47 50:63];% electrodes in deep layers
% WM=[1 17 33 49];

%R6
% superL=[];% electrodes in superficial layers
% inputL=[15:16 15+32:32+16 15+48:16+48 15+16:16+16];%electrodes in input layers
% deepL=[1:14 (1+32):(14+32) (1+48):(14+48) (1+16):(14+16)];% electrodes in deep layers
% WM=[];
%% reading connectome data - https://bbp.epfl.ch/nmc-portal/downloads.html


fileName = 'C:\Users\smei0006\Documents\Experimental_Design\pathways_anatomy_factsheets_simplified.json'; % filename in JSON extension synapse count
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data_L = jsondecode(str); % Using the jsondecode function to parse JSON from string
fn_L = fieldnames(data_L);
from_L_str=[{'1'} {'2'} {'4'} {'5'} {'6'}];
to_L_str=[{'1'} {'2'} {'4'} {'5'} {'6'}];

fileName = 'C:\Users\smei0006\Documents\Experimental_Design\pathways_physiology_factsheets_simplified.json'; % filename in JSON extension excite/inhib
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data_anatomy_L = jsondecode(str); % Using the jsondecode function to parse JSON from string
names_anatomy_L = fieldnames(data_anatomy_L);
firing_rates_EI=[1.09 6.00];%https://www.sciencedirect.com/science/article/pii/S0092867415011915#sec4
%%
it=0;
Numconnections=[];
L_combinations=[];
totalconnect_L=[];
 Numconnections_3groups=[];
totalconnect=0;
for from_L=1:5 %creats structures for storing connection data
    if from_L<4
        from_L_group='s';
    elseif from_L==4
        from_L_group='g';
    elseif from_L>4
        from_L_group='i';
    end
    for to_L=1:5
        %split into original layers
        it=it+1;
        L_combinations{it}=['L' strcat(from_L_str{from_L}, to_L_str{to_L})];
        Numconnections.(L_combinations{it})=0;
        
        %split into 3 groups
        if to_L<4
            to_L_group='s';
        elseif to_L==4
            to_L_group='g';
        elseif to_L>4
            to_L_group='i';
        end
        L_combinationsgroup{it}=strcat(from_L_group, to_L_group);

        Numconnections_3groups.(L_combinationsgroup{it})=0;
        count_inhib_excite.(L_combinationsgroup{it})=[0 0];
    end
    totalconnect_L.(L_combinations{it}(1:2))=0;
end
totalconnect_L.all=0;

 Numconnections_3groups_array=zeros(3,3);
 inhib_excite_arrayFR=zeros(3,3);
 inhib_excite_arrayCOUNT=zeros(3,3);

for it_fn=1:length(fn_L) %calculates number of connections while looping through possible neuron combinations
    underscore_L = strfind(fn_L{it_fn},'___');
    Lconnect=(fn_L{it_fn}([2 underscore_L+4]));
    L_check=['L' Lconnect];
    
    Numconnections.(L_check)=Numconnections.(L_check)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
    
    totalconnect_L.(L_check(1:2))=totalconnect_L.(L_check(1:2))+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
    totalconnect_L.all=totalconnect_L.all+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
    
    %3 group split
    if str2double(fn_L{it_fn}(underscore_L+4))<4
        to_L_group='s';
        to_L_groupnum=1;
    elseif str2double(fn_L{it_fn}(underscore_L+4))==4
        to_L_group='g';
        to_L_groupnum=2;
    elseif str2double(fn_L{it_fn}(underscore_L+4))>4
        to_L_group='i';
        to_L_groupnum=3;
    end
    if str2double(fn_L{it_fn}(2))<4
        from_L_group='s';
        from_L_groupnum=1;
    elseif str2double(fn_L{it_fn}(2))==4
        from_L_group='g';
        from_L_groupnum=2;
    elseif str2double(fn_L{it_fn}(2))>4
        from_L_group='i';
        from_L_groupnum=3;
    end
    
    synapse_type_connection=data_anatomy_L.(names_anatomy_L{it_fn}).synapse_type;

    inhib_1=strfind(synapse_type_connection,'Inhibitory')+1;
    if isempty(inhib_1)
        inhib_1=1;
    end
     L_3layercheck=strcat(from_L_group, to_L_group);
     Numconnections_3groups.(L_3layercheck)=Numconnections_3groups.(L_3layercheck)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection)-data_L.(fn_L{it_fn}).number_of_convergent_neuron_mean+data_L.(fn_L{it_fn}).number_of_divergent_neuron_mean;
     Numconnections_3groups_array(to_L_groupnum,from_L_groupnum)=Numconnections_3groups_array(to_L_groupnum,from_L_groupnum)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection)-data_L.(fn_L{it_fn}).number_of_convergent_neuron_mean+data_L.(fn_L{it_fn}).number_of_divergent_neuron_mean;
     count_inhib_excite.(L_3layercheck)(inhib_1)=count_inhib_excite.(L_3layercheck)(inhib_1)+1;
     inhib_excite_arrayFR(to_L_groupnum,from_L_groupnum)=inhib_excite_arrayFR(to_L_groupnum,from_L_groupnum)+firing_rates_EI(inhib_1);
     inhib_excite_arrayCOUNT(to_L_groupnum,from_L_groupnum)=inhib_excite_arrayCOUNT(to_L_groupnum,from_L_groupnum)+1;
end
layers_depthdefinition=[353+149+165 190+353+149+165 525+700+190+353+149+165];%thickness of microcircuit
layers_thickness=[368 154 718];
percent_connection=layers_thickness./(layers_depthdefinition-[0 667 857]);%microcircuit is thicker than vis cortex
Numconnections_3groups_array=Numconnections_3groups_array.*percent_connection;
%% plot layer connections
cmp=[0.6350, 0.0780, 0.1840; 0.8500, 0.3250, 0.0980;0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.4940, 0.1840, 0.5560]; 
figure(1)
hold on
for from_L=1:5 %calculates number of connections while looping through possible neuron combinations
    for to_L=1:5
        
        L_check=['L' from_L_str{from_L} to_L_str{to_L}];
        
        L_w=plot([1 5], [str2double(from_L_str{from_L}) str2double(to_L_str{to_L})]);
        
        L_w.LineWidth=Numconnections.(L_check).*100./totalconnect_L.all;
        L_w.Color=cmp(str2double(from_L_str{from_L}),:);
    end
end
ylim([0 7])


%% plot 3 layer connections
cmp=[0.9290, 0.6940, 0.1250; 0.4660, 0.6740, 0.1880;0.3010, 0.7450, 0.9330;0.4940, 0.1840, 0.5560]; 
figure(2)
hold on
for from_L=1:5 %calculates number of connections while looping through possible neuron combinations
    if from_L<4
        from_L_group='s';
        from_3layer=1;
    elseif from_L==4
        from_L_group='g';
        from_3layer=2;
    elseif from_L>4
        from_L_group='i';
        from_3layer=3;
    end
    for to_L=1:5
        %split into 3 groups
        if to_L<4
            to_L_group='s';
            to_3layer=1;
        elseif to_L==4
            to_L_group='g';
            to_3layer=2;
        elseif to_L>4
            to_L_group='i';
            to_3layer=3;
        end
        L_check=strcat(from_L_group, to_L_group);
        
        L_w=plot([1 5], [from_3layer to_3layer]);
        
        L_w.LineWidth=Numconnections_3groups.(L_check).*100./totalconnect_L.all;
        L_w.Color=cmp(from_3layer,:);
    end
end
ylim([0 4])

%% finding if superficial or input or deep



%% Initial loop through data to organise it for plotting
%for half amplitude, change poolnormalised activity flag at beginning of
%function. To change to baseline activity, alter Csplit line 216 to basespk_all
%instead of normspk_all. Single animal, just change ratN to the animal - R16
%for 5 elect and R19 for 7,9 elect.
currentdistributions=flipud([1 0; 0.75 0.25; 0.5 0.5; 0.25 0.75; 0 1]);
marmoflag=0;
startpointseconds=2;
secondstoanalyse=8;
numdistsigelect_array=NaN(21,1700);
distsigelect_array=NaN(21,1700);
spkrate_dist=NaN(21,1700);
counterstimchn=0;
counternumberofpairs=0;
for sepdist=-6:-2:-10 %%% Laminar
    AMP=[0 1 2 3 4 5 6 8 10 15];
    sepcheck=['sep' num2str((sepdist+1)*-1)];
    rollingsum.(sepcheck)=nan(16,1000);
    
    altogethersumforavg.(sepcheck)=nan(1,1000);
    for trial=1:5
        check1=['T' num2str(trial)];
        for current = 1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            Csplit.(currcheck).(check1)=[];
            basesplit.(currcheck).(check1)=[];
            Csplit_simulation.(currcheck).(check1)=[];
            CsplitLayers_weighted.(currcheck).(check1)=[];
            CsplitLayers.ss.(currcheck).(check1)=[];
            CsplitLayers.ii.(currcheck).(check1)=[];
            CsplitLayers.gg.(currcheck).(check1)=[];
            CsplitLayers.si.(currcheck).(check1)=[];
            CsplitLayers.sg.(currcheck).(check1)=[];
            CsplitLayers.gi.(currcheck).(check1)=[];
            CsplitLayers.WM.(currcheck).(check1)=[];
            CsplitLayers_model.ss.(currcheck).(check1)=[];
            CsplitLayers_model.ii.(currcheck).(check1)=[];
            CsplitLayers_model.gg.(currcheck).(check1)=[];
            CsplitLayers_model.si.(currcheck).(check1)=[];
            CsplitLayers_model.sg.(currcheck).(check1)=[];
            CsplitLayers_model.gi.(currcheck).(check1)=[];
            CsplitLayers_model.WM.(currcheck).(check1)=[];
            sepdistsig_dual.(sepcheck).(currcheck)=nan(32,1000);
            sepdistsig.(sepcheck).(currcheck)=nan(32,1000);
            CsplitCentroid.(currcheck)=[];
            for shanksep=0:3
                bincheck=['D' num2str(shanksep)];
                bin.(currcheck).(check1).(bincheck)=[];
            end
            for binx=0:5
                binnum=['B' num2str(binx)];
                Csplitbinlayers.(currcheck).(check1).(binnum)=[];
            end
        end
    end
    
    countershank=0;
    stimshank=[];
    AR=[1 2 3 4];%array shanks
    binchecksave='blank';
    
    for ratN=14:20%[6 9 13 21:23]%[6:13 21:23]%14:20 %excluding R7 and R10 and R11
        if ratN<10 && marmoflag==0
            Ratnum=['Rat_00' num2str(ratN)];
        elseif ratN>=10 && marmoflag==0
            Ratnum=['Rat_0' num2str(ratN)];
        else%%marmo
            Ratnum=['Marmo_00' num2str(ratN)];
        end
        AMP=[0 1 2 3 4 5 6 8 10 15];
        
        cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
        D_data=dir;
        if ~any(strcmp({D_data.name}, 'ElectLayerClass.mat'))
            continue
        end
        load('ElectLayerClass.mat','ElectLayerClass')
        layers_depthdefinition=[368 522 1240];
        Depthestimate.s1=zeros(16,1);
        Depthestimate.s2=zeros(16,1);
        Depthestimate.s3=zeros(16,1);
        Depthestimate.s4=zeros(16,1);
        shankorder=[1 4 2 3];
        for shankiterate=0:3
            checkshank=['s' num2str(shankorder(shankiterate+1))];
            if any(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==4)
                firste1=find(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==3,1,'first');
                Depthestimate.(checkshank)(1:firste1-1)=layers_depthdefinition(end)+50*(firste1-1):-50:layers_depthdefinition(end)+50;
                Depthestimate.(checkshank)(firste1:end)=layers_depthdefinition(end):-50:layers_depthdefinition(end)-50*(16-firste1);
            else
                firste1=find(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==2,1,'first');
                Depthestimate.(checkshank)(1:firste1-1)=layers_depthdefinition(2)+50*(firste1-1):-50:layers_depthdefinition(2)+50;
                Depthestimate.(checkshank)(firste1:end)=layers_depthdefinition(2):-50:layers_depthdefinition(2)-50*(16-firste1);
            end
        end
        stimchn_r.(Ratnum)=[];
        savestimshanks=zeros(length(D_data),1);
        for k = 3:length(D_data) % avoid using the first ones
            currD = D_data(k).name; % Get the current subdirectory name
            try
                cd([D_data(k).folder filesep currD])
                loadStimChn;
                if (stimChn(1)-stimChn(2))~=sepdist
                    continue
                end
                load('Averagetrialresponse.mat','avgnospT')
            catch
                stop=0;
                continue
            end
            
            [normspk_all,centroidpershank, basespk_all, respspk_all]=PoolNormalisedActivity(avgnospT,startpointseconds, secondstoanalyse);
            if sum(normspk_all.S1_T1,'all')==0 && sum(centroidpershank,'all')==140
                continue
            end
            stimChn_orig=stimChn;
            Stimchnpos_layers=ElectLayerClass(stimChn);



            if all(Stimchnpos_layers==1)
                layers_stim='ss';
            elseif all(Stimchnpos_layers==2)
                layers_stim='gg';
            elseif all(Stimchnpos_layers==3)
                layers_stim='ii';
            elseif any(Stimchnpos_layers==1) && any(Stimchnpos_layers==2)
                layers_stim='sg';
            elseif any(Stimchnpos_layers==2) && any(Stimchnpos_layers==3)
                layers_stim='gi';
            elseif any(Stimchnpos_layers==1) && any(Stimchnpos_layers==3)
                layers_stim='si';
            elseif any(Stimchnpos_layers==4)
                layers_stim='WM';
                warning('Stim electrodes are in white matter!!!!!')
            else
                error('Layers not saved properly')
            end
            
            
            if stimChn(1)<17 %determines the shank with the stimulating electrodes
                shank=1;
            elseif stimChn(1)<33 && stimChn(1)>16
                shank=4;
                stimChn=stimChn-16;
            elseif stimChn(1)<49 && stimChn(1)>32
                shank=2;
                stimChn=stimChn-32;
            else
                shank=3;
                stimChn=stimChn-48;
            end
            L_GI_boundary=intersect(find(diff(ElectLayerClass)==-1),find(ElectLayerClass==3));
            bin_size=abs(sepdist/2);
            if isempty(L_GI_boundary)
                L_IWM_boundary=intersect(find(diff(ElectLayerClass)==-1),find(ElectLayerClass==4));
                stimchn_bins=ceil((stimChn-L_IWM_boundary(1))/bin_size)-ceil(((layers_depthdefinition(3)-layers_depthdefinition(2))/50)/bin_size);%granular layers are 1, first infra are 0
            else
                stimchn_bins=ceil((stimChn-L_GI_boundary(1))/bin_size);%granular layers are 1, first infra are 0
            end
            stimchn_bins=stimchn_bins+5;%Assuming there will not be more than -5 bins, granular layer is 6 and first infra are 5

            
            AMP=loadAMP;
            AMP(AMP==-1)=0;
            counterstimchn=counterstimchn+1;
            clear baslinespikestruct IDstruct
            load('IDstruct.mat')
            Ampcheck=loadAMP;
            trialinfo=loadTrialInfo(0);
            trialinfo(:,3)=[];
            trialinfo=cell2mat(trialinfo);
            for shanknum=1:4
                if shanknum==2 %used to index data as it is stored shank 1,4,2,3
                    shank_orig=3;
                elseif shanknum==3
                    shank_orig=4;
                elseif shanknum==4
                    shank_orig=2;
                else
                    shank_orig=1;
                end
%                 if shank==shanknum
%                     continue
%                 end
                
                countershank=countershank+1;
                trialcheck=0;
                
                for trial=1:5
                    check=['S' num2str(shanknum) '_T' num2str(trial)];
                    % Identifying the number of significantly responding electrodes based on shank, distance from stimulating electrode, dual and single electrode conditions
                    [~,uniquestimchn]=intersect((stimChn.*shank),stimchn_r.(Ratnum));%multiply by shanknum to ensure unique single stimulation channel
                    if (trial==1 && (isempty(uniquestimchn) || all(uniquestimchn~=1)))  || (trial==5 && (isempty(uniquestimchn)|| all(uniquestimchn~=2))) %single electrode trials only
                        if (trial==5)
                            chnstim=2;
                        else
                            chnstim=1;
                        end
                        if any(Ampcheck==6) && exist('baslinespikestruct')==1
                            trialsAmpmatch=find((trialinfo(:,2)==stimChn_orig(chnstim)) & (trialinfo(:,18)==6));
                            [trialsAmpmatchindex]=find(trialinfo(trialsAmpmatch+1,2)==0);
                            trialmatch=trialinfo(trialsAmpmatch(trialsAmpmatchindex),1);
                            ID=trialmatch;
                            checkID=['T' num2str(ID)];
                            E_MAP = Depth(1);
                            IDstructID=IDstruct.(checkID)(E_MAP,:);
                            baslinespikestructID=baslinespikestruct.(checkID)(E_MAP,:); %comparing baseline to stimualtion trials to determine significance
                            for numelect=1:size(basespk_all.(check),1)
                                [psigchan,hsigchan,~] = signrank(baslinespikestructID(numelect*shank_orig,:), IDstructID(numelect*shank_orig,:),'alpha',0.05,'tail','left');
                                rollingsum.(sepcheck)(numelect,(counterstimchn*2)-(chnstim-1))=nansum([rollingsum.(sepcheck)(numelect,(counterstimchn*2)-(chnstim-1)) hsigchan]); % all significant channels
                                distancerecordstim=round(sqrt(((stimChn(chnstim)-numelect)*50)^2+((shank-shanknum)*200)^2)/50); %distance of the electrode from the single stimulating electrode sorted into 50um bins
                                numdistsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1))=nansum([numdistsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1)) 1]);
                                distsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1))=nansum([distsigelect_array(distancerecordstim+1,(counterstimchn*2)-(chnstim-1)) hsigchan]);
                                if hsigchan==1
                                    spkrate_dist(distancerecordstim+1,(counterstimchn*2)-(chnstim-1))=nansum([spkrate_dist(distancerecordstim+1,(counterstimchn*2)-(chnstim-1)) normspk_all.(check)(numelect,Ampcheck==6)]);
                                end
                            end
                            altogethersumforavg.(sepcheck)((counterstimchn*2)-(chnstim-1))=nansum([altogethersumforavg.(sepcheck)((counterstimchn*2)-(chnstim-1)) 1]);
                        elseif any(Ampcheck==6)
                            stop=0;
                        end
                        if exist('baslinespikestruct')==1
                            for current = 1:length(AMP)
                                if AMP(current)==0 || ~any(Ampcheck==AMP(current))
                                    continue
                                end
                                currcheck=['C' num2str(AMP(current))];
                                trialsAmpmatch=find((trialinfo(:,2)==stimChn_orig(chnstim)) & (trialinfo(:,18)==AMP(current)));
                                [trialsAmpmatchindex]=find(trialinfo(trialsAmpmatch+1,2)==0);
                                trialmatch=trialinfo(trialsAmpmatch(trialsAmpmatchindex),1);
                                ID=trialmatch;
                                checkID=['T' num2str(ID)];
                                E_MAP = Depth(1);
                                IDstructID=IDstruct.(checkID)(E_MAP,:);
                                baslinespikestructID=baslinespikestruct.(checkID)(E_MAP,:);
                                for numelect=1:size(basespk_all.(check),1)
                                    [psigchan,hsigchan,~] = signrank(baslinespikestructID(numelect*shank_orig,:), IDstructID(numelect*shank_orig,:),'alpha',0.05,'tail','left');
                                    sepdistsig.(sepcheck).(currcheck)(numelect+(16-stimChn(1)),(counterstimchn*2)-(chnstim-1))=nansum([sepdistsig.(sepcheck).(currcheck)(numelect+(16-stimChn(1)),(counterstimchn*2)-(chnstim-1)) hsigchan]);
                                end
                            end
                        end
                        
                        
                    elseif ((isempty(uniquestimchn) || all(uniquestimchn~=1)) || (isempty(uniquestimchn)|| all(uniquestimchn~=2))) && (trial~=1 && trial~=5)
                        if exist('baslinespikestruct')==1
                            for current = 1:length(AMP)
                                if AMP(current)==0 || ~any(Ampcheck==AMP(current))
                                    continue
                                end
                                currcheck=['C' num2str(AMP(current))];
                                trialsAmpmatch=find((trialinfo(:,2)==stimChn_orig(1)) & (trialinfo(:,18)==(trial-1)*0.25*AMP(current)))+1;
                                trialsAmpmatch2=find((trialinfo(:,2)==stimChn_orig(2)) & (trialinfo(:,18)==(5-trial)*0.25*AMP(current)));
                                trialsAmpmatchindex=intersect(trialsAmpmatch,trialsAmpmatch2);
                                trialmatch=trialinfo(trialsAmpmatchindex,1);
                                ID=trialmatch;
                                checkID=['T' num2str(ID)];
                                E_MAP = Depth(1);
                                IDstructID=IDstruct.(checkID)(E_MAP,:);
                                baslinespikestructID=baslinespikestruct.(checkID)(E_MAP,:);
                                for numelect=1:size(basespk_all.(check),1)
                                    [psigchan,hsigchan,~] = signrank(baslinespikestructID(numelect*shank_orig,:), IDstructID(numelect*shank_orig,:),'alpha',0.05,'tail','left');
                                    chnstim=1;
                                    sepdistsig_dual.(sepcheck).(currcheck)(numelect+(16-stimChn(chnstim)),(counterstimchn*3)-(trial-2))=nansum([sepdistsig_dual.(sepcheck).(currcheck)(numelect+(16-stimChn(chnstim)),(counterstimchn*3)-(trial-2)) hsigchan]);
                                end
                            end
                        end
                    end
                    
                    % spike rate and centroid calculation
                    check1=['T' num2str(trial)];
                    for current = 1:length(AMP)
                        currcheck=['C' num2str(AMP(current))];
                        basesplit.(currcheck).(check1)=[basesplit.(currcheck).(check1) [NaN(16-stimChn(1),1); basespk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        Csplit.(currcheck).(check1)=[Csplit.(currcheck).(check1) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        CsplitLayers.(layers_stim).(currcheck).(check1)=[CsplitLayers.(layers_stim).(currcheck).(check1) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        binnum=['B' num2str(stimchn_bins(1))];
                        Csplitbinlayers.(currcheck).(check1).(binnum)=[Csplitbinlayers.(currcheck).(check1).(binnum) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        
                        if shanknum==1
                            %%%%%%%%%calculate theoretical weightings

                            %1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
                            minimum_current=1;%uA if electrode is touching the axon
                            %note that current dissipates and E-fields never
                            %interact; 34um max dist until current dissipates
                            stimchn1_current=currentdistributions(trial,1)*AMP(current);
                            stimchn2_current=currentdistributions(trial,2)*AMP(current);

                            %model activity

                            stimchn1depth=Depthestimate.(['s' num2str(shanknum)])(stimChn(1));%um
                            stimchn2depth=Depthestimate.(['s' num2str(shanknum)])(stimChn(2));%um


                            %1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758

                            layers_depthdefinition=[368 154+368 718+154+368];
                            Neuron_densities=[87533.33, 177300, 107700]./10^9;%microcircuit
                            layerstart=[0 368 368+154];
                            middlelayerpoint=((layers_depthdefinition-layerstart)/2)+layerstart;
                            layer_thickness=(layers_depthdefinition-layerstart);
                            %note that current dissipates and E-fields never
                            %interact; 34um max dist until current dissipates
                            percentageactivated=0.01:0.01:1;%0 to 100 for use later
                            ALL_k_const=interp1([1 50 100],[2100 8850 27500],1:100,'spline');% 0 25 50 75 100% of neurons activated 
                            k_const=ALL_k_const(1);%constant ua/mm^2, const at 0-1%neurons activated will multiply primary activation by ALL-k-const probability function
                            if stimchn1_current>=minimum_current
                                radius_activated1=((stimchn1_current-minimum_current)/k_const)^0.5;
                                radius_activated1=radius_activated1*10.^3;
                                radius_probability1=(((stimchn1_current-minimum_current)./ALL_k_const).^0.5).*10^3;
                            else
                                radius_activated1=0;
                                radius_probability1=0;
                            end
                            if stimchn2_current>=minimum_current
                                radius_activated2=((stimchn2_current-minimum_current)/k_const)^0.5;
                                radius_activated2=radius_activated2*10.^3;
                                radius_probability2=(((stimchn2_current-minimum_current)./ALL_k_const).^0.5).*10^3;
                            else
                                radius_activated2=0;
                                radius_probability2=0;
                            end
                            %2. Is radius confined to one group of layers? -
                            %use stim elect positions
                            %need to define distance based on histology and
                            %microdrive depth and LFP

                            UpperLstimchn1=stimchn1depth-radius_activated1;
                            LowerLstimchn1=stimchn1depth+radius_activated1;
                            UpperLstimchn2=stimchn2depth-radius_activated2;
                            LowerLstimchn2=stimchn2depth+radius_activated2;



                            percentagelayers1u=[1-(layers_depthdefinition-UpperLstimchn1)./(layers_depthdefinition-layerstart) -1];
                            percentagelayers1l=[1-(layers_depthdefinition-LowerLstimchn1)./(layers_depthdefinition-layerstart) -1];
                            percentagelayers2u=[1-(layers_depthdefinition-UpperLstimchn2)./(layers_depthdefinition-layerstart) -1];
                            percentagelayers2l=[1-(layers_depthdefinition-LowerLstimchn2)./(layers_depthdefinition-layerstart) -1];


                            firstupper1=find(percentagelayers1u<0,1, 'first')-1;
                            firstlower1=find(percentagelayers1l<0,1, 'first')-1;
                            firstupper2=find(percentagelayers2u<0,1, 'first')-1;
                            firstlower2=find(percentagelayers2l<0,1, 'first')-1;



                            %3. what is the percentage of each activated layer in
                            %the microcircuit column? i.e. neurons activated
                            %only through stimulation - primary
                            %0.29mm^3 volume in the circuit 2082um long
                            radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
                            area_activated=[0 0 0];
                            volumelayers= (layers_depthdefinition-layerstart).* (pi.*(radius_circuit).^2);
                            volume_upper=[];
                            volume_lower=[];
%                             %primary connections activated
%                             if firstupper1~=firstlower1 %cut by layer boundary
%                                 heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
%                                 if heightsegment>radius_activated1 %reverse cap calc
%                                     volume_lower=((pi*(2*radius_activated1-heightsegment)^2)/3)*(3*radius_activated1-(2*radius_activated1-heightsegment));
%                                     volume_upper=((4/3) * pi*(radius_activated1^3))-volume_lower;
%                                 else % calc cap
%                                     volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
%                                     volume_lower=((4/3) * pi*(radius_activated1^3))-volume_upper;
%                                 end
%                                 elect1_primaryactivated=volume_upper*Neuron_densities(firstupper1);
%                                 elect1_primaryactivated=elect1_primaryactivated+volume_lower*Neuron_densities(firstlower1);
%                                 area_activated(firstupper1)=volume_upper/volumelayers(firstupper1);%percentage volume activated;
%                                 area_activated(firstlower1)=volume_lower/volumelayers(firstlower1);%percentage volume activated;
%                             else %whole volume in one layer group
%                                 area_activated(firstupper1)=((4/3) * pi*(radius_activated1^3))/volumelayers(firstupper1);%percentage volume activated
%                                 elect1_primaryactivated=((4/3) * pi*(radius_activated1^3))*Neuron_densities(firstupper1);
%                             end

                            if firstupper1~=firstlower1 %cut by layer boundary
                                heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
                                if heightsegment>radius_activated1 %reverse cap calc
                                    for iterate_sphereprobabilitydist=1:length(radius_probability1)
                                        heightsegment=layers_depthdefinition(firstupper1)-stimchn1depth+radius_probability1(iterate_sphereprobabilitydist);
                                        if heightsegment<(1*radius_probability1(iterate_sphereprobabilitydist))
                                            volume_lower(iterate_sphereprobabilitydist)=((pi*((2*radius_probability1(iterate_sphereprobabilitydist)-heightsegment)^2)/3)*(3*radius_probability1(iterate_sphereprobabilitydist)-(2*radius_probability1(iterate_sphereprobabilitydist)-heightsegment)));
                                            volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3))-volume_lower(iterate_sphereprobabilitydist);
                                        else
                                            volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3));
                                        end
                                    end
                                else % calc cap
                                    for iterate_sphereprobabilitydist=1:length(radius_probability1)
                                        heightsegment=layers_depthdefinition(firstupper1)-stimchn1depth+radius_probability1(iterate_sphereprobabilitydist);
                                        if heightsegment>0
                                            volume_upper(iterate_sphereprobabilitydist)=((pi*heightsegment^2)/3)*(3*radius_probability1(iterate_sphereprobabilitydist)-heightsegment);
                                            volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3))-volume_upper(iterate_sphereprobabilitydist);
                                        else
                                            volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3));
                                        end
                                    end
                                end
                                volume_upper_percentage=(diff([volume_upper 0]).*-1).*percentageactivated(1:length(volume_upper));%volume activated 
                                volume_lower_percentage=(diff([volume_lower 0]).*-1).*percentageactivated(1:length(volume_lower));%volume activated 
                                 elect1_primaryactivated=sum(volume_upper_percentage)*Neuron_densities(firstupper1);
                                 elect1_primaryactivated=elect1_primaryactivated+sum(volume_lower_percentage)*Neuron_densities(firstlower1);
                                area_activated(firstupper1)=area_activated(firstupper1)+(sum(volume_upper_percentage)/volumelayers(firstupper1));%percentage volume activated;;
                                area_activated(firstlower1)=area_activated(firstlower1)+(sum(volume_lower_percentage)/volumelayers(firstlower1));%percentage volume activated;;
                            else %whole volume in one layer group
                                volume_whole_layer=sum(diff([(4/3).* pi.*(radius_probability1.^3) 0]).*-1.*percentageactivated);
                                area_activated(firstupper1)=area_activated(firstupper1)+((volume_whole_layer)/volumelayers(firstupper1));%percentage volume activated;;
                                elect1_primaryactivated=volume_whole_layer*Neuron_densities(firstupper1);
                            end

                            volume_upper=[];
                            volume_lower=[];
                            %primary connections activated
                            if firstupper2~=firstlower2 %cut by layer boundary
                                heightsegment=layers_depthdefinition(firstupper2)-UpperLstimchn2;
                                if heightsegment>radius_activated2 %reverse cap calc
                                    for iterate_sphereprobabilitydist=1:length(radius_probability2)
                                        heightsegment=layers_depthdefinition(firstupper2)-stimchn2depth+radius_probability2(iterate_sphereprobabilitydist);
                                        if heightsegment<(2*radius_probability2(iterate_sphereprobabilitydist))
                                            volume_lower(iterate_sphereprobabilitydist)=((pi*((2*radius_probability2(iterate_sphereprobabilitydist)-heightsegment)^2)/3)*(3*radius_probability2(iterate_sphereprobabilitydist)-(2*radius_probability2(iterate_sphereprobabilitydist)-heightsegment)));
                                            volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3))-volume_lower(iterate_sphereprobabilitydist);
                                        else
                                            volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3));
                                        end
                                    end
                                else % calc cap
                                    for iterate_sphereprobabilitydist=1:length(radius_probability2)
                                        heightsegment=layers_depthdefinition(firstupper2)-stimchn2depth+radius_probability2(iterate_sphereprobabilitydist);
                                        if heightsegment>0
                                            volume_upper(iterate_sphereprobabilitydist)=((pi*heightsegment^2)/3)*(3*radius_probability2(iterate_sphereprobabilitydist)-heightsegment);
                                            volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3))-volume_upper(iterate_sphereprobabilitydist);
                                        else
                                            volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3));
                                        end
                                    end
                                end
                                
                                volume_upper_percentage=(diff([volume_upper 0]).*-1).*percentageactivated(1:length(volume_upper));%volume activated 
                                volume_lower_percentage=(diff([volume_lower 0]).*-1).*percentageactivated(1:length(volume_lower));%volume activated 
                                 elect2_primaryactivated=sum(volume_upper_percentage)*Neuron_densities(firstupper2);
                                 elect2_primaryactivated=elect2_primaryactivated+sum(volume_lower_percentage)*Neuron_densities(firstlower2);
                                area_activated(firstupper2)=area_activated(firstupper2)+(sum(volume_upper_percentage)/volumelayers(firstupper2));%percentage volume activated;;
                                area_activated(firstlower2)=area_activated(firstlower2)+(sum(volume_lower_percentage)/volumelayers(firstlower2));%percentage volume activated;;
                            else %whole volume in one layer group
                                volume_whole_layer=sum(diff([(4/3).* pi.*(radius_probability2.^3) 0]).*-1.*percentageactivated);
                                area_activated(firstupper2)=area_activated(firstupper2)+((volume_whole_layer)/volumelayers(firstupper2));%percentage volume activated;;
                                elect2_primaryactivated=volume_whole_layer*Neuron_densities(firstupper2);
                            end

                            %4. Based on the percentage of primary, what is the
                            %percentage activated through the secondary
                            %connections?
                            %overlappping or non-overlappping activation?
                            %second

                            connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
                            Secondary_neurons_activated=ceil(connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
                            baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);
                            tertiary_connections=sum((Secondary_neurons_activated./(Neuron_densities'.*volumelayers'))'.*Numconnections_3groups_array,2);
                            tertiary_neurons_activated=ceil(tertiary_connections.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
                            SecondaryTertiary_neurons_activated=Secondary_neurons_activated;%+tertiary_neurons_activated;

                            %5. weightings based primary and
                            %secondary or secondary only
                            %Neuron_densities=[14200, 83800, 164600, 83900,
                            %177300, 131500]; split into layers
                            %Neuron_number_groupedlayers=[338+7524 4656 6114+12651].*area_activated;
                            Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);

                            neurons_distributed_depth=zeros(ceil(layers_depthdefinition(end)/50),1);
%                             neurons_distributed_depth(round(stimchn1depth/50)-ceil(radius_activated1/50):round(stimchn1depth/50)+ceil(radius_activated1/50))=ceil(elect1_primaryactivated/(ceil(radius_activated1*2/50)+1));
%                             neurons_distributed_depth(round(stimchn2depth/50)-ceil(radius_activated2/50):round(stimchn2depth/50)+ceil(radius_activated2/50))=ceil(elect2_primaryactivated/(ceil(radius_activated2*2/50)+1));

                            for layertypes=1:3
                                midpoint=(middlelayerpoint(layertypes));
                                numneurons=SecondaryTertiary_neurons_activated(layertypes);
                                if numneurons<2
                                    neurons_distributed_depth(ceil(midpoint/50))=neurons_distributed_depth(ceil(midpoint/50))+SecondaryTertiary_neurons_activated(layertypes);
                                    neuronposition=ceil(midpoint/50);
                                elseif rem(numneurons,2)%odd
                                    distbetweenneurons=(layer_thickness(layertypes)/numneurons);
                                    neuronposition=(midpoint-(distbetweenneurons)*((numneurons-1)/2)):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons-1)/2);
                                    neuronposition=ceil(neuronposition/50);
                                    neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
                                elseif rem(numneurons+1,2)%even
                                    distbetweenneurons=(layer_thickness(layertypes)/numneurons);
                                    neuronposition=(midpoint-(distbetweenneurons)*((numneurons)/2)+distbetweenneurons/2):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons)/2);
                                    neuronposition=ceil(neuronposition/50);
                                    neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
                                end

                                %     if length(unique(neuronposition))~=length(neuronposition)
                                %         error('sampling resolution too large')
                                %     end
                            end



                            weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);

                            %bursting activity due to extracellular stim https://www.sciencedirect.com/science/article/pii/S1935861X09000424#app1

                            burstprob=1./8;%assume cell might fire twice in 8ms window and record for 6ms(2-8)
                            firingrate=neurons_distributed_depth.*burstprob.*1000;
                            filtmov=1/5*ones(5,1);
                            firingrate=filtfilt(filtmov,1,firingrate);
                            [~,shankposlayers]=min(abs([50:50:1250]-Depthestimate.s2(end)));
                            simulationFR=NaN(16,1);
                            simulationFR(1:length(firingrate(shankposlayers:end)))=firingrate(shankposlayers:end);
                            simulationFR=flip(simulationFR(1:16));

                            Csplit_simulation.(currcheck).(check1)=[Csplit_simulation.(currcheck).(check1) [NaN(16-stimChn(1),1); simulationFR; NaN(32-(16-stimChn(1))-16,1)]];
                             CsplitLayers_model.(layers_stim).(currcheck).(check1)=[CsplitLayers_model.(layers_stim).(currcheck).(check1) [NaN(16-stimChn(1),1); simulationFR; NaN(32-(16-stimChn(1))-16,1)]];
                        end



                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
%                         k_const=8850;%constant ua/mm^2
%                         if stimchn1_current>=minimum_current
%                             radius_activated1=((stimchn1_current-minimum_current)/k_const)^0.5;
%                             radius_activated1=radius_activated1*10.^3;
%                         else
%                             radius_activated1=0;
%                         end
%                         if stimchn2_current>=minimum_current
%                             radius_activated2=((stimchn2_current-minimum_current)/k_const)^0.5;
%                             radius_activated2=radius_activated2*10.^3;
%                         else
%                             radius_activated2=0;
%                         end
%                         
%                         
%                         %2. Is radius confined to one group of layers? -
%                         %use stim elect positions
%                         %need to define distance based on histology and
%                         %microdrive depth and LFP
%                         Shankdepthcheck=['s' num2str(shank)];
%                         UpperLstimchn1=Depthestimate.(Shankdepthcheck)(stimChn(1))-radius_activated1;
%                         LowerLstimchn1=Depthestimate.(Shankdepthcheck)(stimChn(1))+radius_activated1;
%                         UpperLstimchn2=Depthestimate.(Shankdepthcheck)(stimChn(2))-radius_activated2;
%                         LowerLstimchn2=Depthestimate.(Shankdepthcheck)(stimChn(2))+radius_activated2;
% 
%                         
%                         layers_thickness=[368 154 718];
% %                         layerstim1=ElectLayerClass(stimChn_orig(1));
% %                         layerstim2=ElectLayerClass(stimChn_orig(2));
% %                         percentagelayers1u=[1-(layers_depthdefinition-UpperLstimchn1)./(layers_depthdefinition-[0 667 857]) -1];
% %                         percentagelayers1l=[1-(layers_depthdefinition-LowerLstimchn1)./(layers_depthdefinition-[0 667 857]) -1];
% %                         percentagelayers2u=[1-(layers_depthdefinition-UpperLstimchn2)./(layers_depthdefinition-[0 667 857]) -1];
% %                         percentagelayers2l=[1-(layers_depthdefinition-LowerLstimchn2)./(layers_depthdefinition-[0 667 857]) -1];
% % 
% %                         firstupper1=find(percentagelayers1u<0,1, 'first')-1;
% %                         firstlower1=find(percentagelayers1l<0,1, 'first')-1;
% %                         firstupper2=find(percentagelayers2u<0,1, 'first')-1;
% %                         firstlower2=find(percentagelayers2l<0,1, 'first')-1;
% 
%                         firstupper1=ElectLayerClass(stimChn_orig(1));
%                         firstlower1=ElectLayerClass(stimChn_orig(1));
%                         firstupper2=ElectLayerClass(stimChn_orig(2));
%                         firstlower2=ElectLayerClass(stimChn_orig(2));
% 
%                         
%                     
%                         
%                         
%                         %3. what is the percentage of each activated layer in
%                         %the microcircuit column? i.e. neurons activated
%                         %only through stimulation - primary
%                         %0.29mm^3 volume in the circuit 2082um long
%                         radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
% %                         if radius_activated>radius_circuit
% %                             error('radius bigger than microcircuit')
% %                         elseif radius_activated>190%layer4size
% %                             error('radius bigger than L4')
% %                         end
% 
%                         area_activated=[0 0 0];
%                         volumelayers= (layers_thickness).* (pi.*(radius_circuit).^2);
%                         
%                         %primary connections activated
%                         if firstupper1~=firstlower1 %cut by layer boundary
%                             error('should not reach')
%                             heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
%                             if heightsegment>radius_activated1 %reverse cap calc
%                                 volume_lower=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
%                                 volume_upper=((4/3) * pi*(radius_activated1^3))-volume_lower;
%                             else % calc cap
%                                 volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
%                                 volume_lower=((4/3) * pi*(radius_activated1^3))-volume_upper;
%                             end
%                             area_activated(firstupper1)=volume_upper/volumelayers(firstupper1);%percentage volume activated;
%                             area_activated(firstlower1)=volume_lower/volumelayers(firstlower1);%percentage volume activated;
%                         else %whole volume in one layer group
%                             area_activated(firstupper1)=((4/3) * pi*(radius_activated1^3))/volumelayers(firstupper1);%percentage volume activated
%                         end
% 
%                         %primary connections activated
%                         if firstupper2~=firstlower2 %cut by layer boundary
%                             error('should not reach')
%                             heightsegment=layers_depthdefinition(firstupper2)-UpperLstimchn2;
%                             if heightsegment>radius_activated2 %reverse cap calc
%                                 volume_lower=((pi*heightsegment^2)/3)*(3*radius_activated2-heightsegment);
%                                 volume_upper=((4/3) * pi*(radius_activated2^3))-volume_lower;
%                                 
%                             else % calc cap
%                                 volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated2-heightsegment);
%                                 volume_lower=((4/3) * pi*(radius_activated2^3))-volume_upper;
%                             end
%                             area_activated(firstupper2)=area_activated(firstupper2)+volume_upper/volumelayers(firstupper2);%percentage volume activated;;
%                             area_activated(firstlower2)=area_activated(firstlower2)+volume_lower/volumelayers(firstlower2);%percentage volume activated;;
%                         else %whole volume in one layer group
%                             area_activated(firstupper2)=area_activated(firstupper2)+((4/3) * pi*(radius_activated2^3))/volumelayers(firstupper2);%percentage volume activated;;
%                         end
%                  
%                         %4. Based on the percentage of primary, what is the
%                         %percentage activated through the secondary
%                         %connections?
%                         %overlappping or non-overlappping activation?
%                         %second
%                         layers_depthdefinition=[353+149+165 190+353+149+165 525+700+190+353+149+165];%thickness of microcircuit
%                        
%                         Neuron_densities=[87533.33, 177300, 107700]./10^9;
%                         connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
%                         Secondary_neurons_activated=connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2);
%                         baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);
%                         tertiary_connections=sum((Secondary_neurons_activated./(Neuron_densities'.*volumelayers'))'.*Numconnections_3groups_array,2);
%                         tertiary_neurons_activated=tertiary_connections.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2);
%                         %5. weightings based primary and
%                         %secondary or secondary only
%                            %Neuron_densities=[14200, 83800, 164600, 83900,
%                            %177300, 131500]; split into layers
%                            %Neuron_number_groupedlayers=[338+7524 4656 6114+12651].*area_activated;
%                            Primary_neurons_activated=Neuron_densities.*area_activated.*volumelayers;
%                            
%                            
%                            weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);
%                            
% 
%                            %6. weighted results
%                            CsplitLayers_weighted.(currcheck).(check1)=[CsplitLayers.(layers_stim).(currcheck).(check1)...
%                                [NaN(16-stimChn(1),1); normspk_all.(check)(:,current).*weightings(ElectLayerClass((shank_orig-1)*16+1:(shank_orig)*16))'; NaN(32-(16-stimChn(1))-16,1)]];
% 
%                            
%                         %%%%%%%%%model
%                         recorddist=50;%um max radius dist electrode can record a neuron - need to adjust
%                          check=['S' num2str(shanknum) '_T' num2str(trial)];
%                          %could use this as a ratio but the neurons
%                          %activated will preferentially be around the
%                          %electrode
%                          %volume_recordedfrom=(4/3)*pi*recorddist*16*4;
%                          %volume_notrecordedfrom=(0.29*10^9)-volume_recordedfrom;
%                          %
%                         neurons_activatedperUM=(Primary_neurons_activated+Secondary_neurons_activated')./(layers_thickness);
%                         spreadrecordedvertical=[sum(ElectLayerClass==1)./4.*recorddist*2 sum(ElectLayerClass==2)./4.*recorddist*2 sum(ElectLayerClass==3)./4.*recorddist*2];
%                         NeuronsRecordedFrom_layer=neurons_activatedperUM.*spreadrecordedvertical;

                        
                        
                        if trialcheck~=1
                            CsplitCentroid.(currcheck)=[CsplitCentroid.(currcheck) centroidpershank((shanknum-1)*5+1:(shanknum*5),current)-stimChn(1)];
                        end
                        % plot resutls of 200um away, 400um away and 600um away
                        Shanksep=shanknum-shank;
                        bincheck=['D' num2str(abs(Shanksep))];
                        bin.(currcheck).(check1).(bincheck)=[bin.(currcheck).(check1).(bincheck) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        %                         if strcmp(bincheck,binchecksave)
                        %                             bin.(currcheck).(check1).(bincheck)(:,end)=mean([bin.(currcheck).(check1).(bincheck)(:,end) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]],2);
                        %                         else
                        %                             bin.(currcheck).(check1).(bincheck)=[bin.(currcheck).(check1).(bincheck) [NaN(16-stimChn(1),1); normspk_all.(check)(:,current); NaN(32-(16-stimChn(1))-16,1)]];
                        %                         end
                    end
                    trialcheck=1;
                end
                if ~strcmp(bincheck,'D0')
                    binchecksave=bincheck;
                end
            end
            binchecksave='blank';
            counternumberofpairs=counternumberofpairs+1;
            countershank=0;
            stimchn_r.(Ratnum)=[stimchn_r.(Ratnum) (stimChn.*shank)];
        end
        savestimshanks(savestimshanks==0)=[];
        stimshank=[stimshank; savestimshanks];
    end
    check=['sep' num2str(-1*sepdist-1)];
    saveCplit.(check)=Csplit;
    savebasesplit.(check)=basesplit;
    saveCplit_weightedlayers.(check)=CsplitLayers_weighted;
    saveCplitlayers.(check)=CsplitLayers;
    saveCsplitcentroid.(check)=CsplitCentroid;
    saveCsplitbinlayers.(check)=Csplitbinlayers;
    saveshanksepdist.(check)=bin;
    saveCplit_simulation.(check)=Csplit_simulation;
    saveCplitlayers_simulation.(check)=CsplitLayers_model;
end


%% N's for everything



%% undo weighting of the responses based on stim layer

%Numconnections_3groups.(L_check)./totalconnect_L.all
 CsplitLayers
 %needs to go into loop^^ split into primary, primary&secondary, secondary
 %- need to determine area activated within the mmicrocircuit layers -
 % % activated -> density -> % connections ->weights 
 
 %connection probability??
 
 %https://www.sciencedirect.com/science/article/pii/0165027096000659 -
 %8280uA/mm^2 activation average of intrinsic both above and below
 I=5;%current in ua
 k=8280;%constant ua/mm^2
 r=(I/k)^0.5;
 
 %% normalise data by max

AMP=[0 1 2 3 4 6 8 10];
 for sepdist=5:2:9
     sepcheck=['sep' num2str(sepdist)];
     for current=1:length(AMP)
         currcheck=['C' num2str(AMP(current))];
%         all_sim_curr= cat(3,saveCplit_simulation.(sepcheck).(currcheck).T1, saveCplit_simulation.(sepcheck).(currcheck).T2, saveCplit_simulation.(sepcheck).(currcheck).T3, saveCplit_simulation.(sepcheck).(currcheck).T4, saveCplit_simulation.(sepcheck).(currcheck).T5);
      %   max_sim=max(all_sim_curr,[],3);
         all_dat_curr=cat(3,saveCplit.(sepcheck).(currcheck).T1, saveCplit.(sepcheck).(currcheck).T2, saveCplit.(sepcheck).(currcheck).T3, saveCplit.(sepcheck).(currcheck).T4, saveCplit.(sepcheck).(currcheck).T5);
         max_data=max(all_dat_curr,[],3);
       %  D1_data=cat(3,saveshanksepdist.(sepcheck).(currcheck).T1.D1,saveshanksepdist.(sepcheck).(currcheck).T2.D1,saveshanksepdist.(sepcheck).(currcheck).T3.D1,saveshanksepdist.(sepcheck).(currcheck).T4.D1,saveshanksepdist.(sepcheck).(currcheck).T5.D1);
      %   max_D1_data=max(D1_data,[],3);
         for trial=trials
             check=['T' num2str(trial)];
          %   normalise_Csplit_simulation.(sepcheck).(currcheck).(check)=saveCplit_simulation.(sepcheck).(currcheck).(check)./max_sim;
             normalise_Csplit.(sepcheck).(currcheck).(check)=saveCplit.(sepcheck).(currcheck).(check)./max_data;
          %   normalise_Csplit_Ddat.(sepcheck).(currcheck).(check)=saveshanksepdist.(sepcheck).(currcheck).(check).D1./max_D1_data;
         end
     end
 end

 %% normalise data by baseline

AMP=[0 1 2 3 4 6 8 10];
 for sepdist=5:2:9
     sepcheck=['sep' num2str(sepdist)];
     for current=1:length(AMP)
         currcheck=['C' num2str(AMP(current))];
         all_dat_curr=cat(3,savebasesplit.(sepcheck).(currcheck).T1, savebasesplit.(sepcheck).(currcheck).T2, savebasesplit.(sepcheck).(currcheck).T3, savebasesplit.(sepcheck).(currcheck).T4, savebasesplit.(sepcheck).(currcheck).T5);
         mean_data=mean(all_dat_curr,3,'omitnan');
         for trial=trials
             check=['T' num2str(trial)];
             normalise_Csplit.(sepcheck).(currcheck).(check)=saveCplit.(sepcheck).(currcheck).(check)./mean_data;
         end
     end
 end

%% normalise data to 50:50 condition
AMP=[0 1 2 3 4 6 8 10];
 for sepdist=5:2:9
     sepcheck=['sep' num2str(sepdist)];
     for current=1:length(AMP)
         currcheck=['C' num2str(AMP(current))];
         sim_5050=saveCplit_simulation.(sepcheck).(currcheck).T3;
         data_5050=saveCplit.(sepcheck).(currcheck).T3;
         D1_data_5050=saveshanksepdist.(sepcheck).(currcheck).T3.D1;
         data_5050(data_5050<0.001)=nan;
         for trial=trials
             check=['T' num2str(trial)];
             normalise_Csplit_simulation.(sepcheck).(currcheck).(check)=saveCplit_simulation.(sepcheck).(currcheck).(check)./sim_5050;
             normalise_Csplit.(sepcheck).(currcheck).(check)=saveCplit.(sepcheck).(currcheck).(check)./data_5050;
             normalise_Csplit_Ddat.(sepcheck).(currcheck).(check)=saveshanksepdist.(sepcheck).(currcheck).(check).D1./D1_data_5050;
         end
     end
 end

 %% plot split bin layers
singleCurrent=6; %plot 6uA results
numshanksToAvg=1;
N=128;
vec = [100;80;50;30;15;0];
     hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
     raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
     map = interp1(vec,raw,linspace(100,0,N),'pchip');
     cmap=colormap(map);
 if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end
 
 for sepdist=5:2:9
     sepcheck=['sep' num2str(sepdist)];

    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for binx=0:5
            binnum=['B' num2str(binx)];
            figure;
            for trial=1:5
                check=['T' num2str(trial)];
                dat=saveCsplitbinlayers.(sepcheck).(currcheck).(check).(binnum);
                if isempty(dat)
                    continue
                end
                for paircount=1:size(dat,2)/numshanksToAvg
                    pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*numshanksToAvg+1:numshanksToAvg* paircount),2);
                end
                hold on
                stdshade(pairavgcur',0.2,cmap(trial*floor((length(cmap))/5),:));
                title(binnum)
            end
            if isempty(dat)
                close;
            else 
                 ylim([0, round(max(pairavgcur,[],'all')/50)*50+50])
                %ylim([0 1])
%                 if sepdist==5
%                     xlim([16-3 16+sepdist+1+3])
%                 elseif sepdist==7
%                     xlim([16-2 16+sepdist+1+2])
%                 elseif sepdist==9
%                     xlim([16-1 16+sepdist+1+1])
%                 end
                xline(16,'r')
                xline(16+sepdist+1,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                yt = yticks;
                yticklabels(yt(1:end-1));
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
                set(gca,'TickDir','out');
                ax1=gca;
                legend('0:100','25:75','50:50','75:25','100:0')
            end
        end
    end
end
%% Plotting main result - single shank or multiple
plotsimulation=0;%simulation
plotnormalised=1;%normalise data^^^
plotAllShanks=1; %1=yes, 0=no. 0 only plots shanks directly next to stim shank
numshanksToAvg=4;%if ^^ is yes. if you take out stim shank ==3. if all shanks ==4
downsampleYN=1; %do you want to downsample 1=yes
singleCurrent=6; %plot 6uA results
ploty=1; %Do the plotting
splitlayersYN=0;%do you want to split into layers

if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end
vec = [100;80;50;30;15;0];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
centroidpos_all=[];
peak_all=[];
samples_rand=[];
singleelectspread=[];
dualelectspread=[];
trials=[1 2 3 4 5];
if splitlayersYN==1
    layer_labels={'ss' 'gg' 'ii' 'sg' 'gi' 'si'};
    numlayers=6;
else
    layer_labels={''};
    numlayers=1;
end

clear spread_electno spread_electnomaxall sepdistsigdwnsampe sepdistsigdwnsampe_dual pairavgdwnsample_all_nomean
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    dualelectspreadsepdist.(sepcheck)=[];
    
    hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
    raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
    map = interp1(vec,raw,linspace(100,0,N),'pchip');
    cmap=colormap(map);
    for layernum=1:numlayers
        layercheck=layer_labels{layernum};
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            if length(AMP)>1 && ploty==1
                figure(1+sepdist+layernum*10)
                subplot(2,4,current)
                hold on
            end
            savdwnsamplespkrate=zeros(32,5);
            if ploty==1 && length(AMP)==1
                figure(1+sepdist+layernum*10)
                axes('Position',[0.13         0.112396822842341                     0.775         0.62])
                hold on
                set(gca,'FontSize',12)
            end
            for trial=trials
                check=['T' num2str(trial)];
                if plotsimulation==1 && plotnormalised==1 
                    DataInput2Plot=normalise_Csplit_simulation.(sepcheck).(currcheck).(check);
                    downsampleYN=0;
                elseif plotnormalised==1 &&  splitlayersYN==0 && (plotAllShanks==1) && plotsimulation==0
                    DataInput2Plot=normalise_Csplit.(sepcheck).(currcheck).(check);%saveCplitlayers.(sepcheck).(layercheck).(currcheck).(check);% single shank - saveshanksepdist.(sepcheck).(currcheck).(check).D1; all shanks - saveCplit.(sepcheck)
                    lengthneeded=size(normalise_Csplit.(sepcheck).C1.T1,2)./numshanksToAvg;%saveCplitlayers.(sepcheck).(layercheck).C1.T1 find the smallest number of pairs for any current - we need to down-sample the others to match
                elseif plotnormalised==1 &&  splitlayersYN==0 && (plotAllShanks==0) && plotsimulation==0
                    DataInput2Plot=normalise_Csplit_Ddat.(sepcheck).(currcheck).(check);%saveCplitlayers.(sepcheck).(layercheck).(currcheck).(check);% single shank - saveshanksepdist.(sepcheck).(currcheck).(check).D1; all shanks - saveCplit.(sepcheck)
                    lengthneeded=size(normalise_Csplit_Ddat.(sepcheck).C1.T1,2)./numshanksToAvg;%saveCplitlayers.(sepcheck).(layercheck).C1.T1 find the smallest number of pairs for any current - we need to down-sample the others to match
                elseif plotsimulation==1 && splitlayersYN==0
                    DataInput2Plot=saveCplit_simulation.(sepcheck).(currcheck).(check);
                    downsampleYN=0;
                elseif plotsimulation==1 && splitlayersYN==1
                    DataInput2Plot=saveCplitlayers_simulation.(sepcheck).(layercheck).(currcheck).(check);
                    downsampleYN=0;
                elseif (plotAllShanks==1) &&  splitlayersYN==1 %are we splitting based on LFP
                    DataInput2Plot=saveCplitlayers.(sepcheck).(layercheck).(currcheck).(check);% single shank - saveshanksepdist.(sepcheck).(currcheck).(check).D1; all shanks - saveCplit.(sepcheck)
                    lengthneeded=size(saveCplitlayers.(sepcheck).(layercheck).C1.T1,2)./numshanksToAvg;% find the smallest number of pairs for any current - we need to down-sample the others to match
                elseif (plotAllShanks==1) &&  splitlayersYN==0
                    DataInput2Plot=saveCplit.(sepcheck).(currcheck).(check);%saveCplitlayers.(sepcheck).(layercheck).(currcheck).(check);% single shank - saveshanksepdist.(sepcheck).(currcheck).(check).D1; all shanks - saveCplit.(sepcheck)
                    lengthneeded=size(saveCplit.(sepcheck).C1.T1,2)./numshanksToAvg;%saveCplitlayers.(sepcheck).(layercheck).C1.T1 find the smallest number of pairs for any current - we need to down-sample the others to match

                else
                    DataInput2Plot=saveshanksepdist.(sepcheck).(currcheck).(check).D1;
                    lengthneeded=size(saveshanksepdist.(sepcheck).C1.T1.D1,2);% find the smallest number of pairs for any current - we need to down-sample the others to match
                end
                dat=DataInput2Plot;
                if isempty(dat)
                    continue
                end
                dat(isinf(dat))=nan;
                if plotAllShanks==1
                    clear pairavgcur erpairavgcur
                    for paircount=1:size(dat,2)/numshanksToAvg
                        pairavgcur(:,paircount)=nanmean(dat(:,( paircount-1)*numshanksToAvg+1:numshanksToAvg* paircount),2);
                    end
                else
                    pairavgcur=dat;
                end
                
                if downsampleYN==1
                    if trial==1
                        [samples_rand1]=DownSample(pairavgcur,lengthneeded,s,seedpoint);
                        samples_rand.(sepcheck)=samples_rand1;
                    end
                    pairavg_dwnsample=pairavgcur(:, samples_rand.(sepcheck));% downsample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else
                    pairavg_dwnsample=pairavgcur;
                end
                % pairavg_dwnsample(sum(~isnan(pairavg_dwnsample),2)<10,:)=nan; %remove rows with less than 10 pairs in the average
                pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check)=pairavg_dwnsample;
                pairavgdwnsample_all.(sepcheck).(currcheck)(:,trial)=nanmean(pairavg_dwnsample,2);
                %             pairavgdwnsample_all.(sepcheck).(currcheck)(:,trial)=nanmean(dat,2);
                %             pairavg_dwnsample=dat;
                if ploty==1
                    stdshade(pairavg_dwnsample',0.2,cmap(trial*floor((length(cmap))/5),:));
                end
                check=['T' num2str(trial)];
            end
            
            if isempty(dat)
                continue
            end
            if ploty==1
                ylim([0, round(max(pairavgdwnsample_all.(sepcheck).(currcheck),[],'all')/50)*50+50])
                %ylim([0 1])
%                 if sepdist==5
%                     xlim([16-3 16+sepdist+1+3])
%                 elseif sepdist==7
%                     xlim([16-2 16+sepdist+1+2])
%                 elseif sepdist==9
%                     xlim([16-1 16+sepdist+1+1])
%                 end
                xline(16,'r')
                xline(16+sepdist+1,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                yt = yticks;
                yticklabels(yt(1:end-1));
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
                set(gca,'TickDir','out');
                ax1=gca;
                legend('0:100','25:75','50:50','75:25','100:0')
            end
        end
   
    if ploty==1
        set(gca,'TickDir','out');
        hex = ['#1b0c36';'#532e9e';'#147df5';'#02a612';'#fffb7d';'#ffe26e'];
        raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
        map = interp1(vec,raw,linspace(100,0,N),'pchip');
        cmap=colormap(map);
    end
    % centroid per pair %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pcur=[];
            if isempty(dat)
                continue
            end
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            %         for trial=trials
            %             check=['T' num2str(trial)];
            %             dat=pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check);
            %             dat(isinf(dat))=nan;
            %             dat(sum(~isnan(dat),2)<3,:)=nan;
            %             pairavgcur=zeros(size(dat,1),size(dat,2));
            %             clear pairavgcur erpairavgcur
            %             for paircount=1:size(dat,2)
            %                 pairavgcur(:,paircount)=dat(:,paircount);
            %                 [electrodecentroid, mincentroidstart]=centroidFspiking(pairavgcur(:,paircount),sepdist);
            %                 centroidpos_all.(sepcheck).(currcheck)(trial,paircount)=(electrodecentroid-1-mincentroidstart).*50;%(centroidpos-16).*50;%(electrodecentroid-1).*50;%
            %             end
            %
            %         end
            
            %%%%%%%%%%%%%%%%%smoothed peaks
            for trial=trials
                check=['T' num2str(trial)];
                dat=pairavgdwnsample_all_nomean.(sepcheck).(currcheck).(check);
                dat(isinf(dat))=nan;
                %dat(sum(~isnan(dat),2)<3,:)=nan;
                clear pairavgcur erpairavgcur
                for paircount=1:size(dat,2)
                    datnonan=dat(:,paircount);
                    %                     firstnonan=find(~isnan(datnonan),1,'first');
                    %                     lastnonan=find(~isnan(datnonan),1,'last');
                    endnonan=find(diff(isnan(datnonan))==1);
                    startnonan=find(diff(isnan(datnonan))==-1);
                    [~,p]=max(endnonan-startnonan);
                    lastnonan=endnonan(p);
                    firstnonan=startnonan(p)+1;
                    datnonan=datnonan(firstnonan:lastnonan);
                    if length(datnonan)<7
                        continue
                    end
                    FilterLength=3;
                    b=ones(FilterLength,1)./FilterLength;
                    smoothedenvelope=filtfilt(b,1,datnonan);
                    pairavgcur(:,paircount)=[nan(firstnonan-1,1); smoothedenvelope; nan(32-lastnonan,1)];
                end
                if isempty(datnonan)
                    continue
                end
                pairavgcur(:,sum(pairavgcur>0,1)==0)=nan;
                [peak,c]=find(pairavgcur==max(pairavgcur));
                peak_all.(sepcheck).(currcheck)(trial,c)=(peak'-16).*50;
            end
            
            if ploty==1
%                 figure(1+sepdist+layernum*10)
%                 if sepdist==5
%                     xlim([16-3 16+sepdist+1+3])
%                 elseif sepdist==7
%                     xlim([16-2 16+sepdist+1+2])
%                 elseif sepdist==9
%                     xlim([16-1 16+sepdist+1+1])
%                 end
                xline(16,'r')
                xline(16+sepdist+1,'r')
                xt = xticks;
                xtl=(xt-16)*50;
                xticklabels(xtl)
                ylabel('Sp/s')
                xlabel('Distance from deepest stim elect (\mum)')
                set(gca,'TickDir','out');
                legend('0:100','25:75','50:50','75:25','100:0')
            end
            %%%%%%%%%%%%%%%%%%%%plotting
            %         avgpcurrc=mean(centroidpos_all.(sepcheck).(currcheck),2);
            %         stdavgc=std(centroidpos_all.(sepcheck).(currcheck),[],2)./sqrt(size(centroidpos_all.(sepcheck).(currcheck),2));
            avgpcurr=nanmean(peak_all.(sepcheck).(currcheck),2);
            stdavg=nanstd(peak_all.(sepcheck).(currcheck),[],2)./sqrt(sum(~isnan(peak_all.(sepcheck).(currcheck)),2));
            
            if ploty==1
                avgpcurr=avgpcurr(trials);
                if length(AMP)==1
                    figure(1+sepdist+layernum*10)
                    axes('Position',[ax1.Position(1) .73 ax1.Position(3) .2])
                    ax2=gca;
                    hold on
                    set(gca,'FontSize',12)
                    for colorbar=1:length(avgpcurr)
                        color_current=cmap(trials(colorbar)*floor((length(cmap))/5),:);
                        er = errorbar(avgpcurr(colorbar),colorbar,stdavg(trials(colorbar)),stdavg(trials(colorbar)),'horizontal');er.Color = [0, 0, 0]; er.LineStyle = 'none';
                        scatter(avgpcurr(colorbar), colorbar, [], 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', color_current)
                    end
                    yticks([0 1 2 3 4 5 6 7 8])
                    labels=fliplr({'0:100','25:75','50:50','75:25','100:0'});
                    labels=labels(trials);
                    yticklabels([{''} labels {''}])
                    ylim([0.5 length(trials)+0.5])
                    xlim((ax1.XLim-16).*50)
                    set(gca,'TickDir','out');
                    set(gca,'xtick',[])
                    set(gca,'xticklabel',[])
                    set(gca,'ytick',[])
                    set(gca,'yticklabel',[])
                    xline(0,'r')
                    xline((sepdist+1)*50,'r')
                    title(layercheck)
                end
            end
        end
    end
end





%% single model

%model activity
current=10;%uA
stimchn1depth=200;%um
stimchn2depth=500;%um
trial=3;

%1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
minimum_current=1;%uA if electrode is touching the axon
currentdistributions=flipud([1 0; 0.75 0.25; 0.5 0.5; 0.25 0.75; 0 1]);
%layers_depthdefinition=[353+149+165 190+353+149+165 525+700+190+353+149+165];
layers_depthdefinition=[368 154+368 718+154+368];
Neuron_densities=[87533.33, 177300, 107700]./10^9;%microcircuit
layerstart=[0 368 368+154];
middlelayerpoint=((layers_depthdefinition-layerstart)/2)+layerstart;
layer_thickness=(layers_depthdefinition-layerstart);
%note that current dissipates and E-fields never
%interact; 34um max dist until current dissipates
stimchn1_current=currentdistributions(trial,1)*current;
stimchn2_current=currentdistributions(trial,2)*current;
k_const=8850;%constant ua/mm^2
if stimchn1_current>=minimum_current
    radius_activated1=((stimchn1_current-minimum_current)/k_const)^0.5;
    radius_activated1=radius_activated1*10.^3;
else
    radius_activated1=0;
end
if stimchn2_current>=minimum_current
    radius_activated2=((stimchn2_current-minimum_current)/k_const)^0.5;
    radius_activated2=radius_activated2*10.^3;
else
    radius_activated2=0;
end
%2. Is radius confined to one group of layers? -
%use stim elect positions
%need to define distance based on histology and
%microdrive depth and LFP

UpperLstimchn1=stimchn1depth-radius_activated1;
LowerLstimchn1=stimchn1depth+radius_activated1;
UpperLstimchn2=stimchn2depth-radius_activated2;
LowerLstimchn2=stimchn2depth+radius_activated2;



percentagelayers1u=[1-(layers_depthdefinition-UpperLstimchn1)./(layers_depthdefinition-layerstart) -1];
percentagelayers1l=[1-(layers_depthdefinition-LowerLstimchn1)./(layers_depthdefinition-layerstart) -1];
percentagelayers2u=[1-(layers_depthdefinition-UpperLstimchn2)./(layers_depthdefinition-layerstart) -1];
percentagelayers2l=[1-(layers_depthdefinition-LowerLstimchn2)./(layers_depthdefinition-layerstart) -1];


firstupper1=find(percentagelayers1u<0,1, 'first')-1;
firstlower1=find(percentagelayers1l<0,1, 'first')-1;
firstupper2=find(percentagelayers2u<0,1, 'first')-1;
firstlower2=find(percentagelayers2l<0,1, 'first')-1;



%3. what is the percentage of each activated layer in
%the microcircuit column? i.e. neurons activated
%only through stimulation - primary
%0.29mm^3 volume in the circuit 2082um long
radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
area_activated=[0 0 0];
volumelayers= (layers_depthdefinition-layerstart).* (pi.*(radius_circuit).^2);

%primary connections activated
if firstupper1~=firstlower1 %cut by layer boundary
    heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
    if heightsegment>radius_activated1 %reverse cap calc
        volume_lower=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
        volume_upper=((4/3) * pi*(radius_activated1^3))-volume_lower;
    else % calc cap
        volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
        volume_lower=((4/3) * pi*(radius_activated1^3))-volume_upper;
    end
    elect1_primaryactivated=volume_upper*Neuron_densities(firstupper1);
    elect1_primaryactivated=elect1_primaryactivated+volume_lower*Neuron_densities(firstlower1);
    area_activated(firstupper1)=volume_upper/volumelayers(firstupper1);%percentage volume activated;
    area_activated(firstlower1)=volume_lower/volumelayers(firstlower1);%percentage volume activated;
else %whole volume in one layer group
    area_activated(firstupper1)=((4/3) * pi*(radius_activated1^3))/volumelayers(firstupper1);%percentage volume activated
    elect1_primaryactivated=((4/3) * pi*(radius_activated1^3))*Neuron_densities(firstupper1);
end

%primary connections activated
if firstupper2~=firstlower2 %cut by layer boundary
    heightsegment=layers_depthdefinition(firstupper2)-UpperLstimchn2;
    if heightsegment>radius_activated2 %reverse cap calc
        volume_lower=((pi*heightsegment^2)/3)*(3*radius_activated2-heightsegment);
        volume_upper=((4/3) * pi*(radius_activated2^3))-volume_lower;
        
    else % calc cap
        volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated2-heightsegment);
        volume_lower=((4/3) * pi*(radius_activated2^3))-volume_upper;
    end
    elect2_primaryactivated=volume_upper*Neuron_densities(firstupper1);
    elect2_primaryactivated=elect2_primaryactivated+volume_lower*Neuron_densities(firstlower2);
    area_activated(firstupper2)=area_activated(firstupper2)+(volume_upper/volumelayers(firstupper2));%percentage volume activated;;
    area_activated(firstlower2)=area_activated(firstlower2)+(volume_lower/volumelayers(firstlower2));%percentage volume activated;;
else %whole volume in one layer group
    area_activated(firstupper2)=area_activated(firstupper2)+(((4/3) * pi*(radius_activated2^3))/volumelayers(firstupper2));%percentage volume activated;;
     elect2_primaryactivated=((4/3) * pi*(radius_activated2^3))*Neuron_densities(firstupper2);
end

%4. Based on the percentage of primary, what is the
%percentage activated through the secondary
%connections?
%overlappping or non-overlappping activation?
%second

connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
Secondary_neurons_activated=ceil(connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);
tertiary_connections=sum((Secondary_neurons_activated./(Neuron_densities'.*volumelayers'))'.*Numconnections_3groups_array,2);
tertiary_neurons_activated=ceil(tertiary_connections.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
SecondaryTertiary_neurons_activated=Secondary_neurons_activated;%+tertiary_neurons_activated;

%5. weightings based primary and
%secondary or secondary only
%Neuron_densities=[14200, 83800, 164600, 83900,
%177300, 131500]; split into layers
%Neuron_number_groupedlayers=[338+7524 4656 6114+12651].*area_activated;
Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);

neurons_distributed_depth=zeros(ceil(layers_depthdefinition(end)/50),1);
neurons_distributed_depth(round(stimchn1depth/50)-ceil(radius_activated1/50):round(stimchn1depth/50)+ceil(radius_activated1/50))=ceil(elect1_primaryactivated/(ceil(radius_activated1*2/50)+1));
neurons_distributed_depth(round(stimchn2depth/50)-ceil(radius_activated2/50):round(stimchn2depth/50)+ceil(radius_activated2/50))=ceil(elect2_primaryactivated/(ceil(radius_activated2*2/50)+1));

for layertypes=1:3
    midpoint=(middlelayerpoint(layertypes));
    numneurons=SecondaryTertiary_neurons_activated(layertypes);
    if numneurons<2 
        neurons_distributed_depth(ceil(midpoint/50))=neurons_distributed_depth(ceil(midpoint/50))+SecondaryTertiary_neurons_activated(layertypes);
        neuronposition=ceil(midpoint/50);
    elseif rem(numneurons,2)%odd
        distbetweenneurons=(layer_thickness(layertypes)/numneurons);
        neuronposition=(midpoint-(distbetweenneurons)*((numneurons-1)/2)):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons-1)/2);
        neuronposition=ceil(neuronposition/50);
        neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
    elseif rem(numneurons+1,2)%even
        distbetweenneurons=(layer_thickness(layertypes)/numneurons);
        neuronposition=(midpoint-(distbetweenneurons)*((numneurons)/2)+distbetweenneurons/2):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons)/2);
        neuronposition=ceil(neuronposition/50);
        neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
     end

%     if length(unique(neuronposition))~=length(neuronposition)
%         error('sampling resolution too large')
%     end
end



weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);

%bursting activity due to extracellular stim https://www.sciencedirect.com/science/article/pii/S1935861X09000424#app1

burstprob=1./8;%assume cell might fire twice in 8ms window and record for 6ms(2-8)
firingrate=neurons_distributed_depth.*burstprob.*1000;
filtmov=1/3*ones(3,1);
firingrate=filtfilt(filtmov,1,firingrate);
shankposlayers=find(0:50:1240==Depthestimate.s1(end));
simulationFR=firingrate(shankposlayers:end);

figure
plot(0:50:1240,firingrate)
xlabel('Distance from cortical surface (\mum)')
ylabel('Firing rate (Sp/s)')

hold on
plot(0:50:layers_depthdefinition(end),(neurons_distributed_depth.*burstprob.*1000))


%% Estimate electrode depths from model



%1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
minimum_current=1;%uA if electrode is touching the axon
currentdistributions=flipud([1 0; 0.75 0.25; 0.5 0.5; 0.25 0.75; 0 1]);
layers_depthdefinition=[353+149+165 190+353+149+165 525+700+190+353+149+165];
layers_thickness=[368 154 718];
Neuron_densities=[87533.33, 177300, 107700]./10^9;%microcircuit
%layerstart=[0 667 857];
layerstart=[0 369 857];
middlelayerpoint=((layers_thickness)/2)+layerstart;
layer_thickness=(layers_depthdefinition-layerstart);
k_const=8850;%constant ua/mm^2

%model activity
current=6;%uA
stimchn1depth=400;%um
stimchn2depth=700;%um
trial=3;
sepdist=5;
sepcheck=['sep' num2str(sepdist)];
trialcheck='T3';
for current=6%[1 2 3 4 6 8 10] %iterate current
    currentcheck=['C' num2str(current)];
    stimchn1_current=currentdistributions(trial,1)*current;
    stimchn2_current=currentdistributions(trial,2)*current;
    for iterateshank=1:size(saveCplit.(sepcheck).(currentcheck).(trialcheck),2) %iterate shank recording
        data=saveCplit.(sepcheck).(currentcheck).(trialcheck)(:,iterateshank);
        MSE_old=10000000000;
        for stimchniterate=100:50:2100-(sepdist+1)*50
            stimchn1depth=stimchniterate;
            stimchn2depth=stimchniterate+(sepdist+1)*50;
        %note that current dissipates and E-fields never
        %interact; 34um max dist until current dissipates


        if stimchn1_current>=minimum_current
            radius_activated1=((stimchn1_current-minimum_current)/k_const)^0.5;
            radius_activated1=radius_activated1*10.^3;
        else
            radius_activated1=0;
        end
        if stimchn2_current>=minimum_current
            radius_activated2=((stimchn2_current-minimum_current)/k_const)^0.5;
            radius_activated2=radius_activated2*10.^3;
        else
            radius_activated2=0;
        end
        %2. Is radius confined to one group of layers? -
        %use stim elect positions
        %need to define distance based on histology and
        %microdrive depth and LFP

        UpperLstimchn1=stimchn1depth-radius_activated1;
        LowerLstimchn1=stimchn1depth+radius_activated1;
        UpperLstimchn2=stimchn2depth-radius_activated2;
        LowerLstimchn2=stimchn2depth+radius_activated2;



        percentagelayers1u=[1-(layers_depthdefinition-UpperLstimchn1)./(layers_depthdefinition-[0 667 857]) -1];
        percentagelayers1l=[1-(layers_depthdefinition-LowerLstimchn1)./(layers_depthdefinition-[0 667 857]) -1];
        percentagelayers2u=[1-(layers_depthdefinition-UpperLstimchn2)./(layers_depthdefinition-[0 667 857]) -1];
        percentagelayers2l=[1-(layers_depthdefinition-LowerLstimchn2)./(layers_depthdefinition-[0 667 857]) -1];


        firstupper1=find(percentagelayers1u<0,1, 'first')-1;
        firstlower1=find(percentagelayers1l<0,1, 'first')-1;
        firstupper2=find(percentagelayers2u<0,1, 'first')-1;
        firstlower2=find(percentagelayers2l<0,1, 'first')-1;



        %3. what is the percentage of each activated layer in
        %the microcircuit column? i.e. neurons activated
        %only through stimulation - primary
        %0.29mm^3 volume in the circuit 2082um long
        radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
        area_activated=[0 0 0];
        volumelayers= (layers_depthdefinition-[0 667 857]).* (pi.*(radius_circuit).^2);

        %primary connections activated
        if firstupper1~=firstlower1 %cut by layer boundary
            heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
            if heightsegment>radius_activated1 %reverse cap calc
                volume_lower=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
                volume_upper=((4/3) * pi*(radius_activated1^3))-volume_lower;
            else % calc cap
                volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
                volume_lower=((4/3) * pi*(radius_activated1^3))-volume_upper;
            end
            elect1_primaryactivated=volume_upper*Neuron_densities(firstupper1);
            elect1_primaryactivated=elect1_primaryactivated+volume_lower*Neuron_densities(firstlower1);
            area_activated(firstupper1)=volume_upper/volumelayers(firstupper1);%percentage volume activated;
            area_activated(firstlower1)=volume_lower/volumelayers(firstlower1);%percentage volume activated;
        else %whole volume in one layer group
            area_activated(firstupper1)=((4/3) * pi*(radius_activated1^3))/volumelayers(firstupper1);%percentage volume activated
            elect1_primaryactivated=((4/3) * pi*(radius_activated1^3))*Neuron_densities(firstupper1);
        end

        %primary connections activated
        if firstupper2~=firstlower2 %cut by layer boundary
            heightsegment=layers_depthdefinition(firstupper2)-UpperLstimchn2;
            if heightsegment>radius_activated2 %reverse cap calc
                volume_lower=((pi*heightsegment^2)/3)*(3*radius_activated2-heightsegment);
                volume_upper=((4/3) * pi*(radius_activated2^3))-volume_lower;

            else % calc cap
                volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated2-heightsegment);
                volume_lower=((4/3) * pi*(radius_activated2^3))-volume_upper;
            end
            elect2_primaryactivated=volume_upper*Neuron_densities(firstupper1);
            elect2_primaryactivated=elect2_primaryactivated+volume_lower*Neuron_densities(firstlower2);
            area_activated(firstupper2)=area_activated(firstupper2)+(volume_upper/volumelayers(firstupper2));%percentage volume activated;;
            area_activated(firstlower2)=area_activated(firstlower2)+(volume_lower/volumelayers(firstlower2));%percentage volume activated;;
        else %whole volume in one layer group
            area_activated(firstupper2)=area_activated(firstupper2)+(((4/3) * pi*(radius_activated2^3))/volumelayers(firstupper2));%percentage volume activated;;
            elect2_primaryactivated=((4/3) * pi*(radius_activated2^3))*Neuron_densities(firstupper2);
        end

        %4. Based on the percentage of primary, what is the
        %percentage activated through the secondary
        %connections?
        %overlappping or non-overlappping activation?
        %second

        connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
        Secondary_neurons_activated=ceil(connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
        baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);

        %5. weightings based primary and
        %secondary or secondary only
        %Neuron_densities=[14200, 83800, 164600, 83900,
        %177300, 131500]; split into layers
        %Neuron_number_groupedlayers=[338+7524 4656 6114+12651].*area_activated;
        Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);
        neurons_distributed_depth=zeros(ceil(layers_depthdefinition(end)/50),1);
        neurons_distributed_depth(round(stimchn1depth/50)-ceil(radius_activated1/50):round(stimchn1depth/50)+ceil(radius_activated1/50))=ceil(elect1_primaryactivated/(ceil(radius_activated1*2/50)+1));
        neurons_distributed_depth(round(stimchn2depth/50)-ceil(radius_activated2/50):round(stimchn2depth/50)+ceil(radius_activated2/50))=ceil(elect2_primaryactivated/(ceil(radius_activated2*2/50)+1));

        for layertypes=1:3
            midpoint=(middlelayerpoint(layertypes));
            numneurons=Secondary_neurons_activated(layertypes);
            if numneurons<2
                neurons_distributed_depth(ceil(midpoint/50))=neurons_distributed_depth(ceil(midpoint/50))+Secondary_neurons_activated(layertypes);
                neuronposition=ceil(midpoint/50);
            elseif rem(numneurons,2)%odd
                distbetweenneurons=(layer_thickness(layertypes)/numneurons);
                neuronposition=(midpoint-(distbetweenneurons)*((numneurons-1)/2)):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons-1)/2);
                neuronposition=ceil(neuronposition/50);
                neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
            elseif rem(numneurons+1,2)%even
                distbetweenneurons=(layer_thickness(layertypes)/numneurons);
                neuronposition=(midpoint-(distbetweenneurons)*((numneurons)/2)+distbetweenneurons/2):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons)/2);
                neuronposition=ceil(neuronposition/50);
                neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
            end

            if length(unique(neuronposition))~=length(neuronposition)
                error('sampling resolution too large')
            end
        end



        weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);

        %bursting activity due to extracellular stim https://www.sciencedirect.com/science/article/pii/S1935861X09000424#app1

        burstprob=2./8;%assume cell might fire twice in 8ms window and record for 6ms(2-8)
        firingrate=neurons_distributed_depth.*burstprob.*1000;
        filtmov=1/5*ones(5,1);
        firingrate=filtfilt(filtmov,1,firingrate);
        
        for lengthCortex=1:42-15
            if lengthCortex*50>stimchn1depth || stimchn2depth>lengthCortex*50+15*50 %ensure electrodes are in the estimated recording region    
                continue
            end
            firingrate_window=firingrate(lengthCortex:lengthCortex+15);
            datanonan=flipud(data(~isnan(data)));%remove nan and flip so shallow is at top
            MSE=(sum((datanonan-firingrate_window).^2,'all','omitnan'))/sum(~isnan(datanonan),'all','omitnan');
            if MSE<MSE_old
                MSE_old=MSE;
                electrodepos.(sepcheck).(currentcheck).(trialcheck)(iterateshank)=stimchn1depth;
                depthrecordingfirst.(sepcheck).(currentcheck).(trialcheck)(iterateshank)=lengthCortex; %begin sample point for recording electrodes  
                FR=firingrate_window;
                FR_raw=neurons_distributed_depth(lengthCortex:lengthCortex+15).*burstprob.*1000;
            end
        end

        end
    end
end


figure
plot(2050:-50:0,flipud(firingrate))
xlabel('Distance from cortical surface (\mum)')
ylabel('Firing rate (Sp/s)')

hold on
plot(0:50:layers_depthdefinition(end),flipud(neurons_distributed_depth.*burstprob.*1000))

