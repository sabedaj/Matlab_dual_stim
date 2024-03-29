%% flash raster
%D:\DATA\Rat_023\flash_001_211011_113004
chn=25;
count_sp_tr=Flash_raster(chn);
%% stim raster
%D:\DATA\Rat_023\S3E5_9elect_001_211011_123736
chn=25;
amp=6; %check trial info - ID 26:30 is 6uA
trig=loadTrig(0);
TP = loadTrialParams;
c = ismember(cell2mat(TP(:,2)), [26:30]);%26:30ID 
indexes = find(c);
tID=indexes(2:2:end)./2;
trig(trig==-500)=[];
trig=trig(tID);
theseTrig = trig./30;
nT=length(trig);
%run Flash_raster(chn) - stop line 40 and load in these trig
%% single pairs R19 E:\DATA\Rat_019\S3E2_9elect_001_210511_104603 and R23 for reviewer
SinglePairWErrorBars
%% save array layer classification
ratN='Rat_020';
cd(['E:\DATA\' ratN])

superL=[];% electrodes in superficial layers
inputL=[16 32 48 64];%electrodes in input layers
deepL=[1:15 (1:15)+32 (1:15)+48 (1:15)+16];% electrodes in deep layers
WM=[];

superL=[10:16 (10:16)+32 (11:16)+48 (11:16)+16];% electrodes in superficial layers
inputL=[7:9 (7:9)+32 (8:10)+48 (8:10)+16];%electrodes in input layers
deepL=[1:6 (1+32):(6+32) (1+48):(7+48) (1+16):(7+16)];% electrodes in deep layers
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
save(['E:\DATA\' ratN filesep 'ElectLayerClass.mat'],'ElectLayerClass')
%%
%r23
% superL=[14:16 (13+32):(16+32) (14+48):(16+48) (15:16)+16];% electrodes in superficial layers
% inputL=[11:13 (10+32):(12+32) (11+48):(13+48) (12+16):(14+16)];%electrodes in input layers
% deepL=[1:10 (1+32):(9+32) (1+48):(10+48) (1+16):(11+16)];% electrodes in deep layers
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
%hist and CSD
% superL=[12:16 (12+32):(16+32) (14+48):(16+48)];% electrodes in superficial layers
% inputL=[8:11 (8+32):(11+32) (10+48):(13+48) (13:16)+16];%electrodes in input layers
% deepL=[1:7 (1+32):(7+32) (1+48):(9+48) (1+16):(12+16)];% electrodes in deep layers
% WM=[];




%R21
% superL=[];% electrodes in superficial layers
% inputL=[];%electrodes in input layers
% deepL=[2:16 (2+32):(16+32) (3+48):(16+48) (3+16):(16+16)];% electrodes in deep layers
% WM=[1 (1+32) (1+48):(2+48) (1+16):(2+16)];

%R13 - bad flash resp ignore
% superL=[];% electrodes in superficial layers
% inputL=[16 31:32 47:48 63:64];%electrodes in input layers
% deepL=[1:15 (1:14)+32 (1:14)+48 (1:14)+16];% electrodes in deep layers
% WM=[];

%R12
%LFP
% superL=[10:16 (10:16)+32 (11:16)+48 (11:16)+16];% electrodes in superficial layers
% inputL=[7:9 (7:9)+32 (8:10)+48 (8:10)+16];%electrodes in input layers
% deepL=[1:6 (1+32):(6+32) (1+48):(7+48) (1+16):(7+16)];% electrodes in deep layers
% WM=[];


%R9
% superL=[];% electrodes in superficial layers
% inputL=[];%electrodes in input layers
% deepL=[7:16 (7:16)+32 (6:16)+48 (6:16)+16];% electrodes in deep layers
% WM=[1:6 (1:6)+32 (1:5)+48 (1:5)+16];


%R6
% superL=[];% electrodes in superficial layers
% inputL=[13:16 (13:16)+32 (13:16)+48 (13:16)+16];%electrodes in input layers
% deepL=[1:12 (1+32):(12+32) (1+48):(12+48) (1+16):(12+16)];% electrodes in deep layers
% WM=[];


%% sort laminar current steering data
skip_stimshank=0;

startpointseconds=2;
secondstoanalyse=8;
AMP=[0 1 2 3 4 6 8 10];
layer_types=['ss';'gg';'ii';'si';'sg';'gi';'WM'];
%setup data structures
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
         psingletrials.(sepcheck).(currcheck)=[];
        for trial=1:5
            trialcheck=['T' num2str(trial)];
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                Depthestimate.(shanksepcheck)=zeros(16,1);
                for layername=1:length(layer_types)
                    saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
                    savecsplit_avgpairs.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)=[];
                    savelayeraligned.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
                end
                for i=1:32
                    pos=['P' num2str(i)];
                    saveStimalignedLayeraligned.(pos).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
                end
            end
        end
    end
     Npairs.(sepcheck)=0;
end
savalllayers=[];
for ratN=14:20%[6 9 12 13 21:23]%loop through animals 14:20%[6 9 12 13 21:23] %23%
    %load data
    if ratN<10
        Ratnum=['Rat_00' num2str(ratN)];
    elseif ratN>=10
        Ratnum=['Rat_0' num2str(ratN)];
    end

    cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
    D_data=dir;
    if ~any(strcmp({D_data.name}, 'ElectLayerClass.mat')) %check that we were able to correctly identify layers
        continue
    end
    load('ElectLayerClass.mat','ElectLayerClass')
    if length(ElectLayerClass)~=64
             error('Layers not saved properly')
    end
    for k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
        currD = D_data(k).name; % Get the current subdirectory name
        try
            cd([D_data(k).folder filesep currD])
            load('Averagetrialresponse.mat','avgnospT') % load sorted neural activity
        catch
            continue;
        end

        loadStimChn;
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
        AMP=loadAMP;
        if stimChn(1)<17
        stimChn_NS=stimChn;
        elseif stimChn(1)<33
            stimChn_NS=stimChn-16;
        elseif stimChn(1)<49
            stimChn_NS=stimChn-32;
        elseif stimChn(1)<65
            stimChn_NS=stimChn-48;
        end
        for sepdist=5:2:9
            if length(stimChn)<2||(stimChn(1)-stimChn(2))~=(-1*sepdist)-1 % check stimchn sepdist between these electrodes
                continue
            end
%             for i=1:64
%                 e1=avgnospT(i,1:5:(length(AMP)*5)).*1000/6;
%                 e2=avgnospT(i,5:5:length(AMP)*5).*1000/6;
%                 if e1(end-1)<10 || e1(end)<10 || e2(end-1)<10 || e2(end)<10
%                     avgnospT(i,:)=nan;
%                 end
%             end
            [Csplit_shankdist, Csplit_depthsep]=PoolNormalisedActivity_refactored(avgnospT,startpointseconds, secondstoanalyse, ElectLayerClass);
            sepcheck=['sep' num2str(sepdist)];
            AMP(AMP==-1)=0;
            for current=1:length(AMP)
                currcheck=['C' num2str(AMP(current))];
                stimpos_layer=find(~isnan(Csplit_depthsep.(currcheck).T1.D0),stimChn_NS(1),'first');
                pos=['P' num2str(stimpos_layer(stimChn_NS(1)))];
%                 if AMP(current)~=0
%                     psingle=SinglePairWErrorBars(AMP(current),0);
%                     psingletrials.(sepcheck).(currcheck)=[psingletrials.(sepcheck).(currcheck) psingle];
%                 end
                
                for trial=1:5
                    trialcheck=['T' num2str(trial)];
                    for shanksep=0:3
                        shanksepcheck=['D' num2str(shanksep)];
                        saveCplit.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[saveCplit.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck) Csplit_shankdist.(currcheck).(trialcheck).(shanksepcheck)];
                        savelayeraligned.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[savelayeraligned.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck) Csplit_depthsep.(currcheck).(trialcheck).(shanksepcheck)];
                        
                        saveStimalignedLayeraligned.(pos).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[saveStimalignedLayeraligned.(pos).(sepcheck).(currcheck).(trialcheck).(shanksepcheck) Csplit_depthsep.(currcheck).(trialcheck).(shanksepcheck)];
                        if AMP(current)==6 && sepdist==5
                            savalllayers=[savalllayers Csplit_depthsep.(currcheck).(trialcheck).(shanksepcheck)];
                        end
                        
                    end
                    
                    if skip_stimshank==1
                        tmp=nanmean([Csplit_shankdist.(currcheck).(trialcheck).D1 Csplit_shankdist.(currcheck).(trialcheck).D2 Csplit_shankdist.(currcheck).(trialcheck).D3],2);
                    else
                        tmp=nanmean([Csplit_shankdist.(currcheck).(trialcheck).D0 Csplit_shankdist.(currcheck).(trialcheck).D1 Csplit_shankdist.(currcheck).(trialcheck).D2 Csplit_shankdist.(currcheck).(trialcheck).D3],2);
                    end
                    savecsplit_avgpairs.(layers_stim).(sepcheck).(currcheck).(trialcheck)=[savecsplit_avgpairs.(layers_stim).(sepcheck).(currcheck).(trialcheck) tmp]; %already averaged pairs
                end
            end
            Npairs.(sepcheck)=Npairs.(sepcheck)+1;
        end
    end
end
saveCplit_raw=saveCplit;
%% simulation
AMP=[0 1 2 3 4 6 8 10];
layer_types=['ss';'gg';'ii';'si';'sg';'gi';'WM'];
%setup data structures
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for trial=1:5
            trialcheck=['T' num2str(trial)];
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                for layername=1:length(layer_types)
                    Csplit_simulation.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
                end
            end
        end
    end
end

% reading connectome data - https://bbp.epfl.ch/nmc-portal/downloads.html


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


layers_depthdefinition=[368 522 1240];
Depthestimate.s1=zeros(16,1);
Depthestimate.s2=zeros(16,1);
Depthestimate.s3=zeros(16,1);
Depthestimate.s4=zeros(16,1);
shankorder=[1 4 2 3];
currentdistributions=flipud([1 0; 0.75 0.25; 0.5 0.5; 0.25 0.75; 0 1]);

for ratN=[6 9 13 21:23]%[6:13 21:23] %loop through animals
    %load data
    if ratN<10
        Ratnum=['Rat_00' num2str(ratN)];
    elseif ratN>=10
        Ratnum=['Rat_0' num2str(ratN)];
    end

    cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
    D_data=dir;
    if ~any(strcmp({D_data.name}, 'ElectLayerClass.mat')) %check that we were able to correctly identify layers
        continue
    end
    load('ElectLayerClass.mat','ElectLayerClass')
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
    for k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
        currD = D_data(k).name; % Get the current subdirectory name
        try
            cd([D_data(k).folder filesep currD])
            load('Averagetrialresponse.mat','avgnospT') % load sorted neural activity
        catch
            continue;
        end
        loadStimChn;
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
        for sepdist=5:2:9
            sepcheck=['sep' num2str(sepdist)];
            if (stimChn(1)-stimChn(2))~=(-1*sepdist)-1 % check stimchn sepdist between these electrodes
                continue
            end
            % spike rate and centroid calculation
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
            for current = 1:length(AMP)
                currcheck=['C' num2str(AMP(current))];

                for trial=1:5
                    trialcheck=['T' num2str(trial)];
                    for shanksep=0:3
                        shanksepcheck=['D' num2str(shanksep)];
                        %%%%%%%%%calculate theoretical weightings

                        %1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
                        minimum_current=1;%uA if electrode is touching the axon
                        %note that current dissipates and E-fields never
                        %interact; 34um max dist until current dissipates
                        stimchn1_current=currentdistributions(trial,1)*AMP(current);
                        stimchn2_current=currentdistributions(trial,2)*AMP(current);

                        %model activity
                        stimshankcheck=['s' num2str(shank)];
                        stimchn1depth=Depthestimate.(stimshankcheck)(stimChn(1));%um
                        stimchn2depth=Depthestimate.(stimshankcheck)(stimChn(2));%um


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
                        firstupper1=find(percentagelayers1u<1,1, 'first');
                        firstlower1=find(percentagelayers1l<1,1, 'first');
                        firstupper2=find(percentagelayers2u<1,1, 'first');
                        firstlower2=find(percentagelayers2l<1,1, 'first');

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
                        if firstupper1~=4 || firstlower1~=4
                            if firstupper1~=firstlower1 && firstlower1==4
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
                                elect1_primaryactivated=sum(volume_upper_percentage)*Neuron_densities(firstupper1);
                                area_activated(firstupper1)=area_activated(firstupper1)+(sum(volume_upper_percentage)/volumelayers(firstupper1));%percentage volume activated;;
                            elseif firstupper1~=firstlower1 && firstlower1~=4%cut by layer boundary
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
                        else
                            elect1_primaryactivated=0;
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
                        SecondaryTertiary_neurons_activated=Secondary_neurons_activated+tertiary_neurons_activated;

                        %5. weightings based primary and
                        %secondary or secondary only
                        %Neuron_densities=[14200, 83800, 164600, 83900,
                        %177300, 131500]; split into layers
                        %Neuron_number_groupedlayers=[338+7524 4656 6114+12651].*area_activated;
                        Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);

                        neurons_distributed_depth=zeros(ceil(layers_depthdefinition(end)/50),1);
                        %                                                     neurons_distributed_depth(round(stimchn1depth/50)-ceil(radius_activated1/50):round(stimchn1depth/50)+ceil(radius_activated1/50))=ceil(elect1_primaryactivated/(ceil(radius_activated1*2/50)+1));
                        %                                                     neurons_distributed_depth(round(stimchn2depth/50)-ceil(radius_activated2/50):round(stimchn2depth/50)+ceil(radius_activated2/50))=ceil(elect2_primaryactivated/(ceil(radius_activated2*2/50)+1));

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
                        Csplit_simulation.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[Csplit_simulation.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck) [NaN(16-stimChn(1),1); simulationFR; NaN(32-(16-stimChn(1))-16,1)]];
                    end
                end
            end
        end


    end
end


%% normalise all shanks
clear saveCplit_norm
%saveCplit=savelayeraligned;
saveCplit=saveCplit_raw;
%saveCplit=Csplit_simulation;
AMP=[0 1 2 3 4 6 8 10];
FilterLength=3;
totalrespelect=0;
totalelect=0;
numstimelect=0;
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)
            for shanksep=0:3
                if shanksep==0 && skip_stimshank==1
                    continue
                end
                shanksepcheck=['D' num2str(shanksep)];
                Data_allcols=cell(5,1);
                max_all=cell(5,1);
                for trial=1:5
                    tname=['T' num2str(trial)];
                    Data_allcols{trial}=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname).(shanksepcheck);
                    Data_allcols{trial}=smoothData(Data_allcols{trial},FilterLength);
                    Data_allcols{trial}(Data_allcols{trial}<0)=0;
                    max_all{trial}=max(Data_allcols{trial});
                    % How many electrodes responded to stimulation (FR>basline)
                    if (trial==1 || trial==5) && ~isempty(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname).(shanksepcheck)) && strcmp(currcheck,'C6')
                        totalrespelect=totalrespelect+sum((saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname).(shanksepcheck)>0),'all');
                        totalelect=totalelect+(size(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname).(shanksepcheck),2)*16);
                    end
                end
                
%                 max_all=zeros(32,size(Data_allcols{1},2));
%                 for i=1:size(Data_allcols{1},2)
%                     max_all(:,i)=max([Data_allcols{1}(:,i), Data_allcols{2}(:,i), Data_allcols{3}(:,i), Data_allcols{4}(:,i), Data_allcols{5}(:,i)],[],2);
%                 end
                if ~isempty(max_all)
                for trial=1:5
                    tname=['T' num2str(trial)];
                    saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname).(shanksepcheck)=Data_allcols{trial}./max_all{trial};
                    saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname).(shanksepcheck)(:,max_all{trial}==0)=Data_allcols{trial}(:,max_all{trial}==0);
                end
                end
            end
        end
    end
end

%% normalise pooled pairs %%not this
% saveCplit=savecsplit_avgpairs;
% AMP=[0 1 2 3 4 6 8 10];
% FilterLength=3;
% for sepdist=5:2:9
%     sepcheck=['sep' num2str(sepdist)];
%     for current=1:length(AMP)
%         currcheck=['C' num2str(AMP(current))];
%         for layername=1:length(layer_types)
%             Data_allcols=cell(5,1);
%             max_all=cell(5,1);
%             for trial=1:5
%                 tname=['T' num2str(trial)];
%                 Data_allcols{trial}=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname);
%                 Data_allcols{trial}=smoothData(Data_allcols{trial},FilterLength);
%                 Data_allcols{trial}(Data_allcols{trial}<0)=0;
%                 max_all{trial}=max(Data_allcols{trial});
%             end
% 
%             %                 max_all=zeros(32,size(Data_allcols{1},2));
%             %                 for i=1:size(Data_allcols{1},2)
%             %                     max_all(:,i)=max([Data_allcols{1}(:,i), Data_allcols{2}(:,i), Data_allcols{3}(:,i), Data_allcols{4}(:,i), Data_allcols{5}(:,i)],[],2);
%             %                 end
%             if ~isempty(max_all)
%                 for trial=1:5
%                     tname=['T' num2str(trial)];
%                     saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname)=Data_allcols{trial}./max_all{trial};
%                     saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(tname)(:,max_all{trial}==0)=Data_allcols{trial}(:,max_all{trial}==0);
%                 end
%             end
%         end
%     end
% end

% %% normalise
% AMP=[0 1 2 3 4 6 8 10];
% for sepdist=5:2:9
%     sepcheck=['sep' num2str(sepdist)];
%     dat2plot=[];
%     for current=1:length(AMP)
%         currcheck=['C' num2str(AMP(current))];
%         for layername=1:length(layer_types)
%             for shanksep=0:3
%                 if shanksep==0 && skip_stimshank==1
%                     continue
%                 end
%                 shanksepcheck=['D' num2str(shanksep)];
%                 Data_allcols1=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T1').(shanksepcheck);
%                 Data_allcols2=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T2').(shanksepcheck);
%                 Data_allcols3=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T3').(shanksepcheck);
%                 Data_allcols4=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T4').(shanksepcheck);
%                 Data_allcols5=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T5').(shanksepcheck);
%                 max_all=max([Data_allcols1, Data_allcols2, Data_allcols3, Data_allcols4, Data_allcols5],[],2);
%                 saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T1').(shanksepcheck)=Data_allcols1./max(Data_allcols1);
%                 saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T2').(shanksepcheck)=Data_allcols2./max(Data_allcols2);%max_all;
%                 saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T3').(shanksepcheck)=Data_allcols3./max(Data_allcols3);
%                 saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T4').(shanksepcheck)=Data_allcols4./max(Data_allcols4);
%                 saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T5').(shanksepcheck)=Data_allcols5./max(Data_allcols5);
%             end
%         end
%     end
% end

%% plotting


%simulation or normalise? pick 1 or none to plot normal data
plotlayeraligned=0;
plotsim=0;
normalise_dat=0;

%YOU CAN ONLY PICK ONE OF THESE
splitshanks=0; %Split into stim, 1 shank away and 2 shanks away etc
splitlayers=0; %split into layers
%AND IF YOU PICK ONE^^^, YOU CAN'T PICK ALL CURRENTS
singleCurrent=6; %plot 6uA results
FilterLength=1;%smoothing parameters %%should be 1 if you have normalised the data above.
%plot laminar
Lam_Acc='across';%options 'across' or 'laminar'
%what do you want to class as indipendent samples
avgpairsamplesplot=0;
if normalise_dat==1
    saveCplit=saveCplit_norm;
elseif plotsim==1
    saveCplit=Csplit_simulation;
elseif plotlayeraligned==1
    saveCplit=savelayeraligned;
else
    saveCplit=saveCplit_raw;
end


if splitshanks==1 || splitlayers==1
    downsampleYN=0; %do you want to downsample 1=yes %%%currently not implemented
    numshanksToAvg=1;
elseif skip_stimshank==1
    numshanksToAvg=3;
    downsampleYN=1; %do you want to downsample 1=yes
else
    numshanksToAvg=6;
    downsampleYN=1; %do you want to downsample 1=yes
end
downsampleYN=1;
if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2
shanknames={'Stimshank'; '1 shank away'; '2 shanks away'; '3 shanks away'};

for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    figure
    if ~splitshanks && ~splitlayers
        %axes('Position',[0.13         0.112396822842341                     0.775         0.62])
        if strcmp(Lam_Acc,'laminar')
            hold on
            ax1=gca;
            hold(ax1, 'on');
            axes('Position',[ax1.Position(1)+ax1.Position(3)*3/4 ax1.Position(2) ax1.Position(3)/4 ax1.Position(4)])
            ax2=gca;
            hold(ax2,'off')
            ax1.Position=[ax1.Position(1) ax1.Position(2) ax1.Position(3)*3/4 ax1.Position(4)];
        else
            hold on
            ax2=gca;
            hold(ax2, 'on');
            axes('Position',[ax2.Position(1) ax2.Position(2) ax2.Position(3) ax2.Position(4)*3/4])
            ax1=gca;
            ax2.Position=[ax2.Position(1) ax2.Position(2)+ax2.Position(4)*3/4 ax2.Position(3) ax2.Position(4)/4];
            hold(ax1,'off')
        end
    end
    clear peak_all
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        dat2plot=[];
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            for layername=1:length(layer_types)
                for shanksep=0:3
                    if shanksep==0 && skip_stimshank==1
                        continue
                    end
                    shanksepcheck=['D' num2str(shanksep)];
                    if splitlayers==1 && avgpairsamplesplot==0
                        dat2plot=[dat2plot saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                        continue
                    elseif splitlayers==1 && avgpairsamplesplot==1 && ((shanksep==0 && skip_stimshank==0) || (shanksep==1 && skip_stimshank==1))
                        if normalise_dat==1 %% make sure you normalised this above^^
                            dat2plot=[dat2plot saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)];
                        else
                            dat2plot=[dat2plot savecsplit_avgpairs.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)];
                        end
                        continue
                    elseif splitshanks==1 && layername==1
                        dat2plot.(shanksepcheck)=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck);
                    elseif splitshanks==1 && layername==length(layer_types)
                        %plot
                        if trial==1
                            %subplot(1,4,shanksep+1)
                            figure
                            ax1.(shanksepcheck)=gca;
                            hold(ax1.(shanksepcheck), 'on');
                            axes('Position',[ax1.(shanksepcheck).Position(1) ax1.(shanksepcheck).Position(2)+ax1.(shanksepcheck).Position(4)*3/4 ax1.(shanksepcheck).Position(3) ax1.(shanksepcheck).Position(4)/4])
                            ax2.(shanksepcheck)=gca;
                            title(shanknames{shanksep+1})
                            hold(ax2.(shanksepcheck),'off')
                            ax1.(shanksepcheck).Position=[ax1.(shanksepcheck).Position(1) ax1.(shanksepcheck).Position(2) ax1.(shanksepcheck).Position(3) ax1.(shanksepcheck).Position(4)*3/4];
                        end
                        dat2plot.(shanksepcheck)=[dat2plot.(shanksepcheck) saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                        dat2plot.(shanksepcheck)((sum(~isnan(dat2plot.(shanksepcheck)),2)<size(dat2plot.(shanksepcheck),2)/2),:)=nan;
                        smoothedenvelope=smoothData(dat2plot.(shanksepcheck),FilterLength);
                        [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1.(shanksepcheck),ax2.(shanksepcheck),Lam_Acc);
                        peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=peak;
                        if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck))
                            [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T1.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T2.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T3.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T4.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T5.(shanksepcheck)(:)],1,'off');
                            %text(ax2.(shanksepcheck),find(~isnan(lineOut.YData),1,'first')+1,DATApos1+9.5,1,['p=' num2str(p)])
                            N=size(peak,1)
                            p
                        end
                    elseif splitshanks==1 && layername>1 && layername<length(layer_types)
                        dat2plot.(shanksepcheck)=[dat2plot.(shanksepcheck) saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                    elseif avgpairsamplesplot==0
                        dat2plot=[dat2plot saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                        continue
                    elseif avgpairsamplesplot==1 && ((shanksep==0 && skip_stimshank==0) || (shanksep==1 && skip_stimshank==1))
                        if normalise_dat==1 %% make sure you normalised this above^^
                            dat2plot=[dat2plot saveCplit_norm.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)];
                        else
                            dat2plot=[dat2plot savecsplit_avgpairs.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)];
                        end

                    end
                end
                if splitlayers==0
                    continue
                end
                layernameFULL=[layer_types(layername),layer_types(layername+length(layer_types))];
                if trial==1
                    figure%subplot(2,4,layername)
                    hold on
                    if strcmp(Lam_Acc,'laminar')
                        hold on
                        ax1.(layernameFULL)=gca;
                        hold(ax1.(layernameFULL), 'on');
                        axes('Position',[ax1.(layernameFULL).Position(1)+ax1.(layernameFULL).Position(3)*3/4 ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3)/4 ax1.(layernameFULL).Position(4)])
                        ax2.(layernameFULL)=gca;
                        title(layernameFULL)
                        hold(ax2.(layernameFULL),'off')
                        ax1.Position=[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3)*3/4 ax1.(layernameFULL).Position(4)];
                    else
                        hold on
                        ax2.(layernameFULL)=gca;
                        hold(ax2.(layernameFULL), 'on');
                        axes('Position',[ax2.(layernameFULL).Position(1) ax2.(layernameFULL).Position(2) ax2.(layernameFULL).Position(3) ax2.(layernameFULL).Position(4)*3/4])
                        ax1.(layernameFULL)=gca;
                        title(layernameFULL)
                        hold(ax1.(layernameFULL),'off')
                        ax2.(layernameFULL).Position=[ax2.(layernameFULL).Position(1) ax2.(layernameFULL).Position(2)+ax2.(layernameFULL).Position(4)*3/4 ax2.(layernameFULL).Position(3) ax2.(layernameFULL).Position(4)/4];
                    end

                end
                if isempty(dat2plot)
                    continue
                end
                dat2plot((sum(~isnan(dat2plot),2)<size(dat2plot,2)/2),:)=nan;
                smoothedenvelope=smoothData(dat2plot,FilterLength);
                [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1.(layernameFULL),ax2.(layernameFULL),Lam_Acc);
                peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=peak;
                if trial==1
                    senv.([layer_types(layername),layer_types(layername+length(layer_types))])=smoothedenvelope;
                end
                if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1)
                    [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T2(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T3(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T4(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T5(:)],1,'off');
                    text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),16,1,['p=' num2str(p)])
                    N=size(peak,1)
                    text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),16,2,['N=' num2str(N)])
                    p
                     psingle=signrank(nanmean(senv.([layer_types(layername),layer_types(layername+length(layer_types))]),2),nanmean(smoothedenvelope,2));
                end
                dat2plot=[];
            end
            if splitlayers==1 || splitshanks==1
                continue
            elseif singleCurrent==0
                subplot(2,4,current)
            end
            hold on
            dat2plot((sum(~isnan(dat2plot),2)<size(dat2plot,2)/2),:)=nan;
            if downsampleYN==1%%%%%%%%%%%%%%%%%%%
                if trial==1
                    lengthneeded=60;%60;
                    [samples_rand1]=DownSample(dat2plot,lengthneeded,s,seedpoint);
                    samples_rand.(sepcheck)=samples_rand1;
                end
                dat2plot=dat2plot(:, samples_rand.(sepcheck));% downsample%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            end
            smoothedenvelope=smoothData(dat2plot,FilterLength);
            [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1,ax2,Lam_Acc);
            if trial==3
                curve=nanmean(smoothedenvelope,2);
                maxsepdist(sepdist)=max(curve);
                FWHMsepdist(sepdist)=sum((curve>=max(curve)/2))./sum(~isnan(curve));
            end
            if trial==1
                senv=smoothedenvelope;
            end
            if trial==3
                senv_mid=smoothedenvelope;
            end
            peak_all.(currcheck).(trialcheck)=peak;
            Peakpos_lam.(sepcheck).(currcheck).(trialcheck)=peak;%Peakpos_lam
            if trial==5 && exist('peak_all') && ~isempty(peak_all.(currcheck).(trialcheck))
                peaks=[peak_all.(currcheck).T1(:),peak_all.(currcheck).T2(:),peak_all.(currcheck).T3(:),peak_all.(currcheck).T4(:),peak_all.(currcheck).T5(:)];
                pgroup={'a','a','a','a','a'};
                %[p,tbl,stats]=friedman(peaks,1,'off');
                p=sigtestrcontinousvar(peaks(:,2:4),1000,'left');
                %[p,tbl]=anova1(peaks,pgroup,'off');
                text(ax2,DATApos1+9.5,1,['p=' num2str(p)])
                N=size(peaks,1)
                p
                %test whether single electrode conditions are sig diff
                psingle=signrank(nanmean(senv,2),nanmean(smoothedenvelope,2));
                %testing sequential
                sequentialpeaks=diff(peaks');
                fprintf(['Percentage sequentially shifted: ' num2str(sum(sum(sequentialpeaks,1)<0 & all(sequentialpeaks<1))*100/length(sequentialpeaks)) '\n'])
                
            end
            dat2plot=[];
        end

    end
end
%% test significance of peak shift with nic's method


%% single electrode stim - separation based on layer
distarray1=[(15:-1:1) 0 1:16]*50;
distgroup=cell(3,1);
distgroup(1:3)={cell(41,1)};
individuallayers=['s';'g';'i'];
currcheck='C6';
layer_types=['ss';'gg';'ii';'si';'sg';'gi';'WM'];
trials=[1 5];
singleElect=cell(6,1);
for layerIT=1:length(layer_types)-1
    layercheck=layer_types(layerIT,:);
    for sepdist=5:2:9
        sepcheck=['sep' num2str(sepdist)];
        for trialIT=1:length(trials)
            trialcheck=['T' num2str(trials(trialIT))];
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                if ~isempty(saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck))
                    layer=find(layercheck(trialIT)==individuallayers);
                    if trials(trialIT)==5
                        singleElect{layer}=[singleElect{layer} [saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck); nan(4,size(saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck),2))]];
                        
                        
                        for row=1:32
                            distanceofgroup=sqrt(distarray1(row)^2 + (200*shanksep)^2);
                            [~,group]=min(abs((0:50:2000)-distanceofgroup));
                            distgroup{layer}{group}=[distgroup{layer}{group} saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)(row,:)];
                        end
                    else
                        differenceinelectdist=9-sepdist;
                        nantoadd=nan(differenceinelectdist,size(saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck),2));
                        nanbottom=nan(4-differenceinelectdist,size(saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck),2));
                        singleElect{layer+3}=[singleElect{layer+3} [nantoadd; saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck); nanbottom]];
                        
                        
                        distarray2=[(16+sepdist:-1:1) 0 1:16-sepdist]*50;
                        for row=1:32
                            distanceofgroup=sqrt(distarray2(row)^2 + (200*shanksep)^2);
                            [~,group]=min(abs((0:50:2000)-distanceofgroup));
                            distgroup{layer}{group}=[distgroup{layer}{group} saveCplit_raw.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)(row,:)];
                        end
                        
                    end
                end
                
            end
        end
    end
end



%%
figure; hold on
colours=[1,0,0;0,1,0;0,0,1];
for layers=1:3
    difftoadd=nan(10,size(singleElect{layers},2));
    catarray1=[difftoadd; singleElect{layers}];
    difftoadd2=nan(size(catarray1,1)-size(singleElect{layers+3},1),size(singleElect{layers+3},2));
    catarray2=[[singleElect{layers+3};difftoadd2], catarray1];
    stdshade(catarray2',0.2,colours(layers,:))
end
legend('s','g','i')
figure; hold on
colours=[1,0,0;0,1,0;0,0,1;1,1,0;0,1,1;1,0,1];
for layers=2:6
    
    stdshade(singleElect{layers}',0.2,colours(layers,:))
end
legend('g','i','s','g','i')
%% distance
figure; hold on
for layers=1:3
    plot(cellfun(@nanmean, distgroup{layers}))
end
%% accross vs laminar peak plots
load('Peakpos_acc.mat','Peakpos_acc')
load('Peakpos_lam.mat','Peakpos_lam')
color_current=cmap([1:5]*floor((length(cmap))/5),:);
mean_acc=zeros(5,1);
mean_lam=zeros(5,1);
std_acc=zeros(5,1);
std_lam=zeros(5,1);
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for trial=1:5
            trialcheck=['T' num2str(trial)];
            mean_acc(trial,1)=(mean(Peakpos_acc.(sepcheck).(currcheck).(trialcheck))-16)*50;
            mean_lam(trial,1)=(mean(Peakpos_lam.(sepcheck).(currcheck).(trialcheck))-16)*50;
            std_acc(trial,1)=(std(Peakpos_acc.(sepcheck).(currcheck).(trialcheck))-16)*50/sqrt(length(Peakpos_acc.(sepcheck).(currcheck).(trialcheck)));
            std_lam(trial,1)=(std(Peakpos_lam.(sepcheck).(currcheck).(trialcheck))-16)*50/sqrt(length(Peakpos_lam.(sepcheck).(currcheck).(trialcheck)));
        end
    end
    figure
    hold on
%     errorbar(mean_acc,mean_lam,std_lam,std_lam,std_acc,std_acc,'k.')
%     scatter(mean_acc,mean_lam,60,color_current,'filled')
    errorbar(mean_lam,mean_acc,std_acc,std_acc,std_lam,std_lam,'k.')
    scatter(mean_lam,mean_acc,60,color_current,'filled')
    xlim([-50 (sepdist+1)*50])
    ylim([-50 (sepdist+1)*50])
    axis square
    set(gca,'TickDir','out');
    plot([-50 (sepdist+1)*50],[-50 (sepdist+1)*50],'--')
    title([num2str((sepdist+1)*50) '\mum Separation between stimulating electrodes'])
    ylabel('Across peak distance from stim electrode closest to tip (\mum)')
    xlabel('Laminar peak distance from stim electrode closest to tip (\mum)')
end

%% firing rate comparison laminar and across
meanFR_acc=nan(9,5);
meanFR_lam=nan(9,5);
stdFR_acc=nan(9,5);
stdFR_lam=nan(9,5);

for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        dat_acc=[];
        dat_lam=[];
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            for layername=1:length(layer_types)
                layercheck=[layer_types(layername),layer_types(layername+length(layer_types))];
                for shanksep=0:3
                    shanksepcheck=['D' num2str(shanksep)];
                    dat_acc=[dat_acc saveCplit_raw_acc.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                    dat_lam=[dat_lam saveCplit_raw_lam.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                end
            end
        end

        meanFR_acc(sepdist,trial)=mean(dat_acc(16:16+sepdist+1,:),'all');
        stdFR_acc(sepdist,trial)=(std(dat_acc(16:16+sepdist+1,:),0 ,'all'))/sqrt(size(dat_acc(16:16+sepdist+1,:),2)*(sepdist+2));
        meanFR_lam(sepdist,trial)=mean(dat_lam(16:16+sepdist+1,:),'all');
        stdFR_lam(sepdist,trial)=(std(dat_lam(16:16+sepdist+1,:),0 ,'all','omitnan'))/sqrt(size(dat_lam(16:16+sepdist+1,:),2)*(sepdist+2));

    end
end
%%
cmap=colormap(hot);
meanFR_acc(isnan(meanFR_acc(:,1)),:)=[];
meanFR_lam(isnan(meanFR_lam(:,1)),:)=[];
stdFR_lam(isnan(stdFR_lam(:,1)),:)=[];
stdFR_acc(isnan(stdFR_acc(:,1)),:)=[];
figure
        hold on
        %errorbar(meanFR_acc,meanFR_lam, stdFR_lam, stdFR_lam,stdFR_acc,stdFR_acc,'k.')
        %color_current=cmap(trial*floor((length(cmap))/5),:);
        %scatter(meanFR_acc,meanFR_lam,60,color_current,'filled')
        color_current=cmap(([1 3 5])*floor((length(cmap))/7),:);
        for group=1:3
            gnum=['G' num2str(group)];
            errorbar(meanFR_lam(group,:), meanFR_acc(group,:),stdFR_acc(group,:),stdFR_acc(group,:),stdFR_lam(group,:), stdFR_lam(group,:),'k.')
            l.(gnum)=scatter(meanFR_lam(group,:),meanFR_acc(group,:),60,color_current(group,:),'filled');
%             errorbar(meanFR_acc(group,:),meanFR_lam(group,:), stdFR_lam(group,:), stdFR_lam(group,:),stdFR_acc(group,:),stdFR_acc(group,:),'k.')
%             l.(gnum)=scatter(meanFR_acc(group,:),meanFR_lam(group,:),60,color_current(group,:),'filled');

        end
        P = polyfit(meanFR_lam,meanFR_acc,1);
    yfit = P(1)*[0 25 50 75 100]+P(2);
    hold on;
    plot([0 25 50 75 100],yfit,'r-.');
        ylim([50 200])
        xlim([0 100])
        ylabel('Across firing rate (Sp/s)')
        xlabel('Laminar firing rate (Sp/s)')
        txt = ['y=' num2str(P(1)) 'x+' num2str(P(2))];
text(10,150,txt)
legend([l.G1 l.G2 l.G3],'300','400','500')
%% firing rate comparison laminar and across
meanFR_acc=nan(1,5);
meanFR_lam=nan(1,5);
stdFR_acc=nan(1,5);
stdFR_lam=nan(1,5);
for trial=1:5
    trialcheck=['T' num2str(trial)];
    for sepdist=5:2:9
        sepcheck=['sep' num2str(sepdist)];

        dat_acc=[];
        dat_lam=[];
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            for layername=1:length(layer_types)
                layercheck=[layer_types(layername),layer_types(layername+length(layer_types))];
                for shanksep=0:3
                    shanksepcheck=['D' num2str(shanksep)];
                    dat_acc=[dat_acc saveCplit_raw_acc.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                    dat_lam=[dat_lam saveCplit_raw_lam.(layercheck).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                end
            end
        end

       
    end
     meanFR_acc(trial)=mean(dat_acc(16:16+sepdist+1,:),'all');
        stdFR_acc(trial)=(std(dat_acc(16:16+sepdist+1,:),0 ,'all'))/sqrt(size(dat_acc(16:16+sepdist+1,:),2)*(sepdist+2));
        meanFR_lam(trial)=mean(dat_lam(16:16+sepdist+1,:),'all');
        stdFR_lam(trial)=(std(dat_lam(16:16+sepdist+1,:),0 ,'all','omitnan'))/sqrt(size(dat_lam(16:16+sepdist+1,:),2)*(sepdist+2));

end
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
meanFR_acc(isnan(meanFR_acc(:,1)),:)=[];
meanFR_lam(isnan(meanFR_lam(:,1)),:)=[];
stdFR_lam(isnan(stdFR_lam(:,1)),:)=[];
stdFR_acc(isnan(stdFR_acc(:,1)),:)=[];
figure
        hold on
        for trial=1:5
        errorbar(meanFR_acc(trial),meanFR_lam(trial), stdFR_lam(trial), stdFR_lam(trial),stdFR_acc(trial),stdFR_acc(trial),'k.')
        color_current=cmap(trial*floor((length(cmap))/5),:);
        scatter(meanFR_acc(trial),meanFR_lam(trial),60,color_current,'filled')
        end
%         color_current=cmap(([1 3 5])*floor((length(cmap))/7),:);
%         for group=1:3
%             gnum=['G' num2str(group)];
%             errorbar(meanFR_lam(group,:), meanFR_acc(group,:),stdFR_acc(group,:),stdFR_acc(group,:),stdFR_lam(group,:), stdFR_lam(group,:),'k.')
%             l.(gnum)=scatter(meanFR_lam(group,:),meanFR_acc(group,:),60,color_current(group,:),'filled');
% %             errorbar(meanFR_acc(group,:),meanFR_lam(group,:), stdFR_lam(group,:), stdFR_lam(group,:),stdFR_acc(group,:),stdFR_acc(group,:),'k.')
% %             l.(gnum)=scatter(meanFR_acc(group,:),meanFR_lam(group,:),60,color_current(group,:),'filled');
% 
%         end
%         P = polyfit(meanFR_lam,meanFR_acc,1);
%     yfit = P(1)*[0 25 50 75 100]+P(2);
%     hold on;
%     plot([0 25 50 75 100],yfit,'r-.');
%         ylim([50 200])
%         xlim([0 100])
%         ylabel('Across firing rate (Sp/s)')
%         xlabel('Laminar firing rate (Sp/s)')
%         txt = ['y=' num2str(P(1)) 'x+' num2str(P(2))];
% text(10,150,txt)
% legend([l.G1 l.G2 l.G3],'300','400','500')




%% plot layer separation
normalise_dat=0;
if normalise_dat==1
    saveCplit=saveCplit_norm;
else
    saveCplit=savelayeraligned;
end

current_all=[1 6 10]; %plot 6uA results
for itcur=1:length(current_all)
    singleCurrent=current_all(itcur);
if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end
%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
FilterLength=3;%smoothing parameters
splitsinglestimlayers.s=[];
splitsinglestimlayers.g=[];
splitsinglestimlayers.i=[];
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)-1
            l1=layer_types(layername);
            l2=layer_types(layername+length(layer_types));
            for shanksep=1:3
                shanksepcheck=['D' num2str(shanksep)];
                splitsinglestimlayers.(l1)=[splitsinglestimlayers.(l1) saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T1.(shanksepcheck)];
                splitsinglestimlayers.(l2)=[splitsinglestimlayers.(l2) saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T5.(shanksepcheck)];
            end
        end
    end
end
% splitsinglestimlayers.s(splitsinglestimlayers.s<0)=0;
% splitsinglestimlayers.g(splitsinglestimlayers.g<0)=0;
% splitsinglestimlayers.i(splitsinglestimlayers.i<0)=0;
lim1=0;
lim2=180;
figure(1)
hold on
plot(nanmean(splitsinglestimlayers.s,2),1550:-50:0,'Color', [(length(current_all)-itcur)/length(current_all) 0 0])
if singleCurrent==10
YLims_fig=get(gca, 'XLim');
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[0 0 345 345],'r','FaceAlpha', 0.1,'linestyle','none')
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[345 345 495 495],'b','FaceAlpha', 0.1,'linestyle','none')
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[495 495 1175 1175],'k','FaceAlpha', 0.1,'linestyle','none')
ylabel('Distance from cortical surface (\mum)')
xlabel('Firing rate (sp/s)')
title('Supra stim')
ylim([0 1240])
set(gca, 'YDir','reverse')
end
plot(nanmean(splitsinglestimlayers.s,2),1550:-50:0,'Color', [(length(current_all)-itcur)/length(current_all) 0 0])
xlim([lim1 lim2])

figure(2)
hold on
plot(nanmean(splitsinglestimlayers.g,2),1550:-50:0,'Color', [0 0 (length(current_all)-itcur)/length(current_all)])
if singleCurrent==10
YLims_fig=get(gca, 'XLim');
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[0 0 345 345],'r','FaceAlpha', 0.1,'linestyle','none')
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[345 345 495 495],'b','FaceAlpha', 0.1,'linestyle','none')
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[495 495 1175 1175],'k','FaceAlpha', 0.1,'linestyle','none')
ylabel('Distance from cortical surface (\mum)')
xlabel('Firing rate (sp/s)')
title('Gran stim')
ylim([0 1240])
set(gca, 'YDir','reverse')
end
plot(nanmean(splitsinglestimlayers.g,2),1550:-50:0,'Color', [0 0 (length(current_all)-itcur)/length(current_all)])
xlim([lim1 lim2])
figure(3)
hold on
plot(nanmean(splitsinglestimlayers.i,2),1550:-50:0,'Color', [(length(current_all)-itcur)/length(current_all) (length(current_all)-itcur)/length(current_all) (length(current_all)-itcur)/length(current_all)])
if singleCurrent==10
YLims_fig=get(gca, 'XLim');
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[0 0 345 345],'r','FaceAlpha', 0.1,'linestyle','none')
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[345 345 495 495],'b','FaceAlpha', 0.1,'linestyle','none')
fill([YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],[495 495 1175 1175],'k','FaceAlpha', 0.1,'linestyle','none')
ylabel('Distance from cortical surface (\mum)')
xlabel('Firing rate (sp/s)')
title('Infra stim')
ylim([0 1240])
set(gca, 'YDir','reverse')
end
plot(nanmean(splitsinglestimlayers.i,2),1550:-50:0,'Color', [(length(current_all)-itcur)/length(current_all) (length(current_all)-itcur)/length(current_all) (length(current_all)-itcur)/length(current_all)])
xlim([lim1 lim2])
end
%% layer and stim elect align

saveCplit=saveStimalignedLayeraligned;
singleCurrent=[6 8 10];
Trialcheck='T3';

if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end

for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    counter=0;
        figure(sepdist)
    for Pl=1:32
        pos=['P' num2str(Pl)];
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            %for shanksep=0:3
               % shanksepcheck=['D' num2str(shanksep)];
                %dat2plot=nanmean(saveCplit.(pos).(sepcheck).(currcheck).(Trialcheck).(shanksepcheck),2);
                dat=[saveCplit.(pos).(sepcheck).(currcheck).(Trialcheck).D0];%,saveCplit.(pos).(sepcheck).(currcheck).(Trialcheck).D1,saveCplit.(pos).(sepcheck).(currcheck).(Trialcheck).D2,saveCplit.(pos).(sepcheck).(currcheck).(Trialcheck).D3];
                max_all=max(dat);
                dat2plot=dat;%nanmean(dat,2);
                if ~isempty(dat2plot) && any(~isnan(dat2plot),'all')
                    for i=1:size(dat2plot,2)
                    hold on
                    cm=colormap;
                    plot(dat2plot(:,i)+counter,'Color',cm(round((counter)/150+1),:))
                    scatter([Pl Pl+sepdist+1],[counter counter],'MarkerEdgeColor',cm(round((counter)/150+1),:))
                    counter=counter+400;
                    yline(counter,'--')
                    end
                end
            %end
        end
    end
    YLims_fig=get(gca, 'YLim');
    fill([0 0 8 8],[YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],'g','FaceAlpha', 0.1,'linestyle','none')
    fill([8 8 21 21],[YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],'r','FaceAlpha', 0.1,'linestyle','none')
    fill([21 21 24 24],[YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],'b','FaceAlpha', 0.1,'linestyle','none')
    fill([24 24 32 32],[YLims_fig(2) YLims_fig(1) YLims_fig(1) YLims_fig(2)],'k','FaceAlpha', 0.1,'linestyle','none')
end

%% centre layer stim electrodes to central electrode for layers
normalise_dat=1;%normalise?
%setup data structures

for current=1:length(AMP)
    currcheck=['C' num2str(AMP(current))];
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        for layername=1:length(layer_types)
            saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=[];
        end
    end
end

if normalise_dat==1
    saveCplit=saveCplit_norm;
else
    saveCplit=saveCplit_raw;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2

%sort data and merge sep dists by aligning to central electrode
for trial=1:5
    trialcheck=['T' num2str(trial)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                sep5dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep5').(currcheck).(trialcheck).(shanksepcheck);
                sep7dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep7').(currcheck).(trialcheck).(shanksepcheck);
                sep9dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep9').(currcheck).(trialcheck).(shanksepcheck);
                if ~isempty(sep5dat) && ~isempty(sep7dat) && ~isempty(sep9dat)
                    saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=[saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck) sep5dat [sep7dat(2:end,:); nan(1,size(sep7dat,2))] [sep9dat(3:end,:); nan(2,size(sep9dat,2))]];
                end
            end
        end
    end
end


for current=1:length(AMP)
    currcheck=['C' num2str(AMP(current))];
    figure
    for layername=1:length(layer_types)
        for trial=1:5
            trialcheck=['T' num2str(trial)];
            layernameFULL=[layer_types(layername),layer_types(layername+length(layer_types))];
           smoothedenvelope=smoothData(saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck),FilterLength);
            if trial==1 && ~isempty(smoothedenvelope)
                figure %subplot(2,4,layername)
                hold on
                ax1.(layernameFULL)=gca;
                hold(ax1.(layernameFULL), 'on');
                axes('Position',[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2)+ax1.(layernameFULL).Position(4)*3/4 ax1.(layernameFULL).Position(3) ax1.(layernameFULL).Position(4)/4])
                ax2.(layernameFULL)=gca;
                title(layernameFULL)
                hold(ax2.(layernameFULL),'off')
                ax1.(layernameFULL).Position=[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3) ax1.(layernameFULL).Position(4)*3/4];
            end
            if ~isempty(smoothedenvelope)
                [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,0,trial,cmap,ax1.(layernameFULL),ax2.(layernameFULL));
                peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=peak;
                if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1)
                    [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T2(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T3(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T4(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T5(:)],1,'off');
                    text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),DATApos1+1,1,['p=' num2str(p)]) 
                end
            end
        end
    end
end

%% centre layer stim electrodes for layers split into two graphs
%stim1 centralise and stim2 centralise
normalise_dat=1;%normalise?
%setup data structures
clear saveCplit_centralisedallsepdist
for current=1:length(AMP)
    currcheck=['C' num2str(AMP(current))];
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        for layername=1:length(layer_types)
            saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).('stim1')=[];
            saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).('stim2')=[];
        end
    end
end

if normalise_dat==1
    saveCplit=saveCplit_norm;
else
    saveCplit=saveCplit_raw;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2
FilterLength=3;
%sort data and merge sep dists by aligning to central electrode
for trial=1:5
    trialcheck=['T' num2str(trial)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                sep5dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep5').(currcheck).(trialcheck).(shanksepcheck);
                sep7dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep7').(currcheck).(trialcheck).(shanksepcheck);
                sep9dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep9').(currcheck).(trialcheck).(shanksepcheck);
                if ~isempty(sep5dat) && ~isempty(sep7dat) && ~isempty(sep9dat)
                    sep5dat=smoothData(sep5dat,FilterLength);
                    sep7dat=smoothData(sep7dat,FilterLength);
                    sep9dat=smoothData(sep9dat,FilterLength);

                    sep5dat1=[sep5dat(1:19,:); nan(13,size(sep5dat,2))];
                    sep7dat1=[sep7dat(1:20,:); nan(12,size(sep7dat,2))];
                    sep9dat1=[sep9dat(1:21,:); nan(11,size(sep9dat,2))];
                    sep5dat2=[nan(19,size(sep5dat,2)); sep5dat(20:end,:)];
                    sep7dat2=[nan(18,size(sep7dat,2)); sep7dat(21:end,:); nan(2,size(sep7dat,2))];
                    sep9dat2=[nan(17,size(sep9dat,2)); sep9dat(22:end,:); nan(4,size(sep9dat,2))];

                    saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).('stim1')=[saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).('stim1') sep5dat1 sep7dat1 sep9dat1];
                    saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).('stim2')=[saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).('stim2') sep5dat2 sep7dat2 sep9dat2];
                end
            end
        end
    end
end

for stime=1:2
    stimcheck=['stim' num2str(stime)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)
            for trial=1:5
                trialcheck=['T' num2str(trial)];
                layernameFULL=[layer_types(layername),layer_types(layername+length(layer_types))];
                dat2plot=saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck).(stimcheck);
                if trial==1 && ~isempty(dat2plot)
                    figure %subplot(2,4,layername)
                    hold on
                    ax1.(layernameFULL)=gca;
                    hold(ax1.(layernameFULL), 'on');
                    axes('Position',[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2)+ax1.(layernameFULL).Position(4)*3/4 ax1.(layernameFULL).Position(3) ax1.(layernameFULL).Position(4)/4])
                    ax2.(layernameFULL)=gca;
                    title(layernameFULL)
                    hold(ax2.(layernameFULL),'off')
                    ax1.(layernameFULL).Position=[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3) ax1.(layernameFULL).Position(4)*3/4];
                end
                if ~isempty(dat2plot)
                    [peak,DATApos1]=plotActivityAndPeak(dat2plot,5,trial,cmap,ax1.(layernameFULL),ax2.(layernameFULL));
                    peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=peak;
                    if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1)
                        [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T2(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T3(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T4(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T5(:)],1,'off');
                        if stime==1
                        text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),18,1,['p=' num2str(p)])
                        else
                            text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),24,1,['p=' num2str(p)])
                        end
                    end
                end
            end
        end
    end
end

%% centre layer stim electrodes and remove electrodes either random, central or outer
%stim1 centralise and stim2 centralise
normalise_dat=1;%normalise?
singleCurrent=6; %plot 6uA results
if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end
Lam_Acc='across';
%setup data structures
clear saveCplit_centralisedallsepdist
for current=1:length(AMP)
    currcheck=['C' num2str(AMP(current))];
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        for layername=1:length(layer_types)
            saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=[];
        end
    end
end

if normalise_dat==1
    saveCplit=saveCplit_norm;
else
    saveCplit=saveCplit_raw;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
seedpoint=20;
s = RandStream('mlfg6331_64','Seed',seedpoint);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2
FilterLength=3;
%sort data and merge sep dists by aligning to central electrode
for trial=1:5
    trialcheck=['T' num2str(trial)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                sep5dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep5').(currcheck).(trialcheck).(shanksepcheck);
                sep7dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep7').(currcheck).(trialcheck).(shanksepcheck);
                sep9dat=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).('sep9').(currcheck).(trialcheck).(shanksepcheck);
                if ~isempty(sep5dat) && ~isempty(sep7dat) && ~isempty(sep9dat)
%                     sep5dat=smoothData(sep5dat,FilterLength);
%                     sep7dat=smoothData(sep7dat,FilterLength);
%                     sep9dat=smoothData(sep9dat,FilterLength);
                    
                    %remove centre electrodes
%                     sep7dat=[sep7dat(1:19,:); sep7dat(22:end,:); nan(2,size(sep7dat,2))];
%                     sep9dat=[sep9dat(1:19,:); sep9dat(24:end,:); nan(4,size(sep9dat,2))];
%                     
                    %remove outer electrodes (in centre closest to the stim
                    %positions
%                     sep7dat=[sep7dat(1:16,:); sep7dat(18:22,:); sep7dat(24:end,:); nan(2,size(sep7dat,2))];
%                     sep9dat=[sep9dat(1:16,:); sep9dat(19:23,:); sep9dat(26:end,:); nan(4,size(sep9dat,2))];

                    %remove random electrodes between stim electrodes
                    sep7dat(randsample(s,7,2,'false')+16,:)=[];
                    sep9dat(randsample(s,9,4,'false')+16,:)=[];
                    sep7dat=[sep7dat; nan(2,size(sep7dat,2))];
                    sep9dat=[sep9dat; nan(4,size(sep9dat,2))];
               
                    saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=[saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck) sep5dat sep7dat sep9dat];
                end
            end
        end
    end
end
for current=1:length(AMP)
    currcheck=['C' num2str(AMP(current))];
    figure
    for layername=1:length(layer_types)
        for trial=1:5
            trialcheck=['T' num2str(trial)];
            layernameFULL=[layer_types(layername),layer_types(layername+length(layer_types))];
           smoothedenvelope=smoothData(saveCplit_centralisedallsepdist.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck),FilterLength);
            if trial==1 && ~isempty(smoothedenvelope)
               figure%subplot(2,4,layername)
                    hold on
                    if strcmp(Lam_Acc,'laminar')
                        hold on
                        ax1.(layernameFULL)=gca;
                        hold(ax1.(layernameFULL), 'on');
                        axes('Position',[ax1.(layernameFULL).Position(1)+ax1.(layernameFULL).Position(3)*3/4 ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3)/4 ax1.(layernameFULL).Position(4)])
                        ax2.(layernameFULL)=gca;
                        title(layernameFULL)
                        hold(ax2.(layernameFULL),'off')
                        ax1.Position=[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3)*3/4 ax1.(layernameFULL).Position(4)];
                    else
                        hold on
                        ax2.(layernameFULL)=gca;
                        hold(ax2.(layernameFULL), 'on');
                        axes('Position',[ax2.(layernameFULL).Position(1) ax2.(layernameFULL).Position(2) ax2.(layernameFULL).Position(3) ax2.(layernameFULL).Position(4)*3/4])
                        ax1.(layernameFULL)=gca;
                        title(layernameFULL)
                        hold(ax1.(layernameFULL),'off')
                        ax2.(layernameFULL).Position=[ax2.(layernameFULL).Position(1) ax2.(layernameFULL).Position(2)+ax2.(layernameFULL).Position(4)*3/4 ax2.(layernameFULL).Position(3) ax2.(layernameFULL).Position(4)/4];
                    end            
            end
            if ~isempty(smoothedenvelope)
                [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,5,trial,cmap,ax1.(layernameFULL),ax2.(layernameFULL),Lam_Acc);
                peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=peak;
                if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1)
                    [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T2(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T3(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T4(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T5(:)],1,'off');
                    text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),DATApos1+9.5,1,['p=' num2str(p)]) 
                    N=size(peak,1)
                end
            end
        end
    end
end
%% plot single lines
%YOU CAN ONLY PICK ONE OF THESE
splitshanks=0; %Split into stim, 1 shank away and 2 shanks away etc
splitlayers=0; %split into layers
%AND IF YOU PICK ONE^^^, YOU CAN'T PICK ALL CURRENTS
singleCurrent=6; %plot 6uA results
%what do you want to class as indipendent samples
avgpairsamplesplot=0;

if splitshanks==1 || splitlayers==1
    downsampleYN=0; %do you want to downsample 1=yes %%%currently not implemented
    numshanksToAvg=1;
elseif skip_stimshank==1
    numshanksToAvg=3;
    downsampleYN=1; %do you want to downsample 1=yes
else
    numshanksToAvg=4;
    downsampleYN=1; %do you want to downsample 1=yes
end

if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2
shanknames={'Stimshank'; '1 shank away'; '2 shanks away'; '3 shanks away'};
FilterLength=3;%smoothing parameters
num=0;
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    dat2plot=[];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for layername=1:length(layer_types)
            for shanksep=0:3
                if shanksep==0 && skip_stimshank==1
                    continue
                end
                shanksepcheck=['D' num2str(shanksep)];

                Data_allcols1=smoothData(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T1').(shanksepcheck),FilterLength);
                Data_allcols2=smoothData(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T2').(shanksepcheck),FilterLength);
                Data_allcols3=smoothData(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T3').(shanksepcheck),FilterLength);
                Data_allcols4=smoothData(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T4').(shanksepcheck),FilterLength);
                Data_allcols5=smoothData(saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).('T5').(shanksepcheck),FilterLength);
                max_all=max([Data_allcols1; Data_allcols2; Data_allcols3; Data_allcols4; Data_allcols5]);
                if isempty(Data_allcols1)
                    continue
                end
                for i=1:size(Data_allcols1,2)
                    figure
                    %axes('Position',[0.13         0.112396822842341                     0.775         0.62])
                    ax1=gca;
                    hold(ax1, 'on');
                    axes('Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)*3/4 ax1.Position(3) ax1.Position(4)/4])
                    ax2=gca;
                    hold(ax2,'off')
                    ax1.Position=[ax1.Position(1) ax1.Position(2) ax1.Position(3) ax1.Position(4)*3/4];
                    mall=max([Data_allcols1(:,i); Data_allcols2(:,i); Data_allcols3(:,i); Data_allcols4(:,i); Data_allcols5(:,i)]);
                    plotActivityAndPeak(Data_allcols1(:,i)./mall,sepdist,1,cmap,ax1,ax2)
                    plotActivityAndPeak(Data_allcols2(:,i)./mall,sepdist,2,cmap,ax1,ax2)
                    plotActivityAndPeak(Data_allcols3(:,i)./mall,sepdist,3,cmap,ax1,ax2)
                    plotActivityAndPeak(Data_allcols4(:,i)./mall,sepdist,4,cmap,ax1,ax2)
                    plotActivityAndPeak(Data_allcols5(:,i)./mall,sepdist,5,cmap,ax1,ax2)
                    num=num+1;
                    savename=['fig' num2str(num) '.png'];
                    saveas(gcf,savename)
                end
            end
        end
    end
end

%% sort across current steering data
skip_stimshank=0;

startpointseconds=2;
secondstoanalyse=8;
AMP=[0 1 2 3 4 6 8 10];
layer_types=['ss';'gg';'ii';'si';'sg';'gi';'WM'];
%setup data structures
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        for trial=1:5
            trialcheck=['T' num2str(trial)];
            for shanksep=0:3
                shanksepcheck=['D' num2str(shanksep)];
                saveCplit.(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
                savecsplit_avgpairs.(sepcheck).(currcheck).(trialcheck)=[];
            end
        end
    end
end

for ratN=14:20%[6:13 21:23] %loop through animals
    %load data
    if ratN<10
        Ratnum=['Rat_00' num2str(ratN)];
    elseif ratN>=10
        Ratnum=['Rat_0' num2str(ratN)];
    end

    cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
    D_data=dir;

    for k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
        currD = D_data(k).name; % Get the current subdirectory name
        try
            cd([D_data(k).folder filesep currD])
            load('Averagetrialresponse.mat','avgnospT') % load sorted neural activity
        catch
            continue;
        end
        loadStimChn;
        for sepdist=5:2:9
            if (stimChn(1)-stimChn(2))~=(-1*sepdist)-1 % check stimchn sepdist between these electrodes
                continue
            end
            [Csplit_shankdist]=PoolNormalisedActivity_refactored(avgnospT,startpointseconds, secondstoanalyse);
            sepcheck=['sep' num2str(sepdist)];
            AMP=loadAMP;
            AMP(AMP==-1)=0;
            for current=1:length(AMP)
                currcheck=['C' num2str(AMP(current))];
                for trial=1:5
                    trialcheck=['T' num2str(trial)];
                    for shanksep=0:3
                        shanksepcheck=['D' num2str(shanksep)];
                        saveCplit.(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[saveCplit.(sepcheck).(currcheck).(trialcheck).(shanksepcheck) Csplit_shankdist.(currcheck).(trialcheck).(shanksepcheck)];
                    end
                    if skip_stimshank==1
                        tmp=nanmean([Csplit_shankdist.(currcheck).(trialcheck).D1 Csplit_shankdist.(currcheck).(trialcheck).D2 Csplit_shankdist.(currcheck).(trialcheck).D3],2);
                    else
                        tmp=nanmean([Csplit_shankdist.(currcheck).(trialcheck).D0 Csplit_shankdist.(currcheck).(trialcheck).D1 Csplit_shankdist.(currcheck).(trialcheck).D2 Csplit_shankdist.(currcheck).(trialcheck).D3],2);
                    end
                    savecsplit_avgpairs.(sepcheck).(currcheck).(trialcheck)=[savecsplit_avgpairs.(sepcheck).(currcheck).(trialcheck) tmp]; %already averaged pairs
                end
            end
        end
    end
end
saveCplit_raw_across=saveCplit;
%% normalise across
saveCplit=saveCplit_raw_across;
AMP=[0 1 2 3 4 6 8 10];
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    dat2plot=[];
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];

        for shanksep=0:3
            if shanksep==0 && skip_stimshank==1
                continue
            end
            shanksepcheck=['D' num2str(shanksep)];
            Data_allcols1=saveCplit.(sepcheck).(currcheck).('T1').(shanksepcheck);
            Data_allcols2=saveCplit.(sepcheck).(currcheck).('T2').(shanksepcheck);
            Data_allcols3=saveCplit.(sepcheck).(currcheck).('T3').(shanksepcheck);
            Data_allcols4=saveCplit.(sepcheck).(currcheck).('T4').(shanksepcheck);
            Data_allcols5=saveCplit.(sepcheck).(currcheck).('T5').(shanksepcheck);
            max_all=max([Data_allcols1; Data_allcols2; Data_allcols3; Data_allcols4; Data_allcols5]);
            saveCplit_norm.(sepcheck).(currcheck).('T1').(shanksepcheck)=Data_allcols1./max_all;
            saveCplit_norm.(sepcheck).(currcheck).('T2').(shanksepcheck)=Data_allcols2./max_all;
            saveCplit_norm.(sepcheck).(currcheck).('T3').(shanksepcheck)=Data_allcols3./max_all;
            saveCplit_norm.(sepcheck).(currcheck).('T4').(shanksepcheck)=Data_allcols4./max_all;
            saveCplit_norm.(sepcheck).(currcheck).('T5').(shanksepcheck)=Data_allcols5./max_all;
        end
    end

end


%% plotting acrosss
%normalise? pick 1 or none to plot normal data
normalise_dat=1;

splitshanks=0; %Split into stim, 1 shank away and 2 shanks away etc
%AND IF YOU PICK THIS^^^, YOU CAN'T PICK ALL CURRENTS
singleCurrent=6; %plot 6uA results
%what do you want to class as indipendent samples
avgpairsamplesplot=0;
if normalise_dat==1
    saveCplit=saveCplit_norm;
else
    saveCplit=saveCplit_raw_across;
end

if splitshanks==1
    downsampleYN=0; %do you want to downsample 1=yes %%%currently not implemented
    numshanksToAvg=1;
elseif skip_stimshank==1
    numshanksToAvg=3;
    downsampleYN=1; %do you want to downsample 1=yes
else
    numshanksToAvg=4;
    downsampleYN=1; %do you want to downsample 1=yes
end

if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2
shanknames={'Stimshank'; '1 shank away'; '2 shanks away'; '3 shanks away'};
FilterLength=3;%smoothing parameters
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    figure
    if ~splitshanks
        %axes('Position',[0.13         0.112396822842341                     0.775         0.62])
        ax1=gca;
        hold(ax1, 'on');
        axes('Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)*3/4 ax1.Position(3) ax1.Position(4)/4])
        ax2=gca;
        hold(ax2,'off')
        ax1.Position=[ax1.Position(1) ax1.Position(2) ax1.Position(3) ax1.Position(4)*3/4];
    end
    clear peak_all
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        dat2plot=[];
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            for shanksep=0:3
                if shanksep==0 && skip_stimshank==1
                    continue
                end
                shanksepcheck=['D' num2str(shanksep)];
                if splitshanks==1
                    dat2plot.(shanksepcheck)=saveCplit.(sepcheck).(currcheck).(trialcheck).(shanksepcheck);
                    if trial==1
                        subplot(1,4,shanksep+1)
                        ax1.(shanksepcheck)=gca;
                        hold(ax1.(shanksepcheck), 'on');
                        axes('Position',[ax1.(shanksepcheck).Position(1) ax1.(shanksepcheck).Position(2)+ax1.(shanksepcheck).Position(4)*3/4 ax1.(shanksepcheck).Position(3) ax1.(shanksepcheck).Position(4)/4])
                        ax2.(shanksepcheck)=gca;
                        title(shanknames{shanksep+1})
                        hold(ax2.(shanksepcheck),'off')
                        ax1.(shanksepcheck).Position=[ax1.(shanksepcheck).Position(1) ax1.(shanksepcheck).Position(2) ax1.(shanksepcheck).Position(3) ax1.(shanksepcheck).Position(4)*3/4];
                    end
                    dat2plot.(shanksepcheck)((sum(~isnan(dat2plot.(shanksepcheck)),2)<size(dat2plot.(shanksepcheck),2)/2),:)=nan;
                    smoothedenvelope=smoothData(dat2plot.(shanksepcheck),FilterLength);
                    [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1.(shanksepcheck),ax2.(shanksepcheck));
                    peak_all.(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=peak;
                    if trial==5 && exist('peak_all') && ~isempty(peak_all.(sepcheck).(currcheck).(trialcheck).(shanksepcheck))
                        [p,tbl,stats]=friedman([peak_all.(sepcheck).(currcheck).(trialcheck).T1(:),peak_all.(sepcheck).(currcheck).(trialcheck).T2(:),peak_all.(sepcheck).(currcheck).(trialcheck).T3(:),peak_all.(sepcheck).(currcheck).(trialcheck).T4(:),peak_all.(sepcheck).(currcheck).(trialcheck).T5(:)],1,'off');
                        text(ax2.(shanksepcheck),find(~isnan(lineOut.YData),1,'first')+1,DATApos1+1,1,['p=' num2str(p)])
                    end
                elseif avgpairsamplesplot==0
                    dat2plot=[dat2plot saveCplit.(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                    continue
                elseif avgpairsamplesplot==1 && ((shanksep==0 && skip_stimshank==0) || (shanksep==1 && skip_stimshank==1))
                    dat2plot=[dat2plot savecsplit_avgpairs.(sepcheck).(currcheck).(trialcheck)];
                end
            end
        end
        if splitshanks==1
            continue
        elseif singleCurrent==0
            subplot(2,4,current)
        end
        hold on
        dat2plot((sum(~isnan(dat2plot),2)<size(dat2plot,2)/2),:)=nan;
        smoothedenvelope=smoothData(dat2plot,FilterLength);
        [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1,ax2);
        peak_all.(currcheck).(trialcheck)=peak;
        if trial==5 && exist('peak_all') && ~isempty(peak_all.(currcheck).(trialcheck))
            [p,tbl,stats]=friedman([peak_all.(currcheck).T1(:),peak_all.(currcheck).T2(:),peak_all.(currcheck).T3(:),peak_all.(currcheck).T4(:),peak_all.(currcheck).T5(:)],1,'off');
            text(ax2,DATApos1+9.5,1,['p=' num2str(p)])
        end
        dat2plot=[];
    end
end


%% plot laminar horizontal for comparison across
%%need to rewrite plotActivityAndPeak function for rotated axes. then
%%adjust ax1 and ax2 and text options

%simulation or normalise? pick 1 or none to plot normal data
plotlayeraligned=0;
plotsim=0;
normalise_dat=1;

%YOU CAN ONLY PICK ONE OF THESE
splitshanks=0; %Split into stim, 1 shank away and 2 shanks away etc
splitlayers=0; %split into layers
%AND IF YOU PICK ONE^^^, YOU CAN'T PICK ALL CURRENTS
singleCurrent=6; %plot 6uA results
FilterLength=3;%smoothing parameters
%what do you want to class as indipendent samples
avgpairsamplesplot=0;
if normalise_dat==1
    saveCplit=saveCplit_norm;
elseif plotsim==1
    saveCplit=Csplit_simulation;
elseif plotlayeraligned==1
    saveCplit=savelayeraligned;
else
    saveCplit=saveCplit_raw;
end


if splitshanks==1 || splitlayers==1
    downsampleYN=0; %do you want to downsample 1=yes %%%currently not implemented
    numshanksToAvg=1;
elseif skip_stimshank==1
    numshanksToAvg=3;
    downsampleYN=1; %do you want to downsample 1=yes
else
    numshanksToAvg=4;
    downsampleYN=1; %do you want to downsample 1=yes
end

if singleCurrent==0
    AMP=[0 1 2 3 4 6 8 10];
else
    AMP=singleCurrent;
end

%setup colours and downsample
vec = [100;80;50;30;15;0];
hex = ['#2E1658';'#884DFF';'#14AAF5';'#00C213';'#FFF700';'#FFCC00'];
N = 128;
seedpoint=65;
s = RandStream('mlfg6331_64','Seed',seedpoint);
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
map = interp1(vec,raw,linspace(100,0,N),'pchip');
cmap=colormap(map);
clear ax1 ax2
shanknames={'Stimshank'; '1 shank away'; '2 shanks away'; '3 shanks away'};

for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    figure
    if ~splitshanks && ~splitlayers
        %axes('Position',[0.13         0.112396822842341                     0.775         0.62])
        hold on
        ax1=gca;
        hold(ax1, 'on');
        axes('Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)*3/4 ax1.Position(3) ax1.Position(4)/4])
        ax2=gca;
        hold(ax2,'off')
        ax1.Position=[ax1.Position(1) ax1.Position(2) ax1.Position(3) ax1.Position(4)*3/4];
    end
    clear peak_all
    for trial=1:5
        trialcheck=['T' num2str(trial)];
        dat2plot=[];
        for current=1:length(AMP)
            currcheck=['C' num2str(AMP(current))];
            for layername=1:length(layer_types)
                for shanksep=0:3
                    if shanksep==0 && skip_stimshank==1
                        continue
                    end
                    shanksepcheck=['D' num2str(shanksep)];
                    if splitlayers==1 && avgpairsamplesplot==0
                        dat2plot=[dat2plot saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                        continue
                    elseif splitlayers==1 && avgpairsamplesplot==1 && ((shanksep==0 && skip_stimshank==0) || (shanksep==1 && skip_stimshank==1))
                        dat2plot=[dat2plot savecsplit_avgpairs.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)];
                        continue
                    elseif splitshanks==1 && layername==1
                        dat2plot.(shanksepcheck)=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck);
                    elseif splitshanks==1 && layername==length(layer_types)
                        %plot
                        if trial==1
                            subplot(1,4,shanksep+1)
                            ax1.(shanksepcheck)=gca;
                            hold(ax1.(shanksepcheck), 'on');
                            axes('Position',[ax1.(shanksepcheck).Position(1) ax1.(shanksepcheck).Position(2)+ax1.(shanksepcheck).Position(4)*3/4 ax1.(shanksepcheck).Position(3) ax1.(shanksepcheck).Position(4)/4])
                            ax2.(shanksepcheck)=gca;
                            title(shanknames{shanksep+1})
                            hold(ax2.(shanksepcheck),'off')
                            ax1.(shanksepcheck).Position=[ax1.(shanksepcheck).Position(1) ax1.(shanksepcheck).Position(2) ax1.(shanksepcheck).Position(3) ax1.(shanksepcheck).Position(4)*3/4];
                        end
                        dat2plot.(shanksepcheck)((sum(~isnan(dat2plot.(shanksepcheck)),2)<size(dat2plot.(shanksepcheck),2)/2),:)=nan;
                        smoothedenvelope=smoothData(dat2plot.(shanksepcheck),FilterLength);
                        [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1.(shanksepcheck),ax2.(shanksepcheck));
                        peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=peak;
                        if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck))
                            [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T1.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T2.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T3.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T4.(shanksepcheck)(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).T5.(shanksepcheck)(:)],1,'off');
                            text(ax2.(shanksepcheck),find(~isnan(lineOut.YData),1,'first')+1,DATApos1+9.5,1,['p=' num2str(p)])
                        end
                    elseif splitshanks==1 && layername>1 && layername<length(layer_types)
                        dat2plot.(shanksepcheck)=[dat2plot.(shanksepcheck) saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                    elseif avgpairsamplesplot==0
                        dat2plot=[dat2plot saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)];
                        continue
                    elseif avgpairsamplesplot==1 && ((shanksep==0 && skip_stimshank==0) || (shanksep==1 && skip_stimshank==1))
                        dat2plot=[dat2plot savecsplit_avgpairs.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck)];
                    end
                end
                if splitlayers==0
                    continue
                end
                layernameFULL=[layer_types(layername),layer_types(layername+length(layer_types))];
                if trial==1
                    subplot(2,4,layername)
                    hold on
                    ax1.(layernameFULL)=gca;
                    hold(ax1.(layernameFULL), 'on');
                    axes('Position',[ax1.(layernameFULL).Position(1)+ax1.(layernameFULL).Position(3)*3/4 ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3)/4 ax1.(layernameFULL).Position(4)])
                    ax2.(layernameFULL)=gca;
                    title(layernameFULL)
                    hold(ax2.(layernameFULL),'off')
                    ax1.(layernameFULL).Position=[ax1.(layernameFULL).Position(1) ax1.(layernameFULL).Position(2) ax1.(layernameFULL).Position(3)*3/4 ax1.(layernameFULL).Position(4)];
                end
                if isempty(dat2plot)
                    continue
                end
                dat2plot((sum(~isnan(dat2plot),2)<size(dat2plot,2)/2),:)=nan;
                smoothedenvelope=smoothData(dat2plot,FilterLength);
                [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1.(layernameFULL),ax2.(layernameFULL));
                peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).(trialcheck)=peak;
                if trial==5 && exist('peak_all') && ~isempty(peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1)
                    [p,tbl,stats]=friedman([peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T1(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T2(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T3(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T4(:),peak_all.([layer_types(layername),layer_types(layername+length(layer_types))]).(currcheck).T5(:)],1,'off');
                    %text(ax2.([layer_types(layername),layer_types(layername+length(layer_types))]),1,DATApos1+13.5,['p=' num2str(p)])
                end
                dat2plot=[];
            end
            if splitlayers==1 || splitshanks==1
                continue
            elseif singleCurrent==0
                subplot(2,4,current)
            end
            hold on
            dat2plot((sum(~isnan(dat2plot),2)<size(dat2plot,2)/2),:)=nan;
            smoothedenvelope=smoothData(dat2plot,FilterLength);
            [peak,DATApos1]=plotActivityAndPeak(smoothedenvelope,sepdist,trial,cmap,ax1,ax2);
            peak_all.(currcheck).(trialcheck)=peak;
            if trial==5 && exist('peak_all') && ~isempty(peak_all.(currcheck).(trialcheck))
                [p,tbl,stats]=friedman([peak_all.(currcheck).T1(:),peak_all.(currcheck).T2(:),peak_all.(currcheck).T3(:),peak_all.(currcheck).T4(:),peak_all.(currcheck).T5(:)],1,'off');
                text(ax2,DATApos1+9.5,1,['p=' num2str(p)])
                p
            end
            dat2plot=[];
        end

    end
end


%%

%%%
%useless crap
%saveshanksepdist.(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
%saveCplit.(sepcheck).(currcheck).(trialcheck)=[];
% %saveCplit.(sepcheck).(currcheck).(trialcheck)=[saveCplit.(sepcheck).(currcheck).(trialcheck) Csplit.(currcheck).(trialcheck)];
% %saveCplit_layers.(layers_stim).(sepcheck).(currcheck).(trialcheck)=[];
%                     %%%%%%%%%
%                      Data=saveCplit.([layer_types(layername),layer_types(layername+length(layer_types))]).(sepcheck).(currcheck).(trialcheck).(shanksepcheck);
% if isempty(Data)
%     continue
% end
%                     for i=1:size(Data,2)
%                     figure
%                     if ~splitshanks && ~splitlayers
%                         %axes('Position',[0.13         0.112396822842341                     0.775         0.62])
%                         ax1=gca;
%                         hold(ax1, 'on');
%                         axes('Position',[ax1.Position(1) ax1.Position(2)+ax1.Position(4)*3/4 ax1.Position(3) ax1.Position(4)/4])
%                         ax2=gca;
%                         hold(ax2,'off')
%                         ax1.Position=[ax1.Position(1) ax1.Position(2) ax1.Position(3) ax1.Position(4)*3/4];
%                     end
%                     plotActivityAndPeak(Data(:,i),sepdist,trial,cmap,ax1,ax2)
%                     end
%                     continue
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%