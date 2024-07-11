%% need to use the simulation portion of LaminarCurrentSteering_refactored with this code, only one layer and needs to be capable of multiple pulses as well as multiple current levels. Only one electrode stimulating
%%simulation
AMP=[0 1 2 3 4 6 8 10];

% reading connectome data - https://bbp.epfl.ch/nmc-portal/downloads.html


fileName = 'C:\Matlab_dual_stim\Matlab_dual_stim\pathways_anatomy_factsheets_simplified.json';%'C:\Users\smei0006\Documents\Experimental_Design\pathways_anatomy_factsheets_simplified.json'; % filename in JSON extension synapse count
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data_L = jsondecode(str); % Using the jsondecode function to parse JSON from string
fn_L = fieldnames(data_L);
from_L_str=[{'1'} {'2'} {'4'} {'5'} {'6'}];
to_L_str=[{'1'} {'2'} {'4'} {'5'} {'6'}];

fileName = 'C:\Matlab_dual_stim\Matlab_dual_stim\pathways_physiology_factsheets_simplified.json'; %C:\Users\smei0006\Documents\Experimental_Design\pathways_physiology_factsheets_simplified.json'; % filename in JSON extension excite/inhib
fid = fopen(fileName); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data_anatomy_L = jsondecode(str); % Using the jsondecode function to parse JSON from string
names_anatomy_L = fieldnames(data_anatomy_L);
firing_rates_EI=[1.09 6.00];%https://www.sciencedirect.com/science/article/pii/S0092867415011915#sec4

filename = 'C:\Matlab_dual_stim\Matlab_dual_stim\layer_download.json'; %C:\Users\smei0006\Documents\Experimental_Design\layer_download.json'; % filename in JSON extension
fid = fopen(filename); % Opening the file
raw = fread(fid,inf); % Reading the contents
str = char(raw'); % Transformation
fclose(fid); % Closing the file
data_layer = jsondecode(str); % Using the jsondecode function to parse JSON from string
fn_layer = fieldnames(data_layer);

%need to know the total number of neurons in all layers of the microcircuit
%number of each type of neuron is in data_layer e.g. data_layer.L1.No_OfNeuronsPerMorphologicalTypes.L1_DAC is the number of L1_DAC neurons while L1 indicates layer 1. There are 5 layers as layer 2 and 3 are combined into L23
%and varying numbers of neurons in each layer
total_neurons=0;
layernames={'L1','L23','L4','L5','L6'};
layernumbers=[1 23 4 5 6];
neurontypenumber=[];
neuron_types={};
for i=layernumbers
    layer=['L' num2str(i)];
    fn_layer2=fieldnames(data_layer.(layer).No_OfNeuronsPerMorphologicalTypes);
    for j=1:length(fn_layer2)
        total_neurons=total_neurons+data_layer.(layer).No_OfNeuronsPerMorphologicalTypes.(fn_layer2{j});
        neurontypenumber=[neurontypenumber data_layer.(layer).No_OfNeuronsPerMorphologicalTypes.(fn_layer2{j})]; % this stores the numbers of each neuron so that they can be evenly distributed in the cube space
        neuron_types=[neuron_types fn_layer2{j}];
    end
end

% Generate a list of neuron indices and shuffle it
neuron_list = [];
for i = 1:length(neuron_types)
    neuron_list = [neuron_list; repmat(neuron_types(i), neurontypenumber(i), 1)];
end
neuron_list = neuron_list(randperm(length(neuron_list)));

%now need to loop through the connection information and create a cubic array of coordinates pointing to connected neurons e.g. if neuron 1 is connected to neuron 2, then cubic_array{1,1,1} will contain the coordinates of neuron 2
% connection information is in data_L and data_anatomy_L
%need to consider neighbor bias, connection probability, divergence and convergence
%neighbor bias is under common_neighbor_bias and  is the probability of a neuron connecting to a neighboring neuron and will be something like 2.4 which indicates it is 2.4x more likely to connect to a neighboring neuron than a non-neighboring neuron, something close to 1 indicates no bias
%connection probability is under connection_probability and is the probability of a neuron connecting to another neuron and will be something like 3.1 which indicates a 3.1% chance of connecting to another neuron
%divergence is under number_of_divergent_neuron_mean which will mean if the pre-synaptic neuron does make a connection based on the connection probability, it will connect to this many post-synaptic neurons of that type
%convergence is under number_of_convergent_neuron_mean which will mean if the post-synaptic neuron does make a connection based on the connection probability, it will connect to this many pre-synaptic neurons of that type
%each property is in data_L e.g. data_L.L4_PC___L6_TPC_L1.common_neighbor_bias is the neighbor bias between L4_PC and L6_TPC
connection_array = cell(length(neuron_list), 1);
Synapse_types_array = cell(length(neuron_list), 1);
timeIEDF=zeros(2,2); % this is the time that the connection is active for inhibitory and excitatory neurons, depressing and facilitating. Depressing means the neuron is less likely to fire, facilitating means the neuron is more likely to fire if there is a consecutive pulse
totalIEDF=zeros(2,3); % this is the total number of connections for inhibitory and excitatory neurons, depressing and facilitating
timedecay=zeros(2,3); % this is the time that it takes the post-synaptic potential to decay from peak to baseline
timelatency=zeros(2,3); %the onset time measured as the difference between the time to peak of the presynaptic AP and time taken to reach 5% of peak PSP amplitude i.e. latency between pre and post synaptic neuron
countsynapsetype.ID=0;
countsynapsetype.ED=0;
countsynapsetype.IF=0;
countsynapsetype.EF=0;
countsynapsetype.IP=0;
countsynapsetype.EP=0;
for neuronconnectiontypiteration=1:length(fn_L)
    connection_probability=data_L.(fn_L{neuronconnectiontypiteration}).connection_probability;
    presynaptic_neuron=fn_L{neuronconnectiontypiteration}(1:strfind(fn_L{neuronconnectiontypiteration},'___')-1);
    postsynaptic_neuron=fn_L{neuronconnectiontypiteration}(strfind(fn_L{neuronconnectiontypiteration},'___')+3:end);
    Synapse_type=data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).synapse_type;
    %need to save mean time of all synapse types d_mean if depressing, f_mean if facilitating, and nothing is pseudo-linear
    %need to save the latency and decay time of the synapse types including pseudo-linear
    if contains(Synapse_type,'depressing') && contains(Synapse_type,'Inhibitory')
        timeIEDF(1,1)=timeIEDF(1,1)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).d_mean; %Inhibitory is row 1, Depressing is column 1, Excitatory is row 2, Facilitating is column 2
        totalIEDF(1,1)=totalIEDF(1,1)+1;
        timedecay(1,1)=timedecay(1,1)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).decay_mean;
        timelatency(1,1)=timelatency(1,1)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).latency_mean;
        Stype='ID';
    elseif contains(Synapse_type,'facilitating') && contains(Synapse_type,'Inhibitory')
        timeIEDF(1,2)=timeIEDF(1,2)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).f_mean;
        totalIEDF(1,2)=totalIEDF(1,2)+1;
        timedecay(1,2)=timedecay(1,2)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).decay_mean;
        timelatency(1,2)=timelatency(1,2)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).latency_mean;
        Stype='IF';
    elseif contains(Synapse_type,'depressing') && contains(Synapse_type,'Excitatory')
        timeIEDF(2,1)=timeIEDF(2,1)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).d_mean;
        totalIEDF(2,1)=totalIEDF(2,1)+1;
        timedecay(2,1)=timedecay(2,1)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).decay_mean;
        timelatency(2,1)=timelatency(2,1)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).latency_mean;
        Stype='ED';
    elseif contains(Synapse_type,'facilitating') && contains(Synapse_type,'Excitatory')
        timeIEDF(2,2)=timeIEDF(2,2)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).f_mean;
        totalIEDF(2,2)=totalIEDF(2,2)+1;
        timedecay(2,2)=timedecay(2,2)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).decay_mean;
        timelatency(2,2)=timelatency(2,2)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).latency_mean;
        Stype='EF';
    elseif contains(Synapse_type,'pseudo-linear') && contains(Synapse_type,'Inhibitory')
        timedecay(1,3)=timedecay(2,3)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).decay_mean;
        timelatency(1,3)=timelatency(2,3)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).latency_mean;
        totalIEDF(1,3)=totalIEDF(1,3)+1;
        Stype='IP';
    elseif contains(Synapse_type,'pseudo-linear') && contains(Synapse_type,'Excitatory')
        timedecay(2,3)=timedecay(2,3)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).decay_mean;
        timelatency(2,3)=timelatency(2,3)+data_anatomy_L.(names_anatomy_L{neuronconnectiontypiteration}).latency_mean;
        totalIEDF(2,3)=totalIEDF(2,3)+1;
        Stype='EP';
    end
    
    
    %determin number of possible connections
    numberpresynaptic=neurontypenumber(strcmp(neuron_types,presynaptic_neuron));
    numberpostsynaptic=neurontypenumber(strcmp(neuron_types,postsynaptic_neuron));
    numberpossibleconnections=numberpresynaptic*numberpostsynaptic;
    
    %number probable connections
    numberprobableconnections=round(numberpossibleconnections*connection_probability/100);
    
    
    %find locations presynaptic neurons in neuron_list
    %find locations postsynaptic neurons in neuron_list
    xpre=find(cellfun(@(x) strcmp(x, presynaptic_neuron), neuron_list));
    xpost=find(cellfun(@(x) strcmp(x, postsynaptic_neuron), neuron_list));
    selectrandompreneuron = randi(numberpresynaptic, [numberprobableconnections, 1]);
    selectrandompostneuron = randi(numberpostsynaptic, [numberprobableconnections, 1]);
    %for each selected pre-synaptic neuron, connect to a random post-synaptic neuron, adding this connection to the exising ones if they exist and use a vectorised method to do this
    %use linear neuron_list for this
    % Create the logical index arrays
    pre_idx = xpre(selectrandompreneuron);
    post_idx = xpost(selectrandompostneuron);
    grouped_updates = accumarray(pre_idx(:), post_idx(:), [], @(x) {x});
    grouped_updates=[grouped_updates; cell(length(connection_array)-length(grouped_updates),1)];
    % Use cellfun to update the connection array
    connection_array = cellfun(@(x, y) [x;y], connection_array, grouped_updates, 'UniformOutput', false);
    %need to save Synapse_type in array of the same structure as connection_array
    
    newsynapsetypes = repmat({Stype}, [numberprobableconnections, 1]);
    % Use accumarray to group new values and synapse types by idstoupdate
    grouped_types = accumarray(pre_idx(:), (1:length(newsynapsetypes))', [], @(x) {newsynapsetypes(x)});
    grouped_types=[grouped_types; cell(length(Synapse_types_array)-length(grouped_types),1)];
    Synapse_types_array = cellfun(@(x, y) [x; y], Synapse_types_array, grouped_types, 'UniformOutput', false);
    countsynapsetype.(Stype)=countsynapsetype.(Stype)+length(newsynapsetypes);
    
end

timelatency=round(timelatency./totalIEDF);
timelatency(timelatency==0)=1;%can't have 0ms delay -rounded down
timedecay=round(timedecay./totalIEDF);
timeIEDF=timeIEDF./totalIEDF(1:2,1:2);
%now that we knokw the total number of neurons, they need to be distributed evenly within a cube
% Distribute the neurons in the cubic array
% Determine the size of the cubic array
cube_size = ceil(total_neurons^(1/3));

% Create the cubic array
cubic_array = nan(cube_size, cube_size, cube_size); %distribute connection_array into cube

index = 1;
for x = 1:cube_size
    for y = 1:cube_size
        for z = 1:cube_size
            cubic_array(x, y, z) = index; %index to connection_array
            index = index + 1;
        end
    end
end

%find centre of cube
centrecube = ceil(cube_size/2);

stimchn1_current=10;%uA
pulses=1;
timebetweenpulse=1000/300;%ms

%1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
minimum_current=1;%uA if electrode is touching the axon
%stimulation current levels and determine the activated area
%also number of pulses and whether the connections are fascillitating or depressing with each subsequent pulse

%1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
Neuron_densities=[87533.33, 177300, 107700]./10^9;%microcircuit
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
%homogenise layer properties and then just use activated area from the average of the layer properties

volume_whole_layer=sum(diff([(4/3).* pi.*(radius_probability1.^3) 0]).*-1.*percentageactivated);
elect1_primaryactivated=round(volume_whole_layer*mean(Neuron_densities));

%the nubmer of neurons in elect1_primaryactivated will be activated by the primary stimulation which will be at the centre of the cube

%number of neurons from the centre in each direction
neurons_distributed_depth = round((elect1_primaryactivated^(1/3))/2);

%find the coordinates of the neurons in the cube activated
neuronindex=cubic_array(centrecube-neurons_distributed_depth+1:centrecube+neurons_distributed_depth,centrecube-neurons_distributed_depth+1:centrecube+neurons_distributed_depth,centrecube-neurons_distributed_depth+1:centrecube+neurons_distributed_depth);
neuronindex=neuronindex(:);
%remove random neurons from the outer ring to make the number of neurons activated equal to elect1_primaryactivated
numneuronsremove=length(neuronindex)-elect1_primaryactivated;
%neuronindex needs those removed from the outer ring
indextoremove=randsample(1:length(neuronindex),[numneuronsremove]);
neuronindex(indextoremove)=[];

%%
% remove any duplicate connections, need to keep the same order of connections
for i = 1:length(connection_array)
    [connection_array{i},orderar] = unique(connection_array{i}, 'stable');
    %fix Synapse_types_array too
    Synapse_types_array{i} = Synapse_types_array{i}(orderar);
end

% Define the anonymous function to subtract 1 from elements greater than 1
subtractFunc = @(x) x - (x > 1);
pre_fd = connection_array; % this will hold the time for facilitation or depression in ms
%set all values to 0 using logical indexing
pre_fd = cellfun(@(x) zeros(size(x)), pre_fd, 'UniformOutput', false);
maxtime=100;%ms
%need to hold furture action potentials due to latency for max(timelatency) ms
maxfuturefiring=repmat({zeros(1,round(max(timelatency(:))))},[length(connection_array),1]);%this will hold the future action potentials due to latency
activeneurons=zeros(length(connection_array),maxtime); %this will be updated each iteration to show which neurons are activated
refactorytimeNeuron=randi([1,3],[size(activeneurons,1),1]);
lambdap = 10/1000; % Mean of the Poisson distribution
n = size(activeneurons,1); % Number of random numbers to generate
randomNumbers = poissrnd(lambdap, n, maxtime);
randomNumbers(randomNumbers>1)=1;
activeneurons=randomNumbers;
%each action potential lasts 1ms, need to iterate through 100 ms time to see what happens in the network
for timeiterate=1:maxtime-max(timelatency,[],'all')
    %need to work out if this is a stimulus time and add those neurons into the activated neurons    
    % absolute refactory period of 2ms means any stimulus that was positive in the last two ms cannot be positive again
    activeneurons(:,timeiterate)%%%%%% need to put in a poissson distribution of random neurons here for baseline
    if timeiterate>3
        %activeneurons(activeneurons(:,timeiterate-2)>=1,timeiterate)=0;
        for neuronit=1:length(refactorytimeNeuron)
            if any(activeneurons(neuronit,timeiterate-refactorytimeNeuron(neuronit):timeiterate-1)>=1) && activeneurons(neuronit,timeiterate)>=1
                activeneurons(neuronit,timeiterate)=0;
            end
        end
    end
    if any(timeiterate==[1 round([1:pulses-1].*timebetweenpulse)]) && timeiterate<=timebetweenpulse*pulses
        activeneurons(neuronindex,timeiterate)=activeneurons(neuronindex,timeiterate)+1;
        activeneurons(neuronindex(activeneurons(neuronindex,timeiterate)<=0),timeiterate)=1;
        %activated neurons will fire immediately, the connections listed in connection_array will become active giving +1 for excitatory and -1 for inhibitory on the secondary neurons
        %latency of the reaction will be added to the maxfuturefiring array
    end

    % step 1: work out the types of connections that are active
    active_connections = connection_array(activeneurons(:,timeiterate) >= 1);
    Synapse_types_connections = Synapse_types_array(activeneurons(:,timeiterate) >= 1);
    index_time=find(activeneurons(:,timeiterate) >= 1);
    % step 2: based on type of connection, add a +1 or -1 to the iteration of maxfuturefiring specified by the latency in timelatency
    for activeneuronnum = 1:length(active_connections)
        neuronIndividual=active_connections{activeneuronnum};
        for individual_connection = 1:length(neuronIndividual)
            connection_type = Synapse_types_connections{activeneuronnum}{individual_connection};
            if contains(connection_type, 'E') && contains(connection_type, 'F')
                pre_fd{index_time(activeneuronnum)}(individual_connection) = pre_fd{index_time(activeneuronnum)}(individual_connection)+timeIEDF(2,2);
                latency = timelatency(2,2);
                timetof=timeiterate + latency + timedecay(2,2);
                lambda=0.02*300/timedecay(2,2);
                activeneuronsadd=round(pre_fd{index_time(activeneuronnum)}(individual_connection)/timeIEDF(2,2)):-(round(pre_fd{index_time(activeneuronnum)}(individual_connection)/timeIEDF(2,2))/length(timeiterate + latency:timetof)):0;
                timev=1:length(activeneuronsadd);
                activeneuronsadd=max(activeneuronsadd).*exp(-lambda.*timev)./max(exp(-lambda.*timev));
                timetof(timetof>maxtime)=maxtime;
                activeneuronsadd(length(timeiterate + latency:timetof)+1:end)=[];
                activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) =  activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof)+activeneuronsadd;
            elseif contains(connection_type, 'E') && contains(connection_type, 'D')
                latency = timelatency(2,1);
                timetof=timeiterate + latency + timedecay(2,1);
                lambda=0.02*300/timedecay(2,1);
                activeneuronsadd=1:-1/length(timeiterate + latency:timetof):0;
                timev=1:length(activeneuronsadd);
                activeneuronsadd=exp(-lambda.*timev)./max(exp(-lambda.*timev));
                timetof(timetof>maxtime)=maxtime;
                activeneuronsadd(length(timeiterate + latency:timetof)+1:end)=[];
                pre_fd{index_time(activeneuronnum)}(individual_connection) = pre_fd{index_time(activeneuronnum)}(individual_connection)+timeIEDF(2,1);
                activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) =  activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) + (1/((pre_fd{index_time(activeneuronnum)}(individual_connection))/timeIEDF(2,1))).*activeneuronsadd;
            elseif contains(connection_type, 'E') && contains(connection_type, 'P')
                latency = timelatency(2,3);
                timetof=timeiterate + latency + timedecay(2,3);
                lambda=0.02*300/timedecay(2,3);
                activeneuronsadd=1:-1/length(timeiterate + latency:timetof):0;
                 timev=1:length(activeneuronsadd);
                activeneuronsadd=exp(-lambda.*timev)./max(exp(-lambda.*timev));
                timetof(timetof>maxtime)=maxtime;
                activeneuronsadd(length(timeiterate + latency:timetof)+1:end)=[];
                activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) =  activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) + activeneuronsadd;
            elseif contains(connection_type, 'I') && contains(connection_type, 'F')
                latency = timelatency(1,2);
                timetof=timeiterate + latency + timedecay(1,2);
                lambda=0.02*300/timedecay(1,2);
                pre_fd{index_time(activeneuronnum)}(individual_connection) = pre_fd{index_time(activeneuronnum)}(individual_connection)+timeIEDF(1,2);
                activeneuronsadd=-round(pre_fd{index_time(activeneuronnum)}(individual_connection)/timeIEDF(1,2)):(round(pre_fd{index_time(activeneuronnum)}(individual_connection)/timeIEDF(1,2))/length(timeiterate + latency:timetof)):0;
                timev=1:length(activeneuronsadd);
                activeneuronsadd=min(activeneuronsadd).*exp(-lambda.*timev)./max(exp(-lambda.*timev));
                timetof(timetof>maxtime)=maxtime;
                activeneuronsadd(length(timeiterate + latency:timetof)+1:end)=[];
                activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) =  activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) + activeneuronsadd;
            elseif contains(connection_type, 'I') && contains(connection_type, 'D')
                latency = timelatency(1,1);
                timetof=timeiterate + latency + timedecay(1,1);
                lambda=0.02*300/timedecay(1,1);
                activeneuronsadd=-1:1/length(timeiterate + latency:timetof):0;
                timev=1:length(activeneuronsadd);
                activeneuronsadd=min(activeneuronsadd).*exp(-lambda.*timev)./max(exp(-lambda.*timev));
                timetof(timetof>maxtime)=maxtime;
                activeneuronsadd(length(timeiterate + latency:timetof)+1:end)=[];
                pre_fd{index_time(activeneuronnum)}(individual_connection) = pre_fd{index_time(activeneuronnum)}(individual_connection)+timeIEDF(1,1);
                activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) =  activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) + (1/((pre_fd{index_time(activeneuronnum)}(individual_connection))/timeIEDF(1,1))).*activeneuronsadd;
            elseif contains(connection_type, 'I') && contains(connection_type, 'P')
                latency = timelatency(1,3);
                timetof=timeiterate + latency + timedecay(1,3);
                lambda=0.02*300/timedecay(1,3);
                activeneuronsadd=-1:1/length(timeiterate + latency:timetof):0;
                timev=1:length(activeneuronsadd);
                activeneuronsadd=min(activeneuronsadd).*exp(-lambda.*timev)./max(exp(-lambda.*timev));
                timetof(timetof>maxtime)=maxtime;
                activeneuronsadd(length(timeiterate + latency:timetof)+1:end)=[];
                activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) =  activeneurons(active_connections{activeneuronnum}(individual_connection),timeiterate + latency:timetof) + activeneuronsadd;
            end
        end
    end
     pre_fd=cellfun(subtractFunc, pre_fd, 'UniformOutput', false);
end
%%
%smooth response in time for each neuron
activeneurons_filtered=filtfilt(1/3*ones(3,1),1,activeneurons')';
%add nan values to bottom of activeneurons to make it the same size as the cube
activeneurons_filtered=[activeneurons_filtered; nan(numel(cubic_array)-size(activeneurons_filtered,1),size(activeneurons_filtered,2))];
%%
%create a figure that iterates through time in the cube and shows the activated neurons
filterSize = 3; % Size of the filter (e.g., 5x5)
sigma = 1.0;   % Standard deviation of the Gaussian function
electrodeposition=[centrecube,centrecube,centrecube];
firing_rate=zeros(maxtime,1);
% Create the 3D Gaussian filter
gaussianFilter = fspecial3('gaussian', filterSize, sigma);
for timeit=1:maxtime
    neuronstoplottime=activeneurons_filtered(:,timeit);
    neuronstoplottime=neuronstoplottime(cubic_array);
    % Apply the Gaussian filter to the data using imfilter
%filteredData = imfilter(neuronstoplottime, gaussianFilter, 'symmetric');
filteredData=neuronstoplottime;
filteredData(filteredData==0)=nan;

    [x, y, z] = ind2sub(size(filteredData), find(~isnan(filteredData)));
    values = filteredData(~isnan(filteredData));
    figure(1);
    %imagesc(values)
    scatter3(x, y, z, 5, values);  % Use '100' to scale the marker size and 'values' to color the markers
    %view(2)
    colormap(parula); % Use a linear colormap
    colorbar;
    title(['Time: ' num2str(timeit) 'ms']);
    clim([-1 1]);
    %zlim([15 30])

    figure(2)
    %values=squeeze(mean(neuronstoplottime,1,'omitnan'));
    values=squeeze(neuronstoplottime(1:16,1:16,16));
    imagesc(values)
       colormap(parula); % Use a linear colormap
    colorbar;
    title(['Time: ' num2str(timeit) 'ms']);
    clim([-1 1]);
 

    %need to model an electrode that can record a 5x5x5 cube of neurons from the position specified in electrodeposition
    excitedandinhibneurons=(neuronstoplottime<=-1).*-1;
    excitedandinhibneurons=[excitedandinhibneurons+(neuronstoplottime>=1)];
    firing_rate(timeit)=sum(excitedandinhibneurons(electrodeposition(1)-2:electrodeposition(1)+2,electrodeposition(2)-2:electrodeposition(2)+2,electrodeposition(3)-2:electrodeposition(3)+2),"all");
   firing_rate(timeit)=firing_rate(timeit)*1000/(5*5*5);
    %pause(0.7);
   
end
    figure(3)
    scatter(1:maxtime,firing_rate);
    hold on
    xlabel('Time (ms)');
    ylabel('Firing rate (Hz)');
    %fit line to data
    p = polyfit(1:maxtime-3,firing_rate(1:end-3)',2);
    yfit = polyval(p,1:maxtime-3);
    plot(1:maxtime-3,yfit,'r-');
    %%
    for i=1:27000
    interArrivalTimes = exprnd(1/lambda, 1, 1000);
% Calculate the arrival times by taking the cumulative sum
arrivalTimes = cumsum(interArrivalTimes);
% Keep only the arrival times within the specified time window
arrivalTimes = arrivalTimes(arrivalTimes <= T);
if ~isempty(arrivalTimes)
figure;
stem(arrivalTimes, ones(size(arrivalTimes)), 'filled');
xlabel('Time');
ylabel('Events');
title('Poisson Process');
end
    end
%%
%bursting activity due to extracellular stim https://www.sciencedirect.com/science/article/pii/S1935861X09000424#app1
burstprob=1./8;%assume cell might fire twice in 8ms window and record for 6ms(2-8)
firingrate=neurons_distributed_depth.*burstprob.*1000;

% %%
%
%
% % layerstart=[0 368 368+154];
% %middlelayerpoint=((layers_depthdefinition-layerstart)/2)+layerstart;
% %layer_thickness=(layers_depthdefinition-layerstart);
% %note that current dissipates and E-fields never
% %interact; 34um max dist until current dissipates
% %layers_depthdefinition=[368 154+368 718+154+368];
%
% radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
% filtmov=1/5*ones(5,1);
% firingrate=filtfilt(filtmov,1,firingrate);
% %based on the primary activation, how many secondary neurons are activated consider the probability of connection, and the weight of
% %the connection
% volumelayers= (layers_depthdefinition-layerstart).* (pi.*(radius_circuit).^2);
% area_activated=((volume_whole_layer)./sum(volumelayers));%percentage volume activated;;
% sum(Neuron_densities.*volumelayers)
% connectionsperneuron=totalconnections./sum(Neuron_densities.*volumelayers);
% Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);
% connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
% Secondary_neurons_activated=ceil(connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
% baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);
% tertiary_connections=sum((Secondary_neurons_activated./(Neuron_densities'.*volumelayers'))'.*Numconnections_3groups_array,2);
% tertiary_neurons_activated=ceil(tertiary_connections.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
% SecondaryTertiary_neurons_activated=Secondary_neurons_activated+tertiary_neurons_activated;
%
% weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);
% %     %should convergent and divergent be considered (are there enough connections?)
% %     if numberprobableconnections/numberpresynaptic<number_of_divergent_neuron_mean && numberprobableconnections/numberpostsynaptic<number_of_convergent_neuron_mean
% %         ignoreprepost=[1,1];
% %     elseif numberprobableconnections/numberpresynaptic<number_of_divergent_neuron_mean
% %         ignoreprepost=[1,0];
% %     elseif numberprobableconnections/numberpostsynaptic<number_of_convergent_neuron_mean
% %         ignoreprepost=[0,1];
% %     else
% %         ignoreprepost=[0,0];
% %     end
%
%     %check this is higher than the number of divergent neurons and convergent neurons
% %     if numberprobableconnections<number_of_divergent_neuron_mean || numberprobableconnections<number_of_convergent_neuron_mean
% %         error('Number of probable connections is less than the number of divergent or convergent neurons')
% %     end
%
%
% %         if ignoreprepost(1)==1 && ignoreprepost(2)==1
% %                     elseif ignoreprepost(1)==1 && ignoreprepost(2)==0
% %             %connect to random neuron
% %             if flagpost==0
% %                 selectrandompreneuron=randi(numberpresynaptic);
% %                 selectrandompostneuron=randi(numberpostsynaptic);
% %                 connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}=[connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}; [xpost(selectrandompostneuron),ypost(selectrandompostneuron),zpost(selectrandompostneuron)]];
% %                 flagpost=round(postsynaptic_neuron)-1;
% %             elseif flagpost>0
% %                 selectrandompreneuron=randi(numberpresynaptic);
% %                 connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}=[connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}; [xpost(selectrandompostneuron),ypost(selectrandompostneuron),zpost(selectrandompostneuron)]];
% %                 flagpost=flagpost-1;
% %             end
% %         elseif ignoreprepost(1)==0 && ignoreprepost(2)==1
% %             if flagpre==0
% %                 selectrandompreneuron=randi(numberpresynaptic);
% %                 selectrandompostneuron=randi(numberpostsynaptic);
% %                 connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}=[connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}; [xpost(selectrandompostneuron),ypost(selectrandompostneuron),zpost(selectrandompostneuron)]];
% %                 flagpre=round(presynaptic_neuron)-1;
% %             elseif flagpre>0
% %                 selectrandompostneuron=randi(numberpostsynaptic);
% %                 connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}=[connectioncubicarray{xpre(selectrandompreneuron),ypre(selectrandompreneuron),zpre(selectrandompreneuron)}; [xpost(selectrandompostneuron),ypost(selectrandompostneuron),zpost(selectrandompostneuron)]];
% %                 flagpre=flagpre-1;
% %             end
% %
% %         end
% %setup data structures
% for sepdist=5:2:9
%     sepcheck=['sep' num2str(sepdist)];
%     for current=1:length(AMP)
%         currcheck=['C' num2str(AMP(current))];
%         for trial=1:5
%             trialcheck=['T' num2str(trial)];
%             for shanksep=0:3
%                 shanksepcheck=['D' num2str(shanksep)];
%                 Csplit_simulation.(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[];
%             end
%         end
%     end
% end
%
% % need to just take an average of all connection data related to layers, average accross s, g and i
% it=0;
% Numconnections=[];
% L_combinations=[];
% totalconnect_L=[];
% Numconnections_3groups=[];
% totalconnect=0;
% for from_L=1:5 %creats structures for storing connection data
%     if from_L<4
%         from_L_group='s';
%     elseif from_L==4
%         from_L_group='g';
%     elseif from_L>4
%         from_L_group='i';
%     end
%     for to_L=1:5
%         %split into original layers
%         it=it+1;
%         L_combinations{it}=['L' strcat(from_L_str{from_L}, to_L_str{to_L})];
%         Numconnections.(L_combinations{it})=0;
%
%         %split into 3 groups
%         if to_L<4
%             to_L_group='s';
%         elseif to_L==4
%             to_L_group='g';
%         elseif to_L>4
%             to_L_group='i';
%         end
%         L_combinationsgroup{it}=strcat(from_L_group, to_L_group);
%
%         Numconnections_3groups.(L_combinationsgroup{it})=0;
%         count_inhib_excite.(L_combinationsgroup{it})=[0 0];
%     end
%     totalconnect_L.(L_combinations{it}(1:2))=0;
% end
% totalconnect_L.all=0;
%
% Numconnections_3groups_array=zeros(3,3);
% inhib_excite_arrayFR=zeros(3,3);
% inhib_excite_arrayCOUNT=zeros(3,3);
% Numconnections_inhib_excite=[0 0];
% totalconnections=[0 0];
% ratioEI_prepostconnection=[0 0];
%
%
% for it_fn=1:length(fn_L) %calculates number of connections while looping through possible neuron combinations
%     underscore_L = strfind(fn_L{it_fn},'___');
%     Lconnect=(fn_L{it_fn}([2 underscore_L+4]));
%     L_check=['L' Lconnect];
%
%     Numconnections.(L_check)=Numconnections.(L_check)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
%
%     totalconnect_L.(L_check(1:2))=totalconnect_L.(L_check(1:2))+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
%     totalconnect_L.all=totalconnect_L.all+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
%
%     %3 group split
%     if str2double(fn_L{it_fn}(underscore_L+4))<4
%         to_L_group='s';
%         to_L_groupnum=1;
%     elseif str2double(fn_L{it_fn}(underscore_L+4))==4
%         to_L_group='g';
%         to_L_groupnum=2;
%     elseif str2double(fn_L{it_fn}(underscore_L+4))>4
%         to_L_group='i';
%         to_L_groupnum=3;
%     end
%     if str2double(fn_L{it_fn}(2))<4
%         from_L_group='s';
%         from_L_groupnum=1;
%     elseif str2double(fn_L{it_fn}(2))==4
%         from_L_group='g';
%         from_L_groupnum=2;
%     elseif str2double(fn_L{it_fn}(2))>4
%         from_L_group='i';
%         from_L_groupnum=3;
%     end
%
%     synapse_type_connection=data_anatomy_L.(names_anatomy_L{it_fn}).synapse_type;
%
%     inhib_1=strfind(synapse_type_connection,'Inhibitory')+1;
%     if isempty(inhib_1)
%         inhib_1=1;
%     end
%     L_3layercheck=strcat(from_L_group, to_L_group);
%     totalconnections(inhib_1)=totalconnections(inhib_1)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
%     ratioEI_prepostconnection(inhib_1)=ratioEI_prepostconnection(inhib_1)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection)*(data_L.(fn_L{it_fn}).number_of_divergent_neuron_mean/data_L.(fn_L{it_fn}).number_of_convergent_neuron_mean);
%     Numconnections_3groups.(L_3layercheck)=Numconnections_3groups.(L_3layercheck)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection)-data_L.(fn_L{it_fn}).number_of_convergent_neuron_mean+data_L.(fn_L{it_fn}).number_of_divergent_neuron_mean;
%     Numconnections_3groups_array(to_L_groupnum,from_L_groupnum)=Numconnections_3groups_array(to_L_groupnum,from_L_groupnum)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection)-data_L.(fn_L{it_fn}).number_of_convergent_neuron_mean+data_L.(fn_L{it_fn}).number_of_divergent_neuron_mean;
%     count_inhib_excite.(L_3layercheck)(inhib_1)=count_inhib_excite.(L_3layercheck)(inhib_1)+1;
%     inhib_excite_arrayFR(to_L_groupnum,from_L_groupnum)=inhib_excite_arrayFR(to_L_groupnum,from_L_groupnum)+firing_rates_EI(inhib_1);
%     inhib_excite_arrayCOUNT(to_L_groupnum,from_L_groupnum)=inhib_excite_arrayCOUNT(to_L_groupnum,from_L_groupnum)+1;
%     Numconnections_inhib_excite(inhib_1)=Numconnections_inhib_excite(inhib_1)+(data_L.(fn_L{it_fn}).total_synapse_count./data_L.(fn_L{it_fn}).mean_number_of_synapse_per_connection);
% end
% layers_depthdefinition=[353+149+165 190+353+149+165 525+700+190+353+149+165];%thickness of microcircuit
% layers_thickness=[368 154 718];
% percent_connection=layers_thickness./(layers_depthdefinition-[0 667 857]);%microcircuit is thicker than vis cortex
% Numconnections_3groups_array=Numconnections_3groups_array.*percent_connection;
% prepost=ratioEI_prepostconnection./totalconnections;
%
% layers_depthdefinition=[368 522 1240];
% %make homogenous
%
%
% %1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
% minimum_current=1;%uA if electrode is touching the axon
% %note that current dissipates and E-fields never
% %interact; 34um max dist until current dissipates
% stimchn1_current=AMP(current);
% %1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
%
% layers_depthdefinition=[368 154+368 718+154+368];
% Neuron_densities=[87533.33, 177300, 107700]./10^9;%microcircuit
% layerstart=[0 368 368+154];
% middlelayerpoint=((layers_depthdefinition-layerstart)/2)+layerstart;
% layer_thickness=(layers_depthdefinition-layerstart);
% %note that current dissipates and E-fields never
% %interact; 34um max dist until current dissipates
% percentageactivated=0.01:0.01:1;%0 to 100 for use later
% ALL_k_const=interp1([1 50 100],[2100 8850 27500],1:100,'spline');% 0 25 50 75 100% of neurons activated
% k_const=ALL_k_const(1);%constant ua/mm^2, const at 0-1%neurons activated will multiply primary activation by ALL-k-const probability function
% if stimchn1_current>=minimum_current
%     radius_activated1=((stimchn1_current-minimum_current)/k_const)^0.5;
%     radius_activated1=radius_activated1*10.^3;
%     radius_probability1=(((stimchn1_current-minimum_current)./ALL_k_const).^0.5).*10^3;
% else
%     radius_activated1=0;
%     radius_probability1=0;
% end
% %homogenise layer properties and then just use activated area from the average of the layer properties
% radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
% volume_whole_layer=sum(diff([(4/3).* pi.*(radius_probability1.^3) 0]).*-1.*percentageactivated);
% elect1_primaryactivated=volume_whole_layer*mean(Neuron_densities);
% %based on the primary activation, how many secondary neurons are activated consider the probability of connection, and the weight of
% %the connection
% volumelayers= (layers_depthdefinition-layerstart).* (pi.*(radius_circuit).^2);
% area_activated=((volume_whole_layer)./sum(volumelayers));%percentage volume activated;;
% sum(Neuron_densities.*volumelayers)
% connectionsperneuron=totalconnections./sum(Neuron_densities.*volumelayers);
% Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);
% connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
% Secondary_neurons_activated=ceil(connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
% baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);
% tertiary_connections=sum((Secondary_neurons_activated./(Neuron_densities'.*volumelayers'))'.*Numconnections_3groups_array,2);
% tertiary_neurons_activated=ceil(tertiary_connections.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
% SecondaryTertiary_neurons_activated=Secondary_neurons_activated+tertiary_neurons_activated;
%
% weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);
%
% %bursting activity due to extracellular stim https://www.sciencedirect.com/science/article/pii/S1935861X09000424#app1
%
% burstprob=1./8;%assume cell might fire twice in 8ms window and record for 6ms(2-8)
% firingrate=neurons_distributed_depth.*burstprob.*1000;
% filtmov=1/5*ones(5,1);
% firingrate=filtfilt(filtmov,1,firingrate);
% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for ratN=[6 9 13 21:23]%[6:13 21:23] %loop through animals
%     %load data
%     if ratN<10
%         Ratnum=['Rat_00' num2str(ratN)];
%     elseif ratN>=10
%         Ratnum=['Rat_0' num2str(ratN)];
%     end
%
%     cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
%     D_data=dir;
%     if ~any(strcmp({D_data.name}, 'ElectLayerClass.mat')) %check that we were able to correctly identify layers
%         continue
%     end
%     load('ElectLayerClass.mat','ElectLayerClass')
%     for shankiterate=0:3
%         checkshank=['s' num2str(shankorder(shankiterate+1))];
%         if any(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==4)
%             firste1=find(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==3,1,'first');
%             Depthestimate.(checkshank)(1:firste1-1)=layers_depthdefinition(end)+50*(firste1-1):-50:layers_depthdefinition(end)+50;
%             Depthestimate.(checkshank)(firste1:end)=layers_depthdefinition(end):-50:layers_depthdefinition(end)-50*(16-firste1);
%         else
%             firste1=find(ElectLayerClass(1+(shankiterate*16):16+(shankiterate*16))==2,1,'first');
%             Depthestimate.(checkshank)(1:firste1-1)=layers_depthdefinition(2)+50*(firste1-1):-50:layers_depthdefinition(2)+50;
%             Depthestimate.(checkshank)(firste1:end)=layers_depthdefinition(2):-50:layers_depthdefinition(2)-50*(16-firste1);
%         end
%     end
%     for k = 3:length(D_data) % loop through the stimulation pairs. Avoid using the first ones
%         currD = D_data(k).name; % Get the current subdirectory name
%         try
%             cd([D_data(k).folder filesep currD])
%             load('Averagetrialresponse.mat','avgnospT') % load sorted neural activity
%         catch
%             continue;
%         end
%         loadStimChn;
%         Stimchnpos_layers=ElectLayerClass(stimChn);
%         if all(Stimchnpos_layers==1)
%             layers_stim='ss';
%         elseif all(Stimchnpos_layers==2)
%             layers_stim='gg';
%         elseif all(Stimchnpos_layers==3)
%             layers_stim='ii';
%         elseif any(Stimchnpos_layers==1) && any(Stimchnpos_layers==2)
%             layers_stim='sg';
%         elseif any(Stimchnpos_layers==2) && any(Stimchnpos_layers==3)
%             layers_stim='gi';
%         elseif any(Stimchnpos_layers==1) && any(Stimchnpos_layers==3)
%             layers_stim='si';
%         elseif any(Stimchnpos_layers==4)
%             layers_stim='WM';
%             warning('Stim electrodes are in white matter!!!!!')
%         else
%             error('Layers not saved properly')
%         end
%         for sepdist=5:2:9
%             sepcheck=['sep' num2str(sepdist)];
%             if (stimChn(1)-stimChn(2))~=(-1*sepdist)-1 % check stimchn sepdist between these electrodes
%                 continue
%             end
%             % spike rate and centroid calculation
%             if stimChn(1)<17 %determines the shank with the stimulating electrodes
%                 shank=1;
%             elseif stimChn(1)<33 && stimChn(1)>16
%                 shank=4;
%                 stimChn=stimChn-16;
%             elseif stimChn(1)<49 && stimChn(1)>32
%                 shank=2;
%                 stimChn=stimChn-32;
%             else
%                 shank=3;
%                 stimChn=stimChn-48;
%             end
%             for current = 1:length(AMP)
%                 currcheck=['C' num2str(AMP(current))];
%
%                 for trial=1:5
%                     trialcheck=['T' num2str(trial)];
%                     for shanksep=0:3
%                         shanksepcheck=['D' num2str(shanksep)];
%                         %%%%%%%%%calculate theoretical weightings
%
%
%
%                         %model activity
%                         stimshankcheck=['s' num2str(shank)];
%                         stimchn1depth=Depthestimate.(stimshankcheck)(stimChn(1));%um
%                         stimchn2depth=Depthestimate.(stimshankcheck)(stimChn(2));%um
%
%
%                         %1. area activated % https://reader.elsevier.com/reader/sd/pii/0165027096000659?token=020F942506A1B062F58C5B4AA93E50E7C1A13CAF041088D85EC6EFDE6E6F19188ADAEF399D35974672C0963CACA62196&originRegion=us-east-1&originCreation=20211022041758
%
%                         layers_depthdefinition=[368 154+368 718+154+368];
%                         Neuron_densities=[87533.33, 177300, 107700]./10^9;%microcircuit
%                         layerstart=[0 368 368+154];
%                         middlelayerpoint=((layers_depthdefinition-layerstart)/2)+layerstart;
%                         layer_thickness=(layers_depthdefinition-layerstart);
%                         %note that current dissipates and E-fields never
%                         %interact; 34um max dist until current dissipates
%                         percentageactivated=0.01:0.01:1;%0 to 100 for use later
%                         ALL_k_const=interp1([1 50 100],[2100 8850 27500],1:100,'spline');% 0 25 50 75 100% of neurons activated
%                         k_const=ALL_k_const(1);%constant ua/mm^2, const at 0-1%neurons activated will multiply primary activation by ALL-k-const probability function
%                         if stimchn1_current>=minimum_current
%                             radius_activated1=((stimchn1_current-minimum_current)/k_const)^0.5;
%                             radius_activated1=radius_activated1*10.^3;
%                             radius_probability1=(((stimchn1_current-minimum_current)./ALL_k_const).^0.5).*10^3;
%                         else
%                             radius_activated1=0;
%                             radius_probability1=0;
%                         end
%                         if stimchn2_current>=minimum_current
%                             radius_activated2=((stimchn2_current-minimum_current)/k_const)^0.5;
%                             radius_activated2=radius_activated2*10.^3;
%                             radius_probability2=(((stimchn2_current-minimum_current)./ALL_k_const).^0.5).*10^3;
%                         else
%                             radius_activated2=0;
%                             radius_probability2=0;
%                         end
%                         %2. Is radius confined to one group of layers? -
%                         %use stim elect positions
%                         %need to define distance based on histology and
%                         %microdrive depth and LFP
%
%                         UpperLstimchn1=stimchn1depth-radius_activated1;
%                         LowerLstimchn1=stimchn1depth+radius_activated1;
%                         UpperLstimchn2=stimchn2depth-radius_activated2;
%                         LowerLstimchn2=stimchn2depth+radius_activated2;
%                         percentagelayers1u=[1-(layers_depthdefinition-UpperLstimchn1)./(layers_depthdefinition-layerstart) -1];
%                         percentagelayers1l=[1-(layers_depthdefinition-LowerLstimchn1)./(layers_depthdefinition-layerstart) -1];
%                         percentagelayers2u=[1-(layers_depthdefinition-UpperLstimchn2)./(layers_depthdefinition-layerstart) -1];
%                         percentagelayers2l=[1-(layers_depthdefinition-LowerLstimchn2)./(layers_depthdefinition-layerstart) -1];
%                         firstupper1=find(percentagelayers1u<1,1, 'first');
%                         firstlower1=find(percentagelayers1l<1,1, 'first');
%                         firstupper2=find(percentagelayers2u<1,1, 'first');
%                         firstlower2=find(percentagelayers2l<1,1, 'first');
%
%                         %3. what is the percentage of each activated layer in
%                         %the microcircuit column? i.e. neurons activated
%                         %only through stimulation - primary
%                         %0.29mm^3 volume in the circuit 2082um long
%                         radius_circuit=sqrt(0.29/(2.082*pi))*10^3;
%                         area_activated=[0 0 0];
%                         volumelayers= (layers_depthdefinition-layerstart).* (pi.*(radius_circuit).^2);
%                         volume_upper=[];
%                         volume_lower=[];
%                         %                             %primary connections activated
%                         %                             if firstupper1~=firstlower1 %cut by layer boundary
%                         %                                 heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
%                         %                                 if heightsegment>radius_activated1 %reverse cap calc
%                         %                                     volume_lower=((pi*(2*radius_activated1-heightsegment)^2)/3)*(3*radius_activated1-(2*radius_activated1-heightsegment));
%                         %                                     volume_upper=((4/3) * pi*(radius_activated1^3))-volume_lower;
%                         %                                 else % calc cap
%                         %                                     volume_upper=((pi*heightsegment^2)/3)*(3*radius_activated1-heightsegment);
%                         %                                     volume_lower=((4/3) * pi*(radius_activated1^3))-volume_upper;
%                         %                                 end
%                         %                                 elect1_primaryactivated=volume_upper*Neuron_densities(firstupper1);
%                         %                                 elect1_primaryactivated=elect1_primaryactivated+volume_lower*Neuron_densities(firstlower1);
%                         %                                 area_activated(firstupper1)=volume_upper/volumelayers(firstupper1);%percentage volume activated;
%                         %                                 area_activated(firstlower1)=volume_lower/volumelayers(firstlower1);%percentage volume activated;
%                         %                             else %whole volume in one layer group
%                         %                                 area_activated(firstupper1)=((4/3) * pi*(radius_activated1^3))/volumelayers(firstupper1);%percentage volume activated
%                         %                                 elect1_primaryactivated=((4/3) * pi*(radius_activated1^3))*Neuron_densities(firstupper1);
%                         %                             end
%                         if firstupper1~=4 || firstlower1~=4
%                             if firstupper1~=firstlower1 && firstlower1==4
%                                 heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
%                                 if heightsegment>radius_activated1 %reverse cap calc
%                                     for iterate_sphereprobabilitydist=1:length(radius_probability1)
%                                         heightsegment=layers_depthdefinition(firstupper1)-stimchn1depth+radius_probability1(iterate_sphereprobabilitydist);
%                                         if heightsegment<(1*radius_probability1(iterate_sphereprobabilitydist))
%                                             volume_lower(iterate_sphereprobabilitydist)=((pi*((2*radius_probability1(iterate_sphereprobabilitydist)-heightsegment)^2)/3)*(3*radius_probability1(iterate_sphereprobabilitydist)-(2*radius_probability1(iterate_sphereprobabilitydist)-heightsegment)));
%                                             volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3))-volume_lower(iterate_sphereprobabilitydist);
%                                         else
%                                             volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3));
%                                         end
%                                     end
%                                 else % calc cap
%                                     for iterate_sphereprobabilitydist=1:length(radius_probability1)
%                                         heightsegment=layers_depthdefinition(firstupper1)-stimchn1depth+radius_probability1(iterate_sphereprobabilitydist);
%                                         if heightsegment>0
%                                             volume_upper(iterate_sphereprobabilitydist)=((pi*heightsegment^2)/3)*(3*radius_probability1(iterate_sphereprobabilitydist)-heightsegment);
%                                             volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3))-volume_upper(iterate_sphereprobabilitydist);
%                                         else
%                                             volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3));
%                                         end
%                                     end
%                                 end
%                                 volume_upper_percentage=(diff([volume_upper 0]).*-1).*percentageactivated(1:length(volume_upper));%volume activated
%                                 elect1_primaryactivated=sum(volume_upper_percentage)*Neuron_densities(firstupper1);
%                                 area_activated(firstupper1)=area_activated(firstupper1)+(sum(volume_upper_percentage)/volumelayers(firstupper1));%percentage volume activated;;
%                             elseif firstupper1~=firstlower1 && firstlower1~=4%cut by layer boundary
%                                 heightsegment=layers_depthdefinition(firstupper1)-UpperLstimchn1;
%                                 if heightsegment>radius_activated1 %reverse cap calc
%                                     for iterate_sphereprobabilitydist=1:length(radius_probability1)
%                                         heightsegment=layers_depthdefinition(firstupper1)-stimchn1depth+radius_probability1(iterate_sphereprobabilitydist);
%                                         if heightsegment<(1*radius_probability1(iterate_sphereprobabilitydist))
%                                             volume_lower(iterate_sphereprobabilitydist)=((pi*((2*radius_probability1(iterate_sphereprobabilitydist)-heightsegment)^2)/3)*(3*radius_probability1(iterate_sphereprobabilitydist)-(2*radius_probability1(iterate_sphereprobabilitydist)-heightsegment)));
%                                             volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3))-volume_lower(iterate_sphereprobabilitydist);
%                                         else
%                                             volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3));
%                                         end
%                                     end
%                                 else % calc cap
%                                     for iterate_sphereprobabilitydist=1:length(radius_probability1)
%                                         heightsegment=layers_depthdefinition(firstupper1)-stimchn1depth+radius_probability1(iterate_sphereprobabilitydist);
%                                         if heightsegment>0
%                                             volume_upper(iterate_sphereprobabilitydist)=((pi*heightsegment^2)/3)*(3*radius_probability1(iterate_sphereprobabilitydist)-heightsegment);
%                                             volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3))-volume_upper(iterate_sphereprobabilitydist);
%                                         else
%                                             volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability1(iterate_sphereprobabilitydist)^3));
%                                         end
%                                     end
%                                 end
%                                 volume_upper_percentage=(diff([volume_upper 0]).*-1).*percentageactivated(1:length(volume_upper));%volume activated
%                                 volume_lower_percentage=(diff([volume_lower 0]).*-1).*percentageactivated(1:length(volume_lower));%volume activated
%                                 elect1_primaryactivated=sum(volume_upper_percentage)*Neuron_densities(firstupper1);
%                                 elect1_primaryactivated=elect1_primaryactivated+sum(volume_lower_percentage)*Neuron_densities(firstlower1);
%                                 area_activated(firstupper1)=area_activated(firstupper1)+(sum(volume_upper_percentage)/volumelayers(firstupper1));%percentage volume activated;;
%                                 area_activated(firstlower1)=area_activated(firstlower1)+(sum(volume_lower_percentage)/volumelayers(firstlower1));%percentage volume activated;;
%                             else %whole volume in one layer group
%                                 volume_whole_layer=sum(diff([(4/3).* pi.*(radius_probability1.^3) 0]).*-1.*percentageactivated);
%                                 area_activated(firstupper1)=area_activated(firstupper1)+((volume_whole_layer)/volumelayers(firstupper1));%percentage volume activated;;
%                                 elect1_primaryactivated=volume_whole_layer*Neuron_densities(firstupper1);
%                             end
%                         else
%                             elect1_primaryactivated=0;
%                         end
%
%                         volume_upper=[];
%                         volume_lower=[];
%                         %primary connections activated
%                         if firstupper2~=firstlower2 %cut by layer boundary
%                             heightsegment=layers_depthdefinition(firstupper2)-UpperLstimchn2;
%                             if heightsegment>radius_activated2 %reverse cap calc
%                                 for iterate_sphereprobabilitydist=1:length(radius_probability2)
%                                     heightsegment=layers_depthdefinition(firstupper2)-stimchn2depth+radius_probability2(iterate_sphereprobabilitydist);
%                                     if heightsegment<(2*radius_probability2(iterate_sphereprobabilitydist))
%                                         volume_lower(iterate_sphereprobabilitydist)=((pi*((2*radius_probability2(iterate_sphereprobabilitydist)-heightsegment)^2)/3)*(3*radius_probability2(iterate_sphereprobabilitydist)-(2*radius_probability2(iterate_sphereprobabilitydist)-heightsegment)));
%                                         volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3))-volume_lower(iterate_sphereprobabilitydist);
%                                     else
%                                         volume_upper(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3));
%                                     end
%                                 end
%                             else % calc cap
%                                 for iterate_sphereprobabilitydist=1:length(radius_probability2)
%                                     heightsegment=layers_depthdefinition(firstupper2)-stimchn2depth+radius_probability2(iterate_sphereprobabilitydist);
%                                     if heightsegment>0
%                                         volume_upper(iterate_sphereprobabilitydist)=((pi*heightsegment^2)/3)*(3*radius_probability2(iterate_sphereprobabilitydist)-heightsegment);
%                                         volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3))-volume_upper(iterate_sphereprobabilitydist);
%                                     else
%                                         volume_lower(iterate_sphereprobabilitydist)=((4/3) * pi*(radius_probability2(iterate_sphereprobabilitydist)^3));
%                                     end
%                                 end
%                             end
%
%                             volume_upper_percentage=(diff([volume_upper 0]).*-1).*percentageactivated(1:length(volume_upper));%volume activated
%                             volume_lower_percentage=(diff([volume_lower 0]).*-1).*percentageactivated(1:length(volume_lower));%volume activated
%                             elect2_primaryactivated=sum(volume_upper_percentage)*Neuron_densities(firstupper2);
%                             elect2_primaryactivated=elect2_primaryactivated+sum(volume_lower_percentage)*Neuron_densities(firstlower2);
%                             area_activated(firstupper2)=area_activated(firstupper2)+(sum(volume_upper_percentage)/volumelayers(firstupper2));%percentage volume activated;;
%                             area_activated(firstlower2)=area_activated(firstlower2)+(sum(volume_lower_percentage)/volumelayers(firstlower2));%percentage volume activated;;
%                         else %whole volume in one layer group
%                             volume_whole_layer=sum(diff([(4/3).* pi.*(radius_probability2.^3) 0]).*-1.*percentageactivated);
%                             area_activated(firstupper2)=area_activated(firstupper2)+((volume_whole_layer)/volumelayers(firstupper2));%percentage volume activated;;
%                             elect2_primaryactivated=volume_whole_layer*Neuron_densities(firstupper2);
%                         end
%
%                         %4. Based on the percentage of primary, what is the
%                         %percentage activated through the secondary
%                         %connections?
%                         %overlappping or non-overlappping activation?
%                         %second
%
%                         connections_activated_secondaryonly=sum(Numconnections_3groups_array.*area_activated,2);
%                         Secondary_neurons_activated=ceil(connections_activated_secondaryonly.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
%                         baselineFR_connections=mean(inhib_excite_arrayFR./inhib_excite_arrayCOUNT);
%                         tertiary_connections=sum((Secondary_neurons_activated./(Neuron_densities'.*volumelayers'))'.*Numconnections_3groups_array,2);
%                         tertiary_neurons_activated=ceil(tertiary_connections.*Neuron_densities'.*volumelayers'./sum(Numconnections_3groups_array,2));
%                         SecondaryTertiary_neurons_activated=Secondary_neurons_activated+tertiary_neurons_activated;
%
%                         %5. weightings based primary and
%                         %secondary or secondary only
%                         %Neuron_densities=[14200, 83800, 164600, 83900,
%                         %177300, 131500]; split into layers
%                         %Neuron_number_groupedlayers=[338+7524 4656 6114+12651].*area_activated;
%                         Primary_neurons_activated=ceil(Neuron_densities.*area_activated.*volumelayers);
%
%                         neurons_distributed_depth=zeros(ceil(layers_depthdefinition(end)/50),1);
%                         %                                                     neurons_distributed_depth(round(stimchn1depth/50)-ceil(radius_activated1/50):round(stimchn1depth/50)+ceil(radius_activated1/50))=ceil(elect1_primaryactivated/(ceil(radius_activated1*2/50)+1));
%                         %                                                     neurons_distributed_depth(round(stimchn2depth/50)-ceil(radius_activated2/50):round(stimchn2depth/50)+ceil(radius_activated2/50))=ceil(elect2_primaryactivated/(ceil(radius_activated2*2/50)+1));
%
%                         for layertypes=1:3
%                             midpoint=(middlelayerpoint(layertypes));
%                             numneurons=SecondaryTertiary_neurons_activated(layertypes);
%                             if numneurons<2
%                                 neurons_distributed_depth(ceil(midpoint/50))=neurons_distributed_depth(ceil(midpoint/50))+SecondaryTertiary_neurons_activated(layertypes);
%                                 neuronposition=ceil(midpoint/50);
%                             elseif rem(numneurons,2)%odd
%                                 distbetweenneurons=(layer_thickness(layertypes)/numneurons);
%                                 neuronposition=(midpoint-(distbetweenneurons)*((numneurons-1)/2)):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons-1)/2);
%                                 neuronposition=ceil(neuronposition/50);
%                                 neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
%                             elseif rem(numneurons+1,2)%even
%                                 distbetweenneurons=(layer_thickness(layertypes)/numneurons);
%                                 neuronposition=(midpoint-(distbetweenneurons)*((numneurons)/2)+distbetweenneurons/2):distbetweenneurons:(midpoint+(distbetweenneurons)*(numneurons)/2);
%                                 neuronposition=ceil(neuronposition/50);
%                                 neurons_distributed_depth(neuronposition)=neurons_distributed_depth(neuronposition)+1;
%                             end
%
%                             %     if length(unique(neuronposition))~=length(neuronposition)
%                             %         error('sampling resolution too large')
%                             %     end
%                         end
%
%
%
%                         weightings=(Primary_neurons_activated)./(Secondary_neurons_activated'+Primary_neurons_activated);
%
%                         %bursting activity due to extracellular stim https://www.sciencedirect.com/science/article/pii/S1935861X09000424#app1
%
%                         burstprob=1./8;%assume cell might fire twice in 8ms window and record for 6ms(2-8)
%                         firingrate=neurons_distributed_depth.*burstprob.*1000;
%                         filtmov=1/5*ones(5,1);
%                         firingrate=filtfilt(filtmov,1,firingrate);
%                         [~,shankposlayers]=min(abs([50:50:1250]-Depthestimate.s2(end)));
%                         simulationFR=NaN(16,1);
%                         simulationFR(1:length(firingrate(shankposlayers:end)))=firingrate(shankposlayers:end);
%                         simulationFR=flip(simulationFR(1:16));
%                         Csplit_simulation.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck)=[Csplit_simulation.(layers_stim).(sepcheck).(currcheck).(trialcheck).(shanksepcheck) [NaN(16-stimChn(1),1); simulationFR; NaN(32-(16-stimChn(1))-16,1)]];
%                     end
%                 end
%             end
%         end
%
%
%     end
% end
%
%
%
%
%
%
% %
% %
% % % Parameters
% % num_exc = 80;
% % num_inh = 20;
% % total_neurons = num_exc + num_inh;
% % duration = 1000; % ms
% % dt = 0.1; % time step in ms
% % time = 0:dt:duration;
% % stim_amplitude = 100; % amplitude of the stimulus in arbitrary units
% %
% % % Define pulse parameters
% % num_pulses = 10; % number of pulses
% % pulse_duration = 0.2; % duration of each pulse in ms
% % inter_pulse_interval = 1000/300; % interval between pulses in ms
% %
% % % Generate pulse train
% % stimulus = zeros(1, length(time));
% % for pulse = 1:num_pulses
% %     start_time = (pulse - 1) * (pulse_duration + inter_pulse_interval);
% %     end_time = start_time + pulse_duration;
% %     stimulus(time >= start_time & time < end_time) = stim_amplitude;
% % end
% %
% % % Neuron properties
% % V_rest = -70; % resting potential in mV
% % V_thresh = -50; % threshold potential in mV
% % V_reset = -70; % reset potential in mV
% %
% % % Synaptic weights (more realistic values)
% % exc_weight = 0.05; % excitatory synaptic weight (arbitrary units)
% % inh_weight = -0.2; % inhibitory synaptic weight (arbitrary units)
% %
% % % Connectivity matrix
% % connectivity = rand(total_neurons) < 0.1; % 10% connectivity
% %
% % % Weight matrix
% % weights = exc_weight * connectivity;
% % weights(:, num_exc+1:end) = inh_weight * connectivity(:, num_exc+1:end);
% %
% % % Membrane potentials
% % V = V_rest * ones(total_neurons, length(time));
% % spikes = zeros(total_neurons, length(time));
% %
% % % Simulation loop
% % for t = 2:length(time)
% %     I_syn = weights * spikes(:, t-1); % synaptic input current
% %     dV = (V_rest - V(:, t-1) + I_syn) * dt; % membrane potential change
% %     V(:, t) = V(:, t-1) + dV; % update membrane potential
% %
% %     % Apply stimulus to excitatory neurons and inhibitory neurons
% %     V(1:num_exc, t) = V(1:num_exc, t) + stimulus(t);
% %     V(num_exc+1:end, t) = V(num_exc+1:end, t) + stimulus(t);
% %
% %     % Spike generation
% %     spikes(:, t) = V(:, t) >= V_thresh;
% %     if sum(spikes(:, t))>0
% %     V(spikes(:, t)==1, t) = V_reset; % reset potential for spiking neurons
% %     end
% % end
% %
% % % Plot results
% % figure;
% % subplot(2, 1, 1);
% % hold on;
% % for i = 1:total_neurons
% %     spike_times = time(spikes(i, :) > 0);
% %     plot(spike_times, i * ones(size(spike_times)), 'k.');
% % end
% % xlabel('Time (ms)');
% % ylabel('Neuron index');
% % title('Spike Raster Plot');
% %
% % subplot(2, 1, 2);
% % plot(time, V(1, :));
% % xlabel('Time (ms)');
% % ylabel('Membrane Potential (mV)');
% % title('Membrane Potential of Neuron 1');
% %
% % %%
% % % % Parameters
% % % num_exc = 80;
% % % num_inh = 20;
% % % total_neurons = num_exc + num_inh;
% % % duration = 1000; % ms
% % % dt = 0.1; % time step in ms
% % % time = 0:dt:duration;
% % % stim_amplitude = 10; % amplitude of the stimulus in arbitrary units
% % %
% % % % Neuron properties
% % % V_rest = -70; % resting potential in mV
% % % V_thresh = -50; % threshold potential in mV
% % % V_reset = -70; % reset potential in mV
% % %
% % % % Synaptic weights
% % % exc_weight = 0.1;
% % % inh_weight = -0.4;
% % %
% % % % Connectivity matrix
% % % connectivity = rand(total_neurons) < 0.1; % 10% connectivity
% % % weights = exc_weight * ones(total_neurons) + inh_weight * (1:total_neurons > num_exc);
% % %
% % % % Membrane potentials
% % % V = V_rest * ones(total_neurons, length(time));
% % % spikes = zeros(total_neurons, length(time));
% % %
% % % % Stimulation protocol (varying pulse trains)
% % % stimulus = stim_amplitude * (mod(time, 200) < 100); % 100 ms on, 100 ms off
% % %
% % % % Simulation loop
% % % for t = 2:length(time)
% % %     I_syn = weights * spikes(:, t-1); % synaptic input current
% % %     dV = (V_rest - V(:, t-1) + I_syn) * dt; % membrane potential change
% % %     V(:, t) = V(:, t-1) + dV; % update membrane potential
% % %
% % %     % Apply stimulus to excitatory neurons
% % %     V(1:num_exc, t) = V(1:num_exc, t) + stimulus(t);
% % %
% % %     % Spike generation
% % %     spikes(:, t) = V(:, t) >= V_thresh;
% % %     V(spikes(:, t), t) = V_reset; % reset potential for spiking neurons
% % end
% %
% % % Plot results
% % figure;
% % subplot(2, 1, 1);
% % hold on;
% % for i = 1:total_neurons
% %     spike_times = time(spikes(i, :) > 0);
% %     plot(spike_times, i * ones(size(spike_times)), 'k.');
% % end
% % xlabel('Time (ms)');
% % ylabel('Neuron index');
% % title('Spike Raster Plot');
% %
% % subplot(2, 1, 2);
% % plot(time, V(1, :));
% % xlabel('Time (ms)');
% % ylabel('Membrane Potential (mV)');
% % title('Membrane Potential of Neuron 1');
%
% %%
% % % Hodgkin-Huxley Neuron Population Model with Synaptic Interactions
% %
% % % Parameters
% % C_m = 1.0;  % �F/cm^2
% % g_Na = 120.0;  % mS/cm^2
% % g_K = 36.0;  % mS/cm^2
% % g_L = 0.3;  % mS/cm^2
% % E_Na = 50.0;  % mV
% % E_K = -77.0;  % mV
% % E_L = -54.387;  % mV
% % E_syn_exc = 0;  % mV, excitatory synapse reversal potential
% % E_syn_inh = -80;  % mV, inhibitory synapse reversal potential
% %
% % % Transition rates
% % alpha_m = @(V) 0.1 * (25 - V) / (exp((25 - V) / 10) - 1);
% % beta_m = @(V) 4.0 * exp(-V / 18);
% % alpha_h = @(V) 0.07 * exp(-V / 20);
% % beta_h = @(V) 1.0 / (exp((30 - V) / 10) + 1);
% % alpha_n = @(V) 0.01 * (10 - V) / (exp((10 - V) / 10) - 1);
% % beta_n = @(V) 0.125 * exp(-V / 80);
% %
% % % Neuron population parameters
% % N = 10;  % Number of neurons
% % exc_ratio = 0.5;  % Proportion of excitatory neurons
% % g_syn = 0.1;  % mS/cm^2, synaptic conductance
% %
% % % Time vector
% % tspan = [0 50];  % ms
% %
% % % Initial conditions
% % V0 = -65.0;  % mV
% % m0 = alpha_m(V0) / (alpha_m(V0) + beta_m(V0));
% % h0 = alpha_h(V0) / (alpha_h(V0) + beta_h(V0));
% % n0 = alpha_n(V0) / (alpha_n(V0) + beta_n(V0));
% % X0 = repmat([V0; m0; h0; n0], N, 1);
% %
% % % Parameters for stimulation
% % pulse_amplitude = 10.0;  % �A/cm^2
% % pulse_duration = 0.2;  % ms
% % pulse_interval = 1000/300;  % ms
% % num_pulses = 3;
% %
% % % Synaptic connectivity matrix
% % syn_conn = rand(N) < 0.2;  % Random sparse connectivity
% % syn_types = rand(N, N) < exc_ratio;  % Excitatory or inhibitory
% %
% %
% %
% % % Solve ODEs for the population
% % [t, X] = ode45(@(t, X) HH_population(t, X, pulse_amplitude, pulse_duration, pulse_interval, num_pulses, C_m, g_Na, E_Na, g_K, E_K, g_L, E_L, alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n, g_syn, E_syn_exc, E_syn_inh, syn_conn, syn_types), tspan, X0);
% %
% % % Plot results for all neurons
% % figure;
% % hold on;
% % for i = 1:N
% %     V = X(:, 4*(i-1) + 1);
% %     plot(t, V);
% % end
% % xlabel('Time (ms)');
% % ylabel('Membrane Potential (mV)');
% % title('Neuron Population Response to Electrical Stimulation');
% % grid on;
% % legend(arrayfun(@(i) ['Neuron ' num2str(i)], 1:N, 'UniformOutput', false));
% % hold off;
% % % ODE function for the population
% % function dXdt = HH_population(t, X, pulse_amplitude, pulse_duration, pulse_interval, num_pulses, C_m, g_Na, E_Na, g_K, E_K, g_L, E_L, alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n, g_syn, E_syn_exc, E_syn_inh, syn_conn, syn_types)
% %     N = length(X) / 4;
% %     dXdt = zeros(size(X));
% %
% %     for i = 1:N
% %         V = X(4*(i-1) + 1);
% %         m = X(4*(i-1) + 2);
% %         h = X(4*(i-1) + 3);
% %         n = X(4*(i-1) + 4);
% %
% %         % External current (stimulation)
% %         I_stim = @(t) pulse_amplitude * any(arrayfun(@(i) t >= (i-1) * pulse_interval && t < (i-1) * pulse_interval + pulse_duration, 1:num_pulses));
% %
% %         % Synaptic current
% %         I_syn = 0;
% %         for j = 1:N
% %             if syn_conn(i, j)
% %                 V_pre = X(4*(j-1) + 1);
% %                 if syn_types(i, j)
% %                     I_syn = I_syn + g_syn * (V - E_syn_exc) * (V_pre > -20);  % Excitatory
% %                 else
% %                     I_syn = I_syn + g_syn * (V - E_syn_inh) * (V_pre > -20);  % Inhibitory
% %                 end
% %             end
% %         end
% %
% %         I_Na = g_Na * m^3 * h * (V - E_Na);
% %         I_K = g_K * n^4 * (V - E_K);
% %         I_L = g_L * (V - E_L);
% %         I_ion = I_Na + I_K + I_L + I_syn;
% %
% %         dVdt = (I_stim(t) - I_ion) / C_m;
% %         dmdt = alpha_m(V) * (1 - m) - beta_m(V) * m;
% %         dhdt = alpha_h(V) * (1 - h) - beta_h(V) * h;
% %         dndt = alpha_n(V) * (1 - n) - beta_n(V) * n;
% %
% %         dXdt(4*(i-1) + 1) = dVdt;
% %         dXdt(4*(i-1) + 2) = dmdt;
% %         dXdt(4*(i-1) + 3) = dhdt;
% %         dXdt(4*(i-1) + 4) = dndt;
% %     end
% % end