function [Sp,Spktimes,Thresh]=OnlineSpkExtract(Ampsignal,Thresh)
%% "online" spike analysis
fs=30000;%Hz
DT = (fs/1e3)*1.6;
%high pass filter
[FilterCoef,~] = generate_Filters;
sig_filt=FiltFiltM(FilterCoef,1,Ampsignal);
%threshold spikes
if Thresh==1
SD=std(sig_filt(1:fs*19));
Thresh=-4.5*SD;
end
art_reject=-1000;
tmu= zeros(1,size(sig_filt,2));
tmu((sig_filt < Thresh) & (sig_filt > art_reject))=1;%get rid of artefact triggers
Spktimes = find(diff(tmu)>0)';
n = find(Spktimes > DT & Spktimes < size(sig_filt,2)-DT);
Spktimes = Spktimes(n);
nsp = length(n);


%find spk wvfms
if nsp>0
    base = ones(nsp,1)*[-DT:DT];
    align = Spktimes(:,ones(1,size(base,2)));
    intspike = resample(reshape(sig_filt(base+align),size(align))',4,1);  % Upsampling
    neg = find(intspike(4*DT,:)<0);
    [~,negind] = min(intspike(3*DT:5*DT,neg));    %  Negative threshold
    pos = find(intspike(4*DT,:)>0);
    [~,posind] = max(intspike(3*DT:5*DT,pos));    %  Positive threshold
    [~,interleaveind]=sort([neg pos],'ascend');
    negposind = [negind posind];
    negposind = negposind(interleaveind);
    ind = negposind + 3*DT-1;
        ind = ind';
    % Sp = intspike;
    %  The following is to avoid indexing spikes at edge of spike extraction
    %  window
    n = find(ind > DT+1 & ind < size(intspike,1)-2*DT);  %  What are the right numbers here for 8x oversampling?
    %disp(['  ' num2str(length(n)) ' spikes']);
    if n > 0
        base = ones(length(n),1)*[-DT:3*DT];
        x = 1:length(n);
        xx = x(ones(1,size(base,2)),:)*size(intspike,1);
        xx = xx'-size(intspike,1);
        align = ind(n,ones(1,size(base,2)));
        index = base+align+xx;
        index = index';
        intspike = intspike(:,n);
        Sp = intspike(index(1:4:end,:));            %  Downsampling
        % Adjustment to how this works
        % We know the minimum of each wave - this should be centered
        % And spktime should reflect the minima of each wave
        % This will allow for removal of duplicate threshold crossings        
        % spktimes are taken from index DT*4+1 of intspike
        % ind has the intspike minimum
        Spktimes = Spktimes + (ind - (DT*4+1))./4;
        Spktimes = Spktimes(n')./fs*1000;
        if nargin ==5
            [~,c]=find(Sp>artefact_high);
            Sp(:,c)=[];
            Spktimes(c)=[];
        end
    else
        Sp = zeros(1,DT+1);
        Spktimes = zeros(1,1);
    end
else
    %Sp = zeros(1,DT+1);
    Sp = [];
    Spktimes = zeros(1,1);
end
Sp =Sp';
% try to remove artefact spikes
%remove any spike whose min is not at sample 13, if the sample is positive
%at 15, if the beginning and end of the spike file are more than 100uV
%apart
if ~isempty(Sp)
Spktimes(abs(Sp(:,1)-Sp(:,end))>100)=[];
Sp(abs(Sp(:,1)-Sp(:,end))>100,:)=[];
[r,c]=find(Sp==min(Sp,[],2));
Spktimes(r(c~=13))=[];
Sp(r(c~=13),:)=[];
Spktimes(Sp(:,14)>0)=[];
Sp(Sp(:,14)>0,:)=[];
end
% %%plotting
% figure;hold on;
% for i=1:size(Sp,1)
%      plot(Sp(i,2:end))
% end
%%
end

