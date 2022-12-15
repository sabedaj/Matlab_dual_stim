function [maxcor,lagcor]=ratecorrelations(typecor,ID)
%typecor - 'noise','stim'
sp=loadSpikes;
%load('ElectLayerClass.mat');
%%
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

order = Depth(1);
spsort=sp(order);
maxcor=zeros(64,64);
lagcor=zeros(64,64);
%enter trialnums 
%channels

if strcmp(typecor,'noise')
    BIN = [-30 -10]; %time ms to look

elseif strcmp(typecor,'stim')
    BIN = [0 20]; %time ms to look


end

trig = loadTrig(0);
TP = loadTrialParams;

tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
chns=unique(cell2mat(TP(cell2mat(TP(:,2)) == ID,3)));
theseTrig = trig(tID)./30;
nT=length(theseTrig);
fffb=cell(4,4);
combsave=[0 0];
for chn1=1:nChn
    spchn1=spsort{chn1};
    xdata=[];
    for tr=1:nT
        theseSp = (spchn1(spchn1(:,1) > theseTrig(tr)+BIN(1) & spchn1(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
        for i = 1:length(theseSp)
            xdata = [xdata, (theseSp(i))];
        end
    end
    ratechn1 = hist(xdata,BIN(1):2:BIN(2)); %raw no smoothing

    for chn2=1:nChn
        xdata=[];
        spchn2=spsort{chn2};
        for tr=1:nT
            theseSp = (spchn2(spchn2(:,1) > theseTrig(tr)+BIN(1) & spchn2(:,1) < theseTrig(tr)+BIN(2)) - theseTrig(tr));
            for i = 1:length(theseSp)
                xdata = [xdata, (theseSp(i))];
            end
        end
        ratechn2 = hist(xdata,BIN(1):2:BIN(2)); %raw no smoothing

        [xc,~] = xcorr(ratechn1,ratechn2,'coeff');

%         if chn1~=chn2 && ~any(ismember(combsave,sort([chn1,chn2]),'rows'))
%                 index=sort([ElectLayerClass(chn2),ElectLayerClass(chn1)]);
%                 fffb{index(1),index(2)}=[fffb{index(1),index(2)} max(xc)];
% 
%         end
%         combsave=[combsave; sort([chn1,chn2])];
        if ~all(isnan(xc))
        maxcor(chn1,chn2)=max(xc);
        lagcor(chn1,chn2)=find(max(xc)==xc,1,'first');
        end
    end

end



% figure
% surf(maxcor)
% title(['Maximum correlation ' typecor ' ID:' num2str(ID) ' Channel: ' num2str(chns(2))])
% figure
% surf(lagcor)
% title(['Lag of maximum correlation ' typecor ' ID:' num2str(ID) ' Channel: ' num2str(chns(2))])
% 


%try smoothing or different time windows
%
%[nanmean(fffb{1,1}),nanmean(fffb{1,2}),nanmean(fffb{1,3}); nanmean(fffb{2,1}),nanmean(fffb{2,2}),nanmean(fffb{2,3});nanmean(fffb{3,1}),nanmean(fffb{3,2}),nanmean(fffb{3,3})]


end