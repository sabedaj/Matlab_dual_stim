function maxcor=noisecorrelations(ElectLayerClass)

sp=loadSpikes;
%load('ElectLayerClass.mat');
%%
[amplifier_channels,frequency_parameters]=read_Intan_RHS2000_file;
nChn=size(amplifier_channels,2);
FS=frequency_parameters.amplifier_sample_rate;

order = Depth(1);
spsort=sp(order);
maxcor=nan(64,64);
lagcor=zeros(64,64);
%enter trialnums
%channels
maxtimecor=10000;%ms


fffb=cell(4,4);
combsave=[0 0];
for chn1=1:nChn
    spchn1=spsort{chn1};
    xdata=[];
    
    theseSp = (spchn1(spchn1(:,1) > 0 & spchn1(:,1) < maxtimecor));
    for i = 1:length(theseSp)
        xdata = [xdata, (theseSp(i))];
    end
    ratechn1 = hist(xdata,0:2:maxtimecor); %raw no smoothing
    
    for chn2=1:nChn
        xdata=[];
        spchn2=spsort{chn2};
        
        theseSp = (spchn2(spchn2(:,1) > 0 & spchn2(:,1) < maxtimecor));
        for i = 1:length(theseSp)
            xdata = [xdata, (theseSp(i))];
        end
        
        ratechn2 = hist(xdata,0:2:maxtimecor); %raw no smoothing
        
        [xc,~] = xcorr(ratechn1,ratechn2,'coeff');
%         
%                 if chn1~=chn2 && ~any(ismember(combsave,sort([chn1,chn2]),'rows'))
%                         index=sort([ElectLayerClass(chn2),ElectLayerClass(chn1)]);
%                         fffb{index(1),index(2)}=[fffb{index(1),index(2)} max(xc)];
%         
%                 end
%                 combsave=[combsave; sort([chn1,chn2])];
        if ~all(isnan(xc))
            maxcor(chn1,chn2)=max(xc);
            lagcor(chn1,chn2)=find(max(xc)==xc,1,'first');
        end
        
        
    end
    
    
    
    
    
    
    
end
end