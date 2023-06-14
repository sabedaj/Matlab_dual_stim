function plotV1V2matrix(Data,trial)

E_map=Depth;
V2elect=reshape(E_map(1:64),16,4);
V1elect=reshape(E_map(65:128),16,4);
colorder=[1 3 4 2];
V1elect=V1elect(:,colorder);
V2elect=V2elect(:,colorder);
close all
for chn=1:128
    tmpdat=mean(squeeze(Data.(['T' num2str(trial)])(chn,1:181,:)),2,'omitnan');
             window = normpdf(-3:3,0,1);
    rate = (1000)*conv(tmpdat,window);
    rate = rate(3+1:end-3);
    if chn<65
        x=find(V2elect==chn);
        figure(2)
        hold on
        subplot(4,16,x)
        plot(-90:90,rate)
        xlim([-85,85])
    else
        x=find(V1elect==chn);
        figure(1)
        hold on
        subplot(4,16,x)
        plot(-90:90,rate)
         xlim([-85,85])
    end
    
    
end


end