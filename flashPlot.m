
%% PLOTTING
colorsplotlayer=[1 0 0; 0 1 0; 0 0 1; 1 1 0];
E_Mapnumber=1;
E_MAP = Depth(E_Mapnumber);
load('meanlfp.mat','meanlfpstruct')
%load('meanresp.mat','meanrespstruct')
avgnospT=meanlfpstruct(E_MAP,:);
%avgnospT=avgnospT(1:16,30*150:30*550);

filefold=pwd;
parentDirectory = fileparts(cd);
cd(parentDirectory)
try
    load('ElectLayerClass.mat','ElectLayerClass')
catch
    ElectLayerClass=ones(64);
end
cd(filefold)
figure
subplot(1,4,1)
hold on
title('Shank 1')
xlabel('Time (ms)')
%ylabel('uV')
sepdist=75;
xvalAX=-250:1500-1;%-250:1/30:500-1/30;
for i=1:1:16
plot(xvalAX,avgnospT(i,:)+i*sepdist,'color',colorsplotlayer(ElectLayerClass(i),:))
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-100 300])
ylim ([sepdist-sepdist-25 sepdist*16+sepdist])
set(gca,'TickDir','out');
x=ones(1,101).*100;
plot(x, -150:-50, 'k')
subplot(1,4,4)
hold on
title('Shank 4')
xlabel('Time (ms)')
%ylabel('uV')
for i=17:1:32
plot(xvalAX,avgnospT(i,:)+i*sepdist,'color',colorsplotlayer(ElectLayerClass(i),:))%plot(-250:1500-1,avgnospT(i,:)+i*sepdist,'color',colorsplotlayer(ElectLayerClass(i),:))
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-100 300])
ylim ([sepdist*17-sepdist-25 sepdist*32+sepdist])
set(gca,'TickDir','out');
subplot(1,4,2)
hold on
title('Shank 2')
for i=33:48
plot(xvalAX,avgnospT(i,:)+i*sepdist,'color',colorsplotlayer(ElectLayerClass(i),:))
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xlabel('Time (ms)')
%ylabel('uV')
xline(1111.11,'r')
xlim([-100 300])
ylim ([sepdist*33-sepdist-25 sepdist*48+sepdist])
set(gca,'TickDir','out');
subplot(1,4,3)
hold on
title('Shank 3')
xlabel('Time (ms)')
%ylabel('uV')
for i=49:64
plot(xvalAX,avgnospT(i,:)+i*sepdist,'color',colorsplotlayer(ElectLayerClass(i),:))
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-100 300])
ylim ([sepdist*49-sepdist-25 sepdist*64+sepdist])
set(gca,'TickDir','out');


plot(x,1200:1300, 'k')

% figure
% for chsp=1:1:nChn
%             check=['T', num2str(chsp)];
%             
%             hold on
%             plot(-250:1000/30000:500-1000/30000,LFPstruct.(check)(1:200:end,:),'k')
% end
% figure
% for chsp=1:1:nChn
%             check=['T', num2str(chsp)];
%             
%             hold on
%             plot(-250:1000/30000:500-1000/30000,meanlfpstruct.(check),'k')
% end
% xlabel('Time (ms)')
% ylabel('uV')
% title('Average LFP waveform of each electrode timelocked to flash stimuli n_r_e_p=900')
% xline(0,'r')