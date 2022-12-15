
%% PLOTTING
E_Mapnumber=1;
E_MAP = Depth(E_Mapnumber);
load('meanlfp.mat','meanlfpstruct')
%load('meanresp.mat','meanrespstruct')
avgnospT=meanlfpstruct(E_MAP,:);
figure
subplot(1,4,1)
hold on
title('Shank 1')
xlabel('Time (ms)')
%ylabel('uV')
sepdist=150;
for i=16:-1:1
plot(-250:1500-1,avgnospT(i,:)-i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-50 1500])
ylim ([-sepdist*16-sepdist -sepdist+sepdist])
x=ones(1,101).*100;
plot(x, -150:-50, 'k')
subplot(1,4,4)
hold on
title('Shank 4')
xlabel('Time (ms)')
%ylabel('uV')
for i=32:-1:17
plot(-250:1500-1,avgnospT(i,:)-i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-50 1500])
ylim ([-sepdist*32-sepdist -sepdist*17+sepdist])

subplot(1,4,2)
hold on
title('Shank 2')
for i=48:-1:33
plot(-250:1500-1,avgnospT(i,:)-i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xlabel('Time (ms)')
%ylabel('uV')
xline(1111.11,'r')
xlim([-50 1500])
ylim ([-sepdist*48-sepdist -sepdist*33+sepdist])

subplot(1,4,3)
hold on
title('Shank 3')
xlabel('Time (ms)')
%ylabel('uV')
for i=64:-1:49
plot(-250:1500-1,avgnospT(i,:)-i*sepdist,'k')
end
set(gca,'YTickLabel',[]);
xline(0,'r')
xline(1111.11,'r')
xlim([-50 1500])
ylim ([-sepdist*64-sepdist -sepdist*49+sepdist])



plot(x,-4950:-4850, 'k')

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