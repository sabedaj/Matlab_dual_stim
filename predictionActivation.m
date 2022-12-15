
%potential at r distance from single current source 
%  https://books.google.com.au/books?id=D8DzBJWvRDcC&pg=PA261&lpg=PA261&dq=current+spread+calculation+electrode+poisson&source=bl&ots=r-Q_p8vMku&sig=ACfU3U3RCQKiGwNKPVxGpb957WMFPedM_w&hl=en&sa=X&ved=2ahUKEwjZt5iG2bHoAhV5zjgGHRTjBAMQ6AEwCnoECAoQAQ#v=onepage&q=current%20spread%20calculation%20electrode%20poisson&f=false
%% 23 pairs of stimulating electrodes 
number_pair_elect=23;
a = 200;
b = 500;
%depth_rand = ceil(a.*randn(number_pair_elect,1) + b);

depth_rand = ceil(1 + (1199-1) .* rand(number_pair_elect,1));
depth_rand(depth_rand<0)=1;
%rho=556*100;%resistivity of tissue (ohms per m)
conduct=0.3; %conductivity of tissue

maxdistance=1500e-6; %distance from electrode (radius assuming a homogenous tissue composition and infinite) (m)
InfPot=0; %potential at infinity (V)
a=177e-12; %area at electrode (m^2)
increment=1e-6;%increment of radius (m)
Interelecdist=300e-6; %distance between two stimulating electrodes


% % stimulation supergranular
% sp_s_sup=175;%firing rate supragranular with supergranular stim
% sp_s_gran=250;%granular firing rate
% sp_s_inf=25;%infragranular firing rate
% 
% % stimulation gran
% sp_g_sup=175;%firing rate supragranular
% sp_g_gran=325;%granular firing rate
% sp_g_inf=150;%infragranular firing rate
% 
% % stimulation infragranular
% sp_i_sup=50;%firing rate supragranular
% sp_i_gran=25;%granular firing rate
% sp_i_inf=100;%infragranular firing rate

% stimulation supergranular
sp_s_sup=225;%firing rate supragranular with supergranular stim
sp_s_gran=400;%granular firing rate
sp_s_inf=30;%infragranular firing rate

% stimulation gran
sp_g_sup=300;%firing rate supragranular
sp_g_gran=400;%granular firing rate
sp_g_inf=225;%infragranular firing rate

% stimulation infragranular
sp_i_sup=50;%firing rate supragranular
sp_i_gran=25;%granular firing rate
sp_i_inf=100;%infragranular firing rate

stim_layer=[sp_s_sup sp_g_sup sp_i_sup;sp_s_gran sp_g_gran sp_i_gran;sp_s_inf sp_g_inf sp_g_inf];

%depth layers
s_depth=580;%depth of supragranular layers
g_depth=581:800;%depth of granular layer
i_depth=801:1500;%depth of infragranular layer - 1490 but for simplicity 1500

%Layer thresholds
%20mV avg
% sup_thresh=[22 18 17 14 14 13 13 15]*10^-3;%supragranular threshold for MUA sup,gran and infra stimulation
% gran_thresh=[14 17.5 22 28 27]*10^-3;%granular threshold for MUA sup,gran and infra stimulation
% infra_thresh=[18 17 32 35 29]*10^-3;%infragranular threshold for MUA sup,gran and infra stimulation

%15mV avg
sup_thresh=[16.5 13.5 12.75 10.5 10.5 9.75 9.75 11.25]*10^-3;%supragranular threshold for MUA sup,gran and infra stimulation
gran_thresh=[10.5 13.125 16.5 21 20.25]*10^-3;%granular threshold for MUA sup,gran and infra stimulation
infra_thresh=[13.5 12.75 24 26.25 21.75]*10^-3;%infragranular threshold for MUA sup,gran and infra stimulation

sup_thresh=repelem(sup_thresh,ceil(s_depth/length(sup_thresh)));
sup_thresh(end-abs(length(sup_thresh)-s_depth)+1:end)=[];
gran_thresh=repelem(gran_thresh,ceil((g_depth(end)-g_depth(1))/length(gran_thresh)));
gran_thresh(end-abs(length(gran_thresh)-length(g_depth))+1:end)=[];
infra_thresh=repelem(infra_thresh,ceil((i_depth(end)-i_depth(1))/length(infra_thresh)));
infra_thresh(end-abs(length(infra_thresh)-length(i_depth))+1:end)=[];
Threshold_all=[sup_thresh gran_thresh infra_thresh];
filtmov=1/100*ones(100,1);
Threshold_all=filtfilt(filtmov,1,Threshold_all);
%Threshold_all=0.015*ones(1,1500);

current=[1 2 3 4 6 8 10]./2;

activation_perpair=zeros(number_pair_elect*length(current),i_depth(end));%to make an average over a number of penetrations
for current_it=1:length(current)
    amp=current(current_it);
    I1=amp*1e-6; %current injected by electrode (A)I1
I2=amp*1e-6; %current injected by electrode (A)I2
for e_pos=1:length(depth_rand)
electrode1 =depth_rand(e_pos)*1e-6;%where electrode 1 is positioned (m) 
VrC=zeros(1,round((maxdistance/increment),0));

for r=1:round(((maxdistance*2+Interelecdist)/increment),0)
    for i=1:round((maxdistance*2/increment),0)
            VrC(r,i)=(I1/(4*pi*conduct*(sqrt((maxdistance+Interelecdist+increment-r*increment)^2+(maxdistance+increment-i*increment)^2))))+(I2/(4*pi*conduct*(sqrt((maxdistance+increment-r*increment)^2+(maxdistance+increment-i*increment)^2)))); %voltage at distance r from electrode
    end
end

r=1:round(((maxdistance*2+Interelecdist)/increment),0);
i=1:round((maxdistance*2/increment),0);
VrC((round((maxdistance+Interelecdist)/increment)+1),(round(maxdistance/increment))+1)=VrC((round((maxdistance+Interelecdist)/increment)),(round(maxdistance/increment)));
VrC((round((maxdistance)/increment)+1),(round(maxdistance/increment))+1)=VrC((round((maxdistance)/increment)),(round(maxdistance/increment)));
% 
 figure(1)
mesh(i,r,VrC)
set(gca,'colorscale','log')
xlabel('Distance (\mum)');
ylabel('Distance (\mum)');
zlabel('Voltage (V)');
title ('Voltage distribution from two point source electrodes')

figure (2)
plot(VrC((1:round(((maxdistance*2+Interelecdist)/increment))),round(maxdistance/increment)+1));
xlabel('Distance (\mum)');
ylabel('Voltage (V)');
title ('Voltage distribution from two point source electrodes')

VrC_act=VrC;
% points=100;
% filtmov=1/points*ones(points,1);
% for filtloop=1:size(VrC_act,1)
% VrC_act(filtloop,:)=filtfilt(filtmov,1,VrC_act(filtloop,:));
% end
% for filtloop=1:size(VrC_act,2)
% VrC_act(:,filtloop)=filtfilt(filtmov,1,VrC_act(:,filtloop));
% end
% figure ; mesh(VrC_act)

VrC_act(1:round((electrode1-maxdistance)*-1*10^6),:)=[];
VrC_act(end-(size(VrC_act,1)-(i_depth(end))-1):end,:)=[];

figure (3)
plot(VrC_act(:,round(maxdistance/increment)+1));
xlabel('Distance (\mum)');
ylabel('Voltage (V)');
title ('Voltage distribution from two point source electrodes')

checkanyactivate=0;
for row=1:size(VrC_act,1)
    for column=1:size(VrC_act,2)
        if VrC_act(row,column)>Threshold_all(row)%threshold for activation is approximately 20mV change in extracellular potential
            VrC_act(row,column)=Threshold_all(row);
        else 
            checkanyactivate=checkanyactivate+1;
        end
                VrC_act(row,column)=VrC_act(row,column)/Threshold_all(row); %normalise voltage to convolve with likelihood activation
    end
end
if checkanyactivate==size(VrC_act,2)*size(VrC_act,1)
    error()
end
% tmin=0;
%         tmax=1;
%         rmax=max(VrC_act,[],'all');
%         rmin=min(VrC_act,[],'all');
%         VrC_act=(VrC_act-rmin)./(rmax-rmin).*(tmax-tmin)+tmin; %normalise voltage to convolve with likelihood activation
% 


figure (4)
plot(VrC_act(:,round(maxdistance/increment)+1));
xlabel('Distance (\mum)');
ylabel('Probability of evoking neural activity');
title ('Voltage distribution from two point source electrodes')




layers=zeros(1,3);
% 
% if electrode1<s_depth*1e-6%E1 in supragran layers 
%     layers(1)=1;
% elseif electrode1<g_depth(end)*1e-6%E1 in gran layers
%     layers(2)=1;
% elseif electrode1<i_depth(end)*1e-6 %E1 in infragran layers
%     layers(3)=1;
% else %E1 is at end of cortex so E2 in white matter
%     
% end
[r,c]=find(VrC_act==1);
if any(r<s_depth) %E2 in supragran layers
    layers(1)=1;
    r(r<s_depth)=1700;
end
if any(r<g_depth(end)) %E2 in gran layer
    layers(2)=1;
      r(r<g_depth(end))=1700;
end
if any(r<i_depth(end)) %E2 in infragran layers
    layers(3)=1;
end
activatedlayer=stim_layer(:,logical(layers));
activatedlayer=max(activatedlayer,[],2);


laminar_sp_max=zeros(size(VrC_act,1),size(VrC_act,2));
laminar_sp_max(1:s_depth,:)=activatedlayer(1);
laminar_sp_max(g_depth(1):g_depth(end),:)=activatedlayer(2);
laminar_sp_max(i_depth(1):i_depth(end),:)=activatedlayer(3);
filtmov=1/200*ones(200,1);
out=filtfilt(filtmov,1,laminar_sp_max(:,1));
for row=1:size(laminar_sp_max,2)
laminar_sp_max(:,row)=out;
end
 Laminar_activated_overall=laminar_sp_max.*VrC_act;
 filtmov=1/100*ones(100,1);
out=filtfilt(filtmov,1,laminar_sp_max(:,1));
figure(5)
mesh(i,1:i_depth(end),laminar_sp_max)
set(gca,'colorscale','log')
xlabel('Distance (\mum)');
ylabel('Distance (\mum)');
zlabel('Sp/s');
title (['Average spike rate from stimulation. Electrode depth ' num2str((electrode1)*1e-6) '\mum & ' num2str((electrode1+Interelecdist)*1e-6) '\mum'])

% 
figure(6)
mesh(i,1:i_depth(end),VrC_act)
set(gca,'colorscale','log')
xlabel('Distance (\mum)');
ylabel('Distance (\mum)');
zlabel('Probability of eliciting spiking');
title ('Two point source electrodes')
% 
% 
figure(7)
mesh(i,1:i_depth(end),Laminar_activated_overall)
set(gca,'colorscale','log')
xlabel('Distance (\mum)');
ylabel('Distance (\mum)');
zlabel('Sp/s');
title ('Neural activation laminar direction from two point source electrodes')

figure(8)
plot(Laminar_activated_overall(:,maxdistance*1e6))
xlabel('Distance (\mum)');
ylabel('Sp/s');
title ('Neural activation laminar direction from two point source electrodes')


activation_perpair(e_pos+(current_it-1)*number_pair_elect,:)=Laminar_activated_overall(:,maxdistance*1e6);

end
end
%% 23 pairs of stimulating electrodes 
adepth_adjust=zeros(number_pair_elect*length(current),3000);
loopnum=0;
    for i=1:number_pair_elect*length(current)
        loopnum=loopnum+1;
        loopnum(loopnum>number_pair_elect)=1;
        depthincrement=loopnum;
        adepth_adjust(i,:)=[zeros(1,1200-depth_rand(loopnum))  activation_perpair(i,:) zeros(1,3000-length(activation_perpair(i,:))-(1200-depth_rand(loopnum)))];
    end
%%
    removezeros=adepth_adjust;%adepth_adjust;
removezeros(removezeros==0)=NaN;
    for i=1:size(removezeros,2)
        testlength=removezeros(:,i);
        nonan_len(i) = length(testlength(~isnan(testlength)));
        leng=100; %set to number of trials to avoid edge effects
        nonan_len(nonan_len<leng)=0;
         if nonan_len(i)==0
             removezeros(:,i)=NaN;
            continue
         end
          input=find(~isnan(testlength));
         s = RandStream('mlfg6331_64','seed','shuffle');
        [y1,idx] = datasample(s,input,(nonan_len(i)-leng),'Replace',false);
        removezeros(y1,i)=NaN;
    end
%%

figure
plot(nanmean(removezeros,1))
figure
plot(removezeros(161,:))
%%
x = [0:1:1490];
current = [1,0.75,0.5,0.25,0];
spread_range = [5,5,5,5];
plot_type = 1; % 1= current, 2 = spread;
cond=zeros(length(current),length(x)*2);
electrode1 =11*50;
electrode2 =18*50;
for i = 1:length(current)

if(plot_type == 1)
    current1 = current(i);
    current2 = 1-current(i);
    spread = 150; %random
else
    electrode1 = 2+i-1;
    current1 = 1;
    current2 = 1;
    spread = spread_range(i);
end
y1 = current1*pdf('Normal',x,electrode1,spread);%.*x./801;
y2 = current2*pdf('Normal',x,electrode2,spread);%.*x./801;
cond(i,1:length(x)*2)=[y1,y2];
end
y1=cond(1:length(current),1:length(x));
y2=cond(1:length(current),length(x)+1:length(x)*2);
y1y2=y1+y2;
figure
for i = 1:length(current)
subplot(2,3,i)
box off
hold
axis square
plot(x,y1(i,:)./max([y1;y2;y1y2]),'k--')
plot(x,y2(i,:)./max([y1;y2;y1y2]),'k:')
plot(x,y1y2(i,:)./max([y1;y2;y1y2]),'k')
xlabel('Distance')
ylabel('Normalised spike rate')
xline(electrode1,'r')
xline(electrode2,'r')
end
figure
for i = 1:length(current)
subplot(2,3,i)
box off
hold
axis square
plot(x,y1(i,:),'k--')
plot(x,y2(i,:),'k:')
plot(x,y1y2(i,:),'k')
xlabel('Distance')
ylabel('Normalised spike rate')
xline(electrode1,'r')
xline(electrode2,'r')
end
%% bias
for i = 1:length(current)
if(plot_type == 1)
    current1 = current(i);
    current2 = 1-current(i);
    spread = 150; %random
else
    electrode1 = 2+i-1;
    current1 = 1;
    current2 = 1;
    spread = spread_range(i);
end
if electrode1<690
    coefficients = polyfit([601, 700], [1/8*current1*pdf('Normal',x(600),electrode1,spread) 0], 1);
    y1 = [1/8*current1*pdf('Normal',x(1:600),electrode1,spread) coefficients(1).*x(601:700)+coefficients(2) current1*1/20*pdf('Normal',x(701:end),1200,spread)];%.*x./801;
else
    coefficients = polyfit([601, 700], [0, current1*pdf('Normal',x(701),electrode1,spread*3)], 1);
    y1 = [zeros(1,600) coefficients(1).*x(601:700)+coefficients(2) current1*pdf('Normal',x(599:end),1200,spread*3)];%.*x./801;
end
if electrode2<690
    coefficients = polyfit([601, 700], [1/8*current1*pdf('Normal',x(600),electrode2,spread) 0], 1);
    y2 = [1/8*current2*pdf('Normal',x(1:600),electrode2,spread) coefficients(1).*x(601:700)+coefficients(2) current2*1/20*pdf('Normal',x(701:end),1200,spread)];%.*x./801;
else
    coefficients = polyfit([601, 700], [0, current2*pdf('Normal',x(701),electrode2,spread*3)], 1);
    y2 = [zeros(1,600) coefficients(1).*x(601:700)+coefficients(2) current2*pdf('Normal',x(701:end),electrode2,spread*3)];%.*x./801;
end
cond(i,1:length(x)*2)=[y1,y2];
end
y1=cond(1:length(current),1:length(x));
y2=cond(1:length(current),length(x)+1:length(x)*2);
y1y2=y1+y2;
    
    
    
biasactivation=[0:1/450:1 ones(1,450)];
figure
for i = 1:length(current)
subplot(2,3,i)
box off
hold
axis square

plot(x,y1(i,:),'k--')
plot(x,y2(i,:),'k:')
plot(x,y1y2(i,:),'k')
xlabel('Distance')
ylabel('Normalised spike rate')
xline(electrode1,'r')
xline(electrode2,'r')
ylim([0 2*10^-3])
end



figure
for i = 1:length(current)
subplot(2,3,i)
box off
hold
axis square
y1y2filt=smooth(y1y2(i,:)./max([y1;y2;y1y2]));
plot(x,y1(i,:)./max([y1;y2;y1y2]),'k--')
plot(x,y2(i,:)./max([y1;y2;y1y2]),'k:')
plot(x,y1y2filt,'k')
xlabel('Distance')
ylabel('Normalised spike rate')
xline(electrode1,'r')
xline(electrode2,'r')

end
