%setup the structure to pool the data
endloop=0.155;
step=0.01;
rsclass=[];
for i=0:step:endloop
    
    ncheck=['T' num2str(i^2)];
    old = ".";
    new = "_";
    check = replace(ncheck,old,new);
    rsclass.(check)=[];
end


% go through each animal
for ratN=6:11
if ratN<10
    Ratnum=['Rat_00' num2str(ratN)];
else
    Ratnum=['Rat_0' num2str(ratN)];
end
cd(['E:\DATA' filesep Ratnum])%change working directory to where data is stored - currently manually input
D_data=dir;
load('Ratio_all.mat')
 %load data
currentvariation50_all.(Ratnum)=currentvariation50;
currentvariation75_all.(Ratnum)=currentvariation75;
currentvariation25_all.(Ratnum)=currentvariation25;
stimchnpair_all.(Ratnum)=Stimchnall;

% find laminar channel pairs
for chnpair=1:length(Stimchnall)
    if Stimchnall(chnpair,1)<17&&Stimchnall(chnpair,2)<17
    elseif Stimchnall(chnpair,1)<33&&Stimchnall(chnpair,2)<33 && Stimchnall(chnpair,1)>16&&Stimchnall(chnpair,2)>16
    elseif Stimchnall(chnpair,1)<49&&Stimchnall(chnpair,2)<49 && Stimchnall(chnpair,1)>32&&Stimchnall(chnpair,2)>32
    elseif Stimchnall(chnpair,1)<65&&Stimchnall(chnpair,2)<65 && Stimchnall(chnpair,1)>48&&Stimchnall(chnpair,2)>48
    else
        %across
        Stimchnall(chnpair,:)=[0,0];
    end
end
Stimchnall(Stimchnall(:,1)-Stimchnall(:,2)~=-6,:)=0; %-6 is a 5 electrode spacing
StimchnNoac_all.(Ratnum)=Stimchnall;

%Determine voltage based on distance for each electrode from the stimulating electrodes
distpair=[];
stimchns=StimchnNoac_all.(Ratnum);%change accordingly
AMP=[1 2 3 4 6 8 10].*10^-6;
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_all.(Ratnum)=distpair;

%take care of input and make it reasonable. Remove null placeholders
rinput=currentvariation50_all.(Ratnum);
rinput(all(rinput == 0,2),:)=-500;
rinput(abs(rinput)>10000)=-500;


%put into a common structure for all rodents
D=[];
stimchns=StimchnNoac_all.(Ratnum);%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); rinput(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end

end









%% assuming 5 electrode separation
% StimchnNoac_r5(StimchnNoac_r5(:,1)-StimchnNoac_r5(:,2)~=-6,:)=0;
% StimchnNoac_r6(StimchnNoac_r6(:,1)-StimchnNoac_r6(:,2)~=-6,:)=0;
% StimchnNoac_r7(StimchnNoac_r7(:,1)-StimchnNoac_r7(:,2)~=-6,:)=0;
% StimchnNoac_r8(StimchnNoac_r8(:,1)-StimchnNoac_r8(:,2)~=-6,:)=0;
% StimchnNoac_r9(StimchnNoac_r9(:,1)-StimchnNoac_r9(:,2)~=-6,:)=0;
% 
% StimchnNoac_r11(StimchnNoac_r11(:,1)-StimchnNoac_r11(:,2)~=-6,:)=0;

%% four shank distance from stim electrodes
distpair=[];
stimchns=StimchnNoac_all.Rat_006;%change accordingly
AMP=[1 2 3 4 6 8 10].*10^-6;
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_007;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r7=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_008;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_009;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r9=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_010;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r10=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_011;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r11=distpair;
%% single shank distance from stimulating electrodes
% distpair=[];
% AMP=[1 2 3 4 6 8]*10^-6;
% stimchns=StimchnNoac_r5;%change accordingly
% for currentsi=1:length(AMP)
%     currents=AMP(currentsi);
%     for i=1:size(stimchns,1)
%         distanceEd=ones(32,1).*0.00000000001;
%         distanceEs=ones(32,1).*0.00000000001;
%         if stimchns(i,1)~=0
%             for j=1:stimchns(i,1)-1
%                 distanceEd(stimchns(i,1)-j,1)=50*j;
%             end
%             
%             for j=1:32-stimchns(i,1)
%                 distanceEd(stimchns(i,1)+j,1)=50*j;
%             end
%             for j=1:stimchns(i,2)-1
%                 distanceEs(stimchns(i,2)-j,1)=50*j;
%             end
%             
%             for j=1:32-stimchns(i,2)
%                 distanceEs(stimchns(i,2)+j,1)=50*j;
%             end
%             distanceAdditive=(((currents*0.5)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.5)./(distanceEd*10^-6*0.3*2*pi)));
%             check=['T' num2str(i) num2str(currents/(10^-6))];
%             distpair.(check)=distanceAdditive;
%         end
%     end
% end
% distpair_r5=distpair;
%% Put together

% r5input=currentvariation50_r5;
% r5input(all(r5input == 0,2),:)=-500;
% r5input(abs(r5input)>10000)=-500;
r6input=currentvariation50_all.Rat_006;
r6input(all(r6input == 0,2),:)=-500;
r6input(abs(r6input)>10000)=-500;
r7input=currentvariation50_all.Rat_007;
r7input(all(r7input == 0,2),:)=-500;
r7input(abs(r7input)>10000)=-500;
r8input=currentvariation50_all.Rat_008;
r8input(all(r8input == 0,2),:)=-500;
r8input(abs(r8input)>10000)=-500;
r9input=currentvariation50_all.Rat_009;
r9input(all(r9input == 0,2),:)=-500;
r9input(abs(r9input)>10000)=-500;
r10input=currentvariation50_all.Rat_010;
r10input(all(r10input == 0,2),:)=-500;
r10input(abs(r10input)>10000)=-500;
r11input=currentvariation50_all.Rat_011;
r11input(all(r11input == 0,2),:)=-500;
r11input(abs(r11input)>10000)=-500;
D=[];
% stimchns=StimchnNoac_r5;%change accordingly
% if (32-max(stimchns,[],'all'))>(min(stimchns(stimchns~=0),[],'all'))
%     [r,c]=find(min(stimchns(stimchns~=0),[],'all')==stimchns);
% else
%     [r,c]=find(max(stimchns,[],'all')==stimchns);
% end
% distinterest=distpair_r5.(['T' num2str(r)]);

endloop=0.155;
step=0.01;
rsclass=[];
for i=0:step:endloop
    
    ncheck=['T' num2str(i^2)];
    old = ".";
    new = "_";
    check = replace(ncheck,old,new);
    rsclass.(check)=[];
end
%% single
% D=[];
% AMP=[1 2 3 4 6 8]*10^-6;
% for currentsi=1:length(AMP)
%     currents=AMP(currentsi);
%     for stimchncounter=1:size(stimchns,1)
%         check=['T' num2str(stimchncounter) num2str(currents*10^6)];
%         try
%             distinterest=distpair_r5.(check);
%             
%             
%             for i=0:step:endloop
%                 ncheck=['T' num2str(i^2)];
%                 old = ".";
%                 new = "_";
%                 check = replace(ncheck,old,new);
%                 if i==endloop
%                     D.(check)=find((distinterest>=(((i-step)^2))));
%                 elseif i==step
%                     D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
%                 elseif i~=endloop
%                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
%                 end
%                 if D.(check)~=0
%                     row=D.(check)+(stimchncounter-1)*32;
%                     rsclass.(check)=[rsclass.(check); r5input(row,currentsi+1)];
%                 end
%             end
%             
%         catch
%             %not laminar or suitable sep dist
%         end
%         D=[];
%     end
% end

%% four shank
%6
AMP=[1 2 3 4 6 8 10].*10^-6;
D=[];
stimchns=StimchnNoac_all.Rat_006;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair_r6.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); r6input(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end
%7
D=[];
stimchns=StimchnNoac_all.Rat_007;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair_r7.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                    D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); r7input(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end

%8
D=[];
stimchns=StimchnNoac_all.Rat_008;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair_r8.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); r8input(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end

%9
D=[];
stimchns=StimchnNoac_all.Rat_009;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair_r9.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); r9input(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end

%9
D=[];
stimchns=StimchnNoac_all.Rat_010;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair_r10.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); r10input(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end

%11
D=[];
stimchns=StimchnNoac_all.Rat_011;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for stimchncounter=1:size(stimchns,1)
        check=['T' num2str(stimchncounter) num2str(currents*10^6)];
        try
            distinterest=distpair_r11.(check);%change accordingly
            
            for i=0:step:endloop
                ncheck=['T' num2str(i^2)];
                old = ".";
                new = "_";
                check = replace(ncheck,old,new);
                if i==endloop
                    D.(check)=find((distinterest>=(((i-step)^2))));
                elseif i==step
                    D.(check)=find((distinterest<i^2)&(distinterest>=((i^2-(step)^2))));
                elseif i~=endloop
                     D.(check)=find((distinterest<i^2)&(distinterest>=(((i-step)^2))));
                end
                if D.(check)~=0
                    row=D.(check)+(stimchncounter-1)*64;
                    rsclass.(check)=[rsclass.(check); r11input(row,currentsi+1)];%change accordingly
                end
            end
            
        catch
            %not laminar or suitable sep dist
        end
        D=[];
    end
end
%% Results calc without currents
resultind=zeros(4,length(fieldnames(rsclass)));
currents=1;
for i=1:length(fieldnames(rsclass))
    str=fieldnames(rsclass);
    str=char(str(i));
    meansig50array=rsclass.(str);
    if ~isempty(meansig50array)
        %meansig50array=meansig50array(:,2:end);
        meansig50array(meansig50array==-500)=[];
        meansig50array(isnan(meansig50array))=-500;
        meansig50array(isinf(meansig50array))=-500;
        meansig50array(abs(meansig50array)>100)=-500;
        if ~isempty(meansig50array)
            meansig50array(all(meansig50array == 0,2),:)=[];
            stder=std(meansig50array(meansig50array~=-500))/sqrt(length(meansig50array(meansig50array~=-500)));
            alllength=length(meansig50array(meansig50array~=-500));
            %norm=jbtest(meansig50array(meansig50array~=-500));
            [h,p1]=ttest(meansig50array(meansig50array~=-500),0);
            meanaltogether=  mean(meansig50array(meansig50array~=-500));
            resultind((currents-1)*4+1:currents*4,i)=[meanaltogether; stder; p1; alllength];
        end
    end
end
%% result calc with currents
resultind=zeros(4,length(fieldnames(rsclass)));
for currents=1:8
    for i=1:length(fieldnames(rsclass))
        str=fieldnames(rsclass);
        str=char(str(i));
        meansig50array=rsclass.(str);
        if ~isempty(meansig50array)
            meansig50array=meansig50array(:,currents);
            meansig50array(meansig50array==-500)=[];
            meansig50array(isnan(meansig50array))=-500;
            meansig50array(isinf(meansig50array))=-500;
            meansig50array(abs(meansig50array)>100)=-500;
            if ~isempty(meansig50array)&&length(meansig50array)>1
                meansig50array(all(meansig50array == 0,2),:)=[];
                stder=std(meansig50array(meansig50array~=-500))/sqrt(length(meansig50array(meansig50array~=-500)));
                alllength=length(meansig50array(meansig50array~=-500));
                norm=jbtest(meansig50array(meansig50array~=-500));
                [h,p1]=ttest(meansig50array(meansig50array~=-500),0);
                meanaltogether=  mean(meansig50array(meansig50array~=-500));
                resultind((currents-1)*4+1:currents*4,i)=[meanaltogether; stder; p1; alllength];
            end
        end
    end
end

%% Plotting
figure
model_series=[];
model_error=[];
pval=[];
nsamp=[];

if size(resultind,1)/4>1
    
    resultind(:,all(resultind == 0,1))=[];
    for i=2:size(resultind,1)/4
        model_series=([model_series;resultind((4*(i-1))+1,:).*100;]);
        model_error = ([model_error; resultind((4*(i-1))+2,:).*100;]);%[0.143651201; 	0.122944593; 	0.197613506].*100;
        pval=[pval;resultind((4*(i-1))+3,:);];%[9.54*10^(-08),	1.93*10^(-06),	0.007309218];
        nsamp=[nsamp;resultind((4*(i-1))+4,:);];
    end
    
else
    resultind(:,all(resultind == 0,1))=[];
    model_series=resultind(1,:).*100;
    model_error = resultind(2,:).*100;%[0.143651201; 	0.122944593; 	0.197613506].*100;
    pval=resultind(3,:);%[9.54*10^(-08),	1.93*10^(-06),	0.007309218];
    nsamp=resultind(4,:);
end

% model_series=flipud(model_series');
% model_error=flipud(model_error');
% pval=fliplr(pval);
% nsamp=fliplr(nsamp);

model_series=(model_series');
model_error=(model_error');
pval=(pval);
nsamp=(nsamp);

% model_series=(model_series);
% model_error=(model_error);
% pval=(pval');
% nsamp=(nsamp');

cmap=colormap('gray');

b=bar(model_series);
for colorchange=1:length(b)
    b(length(b)-colorchange+1).FaceColor = cmap(floor(size(cmap,1)*colorchange/(length(b))),:);
    
end
if length(b)==1
    b(1).FaceColor=cmap(1,:);
end
hold on
% Find the number of groups and the number of bars in each group
ngroups = size(model_series, 1);
nbars = size(model_series, 2);
% Calculate the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
% Set the position of each error bar in the centre of the main bar
% Based on barweb.m by Bolu Ajiboye from MATLAB File Exchange
for i = 1:nbars
    %    Calculate center of each bar
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, model_series(:,i), model_error(:,i), 'k', 'linestyle', 'none');
end

hold off
%%For MATLAB 2019b or later releases
hold on
% Calculate the number of bars in each group
nbars = size(model_series, 2);
% Get the x coordinate of the bars
x = [];
for i = 1:nbars
    x = [x ; b(i).XEndPoints];
end
% Plot the errorbars
errorbar(x',model_series,model_error,'k','linestyle','none');
% plot sig stars
x2s=0.05;
x3s=0.075;
pvall01=pval<0.05;
pvall001=pval<0.01;
pvall01(pvall001)=0;
pvall0001=pval<0.001;
pvall01(pvall0001)=0;
pvall001(pvall0001)=0;
for i=1:size(pvall01,2)
    if pvall01(:,i).*x(:,i)==0
        continue
    end
    ypos=(model_series(i,:)+model_error(i,:)+0.15.*100);
    ypos(ypos<10)=10+15;
    xpos=x(:,i);
    
    scatter(pvall01(:,i).*xpos, pvall01(:,i)'.*ypos,'*','k')%delete apostrophe
    
end
for i=1:size(pvall001,2)
    if pvall001(:,i).*x(:,i)==0
        continue
    end
    ypos=(model_series(i,:)+model_error(i,:)+0.15.*100);
    ypos(ypos<10)=10+15;
    xpos=x(:,i);
    scatter(pvall001(:,i).*xpos-x2s, pvall001(:,i)'.*ypos,'*','k')%delete apostrophe-0.075
    scatter(pvall001(:,i).*xpos+x2s, pvall001(:,i)'.*ypos,'*','k')%delete apostrophe 0.075
    
    
end
for i=1:size(pvall0001,2)
    if pvall0001(:,i).*x(:,i)==0
        continue
    end
    xpos=x(:,i);
    ypos=(model_series(i,:)+model_error(i,:)+0.15.*100);
    ypos(ypos<10)=10+15;
    scatter(pvall0001(:,i).*xpos-x3s, pvall0001(:,i)'.*ypos,'*','k')%delete apostrophe
    scatter(pvall0001(:,i).*xpos+x3s, pvall0001(:,i)'.*ypos,'*','k')%delete apostrophe
    scatter(pvall0001(:,i).*xpos, pvall0001(:,i)'.*ypos,'*','k')%delete apostrophe
end
hold off
%xticklabels({'75%', '50%', '25%'})
%xticklabels({'R5', 'R6', 'R7', 'R8', 'R9'})
%xticklabels({'1', '2', '3', '4', '6', '8', '10'})
xticklabels(({'20^2\leqv<30^2','30^2\leqv<40^2','40^2\leqv<50^2','50^2\leqv<60^2','60^2\leqv<70^2','70^2\leqv<80^2','80^2\leqv<90^2','90^2\leqv<100^2','100^2\leqv<110^2','110^2\leqv<120^2','120^2\leqv<130^2','130^2\leqv<140^2','140^2\leqv'}))%xlabel('Ratio of current towards deep electrode')
xlabel('v=Voltage(mV)')
ylabel('Additivity Index (%)')
%title('Ratio 25% towards favoured electrode');
ymax=150;
ymin=-50;
ytotal=abs(ymin)+ymax;
ylim([ymin ymax])
%xlim([0.5 7.5])
%xlim([1.5 8.5])
%title(['Across: Comparison of current levels in eliciting spike rates at 50% @ shank 1 and shank 4']);
%title(['Laminar: Comparison of current levels in eliciting spike rates for all 50%/50% condition']);
title(['Electrode grouping based on estimated voltage at each electrode with 50%/50% condition']);

for i=1:size(nsamp,2)
    xpos=x(:,i);
    ypos=(model_series(i,:)+model_error(i,:)+0.07.*100);
    ypos(ypos<10)=10;
    th=text(xpos, ypos,num2str(nsamp(:,i)),'HorizontalAlignment', 'center');%delete apostrophe
end
if length(b)~=1
    leg=legend('1\muA','2\muA','3\muA','4\muA','6\muA','8\muA','10\muA');
    title(leg,'Current total')
end

% b.CData(7,:)=cmap(floor(size(cmap,1)*1/(7)),:);
% b.CData(6,:)=cmap(floor(size(cmap,1)*2/(7)),:);
% b.CData(5,:)=cmap(floor(size(cmap,1)*3/(7)),:);
% b.CData(4,:)=cmap(floor(size(cmap,1)*4/(7)),:);
% b.CData(3,:)=cmap(floor(size(cmap,1)*5/(7)),:);
% b.CData(2,:)=cmap(floor(size(cmap,1)*6/(7)),:);
% b.CData(1,:)=cmap(floor(size(cmap,1)*7/(7)),:);
b.FaceColor = 'flat';
endj=13;
for j=1:endj
b.CData(endj-j+1,:)=cmap(floor(size(cmap,1)*j/(endj)),:);
end


%% 75
%% four shank orig
%% four shank distance from stim electrodes
distpair=[];
stimchns=StimchnNoac_all.Rat_006;%change accordingly
AMP=[1 2 3 4 6 8 10].*10^-6;
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.25)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.75)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_007;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.25)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.75)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r7=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_008;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.25)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.75)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_009;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.25)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.75)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r9=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_010;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.25)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.75)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r10=distpair;

distpair=[];
stimchns=StimchnNoac_all.Rat_011;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents*0.25)./(distanceEs*10^-6*0.3*2*pi))+((currents*0.75)./(distanceEd*10^-6*0.3*2*pi)));
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r11=distpair;

%% new 75
distpair=[];
stimchns=StimchnNoac_r6;%change accordingly
AMP=[1 2 3 4 6 8 10].*10^-6;
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=((currents./(distanceEs*0.3*2*pi))+(currents./(distanceEd*0.3*2*pi))).*10^9;
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_r7;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=((currents./(distanceEs*0.3*2*pi))+(currents./(distanceEd*0.3*2*pi))).*10^9;
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r7=distpair;

distpair=[];
stimchns=StimchnNoac_r8;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=((currents./(distanceEs*0.3*2*pi))+(currents./(distanceEd*0.3*2*pi))).*10^9;
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_r9;%change accordingly
for currentsi=1:length(AMP)
    currents=AMP(currentsi);
    for i=1:size(stimchns,1)
        if stimchns(i,1)~=0
            if stimchns(i,1)<17&&stimchns(i,2)<17
                shank=1;
            elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
                shank=4;
                stimchns(i,:)=stimchns(i,:)-16;
            elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
                shank=2;
                stimchns(i,:)=stimchns(i,:)-32;
            elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
                shank=3;
                stimchns(i,:)=stimchns(i,:)-48;
            end
            distanceEd=ones(16,4).*0.00000000001;
            distanceEs=ones(16,4).*0.00000000001;
            for j=1:stimchns(i,1)-1
                distanceEd(stimchns(i,1)-j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=400;
                elseif shank==3
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=200;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=400;
                elseif shank==4
                    distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=200;
                    distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=400;
                    distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),1)=600;
                elseif shank==1
                    distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),2)=200;
                    distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),3)=400;
                    distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEd(stimchns(i,1),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,1)
                distanceEd(stimchns(i,1)+j,shank)=50*j;
                if shank==2
                    distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            for j=1:stimchns(i,2)-1
                distanceEs(stimchns(i,2)-j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=400;
                elseif shank==3
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=200;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=400;
                elseif shank==4
                    distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=200;
                    distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=400;
                    distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),1)=600;
                elseif shank==1
                    distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),2)=200;
                    distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),3)=400;
                    distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                    distanceEs(stimchns(i,2),4)=600;
                end
            end
            
            for j=1:16-stimchns(i,2)
                distanceEs(stimchns(i,2)+j,shank)=50*j;
                if shank==2
                    distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
                elseif shank==3
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
                elseif shank==4
                    distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
                elseif shank==1
                    distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                    distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
                end
            end
            distanceAdditive=(((currents.*0.25)./(distanceEs*0.3*2*pi))+((currents.*0.75)./(distanceEd*0.3*2*pi))).*10^9;
            distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
            check=['T' num2str(i) num2str(currents/(10^-6))];
            distpair.(check)=distanceAdditive;
        end
    end
end
distpair_r9=distpair;


%% single shank
distpair=[];

stimchns=StimchnNoac_r5;%change accordingly
for i=1:size(stimchns,1)
    distanceEd=ones(32,1).*0.00000000001;
    distanceEs=ones(32,1).*0.00000000001;
    if stimchns(i,1)~=0
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,1)=50*j;
        end
        
        for j=1:32-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,1)=50*j;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,1)=50*j;
        end
        
        for j=1:32-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,1)=50*j;
        end
        distanceAdditive=((0.25./(distanceEs))+(0.75./(distanceEd))).*10e3;
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end

distpair_r5=distpair;
%% Put together

% r5input=currentvariation75_r5;
% r5input(all(r5input == 0,2),:)=-500;
r6input=currentvariation75_r6;
r6input(all(r6input == 0,2),:)=-500;
r7input=currentvariation50_r7;
r7input(all(r7input == 0,2),:)=-500;
r8input=currentvariation75_r8;
r8input(all(r8input == 0,2),:)=-500;
r9input=currentvariation75_r9;
r9input(all(r8input == 0,2),:)=-500;
% D=[];
% stimchns=StimchnNoac_r5;%change accordingly
% if (32-max(stimchns,[],'all'))>(min(stimchns(stimchns~=0),[],'all'))
%     [r,c]=find(min(stimchns(stimchns~=0),[],'all')==stimchns);
% else
%     [r,c]=find(max(stimchns,[],'all')==stimchns);
% end
% distinterest=distpair_r5.(['T' num2str(r)]);
endloop=150;
step=50;
rsclass=[];
for i=0:step:endloop
    check=['T' num2str(i)];
    rsclass.(check)=[];
end
%single electrode
% D=[];
% for stimchncounter=1:size(stimchns,1)
%     check=['T' num2str(stimchncounter)];
%     try
%         distinterest=distpair_r5.(check);
%
%         for i=0:step:endloop
%             check=['T' num2str(i)];
%
%             if i==endloop
%                 D.(check)=find((distinterest>=((i-step))));
%             elseif i~=endloop
%                 D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
%             end
%             if D.(check)~=0
%                 row=D.(check)+(stimchncounter-1)*32;
%                 rsclass.(check)=[rsclass.(check); r5input(row,:)];
%             end
%         end
%
%     catch
%         %not laminar or suitable sep dist
%     end
%     D=[];
% end


%% four shank
%6
D=[];
stimchns=StimchnNoac_r6;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r6.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r6input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%7
D=[];
stimchns=StimchnNoac_r7;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r7.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r7input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%8
D=[];
stimchns=StimchnNoac_r8;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r8.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r8input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%9
D=[];
stimchns=StimchnNoac_r9;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r9.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r9input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%% 25
distpair=[];
stimchns=StimchnNoac_r6;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17&&stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
            shank=4;
            stimchns(i,:)=stimchns(i,:)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
            shank=2;
            stimchns(i,:)=stimchns(i,:)-32;
        elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
            shank=3;
            stimchns(i,:)=stimchns(i,:)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_r7;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17&&stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
            shank=4;
            stimchns(i,:)=stimchns(i,:)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
            shank=2;
            stimchns(i,:)=stimchns(i,:)-32;
        elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
            shank=3;
            stimchns(i,:)=stimchns(i,:)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r7=distpair;

distpair=[];
stimchns=StimchnNoac_r8;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17&&stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
            shank=4;
            stimchns(i,:)=stimchns(i,:)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
            shank=2;
            stimchns(i,:)=stimchns(i,:)-32;
        elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
            shank=3;
            stimchns(i,:)=stimchns(i,:)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_r9;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17&&stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,1)<33&&stimchns(i,2)<33 && stimchns(i,1)>16&&stimchns(i,2)>16
            shank=4;
            stimchns(i,:)=stimchns(i,:)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,1)>32&&stimchns(i,2)>32
            shank=2;
            stimchns(i,:)=stimchns(i,:)-32;
        elseif stimchns(i,1)<65&&stimchns(i,2)<65 && stimchns(i,1)>48&&stimchns(i,2)>48
            shank=3;
            stimchns(i,:)=stimchns(i,:)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r9=distpair;






%% Across
%% Across
Stimchnall=Stimchnall_r6;%change accordingly
for chnpair=1:length(Stimchnall)
    if Stimchnall(chnpair,1)<17&&Stimchnall(chnpair,2)<17
        Stimchnall(chnpair,:)=[0,0];
    elseif Stimchnall(chnpair,1)<33&&Stimchnall(chnpair,2)<33 && Stimchnall(chnpair,1)>16&&Stimchnall(chnpair,2)>16
        Stimchnall(chnpair,:)=[0,0];
    elseif Stimchnall(chnpair,1)<49&&Stimchnall(chnpair,2)<49 && Stimchnall(chnpair,1)>32&&Stimchnall(chnpair,2)>32
        Stimchnall(chnpair,:)=[0,0];
    elseif Stimchnall(chnpair,1)<65&&Stimchnall(chnpair,2)<65 && Stimchnall(chnpair,1)>48&&Stimchnall(chnpair,2)>48
        Stimchnall(chnpair,:)=[0,0];
    else
        %across
        
    end
end
StimchnNoac_r6=Stimchnall;%change accordingly

%% four shank calcualte distance
distpair=[];
stimchns=StimchnNoac_r6;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1./(distanceEs))+(1./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_r8;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1./(distanceEs))+(1./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_r9;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1./(distanceEs))+(1./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r9=distpair;

%% put together
r6input=currentvariation50_r6;
r6input(all(r6input == 0,2),:)=-500;
r8input=currentvariation50_r8;
r8input(all(r8input == 0,2),:)=-500;
r9input=currentvariation50_r9;
r9input(all(r8input == 0,2),:)=-500;

endloop=150;
step=50;
rsclass=[];
for i=0:step:endloop
    check=['T' num2str(i)];
    rsclass.(check)=[];
end

%% Group ratios
%6
D=[];
stimchns=StimchnNoac_r6;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r6.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r6input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%8
D=[];
stimchns=StimchnNoac_r8;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r8.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r8input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%9
D=[];
stimchns=StimchnNoac_r9;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r9.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r9input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%% result calc with currents
resultind=zeros(4,length(fieldnames(rsclass)));
for currents=1:8
    for i=1:length(fieldnames(rsclass))
        str=fieldnames(rsclass);
        str=char(str(i));
        meansig50array=rsclass.(str);
        if ~isempty(meansig50array)
            meansig50array=meansig50array(:,currents);
            meansig50array(meansig50array==-500)=[];
            meansig50array(isnan(meansig50array))=-500;
            meansig50array(isinf(meansig50array))=-500;
            meansig50array(abs(meansig50array)>100)=-500;
            if ~isempty(meansig50array)&&length(meansig50array)>1
                meansig50array(all(meansig50array == 0,2),:)=[];
                stder=std(meansig50array(meansig50array~=-500))/sqrt(length(meansig50array(meansig50array~=-500)));
                alllength=length(meansig50array(meansig50array~=-500));
                norm=jbtest(meansig50array(meansig50array~=-500));
                [h,p1]=ttest(meansig50array(meansig50array~=-500),0);
                meanaltogether=  mean(meansig50array(meansig50array~=-500));
                resultind((currents-1)*4+1:currents*4,i)=[meanaltogether; stder; p1; alllength];
            end
        end
    end
end


%% four shank calcualte distance 75
distpair=[];
stimchns=StimchnNoac_r6;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((0.5./(distanceEs))+(1.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_r8;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((0.5./(distanceEs))+(1.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_r9;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((0.5./(distanceEs))+(1.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r9=distpair;


%% put together
r6input=currentvariation75_r6;
r6input(all(r6input == 0,2),:)=-500;
r8input=currentvariation75_r8;
r8input(all(r8input == 0,2),:)=-500;
r9input=currentvariation75_r9;
r9input(all(r8input == 0,2),:)=-500;

endloop=150;
step=50;
rsclass=[];
for i=0:step:endloop
    check=['T' num2str(i)];
    rsclass.(check)=[];
end

%% Group ratios
%6
D=[];
stimchns=StimchnNoac_r6;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r6.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r6input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%8
D=[];
stimchns=StimchnNoac_r8;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r8.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r8input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%9
D=[];
stimchns=StimchnNoac_r9;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r9.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r9input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%% four shank calcualte distance 25
distpair=[];
stimchns=StimchnNoac_r6;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r6=distpair;

distpair=[];
stimchns=StimchnNoac_r8;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r8=distpair;

distpair=[];
stimchns=StimchnNoac_r9;%change accordingly
for i=1:size(stimchns,1)
    if stimchns(i,1)~=0
        if stimchns(i,1)<17
            shank=1;
        elseif stimchns(i,1)<33 && stimchns(i,1)>16
            shank=4;
            stimchns(i,1)=stimchns(i,1)-16;
        elseif stimchns(i,1)<49 && stimchns(i,1)>32
            shank=2;
            stimchns(i,1)=stimchns(i,1)-32;
        elseif stimchns(i,1)<65 && stimchns(i,1)>48
            shank=3;
            stimchns(i,1)=stimchns(i,1)-48;
        end
        distanceEd=ones(16,4).*0.00000000001;
        distanceEs=ones(16,4).*0.00000000001;
        for j=1:stimchns(i,1)-1
            distanceEd(stimchns(i,1)-j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=400;
            elseif shank==3
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=200;
                distanceEd(stimchns(i,1)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=400;
            elseif shank==4
                distanceEd(stimchns(i,1)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=200;
                distanceEd(stimchns(i,1)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=400;
                distanceEd(stimchns(i,1)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),1)=600;
            elseif shank==1
                distanceEd(stimchns(i,1)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1),2)=200;
                distanceEd(stimchns(i,1)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1),3)=400;
                distanceEd(stimchns(i,1)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEd(stimchns(i,1),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,1)
            distanceEd(stimchns(i,1)+j,shank)=50*j;
            if shank==2
                distanceEd(stimchns(i,1)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEd(stimchns(i,1)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEd(stimchns(i,1)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEd(stimchns(i,1)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        if stimchns(i,2)<17
            shank=1;
        elseif stimchns(i,2)<33 &&stimchns(i,2)>16
            shank=4;
            stimchns(i,2)=stimchns(i,2)-16;
        elseif stimchns(i,1)<49&&stimchns(i,2)<49 && stimchns(i,2)>32
            shank=2;
            stimchns(i,2)=stimchns(i,2)-32;
        elseif stimchns(i,2)<65 && stimchns(i,2)>48
            shank=3;
            stimchns(i,2)=stimchns(i,2)-48;
        end
        for j=1:stimchns(i,2)-1
            distanceEs(stimchns(i,2)-j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)-j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=400;
            elseif shank==3
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=200;
                distanceEs(stimchns(i,2)-j,1)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=400;
            elseif shank==4
                distanceEs(stimchns(i,2)-j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=200;
                distanceEs(stimchns(i,2)-j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=400;
                distanceEs(stimchns(i,2)-j,1)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),1)=600;
            elseif shank==1
                distanceEs(stimchns(i,2)-j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2),2)=200;
                distanceEs(stimchns(i,2)-j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2),3)=400;
                distanceEs(stimchns(i,2)-j,4)=sqrt((600^2)+(50*j)^2);
                distanceEs(stimchns(i,2),4)=600;
            end
        end
        
        for j=1:16-stimchns(i,2)
            distanceEs(stimchns(i,2)+j,shank)=50*j;
            if shank==2
                distanceEs(stimchns(i,2)+j,1)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((400^2)+(50*j)^2);
            elseif shank==3
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((400^2)+(50*j)^2);
            elseif shank==4
                distanceEs(stimchns(i,2)+j,3)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,2)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,1)=sqrt((600^2)+(50*j)^2);
            elseif shank==1
                distanceEs(stimchns(i,2)+j,2)=sqrt((200^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,3)=sqrt((400^2)+(50*j)^2);
                distanceEs(stimchns(i,2)+j,4)=sqrt((600^2)+(50*j)^2);
            end
        end
        distanceAdditive=((1.5./(distanceEs))+(0.5./(distanceEd))).*10e3;
        distanceAdditive=[distanceAdditive(:,1),distanceAdditive(:,4),distanceAdditive(:,2),distanceAdditive(:,3)];%change to electrode mapping
        check=['T' num2str(i)];
        distpair.(check)=distanceAdditive;
    end
end
distpair_r9=distpair;


%% put together
r6input=currentvariation25_r6;
r6input(all(r6input == 0,2),:)=-500;
r8input=currentvariation25_r8;
r8input(all(r8input == 0,2),:)=-500;
r9input=currentvariation25_r9;
r9input(all(r8input == 0,2),:)=-500;

%% Group ratios 25
%6
D=[];
stimchns=StimchnNoac_r6;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r6.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r6input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%8
D=[];
stimchns=StimchnNoac_r8;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r8.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r8input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%9
D=[];
stimchns=StimchnNoac_r9;%change accordingly
for stimchncounter=1:size(stimchns,1)
    check=['T' num2str(stimchncounter)];
    try
        distinterest=distpair_r9.(check);%change accordingly
        
        for i=0:step:endloop
            check=['T' num2str(i)];
            if i==endloop
                D.(check)=find((distinterest>=((i-step))));
            elseif i~=endloop
                D.(check)=find((distinterest<i)&(distinterest>=((i-step))));
            end
            if D.(check)~=0
                row=D.(check)+(stimchncounter-1)*64;
                rsclass.(check)=[rsclass.(check); r9input(row,:)];%change accordingly
            end
        end
        
    catch
        %not laminar or suitable sep dist
    end
    D=[];
end

%% result calc with currents
resultind=zeros(4,length(fieldnames(rsclass)));
for currents=1:8
    for i=1:length(fieldnames(rsclass))
        str=fieldnames(rsclass);
        str=char(str(i));
        meansig50array=rsclass.(str);
        if ~isempty(meansig50array)
            meansig50array=meansig50array(:,currents);
            meansig50array(meansig50array==-500)=[];
            meansig50array(isnan(meansig50array))=-500;
            meansig50array(isinf(meansig50array))=-500;
            meansig50array(abs(meansig50array)>100)=-500;
            if ~isempty(meansig50array)&&length(meansig50array)>1
                meansig50array(all(meansig50array == 0,2),:)=[];
                stder=std(meansig50array(meansig50array~=-500))/sqrt(length(meansig50array(meansig50array~=-500)));
                alllength=length(meansig50array(meansig50array~=-500));
                norm=jbtest(meansig50array(meansig50array~=-500));
                [h,p1]=ttest(meansig50array(meansig50array~=-500),0);
                meanaltogether=  mean(meansig50array(meansig50array~=-500));
                resultind((currents-1)*4+1:currents*4,i)=[meanaltogether; stder; p1; alllength];
            end
        end
    end
end





%% 7 electrode distance matrix laminar
Stimchnall=Stimchnall_r9;%change accordingly
for chnpair=1:length(Stimchnall)
    if Stimchnall(chnpair,1)<17&&Stimchnall(chnpair,2)<17
    elseif Stimchnall(chnpair,1)<33&&Stimchnall(chnpair,2)<33 && Stimchnall(chnpair,1)>16&&Stimchnall(chnpair,2)>16
    elseif Stimchnall(chnpair,1)<49&&Stimchnall(chnpair,2)<49 && Stimchnall(chnpair,1)>32&&Stimchnall(chnpair,2)>32
    elseif Stimchnall(chnpair,1)<65&&Stimchnall(chnpair,2)<65 && Stimchnall(chnpair,1)>48&&Stimchnall(chnpair,2)>48
    else
        %across
        Stimchnall(chnpair,:)=[0,0];
    end
end
StimchnNoac_r9=Stimchnall;
%% assuming 5 electrode separation
StimchnNoac_r5(StimchnNoac_r5(:,1)-StimchnNoac_r5(:,2)~=-10,:)=0;
StimchnNoac_r6(StimchnNoac_r6(:,1)-StimchnNoac_r6(:,2)~=-10,:)=0;
StimchnNoac_r7(StimchnNoac_r7(:,1)-StimchnNoac_r7(:,2)~=-10,:)=0;
StimchnNoac_r8(StimchnNoac_r8(:,1)-StimchnNoac_r8(:,2)~=-10,:)=0;
StimchnNoac_r9(StimchnNoac_r9(:,1)-StimchnNoac_r9(:,2)~=-10,:)=0;

%% population level response
%rat sep
rsclass=[];
rsclass.T1=r5input;
rsclass.T2=r6input;
rsclass.T3=r7input;
rsclass.T4=r8input;
rsclass.T5=r9input;
rsclass.T1=[rsclass.T1,-500.*ones(length(rsclass.T1),1)];



%% total pop
rsclass=[];
rsclass.T1=[[r5input,-500.*ones(length(r5input),1)]; r6input;r7input;r8input;r9input;];
