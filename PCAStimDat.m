%PCA stim data
%Pick electrode or use centre and separate electrode

%E:\DATA\Rat_012\S1E3_001_201127_100004
%Set up data and get ptsh
sp=loadSpikes;
order=Depth(1);
spsort=sp(order);
%%
ID=26;
clear rate
BIN = [-100 100]; 
trig = loadTrig(0);
TP = loadTrialParams;
tID = cell2mat(TP(cell2mat(TP(:,2)) == ID,1));
theseTrig = trig(tID)./30;
nT=length(theseTrig);
MAX = 400;
count=0;
for chns=3:2:7
spchn=spsort{chns};
xdata = [];
ydata = [];
for tr = 1:nT
    theseSp = (spchn(spchn > theseTrig(tr)+BIN(1) & spchn < theseTrig(tr)+BIN(2)) - theseTrig(tr));
    for i = 1:length(theseSp)
        xdata = [xdata, (theseSp(i)+ abs(BIN(1)))]; %#ok<*AGROW>
        ydata = [ydata, tr*(MAX/nT)];
    end
end
Z = hist(xdata,0:abs(BIN(1))+abs(BIN(2))); %#ok<HIST>
count=count+1;
rate(count,:) = (1000/nT)*Z;%conv(Z,window);
end



%%

dat=rate';
%dat=dat(:,1:2);
%doing this manually so I understand what is going on
% varx = mean(dat(:, 1).^2) - mean(dat(:, 1))^2;		% <xx> - <x><x>
% vary = mean(dat(:, 2).^2) - mean(dat(:, 2))^2;		% <yy> - <y><y>
% covarxy = mean(dat(:, 1) .* dat(:, 2)) - mean(dat(:, 1)) * mean(dat(:, 2));		% <xy> - <x><y>
% 
% Covar = [varx 		covarxy
% 		 covarxy 	vary];
Covar=cov(dat);
[EigVec, EigVal] = eigs(Covar);


figure(1)
clf
subplot(1, 2, 1)
hold on
plot(dat(:, 1), dat(:, 2), '.');
axis([-3 3 -3 3])
axis square
xlabel('x')
ylabel('y')
PlotVector(EigVec(:, 1), 'r');
PlotVector(EigVec(:, 2), 'g');
hold off

Proj1 = dat * EigVec(:, 1);
Proj2 = dat * EigVec(:, 2);
Proj3= dat * EigVec(:,3);
figure; scatter3(Proj1(1:100),Proj2(1:100),Proj3(1:100),'b')
hold on
scatter3(Proj1(101:110),Proj2(101:110),Proj3(101:110),'k')
scatter3(Proj1(111:200),Proj2(111:200),Proj3(111:200),'r')


