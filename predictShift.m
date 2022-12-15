AMP=[1 2 3 4 6 8 10];
avgpcurr=zeros(5,length(AMP));
stdavg=zeros(5,length(AMP));
 s = RandStream('mlfg6331_64','Seed','shuffle'); 
for sepdist=5:2:9
    sepcheck=['sep' num2str(sepdist)];
    CsplitCentroid= centroidpos_all.(sepcheck);
    clear pairavg stdpairavg cavgspk
    for current=1:length(AMP)
        currcheck=['C' num2str(AMP(current))];
        
        % Define number of columns to average
        AVG_COLS = 4;
        % Dimension over which to average
        DIM = 2; % Columns
        % Use filter to calculate the moving average across EVERY combination of columns
        pairavg.(currcheck) = filter(ones(1,AVG_COLS)/AVG_COLS,1,CsplitCentroid.(currcheck),[],DIM);
        % Grab only the column averages that were actually wanted
        pairavg.(currcheck) = pairavg.(currcheck)(:,AVG_COLS:AVG_COLS:end);
        
        %downsample
%         lengthneeded=length(pairavg.C1); 
%         samplesneeded = datasample(s,1:size(pairavg.(currcheck),2),lengthneeded,'Replace',false);
%         
%         pairavg.(currcheck) = pairavg.(currcheck)(:,samplesneeded);
        
        avgpcurr(:,current)=mean(pairavg.(currcheck),2);%mean(pairavg(:,current:length(AMP):counternumberofpairs*length(AMP)),2);
        stdavg(:,current)=std(pairavg.(currcheck),[],2)./sqrt(size(pairavg.(currcheck),2));
    end
    avgshift.(sepcheck)=avgpcurr;
    stdshift.(sepcheck)=stdavg;
end
%is the change in centroid significantly less for 5 sep than 9 sep
sep5shiftdual=centroidpos_all.sep5.C6(2:4,:);%avgshift.sep5(2:4,:); 
sep9shiftdual=centroidpos_all.sep9.C6(2:4,:);%avgshift.sep9(2:4,:);
[p59_centroid,h,stats] = ranksum(sep5shiftdual(:), sep9shiftdual(:),'tail','left');

%is the change in centroid significantly less for 5 sep than 7 sep
sep7shiftdual=centroidpos_all.sep7.C6(2:4,:);%avgshift.sep7(2:4,:);
[p57_centroid,h,stats] = ranksum(sep5shiftdual(:), sep7shiftdual(:),'tail','left');

%is the change in centroid significantly less for 5 sep than 7 sep
[p79_centroid,h,stats] = ranksum(sep7shiftdual(:), sep9shiftdual(:),'tail','left');

avgshiftarray=[diff(avgshift.sep5); diff(avgshift.sep7); diff(avgshift.sep9)];
dualvsingle=[ones(1,length(AMP)).*2; ones(2,length(AMP)); ones(1,length(AMP)).*2; ones(1,length(AMP)).*2; ones(2,length(AMP)); ones(1,length(AMP)).*2; ones(1,length(AMP)).*2; ones(2,length(AMP)); ones(1,length(AMP)).*2;];
sepdistarray=[ones(4,length(AMP)).*300; ones(4,length(AMP)).*400; ones(4,length(AMP)).*500];
currentarray=repelem(AMP,12,1);

%current and sepdist
fun = @(b,X) (X(:,1).*b(1)+X(:,2).*b(2));
b0=[1 1];
mdl1 = fitnlm([sepdistarray(:) currentarray(:)], avgshiftarray(:),fun,b0);
MSE1=mdl1.MSE;
rsquared1=mdl1.Rsquared.Adjusted;
AIC1=mdl1.ModelCriterion.AIC;

%only sepdist
fun = @(b,X) (X(:,1).*b(1));
b0=[1];
mdl = fitnlm([sepdistarray(:)], avgshiftarray(:),fun,b0);

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;

%sepdist and trial
fun = @(b,X) (X(:,1).*b(1)+X(:,2)*b(2));
b0=[1 1];
mdl = fitnlm([sepdistarray(:) dualvsingle(:)], avgshiftarray(:),fun,b0);

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;

%sepdist and trial and current
fun = @(b,X) (X(:,1).*b(1)+X(:,2)*b(2)+X(:,3)*b(3));
b0=[1 1 1];
mdl = fitnlm([sepdistarray(:) dualvsingle(:) currentarray(:)], avgshiftarray(:),fun,b0);

MSE=mdl.MSE;
rsquared=mdl.Rsquared.Adjusted;
AIC=mdl.ModelCriterion.AIC;
     