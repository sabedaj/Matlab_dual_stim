function [ordershapearrayV1,ordershapearrayV2]=sortCJelectrodeOrder_V1V2

order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearrayV1=ordershapearray(:,5:8);
ordershapearrayV2=ordershapearray(:,1:4);

%sorting dataitself
% x=fieldnames(Spike_trialstruct);
% startchn = str2double(regexp(x{1},'\d*','Match'));
% endchn = str2double(regexp(x{end},'\d*','Match'));
% V2array=nan(16,4);
% V1array=nan(16,4);
% for recchn=startchn:endchn
%  chnname=['Chn_' num2str(recchn)];
% 
%  if recchn<65
%      [rrecchn,crecchn]=find(ordershapearrayV2==recchn);
%     V2array(rrecchn,crecchn)=Spike_trialstruct.(chnname);
%  else
%     [rrecchn,crecchn]=find(ordershapearrayV1==recchn);
%     V1array(rrecchn,crecchn)=Spike_trialstruct.(chnname);
%  end
% end

end