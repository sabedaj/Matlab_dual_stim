function heatmap_array=plotheatmapdistdata(data, stimchn)
%% sorting heatmap from sigmoiddata
stimVA=1; %visual area being stimulated
VA=1;
columns=[1,3,4,2];
chnnrange=1:64;
order=Depth(1);
ordershapearray=reshape(order,16,8);
ordershapearray=ordershapearray(:,columns);

    stimchnarray=[1:16;33:48;49:64;17:32]';
        schn=regexp((stimchn),'\d*','Match');
        schn=str2double(schn{1});
        schn=schn-64;
heatmap_array=nan(31,7);
for recchn=chnnrange(1):chnnrange(end)
    [rrecchn,crecchn]=find(ordershapearray==recchn);
        %if c~=crecchn %stim shank only
        if VA==stimVA
            [r,c]=find(schn==stimchnarray);
            if isempty(r)
                continue
            end
            rplus=16-r;
            cplus=4-c;
        end
            heatmap_array(rrecchn+rplus(1),crecchn+cplus(1))=data(recchn);%iterates forwards
end