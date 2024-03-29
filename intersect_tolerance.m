function sp=intersect_tolerance(sp)
tolerance=0.5;%in ms

for headstage=1:4 % go through headstages
    for elect_overall=1:32 % compare these electrodes against ...
        sp1=sp{elect_overall+(headstage-1)}(:,1);
        for spiketocheck=1:size(sp1,1) % spikes from E1^
            spiketime1=sp1(spiketocheck);
            spiketimemax=spiketime1+tolerance;
            spiketimemin=spiketime1-tolerance;
            overallnumspikesintersect=0;
            for elect_check=1:32 % ... these electrodes
                if elect_check==elect_overall
                    continue
                end
                sp2=sp{elect_check+(headstage-1)}(:,1);
                spikeintersect=sp2(sp2>spiketimemin & sp2<spiketimemax);
                if ~isempty(spikeintersect)
                overallnumspikesintersect=overallnumspikesintersect+1;
                end
                if overallnumspikesintersect>16
                    for e=1:32
                        sp3=sp{e+(headstage-1)};
                        sp3(sp3(:,1)>spiketimemin & sp3(:,1)<spiketimemax,:)=[];
                        sp{e+(headstage-1)}=sp3;
                    end
                   break
                end
            end
        end
 
    end
end

end