function [template,r2t] = loadSpikeTemplate
filepath = pwd;
tmp = dir([filepath '\*.sp.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    t = load(tmp,'template');
    template = t.template;
    r = load(tmp,'r2t');
    r2t = r.r2t;
end
end