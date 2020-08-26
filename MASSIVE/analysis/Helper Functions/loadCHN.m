filepath = pwd;
tmp = dir([filepath '\*exp_datafile_*.mat']);
if ~isempty(tmp)
    tmp = [filepath '\' tmp.name];
    CHN = load(tmp,'CHN');
    CHN = CHN.CHN;
end