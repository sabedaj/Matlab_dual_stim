function depth = Depth(varargin)
filepath = pwd;
NN_cutoff = datetime('06-Nov-2018 00:00:00');
Flex_cutoff = datetime('14-Jul-2020 00:00:00');
fileinfo = dir([filepath filesep 'info.rhs']);
nChn=32;
if datetime(fileinfo.date) < NN_cutoff
    NN_order = 0;
else
    NN_order = 1;
end
if datetime(fileinfo.date) > datetime('28-Nov-2018 00:00:00')
    NN_I = 1;
else 
    NN_I = 0;
end
if ~(NN_order)
    p = 1;
else
    p = 2;
end
if (NN_I)
    p = 3;
end
if datetime(fileinfo.date) > Flex_cutoff
    p = 5;
end
if nargin==1
    E_Mapnumber=cell2mat(varargin(1));
    if E_Mapnumber>0
        nChn=64;
        p=6;
    else
        nChn=32;
        p=5;
    end
end
E_MAP = ProbeMAP;
depth = zeros(1,nChn);
for n = 2:nChn+1
    str = E_MAP{n,p};
    str = str(4:5);
    str = str2double(str);
    depth(n-1) = str;
end
depth = depth';
if ~isempty(depth == 0)
    depth = depth + 1;
end
end