 function stdStruct=stderrorstruct(stdStruct, ID, IDstruct)
check=['T', num2str(ID)];
N=length(IDstruct.(check));
stdStruct=[stdStruct, std(IDstruct.(check),[],2)/sqrt(N)];
 end