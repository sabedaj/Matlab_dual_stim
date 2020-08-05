 function avgStruct=meanstruct(avgStruct, ID, IDstruct)
check=['T', num2str(ID)];
avgStruct=[avgStruct, mean(IDstruct.(check),2)];
 end