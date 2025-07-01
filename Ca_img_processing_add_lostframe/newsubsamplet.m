function list=newsubsamplet(nsamples,zi,process)
if isempty(process.T_lessframe)
    list=(1:nsamples);
else
    list=(1:nsamples+length(find(process.T_lessframe<=nsamples)));
    for i=1:length(process.T_lessframe)
        ind=find(process.T_lessframe(i)==(1:nsamples+length(find(process.T_lessframe<=nsamples))));
        if ~isempty(ind)
                list(find(list==process.T_lessframe(i)))=0;
        end
    end
end
list(list==0)=[];
list=list(1:nsamples);