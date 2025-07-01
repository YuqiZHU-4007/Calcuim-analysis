function repmap=repclrmap(clrmap,gIX)
 repmap=zeros(size(gIX,1),size(clrmap,2));
 for ii=1:length(gIX)
     ind=gIX(ii);
     repmap(ii,:)=clrmap(ind,:);
 end
end