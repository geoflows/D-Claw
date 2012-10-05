
mcolors=60;
zcbar(:,1)= linspace(-TsAmp,TsAmp,mcolors)';
zbarcolors=mapzcolors(zcbar,flow_colormap);
zcmap=zeros(mcolors,3);
for i=1:mcolors
    zcmap(i,:)=zbarcolors(i,1,:);
end
colormap(zcmap);
caxis([-TsAmp TsAmp]);
colorbar('peer',gca);