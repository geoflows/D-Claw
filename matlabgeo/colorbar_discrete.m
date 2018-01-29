function hcbar = colorbar_discrete(flow_colormap,axhandle)
%
colormap(flow_colormap(:,2:4));
tvalues = [flow_colormap(:,1)];
[ncolors,junk] = size(tvalues);
ticks = linspace(0,1,ncolors+1);
dtick = ticks(2)-ticks(1);
cticks = ticks(1:end-1)+0.5*dtick;
strtvalues = cellstr(num2str(tvalues));
celllabels = cellstr(strtvalues);

celllabels(1) = cellstr(strcat('< ',strtvalues(1)));
for jl = 2:ncolors-1
    celllabels(jl) = cellstr(strcat(strtvalues(jl-1),'-',strtvalues(jl)));
end
celllabels(ncolors) = cellstr(strcat('> ',strtvalues(ncolors-1)));
hcbar=colorbar('peer',axhandle,'Ticks',cticks,'TickLabels',celllabels);
