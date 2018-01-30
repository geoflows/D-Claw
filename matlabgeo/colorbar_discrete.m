function hcbar = colorbar_discrete(flow_colormap,axhandle)
%
colormap(flow_colormap(:,2:4));
tvalues = [flow_colormap(:,1)];
[ncolors,junk] = size(tvalues);
ticks = linspace(0,1,ncolors+1);
dtick = ticks(2)-ticks(1);
cticks = ticks(1:end-1)+0.5*dtick;
%strtvalues = cellstr(num2str(tvalues));
%celllabels = cellstr(strtvalues);
celllabels = {};

celllabels{1} = ['< ',num2str(tvalues(1))];
for jl = 2:ncolors-1
    celllabels{jl} = [num2str(tvalues(jl-1)),' to ',num2str(tvalues(jl))];
end
celllabels{ncolors} = ['> ',num2str(tvalues(ncolors-1))];
hcbar=colorbar('peer',axhandle,'Ticks',cticks,'TickLabels',celllabels);
