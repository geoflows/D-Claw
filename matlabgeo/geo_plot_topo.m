
% geo_plot_topo.m plots the topography
%    topo = aux(:,:,1)
%    topo2 = augmented version of size (mx+1) by (my+1) for call to surf.



% mask points in topo2 covered by finer grids:
%erasegrid;
%topo2masked = topo2 .* covered_ind;


if ~exist('topo_colormap')
  disp('*** You must define topo_colormap, e.g. to one of the maps')
  disp('    defined in geo_setzcolormaps.m')
  break
  end

topo2colors = mapzcolors(topo2,topo_colormap);

% set color to NaN in cells covered by finer grids:
set_covered_ind
topo2colors(:,:,1) = topo2colors(:,:,1) .* covered_ind;

%cw=surf(X,Y,topo2,topo2colors);
cw=surf(X,Y,topo2.*covered_ind,topo2colors);

if (PlotGrid(level)==1)
    set(cw,'FaceColor','interp','EdgeColor',[0 0 0]);
else
    set(cw,'FaceColor','interp','EdgeColor','none');
end
set(cw,'FaceLighting','flat','SpecularColorReflectance',0.0,'SpecularStrength',1)

if (PlotGridEdges(level)==1)
    l1=line(xedge,0*xedge+yedge(1),topo2(:,1)+1000,'Color','k','LineWidth',1);
    l2=line(xedge,0*xedge+yedge(end),topo2(:,end)+1000,'Color','k','LineWidth',1);
    l3=line(0*yedge+xedge(1),yedge,topo2(1,:)+1000,'Color','k','LineWidth',1);
    l4=line(0*yedge+xedge(end),yedge,topo2(end,:)+1000,'Color','k','LineWidth',1);
    
end

%ylabel('Latitude','Fontsize',12)
%xlabel('Longitude','Fontsize',12)
