
% geo_plot_topo.m plots the topography
%    topo = aux(:,:,1)
%    topo2 =topo




if ~exist('topo_colormap')
  disp('*** You must define topo_colormap, e.g. to one of the maps')
  disp('    defined in geo_setzcolormaps.m')
  break
  end

topo2colors = mapzcolors(topo2,topo_colormap);


cw=surf(X,Y,topo2,topo2colors);

if (PlotGrid(1)==1)
    set(cw,'FaceColor','flat','EdgeColor',[0 0 0]);
else
    set(cw,'FaceColor','flat','EdgeColor','none');
end


ylabel('Latitude','Fontsize',12)
xlabel('Longitude','Fontsize',12)
