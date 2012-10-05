function [fid]= plottopotype(fname,topotype)

%plottopotype3(fname,topotype). 
%Plot the file fname using standard GeoClaw plotting features

Grid=topofile2Grid(fname,topotype);
 
PlotType = 11;     % = 11 for colored surface plot
                    % = 12 for contour plot

PlotFlow = 1;      % plot the surface of the flow
PlotTopo = 1;      % plot the topography
ContourValues = linspace(-1,1,21);
topoContourValues = 30;   % Contour levels for topo. 
                          % Set to either a scalar or vector

                        
geo_setzcolormaps;    % set up some useful default colormaps for land, water
flow_colormap = zTsunamiColors;
topo_colormap = zBathyTopo;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('topo_colormap')
 disp('*** You must define topo_colormap, e.g. to one of the maps')
 disp('    defined in geo_setzcolormaps.m')
end
topo2colors = mapzcolors(Grid.Z,topo_colormap);
 
cw=surf(Grid.X,Grid.Y,Grid.Z,topo2colors);
set(cw,'FaceColor','flat','EdgeColor','none');
 
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (PlotFlow)
if ~exist('flow_colormap')
  disp('*** You must define flow_colormap, e.g. to one of the maps')
  disp('     defined in geo_setzcolormaps.m')
end

eta2=zeros(size(Grid.Z));
eta2colors = mapzcolors(eta2,flow_colormap);

h=eta2-Grid.Z;

h(find(h<0))=0;
 
% set color to NaN in dry cells:
cutoff=.99e-3;
dry_ind=ones(size(Grid.X));
dry_ind(find((h./cutoff)<1))=NaN;
eta2colors(:,:,1) = eta2colors(:,:,1) .* dry_ind;

cw=surf(Grid.X,Grid.Y,eta2,eta2colors);
 %cw=surf(Grid.X,Grid.Y,eta2.*dry_ind,eta2colors);
set(cw,'FaceColor','flat','EdgeColor','none');

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%daspect([1 1 40000]) 
ylabel('Latitude','Fontsize',12)
xlabel('Longitude','Fontsize',12)
view([25 25])
hold off
 
return