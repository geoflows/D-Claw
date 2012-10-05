
% geo_plot_surface.m plots the surface of the flow  
%    eta = h + topo = q(:,:,1) + aux(:,:,1)
%    eta2 = augmented version of size (mx+1) by (my+1) for call to surf.



% mask points in eta2 covered by finer grids:
%erasegrid;
%eta2masked = eta2 .* covered_ind;

% mask points in eta2 where depth < tol (dry land)
%erasedry;
%eta2masked = eta2masked .* dry_ind;


if ~exist('flow_colormap')
  disp('*** You must define flow_colormap, e.g. to one of the maps')
  disp('     defined in geo_setzcolormaps.m')
  break
  end


eta2colors = mapzcolors(eta2,flow_colormap);


% set color to NaN in dry cells:
geo_set_dry_ind
eta2colors(:,:,1) = eta2colors(:,:,1) .* dry_ind;

%cw=surf(X,Y,eta2,eta2colors);
cw=surf(X,Y,eta2.*dry_ind,eta2colors); % Matlab has a bug regarding plotting edges--1/17/08 dlg

if (PlotGrid(1)==1)
    set(cw,'FaceColor','flat','EdgeColor',[0 0 0]);
else
    set(cw,'FaceColor','flat','EdgeColor','none');
end


ylabel('Latitude','Fontsize',12)
xlabel('Longitude','Fontsize',12)
