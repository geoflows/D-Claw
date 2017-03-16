
% dig_plot_surface.m plots the surface of the flow  
%    eta = h + topo = q(:,:,1) + aux(:,:,1)
%    eta2 = augmented version of size (mx+1) by (my+1) for call to surf.



% mask points in eta2 covered by finer grids:
%erasegrid;
%eta2masked = eta2 .* covered_ind;

% mask points in eta2 where depth < tol (dry land)
%erasedry;
%eta2masked = eta2masked .* dry_ind;

%eta is colored by pressure


while ~exist('flow_colormap')
  disp('*** You must define flow_colormap, e.g. to one of the maps')
  disp('     defined in geo_setzcolormaps.m')
  break
  end


eta2colors = mapzcolors(eta2color,flow_colormap);

% set color to NaN in cells covered by finer grids:
set_covered_ind
eta2colors(:,:,1) = eta2colors(:,:,1) .* covered_ind;

% set color to NaN in dry cells:
geo_set_dry_ind
eta2colors(:,:,1) = eta2colors(:,:,1) .* dry_ind;

%cw=surf(X,Y,eta2,eta2colors);
cw=surf(X,Y,eta2.*dry_ind.*covered_ind,eta2colors); % Matlab has a bug regarding plotting edges--1/17/08 dlg

if (PlotGrid(level)==1)
    set(cw,'FaceColor','flat','EdgeColor',[0 0 0]);
else
    set(cw,'FaceColor','flat','EdgeColor','none');
end
%cw.AlphaData = covered_ind;

if (PlotGridEdges(level)==1)
    l1=line(xedge,0*xedge+yedge(1),eta2(:,1)+1000,'Color','k','LineWidth',1);
    l2=line(xedge,0*xedge+yedge(end),eta2(:,end)+1000,'Color','k','LineWidth',1);
    l3=line(0*yedge+xedge(1),yedge,eta2(1,:)+1000,'Color','k','LineWidth',1);
    l4=line(0*yedge+xedge(end),yedge,eta2(end,:)+1000,'Color','k','LineWidth',1);
    
end

%ylabel('Latitude','Fontsize',12)
%xlabel('Longitude','Fontsize',12)
%hold on
if (level>1&quiverplot)
    sq=10;
    xq = X(1:sq:end,1:sq:end);
    yq = Y(1:sq:end,1:sq:end);
    zq = eta2(1:sq:end,1:sq:end) + 10.;
    uq = hu2(1:sq:end,1:sq:end)./h2(1:sq:end,1:sq:end);
    vq = hv2(1:sq:end,1:sq:end)./h2(1:sq:end,1:sq:end);
    wq = 0.0*vq;
    quiver3(xq,yq,zq,uq,vq,wq,0,'r')
end