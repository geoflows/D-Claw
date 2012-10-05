
% Called from plotframe2 if 
%   PlotType == 11 (colored surface plot of topography and/or fluid)
%   PlotType == 12 (contour plot of topography and/or fluid)
% Set
%   PlotTopo == 1 to plot topography, 0 to suppress
%   PlotFlow == 1 to plot flow, 0 to suppress

if ~exist('PlotTopo')
  PlotTopo = 1;
  disp('*** setting PlotTopo = 1.  Set to 0 to suppress topography')
  end
if ~exist('PlotFlow')
  PlotFlow = 1;
  disp('*** setting PlotFlow = 1.  Set to 0 to suppress plotting flow')
  end

h = fgdata(ng).h;                % depth
hu = fgdata(ng).hu;               % momentum
hv = fgdata(ng).hv;
eta = fgdata(ng).eta;              % surface
topo = fgdata(ng).b;   % topography

my=fginfo(ng).my;
mx=fginfo(ng).my;
xlow=fginfo(ng).xlow;
xhi = fginfo(ng).xhi;
ylow=fginfo(ng).ylow;
yhi = fginfo(ng).yhi;
dx= fginfo(ng).dx;
dy= fginfo(ng).dy;

xedge=linspace(xlow,xhi,mx);
yedge=linspace(ylow,yhi,my);
xcenter=linspace(xlow+0.5*dx , xhi-0.5*dx, mx-1);
ycenter=linspace(ylow+0.5*dy , yhi-0.5*dy, my-1);

if PlotType == 11

    [X,Y]=ndgrid(xedge,yedge);  

    h2 = h; %  h2(:,my+1) = h2(:,my);   h2(mx+1,:) = h2(mx,:);
    hu2 = hu;%   hu2(:,my+1) = hu2(:,my);   hu2(mx+1,:) = hu2(mx,:);
    hv2 = hv; %  hv2(:,my+1) = hv2(:,my);   hv2(mx+1,:) = hv2(mx,:);
    eta2 = eta;%   eta2(:,my+1) = eta2(:,my);   eta2(mx+1,:) = eta2(mx,:);
    topo2 = topo;%   topo2(:,my+1) = topo2(:,my);   topo2(mx+1,:) = topo2(mx,:);
        
    if PlotTopo
      % plot topography:
      geo_plot_topofg   
      end
    
    if PlotFlow
      % plot surface of flow:
      geo_plot_surfacefg
      end
          
    end

%----------------------------------

if PlotType == 12

    % contour plots

    if ~exist('topoContourValues')
      topoContourValues = 30;
      disp(sprintf('*** Using topoContourValues = %g',topocontourValues))
      end
  
    if PlotTopo
      % plot topography:
      contour(xcenter,ycenter,topo,topoContourValues,'g')
      end
    
    if PlotFlow
      % plot surface of flow:
      contour(xcenter,ycenter,eta,ContourValues,'r')
      end

    end
