
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

h = reshape(data(:,1),mx,my);                % depth
hu = reshape(data(:,2),mx,my);               % momentum
hv = reshape(data(:,3),mx,my);
eta = reshape(data(:,4),mx,my);              % surface
topo = reshape(data(:,4)-data(:,1),mx,my);   % topography

if PlotType == 11

    [X,Y]=ndgrid(xedge,yedge);   % (mx+1) by (my+1) arrays

    % augmented versions for surface commands, which ignore last row and col:
    %h2 = nan(mx+1,my+1);     h2(1:mx,1:my) = h;
    %hu2 = nan(mx+1,my+1);    hu2(1:mx,1:my) = hu;
    %hv2 = nan(mx+1,my+1);    hv2(1:mx,1:my) = hv;
    %eta2 = nan(mx+1,my+1);   eta2(1:mx,1:my) = eta;
    %topo2 = nan(mx+1,my+1);  topo2(1:mx,1:my) = topo;

    h2 = h;   h2(:,my+1) = h2(:,my);   h2(mx+1,:) = h2(mx,:);
    hu2 = hu;   hu2(:,my+1) = hu2(:,my);   hu2(mx+1,:) = hu2(mx,:);
    hv2 = hv;   hv2(:,my+1) = hv2(:,my);   hv2(mx+1,:) = hv2(mx,:);
    eta2 = eta;   eta2(:,my+1) = eta2(:,my);   eta2(mx+1,:) = eta2(mx,:);
    topo2 = topo;   topo2(:,my+1) = topo2(:,my);   topo2(mx+1,:) = topo2(mx,:);
        
    if PlotTopo
      % plot topography:
      geo_plot_topo   
      end
    
    if PlotFlow
      % plot surface of flow:
      geo_plot_surface
      end
          
    view(2)  % top view
    axis tight
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
