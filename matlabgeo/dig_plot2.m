
% Called from plotframe2 if
%   PlotType == 13 (colored surface plot of topography and/or fluid)
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

meqnS = size(data);
meqn = meqnS(2);
neta = meqn;
h = reshape(data(:,1),mx,my);                % depth
hu = reshape(data(:,2),mx,my);               % momentum
hv = reshape(data(:,3),mx,my);
hm = reshape(data(:,4),mx,my);
p = reshape(data(:,5),mx,my);
trac = reshape(data(:,6),mx,my);
erode = reshape(data(:,7),mx,my);
eta = reshape(data(:,neta),mx,my);              % surface
topo = reshape(data(:,neta)-data(:,1),mx,my);   % topography

if PlotType == 13

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
    hm2 = hm;   hm2(:,my+1) = hm2(:,my);   hm2(mx+1,:) = hm2(mx,:);
    p2 = p;   p2(:,my+1) = p2(:,my);   p2(mx+1,:) = p2(mx,:);
    trac2 = trac;   trac2(:,my+1) = trac2(:,my);   trac2(mx+1,:) = trac2(mx,:);
    erode2 = erode;   erode2(:,my+1) = erode2(:,my);   erode2(mx+1,:) = erode2(mx,:);
    eta2 = eta;   eta2(:,my+1) = eta2(:,my);   eta2(mx+1,:) = eta2(mx,:);
    topo2 = topo;   topo2(:,my+1) = topo2(:,my);   topo2(mx+1,:) = topo2(mx,:);

    sv = hm2./h2;
    rho = 2700.0*sv + 1000.0*(1.-sv);
    theta = 9.81*ones(size(rho));
    deg2rad = 3.14159/180.0;
    flumelen = 78.0;
    flumerad = 10.0;
    theta1 = 31.0;
    theta2 = 3.0;
    D2 = flumelen + flumerad*(theta1 - theta2)*deg2rad;
    theta(X(:,1)<=flumelen,:) = theta1;
    theta(X(:,1)>=D2,:) = theta2;
    theta(X(:,1)>flumelen&X(:,1)<D2,:) = theta1 - (X(X(:,1)>flumelen&X(:,1)<D2,:) - flumelen)/(deg2rad*flumerad);
    theta = deg2rad*theta;
    if mq == 5
        eta2color = p2./(cos(theta).*9.81.*rho.*h2);
        %eta2color = (p2 - 9.81*1000.0*h2)./(9.81*rho.*h2 - 9.81*1000.0.*h2);
        %eta2color = (9.81*rho.*h2 - p2)./(9.81*rho.*h2 - 9.81*1000.0.*h2);
    elseif mq==1
        eta2color = h2;
    elseif mq==2
        eta2color = hu2./h2;
    elseif mq==3
        eta2color = hv2./h2;
    elseif mq==4
        eta2color = sv;
    elseif mq==6
        eta2color = trac2./h2;
    elseif mq==7
        eta2color = erode2;
    elseif mq==neta
        eta2color=eta2;
    elseif mq==neta+1
        eta2color = sqrt((hu2./h2).^2 + (hv2./h2).^2);
    end


    if PlotTopo
      % plot topography:
      geo_plot_topo
      end

    if PlotFlow
      % plot surface of flow:
      dig_plot_surface
      end

    view(2)  % top view
    axis tight

end



%----------------------------------


