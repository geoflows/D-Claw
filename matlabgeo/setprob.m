
% This routine is called after setplot2.m by plotclaw2.
%
% Set some additional things for plotting GeoClaw output.
%

PlotType = 11;     % = 11 for colored surface plot
                   % = 12 for contour plot

PlotFlow = 1;      % plot the surface of the flow
PlotTopo = 1;      % plot the topography
ContourValues = linspace(-1,1,21);
topoContourValues = 30;   % Contour levels for topo. 
                          % Set to either a scalar or vector

TsAmp = .6;   %max amplitude for coloring
                        


geo_setzcolormaps;    % set up some useful default colormaps for land, water
flow_colormap = zTsunamiColors;
topo_colormap = zLandColors2;

% or for non default colormaps for the water 
% set flowcolormatrix to any colormap desired (ie any m by 3 rgb matrix) 
%[flowcolormatrix,ncolors]=deacolor; 
%[ncolors,n]=size(flowcolormatrix);
%flow_colormap=[linspace(-TsAmp,TsAmp,ncolors)',flowcolormatrix];

