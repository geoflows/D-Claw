function zcolors = mapzcolors(z,zcolormap)
%
% maps a 2d array z of z values to RGB triples that can be used in 
%   surf(x,y,z,zcolors)
% The triples are determined by the array zcolormap, with lines of the form
%   zval Rval Gval Bval
% indicating that value z=zval should be mapped to this color.
% For z values in between values in this array, linear interpolation is used.
% z values outside the range [min(zval) max(zval)] are mapped to the
% corresponding endpoint.
%
% If z is m by n then zcolors is m by n by 3.
%
% This is more useful than the standard matlab colormap function since it
% allows different color maps to be used for different data on the same plot.
%
% To plot a color bar for zcolormap:
%   zcolorbar(zcolormap,ifig)
% will do figure(ifig) and plot a colorbar.


zmin = min(zcolormap(:,1));
zmax = max(zcolormap(:,1));
Z = z;

% restrict to range of zcolormap, leaving NaN's unchanged
ij = find(~isnan(Z));
Z(ij) = max(Z(ij),zmin);
Z(ij) = min(Z(ij),zmax);

% interpolate each value of Z array into zcolormap:

if (exist('discrete_colormap')&(discrete_colormap==1))
	ppR = griddedInterpolant(zcolormap(:,1),zcolormap(:,2),'next');
	Rvalues = ppR(Z);
	ppG = griddedInterpolant(zcolormap(:,1),zcolormap(:,3),'next');
	Gvalues = ppG(Z);
	ppB = griddedInterpolant(zcolormap(:,1),zcolormap(:,4),'next');
	Bvalues = ppB(Z);
else
    ppR = interp1(zcolormap(:,1),zcolormap(:,2),'linear','pp');
    Rvalues = ppval(ppR, Z);
    ppG = interp1(zcolormap(:,1),zcolormap(:,3),'linear','pp');
    Gvalues = ppval(ppG, Z);
    ppB = interp1(zcolormap(:,1),zcolormap(:,4),'linear','pp');
    Bvalues = ppval(ppB, Z);
end

zcolors = nan([size(z) 3]);
zcolors(:,:,1) = Rvalues;
zcolors(:,:,2) = Gvalues;
zcolors(:,:,3) = Bvalues;

%keyboard
