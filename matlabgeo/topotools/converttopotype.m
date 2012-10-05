function [fidout]= converttopotype(fnamein,fnameout,topotypein,topotypeout,nodata_valueout,nodata_valuein)

%function [fidout]= converttopotype(fnamein,fnameout,topotypein,topotypeout)
%
% converts topography files of one topotype to files of another topotype.
% nodata_valueout is optional. If not set it is the same value as
% nodata_value of the input file if topotype is 2 or 3. 
% If input file is of topotype 1,
% the default nodata_valueout is NaN, which is not recommended.
%
%function [fidout]=
%converttopotype(fnamein,fnameout,topotypein,topotypeout,nodata_valueout,notdata_valuein)

if nargin<6
  nodata_valuein=[];
end
if nargin<5
    nodata_valueout=[];
end


if isempty(nodata_valuein)
  Grid=topofile2Grid(fnamein,topotypein);
else
  Grid=topofile2Grid(fnamein,topotypein,nodata_valuein)
end

if isempty(nodata_valueout)
    nodata_valueout= Grid.nodata;
end
     
if isempty(nodata_valueout)
    fid=Grid2topofile(Grid,fnameout,topotypeout);
else
    fid=Grid2topofile(Grid,fnameout,topotypeout,nodata_valueout);
end
return