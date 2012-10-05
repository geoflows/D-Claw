function Grid=topofile2Grid(fname,topotype,nodata_value)

% Grid=topofile2Grid(fname,topotype)
% topofile2Grid takes a file of topotype1,2 or 3, and returns the object Grid, 
% containing matrices Grid.Z, Grid.X Grid.Y, 
% and parameters, Grid.mx, Grid.my, Grid.xlow, Grid.ylow,
% Grid.dx, Grid.dy, Grid.xend, Grid.yend and Grid.nodata
%
% Grid=topofile2Grid(fname,topotype,nodata_value)
% use nodata_value to specify a nodata_value for the Grid object.
    
if nargin<3
    nodata_value=[];
end

if (topotype~=1&topotype~=2&topotype~=3)
    display(['topofile2Grid ERROR: topotype must be 1,2 or 3'])
    return
end

% read data 
fid=fopen(fname);

if (topotype==1)
    if nargin<3 | isempty(nodata_value)
         nodata_value= NaN;
    end
    Grid.nodata=nodata_value;
    data=fscanf(fid,'%g',[3,inf]);

    data=data';
    x=data(:,1);
    y=data(:,2);
    z=data(:,3);
    clear data;

    Grid.xlow=x(1);
    Grid.ylow=y(end);

    Grid.dx=max(diff(x));
    Grid.dy=-min(diff(y));

%    inddiffy=find(diff(y)<0);
    inddiffx=find(diff(x)<0);

    Grid.mx=inddiffx(1);
    Grid.my= length(x)/Grid.mx;

    if(min(Grid.dx)<0)
       display(['xyzfile2Grid ERROR: dx<0: data does not advance as expected.'])
       return
    end
    Grid.Z=reshape(z,Grid.mx,Grid.my)';
    Grid.X=reshape(x,Grid.mx,Grid.my)';
    Grid.Y=reshape(y,Grid.mx,Grid.my)';
    Grid.xend= Grid.xlow + (Grid.mx-1)*Grid.dx;
    Grid.yend= Grid.ylow + (Grid.my-1)*Grid.dy;
    clear x y z inddiffx
    
else
    Grid.mx=fscanf(fid,'%g',1);      fscanf(fid,'%s',1); 
    Grid.my=fscanf(fid,'%g',1);      fscanf(fid,'%s',1); 
    Grid.xlow=fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
    Grid.ylow=fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
    Grid.dx=fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
    Grid.nodata=fscanf(fid,'%g',1);  fscanf(fid,'%s',1);
    if nargin==3 & ~isempty(nodata_value)
         Grid.nodata=nodata_value;
    end
    
    Grid.dy=Grid.dx;

    Grid.Z=fscanf(fid,'%g',[Grid.mx,Grid.my]);
    
    if topotype==2
      Grid.Z=reshape(Grid.Z,Grid.mx,Grid.my);
    end
    Grid.Z=Grid.Z';
    Grid.xend= Grid.xlow + (Grid.mx-1)*Grid.dx;
    Grid.yend= Grid.ylow + (Grid.my-1)*Grid.dy;
    [Grid.X,Grid.Y]= meshgrid(linspace(Grid.xlow,Grid.xend,Grid.mx),...%
                           linspace(Grid.ylow,Grid.yend,Grid.my));
    Grid.Y = flipud(Grid.Y);
end

fclose(fid);

Grid;

return