function fid=Grid2topofile(Grid,fnameout,topotypeout,nodata_valueout)

%function fid=Grid2topofile(Grid,fname,topotypeout)
%
%takes a Grid object and prints out a topo file, of topotype1,2 or 3.
%incoming Grid object should have Nans where no data is present. The output
%file has nodata_valueout for those values.
%
%function fid=Grid2topofile(Grid,fname,topotypeout,nodata_valueout)

if  nargin<4
    nodata_valueout=[];
end

nanind=find(Grid.Z==NaN);
nodataind=find(Grid.Z==Grid.nodata);

if (~isempty(nodata_valueout))
   Grid.Z(nanind)=nodata_valueout;
   Grid.Z(nodataind)=nodata_valueout;
else
   Grid.Z(nodataind)=Grid.nodata;
   Grid.Z(nanind)=Grid.nodata;
    
end

fid=fopen(fnameout,'wt');

if (topotypeout==1)
   
   Grid.Z=Grid.Z';
   Grid.X=Grid.X';
   Grid.Y=Grid.Y';

   Z=reshape(Grid.Z,1,Grid.mx*Grid.my);
   X=reshape(Grid.X,1,Grid.mx*Grid.my);
   Y=reshape(Grid.Y,1,Grid.mx*Grid.my);

   XYZ=[X;Y;Z];
  

   fprintf(fid,'%10.5f %10.5f %10.5f\n',XYZ);
   
elseif (topotypeout>1)
   
   fprintf(fid,'%-10g',Grid.mx);           fprintf(fid,'%s\n',[' ncols']);  
   fprintf(fid,'%-10g',Grid.my);           fprintf(fid,'%s\n',[' nrows']); 
   fprintf(fid,'%-10.5f',Grid.xlow);       fprintf(fid,'%s\n',[' xll']); 
   fprintf(fid,'%-10.5f',Grid.ylow);       fprintf(fid,'%s\n',[' yll']); 
   fprintf(fid,'%-10.5f',Grid.dx);         fprintf(fid,'%s\n',[' cellsize']);
   if nargin<4
      fprintf(fid,'%-10g',Grid.nodata); fprintf(fid,'%s\n',[' nodata_value']); 
   else
      fprintf(fid,'%-10g',nodata_valueout); fprintf(fid,'%s\n',[' nodata_value']);
   end
   
   if (topotypeout==2)
      fprintf(fid,'%g\n',Grid.Z');
   elseif (topotypeout==3)
      
       for j=1:Grid.my
       for i=1:Grid.mx
           fprintf(fid,'%10g',Grid.Z(j,i));
       end
           fprintf(fid,'\n');
       end
   end
end

fclose(fid);
return