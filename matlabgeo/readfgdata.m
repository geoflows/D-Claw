function [fginfo,fgdata]= readfgdata(fgfilename) 

   fid=fopen(fgfilename);
     fgdata.t =fscanf(fid,'%g',1);      fscanf(fid,'%s',1);
     fginfo.mx=fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
     fginfo.my=fscanf(fid,'%d',1);      fscanf(fid,'%s',1);
     fginfo.xlow=fscanf(fid,'%g',1);    fscanf(fid,'%s',1);
     fginfo.ylow=fscanf(fid,'%g',1);    fscanf(fid,'%s',1);
     fginfo.xhi=fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
     fginfo.yhi=fscanf(fid,'%g',1);     fscanf(fid,'%s',1);
     fginfo.columns=fscanf(fid,'%d',1); fscanf(fid,'%s',1);
     
     size=[fginfo.columns,fginfo.mx*fginfo.my];
     data=fscanf(fid,'%g',size);
     data=data';
     
     fgdata.h   = reshape(data(:,1),fginfo.mx,fginfo.my);
     fgdata.hu  = reshape(data(:,2),fginfo.mx,fginfo.my);
     fgdata.hv  = reshape(data(:,3),fginfo.mx,fginfo.my);
     fgdata.b   = reshape(data(:,4),fginfo.mx,fginfo.my);
     fgdata.eta = reshape(data(:,5),fginfo.mx,fginfo.my);
     
     if (fginfo.columns>5)
        fgdata.etamin  = reshape(data(:,6),fginfo.mx,fginfo.my);
        fgdata.etamax  = reshape(data(:,7),fginfo.mx,fginfo.my);
     end
     
     
     if fginfo.mx>1
         fginfo.dx=(fginfo.xhi-fginfo.xlow)/(fginfo.mx-1);
     else
         fginfo.dx=0;
     end
     if fginfo.my>1
         fginfo.dy=(fginfo.yhi-fginfo.ylow)/(fginfo.my-1);
     else
         fginfo.dy=0;
     end 
   fclose(fid);
   
   
   
   