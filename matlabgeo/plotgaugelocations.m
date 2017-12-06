fid = fopen('setgauges.data');
for i=1:6
    fgetl(fid);
end
[mgauges,count]=fscanf(fid,'%g',[1]);
junk = fscanf(fid,'%s',2);

[gaugelocs,count2] = fscanf(fid,'%g',[5,inf]);
gaugelocs=gaugelocs';
status = fclose(fid);

hold on
for i=1:mgauges
  plot3(gaugelocs(i,2),gaugelocs(i,3),9000,'o','MarkerSize',30)
  plot3(gaugelocs(i,2),gaugelocs(i,3),9000,'+','MarkerSize',25)
  text(gaugelocs(i,2),gaugelocs(i,3),9000,num2str(gaugelocs(i,1)),'FontSize',20)
  end
hold off
