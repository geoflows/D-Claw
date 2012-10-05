fid = fopen('setgauges.data');
mgauges=fscanf(fid,'%g',[1]);
gaugelocs = fscanf(fid,'%g',[3,inf])';
status = fclose(fid);

hold on
for i=1:mgauges
  plot3(gaugelocs(i,2),gaugelocs(i,3),100,'o','MarkerSize',30)
  plot3(gaugelocs(i,2),gaugelocs(i,3),100,'+','MarkerSize',25)
  text(gaugelocs(i,2),gaugelocs(i,3),100,num2str(gaugelocs(i,1)),'FontSize',20)
  end
hold off
