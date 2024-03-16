% SAVEFRAMESVN produces a .svg file of the current frame.
%
%
%
%     See also PRINTGIF, PRINTJPG, PRINTPNG
%

framest = num2str(Frame);
while(size(framest,2))<5,
  framest = ['0',framest];
  end;
fname = ['frame' framest];
saveas(gcf,fname,'svg');