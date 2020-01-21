% MAKEFRAMEPNG produces a .png file of the current frame.
%
%     Produce a png file of this frame with numbering so that gifmerge
%     can be used to make an animation:  frame00000.png, frame00001.png, etc.
%
%     By putting the command makeframepng in the afterframe.m file,
%
%     Then catenate all these files together into an animation.
%     There are various ways to do this in unix or linux depending on what
%     software is on your system, e.g.
%         convert -delay 20 frame*.png movie.gif
%
%
%     See also PRINTGIF, PRINTJPG, PRINTPNG
%

framest = num2str(Frame);
while(size(framest,2))<5,
  framest = ['0',framest];
  end;
fname = ['frame' framest];
printpng(fname);
