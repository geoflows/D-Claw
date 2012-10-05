

if exist('NoQuery')
  if NoQuery == 1
    % set NoQuery=1 if you want the plots to be produced with no
    % query in between.  Particularly useful if you want lots of frames to
    % be printed out for an animation (put a command like makeframegif
    % in afterframe.m and set NoQuery=1)
    pause(1)
    Frame = Frame + 1;
    return
  end
end

inp = 'k';
while strcmp(inp,'k')

  inp = input(...
      ['Hit <return> for next plot, or type k, r, rr, j, i, q, or ?  '],'s');

  if strcmp(inp,'?')
    disp('  k  -- keyboard input.  Type any commands and then "return"')
    disp('  r  -- redraw current frame, without re-reading data')
    disp('  rr -- re-read current file,and redraw frame');
    disp('  j  -- jump to a particular frame')
    disp('  i  -- info about parameters and solution')
    disp('  q  -- quit')
  elseif strcmp(inp,'k')
    keyboard
  elseif strcmp(inp,'r')
    % redraw:  leave Frame counter alone
    if Frame==-1
      disp('Cannot redraw yet')
      inp = 'k';
    end
  elseif strcmp(inp,'rr')
    % redraw frame AND re-read data
    fgdata(ng) = [];
  elseif strcmp(inp,'j')
    Frame = input('Frame to jump to? ');
  elseif strcmp(inp,'i')
      disp(' ')
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      disp(['Grid Number = ',num2str(ng)]);
      disp(['Time = ', num2str(fgdata.t)]);
      disp(['fginfo(',num2str(ng),'):']); fginfo(ng)
      
      inp = 'k';
  elseif isempty(inp)
    % go to next frame
    Frame = Frame + 1;
  elseif (~strcmp(inp,'q'))
    % quit handled separately below.
    % Otherwise unrecognized input, go back and try again
    inp = 'k';
  end % if strcmp
end % while strcmp

if strcmp(inp,'q')
  % quit now
  break
end