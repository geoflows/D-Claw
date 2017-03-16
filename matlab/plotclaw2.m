%
% PLOTCLAW2 is the main driver routine for plotting 2d graphs for Clawpack.
%
%    PLOTCLAW2 is the main routine that the user calls to step through a
%    series of frames (i.e. fort.tXXXX and fort.qXXXX files).
%
%    Various parameters can be set in SETPLOT2.  The default version in
%    claw/matlab/setplot2.m can be copied to your directory and modifed to
%    set things up differently, or type 'k' at the prompt to get keyboard
%    control and change a value.
%
%    See also SETPLOT, PLOTFRAME2.


%
% plotclaw2.m
%
% generic plotting routine for clawpack and amrclaw output in matlab
% R. J. LeVeque, 1999
%
% Various parameters are set in setplot2.m
% The default version in claw/matlab/setplot2.m can be copied to your
% directory and modified to set things up differently, or type 'k'
% at the prompt to get keyboard control and change a value.
%
%---------------------------------------------------------------------

clawdim = 2;

disp(' ')
disp('plotclaw2  plots 2d results from clawpack or amrclaw')

% set plotting parameters:
whichfile = which('setplot2');
if strcmp(whichfile,'')
  disp('*** No setplot2 file found')
else
  inp = input(['Execute setplot2 (default = no)? '],'s');
  inpd = findstr('y',lower(inp));
  if (inpd == 1)
    setplot2;
    disp(['Executing m-script ' whichfile])
    disp(' ')
  end
end
disp(' ')

% the file setprob.m can be used to set up any necessary physical parameters
% or desired values of plotting parameters for this particular problem.

whichfile = which('setprob');
if strcmp(whichfile,'')
  %disp('*** No setprob file found')
else
  disp(['Executing m-script ' whichfile])
  disp(' ')
  setprob
end

%=============================================
% MAIN LOOP ON FRAMES:
%=============================================

Frame = -1;  % initialize frame counter

if ~exist('MaxFrames')
   disp('MaxFrames parameter not set... you may need to execute setplot2')
   return
   end

set_value('outputdir','OutputDir','./');
set_value('outputflag','OutputFlag','ascii');

amrdata = [];
while Frame <= MaxFrames

  % pause for input from user to determine if we go to next frame,
  % look at data, or skip around.  This may reset Frame counter.

  Frame_old = Frame;
  %queryframe;  % this sets Frame
  %  QUERYFRAME used by PLOTCLAW1, PLOTCLAW2 and PLOTCLAW3 to loop through data.
%
%      QUERYFRAME is normally called directly from PLOWCLAW<N>, and so is
%      not typically called by the user.  However, the user can supress the
%      query option by setting the 'NoQuery' parameter to 1.  This is used
%      for example, in making movies.
%
%      See also PLOTCLAW1, PLOTCLAW2, PLOTCLAW3.

if exist('NoQuery')
  if NoQuery == 1
    % set NoQuery=1 if you want the plots to be produced with no
    % query in between.  Particularly useful if you want lots of frames to
    % be printed out for an animation (put a command like makeframegif
    % in afterframe.m and set NoQuery=1)
    pause(1)
    Frame = Frame + 1;
    if Frame > MaxFrames
      return;   % break out of plotclawN after last frame
    end
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
    amrdata = [];
  elseif strcmp(inp,'j')
    Frame = input('Frame to jump to? ');
  elseif strcmp(inp,'i')
    if clawdim == 1
      infoplot1
      disp(' ')
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      infoframe1
    end
    if clawdim == 2
      infoplot2
      disp(' ')
      disp('hit <return> for information about this frame')
      pause
      infoframe2
    end
    if clawdim == 3
      infoplot3
      disp(' ');
      disp('hit <return> for information about this frame');
      pause
      infoframe3
    end;
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

  if (Frame ~= Frame_old | isempty(amrdata))
    [amrdata,t] = readamrdata(clawdim,Frame,outputdir,outputflag);
  end;

  plotframe2;

end % main loop on frames
