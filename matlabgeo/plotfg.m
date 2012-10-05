

setprob;
setplot2;

ng = input('Fixed grid number to plot? ');
while (isempty(ng))
    ng = input('Enter a Grid Number? ');
end
ngs='00';
ngs(1)=num2str((ng-mod(ng,10))/10);
ngs(2)=num2str(mod(ng,10));
fgfilename1 = ['fort.fg',ngs,'_'];

if (exist([fgfilename1,'0001'],'file'));
    cc=1;
else
   disp([' ']);
   disp([fgfilename1,'0001', ' does not exist.']);
   disp([' ']);
   cc=0;
 end

Frame=0;

while (cc)
 queryfg;
 framename='xxxx';
 framenumber=Frame;
 for ind=4:-1:1 % fort.fgYY_XXXX
   framename(ind) = num2str(mod(framenumber,10));
   framenumber=floor(framenumber/10);
 end

 fgfilename=[fgfilename1,framename];

 if (exist(fgfilename,'file'))
     [fginfo(ng),fgdata(ng)]= readfgdata(fgfilename);
     %fgdata(ng).eta= fgdata(ng).h + fgdata(ng).b;
     disp(['Reading Data from ./',fgfilename]);
     disp([' ']);
     disp(['Frame ', num2str(Frame), ' at time t = ', num2str(fgdata(ng).t)]);
     disp([' ']);
     clf
     hold on
     geo_plot2fg;
     %afterframe;
     axis([fginfo(ng).xlow,fginfo(ng).xhi,fginfo(ng).ylow,fginfo(ng).yhi]);
     title(['t = ',num2str(fgdata(ng).t)],'FontSize', 12)
     hold off
     
 else
   disp([' ']);
   disp([fgfilename, ' does not exist.']);
   disp([' ']);
   cc=0;
 end
 hold off
end
hold off