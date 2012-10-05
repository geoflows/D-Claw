function plotgauges(gaugenumber)



fid = fopen('fort.gauge');
gaugedata = fscanf(fid,'%g',[6,inf])';
fclose(fid);



ngauges = max(gaugedata(:,1));
tmax = max(gaugedata(:,2));
hmax = max(gaugedata(:,3));
etamax = max(gaugedata(:,6));

for n=1:ngauges
   in = find(gaugedata(:,1) == n);
   t  = gaugedata(in,2);
   h  = gaugedata(in,3);
   hu = gaugedata(in,4);
   hv = gaugedata(in,5);
   eta  = gaugedata(in,6);

   plot(t,eta,'r',t,0*t,'k','LineWidth',2)
   %axis([0 tmax -1.1*etamax 1.1*etamax])
   title(['Surface elevation at Gauge ' num2str(n)],'FontSize',20)
   fname = ['gauge' num2str(n)];
   
   an=input('Do you want to print the figure as a .jpg? (y/n)\n','s');
   if (strcmp(an,'y'))
     eval('printjpg(fname)');
   end
   input('Hit return for next gauge...');
   end
