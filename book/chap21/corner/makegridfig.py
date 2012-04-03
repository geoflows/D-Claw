
# make the figure for the book showing the interface.

from pylab import *

figure(5)
clf()
plot([0,0,1],[-1,0,0.55],'g',linewidth=2)
text(.4,-.6,'right',fontsize=15)
text(-.6,0,'left',fontsize=15)
axis('scaled')
xlim([-1,1])
ylim([-1,1])
savefig('interface.png')
print "Created interface.png"

