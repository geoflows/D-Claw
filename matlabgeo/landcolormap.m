% landcolormap.m builds a colormap matrix 'land' that gives a land topographical look

land=zeros(64,3);
%shades of green
land([1:16],1)=linspace(0,0,16)';
land([1:16],2)=linspace(.5,1,16)';
land([1:16],3)=linspace(0,0,16)';

land([17:44],1)=linspace(0,.8,28)';
land([17:44],2)=linspace(1,1,28)';
land([17:44],3)=linspace(0,.45,28)';
%shades of brown
land([45:60],1)=linspace(.8,.8,16)';
land([45:60],2)=linspace(1,.25,16)';
land([45:60],3)=linspace(.45,.15,16)';
%snow
land([61:64],1)=linspace(.9,1,4)';
land([61:64],2)=linspace(.8,1,4)';
land([61:64],3)=linspace(.5,1,4)';

