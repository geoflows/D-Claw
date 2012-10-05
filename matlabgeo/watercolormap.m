% watercolormap.m builds a colormap matrix 'water' that gives a scale from
% blue to red, without any green, so that land is easily distinguished from
% water.

water=zeros(64,3);
%shades of dark blue
water([1:31],1)=linspace(0,0,31)';
water([1:31],2)=linspace(0,1,31)';
water([1:31],3)=linspace(.5625,1,31)';

% light blue (cyan) zero value
water([32:33],1)=linspace(0,0,2)';
water([32:33],2)=linspace(1,1,2)';
water([32:33],3)=linspace(1,1,2)';

% cyan to yellow
water([34:45],1)=linspace(0,1,12)';
water([34:45],2)=linspace(1,1,12)';
water([34:45],3)=linspace(1,0,12)';
%shades of orange
water([46:57],1)=linspace(1,1,12)';
water([46:57],2)=linspace(1,0,12)';
water([46:57],3)=linspace(0,0,12)';

%shades of red
water([58:64],1)=linspace(1,.5,7)';
water([58:64],2)=linspace(0,0,7)';
water([58:64],3)=linspace(0,0,7)';






