
% set some useful colors as RGB vectors and zcolormaps:
% to view, try e.g.  showcolor(DarkGreen)
%              or    zcolorbar(zTsunamiColors)
% see also the sample code at the end of this file.


Black         = [0.0  0.0  0.0];
White         = [1.0  1.0  1.0];
Red           = [1.0  0.0  0.0];
Green         = [0.0  1.0  0.0];
DarkGreen     = [0.1  0.4  0.0];
LightGreen    = [0.8  1.0  0.5];
Blue          = [0.0  0.0  1.0];
DarkBlue      = [0.2  0.2  0.7];
LightBlue     = [0.5  0.5  1.0];
BlueGreen     = [0.0  1.0  1.0];
Tan           = [0.9  0.8  0.2];
Tan           = [0.8  0.5  0.2];
Brown         = [0.9  0.8  0.2];
Gray8         = [0.8  0.8  0.8];
Purple        = [0.8  0.3  0.8];

%-------------------------------------------------
if(~exist('TsAmp'))
TsAmp = 0.6;   % default max amplitude for coloring
end

zTsunamiColors = [-TsAmp  Blue;
                       0  BlueGreen;
                   TsAmp  Red];

%-------------------------------------------------

zLandColors1 = [   0  DarkGreen;
                1000  Green;
                2000  LightGreen;
                3000  Tan;
                4000  White];

%-------------------------------------------------

%-------------------------------------------------

zLandColors2 = [   0  DarkGreen;
                 50  Green;
                 100  LightGreen;
                 200  Tan];

%-------------------------------------------------

zWaterLand = [-1000  DarkBlue;
               -500  Blue;
                  0  LightBlue;
                 .1  Tan;                  % beach region
                  5  Tan;
                  6  DarkGreen;
               1000  Green;
               2000  LightGreen;
               4000  Tan];

%-------------------------------------------------

zBathyTopo = [-1000  Brown;
                  0  Tan;
                 .1  DarkGreen;
               1000  Green;
               2000  LightGreen];
           
zBathyTopo = [-1000  Brown;
               -100  Tan;
                  0  DarkGreen;
                 .1  DarkGreen;
               1000  Green;
               2000  LightGreen];

%-------------------------------------------------

zAllRed = [0  Red;
           1  Red];

zAllBlue = [0  Blue;
            1  Blue];
        
        
zRedBlue = [0  Red;
            1  Blue];
                
zBlueRed = [0  Blue;
            1  Red];
        
zDarkRedBlue = [0  Red;
            1  DarkBlue];
                
zDarkBlueRed = [0  DarkBlue;
            1  Red];
        
zWhiteBlue = [0 White;
              1 Blue];
          

        
        

%-------------------------------------------------

  return

%-------------------------------------------------------------------------

% code to demo one of these maps:

zcolormap = zBathyTopo;

x = linspace(-1,1,101);
y = linspace(-1,1,101);
[X,Y] = meshgrid(x,y);
h = 1000*(X.^2 - Y.^2);

hcolors = mapzcolors(h,zcolormap);

figure(1)
surf(X,Y,h,hcolors,'EdgeColor','none');

zcolorbar(zcolormap)

