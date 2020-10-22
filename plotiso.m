function p = plotiso(x, y, z, v, thresh, color);
% function p = plotiso(x, y, z, v, thresh, color, transparency);
% Plot an isosurface of the data v, specified on the 3-d arrays, x, y, z;
% with thresh as the threshold.  Makes a nice looking plot
% Adam Cohen, 5 Oct. 2011

pdata = isosurface(x, y, z, v, thresh);
p = patch(pdata);
set(p, 'FaceColor',color,'EdgeColor','none');
alpha(.5);
daspect([1 1 1])
axis off
lightangle(45,30); lighting phong; 
% set(p,'AmbientStrength',.5)
% set(p,'SpecularColorReflectance',.8,'SpecularExponent',50);