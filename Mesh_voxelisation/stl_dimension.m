function [minx, maxx, miny, maxy, minz, maxz] = stl_dimension(filename)
[~,gridCOx,~,~] = VOXELISE(1024*5,1,1,filename,'xyz');
[~,~,gridCOy,~] = VOXELISE(1,1024*5,1,filename,'xyz');
[~,~,~,gridCOz] = VOXELISE(1,1,1024*5,filename,'xyz');
% [~,gridCOx,gridCOy,gridCOz] = VOXELISE(512,512,256,filename,'xyz');
minx = min(gridCOx(:));  maxx = max(gridCOx(:));
miny = min(gridCOy(:));  maxy = max(gridCOy(:));
minz = min(gridCOz(:));  maxz = max(gridCOz(:));
disp([maxx - minx, maxy - miny, maxz - minz])