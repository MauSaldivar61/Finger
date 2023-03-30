function [G_S, G_H, G_F, G_F1] = import_SoftHard_design_v2(filename_S, filename_H, filename_F, filename_F1, RVE_xy, RVE_z)

% We first need to obtain the dimensions of the devices:
% The hard part:
[minx_H,maxx_H,miny_H,maxy_H,minz_H,maxz_H] = stl_dimension(filename_H);

% The soft part:
[minx_S,maxx_S,miny_S,maxy_S,minz_S,maxz_S] = stl_dimension(filename_S);

% The Fluid TPMS part:
[minx_T,maxx_T,miny_T,maxy_T,minz_T,maxz_T] = stl_dimension(filename_F);

% The Fluid part:
[minx_F,maxx_F,miny_F,maxy_F,minz_F,maxz_F] = stl_dimension(filename_F1);

% Obtain the true maximum
minx = min([minx_H, minx_S, minx_T, minx_F]);
miny = min([miny_H, miny_S, miny_T, miny_F]);
minz = min([minz_H, minz_S, minz_T, minz_F]);
maxx = max([maxx_H, maxx_S, maxx_T, maxx_F]);
maxy = max([maxy_H, maxy_S, maxy_T, maxy_F]);
maxz = max([maxz_H, maxz_S, maxz_T, maxz_F]);
% [minz maxz]
% Obtain the array of RVEs required
RVE_xg = linspace( minx + RVE_xy/2, maxx, round((maxx-minx)./(RVE_xy)));
RVE_yg = linspace( miny + RVE_xy/2, maxy, round((maxy-miny)./(RVE_xy)));
RVE_zg = linspace( minz + RVE_xy/2, maxz, abs(round((maxz-minz)./(RVE_z))));

% And now you can voxelize the designs
[G_H] = double(VOXELISE(RVE_xg,RVE_yg,RVE_zg,filename_H,'xyz')); %Voxelize HARD design to given resolution
[G_S] = double(VOXELISE(RVE_xg,RVE_yg,RVE_zg,filename_S,'xyz')); %Voxelize SOFT design to given resolution
[G_F] = double(VOXELISE(RVE_xg,RVE_yg,RVE_zg,filename_F,'xyz')); %Voxelize SOFT design to given resolution
[G_F1] = double(VOXELISE(RVE_xg,RVE_yg,RVE_zg,filename_F1,'xyz')); %Voxelize SOFT design to given resolution

% We rearrange the designs so the long side is in the column direction:
[~,iM] = max(size(G_H));
if iM == 1
    G_H = permute(G_H,[ 2, 1, 3]);
    G_S = permute(G_S,[ 2, 1, 3]);
    G_F = permute(G_F,[ 2, 1, 3]);
    G_F1 = permute(G_F1,[ 2, 1, 3]);
end
% The rows always are inverted, so we invert the design across the rows:
G_H = G_H(end:-1:1,:,:);
G_S = G_S(end:-1:1,:,:);
G_F = G_F(end:-1:1,:,:);
G_F1 = G_F1(end:-1:1,:,:);
