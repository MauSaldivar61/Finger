function [G_S, G_H] = import_SoftHard_design(filename_S, filename_H, RVE_xy, RVE_z)

% We first need to obtain the dimensions of the devices:
% The hard part:
[minx_H,maxx_H,miny_H,maxy_H,minz_H,maxz_H] = stl_dimension(filename_H);

% The soft part:
[minx_S,maxx_S,miny_S,maxy_S,minz_S,maxz_S] = stl_dimension(filename_S);

% Obtain the true maximum
minx = min([minx_H, minx_S]);
miny = min([miny_H, miny_S]);
minz = min([minz_H, minz_S]);
maxx = max([maxx_H, maxx_S]);
maxy = max([maxy_H, maxy_S]);
maxz = max([maxz_H, maxz_S]);
% [minz maxz]
% Obtain the array of RVEs required
RVE_xg = linspace( minx, maxx, round((maxx-minx)./(RVE_xy)));
RVE_yg = linspace( miny, maxy, round((maxy-miny)./(RVE_xy)));
RVE_zg = linspace( minz, maxz, abs(round((maxz-minz)./(RVE_z))));

% And now you can voxelize the designs
[G_H] = double(VOXELISE(RVE_xg,RVE_yg,RVE_zg,filename_H,'xyz')); %Voxelize HARD design to given resolution
[G_S] = double(VOXELISE(RVE_xg,RVE_yg,RVE_zg,filename_S,'xyz')); %Voxelize SOFT design to given resolution

% We rearrange the designs so the long side is in the column direction:
[~,iM] = max(size(G_H));
if iM == 1
    G_H = permute(G_H,[ 2, 1, 3]);
    G_S = permute(G_S,[ 2, 1, 3]);
end
% The rows always are inverted, so we invert the design across the rows:
G_H = G_H(end:-1:1,:,:);
G_S = G_S(end:-1:1,:,:);