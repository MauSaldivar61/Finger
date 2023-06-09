function [vol_handle,FV]=VoxelPlotter(VoxelMat,Vox_Size)
%detect the external voxels and faces
vol_handle=0;
if nargin==1
Vox_Size=1;
end
    FV=FindExternalVoxels(VoxelMat,Vox_Size);
%plot only external faces of external voxels
cla;
if size(FV.vertices,1)==0
    cla;
else
vol_handle=patch(FV,'FaceColor',[0 102 162]/255,'EdgeColor','k','LineWidth',0.5);
%use patchslim here for better results
end
end

