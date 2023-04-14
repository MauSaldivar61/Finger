% LOADING DATA
arr = load('example.mat');
Data = arr.Data(:,:,1:end-1); 

vox_xyz = 1;
gridX = linspace(vox_xyz/2,size(Data,2)*vox_xyz,size(Data,2)*vox_xyz);
gridY = linspace(vox_xyz/2,size(Data,1)*vox_xyz,size(Data,1)*vox_xyz);
gridZ = linspace(vox_xyz/2,size(Data,3)*vox_xyz,size(Data,3)*vox_xyz);

%               gridDATA  - 3D logical array of size (P,Q,R) - Voxelised data
%                                     1 => Inside the object
%                                     0 => Outside the object
%
%               gridX     - A 1xP array       - List of the X axis coordinates.
%               gridY     - A 1xQ array       - List of the Y axis coordinates.
%               gridZ     - A 1xR array       - List of the Z axis coordinates.

CONVERT_voxels_to_stl('testl_H.stl',logical(Data),gridX,gridY,gridZ)
CONVERT_voxels_to_stl('testl_S.stl',~logical(Data),gridX,gridY,gridZ)

