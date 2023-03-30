clear all
close all
clc

% Go to the folder of this file:
foldername = fileparts(matlab.desktop.editor.getActiveFilename);
if(~isdeployed)
    cd(foldername);
end
% Add the tools:
addpath('Mesh_voxelisation')

%% Color Formating:
TU_cyan = [0,102,162]/255; % Used for Hard
TU_white = [249,249,255]/255; % Used for soft
TU_red = [195,49,47]/255;

%% DECISIONS TO TAKE:
% 0. What's the size of each RVE in the xy direction (how 'chuncky' will the gradient be)
U_xyz = 8; % Edge size (in voxels) of the greyscale cube
% 1. What's the thickness (in mm!) of the gradient? (I'll correct that to RVEs in the lines after)
W_G = 4.2;
% 2. Over which direction do you want to make the gradient?
D_G = 1; % (0) Made over the soft material. (1) Made over the hard material. (2) Made over both materials.
% 3. What is the type of interface function?
fun_shape = 'Cos';%  'Sig' for sigmoid, 'Lin' for linear, 'Cos' for cosine

design_num = 1:1; % How many variations of the design you wanna test?

for i = 1 : length(design_num)
    
    % What's the name of the hard (H) and soft (S) files
    filename_H = [foldername,'\Finger_hard_v',num2str(i),'.stl'];
    filename_S = [foldername,'\Finger_soft_v',num2str(i),'.stl'];
    
    
    %% Voxel and RVE sizes: (You shouldn't touch this part)
    % The xy direction (For simplicity, the voxels will be squared plates)
    vox_xyz = 25.4/300; % Edge size of the cubic voxel
    RVE_xyz = vox_xyz*U_xyz; % Edge size of the RVE
    W_G = round(W_G/RVE_xyz);
    
    %% Importing the files:
    [G_S, G_H] = import_SoftHard_design(filename_S, filename_H, RVE_xyz, RVE_xyz); %G_S is the soft material 2D image, G_H is the hard one
    % G_S(end-5:end,end-5:end,end-5:end) = 0;
    
    % We plot the 2D designs
    G_color = grs2rgb3D(G_H + 0.001*G_S, TU_white, TU_cyan);
    figure; imshow3D(G_color); % Cyan is hard material, white is soft material, black is no material
    
    % CREATE GRADED INTERFACE
    [G_greyscale, rho_e] = create_design_layers(G_S, G_H, W_G, D_G, fun_shape);
    
    %% Plotting things:
    % First the function:
    X_e = [0:length(rho_e)-1].*RVE_xyz;
    figure; hold on;
    pbaspect([33  10 1]); hold on;  set(gca,'fontsize', 21)
    box off; ax = gca;ax.LineWidth = 2;
    grid off
    plot(X_e,100*rho_e,'-','Color',TU_red,'LineWidth',2); hold on;
    xlabel('X (mm)'); %Position
    ylabel('\rho (%)'); %Hard/soft ratio
    xlim([0,max(X_e)]);
    ylim([0,100]);
    
    % Now the final interface
    G_color = grs2rgb3D(G_greyscale, TU_white, TU_cyan);
    figure; imshow3D(G_color); % Cyan is hard material, white is soft material, black is no material
    
    %% Store the file for Abaquss
    G_greyscale(G_greyscale>0&G_greyscale<0.01) = 0.01;
    folder_abbas = [foldername,'\abbas'];
    name_struct = ['Finger_v',num2str(i),'_grey'];
    TPMS_write_AbaqusInput_v2(folder_abbas,name_struct,double(G_greyscale),...
        vox_xyz, vox_xyz, vox_xyz, 1, 1, 1);
    disp(['The FEM model will have ',num2str(round(sum(ceil(G_greyscale(G_greyscale>0))))), ' elements!'])
    
    %% Making it bitmap
    [G_S_bits, G_H_bits] = make_bitmap(G_greyscale, rho_e, G_S, G_H, U_xyz, U_xyz);
    
    G_color = grs2rgb3D(G_H_bits+G_S_bits*0.001, TU_white, TU_cyan);
    figure; imshow3D(G_color); % Cyan is hard material, white is soft material, black is no material
    % Storing example images of the top layer:
    figure; imshow(squeeze(G_color(:,:,1,:))); % Cyan is hard material, white is soft material, black is no material
    plot2svg(['C:\Scratch\Finger\',name_struct,'_col.svg']);
    
    
    %% Extending the designs for printing and making the files:
    % The printing sizes:
    vox_x = 25.4/300; % vox_x (along the rows) is 84 um
    vox_y = 25.4/600; % vox_y (along the rows) is 42 um
    vox_z = 0.027; % vox_z (along each layer) is 27 um
    
    % We need to repeat the elements than
    G_S_bits = repelem(G_S_bits,round(vox_xyz./vox_x),round(vox_xyz./vox_y),round(vox_xyz./vox_z));
    G_H_bits = repelem(G_H_bits,round(vox_xyz./vox_x),round(vox_xyz./vox_y),round(vox_xyz./vox_z));
    % We check the design. NOTE: if vox_xyz is cubic, this new image should
    % look "extended" towards the columns:
    figure; imshow( double(G_H_bits(:,:,1))*1.00 + double(G_S_bits(:,:,1))*0.25); %Just show the first image
    
    %% Making the print files:
    disp(name_struct)
    
    toprint = 1; % Set this constant to 1 if you want to generate printing files (it may take a while)
    if toprint == 1
        folder_prints = 'C:\Scratch\Finger\Prints'; % The folder that stores all the files:
        Mat1 = 'VeroCY-V';  Mat2 = 'VeroClear';  Mat3 = 'VeroPureWht';
        Mat4 = 'M.Cleanser';    Mat5 = 'VeroMGT-V';  Mat6 = 'Agilus30Clr';
        G_F_bits = 0*G_S_bits; % This is a 'fluid' file, right now we haven't added it so it's just empty
        
        make_print_files(folder_prints, name_struct, G_H_bits, G_S_bits, G_F_bits,...
            Mat1,Mat2,Mat3,Mat4,Mat5,Mat6);
        %CAREFUL, FOR THIS CODE TO WORK AS IT ITS, THE FOLLOWING SETTING NEEDS TO BE IN THE PRINTER:
        % Mat1 has to be M.Cleanser (FLUID)
        % Mat3 has to be Agilus30Clr (SOFT);
        % Mat6 has to be VeroCY-V (HARD);
        % ALL THE OTHER NAMES HAVE TO BE CORRECT IN THE FILE
        % ANY DEVIATION WILL NEED ADJUSTMENTS
        % ALSO, THE OUTPUT FILES NEED TO BE ADJUSTED TO THE CORRECT PATH IN THE 3D-PRINTER COMPUTER
        
    end
end
% All done:
disp('How nice, to feel nothing, and still get full credit for being alive.');

