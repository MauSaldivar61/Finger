clear all
close all
clc

% Go to the folder of this file:
foldername = fileparts(matlab.desktop.editor.getActiveFilename);
if(~isdeployed)
    cd(foldername);
end
foldername_files = 'C:\Scratch\Finger';
folder_prints = [foldername_files,'\Prints'];
folder_abbas = [foldername_files,'\abbas'];
% Add the tools:
addpath('Mesh_voxelisation')

%% Color Formating:
TU_cyan = [0,102,162]/255; % Used for Hard
TU_white = [249,249,255]/255; % Used for soft
TU_red = [195,49,47]/255;

%% DECISIONS TO TAKE:
% 0. What's the size of each RVE in the xy direction (how 'chuncky' will the gradient be)
U_xyz = 3; % Edge size (in voxels) of the greyscale cube
% 1. What's the thickness (in mm!) of the gradient? (I'll correct that to RVEs in the lines after)
W_G = 1;
% 2. Over which direction do you want to make the gradient?
D_G = 1; % (0) Made over the soft material. (1) Made over the hard material. (2) Made over both materials.
% 3. What is the type of interface function?
fun_shape = 'Cos';%  'Sig' for sigmoid, 'Lin' for linear, 'Cos' for cosine

% TPMS-FLUID options
fun_tpms = 'Gyroid'; % Gyroid, Diamond, Iwp, Bcc, Octo, Neovius, Schwarz

%NOTE: There's a density tpms loop integrated now. For each design you
%send, it'll create verssions for each of the following densities
rho_tpms_v = [0.33,0.46,0.58]; % Percentage of solid material within the TPMS   use 0.33, 0.46, 0.58

sheet_q = 1; % Is the equation for sheet- (1) or beam- (0) based TPMS?
tpms_UC_vox = 60; %number of voxels for each unit cell (UC) of the TPMS cell (e.g., ideally, keep it in multiples of 12)

design_num = 1:2; % How many variations of the design you wanna test?

for i = design_num
    % What's the name of the hard (H), soft (S), and Fluid (F, F1) files
    %   filename_H = [foldername_files,'\Finger_hard_v',num2str(i),'.stl'];
    filename_H = [foldername_files,'\Finger_senseHard_v',num2str(i),'.stl'];
    %   filename_S = [foldername_files,'\Finger_soft_v',num2str(i),'.stl'];
    filename_S = [foldername_files,'\Finger_senseBoundary_v',num2str(i),'.stl'];
    %   filename_F = [foldername_files,'\Finger_TnF_v',num2str(i),'.stl']; % to generate Soft-fluid part with porous
    filename_F = [foldername_files,'\Finger_senseTPMS_v',num2str(i),'.stl'];
    %   filename_F1 = [foldername_files,'\Finger_F_v',num2str(i),'.stl'];  % to generate 100% fluid
    filename_F1 = [foldername_files,'\Finger_senseLiquid_v',num2str(i),'.stl'];
    
    %% Voxel and RVE sizes: (You shouldn't touch this part)
    % The xy direction (For simplicity, the voxels will be squared plates)
    vox_xyz = 25.4/300; % Edge size of the cubic voxel
    RVE_xyz = vox_xyz*U_xyz; % Edge size of the RVE
    W_G = round(W_G/RVE_xyz);
    
    %% Importing the files:
    [G_S, G_H, G_F, G_F1] = import_SoftHard_design_v2(filename_S, filename_H, filename_F, filename_F1, RVE_xyz, RVE_xyz); %G_S is the soft material 2D image, G_H is the hard one
    % G_S(end-5:end,end-5:end,end-5:end) = 0;
    
    % We plot the designs
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
    name_struct = ['Finger_v',num2str(i),'_grey'];
    TPMS_write_AbaqusInput_v2(folder_abbas,name_struct,double(G_greyscale),...
        RVE_xyz, RVE_xyz, RVE_xyz, 1, 1, 1);
    disp(['The FEM model will have ',num2str(round(sum(ceil(G_greyscale(G_greyscale>0))))), ' elements!'])
    
    % Since the next step is making the bitmap and fluid, we add the
    % TPMS density loop here to make the code more efficient:
    for jj = 1 : length(rho_tpms_v)
        rho_tpms = rho_tpms_v(jj);
        name_struct = ['Finger_v',num2str(i),'_grey','_',fun_tpms,num2str(rho_tpms*100)];
        close all
        % Making it bitmap
        [G_S_bits, G_H_bits] = make_bitmap(G_greyscale, rho_e, G_S, G_H, U_xyz, U_xyz);
        G_F_bits = repelem( G_F , size(G_S_bits,1)./size(G_F,1), size(G_S_bits,2)./size(G_F,2), size(G_S_bits,3)./size(G_F,3) );
        G_F1_bits = repelem( G_F1 , size(G_S_bits,1)./size(G_F1,1), size(G_S_bits,2)./size(G_F1,2), size(G_S_bits,3)./size(G_F1,3) );
        
        % CREATE TPMS-FLUID PART
        rho_fun = rho_tpms*ones(1,size(G_F_bits,2)); % Get the density function
        W_t = ceil(size(G_F_bits,2)/tpms_UC_vox)*tpms_UC_vox*vox_xyz; % Get a dimension that fits the TPMS on the size of the entire finger
        H_t = ceil(size(G_F_bits,1)/tpms_UC_vox)*tpms_UC_vox*vox_xyz; % Get a dimension that fits the TPMS on the size of the entire finger
        L_t = ceil(size(G_F_bits,3)/tpms_UC_vox)*tpms_UC_vox*vox_xyz; % Get a dimension that fits the TPMS on the size of the entire finger
        UC_tpms = tpms_UC_vox*vox_xyz;
        %   col,     row,     lay,
        G_T_bits = permute( TPMS_generate_Voxels(    H_t,     W_t,     L_t,...
            UC_tpms, UC_tpms, UC_tpms,...
            vox_xyz, vox_xyz, vox_xyz, fun_tpms, rho_fun, sheet_q), [2,1,3]) ; % Generating the TPMS and arranging the orientation
        G_T_bits = G_T_bits(1:size(G_F_bits,1),1:size(G_F_bits,2),1:size(G_F_bits,3)); % Cropping the design to fit the finger dimensions
        
        G_S_bits =  G_S_bits |(G_T_bits&G_F_bits); % The Soft material is what it was already defined and the TPMS function within the block
        G_F_bits =(~G_T_bits & G_F_bits)|G_F1_bits; % The fluid is everything is not in the TPMS and is within the block
        
        % Making the colored figures for the show!
        G_color = grs2rgb3D(G_H_bits+G_S_bits*0.001, TU_white, TU_cyan);
        G_color_fluid = grs2rgb3D(G_F_bits, TU_white, TU_red) + G_color;
        figure; imshow3D(G_color); % Cyan is hard material, white is soft material, black is no material
        figure; imshow3D(G_color_fluid); % Cyan is hard material, white is soft material, black is no material
        % Storing example images of the top layer:
        figure; imshow(squeeze(G_color_fluid(:,:,1,:))); % Cyan is hard material, white is soft material, black is no material
        plot2svg([foldername_files,name_struct,'_F_col.svg']);
        
        % MAKING A VIDEO
        makeVideosSIGMA( repelem( G_color_fluid, ceil(720/size(G_color_fluid,2)),ceil(720/size(G_color_fluid,2)))  ,[foldername_files,name_struct,'_F_col'],12);
        
        %% Extending the designs for printing and making the files:
        % The printing sizes:
        vox_x = 25.4/300; % vox_x (along the rows) is 84 um
        vox_y = 25.4/600; % vox_y (along the rows) is 42 um
        vox_z = 0.027; % vox_z (along each layer) is 27 um
        
        % We need to repeat the elements than
        G_S_bits = repelem(G_S_bits,round(vox_xyz./vox_x),round(vox_xyz./vox_y),round(vox_xyz./vox_z));
        G_H_bits = repelem(G_H_bits,round(vox_xyz./vox_x),round(vox_xyz./vox_y),round(vox_xyz./vox_z));
        G_F_bits = repelem(G_F_bits,round(vox_xyz./vox_x),round(vox_xyz./vox_y),round(vox_xyz./vox_z));
        
        G_S_bits(G_S_bits&G_H_bits) = 0;
        G_F_bits(G_S_bits&G_F_bits) = 0;
        G_F_bits(G_H_bits&G_F_bits) = 0;
        
        % We check the design. NOTE: if vox_xyz is cubic, this new image should
        % look "extended" towards the columns:
        figure; imshow( double(G_H_bits(:,:,1))*1.00 + double(G_S_bits(:,:,1))*0.1 + double(G_F_bits(:,:,1))*0.5); %Just show the first image
        
        %% Making the print files:
        disp(name_struct)
        
        toprint = 1; % Set this constant to 1 if you want to generate printing files (it may take a while)
        if toprint == 1
            Mat1 = 'VeroCY-V';  Mat2 = 'VeroClear';  Mat3 = 'VeroPureWht';
            Mat4 = 'M.Cleanser';    Mat5 = 'VeroMGT-V';  Mat6 = 'Agilus30Clr';
            
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
end
% All done:
disp('How nice, to feel nothing, and still get full credit for being alive.');


