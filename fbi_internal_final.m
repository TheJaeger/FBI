% Perform fiber ball imaging (FBI) and FBI white matter model (FBWM) analyses
% Author(s): Russell Glenn, Emilie McKinnon and Hunter Moss
% Medical University of South Carolina (MUSC)
% Modified further by Hunter Moss and Emilie McKinnon (2016, 2017, 2018, 2019, 2020)
%--------------------------------------------------------------------------
%
%   REFERENCES
%   
%   Moss, H. G., et al. (2019) Optimization of data acquisition and analysis for fiber ball imaging. NeuroImage, 200, 670-703.
%
%   McKinnon, E. T., Helpern, J. A., & Jensen, J. H. (2018). Modeling white matter microstructure with fiber ball imaging. NeuroImage, 176, 11-21.
%
%   Jensen, J. H., Glenn, G. R., & Helpern, J. A. (2016). Fiber ball imaging. Neuroimage, 124, 824-833.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: for the output to make sense, currently need a whole-brain mask to work with. 
%       If the input is some subset (e.g., a WM binary mask), the FBWM metrics are wrong.
%       The reason for this is currently unclear to me. 
%
% NOTE: The make_fib flag is under construction still...keep it set to 0 for OFF.
%
%
% NOTE: for fODF rectification: re-expansion needs to be the same degree (lmax) as the initial DWI spherical harmonic expansion.
%       If it is not, the cost function does strange things and the outputs are clearly incorrect. Again, I am unsure why this \
%       is the case.
%
%
%       -Hunter Moss (04/30/2020)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: USER INPUT   
%
FBI_struc.subj = 'IAM_1116'; % subject name
FBI_struc.study = 'fbi_internal'; % study folder (i.e, TDE, HARDI, IAM, etc.)
FBI_struc.root = fullfile('/Volumes','Sandy',FBI_struc.study,FBI_struc.subj,'fbwm'); % main directory where data is held
FBI_struc.bval = 6;     %b-value ms/microm^2 you want to use; must be in b-vals table; usually 4, 5 or 6
FBI_struc.degree = 6; % can be an even degree up to and including 20; l = 6 is defualt, could possibly go up to l = 8.
FBI_struc.degree_rec = FBI_struc.degree; % can be an even degree up to and including 20; l = 6 is defualt. (Should be the same as degree currently)
FBI_struc.fn_mask = fullfile(FBI_struc.root,'b0_brain_mask_eroF.nii'); % path to mask of voxels to calculate FBI parameters 
FBI_struc.path_data = fullfile(FBI_struc.root,'fbwm.nii'); % path to the pre-processed 4D data (i.e., denoised, GRC, Rician corrected, Gaussian smoothed, Eddy and motion corrected, etc.)
FBI_struc.path_gradient = fullfile(FBI_struc.root,'fbwm.bvec'); % path to the gradient text file
FBI_struc.path_bval = fullfile(FBI_struc.root,'fbwm.bval'); % path to the b-value text file
FBI_struc.outdir = fullfile(FBI_struc.root,'output'); % path to directory to place the output images and files
FBI_struc.pre_name = ''; % pre-pend string to output images, usually subject name (i.e., sprintf('%s_',subj))
FBI_struc.post_name = ''; % post-pend string to output images, usually b-value and degree (i.e., sprintf('_b%d000_l%d',bval,degree))
FBI_struc.fbwm_flag = 1; % perform (1) or not (0) the FBI white matter (WM) model (FBWM)
FBI_struc.rectification = 1; % perform (1) or not (0) the fODF rectification
FBI_struc.make_fib = 0; % Write .fib output for dsi studio (1 := yes, 0:= no)
FBI_struc.path_dke_parameters = '/Volumes/Sandy/fbi_internal/DKEParameters.txt'; % path to DKE parameters file.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: WMFT FIELDS
FBI_struc.ft_flag = 1; % perform (1) or not (0) WMFT with EULER FT (deterministic).
FBI_struc.fa_threshold = 0.1; % if DKI derived FA then 0.1, if FBI derived FAA then 0.3 should be okay.
FBI_struc.angle_threshold = 35; % angular threshold between steps to terminate tracking.
FBI_struc.trk_length = 30; % in mm
FBI_struc.step_size = 0.1; % in mm
FBI_struc.trk_mask = FBI_struc.fn_mask; % binary mask to allow fiber tracking to traverse; typically, brain mask.
FBI_struc.seed_flag = 2; % binary mask where FT seeding will occur (randomly), (1) whole-brain, (2) FA > 0.2 or (3) FAA > 0.3
FBI_struc.name = sprintf('%s_FBI_b%d000',FBI_struc.subj,FBI_struc.bval); % name to define the structure that is output (e.g., sprintf('%s_FBI_b%d000',subj,bval))
FBI_struc.seednum = 50000; % number of random seed points to generate; typically, 100k-250k; more takes longer.
FBI_struc.hdr = fullfile(FBI_struc.outdir,[FBI_struc.pre_name 'FAA' FBI_struc.post_name '.nii']); % typicall FA or FAA image hdr to be read by EULER_DKE for tracking; should match whatever fa_threshold will be using for tracking.
% END: WMFT FIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           ----- CAUTION --- WARNING --- CAUTION --- 
%
%       careful consideration should be taken 
%      before modifiying below even these following line.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: OPTIONAL FIELDS: DEFAULTS ARE RECOMMENDED CURRENTLY -----
FBI_struc.odf_size = 20000; %Note: length(mask_idx)<=odf_size*11 per DSI studio fib file creation
FBI_struc.find_orientations = 1; %Do you care about the directions of local maxima? It takes lots of extra time if you do (default := 1)
FBI_struc.lambda = 0; %Laplace-Beltrami regularization parameter; eg 0.006 (never changed from 0, not sure what it affects: HM)
FBI_struc.Dn = 1.5; %Estimate of intra-neurite diffusivity to estimate fns (um^2/ms)
FBI_struc.D0 = 3.0;  %Diffusivity of free water (um^2/ms) at body temp (37 degree C), FYI: this is the upper bound on Da
FBI_struc.scale_odf = 0.1; % only affects DSI studio visualization
FBI_struc.fbvc = 1; % apply finite b-value correction (default := 1), usually want to do this but we see no significant difference if not performed, but it doesn't take much extra time.
FBI_struc.img_orientation = 'LAS'; % Orientation of image volumes
FBI_struc.odf_orientation = 'LAS'; % Orientation of gradient table
FBI_struc.peak_threshold = 0.4; % threshold to include peaks for WMFT; found 0.4 of max peak per voxel was sufficient; CSD for instance uses 0.1;
% END: OPTIONAL FIELDS 
% END: USER INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------ DO NOT MODIFY BELOW THIS LINE --------
%
subj = FBI_struc.subj;
study = FBI_struc.study;
bval = FBI_struc.bval;
root = FBI_struc.root;
dataIn = FBI_struc.path_data;
btIn = FBI_struc.path_bval;
gtIn = FBI_struc.path_gradient;    
odf_size = FBI_struc.odf_size;
pre_name = FBI_struc.pre_name;
post_name = FBI_struc.post_name;
degree = FBI_struc.degree;
degree_rec = FBI_struc.degree_rec;
find_orientations = FBI_struc.find_orientations; 
lambda = FBI_struc.lambda;             
Dn = FBI_struc.Dn;               
D0 = FBI_struc.D0;                
scale_odf = FBI_struc.scale_odf;        
apply_fbvc = FBI_struc.fbvc;
fn_mask = FBI_struc.fn_mask;
img_orientation = FBI_struc.img_orientation;
odf_orientation = FBI_struc.odf_orientation;
make_fib = FBI_struc.make_fib;
peak_threshold = FBI_struc.peak_threshold;
outdir = FBI_struc.outdir;
%
%--------------------------------------------------------------------------
% Image orientation stuff...don't change
img_orientation = lower(img_orientation); 
odf_orientation = lower(odf_orientation);

permute_img = [[strfind(img_orientation,'l') strfind(img_orientation,'r')],[strfind(img_orientation,'a') strfind(img_orientation,'p')],[strfind(img_orientation,'s') strfind(img_orientation,'i')]]; 
img_orientation = img_orientation(permute_img); 
invert_img = 2*[strcmpi(img_orientation(1),'l'), strcmpi(img_orientation(2),'p'), strcmpi(img_orientation(3),'s')]-1; 

permute_odf = [[strfind(odf_orientation,'l') strfind(odf_orientation,'r')],[strfind(odf_orientation,'a') strfind(odf_orientation,'p')],[strfind(odf_orientation,'s') strfind(odf_orientation,'i')]]; 
img_orientation = img_orientation(permute_odf); 
invert_odf = 2*[strcmpi(odf_orientation(1),'l'), strcmpi(odf_orientation(2),'p'), strcmpi(odf_orientation(3),'s')]-1; 


fprintf('\n Fiber ball imaging (FBI) analysis \n')

% read in b-values and gradient table 
fprintf('\t Loading b-value and gradient tables...\n')
bt = load(btIn); % b-value table
gt = load(gtIn); % gradient table (should be N X 3)

bt_unique = unique(bt); % check if there is multiple unique b-values

% Separate b0 from DWI's
% idx = find(bt ~= 0); % b0 indexing
idx = find(bt == bval*10^3); % FBI b-value indexing
idx_b0 = find(bt == 0);  % b0 indexing

% make output directory if it doesn't already exist
if~isdir(outdir); mkdir(outdir); end

% Read in data
fprintf('\n \t Reading in data...\n')
img = spm_read_vols(spm_vol(dataIn)); % data to be processed

%Get image information
fprintf('\t Getting image information...\n')
hdr = spm_vol(fn_mask); 
dimension = hdr.dim; % image dimensions
voxel_size = sqrt(sum(hdr.mat(1:3).^2)); % voxel size 

fprintf('\n \t subject: %s \n \t b-value: %d000 s/mm^2 \n \t ndirs: %d \n',subj,bval,length(idx))

% determine if the SH expansion degree selected is even and throw error if not 
if mod(degree,2); error('Only even degree SH expansions are allowed'); end 

fprintf('\t SH expansion: l = %d \n',degree)

% Define what spherical harmonics to grab depending on the degree set by the user
% Only even degrees are allowed...
l_num = 2*[0:2:degree] + 1;
harmonics = [];
sh_end = 1;

for h = 1:length(0:2:degree)
    
    sh_start = sh_end + l_num(h) - 1;
    sh_end = sh_start + l_num(h) - 1;  
    harmonics = [harmonics sh_start:sh_end];
        
end

if FBI_struc.rectification == 1
    
    l_num = 2*[0:2:degree_rec] + 1;
    harmonics_rec = [];
    sh_end = 1;

    for h = 1:length(0:2:degree_rec)

        sh_start = sh_end + l_num(h) - 1;
        sh_end = sh_start + l_num(h) - 1;  
        harmonics_rec = [harmonics_rec sh_start:sh_end];

    end

end

% Separate the non-DWI's and the DWI's
b0 = img(:,:,:,idx_b0); % non-DWI's 
% average the b0's if needed (doesn't cause problems if not needed)
if length(idx_b0) > 1; b0 = mean(b0,4); end

% Load brain mask in 
fprintf('\n \t Loading binary brain mask...\n')
mask = spm_read_vols(hdr);
% mask_idx = reshape(1:prod(dimension),dimension).*mask; 
% mask_idx = mask_idx(mask_idx > 0);
mask_reshape = reshape(mask,[1,prod(dimension)]);

if length(bt_unique) == 4 && FBI_struc.fbwm_flag == 1 % if there are, and FBWM flag is set, then find the images
    
    if ~isdir(fullfile(outdir,'DKE')); mkdir(fullfile(outdir,'DKE')); end

    fprintf('\n\t FBI white matter (FBWM) model is selected...\n')
    fprintf('\t Getting files organized for FBWM...\n')
    
    if (bt_unique(2)/1000) < 1; else bt_unique = bt_unique/1000; end % check and put b-values into ms/um^2

    idx1 = find(bt == bt_unique(2)*10^3); 
    idx2 = find(bt == bt_unique(3)*10^3);

    % Get gradient vectors and b-value tables for each b-value 
    gt1 = gt(idx1,:); % b1000 gradient table
    gt2 = gt(idx2,:); % b2000
    gt = gt(idx,:); 
    
    % I think EM did this as a normalization ?? Need to ask her.
    gt1 = gt1./sqrt(gt1(:,1).^2 + gt1(:,2).^2 + gt1(:,3).^2);
    gt2 = gt2./sqrt(gt2(:,1).^2 + gt2(:,2).^2 + gt2(:,3).^2);
    gt = gt./sqrt(gt(:,1).^2 + gt(:,2).^2 + gt(:,3).^2);
    
    bt1 = bt(idx1); % b1000 b-value table
    bt2 = bt(idx2); % b2000 b-value table
    
    fprintf('\t Writing gradient files...\n')
    % write gradient files for DKE 
    fn = fullfile(outdir,'DKE','gradient_b1.txt'); 
    fid = fopen(fn,'w');   
    fprintf(fid,'%18.15f\t%18.15f\t%18.15f\n',gt1'); 
    fclose(fid);     
    
    fn = fullfile(outdir,'DKE','gradient_b2.txt'); 
    fid = fopen(fn,'w');   
    fprintf(fid,'%18.15f\t%18.15f\t%18.15f\n',gt2'); 
    fclose(fid); 
       
    fprintf('\n\t Writing DKI data: 4D_DKI.nii \n\n')
    % DKI image data
    img1 = img(:,:,:,idx1); % b = 1000 s/mm^2 DWI's
    img2 = img(:,:,:,idx2); % b = 2000 s/mm^2 DWI's
    
    img_dki = cat(4,b0,img(:,:,:,[idx1 idx2])); % make the DKI 4D image file by combining the b0, b1000 and b2000 data in the 4th dimension
    hdr.fname = fullfile(outdir,'DKE','4D_dki.nii'); % create header for DKI 4D image
    hdr_dki = repmat(hdr,[1,length(idx1)+length(idx2)+1]); % make one for each volume to be written
    
    % Write out the 4D DKI image
    for ii = 1:(length(idx1)+length(idx2)+1)
        hdr_dki(ii).n = [ii 1];
        spm_write_vol(hdr_dki(ii),img_dki(:,:,:,ii).*mask);
    end
    
    img1_temp = zeros(length(idx1),prod(dimension));
    img2_temp = zeros(length(idx2),prod(dimension));
    
    for i = 1:length(idx1); img1_temp(i,:) = reshape(img1(:,:,:,i),[1,prod(dimension)]); end
    for i = 1:length(idx2); img2_temp(i,:) = reshape(img2(:,:,:,i),[1,prod(dimension)]); end
    
    img1 = img1_temp; % reshaped b = 1000 s/mm^2 data
    img2 = img2_temp; % reshaped b = 2000 s/mm^2 data
    
    clear img1_temp img2_temp
          
    fprintf('\n\t Starting Diffusional Kurtosis Estimator (DKE)...\n\n')
    
    fid = fopen(FBI_struc.path_dke_parameters); % original DKE paramters file that will be modified
    fout = fullfile(outdir,'DKE','DKEParameters.txt'); % modified DKE parameters file 
    fidout = fopen(fout,'w');
    
    while(~feof(fid))
        
        s = fgetl(fid);
        s = strrep(s,'''/Users/Example/output/gradient.txt''',['{' '''' outdir '/DKE/gradient_b1.txt'',''' outdir '/DKE/gradient_b2.txt''}']);
        s = strrep(s,'/Users/Example/output/DKE',[outdir '/DKE']);
        s = strrep(s,'ndir = 30',['ndir = ' num2str(length(idx1))]);
        s = strrep(s,'bval = [0 1000 2000]',['bval = [0 ' num2str(bt_unique(2)*10^3) ' ' num2str(bt_unique(3)*10^3) ']']);

        fprintf(fidout,'%s\n',s);
%         disp(s)
    end
    
    fclose(fid);
    fclose(fidout);
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    dke([outdir '/DKE/DKEParameters.txt'])
    delete(poolobj);
    
    fprintf('\n\t Loading in diffusion tensor (DT)...\n')
    load([outdir '/DKE/DT.mat'])
    
    fprintf('\t Getting SH basis functions for DKI...\n')

    if FBI_struc.rectification == 1
        
        % Generate SH basis functions for the DKI b-values 
        B1 = getSH(degree_rec,[atan2(gt1(:,2),gt1(:,1)) acos(gt1(:,3))],'complex'); 
        B1 = B1(:,harmonics_rec);
        B1(:,~harmonics) = 0;
        
        
        B2 = getSH(degree_rec,[atan2(gt2(:,2),gt2(:,1)) acos(gt2(:,3))],'complex'); 
        B2 = B2(:,harmonics_rec);
        B2(:,~harmonics) = 0;
        
    else
        
        B1 = getSH(degree,[atan2(gt1(:,2),gt1(:,1)) acos(gt1(:,3))],'complex'); 
        B1 = B1(:,harmonics);
        
        B2 = getSH(degree,[atan2(gt1(:,2),gt1(:,1)) acos(gt1(:,3))],'complex'); 
        B2 = B2(:,harmonics);
        
    end
    
    GT = {gt1, gt2, gt};
    BT = bt_unique(2:end);
    ndir = [length(idx1) length(idx2) length(idx)];
    
else
    
    gt = gt(idx,:);
    gt = gt./sqrt(gt(:,1).^2 + gt(:,2).^2 + gt(:,3).^2);
     
end

% mask_reshape = mask_reshape.*wm;

b0 = reshape(b0,[1,prod(dimension)]);

img = img(:,:,:,idx); % high b-value DWI's
img_temp = zeros(length(idx),prod(dimension));
for i = 1:length(idx); img_temp(i,:) = reshape(img(:,:,:,i),[1,prod(dimension)]); end
img = img_temp; % reshaped high b-value data
clear img_temp

%Initialize Things---------------------------------------------------------
[odf_vertices odf_faces] = odf_vertices_faces_DSIStudio_8; %Fib Output
vidx = 1:size(odf_vertices,2)/2; 
V = [acos(odf_vertices(3,vidx))' atan2(odf_vertices(2,vidx),odf_vertices(1,vidx))']; 
    
[S,IDX,~,AREA] = sphericalgrid5;  % Peak Detection grid
R = [sin(S(:,1)).*cos(S(:,2)),sin(S(:,1)).*sin(S(:,2)),cos(S(:,1))];  % Cartesian Coordinates
    
% Calculate the SH basis functions
% For FBI, these are complex valued
fprintf('\t Calculating SH basis functions...\n')
B = getSH(degree,[atan2(gt(:,2),gt(:,1)) acos(gt(:,3))],'complex'); 
B = B(:,harmonics);

H = getSH(degree,[S(:,2) S(:,1)],'complex');
H = H(:,harmonics);

HV = getSH(degree,[atan2(odf_vertices(2,vidx),odf_vertices(1,vidx))' acos(odf_vertices(3,vidx))'],'complex');
HV = HV(:,harmonics);

if FBI_struc.rectification == 1
    
    Hr = getSH(degree_rec,[S(:,2) S(:,1)],'complex');
    Hr = Hr(:,harmonics_rec);

    Br = getSH(degree_rec,[atan2(gt(:,2),gt(:,1)) acos(gt(:,3))],'complex'); 
    Br = Br(:,harmonics_rec);

end

% ADDITIONAL MATRICES---------------------------------------------------

if FBI_struc.rectification == 1
  
    degree_term = degree_rec;
    
else
    
    degree_term = degree;
    
end

idx_Y = 0;

 for l = 0:2:degree_term
    
    Pl0(idx_Y+1:idx_Y+(2*l+1),:) = ((-1)^(l/2) * factorial(l))/(4^(l/2)*factorial(l/2)^2)*ones(2*l+1,1); % Jensen (2016) Eq. 12
    gl(idx_Y+1:idx_Y+(2*l+1),:) = (factorial(l/2)*(bval*D0)^((l+1)/2)/gamma(l+3/2)*hypergeom(((l+1)/2),l+3/2,-bval*D0))*ones(2*l+1,1); % Jensen (2016) Eq. 18 
    idx_Y = idx_Y + 2*l+1;
    
 end

%-------------------------------------------------------------------------- 

% Initialize output parameters
fprintf('\t Initializing output arrays...\n')
odf_f = cell(1,prod(dimension)); % fODF directions
odf_fp = cell(1,prod(dimension)); % fODF peak amplitudes

SH = zeros(length(harmonics),prod(dimension)); % spherical harmonic coefficients
SH_reshape = zeros([dimension,length(harmonics)]);
zeta = zeros(1,prod(dimension)); % tissue prameter defined as f/sqrt(Da), see Jensen (2016)
FAA = zeros(1,prod(dimension)); % intra-axonal fractional anisotropy, see McKinnon (2018)

if FBI_struc.fbwm_flag == 1
    
    De = zeros(3,3,prod(prod(dimension)));
    aDT = zeros(6,prod(dimension)); % inta-axonal DT, linear
    cost_fn = zeros(100,prod(dimension));

    iDT_img = zeros(3,3,prod(prod(dimension(1:3))));
    iaDT_img = zeros(3,3,prod(prod(dimension(1:3))));

    De_mean=zeros(1,prod(dimension(1:3)));
    De_ax=zeros(1,prod(dimension(1:3)));
    De_rad=zeros(1,prod(dimension(1:3)));
    De_fa=zeros(1,prod(dimension(1:3)));

end

if FBI_struc.rectification == 1
    
    SH_rec = zeros(length(harmonics_rec),prod(dimension));

end

%PROCESS EVERYTHING--------------------------------------------------------
%Go through all voxels in the brain masks in steps based on the odf_size
fprintf('\t Processing...\n')
tic
for vox = 2174 %1:prod(dimension)

    b0n = b0(vox); % b0 image data subset
    imgn = img(:,vox); % high b-value data subset
    
    if mask_reshape(vox) == 1 && ~isequal(b0n,0)
        
        % Get the SH coefficients
        alm = (B'*B)^-1*B'*imgn/b0n;    %SH fit to normalized signal attenuation (S/S0)
        alm(isnan(alm)) = 0; % find any NaN values (from denoising usually) and set them to zero
        a00 = alm(1,:); % for normalization
    
        clm = alm*gl(1).*(sqrt(4*pi)*Pl0*alm(1).*gl).^-1; % equation 8 FBWM paper (McKinnon 2018) and Eq. 21, Jensen (2016)
        c00 = clm(1,:); % first fODF coefficient
        Y00 = 1/(sqrt(4*pi)); % first SH basis value
    
        % Normalize the fODF coefficients
        clm = bsxfun(@rdivide,clm,c00); % first by C00
        clm = bsxfun(@times,clm,Y00); % now to Y00
        SH(:,vox) = clm; % store the clm's for use in MRTrix3
    
        %Get fODFs
        ODF = H*clm; % for peak detection (and potentially, if flagged, fODF rectification) 
 
        if FBI_struc.rectification == 1
            
            % fODF rectification
            % Here we perform a bisection search for the root f(c) == 0, 
            % to determine the epsilon that will get rid of negative fODF peaks.

            fODF = real(ODF); % the fODF 
            fODF(isnan(fODF)) = 0;  
            Fmax = max(fODF); % get the max value that will bound the fODF between [0 and Fmax]

            lB = 0; % initial lower bound
            uB = Fmax; % initial upper bound

            M = 1; % initialized number of iterations
            Mmax = 100; % max number of iterations

            if ~isequal(Fmax,0)
                while M <= Mmax 

                    midpt = (lB + uB)/2; % the new midpoint for the search
                    fODF_lB = sum((abs(fODF - lB) - fODF - lB).*AREA,1); % calculate f(a)
                    fODF_midpt = sum((abs(fODF - midpt) - fODF - midpt).*AREA,1); % calculate f(c)

                    if fODF_midpt == 0 || (uB - lB)/2 < 10^-8 % tolerance criteria

                        EPS = midpt; % output epsilon as EPS

                        break % break the while loop

                    else

                        M = M + 1; % increase iteration by 1 step

                        if isequal(sign(fODF_midpt),sign(fODF_lB)) % determine sign since a and b must be of different sign for this to work

                            lB = midpt; % if sign(f(midpt)) == sign(f(lB))

                        else 

                            uB = midpt;

                        end

                    end                           
                end

                % Apply the rectification epsilon value to the ODF    
                % the rectification does not affect our previous normalziation as it affects all points/peaks equally
                % for this verison of the fODF, calling sum(AREA.*ODF) will give a value of approximately 0.50
                % because the points that create ODF are sampled from only a half-shpere 

                ODF = (1/2)*(abs(ODF - EPS) + ODF - EPS);

            end

            % Since the bi-section algorithm is only going to get epsilon to be very close to making the integral zero, we sort the rest by zeroing out what is left 
            % this is a tiny tiny fractional amount that we are correcting here...
            
            ODF(ODF > -10^-8 & ODF < 0) = 0;
            ODF(ODF < 10^-8 & ODF > 0) = 0;

            % re-expand the fODF coefficients 
            clm_rec = (AREA.*ODF)'*conj(Hr);
            clm_rec = permute(clm_rec,[2 1]);
            
            clm = clm_rec; % redefine the coefficient variable from clm_rec to clm for continuity
            c00 = clm(1,:); % for normalization
            clm = bsxfun(@rdivide,clm,c00); % normalize the first SH coefficient
            clm = bsxfun(@times,clm,Y00); % scales it now to 1/sqrt(4pi)

            SH_rec(:,vox) = clm; % store the SH coefficients of the fODF for saving in a 4D nifti for MRTrix3
    
        end
        
        aDT(:,vox) = axonal_DT(clm); % calculate the intra-axonal DT

        % Rotationally invariant tissue paramter calculations
        zeta(vox) = a00.*sqrt(bval)./pi; % Jensen (2016), Eq. 14 
        FAA(vox) = sqrt(3*sum(abs(clm(2:6).^2))./(5*abs(clm(1)).^2 + 2*sum(abs(clm(2:6).^2)))); % intra-axonal FA using complex SH coefficients, see McKinnon (2018)
     
        ODF = real(ODF);
        % Peak selection for WMFT
        peaks = zeros(size(IDX,1),size(ODF,2));
      
        if find_orientations            
            for j = 1:size(IDX,1)

                x = ODF(IDX(j,:));  
                peaks(j,:) = x(1,:) == max(x);

            end

            for j = 1:size(peaks,2)

                idxp = find(peaks(:,j)); % Get peak index

                if ~isempty(idxp) % added since there were zero peaks being indexed here for some reason

                    p = ODF(IDX(idxp,1));          % Get peak magnitudes
                    [p,idxj] = sort(p,'descend');   % sort the peaks in descending order

                    % This will grab only the top 3 peaks for WMFT
                    if size(p,1) > 3

                        p = p(1:3,1);
                        idxj = idxj(1:3,1);

                    end

                    odf_f{vox} = R(idxp(idxj),:)'; % Save sorted peaks
                    odf_fp{vox} = [p' min(ODF)];  % Save Peak magnitudes

                end

            end
        end
        
    if FBI_struc.fbwm_flag == 1
        
        if FBI_struc.rectification == 1
            
            shB = {B1,B2,Br};
            
            g2l_fa_R_b = zeros(length(bt_unique(2:end)),100, length(harmonics_rec));
            g2l_fa_R = zeros(length(harmonics_rec),100);
            g2l_fa_R_large = zeros(length(harmonics_rec),100);
            
            degree_term = degree;
            
        else

            shB = {B1,B2,B}; 
            
            g2l_fa_R_b = zeros(length(bt_unique(2:end)),100, length(harmonics));
            g2l_fa_R = zeros(length(harmonics),100);
            g2l_fa_R_large = zeros(length(harmonics),100);
            
            degree_term = degree_rec;

        end
            
        imgn1 = img1(:,vox); % b1000 image data subset
        imgn2 = img2(:,vox); % b2000 image data subset 

        IMG = {imgn1, imgn2, imgn};
        
        f_grid = linspace(0,1,100);

        iaDT = [aDT(1,vox), aDT(4,vox), aDT(5,vox);
                aDT(4,vox), aDT(2,vox), aDT(6,vox);
                aDT(5,vox), aDT(6,vox), aDT(3,vox)];

         iDT = [DT(1,vox), DT(4,vox), DT(5,vox);
                DT(4,vox), DT(2,vox), DT(6,vox);
                DT(5,vox), DT(6,vox), DT(3,vox)];

        for b = 1:length(BT)

            idx_hyper = BT(b)*(f_grid.^2)*zeta(vox)^-2 < 20; % when bda becomes large use equation 20 FBI paper                
            idx_Y = 0;

            for l = 0:2:degree_term

                hypergeom_opt = sum((gamma((l+1)/2+(0:99)).*gamma(l+3/2).*((-BT(b)*(f_grid(idx_hyper).^2*zeta(vox)^-2))'.^(0:99))./(factorial((0:99)).*gamma(l+3/2+(0:99)).*gamma((l+1)/2))),2)';
%                 hypergeom_opt = sum((gamma((l+1)/2+(0:100)).*gamma(l+3/2).*((-BT(b)*(f_grid(idx_hyper).^2*zeta(vox)^-2))'.^(0:100))./(factorial((0:100)).*gamma(l+3/2+(0:100)).*gamma((l+1)/2))),2)';
                g2l_fa_R(idx_Y+1:idx_Y+(2*l+1), idx_hyper) = repmat((factorial(l/2)*(BT(b).*(f_grid(idx_hyper).^2*zeta(vox)^-2)).^((l+1)/2)/gamma(l+3/2).*hypergeom_opt),2*l+1,1); % equation 9 FBWM paper 
                idx_Y = idx_Y + 2*l+1;

            end

            g2l_fa_R_b(b,idx_hyper,:) = g2l_fa_R(:, idx_hyper)';
            idx_Y = 0;

            for l = 0:2:degree_term

                g2l_fa_R_large(idx_Y+1:idx_Y+(2*l+1), ~idx_hyper) = repmat((exp(-l/2*(l+1)./((2*BT(b).*(f_grid(~idx_hyper).^2*zeta(vox)^-2))))),2*l+1,1); %equation 20 FBI paper 
                idx_Y = idx_Y + 2*l+1;

            end

            g2l_fa_R_b(b,~idx_hyper,:) = g2l_fa_R_large(:, ~idx_hyper)';

        end

        
        % grid search 
        for grid = 1:100 
             for b = 1:length(BT)

                awf = f_grid(grid);
                Se = (b0n.*exp((-BT(b)*(1-awf).^-1).*diag((GT{b}*(iDT-(awf.^3*zeta(vox).^-2)*iaDT)*GT{b}'))))*(1-awf); %equation 3 FBWM
                Sa = (ones(length(GT{b}),1).*(2*pi*b0n*zeta(vox)*sqrt(pi./BT(b)))).*(shB{b}*(Pl0.*squeeze(g2l_fa_R_b(b,grid,:)).*clm)); % equation 15 FBWM
                
                cost_fn(grid,vox) = cost_fn(grid,vox) + ndir(b)^-1*sum((IMG{b}-Se-Sa).^2);

             end


             cost_fn(grid,vox)=b0n^-1*sqrt(length(BT)^-1*cost_fn(grid,vox)); % equation 21


             iDT_img(:,:,vox) = iDT;
             iaDT_img(:,:,vox) = iaDT;

        end
         
    end
    
    end
end


if FBI_struc.fbwm_flag == 1
    
    L_1 = zeros(1,prod(dimension));
    L_2 = zeros(1,prod(dimension));
    L_3 = zeros(1,prod(dimension));
    
    % Calculate the tissue microstructure parameters
    [min_cost_fn,min_cost_fn_idx] = sort(cost_fn); % find index where cost function is minimal
    awf_grid = linspace(0,1,100);
    min_awf = awf_grid(min_cost_fn_idx(1,:)); % which grid point is this?
    Da = min_awf.^2./zeta.^2; % Equation 22 FBWM paper (McKinnon 2018)

    for vox = 1:prod(dimension)
        if mask_reshape(vox) == 1 && ~isequal(b0n,0)


            De(:,:,vox) = (iDT_img(:,:,vox)-(min_awf(vox).^3.*zeta(vox).^-2).*iaDT_img(:,:,vox))./(1-min_awf(vox)); %equation 19 FBWM

            int = De(:,:,vox);
            int(isnan(De(:,:,vox))) = 0; %get rid of NaNs
            int(isinf(De(:,:,vox))) = 0; %get rid of INFs

            [~, L] = eig(int);
            L = diag(L);
            L = sort(L, 'descend');
            N = 1;

               while (L(1)<0) || (L(2)<0) || (L(3)<0) % Find different AWF value if the eigenvalues are < 0 
                    N = N+1;
                    if N < 100

                        min_awf(vox) = awf_grid(min_cost_fn_idx(N,vox));      

                    else % if all awf values result in negative eigenvalues set AWF to 0

                        min_awf(vox) = 0;
                        De(:,:,vox) = (iDT_img(:,:,vox)-(min_awf(vox).^3.*zeta(vox).^-2).*iaDT_img(:,:,vox))./(1-min_awf(vox));
                        Da(vox) = min_awf(vox).^2./zeta(vox).^2;

                        break
                    end

                De(:,:,vox) = (iDT_img(:,:,vox)-(min_awf(vox).^3.*zeta(vox).^-2).*iaDT_img(:,:,vox))./(1-min_awf(vox));
                Da(vox) = min_awf(vox).^2./zeta(vox).^2;

               end

            % calculate eigenvalues again with the correct AWF value 

            int = De(:,:,vox);
            int(isnan(De(:,:,vox))) = 0; %get rid of NaNs
            int(isinf(De(:,:,vox))) = 0; %get rid of INFs

            [~, L] = eig(int);  
            L = diag(L);
            L = sort(L, 'descend');   
            De_ax(vox) = L(1); % equation 24 FBWM
            De_rad(vox) = sum(L(2:3))*0.5; %equation 25 FBWM
            De_fa(vox) = sqrt(((L(1) - L(2)) ^ 2 + (L(1) - L(3)) ^ 2 + (L(2) - L(3)) ^ 2) / (2 * sum(L .^ 2))); 
            De_mean(vox) = 1/3*(2*De_rad(vox)+De_ax(vox)); % equation 23 FBWM ( sort of ) 
            
            L_1(i)=L(1);
            L_2(i)=L(2);
            L_3(i)=L(3); 

        end
    end
    
    De(isnan(De)) = 0;
    De(isnan(De)) = 0;

end

toc

fprintf('\t Writing output images...\n')
     
odf_f = reshape(odf_f,dimension);
odf_fp = reshape(odf_fp,dimension);
zeta = reshape(zeta,dimension); 
zeta(zeta < 0) = 0; % contrain zeta to be greater than 0 since if it is < 0 it is unphysical
zeta(zeta > 1) = 1; % contrain zeta to be less than 1 since if > 1 it is unphysical
FAA = reshape(FAA,dimension);

if FBI_struc.fbwm_flag == 1

    min_awf = reshape(min_awf,dimension);
    Da = reshape(Da,dimension);
    De_mean = reshape(De_mean,dimension);
    De_ax = reshape(De_ax,dimension);
    De_rad = reshape(De_rad,dimension);
    De_fa = reshape(De_fa,dimension);
    min_cost_fn = reshape(min_cost_fn(1,:),dimension);

end

%Save Things---------------------------------------------------------------

hdr.dt = [16 0]; % needed by default for DWI image header
hdr.fname = fullfile(outdir,[pre_name 'zeta' post_name '.nii']); spm_write_vol(hdr,zeta);  
hdr.fname = fullfile(outdir,[pre_name 'FAA' post_name '.nii']); spm_write_vol(hdr,FAA);

if FBI_struc.fbwm_flag == 1
    
    hdr.fname = fullfile(outdir,[pre_name 'awf_fbwm' post_name '.nii']); spm_write_vol(hdr,min_awf);
    hdr.fname = fullfile(outdir,[pre_name 'Da_fbwm' post_name '.nii']); spm_write_vol(hdr,Da);
    hdr.fname = fullfile(outdir,[pre_name 'De_ax_fbwm' post_name '.nii']); spm_write_vol(hdr,De_ax);
    hdr.fname = fullfile(outdir,[pre_name 'De_rad_fbwm' post_name '.nii']); spm_write_vol(hdr,De_rad);
    hdr.fname = fullfile(outdir,[pre_name 'De_mean_fbwm' post_name '.nii']); spm_write_vol(hdr,De_mean);
    hdr.fname = fullfile(outdir,[pre_name 'De_FA_fbwm' post_name '.nii']); spm_write_vol(hdr,De_fa);
    hdr.fname = fullfile(outdir,[pre_name 'min_cost_fn_fbwm' post_name '.nii']); spm_write_vol(hdr,min_cost_fn);

    save(fullfile(outdir,[pre_name 'De.mat']),'De');
    save(fullfile(outdir,[pre_name 'iaDT.mat']),'iaDT_img');
    save(fullfile(outdir,[pre_name 'iDT.mat']),'iDT_img');
    save(fullfile(outdir,[pre_name 'cost.mat']),'cost_fn');
     
end

save(fullfile(outdir,[pre_name 'fODF_SH' post_name '.mat']),'SH');
save(fullfile(outdir,[pre_name 'odf_f' post_name '.mat']),'odf_f');
save(fullfile(outdir,[pre_name 'odf_fp' post_name '.mat']),'odf_fp');


% This will write-out the SH in the order (convention) of MRTrix3 so that they an be visualized or used for FT with CSD.
% Haven't yet worked out what the SH_idx is for l >= 8 yet.

if degree <= 8
    
    if isequal(degree,6)
        
        SH_idx = [3 5 8 10 12 14 17 19 21 23 25 27]; 
        
     elseif isequal(degree,8)
         
        SH_idx = [3 5 8 10 12 14 17 19 21 23 25 27 30 32 34 36 38 40 42 44];
        
     end

     SH_reshape = zeros([dimension,length(harmonics)]);


     for i=1:length(harmonics)

        SH_reshape(:,:,:,i) = reshape(SH(i,:),dimension);

     end
     
     
    hdr_sh = spm_vol(fullfile(outdir,[pre_name 'zeta' post_name '.nii'])); 
    sh_hdr = repmat(hdr_sh,[length(harmonics),1]);

    SH_reshape(:,:,:,SH_idx) = SH_reshape(:,:,:,SH_idx) * -1; % this is to account for the Condon-Shortley phase that is missing in MRTrix3

    hdr.fname = fullfile(outdir,[pre_name 'fODF_SH' post_name '.nii']);
    hdr_4D=repmat(hdr,[1,length(harmonics)]);

    for ii=1:length(harmonics)
     hdr_4D(ii).n=[ii 1];
     spm_write_vol(hdr_4D(ii),SH_reshape(:,:,:,ii));
    end
    
%    % This is to convert the DKE DT.mat to DT.nii for pyFBI use and testing...
%     DTre = zeros([dimension,6]);
%     
%     for i = 1:6
%         DTre(:,:,:,i) = reshape(DT(i,:),dimension);
%     end
%     
%     hdr.fname = fullfile(outdir,'DKE',[pre_name 'DT' post_name '.nii']);
%     hdr_4D=repmat(hdr,[1,size(DTre,4)]);
% 
%     for ii=1:size(DTre,4)
%      hdr_4D(ii).n=[ii 1];
%      spm_write_vol(hdr_4D(ii),DTre(:,:,:,ii));
%     end

end

if FBI_struc.ft_flag == 1
    
    fprintf('\t Performing WMFT...\n\n\n')
    % Deterministic white matter fiber tractography (dWMFT)
    FT_struc.fa_threshold = FBI_struc.fa_threshold; 
    FT_struc.angle_threshold = FBI_struc.angle_threshold; 
    FT_struc.trk_length = FBI_struc.trk_length; % in mm
    FT_struc.step_size = FBI_struc.step_size; % in mm
    FT_struc.trk_mask = FBI_struc.trk_mask; 
    
    if FBI_struc.seed_flag == 1
        
        FT_struc.seed_mask = FBI_struc.fn_mask;
        
    elseif FBI_struc.seed_flag == 2
        
        fa = spm_read_vols(spm_vol(fullfile(outdir,'DKE','fa.nii')));
        fa_mask = fa > 0.2;
        hdr.fname = fullfile(outdir,'DKE',[pre_name 'fa_mask' post_name '.nii']); spm_write_vol(hdr,fa_mask);
        eval(['!fslmaths ' hdr.fname ' -fillh26 ' hdr.fname]);
        
        FT_struc.seed_mask = hdr.fname;
        
    elseif FBI_struc.seed_flag == 3
        
        FAA_mask = FAA > 0.5;
        hdr.fname = fullfile(outdir,[pre_name 'FAA_mask' post_name '.nii']); spm_write_vol(hdr,FAA_mask);
        eval(['!fslmaths ' hdr.fname ' -fillh26 ' hdr.fname]);
        FT_struc.seed_mask = fullfile(outdir,'FAA_mask.nii');
        
    end
    
    FT_struc.shift = 0.5; 
    FT_struc.name = FBI_struc.name; 
    FT_struc.seednum = FBI_struc.seednum; 
    FT_struc.permute_odf = permute_odf;  
    FT_struc.invert_odf = invert_odf; 
    FT_struc.permute_img = permute_img; 
    FT_struc.invert_img = invert_img; 
    FT_struc.hdr = spm_vol(FBI_struc.hdr);
    FT_struc.outdir = outdir; 
    FT_struc.reset_memory = 1; 
    FT_struc.output_DTI_trks = 0;  
    FT_struc.pre_name = pre_name; 
    FT_struc.post_name =  post_name; 
    SEED = random_seed_FT(FT_struc);
    FT_struc.SEED = SEED;
    FT_struc.vox_to_ras = hdr.mat([FT_struc.permute_img 4],:)*[diag(FT_struc.invert_img) double(FT_struc.invert_img==-1)'.*FT_struc.hdr.dim(FT_struc.permute_img)'-0.5*[1;1;1];0 0 0 1];
    save(fullfile(outdir,sprintf('%s_FBI_FT_struc_b%d000.mat',subj,bval)),'FT_struc');

    tic

    EULER_DKE(FT_struc,odf_f,'FBI_FT')

    toc

end
