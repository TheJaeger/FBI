%% "FBWM": fiber ball white matter model     
    %       input:
    %           - path_data: path to 4D nifti 
    %           - path_gradient:   path to gradient directions in same order as data
    %           - path_bval:  path to bvalues in same order as data
    %           - degree: degree of spherical harmonics used for fODF 
    %           - path_brain_mask: path to binary mask 
    %           - path_output: path to output data
    %           - path_dke_parameters= path to template file DKEParameters.txt

    %   Authors: Emilie McKinnon (mckinnon@musc.edu), August 31, 2018
    %   Cleaner: Hunter Moss (mossh@musc.edu), November 21, 2019
    %   Copyright (c) 2019 medical university of south carolina (MUSC)

%--------------------------------------------------------------------------
%   REFERENCES
%   
%   Moss, H. G., et al. (2019) Optimization of data acquisition and analysis for fiber ball imaging. NeuroImage, 200, 670-703.
%
%   McKinnon, E. T., Helpern, J. A., & Jensen, J. H. (2018). Modeling white matter microstructure with fiber ball imaging. NeuroImage, 176, 11-21.
%
%   Jensen, J. H., Glenn, G. R., & Helpern, J. A. (2016). Fiber ball imaging. Neuroimage, 124, 824-833.
%
%

%% input 
% study = FBI_struc.study;
% subj = FBI_struc.subj;
% root = FBI_struc.root;
% path_data = FBI_struc.path_data;
% path_gradient = FBI_struc.path_gradient;
% path_brain_mask = FBI_struc.fn_mask;
% degree = FBI_struc.degree;
% path_bval = FBI_struc.path_bval;
% path_output = FBI_struc.path_output ;
% path_dke_parameters = fullfile('/Users','gonzo','Documents','MATLAB','dMRI_updated','DKEParameters.txt');

% For testing purposes...
FBI_struc.subj = 'IAM_1116'; % subject name
FBI_struc.study = 'fbi_internal'; % study folder (i.e, TDE, HARDI, IAM, etc.)
FBI_struc.root = fullfile('/Volumes','Sandy',FBI_struc.study,FBI_struc.subj,'fbwm'); % main directory where data is held
FBI_struc.bval = 6;     %b-value ms/microm^2 you want to use; must be in b-vals table; usually 4, 5 or 6
FBI_struc.degree = 6; % can be an even degree up to and including 20; l = 6 is defualt.
FBI_struc.degree_rec = 14; % can be an even degree up to and including 20; l = 6 is defualt.
FBI_struc.fn_mask = fullfile(FBI_struc.root,'b0_brain_mask_eroF.nii'); % path to mask of voxels to calculate FBI parameters 
FBI_struc.path_data = fullfile(FBI_struc.root,'fbwm.nii'); % path to the pre-processed 4D data (i.e., denoised, GRC, Rician corrected, Gaussian smoothed, Eddy and motion corrected, etc.)
FBI_struc.path_gradient = fullfile(FBI_struc.root,'fbwm.bvec'); % path to the gradient text file
FBI_struc.path_bval = fullfile(FBI_struc.root,'fbwm.bval'); % path to the b-value text file
FBI_struc.outdir = fullfile(FBI_struc.root,'output'); % path to directory to place the output images and files
FBI_struc.pre_name = ''; % pre-pend string to output images, usually subject name (i.e., sprintf('%s_',subj))
FBI_struc.post_name = ''; % post-pend string to output images, usually b-value and degree (i.e., sprintf('_b%d000_l%d',bval,degree))
FBI_struc.fbwm_flag = 1; % perform (1) or not (0) the FBI white matter (WM) model (FBWM)
FBI_struc.path_dke_parameters = '/Volumes/Sandy/fbi_internal/DKEParameters.txt'; % path to DKE parameters file.
path_dke_parameters = FBI_struc.path_dke_parameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BEGIN: WMFT FIELDS
FBI_struc.fa_threshold = 0.1; % if DKI derived FA then 0.1, if FBI derived FAA then 0.3 is okay.
FBI_struc.angle_threshold = 35; % angular threshold between steps to terminate tracking.
FBI_struc.trk_length = 30; % in mm
FBI_struc.step_size = 0.1; % in mm
FBI_struc.trk_mask = FBI_struc.fn_mask; % binary mask to allow fiber tracking to traverse; typically, brain mask.
FBI_struc.seed_mask = FBI_struc.fn_mask; % binary mask where FT seeding will occur (randomly), usually either FA > 0.2 mask or brain mask.
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
FBI_struc.odf_orientation = 'LPS'; % Orientation of gradient table
FBI_struc.peak_threshold = 0.4; % threshold to include peaks for WMFT; found 0.4 of max peak per voxel was sufficient; CSD for instance uses 0.1;
% END: OPTIONAL FIELDS 
% END: USER INPUT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------ DO NOT MODIFY BELOW THIS LINE --------
%
subj = FBI_struc.subj;
study = FBI_struc.study;
% bval = FBI_struc.bval;
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

peak_threshold = FBI_struc.peak_threshold;
outdir = [FBI_struc.outdir '/Ems_code'];
path_output = outdir;
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
%%
% read in bvalues and gradient table 
bt= textscan(fopen(fullfile(btIn)),'%f');
gt= textscan(fopen(fullfile(gtIn)),'%f%f%f');

% define bvalues used
bval=unique(bt{1});
bval_fbi=bval(end);
bval_dki=bval(2:end-1); 

index_b1= find(bt{1}==bval_dki(1));
index_b2= find(bt{1}==bval_dki(2));
index_fbi= find(bt{1}==bval_fbi);

% number of gradient directions for each b-value 
ndir_b1=length(index_b1) ;
ndir_b2=length(index_b2) ;
ndir_fbi=length(index_fbi);
ndir=[ndir_b1 ndir_b2 ndir_fbi];

% check units of bvalues and put it in ms/um2 
if (bval_dki/1000)<1
else
bval_dki=bval_dki/1000;
bval_fbi=bval_fbi/1000; 
end
bval=[bval_dki' bval_fbi'];

%Separate b0 from DWIs
idx_b0 = find(bt{1}==0); 

% Read in data
S=spm_read_vols(spm_vol(fullfile(dataIn)));
dim=size(S);

S_b1=zeros(1,prod(dim(1:3)));
S_b2=zeros(1,prod(dim(1:3)));
Sb6_reshape=zeros(1,prod(dim(1:3)));

S_b1_4D=S(:,:,:,index_b1);
S_b2_4D=S(:,:,:,index_b2);

for i=1:ndir_b1 
    S_b1(i,:)=reshape(S_b1_4D(:,:,:,i),[1,prod(dim(1:3))]);
end

for i=1:ndir_b2
    S_b2(i,:)=reshape(S_b2_4D(:,:,:,i),[1,prod(dim(1:3))]);
end
    
S_b1_b6=S(:,:,:,index_fbi);
for i=1:length(index_fbi)
    Sb6_reshape(i,:)=reshape(S_b1_b6(:,:,:,i),[1,prod(dim(1:3))]);
end

S_b0=S(:,:,:,idx_b0);
S_b0_mean=mean(S_b0,4); 
s0=reshape(S_b0_mean,[1, prod(dim(1:3))]);

% brain mask 
hdr=spm_vol(fn_mask);
brainmask=spm_read_vols(hdr);
brain_mask_reshape=reshape(brainmask,[1,prod(dim(1:3))]);

% seperate gradient tables 
GT_fbi = [gt{1}(index_fbi) gt{2}(index_fbi) gt{3}(index_fbi)];  
GT_fbi=GT_fbi./sqrt(GT_fbi(:,1).^2+GT_fbi(:,2).^2+GT_fbi(:,3).^2);

GT_b1 = [gt{1}(index_b1) gt{2}(index_b1) gt{3}(index_b1)];  
GT_b1=GT_b1./sqrt(GT_b1(:,1).^2+GT_b1(:,2).^2+GT_b1(:,3).^2);
GT_b1=GT_b1(1:ndir_b1,:);

GT_b2 = [gt{1}(index_b2) gt{2}(index_b2) gt{3}(index_b2)];  
GT_b2=GT_b2./sqrt(GT_b2(:,1).^2+GT_b2(:,2).^2+GT_b2(:,3).^2);
GT_b2=GT_b2(1:ndir_b2,:);

GT={GT_b1,GT_b2,GT_fbi};

%GET SH BASIS FUNCTIONS
B_fbi=getSH(degree,[ atan2(GT_fbi(:,2),GT_fbi(:,1)) acos(GT_fbi(:,3))],'complex');
B_b1=getSH(degree,[ atan2(GT_b1(:,2),GT_b1(:,1)) acos(GT_b1(:,3))],'complex');
B_b2=getSH(degree,[ atan2(GT_b2(:,2),GT_b2(:,1)) acos(GT_b2(:,3))],'complex');

Harm_id={1,5:9,17:25,37:49,65:81}; % max degree = 8 
B_fbi=B_fbi(:,cell2mat(Harm_id(1:degree/2+1)));
B_b1=B_b1(:,cell2mat(Harm_id(1:degree/2+1)));
B_b2=B_b2(:,cell2mat(Harm_id(1:degree/2+1)));

B={B_b1,B_b2,B_fbi};


%% Calculate diffusion tensor 
mkdir([outdir '/DKE/'])
S_DKE=cat(4,S_b0_mean,S_b1_4D,S_b2_4D);
hdr.fname=[path_output '/DKE/4D_dki.nii'];
hdr_4D=repmat(hdr,[1,ndir_b1+ndir_b2+1]);

% write the 4D nifti file for DKE 

for ii=1:(ndir_b1+ndir_b2+1)
   hdr_4D(ii).n=[ii 1];
   spm_write_vol(hdr_4D(ii),S_DKE(:,:,:,ii).*brainmask);
end

% write gradient file for DKE 
    fn = fullfile(path_output,'gradient_b1.txt'); 
    fid = fopen(fn,'w');   
    fprintf(fid,'%18.15f\t%18.15f\t%18.15f\n',GT_b1'); 
    fclose(fid);     
    
    fn = fullfile(path_output,'gradient_b2.txt'); 
    fid = fopen(fn,'w');   
    fprintf(fid,'%18.15f\t%18.15f\t%18.15f\n',GT_b2'); 
    fclose(fid);     
    
% adjust dkeparameters template file    
fid=fopen(path_dke_parameters); % original DKE_parameters
fout=fullfile(path_output,'DKE','DKEParameters.txt'); % new ft_parameters 
fidout=fopen(fout,'w');

while(~feof(fid))
    s=fgetl(fid);
    s=strrep(s,'''/Users/Example/output/gradient.txt''',['{' '''' path_output '/gradient_b1.txt'',''' path_output '/gradient_b2.txt''}']);
    s=strrep(s,'/Users/Example/output/',[path_output '/']);
    s=strrep(s,'ndir = 30',['ndir=[' num2str(ndir_b1) ',' num2str(ndir_b2) ']']);
    s=strrep(s,'bval=[0 1000 2000]',['bval=[0' num2str(bval_dki'*1000)]);

    fprintf(fidout,'%s\n',s);
    disp(s)
end

fclose(fid);
fclose(fidout);

poolobj = gcp('nocreate');
delete(poolobj);
dke([path_output '/DKE/DKEParameters.txt'])
load([path_output '/DKE/DT.mat'])
% load([root '/DT.mat']);


%% Core calculations (can take about 10 min) 
% initialisation 
tic 
De=zeros(3,3,prod(prod(dim(1:3))));
intermediate_A_img=zeros(3,3,prod(prod(dim(1:3))));
Cost=zeros(100,prod(dim(1:3)));
intermediate_DT_img=zeros(3,3,prod(prod(dim(1:3))));
zeta_img=zeros(1,prod(dim(1:3)));
faa_f=zeros(1,prod(dim(1:3)));

% calculate legendre functions and g2l 
D0=3;
idx_Y = 0;
    for l = 0:degree/2
        P2l0(idx_Y+1:idx_Y+(4*l+1), :) =  ((-1)^(l) * factorial(2*l))/(4^(l)*factorial(l)^2)*ones(4*l+1,1); %equation 12 FBI paper
        G2l(idx_Y+1:idx_Y+(4*l+1), :) = (factorial(l)*(bval_fbi*D0)^((l+0.5))/gamma(2*l+3/2)*hypergeom((l+0.5),2*l+3/2,-bval_fbi*D0))*ones(4*l+1,1); %equation18 FBI paper
        idx_Y = idx_Y + 4*l+1;
    end
    
% voxel-wise estimation of other parameters 
for vox=1:prod(dim(1:3)) % turned parfor off bc of error being thrown about the inaccesable 'DT'?

G2l_fa_R_b=zeros(length(bval),100, length(cell2mat(Harm_id(1:degree/2+1))));
G2l_fa_R=zeros(length(cell2mat(Harm_id(1:degree/2+1))),100);
G2l_fa_R_large=zeros(length(cell2mat(Harm_id(1:degree/2+1))),100);

%voxelwise signal values
S0_vox=s0(vox);
Sfbi_vox=Sb6_reshape(:,vox);
Sb1_vox=squeeze(S_b1(:,vox));
Sb2_vox=squeeze(S_b2(:,vox));
S={Sb1_vox, Sb2_vox, Sfbi_vox};
    if brain_mask_reshape(vox)==1

     
A2l=(B_fbi'*B_fbi)^-1*B_fbi'*Sfbi_vox/s0(vox); %equation 4 FBWM paper 
C2l=A2l*G2l(1).*(sqrt(4*pi)*P2l0*A2l(1).*G2l).^-1; %equation 8 FBWM paper 

C2l = bsxfun(@rdivide,C2l,C2l(1)); %normalization
C2l= bsxfun(@times,C2l,1/(sqrt(4*pi))); %normalization
zeta=A2l(1)*sqrt(bval_fbi)./pi;

scale=repmat((((C2l(1)*sqrt(30))).^-1),[1 6]); % scale factor A 

A=scale.*[sqrt(30)/3*C2l(1)-sqrt(6)/3*C2l(4)+C2l(6)+C2l(2),...
            sqrt(30)/3*C2l(1)-sqrt(6)/3*C2l(4)-C2l(6)-C2l(2),...
            sqrt(30)/3*C2l(1)+2*sqrt(6)/3*C2l(4),...
            1i*C2l(6)-1i*C2l(2),...
            -C2l(5)+C2l(3),...
            -1i*C2l(5)-1i*C2l(3)]; %equation 12 FBWM paper 



intermediate_DT=[DT(1,vox), DT(4,vox), DT(5,vox);
                 DT(4,vox), DT(2,vox), DT(6,vox);
                 DT(5,vox), DT(6,vox), DT(3,vox)];

intermediate_A=[A(1), A(4), A(5);
                A(4), A(2), A(6);
                A(5), A(6), A(3)];

faa_f(vox)= sqrt(3*sum(abs(C2l(2:6)').^2)./(5*abs(C2l(1)').^2 + 2*sum(abs(C2l(2:6)').^2))); % equation 13 FBWM 

awf_grid=linspace(0,1,100);
        
            for b=1:length(bval)
                idx_hyper=bval(b)*(awf_grid.^2)* zeta^-2 <20; % when bda becomes large use equation 20 FBI paper                
                idx_Y = 0;
                for l = 0:degree/2
                hypergeom_opt=sum((gamma((l+0.5)+(0:100)).*gamma(2*l+3/2).*((-bval(b)*(awf_grid(idx_hyper).^2*zeta^-2))'.^(0:100))./(factorial((0:100)).*gamma(2*l+3/2+(0:100)).*gamma((l+0.5)))),2)';
                G2l_fa_R(idx_Y+1:idx_Y+(4*l+1), idx_hyper)=repmat((factorial(l)*(bval(b).*(awf_grid(idx_hyper).^2*zeta^-2)).^(l+0.5)/gamma(2*l+3/2).*hypergeom_opt),4*l+1,1); % equation 9 FBWM paper 
                idx_Y = idx_Y + 4*l+1;
                end
               
                G2l_fa_R_b(b,idx_hyper,:)=G2l_fa_R(:, idx_hyper)';

                idx_Y = 0;
                for l = 0:degree/2
                G2l_fa_R_large(idx_Y+1:idx_Y+(4*l+1), ~idx_hyper)=repmat((exp(-l*(2*l+1)./((2*bval(b).*(awf_grid(~idx_hyper).^2*zeta^-2))))),4*l+1,1); %equation 20 FBI paper 
                idx_Y = idx_Y + 4*l+1;
                end
                G2l_fa_R_b(b,~idx_hyper,:)=G2l_fa_R_large(:, ~idx_hyper)';
            end
 
% grid search 
for grid=1:100 
     for b=1:length(bval)
            awf=awf_grid(grid);
            Se=(S0_vox.*exp((-bval(b)*(1-awf).^-1).*diag((GT{b}*(intermediate_DT-(awf.^3*zeta.^-2)*intermediate_A)*GT{b}'))))*(1-awf); %equation 3 FBWM
            Sa=(ones(length(GT{b}),1).*(2*pi*S0_vox*zeta*sqrt(pi./bval(b)))).*(B{b}*(P2l0.*squeeze(G2l_fa_R_b(b,grid,:)).*C2l));       %equation 15 FBWM
            Cost(grid,vox)= Cost(grid,vox)+ndir(b)^-1*sum((S{b}-Se-Sa).^2);
     end
      Cost(grid,vox)=S0_vox^-1*sqrt(length(bval)^-1*Cost(grid,vox)); % equation 21
end

intermediate_DT_img(:,:,vox)=intermediate_DT;
zeta_img(vox)=zeta;
intermediate_A_img(:,:,vox)=intermediate_A;


    end
    
end
toc
%% Calculate other tissue microstructure parameters 
De_par=zeros(1,prod(dim(1:3)));
De_rad=zeros(1,prod(dim(1:3)));
De_mean=zeros(1,prod(dim(1:3)));
De_fa=zeros(1,prod(dim(1:3)));
De(isnan(De))=0;
De(isinf(De))=0;
L_1=zeros(1,prod(dim(1:3)));
L_2=zeros(1,prod(dim(1:3)));
L_3=zeros(1,prod(dim(1:3)));

[min_cost,idx]=sort(Cost); % find index where cost function is minimal 
awf_grid=linspace(0,1,100);
min_awf=awf_grid(idx(1,:)); % which grid point is this?
Da=min_awf.^2./zeta_img.^2; % Equation 22 FBWM 

for i=1:prod(dim(1:3))  
    if brain_mask_reshape(i)==1
    De(:,:,i)=(intermediate_DT_img(:,:,i)-(min_awf(i).^3.*zeta_img(i).^-2).*intermediate_A_img(:,:,i))./(1-min_awf(i)); %equation 19 FBWM
    
    int=De(:,:,i);
    int(isnan(De(:,:,i)))=0; %get rid of NaNs
    int(isinf(De(:,:,i)))=0; %get rid of INFs
    
    [~, L] = eig(int);
    L = diag(L);
    L = sort(L, 'descend');
    n=1;
    
       while   (L(1)<0) || (L(2)<0) || (L(3)<0) % Find different AWF value if the eigenvalues are < 0 
            n=n+1;
            if n<100
            min_awf(i)=awf_grid(idx(n,i));      
            else % if all awf values result in negative eigenvalues set AWF to 0
            min_awf(i)=0;
            De(:,:,i)=(intermediate_DT_img(:,:,i)-(min_awf(i).^3.*zeta_img(i).^-2).*intermediate_A_img(:,:,i))./(1-min_awf(i));
            Da(i)=min_awf(i).^2./zeta_img(i).^2;
            break
            end
        De(:,:,i)=(intermediate_DT_img(:,:,i)-(min_awf(i).^3.*zeta_img(i).^-2).*intermediate_A_img(:,:,i))./(1-min_awf(i));
        Da(i)=min_awf(i).^2./zeta_img(i).^2;
       end
       
    % calculate eigenvalues again with the correct AWF value 
 
    int=De(:,:,i);
    int(isnan(De(:,:,i)))=0; %get rid of NaNs 
    int(isinf(De(:,:,i)))=0; %get rid of INFs
    
    [~, L] = eig(int);  
    L = diag(L);
    L = sort(L, 'descend');   
    De_par(i)=L(1); % equation 24 FBWM
    De_rad(i)=sum(L(2:3))*0.5; %equation 25 FBWM
    De_fa(i)=sqrt(((L(1) - L(2)) ^ 2 + (L(1) - L(3)) ^ 2 + (L(2) - L(3)) ^ 2) / (2 * sum(L .^ 2))); 
    De_mean(i)=1/3*(2*De_rad(i)+De_par(i)); % equation 23 FBWM ( sort of ) 
    
    L_1(i)=L(1);
    L_2(i)=L(2);
    L_3(i)=L(3);    
          
    end
       
end
 
%% write outputs 
mkdir([path_output '/FBWM/']);
hdr.dt=[16 0];
zeta_reshape=reshape(zeta_img,dim(1:3));hdr.fname=[path_output '/FBWM/zeta_fbwm.nii'];spm_write_vol(hdr,zeta_reshape);
min_awf_reshape=reshape(min_awf,dim(1:3));hdr.fname=[path_output '/FBWM/awf_fbwm.nii'];spm_write_vol(hdr,min_awf_reshape);
Da_reshape=reshape(Da,dim(1:3));hdr.fname=[path_output '/FBWM/da_fbwm.nii']; spm_write_vol(hdr,Da_reshape);
De_par_reshape=reshape(De_par,dim(1:3));hdr.fname=[path_output '/FBWM/De_par_fbwm.nii']; spm_write_vol(hdr,De_par_reshape);
De_rad_reshape=reshape(De_rad,dim(1:3));hdr.fname=[path_output '/FBWM/De_rad_fbwm.nii']; spm_write_vol(hdr,De_rad_reshape);
De_fa_reshape=reshape(De_fa,dim(1:3));hdr.fname=[path_output '/FBWM/De_fa_fbwm.nii']; spm_write_vol(hdr,De_fa_reshape);
De_mean_reshape=reshape(De_mean,dim(1:3));hdr.fname=[path_output '/FBWM/De_mean_fbwm.nii']; spm_write_vol(hdr,De_mean_reshape);

min_cost_reshape=reshape(min_cost(1,:),dim(1:3));hdr.fname=[path_output '/FBWM/min_cost_fbwm.nii'];spm_write_vol(hdr,min_cost_reshape);
faa_f_reshape=reshape(faa_f,dim(1:3));hdr.fname=[path_output '/FBWM/faa_f_fbwm.nii'];spm_write_vol(hdr,faa_f_reshape);

L_1_reshape=reshape(L_1,dim(1:3));hdr.fname=[path_output '/FBWM/L_1_De_fbwm.nii'];spm_write_vol(hdr,L_1_reshape);
L_2_reshape=reshape(L_2,dim(1:3));hdr.fname=[path_output '/FBWM/L_2_De_fbwm.nii'];spm_write_vol(hdr,L_2_reshape);
L_3_reshape=reshape(L_3,dim(1:3));hdr.fname=[path_output '/FBWM/L_3_De_fbwm.nii'];spm_write_vol(hdr,L_3_reshape);

save([path_output '/FBWM/de.mat'],'De');
save([path_output '/FBWM/A.mat'],'intermediate_A_img');
save([path_output '/FBWM/Cost.mat'],'Cost');
save([path_output '/FBWM/intermediate_DT.mat'],'intermediate_DT_img');
save([path_output '/FBWM/intermediate_A_img.mat'],'intermediate_A_img');


