% This will combine the DKI and FBI data into one 4-D file that can then be preprocessed (denoising, etc.)
% Also output is the gradient and bvals file that you will need.
%
% for pre-processing you will need MRTrix3.
% FSL implemented within MATLAB is also a requirement.

root = fullfile('/Volumes','Sandy','fbi_internal'); % root is the parent directory holding the subject
b = 6000; % b is the large b-value of the scan, in s/mm^2, usually b > 4000 s/mm^2
subjectList = {'IAM_1116'}; % list of the subjects in the study to be used


for s = 1:size(subjectList,2) 
    
    subject = subjectList{s}; % subject is the name of the folder and prefix of the files

    fprintf('%s...\n',subject)
    
    dataPath = fullfile(root,subject,'nifti','dwi');
    fbwmPath = fullfile(root,subject,'fbwm');
    
    if ~isdir(fbwmPath); mkdir(fbwmPath); end
    
    cd(dataPath);

    eval(['!fslmerge -t fbwm.nii ' sprintf('%s_18_DKI_30_Dir.nii',subject) ' ' sprintf('%s_17_FBI_b%d_128.nii',subject,b)]);
    movefile('fbwm.nii',fullfile(fbwmPath,'fbwm.nii'));
    
    % dki_img = spm_read_vols(spm_vol(fullfile(root,subject,sprintf('%s_dki.nii',subject))));
    dki_bval = load(sprintf('%s_18_DKI_30_Dir.bval',subject));
    dki_bvec = load(sprintf('%s_18_DKI_30_Dir.bvec',subject));

    % fbi_img = spm_read_vols(spm_vol(fullfile(root,subject,sprintf('%s_b%d.nii',subject,b))));
    fbi_bval = load(sprintf('%s_17_FBI_b%d_128.bval',subject,b));
    fbi_bvec = load(sprintf('%s_17_FBI_b%d_128.bvec',subject,b));

    BVAL = [dki_bval fbi_bval]';
    BVEC = [dki_bvec fbi_bvec]';

    cd(fbwmPath);
    
    G = 'fbwm.bvec'; % filename for new coregistered gradient table
    fid = fopen(G,'w+'); % open grad_fn in write-mode
    fprintf(fid,'%18.15f\t%18.15f\t%18.15f\n',BVEC'); % write new gradient directions to file GT
    fclose(fid); % close GT file

    B = 'fbwm.bval'; % filename for new coregistered gradient table
    fid = fopen(B,'w+'); % open grad_fn in write-mode
    fprintf(fid,'%d\n',BVAL); % write new gradient directions to file BT
    fclose(fid); % close BT file
    
    eval('!fslroi fbwm.nii b0.nii 0 1');
    eval('!bet b0.nii b0_brain.nii -f 0.3 -g 0.2 -m');
    eval('!fslmaths b0_brain_mask.nii -eroF b0_brain_mask_eroF.nii');
    
end
