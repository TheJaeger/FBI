#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
% Perform fiber ball imaging (FBI) and FBI white matter model (FBWM) analyses
% Author(s): Russell Glenn, Emilie McKinnon and Hunter Moss (MATLAB code)
%            Hunter Moss (1st python 3 iteration, 05/05/2020)
%
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
%   Under review:  Moss, H.G. and Jensen, J. H., Optimized rectification of fiber orientation density fucntion (2020)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE: for the output to make sense, currently need a whole-brain mask to work with.
%       If the input is some subset (e.g., a WM binary mask), the FBWM metrics are wrong.
%       The reason for this is currently unclear to me.
%
% NOTE: The make_fib flag is under construction still...keep it set to 0 for OFF.
%       (ODF outputs are not supported just yet in python version, WIP, 05/05/2020)
%
% NOTE: for fODF rectification: re-expansion needs to be the same degree (lmax) as the initial DWI spherical harmonic expansion.
%       If it is not, the cost function does strange things and the outputs are clearly incorrect. Again, I am unsure why this \
%       is the case.
%
#
% FOR PYTHON VERSION ONLY:
#
% NOTE: This assumes that DKI analysis has already been run
%       and the DT is findable in a directory called 'DKE' in the parent folder
%       - Should be modified to fit with PyDesigner DKI output!!!!
%       - See lines 93-97
%
% NOTE: Need to add peak detection for WMFT
%
% NOTE: Need to add SH output for used with MRTrix3 WMFT (CSD algorithm) - python only (MATLAB version already does this)
%
% NOTE: automate brain extraction and masking? Currently must supply a mask
%
% NOTE: ODF outputs need to be added.
%
% NOTE: Need to automate dimension extraction (currenlty set at [74,74,42] for writing out images, see end of program)
%
%       -Hunter Moss (04/30/2020) - MATLAB complete iteration finished
%       -Hunter Moss (05/05/2020) - python 1st iteration finished
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
'''
#---------------------------------------------------------------------
# Package Management
#---------------------------------------------------------------------
import sys as sys
import os # mkdir
import os.path as op # path
import numpy as np # array, ndarray
import numpy.matlib as npm
import nibabel as nib
import scipy as sp
from scipy.io import loadmat
import math

# Define study/dataset specifics
subj = 'IAM_1116'
study = 'fbi_internal'
fn = 'fbwm'
root = op.join('/Volumes/Sandy/',study,subj,'pyFBI')

# assumes bval, bvec and image files have same leading name (variable: fn [see above])
dataIn = op.join(root, fn + '.nii')
btIn = op.join(root, fn + '.bval')
gtIn = op.join(root, fn + '.bvec')

# get a hdr for a template to write out images later
# need to automate this as well...
hdr = nib.load(op.join(root,'b0.nii'))

fn_mask = op.join(root,'b0_brain_mask_eroF.nii') # must provide brain mask currently (should automate)
outdir = op.join(root)

bval = 6
degree = 6
degree_rec = 6
Dn = 1.5
D0 = 3.0
pkT = 0.4
voxSize = 3.0

rectification_flag = 1 # Do (1) or don't do (0) fODF rectification
fbwm_flag = 1 # Do (1) or don't do (0) FBI WM model (FBWM)

if fbwm_flag == 1:

    tmp = nib.load(op.join(root,'DKE','DT.nii'))
    DT = np.array(tmp.dataobj)
    DT = np.reshape(DT,(np.prod(DT[:,:,:,0].shape),6), order = 'F')

print('\n Fiber ball imaging (FBI) analysis \n')
print('\t SUBJECT: ' + subj + '\n')

print('\n Loading b-value table...\n')
bt = np.rint(np.loadtxt(btIn)) # load the b-val table
bt_unique = np.unique(bt).astype('int') # get the unique b-values

print('\n Loading gradient table...\n')
gt = np.loadtxt(gtIn) # load in the gradient file

if fbwm_flag == 1:

    # if FBWM is slectted, then split the gradients for DKI out
    gt1 = gt[bt == 1000,:]
    gt2 = gt[bt == 2000,:]

gt = gt[bt == bval*10**3,:] # FBI gradients

print('\n Loading binary brain mask...')
# load in binary brain mask
tmp = nib.load(fn_mask)
mask = np.array(tmp.dataobj).astype(bool)
mask = np.reshape(mask,(np.prod(mask.shape)), order = 'F') # reshape it (vecotrize it?)

print('\n Loading image data...')
# Load in the DWI image data
tmp = nib.load(op.join(root,fn + '.nii'))
img = np.array(tmp.dataobj)

# split out the b0 images (should maybe flag)
b0 = img[:,:,:, bt == 0]
b0 = np.mean(b0,axis = 3)

if fbwm_flag == 1:

    # split out DKI shells if FBWM is selected

    img1 = img[:,:,:, bt == 1000]
    img2 = img[:,:,:, bt == 2000]

    img1 = np.transpose(np.reshape(img1,(np.prod(img1[:,:,:,0].shape),img1.shape[3]), order = 'F'))
    img2 = np.transpose(np.reshape(img2,(np.prod(img2[:,:,:,0].shape),img2.shape[3]), order = 'F'))

img = img[:,:,:, bt == bval*1000] # FBI image data

# reshape (vectorize) the b0 and DWI image data
b0 = np.reshape(b0,(np.prod(b0.shape)), order = 'F')
img = np.transpose(np.reshape(img,(np.prod(img[:,:,:,0].shape),img.shape[3]), order = 'F'))

# HERE will be the SH basis set calculation...
# get the harmonics we will need
degs = np.arange(degree + 1)

l_tot = 2*degs + 1 # total harmonics in the degree
l_num = 2 * degs[::2] + 1 # how many per degree (evens only)
harmonics = []
sh_end = 1 # initialize the SH set for indexing

for h in range(0,len(degs[::2])):
    sh_start = sh_end + l_num[h] - 1
    sh_end = sh_start + l_num[h] - 1
    harmonics.extend(np.arange(sh_start,sh_end+1))

# Currently, the rectificaiton and original MUST have same degree,
# not sure why, will have to dig deeper...
if rectification_flag == 1:

    degs = np.arange(degree_rec + 1)

    l_tot = 2*degs + 1
    l_num = 2 * degs[::2] + 1
    harmonics_rec = []
    sh_end = 1

    for h in range(0,len(degs[::2])):
        sh_start = sh_end + l_num[h] - 1
        sh_end = sh_start + l_num[h] - 1
        harmonics_rec.extend(np.arange(sh_start,sh_end+1))

# Define the azimuthal (phi) and polar(theta) angles for our spherical exapnsion
# using the experimentally defined gradients from the scanner
theta = np.arccos(gt[:,2])
phi = np.arctan2(gt[:,1],gt[:,0])

# This was from Russel (I made his mat file into a text file here)
# this would be used for the ODF exapnsion which hasnt' been implemented yet
spherical_grid = np.loadtxt(op.join(root,'spherical_grid.txt')) # this is only HALF-SPHERE!
S1 = spherical_grid[:,0] # theta, i think
S2 = spherical_grid[:,1] # phi, i think
AREA = spherical_grid[:,2] # need the area since it is impossible to get exact isotropic (uniform) sampling

# initialze the SH basis set variable
B = np.zeros((len(gt),np.sum(l_num)),dtype=np.complex,order = 'F') # This is for the dMRI SH expasion
H = np.zeros((len(S1),np.sum(l_num)),dtype=np.complex,order = 'F') # This is for ODF creation (once it happens)

if fbwm_flag == 1:

    # SH basis set for the DKI data, just in case...
    B1 = np.zeros((len(gt1),np.sum(l_num)),dtype=np.complex,order='F')
    B2 = np.zeros((len(gt2),np.sum(l_num)),dtype=np.complex,order='F')

cnt = 0

for n in degs:
    for m in range(-n,n+1):

        if (n % 2) == 0:

            # FBI SH basis set generation
            B[:,cnt] = sp.special.sph_harm(m,n,phi,theta)
            H[:,cnt] = sp.special.sph_harm(m,n,S2,S1)


            if fbwm_flag == 1: # DKI SH basis sets

                theta1 = np.arccos(gt1[:,2])
                phi1 =  np.arctan2(gt1[:,1],gt1[:,0])

                theta2 = np.arccos(gt2[:,2])
                phi2 =  np.arctan2(gt2[:,1],gt2[:,0])

                B1[:,cnt] = sp.special.sph_harm(m,n,phi1,theta1)
                B2[:,cnt] = sp.special.sph_harm(m,n,phi2,theta2)

            cnt = cnt + 1

idx_Y = 0
Pl0 = np.zeros((len(harmonics),1),order = 'F') # need Legendre polynomial Pl0
gl = np.zeros((len(harmonics),1),order = 'F') # calculate correction factor (see original FBI paper, Jensne 2016)

for l in degs[::2]:

    Pl0[idx_Y:idx_Y+(2*l+1),:] = (np.power(-1,l/2)* np.math.factorial(l)) / (np.power(4,l/2)*np.power(np.math.factorial(l/2),2))*np.ones((2*l+1,1))
    gl[idx_Y:idx_Y+(2*l+1),:] = (np.math.factorial(l/2)*np.power(bval*D0,(l+1)/2))/sp.special.gamma(l+3/2)*sp.special.hyp1f1((l+1)/2,l+3/2,-bval*D0)*np.ones((2*l+1,1))
    idx_Y = idx_Y + (2*l+1)

Pl0 = np.squeeze(Pl0)
gl = np.squeeze(gl)

# initialize the outputs
SH = np.zeros((len(harmonics),img.shape[1]), dtype = 'complex', order = 'F') # will hold the clm (SH coefficients)
#SH_reshape = np.zeros((img.shape[1]),len(harmonics)), dtype = 'complex',order = 'F') # This would be to read out SH image to load into MRTrix3
zeta = np.zeros(img.shape[1],order = 'F') # zeta (see original FBI paper)
faa = np.zeros(img.shape[1],order = 'F') # FAA (see McKinnon 2018, Moss 2019)

if rectification_flag == 1:

    # If rectificaion is selected: SH_rec will hold the rectified clm's (SH coefficients corected)
    SH_rec = np.zeros((len(harmonics_rec),img.shape[1]), dtype = 'complex',order = 'F')

if fbwm_flag == 1: # initialize FBWM metrics, if selected

    De = np.zeros((3,3,img.shape[1]), order = 'F')
    aDT = np.zeros((6,img.shape[1]), order = 'F')
    cost_fn = np.zeros((100,img.shape[1]), order = 'F')

    iDT_img = np.zeros((3,3,img.shape[1]), order = 'F')
    iaDT_img = np.zeros((3,3,img.shape[1]), order = 'F')

    De_mean = np.zeros(img.shape[1], order = 'F')
    De_ax = np.zeros(img.shape[1], order = 'F')
    De_rad = np.zeros(img.shape[1], order = 'F')
    De_fa = np.zeros(img.shape[1], order = 'F')

    # these three are used in the cost function of FBWM
    BT = bt_unique[1:]/1000; # holds the b-vlaues (usually, 1, 2 and usually either 4, 5 or 6)
    GT = [gt1,gt2,gt] # holds all b-vale gradient vectors
    ndir = [len(B1), len(B2),len(B)] # holds the number of directions for each b-values

# get the entire set of voxel integers as an array
voxList = np.arange(0,img.shape[1], dtype = int)

# loop over those voxels
for vox in voxList:

    b0n = b0[vox] # get b0 (non-diffusion weighted) value (i.e., starting point for decay)
    imgn = img[:,vox] # get the DWI value for each gradient direction

    # grab only brain voxels (need to figure out how to only do WM, the FBWM does not like that though currenlty)
    if mask[vox] == 1 and b0n > 0:

        # For references to alm and clm see FBI papers, they (alm and clm) are defined in all of them
        alm = np.dot(np.linalg.pinv(B),(imgn/b0n)) # DWI signal SH coefficients (these are complex)
        alm[np.isnan(alm)] = 0
        a00 = alm[0].real # the imaginary part is on the order of 10^-18 (this is for zeta)

        clm = alm*gl[0]*np.power(np.sqrt(4*np.pi)*alm[0]*Pl0*gl,-1) # fODF SH coefficients (these are complex)
        c00 = clm[0]
        clm = clm/c00
        clm = clm*(1/np.sqrt(4*np.pi))

        SH[:,vox] = clm

        # need to figure out how to do peak detection (on this variable and then read out odf structures like in MATLAB code)
        # only the real part would be read out but that would need to be done later on after the rectification process below
        ODF = np.matmul(H,clm)

        # Here is where the fODF rectifcation begins
        # It is simple and uses a bi-section algorithm to find the ONLY root
        # There is only 1 root (see Optimized Rectificaion paper, Moss/Jensen 2020)
        # This will eliminate all negative peaks (and some small spurious peaks that are noise induced)
        # Makes the fODF completely positive

        if rectification_flag == 1:

            # fODF rectification
            fODF = ODF.real # grab real part of the fODF
            fODF[np.isnan(fODF)] = 0
            Fmax = np.max(fODF) # get the max peak value of the ODF

            lB = 0 # initial lower bound
            uB = Fmax # initial upper bound

            M = 1 # initialze iteration counter
            Mmax = 1000 # max iterations (could prbably be 100 too)

            if Fmax > 0:
                while M <= Mmax:

                    # BEGIN: bi-section algorithm
                    midpt = (lB + uB)/2
                    fODF_lB = np.sum((np.abs(fODF - lB) - fODF - lB)*AREA,0)
                    fODF_midpt = np.sum((np.abs(fODF - midpt) - fODF - midpt)*AREA,0)

                    if fODF_midpt == 0 or (uB - lB)/2 < 10**-8:

                        EPS = midpt

                        break

                    else:

                        M = M + 1

                        if np.sign(fODF_midpt) == np.sign(fODF_lB):

                            lB = midpt

                        else:

                            uB = midpt
                    # END: bi-section algorithm

                # Subract solution from each ODF point
                ODF = (1/2)*(np.abs(ODF - EPS) + ODF - EPS)
                ODF = ODF.real

                odf_pts = np.arange(ODF.size)

                # due to numerical error, we manually set
                # very very very tiny peaks to zero afte the fact...
                for p in odf_pts:
                    if ODF[p] > -10**-8 and ODF[p] < 0:
                        ODF[p] = 0
                    elif ODF[p] < 10**-8 and ODF[p] > 0:
                        ODF[p] = 0

            # Re-expand the rectified fODF into SH's
            clm_rec = np.matmul(AREA*ODF,np.conj(H))
            clm = clm_rec
            c00 = clm[0]
            clm = clm/c00
            clm = clm*(1/np.sqrt(4*np.pi))

            SH_rec[:,vox] = clm

        # zeta and FAA calculations
        # NOTE: zeta is not affected by the rectification, only FAA
        zeta[vox] = a00*np.sqrt(bval)/np.pi
        faa[vox] = np.sqrt(3*np.sum(np.abs(clm[1:6]**2))/(5*np.abs(clm[0])**2 + 2 * np.sum(np.abs(clm[1:6]**2))))

        # BEGIN: construct axonal DT (aDT)

        c00 = clm[0]
        c2_2 = clm[1]
        c2_1 = clm[2]
        c20 = clm[3]
        c21 = clm[4]
        c22 = clm[5]

        A11 = ((np.sqrt(30)/3)*c00 - (np.sqrt(6)/3)*c20 + c22 + c2_2)
        A22 = ((np.sqrt(30)/3)*c00 - (np.sqrt(6)/3)*c20 - c22 - c2_2)
        A33 = ((np.sqrt(30)/3)*c00 + (2*np.sqrt(6)/3)*c20)

        A12 = (1j*(c22 - c2_2))
        A13 = ((-c21 + c2_1))
        A23 = (1j*(-c21 - c2_1))

        aDT = np.array([A11, A12, A13, A12, A22, A23, A13, A23, A33]).real
        aDT = 1/(c00*np.sqrt(30))*aDT
        iaDT = np.reshape(aDT,(3,3)).real

        # END: construct axonal DT (aDT)

        # BEGIN: FBWM portion...
        if fbwm_flag == 1:

            f_grid = np.linspace(0,1,100) # define AWF grid (100 pts evenly spaced between 0 (min) and 1 (max))
            f_grid = f_grid * np.ones((1,100)) # makes it in to a proper array...weird but it works

            int_grid = np.linspace(0,99,100) # define grid points to iterate over (100 of them)
            int_grid = int_grid * np.ones((1,100)) # same as above, makes it a proper array object...?

            # This holds the SH basis sets for each b-value shell
            shB = [B1,B2,B] # list object: to access, shB[0] = B1 (for example)

            imgn1 = img1[:,vox] # b1000 DWI images
            imgn2 = img2[:,vox] # b2000 DWI images

            # This hold all DWI volumes for each b-vlaue shell
            IMG = [imgn1, imgn2, imgn] # list object: to access, IMG[0] = imgn1 (for example)

            # BEGIN: DT construction (should be modified to fit with PyDesigner output)
            # This is based on DKE DT output as of 05/05/2020
            iDT = np.array([DT[vox,0],DT[vox,3],DT[vox,4],DT[vox,3],DT[vox,1],DT[vox,5],DT[vox,4],DT[vox,5],DT[vox,2]])
            iDT = np.reshape(iDT,(3,3))
            # END: DT construction

            # initialze correction factor elements that will be looped over and filled accordingly...
            g2l_fa_R = np.zeros((len(harmonics),f_grid.shape[1]), order = 'F')
            g2l_fa_R_b = np.zeros((len(BT),f_grid.shape[1],len(harmonics)), order = 'F')
            g2l_fa_R_large = np.zeros((len(harmonics),f_grid.shape[1]), order = 'F')



            # BEGIN: cost function
            # Not many comments here, See McKinnon 2018 FBWm paper for details
            for b in range(0,len(BT)):

                idx_hyper = BT[b] * np.power(f_grid,2) * np.power(zeta[vox],-2) < 20 # when should hypergeometric function be implemented? When b*D is small
                idx_Y = 0

                for l in degs[::2]:

                    hypergeom_opt = np.sum((sp.special.gamma((l+1)/2 + int_grid) * sp.special.gamma(l+(3/2)) * ((-BT[b] * f_grid[idx_hyper]**2 * zeta[vox]**-2)*np.ones((1,len(f_grid[idx_hyper])))).T ** int_grid / (sp.special.factorial(int_grid) * sp.special.gamma(l+(3/2) + int_grid) * sp.special.gamma((l+1)/2))),1)*np.ones((1,len(f_grid[idx_hyper])))
                    g2l_fa_R[idx_Y:idx_Y+(2*l+1),np.squeeze(idx_hyper)] = npm.repmat((sp.special.factorial(l/2) * (BT[b] * f_grid[idx_hyper]**2 * zeta[vox]**-2) ** ((l+1)/2) / sp.special.gamma(l+(3/2)) * hypergeom_opt),(2*l+1),1) # Eq. 9 FBWM paper
                    idx_Y = idx_Y + (2*l+1)

                g2l_fa_R_b[b,np.squeeze(idx_hyper),:] = g2l_fa_R[:,np.squeeze(idx_hyper)].T

                idx_Y = 0

                for l in degs[::2]:

                    g2l_fa_R_large[idx_Y:idx_Y+(2*l+1), np.squeeze(~idx_hyper)] = npm.repmat((np.exp(-l/2 * (l+1) / ((2*BT[b] * (f_grid[~idx_hyper]**2 * zeta[vox]**-2))))),(2*l+1),1) # Eq. 20 FBI paper
                    idx_Y = idx_Y + (2*l+1)

                g2l_fa_R_b[b,np.squeeze(~idx_hyper),:] = g2l_fa_R_large[:,np.squeeze(~idx_hyper)].T

            # here is the core piece of the cost function:
            for grid in np.squeeze(int_grid.astype('int')):
                for b in range(0,len(BT)):

                    awf = f_grid[:,grid] # define AWF grid

                    # Se and Sa are the theoretical extra-axonal and intra-axonal signals that will be compared with IMG[:] DWI values for each voxel element
                    Se = (b0n * np.exp((-BT[b] * (1-awf)**-1) * np.diag((GT[b].dot((iDT - (awf**3 * zeta[vox]**-2) * iaDT).dot(GT[b].T)))))) * (1 - awf) # Eq. 3 FBWM paper
                    Sa = (2*np.pi*b0n*zeta[vox]*np.sqrt(np.pi/BT[b])) * (shB[b].dot((Pl0 * np.squeeze(g2l_fa_R_b[b,grid,:])*clm))) # Eq. 4 FBM paper

                    cost_fn[grid,vox] = cost_fn[grid,vox] + ndir[b]**-1 * np.sum((IMG[b] - Se - Sa)**2)

                cost_fn[grid,vox] = b0n**-1 * np.sqrt(len(BT)**-1 * cost_fn[grid,vox]) # Eq. 21 FBWM paper

                iDT_img[:,:,vox] = iDT
                iaDT_img[:,:,vox] = iaDT

if fbwm_flag == 1:

    min_cost_fn_idx = np.argsort(cost_fn, axis = 0) # find the indexes of the sorted cost_fn values
    min_cost_fn = np.take_along_axis(cost_fn,min_cost_fn_idx,axis=0) # sort those values

    awf_grid = np.linspace(0,1,100) # another AWF grid

    min_awf = awf_grid[min_cost_fn_idx[0,:]] # grad the minimum AWF value based on the cost_fn sorting done immeidately prior to this...

    Da = min_awf**2 / zeta**2 # Eq. 22 McKinnon (2018). intrinsic intra-axonal diffusivity

    # loop over the voxels to get extra-axonal diffusion tensor (De)...
    for vox in voxList:
        if mask[vox] == 1 and b0n > 0:

            De[:,:,vox] = (iDT_img[:,:,vox] - (min_awf[vox]**3 * zeta[vox]**-2) * iaDT_img[:,:,vox]) / (1 - min_awf[vox])

            iDe = De[:,:,vox] # intermeidate De
            iDe[np.isnan(iDe)] = 0
            iDe[np.isinf(iDe)] = 0
            L,V = np.linalg.eig(iDe) # L : eigVals and V: eigVecs
            L = np.sort(L) # sort them (this is ascending)
            L = L[::-1] # reverse the order so they are descending (high -> low)

            N = 1 # initialize counter

            while L[0] < 0 or L[1] < 0 or L[2] < 0: # find new AWF values if L's are < 0

                N = N + 1

                if N < 100:

                    min_awf[vox] = awf_grid[min_cost_fn_idx[N,vox]]

                else:

                    min_awf[vox] = 0
                    De[:,:,vox] = (iDT_img[:,:,vox] - (min_awf[vox]**3 * zeta[vox]**-2) * iaDT_img[:,:,vox]) / (1 - min_awf[vox])
                    Da[vox] = min_awf[vox]**2 / zeta[vox]**2

                    break

                # update De here...
                De[:,:,vox] = (iDT_img[:,:,vox] - (min_awf[vox]**3 * zeta[vox]**2) * iaDT_img[:,:,vox]) / (1 - min_awf[vox])
                Da[vox] = min_awf[vox]**2 / zeta[vox]**2 # recalculate Da too...

            # Now recalculate eigVals again with correct AWF values

            iDe = De[:,:,vox]
            iDe[np.isnan(iDe)] = 0
            iDe[np.isinf(iDe)] = 0
            L,V = np.linalg.eig(iDe) # L : eigVals and V: eigVecs
            L = np.sort(L) # again, ascending
            L = L[::-1] # now, descending

            De_ax[vox] = L[0] # Eq. 24 FBWM paper, axial extra-axonal diffusivity
            De_rad[vox] = (L[1] + L[2])/2 # radial De
            De_fa[vox] = np.sqrt(((L[0] - L[1]) ** 2 + (L[0] - L[2]) ** 2 + (L[1] - L[2]) ** 2 ) / (2 * np.sum(L ** 2))) # extra-axonal FA
            De_mean[vox] = (1/3) * (2 * De_rad[vox] + De_ax[vox]) # average De

    De[np.isnan(De)] = 0
    De[np.isinf(De)] = 0


# HERE IS WHERE IMAGES ARE READ OUT....
# need to automate dimension pull and reshape...

zeta = np.reshape(zeta,(74,74,42), order = 'F')
faa = np.reshape(faa,(74,74,42), order = 'F')

zetaImg = nib.Nifti1Image(zeta,hdr.affine,hdr.header)
faaImg = nib.Nifti1Image(faa,hdr.affine,hdr.header)

zetaImg.to_filename(op.join(root,'zeta.nii'))
faaImg.to_filename(op.join(root,'faa.nii'))

if fbwm_flag == 1:

    min_awf = np.reshape(min_awf,(74,74,42), order = 'F')
    Da = np.reshape(Da,(74,74,42), order = 'F')
    De_mean = np.reshape(De_mean,(74,74,42), order = 'F')
    De_ax = np.reshape(De_ax,(74,74,42), order = 'F')
    De_rad = np.reshape(De_rad,(74,74,42), order = 'F')
    De_fa = np.reshape(De_fa,(74,74,42), order = 'F')
    min_cost_fn = np.reshape(min_cost_fn[0,:],(74,74,42), order = 'F')

    awfImg = nib.Nifti1Image(min_awf,hdr.affine,hdr.header)
    DaImg = nib.Nifti1Image(Da,hdr.affine,hdr.header)
    De_meanImg = nib.Nifti1Image(De_mean,hdr.affine,hdr.header)
    De_axImg = nib.Nifti1Image(De_ax,hdr.affine,hdr.header)
    De_radImg = nib.Nifti1Image(De_rad,hdr.affine,hdr.header)
    De_faImg = nib.Nifti1Image(De_fa,hdr.affine,hdr.header)
    minCostImg = nib.Nifti1Image(min_cost_fn,hdr.affine,hdr.header)

    awfImg.to_filename(op.join(root,'awf.nii'))
    DaImg.to_filename(op.join(root,'da.nii'))
    De_meanImg.to_filename(op.join(root,'de_mean.nii'))
    De_axImg.to_filename(op.join(root,'de_ax.nii'))
    De_radImg.to_filename(op.join(root,'de_rad.nii'))
    De_faImg.to_filename(op.join(root,'fae.nii'))
    minCostImg.to_filename(op.join(root,'minCost.nii'))
