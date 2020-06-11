#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Perform fiber ball imaging (FBI) and FBI white matter model (FBWM)
analyses
Author(s): Russell Glenn, Emilie McKinnon and Hunter Moss (MATLAB code)
           Hunter Moss (1st python 3 iteration, 05/05/2020)
Medical University of South Carolina (MUSC)
Modified further by Hunter Moss, Emilie McKinnon and Siddhartha Dhiman
(2016, 2017, 2018, 2019, 2020)
'''

import sys as sys
import os # mkdir
import os.path as op # path
import numpy as np # array, ndarray
import numpy.matlib as npm
import nibabel as nib
import scipy as sp
from scipy.io import loadmat
import math

def loaddwi(dwi, bval, bvec):
    '''
    Loads DWI file into memory, along with it's bvals and bvecs

    Parameters
    ----------
    dwi : str
        Path to NifTi file
    bval : str
        Path to bval file
    bvec : str
        Path to bvec file

    Returns
    -------
    img : (x, y, z, n) ndarray
        Numpy array object containing DWI volume
    bvals : (n,) ndarray
        Numpy array object containg bvals
    bvecs : (n, 3) ndarray
        Numpy object containing bvecs
    '''
    if not op.exists(dwi):
        raise FileNotFoundError('Input volume does not exist: {}'
        ''.format(dwi))
    if not op.exists(bval):
        raise FileNotFoundError('Input bval file does not exist: {}'
        ''.format(bval))
    if not op.exists(bval):
        raise FileNotFoundError('Input bvec file does not exist: {}'
        ''.format(bvec))
    img = np.array((nib.load(dwi).dataobj))
    bvals = np.rint(np.loadtxt(bval))
    bvecs = np.loadtxt(bvec)
    bvecs = np.reshape(bvecs, (max(bvecs.shape), 3))
    return img, bvals, bvecs

def loaddwivols(dwi, bval, bvec):
    '''
    Loads multiple 4D DWIs into single 4D numpy array

    Parameters
    ----------
    dwi : list of str
        Paths of 4D DWI files to read
    bval : list of str
        Paths of bval files
    bvec : list of str
        Path of bvec files

    Returns
    -------
    img : (x, y, z, n) ndarray
        Numpy array object containing DWI volumes
    bvals : (n,) ndarray
        Numpy array object containg bvals
    bvecs : (n, 3) ndarray
        Numpy object containing bvecs
    '''
    if not (isinstance(dwi, list) and \
        isinstance(bval, list) and \
        isinstance(bvec, list)):
        raise TypeError('Inputs must be provided as a list')
    exist_idx = [op.exists(x) for x in dwi]
    if not any(exist_idx):
        raise FileNotFoundError('Input volume(s) do not exist: '
        '\n{}'.format(dwi[exist_idx]))
    exist_idx = [op.exists(x) for x in bvec]
    if not any(exist_idx):
        raise FileNotFoundError('Input bvals(s) do not exist: '
        '\n{}'.format(bvec[exist_idx]))
    exist_idx = [op.exists(x) for x in bval]
    if not any(exist_idx):
        raise FileNotFoundError('Input bvecs(s) do not exist: '
        '\n{}'.format(bval[exist_idx]))
    img = np.array[np.concatenate(np.array((nib.load(x).dataobj))) \
        for x in dwis]
    img = np.array((nib.load(dwi[0]).dataobj))
    bvals = np.rint(np.loadtxt(bval[0]))
    bvecs = np.loadtxt(bvec[0])
    bvecs = np.reshape(bvecs, (max(bvecs.shape), 3))
    if len(dwi) > 1:
        for i in dwi[:, :, :, 1:]:
            img = np.concatenate((img, np.array(nib.load(i).dataobj)),
            axis=3)
        for i in bval[1:]:
            bvals = np.concatenate((bvals, np.loadtxt(i)), axis=0)
        for i in bvec[1:]:
            tmp = np.loadtxt(i)
            tmp = np.reshape(tmp, (max(tmp.shape), 3))
            bvecs = np.concatenate((bvecs, tmp), axis=0)
    return img, bvals, bvecs

def loadnifti(path):
    '''
    Loads single NifTi file into memory

    Parameters
    ----------
    path : str
        Path to NifTi file
    
    Returns
    -------
    out : (x, y, z) ndarray
        Numpy array object containing volume
    hdr : Nibabel header
        NifTi file header information
    '''
    tmp = nib.load(path)
    img = tmp.dataobj
    hdr = tmp.header
    return img, hdr
    
    
    

