
# Python function that reads in MaStar SSP models and interpolates
# Daniel Thomas, 14/11/2019

import numpy as np
from astropy.io import fits


def t(ver):
	hdul = fits.open('MaStar_SSP_'+ver+'.fits.gz')
	t=hdul[1].data[:,0,0,0] #parameter matrix
	hdul.close()
	return t

def Z(ver):
	hdul = fits.open('MaStar_SSP_'+ver+'.fits.gz')
	Z=hdul[1].data[0,:,0,1] #parameter matrix
	hdul.close()
	return Z

def s(ver):
	hdul = fits.open('MaStar_SSP_'+ver+'.fits.gz')
	s=hdul[1].data[0,0,:,2] #parameter matrix
	hdul.close()
	return s

def wave(ver):
	hdul = fits.open('MaStar_SSP_'+ver+'.fits.gz')
	wave=hdul[2].data[0,:] #lambda array	
	hdul.close()
	return wave

def flux(ver,lib):
	hdul = fits.open('MaStar_SSP_'+ver+'.fits.gz')
	if (lib=='Th'):
		fluxgrid=hdul[3].data #fluxgrid matrix
	if (lib=='E'):
		fluxgrid=hdul[4].data #fluxgrid matrix
	if (lib=='gold'):
		fluxgrid=hdul[3].data #fluxgrid matrix
	hdul.close()
	return fluxgrid


def inter(t,Z,s,fluxgrid,tin,Zin,sin):

	tfluxgrid=np.ndarray((len(Z),len(s),len(fluxgrid[0,0,0,:])))
	tZfluxgrid=np.ndarray((len(s),len(fluxgrid[0,0,0,:])))
	flux=np.ndarray(len(fluxgrid[0,0,0,:]))

	if (Zin<-1.35 and tin<1):
		for i in range(len(fluxgrid[0,0,0,:])):
			flux[i]=-99
	else:
		for i in range(1,len(t)):
			if (tin<=t[i]):
				t0=t[i-1]
				t1=t[i]
				f0=fluxgrid[i-1,:,:,:]
				f1=fluxgrid[i,:,:,:]
				tfluxgrid[:,:,:]=f0+(f1-f0)/(t1-t0)*(tin-t0)
				break
			t0=t[i-1]
			t1=t[i]
			f0=fluxgrid[i-1,:,:,:]
			f1=fluxgrid[i,:,:,:]
			tfluxgrid[:,:,:]=f0+(f1-f0)/(t1-t0)*(tin-t0)
		for j in range(len(Z)):
			if (Zin<=Z[j]):
				Z0=Z[j-1]
				Z1=Z[j]
				f0=tfluxgrid[j-1,:,:]
				f1=tfluxgrid[j,:,:]
				tZfluxgrid[:,:]=f0+(f1-f0)/(Z1-Z0)*(Zin-Z0)
				break
			Z0=Z[j-1]
			Z1=Z[j]
			f0=tfluxgrid[j-1,:,:]
			f1=tfluxgrid[j,:,:]
			tZfluxgrid[:,:]=f0+(f1-f0)/(Z1-Z0)*(Zin-Z0)
		for k in range(len(s)):
			if (sin<=s[k]):
				s0=s[k-1]
				s1=s[k]
				f0=tZfluxgrid[k-1,:]
				f1=tZfluxgrid[k,:]
				flux[:]=f0+(f1-f0)/(s1-s0)*(sin-s0)
				break
			s0=s[k-1]
			s1=s[k]
			f0=tZfluxgrid[k-1,:]
			f1=tZfluxgrid[k,:]
			flux[:]=f0+(f1-f0)/(s1-s0)*(sin-s0)
			
	return flux

