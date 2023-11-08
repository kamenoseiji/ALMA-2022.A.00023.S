import numpy as np
import matplotlib.pyplot as plt
prefix =  'uid___A002_X10a7a20_X65a8'
from almahelpers_localcopy import tsysspwmap
#Fields: 4
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0238+1636          02:38:38.930108 +16.36.59.27461 ICRS    0       19442150
#  1    none J0235-0402          02:35:07.341781 -04.02.05.31655 ICRS    1       17905550
#  2    none NGC_1052            02:41:04.798500 -08.15.20.75195 ICRS    2       25575800
#  3    none J0215-0343          02:15:11.506503 -03.43.07.89436 ICRS    3        7670950
BPCAL = 'J0238+1636'
PHCAL= 'J0235-0402'
TARGET= 'NGC_1052'
CHECK= 'J0215-0343'
REFANT = 'DA44'
WVRSPW = '4'
'''
os.system('rm -rf ' + prefix + '_flagonline.txt')
importasdm(prefix)
plotants(vis=prefix + '.ms', figfile=prefix + '_plotants.png')
#browsetable(tablename = prefix +'.ms')
listobs(vis=prefix+'.ms', scan='', spw='', verbose=True, listfile=prefix+'.listobs')
msfile = prefix + '.ms'
flagdata(vis=msfile,mode='manual', autocorr=True, flagbackup=False)
flagdata(vis=msfile,mode='manual', intent='*POINTING*,*ATMOSPHERE*', flagbackup=False)
flagdata(vis=msfile,mode='manual', antenna='DA54', timerange='09:46:10~10:05:41', flagbackup=False) # PRTSPR-67128 DA54 EL brake engaged
flagmanager(vis=msfile,mode='save',versionname='Original')
#---- Tsys table
os.system('rm -rf *T0.tsys')
gencal(vis = msfile, caltable = prefix + 'T0.tsys', caltype = 'tsys')
plotms(prefix + 'T0.tsys', xaxis='freq', yaxis='tsys', spw='17:4~123,19:4~123,21:4~123,23:4~123', iteraxis='antenna', coloraxis='corr', gridrows=4, gridcols=4)
#---- WVR
os.system('rm -rf ' + prefix + '.WVR*')
wvrgcal(segsource=True, caltable=prefix + '.WVR', vis=prefix+'.ms', wvrspw=[4], toffset=0, wvrflag=[], statsource=TARGET, minnumants=1, maxdistm=1500.0, mingoodfrac=0.0)
plotms(prefix + '.WVR', xaxis='time', yaxis='phase', spw='17', iteraxis='antenna')
#---- Apply Tsys and WVR
tsysmap = tsysspwmap(vis=prefix+'.ms', tsystable=prefix + 'T0.tsys', tsysChanTol=1)
for fields in list(set([BPCAL, PHCAL, TARGET])):
    applycal(vis=prefix+'.ms', field=fields, flagbackup=False, spw='17,19,21,25', interp='linear', gaintable=[prefix+'T0.tsys', prefix + '.WVR'], gainfield=fields, spwmap=[tsysmap, []], calwt=True)
applycal(vis=prefix+'.ms', field=CHECK, flagbackup=False, spw='17,19,21,25', interp='linear', gaintable=[prefix+'T0.tsys', prefix + '.WVR'], gainfield=['',''], spwmap=[tsysmap, []], calwt=True)
os.system('rm -rf ' + prefix+'_Cal.ms*'); split(vis=prefix+'.ms', outputvis=prefix+'_Cal.ms', datacolumn='corrected', spw='17,19,21,25')
#plotms(prefix + '_Cal.ms', scan='3', xaxis='time', yaxis='phase', spw='3', avgchannel='1920', coloraxis='corr', antenna='*&', iteraxis='baseline', gridrows=4, gridcols=4)
#-------- Phase Cal for bandpass
os.system('rm -rf P0*')
gaincal(vis=prefix+'_Cal.ms', caltable='P0', spw='0:4~123,1:4~123,2:4~123,3', field= BPCAL, scan='3', selectdata=True, solint='int', refant=REFANT, refantmode='strict', calmode='p')
plotms('P0', xaxis='time', yaxis='phase', spw='*', iteraxis='antenna', coloraxis='corr', gridrows=4, gridcols=4)
os.system('rm -rf B0*')
bandpass(vis = prefix + '_Cal.ms', caltable = 'B0', gaintable='P0', spw='0:4~123,1:4~123,2:4~123,3', field=BPCAL, scan='3', minblperant=5, minsnr=4, solint='inf', bandtype='B', fillgaps=1, refant = REFANT, solnorm = True)
plotms('B0', xaxis='frequency', yaxis='amp', spw='*', iteraxis='antenna', coloraxis='corr', gridrows=4, gridcols=4, plotrange=[0,0,0,1.5])
execfile('/home/skameno/CASAimaging/2022.A.00023.S/smBP.py') # 9-channel spline smoothing, generating B0.smooth
plotms('B0.smooth', xaxis='frequency', yaxis='amp', spw='*', iteraxis='antenna', coloraxis='corr', gridrows=4, gridcols=4, plotrange=[0,0,0,1.5])
#-------- Phase Cal to align SPW 
os.system('rm -rf P1')
gaincal(vis=prefix+'_Cal.ms', caltable='P1', spw='0:4~123,1:4~123,2:4~123,3', timerange='09:09:00~09:09:20', gaintable = 'B0.smooth', field=BPCAL, selectdata=True, solint='inf', refant=REFANT, gaintype='G', calmode='p', minsnr=7, spwmap=[0,1,2,3])
plotms('P1', xaxis='spw', yaxis='phase', spw='*', iteraxis='antenna', coloraxis='corr')
#-------- Phase Cal for all
os.system('rm -rf P2')
gaincal(vis=prefix+'_Cal.ms', caltable='P2', spw='0:4~123,1:4~123,2:4~123,3', gaintable = ['B0.smooth','P1'], field='0,1,2,3', selectdata=True, solint='int', refant=REFANT, gaintype='T', combine='spw', calmode='p', minsnr=3, spwmap=[[0,1,2,3],[0,1,2,3]])
plotms('P2', xaxis='time', yaxis='phase', spw='*', iteraxis='antenna', coloraxis='field')
#-------- Flux cal
IQUV = [1.304, 0.0, 0.0, 0.0]
setjy(vis=prefix + '_Cal.ms', field=BPCAL, spw='', standard='manual', fluxdensity=IQUV, spix = [-0.7,0], reffreq = '320.0GHz', usescratch=False)
os.system('rm -rf G0*')
gaincal(vis=prefix + '_Cal.ms', caltable = 'G0', spw ='0:4~123,1:4~123,2:4~123,3', field='0,1,2,3', minsnr=5.0, solint='inf', selectdata=True, solnorm=False, refant = REFANT, refantmode='strict', gaintable = ['B0.smooth', 'P1', 'P2'], calmode = 'ap', gaintype='G', spwmap=[[0,1,2,3],[0,1,2,3],[0,0,0,0]], interp=['nearest', 'nearest','linear'])
plotms('G0', xaxis = 'time', yaxis = 'amp', spw='*', iteraxis = 'antenna', coloraxis='field')
'''
fluxscale(vis=prefix + '_Cal.ms', caltable= 'G0', fluxtable='G0.flux', reference=BPCAL, transfer=[PHCAL,TARGET], refspwmap=[0,1,2,3])
#2023-08-26 23:18:36 INFO fluxscale	##########################################
#2023-08-26 23:18:36 INFO fluxscale	##### Begin Task: fluxscale          #####
#2023-08-26 23:18:36 INFO fluxscale	fluxscale( vis='uid___A002_X10a7a20_X65a8_Cal.ms', caltable='G0', fluxtable='G0.flux', reference=['J0238+1636'], transfer=['J0235-0402', 'NGC_1052'], listfile='', append=False, refspwmap=[0, 1, 2, 3], gainthreshold=-1.0, antenna='', timerange='', scan='', incremental=False, fitorder=1, display=False )
#2023-08-26 23:18:36 INFO fluxscale	****Using NEW VI2-driven calibrater tool****
#2023-08-26 23:18:36 INFO fluxscale	Opening MS: uid___A002_X10a7a20_X65a8_Cal.ms for calibration.
#2023-08-26 23:18:36 INFO fluxscale	Initializing nominal selection to the whole MS.
#2023-08-26 23:18:36 INFO fluxscale	Beginning fluxscale--(MSSelection version)-------
#2023-08-26 23:18:37 INFO fluxscale	 Found reference field(s): J0238+1636
#2023-08-26 23:18:37 INFO fluxscale	 Found transfer field(s):  J0235-0402 NGC_1052
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for J0235-0402 in SpW=0 (freq=3.05705e+11 Hz) is: 0.101811 +/- 0.00182296 (SNR = 55.8495, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for J0235-0402 in SpW=1 (freq=3.07588e+11 Hz) is: 0.101953 +/- 0.00178749 (SNR = 57.0368, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for J0235-0402 in SpW=2 (freq=3.17755e+11 Hz) is: 0.0980689 +/- 0.00184071 (SNR = 53.2779, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for J0235-0402 in SpW=3 (freq=3.19637e+11 Hz) is: 0.0977055 +/- 0.00193473 (SNR = 50.5008, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for NGC_1052 in SpW=0 (freq=3.05705e+11 Hz) is: 0.701114 +/- 0.00382401 (SNR = 183.345, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for NGC_1052 in SpW=1 (freq=3.07588e+11 Hz) is: 0.697079 +/- 0.00381092 (SNR = 182.916, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for NGC_1052 in SpW=2 (freq=3.17755e+11 Hz) is: 0.670501 +/- 0.00386092 (SNR = 173.664, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Flux density for NGC_1052 in SpW=3 (freq=3.19637e+11 Hz) is: 0.672436 +/- 0.00393823 (SNR = 170.746, N = 100)
#2023-08-26 23:18:39 INFO fluxscale	 Fitted spectrum for J0235-0402 with fitorder=1: Flux density = 0.0998689 +/- 0.000209098 (freq=312.612 GHz) spidx: a_1 (spectral index) =-1.01994 +/- 0.107927 covariance matrix for the fit:  covar(0,0)=8.5644e-05 covar(0,1)=0.000864013 covar(1,0)=0.000864013 covar(1,1)=1.20656
#2023-08-26 23:18:39 INFO fluxscale	 Fitted spectrum for NGC_1052 with fitorder=1: Flux density = 0.68513 +/- 0.00150322 (freq=312.612 GHz) spidx: a_1 (spectral index) =-1.03001 +/- 0.112596 covariance matrix for the fit:  covar(0,0)=7.94287e-06 covar(0,1)=5.83105e-05 covar(1,0)=5.83105e-05 covar(1,1)=0.110906
#2023-08-26 23:18:39 INFO fluxscale	Storing result in G0.flux
#2023-08-26 23:18:39 INFO fluxscale	Writing solutions to table: G0.flux
#2023-08-26 23:18:39 INFO fluxscale	Task fluxscale complete. Start time: 2023-08-27 08:18:36.337003 End time: 2023-08-27 08:18:39.151158
#2023-08-26 23:18:39 INFO fluxscale	##### End Task: fluxscale            #####
#2023-08-26 23:18:39 INFO fluxscale	##########################################
applycal(vis=prefix + '_Cal.ms', field=BPCAL, flagbackup=False, spw='0:4~123,1:4~123,2:4~123,3', interp=['nearest', 'nearest', 'linear', 'nearest'], gaintable=['B0.smooth', 'P1', 'P2','G0.flux'], gainfield=[BPCAL, BPCAL, BPCAL, BPCAL], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]], calwt=True)
applycal(vis=prefix + '_Cal.ms', field=PHCAL, flagbackup=False, spw='0:4~123,1:4~123,2:4~123,3', interp=['nearest', 'nearest', 'linear', 'nearest'], gaintable=['B0.smooth', 'P1', 'P2','G0.flux'], gainfield=[BPCAL, BPCAL, TARGET, PHCAL], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]], calwt=True)
applycal(vis=prefix + '_Cal.ms', field=TARGET, flagbackup=False, spw='0:4~123,1:4~123,2:4~123,3', interp=['nearest', 'nearest', 'linear', 'nearest'], gaintable=['B0.smooth', 'P1', 'P2','G0.flux'], gainfield=[BPCAL, BPCAL, TARGET, TARGET], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]], calwt=True)
applycal(vis=prefix + '_Cal.ms', field=CHECK, flagbackup=False, spw='0:4~123,1:4~123,2:4~123,3', interp=['nearest', 'nearest', 'linear', 'nearest'], gaintable=['B0.smooth', 'P1', 'P2','G0.flux'], gainfield=[BPCAL, BPCAL, TARGET, TARGET], spwmap=[[0,1,2,3], [0,1,2,3], [0,0,0,0], [0,1,2,3]], calwt=True)
for sourceName in [BPCAL, PHCAL, CHECK]:
    os.system('rm -rf ' + sourceName + '.ms')
    split(prefix + '_Cal.ms', spw='0,1,2,3', field=sourceName, outputvis=sourceName + '.ms', datacolumn='corrected', timebin='18.2s')
#
sourceName = TARGET
os.system('rm -rf ' + sourceName + '.*.ms')
split(prefix + '_Cal.ms', spw='0,1,2', field=sourceName, outputvis=sourceName + '.SPW012.ms', datacolumn='corrected', timebin='18.2s')
split(prefix + '_Cal.ms', spw='3', field=sourceName, outputvis=sourceName + '.SPW3.ms', datacolumn='corrected', timebin='18.2s')
plotms(TARGET + '.SPW3.ms', spw='*', xaxis='frequency', yaxis='amp', ydatacolumn='corrected', antenna='*&', avgtime='1e6', avgbaseline=True, avgscan=True, coloraxis='corr')
'''
#-------- UVCONTSUB for NGC 1052
os.system('rm -rf ' + TARGET + '.SPW3.G0*') 
gaincal(vis=TARGET+'.SPW3.ms', caltable=TARGET + '.SPW3.G0', spw='*', selectdata=True, solint=['int','inf'], refant=REFANT, refantmode='strict', gaintype='G', calmode='p', minsnr=5)
plotms(TARGET + '.SPW3.G0', spw='*', xaxis='time', yaxis='phase', antenna='*', iteraxis='antenna', coloraxis='corr', gridrows=4, gridcols=4, plotrange=[0,0,-180,180])
applycal(vis=TARGET+'.SPW3.ms', flagbackup=False, spw='*', interp='nearest', gaintable=TARGET + '.SPW3.G0', calwt=True)
os.system('rm -rf ' + TARGET + '.SPW3.cont*') 
uvcontsub(TARGET + '.SPW3.ms', fitspw='0:1~699;1360~1919', fitorder=1, spw='0', want_cont=True) 
plotms(TARGET + '.SPW3.ms.cont', spw='*', xaxis='frequency', yaxis='phase', antenna='*&', avgtime='1e6', iteraxis='baseline', avgscan=True, coloraxis='corr', gridrows=4, gridcols=4, plotrange=[0,0,-30,30])
#-------- Bandpass delay fine adjustment using continuum data
#scanList = [ 9 , 11, 16, 18, 22, 27, 29, 31, 38, 40, 43, 49, 52, 56, 58, 63, 67, 69, 76, 78, 81, 87, 90]
os.system('rm -rf B.cont*')
bandpass(vis = TARGET + '.SPW3.ms.cont', caltable = 'B.cont', spw='*', field=  TARGET, scan='9 , 11, 16, 18, 22, 27, 29, 31, 38, 40, 43, 49, 52, 56, 58, 63, 67, 69, 76, 78, 81, 87, 90', combine='', minblperant=5, minsnr=3, solint='inf', bandtype='B', fillgaps=1, refant = REFANT, solnorm = True)
plotms('B.cont', xaxis='frequency', yaxis='phase', spw='*', iteraxis='antenna', coloraxis='corr', gridrows=4, gridcols=4, plotrange=[0,0,-180,180])
applycal(vis=TARGET + '.SPW3.ms', field=TARGET, flagbackup=False, spw='0', interp='nearest', gaintable='B.cont', spwmap=[0], calwt=True)
plotms(TARGET + '.SPW3.ms', xaxis='frequency', yaxis='phase', spw='0', iteraxis='baseline', avgtime='1e6', avgscan=True, coloraxis='corr', gridrows=4, gridcols=4, plotrange=[0,0,-10,10])
os.system('rm -rf ' + TARGET + '*.cont*')
uvcontsub(TARGET + '.SPW3.ms', fitspw='0:1~699;1360~1919', fitorder=1, spw='0', want_cont=True)
plotms(TARGET + '.SPW3.ms.cont', spw='0', xaxis='frequency', yaxis='amp', ydatacolumn='corrected', antenna='*&', avgtime='1e6', avgbaseline=True, avgscan=True, coloraxis='corr')
plotms(TARGET + '.ms.cont', spw='*', xaxis='frequency', yaxis='phase', ydatacolumn='corrected', antenna='*&', avgtime='1e6', coloraxis='corr', iteraxis='baseline', gridrows=4, gridcols=4)
# Now, we have calibrated visibilities:
# - NGC_1052.SPW3.ms.cont    : phase-slope corrected continuum 
# - NGC_1052.SPW3.ms.contsub : phase-slope corrected line 
#-------- Maser CLEAN
tclean(vis= TARGET + '.SPW3.ms.contsub', imagename=TARGET + '.SPW3.contsub.NatClean', field=TARGET, cell=['0.001arcsec'], imsize=[512,512], spw='', specmode='cube', start=1, nchan=1919, width=1, outframe='LSRK', veltype='radio', restfreq='321.2256760GHz', deconvolver='hogbom', weighting='natural', interactive=True, niter=100000, mosweight=True, pbcor=False, savemodel='modelcolumn')
#-------- UVFITS for difmap imaging in continuum-only SPWS
os.system('rm -rf ' + TARGET + '.SPW012.CHAV.ms*')
mstransform(TARGET + '.SPW012.ms', outputvis=TARGET + '.SPW012.CHAV.ms', spw='*:4~123', datacolumn='data', chanaverage = True, chanbin=120)
os.system('rm -rf ' + TARGET + '.SPW*.fits*')
exportuvfits(TARGET + '.SPW012.CHAV.ms', fitsfile=TARGET + '.SPW012.fits', spw='0,1,2', multisource=False, field=TARGET, combinespw = True, writestation=False, padwithflags=True) 
#-------- UVFITS for maser blue component
velocFile = open('X10a7a20_NGC_1052.SPW3.contsub.NatClean.veloc.txt','r') # Read spectral data
velocLines = velocFile.readlines()
velocFile.close()
veloc, flux = [], []
for lineElement in velocLines:
    if lineElement[0] == '#': continue
    if lineElement[0] == '\n': continue
    element = lineElement.rstrip().split(' ')
    veloc = veloc + [float(element[0])]
    flux  = flux + [float(element[1])]
#
veloc, flux = np.array(veloc), np.array(flux)
plt.step( veloc, flux, where='mid')
os.system('rm -rf ' + TARGET + '.maser.*.ms*')
chList = [index for index, v in enumerate(veloc) if v > 1360 and v < 1400]
spwch = '0:%d~%d' % (chList[0], chList[-1])
mstransform(TARGET + '.SPW3.ms.contsub', outputvis=TARGET + '.maser.blue.ms', spw=spwch, datacolumn='corrected', chanaverage = True, chanbin=abs(chList[-1] - chList[0])+1)
exportuvfits(TARGET + '.maser.blue.ms', fitsfile=TARGET + '.maser.blue.fits', spw='0', multisource=False, field=TARGET, combinespw = True, writestation=False, padwithflags=False) 
#-------- UVFITS for maser red component
chList = [index for index, v in enumerate(veloc) if v > 1500 and v < 1625]
spwch = '0:%d~%d' % (chList[0], chList[-1])
mstransform(TARGET + '.SPW3.ms.contsub', outputvis=TARGET + '.maser.red.ms', spw=spwch, datacolumn='corrected', chanaverage = True, chanbin=abs(chList[-1] - chList[0])+1)
exportuvfits(TARGET + '.maser.red.ms', fitsfile=TARGET + '.maser.red.fits', spw='0', multisource=False, field=TARGET, combinespw = True, writestation=False, padwithflags=False) 
'''
