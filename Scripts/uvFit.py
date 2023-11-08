SCR_DIR = '/home/skameno/ALMA_Py3/'
exec(open(SCR_DIR + 'interferometry.py').read())
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#-------- Gaussian visibility model 
def GaussVis(uvwave, theta):    # 1-dimensional Gaussian visibility (complex)
    # uvwave : uv distance in the unit of wavelength
    # theta  : FWHM in arcsec
    Factor = np.pi* (2.0 + 0.0j) / 206264.80624709636
    return np.exp( -(Factor* theta* uvwave)**2 )
#
def GaussResid(vis, flux, uvwave, theta): # Calculate residual = observed visisbility - Gaussian model
    residVis = vis - flux* GaussVis(uvwave, theta)
    return  (residVis.dot(residVis.conjugate())).real
#
def quadMin(x, y):
    a = y[0]/(x[0] - x[1])/(x[0] - x[2]) + y[1]/(x[1] - x[2])/(x[1] - x[0]) + y[2]/(x[2] - x[0])/(x[2] - x[1])
    b = y[0]* (x[1] + x[2])/(x[0] - x[1])/(x[0] - x[2]) + y[1]* (x[2] + x[0])/(x[1] - x[2])/(x[1] - x[0]) + y[2]* (x[0] + x[1])/(x[2] - x[0])/(x[2] - x[1])
    return 0.5*b/a
#
def visFit(spatialFreq, Vis):
    flux, size, step = abs(np.mean(Vis)), 0.001, 0.2
    residImag = 2.0* abs(Vis.imag.dot(Vis.imag))
    dof       = len(Vis) - 2
    while( step > 0.0001):
        print('flux=%.3f size=%e resid=%.1f/%d ' % (flux, size, dof*GaussResid(Vis, flux, spatialFreq, size)/residImag,dof))
        fluxTemp = np.array([(1.0+step)*flux, flux, flux/(1.0+step)])
        sizeTemp = np.array([(1.0+step)*size, size, size/(1.0+step)])
        resid = [GaussResid( Vis, flux, spatialFreq, sizeTemp[0]),
                 GaussResid( Vis, flux, spatialFreq, sizeTemp[1]),
                 GaussResid( Vis, flux, spatialFreq, sizeTemp[2])]
        if np.std(resid) < 0.0001: break
        size = quadMin(sizeTemp,resid)
        resid = [GaussResid( Vis, fluxTemp[0], spatialFreq, size),
                 GaussResid( Vis, fluxTemp[1], spatialFreq, size),
                 GaussResid( Vis, fluxTemp[2], spatialFreq, size)]
        if np.std(resid) < 0.0001: break
        flux = quadMin(fluxTemp,resid)
        step = 0.75* step
    #
    bestSize, bestFlux, bestResid = size, flux, dof*GaussResid(Vis, flux, spatialFreq, size)/residImag
    step = bestSize * 0.001
    minSize, maxSize = 0.0, bestSize
    minResid, maxResid = dof*GaussResid(Vis, bestFlux, spatialFreq, minSize)/residImag, bestResid
    #---- minimum acceptable size
    while minResid > bestResid + 2.57: # 99% confidence level
        minSize = minSize + step
        minResid = dof*GaussResid(Vis, bestFlux, spatialFreq, minSize)/residImag
    #---- maximum acceptable size
    while maxResid < bestResid + 2.57: # 99% confidence level
        maxSize = maxSize + step
        maxResid = dof*GaussResid(Vis, bestFlux, spatialFreq, maxSize)/residImag
    #
    return bestFlux, bestSize, minSize, maxSize
#
#-------- Maser and Cont data
maserBlueMS = 'NGC_1052.maser.blue.ms'
maserRedMS  = 'NGC_1052.maser.red.ms'
contMS  = 'NGC_1052.SPW012.CHAV.ms'
scanList = [10, 12, 17, 19, 23, 28, 30, 32, 39, 41, 44, 50, 53, 57, 59, 66, 68, 70, 77, 79, 84, 86, 93] # for X109d26e
#scanList = [9,11,16,18,22,27,29,31,38,40,43,49,52,56,58,63,67,69,76,78,81,87,90]   # for X10a7a20
chNum, chWid, freq = GetChNum(maserBlueMS, 0)
maserBlueVisList, maserRedVisList, contVisList, uvList = [], [], [], []
for scan in scanList:
    timeStamp, uvw = GetUVW(maserBlueMS, 0, scan)
    timeStamp, maserPS, maserBlueVis=  GetVisAllBL(maserBlueMS, 0, scan)
    timeStamp, maserPS, maserRedVis =  GetVisAllBL(maserRedMS,  0, scan)
    timeStamp, contPS,  contVis  =  GetVisAllBL(contMS, 0, scan)
    timeNum = len(timeStamp)
    for time_index in list(range(timeNum)):
        uvList = uvList + [np.sqrt(uvw[0,:,time_index]**2 + uvw[1,:,time_index]**2)]
        maserRedVisList  = maserRedVisList  + [np.mean(maserRedVis, axis=(0,1))[:,time_index]]
        maserBlueVisList = maserBlueVisList + [np.mean(maserBlueVis, axis=(0,1))[:,time_index]]
        contVisList  = contVisList  + [np.mean(contVis,  axis=(0,1))[:,time_index]]
    #
#
timeNum = len(uvList)
blNum = uvw.shape[1]
#-------- Visibility vs baseline length
uvWaveArray = np.array(uvList)* freq[0] / 299792458   # uv distance in unit of wavelength, uvWave[time, bl]
maserRedVisArray = np.array(maserRedVisList)   # maserRedVisArray[time, bl]
maserBlueVisArray = np.array(maserBlueVisList) # maserBlueVisArray[time, bl]
contVisArray  = np.array(contVisList)    # maserVis[time, bl]
uvWave, maserRedVis, maserBlueVis, contVis = np.zeros(blNum), np.zeros(blNum, dtype=complex), np.zeros(blNum, dtype=complex), np.zeros(blNum, dtype=complex)
useBlList = []
for bl_index in list(range(blNum)):
    flagIndex = np.where(abs(contVisArray[:,bl_index]) > 0.1* np.median(abs(contVisArray)))[0]
    uvWave[bl_index]   = np.mean( uvWaveArray[flagIndex, bl_index], axis=0 )
    #print('%d / %d' % (len(flagIndex), timeNum))
    if len(flagIndex) > 2:
        maserRedVis[bl_index]  = np.mean( maserRedVisArray[flagIndex, bl_index], axis=0 )
        maserBlueVis[bl_index] = np.mean( maserBlueVisArray[flagIndex, bl_index], axis=0 )
        contVis[bl_index]  = np.mean( contVisArray[flagIndex, bl_index], axis=0 )
        useBlList = useBlList + [bl_index]
    #
#
uvWave, contVis, maserBlueVis, maserRedVis = uvWave[useBlList], contVis[useBlList], maserBlueVis[useBlList], maserRedVis[useBlList]
MEGA = 1.0e6
plt.plot(uvWave/MEGA, abs(contVis),  'k.', label='continuum')
plt.plot(uvWave/MEGA, abs(maserBlueVis), 'b*', markersize=3, label='H$_2$O maser (blue)')
plt.plot(uvWave/MEGA, abs(maserRedVis),  'r+', markersize=3, label='H$_2$O maser (red)')
plt.yscale('log')
# plt.title('uid://A002/X10a7a20/X65a8')
plt.title('uid://A002/X109d26e/X12c34')
plt.xlabel('Spatial frequency [Mλ]')
plt.ylabel('Visibility amplitude [Jy]')
#-------- Continuum Size
contFlux, contSize, minContSize, maxContSize = visFit( uvWave, contVis )
maserBlueFlux, maserBlueSize, minMaserBlueSize, maxMaserBlueSize = visFit( uvWave, maserBlueVis )
maserRedFlux, maserRedSize, minMaserRedSize, maxMaserRedSize = visFit( uvWave, maserRedVis )
modelUVwave = np.arange(0.0, max(uvWave), 100)
plt.fill_between(modelUVwave/MEGA, contFlux* abs(GaussVis(modelUVwave, minContSize)), contFlux* abs(GaussVis(modelUVwave, maxContSize)), facecolor='black', alpha=0.25)
plt.fill_between(modelUVwave/MEGA, maserBlueFlux* abs(GaussVis(modelUVwave, minMaserBlueSize)), maserBlueFlux* abs(GaussVis(modelUVwave, maxMaserBlueSize)), facecolor='navy', alpha=0.5)
plt.fill_between(modelUVwave/MEGA, maserRedFlux* abs(GaussVis(modelUVwave, minMaserRedSize)), maserRedFlux* abs(GaussVis(modelUVwave, maxMaserRedSize)), facecolor='Red', alpha=0.5)
print('Continuum: flux=%.3f bestSize=%e maxSize=%e' % (contFlux, contSize, maxContSize))
print('H$_2$O (Blue): flux=%.3f bestSize=%e maxSize=%e' % (maserBlueFlux, maserBlueSize, maxMaserBlueSize))
print('H$_2$O (Red) : flux=%.3f bestSize=%e maxSize=%e' % (maserRedFlux, maserRedSize, maxMaserRedSize))
plt.legend()
plt.savefig('UVplot.pdf')
plt.close('all')
