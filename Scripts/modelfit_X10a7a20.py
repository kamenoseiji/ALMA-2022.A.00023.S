import numpy as np
import matplotlib.pyplot as plt
BPCAL = 'J0238+1636'
PHCAL= 'J0235-0402'
TARGET= 'NGC_1052'
CHECK= 'J0215-0343'
REFANT = 'DA44'
#-------- Read spectral data
velocFile = open('X10a7a20_NGC_1052.SPW3.contsub.NatClean.veloc.txt','r')
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
#-------- modelfit-size for blue component
'''
chList = [index for index, v in enumerate(veloc) if v > 1375 and v < 1380]
spwch = '0:%d~%d' % (chList[0], chList[-1])
uvmodelfit(vis= TARGET + '.SPW3.ms.contsub', field=TARGET, spw=spwch, comptype='G', sourcepar=[np.mean(flux[chList]), 0.0, 0.0, 1.0e-3, 0.8, 25.6], varypar=[False, False, False, True, True, True], outfile='cl.blue')
#-------- modelfit-size for red component
#chList = [index for index, v in enumerate(veloc) if v > 1500 and v < 1625]
'''
#
chWidth = 5
chExt   = int(chWidth/2)
chList = [index for index, v in enumerate(veloc) if v > 1250 and v < 1750]
for ch in chList:
    if ch % chWidth != 0: continue
    spwch = '0:%d~%d' % (ch-chExt, ch+chExt)
    meanFlux = np.mean(flux[(ch-chExt):(ch+chExt)])
    clfile= 'cl.%d' % (ch)
    uvmodelfit(vis= TARGET + '.SPW3.ms.contsub', field=TARGET, spw=spwch, sourcepar=[flux[ch], 0.0, 0.0], varypar=[False, True, True], outfile=clfile)
#
