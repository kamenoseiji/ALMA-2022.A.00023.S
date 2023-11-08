import numpy as np
REFANT = 'DA64'
#Fields: 5
#  ID   Code Name                RA               Decl           Epoch   SrcId      nRows
#  0    none J0423-0120          04:23:15.800720 -01.20.33.06550 ICRS    0       16342348
#  1    none J0239-0234          02:39:45.472276 -02.34.40.91451 ICRS    1        1706628
#  2    none J0235-0402          02:35:07.341781 -04.02.05.31655 ICRS    2       15286084
#  3    none NGC_1052            02:41:04.798500 -08.15.20.75195 ICRS    3       21737232
#  4    none J0215-0343          02:15:11.506503 -03.43.07.89436 ICRS    4        6451148
BPCAL = 'J0423-0120'
PHCAL= 'J0235-0402'
TARGET= 'NGC_1052'
CHECK= 'J0215-0343'
#-------- Read spectral data
velocFile = open('X109d26e_NGC_1052.SPW3.contsub.NatClean.veloc.txt', 'r')
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
#plt.step( veloc, flux, where='mid')
#-------- modelfit-size for blue component
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
