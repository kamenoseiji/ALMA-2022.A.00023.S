SCR_DIR = '/home/skameno/ALMA_Py3/'  # Path to the directory where you deployed git clone
wd = './'
BPFILE = 'B0'
smoothCH = [0,0,0,9]    # smoothing channels for each SPW. 0: do not smooth
exec(open(SCR_DIR + 'interferometry.py').read())
#-------- Save original and copy for new
os.system('rm -rf ' + BPFILE + '.smooth')
os.system('cp -r ' + BPFILE + ' ' + BPFILE + '.smooth')
#-------- Open and read Bandpass Table
tb.open(BPFILE + '.smooth', nomodify=False)
colnameList = tb.colnames()
#-------- List of antennas, spws
refantIndex = unique(tb.getcol('ANTENNA2').tolist())[0]
antList = unique(tb.getcol('ANTENNA1').tolist()).tolist()
spwList = unique(tb.getcol('SPECTRAL_WINDOW_ID').tolist()).tolist()
#-------- For each antenna and SPW
for spw in spwList:
    if smoothCH[spw] < 1: continue
    SPWQ = tb.query('SPECTRAL_WINDOW_ID == %d'%(spw))
    BP = SPWQ.getcol('CPARAM')
    polNum, chNum = BP.shape[0], BP.shape[1]
    for ant in antList:
        for pol_index in list(range(polNum)):
            smBP = splineComplex(arange(chNum), BP[pol_index][:,ant], smoothCH[spw])
            BP[pol_index][:,ant] = smBP
        #
    #
    SPWQ.putcol('CPARAM', BP)
    SPWQ.close()
#
tb.close()
