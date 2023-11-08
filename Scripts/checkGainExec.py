SCR_DIR = '/home/skameno/ALMA_Py3/'
wd = './'
DELAYCAL = False
BPPLOT   = True
antFlag = []
PLOTFMT  = 'pdf'
polName = ['XX', 'YY']
#---- DelayCal Scan
#prefix = 'uid___A002_X109d26e_X12c34'; scanList = [3, 10, 12, 17, 19, 23, 28, 30, 32, 39, 41, 44, 50, 53, 57, 59, 66, 68, 70, 77, 79, 84, 86, 93]; spwList = [17]
#prefix = 'uid___A002_X10a7a20_X65a8'; scanList = [3 , 9 , 11, 16, 18, 22, 27, 29, 31, 38, 40, 43, 49, 52, 56, 58, 63, 67, 69, 76, 78, 81, 87, 90]; spwList = [17]
prefix = 'uid___A002_X10a7a20_X65a8'; scanList = [56, 57, 58, 59, 60, 62, 63, 64, 66, 67, 68, 69]; spwList = [17]
for spw in spwList:
    exec(open(SCR_DIR + 'checkGain.py').read())
#


