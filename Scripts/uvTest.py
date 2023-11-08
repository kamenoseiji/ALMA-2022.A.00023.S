prefix = 'uid___A002_X109d26e_X12c34'
BPCAL = 'J0423-0120'
IQUV = [2.104, 0.0, 0.0, 0.0]
#-------- Split channel-averaged data
'''
split(prefix + '.ms', spw='18', field='0', scan='3', outputvis='J0423-0120.chav.ms', datacolumn='corrected')
'''
msfile='J0423-0120.chav.ms'
REFANT = 'DA64'
'''
plotms(msfile, xaxis='time', yaxis='amp', antenna='*&', spw='0', iteraxis='baseline', coloraxis='corr', gridrows=4, gridcols=4)
setjy(vis=msfile, field=BPCAL, spw='', standard='manual', fluxdensity=IQUV, spix = [-0.7,0], reffreq = '343.4GHz', usescratch=False)
os.system('rm -rf TG0*')
gaincal(vis=msfile, caltable = 'TG0', spw ='0', field='0', minsnr=5.0, solint='int', selectdata=True, solnorm=False, refant = REFANT, refantmode='strict', calmode = 'ap', gaintype='G')
plotms('TG0', xaxis='time', yaxis='phase', spw='*', iteraxis='antenna')
applycal(vis=msfile, flagbackup=False, spw='0', interp='nearest', gaintable='TG0', calwt=True)
plotms(msfile, xaxis='time', yaxis='amp', ydatacolumn='corrected', antenna='*&', spw='0', iteraxis='baseline', coloraxis='corr', gridrows=4, gridcols=4)
os.system('rm -rf ' + BPCAL + '.testclean*')
tclean(vis=msfile,  imagename=BPCAL + '.testclean', field=BPCAL, cell=['0.001arcsec'], imsize=[512,512], spw='0', deconvolver='hogbom', weighting='natural', interactive=True, niter=1000, mosweight=True, pbcor=False)
                                                         
phaseshift(msfile, outputvis=BPCAL + '.shift.ms', datacolumn='corrected', phasecenter='J2000 04:23:15.797 -01.20.33.00')
plotms(BPCAL + '.shift.ms', xaxis='time', yaxis='phase', ydatacolumn='corrected', antenna='*&', spw='0', iteraxis='baseline', coloraxis='corr', gridrows=4, gridcols=4)
os.system('rm -rf ' + BPCAL + '.shiftclean*')
tclean(vis=BPCAL + '.shift.ms',  imagename=BPCAL + '.shiftclean', field=BPCAL, cell=['0.001arcsec'], imsize=[512,512], spw='0', deconvolver='hogbom', weighting='natural', interactive=True, niter=1000, mosweight=True, pbcor=False)
'''
os.system('rm -rf cl.test*')
uvmodelfit(vis= BPCAL + '.shift.ms', spw='0', sourcepar=[2.20, 0.071, -0.060], varypar=[True, True, True], outfile='cl.test')


