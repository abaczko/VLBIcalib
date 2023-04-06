#!/usr/bin/env ParselTongue

###########################################################################
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###########################################################################
# Original author: Anne-Kathrin Baczko (MPIfR) baczko(AT)mpifr-bonn.mpg.de
###########################################################################

'''
This is the main python file for semi-automatic calibration.

Everything that is done is written into two log files: 
PT_inputs_xxx: All inputs to the tasks are saved and some additional information
AIPS_messages: everything printed to the aips message server

There is the option to run with AIPSLite() if no local AIPS installation is 
present.

The main script that has to be changed for each observation is 
 observation_parameters.py
it should reside in the same direcetory as this script is saved.

In addition the following private modules must be called, either by adding their
location to the local Pythonpath or by putting them also in the same directory
as this file:
 functions
 aips_tasks

The scripts will automatically create a directory 'aips_out', if not existing
where all output files are saved.
'''
import sys,os,json,logging,logging.config
import observation_parameters as OP
import numpy as np
from astropy.table import Table
#
# Load personal functions
#
import VLBIcalib.modules.functions as AF
import VLBIcalib.modules.hrk as hrk
import VLBIcalib.modules.aips_tasks as AT
#
##################################################################
# Setup logging ##
##################
logfile = open(OP.tasklog,'a')
#
log_r = open(OP.tasklog,'r')
if len(log_r.read())==0:
	log_r.close()
	logfile.write('########################################################\n')
	logfile.write('# This is a ParselTongue Log file to keep track of the # \n')
	logfile.write('# inputs to AIPS Tasks in ParselTongue. 								#\n')
	logfile.write('########################################################\n')
	logfile.write('#\n'*2)
else:
	log_r.close()
#
logfile.close()
with open('logging_config.json','r') as f:
  conf = json.load(f)

conf["handlers"]["info_file_handler"]["filename"] = OP.tasklog
logging.config.dictConfig(conf)
logger = logging.getLogger('__name__')

#######################################
aips_out = OP.aips_out+'/'
###################################################################
# Initialise AIPS disk, make the AIPS directories, imports needed ##
####################################################################
#
# Import python AIPS libraries 
#
if OP.AIPSLite:
	logger.info('Will use AIPSLite')
	import AIPSLite
	AIPSLite.setup()
else:
	logger.info('Will use local AIPS installation')
#
from AIPS import AIPS#, AIPSDisk # I am not sure if I need AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage,AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
from AIPSTV import AIPSTV

###############################################
# Starting up AIPS ############################
###############################################
#
aips_out=OP.aips_out+'/'
sys.stdout.write('Will write AIPS messages to file \n{0} \n and task inputs to file \n{1}\n'.format(OP.aipslog,OP.tasklog))
AIPS.log  = open(OP.aipslog,'a')
#
AIPS.userno	= OP.aipsid
sys.stdout.write('Printing AIPSCat.\n')
print (AIPSCat(OP.outdisk)) #list the catalogs on the desired disk
#input('Press Enter to continue if you are happy with what you see ...')
#
#logfile.close()
###############################################
# Calibration start ###########################
###############################################
if OP.load_data:
	if len(OP.uvrawdatafile)>1:
		for files in OP.uvrawdatafile:
			uvdata_catalog = hrk.loadata(files,OP.outname,OP.outclass,OP.outdisk,OP.outseq,clint=OP.cl_table_interval,digicor=OP.digicor,doconcat=1)
	else:
		uvdata_catalog = hrk.loadata(OP.uvrawdatafile[0],OP.outname,OP.outclass,OP.outdisk,OP.outseq,clint=OP.cl_table_interval,digicor=OP.digicor)
	uvdata = AIPSUVData(*uvdata_catalog)
	AF.data_exists(uvdata)
elif not OP.load_data and OP.dataPreProcessing:
	logger.info('Use existing data\n')
	uvdata_catalog = [OP.outname,OP.outclass,OP.outdisk,OP.outseq]
	print (uvdata_catalog)
	uvdata = AIPSUVData(*uvdata_catalog)
	AF.data_exists(uvdata)
	uvdata.clrstat()
#
if OP.dataPreProcessing:
	if OP.swpol:
		swpoldata = AT.swpol(uvdata,antennas=OP.swpol_antennas)
		if swpoldata.exists() == True:
			sys.stdout.write("SWPOL data exists.  Will set local variable 'uvdata' to point towards swpol data.\n")
			logger.info('SWPOL data exists.  Will set local variable uvdata to point towards swpol data.\n')
			uvdata = swpoldata
		else:
			sys.stdout.write("Something went wrong with Task('SWPOL')")
			logger.error("Something went wrong with Task('SWPOL')")
			sys.exit()
	msortdata = AT.msort(uvdata)
	msortdata = [OP.outname,'MSORT',OP.outdisk,OP.outseq]
#	if msortdata[3]<uvdata.seq:
#		logger.error('Something went wrong with Task MSORT \n')
#		sys.exit()
	msortuv = AIPSUVData(*msortdata)
	if msortuv.exists()==True:
		logger.info('Msort data exists. Assume it is already TB sorted (INDXR has been run).  Will set local variable uvdata to point towards msorted data.\n')
		uvdata = msortuv
	else:
		logger.error('Something went wrong with Task MSORT \n')
		sys.exit()
	AT.indxr(uvdata)
else:
	logger.info('No data Pre Processing (SWPOL,TABED)')
	uvdata_catalog = [OP.outname,OP.outclass,OP.outdisk,OP.outseq]
	uvdata = AIPSUVData(*uvdata_catalog)
#	AT.indxr(msortuv)

if uvdata.exists()==True:
	logger.info('Msort data exists. Assume it is already TB sorted (INDXR has been run).  Will set local variable uvdata to point towards msorted data.\n')
	#uvdata = msortuv
	uvdata.clrstat()
else:
  logger.error('Something went wrong with loading the catalog. Please check the input parameters for outname, outclass, outdisk, outseq to ensure that these point to the correct catalog. \n')
  sys.exit()

logger.info('Data Loaded: (%s, %s, %d, %d)\n',uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq)
if OP.tabed:
	if type(OP.tabed)==int:
		OP.tabed_ine			= [OP.tabed_ine]			
		OP.tabed_optype		= [OP.tabed_optype]		
		OP.tabed_keyvalue = [OP.tabed_keyvalue] 
		OP.tabed_aparm		= [OP.tabed_aparm]
		OP.tabed_inv			= [OP.tabed_inv]	
		OP.tabed_outv			= [OP.tabed_outv]	
	for i,tabed in enumerate(OP.tabed_ine): 
		AT.tabed(uvdata,ine=OP.tabed_ine[i],optype=OP.tabed_optype[i],keyvalue=OP.tabed_keyvalue[i],aparm=OP.tabed_aparm[i],inv=OP.tabed_inv[i],outv=OP.tabed_outv[i])

if OP.printBasicInfo:
	'''
	Runs: listr,imhead,possm,prtan,snplt
	'''
	AF.print_basic_infos(uvdata)

#
if OP.delete_cal_tables:
	uvdata.zap_table('SN',-1)
	clh=uvdata.table_highver('CL')
	for clv in np.arange(2,clh+1):
		uvdata.zap_table('CL',clv)
	uvdata.zap_table('FG',-1)

if OP.load_antab:
	uvdata.zap_table('TY',-1)
	uvdata.zap_table('GC',-1)
	#if len(OP.antabfile)==1:
	#	OP.antabfile=[OP.antabfile]
	for antabf in OP.antabfile:
		AT.antab(uvdata,infile=antabf)
	AT.snplt(uvdata,ine='TY',inv=1,optype='TSYS',stokes='',sources=[],opcode='ALIF')
#
if OP.get_vlba_cal:
	AT.vlog(uvdata,OP.vlbacal,OP.vlog_out)
# 

if OP.usuba:
	AT.usuba(uvdata,opcode='AUTO')
	AT.indxr(uvdata) #To be on the safe side.
	AT.snplt(uvdata,ine='CL',inv=1,optype='PHAS',stokes='',sources=[],opcode='ALST',suba=-1,plotname='after_usuba+in')
	AT.listr(uvdata,outfile='MB005_w_2_list_after_usuba.txt')

if OP.flag_outer_ch:
	uvdata.zap_table('FG',-1)
	AT.uvflg(uvdata,bif=1,eif=0,bchan=1,echan=3,outfgv=1)
	AT.uvflg(uvdata,bif=1,eif=0,bchan=AF.max_ch(uvdata)-3,echan=AF.max_ch(uvdata),outfgv=1)
#
if OP.load_flag_info:
	try:
		for flagf in OP.flagfile:
			AT.uvflg_flagfile(uvdata,intext=flagf,outfgv =1)
	except:
		print ('Something went wrong. Flag table not loaded. Does it exist?')
#
if OP.uvflg_additional:
	if type(OP.fg_bif)==int:
		OP.fg_bif=[OP.fg_bif]
		OP.fg_eif=[OP.fg_eif]
		OP.fg_timer=[OP.fg_timer]
		OP.fg_antennas=[OP.fg_antennas]
		OP.fg_stokes=[OP.fg_stokes]
		OP.fg_reason=[OP.fg_reason]
		OP.fg_bchn=[OP.fg_bchan]
		OP.fg_echan=[OP.fg_echan]
	for i in range(len(OP.fg_reason)):
		AT.uvflg(uvdata,antennas=[OP.fg_antennas[i]],timer=OP.fg_timer[i],bif=OP.fg_bif[i],eif=OP.fg_eif[i],bchan=OP.fg_bchan[i],echan=OP.fg_echan[i],outfgv=1,stokes=OP.fg_stokes[i],reason=OP.fg_reason[i])

#

if any([OP.load_flag_info,OP.flag_outer_ch,OP.uvflg_additional]): 
	logger.info('Printing FG table to file\n')
	for i in range(uvdata.table_highver('FG')):
		AT.prtab(uvdata,'FG',i+1)
else:
	logger.info('Seems no FG table has been written now or something else went wrong.\n')

if OP.runaccor:
	#uvdata.zap_table('SN',-1)
	AT.accor(uvdata)
	sn_hv = uvdata.table_highver('SN')
	
	AT.snsmo(uvdata,samptype='MWF',smotype='AMPL',bparm=[5.0,0],dobtween=1,doblank=-1,inv=sn_hv,outv=sn_hv+1,refant=AF.get_antenna_number(uvdata,OP.refant_global))
	AT.extd(uvdata,'SN',sn_hv)
	AT.tacop(uvdata,ine='SN',inv=sn_hv+1,outv=sn_hv)
	AT.extd(uvdata,'SN',sn_hv+1)
	AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='AMP',stokes='',sources=[],opcode='ALST',plotname='accor_smoothed')

if OP.correct_pang:
	if OP.usuba:
		maxan = uvdata.table_highver('AN')
		gainv = uvdata.table_highver('CL')
		gainu = gainv+1
		for i in range(maxan):
			AT.pang(uvdata,gainv=gainv,gainu=gainu,suba=i+1)
			gainv = gainu
	else:
		AT.pang(uvdata)
	AT.snplt(uvdata,ine='CL',inv=uvdata.table_highver('CL')-1,optype='PHAS',stokes='',sources=[],opcode='ALIF',plotname='before_pang')
	AT.snplt(uvdata,ine='CL',inv=uvdata.table_highver('CL'),optype='PHAS',stokes='',sources=[],opcode='ALIF',plotname='after_pang')


#
if OP.correct_eop:
	AT.eops(uvdata)
#
if OP.runaccor:
	AT.clcal(uvdata,gainv=0,gainu=0,snv=uvdata.table_highver('SN'),refant=AF.get_antenna_number(uvdata,OP.refant_global),doblank=1,dobtween=1,interpol='2PT')
	AT.snplt(uvdata,ine='CL',inv=uvdata.table_highver('CL'),optype='AMP',stokes='',sources=[],opcode='ALIF',plotname='ACCOR')
	cl_first_cal_steps = uvdata.table_highver('CL')
	sn_first_cal_steps = uvdata.table_highver('SN')
	AT.possm(uvdata,sources=[],timer=[0,0,0,0],solint=-1,gainu= cl_first_cal_steps-1,aparm=[0,1,0,0,-200,200,0,1,1,0],plotname='AC')
	AT.possm(uvdata,sources=[],timer=[0,0,0,0],solint=-1,gainu= cl_first_cal_steps,aparm=[0,1,0,0,-200,200,0,1,1,0],plotname='AC')
#

if_exit=False
if OP.cl_cal_before_mp and OP.sn_cal_before_mp:
	cl_first_cal_steps = int(OP.cl_cal_before_mp)
	sn_first_cal_steps = int(OP.sn_cal_before_mp)

elif not OP.cl_cal_before_mp and not OP.sn_cal_before_mp and uvdata.table_highver('SN')>=1:
	logging.info('Searching for first appearance of CLCAL in history file to get CL and SN resulting from ACCOR, last task before Manual phasecal.')
	history=AF.print_history(uvdata)
	for r in history:
		if r.startswith('CLCAL GAINV'):
			row_cl = r
			clsplit=row_cl.split()
		if r.startswith('CLCAL SNV'):
			row_sn = r
			snsplit=row_sn.split()
			break

if 'clsplit' in locals() and 'snsplit' in locals():
	cl_first_cal_steps= int([clsplit[i+2] for i in range(len(clsplit)) if clsplit[i].startswith('GAINU')][0])
	sn_first_cal_steps= int([snsplit[i+2] for i in range(len(snsplit)) if snsplit[i].startswith('SNV')][0])
#
	logger.info('#'*20+'\n')
	logger.info('''Calibration Tables are loaded: ANTAB, PC, FG (if specified)
	Basic correction have been applied: SWPOL, TABED, FLAG (outer channels) (if specified)
	Initial calibration steps ACCOR, PANG, EOP (if specified) are finalized.
	Please check output tables and files in 'aips_out' folder.
	Current highest CL=%d and SN=%d.
	Now instrumental phase calibration will start, either with Phascals or manual.
	#####
	The AIPS Catalog file [%s, %s, %d, %d] now has the following tables attached:\n
	''',cl_first_cal_steps,sn_first_cal_steps,uvdata.name,uvdata.klass,uvdata.disk,uvdata.seq)
else:
	if_exit=True

if OP.cl_cal_before_mp:
	if_exit=False

for row in uvdata.tables:
	logger.info(row)
logger.info('#'*20)
logger.info('\n')
#
#if runpccor:
#
if OP.get_some_listings:
	AF.get_several_listgings(uvdata)
	MATX=AF.get_matx(uvdata)
#
if OP.clear_manual_phasecal and not if_exit:
	logger.info('Will delete all calibration steps starting with manual phasecal\n')
	logger.info('To start instrumental Phase Correction over again\nDelete CL>%d and SN>%d\n',cl_first_cal_steps,sn_first_cal_steps)
	input('Press Enter if that is ok')
	logger.info('#'*20+'\n')
	for n in range(cl_first_cal_steps,uvdata.table_highver('AIPS CL')):
		uvdata.zap_table('CL',n+1)
	for n in range(sn_first_cal_steps,uvdata.table_highver('AIPS SN')):
		uvdata.zap_table('AIPS SN',n+1)
elif OP.clear_manual_phasecal and if_exit:
	logger.info('As no further calibration steps had been applied before, there is no need to delete tables.')
	input('If that is ok press enter')
#else:
#	logger.info('No deletion of tables before manual phasecal')
#
if OP.runpccor:
	AF.do_pccor(uvdata,antennas=OP.pcal_antennas,calibrator=OP.pcal_calibrator,timer=OP.pcal_timer,refant=OP.refant_global)
	logger.info('Pccor was run. The following antennas did not had PCal information:%s'%', '.join(antennas_missing))
#
if OP.manual_phase_cal:
	AF.do_manual_phasecal(uvdata,OP.mpc_calibrator,OP.mpc_antennas_fring,OP.mpc_antennas_clcal,OP.mpc_refant,OP.mpc_timer,OP.mpc_aparm,OP.mpc_dparm,OP.mpc_solint,uvdata.table_highver('FG'),OP.mpc_suba)
	AT.possm(uvdata,sources=OP.calibrators,gainu=uvdata.table_highver('CL'),aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal_After_mp',antennas=OP.possm_antennas)
	AT.possm(uvdata,sources=OP.target,gainu=uvdata.table_highver('CL'),aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllTargets_After_mp',antennas=OP.possm_antennas)
#
if OP.setjy:
	AT.setjy(uvdata)

if OP.bandpass_correction:
	uvdata.zap_table('BP',-1)
	cl_hv = uvdata.table_highver('CL')
	if not OP.bp_ichansel:
		mch=AF.max_ch(uvdata)
		ichansel=[]
		for i in range(8):
			ichansel.append([3,mch-3,1,i+1])
	AT.bpass(uvdata,antennas=OP.antennas,gainu=cl_hv,solint=0,refant=OP.bp_refan,cals=OP.bp_cal,ichansel=ichansel,bpassprm=OP.bp_bpassprm,suba=OP.bp_suba,fgv=0)
	AT.possm(uvdata,sources=OP.bp_cal,bpv=1,doband=2,fgv=0,gainu=cl_hv,aparm=[0,0,0,0,0,0,0,2,1],plotname='BP_Bandpass')
#	AT.possm(uvdata,sources=OP.calibrators,bpv=1,doband=2,timer=OP.bp_timer,fgv=0,solint=-1,gainu= cl_hv,plotname='BP_ACal')

#
if OP.runacscl:
    cl_hv = uvdata.table_highver('CL')
    if OP.use_usuba:
        sn_aa = []
        maxan = uvdata.table_highver('AN')
        for i in range(maxan):
            AT.acscl(uvdata,doband=1,bpv=1,suba=i+1)
            sn_hv = uvdata.table_highver('SN')
            AT.snsmo(uvdata,samptype='MWF',smotype='BOTH',doblank=1,dobtween=1,inv=sn_hv,outv=sn_hv+1,refant=OP.bp_refan,suba=i+1)
            sn_aa.append(uvdata.table_highver('SN'))
            AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='AMP',stokes='',sources=[],opcode='ALST',plotname='acscl_'+str(i+1))

        AT.clcal(uvdata,gainv=cl_hv,gainu=cl_hv+1,snv=sn_aa[0],inv=sn_aa[-1],refant=AF.get_antenna_number(uvdata,OP.refant_global),doblank=1,dobtween=1,interpol='self',suba=-32000)
    else:
        AT.acscl(uvdata,doband=2,bpv=1)
        sn_hv = uvdata.table_highver('SN')
        AT.snsmo(uvdata,samptype='MWF',smotype='BOTH',doblank=1,dobtween=1,inv=sn_hv,outv=sn_hv+1,refant=OP.bp_refan)
        sn_aa.append(uvdata.table_highver('SN'))
        AT.clcal(uvdata,gainv=cl_hv,gainu=cl_hv+1,snv=uvdata.table_highver('SN'),refant=AF.get_antenna_number(uvdata,OP.refant_global),doblank=1,dobtween=1,interpol='self')
        AT.snplt(uvdata,ine='SN',inv=uvdata.table_highver('SN'),optype='AMP',stokes='',sources=[],opcode='ALST',plotname='acscl')
#

if OP.runapcal:
	gainv = uvdata.table_highver('CL')
	gainu = gainv+1
	if OP.use_usuba:
		sn_aa = []
		maxan = uvdata.table_highver('AN')
		for i in range(maxan):
			if OP.opacity:
				AF.do_apcal(uvdata,aparm=OP.apcal_aparm,tyv=OP.apcal_tyv,dofit=OP.apcal_dofit,tau0=OP.apcal_tau0,trecvr=OP.apcal_trecvr,opcode=OP.apcal_opcode,calin=OP.apcal_calin,solint=OP.apcal_solint,suba=i+1,savelog=True)
			else:
				AT.apcal(uvdata,suba=i+1)
			sn_hv = uvdata.table_highver('SN')
			if OP.smooth_apcal_sn:
				AT.snsmo(uvdata,antennas=OP.antennas,smotype='BOTH',samptype='BOX',bparm=[1,0],doblank=1,dobtween=1,inv=sn_hv,outv=0,refant=AF.get_antenna_number(uvdata,OP.refant_global),suba=i+1)
			sn_aa.append(uvdata.table_highver('SN'))
		print('should run clcal now')
		AT.clcal(uvdata,gainv=gainv,gainu=gainu,snv=sn_aa[0],inv=sn_aa[-1],refant=AF.get_antenna_number(uvdata,OP.refant_global),doblank=1,dobtween=1,interpol='SELF',suba=-32000)
		AT.snplt(uvdata,ine='SN',inv=sn_aa[-1],optype='AMP',stokes='HALF',opcode='ALSI',plotname='apcal')
		logger.info('APCAL finalized. Please check output files:\nSN%d, CL%d',sn_aa[-1],uvdata.table_highver('CL'))
	else:
		cl_aa=[]
		if OP.opacity:
			AF.do_apcal(uvdata,aparm=OP.apcal_aparm,tyv=OP.apcal_tyv,dofit=OP.apcal_dofit,tau0=OP.apcal_tau0,trecvr=OP.apcal_trecvr,opcode=OP.apcal_opcode,calin=OP.apcal_calin,solint=OP.apcal_solint,savelog=True)
		else:
			AT.apcal(uvdata)
		sn_hv = uvdata.table_highver('SN')
		if OP.smooth_apcal_sn:
			AT.snsmo(uvdata,antennas=OP.antennas,smotype='BOTH',samptype='BOX',bparm=[1,0],doblank=1,dobtween=1,inv=sn_hv,outv=0,refant=AF.get_antenna_number(uvdata,OP.refant_global))
		sn_hv = uvdata.table_highver('SN')
		AT.clcal(uvdata,gainv=gainv,gainu=gainu,snv=sn_hv,refant=AF.get_antenna_number(uvdata,OP.refant_global),doblank=1,interpol='SELF')
		AT.snplt(uvdata,ine='SN',inv=sn_hv,optype='AMP',stokes='HALF',opcode='ALSI',plotname='apcal')
		logger.info('APCAL finalized. Please check output files:\nSN%d, CL%d',sn_hv,uvdata.table_highver('CL'))
#	AT.possm(uvdata,gainu=uvdata.table_highver('CL'),aparm=[0,1,0,5,-200,200,0,0,1,0])


if OP.get_best_solint:
	if OP.st_gainu:
		cl_hv=OP.st_gainu
	else:
		cl_hv= uvdata.table_highver('CL')
	st_timer=[]
	for i in OP.st_scan:
		st_timer.append(AF.scantime(uvdata,i))

	for i in range(len(OP.st_refant)):
		AF.derive_solint(uvdata,timer=st_timer[i],refant=OP.st_refant[i],gainu=cl_hv,plotname=OP.st_plotname[i],antennas=OP.st_antennas,aparm=OP.st_aparm,suba=OP.st_suba,solint=OP.st_solint)

if OP.global_fring:
	get2n = False
	if OP.gf_scan:
		OP.gf_timer=[]
		for i in OP.gf_scan:
			OP.gf_timer=AF.scantime(uvdata,i)
	if OP.get2n:
	#	imfile=AIPSImage(OP.cmap_name,'CMAP',uvdata.disk,1)
	#	if imfile.exists():
	#		get2n=[OP.cmap_name,'CMAP',uvdata.disk,1]
	#		logger.info('Using clean image')
	#	else:
		catnr=AT.imlod(OP.cmap_name,uvdata.disk,OP.cmap_file)
		get2n=AT.getndata(uvdata.disk,catnr)
	if OP.use_bp:
		if OP.smooth_gf_sn:
			cl_hv_af=AF.do_global_fring(uvdata,cals=OP.gf_cals,sources=OP.gf_sources,bpv=1,doband=2,timer=OP.gf_timer,dofit=OP.gf_dofit,antennas=OP.gf_antennas,aparm=OP.gf_aparm,dparm=OP.gf_dparm,refant=OP.gf_refant,interpol=OP.gf_interpol,search=OP.gf_search,solint=OP.gf_solint,solsub=OP.gf_solsub,gainu=OP.gf_gainu,dosnsmo=OP.smooth_gf_sn,smopa=[OP.smooth_gf_doblank,OP.smooth_gf_dobtween,OP.smooth_gf_bparm,OP.smooth_gf_cparm],get2n=get2n)
		else:
			cl_hv_af=AF.do_global_fring(uvdata,cals=OP.gf_cals,sources=OP.gf_sources,bpv=1,doband=2,timer=OP.gf_timer,dofit=OP.gf_dofit,antennas=OP.gf_antennas,aparm=OP.gf_aparm,dparm=OP.gf_dparm,refant=OP.gf_refant,interpol=OP.gf_interpol,search=OP.gf_search,solint=OP.gf_solint,solsub=OP.gf_solsub,gainu=OP.gf_gainu,get2n=get2n)

	else:
		cl_hv_af=AF.do_global_fring(uvdata,cals=OP.gf_cals,sources=OP.gf_sources,timer=OP.gf_timer,dofit=OP.gf_dofit,antennas=OP.gf_antennas,aparm=OP.gf_aparm,dparm=OP.gf_dparm,refant=OP.gf_refant,interpol=OP.gf_interpol,search=OP.gf_search,solint=OP.gf_solint,solsub=OP.gf_solsub,solmin=OP.gf_solmin,gainu=OP.gf_gainu,dosnsmo=OP.smooth_gf_sn,get2n=get2n,dobtween=OP.gf_dobtween,doblank=OP.gf_doblank,suba=OP.gf_suba,smopa=[OP.smooth_gf_doblank,OP.smooth_gf_dobtween,OP.smooth_gf_bparm,OP.smooth_gf_cparm])

	try:
		if OP.do_ad_mpc:
			for i,admp in enumerate(OP.ad_mpc_cals):
				logger.info('#'*20)
				cl_hv=uvdata.table_highver('CL')
				logger.info('\nStarting Fring for additional manual-phase-cal #%d\n',i+1)
				AT.fring_instr(uvdata,cals=[admp],antennas=OP.ad_mpc_antennas_fring[i],timer=OP.ad_mpc_timer[i],refant=OP.ad_mpc_refant[i],aparm=OP.ad_mpc_aparm[i],dparm=OP.ad_mpc_dparm[i],solint=OP.ad_mpc_solint[i],suba=OP.ad_mpc_suba[i],fgv=0,gainu=cl_hv,snv=0,bpv=0,doband=-1)
				AT.clcal(uvdata,sources=OP.ad_mpc_sources[i],cals=[admp],antennas=OP.ad_mpc_antennas_clcal[i],gainv=cl_hv,gainu=cl_hv+1,snv=uvdata.table_highver('SN'),refant=OP.ad_mpc_refant[i],interpol='2PT',timer=OP.ad_mpc_timer_clcal[i])
				logger.info('Manual-phase-cal #%d done\nCurrent tables: CL%d, SN%d',i+1,uvdata.table_highver('CL'),uvdata.table_highver('SN'))
	except:
		logger.warning('No additional manual phasecal done.\n Process without.')

	#logger.info('CL tables produced in global fring are:{}'.format(cl_hv_af))
	if OP.make_fring_tests:
		cl_af = cl_hv_af
	#AT.possm(uvdata,sources=[],gainu=uvdata.table_highver('CL'),aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All',antennas=OP.possm_antennas)

#

if OP.complexbp:
	cl_hv = uvdata.table_highver('CL')

	mch=AT.max_ch(uvdata)
	ichansel=[]
	for i in range(8):
		ichansel.append([3,mch-3,1,i+1])
	AT.bpass(uvdata,antennas=OP.antennas,gainu=cl_hv,solint=0,refant=OP.bp_refan,cals=[OP.bp_cal],ichansel=ichansel,bpassprm=OP.bp_bpassprm,fgv=0)

#	AT.bpass(uvdata,antennas=OP.antennas,gainu=cl_hv,solint=0,refant=OP.bp_refan,cals=[OP.bp_cal],ichansel=ichansel,bpassprm=[0,0,0,0,0,0,0,0,1,3,0],fgv=0)
	AT.possm(uvdata,sources=[OP.bp_cal],bpv=uvdata.table_highver('BP'),doband=2,fgv=0,gainu=cl_hv,aparm=[0,0,0,0,0,0,0,2,1],plotname='BP_new_Bandpass')
	AT.possm(uvdata,sources=OP.calibrators,bpv=uvdata.table_highver('BP'),doband=2,fgv=1,solint=-1,gainu= cl_hv,aparm=[0,0,0,0,-200,200,0,0,1,0],plotname='BP_new_ACal')


if OP.finalize:
	if OP.make_fring_tests:
		cl_hv = cl_aa
	else:
		cl_hv = [uvdata.table_highver('CL')]
	for clh in cl_hv:
		if OP.use_bp:
			bpv=uvdata.table_highver('BP')
			AT.split(uvdata,bpv=bpv,doband=2,sources=[],gainu=clh,fgv=1,antennas=OP.split_antennas,outd=OP.outdisk,aparm=OP.split_aparm)
		else: 
			#I added the loop afterwards and commented the else. to change it back to before, the split is under the else, the rest below has to be one tab removed.
		#for jj in [15,clh]:
			AT.split(uvdata,sources=[],gainu=clh,fgv=1,antennas=OP.split_antennas,outd=OP.outdisk,aparm=OP.split_aparm)
		data_rows=[ac for ac in AIPSCat(OP.outdisk)[OP.outdisk]]
		t=Table(rows=data_rows)
		sources=[]
		if sources==[]:
			sources = uvdata.sources
		for source in sources:
			imcat=t[t['name']==source]
			if len(imcat)>0:
				imcat=imcat[imcat['seq']==np.max(imcat['seq'])]
				imdata=[str(imcat['name'][0]),str(imcat['klass'][0]),OP.outdisk,int(imcat['seq'][0])]
				AT.fittp(imdata,outfile=imdata[0]+'_'+OP.date+'_CL'+str(clh)+'.uvfits')	

if OP.make_final_possm_plots:
	cl_hv=uvdata.table_highver('CL')
	if OP.use_bp:
		bpv=uvdata.table_highver('BP')
		#bpv=1
		AT.possm(uvdata,bpv=bpv,doband=2,sources=[],gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All',antennas=OP.possm_antennas)
		#AT.possm(uvdata,bpv=1,doband=2,sources=OP.calibrators,gainu=cl_hv-6,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal',antennas=OP.possm_antennas)
		#AT.possm(uvdata,bpv=1,doband=2,sources=OP.calibrators,gainu=cl_hv-7,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal',antennas=OP.possm_antennas)
#		AT.possm(uvdata,bpv=1,doband=2,sources=OP.calibrators,gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal_bp1',antennas=OP.possm_antennas)
	else:
#		AT.possm(uvdata,sources=[],gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All_suba1',antennas=OP.possm_antennas,suba=1)
#		AT.possm(uvdata,sources=[],gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All_suba2_solint2',antennas=OP.possm_antennas,suba=2,solint=1,timer=[0,17,0,0,0,19,30,0])
#		AT.possm(uvdata,sources=[],gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All_suba2_solint1',antennas=OP.possm_antennas,suba=1,solint=1,timer=[0,17,0,0,0,19,30,0])

#		AT.possm(uvdata,sources=[],gainu=cl_hv,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='ALL',antennas=OP.possm_antennas)
		AT.possm(uvdata,sources=['0954+658','1156+295','M84'],gainu=1,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='ALLpol',antennas=OP.possm_antennas,stokes='FULL')
#		AT.possm(uvdata,sources=[],gainu=cl_hv-7,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All_suba2',antennas=OP.possm_antennas,suba=2)
#		AT.possm(uvdata,sources=[],gainu=cl_hv-11,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='All_suba2',antennas=OP.possm_antennas,suba=2)


	#	AT.possm(uvdata,bpv=0,doband=-1,sources=OP.calibrators,gainu=cl_hv-2,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal',antennas=OP.possm_antennas)
	#	AT.possm(uvdata,bpv=0,doband=-1,sources=OP.calibrators,gainu=cl_hv-3,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal',antennas=OP.possm_antennas)
	#AT.possm(uvdata,sources=OP.calibrators,gainu=3,aparm=[0,1,0,0,-200,200,0,0,1,0],plotname='AllCal',antennas=OP.possm_antennas)
	#AT.possm(uvdata,sources=[],timer=[0,0,0,0],solint=-1,gainu= cl_first_cal_steps-1,aparm=[0,1,0,0,-200,200,0,1,1,0],plotname='AC')
	#AT.possm(uvdata,sources=[],timer=[0,0,0,0],solint=-1,gainu= cl_first_cal_steps,aparm=[0,1,0,0,-200,200,0,1,1,0],plotname='AC')


AIPS.log.close()
