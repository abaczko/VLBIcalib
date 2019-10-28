import os,time,datetime,sys
import logging
logger = logging.getLogger(__name__)
###############################################
# SETTINGS - aipsid,disks,outpaths,datainfo  ##
###############################################
#
AIPSLite	= False
aipsid		= xxx
aips_cf		= '' 	# location/directory of calibration files as antab,flagtables
aips_if		= '' # location of input fits files
#
######################
# Info on experiment #
######################
uvrawdatafile	= [aips_if+'XX.fits'] #
experiment		= 'XX' #name of file in aipscat
indata 				= 1
outname				= experiment
outclass			= 'UVDATA' #outclass
outdata 			= indata
outdisk				= 1
outseq				= 1
##############################################################
# Initialise Log #############################################
# All Aips messages will be printed to the specific files    #
# All inputs for tasks in ParselTongue will be saved in      #
# another file                                               #
##############################################################
#####################
# Info for logs #####
local_dir = os.getcwd()
now = datetime.datetime.now()
date =str(now.day)+'_'+str(now.month)+'_'+str(now.year)
#

aipslog = local_dir+'/AIPS_messages_'+experiment+'_'+date+'.log'
tasklog = local_dir+'/PT_inputs_'+experiment+'_'+date+'.log'

#######################################
# Create aips output folder ###########
#######################################
aips_out=local_dir+'/aips_out'
if not os.path.isdir(aips_out):
  os.makedirs(aips_out)
###############################################
## Settings for the calibration ###############
###############################################
#
solint	= float(120.)		# main solution interval in seconds
t_int		= float(0.5)		# correlator integration time in seconds.
cl_table_interval	= 0.1	# desired CL table interval in minutes
digicor = 1
###############
# Calibrators #
###############
calibrators = ['3C84','0224+069']	# specify all calibrators	
target			= 'NGC1052' # target source
cal_instr		= '3C84'	# phase-cal calibrator
cal_bp			= '' 	# bandpass calibrator
refant_global	= 'GB' 	# typical reference antenna
#
#####################
# Calibration files #
#####################
antabfile	=[aips_cf+'xxx.antab']

#in case vlog should be used to extract info about VLBA antennas
vlbacal='xxxcal.vlba' #vlbacal file for use with vlog 
vlba_cal_out=aips_cf+'vlog/' # output dir for vlog
if not os.path.isdir(vlba_cal_out):
	os.makedirs(vlba_cal_out)
vlog_out = vlba_cal_out+experiment #name for vlog output files 

#Pulsecal calibration parameters
pcfile		=aips_cf+'vlog/'+experiment+'.PCAL'
pcal_calibrator = cal_instr #Specify a calibrator for Pulsecal if desired
pcal_timer = []
pcal_antennas=[]
#
flagfile= [vlba_cal_out+experiment+'.FLAG'] #list of flag files to be loaded

####################
# additional flags #
####################
#Flags that should by allied using uvflg.
fg_bif		= [1,1,1,8]
fg_eif		= [0,0,0,8]
fg_bchan	= [1,1,1,1]
fg_echan  = [0,0,0,0]
fg_timer	= [[0],[0],[0,10,15,0,0,10,22,20],[0,16,40,0,0,16,50,0]]
fg_antennas=[3,15,11,11]
fg_stokes	= ['','1001','','RR']
fg_reason	= ['flag dbbc2 ef','YS only LL','bad scan','bad if']
mh_timer=[[0,10,15,0,0,10,22,20],[0,10,30,0,0,10,37,20]]
for t in mh_timer:
	fg_bchan.append(1)
	fg_echan.append(0)
	fg_bif.append(8)
	fg_eif.append(8)
	fg_timer.append(t)
	fg_antennas.append(8)
	fg_stokes.append('RR')
	fg_reason.append('bad if')

########################
# Parameters for swpol #
########################
swpol_antennas=[0]
########################
# Parameters for tabed #
########################
# here as example to correct mounting for PV, which was antennas 14
tabed_ine='AN'
tabed_optype='REPL'
tabed_keyvalue=[5,0]
tabed_aparm=[5,0,0,4,4,14,0]
tabed_inv=1
tabed_outv=1

#########################################
# Parameters for manual phase cal steps #
#########################################
mpc_calibrator		= ['3C84','3C84','0224+069']
mpc_antennas_fring= [[2,11,14,15],[14,8],[14,5]]
mpc_antennas_clcal= [[2,11,14,15],[8],[5]]
mpc_refant				= [14,14,14]
mpc_timer					= [[0,9,50,30,0,9,51,0],[0,9,51,0,0,9,52,0],[0,16,43,0,0,16,43,30]]
mpc_aparm					= [[2,0,0,0,0,2,5,0],[2,0,0,0,0,2,4.5,0],[2,0,0,0,0,2,4.7,0]]
mpc_dparm					= [[1,400,0,0,0,0,0,1,1,0],[1,50,0,0,0,0,0,1,1,0],[1,200,0,0,0,0,0,1,1,0]]
mpc_solint				= [0.5,1.0,0.5,4,1.0,1.5,1.5,4.0]#0.5]
#suba					= [1,1,1]
if len(mpc_calibrator)==len(mpc_antennas_fring)==len(mpc_antennas_clcal)==len(mpc_refant)==len(mpc_timer)==len(mpc_aparm)==len(mpc_dparm)==len(mpc_solint):
	sys.stdout.write('Parameters for manual phase cal set correctly\n')
else:
	sys.stdout.write('Setting the Parameters for manual phase cal you forgot something.\nPlease check.\n')
	sys.exit()
#
###############
# solint test #
'''
Parameters for finding the best solution interval. 
If this is run plot files are created. There can be made tests on several scans in a row, therefore please give parameters as arrays.
run by setting get_best_solint=True at the end of the file.
'''
st_refant = [14,14,14,5,5,5]
st_gainu = 11
st_scan=[1,3,6,24,33,50]
st_plotname=['_scan1_ref14','_scan3_ref14','_scan6_ref14','_scan24_ref5','_scan33_ref5','_scan50_ref5']
st_solint=[0.05,4]
st_snr_cut=4.1
#
################
# Global Fring #
################
'''
There is the possibility to run several global fringes after each other. E.g to do a test. 
Therefore the structure of the parameters has to be an array as seen below.
'''
gf_scan = False
gf_timer		=[[0,0,0,0]]
gf_aparm		= [[2,0,0,0,1,2,4.0,0,1,0]]
gf_dparm		= [[1,400,400,1,0,0,0,0]]
gf_antennas=[[]]
gf_refant	= ['GB']
gf_dofit = [[]]
gf_search  =[['GB','PV','FD','EB','LA','ON']]
gf_interpol= ['2PT']
gf_solint		= [0.15]
gf_gainu=[0]
gf_gainu=[11]

if len(gf_timer)==len(gf_aparm)==len(gf_dparm)==len(gf_refant)==len(gf_interpol)==len(gf_search)==len(gf_solint):
	pass
else:
	sys.stdout.write('Setting the Parameters for global fring fit you forgot something.\nPlease check.\n')
	sys.exit()

###########################
# Parameters for bandpass #
###########################
bp_refan		= 5
bp_cal			= cal_instr
bp_ichansel =	[]
for i in range(8):
	bp_ichansel.append([3,29,1,i])
bp_bpassprm = [1,2,0,0,1,0,0,0,1,6,0]
#
###########################
# Parameters for APCAL	  #
###########################
apcal_dofit	= 15*[1] 
apcal_tyv=2
apcal_tau0	= 15*[0.1]
apcal_trecvr= 30*[100]
apcal_calin	= aips_cf+'xxx.wx' #as an example
apcal_aparm	=[1,0,0,1,4,1,3]

##############################
# Settings for finalize=True #
##############################
split_antennas=[] # any antennas that should not be used in SPLIT
####################
# inputs for POSSM #
####################
possm_antennas=[]

###############################################
# SETTINGS - Which calibration to be run ### ##
###############################################
#
load_data						= False # loading the data from the correlator
swpol								= False # applying swpol 
tabed								= False # using tabed, at the moment to correct mounting for PV
dataPreProcessing		= False # doing some processing after loading: if specified swpol,tabed. Then MSORT and INDXR.
printBasicInfo			= False # print basic info as prtan, listr(scan), possm plot
delete_cal_tables		= False # if the calibration should be started over again, set this parameter to True. All SN tables, CL>1 and FG tables will be deleted
get_vlba_cal				= False # use vlog to print info from xxxvlba.cal file to individual files
load_antab					= False # load antab table(s)
load_flag_info			= False # load flag fable(s)
flag_outer_ch				= False # flag the first and last 3 channels
uvflg_additional		= False # apply uvflg according to fg_xxx parameters
usuba								= False # run usuba with opcode='auto'
correct_pang				= False # correct for PANG
correct_eop					= False # correct for EOP. The file will automatically be downloaded. Currently the Goddard server is down, therefore a copy as linked from aips.nrao is used.
runaccor						= False # to run ACCOR
clear_manual_phasecal= False # to clear all tables after ACCOR to start delay calibration over again. The history file will be searched to find the last SN and CL tables as were produced by ACCOR and deletes all SN and CL tables higher then these found
# Parameters for from which cl and sn table on tables should be deleted
cl_cal_before_mp		= False # if the above assumption of finding the first SN and CL tables produced by manual phasecal does nto work as suspected, the last CL and SN tables that should not be deleted can be specified here
sn_cal_before_mp		= False
#
get_some_listings		= False # currently that is used to print several versions of MATX
runpclod						= False # if True: PC tables are deleted if existent and PC-table is loaded
runpccor						= False # to run PCCOR: function not tested very well, yet
manual_phase_cal		= False # do do manual pahse cal as defined in 'mfc_' parameters
bandpass_correction	= False # apply a bandpass
runacscl						= False # runacscl as recommended in the COOKBOOK after applying a bandpass
get_best_solint			= False # run test function to get the best solution interval
global_fring				= False # do global fring
make_fring_tests		= False # if different settings for global fringing should be tested, set this parameter. If APCAL and Finalize (split and fittp) is 'TRUE' this will be done for each CL table resulting from the different fring settings. Just give a list to the 'gf_' parameters specifying the different settings for each fring run. If make_fring_tests='False' only the highest CL table will be further processed, even if a list of entries is given to 'gf_'
runapcal						= False # run APCAL
opacity							= False # fit for opacity, currently APCAL is run only one time. If you like to use the fitting results as new input parameters, you have to do this for now by hand. Just run APCAL again with specifying the fit values in the 'apcal_' parameters.
finalize						= False # finalize the calibration, by applying SPLIT and FITTP
make_final_possm_plots= False # produce a large set of POSSM plots. 