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
outdisk				= 1 # to write to aips disk 1
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
t_int		= float(0.5)		# correlator integration time in seconds.
cl_table_interval	= 0.1	# desired CL table interval in minutes
digicor = 1
###############
# Calibrators #
###############
calibrators = ['XXX']	# specify all calibrators	
target			= 'XXX' # target source
cal_instr		= 'XXX'	# phase-cal calibrator
cal_bp			= ['XXX'] 	# bandpass calibrator(s)
refant_global	= 'XX' 	# typical reference antenna
antennas    = [0]#If an antennas should not be used for the whole calibration or only a set of them, put it here
#
#####################
# Calibration files #
#####################
antabfile	=[aips_cf+'XXX.antab']
weatherfile =aips_cf+'XXX.wx'

# In case vlog should be used to extract info about VLBA antennas give here the path to the vlba-cal file
vlbacal='xxxcal.vlba' #vlbacal file for use with vlog 
vlba_cal_out=aips_cf+'vlog/' # output dir for vlog
if not os.path.isdir(vlba_cal_out):
	os.makedirs(vlba_cal_out)
vlog_out = vlba_cal_out+experiment #outputname for vlog output files 

# Direcotry to Pulsecal calibration parameters as it would be after VLOG did run.
# If VLOG is not run, give the path to the file here
pcfile= vlba_cal_out+experiment+'.PCAL'
pcal_calibrator = cal_instr #Specify a calibrator for Pulsecal if desired
pcal_timer = []
pcal_antennas=[]
#
####################
# additional flags #
####################
# List of flag files to be loaded.
# More than one can be used. Default is VLBA flag file witten during VLOG
flagfile= [vlba_cal_out+experiment+'.FLAG']
# Flags that should by allied using uvflg. Here are shown a few examples
fg_bif		= [1,1,1,8]
fg_eif		= [0,0,0,8]
fg_bchan	= [1,1,1,1]
fg_echan  = [0,0,0,0]
fg_timer	= [[0],[0],[0,10,15,0,0,10,22,20],[0,16,40,0,0,16,50,0]]
fg_antennas=[3,15,11,11]
fg_stokes	= ['','1001','','RR']
fg_reason	= ['flag dbbc2 ef','YS only LL','bad scan','bad if']
# Another example: here one IF should be flagged for several scans:
#mh_timer=[[0,10,15,0,0,10,22,20],[0,10,30,0,0,10,37,20]]
#for t in mh_timer:
#	fg_bchan.append(1)
#	fg_echan.append(0)
#	fg_bif.append(8)
#	fg_eif.append(8)
#	fg_timer.append(t)
#	fg_antennas.append(8)
#	fg_stokes.append('RR')
#	fg_reason.append('bad if')
#
########################
# Parameters for swpol #
########################
swpol_antennas=[0]
########################
# Parameters for tabed #
########################
# Here as example to correct mounting for PV, which was antennas 14
tabed_ine='AN'
tabed_optype='REPL'
tabed_keyvalue=[5,0]
tabed_aparm=[5,0,0,4,4,14,0]
tabed_inv=1
tabed_outv=1

#########################################
# Parameters for manual phase cal steps #
#########################################
"""
The functions are written in a way to automatically do more than one fring run in a row, in case there is not the one and only good scan to align the phases between IFs. The Information has to be given in Arrays. E.g.
calibrator      = ['3C84','3C84','3C279']
antennas_fring  = [[1,2,3,5,6],[4,5,7,8,10],[9,10,11,12,13,14]
refant          = [5,5,10]
antennas_clcal  = [[1,2,3,5,6],[4,7,8,10],[9,11,12,13,14]]  #Probably you want to exclude refant antennas that had already been corrected before
suba            = [1,1,2]
timer           = [[0,9,49,1,0,9,53,1],[0,9,33,1,0,9,37,0],[0,16,40,0,0,16,50,0]]
aparm           = [[2,0,0,0,0,2,4.5,0],[2,0,0,0,0,2,5,0],[2,0,0,0,0,2,7,0]]
dparm           = [[1,400,400,0,0,0,1,1,1],[1,400,400,0,0,0,1,1],[1,400,400,0,0,0,1,1]
solint          = [-1,-1,4]
"""
mpc_calibrator    = []
mpc_antennas_fring= []
mpc_antennas_clcal= []
mpc_refant        = []
mpc_timer         = []
mpc_aparm         = []
mpc_dparm         = []
mpc_solint        = []
#mpc_suba         = []

if len(mpc_calibrator)==len(mpc_antennas_fring)==len(mpc_antennas_clcal)==len(mpc_refant)==len(mpc_timer)==len(mpc_aparm)==len(mpc_dparm)==len(mpc_solint):
	sys.stdout.write('Parameters for manual phase cal set correctly\n')
else:
	sys.stdout.write('While setting the Parameters for manual phase cal you forgot something.\nPlease check.\n')
	sys.exit()
#
###############
# solint test #
'''
Parameters for finding the best solution interval. 
If this is run plot files are created. There can be made tests on several scans in a row, therefore please give parameters as arrays.
run by setting get_best_solint=True at the end of the file.
This functionality has not been tested very extensively. Be careful
'''
st_refant		= [14,14,14,5,5,5] #reference antenna to be used
st_gainu		= 11	#CL table to be used
st_scan			= [1,3,6,24,33,50] #At the moment it is prefered to give a scan number
st_plotname	= []
for i in range(len(st_refant)):
	    st_plotname.append('_ref'+str(st_refant[i])+'_scan'+str(st_scan[i]))
st_solint		= [0.05,4] #The range of solution intervals which should be tested
st_snr_cut	= 4.1 #SNR cutoff during FRING
#
################
# Global Fring #
################
'''
There is the possibility to run several global fringes after each other. E.g to do a test. 
Therefore the structure of the parameters has to be an array as seen below.
'''
gf_scan			= False # it is also possible to give the scan number and not the time range
gf_timer		= [[0,0,0,0]] # please fill in the time range
gf_cals     = [[]]
gf_sources  = [[]]
gf_aparm		= [[2,0,0,0,1,2,4.0,0,1,0]]
gf_dparm		= [[1,400,400,1,0,0,0,1,0]]
gf_antennas	= [[antennas]]
gf_refant		= ['GB']
gf_dofit		= [[0]] #fit for all antennas
gf_search		= [['GB','PV','FD','EB','LA','ON']]
gf_interpol	= ['2PT']
gf_solint		= [4.0]
gf_gainu		= [0] # to use the highest number CL table during fringe

#####################################
# If a clean image should be loaded #
#####################################
get2n				= False # if an image file should be used
cmap_file = aips_if+''
cmap_name = ''
##################################
# If Fring SN should be smoothed #
##################################
'''
In case the SN table resulting from FRING should be smoothed with SNSMO
Please check input parameters carefully.
'''
smooth_gf_sn= False
smooth_gf_bparm=[0,0.05,0.05,0.05,0]
smooth_gf_cparm=[0,0,0.05,0.05,0,0,0,100,500,0]
smooth_gf_doblank     = -1
smooth_gf_dobtween    = -1

'''
There are many times there are remaining jumps between IFs.
By setting parameter fo ad_mpc an additional run of manual phase calibration can be done.
Below is just an example for three additional runs
'''
ad_mpc_cals           =['','','']
ad_mpc_sources        = [[],[],[]]
ad_mpc_antennas_fring =[[5,2,14,15],[8,11,14],[4,7]]
ad_mpc_antennas_clcal =[[14,2,15],[11,8],[4]]
ad_mpc_refant         =[5,14,7]
ad_mpc_timer          =[[0,16,40,0,0,16,50,0],[0,9,45,0,0,9,52,0],[0,15,40,0,0,15,50,0]]
ad_mpc_timer_clcal  =[[0,0,0,0],[0,0,0,0],[0,0,0,0]]
ad_mpc_aparm          =[[2,0,0,0,0,2,5,0],[2,0,0,0,0,2,5,0],[2,0,0,0,0,2,5,0]]
ad_mpc_dparm          =[[1,0,0,1,0,0,0,1,0],[1,0,0,1,0,0,0,1,0],[1,0,0,1,0,0,0,1,0]]
ad_mpc_solint         =[-1,-1,-1]



if len(gf_timer)==len(gf_aparm)==len(gf_dparm)==len(gf_refant)==len(gf_interpol)==len(gf_search)==len(gf_solint):
	pass
else:
	sys.stdout.write('Setting the Parameters for global fring fit you forgot something.\nPlease check.\n')
	sys.exit()

###########################
# Parameters for bandpass #
###########################
bp_refan		= 5
bp_cal			= cal_bp
bp_ichansel =	False
bp_bpassprm = [1,2,0,0,1,0,0,0,1,6,0]
#
###########################
# Parameters for APCAL	  #
###########################
apcal_dofit	= 15*[1] #asuming 15 antennas
apcal_tyv		= 1
apcal_tau0	= 15*[0.1]
apcal_trecvr= 30*[100]
apcal_calin	= weatherfile
apcal_aparm	= [1,0,0,1,4,1,3]
apcal_opcode='GRDR'
##############################
# Settings for finalize=True #
##############################
split_antennas= antennas # any antennas that should not be used in SPLIT
split_aparm		= [2,0,1,1]
####################
# inputs for POSSM #
####################
possm_antennas= antennas

###############################################
# SETTINGS - Which calibration to be run ### ##
###############################################
#
load_data						= False # loading the data from the correlator
swpol								= False # applying swpol 
dataPreProcessing		= False # doing some processing after loading: if specified swpol,tabed. Then MSORT and INDXR.
tabed								= False # using tabed, at the moment to correct mounting for PV
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
# if the above assumption of finding the first SN and CL tables produced by manual phasecal does nto work as suspected, the last CL and SN tables that should not be deleted can be specified here, typically this should not be needed
cl_cal_before_mp		= False 
sn_cal_before_mp		= False
#
get_some_listings		= False # currently that is used to print several versions of MATX
runpclod						= False # if True: PC tables are deleted if existent and PC-table is loaded
runpccor						= False # to run PCCOR: function not tested very well, yet
manual_phase_cal		= False # do do manual pahse cal as defined in 'mfc_' parameters
setjy               = False # run setij?
bandpass_correction	= False # apply a bandpass
use_bp              = False # if the bandpass table should be used
runacscl						= False # runacscl as recommended in the COOKBOOK after applying a bandpass
get_best_solint			= False # run test function to get the best solution interval
runapcal            = False # run apcal
smooth_apcal_sn     = False # smooth SN table resulting from antab
opacity             = False # fit for opacity
complexbp           = False # estimate a complex bandpass
global_fring				= False # do global fring
make_fring_tests		= False # if different settings for global fringing should be tested, set this parameter. If APCAL and Finalize (split and fittp) is 'TRUE' this will be done for each CL table resulting from the different fring settings. Just give a list to the 'gf_' parameters specifying the different settings for each fring run. If make_fring_tests='False' only the highest CL table will be further processed, even if a list of entries is given to 'gf_'
finalize						= False # finalize the calibration, by applying SPLIT and FITTP
make_final_possm_plots= False # produce a large set of POSSM plots. 
