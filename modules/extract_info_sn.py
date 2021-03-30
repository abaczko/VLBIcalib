#!/usr/bin/env python

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

import numpy as np
def extract_info_sn(uvdata,snver,minSNR=20):
	sntab=uvdata.table('SN',snver)
	nif=len(sntab[0].real1)
	nsou = len(uvdata.sources)
	nants=len(uvdata.antennas)
	RIF=3

	# How many polarizations do we have?
	if uvdata.header['naxis'][1]==1:
		NPOL = 1
	else:
		NPOL = 2

	fring_sol=[[] for i in range(nif)]
	fring_sol2=[[] for i in range(nif)]
	
	for row in sntab:
		if NPOL==1:
			if all([x>=minSNR for x in row.weight_1]):
				if all([x==row.refant_1[RIF] for x in row.refant_1]):
					for i in range(nif):
						fring_sol[i].append([row.time,row.time_interval,row.antenna_no,row.refant_1[i],row.delay_1[i],row.rate_1[i],row.mbdelay1,row.weight_1[i]])

		elif NPOL==2:
			if all([x>=minSNR for x in row.weight_1]):
				if all([x==row.refant_1[RIF]==row.refant_2[RIF] for x in row.refant_1]):
					for i in range(nif):
						fring_sol[i].append([row.time,row.time_interval,row.antenna_no,row.refant_1[i],row.delay_1[i],row.rate_1[i],row.mbdelay1,row.weight_1[i]])

			if all([x>=minSNR for x in row.weight_2]):
				if all([x==row.refant_1[RIF] for x in row.refant_2]):
					for i in range(nif):
						fring_sol2[i].append([row.time,row.time_interval,row.antenna_no,row.refant_2[i],row.delay_2[i],row.rate_2[i],row.mbdelay2,row.weight_2[i]])

  # ARRANGE SOLUTIONS BY ANTENNA:

  # Sort all data by time:


	for i in range(nif):
		if len(fring_sol[i])==0:
			print('There is no solution for Pol1 and IF{} and SN{} for the given SNR threshold. Try lowering it.'.format(i,snver))
			return
		fring_sol[i] = np.array(fring_sol[i])
		ind= np.argsort(fring_sol[i][:,0])
		fring_sol[i]=fring_sol[i][ind]

	if NPOL==2:
		for i in range(nif):
			if len(fring_sol2[i])==0:
				print('There is no solution for Pol2 and IF{} and SN{} for the given SNR threshold. Try lowering it.'.format(i,snver))
				return
			fring_sol2[i] = np.array(fring_sol2[i])
			ind= np.argsort(fring_sol2[i][:,0])
			fring_sol2[i]=fring_sol2[i][ind]
	
	ants=np.unique(fring_sol[RIF][:,2])
	ants=[int(x) for x in ants]
	nants_fit=len(ants)
	
	# Find the most common reference antenna:
	goodrefs=[[] for i in range(nants_fit)]
	bestref=[j for j in range(nants_fit)]
	j=0
	for a in ants:
		goodrefs[j]=[len(fring_sol[RIF][np.logical_and(fring_sol[RIF][:,3]==i+1,fring_sol[RIF][:,2]==a)]) for i in range(nants)]
		bestref[j]=np.argmax(np.array(goodrefs[j]))
		j+=1
	bestset=set(bestref)
	#bestset.remove(0)
	bestall = max(bestset, key=bestref.count)

	if NPOL==2:
		ants2=np.unique(fring_sol2[RIF][:,2])
		ants2=[int(x) for x in ants2]
		nants_fit2=len(ants2)
	
		goodrefs2=[[] for i in range(nants_fit2)]
		bestref2=[j for j in range(nants_fit2)]
		j=0
		for a in ants2:
			goodrefs2[j]=[len(fring_sol2[RIF][np.logical_and(fring_sol2[RIF][:,3]==i+1,fring_sol2[RIF][:,2]==a)]) for i in range(nants)]
			bestref2[j]=np.argmax(np.array(goodrefs2[j]))
			j+=1
		bestset=set(bestref2)
		#bestset.remove(0)
		bestall2 = max(bestset, key=bestref2.count)
	
	#refants = np.where(np.array(goodrefs))
	time		= [[[] for r in range(nants)] for j in range(nants_fit)]
	time_iv	=	[[[] for r in range(nants)] for j in range(nants_fit)]
	weight	=	[[[] for r in range(nants)] for j in range(nants_fit)] 
	dela		=	[[[] for r in range(nants)] for j in range(nants_fit)] 
	rate		=	[[[] for r in range(nants)] for j in range(nants_fit)] 
	mbdelay	=	[[[] for r in range(nants)] for j in range(nants_fit)] 
	antenna	=	[[[] for r in range(nants)] for j in range(nants_fit)] 
	refant	=	[[[] for r in range(nants)] for j in range(nants_fit)] 

	#j is antenna; r is refant
	j=0
	for a in ants:
		for r in range(nants):
			for i in range(nif):
				cond=np.logical_and(fring_sol[i][:,2] == a, fring_sol[i][:,3]==r+1)
				temp1 = np.array(fring_sol[i][cond,:])
				time[j][r].append(temp1[:,0])
				time_iv[j][r].append(temp1[:,1])
				antenna[j][r].append(temp1[:,2])
				refant[j][r].append(temp1[:,3])
				dela[j][r].append(temp1[:,4])
				rate[j][r].append(temp1[:,5])
				mbdelay[j][r].append(temp1[:,6])
				weight[j][r].append(temp1[:,7])
		j+=1

	if NPOL==2:
		time2		= [[[] for k in range(nants)] for i in range(nants_fit2)]
		time_iv2=	[[[] for k in range(nants)] for i in range(nants_fit2)]
		weight2	=	[[[] for k in range(nants)] for i in range(nants_fit2)] 
		dela2		=	[[[] for k in range(nants)] for i in range(nants_fit2)] 
		rate2		=	[[[] for k in range(nants)] for i in range(nants_fit2)] 
		mbdelay2=	[[[] for k in range(nants)] for i in range(nants_fit2)] 
		antenna2=	[[[] for k in range(nants)] for i in range(nants_fit2)] 
		refant2	=	[[[] for k in range(nants)] for i in range(nants_fit2)] 
		j=0
		#j is antenna; r is refant
		for a in ants2:
			for r in range(nants):
				for i in range(nif):
					cond=np.logical_and(fring_sol2[i][:,2] == a, fring_sol2[i][:,3]==r+1)
					temp2 = np.array(fring_sol2[i][cond,:])
					time2[j][r].append(temp2[:,0])
					time_iv2[j][r].append(temp2[:,1])
					antenna2[j][r].append(temp2[:,2])
					refant2[j][r].append(temp2[:,3])
					dela2[j][r].append(temp2[:,4])
					rate2[j][r].append(temp2[:,5])
					mbdelay2[j][r].append(temp2[:,6])
					weight2[j][r].append(temp2[:,7])
			j+=1

	return{'Pol1':{'time':time,'time_iv':time_iv,'weight':weight,'dela':dela,'rate':rate,'mbdelay':mbdelay,'antenna':antenna,'refant':refant,'goodrefs':goodrefs,'bestref':bestref,'bestall':bestall,'antennas':ants},'Pol2':{'time':time2,'time_iv':time_iv2,'weight':weight2,'dela':dela2,'rate':rate2,'mbdelay':mbdelay2,'antenna':antenna2,'refant':refant2,'goodrefs':goodrefs2,'bestref':bestref2,'bestall':bestall2,'antennas':ants2}}
