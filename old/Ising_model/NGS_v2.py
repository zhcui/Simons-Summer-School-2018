#!/usr/bin/env python
#Fsqt_v3.py    A script for F_s(q,t) calculation, should be combined with get_pos.sh and LAMMPS
#              You should specify the # of particles, lattice constant, first peak of RDF, begin and end index for calculation and the temp you want calculate.
#              get_pos.py is to treat the LAMMPS dummp file to get the posistion informations. You should specify the lines you want to hold.
#              This script will generate the F_s(q,t)-vs-logt data in plot_fsqt file.

import os
import numpy
import scipy
import string
import math
import cmath
from math import floor
from numpy import linspace

### FUNCTION ###

def Nint(x):
	if x<0:
		y=-floor(0.5-x)
	else:
		y=floor(0.5+x)
	
	return y

def read_pos(index):
	xvec=[]
	yvec=[]
	zvec=[]
	file1=open("pos_%s.dat" % index, "r")
	while 1:
		line=file1.readline()
		if not line:
			break
		line=line.split()
		xvec.append(string.atof(line[0]))
		yvec.append(string.atof(line[1]))
		zvec.append(string.atof(line[2]))
	file1.close()
	return [xvec,yvec,zvec]
	
def fsqt(q,start,end,length,step,latvec):
	if end+length>len(q):
		print "end+length OUT OF RANGE!!"
		return 0
	disp=[]
	for i in range(start,end+1,step):
		rj=q[i]
		r0=q[i]
		for j in range(1,length+1):
			p2=q[i+j]
			rj=p2+Nint((rj-p2)/float(latvec))*latvec
		disp.append(rj-r0)
	return disp


### MAIN ###
	
os.system('rm -f plot_NGS_v2')
os.system('rm -f plot_MSD_v1')
qq=2*math.pi/1.0
q0=[qq,qq,qq]
latconst=6.50508
begin=0
end=4000
th=15000
step=100
numpart=256
timestep=10

linsp=linspace(0,math.log10(th),num=40)
linsp=list(linsp)
for w in range(0,len(linsp)):
	linsp[w]=int(floor(10**(linsp[w])))


ar_q = [read_pos(index) for index in range(1, numpart+1)]

for time in linsp:
#	SUM=0
	SUM2=0
	SUM4=0
	for index in range(1,numpart+1):
		q=ar_q[index-1]
	#time=30
#	print len(q[0])
		zx=fsqt(q[0],begin,end,time,step,latconst)
		zy=fsqt(q[1],begin,end,time,step,latconst)
		zz=fsqt(q[2],begin,end,time,step,latconst)
#	print(zx[:10])
	#zy=fsqt(q[1],len(q[0])-time-52,len(q[0])-time-2,time,18.05)
	#zz=fsqt(q[2],len(q[0])-time-52,len(q[0])-time-2,time,18.05)
#	print len(zx)
#	print "\n"
#	print zx
#	print "\n"
#	print zy
#	print "\n"
#	print zz
#	print "\n"
		im=1j
#		sum_phase=0
		sum_r4=0
		sum_r2=0
		for k in range(0,len(zx)):
			vecr=[zx[k],zy[k],zz[k]]
			normsq=numpy.linalg.norm(vecr)**2
			normsq2=numpy.linalg.norm(vecr)**4
#			ph=cmath.exp(numpy.dot(vecr,q0)*im)
#			print time
#			print " "
#			print normsq
#			print " \n"
#			sum_phase=sum_phase+ph
#			sum_phase=sum_phase+normsq
			sum_r2=sum_r2+normsq
			sum_r4=sum_r4+normsq2
#		ave=sum_phase/(len(zx)/((int)(max(1,math.floor(len(zx)/20)))))
		ave2=sum_r2/(float(len(zx)))
		ave4=sum_r4/(float(len(zx)))
		SUM2=SUM2+ave2
		SUM4=SUM4+ave4
	AVE2=SUM2/float(numpart)
	AVE4=SUM4/float(numpart)
#	print time
#	print " "
#	print AVE4
#	print " "
#	print AVE2
#	print " \n"
	NGS=0.6*(AVE4)/(AVE2**2)-1
	file2=open('plot_NGS_v2','a')
	file2.write(str(math.log10(time*timestep)))
	file2.write(" ")
#	file2.write(str(math.log10(SUM)))
	file2.write(str(NGS))
	file2.write("\n")
	file2.close()
	
	file3=open('plot_MSD_v1','a')
	file3.write(str(math.log10(time*timestep)))
	file3.write(" ")
	file3.write(str(math.log10(AVE2)))
	file3.write("\n")
	file3.close()


