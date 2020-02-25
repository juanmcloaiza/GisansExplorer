import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import sys
import os
import tkFileDialog

print ""
print "GISANS viever for MARIA data - Jan 2019"
print "for questions contact a.koutsioumpas@fz_juelich.de"
print "\n"

minmax = raw_input("enter min/max (yes/no): ");
#minmax='no'
if minmax=='yes':
	gmin=input('Give min value:');
	gmax=input('Give max value:');

#qzvalue=0.0
#dqzvalue=0.05
#qzvalue=input("Give qz: ");
#dqzvalue=input("Give delta qz: ");

# Open file dialog
ftypes = [
    ('Measurement dat dile', '*.dat'), 
]

file_path_string = tkFileDialog.askopenfilename(filetypes=ftypes)

basename = file_path_string[:-4]

# Open and read the dat file
file = open(file_path_string, 'r') 

for line in file:
	if line.find('omega_value')>0:
		omega_line_list = line.split()
		omega=omega_line_list[3]
		print 'Angle of incidence (degrees): '+str(omega)
	if line.find('selector_lambda_value')>0:
		lambda_line_list = line.split()
		selector_lambda=lambda_line_list[3]
		print 'Neutron wavelength (Angtrom): '+str(selector_lambda)
	if line.find('.gz')>0:
		line_list = line.split()
		for word in line_list:
			if word.find('.gz')>0:
				file_path_string = word
	if line.find('.yaml')>0:
		line_list = line.split()
		for word in line_list:
			if word.find('.yaml')>0:
				file_path_stringB = word

file.close()
# Open and read the yaml file
line1=''
line2=''
line3=''
line4=''
file = open(file_path_stringB, 'r') 
for line in file:
	if line4.find('name: mon1')>0:
		line_list=line1.split()
		monitor=line_list[1]
		print 'Monitor conuts: '+str(monitor)+'\n'
	line4=line3
	line3=line2
	line2=line1
	line1=line

file.close()
# Define the position of the direct beam on the detector
yc=528 #530
zc=211 #220
#zc=512

y=[]
z=[]
I=[]

print 'Loading data and applying corrections... please wait..'+'\n'

# Unzip the file
os.system('cp '+file_path_string+' temp.gz')
os.system('gzip -d -q temp.gz')


sens = np.loadtxt("sensitivity_map", dtype='i', delimiter=' ')
#sens = np.loadtxt("sens.map", dtype='i', delimiter=' ')

num=0.0
meansens=0.0
for i in xrange(1,1024):
 	for j in xrange(1,1024):
 		if i>=162 and i<=862:
 			if j>=162 and j<=862:
 				meansens=meansens+float(sens[i,j])
 				num=num+1.0

meansens=meansens/num

inputd = np.loadtxt("temp", dtype='i', delimiter=' ')
os.system('rm temp')
#os.system('rm temp.gz')

#file = open('output.dat','w') 

# Addition important... qz=0 at the horizon of the sample

print(zc)
zc=zc+int((1990.0*np.tan(np.pi*float(omega)/180.0))/0.5755)
print(zc)

file = open(basename+"_GISANS.map", "w")

for i in xrange(1,1024):
 	for j in xrange(1,1024):
 		if i>=162 and i<=862:
 			if j>=162 and j<=862:
 				y.append((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))
 				z.append((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
 				if float(sens[i,j]) > 0.0:
 					I.append((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))
 					file.write(str((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))+'   '+str((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))+'   '+str((meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor))+'\n')
 				#file.write(str((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/162.0)*((yc-i)/np.sqrt(1990*1990+(0.6*(yc-i))*(0.6*(zc-i)))))+' '+str((-2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/162.0)+((zc-j)/np.sqrt(1990*1990+(0.6*(zc-j))*(0.6*(zc-j))))))+' '+str(float(input[i,j])/float(monitor))+'\n')

file.close()

min_value=min(I)+float(1)/float(monitor)


max_value=max(I)

if minmax=='yes':
	min_value=gmin
	max_value=gmax


plt.scatter(z,y,edgecolors='none',c=I, norm = LogNorm(vmin=min_value, vmax=max_value), marker = 's')
cbar = plt.colorbar()
plt.ylabel(r'$Q_{y}(\AA^{-1})$')
plt.xlabel(r'$Q_{z}(\AA^{-1})$') 
plt.show()


linecut = raw_input("Line cut (yes/no): ");

qdirection=' '
if linecut=='yes':
	qdirection=raw_input('qy or qz? ');


if qdirection=='qy' and linecut=='yes':
	qzvalue=input("Give qz: ");
	dqzvalue=input("Give delta qz: ");

	cut_qy=[0] * len(z)
	cut_qy_xaxis=[0] * len(z)
	for i in xrange(1,1024):
		for j in xrange(1,1024):
			if i>=162 and i<=862:
				if j>=162 and j<=862:
					cut_qy_xaxis[i]=((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))
					cqz=((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
					if cqz>=qzvalue-dqzvalue and cqz<=qzvalue+dqzvalue:
						cut_qy[i]=cut_qy[i]+(meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor)
	plt.loglog(cut_qy_xaxis,cut_qy,label="cut at "+r'$Q_{z}$= '+'{:06.4f}'.format(qzvalue)+'$ \AA^{-1})$'+'\n'+r'$\Delta Q_{z}$= '+'{:06.4f}'.format(dqzvalue)+'$ \AA^{-1})$')
	plt.legend()
	plt.ylabel(r'$I(A.U.)$')
	plt.xlabel(r'$Q_{y}(\AA^{-1})$')
	plt.show()


	file = open(basename+"_line_cut_qy.out", "w")
	for i in xrange(1,1024):
		if i>=162 and i<=862:
			file.write(str(cut_qy_xaxis[i])+'   '+str(cut_qy[i])+"\n")
	file.close()

if qdirection=='qz' and linecut=='yes':
	qyvalue=input("Give qy: ");
	dqyvalue=input("Give delta qy: ");

	cut_qz=[0] * len(y)
	cut_qz_xaxis=[0] * len(y)
	for i in xrange(1,1024):
		for j in xrange(1,1024):
			if i>=162 and i<=862:
				if j>=162 and j<=862:
					cut_qz_xaxis[j]=((2*np.pi/float(selector_lambda))*(np.sin(np.pi*float(omega)/180.0)+(0.5755*(j-zc)/np.sqrt(1990*1990+(0.5755*(zc-j))*(0.5755*(zc-j))))))
					cqy=((2*np.pi/float(selector_lambda))*np.cos(np.pi*float(omega)/180.0)*(0.5755*(yc-i)/np.sqrt(1990*1990+(0.5755*(yc-i))*(0.5755*(yc-i)))))
					if cqy>=qyvalue-dqyvalue and cqy<=qyvalue+dqyvalue:
						cut_qz[j]=cut_qz[j]+(meansens*float(inputd[i,j])/float(sens[i,j]))/float(monitor)
	plt.loglog(cut_qz_xaxis,cut_qz,label="cut at "+r'$Q_{y}$= '+'{:06.4f}'.format(qyvalue)+'$ \AA^{-1})$'+'\n'+r'$\Delta Q_{y}$= '+'{:06.4f}'.format(dqyvalue)+'$ \AA^{-1})$')
	plt.legend()
	plt.ylabel(r'$I(A.U.)$')
	plt.xlabel(r'$Q_{z}(\AA^{-1})$')
	plt.show()


	file = open(basename+"_line_cut_qz.out", "w")
	for i in xrange(1,1024):
		if i>=162 and i<=862:
			omega_value=cut_qz_xaxis[i]*float(selector_lambda)*180.0/(4.0*np.pi*np.pi)
			file.write(str(omega_value)+'   '+str(cut_qz_xaxis[i])+'   '+str(cut_qz[i])+"\n")
	file.close()
