#!/usr/bin/env python
import sys 
import os
import gzip
import operator
import yaml
import cProfile
import pstats
import numpy as np
import subprocess
import tempfile
from string import *
from glob import *
import datetime
from math import *
from copy import  copy,deepcopy
import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import platform

NICOS=1
#replaced DejaVuSans by Courier due to segmentation faults when typing fi or simelar ...
#global variables for MCS6
H_MIN_CHANNEL=165
H_MAX_CHANNEL=860
V_MIN_CHANNEL=165
V_MAX_CHANNEL=860
ALL_CHANNELS=1024
#S1_left=20.907
#S1_right=26.87
#S2_left=22.95
#S2_right=29.15
S1_left=0.0
S1_right=0.0
S2_left=0.0
S2_right=0.0
DETECTORWIDTH=400.0
DIST_SAMP_DET=2093.0 #!! Gemessen am 15.02.2015
DIST_SAMP_DET=1910.0 #!! Gemessen am 01.05.2017 withnew Argon Box
DIST_SAMP_DET=1930.0
CHANNELWIDTH=DETECTORWIDTH/(H_MAX_CHANNEL-H_MIN_CHANNEL) #0.575
L3=4000 #Laenge der Kollimation
C=450  #Abstand Probe-S2
RES=0.001
GNUFONT="Courier"
GNUFONT="DejaVuSans"

sname2={0:"pf0af0", 1:"pf1af0", 2:"pf0af1", 3:"pf1af1"}
if NICOS==1:
  FS_DICT={0:"uu", 1:"du",2:"ud",3:"dd" }
  Pflipperdic={'up':1,'down':0}
  Aflipperdic={'up':0,'down':1}
  Analyzerdic={'in':1,'out':0,'none':3}
  #NIST_DICT={0:"# \"polarization\": \"--\"", 1:"# \"polarization\": \"++\"",2:"# \"polarization\": \"--\"",3:"# \"polarization\": \"++\""}
  NIST_DICT_POL1={0:"# \"polarization\": \"++\"", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"# \"polarization\": \"--\""}
  NIST_DICT_POL2={0:"# \"polarization\": \"--\"", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"# \"polarization\": \"++\""}
  NIST_DICT_POL3={0:"# \"unpolarized\":",1:"# \"unpolarized\":",2:"# \"unpolarized\":",3:"# \"unpolarized\":"}
  NIST_DICT_POL4={0:"# \"polarization\": \"+-\"", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"# \"polarization\": \"-+\""}
  NIST_DICT_POL5={0:"# \"polarization\": \"+-\"", 1:"# not existing with pol5 SF",2:"# not existing with pol5 SF",3:"# \"polarization\": \"-+\""}
  NIST_DICT_POL6={0:"# not existing with pol6 NSF", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"# not existing with pol6 NSF"}
  NIST_DICT_POL7={0:"# \"polarization\": \"+-\"", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"#not existing with pol7 full pol wizth one SF"}
  NIST_DICT_POL8={0:"#not existing with pol8 full pol wizth one SF", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"# \"polarization\": \"-+\""}
  NIST_DICT_POL9={0:"# \"polarization\": \"+-\"", 1:"# \"polarization\": \"--\"",2:"# \"polarization\": \"++\"",3:"# \"polarization\": \"-+\""}
  NIST_DICT=[NIST_DICT_POL1,NIST_DICT_POL2,NIST_DICT_POL3,NIST_DICT_POL4,NIST_DICT_POL5,NIST_DICT_POL6,NIST_DICT_POL7,NIST_DICT_POL8,NIST_DICT_POL9]
  # "polarization": "++"
else:
  FS_DICT={0:"ud", 1:"dd",2:"uu",3:"du" }
  Pflipperdic={'up':1,'down':0}
  Aflipperdic={'up':0,'down':1}
  Analyzerdic={'in':1,'out':0,'none':3}

SUMUP_ARRAY=[[[0 for i in range(ALL_CHANNELS+1)] for j in range(ALL_CHANNELS+1)] for fs in range(4)]


def usage():
  global DEBUG
  print( "usage: images.py [options] -base name -seq list")
  print()
  print( "options:")
  print( "          -help            Print this" )
  print( "          -base            Give the basename of the file" )
  print( "          -seq             List of number (space separated) of basename" )
  print( "          -div # # #...    List of numbers (space separated) as dividers between sequences" )
  print( "          -fread           Force reading full data set again" )
  print( "          -fpng            Force writing all png files" )
  print( "          -favi            Force writing movie file")
  print( "          -faiaf           Force writing aiaf-map" )
  print( "          -nopngcairo           Don't create jpeg images of detector images" )
  print( "          -noavi           Don't create avi-movie of detector images (pngcairo have to exist)" )
  print( "          -noaiaf          Don't create alpha_i/alpha_f map" )
  print( "          -aiafmm # #      Minimum und maximum value for aiaf plot" )
  print( "          -bglin #         Use linear background correction and # of points beside ROI, (0 switches of)")
  print( "          -roffset #       Offset of to the center channel at 512 (def: 0)")
  print( "          -roi # #         Give the width and height for the roi (def: 24 440)")
  print( "                           If only one number is give it's just width"  )
  print( "          -nmon #          Normalize reflectivity, aiaf-map: 0=Time 1=Mon1 2=Mon2 -1=Nothing counts/s")
  print( "                           Mon1 takes s1 and s2 into account, Mon2 only s2 and Time nothing else")
  print( "                           In the case of counts/s, normalize to the direct beam"  )
  print( "          -scale #         Scale reflectivity by given value (default value 1, -1 will try to find best scale value)")
  print( "          -sfc             Simple footprint( correction (divide by omega) (default off)" )
  print( "          -fc #            Testing!!! Footprint( correction with sample length in mm (defa))ult off)" )
  print( "          -HMIN            Set HMIN-Channel for measurement" )
  print( "          -mansel #        Manually set wavelength of selektor" )
  print( "          -divdet          Divide 2 det images by each other (0/1 2/3 4/5 etc.)" )
  print( "          -subdet          Substract two det images by each other (0/1 2/3 4/5 etc.)" )
  print( "          -rev             Revers the filenames om divdet and subdet" )
  print( "          -sumup           Sum all detector images up (for splitted gisans/timescan measurements)" )
  print( "          -noref           No reflectivity curve" )
  print( "          -maxi #          Maximum intensity for gnuplot detector images" )
  print( "          -mini #          Maximum intensity for gnuplot detector images" )
  print( "          -smaxi #         Maximum intensity for gnuplot detector sumup images" )
  print( "          -smini #         Maximum intensity for gnuplot detector sumup images" )
  print( "          -submaxi #       Maximum intensity for gnuplot detector subdet images" )
  print( "          -submini #       Minimum intensity for gnuplot detector subdet images" )
  print( "          -divmaxi #       Maximum intensity for gnuplot detector divdet images" )
  print( "          -divmini #       Minimum intensity for gnuplot detector divdet images" )
  print( "          -divsmaxi #      Maximum intensity for gnuplot detector divdet sumup images" )
  print( "          -divsmini #      Minimum intensity for gnuplot detector divdet sumup images" )
  print( "          -sens #          Read detector image as sesnitivity file (default: ~/bin/det_sens_img.gz)" )
  print( "          -nosens          Switch off sensitivity map correction of detector images." )
  print( "          -gisans          Switch to gisans mode." )
  print( "          -cut x1 y1 x2 y2 Export cut of sumup.img to gnuplot (x1,y1,x2,y2 in detector channels)" )
  print( "          -leak #          Substract from the spin flip map the leakage of the non spin flip curve" )
  print( "          -noq  \"string\"   Instead plotting reflectivity intensity versus q but 'parameter' e.g. magnet" )
  print( "          -genx            Export the data for genx data format instead parat" )
  print( "          -coh             Analyse the data of a soft matter experiment with polarizer and analyzer")
  print( "                           to measure the incoherent background with SF and NSF channel (Highly experimental)")
  print( "          -nicos           switch to NICOS mode")
  print( "          -pyfrid          switch to PYFRID mode (default)")
  print( "          -specfromaiaf    Get specular reflectivity from aiaf plot")
  print( "          -vert            Measure specular by rotating rx instead of omega")
  print()
  print()
  
  sys.exit(1)



def get_switches ():
  global DEBUG, H_MIN_CHANNEL, H_MAX_CHANNEL,NICOS
  global FS_DICT, Pflipperdic, Aflipperdic, Analyzerdic

  base=""
  sequence=[]
  divisors=[]
  window=[]
  logz=0
  nopng=0
  noavi=0
  noaiaf=0
  noref=0
  aiafmm=[0,0]
  bglin=-1
  roffset=[]
  roi= [24,440]
  fread=0
  fpng=0
  favi=0
  faiaf=0
  nmon=1
  scale=1
  sfc=0
  fc=0
  mansel=0
  divdet=0
  subdet=0
  sumup=0
  maxi=0
  mini=0
  smaxi=0
  smini=0
  submaxi=10
  submini=-10
  divmaxi=2
  divmini=.5
  divsmaxi=30
  divsmini=0
  sdet=1
  gisans=0
  rev=0
  leak=0
  noq=""
  genx=0
  coh=0
  specfromaiaf=0
  vert=0

  if(os.path.isfile("/home/maria/bin/sens_det_image.gz")==True):
    sens_det_image="/home/maria/bin/sens_det_image.gz"
    print( "Using default sensitivity map: /home/maria/bin/sens_det_image.gz")
  elif (os.path.isfile("/Users/mattauch/bin/sens_det_image.gz")==True):
    sens_det_image="/Users/mattauch/bin/sens_det_image.gz"
    print( "Using default sensitivity map: /Users/mattauch/bin/sens_det_image.gz")
  else:
    print( "Cannot find the default sens_det_image switching to -sdet 0")
    sens_det_image=""
    sdet=0
  if len(sys.argv) >= 1:
    for i in range(1,len(sys.argv)):
      if sys.argv[i] == "-base":
        base=sys.argv[i+1]
      if sys.argv[i] == "-seq":
        while (i+1 < len(sys.argv)) and sys.argv[i+1].isdigit():
          sequence.append(int(sys.argv[i+1]))
          i=i+1
      if sys.argv[i] == "-div":
        while (i+1 < len(sys.argv)) and (is_number(sys.argv[i+1])==True):
          divisors.append(float(sys.argv[i+1]))
          i=i+1
      if sys.argv[i] == "-cut":
        while (i+1 < len(sys.argv)) and (is_number(sys.argv[i+1])==True):
          window.append(float(sys.argv[i+1]))
          i=i+1
        print( "window:",window)
      if sys.argv[i] == "-logz":
        logz=1
      if sys.argv[i] == "-roffset":
#        roffset=int(sys.argv[i+1])
        while (i+1 < len(sys.argv)) and (is_number(sys.argv[i+1])==True):
          roffset.append(float(sys.argv[i+1]))
          i=i+1
      if sys.argv[i] == "-roi":
	roi[0]=int(sys.argv[i+1])
        if sys.argv[i+2].isdigit():
          roi[1]=int(sys.argv[i+2])
          i=i+2
        else:
          roi[1]=440
          i=i+1
      if sys.argv[i] == "-genx":
	genx=1
      if sys.argv[i] == "-fread":
	fread=1
      if sys.argv[i] == "-fpng":
	fpng=1
      if sys.argv[i] == "-favi":
	favi=1
      if sys.argv[i] == "-faiaf":
	faiaf=1
      if sys.argv[i] == "-nopng":
	nopng=1
      if sys.argv[i] == "-noref":
	noref=1
      if sys.argv[i] == "-noavi":
	noavi=1
      if sys.argv[i] == "-noaiaf":
	noaiaf=1
      if sys.argv[i] == "-sfc":
	sfc=float(sys.argv[i+1])
      if sys.argv[i] == "-divdet":
	divdet=1
      if sys.argv[i] == "-subdet":
	subdet=1
      if sys.argv[i] == "-rev":
        rev=1
      if sys.argv[i] == "-sumup":
	sumup=1
      if sys.argv[i] == "-vert":
	vert=float(sys.argv[i+1])
      if sys.argv[i] == "-specfromaiaf":
	specfromaiaf=1
      if sys.argv[i] == "-gisans":
	gisans=1
        noaiaf=1
        noavi=1
        noref=1
      if sys.argv[i] == "-nicos":
        FS_DICT={0:"uu", 1:"du",2:"ud",3:"dd" }
        Pflipperdic={'up':1,'down':0}
        Aflipperdic={'up':0,'down':1}
        Analyzerdic={'in':1,'out':0,'none':3}
        NICOS=1
      if sys.argv[i] == "-pyfrid":
        FS_DICT={0:"ud", 1:"dd",2:"uu",3:"du" }
        Pflipperdic={'up':1,'down':0}
        Aflipperdic={'up':0,'down':1}
        Analyzerdic={'in':1,'out':0,'none':3}
        NICOS=0
      if sys.argv[i] == "-leak":
        leak=float(sys.argv[i+1])
      if sys.argv[i] == "-maxi":
        maxi=float(sys.argv[i+1])
      if sys.argv[i] == "-mini":
        mini=float(sys.argv[i+1])
      if sys.argv[i] == "-smaxi":
        smaxi=float(sys.argv[i+1])
      if sys.argv[i] == "-smini":
        smini=float(sys.argv[i+1])
      if sys.argv[i] == "-submaxi":
        submaxi=float(sys.argv[i+1])
      if sys.argv[i] == "-submini":
        submini=float(sys.argv[i+1])
      if sys.argv[i] == "-divmaxi":
        divmaxi=float(sys.argv[i+1])
      if sys.argv[i] == "-divmini":
        divmini=float(sys.argv[i+1])
      if sys.argv[i] == "-divsmaxi":
        divsmaxi=float(sys.argv[i+1])
      if sys.argv[i] == "-divsmini":
        divsmini=float(sys.argv[i+1])
      if sys.argv[i] == "-fc":
        fc=float(sys.argv[i+1])
        sfc=0
      if sys.argv[i] == "-aiafmm":
        aiafmm[0]=float(sys.argv[i+1])
        aiafmm[1]=float(sys.argv[i+2])
        aiafmm.sort()
      if sys.argv[i] == "-sens":
        sdet=1
        if sys.argv[i+1][0]=="-":
          sens_det_image=os.path.expanduser(sens_det_image)
        else:
          sens_det_image=os.path.expanduser(sys.argv[i+1])
      if sys.argv[i] == "-nosens":
        sdet=0
      if sys.argv[i] == "-coh":
        coh=1
      if sys.argv[i] == "-scale":
        scale=float(sys.argv[i+1])
      if sys.argv[i] == "-bglin":
        bglin=int(sys.argv[i+1])
      if sys.argv[i] == "-nmon":
        nmon=int(sys.argv[i+1])
      if sys.argv[i] == "-mansel":
        mansel=float(sys.argv[i+1])
      if sys.argv[i] == "-noq":
        noq=sys.argv[i+1]
      if sys.argv[i] == "-HMIN":
        H_MIN_CHANNEL=int(sys.argv[i+1])
        if nmon < -1 or nmon > 2:
          print( "Wrong value (%d) given to flag -nmon!" %(nmon))
      if (sys.argv[i] == "-help" or sys.argv[i] == "-h"):
        usage()
      i+=1
  

  if len(roffset)==0:
    roffset.append(0)

  if len(sequence)==0:
    if len(base)==0:
      usage()
    else:
      if base[base.rfind("_")+1:].isdigit():
        sequence.append(int(base[base.rfind("_")+1:]))
        base=base[:base.rfind("_")+1]
      print( "base=%s  seq=%d" %(base,sequence[0]))

  if len(divisors) > 0:
    if len(sequence)!=len(divisors):
      print( "The number of values given for divisors (%d) differs from sequence (%d)!" %(len(divisors),len(sequence)))
      print( "I have to exit now")
      sys.exit(-1)

  if len(roffset) > 0 and len(roffset) != 1:
    if len(sequence)!=len(roffset):
      print( "The number of values given for roffset (%d) differs from the one for sequence (%d) ansd it is not 1!" %(len(roffset),len(sequence)))
      print( "I have to exit now")
      sys.exit(-1)


  fp=open("go_images_log",'a')
  tmp=""
  for i in range(1,len(sys.argv)):
    tmp+=sys.argv[i]+" "
  tmp+="\n"
  fp.write(tmp)
  fp.close()


  return base,sequence,divisors,logz,nopng,noavi,noaiaf,aiafmm,bglin,roi,roffset,fread,fpng,favi,faiaf,nmon,scale,sfc,fc,mansel,divdet,subdet,noref,maxi,mini,smaxi,smini,submaxi,submini,divmaxi,divmini,divsmaxi,divsmini,sdet,sens_det_image,sumup,gisans,window,rev,leak,noq,genx,coh,specfromaiaf,vert



def create_number_from_flipper_states(alldata,seq,number):
  global Aflipperdic,Pflipperdic,NICOS

  af=0
  pf=0

  if NICOS==1:
    if 'aflipper' in alldata[seq][number][5]:
      af=Aflipperdic[alldata[seq][number][5]['aflipper']]
    else:
      af=0
    if 'pflipper' in alldata[seq][number][5]:
      pf=Pflipperdic[alldata[seq][number][5]['pflipper']]
    else:
      pf=0
  else:
    if 'aflipper' in alldata[seq][number][5]:
      if is_number(alldata[seq][number][5]['aflipper']):
        af = int(alldata[seq][number][5]['aflipper'])
      else :
        af=0
    else:
      af=0

    if 'pflipper' in alldata[seq][number][5]:
      if is_number(alldata[seq][number][5]['pflipper']):
        pf = int(alldata[seq][number][5]['pflipper'])
      else :
        pf=0
    else:
      pf=0


#  print( "create_number_from_flipper: pf=",pf," af=",af,"  number=",pf+2*af)
  return pf+2*af


def isFloat(f):
    if not operator.isNumberType(f):
        return 0
    if f % 1:
        return 1
    else:
        return 0


def is_number(s):
  try:
    float(s)
    return True
  except ValueError:
    return False


def bg_fit_neu_bin(npo,bgo,width,niter,sigma):

  delta=(width-1)/2
  deltaf=float(delta)
  sep=1
  sigma = .5
  left_integral=0.0
  right_integral=0.0
  bg=[]
  bgn=[]

  bgn=np.zeros(len(bgo))
  bg=np.zeros(int(len(bgo)/2.0))

  for i in range(0,len(bgo),2):
    bg[i/2]=float((bgo[i]+bgo[i+1])/2.0)
  npl=int(len(bgo)/2.0)

  ## ---- Peak stripping ...                                                                                                                  
  for iter in range(niter):
    modified = 0
    for i in range(delta+1+sep,npl-1-delta-sep):
      # left side                                                                                                                             
      for j in range(i-sep,i-sep-delta,-1):
        left_integral += bg[j]
      # right side                                                                                                                            
      for j in range(i+sep,i+sep+delta):
        right_integral += bg[j]

      left_integral = float(left_integral)/deltaf

      right_integral = float(right_integral)/deltaf

      mean = float((left_integral + right_integral)/2.0)+1 #average                                                                           

      if (mean >= bg[i]):
        continue

      left_diff  = fabs(left_integral - bg[i])
      sigma_left = sigma * sqrt(left_integral)
      if (left_diff > sigma_left):
        bg[i] = mean
        modified = 1
        continue

      right_diff  = fabs(right_integral - bg[i])
      sigma_right = sigma * sqrt(right_integral)
      if (right_diff > sigma_right):
        bg[i] = mean
        modified = 1
        continue

        #      if (left_diff > sigma_left) or (right_diff > sigma_right):                                                                     
#        bg[i] = mean                                                                                                                         
#        modified = 1                                                                                                                         

    # ----  ... Reduce width for stripping during final cycles                                                                                
    if modified==0:
      break
#  print( modified,iter,niter,width                                                                                                            )
#  print( "bgfit Finished after iter=",iter                                                                                                    )

  for i in range(0,len(bg)-2,1):
    bgn[i*2]=bg[i]
    bgn[i*2+1]=float((bg[i]+bg[i+1])/2.0)
  bgn[(i+1)*2]=float(bg[i+1])


  return bgn



def bg_fit(npl,bg,width,niter,sigma):

  delta=(width-1)
  deltaf=float(delta)
  sep=1
  sigma = .5
  left_integral=0.0
  right_integral=0.0
  
  ## ---- Peak stripping ...
  for iter in range(niter):
    modified = 0
    for i in range(delta+1+sep,npl-1-delta-sep):
      # left side
      for j in range(i-sep,i-sep-delta,-1):
        left_integral += bg[j]
      # right side
      for j in range(i+sep,i+sep+delta):
        right_integral += bg[j]

      left_integral = float(left_integral)/deltaf

      right_integral = float(right_integral)/deltaf

      mean = float((left_integral + right_integral)/2.0)+1 #average

      if (mean >= bg[i]):
        continue

      left_diff  = fabs(left_integral - bg[i])
      sigma_left = sigma * sqrt(left_integral)
      if (left_diff > sigma_left):
        bg[i] = mean
        modified = 1
        continue
        
      right_diff  = fabs(right_integral - bg[i])
      sigma_right = sigma * sqrt(right_integral)
      if (right_diff > sigma_right):
        bg[i] = mean
        modified = 1
        continue

    # ----  ... Reduce width for stripping during final cycles
    if modified==0:
      break

  return bg


def bg_fit_orig(np,h,fwhm,niter,sigma):

  width=0.0
  left_integral=0.0
  left_diff=0.0
  sigma_left=0.0
  right_integral=0.0
  right_diff=0.0
  sigma_right=0.0
  mean=0.0
  
  iter=0
  i=0
  j=0
  delta=fwhm-1
  left=0
  right=0
  left_points=0
  right_points=0
  erster_punkt=0
  
  bg=h[:]

  width=fwhm
  sigma = 3.0 / sqrt(sigma)
  sigma = .5

  ## ---- Peak stripping ...
#  print( bg)

  erster_punkt=0 # continuous scan 
  np=np-1
  for iter in range(niter):
    modified = 0
    for i in range(erster_punkt,np):
      left  = max((i-width),erster_punkt)
      right = min((i+width),(np - 1))
      
      # left side
      lmax=max(left,i-width+delta)
      left_points=lmax-left
      if left_points>1 :
        left_integral = 0
#        print( "left range:", lmax, left-1,i,left_points)
        for j in range(lmax,left-1,-1):
          left_integral = left_integral + bg[j]
	left_integral = float(left_integral / left_points)
      else:
        left_integral = bg[min(i-width,0)]
      left_diff  = fabs(left_integral - bg[i])
      sigma_left = sigma * sqrt(left_integral)
  
      # right side
      rmin=min(right,i+width-delta)
      right_points = right-rmin
      if right_points>1:
        right_integral = 0
#        print( "right range:", rmin, right,i)
        for j in range(rmin,right):
          right_integral = right_integral + bg[j]
	right_integral = float(right_integral / right_points)
      else:
        right_integral = bg[rmin]
      right_diff  = fabs(right_integral - bg[i])
      sigma_right = sigma * sqrt(right_integral)
      
      mean = float(0.5*(left_integral + right_integral))+1 #average
      
      if (mean >= bg[i]):
	continue

      if (left_diff > sigma_left) or (right_diff > sigma_right):
#        print( i,bg[i],mean,left_diff,sigma_left,right_diff,sigma_right)
	bg[i] = mean
	modified = 1

    # ----  ... Reduce width for stripping during final cycles
    if modified==0:
      break

#  print( bg)
#  print( niter,iter,modified,i,mean,bg[i])
#  print( bg)
#  sys.exit(0)

  return bg



def get_file_number(name):
  return int(name[(name.rfind("_")+1):(name.find(".gz"))])



def filename_compare(x, y):
  
  file_x = get_file_number(x)
  file_y = get_file_number(y)
  if file_x > file_y:
    return 1
  elif file_x == file_y:
    return 0
  else:
    return -1



def get_file_list(base,sequence):
  
  filelist=[]
  sequence.sort()

  start=0
  
  for entry in sequence:
    for i in range(start,entry):
      filelist.append([])
    start=entry+1
#    name="%s%d.img.*.gz" %(base,entry)
    name="%s%d_*.gz" %(base,entry)
    list=glob(name)    
    list.sort(filename_compare)
    filelist.append(list)

  return filelist



def get_nicos_file_list(base,sequence):
  
  filelist=[]
  sequence.sort()

  start=0
  for entry in sequence:
    for i in range(start,entry):
      filelist.append([])
    start=entry+1

    datname=base+"%08d.dat" %(entry)
    fpdat=open(datname,'r')
    line=fpdat.readline()
    while line.find("Scan data")==-1:
      line=fpdat.readline()

    while line[0]=="#":
      line=fpdat.readline()
    
    lfilelist=[]
    while line[0]!="#":
      tmp=line.split()
      lfilelist.append(tmp[len(tmp)-1])
      line=fpdat.readline()
      if len(line)<1:
        break

    fpdat.close()
    filelist.append(lfilelist)

  return filelist




def calc_illumination(s1,s2):
  global C,L3

  beleuchtung=[]
  tmp=[]
  maxpos=0
  minpos=-10000
  slitres=0.001 #in mm

  s1l=int(round(-s1/2.0/slitres,0))
  s1r=int(round(s1/2.0/slitres,0))
  s2l=int(round(-s2/2.0/slitres,0))
  s2r=int(round(s2/2.0/slitres,0))

  for i in range(s1l,s1r+1,1):
    ai=float(i)
    for j in range(s2l,s2r+1,1):
      aj=float(j)
      pos=ai+(((aj-ai)/L3)*(L3+C))
      if pos>maxpos:
        maxpos=pos

  maxpos=int(round(maxpos,0))
#  print( "s1l=",s1l,"s2l=",s2l,"s1r=",s1r,"s2r=",s2r)
#  print( "maxpos=",maxpos)

  tmp=(0)
  for i in range(0,maxpos+1,1):
    beleuchtung.append(tmp)

#  print( "s1l=",s1l,"s2l=",s2l,"s1r=",s1r,"s2r=",s2r)
#  print( "maxpos=",maxpos)
#  print( "Res:",slitres)

  for i in range(s1l,s1r+1,1):
    ai=float(i)
    for j in range(s2l,s2r+1,1):
      aj=float(j)
      pos=ai+((aj-ai)/(L3))*(L3+C)
      if pos > -.5 :#and pos < maxpos:
        beleuchtung[int(round(pos,0))]+=1


#  for i in range(0,maxpos+1,1):
#    beleuchtung[i]=beleuchtung[i]*sin((pi/(2.0*maxpos))*(maxpos-i))**2
##    beleuchtung[i]=beleuchtung[i]*sin((pi/(2.0*maxpos))*(maxpos-i))**6

  maxint=-1
  for i in range(0,maxpos+1,1):
    if beleuchtung[i]>maxint:
      maxint=beleuchtung[i]

  integral=0
  for i in range(0,maxpos+1,1):
      integral+=beleuchtung[i]
  
#  print( "Integral=",integral,"  maxint=",maxint,"  maxpos=",maxpos)
#  print( beleuchtung)

  return beleuchtung, integral, maxint, maxpos, slitres



def calc_correction_for_illumination(beleuchtung, integral, maxint, maxpos, om, slength, res, offset,wavelength):
  global C,L3,DIST_SAMP_DET,CHANNELWIDTH

  minpos=-maxpos
  kf=1
  wert=1
  
  dom=(offset*CHANNELWIDTH)/DIST_SAMP_DET/2.0
  if om < 0.02:
    comp=1-(.02-om)*(126/wavelength)
    comp=1-(.02-om)*(-2.5*wavelength+40)
    comp=1-(.02-om)*(-1.8*wavelength+32)
  else:
    comp=1

  comp=1
  l=int(round((slength*(om))/res/2.0,0))


  if l>= 0:
    for i in range(0,maxpos+1,1):
      if abs(i) <= l:
        wert+=beleuchtung[i]

    kf= 1.0/(float(wert)/float(integral))*comp

  tmp="  maxpos=%d  res=%.3f  om=%.4f  comp=%.4f  fp/2=%.3f  wert=%d  integral=%d => kf= %.3f" %(maxpos,res,om,comp,l*res,wert,integral,kf)
  
  return kf,tmp



def calc_intensity(hsum,roi,offset,paramdic,nmon,fc):
  global S1_left, S1_right, S2_left, S2_right, NICOS
  
  intensity=0
  sigma=0

  if NICOS==1:
    time=float(paramdic['timer'][0])
    mon1=float(paramdic['mon0'][0])
    mon2=float(paramdic['mon1'][0])
  else:
    time=float(paramdic['time'])
    mon1=float(paramdic['mon1'])
    mon2=float(paramdic['mon2'])

  if check_dic(paramdic,'s1_left')!="NONE":
    s1_left=float(paramdic['s1_left'])
    s1_right=float(paramdic['s1_right'])
    s2_left=float(paramdic['s2_left'])
    s2_right=float(paramdic['s2_right'])
  else:
    print( "Slits not available!!!!")
    s1_left=1
    s1_right=1
    s2_left=1
    s2_right=1

  lroi=int(512+offset-int(roi[0]/2))
  rroi=int(512+offset+int(roi[0]/2))
    
  for i in range(lroi,rroi):
    intensity=intensity+hsum[i][4]

  for i in range(lroi,rroi):
    sigma=sigma+hsum[i][2]
  sigma=sqrt(sigma)

#  if intensity >0:
#    sigma=sqrt(intensity)
#  else:
#    sigma=0
##  if intensity > 0 and  hsum[lroi-1][3] > 0 and hsum[rroi-1][3] > 0:
##    sigma=sqrt(intensity)+sqrt((hsum[lroi-1][3]+hsum[rroi+1][3])/2)
##  else:
##    sigma=1000

  if nmon==1:
#    factor=mon1*(S1_left-s1_left+S1_right-s1_right)*(S2_left-s2_left+S2_right-s2_right)
    if fc>0:
      factor=mon1*(s1_left+s1_right)*(s2_left+s2_right)*(s1_left+s1_right)*(s2_left+s2_right)
    else:
      factor=mon1*(s1_left+s1_right)*(s2_left+s2_right)
    if factor ==0:
      factor=1.1754943508222875e-38
    intensity=intensity/factor
    sigma=sigma/factor
  elif nmon==2:
#    factor=mon2*(S2_left-s2_left+S2_right-s2_right)
    if fc>0:
      factor=mon2*(s2_left+s2_right)*(s2_left+s2_right)
    else:
      factor=mon2*(s2_left+s2_right)
    if factor ==0:
      factor=1.1754943508222875e-38
    intensity=intensity/factor
    sigma=sigma/factor
  elif nmon==0:
    intensity=intensity/time
    sigma=sigma/time
  elif nmon==-1:
    intensity=intensity
    sigma=sigma
  else:
    intensity=intensity/time/nmon
    sigma=sigma/time/nmon
    

  return intensity,sigma


  
def smooth_intensity(hsum,flag):
  global H_MIN_CHANNEL,H_MAX_CHANNEL
  
  bg=[]
  
  for i in range(H_MIN_CHANNEL,H_MAX_CHANNEL+1):
    if flag>0:
      bg.append(int(floor(float(hsum[i-1]+hsum[i]+hsum[i+1])/3.0)))
    elif flag==0:
      bg.append(hsum[i])
    elif flag==-1:
      bg.append(0)
    else:
      print( "unknown flag")
      sys.exit(-1)

  return bg



def find_bg_hsum(hsum,bglin,roi,offset,sdet):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS

  allint=[]
  
  if bglin==-1:
    if sdet == 0:
#      bg=smooth_intensity(hsum,1) # before the rotation of the detector
      bg=smooth_intensity(hsum,0) # now this
    else:
      bg=smooth_intensity(hsum,0)
    intens=bg[:]
    fwhm=25
    fwhm=35
    fwhm=55
    fwhm=roi[0]*2
    sigma=3
    niter=10000
    np=H_MAX_CHANNEL-H_MIN_CHANNEL+1
    bg=bg_fit(np,bg,fwhm,niter,sigma)
  else:
    if bglin==0:
      bg=smooth_intensity(hsum,0)
      intens=bg[:]
      for i in range(0,len(bg)):
        bg[i]=0
    if bglin >0:
      if sdet == 0:
#        bg=smooth_intensity(hsum,1) # before therotation of the detector
        bg=smooth_intensity(hsum,0) # now this
      else:
        bg=smooth_intensity(hsum,0)
      intens=bg[:]

      lroi=512+offset-int(roi[0]/2)
      rroi=512+offset+int(roi[0]/2)

      lbg=0
      rbg=0

      for i in range(0,bglin):
        lbg=lbg+bg[lroi-i-H_MIN_CHANNEL]
        rbg=rbg+bg[rroi+i-H_MIN_CHANNEL]

      lbg=float(lbg)/float(bglin)
      rbg=float(rbg)/float(bglin)
      m=float(lbg-rbg)/float(roi[0])

      for i in range(0,roi[0]):
        bg[lroi+i-H_MIN_CHANNEL]=lbg-i*m

      for i in range(0,bglin):
        bg[lroi-i-H_MIN_CHANNEL]=lbg
        bg[rroi+i-H_MIN_CHANNEL]=rbg

  for i in range(ALL_CHANNELS):
    if i <= H_MIN_CHANNEL or i > H_MAX_CHANNEL:
      allint.append([i,hsum[i],0,0,0])
    elif bglin==0:
      allint.append([i,hsum[i],intens[i-H_MIN_CHANNEL],0,intens[i-H_MIN_CHANNEL]])
    else:
#      print( i,intens[i-H_MIN_CHANNEL],intens[i-H_MIN_CHANNEL]-bg[i-H_MIN_CHANNEL] )
      allint.append([i,hsum[i],intens[i-H_MIN_CHANNEL],bg[i-H_MIN_CHANNEL],intens[i-H_MIN_CHANNEL]-bg[i-H_MIN_CHANNEL]])

  return allint



def get_params(paramfile,filename):

  i=0
#  while(paramfile[i].find("image_file")==-1):
  while(paramfile[i].find("shutter__1")==-1):
    i=i+1

  tmp=paramfile[i].split()
  paramdic={}
#  j=-1
  j=0

  for k in range(len(tmp)):
#    if tmp[k].find("s1_left")!=-1:
    if tmp[k].find("s1_left")!=-1:
      paramdic['s1_left']=j
    elif tmp[k].find("s1_right")!=-1:
      paramdic['s1_right']=j
    elif tmp[k].find("s1_top")!=-1:
      paramdic['s1_top']=j
    elif tmp[k].find("s1_bottom")!=-1:
      paramdic['s1_bottom']=j
    elif tmp[k].find("s2_left")!=-1:
      paramdic['s2_left']=j
    elif tmp[k].find("s2_right")!=-1:
      paramdic['s2_right']=j
    elif tmp[k].find("s2_top")!=-1:
      paramdic['s2_top']=j
    elif tmp[k].find("s2_bottom")!=-1:
      paramdic['s2_bottom']=j
    elif tmp[k].find("pol")!=-1:
      paramdic['pol']=j
    elif tmp[k].find("giref")!=-1:
      paramdic['giref']=j
    elif tmp[k].find("time")!=-1:
      paramdic['time']=j
#    elif tmp[k].find("Mon1")!=-1:
    elif tmp[k].find("monitor1")!=-1:
      paramdic['mon1']=j
#    elif tmp[k].find("Mon2")!=-1:
    elif tmp[k].find("monitor2")!=-1:
      paramdic['mon2']=j
#    elif tmp[k].find("Mon3")!=-1:
    elif tmp[k].find("monitor3")!=-1:
      paramdic['mon3']=j
#    elif tmp[k].find("Mon4")!=-1:
    elif tmp[k].find("monitor4")!=-1:
      paramdic['mon4']=j
    elif tmp[k].find("monitor5")!=-1:
      paramdic['mon5']=j
    elif tmp[k].find("full")!=-1:
      paramdic['full']=j
    elif tmp[k].find("roi1")!=-1:
      paramdic['roi1']=j
    elif tmp[k].find("roi2")!=-1:
      paramdic['roi2']=j
    elif tmp[k].find("roi3")!=-1:
      paramdic['roi3']=j
    elif tmp[k].find("roi4")!=-1:
      paramdic['roi4']=j
    elif tmp[k].find("roi5")!=-1:
      paramdic['roi5']=j
    elif tmp[k].find("roi6")!=-1:
      paramdic['roi6']=j
    elif tmp[k].find("roi7")!=-1:
      paramdic['roi7']=j
    elif tmp[k].find("roi8")!=-1:
      paramdic['roi8']=j
    elif tmp[k].find("roi9")!=-1:
      paramdic['roi9']=j
#    elif tmp[k].find("selector")!=-1:
    elif tmp[k].find("wavelength")!=-1:
      paramdic['selector']=j
    elif tmp[k].find("shutter")!=-1:
      paramdic['shutter']=j
    elif tmp[k].find("tx")!=-1:
      paramdic['tx']=j
    elif tmp[k].find("ty")!=-1:
      paramdic['ty']=j
    elif tmp[k].find("tz")!=-1:
      paramdic['tz']=j
    elif tmp[k].find("rx")!=-1:
      paramdic['rx']=j
    elif tmp[k].find("ry")!=-1:
      paramdic['ry']=j
    elif tmp[k].find("rz")!=-1:
      paramdic['rz']=j
    elif tmp[k].find("omega")!=-1:
      paramdic['omega']=j
    elif tmp[k].find("detarm")!=-1:
      paramdic['detarm']=j
    elif tmp[k].find("filename")!=-1:
      paramdic['image_file']=j
    elif tmp[k].find("pflipper")!=-1:
      paramdic['pflipper']=j
    elif tmp[k].find("aflipper")!=-1:
      paramdic['aflipper']=j
    elif tmp[k].find("pow1curr")!=-1:
      paramdic['pow1curr']=j
    elif tmp[k].find("pow1volt")!=-1:
      paramdic['pow1volt']=j
    elif tmp[k].find("pow3curr1")!=-1:
      paramdic['pow3curr1']=j
    elif tmp[k].find("pow3volt1")!=-1:
      paramdic['pow3volt1']=j
    elif tmp[k].find("pow3curr2")!=-1:
      paramdic['pow3curr2']=j
    elif tmp[k].find("pow3volt2")!=-1:
      paramdic['pow3volt2']=j
    elif tmp[k].find("pow4curr1")!=-1:
      paramdic['pow4curr1']=j
    elif tmp[k].find("pow4volt1")!=-1:
      paramdic['pow4volt1']=j
    elif tmp[k].find("pow4curr2")!=-1:
      paramdic['pow4curr2']=j
    elif tmp[k].find("pow4volt2")!=-1:
      paramdic['pow4volt2']=j
    elif tmp[k].find("powHP7")!=-1:
      paramdic['powHP7']=j
    elif tmp[k].find("powsampcurr1")!=-1:
      paramdic['powsampcurr1']=j
    elif tmp[k].find("powsampvolt1")!=-1:
      paramdic['powsampvolt1']=j
    elif tmp[k].find("powsampcurr2")!=-1:
      paramdic['powsampcurr2']=j
    elif tmp[k].find("powsampvolt2")!=-1:
      paramdic['powsampvolt2']=j
    elif tmp[k].find("bs1_rot")!=-1:
      paramdic['bs1_rot']=j
    elif tmp[k].find("bs1_trans")!=-1:
      paramdic['bs1_trans']=j
    elif tmp[k].find("bs2_rot")!=-1:
      paramdic['bs2_rot']=j
    elif tmp[k].find("bs2_trans")!=-1:
      paramdic['bs2_trans']=j
    elif tmp[k].find("field")!=-1:
      paramdic['field']=j
    elif tmp[k].find("magnet")!=-1:
      paramdic['magnet']=j
    elif tmp[k].find("temperature")!=-1:
      paramdic['temperature']=j
#    elif tmp[k].find("analyzer_height")!=-1:
#      paramdic['analyzer_height']=j
    elif tmp[k].find("analyzer_shift")!=-1:
      paramdic['analyzer_shift']=j
#    elif tmp[k].find("lift1")!=-1:
#      paramdic['analyzer_lift1']=j
#    elif tmp[k].find("lift2")!=-1:
#      paramdic['analyzer_lift2']=j
    elif tmp[k].find("bsd")!=-1:
      paramdic['bsd']=j
    elif tmp[k].find("Atten_A")!=-1:
      paramdic['attena']=j
    elif tmp[k].find("Atten_B")!=-1:
      paramdic['attenb']=j
    elif tmp[k].find("Atten_C")!=-1:
      paramdic['attenc']=j
    elif tmp[k].find("hvstick")!=-1:
      paramdic['hv']=j
    j=j+1
  
#  while(paramfile[i].find(filename)==-1):
#  while(paramfile[i].find("sps5__1")==-1):
#    print( "i=",i)
#    i=i+1
  i=get_file_number(filename)+1

  liste={}
  tmp=paramfile[i]
  tmp=tmp.replace(', ',',')
  tmp=tmp.replace(': ',':')
  tmp=tmp.split()

#  tmp=paramfile[i].split()

  for i in paramdic:
#    print( i,paramdic[i])
    liste[i]=tmp[paramdic[i]]

#  print( liste)

  
  return liste



def get_nicos_params(paramfile):
  d = None
  with open(paramfile, 'r') as fd:
#    d = yaml.load(fd,Loader=yaml.FullLoader)
    d = yaml.load(fd)

  devices = d["measurement"].pop("devices")
  dev = {}
  for itm in devices:
    name = itm.pop("name")
    dev[name] = itm['value']

#  print( dev)

  return dev



def check_dic(paramdic,parameter):
  global Pflipperdic, Aflipperdic, NICOS

  if NICOS==1:
    if parameter in paramdic:
      if parameter=='pflipper' or parameter=='aflipper' or parameter=='analyzer_shift' or parameter=='atten_A' or parameter=='atten_B' or parameter=='atten_C' or parameter=='compressor' or parameter=='shutter':
        tmp=paramdic[parameter]
      elif is_number(paramdic[parameter]):
        tmp=paramdic[parameter]
      else:
        tmp="NONE"
    else:
      tmp="NONE"
  else:
    if parameter in paramdic:
      if is_number(paramdic[parameter]):
        tmp=paramdic[parameter]
      else:
        tmp="NONE"
    else:
      tmp="NONE"


  return tmp


def make_gnuplot_pngfile(base,number,seq,filename,alldata,max,logz,roi,offset,bglin,mansel,maxi,mini,gisans,tfilename):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL,V_MAX_CHANNEL, ALL_CHANNELS, DIST_SAMP_DET, CHANNELWIDTH,GNUFONT,NICOS


  paramdic=alldata[seq][number][5]

  if NICOS==1:
    savetime=get_nicos_savetime(filename[:-2]+"yaml")
    rois=get_nicos_roi_det(filename)
  else:
    savetime=get_savetime(filename[:-2]+"dev")
    rois=get_roi_det(filename)

  dmax=alldata[seq][number][4][0]
  hmax=alldata[seq][number][4][1]
  vmax=alldata[seq][number][4][2]

  tmp="%s.gpl" %(filename)
  fp=open(tmp,'w')
  
  tmp="#!/usr/bin/gnuplot\n"
  fp.write(tmp)
#  tmp="set terminal pngcairo truecolor font \"Courier,36\" size 2048,2048 0xffffff\n"
  tmp="set terminal jpeg font \"%s,36\" size 2048,2048\n" %(GNUFONT)
  tmp="set terminal svg font \"%s,36\" size 2048,2048 fixed\n"%(GNUFONT)
  tmp="set terminal pngcairo truecolor noenhanced font \"%s,36\" size 2048,2048\n" %(GNUFONT)
  fp.write(tmp)
#  tmp="set output \"%s%d.%04d.jpg\"\n" %(base,seq,number)
#  tmp="set output \"%s%d.%04d.svg\"\n" %(base,seq,number)
  tmp="set output \"%s%d.%04d.png\"\n" %(base,seq,number)
  fp.write(tmp)
  tmp="set multiplot\n"
  fp.write(tmp)
  tmp="set size noratio \n"
  fp.write(tmp)
  tmp="set origin 0,0\n"
  fp.write(tmp)
  tmp="set bmargin at screen 0.1\n"
  fp.write(tmp)
  tmp="set tmargin at screen 0.8\n"
  fp.write(tmp)
  tmp="set lmargin at screen 0.13\n"
  fp.write(tmp)
  tmp="set rmargin at screen 0.83\n"
  fp.write(tmp)
#  tmp="set pm3d map interpolate 0,0 corners2color median\n"
  tmp="set pm3d map interpolate 10,10\n"
  tmp="set pm3d map interpolate 0,0 \n"
  fp.write(tmp)
  if gisans==0:
    tmp="set palette file \"-\"\n"
    fp.write(tmp)
    if logz==1:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 0\n"
      fp.write(tmp)
    else:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 .3\n"
      fp.write(tmp)
    tmp="0 1 0\n"
    fp.write(tmp)
    tmp="0 0 1\n"
    fp.write(tmp)
    tmp="1 0 0\n"
    fp.write(tmp)
    tmp="0 1 1\n"
    fp.write(tmp)
    tmp="1 0 1\n"
    fp.write(tmp)
    tmp="1 1 0\n"
    fp.write(tmp)
    tmp="1 1 1\n"
    fp.write(tmp)
    tmp="e\n"
    fp.write(tmp)
  else:
    tmp="set palette file \"-\"\n"
    fp.write(tmp)
    if logz==1:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 1\n"
      fp.write(tmp)
    else:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 1\n"
      fp.write(tmp)
    tmp="1 1 0\n"
    fp.write(tmp)
    tmp="1 0 0\n"
    fp.write(tmp)
    tmp="0 1 1\n"
    fp.write(tmp)
    tmp="1 0 1\n"
    fp.write(tmp)
    tmp="0 1 0\n"
    fp.write(tmp)
    tmp="1 1 1\n"
    fp.write(tmp)
    tmp="e\n"
    fp.write(tmp)
#
#    tmp="set palette rgb 33,13,10\n"
#    fp.write(tmp)

  if logz==1:
    tmp="set logscale cb\n"
    fp.write(tmp)
#    tmp="set logscale z\n"
#    fp.write(tmp)
    if maxi ==0:
#      tmp="set zrange [.01:%d]\n" %(dmax)
#      fp.write(tmp)
      tmp="set cbrange [.01:%d]\n" %(dmax)
      fp.write(tmp)
    else:
      if mini == 0 or mini <0:
        mini=0.01
#      tmp="set zrange [%f:%f]\n" %(mini,maxi)
#      fp.write(tmp)
      tmp="set cbrange [%f:%f]\n" %(mini,maxi)
      fp.write(tmp)
  else:
    if maxi ==0:
#      tmp="set zrange [0:%d]\n" %(dmax)
#      fp.write(tmp)
      tmp="set cbrange [0:%d]\n" %(dmax)
      fp.write(tmp)
    else:
#      tmp="set zrange [%f:%f]\n" %(mini,maxi)
#      fp.write(tmp)
      tmp="set cbrange [%f:%f]\n" %(mini,maxi)
      fp.write(tmp)

    
  for i in range(0,len(rois)):
    tmp="set arrow from %f,%f,1 to %f,%f,1 lw 1 lt -1 lc rgb \"grey\" front nohead\n" %(float(rois[i][1][0]),float(rois[i][1][1]),float(rois[i][1][2]),float(rois[i][1][1]))
    fp.write(tmp)
    tmp="set arrow from %f,%f,1 to %f,%f,1 lw 1 lt -1 lc rgb \"grey\" front nohead\n" %(float(rois[i][1][2]),float(rois[i][1][1]),float(rois[i][1][2]),float(rois[i][1][3]))
    fp.write(tmp)
    tmp="set arrow from %f,%f,1 to %f,%f,1 lw 1 lt -1 lc rgb \"grey\" front nohead\n" %(float(rois[i][1][0]),float(rois[i][1][1]),float(rois[i][1][0]),float(rois[i][1][3]))
    fp.write(tmp)
    tmp="set arrow from %f,%f,1 to %f,%f,1 lw 1 lt -1 lc rgb \"grey\" front nohead\n" %(float(rois[i][1][0]),float(rois[i][1][3]),float(rois[i][1][2]),float(rois[i][1][3]))
    fp.write(tmp)
    tmp="set label  at %d,%d,1  offset .1,.2 front \"%d \" tc rgb \"blue\" font \"%s,14\" \n" %(int(float(rois[i][1][0])),int(float(rois[i][1][1])),rois[i][0],GNUFONT)
    fp.write(tmp)

  tmp="set arrow from 0,512,1 to 1024,512,1 front nohead\n"
  fp.write(tmp)
  tmp="set arrow from 512,0,1 to 512,1024,1 front nohead\n"
  fp.write(tmp)
  #
  lr_x= int((ALL_CHANNELS/2)+offset-(float(roi[0])/2))
  lr_y= int((ALL_CHANNELS/2)-(float(roi[1])/2))
  tmp="set arrow from %d,%d,1 to %d,%d,1 lt 1 front nohead\n" %(lr_x,lr_y,lr_x+roi[0],lr_y)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lt 1 front nohead\n" %(lr_x,lr_y,lr_x,lr_y+roi[1])
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lt 1 front nohead\n" %(lr_x,lr_y+roi[1],lr_x+roi[0],lr_y+roi[1])
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lt 1 front nohead\n" %(lr_x+roi[0],lr_y,lr_x+roi[0],lr_y+roi[1])
  fp.write(tmp)
  #
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 4 lt -1 front nohead\n" %(H_MIN_CHANNEL,V_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 4 lt -1 front nohead\n" %(H_MAX_CHANNEL,V_MIN_CHANNEL,H_MAX_CHANNEL,V_MAX_CHANNEL)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 4 lt -1 front nohead\n" %(H_MIN_CHANNEL,V_MIN_CHANNEL,H_MIN_CHANNEL,V_MAX_CHANNEL)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 4 lt -1 front nohead\n" %(H_MIN_CHANNEL,V_MAX_CHANNEL,H_MAX_CHANNEL,V_MAX_CHANNEL)
  fp.write(tmp)

#  tmp="set label at 0,1250,1 front \"s1: l=%s r=%s t=%.2f b=%.2f   s2: l=%.2f r=%.2f t=%.2f b=%.2f                                        %s\" font \"Courier,18\" \n" %(check_dic(paramdic,"s1_left"),check_dic(paramdic,'s1_right'),check_dic(paramdic,'s1_top'),check_dic(paramdic,'s1_bottom'),check_dic(paramdic,'s2_left'),check_dic(paramdic,'s2_right'),check_dic(paramdic,'s2_top'),check_dic(paramdic,'s2_bottom'),filename)
  fp.write(tmp)
  if NICOS==1:
    tmp="set label at 0,1250,1 front \"s1: l=%.3f r=%.3f t=%.3f b=%.3f   s2: l=%.3f r=%.3f t=%.3f b=%.3f     %s\" font \"%s,16\" \n" %(float(paramdic["s1_left"]),float(paramdic['s1_right']),float(paramdic['s1_top']),float(paramdic['s1_bottom']),float(paramdic['s2_left']),float(paramdic['s2_right']),float(paramdic['s2_top']),float(paramdic['s2_bottom']),filename,GNUFONT)
  else:
    tmp="set label at 0,1250,1 front \"s1: l=%s r=%s t=%s b=%s   s2: l=%s r=%s t=%s b=%s                            %s\" font \"%s,18\" \n" %(check_dic(paramdic,"s1_left"),check_dic(paramdic,'s1_right'),check_dic(paramdic,'s1_top'),check_dic(paramdic,'s1_bottom'),check_dic(paramdic,'s2_left'),check_dic(paramdic,'s2_right'),check_dic(paramdic,'s2_top'),check_dic(paramdic,'s2_bottom'),filename,GNUFONT)
  fp.write(tmp)

  if NICOS==1:
    tmp="set label at 0,1220,1 front \"bs1: t:%.3f  r:%.3f   bs2: t:%.3f  r:%.3f   bsd: %.3f   He3 shift: %s   Attenuators: A:%s  B:%s  C:%s   %s\" font \"%s,16\" \n" %(float(paramdic['bs1_trans']),float(paramdic['bs1_rot']),float(paramdic['bs2_trans']),float(paramdic['bs2_rot']),float(paramdic['bsd']),check_dic(paramdic,'analyzer_shift'),check_dic(paramdic,'atten_A'),check_dic(paramdic,'atten_B'),check_dic(paramdic,'atten_C'),savetime,GNUFONT)
  else:
    tmp="set label at 0,1220,1 front \"bs1: t:%s  r:%s   bs2: t:%s  r:%s   bsd: %s   He3 shift: %s     Attenuators: A%s  B%s  C%s             %s\" font \"%s,18\" \n" %(check_dic(paramdic,'bs1_trans'),check_dic(paramdic,'bs1_rot'),check_dic(paramdic,'bs2_trans'),check_dic(paramdic,'bs2_rot'),check_dic(paramdic,'bsd'),check_dic(paramdic,'analyzer_shift'),check_dic(paramdic,'attena'),check_dic(paramdic,'attenb'),check_dic(paramdic,'attenc'),savetime,GNUFONT)
  fp.write(tmp)

  if NICOS==1:
    tmp="set label at 0,1190,1 front \"pol=%s   giref=%s   lambda=%.2f   time=%.1f   shutter=%s  Flipper: p=%s a=%s\" font \"%s,16\" \n" %(check_dic(paramdic,'pol'),check_dic(paramdic,'giref'),check_dic(paramdic,'wavelength'),float(paramdic['timer'][0]),check_dic(paramdic,'shutter'),check_dic(paramdic,'pflipper'),check_dic(paramdic,'aflipper'),GNUFONT)
    wavelength=check_dic(paramdic,'wavelength')
  else:
    
    if 'selector' in paramdic:
      if is_number(paramdic['selector']):
        selector=float(paramdic['selector'])
      else:
        selector=mansel
    else:
      selector=mansel

    tmp="set label at 0,1190,1 front \"pol=%s   giref=%s   lambda=%.2f   time=%.1f   shutter=%s  Flipper: p=%s a=%s\" font \"%s,18\" \n" %(check_dic(paramdic,'pol'),check_dic(paramdic,'giref'),selector,float(paramdic['time']),check_dic(paramdic,'shutter'),check_dic(paramdic,'pflipper'),check_dic(paramdic,'aflipper'),GNUFONT)
  fp.write(tmp)

  if NICOS==1:
    tmp="set label at 0,1160,1 front \"tx:%.2f ty:%.2f tz:%.2f rx:%.2f ry:%.2f rz:%.2f om:%.2f detarm:%.2f \" font \"%s,16\" \n" %(float(paramdic['tx']),float(paramdic['ty']),float(paramdic['tz']),float(paramdic['rx']),float(paramdic['ry']),float(paramdic['rz']),float(paramdic['omega']),float(paramdic['detarm']),GNUFONT)
  else:
    tmp="set label at 0,1160,1 front \"tx:%.2f ty:%.2f tz:%.2f rx:%.2f ry:%.2f rz:%.2f om:%.2f detarm:%.2f \" font \"%s,18\" \n" %(float(paramdic['tx']),float(paramdic['ty']),float(paramdic['tz']),float(paramdic['rx']),float(paramdic['ry']),float(paramdic['rz']),float(paramdic['omega']),float(paramdic['detarm']),GNUFONT)
  fp.write(tmp)
  if NICOS==1:
    tmp="set label at 0,1130,1 front \"detector: full: %s  roi: 1: %s   2: %s   3: %s   4: %s   5: %s   6: %s\" font \"%s,16\" \n" %(float(paramdic['detimg'][0]),float(paramdic['roi1'][0]),float(paramdic['roi2'][0]),float(paramdic['roi3'][0]),float(paramdic['roi4'][0]),float(paramdic['roi5'][0]),float(paramdic['roi6'][0]),GNUFONT)
  else:
    tmp="set label at 0,1130,1 front \"detector: full: %.7g  roi: 1: %.7g   2: %.7g   3: %.7g   4: %.7g   5: %.7g   6: %.7g\" font \"%s,18\" \n" %(float(paramdic['full']),float(paramdic['roi1']),float(paramdic['roi2']),float(paramdic['roi3']),float(paramdic['roi4']),float(paramdic['roi5']),float(paramdic['roi6']),GNUFONT)
  fp.write(tmp)
  if NICOS==1:
    tmp="set label at 0,1100,1 front \"fpga:  monitor: 1: %.7g   2: %.7g   3: %s   4: %s \" font \"%s,16\" \n" %(float(paramdic['mon0'][0]),float(paramdic['mon1'][0]),check_dic(paramdic,'mon2'[0]),check_dic(paramdic,'mon3'[0]),GNUFONT)
  else:
    tmp="set label at 0,1100,1 front \"fpga:  monitor: 1: %.7g   2: %.7g   3: %.7g   4: %.7g   5: %.7g\" font \"%s,18\" \n" %(float(paramdic['mon1']),float(paramdic['mon2']),float(paramdic['mon3']),float(paramdic['mon4']),float(paramdic['mon5']),GNUFONT)
  fp.write(tmp)

  if NICOS==1:
    tmp="set label at 0,1070,1 front \"power[C/V]: 1:%s/%s  4a:%s/%s  4b:%s/%s  sampa:%s/%s  sampb:%s/%s  powHP:%s\" font \"%s,16\" \n" %(check_dic(paramdic,'pow1curr'),check_dic(paramdic,'pow1volt'),check_dic(paramdic,'pow4curr1'),check_dic(paramdic,'pow4volt1'),check_dic(paramdic,'pow4curr2'),check_dic(paramdic,'pow4volt2'),check_dic(paramdic,'powsampcurr1'),check_dic(paramdic,'powsampvolt1'),check_dic(paramdic,'powsampcurr2'),check_dic(paramdic,'powsampvolt2'),check_dic(paramdic,'powHP7'),GNUFONT)
  else:
    tmp="set label at 0,1070,1 front \"power[C/V]: 1:%s/%s  4a:%s/%s  4b:%s/%s  sampa:%s/%s  sampb:%s/%s  powHP:%s\" font \"%s,18\" \n" %(check_dic(paramdic,'pow1curr'),check_dic(paramdic,'pow1volt'),check_dic(paramdic,'pow4curr1'),check_dic(paramdic,'pow4volt1'),check_dic(paramdic,'pow4curr2'),check_dic(paramdic,'pow4volt2'),check_dic(paramdic,'powsampcurr1'),check_dic(paramdic,'powsampvolt1'),check_dic(paramdic,'powsampcurr2'),check_dic(paramdic,'powsampvolt2'),check_dic(paramdic,'powHP7'),GNUFONT)
  fp.write(tmp)

  if NICOS==1:
    tmp="set label at 0,1040,1 front \"Magnetic Field: %s[Gauss]  Magnet: %s[A]  Temperature: %s[K]   HV: %s[V] \" font \"%s,16\" \n" %(check_dic(paramdic,'field'),check_dic(paramdic,'magnet'),check_dic(paramdic,'Ts'),check_dic(paramdic,'hv'),GNUFONT)
  else:
    tmp="set label at 0,1040,1 front \"Magnetic Field: %s[Gauss]  Magnet: %s[A]  Temperature: %s[K]   HV: %s[V] \" font \"%s,18\" \n" %(check_dic(paramdic,'field'),check_dic(paramdic,'magnet'),check_dic(paramdic,'Ts'),check_dic(paramdic,'hv'),GNUFONT)
  fp.write(tmp)

  command=""
  line=-100
  for k in range(len(sys.argv)):
    command+=sys.argv[k]+" "
    if len(command)>180:
      tmp="set label at 0,%d,1 front \"%s\" font \"%s,14\" \n"%(line,command,GNUFONT)
      fp.write(tmp)
      command=""
      line-=10
  tmp="set label at 0,%d,1 front \"%s\" font \"%s,14\" \n"%(line,command,GNUFONT)
  fp.write(tmp)
        

  tmp="set xrange [0:1024]\n"
  fp.write(tmp)
  tmp="set yrange [0:1024]\n"
  fp.write(tmp)

  if gisans==1:
    omega=float(check_dic(paramdic,'omega'))
    ttheta=float(check_dic(paramdic,'detarm'))
    min_det=ttheta-(atan((512-H_MIN_CHANNEL)*CHANNELWIDTH/DIST_SAMP_DET)*180.0/3.1415)
    max_det=ttheta+(atan((H_MAX_CHANNEL-512)*CHANNELWIDTH/DIST_SAMP_DET)*180.0/3.1415)
#  print( "angles:",omega, ttheta,atan(300.0*0.666/1900.0)*180.0/3.1415,min_det,max_det,sin(min_det/2.0/180.0*3.1415)*4*3.1415/selector,sin(max_det/2.0/180.0*3.1415)*4*3.1415/selector)

    if NICOS==1:
      tmp="set x2range [%f:%f]\n" %(sin(min_det/2.0/180.0*3.1415)*4*3.1415/wavelength,sin(max_det/2.0/180.0*3.1415)*4*3.1415/wavelength)
      fp.write(tmp)
    else:
      tmp="set x2range [%f:%f]\n" %(sin(min_det/2.0/180.0*3.1415)*4*3.1415/selector,sin(max_det/2.0/180.0*3.1415)*4*3.1415/selector)
      fp.write(tmp)


  tmp="set border lw 4 lt -1\n" 
  fp.write(tmp)
#  tmp="splot \"tmp.img\" matrix t\"%s\"\n" %(filename)
#  tmp="splot \"tmp.img\" matrix with image t\"\"\n"
  tmp="splot \"%s\" matrix with image t\"\"\n"%(tfilename)
  fp.write(tmp)
  tmp="unset arrow\n"
  fp.write(tmp)
  tmp="unset label\n"
  fp.write(tmp)
  tmp="unset pm3d\n"
  fp.write(tmp)
  ##
  ##
  ##
  ## zoom of spec region: line plots y integration ##
  tmp="set bmargin at screen 0.69\n"
  fp.write(tmp)
  tmp="set tmargin at screen %f\n" %(V_MIN_CHANNEL/1024.0*0.68+0.69)
  fp.write(tmp)
  tmp="set lmargin at screen %f\n" %(H_MIN_CHANNEL/1024.0*0.7+0.13)
  fp.write(tmp)
  tmp="set rmargin at screen %f\n" %(H_MAX_CHANNEL/1024.0*0.7+0.13)
  fp.write(tmp)
  tmp="set format \"%.1e\"\n" 
  fp.write(tmp)
  tmp="unset xtics\n" 
  fp.write(tmp)
  tmp="unset colorbox\n" 
  fp.write(tmp)
  tmp="set ytics mirror in font \"%s,14\" \n" %(GNUFONT) 
  fp.write(tmp)
  if gisans==1:
    tmp="set x2tics offset 0,-3. rotate by 90 in font \"%s,14\" tc rgbcolor \"black\"\n" %(GNUFONT)
    fp.write(tmp)
  tmp="set yrange [*:*]; set xrange [%d:%d]; set border lw 2 lt 2\n" %(512-roi[0]+int(offset),512+roi[0]+int(offset))
  fp.write(tmp)
  zmin=1E6
  zmax=0
  for i in  range(512-roi[0],512+roi[0]):
    if zmin > alldata[seq][number][2][i][3]: 
      zmin=alldata[seq][number][2][i][3]
    if zmax < alldata[seq][number][2][i][3]: 
      zmax=alldata[seq][number][2][i][3]
    if zmin > alldata[seq][number][2][i][2]: 
      zmin=alldata[seq][number][2][i][2]
    if zmax < alldata[seq][number][2][i][2]: 
      zmax=alldata[seq][number][2][i][2]
  tmp="set arrow from %d,%d to %d,%d lt 4 lw 4 front nohead\n" %(512+offset-roi[0]/2,zmin,512+offset-roi[0]/2,zmax)
  fp.write(tmp)
  tmp="set arrow from %d,%d to %d,%d lt 4 lw 4 front nohead\n" %(512+offset+roi[0]/2,zmin,512+offset+roi[0]/2,zmax)
  fp.write(tmp)
  tmp="set arrow from %d,%d to %d,%d lt 5 lw 4 front nohead\n" %(512,zmin,512,zmax)
  fp.write(tmp)
#  tmp="plot [][:*] \"-\" u 1:2 t \"\" w lp lw 2, \"\" u 1:3 t \"\" w lp lw 2, \"\" u 1:4 t \"\" w lp lw 2\n"
  tmp="plot [][:*] \"-\" u 1:2 t \"\" w lp lw 2, \"\" u 1:3 t \"\" w lp lw 2\n"
  fp.write(tmp)
  for i in  range(512-roi[0]+int(offset),512+roi[0]+int(offset)):
    tmp="%d %d %d %d\n" %(alldata[seq][number][2][i][0],alldata[seq][number][2][i][2],alldata[seq][number][2][i][3],alldata[seq][number][2][i][4])
    fp.write(tmp)
  tmp="e\n"
  fp.write(tmp)
  for i in  range(512-roi[0]+int(offset),512+roi[0]+int(offset)):
    tmp="%d %d %d %d\n" %(alldata[seq][number][2][i][0],alldata[seq][number][2][i][2],alldata[seq][number][2][i][3],alldata[seq][number][2][i][4])
    fp.write(tmp)
  tmp="e\n"
  fp.write(tmp)
#  for i in  range(512-roi[0]+offset,512+roi[0]+offset):
#    tmp="%d %d %d %d\n" %(alldata[seq][number][2][i][0],alldata[seq][number][2][i][2],alldata[seq][number][2][i][3],alldata[seq][number][2][i][4])
#    fp.write(tmp)
#  tmp="e\n"
#  fp.write(tmp)
  tmp="set xtics rotate in font \"%s,14\"\n" %(GNUFONT)
  fp.write(tmp)
  tmp="unset ytics\n"
  fp.write(tmp)
  tmp="unset arrow\n"
  fp.write(tmp)
  tmp="unset label\n"
  fp.write(tmp)
  if gisans==1:
    tmp="unset x2tics\n"
    fp.write(tmp)
 
 ##
  ## 
  ##
  ## line plots y integration ##
  tmp="set bmargin at screen 0.1\n"
  fp.write(tmp)
  tmp="set tmargin at screen %f\n" %(V_MIN_CHANNEL/1024.0*0.7+0.1)
  fp.write(tmp)
  tmp="set lmargin at screen %f\n" %(H_MIN_CHANNEL/1024.0*0.7+0.13)
  fp.write(tmp)
  tmp="set rmargin at screen %f\n" %(H_MAX_CHANNEL/1024.0*0.7+0.13)
  fp.write(tmp)
  tmp="set format \"%.1e\"\n" 
  fp.write(tmp)
  tmp="unset xtics\n" 
  fp.write(tmp)
  tmp="unset colorbox\n" 
  fp.write(tmp)
  tmp="set ytics mirror in font \"%s,14\" \n" %(GNUFONT) 
  fp.write(tmp)
  if gisans==1:
    tmp="set x2tics offset 0,-2.5 rotate by 90 out font \"%s,14\" tc rgbcolor \"black\"\n" %(GNUFONT)
    fp.write(tmp)
#  if bglin !=0:
  if bglin >-2:
    tmp="set arrow from %d,1 to %d,%d lt 4 lw 4 front nohead\n" %(512+offset-roi[0]/2,512+offset-roi[0]/2,hmax/2)
    fp.write(tmp)
    tmp="set arrow from %d,1 to %d,%d lt 4 lw 4 front nohead\n" %(512+offset+roi[0]/2,512+offset+roi[0]/2,hmax/2)
    fp.write(tmp)
  tmp="set yrange [*:*]; set xrange [*:*]; set border lw 2 lt 2\n"
  fp.write(tmp)
#  if logz==1:
#    tmp="set logscale y\n"
#    fp.write(tmp)
  tmp="plot [][:*] \"-\" u 1:2 t \"\" w l lw 2, \"\" u 1:3 t \"\" w l lw 2, \"\" u 1:4 t \"\" w l lw 2\n"
  fp.write(tmp)
  for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
    tmp="%d %d %d %d\n" %(alldata[seq][number][2][i][0],alldata[seq][number][2][i][2],alldata[seq][number][2][i][3],alldata[seq][number][2][i][4])
    fp.write(tmp)
  tmp="e\n"
  fp.write(tmp)
  for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
    tmp="%d %d %d %d\n" %(alldata[seq][number][2][i][0],alldata[seq][number][2][i][2],alldata[seq][number][2][i][3],alldata[seq][number][2][i][4])
    fp.write(tmp)
  tmp="e\n"
  fp.write(tmp)
  for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
    tmp="%d %d %d %d\n" %(alldata[seq][number][2][i][0],alldata[seq][number][2][i][2],alldata[seq][number][2][i][3],alldata[seq][number][2][i][4])
    fp.write(tmp)
  tmp="e\n"
  fp.write(tmp)
  tmp="set xtics rotate in font \"%s,14\"\n" %(GNUFONT)
  fp.write(tmp)
  tmp="unset ytics\n"
  fp.write(tmp)
  if gisans==1:
    tmp="unset x2tics\n"
    fp.write(tmp)
  ##
  ##
  ## line plot x integration ##
  if gisans==1:
    min_det=(-atan((512-V_MIN_CHANNEL)*CHANNELWIDTH/DIST_SAMP_DET))
    max_det=(atan((V_MAX_CHANNEL-512)*CHANNELWIDTH/DIST_SAMP_DET))
    if NICOS==1:
      tmp="set y2range [%f:%f]\n" %(sin(min_det/2.0)*4*3.1415/wavelength,sin(max_det/2.0)*4*3.1415/wavelength)
      fp.write(tmp)
    else:
      tmp="set y2range [%f:%f]\n" %(sin(min_det/2.0)*4*3.1415/selector,sin(max_det/2.0)*4*3.1415/selector)
      fp.write(tmp)
    tmp="set y2tics out mirror offset -13,0 font \"%s,14\" tc rgbcolor \"black\"\n" %(GNUFONT)
    fp.write(tmp)

  tmp="set bmargin at screen %f\n" %(V_MIN_CHANNEL/1024.0*0.7+0.1)
  fp.write(tmp)
  tmp="set tmargin at screen %f\n" %(V_MAX_CHANNEL/1024.0*0.7+0.1) 
  fp.write(tmp)
  tmp="set lmargin at screen %f\n" %(H_MAX_CHANNEL/1024.0*0.7+0.13)
  fp.write(tmp)
  tmp="set rmargin at screen 0.83\n"
  fp.write(tmp)
  tmp="unset logscale y\n"
  fp.write(tmp)
  tmp="unset arrow\n"
  fp.write(tmp)
  tmp="set border lw 2 lt 3\n" 
  fp.write(tmp)

  tmp="plot \"-\" u 2:1 t \"\" w l lw 2\n" 
  fp.write(tmp)
  j=0
#  print( "V_MAX_CHANNEL %d   V_MIN_CHANNEL %d" %(V_MAX_CHANNEL,V_MIN_CHANNEL))
  for i in  range(V_MIN_CHANNEL,V_MAX_CHANNEL):
    tmp="%d %d\n" %(j,alldata[seq][number][3][i])
    fp.write(tmp)
    j=j+1
  tmp="e\n"
  fp.write(tmp)
  ##
  ##
  ## line zoom plot x integration ##
  if gisans==1:
    tmp="unset xtics\n"
    fp.write(tmp)
    tmp="unset ytics\n"
    fp.write(tmp)
    tmp="unset logscale y\n"
    fp.write(tmp)
    tmp="unset arrow\n"
    fp.write(tmp)
    tmp="set bmargin at screen %f\n" %(V_MIN_CHANNEL/1024.0*0.7+0.1)
    fp.write(tmp)
    tmp="set tmargin at screen %f\n" %(V_MAX_CHANNEL/1024.0*0.7+0.1) 
    fp.write(tmp)
    tmp="set lmargin at screen %f\n" %(H_MIN_CHANNEL/1024.0*0.7+0.13)
    fp.write(tmp)
    tmp="set rmargin at screen 0.13\n"
    fp.write(tmp)
    tmp="unset logscale y\n"
    fp.write(tmp)
    tmp="unset arrow\n"
    fp.write(tmp)
    tmp="set border lw 2 lt 3\n" 
    fp.write(tmp)
    tmp="set y2tics out mirror offset 6,0 font \"%s,14\" tc rgbcolor \"black\"\n" %(GNUFONT)
    fp.write(tmp)
    centrepak=507+offset
    tmp="set yrange [*:*]\n"
    fp.write(tmp)
    tmp="set y2range [%d:%d]; set xrange [*:*]; set border lw 2 lt 2\n" %(centrepak-(roi[1]/2.0),centrepak+(roi[1]/2.0))
    fp.write(tmp)

    tmp="plot \"-\" u 2:1 t \"\" w p ps 2\n" 
    fp.write(tmp)
    j=0
    for i in  range(int(centrepak-int(roi[1]/2.0)),int(centrepak+int(roi[1]/2.0))):
      tmp="%d %d\n" %(j,alldata[seq][number][6][i])
      fp.write(tmp)
      j=j+1
    tmp="e\n"
    fp.write(tmp)







  tmp="unset multiplot\n" 
  fp.write(tmp)

  fp.close()
  
  return


def make_sumup_pngfile(logz,smaxi,smini,gisans,window,fs):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL,V_MAX_CHANNEL, ALL_CHANNELS, FS_DICT, GNUFONT

  tmp="%s_%s.gpl" %("sumup",FS_DICT[fs])
  fp=open(tmp,'w')

  tmp="#!/usr/bin/gnuplot\n"
  fp.write(tmp)
  tmp="set terminal pngcairo truecolor font \"%s,36\" size 2048,2048\n" %(GNUFONT)
  fp.write(tmp)
  tmp="set output \"%s_%s.png\"\n" %("sumup",FS_DICT[fs])
  fp.write(tmp)
  tmp="set multiplot\n"
  fp.write(tmp)
  tmp="set size noratio \n"
  fp.write(tmp)
  tmp="set origin 0,0\n"
  fp.write(tmp)
  tmp="set bmargin at screen 0.1\n"
  fp.write(tmp)
  tmp="set tmargin at screen 0.8\n"
  fp.write(tmp)
  tmp="set lmargin at screen 0.13\n"
  fp.write(tmp)
  tmp="set rmargin at screen 0.83\n"
  fp.write(tmp)
  tmp="set pm3d map interpolate 10,10\n"
  tmp="set pm3d map interpolate 0,0\n"
  fp.write(tmp)
  if gisans==0:
    tmp="set palette file \"-\"\n"
    fp.write(tmp)
    if logz==1:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 0\n"
      fp.write(tmp)
    else:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 .3\n"
      fp.write(tmp)
    tmp="0 1 0\n"
    fp.write(tmp)
    tmp="0 0 1\n"
    fp.write(tmp)
    tmp="1 0 0\n"
    fp.write(tmp)
    tmp="0 1 1\n"
    fp.write(tmp)
    tmp="1 0 1\n"
    fp.write(tmp)
    tmp="1 1 0\n"
    fp.write(tmp)
    tmp="1 1 1\n"
    fp.write(tmp)
    tmp="e\n"
    fp.write(tmp)
  else:
    tmp="set palette file \"-\"\n"
    fp.write(tmp)
    if logz==1:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 1\n"
      fp.write(tmp)
    else:
      tmp="1 1 1\n"
      fp.write(tmp)
      tmp="0 0 1\n"
      fp.write(tmp)
    tmp="1 1 0\n"
    fp.write(tmp)
    tmp="1 0 0\n"
    fp.write(tmp)
    tmp="0 1 1\n"
    fp.write(tmp)
    tmp="1 0 1\n"
    fp.write(tmp)
    tmp="0 1 0\n"
    fp.write(tmp)
    tmp="1 1 1\n"
    fp.write(tmp)
    tmp="e\n"
    fp.write(tmp)
#
#    tmp="set palette rgb 33,13,10\n"
#    fp.write(tmp)

  if logz==1:
    tmp="set logscale cb\n"
    fp.write(tmp)
    #tmp="set logscale z\n"
    #fp.write(tmp)
    if smaxi ==0:
      #tmp="set zrange [.01:*]\n"
      #fp.write(tmp)
      tmp="set cbrange [.01:*]\n"
      fp.write(tmp)
    else:
      if smini==0:
        smini=0.01
      #tmp="set zrange [%f:%f]\n" %(smini,smaxi)
      #fp.write(tmp)
      tmp="set cbrange [%f:%f]\n" %(smini,smaxi)
      fp.write(tmp)
  else:
    if smaxi ==0:
      #tmp="set zrange [0:*]\n"
      #fp.write(tmp)
      tmp="set cbrange [0:*]\n"
      fp.write(tmp)
    else:
      #tmp="set zrange [%f:%f]\n" %(smini,smaxi)
      #fp.write(tmp)
      tmp="set cbrange [%f:%f]\n" %(smini,smaxi)
      fp.write(tmp)

    
  tmp="set label at 0,1050,1 front \"Sumup File of the %s channel, counts/s\" font \"%s,36\" \n" %(FS_DICT[fs],GNUFONT)
  fp.write(tmp)
  #
  tmp="set arrow from 0,512,1 to 1024,512,1 front nohead\n"
  fp.write(tmp)
  tmp="set arrow from 512,0,1 to 512,1024,1 front nohead\n"
  fp.write(tmp)
  #
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 8 lt -1 front nohead\n" %(H_MIN_CHANNEL,V_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 8 lt -1 front nohead\n" %(H_MAX_CHANNEL,V_MIN_CHANNEL,H_MAX_CHANNEL,V_MAX_CHANNEL)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 8 lt -1 front nohead\n" %(H_MIN_CHANNEL,V_MIN_CHANNEL,H_MIN_CHANNEL,V_MAX_CHANNEL)
  fp.write(tmp)
  tmp="set arrow from %d,%d,1 to %d,%d,1 lw 8 lt -1 front nohead\n" %(H_MIN_CHANNEL,V_MAX_CHANNEL,H_MAX_CHANNEL,V_MAX_CHANNEL)
  fp.write(tmp)
  tmp="set xrange [0:1024]\n"
  fp.write(tmp)
  tmp="set yrange [0:1024]\n"
  fp.write(tmp)
  tmp="set border lw 4 lt -1\n" 
  fp.write(tmp)
  tmp="splot \"tmp_sumup_%s.img\" matrix with image t\"\"\n" %(FS_DICT[fs])
  fp.write(tmp)
  tmp="unset arrow\n"
  fp.write(tmp)
  tmp="unset label\n"
  fp.write(tmp)
  tmp="unset pm3d\n"
  fp.write(tmp)
  #
  #
  tmp="unset multiplot\n" 
  fp.write(tmp)

  fp.close()
  
  return


def analyse_pol_state(base,sequence,filelist,alldata,mr,mansel,noq):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS
  global S1_left, S1_right, S2_left, S2_right
  global Pflipperdic,Aflipperdic,Analyzerdic,NICOS
  noqmaxomrad=0
  maxomrad=0
  maxdetrad=0
  afmaxrad=0
  afminrad=10
  noqminomrad=10
  minomrad=10
  mindetrad=10
  pfs=0
  afs=0
  tpfs=0
  tafs=0
  pfs1=0
  pfs0=0
  afs1=0
  afs0=0
  ruu=0
  rdd=0
  rud=0
  rdu=0
  a_shift_0=0
  a_shift_1=0
  pol_pos_u=0
  pol_pos_d=0
  pol=3
  count =0
  maxv=0
  
  if NICOS!=1:
    polfilename=base+str(sequence[0])+".set"
    polflag=get_polarization(polfilename)
  #print( "polflag from files ->",polflag,"<- may differ from reality")

  for seq in sequence:
    for filename in filelist[seq]:
      count=count+1
      number=get_file_number(filename)
      if NICOS==1:
        pfs += int(Pflipperdic[alldata[seq][number][5]['pflipper']])
        tpfs=int(Pflipperdic[alldata[seq][number][5]['pflipper']])  
        if 'aflipper' in alldata[seq][number][5]:
          afs += int(Aflipperdic[alldata[seq][number][5]['aflipper']])
          tafs=int(Aflipperdic[alldata[seq][number][5]['aflipper']])
        else:
          afs+=0
          tafs=0
      else:
        pfs += int(alldata[seq][number][5]['pflipper'])
        tpfs=int(alldata[seq][number][5]['pflipper'])
        if 'aflipper' in alldata[seq][number][5]:
          if is_number(alldata[seq][number][5]['aflipper']):
            afs += int(alldata[seq][number][5]['aflipper'])
            tafs=int(alldata[seq][number][5]['aflipper'])
          else:
            afs+=0
            tafs=0
        else:
          afs+=0
          tafs=0
      
      if tafs==0:
        afs0+=1
      if tafs==1:
        afs1+=1
      if tpfs==0:
        pfs0+=1
      if tpfs==1:
        pfs1+=1
      if NICOS==1:
        if int(Analyzerdic[alldata[seq][number][5]['analyzer_shift']])==0:
          a_shift_0+=1
        if int(Analyzerdic[alldata[seq][number][5]['analyzer_shift']])==1 or int(Analyzerdic[alldata[seq][number][5]['analyzer_shift']])==3:
          a_shift_1+=1
      else:
        if int(alldata[seq][number][5]['analyzer_shift'])==0:
          a_shift_0+=1
        if int(alldata[seq][number][5]['analyzer_shift'])==1 or int(alldata[seq][number][5]['analyzer_shift'])==3:
          a_shift_1+=1
      if float(alldata[seq][number][5]['pol'])<500000:
        pol_pos_d+=1
      if float(alldata[seq][number][5]['pol'])>500000 and float(alldata[seq][number][5]['pol'])<800000:
        pol_pos_u+=1
      if tpfs==0 and tafs==0:
        rud+=1
      elif tpfs==0 and tafs==1:
        ruu+=1
      elif tpfs==1 and tafs==0:
        rdd+=1
      elif tpfs==1 and tafs==1:
        rdu+=1

      om  = float(alldata[seq][number][5]['omega'])
      det = float(alldata[seq][number][5]['detarm'])
      if NICOS!=1:
        if mansel==0:
          lam = float(alldata[seq][number][5]['selector'])
        else:
          lam=mansel
      else:
        lam = float(alldata[seq][number][5]['wavelength'])
        
      omrad = om/180*pi*1000  *2*pi/lam
      detrad= det/180*pi*1000 *2*pi/lam
      if omrad>maxomrad:
        maxomrad=omrad
      if omrad < minomrad:
        minomrad=omrad
      if detrad>maxdetrad:
        maxdetrad=detrad
      if detrad < mindetrad:
        mindetrad=detrad
      if maxdetrad -maxomrad + (H_MAX_CHANNEL-512)*mr*1000*2*pi/lam > afmaxrad:
        afmaxrad=maxdetrad -maxomrad + (H_MAX_CHANNEL-512)*mr*1000*2*pi/lam
      if mindetrad -minomrad - (512-H_MIN_CHANNEL)*mr*1000*2*pi/lam < afminrad:
        afminrad=mindetrad -minomrad - (512-H_MIN_CHANNEL)*mr*1000*2*pi/lam

      if len(noq)>1:
        if float(alldata[seq][number][5][noq])>noqmaxomrad:
          noqmaxomrad=float(alldata[seq][number][5][noq])
        if float(alldata[seq][number][5][noq])<noqminomrad:
          noqminomrad=float(alldata[seq][number][5][noq])

          
  if len(noq)<1:
    aimaxrad=maxomrad
    aiminrad=minomrad
  else:
    aimaxrad=noqmaxomrad
    aiminrad=noqminomrad
    aimaxrad=maxomrad
    aiminrad=minomrad
#  afmaxrad=maxdetrad -maxomrad + (H_MAX_CHANNEL-512)*mr*1000
#  afminrad=mindetrad -minomrad - (512-H_MIN_CHANNEL)*mr*1000
  pfs=float(pfs)/float(count)
  afs=float(afs)/float(count)

  print( "pfs=",pfs)
  print( "afs=",afs)
  print( "pfs0=",pfs0, pfs0/float(count))
  print( "pfs1=",pfs1, pfs1/float(count))
  print( "afs0=",afs0, afs0/float(count))
  print( "afs1=",afs1, afs1/float(count))
  print( "a_shift_0=",a_shift_0, a_shift_0/float(count))
  print( "a_shift_1=",a_shift_1, a_shift_1/float(count))
  print( "pol_pos_u=",pol_pos_u, pol_pos_u/float(count))
  print( "pol_pos_d=",pol_pos_d, pol_pos_d/float(count))
  print( "real uu=", ruu, ruu/float(count))
  print( "real dd=", rdd, rdd/float(count))
  print( "real du=", rdu, rdu/float(count))
  print( "real ud=", rud, rud/float(count))

  if pol_pos_u > 0 and pol_pos_d == 0:
    print( "Unpolarized beam: pol=3")
    pol=3
  elif pol_pos_d > 0 and pol_pos_u == 0 and a_shift_1 == 0:
    #no analyzer but polarized beam:
    if pfs0 > 0 and  pfs1 > 0:
      print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
      print( "Flipping with pflipper: pol=1!")
      pol=1
    else:
      print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
      print( "No flipping: pol=3")
      pol=3
  elif pol_pos_d > 0 and a_shift_1 > 0 and pol_pos_u == 0 and a_shift_0 == 0:
    #analyzer and polarized beam:
    if pfs0 > 0 and  pfs1 > 0 and afs0 > 0 and afs1 >0:
      if rdu>0 and rud>0:
        print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
        print( "Full poarization flip: pol=4")
        pol=4
      elif rdu>0 and rud==0:
        print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
        print( "Full poarization flip with one SF(rdu) channel: pol=7")
        pol=7
      elif rud>0 and rdu==0:
        print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
        print( "Full poarization flip with one SF(rud) channel: pol=8")
        pol=8
      elif ruu > 0 and rud  == 0  and rdd > 0 and rdu ==0 :
        print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
        print( "Non Spin flip channels: pol=6")
        pol=6
      elif rdu > 0 and ruu  == 0  and rud > 0 and rdd ==0 :
        print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
        print( "Only spin flip channles: pol=5")
        pol=5
      else:
        print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
        print( "Unknown state" )
        sys.exit(0)
#    elif (pfs0 > 0 or pfs1 > 0)  and afs0 > 0 and afs1 >0 :
    elif (rud == 0 and rdd == 0 and rdu>0 and ruu >0) or (rud > 0 and rdd > 0 and rdu == 0 and ruu == 0): 
      print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
      print( "Flipping with pflipper: pol=1")
      pol=1
#    elif (afs0 > 0 or afs1 > 0)  and pfs0 > 0 and pfs1 >0 :
    elif (rdu == 0 and rdd == 0 and rud>0 and ruu >0) or (rdu > 0 and rdd > 0 and rud == 0 and ruu == 0): 
      print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
      print( "Flipping with aflipper: pol=2")
      pol=2
    else:
      print( "shouldn't happen!1")
      print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
  else:
    print( "shouldn't happen!2")
    print( "ruu:", ruu, "  rdd:", rdd,"  rud:",rud,"  rdu:",rdu)
  
  return pol,maxomrad,aimaxrad,aiminrad,afmaxrad,afminrad



def find_maxvals_for_qmap(base,sequence,filelist,alldata,mr,nmon,divisors,mansel):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS
  global S1_left, S1_right, S2_left, S2_right

  qxmin=10
  qxmax=0
  qzmin=10
  qzmax=0
  counter=0
  maxv=0
  afminrad=1
  afmaxrad=0
  aiminrad=1
  aimaxrad=0

  mr=mr/1000.0

  for seq in sequence:
    if len(divisors)!=0:
      divider=divisors[counter]
    else:
      divider=1.0
    counter+=1
    for filename in filelist[seq]:
      number=get_file_number(filename)

      om  = float(alldata[seq][number][5]['omega'])
      det = float(alldata[seq][number][5]['detarm'])
      if NICOS==1:
        lam = float(alldata[seq][number][5]['wavelength'])
      else:
        if mansel==0:
          lam = float(alldata[seq][number][5]['selector'])
        else:
          lam = mansel
      if NICOS==1:
        time=float(alldata[seq][number][5]['timer'][0])
        mon1=float(alldata[seq][number][5]['mon0'][0])
        mon2=float(alldata[seq][number][5]['mon1'][0])
      else:
        time=float(alldata[seq][number][5]['time'])
        mon1=float(alldata[seq][number][5]['mon1'])
        mon2=float(alldata[seq][number][5]['mon2'])
      s1_left=float(alldata[seq][number][5]['s1_left'])
      s1_right=float(alldata[seq][number][5]['s1_right'])
      s2_left=float(alldata[seq][number][5]['s2_left'])
      s2_right=float(alldata[seq][number][5]['s2_right'])


      omrad=om/180*pi
      detrad=det/180*pi
      afmaxrad=detrad-omrad + (H_MAX_CHANNEL-512)*mr
      afminrad=detrad-omrad - (512-H_MIN_CHANNEL)*mr

      qx=2*pi/lam*(cos(afmaxrad)-cos(omrad))
      if qx>qxmax:
        qxmax=qx
      if qx<qxmin: 
        qxmin=qx
      qx=2*pi/lam*(cos(afminrad)-cos(omrad))
      if qx>qxmax:
        qxmax=qx
      if qx<qxmin: 
        qxmin=qx
      qz=2*pi/lam*(sin(afmaxrad)+sin(omrad))
      if qz>qzmax:
        qzmax=qz
      if qz<qzmin:
        qzmin=qz
      qz=2*pi/lam*(sin(afminrad)+sin(omrad))
      if qz>qzmax:
        qzmax=qz
      if qz<qzmin:
        qzmin=qz

      for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
        intensity=alldata[seq][number][2][i][1]/divider
        if nmon==1:
          factor=mon1*(s1_left+s1_right)*(s2_left+s2_right)
          intensity=intensity/factor
        elif nmon==2:
          factor=mon2*(s2_left+s2_right)
          intensity=intensity/factor
        elif nmon==0:
          intensity=intensity/time
        elif nmon==-1:
          intensity=intensity
        else:
          intensity=intensity/time/nmon
        if intensity > maxv:
          maxv=intensity



  return qxmin,qxmax,qzmin,qzmax,maxv




def write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm):
  global GNUFONT

  print( "afmaxrad=",afmaxrad,"  afminrad=",afminrad,"  aimaxrad=",aimaxrad,"  aiminrad=",aiminrad)
  tmp="set size ratio %f \n" %((afmaxrad-afminrad)/(aimaxrad-aiminrad))
  fp.write(tmp)
  tmp="set yrange [%f:%f] \n" %(afminrad+aiminrad,afmaxrad+aimaxrad)
  fp.write(tmp)
  tmp="set xrange [%f:%f] \n" %(aiminrad-afminrad,aimaxrad-afmaxrad)
  fp.write(tmp)

  tmp="set ytics font \"%s,14\" \n" %(GNUFONT)
  fp.write(tmp)
  tmp="set xtics font \"%s,14\" \n" %(GNUFONT)
  fp.write(tmp)
  tmp="set logscale cb\n"
  fp.write(tmp)
  tmp="set logscale z\n"
  fp.write(tmp)
  tmp="set format z \"%.2g\"\n"
  fp.write(tmp)
  tmp="set format cb \"%.2g\"\n"
  fp.write(tmp)
  tmp="set xlabel \"({/Symbol a}_i - {/Symbol a}_f) 2{/Symbol p}/{/Symbol l} [mrad/\305^-^1]\"\n"
  fp.write(tmp)
  tmp="set ylabel \"({/Symbol a}_i + {/Symbol a}_f) 2{/Symbol p}/{/Symbol l} [mrad/\305^-^1]\"\n"
  fp.write(tmp)

  return



def write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm):
  global GNUFONT

  print( "afmaxrad=",afmaxrad,"  afminrad=",afminrad,"  aimaxrad=",aimaxrad,"  aiminrad=",aiminrad)
  tmp="set size ratio %f \n" %((afmaxrad-afminrad)/(aimaxrad-aiminrad))
  fp.write(tmp)
  tmp="set yrange [%f:%f] \n" %(afminrad,afmaxrad)
  fp.write(tmp)
  tmp="set xrange [%f:%f] \n" %(aiminrad,aimaxrad)
  fp.write(tmp)

  tmp="set ytics font \"%s,14\" \n" %(GNUFONT)
  fp.write(tmp)
  tmp="set xtics font \"%s,14\" \n" %(GNUFONT)
  fp.write(tmp)
  tmp="set logscale cb\n"
  fp.write(tmp)
  tmp="set logscale z\n"
  fp.write(tmp)
  tmp="set format z \"%.2g\"\n"
  fp.write(tmp)
  tmp="set format cb \"%.2g\"\n"
  fp.write(tmp)
  tmp="set xlabel \"{/Symbol a}_i 2{/Symbol p}/{/Symbol l} [mrad/\305^-^1]\"\n"
  fp.write(tmp)
  tmp="set ylabel \"{/Symbol a}_f 2{/Symbol p}/{/Symbol l} [mrad/\305^-^1]\"\n"
  fp.write(tmp)

  return



def write_gnuoplot_header_qmap(fp,qxmin,qxmax,qzmin,qzmax,aiafmm):
  global GNUFONT

  tmp="#!/usr/bin/gnuplot\n"
  fp.write(tmp)
  tmp="set encoding iso_8859_1\n"
  fp.write(tmp)
  tmp="set terminal pngcairo enhanced truecolor font \"%s,36\" size 2048,2048\n" %(GNUFONT)
  fp.write(tmp)
  tmp="set pm3d map interpolate 0,0 \n"
  fp.write(tmp)
  tmp="#qxmax=%g qzmax=%g \n" %(qxmax,qzmax)
  fp.write(tmp)
  tmp="#set size square\n"
  fp.write(tmp)
  print( "qxmax=%f  qxmin=%f qzmax=%f  qzmin=%f" %(qxmax,qxmin,qzmax,qzmin))
  tmp="#set size ratio %f \n" %((qzmax-qzmin)/(qxmax-qxmin))
  fp.write(tmp)
  tmp="#set yrange [%f:%f] \n" %(qzmin,qzmax)
  fp.write(tmp)
  tmp="#set xrange [%f:%f] \n" %(qxmin,qxmax)
  fp.write(tmp)
  tmp="set palette file \"-\"\n"
  fp.write(tmp)
#
  tmp="1 1 1\n"
  fp.write(tmp)
  tmp="0 0 1\n"
  fp.write(tmp)
  tmp="1 1 0\n"
  fp.write(tmp)
  tmp="1 0 0\n"
  fp.write(tmp)
  tmp="0 1 1\n"
  fp.write(tmp)
  tmp="1 0 1\n"
  fp.write(tmp)
  tmp="0 1 0\n"
  fp.write(tmp)
  tmp="1 1 1\n"
  fp.write(tmp)
  tmp="e\n"
  fp.write(tmp)
#
  tmp="set ytics font \"%s,14\" \n" %(GNUFONT)
  fp.write(tmp)
  tmp="set xtics font \"%s,14\" \n" %(GNUFONT) 
  fp.write(tmp)
  tmp="set logscale cb\n"
  fp.write(tmp)
  tmp="set logscale z\n"
  fp.write(tmp)
  tmp="set format z \"%.2g\"\n"
  fp.write(tmp)
  tmp="set format cb \"%.2g\"\n"
  fp.write(tmp)
  tmp="set xlabel \"Qx [\305^-^1]\"\n"
  fp.write(tmp)
  tmp="set ylabel \"Qz [\305^-^1]\"\n"
  fp.write(tmp)
  tmp="set cbrange [*:*]\n" 
  fp.write(tmp)
  if aiafmm[0]!=0 and aiafmm[1]!=0:
    tmp="set cbrange [%g:%g]\n" %(aiafmm[0],aiafmm[1])
    fp.write(tmp)
  else:
    tmp="set cbrange [*:*]\n" 
    fp.write(tmp)

  return




def create_qmap_intlist(sequence,filelist,alldata,nmon,roffset,divisors,sfc,fc,mansel,vert,noq):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS,CHANNELWIDTH,DIST_SAMP_DET
  global S1_left, S1_right, S2_left, S2_right
  global Pflipperdic,Aflipperdic,Analyzerdic,NICOS

  maxv=-1
  maxvbg=-1
  intlist=[]
  kf=1
  factor=1

#  mr=400.0/(H_MAX_CHANNEL-H_MIN_CHANNEL)/1900.0
  mr=CHANNELWIDTH/DIST_SAMP_DET


  counter=0
  intlist=[[],[],[],[]]
  for seq in sequence:
    if len(divisors)!=0:
      divider=divisors[counter]
    else:
      divider=1.0
    if len(roffset)>1:
      offset=roffset[counter]
    else:
      offset=roffset[0]
    counter+=1
    
    for filename in filelist[seq]:
      number=get_file_number(filename)

      if vert==0:
        om  = float(alldata[seq][number][5]['omega'])
      else:
        om  = -1*(vert+float(alldata[seq][number][5]['rx']))

      det = float(alldata[seq][number][5]['detarm'])
      if NICOS==1:
        pf  = int(Pflipperdic[alldata[seq][number][5]['pflipper']])
        if 'aflipper' in alldata[seq][number][5]:
          af  = Aflipperdic[alldata[seq][number][5]['aflipper']]
        else:
          af=0
      else:
        pf  = int(alldata[seq][number][5]['pflipper'])
        if 'aflipper' in alldata[seq][number][5]:
          if is_number(alldata[seq][number][5]['aflipper']):
            af  = int(alldata[seq][number][5]['aflipper'])
          else:
            af=0
        else:
          af=0

      omrad=(om/180.0*pi)+((offset*mr)/2.0)
      detrad=det/180.0*pi

      if NICOS==1:
        time=float(alldata[seq][number][5]['timer'][0])
        mon1=float(alldata[seq][number][5]['mon0'][0])
        mon2=float(alldata[seq][number][5]['mon1'][0])
      else:
        time=float(alldata[seq][number][5]['time'])
        mon1=float(alldata[seq][number][5]['mon1'])
        mon2=float(alldata[seq][number][5]['mon2'])


      if NICOS==1:
          lam=float(alldata[seq][number][5]['wavelength'])
      else:
        if mansel==0:
          lam=float(alldata[seq][number][5]['selector'])
        else:
          lam=mansel
      q=omrad*4*pi/lam
      s1_left=float(alldata[seq][number][5]['s1_left'])
      s1_right=float(alldata[seq][number][5]['s1_right'])
      s2_left=float(alldata[seq][number][5]['s2_left'])
      s2_right=float(alldata[seq][number][5]['s2_right'])
      j=int(pf+af*2)
      detector=[]
      
#      if sfc==1:
#        kf=q
#      elif fc>0:
#        kf=1/float(alldata[seq][number][4][5])

      if sfc>0 or fc > 0:
        kf=1/float(alldata[seq][number][4][5])

      for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
        rawintensity=alldata[seq][number][2][i][1]
        intensity=alldata[seq][number][2][i][1]/divider
        intensitybg=alldata[seq][number][2][i][4]/divider
        if nmon==1:
          factor=mon1*(s1_left+s1_right)*(s2_left+s2_right)*kf
          intensity=intensity/factor
          intensitybg=intensitybg/factor
        elif nmon==2:
          factor=mon2*(s2_left+s2_right)*kf
          intensity=intensity/factor
          intensitybg=intensitybg/factor
        elif nmon==0:
          factor=kf
          intensity=intensity/time/factor
          intensitybg=intensitybg/time/factor
        elif nmon==-1:
          factor=kf
          intensity=intensity/factor
          intensitybg=intensitybg/factor
        else:
          factor=kf
          intensity=intensity/time/nmon/factor
          intensitybg=intensitybg/time/mon/factor
        if intensity > maxv:
          maxv=intensity
        if intensitybg > maxvbg:
          maxvbg=intensitybg

        afrad=detrad+(i-511)*mr-omrad

        qx=2*pi/lam*(cos(afrad)-cos(omrad))
        qz=2*pi/lam*(sin(afrad)+sin(omrad))

        if len(noq)<1:
          value_noq=0
        else:
          value_noq=float(alldata[seq][number][5][noq])

        detector.append([qx,qz,pf,af,intensity,intensitybg,i,seq,number,omrad*1000*2*pi/lam,afrad*1000*2*pi/lam,rawintensity,divider,factor,value_noq])
      intlist[j].append(detector)

  print( "maxv=",maxv,"  maxvbg=",maxvbg)

  return intlist, maxv,maxvbg



def write_qmap_data(filename,sequence,filelist,alldata,aiafmm,nmon,pfs,afs,roffset,leak,divisors,bg,flag,intlist,maxv,maxvbg,noq):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS
  global S1_left, S1_right, S2_left, S2_right
  
  fp=open(filename,"w")

  j=pfs+2*afs
  count_sec=-1
  for detector in intlist[j]:
    count_sec+=1
    count_entry=-1
    for entry in detector:
      count_entry+=1

      if filename.find("qmap")!=-1:
        x=entry[0] #qx
        y=entry[1] #qz
      elif filename.find("taiaf")!=-1:
        xp=entry[9] #alphai
        yp=entry[10] #alphaf
        x=xp-yp
        y=xp+yp
      elif filename.find("aiaf")!=-1:
        if len(noq)<1:
          x=entry[9] #noq value
        else:
          x=entry[14] #alphai
        y=entry[10] #alphaf

      pf=entry[2]
      af=entry[3]
      rawintensity=entry[11]
      divider=entry[12]
      factor=entry[13]

      if bg==0:
        intensity=entry[4]/maxv
      else:
        intensity=entry[5]/maxvbg


      if flag==0:
        if pfs==-1 and afs==-1:
          tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pfs,afs,flag,rawintensity,divider,factor)
          fp.write(tmp)
        elif pf==pfs and af==afs:
          tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pfs,afs,flag,rawintensity,divider,factor)
          fp.write(tmp)
      elif flag==1: # substract like in the case of coherent-incoherent (pf0,afx) - (pf1,afx)
        j2=((pfs+1)%2)+2*afs
        if len(intlist[j2]) > count_sec:
          if leak == 0 and pfs==0:
            if bg==0:
              intensity=intensity - (intlist[j2][count_sec][count_entry][4]/maxv)
            else:
              intensity=intensity - (intlist[j2][count_sec][count_entry][5]/maxvbg)
            tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pf,af,flag,rawintensity,divider,factor)
            fp.write(tmp)
          elif leak > 0 :
            if bg==0:
              intensity=intensity - (intlist[j2][count_sec][count_entry][4]/maxv)*leak
            else:
              intensity=intensity - (intlist[j2][count_sec][count_entry][5]/maxvbg)*leak
            tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,-2,af,flag,leak,j,j2,rawintensity,divider,factor)
            fp.write(tmp)
      elif flag==2:# substract like in the case of full spin flip measurement (pfx,af1) - (pfx,af0)
        j2=pfs+2*((afs+1)%2)
        if len(intlist[j2]) > count_sec:
          if leak == 0 and pfs==0:
            if bg==0:
              intensity=intensity - (intlist[j2][count_sec][count_entry][4]/maxv)
            else:
              intensity=intensity - (intlist[j2][count_sec][count_entry][5]/maxvbg)
            tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pf,af,flag,rawintensity,divider,factor)
            fp.write(tmp)
          elif leak > 0:
            if bg==0:
              intensity=intensity - (intlist[j2][count_sec][count_entry][4]/maxv)*leak
            else:
              intensity=intensity - (intlist[j2][count_sec][count_entry][5]/maxvbg)*leak
            tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\t%g\n" %(x,y,intensity,-2,af,flag,leak,rawintensity,divider,factor)
            fp.write(tmp)
      elif flag==4:
        if pfs==-1:
          tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pfs,afs,flag,rawintensity,divider,factor)
          fp.write(tmp)
        elif pf==pfs:
          tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pfs,afs,flag,rawintensity,divider,factor)
          fp.write(tmp)
      elif flag==5:
        if afs==-1:
          tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\t%g\n" %(x,y,intensity,pfs,afs,flag,rawintensity,divider,factor)
          fp.write(tmp)
        elif af==afs:
          tmp="%f\t%f\t%g\t%d\t%d\t%d\t%g\t%g\n" %(x,y,intensity,pfs,afs,flag,rawintensity,divider,factor)
          fp.write(tmp)

    tmp="\n"
    if pfs==-1 and afs==-1:
      fp.write(tmp)
    elif pf==pfs and af==afs:
      fp.write(tmp)
    elif pfs==-2 and pf==1:
      fp.write(tmp)


  fp.close()

  return




def plot_qx_qz_map(fp,basename,sname,sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,flag,intlist,maxv,maxvbg):


  tmp="set output \"%s.qmap.%s.png\"\n" %(basename,sname)
  fp.write(tmp)
  tmp="splot \"%s.qmap.%s.dat\" u 1:2:3 t\"%s\"\n" %(basename,sname,sname)
  fp.write(tmp)
  write_qmap_data(basename+".qmap."+sname+".dat",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,0,flag,intlist,maxv,maxvbg,"")
#
  if flag!=1:
    tmp="set output \"%s.%s.qmap.bg.png\"\n" %(basename,sname)
    fp.write(tmp)
    tmp="splot \"%s.qmap.%s.bg.dat\" u 1:2:3 t\"%s-bg\"\n" %(basename,sname,sname)
    fp.write(tmp)
    write_qmap_data(basename+".qmap."+sname+".bg.dat",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,1,flag,intlist,maxv,maxvbg,"")

  return


def plot_ai_af_map(fp,basename,sname,sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,flag,intlist,maxv,maxvbg,noq):

  tmp="set output \"%s.aiafmap.%s.png\"\n" %(basename,sname)
  fp.write(tmp)
  tmp="splot \"%s.aiafmap.%s.dat\" u 1:2:3 t\"%s\"\n" %(basename,sname,sname)
  fp.write(tmp)
  write_qmap_data(basename+".aiafmap."+sname+".dat",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,0,flag,intlist,maxv,maxvbg,noq)
#
  if flag!=1:
    tmp="set output \"%s.%s.aiafmap.bg.png\"\n" %(basename,sname)
    fp.write(tmp)
    tmp="splot \"%s.aiafmap.%s.bg.dat\" u 1:2:3 t\"%s-bg\"\n" %(basename,sname,sname)
    fp.write(tmp)
    write_qmap_data(basename+".aiafmap."+sname+".bg.dat",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,1,flag,intlist,maxv,maxvbg,noq)

  return


def plot_taiaf_map(fp,basename,sname,sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,flag,intlist,maxv,maxvbg,roi,specfromaiaf,scale):

  tmp="set output \"%s.taiafmap.%s.png\"\n" %(basename,sname)
  fp.write(tmp)
  tmp="splot \"%s.taiafmap.%s.dat\" u 1:2:3 t\"%s\"\n" %(basename,sname,sname)
  fp.write(tmp)
  write_qmap_data(basename+".taiafmap."+sname+".dat",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,0,flag,intlist,maxv,maxvbg,"")
#
  if(find(sname,"-")==-1 and specfromaiaf==1):
    print( "using:",basename+".taiafmap."+sname+".dat")
    extract_specref_from_taiaf(basename+".taiafmap."+sname+".dat",roi,roffset,scale)  
#
  if flag!=1:
    tmp="set output \"%s.%s.taiafmap.bg.png\"\n" %(basename,sname)
    fp.write(tmp)
    tmp="splot \"%s.taiafmap.%s.bg.dat\" u 1:2:3 t\"%s-bg\"\n" %(basename,sname,sname)
    fp.write(tmp)
    write_qmap_data(basename+".taiafmap."+sname+".bg.dat",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,1,flag,intlist,maxv,maxvbg,"")

  
  return



def get_aflipper(alldata,sequence,pos):

  if pos==0:
#    print( alldata[sequence[0]][0])
    paramdic=alldata[sequence[0]][0][5]
    if 'aflipper' in paramdic:
      if is_number(paramdic['aflipper']):
        af  = int(paramdic['aflipper'])
      else:
        af=0
    else:
      af=0
  elif pos==1:
    paramdic=alldata[sequence[0]][1][5]
    if 'aflipper' in paramdic:
      if is_number(paramdic['aflipper']):
        af  = int(paramdic['aflipper'])
      else:
        af=0
    else:
      af=0
  return af


def get_aflipper_nicos(alldata,sequence,filelist,pos):
  global Aflipperdic,NICOS

  number=get_file_number(filelist[sequence[0]][0])

  if pos==0:
    paramdic=alldata[sequence[0]][number][5]
    if 'aflipper' in paramdic:
      af  = Aflipperdic[paramdic['aflipper']]
    else:
      af=0
  elif pos==1:
    paramdic=alldata[sequence[0]][number+1][5]
    if 'aflipper' in paramdic:
      af  = Aflipperdic[paramdic['aflipper']]
    else:
      af=0
  return af



def make_qx_qz_map(base,sequence,filelist,alldata,gdfontpath,aiafmm,nmon,roffset,leak,divisors,roi,sfc,fc,mansel,specfromaiaf,scale,vert,noq):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS,CHANNELWIDTH,DIST_SAMP_DET
  global S1_left, S1_right, S2_left, S2_right
  global Pflipperdic,Aflipperdic,Analyzerdic,NICOS

  if NICOS==1:
    number=get_file_number(filelist[sequence[0]][0])
    paramdic=alldata[sequence[0]][number][5]
  else:
    paramdic=alldata[sequence[0]][0][5]

#  mr=400.0/(H_MAX_CHANNEL-H_MIN_CHANNEL)/1900.0*1000
  mr=CHANNELWIDTH/DIST_SAMP_DET

  pol,maxomrad,aimaxrad,aiminrad,afmaxrad,afminrad=analyse_pol_state(base,sequence,filelist,alldata,mr,mansel,noq)

  qxmin,qxmax,qzmin,qzmax,maxv=find_maxvals_for_qmap(base,sequence,filelist,alldata,mr,nmon,divisors,mansel)
      
  filename="%s%d-%d.qmap.gpl" %(base,sequence[0],sequence[len(sequence)-1])
  fp=open(filename,"w")

  write_gnuoplot_header_qmap(fp,qxmin,qxmax,qzmin,qzmax,aiafmm)

  basename="%s%d-%d" %(base,sequence[0],sequence[len(sequence)-1])

  if NICOS==1:
    if 'aflipper' in paramdic:
      af=Aflipperdic[paramdic['aflipper']]
    else:
      af=0
    pf  = Pflipperdic[paramdic['pflipper']]
  else:
    af=get_aflipper(alldata,sequence,0)
    pf  = int(paramdic['pflipper'])

  intlist,maxv,maxvbg=create_qmap_intlist(sequence,filelist,alldata,nmon,roffset,divisors,sfc,fc,mansel,vert,noq)

  if pol==1: #Flipping with pflipper
    afa=[] # im Falle das aflipper geflippt wurde aber nichtim Strahl war                                                                             
    if NICOS==1:
      afa.append(get_aflipper_nicos(alldata,sequence,filelist,0))
      afa.append(get_aflipper_nicos(alldata,sequence,filelist,1))
    else:
      afa.append(get_aflipper(alldata,sequence,0))
      afa.append(get_aflipper(alldata,sequence,1))

    plot_qx_qz_map(fp,basename,"pf0",sequence,filelist,alldata,aiafmm,nmon,0,afa[1],roffset,leak,divisors,4,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1",sequence,filelist,alldata,aiafmm,nmon,1,afa[0],roffset,leak,divisors,4,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf0-pf1",sequence,filelist,alldata,aiafmm,nmon,0,af,roffset,leak,divisors,1,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"pf0",sequence,filelist,alldata,aiafmm,nmon,0,afa[1],roffset,leak,divisors,4,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1",sequence,filelist,alldata,aiafmm,nmon,1,afa[0],roffset,leak,divisors,4,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf0-pf1",sequence,filelist,alldata,aiafmm,nmon,0,af,roffset,leak,divisors,1,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"pf0",sequence,filelist,alldata,aiafmm,nmon,0,afa[1],roffset,leak,divisors,4,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1",sequence,filelist,alldata,aiafmm,nmon,1,afa[0],roffset,leak,divisors,4,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf0-pf1",sequence,filelist,alldata,aiafmm,nmon,0,af,roffset,leak,divisors,1,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==2: #Flipping with aflipper
    plot_qx_qz_map(fp,basename,"af0",sequence,filelist,alldata,aiafmm,nmon,pf,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"af1",sequence,filelist,alldata,aiafmm,nmon,pf,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"af0-af1",sequence,filelist,alldata,aiafmm,nmon,pf,0,roffset,leak,divisors,2,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"af0",sequence,filelist,alldata,aiafmm,nmon,pf,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"af1",sequence,filelist,alldata,aiafmm,nmon,pf,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"af0-af1",sequence,filelist,alldata,aiafmm,nmon,pf,0,roffset,leak,divisors,2,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"af0",sequence,filelist,alldata,aiafmm,nmon,pf,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"af1",sequence,filelist,alldata,aiafmm,nmon,pf,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"af0-af1",sequence,filelist,alldata,aiafmm,nmon,pf,0,roffset,leak,divisors,2,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==3: #unpolarised
    plot_qx_qz_map(fp,basename,"unpol",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,0,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)
    plot_ai_af_map(fp,basename,"unpol",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)
    plot_taiaf_map(fp,basename,"unpol",sequence,filelist,alldata,aiafmm,nmon,pf,af,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==4: # Full polarization
    plot_qx_qz_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    if leak > 0:
      plot_qx_qz_map(fp,basename,"pf0af0-leak",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,1,intlist,maxv,maxvbg)
      plot_qx_qz_map(fp,basename,"pf1af1-leak",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,1,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    if leak > 0:
      plot_ai_af_map(fp,basename,"pf0af0-leak",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,1,intlist,maxv,maxvbg,noq)
      plot_ai_af_map(fp,basename,"pf1af1-leak",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,1,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    if leak > 0:
      plot_taiaf_map(fp,basename,"pf0af0-leak",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,1,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
      plot_taiaf_map(fp,basename,"pf1af1-leak",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,1,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==5: # Only spin flip channels
    plot_qx_qz_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==6: # Only non spin flip channels
    plot_qx_qz_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==7: # "Full poarization flip with one SF(ud) channel: pol=7"
    plot_qx_qz_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1af1",sequence,filelist,alldata,aiafmm,nmon,1,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

  elif pol==8: # "Full poarization flip with one SF(du) channel: pol=8"
    plot_qx_qz_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg)
    plot_qx_qz_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg)

    write_gnuoplot_header_aiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_ai_af_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)
    plot_ai_af_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,noq)

    write_gnuoplot_header_taiafmap(fp,aimaxrad,aiminrad,afmaxrad,afminrad,aiafmm)

    plot_taiaf_map(fp,basename,"pf0af0",sequence,filelist,alldata,aiafmm,nmon,0,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf1af0",sequence,filelist,alldata,aiafmm,nmon,1,0,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)
    plot_taiaf_map(fp,basename,"pf0af1",sequence,filelist,alldata,aiafmm,nmon,0,1,roffset,leak,divisors,0,intlist,maxv,maxvbg,roi,specfromaiaf,scale)

      
  fp.close()

  filename="%s%d-%d.qmap.gpl" %(base,sequence[0],sequence[len(sequence)-1])
  tmp="export GDFONTPATH=\"%s\"; gnuplot < %s" %(gdfontpath,filename)
  os.system(tmp)

  return


def find_maxvals_for_specref(base,sequence,filelist,alldata,mr,nmon,divisors,mansel,sfc,noq,vert):
  global Aflipperdic,Pflipperdic,NICOS

  maxv=-1
  maxvcoh=-1
  maxq=-1
  ointensity=0
  
  intlist=[[],[],[],[]]

  counter=0
  for seq in sequence:

    if len(divisors)!=0:
      divider=divisors[counter]
    else:
      divider=1.0

    if NICOS!=1:
      polfilename=base+str(sequence[0])+".set"
      polflag=get_polarization(polfilename)
#      print( "polflag from files ->",polflag,"<- may differ from reality")
    
    for filename in filelist[seq]:
      print( "filename:",filename)
      number=get_file_number(filename)
      paramdic=alldata[seq][number][5]

      if vert==0:
        om = float(paramdic['omega'])/180*pi
      else:
        om = (-1*(vert+float(paramdic['rx'])))/180*pi

      det = float(paramdic['detarm'])/180*pi

      if NICOS==1:
        wl= float(paramdic['wavelength'])
      else:
        if mansel==0:
          wl = float(paramdic['selector'])
        else:
          wl = mansel
        

      if NICOS==1:
        if 'pflipper' in paramdic:
          pf = Pflipperdic[paramdic['pflipper']]
        else :
          pf=0
        if 'aflipper' in paramdic:
          af = Aflipperdic[paramdic['aflipper']]
        else:
          af=0
      else:
        if is_number(paramdic['pflipper']):
          pf = int(paramdic['pflipper'])
        else :
          pf=0

        if 'aflipper' in paramdic:
          if is_number(paramdic['aflipper']):
            af = int(paramdic['aflipper'])
          else :
            af=0
        else:
          af=0

      j=pf+2*af

      q=sin(om)*4*pi/wl

 #     print( "omega=",om," q=",q," lambda=",wl," seq=",seq," number=",number)

      if q>maxq:
        maxq=q

      intensity=float(alldata[seq][number][4][3])/divider
      sigma=float(alldata[seq][number][4][4])/divider

#      if sfc>0 and float(paramdic['omega'])>0:
#        opening=float(paramdic['s2_left'])+float(paramdic['s2_right'])
#        if sfc*om<opening :
##        intensity=intensity/float(paramdic['omega'])
#          intensity=intensity/om
##        sigma=sigma/q


      if len(noq)<1:
        intlist[j].append([q,intensity,sigma,pf,af])
      else:
        intlist[j].append([float(paramdic[noq]),intensity,sigma,pf,af])

      if intensity > maxv:
        maxv=intensity
      if ointensity-0.5*intensity > maxvcoh:
        maxvcoh=ointensity-0.5*intensity
      ointensity=intensity

    counter=counter+1

  if maxv==0:
    maxv=1

  return intlist,maxv,maxvcoh,maxq



def write_gnuoplot_header_specref(fp,base,seq,lenseq,noq):

  tmp="#!/usr/bin/gnuplot\n"
  fp.write(tmp)
  tmp="set encoding iso_8859_1\n"
  fp.write(tmp)
  tmp="set terminal postscript enhanced color landscape\n"
  fp.write(tmp)
  tmp="set output \"%s%d-%d.specref.ps\"\n" %(base,seq,lenseq)
  fp.write(tmp)
  tmp="set logscale cb\n"
  fp.write(tmp)
  tmp="set logscale z\n"
  fp.write(tmp)
  if len(noq)<1:
    tmp="set xlabel \"Q [\305^-^1]\"\n"
  elif find(noq,"field")!=-1:
    tmp="set xlabel \"Magnetic field [Gauss]\"\n"
  elif find(noq,"magnet")!=-1:
    tmp="set xlabel \"Current in Magnet [A]\"\n"
  else:
    tmp="set xlabel \"%s\"\n" %(noq)

  fp.write(tmp)

  if len(noq)<1:
    tmp="set ylabel \"Refl. [normalized to 1]\"\n"
  elif find(noq,"field")!=-1:
    tmp="set ylabel \"Intensity [normalized to 1]\"\n"
  elif find(noq,"magnet")!=-1:
    tmp="set ylabel \"Intensity [normalized to 1]\"\n"
  else:
    tmp="set ylabel \"Intensity [normalized to 1]\"\n"

  fp.write(tmp)

  return



def get_jlist(intlist):

  jlist=[]

  for j in range(0,4):
    if len(intlist[j])>0:
      jlist.append(j)
  
  return jlist



def write_specref_data(basename,sname,pol,leak,intlist,maxv,maxvcoh,scale):
  global sname2, NIST_DICT
  
  jlist=get_jlist(intlist)

  maxcj=0
  maxci=0
  max=0
  cj=0
  ci=0

  if scale<0:
    cj=0
    for j in jlist:
      ci=0
      for entry in intlist[j]:
        if entry[1]/maxv>max:
          max=entry[1]/maxv
          maxci=ci
          maxcj=cj
        ci+=1
      cj+=1
    cj=0
    aver=0
    for j in jlist:
      ci=0
      for entry in intlist[j]:
        if entry[1]/maxv < max*0.9 and ci > maxci:
          break
        elif entry[1] >= max*0.9:
          cj+=1
          aver+=entry[1]/maxv
        ci+=1
    #print( "aver=",aver," aver/cj=",aver/cj)
    #print( "max=",max," cj=",cj," ci=",ci," maxcj=",maxcj," maxci",maxci)
    scale=1.0/(aver/cj)
    print( "estimated scale factor:", scale)



  if leak==0:
    filenamenist="%s%s.refl" %(basename,sname)
    fpnist=open(filenamenist,"w")
    tmpnist="#q\t int_norm\t sigma_norm\n"
    fpnist.write(tmpnist)
    for j in jlist:
#      print( jlist)
      print( "j=",j," pol=",pol)
      if pol!=3:
        filename="%s%s%s.dat" %(basename,sname,sname2[j])
      else:
        filename="%s%s.dat" %(basename,sname)
      
      fp=open(filename,"w")
      tmp="#q\t int_norm\t sigma_norm\t pf\t af\t int_raw\t sigma_raw\t\n"
      fp.write(tmp)
      tmpnist=NIST_DICT[pol-1][j]+"\n"
      fpnist.write(tmpnist)
      for entry in intlist[j]:
        tmp="%f\t%g\t%g\t%d\t%d\t%f\t%f\n" %(entry[0],entry[1]*scale/maxv,entry[2]*scale/maxv,entry[3],entry[4],entry[1],entry[2])
        fp.write(tmp)
        tmpnist="%f\t%g\t%g\n" %(entry[0],entry[1]*scale/maxv,entry[2]*scale/maxv)
        fpnist.write(tmpnist)
      if j!=jlist[len(jlist)-1]:
        tmpnist="\n\n"
        fpnist.write(tmpnist)
      fp.close()

    fpnist.close()

    if pol!=3:
      tmp="paste"
      for j in jlist:
        filename="%s%s%s.dat" %(basename,sname,sname2[j])
        tmp+=" %s" %(filename)
      filename="%s%sall.dat" %(basename,sname)
      tmp+=" > %s" %(filename)
      os.system(tmp)

  else:
    print( "leak on reflectivity is not implemented")
    for j in jlist:
      if pol!=3:
        filename="%s%s%s.dat" %(basename,sname,sname2[j])
      else:
        filename="%s%s.dat" %(basename,sname)
      fp=open(filename,"w")
      tmp="#q\t int_norm\t sigma_norm\t pf\t af\t int_raw\t sigma_raw\t\n"
      fp.write(tmp)
      for entry in intlist[j]:
        tmp="%f\t%g\t%g\t%d\t%d\t%f\t%f\n" %(entry[0],entry[1]*scale/maxv,entry[2]*scale/maxv,entry[3],entry[4],entry[1],entry[2])
        fp.write(tmp)
      fp.close()

    if pol!=3:
      tmp="paste"
      for j in jlist:
        filename="%s%s%s.dat" %(basename,sname,sname2[j])
        tmp+=" %s" %(filename)
      filename="%s%sall.dat" %(basename,sname)
      tmp+=" > %s" %(filename)
      os.system(tmp)


  return



def create_plot_string(basename, sname, pol, jlist, raw,curveflag,maxq):
  global sname2
  
  tmp="plot "
  
  sep_flag=""
  for j in jlist:
    if pol!=3:
      filename="%s%s%s.dat" %(basename,sname,sname2[j])
    else:
#      filename="%s%s.dat" %(basename,sname)
      sname2u={0:"unpol", 1:"unpol",2:"unpol",3:"unpol"}
      filename="%s%s.dat" %(basename,sname2u[j])
    if raw==0:
      tmp+=sep_flag+"\"%s\" u 1:2:3 t \"%s\" with errorlines" %(filename,sname2[j])
    else:
      tmp+=sep_flag+"\"%s\" u 1:6:7 t \"%s\" with errorlines" %(filename,sname2[j])
    sep_flag=","

  tmp=tmp[:-1]

  if curveflag==1:
    if os.path.isfile("curve.dat")==True:
      tmp+=", \"curve.dat\" u 1:2 t \"curve.dat\" w lp\n"
      tmp2=" [:%f] " %(maxq)
      tmp=tmp[0:4]+tmp2+tmp[5:]
    else:
      tmp=""
  tmp+="\n"

  return tmp


  
def plot_specref(fp,basename,sname,pol,leak,intlist,maxv,maxvcoh,maxq):
  global sname2

  jlist=get_jlist(intlist)
  
#  print( "jlist", jlist)
  
  tmp=create_plot_string(basename,sname,pol, jlist,0,0,maxq)
  fp.write(tmp)    

  tmp=create_plot_string(basename,sname,pol, jlist,1,0,maxq)
  fp.write(tmp)    

  tmp="set log y\n"
  fp.write(tmp)

  tmp=create_plot_string(basename,sname,pol, jlist,0,0,maxq)
  fp.write(tmp)    

  tmp=create_plot_string(basename,sname,pol, jlist,1,0,maxq)
  fp.write(tmp)    

  tmp=create_plot_string(basename,sname,pol, jlist,0,1,maxq)
  fp.write(tmp)    

  tmp="unset log y\n"
  fp.write(tmp)

  filename="%s%sall.dat" %(basename,sname)
  if pol==1 or pol==2 or pol==5 or pol==6:
    tmp="plot \"%s\" u 1:(($2-$9)/($2+$9)) t \"asymmetry ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[0]],sname2[jlist[1]])
    fp.write(tmp)
  if pol==4:
    tmp="plot \"%s\" u 1:(($2-$9)/($2+$9)) t \"asymmetry ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[0]],sname2[jlist[1]])
    fp.write(tmp)
    tmp="plot \"%s\" u 1:(($16-$23)/($16+$23)) t \"asymmetry ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[2]],sname2[jlist[3]])
    fp.write(tmp)
  if pol==7:
    tmp="plot \"%s\" u 1:(($16-$9)/($16+$9)) t \"asymmetry ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[1]],sname2[jlist[2]])
    fp.write(tmp)
  if pol==8:
    tmp="plot \"%s\" u 1:(($2-$9)/($2+$9)) t \"asymmetry ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[0]],sname2[jlist[1]])
    fp.write(tmp)

  if pol==1 or pol==2 or pol==5 or pol==6:
    tmp="plot \"%s\" u 1:(($2/$9)) t \"flipping ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[0]],sname2[jlist[1]])
    fp.write(tmp)
    tmp="plot \"%s\" u 1:(($9/$2)) t \"flipping ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[1]],sname2[jlist[0]])
    fp.write(tmp)
  if pol==4:
    tmp="plot \"%s\" u 1:(($2/$9)) t \"flipping ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[0]],sname2[jlist[1]])
    fp.write(tmp)
    tmp="plot \"%s\" u 1:(($16/$23)) t \"flipping ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[2]],sname2[jlist[3]])
    fp.write(tmp)
  if pol==7:
    tmp="plot \"%s\" u 1:(($16/$9)) t \"flipping ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[2]],sname2[jlist[1]])
    fp.write(tmp)
  if pol==8:
    tmp="plot \"%s\" u 1:(($2/$9)) t \"flipping ratio: %s / %s\" with lp \n" %(filename,sname2[jlist[0]],sname2[jlist[1]])
    fp.write(tmp)


  return





def make_specref_file(base,sequence,filelist,alldata,divisors,nmon,scale,sfc,fc,mansel,noq,genx,leak,vert):
  global CHANNELWIDTH, DIST_SAMP_DET

  sname={1:"pf", 2:"af",3:"unpol",4:"pfaf",5:"SF",6:"NSF",7:"2*NSF+du",8:"2*NSF+ud" }

#  mr=400.0/(H_MAX_CHANNEL-H_MIN_CHANNEL)/1900.0*1000
  mr=(CHANNELWIDTH/DIST_SAMP_DET)*1000.0

  print( "Inside make_spec_ref")
  pol,maxomrad,aimaxrad,aiminrad,afmaxrad,afminrad=analyse_pol_state(base,sequence,filelist,alldata,mr,mansel,noq)
  print( "analyse_pol: done")
  intlist,maxv,maxvcoh,maxq=find_maxvals_for_specref(base,sequence,filelist,alldata,mr,nmon,divisors,mansel,sfc,noq,vert)
  print( "find_maxvals_for_specref: done")

  #illumination correction is missing here ???

  filename="%s%d-%d.specref.gpl" %(base,sequence[0],sequence[len(sequence)-1])
  fp=open(filename,"w")
  print( "filename=",filename, "fpopen: done")

  write_gnuoplot_header_specref(fp,base,sequence[0],sequence[len(sequence)-1],noq)
  print( "write_gnuplot_header_specref: done=")

  basename="%s%d-%d.specref." %(base,sequence[0],sequence[len(sequence)-1])
  print( "basename=",basename,"  done")
  
  write_specref_data(basename,sname[pol],pol,leak,intlist,maxv,maxvcoh,scale)
  print( "write_specref_dat: done")
  
  plot_specref(fp,basename,sname[pol],pol,leak,intlist,maxv,maxvcoh,maxq)
  print( "plot_specref: done")

  
  fp.close()

  tmp="gnuplot %s" %(filename)
  result=os.system(tmp)
  print( "result:",result)

  return



def append_analyse_data(base,sequence,filelist,bglin,roi,roffset,alldata,maxv,nmon,sfc,fc,sdet,sens_det_array,mansel):
  global C,L3,RES,DIST_SAMP_DET,CHANNELWIDTH
  hsum=[]
  vsum=[]
  vroisum=[]
  hmaxv=0
  vmaxv=0
  maxint=0
  hmaxint=0
  vmaxint=0
  counter=0
  kf=1


#  print( "len(sequence)",len(sequence) )
#  print( "len(alldata)",len(alldata) )

  start=sequence[0]
  
  for seq in sequence:
    if len(sequence) > len(alldata):
      for i in range(start,seq+1):
        np.append(alldata, [])
      #      alldata.np.append([])
      
    start=seq
    if NICOS!=1:
#      print( base,seq,sequence,filelist)
      filename="%s%s.dat" %(base,seq)
      fp=open(filename)
      paramfile=fp.readlines()
      fp.close()
    
    if len(roffset)>1:
      offset=roffset[counter]
    else:
      offset=roffset[0]

    if fc>0 or sfc >0:
      if NICOS==1:
        paramfile="%s.yaml" %(filelist[seq][0][:-3])
        paramdic=get_nicos_params(paramfile)
      else:
        paramdic=get_params(paramfile,filelist[seq][0])
      s1_left=float(paramdic['s1_left'])
      s1_right=float(paramdic['s1_right'])
      s2_left=float(paramdic['s2_left'])
      s2_right=float(paramdic['s2_right'])
 
      if fc>0:
        beleuchtung, integral, maxint, maxpos, res=calc_illumination(s1_left+s1_right,s2_left+s2_right)


#    print( "seq",seq)
#    print( "len(filelist[seq])",len(filelist[seq]) )
#    print( "len(alldata[seq])",len(alldata[seq]))

    if len(filelist[seq])>len(alldata[seq]):
      c2=0
      for filename in filelist[seq]:
        if c2>len(alldata[seq]):
          if NICOS==1:
            paramfile="%s.yaml" %(filename[:-3])
            paramdic=get_nicos_params(paramfile)
          else:
            paramdic=get_params(paramfile,filename)
          
          number=get_file_number(filename)
        #      alldata[seq].append([number])
          np.append(alldata[seq],[number])
          
          if vert!=0:
            loffset=offset-int((tan(2*radians(vert+float(paramdic['rx'])))*DIST_SAMP_DET/CHANNELWIDTH))
          else:
            loffset=offset

          hsum,hmaxv,vsum,vmaxv,maxv,vroisum=integrate(filename,roi,sdet,sens_det_array,loffset,vert)
          
          hsum=find_bg_hsum(hsum,bglin,roi,loffset,sdet)
          
          intensity,sigma=calc_intensity(hsum,roi,loffset,paramdic,nmon,fc)
          
          if fc>0:
            if vert==0:
              om=float(paramdic['omega'])*pi/180.0
            else:
              om=radians(-1*(vert+float(paramdic['rx'])))
              
            if om>0.01*pi/180.0:
              if mansel==0:
                wavelength=float(paramdic['selector'])
              else:
                wavelength=mansel

              if fc>0:
                kf,tmp=calc_correction_for_illumination(beleuchtung, integral, maxint, maxpos, om, fc, res, loffset,wavelength)
              elif sfc>0:
                tmp=""
                if om*sfc<float(s2_left+s2_right):
                  kf=1/(om/(float(s2_left+s2_right)/sfc))
                else:
                  kf=1
                tmp="simple footprint( correction with sample length=%d and factor %.3f" %(round(sfc,0),kf) 
              else:
                kf=1
                tmp="No footprint( correction"
              intensity=kf*intensity
              sigma=sigma
              print( "filename2=",filename,tmp)
            else:
              print( "reading=",filename)

          np.append(alldata[seq][number],[])
          np.append(alldata[seq][number],hsum[:])
          np.append(alldata[seq][number],vsum[:])
          np.append(alldata[seq][number],[maxv,hmaxv,vmaxv,intensity,sigma,kf])
          np.append(alldata[seq][number],paramdic)
          np.append(alldata[seq][number],vroisum[:])
        else:
          c2+=1

    counter+=1
  print( "Found maximal value:",maxv)

  return alldata,maxv


def integrate(filename,roi,sdet,sens_det_array,offset,vert):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS,NICOS

  hsum=[]
  hroisum=[]
  vsum=[]
  vroisum=[]
  tmp=[]
  hmaxv=1
  vmaxv=1
  maxv=1
  det=[]
  detr=[]
  array=[]

  h_min_channel= int((ALL_CHANNELS/2) - (float(roi[1])/2))
  h_max_channel= int((ALL_CHANNELS/2) + (float(roi[1])/2))
  v_min_channel= int((ALL_CHANNELS/2) + offset - (float(roi[0])/2))
  v_max_channel= int((ALL_CHANNELS/2) + offset + (float(roi[0])/2))

  fp=gzip.open(filename,"r")

  #reading in the detector
  for string in fp.readlines():
    det.append(string)

  hsum=np.zeros(ALL_CHANNELS)
  hroisum=np.zeros(ALL_CHANNELS)
  vsum=np.zeros(ALL_CHANNELS)
  vroisum=np.zeros(ALL_CHANNELS)

  for i in range(0,ALL_CHANNELS):
    #flipping the detector picture up down
    if NICOS!=1:
      detr.append(det[ALL_CHANNELS-i-1])
    else:
      detr.append(det[i])

####
  for line in detr:
    tmp=line.split()
    array.append(tmp)

  if vert!=0:
    tmp = list(reversed(zip(*array)))
  else:
    tmp=array
####



  for zeile in range(h_min_channel+1,h_max_channel):
    b=tmp[zeile]
    c=0
    croi=0
    for spalte in range(H_MIN_CHANNEL+1,H_MAX_CHANNEL):
      sens=float(sens_det_array[zeile][spalte])
      if sens>0 and sdet==1:
        value=float(b[spalte])/sens
      else:
        value=float(b[spalte])

      if value > 0:
        c=c+value
        hsum[spalte]=hsum[spalte]+value
        #          if zeile>h_min_channel and zeile<h_max_channel:
        hroisum[spalte]=hsum[spalte]
        if spalte>v_min_channel and spalte<v_max_channel:
          croi=croi+value
        if value > maxv:
          if sens>0 and sdet==1 :
            maxv=float(b[spalte])/sens
          else:
            maxv=float(b[spalte])
    vsum[zeile]=c
    vroisum[zeile]=croi
    if c > vmaxv:
      vmaxv=c
    
  for i in range(0,ALL_CHANNELS):
    if hsum[i] > hmaxv:
      hmaxv=hsum[i]

  fp.close()
#  print( "hsum\n",hsum)
#  print( "vsum\n",vsum)
#  print( hmaxv, vmaxv, maxv)
  return hsum,hmaxv,vsum,vmaxv,maxv,vroisum


def integrate_orig(filename,roi,sdet,sens_det_array,offset,vert):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS,NICOS

  hsum=[]
  hroisum=[]
  vsum=[]
  vroisum=[]
  tmp=[]
  hmaxv=1
  vmaxv=1
  maxv=1
  det=[]
  detr=[]
  array=[]

#  h_min_channel= int((ALL_CHANNELS/2) + offset - (float(roi[1])/2))
#  h_max_channel= int((ALL_CHANNELS/2) + offset + (float(roi[1])/2))
#  v_min_channel= int((ALL_CHANNELS/2) - (float(roi[0])/2))
#  v_max_channel= int((ALL_CHANNELS/2) + (float(roi[0])/2))
  h_min_channel= int((ALL_CHANNELS/2) - (float(roi[1])/2))
  h_max_channel= int((ALL_CHANNELS/2) + (float(roi[1])/2))
  v_min_channel= int((ALL_CHANNELS/2) + offset - (float(roi[0])/2))
  v_max_channel= int((ALL_CHANNELS/2) + offset + (float(roi[0])/2))

  #print( "h_min_channel=",h_min_channel," h_max_channel=",h_max_channel)
  #print( "v_min_channel=",v_min_channel," v_max_channel=",v_max_channel)

  fp=gzip.open(filename,"r")

  #reading in the detector
  for string in fp.readlines():
    det.append(string)

  for i in range(0,ALL_CHANNELS):
    hsum.append(0)
    hroisum.append(0)
    vsum.append(0)
    vroisum.append(0)
    #flipping the detector picture up down
    if NICOS!=1:
      detr.append(det[ALL_CHANNELS-i-1])
    else:
      detr.append(det[i])

####
  for line in detr:
    tmp=line.split()
    array.append(tmp)

  if vert!=0:
    tmp = list(reversed(zip(*array)))
  else:
    tmp=array
####



  zeile=0
#  for string in detr:
#   b=string.split()
  for i in range(0,ALL_CHANNELS):
    b=tmp[i]
    c=0
    croi=0
    if zeile>h_min_channel and zeile<h_max_channel:
      for spalte in range(len(b)):
        if spalte>H_MIN_CHANNEL and spalte<H_MAX_CHANNEL:
          if float(sens_det_array[zeile][spalte])>0 and sdet==1:
            value=float(b[spalte])/float(sens_det_array[zeile][spalte])
          else:
#            value=int(b[spalte])
            value=float(b[spalte])

          if value > 0:
            c=c+value
            hsum[spalte]=hsum[spalte]+value
            if zeile>h_min_channel and zeile<h_max_channel:
              hroisum[spalte]=hsum[spalte]
              if spalte>v_min_channel and spalte<v_max_channel:
                croi=croi+value
            if value > maxv:
              if float(sens_det_array[zeile][spalte])>0 and sdet==1 :
                maxv=float(b[spalte])/float(sens_det_array[zeile][spalte])
              else:
#                maxv=int(b[spalte])
                maxv=float(b[spalte])
      vsum[zeile]=c
      vroisum[zeile]=croi
      if c > vmaxv:
        vmaxv=c
    zeile=zeile+1
    
  for i in range(0,ALL_CHANNELS):
    if hsum[i] > hmaxv:
      hmaxv=hsum[i]

  fp.close()
#  print( "hsum\n",hsum)
#  print( "vsum\n",vsum)
#  print( hmaxv, vmaxv, maxv)
  return hsum,hmaxv,vsum,vmaxv,maxv,vroisum


def analyse_data(base,sequence,filelist,bglin,roi,roffset,nmon,sfc,fc,sdet,sens_det_array,mansel,vert):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS,C,L3,RES,DIST_SAMP_DET,CHANNELWIDTH

  alldata=[]
  hsum=[]
  vsum=[]
  vroisum=[]
  maxv=0
  hmaxv=0
  vmaxv=0
  maxint=0
  hmaxint=0
  vmaxint=0
  counter=0
  kf=1

  start=0
  
  for seq in sequence:
    for i in range(start,seq+1):
      alldata.append([])
      
    start=seq
    if NICOS!=1:
      filename="%s%s.dat" %(base,seq)
      print( "filename=",filename)
      fp=open(filename)
      paramfile=fp.readlines()
      fp.close()

    if len(roffset)>1:
      offset=roffset[counter]
    else:
      offset=roffset[0]
    
    if fc>0 or sfc>0:
      if NICOS==1:
        paramfile="%s.yaml" %(filelist[seq][0][:-3])
        paramdic=get_nicos_params(paramfile)
      else:
        paramdic=get_params(paramfile,filelist[seq][0])
      s1_left=float(paramdic['s1_left'])
      s1_right=float(paramdic['s1_right'])
      s2_left=float(paramdic['s2_left'])
      s2_right=float(paramdic['s2_right'])
      if fc>0:
        beleuchtung, integral, maxint, maxpos, res=calc_illumination(s1_left+s1_right,s2_left+s2_right)

    if NICOS==1:
      if len(filelist[seq])>0:
        number=get_file_number(filelist[seq][0])
        for i in range(0,number):
          alldata[seq].append([i])

    for filename in filelist[seq]:
      if NICOS==1:
        paramfile="%s.yaml" %(filename[:-3])
        paramdic=get_nicos_params(paramfile)
      else:
        paramdic=get_params(paramfile,filename)

      number=get_file_number(filename)
      alldata[seq].append([number])

      if vert!=0:
        loffset=offset-int((2*tan(radians(vert+float(paramdic['rx'])))*DIST_SAMP_DET/CHANNELWIDTH))
#        print( float(paramdic['rx']),float(paramdic['rx'])+vert,2*radians(vert+float(paramdic['rx'])),tan(2*radians(vert+float(paramdic['rx']))),tan(2*radians(vert+float(paramdic['rx'])))*DIST_SAMP_DET, tan(2*radians(vert+float(paramdic['rx'])))*DIST_SAMP_DET/CHANNELWIDTH)
      else:
        loffset=offset

      hsum,hmaxv,vsum,vmaxv,maxv,vroisum=integrate(filename,roi,sdet,sens_det_array,loffset,vert)

      hsum=find_bg_hsum(hsum,bglin,roi,loffset,sdet)

      intensity,sigma=calc_intensity(hsum,roi,loffset,paramdic,nmon,fc)
      if fc>0 or sfc >0:
        if vert==0:
          om=float(paramdic['omega'])*pi/180.0
        else:
          om=radians(-1*(vert+float(paramdic['rx'])))

        if om>0.01*pi/180.0:
          if NICOS==1:
            wavelength=float(paramdic['wavelength'])
          else:
            if mansel==0:
              wavelength=float(paramdic['selector'])
            else:
              wavelength=mansel
          if fc>0:
            kf,tmp=calc_correction_for_illumination(beleuchtung, integral, maxint, maxpos, om, fc, res,loffset,wavelength)
          elif sfc>0:
            if om*sfc<float(s2_left+s2_right):
              kf=1/(om/(float(s2_left+s2_right)/sfc))
            else:
              kf=1
            tmp="simple footprint( correction with sample length=%d and factor %.3f" %(round(sfc,0),kf) 
          else:
            kf=1
            tmp="No footprint( correction"
          intensity=kf*intensity
          sigma=sigma
          print( "filename2=",filename,tmp)
        else:
          print( "filename2=",filename)

      alldata[seq][number].append([])
      alldata[seq][number].append(hsum[:])
      alldata[seq][number].append(vsum[:])
      alldata[seq][number].append([maxv,hmaxv,vmaxv,intensity,sigma,kf])
      alldata[seq][number].append(paramdic)
      alldata[seq][number].append(vroisum[:])
      
    counter+=1
  print( "Found maximal value: maxv %d  hmaxv %d   vmaxv %d" %(maxv,hmaxv,vmaxv))

  return alldata,maxv


def create_tmp_img(tfilename,sens_det_array,sumup,window,fs,time,vert):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL,V_MAX_CHANNEL, ALL_CHANNELS, SUMUP_ARRAY, NICOS

  det=[]
  detr=[]

  tmpfile=tfilename+".gz"
  fp=gzip.open(tmpfile,'r')
  array=[]

  for line in fp:
    det.append(line)

  for i in range(0,ALL_CHANNELS):
    if NICOS!=1:
      detr.append(det[ALL_CHANNELS-i-1])
    else:
      detr.append(det[i])

  for line in detr:
    tmp=line.split()
    array.append(tmp)

  fp.close()

  if vert!=0:
#    tmp = list(zip(*array[::-1]))
    tmp = list(reversed(zip(*array)))


#  counter1=0
  if vert==0:
    for i in range(H_MIN_CHANNEL,H_MAX_CHANNEL):
      for j in range(V_MIN_CHANNEL,V_MAX_CHANNEL):
        if sens_det_array[i][j]>0:
          array[i][j]=str(float(array[i][j])/(sens_det_array[i][j]))
  else:
    for i in range(H_MIN_CHANNEL,H_MAX_CHANNEL):
      for j in range(V_MIN_CHANNEL,V_MAX_CHANNEL):
        if sens_det_array[i][j]>0:
          array[i][j]=str(float(tmp[i][j])/(sens_det_array[i][j]))
  if sumup==1:
    SUMUP_ARRAY[fs][1024][1024]+=1
    for i in range(0,ALL_CHANNELS):
      for j in range(0,ALL_CHANNELS):
        SUMUP_ARRAY[fs][i][j]+=float(array[i][j])/time

  fpt=open(tfilename,'w')
  tmp=""
#  for i in range(ALL_CHANNELS,0,-1):
  for i in range(0,ALL_CHANNELS):
    tmp+=join(array[i-1])
    tmp+="\n"
  fpt.write(tmp)
  fpt.close()

#  print( "pixels above 0: before=%d    after=%d" %(counter1,counter2))

  if sumup==1:
    tmp="tmp_sumup_%s.img" %(FS_DICT[fs])
    fps=open(tmp,'w')
    tmps=""
#    for i in range(ALL_CHANNELS,0,-1):
    for i in range(0,ALL_CHANNELS):
      for j in range(0,ALL_CHANNELS):
        if SUMUP_ARRAY[fs][i-1][j]==0:          
#          tmps+="NAN  "
          tmps+="0  "
        else:
          tmps+=str(SUMUP_ARRAY[fs][i-1][j])+" "
      tmps+="\n"
    fps.write(tmps)
    fps.close()

  
  return
  


def cut_out_window(window,fs):
  global SUMUP_ARRAY, FS_DICT
  counter =0
  while counter < len(window):
    x0=int(window[counter])
    counter+=1
    y0=int(window[counter])
    counter+=1
    x1=int(window[counter])
    counter+=1
    y1=int(window[counter])
    counter+=1
    data=[]
    for x in range(x0,x1+1):
      tmp=[]
      tmp.append(x)
      tmp.append(0)
      data.append(tmp)
    if y0>y1:
      step=-1
    else:
      step=1
    for y in range(y0,y1+1,step):
      for x in range(x0,x1+1):
        data[x-x0][1]+=SUMUP_ARRAY[fs][y][x]

    filename="cut-%s-%s.dat"%(int(counter-4),FS_DICT[fs])
    fp=open(filename,'w')
    for i in range(len(data)):
      tmp="%d %f\n" %(int(data[i][0]),data[i][1])
      fp.write(tmp)
    fp.close()

  return



def create_png_images(base,sequence,filelist,alldata,maxv,gdfontpath,logz,roi,roffset,fpng,bglin,mansel,maxi,mini,smaxi,smini,sdet,sens_det_array,sumup,gisans,window,vert):
  global  SUMUP_ARRAY,FS_DICT

  counter=0
  for seq in sequence:

    if len(roffset)>1:
      offset=roffset[counter]
    else:
      offset=roffset[0]

    for filename in filelist[seq]:
      print( filename)
#      basename=base+"%d" %(seq)
      number=get_file_number(filename)
      targetfile="%s%d.%04d.png" %(base,seq,number)

      fs=create_number_from_flipper_states(alldata,seq,number)
      if NICOS==1:
        time=float(alldata[seq][number][5]['timer'][0])
      else:
        time=float(alldata[seq][number][5]['time'])

      if vert!=0:
        loffset=int(offset-int((tan(2*radians(vert+float(alldata[seq][number][5]['rx'])))*DIST_SAMP_DET/CHANNELWIDTH)))
      else:
        loffset=offset
      
      if(os.path.isfile(targetfile)==False) or fpng==1:
        (fpt, tfilename) = tempfile.mkstemp()
        os.chmod(tfilename, 0o664)
        os.close(fpt)

        tmp="cp -f %s %s.gz" %(filename,tfilename)
        os.system(tmp)
        
        create_tmp_img(tfilename,sens_det_array,sumup,window,fs,time,vert)
#        if (sdet==0) and sumup==0:
#          tmp="gzip -d tmp.img.gz"
#          os.system(tmp)
#        else:
#          create_tmp_img("tmp.img.gz",sens_det_array,sumup,window,fs,time)

        make_gnuplot_pngfile(base,number,seq,filename,alldata,maxv,logz,roi,loffset,bglin,mansel,maxi,mini,gisans,tfilename)

        tmp="export GDFONTPATH=\"%s\"; gnuplot < %s.gpl" %(gdfontpath,filename)
        os.system(tmp)

        tmp="rm -f %s.gpl" %(filename)
        os.system(tmp)
#        os.remove(tfilename)
#        os.remove(tfilename+".gz")
        if sdet==1:          
#          print( tfilename)
#          print( "Sens_"+filename)
          if(os.path.isfile("Sens_"+filename)==True):
               os.remove("Sens_"+filename)
#          os.rename(tfilename,"Sens_"+filename[:-3])
          tmp="mv %s %s" %(tfilename,"Sens_"+filename[:-3])
          os.system(tmp)
          tmp="gzip %s" %("Sens_"+filename[:-3])
          os.system(tmp)
#          print( tmp)
          os.remove(tfilename+".gz")
        else:
          os.remove(tfilename)
          os.remove(tfilename+".gz")

  
  if sumup==1:
    for fs in range(4):
      print( "fs=",fs," : ",SUMUP_ARRAY[fs][1024][1024]," ",FS_DICT[fs])
      if SUMUP_ARRAY[fs][1024][1024]>0:
        make_sumup_pngfile(logz,smaxi,smini,gisans,window,fs)
        if len(window)>0:
          cut_out_window(window,fs)

        tmp="export GDFONTPATH=\"%s\"; gnuplot < %s_%s.gpl" %(gdfontpath,"sumup",FS_DICT[fs])
        os.system(tmp)
        tmp="gzip -c tmp_sumup_%s.img > sumup_%s.gz" %(FS_DICT[fs],FS_DICT[fs])
        os.system(tmp)
        tmp="rm -f tmp_sumup_%s.img" %(FS_DICT[fs])
        os.system(tmp)
    
  return



def create_avi_movie(base,sequence):

  if platform.system()=='MacOS':
    tmp="mencoder -really-quiet \"mf://%s[%d-%d].*.png\" -mf fps=10 -o %s%d-%d.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2000" %(base,sequence[0],sequence[len(sequence)-1],base,sequence[0],sequence[len(sequence)-1])
    print( tmp)
    os.system(tmp)
  elif platform.system()=='Linux':
    globstring=""
    for i in range(sequence[0],sequence[len(sequence)-1]):
      globstring+="-i '%s%d.*.png' " %(base,i)
    tmp="ffmpeg -pattern_type glob %s' -vcodec mpeg2video -pix_fmt yuv420p %s%d-%d.avi" %(globstring,base,sequence[0],sequence[len(sequence)-1])
    print( tmp)
    tmp="ffmpeg -pattern_type glob -i '%s[%d-%d].*.png' -vcodec mpeg2video -pix_fmt yuv420p %s%d-%d.avi" %(base,sequence[0],sequence[len(sequence)-1],base,sequence[0],sequence[len(sequence)-1])
    print( tmp)
    os.system(tmp)

  return



def get_polarization(filename):
  polflag="sf"
  i=0

  fp=open(filename,'r')
  liste=fp.readlines()
  fp.close()

  while(liste[i].find("polarization")==-1):
    i=i+1

  tmp=liste[i].split()
#  print( "polflag",tmp," : ",tmp[2][:-1])
  polflag=tmp[2]

  return polflag


def get_nicos_savetime(filename):

  i=0

  fp=open(filename,'r')
  liste=fp.readlines()
  fp.close()

  while(liste[i].find("stopped:")==-1):
    i=i+1

  tmp=liste[i].split()
  savedtime=tmp[1][1:-1]

  return savedtime


def get_savetime(filename):

  i=0

  fp=open(filename,'r')
  liste=fp.readlines()
  fp.close()

  while(liste[i].find("saved:")==-1):
    i=i+1

  tmp=liste[i]
#  print( "saved",tmp[6:-1])
  savedtime=tmp[6:-1]

  return savedtime


def get_roi_det(filename):

  rois=[]
  filename=filename[:rfind(filename,"_")]+".set"
  fp=open(filename,'r')
  liste=fp.readlines()
  fp.close()

  i=1
  while(liste[i].find("roi")==-1):
    i=i+1

  for j in range(1,7):
#    print( "roi:",liste[i])
    if liste[i+5].find("true")!=-1 :
#      print( "liste[i+5].find:",liste[i+5])
      tmp=[]
#      print( "number:",liste[i][3:4])
      tmp.append(int(liste[i][3:4]))
      string=liste[i+3]
      string=string.replace("-"," ")
      string=string.replace(","," ")
      string=string.replace("["," ")
      string=string.replace("]"," ")
      tmps=string.split()
      tmp.append(tmps)
      rois.append(tmp)
    i=i+7

#  print( rois)

  return rois



def get_nicos_roi_det(filename):

  rois=[]
  for roi in range(1,7):
    tmp=[]
    tmp.append(roi)
    tmp.append(["1","1","2","2"])
    rois.append(tmp)

  return rois


def write_data_to_file_neu(base,sequence,filelist,alldata,mansel):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,ALL_CHANNELS

  filename="%s%d-%d.dat" %(base,sequence[0],sequence[len(sequence)-1])

  np.savez(filename,alldata)
  return



def get_alldata_from_file_neu(base,sequence,filelist):

  nlist=[]
  rlist=[]
  oamaxv=0
  nseq=sequence[:]

  filename="%s%d-%d.dat" %(base,sequence[0],sequence[len(sequence)-1])
  alldata=np.load(filename+".npz").items()[0][1]

#  print( alldata)

  start =0
  for seq in sequence:
    for i in range(start,seq+1):
      nlist.append([])
      
    start=seq
    for file in filelist[seq]:
      flag=0
      for entry in rlist:
        if file==entry[0]:
          flag=1
          break

      if flag==0:
        nlist[seq].append(file)

    for tmp in alldata[seq]:
      maxv=tmp[4][0]
      if maxv>oamaxv:
        oamaxv=maxv

  return alldata,oamaxv,nlist,nseq



def call_divdet(base,sequence,sumup,filelist,divmaxi,divmini,divsmaxi,divsmini,rev,logz):

  if logz==1:
    slogz="-logz"
  else:
    slogz=""
    
  for seq in sequence:
    counter=0
    for filename in filelist[seq]:
      if counter%2==0:
        ofilename=filename
        counter+=1
      else:
        if rev==1:
          tmp=ofilename
          ofilename=filename
          filename=tmp
        #tmp="divdetimg %s %s" %(ofilename,filename)
        tmp="divdetimg_simple.py %s -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(slogz,divmini,divmaxi,ofilename,filename)
        print( "call:", tmp)
        os.system(tmp)
        print( tmp)
        counter+=1

  if sumup==1:
    print( "SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[1][1024][1024]>0",SUMUP_ARRAY[2][1024][1024],SUMUP_ARRAY[1][1024][1024])
    if SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[1][1024][1024]>0:
      tmp="divdetimg_simple.py %s -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(slogz,divsmini,divsmaxi,"sumup_uu.gz","sumup_dd.gz")
      print( "call:", tmp)
      os.system(tmp)
    print( "SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[3][1024][1024]>0",SUMUP_ARRAY[2][1024][1024],SUMUP_ARRAY[3][1024][1024])
    if SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[3][1024][1024]>0:
      tmp="divdetimg_simple.py %s -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(slogz,divsmini,divsmaxi,"sumup_uu.gz","sumup_du.gz")
      print( "call:", tmp)
      os.system(tmp)
    print( "SUMUP_ARRAY[1][1024][1024]>0 and SUMUP_ARRAY[0][1024][1024]>0:",SUMUP_ARRAY[1][1024][1024], SUMUP_ARRAY[3][1024][1024])
    if SUMUP_ARRAY[1][1024][1024]>0 and SUMUP_ARRAY[0][1024][1024]>0:
      tmp="divdetimg_simple.py %s -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(slogz,divsmini,divsmaxi,"sumup_dd.gz","sumup_ud.gz")
      print( "call:", tmp)
      os.system(tmp)


  return



def call_subdet(base,sequence,sumup,filelist,submaxi,submini,smaxi,rev):

  print( "Inside subdet")
  for seq in sequence:
    counter=0
    for filename in filelist[seq]:
      if counter%2==0:
        ofilename=filename
        counter+=1
      else:
        if rev==1:
          tmp=ofilename
          ofilename=filename
          filename=tmp
        #tmp="subdetimg.py %s %s" %(ofilename,filename)
        tmp="subdetimg_simple.py -bin 4 -mini %d -maxi %d -f1 %s -f2 %s" %(submini,submaxi,ofilename,filename)
        print( "call:", tmp)
        os.system(tmp)
        print( tmp)
        counter+=1

  if sumup==1:
    print( "SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[1][1024][1024]>0",SUMUP_ARRAY[2][1024][1024],SUMUP_ARRAY[1][1024][1024])
    if SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[1][1024][1024]>0:
      tmp="subdetimg_simple.py -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(-smaxi/2.0,smaxi/2.0,"sumup_uu.gz","sumup_dd.gz")
      print( "call:", tmp)
      os.system(tmp)
    print( "SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[3][1024][1024]>0",SUMUP_ARRAY[2][1024][1024],SUMUP_ARRAY[3][1024][1024])
    if SUMUP_ARRAY[2][1024][1024]>0 and SUMUP_ARRAY[3][1024][1024]>0:
      tmp="subdetimg_simple.py -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(-smaxi/2.0,smaxi/2.0,"sumup_uu.gz","sumup_du.gz")
      print( "call:", tmp)
      os.system(tmp)
    print( "SUMUP_ARRAY[1][1024][1024]>0 and SUMUP_ARRAY[0][1024][1024]>0:",SUMUP_ARRAY[1][1024][1024], SUMUP_ARRAY[3][1024][1024])
    if SUMUP_ARRAY[1][1024][1024]>0 and SUMUP_ARRAY[0][1024][1024]>0:
      tmp="subdetimg_simple.py -bin 4 -mini %f -maxi %f -f1 %s -f2 %s" %(-smaxi/2.0,smaxi/2.0,"sumup_dd.gz","sumup_ud.gz")
      print( "call:", tmp)
      os.system(tmp)
      
  return




def get_sens_det_image(sens_det_image):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL,V_MAX_CHANNEL, ALL_CHANNELS, NICOS

  fp=gzip.open(sens_det_image,"r")
  array=[]
  sens=[]
  det=[]
  detr=[]

  for line in fp:
    det.append(line)
    
  for i in range(0,ALL_CHANNELS):
    #at the moment identical, will change ...
    if NICOS!=1:
      detr.append(det[ALL_CHANNELS-i-1])
    else:
      detr.append(det[ALL_CHANNELS-i-1])

  for line in detr:
    tmp=line.split()
    array.append(tmp)

  #find maximum in defined window
  value=0

  print( "Number of channels horizontally:",len(array)," and vertically:",len(array[0]))
  print( "Using from H_MIN_CHANNEL=", H_MIN_CHANNEL," to H_MAX_CHANNEL=",H_MAX_CHANNEL," and from V_MIN_CHANNEL=",V_MIN_CHANNEL," to V_MAX_CHANNEL=",V_MAX_CHANNEL)
  
  for i in range(H_MIN_CHANNEL,H_MAX_CHANNEL):
    for j in range(V_MIN_CHANNEL,V_MAX_CHANNEL):
      value+=float(array[i][j])

  aver=float(value)/float(H_MAX_CHANNEL-H_MIN_CHANNEL)/float(V_MAX_CHANNEL-V_MIN_CHANNEL)
  print( "aver aus get_sens_det_image=",aver)

  for i in range(0,ALL_CHANNELS):
    tmp=[]
    for j in range(0,ALL_CHANNELS):
      tmp.append(1.0)
    sens.append(tmp)

  for i in range(H_MIN_CHANNEL,H_MAX_CHANNEL):
    for j in range(V_MIN_CHANNEL,V_MAX_CHANNEL):
      sens[i][j]=float(array[i][j])/float(aver)


  fp.close()

  return sens



def redo_bg_correction(base,sequence,filelist,alldata,roi,roffset,bglin,nmon,fc):
  global H_MIN_CHANNEL,H_MAX_CHANNEL,V_MIN_CHANNEL,V_MAX_CHANNEL, ALL_CHANNELS,NICOS


  nhsum=[]
  for i in  range(0,ALL_CHANNELS):
    nhsum.append(0)

  count=0
  counter=0
  for seq in sequence:

    if len(roffset)>1:
      offset=roffset[counter]
    else:
      offset=roffset[0]

    for filename in filelist[seq]:
      count=count+1
      number=get_file_number(filename)
      if NICOS==1:
        rois=get_nicos_roi_det(filename)
      else:
        rois=get_roi_det(filename)
      paramdic=alldata[seq][number][5]
      pfs = int(alldata[seq][number][5]['pflipper'])
      if 'aflipper' in alldata[seq][number][5]:
        if is_number(alldata[seq][number][5]['aflipper']):
          afs = int(alldata[seq][number][5]['aflipper'])
        else:
          afs=0
      apos = int(alldata[seq][number][5]['analyzer_shift'])

      if apos==0:
        print( "analyzer was not in the beam, so this correction makes no sense!")
        print( "please rerun it and nad remove the switch -coh")
        sys.exit(-1)
      elif pfs!=afs:
        for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
          nhsum[i]=float(alldata[seq][number][2][i][2]-0.5*alldata[seq][number+1][2][i][2])

        n2hsum=find_bg_hsum(nhsum,bglin,roi,offset,0)
        intensity,sigma=calc_intensity(n2hsum,roi,offset,paramdic,nmon,fc)

        alldata[seq][number][4][3]=intensity
        alldata[seq][number][4][4]=sigma

        for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
          alldata[seq][number][2][i][2]=n2hsum[i][2]
          alldata[seq][number][2][i][3]=n2hsum[i][3]
          alldata[seq][number][2][i][4]=n2hsum[i][4]

      else:
        for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
          nhsum[i]=float(1.5*alldata[seq][number][2][i][2])

        n2hsum=find_bg_hsum(nhsum,bglin,roi,offset,0)
        intensity,sigma=calc_intensity(n2hsum,roi,offset,paramdic,nmon,fc)

        alldata[seq][number][4][3]=intensity
        alldata[seq][number][4][4]=sigma

        for i in  range(H_MIN_CHANNEL,H_MAX_CHANNEL):
          alldata[seq][number][2][i][2]=n2hsum[i][2]
          alldata[seq][number][2][i][3]=n2hsum[i][3]
          alldata[seq][number][2][i][4]=n2hsum[i][4]
        
    counter+=1

  return



def lin(x, x0, m, c):
  return (x-x0)*m+c


def gaus(x, a, x0, sigma1):
  return a*exp(-(x-x0)**2/(2*sigma1**2))


def fit3_row(row,lroi,limit,scale):
  global bg
  
  count=0
  bgl=0
  bgr=0
  bgrc=0
  bglc=0
  cp=0
  cn=0
  
  for i in range(len(row)):
    if fabs(row[i][0])>lroi and fabs(row[i][0])<limit:
      count+=1
      
#    print( "count",count)
  if count==0:
    print( "Find new limit, no background in this range!")
    for i in range(len(row)):
      if fabs(row[i][0]+limit)<(row[cn][0]+limit):
        cn=i
      if fabs(row[i][0]-limit)<(row[cp][0]-limit):
        cp=i
      if fabs(row[cp][0])>fabs(row[cp][0]):
        limit=fabs(row[cp][0])
      else:
        limit=fabs(row[cn][0])
      for i in range(len(row)):
        if fabs(row[i][0])>lroi and fabs(row[i][0])<limit:
          count+=1

  xo=ar(range(count), dtype=np.float32)
  yo=ar(range(count), dtype=np.float32)
  count=0
  for i in range(len(row)):
    if fabs(row[i][0])>lroi and fabs(row[i][0])<limit:
      xo[count]=row[i][0]
      yo[count]=row[i][1]
      count+=1

  for i in range(len(row)):
    if row[i][0]>-limit and row[i][0]<-lroi:
      bgr+=row[i][1]
      bgrc+=1
    if row[i][0]<limit and row[i][0]>lroi:
      bgl+=row[i][1]
      bglc+=1

#    print( "bgrc=",bgrc," bglc=",bglc," bglc+bgrc=",bglc+bgrc)
#    print( "bgr=",bgr," bgl=",bgl)
  if bgrc==0:
    if bglc==0:
      bgr=0
    else:
      bgr=bgl/bglc
  else:
    bgr=bgr/bgrc

  if bglc==0:
    if bgrc==0:
      bgl=0
    else:
      bgl=bgr/bgrc
  else:
    bgl=bgl/bglc
  
#    print( "bgr=",bgr," bgl=",bgl)
  
  m=(bgr-bgl)/(((limit-lroi)/2.0)+lroi)
  if bglc>0 and bgrc>0:
    mean=(bgr+bgl)/2.0
  elif bglc>0:
    mean=bgl
  elif bgrc>0:
    mean=bgr
  else:
    mean=0
  
#    print( )
  
  if bglc==0 or bgrc==0:
    if bglc==0:
      print( "Fix background to left side, not enough points for fit")
      poptb=[0,0,mean]
    elif bgrc==0:
      print( "Fix background to right side, not enough points for fit")
      poptb=[0,0,mean]
    elif bglc==0 and bgrc==0:
      print( "Fix background to 0, no points found")
      poptb=[0,0,0]
  else:
    try:
      poptb,pcovb = curve_fit(lin,xo,yo,p0=[0,m,mean])
#            print( "mean=",mean," m=",m)
    except:
      print( "Background fit didn't work, taking fix points")
      poptb=[0,0,(bgr+bgl)/2.0]
      mean=(bgr+bgl)/2.0

  count=0
  for i in range(len(row)):
    if fabs(row[i][0])<limit:
      count+=1

  xc=ar(range(count), dtype=np.float32)
  yc=ar(range(count), dtype=np.float32)
  sc=ar(range(count), dtype=np.float32)
  max=0
  count=0
  for i in range(len(row)):
    if fabs(row[i][0])<limit:
      xc[count]=row[i][0]
      yc[count]=row[i][1]-((row[i][0]-poptb[0])*poptb[1]+poptb[2])
      sc[count]=scale[i][1]
      if yc[count]>max:
        max=yc[count]
      count+=1

  nc=len(xc)
  if sum(yc)!=0:
    mean=sum(xc*yc)/sum(yc)
    if fabs(mean)>lroi:
      mean=0
    if sum(yc*(xc - mean)**2)/sum(yc)>0:
      sigma=np.sqrt(sum(yc*(xc - mean)**2)/sum(yc))
    else:
      #        print( "sum is negative, sigma=1")
      sigma=1
    cmean=mean
    csigma=sigma
    cmax=max
  else:
    mean=0
    sigma=1
    cmean=mean
    csigma=sigma
    cmax=max

  failed=0
  try:
    popt,pcov = curve_fit(gaus,xc,yc,p0=[max,mean,sigma])
    if popt[2]>20:
      print( "Peak fit failed, sigma too large")
#      print( "amp=",popt[0]," x0=",popt[1]," sigma=",popt[2])
      failed=1
      popt=[cmax,cmean,csigma]            
      intensity,sigma=calc_simple_intensity(xc,yc,sc)
    elif fabs(popt[1])>lroi:
      print( "Peak fit failed center of gaussian fit is outside lroi")
#      print( "amp=",popt[0]," x0=",popt[1]," sigma=",popt[2])
      failed=1
      popt=[cmax,cmean,csigma]            
      intensity,sigma=calc_simple_intensity(xc,yc,sc)
    elif popt[0]<0:
      print( "Peak fit failed amplitude of gaussian fit is negative")
#      print( "amp=",popt[0]," x0=",popt[1]," sigma=",popt[2])
      failed=1
      popt=[cmax,cmean,csigma]            
      intensity,sigma=calc_simple_intensity(xc,yc,sc)
    else:
      intensity,sigma=calc_simple_intensity(xc,yc,sc)
      intensity=0
      sigma=0
      scaver=0
      scc=0
      for i in sc:
        scaver+=i
        scc+=1
      scaver=scaver/scc
      #      for i in range(-lroi,lroi+1,1):
      for i in range(-2*lroi,2*lroi+1,1):
        intensity+=gaus(i,popt[0],popt[1],popt[2])
        sigma+=gaus(i,popt[0],popt[1],popt[2])*scaver
      sigma=sqrt(sigma)/scaver
      failed=0
  except:
    failed=1
    popt=[cmax,cmean,csigma]
    intensity,sigma=calc_simple_intensity(xc,yc,sc)

#  print( "amp=",popt[0]," max=",max," x0=",popt[1]," cmean=",cmean," sigma1=",popt[2]," csigma=",csigma )

  return poptb,popt,failed,intensity,sigma


def calc_simple_intensity(xc,yc,sc):
  intensity=0
  sigma=0
  for i in range(len(xc)):
    if i>0 and i<len(xc)-1:
      width=((xc[i+1]-xc[i])/2.0)+((xc[i]-xc[i-1])/2.0)
    elif i>0:
      width=xc[i]-xc[i-1]
    elif i< len(xc)-1:
      width=xc[i+1]-xc[i]
    else:
      print( "Error in width, stopping")
      exit(0)
    intensity+=yc[i]*width
    sigma+=(yc[i]*width*sc[i])/(sc[i]*sc[i])
    
  if sigma>0:
    sigma=sqrt(sigma)
  else:
    sigma=0

  return intensity,sigma


def plot_row(row,y,result,lroi,limit,filename,fcounter):
  global GNUFONT,bg
  
  zmin=100
  zmax=0
  
  for i in range(len(row)):
    if row[i][1]>zmax:
      zmax=row[i][1]
    if zmin>row[i][1]:
      zmin=row[i][1]
            
  fp=open("plot_it.gpl","w")
  
  tmp="#!/usr/bin/gnuplot\n"
  fp.write(tmp)
  tmp="set encoding iso_8859_1\n"
  fp.write(tmp)
  tmp="set terminal postscript enhanced color landscape\n"
  tmp="set terminal pngcairo truecolor noenhanced font \"%s,36\" size 2048,2048 lw 2\n" %(GNUFONT)
  fp.write(tmp)
  tmp="set output \"%s_%03d_%s.png\"\n" %(filename[:-4],fcounter,str(round(y[0],3)))
  fp.write(tmp)
  tmp="set logscale cb\n"
  fp.write(tmp)
  tmp="set logscale z\n"
  fp.write(tmp)
  tmp="set xlabel \"Q [\305^-^1]\"\n"
  fp.write(tmp)
  tmp="set title  \"%s  %s\" font \"%s,36\" \n" %(str(round(y[0],3)),str(round(y[1],3)),GNUFONT)
  fp.write(tmp)
  tmp="set arrow from %d,%f to %d,%f lt 4 lw 4 front nohead\n" %(lroi,zmin,lroi,zmax)
  fp.write(tmp)
  tmp="set arrow from %d,%f to %d,%f lt 4 lw 4 front nohead\n" %(-lroi,zmin,-lroi,zmax)
  fp.write(tmp)
  tmp="set arrow from %d,%f to %d,%f lt 4 lw 4 front nohead\n" %(limit,zmin,limit,zmax/4.0)
  fp.write(tmp)
  tmp="set arrow from %d,%f to %d,%f lt 4 lw 4 front nohead\n" %(-limit,zmin,-limit,zmax/4.0)
  fp.write(tmp)
  
  tmp="a=%f\n" %(result[0])
  fp.write(tmp)
  tmp="xa0=%f\n" %(result[1])
  fp.write(tmp)
  tmp="sigma1=%f\n" %(result[2])
  fp.write(tmp)
  tmp="xb0=%f\n" %(result[3])
  fp.write(tmp)
  tmp="m=%f\n" %(result[4])
  fp.write(tmp)
  tmp="c=%f\n" %(result[5])
  fp.write(tmp)
  #    tmp="g(x)=a*exp(-(x-x0)**2/(2*sigma1**2))\n"
  #    fp.write(tmp)
  tmp="g(x)=a*exp(-(x-xa0)**2/(2*sigma1**2))\n"
  fp.write(tmp)
  tmp="b(x)=(((x-xb0)*m)+c)\n"
  fp.write(tmp)
  tmp="set samples 1000\n"
  fp.write(tmp)

  tmp="plot \"plot_it.dat\" u 1:2 w lp lw 2, g(x), b(x),g(x)+b(x)\n"
  fp.write(tmp)
  fp.close()
  
  fp=open("plot_it.dat","w")
  for a in range(len(row)):
    tmp="%f %f\n" %(row[a][0],row[a][1])      
    fp.write(tmp)
  fp.close()
  
  tmp="gnuplot plot_it.gpl" 
  os.system(tmp)
  
  print( "plotting: %s_%03d_%s" %(filename[:-4],fcounter,str(round(y[0],3))))
  
  tmp="rm -f plot_it.gpl"
  os.system(tmp)
  
  tmp="rm -f plot_it.dat"
  os.system(tmp)
   
#    time.sleep(.3)
  return


def extract_specular(allrow,allres,lroi,limit,filename,scale):

    filename="%s-specrefnew.dat" %(filename[:-4])
    print( filename)
    print( "scale=",scale)
    fp=open(filename,'w')

    maxv=0
    for i in range(len(allres)):
      if allres[i][6]==0:
        intensity=allres[i][7]
        if intensity>maxv:
          maxv=intensity

    for i in range(len(allres)):
      if allres[i][6]==0:
        q=allres[i][9][0]/1000.0
        intensity=allres[i][7]
        tmp="%f %g %g\n" %(q,intensity/maxv*scale,allres[i][8]/maxv)
        fp.write(tmp)
    fp.close()

    filenamegpl="%s.gpl" %(filename[:-4])
    filenameps="%s.ps" %(filename[:-4])
    print( "filename:", filename)
    print( "filenamegpl:", filenamegpl)
    print( "filenameps:", filenameps)

    fp=open(filenamegpl,'w')
    tmp="#!/usr/bin/gnuplot\n"
    fp.write(tmp)
    tmp="set encoding iso_8859_1\n"
    fp.write(tmp)
    tmp="set terminal postscript noenhanced color \"%s\" 18 landscape lw 2\n" %(GNUFONT)
    fp.write(tmp)
    tmp="set output \"%s\"\n" %(filenameps)
    fp.write(tmp)
    tmp="set logscale y\n"
    fp.write(tmp)
    tmp="set xlabel \"Q [\305^-^1]\"\n"
    fp.write(tmp)
    tmp="plot \"%s\" u 1:2:3 t \"%s\" with errorlines" %(filename,filename)
    fp.write(tmp)
    fp.close()

    if platform.system()=='Linux':
      gdfontpath="/usr/share/fonts/dejavu"
    elif platform.system()=='MacOS':
      gdfontpath="/usr/local/TeX/texmf-dist/fonts/truetype/public/dejavu/"

    tmp="export GDFONTPATH=\"%s\"; gnuplot < %s" %(gdfontpath,filenamegpl)
    os.system(tmp)
    print( tmp)

    return



def extract_specref_from_taiaf(filename,roi,lroffset,overallscale):
  global DIST_SAMP_DET,CHANNELWIDTH

  lroi=((roi[0]*CHANNELWIDTH)/DIST_SAMP_DET)*2.0*1000
  lroi=12
  limit=2*lroi
  failed=0

  map=[]
  allres=[]
  allrow=[]
  allscale=[]

  print( "Filename=",filename)
  fp=open(filename,"r")  
  counter=0
  maxcount=0
  for line in fp.readlines():
    tmp=[]
    tmp.append(counter)
    if len(line.split())>0:
      for v in line.split():
        tmp.append(float(v))
      map.append(tmp)
      counter+=1
    else:
      if counter>maxcount:
        maxcount=counter
      counter=0
  fp.close

#  print( len(map))
  map2=[]
  pos1=[]
  pos2=[]
  diff=[]
  counter==0
  while counter==map[counter][0]:
    pos1.append(map[counter][1])
    counter+=1

  c2=1
  for x in range(counter,len(map)):
    pos1[map[x][0]]+=map[x][1]
    if map[x][0]==0:
      c2+=1

  zpos=0
  for y in range(0,len(pos1)):
    pos1[y]=round(pos1[y]/c2,3)
    if fabs(pos1[y])<pos1[zpos]:
      zpos=y

  deltah=fabs(pos1[zpos]-pos1[zpos-1])

  vert=[]
  for x in range(counter,len(map)):
    if fabs(map[x][1]-pos1[zpos])< deltah/2.0 and map[x][2] >0:
      tmp=[]
      tmp.append(round(map[x][2],3))
      tmp.append(0)
      vert.append(tmp)


  for y in range(0,len(vert)-1):
    vert[y][1]=vert[y+1][0]-vert[y][0]

  deltav=vert[0][1]

  fcounter=0
  for y in vert:
    fcounter+=1
    row=[]
    scale=[]
    for x in range(len(map)):
      tmp=[]
      if fabs(map[x][2]-y[0])< deltav/2.0:
        tmp.append(round(map[x][1],3))
        tmp.append(map[x][3])
#        tmp.append(map[x][7])
        row.append(tmp)
        tmp=[]
        tmp.append(map[x][8])
        tmp.append(map[x][9])
        scale.append(tmp)

    row.sort()
    rb,rp,failed,intensity,sigma=fit3_row(row,lroi,limit,scale)
    result=[]
    result.append(rp[0]) #amplitude
    result.append(rp[1]) #x0
    result.append(rp[2]) #sigma
    result.append(rb[0]) #x0 
    result.append(rb[1]) #m steigung
    result.append(rb[2]) #c constante
    result.append(failed)
    result.append(intensity)
    result.append(sigma)
    result.append(y)
    allres.append(result)
    allrow.append(row)
    allscale.append(scale)
#    print( "plotting:",filename," y=",y)
    plot_row(row,y,result,lroi,limit,filename,fcounter)

  extract_specular(allrow,allres,lroi,limit,filename,overallscale)
  asigma=0.0
  ax0=0.0
  bx0=0.0
  for i in range(len(allres)):
    asigma+=allres[i][2]
    ax0+=allres[i][1]
    bx0+=allres[i][3]
  print( "Values for extractet specular reflectivity from %s:" %(filename))
  print( "Aver sigma=", asigma/float(len(allres)))
  print( "Aver peakx0==", ax0/float(len(allres)))
  print( "Aver backx0==", bx0/float(len(allres)))

  if platform.system()=='MacOS':
    tmp="mencoder -really-quiet \"mf://%s*.*.png\" -mf fps=10 -o %s-specrefnew.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=2000" %(filename[:-4],filename[:-4])
  elif platform.system()=='Linux':
    tmp="ffmpeg -r 1/2 -pattern_type glob -i '%s*.*.png' -vcodec mpeg2video -pix_fmt yuv420p -r 25 %s/-specrefnew.avi" %(filename[:-4],filename[:-4])
    tmp="ffmpeg -pattern_type glob -i '%s*.*.png' -vcodec mpeg2video -pix_fmt yuv420p %s/-specrefnew.avi" %(filename[:-4],filename[:-4])
  print( tmp)
  os.system(tmp)
  tmp="rm -f %s*.*.png" % (filename[:-4])
  os.system(tmp)
  print( tmp)

  return



def main_part():
  a = datetime.datetime.now()
  global ALL_CHANNELS,GNUFONT, NICOS

  max =0
  tmp=""
  cycle=0
  lines=0
  grey=0
  logz=0
  base="*"
  gdfontpath=""
  do_not_write_back=0
  divisors=[]

#  if(os.path.isdir("/usr/share/fonts/truetype/ttf-dejavu")==True):
#    gdfontpath="/usr/share/fonts/truetype/ttf-dejavu"
#    GNUFONT="DejaVuSans"
#  else:
#    print( os.path.isdir("/usr/local/TeX/texmf-dist/fonts/truetype/public/dejavu/"))
#    if(os.path.isdir("/usr/local/TeX/texmf-dist/fonts/truetype/public/dejavu/")==True):
#      gdfontpath="/usr/local/TeX/texmf-dist/fonts/truetype/public/dejavu/"
#      GNUFONT="Courier"
#    else:
#      gdfontpath="/usr/share/fonts/dejavu"
#
#  print( gdfontpath)

  if platform.system()=='Linux':
    print( "Platform Linux")
    gdfontpath="/usr/share/fonts/dejavu"
    GNUFONT="DejaVuSans"
  elif platform.system()=='MacOS':
    print( "Platform MacOS")
    gdfontpath="/usr/local/TeX/texmf-dist/fonts/truetype/public/dejavu/"
    GNUFONT="Courier"
  else:
    print( "Platform unknown")
    GNUFONT="DejaVuSans"
  
  
  base,sequence,divisors,logz,nopng,noavi,noaiaf,aiafmm,bglin,roi,roffset,fread,fpng,favi,faiaf,nmon,scale,sfc,fc,mansel,divdet,subdet,noref,maxi,mini,smaxi,smini,submaxi,submini,divmaxi,divmini,divsmaxi,divsmini,sdet,sens_det_image,sumup,gisans,window,rev,leak,noq,genx,coh,specfromaiaf,vert=get_switches()

  if sdet==1:
    sens_det_array=get_sens_det_image(sens_det_image)
  else:
    sens_det_array=[[0 for i in range(ALL_CHANNELS)] for j in range(ALL_CHANNELS)]

  if NICOS==1:
    filelist=get_nicos_file_list(base,sequence)
  else:
    filelist=get_file_list(base,sequence)
  

#  targetfile="%s%d-%d.dat" %(base,sequence[0],sequence[len(sequence)-1])
  targetfile="%s%d-%d.dat" %(base,sequence[0],sequence[len(sequence)-1])

  if os.path.isfile(targetfile)==True and fread==0:
#    alldata,maxv,short_filelist,short_sequence=get_alldata_from_file(base,sequence,filelist)
    print( "get all data from file neu")
    alldata,maxv,short_filelist,short_sequence=get_alldata_from_file_neu(base,sequence,filelist)
    print( "append data")
    alldata,maxv=append_analyse_data(base,short_sequence,short_filelist,bglin,roi,roffset,alldata,maxv,nmon,sfc,fc,sdet,sens_det_array,mansel,vert)
#    write_data_to_file(base,sequence,filelist,alldata,mansel)
    print( "write data to file neu")
    write_data_to_file_neu(base,sequence,filelist,alldata,mansel)
    print( "finished")
  else:
    alldata,maxv=analyse_data(base,sequence,filelist,bglin,roi,roffset,nmon,sfc,fc,sdet,sens_det_array,mansel,vert)
#    write_data_to_file(base,sequence,filelist,alldata,mansel)
    write_data_to_file_neu(base,sequence,filelist,alldata,mansel)
	

  print( "main part: noref=",noref," coh=",coh)
  if noref==0:
    if coh==1:
      redo_bg_correction(base,sequence,filelist,alldata,roi,roffset,bglin,nmon,fc)
    make_specref_file(base,sequence,filelist,alldata,divisors,nmon,scale,sfc,fc,mansel,noq,genx,leak,vert)

  if nopng==0 or fpng==1:
    create_png_images(base,sequence,filelist,alldata,maxv,gdfontpath,logz,roi,roffset,fpng,bglin,mansel,maxi,mini,smaxi,smini,sdet,sens_det_array,sumup,gisans,window,vert)

  if noavi==0 or favi==1:
    create_avi_movie(base,sequence)
    
  if noaiaf==0 or faiaf==1:
    make_qx_qz_map(base,sequence,filelist,alldata,gdfontpath,aiafmm,nmon,roffset,leak,divisors,roi,sfc,fc,mansel,specfromaiaf,scale,vert,noq)

  if divdet==1:
    call_divdet(base,sequence,sumup,filelist,divmaxi,divmini,divsmaxi,divsmini,rev,logz)

  if subdet==1:
    call_subdet(base,sequence,sumup,filelist,submaxi,submini,smaxi,rev)

  b = datetime.datetime.now()
  c = b-a
  print( "Operation took: "+str(c.seconds)+" seconds")
  
#main_part()
cProfile.run('main_part()', 'fooprof')
p = pstats.Stats('fooprof')
p.sort_stats('time')
#p.print(_stats())
