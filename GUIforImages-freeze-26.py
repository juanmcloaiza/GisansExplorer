#!/usr/local/bin/python
#!/usr/bin/env python
from Tkinter import *
from tkMessageBox import *
from tkFileDialog import *
from string import *
#from os import *
from glob import *
import sys
import os
import platform
import yaml


VERS="Pyfrid"
VERS="NICOS"
fileselection=[]
tmp_string=""

def answer(txt):
    showinfo("Error", txt)

def message(txt):
    showinfo("Message", txt)

def message_ask(txt):
    return askokcancel("Message", txt)

def open_filer_base():
    base=basename.get()
    print( "base:", base)
    print( base[len(base)-1:])
    if base[len(base)-1:]=="_":
        base=base[:base.rfind("/")]
#        print( "base:",base
    result=askdirectory(initialdir=base)
    print( "result:",result)
    if find(platform.system(),"Linux")==-1:
        #MACOS
#    print( "new result:",replace(result,"/tmp_mnt/srv/","/home/maria/")
        print( "MACOS")
        result=replace(result,"/tmp_mnt/srv/","/home/maria/")
    else:
        #LINUX
        print( "LINUX")
        result=replace(result,"/tmp_mnt/srv/","/home/maria/")
    basename.delete(0,END)
    basename.insert(0,result[:result.rfind("_")+1])

    #Updating the sequence
    globname=result[:result.rfind("_")+1]+"*"
    print( "globname:",globname)
    tmp=""
    for i in glob(globname):
        tmp+=i[i.rfind("_")+1:]+" "
    sequence.delete(0,END)
    sequence.insert(0,tmp[:-1])

    #Updateing the destination
    if find(platform.system(),"Linux")==-1:
        #MACOS
        print( "MACOS")
        dest=replace(result[:result.rfind("_")+1],"/home/maria/data/","/home/maria/Data_Treatment/")
    else:
        #LINUX
        print( "LINUX")
        dest=replace(result[:result.rfind("_")+1],"/tmp_mnt/srv/","/home/maria/Data_Treatment/")
    vsequence=str(sequence.get())
    seqlist=[]
    for i in vsequence.split(" "):
        seqlist.append(i)

    for i in sorted(seqlist,key=int):
        dest=dest+i+"-"

    if len(dest)>255: #name gets to long!
        if find(platform.system(),"Linux")==-1:
            #MACOS
            print( "MACOS")
            dest=replace(result[:result.rfind("_")+1],"/home/maria/data/","/home/maria/Data_Treatment/")
        else:
            #LINUX
            print( "LINUX")
            dest=replace(result[:result.rfind("_")+1],"/tmp_mnt/srv/","/home/maria/Data_Treatment/")
        vsequence=str(sequence.get())
        seqlist=[]
        for i in vsequence.split(" "):
            seqlist.append(i)
        j=0
        for i in sorted(seqlist,key=int):
            if j==0:
                dest=dest+i+"-"
            j+=1
        dest=dest+i+"-"

    destname.delete(0,END)
    destname.insert(0,dest[:-1])



def open_filer_base_NICOS():
    global fileselection

    base=basename.get()
    print( "base:",base)
    print( base[len(base)-1:])
    if base[len(base)-1:]=="_":
        base=base[:base.rfind("/")]
    print( "base2:",base)
    fileselection=askopenfilenames(initialdir=base,title='Choose files',filetypes=[('dat files', '*.dat'),('All files','*')])
    print( "fileselection:",fileselection)
    result=fileselection[0]
    if find(platform.system(),"Linux")==-1:
        #MACOS
        print( "MACOS")
        print( "not changing filepath selection",result)
        #result=replace(result,"/Users/mattauch/","/home/maria/")
    else:
        #LINUX
        print( "LINUX",result)
#        result=replace(result,"/tmp_mnt","/data")
        print( result)
    basename.delete(0,END)
    basename.insert(0,result[:result.rfind("_")+1])

    #Updating the sequence
    tmp=""
    for i in fileselection:
        tmp+=str(int(i[i.rfind("_")+1:i.rfind(".dat")]))+" "
    sequence.delete(0,END)
    sequence.insert(0,tmp[:-1])

    #Updating the destination
    result=fileselection[0]
    if find(platform.system(),"Linux")==-1:
        #MACOS
        print( "MACOS")
        #dest=replace(result[:result.rfind("_")+1],"/home/maria/data/","/home/maria/Data_Treatment/")
        print( result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1])
        dest=result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1]
    else:
        #LINUX
        print( "LINUX")
#        dest=replace(result[:result.rfind("_")+1],"/tmp_mnt/data/","/home/maria/Data_Treatment/")
#        dest=replace(result[:result.rfind("_")+1],"/tmp_mnt/data/","/home/jcns/Data_Treatment/")
        print( result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1])
        dest=result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1]
    vsequence=str(sequence.get())
    seqlist=[]
    for i in vsequence.split(" "):
        seqlist.append(i)

    for i in sorted(seqlist,key=int):
        dest=dest+i+"-"

    if len(dest)>255: #name gets to long!
        result=fileselection[0]
        if find(platform.system(),"Linux")==-1:
            #MACOS
            print( "MACOS")
            #dest=replace(result[:result.rfind("_")+1],"/home/maria/data/","/home/maria/Data_Treatment/")
            print( result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1])
            dest=result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1]
        else:
            #LINUX
            print( "LINUX")
#            dest=replace(result[:result.rfind("_")+1],"/tmp_mnt/data/","/home/maria/Data_Treatment/")
#            dest=replace(result[:result.rfind("_")+1],"/tmp_mnt/data/","/home/jcns/Data_Treatment/")
            print( result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1])
            dest=result[:result.rfind("data")]+"Data_Treatment"+result[result.rfind("data")+4:result.rfind("_")+1]
        vsequence=str(sequence.get())
        seqlist=[]
        for i in vsequence.split(" "):
            seqlist.append(i)
        j=0
        for i in sorted(seqlist,key=int):
            if j==0:
                dest=dest+i+"-"
            j+=1
        dest=dest+i+"-"

    destname.delete(0,END)
    destname.insert(0,dest[:-1])



def open_filer_dest():
    tmp=destname.get()
#    result=askokenfilename(initialdir=tmp[:tmp.rfind("/")])
    result=askdirectory(initialdir=destname.get())
    destname.delete(0,END)
    destname.insert(0,result)

def show_entry_fields():
    global fileselection,VERS
    fpng_flag=""
    sfc_flag=""
    roi_flag=""
    roffset_flag=""
    bglin_flag=""
    add_flag=""
    sim_flag=""

    tmp_string="images "
    tmp_string="imagesn2 "

    tmp_string+="-nmon "+str(vmonitors.get())+" "

    print( "gisans_mode.get()",vmode.get())
    if vmode.get()==1:
        tmp_string+="-gisans "

    if vpng.get()==1:
        print( "png flag set")
        tmp_string+="-fpng "
        fpng_flag="-fpng "
    else:
        print( "png flag not set")
        #tmp_string+="-nopng "
        fpng_flag=""

    if vread.get()==1:
        tmp_string+="-fread "
    else:
        tmp_string+=" "

    if vaiaf.get()==1:
        tmp_string+="-faiaf "
    else:
        tmp_string+="-noaiaf "

    if vavi.get()==1:
        tmp_string+="-favi "
    else:
        tmp_string+="-noavi "

    if vsens.get()==0:
        tmp_string+="-nosens "

    if vnoref.get()==1:
        tmp_string+="-noref "

    if vsfc.get()==1:
        if find(fcsl.get(),"Sample")==-1  and len(fcsl.get())>0:
            tmp_string+="-sfc "+str(fcsl.get())+" "
            sfc_flag="-sfc "

    if vfc.get()==1:
        if find(fcsl.get(),"Sample")==-1  and len(fcsl.get())>0:
            tmp_string+="-fc "+str(fcsl.get())+" "

    if scale.get()!=1:
        tmp_string+="-scale "+str(scale.get())+" "

    if find(div.get(),"1.0,x,y,z")==-1 and len(div.get())>0:
        tmp_string+="-div "+str(div.get())+" "

    if find(mansel.get(),"-1")==-1 and len(mansel.get())>0:
        tmp_string+="-mansel "+str(mansel.get())+" "

    if find(maxi.get(),"autoscale")==-1 and len(maxi.get())>0:
        tmp_string+="-maxi "+str(maxi.get())+" "
    if find(mini.get(),"autoscale")==-1 and len(mini.get())>0:
        tmp_string+="-mini "+str(mini.get())+" "

    if find(smaxi.get(),"autoscale")==-1 and len(smaxi.get())>0:
        tmp_string+="-smaxi "+str(smaxi.get())+" "
    if find(smini.get(),"autoscale")==-1 and len(smini.get())>0:
        tmp_string+="-smini "+str(smini.get())+" "

    if find(submaxi.get(),"autoscale")==-1 and len(submaxi.get())>0:
        tmp_string+="-submaxi "+str(submaxi.get())+" "
    if find(submini.get(),"autoscale")==-1 and len(submini.get())>0:
        tmp_string+="-submini "+str(submini.get())+" "

    if find(divmaxi.get(),"autoscale")==-1 and len(divmaxi.get())>0:
        print( "len",len(divmaxi.get()))
        tmp_string+="-divmaxi "+str(divmaxi.get())+" "
    if find(divmini.get(),"autoscale")==-1 and len(divmini.get())>0:
        tmp_string+="-divmini "+str(divmini.get())+" "

    if find(aiafmaxi.get(),"autoscale")==-1 and find(aiafmini.get(),"autoscale")==-1 and len(aiafmaxi.get())>0 and len(aiafmini.get())>0:
        tmp_string+="-aiafmm "+str(aiafmaxi.get())+" "+str(aiafmini.get())+" "

    if find(addstring.get(),"empty")==-1 and len(addstring.get())>0:
        tmp_string+=" "+str(addstring.get())+" "
        add_flag=" "+str(addstring.get())+" "

    if find(simstring.get(),"empty")==-1 and len(simstring.get())>0:
        sim_flag=" -sim "
#    print( "sim_flag=",sim_flag

    tmp_string+="-roffset "+str(roioff.get())+" "
    roffset_flag="-roffset "+str(roioff.get())+" "

    if len(roiw.get())>0 and  len(roih.get()) >0:
        tmp_string+="-roi "+str(roiw.get())+" "+str(roih.get())+" "
        roi_flag="-roi "+str(roiw.get())+" "+str(roih.get())+" "

    if len(bglin.get())>0:
        tmp_string+="-bglin "+str(bglin.get())+" "
        bglin_flag="-bglin "+str(bglin.get())+" "

    if vsumup.get()==1:
        tmp_string+="-sumup "

    if vdivdet.get()==1:
        tmp_string+="-divdet "

    if vsubdet.get()==1:
        tmp_string+="-subdet "

    if vmode.get()==3:
        tmps="kinimage3.py %s %s %s %s %s %s %s " %(simstring.get(),add_flag,roi_flag,roffset_flag,fpng_flag,sfc_flag,bglin_flag)
        vbasename=basename.get()
        tmp_basename=vbasename[vbasename.rfind("/")+1:]
        tmps+=" -base %s" %(tmp_basename)
        vsequence=str(sequence.get())
        print( vsequence.split())
        tmp_string=""
        for i in vsequence.split():
            if VERS=="NICOS":
                tmp_string+="%s%08d;" %(tmps,int(i))
            else:
                tmp_string+="%s%s;" %(tmps,i)

        print( "tmp_string kinimage:",tmp_string)

    elif vmode.get()==4:
        if VERS=="NICOS":
            tmp_string="refspec "
            print( "tmp_string NICOS:",tmp_string)
            if find(div.get(),"1.0,x,y,z")==-1 and len(div.get())>0:
                tmp_div=div.get()
            else:
                tmp_div=""
            print( "fileselection:",fileselection)
#             vbasename=basename.get()
#             vsequence=str(sequence.get())
#             vdestname=destname.get()
#             print( "vbasename:",vbasename
#             print( "vsequence:",vsequence
#             print( "vdestname:",vdestname
#             for i in vsequence.split(" "):
# #                os.chdir(vdestname)
#                 filename="%s%08d" %(vbasename,int(i))
#                 print( "filename:",filename
#                 command="ln -s %s%08d.dat ." %(vbasename,int(i))
#                 #os.system(command)
#                 print( "command:",command
#                 globname="%s*_%08d_*" %(vbasename,int(i))
#                 print( "globname:",globname

            for i in range(len(fileselection)):
                filename=fileselection[i][fileselection[i].rfind("/")+1:]
#                print( "check0:",tmp_div.split()
#                print( "check1:",len(tmp_div.split())
#                print( "check2:",len(fileselection)
                if len(tmp_div.split())==len(fileselection):
#                    print( "check3:",tmp_div.split()[i]
                    tmp_string+=" %s %s " %(filename[:-4],tmp_div.split()[i])
                else:
                    tmp_string+=" %s 1 " %(filename[:-4])

            print( "tmp_string refspec:",tmp_string)
        else:
            tmp_string="refspec "
            if find(div.get(),"1.0,x,y,z")==-1 and len(div.get())>0:
                tmp_div=div.get()
            else:
                tmp_div=""
            vbasename=basename.get()
            vdestname=destname.get()
            tmp_basename=vbasename[vbasename.rfind("/")+1:]
            tmps2="%s" %(tmp_basename)
            vsequence=str(sequence.get())
            for i in range(len(vsequence.split())):
                if len(tmp_div.split())==len(vsequence.split()):
                    tmp_string+=" %s%s %s" %(tmps2,vsequence.split()[i],tmp_div.split()[i])
                else:
                    tmp_string+=" %s%s " %(tmps2,vsequence.split()[i])

            print( "tmp_string refspec:",tmp_string)
    else:
        vbasename=basename.get()
        tmp_basename=vbasename[vbasename.rfind("/")+1:]
        tmp_string+=" -base %s" %(tmp_basename)
        vsequence=str(sequence.get())
        tmp_string+=" -seq %s" %(vsequence)

        print( "tmp_string else::",tmp_string)


    vdestname=destname.get()
    rsync="-u "

    print( "before creating Destination")
    if VERS=="NICOS":
        print( "VERS=",VERS)
        if (os.path.isdir(vdestname)==False):
            tmp="Destination (%s) does not exist\n I am going to create the path for you" %(vdestname)
            print( "->:",tmp,":<-")
            print( "Creating folders:",vdestname)
            original_umask = os.umask(0)
            os.makedirs(vdestname,0o775)
#            os.chmod(vdestname,0o775)
            os.umask(original_umask)
        else:
            tmp="Destination (%s) exist!" %(vdestname)
            print( "->:",tmp,":<-")
    else:
        if (os.path.isdir(vdestname)==False):
            tmp="Destination (%s) does not exist\n I am going to create the path for you" %(vdestname)
            print( "->:",tmp,":<-")
            print( "Creating folders:",vdestname)
            os.makedirs(vdestname)
        else:
            tmp="Destination (%s) exist!" %(vdestname)
            print( "->:",tmp,":<-")

    print( "after creating Destination")

    #Rsync data
    if VERS=="NICOS":
        for i in fileselection:
            filename="%s" %(i)
            if (os.path.isfile(filename)==False):
                tmp="Cannot find the folder (%s)!" %(filename)
                answer(tmp)
                return "",-1
            else:
                rsync+=" %s  %s_* " %(filename,filename[:-4])
        print( "after rsync:",rsync)
    else:
        vsequence=str(sequence.get())
        for i in vsequence.split(" "):
            filename="%s%s" %(vbasename,i)
            if (os.path.isdir(filename)==False):
                tmp="Cannot find the folder (%s)!" %(filename)
                answer(tmp)
                return "",-1
            else:
                rsync+=" %s " %(filename)

        print( "after rsync:",rsync)

    if vmode.get()==2:
        tmp_string="creategpl "
        if find(addstring.get(),"empty")==-1 and len(addstring.get())>0:
            tmp_string+=" "+str(addstring.get())+" "

        print( "tmp_string creategpl:",tmp_string)

    if find(simstring.get(),"empty")==-1 and len(simstring.get())>0:
#        print( "inside simstring"
        #.1 1 ; 0 0 0 ; 9.408E-6 400 5 ; 6.6E-6 0 5 ;
        vdestname=destname.get()
        fname="%s/input.str" %(vdestname)
        print( "simstring:",fname)
        fp2=open(fname,'w')
        line=simstring.get()
        print( "line:",line)
        tmp=line.split()
        for i in tmp:
            if i==";":
                fp2.write("---\n")
            else:
                fp2.write(i)
                fp2.write("\n")
        fp2.close()
        command="pwd"
        os.system(command)
        command="cd %s; rm -f curve.dat;refl_tool_silent"%(vdestname)
        os.system(command)


    print( "execute command:",tmp_string)
    message_string="copy %s to %s\nexecuting: %s \n" %(rsync,vdestname,tmp_string)
    result=message_ask(message_string)

#    return tmp_string,True
    return tmp_string,result


def execute():
    message,result=show_entry_fields()
    print( result)
    print( message)
    if result==True:
        vbasename=basename.get()
        vsequence=str(sequence.get())
        vdestname=destname.get()
        for i in vsequence.split(" "):
            if VERS=="NICOS":
                os.chdir(vdestname)
                filename="%s%08d" %(vbasename,int(i))
#                command="ln -s %s%08d.dat ." %(vbasename,int(i))
#                os.system(command)
#                print( "command:",command
                command="ln -s ../../%s%08d.dat ." %(vbasename[vbasename.rfind("/data/")+1:],int(i))
#                print( "command neu:",command
                os.system(command)
                globname="%s*_%08d_*" %(vbasename,int(i))
                for item in glob(globname):
                    #command="ln -s %s ." %(item)
                    #os.system(command)
                    #print( "command glob:",command
                    command="ln -s ../../%s ." %(item[item.rfind("/data/")+1:])
                    os.system(command)
#                    print( "command glob neu:",command
            else:
                filename="%s%s" %(vbasename,i)
                rsync="rsync -av %s/ %s" %(filename,vdestname)
                print( "execute",rsync)
                os.system(rsync)
        command="cd %s;%s" %(vdestname,message)
        print( "executing command:", command)
        print( "executing system(%s)" %(message))
        os.system(command)
        print( "finished")

    return

master = Tk()

vpng=IntVar()
vpng.set(1)
vaiaf=IntVar()
vaiaf.set(1)
vavi=IntVar()
vavi.set(0)
vread=IntVar()
vread.set(0)
vsens=IntVar()
vsens.set(1)
vsumup=IntVar()
vsumup.set(0)
vdivdet=IntVar()
vdivdet.set(0)
vsubdet=IntVar()
vsubdet.set(0)
vnoref=IntVar()
vnoref.set(0)
vsfc=IntVar()
vsfc.set(0)
vfc=IntVar()
vfc.set(0)

monitors={
    ("Time",0,1),
    ("Mon1",1,2),
    ("Mon2",2,3),
    ("Nothing",-1,4)
}

vmonitors=IntVar()
vmonitors.set(2)

mode={
    ("Reflectivity",0,0),
    ("GISANS",1,1),
    ("CreateGPL",2,2),
    ("KinImage",3,3),
    ("Refspec",4,4)
}
vmode=IntVar()
vmode.set(0)


input_row=0
Label(master, text=""" """).grid(row=input_row,column=0)
input_row+=1
Label(master, text="""basename of the folder""").grid(row=input_row,column=0)
if VERS=="NICOS":
    Button(master, text="Browse", command=open_filer_base_NICOS).grid(row=input_row,column=2)
else:
    Button(master, text="Browse", command=open_filer_base).grid(row=input_row,column=2)
Label(master, text="""sequence numbers space seperated""").grid(row=input_row,column=3)
basename = Entry(master, width=40)
sequence = Entry(master)
#basename.insert(0,"/Users/mattauch/tmp/Z34/test/align_")
if VERS=="NICOS":
    basename.insert(0,"/data/")
else:
    basename.insert(0,"~/data/")
sequence.insert(0,"0 1 2")
basename.grid(row=input_row, column=1)
sequence.grid(row=input_row, column=4)
input_row+=1
Label(master, text="""Save into (empty) folder """).grid(row=input_row,column=0)
Button(master, text="Browse", command=open_filer_dest).grid(row=input_row,column=2)
destname = Entry(master, width=40)
#destname.insert(0,"/Users/mattauch/tmp/bla/")
destname.insert(0,"~/Data_Treatment/")
destname.grid(row=input_row, column=1)


#
input_row+=1
Label(master, text="""Choose the main mode:""", justify="left").grid(row=input_row,column=0)
for txt, val, step in mode:
    Radiobutton(master, text=txt,variable=vmode,value=val).grid(row=input_row,column=1+step)

#
input_row+=1
Label(master, text="""Choose from the following options:""").grid(row=input_row,column=0)
Checkbutton(master, text="""create new png files""",variable=vpng).grid(row=input_row,column=1+0)
Checkbutton(master, text="""create an aiaf map""",variable=vaiaf).grid(row=input_row,column=1+1)
Checkbutton(master, text="""create an avi file""",variable=vavi).grid(row=input_row,column=1+2)
Checkbutton(master, text="""read data again""",variable=vread).grid(row=input_row,column=1+3)
Checkbutton(master, text="""Sensitivity correction""",variable=vsens).grid(row=input_row,column=1+4)
Label(master, text=""" """).grid(row=input_row+1,column=0)

#
input_row+=2
Label(master, text="""Choose your normalisation method:""").grid(row=input_row,column=0)
for txt, val, step in monitors:
    Radiobutton(master, text=txt, variable=vmonitors,value=val).grid(row=input_row,column=0+step)
Label(master, text=""" """).grid(row=input_row+1,column=0)

#
input_row+=2
Label(master, text="""Choose from the Reflectivity options:""").grid(row=input_row,column=0)
Checkbutton(master, text="""Simple footprint( correction""",variable=vsfc).grid(row=input_row,column=1)
Checkbutton(master, text="""Full footprint( correction\ngive sample length below""",variable=vfc).grid(row=input_row,column=2)
#Checkbutton(master, text="""Substract det images 0-1,2-1,..""",variable=vsubdet).grid(row=input_row,column=3)
#Checkbutton(master, text="""Calculate no ref-curve""",variable=vnoref).grid(row=input_row,column=4)
Label(master, text=""" """).grid(row=input_row+1,column=0)

#
input_row+=2
Label(master, text="""Choose from the GISANS options:""").grid(row=input_row,column=0)
Checkbutton(master, text="""sumup GISANS images""",variable=vsumup).grid(row=input_row,column=1)
Checkbutton(master, text="""Divide det images 0/1, 2/3,..""",variable=vdivdet).grid(row=input_row,column=2)
Checkbutton(master, text="""Substract det images 0-1,2-1,..""",variable=vsubdet).grid(row=input_row,column=3)
Checkbutton(master, text="""Calculate no ref-curve""",variable=vnoref).grid(row=input_row,column=4)
Label(master, text=""" """).grid(row=input_row+1,column=0)
#
input_row=20
Label(master, text=""" """).grid(row=input_row,column=0)
input_row+=1
Label(master, text="""roi width in pixels""").grid(row=input_row,column=0)
Label(master, text="""roi height in pixels""").grid(row=input_row,column=2)
roiw = Entry(master)
roih = Entry(master)
roiw.insert(10,"24")
roih.insert(10,"600")
roiw.grid(row=input_row, column=1)
roih.grid(row=input_row, column=3)

input_row+=1
Label(master, text="""roi offset in pixels""").grid(row=input_row,column=0)
roioff= Entry(master)
roioff.insert(10,"0")
roioff.grid(row=input_row, column=1)

input_row+=1
Label(master, text="""background correction:""").grid(row=input_row,column=0)
Label(master, text="""-1=automatic or 0=none or >0 linear\n>0: average over n-points on left and right \nof roi and calculate alinear background""").grid(row=input_row,column=2)
Label(master, text=""" """).grid(row=input_row,column=3)
bglin = Entry(master)
bglin.insert(10,"-1")
bglin.grid(row=input_row, column=1)


input_row+=1
Label(master, text=""" """).grid(row=input_row,column=0)
input_row+=1
Label(master, text="""Experimental footprint( correction""").grid(row=input_row,column=0)
Label(master, text="""Scale factor for reflectivity curve""").grid(row=input_row,column=2)
fcsl = Entry(master)
scale = Entry(master)
fcsl.insert(10,"Sample length in mm")
scale.insert(10,"1.0")
fcsl.grid(row=input_row, column=1)
scale.grid(row=input_row, column=3)

input_row+=1
Label(master, text="""Stitching of reflectivity curve:\n divide sections by a number\n(space seperated string)""").grid(row=input_row,column=0)
Label(master, text="""Set wavelength of selector\n(after a selector fault)""").grid(row=input_row,column=2)
div = Entry(master)
mansel = Entry(master)
div.insert(10,"")
mansel.insert(10,"-1")
div.grid(row=input_row, column=1)
mansel.grid(row=input_row, column=3)

input_row+=1
Label(master, text="""Scale range for plots\n default: autoscale""").grid(row=input_row,column=0)
input_row+=1
Label(master, text="""Max int. for det plots""").grid(row=input_row,column=0)
Label(master, text="""Min int. for det plots""").grid(row=input_row,column=2)
maxi= Entry(master)
mini = Entry(master)
maxi.insert(10,"autoscale")
mini.insert(10,"autoscale")
maxi.grid(row=input_row, column=1)
mini.grid(row=input_row, column=3)
input_row+=1
Label(master, text="""Max int. for sumup plots""").grid(row=input_row,column=0)
Label(master, text="""Min int. for sumup plots""").grid(row=input_row,column=2)
smaxi= Entry(master)
smini = Entry(master)
smaxi.insert(10,"autoscale")
smini.insert(10,"autoscale")
smaxi.grid(row=input_row, column=1)
smini.grid(row=input_row, column=3)
input_row+=1
Label(master, text="""Max int. for subdet plots""").grid(row=input_row,column=0)
Label(master, text="""Min int. for subdet plots""").grid(row=input_row,column=2)
submaxi= Entry(master)
submini = Entry(master)
submaxi.insert(10,"autoscale")
submini.insert(10,"autoscale")
submaxi.grid(row=input_row, column=1)
submini.grid(row=input_row, column=3)
input_row+=1
Label(master, text="""Max int. for divdet plots""").grid(row=input_row,column=0)
Label(master, text="""Min int. for divdet plots""").grid(row=input_row,column=2)
divmaxi= Entry(master)
divmini = Entry(master)
divmaxi.insert(10,"autoscale")
divmini.insert(10,"autoscale")
divmaxi.grid(row=input_row, column=1)
divmini.grid(row=input_row, column=3)
input_row+=1
Label(master, text="""Max int. for aiaf map""").grid(row=input_row,column=0)
Label(master, text="""Min int. for aiaf map""").grid(row=input_row,column=2)
aiafmaxi= Entry(master)
aiafmini = Entry(master)
aiafmaxi.insert(10,"autoscale")
aiafmini.insert(10,"autoscale")
aiafmaxi.grid(row=input_row, column=1)
aiafmini.grid(row=input_row, column=3)
input_row+=1
Label(master, text="""Any additional not here \n listed command\n e.g. -noq 'magnet' """).grid(row=input_row,column=0)
addstring= Entry(master)
addstring.insert(10,"empty")
addstring.grid(row=input_row, column=1)
input_row+=1
#
input_row+=2
Label(master, text="""Layers for simulation e.g.\n  0.1 1 ; 0 0 0; 9.408E-6 400 5 ; 6.6E-6 0 5 ;""").grid(row=input_row,column=0)
simstring= Entry(master)
simstring.insert(10,"empty")
simstring.grid(row=input_row, column=1)
input_row+=2
#
Button(master, text='Quit', command=master.quit).grid(row=40, column=0, sticky=W, pady=4)
#Button(master, text='Try', command=show_entry_fields).grid(row=40, column=1, sticky=W, pady=4)
Button(master, text='Execute', command=execute).grid(row=40, column=2, sticky=W, pady=4)

mainloop( )
