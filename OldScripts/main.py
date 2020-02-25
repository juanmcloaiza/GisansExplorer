class FrozenClass(object):
    __isfrozen = False
    def __setattr__(self, key, value):
        if self.__isfrozen and not hasattr(self, key):
            raise TypeError( "%r is a frozen class" % self )
        object.__setattr__(self, key, value)

    def _freeze(self):
        self.__isfrozen = True


class Settings(FrozenClass):
    import platform
    def __init__(self):
        self.sname2 = None
        self.debug = None
        self.channelwidth = None
        self.h_max_channel = None
        self.h_min_channel = None
        self.c = None
        self.nist_dict_pol9 = None
        self.analyzerdic = None
        self.v_min_channel = None
        self.nist_dict_pol4 = None
        self.pflipperdic = None
        self.res = None
        self.nist_dict = None
        self.aflipperdic = None
        self.detectorwidth = None
        self.gnufont = None
        self.nist_dict_pol3 = None
        self.nist_dict_pol6 = None
        self.s2_left = None
        self.s1_left = None
        self.nist_dict_pol1 = None
        self.dist_samp_det = None
        self.s1_right = None
        self.nist_dict_pol7 = None
        self.nist_dict_pol8 = None
        self.all_channels = None
        self.nist_dict_pol5 = None
        self.fs_dict = None
        self.s2_right = None
        self.v_max_channel = None
        self.nicos = None
        self.nist_dict_pol2 = None
        self.sumup_array = None
        self.bg = None
        self.l3 = None
        self._freeze()
    
        self.set_platform()
        self.set_other_settings()



    def set_platform(self) -> None:
        self.platformstr = platform.system()
        self.message = 'Platform ' + platformstr

        if self.platformstr == 'Linux':
            self.gdfontpath = "/usr/share/fonts/dejavu"
            self.gnufont = "DejaVuSans"
        elif self.platformstr == 'MacOS':
            self.gdfontpath="/usr/local/TeX/texmf-dist/fonts/truetype/public/dejavu/"
            self.gnufont = "Courier"
        else:
            print( "Platform unknown")
            self.gdfontpath = ""
            self.gnufont = "DejaVuSans"
            self.message += ": Unknown" 

        print(self.message)

        return


class Switches(FrozenClass):
    def __init__(self):
        self.base=""
        self.sequence=[]
        self.divisors=[]
        self.logz=0
        self.nopng=0
        self.noavi=0
        self.noaiaf=0
        self.aiafmm=[0,0]
        self.bglin=-1
        self.roi= [24,440]
        self.roffset=[]
        self.fread=0
        self.fpng=0
        self.favi=0
        self.faiaf=0
        self.nmon=1
        self.scale=1
        self.sfc=0
        self.fc=0
        self.mansel=0
        self.divdet=0
        self.subdet=0
        self.noref=0
        self.maxi=0
        self.mini=0
        self.smaxi=0
        self.smini=0
        self.submaxi=10
        self.submini=-10
        self.divmaxi=2
        self.divmini=.5
        self.divsmaxi=30
        self.divsmini=0
        self.sdet=1
        self.sens_det_image=""
        self.sumup=0
        self.gisans=0
        self.window=[]
        self.rev=0
        self.leak=0
        self.noq=""
        self.genx=0
        self.coh=0
        self.specfromaiaf=0
        self.vert=0
        self._freeze()

    base,sequence,divisors,logz,nopng,noavi,noaiaf,aiafmm,bglin,roi,roffset,fread,fpng,favi,faiaf,nmon,scale,sfc,fc,mansel,divdet,subdet,noref,maxi,mini,smaxi,smini,submaxi,submini,divmaxi,divmini,divsmaxi,divsmini,sdet,sens_det_image,sumup,gisans,window,rev,leak,noq,genx,coh,specfromaiaf,vert=get_switches()


def main():
    import datetime
    a = datetime.datetime.now()
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


  if sdet==1:
    sens_det_array=get_sens_det_image(sens_det_image)
  else:
    sens_det_array=[[0 for i in range(settings.all_channels)] for j in range(settings.all_channels)]

  if settings.is_nicos():
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
