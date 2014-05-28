##------------reports---------------------
import os
from Bio import SeqIO
from Bio.Seq import Seq
##import tkFileDialog as tfd
import matplotlib as mpl
mpl.use('agg')
from matplotlib import pyplot as plt
from xlwt import Workbook, XFStyle, Borders, Pattern, Font
##import dictsGf
import shutil
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text',usetex=True)

##class xlbook:
##    def __init__(self,sheetNamesList):
##        self.book = Workbook()
##        self.sheetdict={}
##        for name in sheetNamesList:
##            self.sheetdict[name]=self.book.add_sheet(name)
##        self.row = 0
##        self.column = 0
##        self.currentsheet=sheetNamesList[0]
##        
##    def write(self,data):
##        sheetname=self.currentsheet
##        row=self.row
##        column=self.column
##        self.sheetdict[sheetname].write(row,column,data)
##
##    def save_book(self,path):
##        self.book.save(path)
##
##        
##class ReportBQH:
##    def __init__(self,filename):
##        self.record=SeqIO.read(filename,"genbank")
##        self.dict=dictsGf.dicts['BQH'].copy()
##        vals=self.dict['vrange']+['%'+str(v)+'/total' for v in self.dict['vrange']]
##        self.outvalue={}.fromkeys(vals,0)
##        self.outvalue['total']=0
##        path=filename.rpartition['/'][0]+filename.rpartition['/'][2].rpartition['.gb'][0]
##        os.mkdir(path)
##        self.path=path
##        
##    def getoutval(self):
##        seq=self.record.seq
##        ln=len(seq)
##        step=self.dict['step']
##        for i in xrange(0,ln-ln%step,step):
##            self.outvalue['total']+=1
##            self.outvalue[self.dict[seq[i:i+step]]]+=1
##        for v in self.dict['vrange']:
##            self.outvalue['%'+str(v)+'/total']=round(float(self.outvalue[v])/self.outvalue['/total']*100,2)
##            
##    def xlwrite(outvalue):
##        book=xlbook([self.dict['name']])
##        for key in outvalue:
##            book.write(key)
##            book.column+=1
##        book.column=0
##        book.row+=1
##        for key in outvalue:
##            book.write(outvalue[key])
##            book.column+=1
##        book.save_book(self.path)

class xfname:
    def __init__(self,initfname,compstep,compdeep,comptype):
        self.initfname=initfname
##        print initfname
        self.path=initfname.rpartition('/')[0]+'/'
        self.name=initfname.rpartition('/')[2].rpartition('.')[0]
        self.ext='.'+initfname.rpartition('/')[2].rpartition('.')[2]
        self.compstep=compstep
        self.compdeep=compdeep
        self.comptype=comptype
    def get_name(self):
        return self.path+self.name+'-s'+str(self.compstep)+'-'+str(self.compdeep)+'-'+self.comptype+self.ext
    def get_level(self):
        if self.comptype=='':
            return str(self.compdeep)
        return str(self.compdeep)+'/'+self.comptype
    def copy(self):
        return xfname(self.initfname,self.compstep,self.compdeep,self.comptype)

def transposition(fname,compopt,sotype):
    if sotype==[1]:
        report=Report_F_BQH_P_O(fname,compopt,'')
        report.process()
        return
    path=fname.rpartition('/')[0]+'/'        
    if fname[-3:]=='.gb':
        rec=SeqIO.read(fname,"genbank")
    else:
        rec=SeqIO.read(fname,"fasta")
    
    for n in sotype:
        try:
            os.mkdir(path+str(n)+'-plet_trpos')
        except:
            pass
        shifts=[len(list(set([x/n**(n-j-1)-x/n**(n-j)*n for j in xrange(n)])))==n and
                [str(x/n**(n-j-1)-x/n**(n-j)*n) for j in xrange(n)]for x in xrange(n**n)]
        for shift in shifts:
            if shift:
                subpath=path+str(n)+'-plet_trpos/'+''.join(shift)
                try:
                    os.mkdir(subpath)
                except:
                    pass
                resseq=Seq('',rec.seq.alphabet)
                resname=subpath+'/'+fname.rpartition('/')[2]
                res=open(resname,'w')
                for i in range(0,len(rec.seq)-len(rec.seq)%n,n):
                    for p in shift:
                        resseq+=rec.seq[i:i+n][int(p)]
                rec.seq=resseq
                if fname[-3:]=='.gb':
                    SeqIO.write(rec,res,"genbank")
                else:
                    SeqIO.write(rec,res,"fasta")
                res.close()
                report=Report_F_BQH_P_O(resname,compopt,str(n)+'/'+''.join(shift))
                report.process()

class Report_F_BQH_P_O:        
    def __init__(self,filename,compopt,sotype):
        """compopt={'compdeep':compdeep,'compstep':compstep,'comptype':comptype};
        comptype=[posi,poj...] """
        self.report_name='Report_F_BQH_P_O'
        self.sotype=sotype
        ##self.record=SeqIO.read(filename,"genbank")
        self.dictBQH={'tat': -1, 'tgt': -1, 'tct': 1, 'tag': 1, 'ctt': -1, 'ttt': -1, 'tgc': 1, 'tga': 1,
                      'tgg': 1, 'tac': 1, 'ttc': 1, 'tcg': -1, 'tta': 1, 'ttg': 1, 'tcc': -1, 'tca': -1,
                      'gca': 1, 'gta': 1, 'gcc': 1, 'gtc': 1, 'gcg': 1, 'gtg': 1, 'caa': -1, 'ctg': 1,
                      'gtt': -1, 'gct': -1, 'ggt': -1, 'plot': 1, 'cga': 1, 'cgc': 1, 'gat': 1, 'aag': -1,
                      'cgg': 1, 'act': -1, 'ggg': 1, 'gga': 1, 'ggc': 1, 'gag': -1, 'aaa': -1, 'gac': -1,
                      'cgt': -1, 'gaa': -1, 'name': 'H8', 'acc': 1, 'atg': -1, 'aca': 1, 'acg': 1, 'atc': -1,
                      'aac': -1, 'ata': -1, 'agg': -1, 'cct': -1, 'agc': -1, 'step': 3, 'aga': -1, 'cat': 1,
                      'aat': 1, 'att': 1, 'vrange': [1, -1], 'cta': 1, 'ctc': 1, 'cac': -1, 'ccg': 1, 'agt': 1,
                      'cag': -1, 'cca': 1, 'ccc': 1, 'taa': 1}
        if filename.rpartition('.gb')[0]!='':
            path=filename.rpartition('/')[0]+'/'+filename.rpartition('/')[2].rpartition('.gb')[0]
            inputfile=path+'/'+filename.rpartition('/')[2].rpartition('.gb')[0]+'.gb'
        elif filename.rpartition('.fas')[0]!='':
            path=filename.rpartition('/')[0]+'/'+filename.rpartition('/')[2].rpartition('.fas')[0]
            inputfile=path+'/'+filename.rpartition('/')[2].rpartition('.fas')[0]+'.fas'
        inputfile=xfname(inputfile,compopt['compstep'],0,'')
        try:
            os.mkdir(path)
        except:
            pass
##        print inputfile.get_name()
##        print inputfile.get_name().rpartition('/')[2]
        if inputfile.get_name().rpartition('/')[2] not in os.listdir(path):
            shutil.copy2(filename,inputfile.get_name())
        self.inpufile=inputfile
        self.path=path
        self.compopt=compopt
        self.oligo=range(1,6)
        self.plotopt={}
        for oligo in range(1,6):
            self.plotopt[str(oligo)+'-plet']=1
            if oligo > 1:
                self.plotopt['qua'+str(oligo)+'-plet']=2

    def get_C(self,oligo):
        C={'a':'t','t':'a','c':'g','g':'c'}
        coligo=''
        for l in oligo:
            coligo+=C[l]
        return coligo
    
    def get_CR(self,oligo):
        C={'a':'t','t':'a','c':'g','g':'c'}
        roligo=list(oligo)
        roligo.reverse()
        croligo=''
        for l in roligo:
            croligo+=C[l]
        return croligo

    def plot1(self,oligodata,name):
        for level in oligodata:
            valcheck=[]
            x=[]
            y=[]
            maxval=0
            for val in oligodata[level]:
                if oligodata[level][val]>maxval:
                    maxval-=maxval
                    maxval+=oligodata[level][val]
                crval=self.get_CR(val)
                if val!=crval and (val not in valcheck) and (crval not in valcheck):
                    valcheck.append(val)
                    valcheck.append(crval)
                    x.append(oligodata[level][val])
                    y.append(oligodata[level][crval])
            maxval=maxval+0.05*maxval
##            print 'maxval=',maxval
            plt.plot(x,y,'k+')
            plt.plot([-1,maxval],[-1,maxval],'k')
##            plt.xlabel('x')
##            plt.ylabel('y')
            plt.text(maxval/2,maxval-2*float(maxval)/100,r'$S_{'+ level+'}^{'+self.sotype+'}$',va='top',ha='center',fontsize=20)
##            plt.legend()
            resname=self.inpufile.path+name+'_'+level.rpartition('/')[0]+'-'+level.rpartition('/')[2]+'.pdf'
            plt.savefig(resname)
            plt.cla()
        plt.close()

    def plot2(self,oligodata,name):        
        for level in oligodata:
            valcheck=['z(a)','z(c)','z(g)','z(t)','z(t+a)/z(c+g)']
            x=[]
            y=[]
            maxval=0
            for val in oligodata[level]:
                if (val not in valcheck) and (oligodata[level][val]>maxval):
                    maxval-=maxval
                    maxval+=oligodata[level][val]
                if (val not in valcheck):
                    cval=self.get_C(val[0:int(name[3])-1])
                    cval=cval+'a+'+cval+'c+'+cval+'g+'+cval+'t'
                    if (cval not in valcheck):
                        valcheck.append(val)
                        valcheck.append(cval)
                        x.append(oligodata[level][val])
                        y.append(oligodata[level][cval])
            maxval=maxval+0.05*maxval
##            print 'maxval=',maxval
            plt.plot(x,y,'k+')
            plt.plot([-1,maxval],[-1,maxval],'k')
##            plt.xlabel('x')
##            plt.ylabel('y')
            plt.text(maxval/2,maxval-2*float(maxval)/100,r'$S_{'+ level+'}^{'+self.sotype+'}$',va='top',ha='center',fontsize=20)
##            plt.legend()
            resname=self.inpufile.path+name+'_'+level.rpartition('/')[0]+'-'+level.rpartition('/')[2]+'.pdf'
            plt.savefig(resname)
            plt.cla()
        plt.close()        

    def complete_oligolist(self,oligolist,stroligo,n):
        let='acgt'
        if n==0:
            return
        if n ==1 :
            for l in let:
                oligolist.append(stroligo+l)
        else:
            for l in let:
                self.complete_oligolist(oligolist,stroligo+l,n-1)

    def get_pos_on_level(self,poslist,pos,deep):
        comptype="".join([str(c) for c in self.compopt['comptype']])
        if deep==1:
            for p in comptype:
                poslist.append(pos+p)
        else:
            for p in comptype:
                self.get_pos_on_level(poslist,pos+p,deep-1)


    def form_sheets(self,sheets,filename):
        vals=self.dictBQH['vrange']+['%'+str(v)+'/total' for v in self.dictBQH['vrange']]+['total']
        try:
            sheets['BQH']
        except:
            sheets['BQH']={}
        sheets['BQH'][filename.get_level()]={}.fromkeys(vals,0)
        letters='acgt'
        color='wb'
        vals=[]
        for c in color:
            for l in letters:
                vals.append(c+'('+l+')')
        try:
            sheets['white_and_black_letters']
        except:
            sheets['white_and_black_letters']={}
        sheets['white_and_black_letters'][filename.get_level()]={}.fromkeys(vals,0)

                    
        for oligo in self.oligo:
            vals=[]
            pos=range(oligo)
            vals.append('totp')
            for l in letters:
                vals.append(l)
                for p in pos:
                    vals.append(l+str(p)+'/'+l)
                    vals.append(l+str(p)+'/'+str(p))
            try:
                sheets['pos'+str(oligo)+'-plet']
            except:
                sheets['pos'+str(oligo)+'-plet']={}
            sheets['pos'+str(oligo)+'-plet'][filename.get_level()]={}.fromkeys(vals,0)
            vals=[]
            self.complete_oligolist(vals,'',oligo)
            try:
                sheets[str(oligo)+'-plet']
                sheets['sort'+str(oligo)+'-plet']
            except:
                sheets[str(oligo)+'-plet']={}
                sheets['sort'+str(oligo)+'-plet']={}
                
            sheets[str(oligo)+'-plet'][filename.get_level()]={}.fromkeys(vals,0)
            sheets['sort'+str(oligo)+'-plet'][filename.get_level()]={}.fromkeys(vals,0)
            if oligo>1:
                try:
                    sheets['qua'+str(oligo)+'-plet']
                except:
                    sheets['qua'+str(oligo)+'-plet']={}
                vals=[]
                voligo=[]
                self.complete_oligolist(voligo,'',oligo-1)
                for vo in voligo:
                    vals.append(vo+'a+'+vo+'c+'+vo+'g+'+vo+'t')
                vals=vals+['z(a)','z(c)','z(g)','z(t)','z(t+a)/z(c+g)']
                sheets['qua'+str(oligo)+'-plet'][filename.get_level()]={}.fromkeys(vals,0)
        
    def compress(self,filename,cd,pos):
        filename.compdeep=cd-1
        filename.comptype=pos[0:len(pos)-1]
        if filename.ext=='.gb':
            rec=SeqIO.read(filename.get_name(),"genbank")
            ln=len(rec.seq)
        else:
            rec=SeqIO.read(filename.get_name(),"fasta")
            ln=len(rec.seq)
        filename.compdeep=cd
        filename.comptype=pos
        numpos=int(pos[len(pos)-1])
        compstep=self.compopt['compstep']
        resseq=Seq('',rec.seq.alphabet)
        res=open(filename.get_name(),'w')
        oligolist=[]
        self.complete_oligolist(oligolist,'',compstep)
        for i in xrange(0,ln-ln%compstep,compstep):
            if str(rec.seq[i:i+compstep]).lower() in oligolist:
                resseq+=rec.seq[i:i+compstep][numpos]
        rec.seq=resseq
        if filename.ext=='.gb':
            SeqIO.write(rec,res,"genbank")
        else:
            SeqIO.write(rec,res,"fasta")
        res.close()
        return resseq

    def complete_sheets(self,sheets,filename,data):
        data=str(data).lower()
        ln=len(data)
        letters='acgt'
        for i in range(0,ln-ln%self.dictBQH['step'],self.dictBQH['step']):
            trip=data[i:i+self.dictBQH['step']]
            try:
                sheets['BQH'][filename.get_level()]['total']+=1
                sheets['BQH'][filename.get_level()][self.dictBQH[trip]]+=1
                if self.dictBQH[trip]==-1:
                    for l in trip:
                        sheets['white_and_black_letters'][filename.get_level()]['w('+l+')']+=1
                else:
                    for l in trip:
                        sheets['white_and_black_letters'][filename.get_level()]['b('+l+')']+=1
            except:
                pass
        wa=sheets['white_and_black_letters'][filename.get_level()]['w(a)']
        wc=sheets['white_and_black_letters'][filename.get_level()]['w(c)']
        wg=sheets['white_and_black_letters'][filename.get_level()]['w(g)']
        wt=sheets['white_and_black_letters'][filename.get_level()]['w(t)']
        ba=sheets['white_and_black_letters'][filename.get_level()]['b(a)']
        bc=sheets['white_and_black_letters'][filename.get_level()]['b(c)']
        bg=sheets['white_and_black_letters'][filename.get_level()]['b(g)']
        bt=sheets['white_and_black_letters'][filename.get_level()]['b(t)']
        if (wc+wg) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+t)/w(c+g)']=round(float(wa+wt)/(wc+wg),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+t)/w(c+g)']=str(wa+wt)+'/0'

        if (bc+bg) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['b(a+t)/b(c+g)']=round(float(ba+bt)/(bc+bg),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['b(a+t)/b(c+g)']=str(ba+bt)+'/0'

        if (ba+bt) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+t)/b(a+t)']=round(float(wa+wt)/(ba+bt),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+t)/b(a+t)']=str(wa+wt)+'/0'

        if (bc+bg) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['w(c+g)/b(c+g)']=round(float(wc+wg)/(bc+bg),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['w(c+g)/b(c+g)']=str(wc+wg)+'/0'

        if (ba+bc+bg+bt) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+c+g+t)/b(a+c+g+t)']=round(float(wa+wc+wg+wt)/
                                                                                           (ba+bc+bg+bt),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+c+g+t)/b(a+c+g+t)']=str(wa+wc+wg+wt)+'/0'

        if (bc+bg) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+t)/b(c+g)']=round(float(wa+wt)/(bc+bg),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['w(a+t)/b(c+g)']=str(wa+wt)+'/0'

        if (wc+wg) != 0:
            sheets['white_and_black_letters'][filename.get_level()]['b(a+t)/w(c+g)']=round(float(ba+bt)/(wc+wg),4)
        else:
            sheets['white_and_black_letters'][filename.get_level()]['b(a+t)/w(c+g)']=str(ba+bt)+'/0'                
                    
        for v in self.dictBQH['vrange']:
            val='%'+str(v)+'/total'
            if sheets['BQH'][filename.get_level()]['total']!=0:
                vres=float(sheets['BQH'][filename.get_level()][v])/sheets['BQH'][filename.get_level()]['total']
                sheets['BQH'][filename.get_level()][val]=round(vres,2)
            else:
                sheets['BQH'][filename.get_level()][val]=str(sheets['BQH'][filename.get_level()][v])+'/0'
                
        for step in self.oligo:
            for i in xrange(0,ln-ln%step,step):
                oligo=data[i:i+step]
                try:
                    sheets[str(step)+'-plet'][filename.get_level()][oligo]+=1
                    sheets['sort'+str(step)+'-plet'][filename.get_level()][oligo]+=1
                    sheets['pos'+str(step)+'-plet'][filename.get_level()]['totp']+=1
                    for p in xrange(step):
                        l=oligo[p]
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l]+=1
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+l]+=1
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+str(p)]+=1
                except:
                    pass
            for l in letters:
                for p in xrange(step):
                    
                    if sheets['pos'+str(step)+'-plet'][filename.get_level()]['totp']!=0:
                        pres=sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+str(p)]
                        totp=sheets['pos'+str(step)+'-plet'][filename.get_level()]['totp']
                        ptres=float(pres)/totp
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+str(p)]=round(ptres,4)
                    else:
                        pres=sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+str(p)]
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+str(p)]=str(pres)+'/0'
                        
                    if sheets['pos'+str(step)+'-plet'][filename.get_level()][l]!=0:
                        lres=sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+l]
                        totl=sheets['pos'+str(step)+'-plet'][filename.get_level()][l]
                        ltres=float(lres)/totl
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+l]=round(ltres,4)
                    else:
                        lres=sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+l]
                        sheets['pos'+str(step)+'-plet'][filename.get_level()][l+str(p)+'/'+l]=str(lres)+'/0'
            
            if step>1:
                vals=sheets[str(step)+'-plet'][filename.get_level()].keys()
                vals.sort()
                for i in xrange(0,len(vals),4):
                    v=''
                    x=0
                    for j in xrange(4):
                        x+=sheets[str(step)+'-plet'][filename.get_level()][vals[i+j]]
                        sheets['qua'+str(step)+'-plet'][filename.get_level()]['z('+vals[i+j][0]+')']+=sheets[
                            str(step)+'-plet'][filename.get_level()][vals[i+j]]
                        v+=vals[i+j]
                        if j<3:
                            v+='+'
                    sheets['qua'+str(step)+'-plet'][filename.get_level()][v]=x

                za=sheets['qua'+str(step)+'-plet'][filename.get_level()]['z(a)']
                zt=sheets['qua'+str(step)+'-plet'][filename.get_level()]['z(t)']
                zc=sheets['qua'+str(step)+'-plet'][filename.get_level()]['z(c)']
                zg=sheets['qua'+str(step)+'-plet'][filename.get_level()]['z(g)']
                z=float(za+zt)/(zc+zg)
                sheets['qua'+str(step)+'-plet'][filename.get_level()]['z(t+a)/z(c+g)']=round(z,4)                      
        
                
        

    def xlwrite(self,sheets,filename):
        book=Workbook()
        shnamelist=sheets.keys()
        shnamelist.sort()
        for shname in shnamelist:
            if shname in self.plotopt:
                if self.plotopt[shname]==1:
                    self.plot1(sheets[shname],shname)
                if self.plotopt[shname]==2:
                    self.plot2(sheets[shname],shname)
                
            sheet = book.add_sheet(shname)
            row=0
            column=0
            levellist=sheets[shname].keys()
            levellist.sort()
            for level in levellist:
                sheet.write(row,column,level)
                column+=1
                row-=row
                vallist=sheets[shname][level].keys()
                vallist.sort()
                if shname[0:4]=='sort':
                    for i in xrange(0,len(vallist),2):
                        x=vallist.pop(vallist.index(self.get_CR(vallist[i])))
                        vallist.insert(i+1,x)
                elif shname[0:3]=='qua':
                    for i in xrange(0,4**(int(shname[3])-1),2):
                        root=vallist[i][0:int(shname[3])-1]
                        croot=self.get_C(root)
                        v=croot+'a+'+croot+'c+'+croot+'g+'+croot+'t'
                        x=vallist.pop(vallist.index(v))
                        vallist.insert(i+1,x)
                    
                for val in vallist:
                    sheet.write(row,column,val)
                    row+=1
                column+=1
                row-=row
                for val in vallist:
                    sheet.write(row,column,sheets[shname][level][val])
                    row+=1
                column+=1
                row-=row
                
        filename.name+=self.report_name
        filename.comptype="".join([str(c) for c in self.compopt['comptype']])
        filename.compdeep=self.compopt['compdeep']
        filename.ext='.xls'
        book.save(filename.get_name())                  

    def process(self):
        filename=self.inpufile.copy()
        sheets={}
        self.form_sheets(sheets,filename)
        if filename.ext=='.gb':
            data=SeqIO.read(filename.get_name(),"genbank").seq
        else:
            data=SeqIO.read(filename.get_name(),"fasta").seq
        self.complete_sheets(sheets,filename,data)
        for cd in xrange(1,self.compopt['compdeep']+1):
            poslist=[]
            self.get_pos_on_level(poslist,'',cd)
            for pos in  poslist:
                data=self.compress(filename,cd,pos)
                self.form_sheets(sheets,filename)
                self.complete_sheets(sheets,filename,data)
        self.xlwrite(sheets,filename)


    
            
        
    
