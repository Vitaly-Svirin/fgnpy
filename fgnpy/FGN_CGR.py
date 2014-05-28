# -*- coding: utf-8 -*-
'''
для работы программы нужно сформировать
L - список нуклеотидов и
T - соответсвующий список цифр (направлений)

L можно задать вручную:  L=[180,120,60]*50
или сгенерировать псевдослучайно: L=synth([180,120,60,240,300,0],1000); print L

T можно задать вручную: T='0 2 2 3 2 1 0 0 3 2'
а L получить из него автоматически: L=konvert(T)

можно считывать данные из файла c помощью функций viz1 и viz2 в зависимости от того, что за файл нужно визуализовать
для этого просто раскомментируйте одну из двух строк в самом  низу пограммы с указанием имени файла для визуализации
'''


from __future__ import division
import math
import random
import pylab
import os,sys
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
import HTML

from dict_quaternion import *

class xfname:
    def __init__(self,initfname,compstep,compdeep,comptype):
        self.initfname=initfname
##        print initfname
        self.path=initfname.rpartition('/')[0]+'/'
        name=initfname.rpartition('/')[2].rpartition('.')[0]
        name=name.replace(' ','_')
        self.name=name
        self.ext='.'+initfname.rpartition('/')[2].rpartition('.')[2]
        if self.ext=='.bg':
            self.seqtype="genbank"
        elif self.ext=='.fas':
            self.seqtype="fasta"
        self.compstep=compstep
        self.compdeep=compdeep
        self.comptype=comptype
    def get_name(self):
        return self.path+self.name+'-c'+str(self.compstep)+'-'+str(self.compdeep)+'-'+self.comptype+self.ext
    def get_level(self):
        if self.comptype=='':
            return str(self.compdeep)
        return str(self.compdeep)+'/'+self.comptype
    def copy(self):
        return xfname(self.initfname,self.compstep,self.compdeep,self.comptype)

def complete_oligolist(oligolist,stroligo,n):
    let='acgt'
    if n==0:
        return
    if n ==1 :
        for l in let:
            oligolist.append(stroligo+l)
    else:
        for l in let:
            complete_oligolist(oligolist,stroligo+l,n-1)

def get_pos_on_level(poslist,pos,deep,compopt):
    comptype="".join([str(c) for c in compopt['comptype']])
    if deep==1:
        for p in comptype:
            poslist.append(pos+p)
    else:
        for p in comptype:
            get_pos_on_level(poslist,pos+p,deep-1,compopt)
            
    
def compress(filename,cd,pos,compopt):
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
    compstep=compopt['compstep']
    resseq=Seq('',rec.seq.alphabet)
    res=open(filename.get_name(),'w')
    oligolist=[]
    complete_oligolist(oligolist,'',compstep)
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

def form_linklist(linklist,filename):
    fname=filename.copy()
    fname.ext='.html'
    linklist.append(
                    {'deep':fname.compdeep,'type':fname.comptype,
                     'link':HTML.link(fname.get_level(),'report/'+fname.get_name().rpartition('/')[2])
                     }
                    )


def rad(a):
    'перевод градусов в радианы'
    return (a*math.pi)/180 #print rad(568568)

def f(gradus):
    'перевод из угла в градусах в координаты'
    x=math.cos(rad(gradus))
    y=math.sin(rad(gradus))
    return x,y 

def synth(kinds,length):
    '''синтез случайной последовательности  
     состоящей из length элементов, перечисленных в списке kinds'''
    L=[]
    for j in range(length):
        for i in range(len(kinds)):
            L.append(random.choice(kinds))    
    return L
    
def savefig(name,window,color):
    'сохранение графики в файл'
    Type={0:'DIAGRAM',1:'ROAD',2:'DIAGRAM',3:'ROAD',4:'FRAKTAL',5:'FRAKTAL'} # номера окон и типы графиков
    Algebra={'g':'Q','r':'BQ','b':''} # цвета графиков и типы алгебры 
    try:          os.chdir(DIR_OUT)
    except:       pass
    #pylab.show() # если хотим увидеть картинки помимо их сохранения
    pylab.plt.savefig(name, format='png', transparent=False)
    pylab.close(window)

def graph(L,color,window,name,dictype,level,prefix):
    'вывод графики'
    pylab.figure(window)
    legend='dict. type: '+dictype+'\n'+'conv. level: '+level+'\n'+'order: '+prefix+'\n'
    #ax = pylab.subplot(120)
    X=[0]
    Y=[0]
    for i in L:
        X+=[X[-1]+f(i)[0]]
        Y+=[Y[-1]+f(i)[1]]
    #print X,Y
    Fig=pylab.gcf()
    DefaultSize = Fig.get_size_inches()
    ax=pylab.subplot2grid((5, 1), (1, 0),colspan=4, rowspan=4)
    ax.plot(X,Y,'-'+color+'.',label=legend)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
              ncol=2, mode="expand", borderaxespad=0.)
    Fig.set_size_inches( (DefaultSize[0], DefaultSize[1]/0.8) )
    savefig(name,window,color)
    pylab.close(window)
    
def drawbars (T,frqs,Delta,r,color,window,name,dictype,level,prefix):
    'вывод баров с частотами'
    pylab.figure(window)
    x = pylab.arange(4)

    legend=('dict. type: '+dictype+'\n'+'conv. level: '+level+'\n'+'order: '+prefix+'\n'+
            ' 0:'+str(T.count('0'))+str(Delta[0])+r[0]+';'+
            ' 1:'+str(T.count('1'))+str(Delta[1])+r[1]+';'+'\n'+
            ' 2:'+str(T.count('2'))+str(Delta[2])+r[2]+';'+
            ' 3:'+str(T.count('3'))+str(Delta[3])+r[3]+'.'
            )
    #ax = pylab.subplot(121)
    Fig=pylab.gcf()
    DefaultSize = Fig.get_size_inches()
    ax=pylab.subplot2grid((5, 1), (1, 0),colspan=4, rowspan=4)
    ax.bar(x, frqs, color=color,label=legend)
    pylab.xticks( x + 0.5,  ('0', '1', '2', '3') )
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
              ncol=2, mode="expand", borderaxespad=0.)
    Fig.set_size_inches( (DefaultSize[0], DefaultSize[1]/0.8) )
    savefig(name,window,color)
    
def serpinsky(T,color,window,name,dictype,level,prefix):
    'вывод фрактала'
    legend='dict. type: '+dictype+'\n'+'conv. level: '+level+'\n'+'order: '+prefix+'\n'
    pylab.figure(window)
    x=[100];y=[100];
    for i in T:
        if str(i)=='0': x.append(0.5*x[-1]); y.append(0.5*y[-1]);
        if str(i)=='1': x.append(0.5*x[-1]); y.append(0.5*y[-1]+50);
        if str(i)=='2': x.append(0.5*x[-1]+50);  y.append(0.5*y[-1]);
        if str(i)=='3': x.append(0.5*x[-1]+50);  y.append(0.5*y[-1]+50);
    Fig=pylab.gcf()
    DefaultSize = Fig.get_size_inches()
    ax=pylab.subplot2grid((5, 1), (1, 0),colspan=4, rowspan=4)
    ax.plot(x,y,color+'+',label=legend)
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
       ncol=2, mode="expand", borderaxespad=0.)
    Fig.set_size_inches( (DefaultSize[0], DefaultSize[1]/0.8) )
    savefig(name,window,color)
    
def konvert(T):
    'конвертор из символов в углы c фильтрацией'
    S={'0':180,
       '1':120,
       '2':60,
       '3':0}  
    r=[]
    for i in T:
        try:
            r.append(S[i])
        except:
            pass
    return r

def splitCount(s, count):
        'разделение строки на n-плеты n это - count'
        return [''.join(x) for x in zip(*[list(s[z::count]) for z in range(count)])]

def nplets(s,n):
    'фильтрация строки по 4 буквам и ее конвертация в список триплетов'
    t=''
    for i in s:
        if i in ['a','g','c','t','u','A','G','C','T','`U']:
            t+=i
        else:
            pass
        t.split
    return splitCount(t,n)

def triplets(s):
    return nplets(s,3)

def perestanovki(s):
    '''Всего из исходной секвенции с порядком позиций в 
    триплетах 1-2-3 можно сконструировать пять секвенций 
    с таким измененным порядком позиций в триплетах: 
    2-3-1, 3-1-2, 3-2-1, 2-1-3, 1-3-2.'''

    s123=s
    s231=[]
    s312=[]
    s321=[]
    s213=[]
    s132=[]
    for triplet in s:
        s231.append(triplet[1]+triplet[2]+triplet[0])
        s312.append(triplet[2]+triplet[0]+triplet[1])
        s321.append(triplet[2]+triplet[1]+triplet[0])
        s213.append(triplet[1]+triplet[0]+triplet[2])
        s132.append(triplet[0]+triplet[2]+triplet[1])

    return s123,s231,s312,s321,s213,s132



def delta(B,Q):
    ''' B – количество цифры в бикватернионной последовательности,
      Q – ее количество в кватернионной последовательности. 
      Например, в таблице в первой геносеквенции для цифры 0 указано,
       что в кватернионной последовательности ее количество Q равно 
       7297, а в бикватернионной  последовательности ее количество 
       В равно 7261. 
       Применяя названную формулу, получаем ∆=100*(7261/7297 – 1) = - 0,49%. 
    '''
    d=[]
    for i in range(4): d.append (round(100*(B[i]/Q[i] - 1),2))
    return d

def rel(x0,x1,x2,x3):
    ''' сколько процентов от общего количества этих цифр составляет количество каждой 
    цифры  Например, если секвенция содержит 20 штук цифры 0, 
    50 штук цифры 1, 40 штук цифры 2 и 90 штук цифры 3 (общее их количество равно 200), 
    то процентное содержание этих цифр будет 10%, 25%, 20% и 45%.
    '''
    s=(x0+x1+x2+x3)
    return ('('+str(round(100*x0/s,2))+'%)',
           '('+str(round(100*x1/s,2))+'%)', 
           '('+str(round(100*x2/s,2))+'%)', 
           '('+str(round(100*x3/s,2))+'%)')

def viz_number(filename):
    '''визуализация файла с последовательностью вида 1 2 3 0
    пример использования #viz_number('seq_numbers.txt')'''
    T=readfile(filename) 
    L=konvert(T);
    drawbars(T,'b',0,filename)
    graph(L,'b',1,filename)

def counter(data,inputfile):
    '''запись отчета в файл'''

    def count_letters_and_write_to_file(seq,inputfile):
        '''запись 1 отчета в файл'''
        seq=seq.lower()
##        print seq,len(seq)
        c=seq.count('c');g=seq.count('g');a=seq.count('a');t=seq.count('t') # подсчет букв
        S=float(c+g+a+t)
##        print S
        CG=float(c+g)
        AT=float(a+t)
        outfile=inputfile.copy()
        outfile.ext='.html'
        outf=open(outfile.get_name(),'w')
        def p(l):
            'проценты'
            return ' ('+str(round(l*100/S,2))+' %)'

        print >> outf, HTML.link('Home','../report.html'),'<br><hr>',outfile.name,'<br><hr>' 
        print >> outf, 'conv. level:',outfile.get_level(),'<br><hr>'
##        print >> outf, 'total amount of all letters:',S,'<br>' 
##        print >> outf, 'number of letters c:',c,p(c),'<br>'
##        print >> outf, 'number of letters g:',g,p(g),'<br>'
##        print >> outf, 'number of letters a:',a,p(a),'<br>'
##        print >> outf, 'number of letters t:',t,p(t),'<br>'
##        print >> outf, 'total number of letters  (c+g):', CG, p(CG),'<br>'
##        print >> outf, 'total number of letters  (a+t):', AT, p(AT),'<br>'
##        print >> outf, 'relevant amounts of pairs of letters: (a+t)/(c+g):', round(AT/CG,3),'<br>'   
##        print >> outf, 'pairwise relationship of letters: a/c:', round(a/c,3),'<br>'
##        print >> outf, 'pairwise relationship of letters: a/g:', round(a/g,3),'<br>'   
##        print >> outf, 'pairwise relationship of letters: a/t:', round(a/t,3),'<br>'   
##        print >> outf, 'pairwise relationship of letters: g/c:', round(g/c,3),'<br>'    
##        print >> outf, '<br><hr>'
        return outf
##        outf.close()

##    seq=SeqIO.read(inputfile.get_name(),inputfile.seqtype) # print 'всего: ',len(triplets(seq))
##    s=triplets(data)
##    print s,len(data)
    outf=count_letters_and_write_to_file(data,inputfile)
    viz_symbol(data,inputfile,outf)

def viz_symbol(seq,inputfile,outf):
    '''визуализация файла с последовательностью вида agct
    пример использования viz_symbol('Arabidopsis thaliana 200001.txt')'''

    seq=seq.lower() # print 'всего: ',len(triplets(seq))
    s=triplets(seq)
    s123,s231,s312,s321,s213,s132=perestanovki(s)

    def go(s,prefix,inputfile,outf):
        figname=inputfile.copy()
        


        # зеленые - кватернионы
        TQ=''
        for i in s:
            i=i.lower()
            if   i in SQ_0: TQ+='0 '
            elif i in SQ_1: TQ+='1 '
            elif i in SQ_2: TQ+='2 '
            elif i in SQ_3: TQ+='3 '
            else: print i, 'error'
        LQ=konvert(TQ);
        Q = [TQ.count('0'), TQ.count('1'), TQ.count('2'), TQ.count('3')]

##        print 'quaternions:',  Q, 'all: ', sum(Q)

        Delta=[' ',' ',' ',' ']
        r=rel(Q[0],Q[1],Q[2],Q[3]); #print r

##        figname.ext='_qbar_Tr-'+prefix+'.png'
##        name=unicode(figname.get_name())
##        title='dict. type: '+dictype+'\n'+'conv. level: '+level+'\n'+'order: '+prefix+'\n'
##        print >> outf,'<img src=',name,'>','<br>','<hr>'
##        drawbars(TQ,Q,Delta,r,'g',0,name,"quaternions",figname.get_level(),prefix) ;
##        figname.ext='_qpath_Tr-'+prefix+'.png'
##        name=unicode(figname.get_name())
##        print >> outf,'<img src=',name,'>','<br>','<hr>'
##        #print 'drawing quaternions path ...'; 
##        graph(LQ,'g',1,name,"quaternions",figname.get_level(),prefix)  ;
        figname.ext='_qfractal_Tr-'+prefix+'.png'
        namef=figname.get_name()
        namehtm=figname.get_name().rpartition('/')[2]
        print >> outf,'<img src=',namehtm,'>','<br>','<hr>'
        #print 'drawing quaternions fractall ...'; 
        serpinsky(TQ,'g',4,namef,"quaternions",figname.get_level(),prefix);
       
        
        
        



        # красные - бикватернионы
        TBQ=''
        for i in s:
            i=i.lower()
            if   i in SBQ_0: TBQ+='0 '
            elif i in SBQ_1: TBQ+='1 '
            elif i in SBQ_2: TBQ+='2 '
            elif i in SBQ_3: TBQ+='3 '
            else: print i, 'error'
        LBQ=konvert(TBQ);
        B = [TBQ.count('0'), TBQ.count('1'), TBQ.count('2'), TBQ.count('3')]
        
        D=[i for i in delta(B,Q)]
        Delta=['('+ r'$\Delta$'+'=' +str(i)+'%)' for i in D]  #print 'delta=', Delta
        r=['','','',''] #rel(B[0],B[1],B[2],B[3]); print r

##        print 'biquaternions:',B, 'all: ', sum(B)
##        figname.ext='_bqbar_Tr-'+prefix+'.png'
##        name=unicode(figname.get_name())
##        print >> outf,'<img src=',name,'>','<br>','<hr>'
##        drawbars(TBQ,B,Delta,r,'r',2,name,"biquaternions",figname.get_level(),prefix);
##        figname.ext='_bqpath_Tr-'+prefix+'.png'
##        name=unicode(figname.get_name())
##        print >> outf,'<img src=',name,'>','<br>','<hr>'
##        #print 'drawing biquaternions path ...';
##        graph(LBQ,'r',3,name,"biquaternions",figname.get_level(),prefix)  ;
        figname.ext='_bqfractal_Tr-'+prefix+'.png'
        namef=figname.get_name()
        namehtm=figname.get_name().rpartition('/')[2]
        print >> outf,'<img src=',namehtm,'>','<br>','<hr>'
        #print 'drawing biquaternions fractall ...';
        serpinsky(TBQ,'r',5,namef,"biquaternions",figname.get_level(),prefix);
        
    go(s123,'123',inputfile,outf)   
    go(s231,'231',inputfile,outf) 
    go(s312,'312',inputfile,outf) 
    go(s321,'321',inputfile,outf) 
    go(s213,'213',inputfile,outf) 
    go(s132,'132',inputfile,outf)
    outf.close()


    
def form_report(path,linklist,compopt):
    frep=open(path+'/report.html','w')
    reptab=[['deep','level']]
    for i in xrange(compopt["compdeep"]+1):
        reptab.append([str(i)])
    for l in linklist:
        reptab[l['deep']+1].append(l['link'])
    frep.write(HTML.table(reptab))
    frep.close()


        

def main(filename,compopt={'compdeep':3,'compstep':3,'comptype':[0,1,2]}):
    """
    compopt={'compdeep':int compdeep,'compstep':int compstep,'comptype':int comptype};
    comptype=[int posi,int posj...]
    """
    def takt(data,inputfile,linklist):
        form_linklist(linklist,inputfile)
        counter(data,inputfile)
       

      
    if filename.rpartition('.gb')[0]!='':
        path=filename.rpartition('/')[0]+'/'+filename.rpartition('/')[2].rpartition('.gb')[0]+'_CGR'
        reppath=path+'/report'
        inputfile=reppath+'/'+filename.rpartition('/')[2].rpartition('.gb')[0]+'.gb'
    elif filename.rpartition('.fas')[0]!='':
        path=filename.rpartition('/')[0]+'/'+filename.rpartition('/')[2].rpartition('.fas')[0]+'_CGR'
        reppath=path+'/report'
        inputfile=reppath+'/'+filename.rpartition('/')[2].rpartition('.fas')[0]+'.fas'
    inputfile=xfname(inputfile,compopt['compstep'],0,'')
    try:
        os.mkdir(path)
    except:
        pass
    
    try:
        os.mkdir(reppath)
    except:
        pass
##        print inputfile.get_name()
##        print inputfile.get_name().rpartition('/')[2]
    

    linklist=[]
    
    if inputfile.get_name().rpartition('/')[2] not in os.listdir(reppath):
        shutil.copy2(filename,inputfile.get_name())

    if inputfile.ext=='.gb':
##        print inputfile.get_name()
        data=SeqIO.read(inputfile.get_name(),"genbank").seq
    elif inputfile.ext=='.fas':
        data=SeqIO.read(inputfile.get_name(),"fasta").seq
##    print oligox
    
    takt(data,inputfile,linklist)        
    for cd in xrange(1,compopt['compdeep']+1):
        poslist=[]
        get_pos_on_level(poslist,'',cd,compopt)
        for pos in  poslist:
            data=compress(inputfile,cd,pos,compopt)
            takt(data,inputfile,linklist)
    form_report(path,linklist,compopt)




##length_of_random_seq=100 # задание длины случайной последовательности
##print 'случайная последовательность: ', ''.join(synth(['a','g','c','t'],length_of_random_seq))
##
##DIR_IN='IN-300000vyrez' # папка с входами
##DIR_OUT='OUT-300000vyrez' # папка для сохранения результатов

