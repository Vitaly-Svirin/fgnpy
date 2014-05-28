#!/usr/bin/env python
# -*- coding: cp1251 -*-



from semantix import *
import shutil
import os
from Bio import SeqIO
from Bio.Seq import Seq
##from xlwt import Workbook, XFStyle, Borders, Pattern, Font
import HTML

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
        return str(self.compdeep)+'-'+self.comptype
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

def get_Tez(seq,n_nplets,deep=0,distance=1):
    seq=str(seq)
    seq.lower()
	# deep = 0 # порог частот связей для визуализации
	# distance = 1 # 1 - смортим соседей, 2 - через 1 и тд
	# n_nplets = 2 # 3- смотрим триплеты
    w_pic=1200 # ширина картинки в пикселах
    h_pic=800 # высота картинки в пикселах
	#algorythm_layout = kk # graphopt fr  drl  tree circle 
##    print n_nplets,type(n_nplets)
    sentences, sentences_rev, rest = splitgen(seq,n_nplets)
    g,names,vertex_sizes,vertex_names_sizes,Tez_rev = get_data4graf(sentences_rev,deep,distance) 
    g,names,vertex_sizes,vertex_names_sizes,Tez     = get_data4graf(sentences,    deep,distance)

	 	
	##ramka,signs=get_matrix(Tez,Tez_rev)	#print seq
    return Tez,Tez_rev
	

def get_matrix(Tez,Tez_rev):
	'''Создание матрицы смежности'''
	##print HTML.table(G.get_adjacency())
	#print Tez, BR()
	'создание и сортировка по алфавиту списка всех элементов signs проанализированного текста'
	signs=[]
	for i in Tez.keys():
		signs.append(i[0]); signs.append(i[1])
	signs = list(set(signs)); signs.sort()

	#for i,j in Tez.keys(): print i,j, Tez[i,j]
	'создание матрицы смежности с полями и столбцами из signs'
	matrix=[]
	for i in range(len(signs)):
		x = signs[i]
		row=[]
		for j in range(len(signs)):
			y = signs [j]					
			if (x,y) in Tez.keys():
					uper=Tez[x,y]
			else: 
					uper='0'
			if (x,y) in Tez_rev.keys():
					down=Tez_rev[x,y]
			else: 
					down='0'

			if uper=='0' and down=="0":
					row.append('')
			else:
					row.append(str(uper)+'/'+str(down))
		matrix.append(row)


	'матрица с круговой обводкой signs в виде списка списков'
	ramka=[]	
	ramka.append(['']+signs+[''])
	for i in range(len(signs)):
		ramka.append([signs[i]]+matrix[i]+[]+[signs[i]])
	ramka.append(['']+signs+[''])

	return ramka,signs





def nplets(s,n):
	'фильтрация строки по 4 буквам и ее конвертация в список n-плетов по числу n (3 для триплетов)'
	def splitCount(s, count):
		'разделение строки на n-плеты n это - count'
		# for z in range(count):
		# 	print [list(s[z::count])]
		return [''.join(x) for x in zip(*[list(s[z::count]) for z in range(count)])]
	t=''
	for i in s:
		if i in ['a','g','c','t','u','A','G','C','T','U']:
			t+=i
        else:
            pass
        t.split
	return splitCount(t,n)


def splitgen(seq,n):
	'если нужно разбить ген на n-плеты и представить к графу, то используем это'
	l=nplets(seq,n) 
	rev= l [::-1]

	def splitter(l):
		'конвертирует списко в строку и добавляет символ разделения - запятую с пробелами'
		sentences=''
		for i in l:
			sentences += str(i)+' , '
		sentences = sentences[:-3] # удаление лишней запятой в конце
		return [sentences]

	sentences = splitter(l)
	# print l
	# print rev
	sentences_rev = splitter(rev)
	# print '1:',sentences, '<BR>'
	# print '2:',sentences_rev
	rest='000000'
	return sentences, sentences_rev, rest

def form_linklist(linklist,fname,oligo=0):
    if oligo==0:
        linklist.append(
                        {'linktype':'data','deep':fname.compdeep,'type':fname.comptype,
                         'link':HTML.link(fname.get_level(),fname.get_name())
                         }
                        )
    else:
        linklist.append(
                        {'linktype':'table','deep':fname.compdeep,'type':fname.comptype,
                         'link':HTML.link(fname.get_level(),fname.path+str(oligo)+'-plet_'+fname.get_level()+'.html'),
                         'oligo':oligo
                         }
                        )

def html_write(info,data,oligo,inputfile,linklist):
    tabinfo=HTML.table(info)
    tab=HTML.table(data)
    linkhome=HTML.link('Home','../report.html')
    htmlcode=('<haed><link rel="stylesheet" type="text/css" href="http://matrix.petoukhov.com/css/default.css"'
              + 'media="screen" /></haed>'+linkhome+'<hr><br>'+ tabinfo+'<hr><br>'+tab
              )
    fname=inputfile.path+str(oligo)+'-plet_'+inputfile.get_level()+'.html'
    fres=open(fname,'w')
    fres.write(htmlcode)
    fres.close()
    form_linklist(linklist,inputfile,oligo)

def form_report(path,linklist,oligolist,compopt):
##    linktab=HTML.link('Tables','report/tables.html')
    linkdata=HTML.link('datafiles','report/data.html')
    linkhome=HTML.link('Home','../report.html')
    htmlcode=('<haed><link rel="stylesheet" type="text/css" href="http://matrix.petoukhov.com/css/default.css"'
            +'media="screen" /></haed><br>'+linkdata+'<hr><br>'
            )    
    data=('<haed><link rel="stylesheet" type="text/css" href="http://matrix.petoukhov.com/css/default.css"'
        +'media="screen" /></haed>'
        )
    dtab=[]
    for d in xrange(compopt['compdeep']+1):
            dtab.append(['deep='+str(d)])
    
    tables=''
    oltab={}.fromkeys(oligolist)
    for olt in oltab:
        l=[]
        l.append([str(olt)+'-plet'])
        for d in xrange(compopt['compdeep']+1):
            l.append(['deep='+str(d)])
        oltab[olt]=l
        
    
    for dlink in linklist:
        if dlink['linktype']=='data':
            dtab[dlink['deep']].append(dlink['link'])
            
        else:
            oltab[dlink['oligo']][dlink['deep']+1].append(dlink['link'])
    for ot in oltab:
        tables+=HTML.table(oltab[ot])
        tables+='<hr><br>'
    htmlcode+=tables
    frep=open(path+'/report.html','w')
    frep.write(htmlcode)
    frep.close()
    data+=linkhome+'<hr><br>'+HTML.table(dtab)
    fdata=open(path+'/report/data.html','w')
    fdata.write(data)
    fdata.close()



def process(filename,compopt={'compdeep':3,'compstep':3,'comptype':[0,1,2]},oligox=[1,2,3]):
    """
    compopt={'compdeep':compdeep,'compstep':compstep,'comptype':comptype};
    comptype=[posi,poj...]
    """
      
    if filename.rpartition('.gb')[0]!='':
        path=filename.rpartition('/')[0]+'/'+filename.rpartition('/')[2].rpartition('.gb')[0]+'_semantix'
        reppath=path+'/report'
        inputfile=reppath+'/'+filename.rpartition('/')[2].rpartition('.gb')[0]+'.gb'
    elif filename.rpartition('.fas')[0]!='':
        path=filename.rpartition('/')[0]+'/'+filename.rpartition('/')[2].rpartition('.fas')[0]+'_semantix'
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
    else:
        data=SeqIO.read(inputfile.get_name(),"fasta").seq
##    print oligox
    form_linklist(linklist,inputfile)
    for o in oligox:
##        print o
        Tez,Tez_rev=get_Tez(data,o)
        ramka,signs=get_matrix(Tez,Tez_rev)
        infotab=[['oligonucleotide:',str(o)+'-plet'],['compress level:',inputfile.get_level()]]
        html_write(infotab,ramka,o,inputfile,linklist)
        
    for cd in xrange(1,compopt['compdeep']+1):
        poslist=[]
        get_pos_on_level(poslist,'',cd,compopt)
        for pos in  poslist:
            data=compress(inputfile,cd,pos,compopt)
            form_linklist(linklist,inputfile)
            for o in oligox:
##                print o
                Tez,Tez_rev=get_Tez(data,o)
                ramka,signs=get_matrix(Tez,Tez_rev)
                infotab=[['oligonucleotide:',str(o)+'-plet'],['compress level:',inputfile.get_level()]]
                html_write(infotab,ramka,o,inputfile,linklist)
    form_report(path,linklist,oligox,compopt)
                







