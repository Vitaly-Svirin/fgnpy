#!/usr/bin/env python
# -*- coding: cp1251 -*-
'''
обработчик семантическими алгоритмами 
для подготовки к построению 
семантрических графов по Бодякину
'''
##import cgi
##import cgitb; cgitb.enable()
##from graphs import *


##print "Content-type: text/html"
##print
##print "<!doctype html>"


import codecs
import re
from string import lower


#from stops import *
stop_symbols=''
stop_words = []



def couples(words,distance):
    'выделение пар слов с расстояниемм между ними distance'
    S=[]
    for i in range(len(words)-distance):
        S+=[(words[i],words[i+distance])]
    return S

def counter(s):
    'подсчет частот связей (пар) и формирование словаря'
    T={}
    for i in range(len(s)):
        p=s[i]; #r=list(s[i]);    tuple(r.reverse())
        a=list(p);
        a.reverse();
        q=tuple(a); 
        #print 'p=', p,'<BR>', 'q=',q, '<BR>'

        if p in T: 
            T[p]+=1
        else:
            T[p]=1

        # if q in T:
        #     T[q]+=1
        #     del T[p] - это удаление повторов имеет смысл лишь при дублирующихся рядом словах, в генетике оно важно поэтому закоментирван


    return T

def build_txt_tez(sentences,distance):
    'составление тезауруса'
    s=[]
    for sentence in sentences: 
        #print 'sentence=',sentence
        words = sentence.split()
        s+=couples(words,distance)
        #print 's=',s, '<BR>','counter=',counter(s),'<BR>'
    return counter(s)

def cleaner(string,stop_symbols,stop_words):
    'очистка строки от ненужных символов и слов из stops.py и добавление разделиятеля  * '
    rezult= ' '.join( [x for x in [y.strip(stop_symbols) for y in string.lower().split()] if x and (x not in stop_words)] )
    return rezult+' * '

def ochistka (sentences,stop_symbols,stop_words):
    'очистка срок в цикле'
    cleanded=[]
    try:
        for sentence in sentences:
            sentence= sentence.lower()
            sentence= sentence.encode('utf-8')
            cleanded+=cleaner(sentence,stop_symbols,stop_words)
            clean=",".join(cleanded).replace(',','')

        clean=clean.split(' * ')
        clean=[clean][0][:-1] #print 'clean =',clean
    except:
        print 'error decoding'
        clean='error decoding'
    return clean


def prepare(T):
    'подготова к передаче в граф исходного словаря'
    names=[]
    for i in T:
        if i[0] not in names: names.append(i[0])
        if i[1] not in names: names.append(i[1])
    names = [str(i) for i in names]

    g=[]
    for i in T:
        g.append((names.index(i[0]),names.index(i[1]))); #print g
    return g,names
    

def freq_counter(T):
    'подсчет частот концептов (вершин) и возврат словарем'
    keys=[]
    for para in T.keys():
        for i in para:
            keys.append(i) 
    #print keys
    frqs_vertexes={}
    for i in keys:
        frqs_vertexes[i] = keys.count(i)
    return  frqs_vertexes


def graph_deep(Tez, deep):
    'оставяем в тезаурусе  пары выше порога deep'
    Tez_filtered={}
    for i in Tez:
        #print i, Tez[i], deep
        if Tez[i]>deep:
            Tez_filtered[i]=Tez[i]
    return Tez_filtered


def get_data4graf(sentences,deep,distance):
    '''синтез параметров для визуализации графа
    k1 и k2 - коэффициенты размера вершин и надписей'''
    cleaned=ochistka(sentences,stop_symbols,stop_words)
    Tez=build_txt_tez(cleaned,distance)
    Tez=graph_deep(Tez, deep)
    if len(Tez) == 0:
        print 'ПОЛУЧЕН ПУСТОЙ ТЕЗАУРУС'
        return
    else:
        g,names=prepare(Tez)
        frqs_vertexes=freq_counter(Tez)
        vertex_sizes=[frqs_vertexes[i] for i in names]
        vertex_names_sizes=[frqs_vertexes[i] for i in names]
        return g,names,vertex_sizes,vertex_names_sizes,Tez

if __name__=='__main__':    
    sentences=[u'собаки: и кошки и кошки',u'собаки, мотыли, кошки игрушки',
    u'''Поступки мудрых Поступки мудрых людей продиктованы умом, людей менее сообразительных - опытом, самых невежественных - необходимостью, животных - природою.'''] 
    #sentences=['dogs-: | and cats dogs and cats','dogs motyl games','cats games'] # dlya otladki in english
    im_name='/DB/'+'semantix_example'
    
    deep = 0 # порог частот связей для визуализации
    distance=1

    g,names,vertex_sizes,vertex_names_sizes,Tez = get_data4graf(sentences,deep,distance)
    widths=Tez.values();  #print Tez.values()
    vertex_sizes=[sigma(i,0.5)*30 for i in vertex_sizes]
    vertex_names_sizes=[sigma(i,0.5)*15 for i in vertex_names_sizes]
    widths=[sigma(i,0.5)*3 for i in widths]
    #vertex_sizes=20
    #vertex_names_sizes=10
    G=graph(g,names,vertex_sizes,vertex_names_sizes,im_name,widths,(800, 600), layout='kk')
    #layout= kk  graphopt fr  drl  tree circle
    #processing_viz(Tez)