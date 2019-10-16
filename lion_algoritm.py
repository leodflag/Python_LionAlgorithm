# -*- coding: utf-8 -*-
"""
@author: leodflag
Genetic Algorithm
"""
import random,math
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
# 獅群參數設計
string_length = 8 # 基因長度
pop_lion_num=6 # 母群體數量
best_male_lion=[] # 最佳公獅
best_female_lion=[] # 最佳母獅
itera=50  # 迭代次數
mutation_rate=0.1  # 突變率
male_group=[] # 公獅群母體
female_group=[] # 母獅群母體
male_cubs_group=[] # 公幼獅群
female_cubs_group=[] # 母幼獅群
cubs_group=[] # 交配後的幼獅群
adept=[] #適應值
# 基因範圍判斷
lion_gene1_U=255
lion_gene1_L=0
lion_gene2_U=149
lion_gene2_L=120
# 繪圖參數
pltX=[] # x軸為迭代次數
pltY=[] # y軸為平均適應值
# ---------轉換-----------
# 實際對應數值換算成十進位
def floTurnTen(x):
    xx=(x-(-3.0))*255/15.1
    return round(xx)
# 十進位轉換成字串基因碼
def turnStrGene(x):
    return bin(x)[2:].zfill(string_length)
# 二進位的字串轉回十進位
def bin_Int(x):
    return int(x,2)
# 十進位換算成真實浮點數
def tenTurnflo(x):
    xx=15.1*x/255-3.0
    return round(float(xx),1)
#---------基因的操作-----------
# 指定位置突變
def _invert_at(s, index):
    # "**"會先運算，s為字串
    return bin(int(s,2)^2 ** (index))[2:].zfill(string_length)
# 大到小排序
def order(x):
    temp=0
    for i in range(6):
        for j in range(6):
            if x[j]<x[i]:
                temp=x[j]
                x[j]=x[i]
                x[i]=temp
    return x
"""
# 解壓縮   
def zipReturn(pool):
    x=[]
    adept[:],x[:]=zip(*pool)
    print('adept,x',adept,x)
    print('x',x)
    gene_1[:],gene_2[:]=zip(*x)
    print('gene_1,gene_2',gene_1,gene_2)
    return gene_1,gene_2
"""
# 判斷大小，輸入int
def range_gene_1(x):
    if x>=lion_gene1_L and x<=lion_gene1_U:
        return 1
    else:
        return 0
def range_gene_2(x):
    if x>=lion_gene2_L and x<=lion_gene2_U:
        return 1
    else:
        return 0    
# 刪掉前6個基因
def delList(x):
    if itera!=50:
       del x[0:6]     
# 產生初代獅群
def init_lion_gene(gene_1,gene_2,pop_gene_num):
    for i in range(pop_gene_num*2):
        s1=random.randint(0,255)
        if (i%2) == 0:
            gene_1.append(turnStrGene(s1))
        else:
            gene_2.append(turnStrGene(s1))
    print('male lion group',gene_1)
    print('female lion group',gene_2)
    return gene_1,gene_2 
# 交配：含隨機點交配，基因突變，k_meas分群
def lion_mating(gene_1,gene_2):
    # 隨機單點交配，公母交配產生幼獅
    i=0
    cut=0
    gene1=gene_1 
    gene2=gene_2
    while(cut < pop_lion_num):
        dot=random.randint(0,7) # 隨機取交換點的位置
        gene1[cut]=gene_1[cut][0:dot]+gene_2[cut][dot:8] # 公獅1的前段基因跟母獅1的後段基因結合
        gene1[cut+1]=gene_1[cut+1][0:dot]+gene_2[cut+1][dot:8] # 公獅2的前段基因跟母獅2的後段基因結合
        gene2[cut]=gene_2[cut][0:dot]+gene_1[cut][dot:8] # 母獅1的前段基因跟公獅1的後段基因結合
        gene2[cut+1]=gene_2[cut+1][0:dot]+gene_1[cut+1][dot:8] # 母獅1的前段基因跟公獅1的後段基因結合
        cut+=2

    #隨機單點突變，突變部分基因
    s1=gene_1+gene_2 # 幼獅群
    s2=s1 # 2n隻幼獅複製產生新2n隻幼獅
    for i in range(pop_lion_num*2):
        rand = random.random()
        if rand <= mutation_rate: # 若隨機數小於突變率
            h=random.randint(0,7)#隨機位置突變，輸入為2位元，輸出字串
            s2[i]=_invert_at(s2[i],h)#將基因碼突變
    #np.array(s1,)
    s1=s1+s2 # 突變的也抓回去
   # for i in range(pop_lion_num*2):
        
    
    return s1
#---------適應函數Adaptation function---------
def Adaptation_x1(x1):  
    f1=x1*math.sin(4*math.pi*x1)
    return f1
def Adaptation_x2(x2):  
    f2=x2*math.sin(20*math.pi*x2)
    return f2
#---------K-means function---------
def kmeans_clusters(L):
    X = np.array(L)
    plt.scatter(X[:,0],X[:,1],s=50)
    
#---------1.產生獅群---------
male_group,female_group=init_lion_gene(male_group,female_group,pop_lion_num)
#---------2.繁衍後代---------
cubs_group = lion_mating(male_group,female_group)
print('cubs',cubs_group)