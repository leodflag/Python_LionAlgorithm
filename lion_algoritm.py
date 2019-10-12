# -*- coding: utf-8 -*-
"""
@author: leodflag
Genetic Algorithm
"""
import random,math
import matplotlib.pyplot as plt

string_length = 8 # 基因長度
pop_lion_num=6 # 母群體數量
best_male_lion=[] # 最佳公獅
best_female_lion=[] # 最佳母獅
itera=50  # 迭代次數
mutation_rate=0.1  # 突變率

# 獅群參數設計
male_group=[] # 公獅群母體
female_group=[] # 母獅群母體
male_cubs_group=[] # 公幼獅群
female_cubs_group=[] # 母幼獅群
cubs_group=[] # 交配後的幼獅群

# 基因範圍判斷
lion_gene1_U=255
lion_gene1_L=0
lion_gene2_U=149
lion_gene2_L=120

# 繪圖參數
pltX=[] # x軸為迭代次數
pltY=[] # y軸為平均適應值
"""
轉換
"""
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
"""
基因的操作
"""
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
# 解壓縮   
def zipReturn(pool):
    x=[]
    adept[:],x[:]=zip(*pool)
    print('adept,x',adept,x)
    print('x',x)
    gene_1[:],gene_2[:]=zip(*x)
    print('gene_1,gene_2',gene_1,gene_2)
    return gene_1,gene_2
# 判斷大小，輸入int
def range_gene_1(x):
    if x>=gene_1_L and x<=gene_1_U:
        return 1
    else:
        return 0
def range_gene_2(x):
    if x>=gene_2_L and x<=gene_2_U:
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
"""
適應函數Adaptation function
"""
def Adaptation_x1(x1):  
    f1=x1*math.sin(4*math.pi*x1)
    return f1
def Adaptation_x2(x2):  
    f2=x2*math.sin(20*math.pi*x2)
    return f2

"""
#初代基因組
x=list(zip(gene_1,gene_2))
print('初代基因組',x)
print('初代基因組排序',order(x))
"""
"""
1.產生獅群
"""
male_group,female_group=init_lion_gene(male_group,female_group,pop_lion_num)

"""
2.繁衍後代
"""