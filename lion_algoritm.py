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

#繪圖參數
pltX=[]#x軸為迭代次數
pltY=[]#y軸為平均適應值
"""
轉換
"""
#實際對應數值換算成十進位
def floTurnTen(x):
    xx=(x-(-3.0))*255/15.1
    return round(xx)
#十進位轉換成字串基因碼
def turnStrGene(x):
    return bin(x)[2:].zfill(string_length)
#二進位的字串轉回十進位
def bin_Int(x):
    return int(x,2)
#十進位換算成真實浮點數
def tenTurnflo(x):
    xx=15.1*x/255-3.0
    return round(float(xx),1)
