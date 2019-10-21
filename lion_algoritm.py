# -*- coding: utf-8 -*-
"""
@author: leodflag
Genetic Algorithm
"""
import random,math
import numpy as np
import copy
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
# 獅群參數設計
string_length = 8 # 基因長度
pop_lion_num=6 # 母群體數量
best_male_lion=[] # 最佳公獅
best_female_lion=[] # 最佳母獅
itera=5  # 迭代次數
mutation_rate=0.1  # 突變率
growing_time=3
Bstrength=2
Bcount=0
male_group=[] # 公獅群母體
female_group=[] # 母獅群母體
cubs_group=[] # 交配後的幼獅群
lion_group_ALL=[] # 總獅群組
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
#---------適應函數Adaptation function---------
def Adaptation_x1(x1):  
    f1=x1*math.sin(4*math.pi*x1)
    return f1
def Adaptation_x2(x2):  
    f2=x2*math.sin(20*math.pi*x2)
    return f2
def Adaptation(x1,x2):
    f=x1*math.sin(4*math.pi*x1)+x2*math.sin(20*math.pi*x2)
    return f
# 只取前面的基因 越大越好，所以要從後面開始刪
def delList(x,del_num):
    del x[del_num:len(x)] 
# 大到小排序，會用x[i][0]排
def order(x):
    temp=0
    for i in range(len(x)):
        for j in range(len(x)):
            if x[j]<x[i]:
                temp=x[j]
                x[j]=x[i]
                x[i]=temp
    return x
# 比較大小
def comparison_size(x1,x2):
    dis_size=len(x1)-len(x2)
    end_n=0
    x1_A=[]
    x2_A=[]
    x1_G=[]
    x2_G=[]
    if dis_size>0: # 若x1的基因數多於x2
        for i in range(len(x1)):
            x1_A.append(Adaptation_x1(x1[i])) # 計算適應函數並加到新x1_order列表裡
        x1_order=list(zip(x1_A,x1)) # 綁定適應值與基因的順序
        x1_order=order(x1_order)
        end_n=len(x1)-dis_size # 取得x1應該停在哪個位置，也就是要取幾個
        delList(x1_order,end_n) # 只取前面的基因
        x1_A[:],x1_G[:]=zip(*x1_order) # 解壓縮，分離預測值與基因
        return x1_G,x2
    elif dis_size<0:
        for i in range(len(x2)):
            x2_A.append(Adaptation_x2(x2[i])) # 計算適應函數並加到新x1_order列表裡
        x2_order=list(zip(x2_A,x2)) # 綁定適應值與基因的順序
        x2_order=order(x2_order)
        end_n=len(x2)+dis_size # 取得x1應該停在哪個位置，也就是要取幾個
        delList(x2_order,end_n) # 只取前面的基因
        #print('x2_orderx2_order',x2_order)
        x2_A[:],x2_G[:]=zip(*x2_order) # 解壓縮，分離預測值與基因
        return x1,x2_G       
    else:
        return x1,x2
# 判斷基因區間，輸入int
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
#---------K-means function---------
def kmeans_clusters(X): # 輸入n*1大小array做kmeans
    X.reshape(1,-1) # 包含單個樣本，重整數據
    clf = KMeans(n_clusters=2) # 分2類
    clf.fit(X) # 以原設定訓練
    return clf.labels_ # 返回訓練完成的標籤
# 產生初代獅群
def init_lion_gene(gene_1,gene_2,pop_gene_num):
    for i in range(pop_gene_num*2):
        s1=random.randint(0,255)
        if (i%2) == 0:
            gene_1.append(turnStrGene(s1))
        else:
            gene_2.append(turnStrGene(s1))
    #print('male lion group',gene_1)
    #print('female lion group',gene_2)
    return gene_1,gene_2 
# 交配：含隨機點交配，基因突變，k_meas分群
def lion_mating(gene_1,gene_2):
    # 隨機單點交配，公母交配產生幼獅
    i=0
    cut=0
    gene_ALL=[] # 多獅群組
    while(cut < len(gene_1)): # 假設有6對
        gene1=[] # 交配後代基因
        for i in range(2):
            dot=random.randint(0,7) # 隨機取交換點的位置
            gene1.append(gene_1[cut][0:dot]+gene_2[cut][dot:8]) # 公獅1的前段基因跟母獅1的後段基因結合 一對只交配一次產生2後代
            gene1.append(gene_2[cut][0:dot]+gene_1[cut][dot:8]) # 母獅1的前段基因跟公獅1的後段基因結合
        #隨機單點突變，突變部分基因
        cubs_group=gene1 # 幼獅群
        cubs_group2=copy.deepcopy(cubs_group) # 2n隻幼獅複製產生新2n隻幼獅
        for i in range(len(cubs_group2)):
            rand = random.random()
            if rand <= mutation_rate: # 若隨機數小於突變率
                h=random.randint(0,7) # 隨機位置突變，輸入為2位元，輸出字串
                cubs_group2[i]=_invert_at(cubs_group2[i],h) # 將基因碼突變
        cubs_group=cubs_group+cubs_group2 # 整合成 1 list
        cubs_group2_int=[]
        for i in range(len(cubs_group)): # 為了使用kmeans，做成n*1 list
            cubs_group_int=[]
            cubs_group_int.append(bin_Int(cubs_group[i])) # 換回十位數
            cubs_group2_int.append(cubs_group_int)
        X=np.array(cubs_group2_int) # list 轉 array，逗號會消失，但可直接丟到kmeans function做分類
        X=kmeans_clusters(X) # 輸入至kmeans function做分類，輸出分類結果 
        male_cubs_group=[] # 公幼獅群
        female_cubs_group=[] # 母幼獅群
        for i in range(len(X)):
            if X[i]==1: # 依照分類結果分公幼獅群和母幼獅群
                male_cubs_group.append(tenTurnflo(cubs_group2_int[i][0])) # 換回真實浮點數
            else:
                female_cubs_group.append(tenTurnflo(cubs_group2_int[i][0]))
        male_cubs_group_1,female_cubs_group_1=comparison_size(male_cubs_group,female_cubs_group)  
        gene_group=[] # 一公一母與孩子的獅群
        gene_group.append(gene_1[cut]) # 加入交配的一隻公獅
        gene_group.append(gene_2[cut]) # 加入交配的一隻母獅
        gene_group.append(male_cubs_group_1) # 加入公幼獅群
        gene_group.append(female_cubs_group_1) # 加入母幼獅群
        gene_ALL.append(gene_group) # 成為獅群組裡的其中一個獅群
        cut+=1
    return gene_ALL
# 領土防禦
def territorial_defense(group_all,growing_t):
    group_ALL=copy.deepcopy(group_all) # 複製輸入的獅群組
    male_a=0
    for i in range(len(group_ALL)):
        times=growing_t
        while(times!=0): # 若成長時間不歸零
            nomad=random.randint(0,255) # 產生隨機流浪獅
            nomad_a=Adaptation_x1(tenTurnflo(nomad)) # 計算流浪獅的適應值
            male_a=Adaptation_x1(tenTurnflo(bin_Int(group_ALL[i][0]))) # 計算公獅的適應值
            if nomad_a>male_a: # 若流浪獅適應值高於公獅
                nomad_g=[]
                nomad_g.append(turnStrGene(nomad)) # 轉換十進位流浪獅基因成二進位，
                female_g=[] 
                female_g.append(group_ALL[i][1]) # 獨立出原母獅二進位基因
                new_group=lion_mating(nomad_g,female_g) # 交配後產生新獅群
                #print('==========================',i)
                #print('group_ALL[%d]=',group_ALL[i])
                group_ALL[i]=new_group[0] # 將新獅群放回總獅群裡
                #print('group_ALL[%d]=',group_ALL[i])
                times=growing_t # 新幼獅群重新成長
                #print('times=========if',times)
            else:
                times-=1 # 平安無事，幼獅成長
    return group_ALL

def territorial_takeover(group_all):
    best_male_lion=[]
    best_female_lion=[]
    for i in range(len(group_all)):
        female=[] # 獅群組
        male=[]        
        adept=[]
        male.append(group_all[i][2]) # 公幼獅群
        male_group=male[0] # 加幼公獅群到公獅群
        male.clear()
        male_group.append(tenTurnflo(bin_Int(group_all[i][0]))) # 轉換原公獅加進公獅群
        for m in range(len(male_group)):    
            adept.append(Adaptation_x1(male_group[m])) # 計算出適應值
        male=list(zip(adept,male_group))
        male=order(male)
        best_male_lion.append(male[0])
        #print('-------------',i)
        #print(best_male_lion)
        female.append(group_all[i][3]) # 母幼獅群
        female_group=female[0] # 加母公獅群到母獅群
        female.clear()
        female_group.append(tenTurnflo(bin_Int(group_all[i][1]))) # 轉換原母獅加進母獅群
        adept.clear()
        for m in range(len(female_group)):    
            adept.append(Adaptation_x2(female_group[m])) # 計算出適應值
        female1=[]
        female1=list(zip(adept,female_group))
        female1=order(female1)
        Bcount=0
        best_female_lion.append(female1[0])
        for f in range(len(female1)):
            if best_female_lion[i]==female1[f]:
                Bcount+=1
        if Bcount>Bstrength:
            adept=[]
            female1=[]
            new=[]
            new_f=random.randint(0,255)
            new_f_flo=tenTurnflo(new_f)
            new_adept=0
            new_adept=Adaptation_x2(new_f_flo)           
            if new_f_flo != best_male_lion[i][1] and new_adept>best_female_lion[i][0]:
                adept.append(new_adept)
                female1.append(new_f_flo)
                new=list(zip(adept,female1))    
                best_female_lion[i]=new[0]
        else:
            Bcount=0
    
    adept[:],male[:]=zip(*best_male_lion)
    adept[:],female[:]=zip(*best_female_lion)
    #print('male55555',male)
    return male,female    
        
#---------1.產生獅群---------
male_group,female_group=init_lion_gene(male_group,female_group,pop_lion_num)
while(itera>0):
    lion_group_ALL=[]
    lion_group_ALL_1=[]   
    best_male_lion_1=[]
    best_female_lion_1=[]
    #---------2.繁衍後代---------
    lion_group_ALL=lion_mating(male_group,female_group)
    #---------3.領土防禦---------
    lion_group_ALL_1=territorial_defense(lion_group_ALL,growing_time)
    #---------4.領土爭奪---------
    best_male_lion_1,best_female_lion_1=territorial_takeover(lion_group_ALL_1)
    male_group=[]
    female_group=[]
    for i in range(len(best_male_lion_1)):
        male_group.append(turnStrGene(floTurnTen(best_male_lion_1[i])))
        female_group.append(turnStrGene(floTurnTen(best_female_lion_1[i])))
    pltX.append(50-itera)#圖的X軸為次數
    itera-=1