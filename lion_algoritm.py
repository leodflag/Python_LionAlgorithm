# -*- coding: utf-8 -*-
"""
@author: leodflag
Genetic Algorithm
"""
import random,math
import matplotlib.pyplot as plt

string_length = 8 # �����L��
pop_lion_num=6 # ĸȺ�w����
best_male_lion=[] # ��ѹ��{
best_female_lion=[] # ���ĸ�{
itera=50  # �����Δ�
mutation_rate=0.1  # ͻ׃��

# �{Ⱥ�����OӋ
male_group=[] # ���{Ⱥĸ�w
female_group=[] # ĸ�{Ⱥĸ�w
male_cubs_group=[] # ���ת{Ⱥ
female_cubs_group=[] # ĸ�ת{Ⱥ
cubs_group=[] # ��������ת{Ⱥ

# ���򹠇��Д�
lion_gene1_U=255
lion_gene1_L=0
lion_gene2_U=149
lion_gene2_L=120

# �L�D����
pltX=[] # x�S������Δ�
pltY=[] # y�S��ƽ���m��ֵ
"""
�D�Q
"""
# ���H������ֵ�Q���ʮ�Mλ
def floTurnTen(x):
    xx=(x-(-3.0))*255/15.1
    return round(xx)
# ʮ�Mλ�D�Q���ִ�����a
def turnStrGene(x):
    return bin(x)[2:].zfill(string_length)
# ���Mλ���ִ��D��ʮ�Mλ
def bin_Int(x):
    return int(x,2)
# ʮ�Mλ�Q����挍���c��
def tenTurnflo(x):
    xx=15.1*x/255-3.0
    return round(float(xx),1)
"""
����Ĳ���
"""
# ָ��λ��ͻ׃
def _invert_at(s, index):
    # "**"�����\�㣬s���ִ�
    return bin(int(s,2)^2 ** (index))[2:].zfill(string_length)
# ��С����
def order(x):
    temp=0
    for i in range(6):
        for j in range(6):
            if x[j]<x[i]:
                temp=x[j]
                x[j]=x[i]
                x[i]=temp
    return x
# �≺�s   
def zipReturn(pool):
    x=[]
    adept[:],x[:]=zip(*pool)
    print('adept,x',adept,x)
    print('x',x)
    gene_1[:],gene_2[:]=zip(*x)
    print('gene_1,gene_2',gene_1,gene_2)
    return gene_1,gene_2
# �Д��С��ݔ��int
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
# �h��ǰ6������
def delList(x):
    if itera!=50:
       del x[0:6]     
# �a�������{Ⱥ
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
# ���䣺���S�C�c���䣬����ͻ׃��k_meas��Ⱥ
def lion_mating(gene_1,gene_2):
    # �S�C���c����
    i=0
    cut=0
    gene1=gene_1 # ԭ���{
    gene2=gene_2
    while(cut < pop_gene_num):
        dot=random.randint(0,7) # �S�Cȡ���Q�c��λ��
        gene1[cut]=gene_1[cut][0:dot]+gene_2[cut][dot:8] # ���{1��ǰ�λ����ĸ�{1����λ���Y��
        gene1[cut+1]=gene_1[cut+1][0:dot]+gene_2[cut+1][dot:8] # ���{2��ǰ�λ����ĸ�{2����λ���Y��
        gene2[cut]=gene_2[cut][0:dot]+gene_1[cut][dot:8] # ĸ�{1��ǰ�λ�������{1����λ���Y��
        gene2[cut+1]=gene_2[cut+1][0:dot]+gene_1[cut+1][dot:8] # ĸ�{1��ǰ�λ�������{1����λ���Y��
        cut+=2

    #�S�C���cͻ׃��ͻ׃���ֻ���
    s1=gene_1
    s2=gene_2
    for i in range(pop_gene_num):
        rand = random.random()
        if rand <= mutation_rate:
            #ͻ׃��ݔ���2λԪ��ݔ���ִ�
            h=random.randint(0,7)#�S�Cλ��ͻ׃
            s1[i]=_invert_at(gene_1[i],h)#������aͻ׃
            s2[i]=_invert_at(gene_2[i],h)

    return gene1,gene2

"""
�m������Adaptation function
"""
def Adaptation_x1(x1):  
    f1=x1*math.sin(4*math.pi*x1)
    return f1
def Adaptation_x2(x2):  
    f2=x2*math.sin(20*math.pi*x2)
    return f2
"""
1.�a���{Ⱥ
"""
male_group,female_group=init_lion_gene(male_group,female_group,pop_lion_num)
"""
2.�������
"""
cubs_group=male_group.extend(female_group)
print('cubs',cubs_group)