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

#�L�D����
pltX=[]#x�S������Δ�
pltY=[]#y�S��ƽ���m��ֵ
"""
�D�Q
"""
#���H������ֵ�Q���ʮ�Mλ
def floTurnTen(x):
    xx=(x-(-3.0))*255/15.1
    return round(xx)
#ʮ�Mλ�D�Q���ִ�����a
def turnStrGene(x):
    return bin(x)[2:].zfill(string_length)
#���Mλ���ִ��D��ʮ�Mλ
def bin_Int(x):
    return int(x,2)
#ʮ�Mλ�Q����挍���c��
def tenTurnflo(x):
    xx=15.1*x/255-3.0
    return round(float(xx),1)
