# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 13:46:07 2019

@author: USER
"""

import numpy as np
import matplotlib.pyplot as plt

#產生100筆資料，每筆資料都是2個數字
X = np.random.rand(100,2)

#第一筆長這樣
X[0]

#畫出來看看，想當然是平均的佈滿整個畫面
#然後我們會用KMeans硬把他分類(明明沒意義的100個點……但他就是分的出來)
plt.scatter(X[:,0],X[:,1],s=50)

#接下來匯入KMeans函式庫
from sklearn.cluster import KMeans

#請KMeans分成三類
clf = KMeans(n_clusters=2)

#開始訓練！
clf.fit(X)

#這樣就可以取得預測結果了！
clf.labels_

#最後畫出來看看
#真的分成三類！太神奇了………無意義的資料也能分～
plt.scatter(X[:,0],X[:,1], c=clf.labels_)