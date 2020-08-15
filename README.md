# Python_LionAlgorithm
2012 Lion Algorithm

## 目標
    Max f (x1, x2 ) = 21.5 + x1sin(4πx1) + x2sin(20πx2 )

## 條件限制
    −3.0 ≤ x1 ≤ 12.1,  4.1 ≤ x2 ≤ 5.8

## 虛擬碼
    1 產生獅群
        1. 1產生2n隻獅子並隨機均分成公和母兩種，生成隨機數時採用基因範圍。
    2 繁衍後代
        2.1. 公和母交配產生2隻幼獅。
        2.2.  2n隻幼獅複製產生新2n隻幼獅，透過突變率決定是否突變，沒突變的也不管，直接抓回去變成總共4n隻幼獅。
        2.3.  4n隻以k-means分群法分成小公獅小母獅2群，以數量較少的群數量為準（假設是小公獅），數量較多的群（小母獅）
        去計算內部的全部適應值並排序，去掉差的直到跟少數量群數量（小公獅）一樣。若kmeans分群法將幼獅群全部分成同一類
        ，則重新分群，3次後若仍為同一類，則強制將第一隻幼獅換成另一類。
        2.4. 將一公一母加上其產生的幼獅群當作一整個獅群，因此在第二步驟繁衍後代時會產生共n群。
    3 領土防禦
        3.1 設成長時間為3，以公獅範圍產生一隻流浪獅子對一獅群（一公一母＋幼獅群），只跟公獅比較適應值，若流浪獅子挑戰公獅成功，
        則殺掉所有幼獅群，回到第二步驟跟原本的母獅繁衍後代，如果公獅防禦成功則繼續跟新產生的流浪獅子比較，直到幼獅全部成長成成獅
        才能停止迴圈，此時原本有n群獅群，獅群內公母都是成獅，且公母數量相同，比如說獅群1有公母成獅各12隻……獅群10有公母成獅各10隻。
    4.領土爭奪
        4.1 獅群1為例，計算12隻公獅的適應值做由大到小的排列，取出第一隻公獅代表獅群1。
        4.2計算12隻母獅的適應值做由大到小的排列，取出第一隻母獅當最佳，若有跟最佳母獅適應值一樣的算入Bcout＋＋，若超過設定的Bstrength
        假設為2隻，要以母獅的範圍新生成一隻母獅且不能跟公獅一樣，直到新生成的母獅強於原先選出來的最佳母獅便替換掉牠，變成代表獅群1的母獅。

## 函式
    import copy
    copy.deepcopy() # Deep copy operation on arbitrary Python objects.

    import time
    tStart = time.time() # 計時開始
    tEnd = time.time()#計時結束
    print("It cost %f sec" % (tEnd - tStart)) # 會自動做進位
    print(tEnd - tStart) # 原型長這樣

    from sklearn.cluster import Kmeans
    def kmeans_clusters(X): # 輸入n*1大小array做kmeans
        X.reshape(1,-1) # 包含單個樣本，重整數據
        clf = KMeans(n_clusters=2) # 分2類
        clf.fit(X) # 以原設定訓練
        return clf.labels_ # 返回訓練完成的標籤
        
  ## 結果
  [基因](https://github.com/leodflag/Python_LionAlgorithm/blob/master/gene_result.png)
  [平均適應值](https://github.com/leodflag/Python_LionAlgorithm/blob/master/adaptation_result.png)
