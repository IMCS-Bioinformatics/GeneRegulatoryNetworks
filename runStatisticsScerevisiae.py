import Utils.utilsStatisticsScerevisiae

templStats = { \
    "fileInputLinks": "./Data/ScerevisiaeNetwork.json", \
    "fileMotifRandomized": "./Results/dataRandomized-1.json", \
    "fileInputOhnoloque": "./Data/ScerevisiaeHomologues.json", \
    "fileOutputStats": "./Results/Statistics-rand-pairs-scerevisiae.csv", \
    "ratio": 0.05, \
    "swapPairsRandom": True, \
    "seed": 997 \
}

ms = Utils.utilsStatisticsScerevisiae.utilsMotifs(templStats)
ms.finalStats()
