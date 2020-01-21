import Utils.utilsStatistics

templStats = { \
    "fileInputLinks": "./Data/links.json", \
    "fileMotifRandomized": "./Results/dataRandomized-1.json", \
    "fileInputOhnoloque": "./Data/ohnologues.json", \
    "fileOutputStats": "./Results/statistics-rand-pairs.csv", \
    "ratio": 0.05, \
    "swapPairsRandom": True, \
    "seed": 997 \
}

ms = Utils.utilsStatistics.utilsMotifs(templStats)
ms.finalStats()
