import Utils.utilsStatistics

templStats = { \
    "fileInputLinks": "./Data/links.json", \
    "fileMotifRandomized": "./Results/dataRandomized-1.json", \
    "fileInputOhnoloque": "./Data/ohnologues.json", \
    "fileOutputStats": "./Results/statistics.csv", \
    "ratio": 0.05 \
}

ms = Utils.utilsStatistics.utilsMotifs(templStats)
ms.finalStats()
