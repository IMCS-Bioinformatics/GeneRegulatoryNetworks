import Utils.utilsStatistics

templStats = { \
    "fileInputLinks": "./Data/links.json", \
    "fileMotifRandomized": "./Results/dataRandomized.json", \
    "fileInputOhnoloque": "./Data/ohnologues.json", \
    "fileOutputStats": "./Results/statistics.csv", \
    "ratio": 0.0005 \
}

ms = Utils.utilsStatistics.utilsMotifs(templStats)
ms.finalStats()
