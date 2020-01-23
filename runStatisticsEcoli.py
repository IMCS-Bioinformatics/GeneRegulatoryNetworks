import Utils.utilsStatisticsEcoli

templStats = { \
    "fileInputLinks": "./Data/EcoliNetwork.json", \
    "fileInputGenePairs": "./Data/EcoliPairs.json", \
    "fileOutputStats": "./Results/Statistics-rand-pairs-ecoli-111.csv", \
}

ms = Utils.utilsStatisticsEcoli.utilsGraphlets(templStats)
ms.finalStats()
