import WP3statistic

templStats = { \
    "fileInputLinks": "Data/linksWP3.json", \
    "fileMotifRandomized": "Rez/rand_res.json", \
    "fileInputOhnoloque": "Data/ohnologueWP3.json", \
    "fileOutputStats": "Rez/randomize_stats.csv", \
    "ratio": 0.0005 \
}

ms = WP3statistic.utilsMotifs(templStats)
ms.finalStats()