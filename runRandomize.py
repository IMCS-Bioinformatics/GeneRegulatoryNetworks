import Utils.utilsRandomize

templ = { \
    "fileInput": "./Data/ScervisiaeNetwork.json", \
    "fileInputProteins": "./Data/ScervisiaeTF.json", \
	"fileOutput": "./Results/dataRandomizedScervisiae.json", \
	"randomizeCount": 12, \
	"saveAfterCount": 3, \
}

rm = Utils.utilsRandomize.randomizeMotifs(templ)
rm.executeRandomize()