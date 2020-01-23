import Utils.utilsRandomizeScerevisiae

templ = { \
	"fileInput": "./Data/ScerevisiaeNetwork.json", \
	"fileInputProteins": "./Data/ScerevisiaeTF.json", \
	"fileOutput": "./Results/dataRandomized.json", \
	"randomizeCount": 10000, \
	"saveAfterCount": 100, \
}

rm = Utils.utilsRandomizeScerevisiae.randomizeMotifs(templ)
rm.executeRandomize()
