import Utils.utilsRandomize

templ = { \
	"fileInput": "./Data/links.json", \
	"fileInputProteins": "./Data/proteinList.json", \
	"fileOutput": "./Results/dataRandomized.json", \
	"randomizeCount": 10000, \
	"saveAfterCount": 100, \
}

rm = Utils.utilsRandomize.randomizeMotifs(templ)
rm.executeRandomize()
