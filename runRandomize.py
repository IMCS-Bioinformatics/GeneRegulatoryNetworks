import Utils.utilsRandomize

templ = { \
	"fileInput": "./Data/links.json", \
	"fileInputProteins": "./Data/proteinList.json", \
	"fileOutput": "./Results/dataRandomized.json", \
	"randomizeCount": 12, \
	"saveAfterCount": 3, \
}

rm = Utils.utilsRandomize.randomizeMotifs(templ)
rm.executeRandomize()
