import WP3randomize

templ = { \
    "fileInput": "Data/linksWP3.json", \
    "fileInputProteins": "Data/proteinList.json", \
	"fileOutput": "Rez/rand_res.json", \
	"randomizeCount": 12, \
	"saveAfterCount": 3, \
}

rm = WP3randomize.randomizeMotifs(templ)
rm.executeRandomize()