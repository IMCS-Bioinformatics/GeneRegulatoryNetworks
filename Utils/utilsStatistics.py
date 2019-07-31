import os.path, math, random, time
import Utils.utilsReader

class utilsMotifs:
    def __init__(self, templ):
        self.template = templ
        self.delim = "#"

    def gene(self,o): return o.upper()
    def initStatsData(self):
        self.links = Utils.utilsReader.readJson(self.template["fileInputLinks"])
        self.randomizedMotifs = Utils.utilsReader.readJson(self.template["fileMotifRandomized"])
        self.ohnologues = Utils.utilsReader.readJson(self.template["fileInputOhnoloque"])
        self.count = self.randomizedMotifs["randCount"]
        self.motifCount = self.randomizedMotifs["motifCount"]
        self.proteinList = self.randomizedMotifs["proteinList"]
        self.ratio = self.template["ratio"]
        self.motifNames = ["FFLX", "FFLY", "FFLZ", "FBLDX", "FBLDY", "FBLUY", "FBLUZ", "C2", "C3", "C4"]
        if self.template["swapPairsRandom"]:
            self.swapPairsRandom()
    def swapPairsRandom(self): 
        random.seed(self.template["seed"])
        self.proteinList = []
        for i in range(len(self.ohnologues)):
            r = random.randint(0,1) 
            if r == 1:
                v = self.ohnologues[i][0]
                self.ohnologues[i][0] = self.ohnologues[i][1]
                self.ohnologues[i][1] = v
            self.proteinList.append(self.ohnologues[i][0])
            self.proteinList.append(self.ohnologues[i][1])
    def validMotifValue(self, v): return (v / self.count) <= self.ratio
    def calcProteinMotifs(self):
        self.proteinMotifs = {}
        pSet = set(self.proteinList)
        motifKeys = self.motifCount.keys()        
        for p in pSet:
            self.proteinMotifs[p] = {}
            for k in self.motifNames: 
                self.proteinMotifs[p][k] = set()
        for k in ["C2", "C3", "C4"]:
            for m in self.motifCount[k].keys():
                pList = m.split(self.delim)
                for p in pList:
                    if p in pSet: self.proteinMotifs[p][k].add(m)
        for m in self.motifCount["FFL"].keys():
            p = m.split(self.delim)
            if p[0] in pSet: self.proteinMotifs[p[0]]["FFLX"].add(m)
            if p[1] in pSet: self.proteinMotifs[p[1]]["FFLY"].add(m)
            if p[2] in pSet: self.proteinMotifs[p[2]]["FFLZ"].add(m)
        for m in self.motifCount["FBLD"].keys():
            p = m.split(self.delim)
            if p[0] in pSet: self.proteinMotifs[p[0]]["FBLDX"].add(m)
            if p[1] in pSet: self.proteinMotifs[p[1]]["FBLDY"].add(m)
            if p[2] in pSet: self.proteinMotifs[p[2]]["FBLDY"].add(m)
        for m in self.motifCount["FBLU"].keys():
            p = m.split(self.delim)
            if p[0] in pSet: self.proteinMotifs[p[0]]["FBLUY"].add(m)
            if p[1] in pSet: self.proteinMotifs[p[1]]["FBLUY"].add(m)
            if p[2] in pSet: self.proteinMotifs[p[2]]["FBLUZ"].add(m)
    def getLinkCount(self, oList):
        oSet = set(oList)
        n = 0
        for e in self.links:
            if e[0] in oSet or self.gene(e[1]) in oSet: n += 1
        return n
    def removeFromMotif(self, m, v):
        p = m.split(self.delim)
        p.remove(v)
        return p[0] + self.delim + p[1]
    def removeFromSet(self, v, mSet):
        newSet = set()
        for m in mSet: newSet.add(self.removeFromMotif(m, v))
        return newSet
    def countRemoveTriple(self, v, vSet, w, wSet):
        return [len(vSet), len(wSet), len(self.removeFromSet(v, vSet).intersection(self.removeFromSet(w, wSet)))]
    def countTriple(self, vSet, wSet):
        return [len(vSet), len(wSet), len(vSet.intersection(wSet))]
    def getValidMotifSet(self, kk, v):
        # ipr = 10
        vSet = set()
        if kk[0] != "C": k = kk[:-1]
        else: k = kk
        if kk in self.proteinMotifs[v]:
            for m in self.proteinMotifs[v][kk]:
                if self.validMotifValue(self.motifCount[k][m]): 
                    vSet.add(m)
        return vSet
    def ohnologueStats(self):
        header = ["protein-1", "protein-2", "value", "type", "links"]
        # for k in ["FFLX", "FFLY", "FFLZ", "FBLDX", "FBLDY", "FBLUY", "FBLUZ", "C2", "C3", "C4"]:
        #     for kk in ["", "-R"]:
        #         header += [k + kk + "-1", k + kk + "-2", k + kk + "-12"]

        headerNew = ["protein-1", "protein-2", "value", "type", "links"]
        for k in ["FFL-X", "FFL-Y", "FFL-Z", "DFFLUp-X", "DFFLUp-Y", "DFFLDown-Y", "DFFLDown-Z", "C2", "C3", "C4"]:
            for kk in ["", "-R"]:
                headerNew += [k + "-1" + kk, k + "-2" + kk, k + "-12" + kk]

        # rez = [header]
        rez = [headerNew]
        for o in self.ohnologues:
            if o[0] in self.proteinMotifs and o[1] in self.proteinMotifs:
                values = [o[0], o[1], o[2], o[3], self.getLinkCount(o[:2])]
                for k in ["FFLX", "FFLY", "FFLZ", "FBLDX", "FBLDY", "FBLUY", "FBLUZ"]:
                    values += self.countRemoveTriple(o[0], self.proteinMotifs[o[0]][k], o[1], self.proteinMotifs[o[1]][k])
                    values += self.countRemoveTriple(o[0], self.getValidMotifSet(k, o[0]), o[1], self.getValidMotifSet(k, o[1]))
                for k in ["C2", "C3", "C4"]:
                    values += self.countTriple(self.proteinMotifs[o[0]][k], self.proteinMotifs[o[1]][k])
                    values += self.countTriple(self.getValidMotifSet(k, o[0]), self.getValidMotifSet(k, o[1]))
                rez.append(values)
        return rez
    def simgleProteinStats(self):
        header = ["protein", "links to gene"]
        for k in self.motifNames:
            header.append(k)
            header.append(k + "-R")
        rez = [header]
        for p in self.proteinList:
            ecount = 0
            for e in self.links:
                if e[0] == p or self.gene(e[1]) == p: ecount += 1
            row = [p, ecount]
            for k in self.motifNames:
                kk = k
                if len(kk) > 2: kk = kk[:len(kk) - 1]
                rcount = 0
                for m in self.proteinMotifs[p][k]:
                    # if True: 
                    if self.validMotifValue(self.motifCount[kk][m]): 
                        rcount += 1
                row.append(len(self.proteinMotifs[p][k]))
                row.append(rcount)
            rez.append(row)
        return rez
    def addNormalizationByLinks(self, rez):
        # header = rez[0]
        # ilinks = header.index("links")
        # print (ilinks, header)
        # n = ilinks + 1
        # nn = len(header)
        # newHeader = header[:n]
        # for k in header[n:]:
        #     newHeader += [k, k + "_N"]
        # newRez = [newHeader]
        # for d in rez[1:]:
            # x = float(d[n - 1])
            # values = []
            # for i in range(n): values.append(d[i])
            # for i in range(n, nn): values += [d[i], d[i] / x]
        #     newRez.append(values)
        # CReader.saveCsv(self.template["fileOutputStats"], newRez)
        Utils.utilsReader.saveCsv(self.template["fileOutputStats"], rez)
    def finalStats(self):
        self.initStatsData()
        self.calcProteinMotifs()

        rez = self.ohnologueStats()
        self.addNormalizationByLinks(rez)
    def initSubgraphData(self):
        self.randomizedMotifs = CReader.readJson(self.template["fileMotifRandomized"])
        self.count = self.randomizedMotifs["randCount"]
        self.proteinList = self.randomizedMotifs["proteinList"]
        self.motifCount = self.randomizedMotifs["motifCount"]
        self.testProteins = self.template["proteinsForSubgraph"]
        pInters = set(self.proteinList).intersection(set(self.testProteins))
        if len(pInters) != len(self.testProteins):
            print ("\n### Wrong list of 'proteinsForSubgraph'. Proteins ", \
                list(pInters), \
                "is not in randomization 'proteinList'")
    def getCX(self, motif):
        links = []
        pList = motif.split(self.delim)
        v = pList[len(pList) - 1]
        hasM = True
        for w in pList: 
            links.append[self.gene(v), self.protein(v)]
            links.append[self.protein(v), self.gene(w)]
        return links
    def testFBLD(self):
        for m in self.motifCount["FBLD"].keys():
            if m not in self.motifCount["FFL"]:
                print ("FBLD", "FFL m", m)
            p = m.split(self.delim)
            n = p[0] + self.delim + p[2] + self.delim + p[1]
            if n not in self.motifCount["FFL"]:
                print ("FBLD", "FFL n", n)
    def testFBLU(self):
        i = 0
        for m in self.motifCount["FBLU"].keys():
            if m not in self.motifCount["FFL"]:
                print ("FBLU", "FFL m", m)
            p = m.split(self.delim)
            n = p[1] + self.delim + p[0] + self.delim + p[2]
            if n not in self.motifCount["FFL"]:
                print ("FBLU", "FFL n", n)
                i += 1
        print ("FFL", len(self.motifCount["FFL"].keys()))
        print ("FBLU", len(self.motifCount["FBLU"].keys()))