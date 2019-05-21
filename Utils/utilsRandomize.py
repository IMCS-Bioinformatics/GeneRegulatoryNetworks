import os.path, math, random, time
import Utils.utilsReader

class randomizeMotifs:
    def __init__(self, templ):
        self.template = templ
        self.delim = "#"
        self.motifSets = {}
    def initData(self):
        self.links = Utils.utilsReader.readJson(self.template["fileInput"])
        self.proteinList = Utils.utilsReader.readJson(self.template["fileInputProteins"])
        self.updateCand()
    def isProtein(self,o): return o.islower()
    def isGene(self,o): return not self.isProtein(o)
    def isEdge(self,o): return delim in o
    def protein(self,o): return o.lower()
    def gene(self,o): return o.upper()
    def edgeDId(self,v, w): return str(v) + self.delim + str(w)
    def edgeUId(self,v, w):
        if w > v: return self.edgeDId(v, w) 
        return self.edgeDId(w, v)
    def hasLink(self, g, p): return self.edgeDId(g, p) in self.linkSet
    # uzģenerē saistītos
    def updateCand(self):
        gSet = set()
        pSet = set()
        # rēķina pgSet
        for e in self.links:
            if self.isGene(e[1]):
                gSet.add(self.protein(e[1]))
                pSet.add(e[0])
        self.pgSet = pSet.intersection(gSet)
        self.pgCand = [{}, {}]
        self.pCand = {}
        self.gCand = {}
        for e in self.links:
            if self.isGene(e[1]):
                # rēķina pgCand
                for i in [0, 1]:
                    v = self.protein(e[i])
                    w = self.protein(e[1 - i])
                    if v in self.pgSet and v != w:
                        if v not in self.pgCand[i]: self.pgCand[i][v] = set()
                        self.pgCand[i][v].add(w)
                # rēķina pCand
                v = e[0]
                w = self.protein(e[1])
                if w in self.pgSet: 
                    if v not in self.pCand: self.pCand[v] = set()
                    self.pCand[v].add(w)
                if v in self.pgSet: 
                    if w not in self.gCand: self.gCand[w] = set()
                    self.gCand[w].add(v)
        self.pgLinkCount = len(self.links)
        for p in self.pgSet: self.links.append([self.gene(p), p])
        self.linkSet = set()
        for e in self.links: self.linkSet.add(self.edgeDId(e[0], e[1]))
        self.pgCandSet = set(self.pgCand[0].keys()).intersection(set(self.pgCand[1].keys()))
        # print ("updateCand pCand, gCand", len(self.pCand), len(self.gCand))

    def uId(self, c):
        i = 0
        n = len(c)
        for j in range(n):
            if c[i] > c[j]: i = j
        s = ""
        d = ""
        for j in range(n):
            s += d + c[(i + j) % n]
            d = self.delim
        return s
    def dId(self, c):
        s = ""
        d = ""
        for p in c:
            s += d + p
            d = self.delim
        return s
    def pgC4I(self, p):
        pCand = set()
        gCand = set()
        rezSet = set()
        ci = 0
        j = 0
        if self.hasLink(self.gene(p), p) and p in self.pgCand[0] and p in self.pgCand[1]:
            for pi in self.pgCandSet:
                p01Set = (self.pgCand[0][p].intersection(self.pgCand[1][pi]))
                p10Set = (self.pgCand[1][p].intersection(self.pgCand[0][pi]))
                ptSet = p01Set.intersection(p10Set)
                for p0 in p01Set:
                    for p1 in p10Set:
                        if p != p0 and p != p1 and p != pi and p0 != p1 and p0 != pi and p1 != pi:
                            rezSet.add(self.dId([p0, pi, p1]))
                            j += 1
        return rezSet
    def getC4I(self, pOhn):
    	return self.pgC4I(self.protein(pOhn))
    def pgC3I(self, p):
    	pCand = set()
    	gCand = set()
    	rezSet = set()
    	ci = 0
    	if self.hasLink(self.gene(p), p) and p in self.pgCand[0] and p in self.pgCand[1]:
    		for pi in self.pgCand[0][p]:
    			if pi in self.pgCand[0]:
    				pjSet = self.pgCand[1][p].intersection(self.pgCand[0][pi])
    				for pj in pjSet:
    					if pi != pj and p != pi and p != pj:
    						rezSet.add(self.dId([pi, pj]))
    	return rezSet
    def getC3I(self, pOhn):
    	return self.pgC3I(self.protein(pOhn))
    def pgC2I(self, p):
    	pCand = set()
    	gCand = set()
    	rezSet = set()
    	ci = 0
    	if self.hasLink(self.gene(p), p) and p in self.pgCand[0] and p in self.pgCand[1]: 
    		iSet = set(self.pgCand[0][p].intersection(self.pgCand[1][p]))
    		for pi in iSet:
    			rezSet.add(self.uId([pi]))
    	return rezSet
    def getC2I(self, pOhn):
    	return self.pgC2I(self.protein(pOhn))
    def addCX(self, po, cList):
    	mSet = set()
    	for c in cList:
    		cl = c.split(self.delim)
	    	cl.append(po)
    		vmin = min(cl)
    		i = 0
	    	while cl[i] != vmin: i += 1
    		cf = vmin
    		n = len(cl)
    		for j in range(1, n):
    			cf += self.delim + cl[(i + j) % n]
    		mSet.add(cf)
    	return mSet
    def testCX(self, mSet):
        comSet = set()
        for d in mSet:
            cList = d.split(self.delim)
            v = cList[len(cList) - 1]
            hasM = True
            for w in cList: 
                hasM = hasM and self.hasLink(self.protein(v), self.gene(w))
                v = w
            if hasM: comSet.add(d)
        return comSet
    def testFFL(self, mSet):
        comSet = set()
        for d in mSet:
            p = d.split(self.delim)
            hasM = \
                self.hasLink(self.gene(p[1]), p[1]) and \
                self.hasLink(p[0], self.gene(p[1])) and \
                self.hasLink(p[1], self.gene(p[2])) and \
                self.hasLink(p[0], self.gene(p[2]));
            if hasM: comSet.add(d)
        return comSet
    def testFBLD(self, mSet):
        comSet = set()
        for d in mSet:
            p = d.split(self.delim)
            hasM = \
                self.hasLink(self.gene(p[1]), p[1]) and \
                self.hasLink(self.gene(p[2]), p[2]) and \
                self.hasLink(p[1], self.gene(p[2])) and \
                self.hasLink(p[2], self.gene(p[1])) and \
                self.hasLink(p[0], self.gene(p[1])) and \
                self.hasLink(p[0], self.gene(p[2]));
            if hasM: comSet.add(d)
        return comSet
    def testFBLU(self, mSet):
        comSet = set()
        for d in mSet:
            p = d.split(self.delim)
            hasM = \
                self.hasLink(self.gene(p[0]), p[0]) and \
                self.hasLink(self.gene(p[1]), p[1]) and \
                self.hasLink(p[0], self.gene(p[1])) and \
                self.hasLink(p[1], self.gene(p[0])) and \
                self.hasLink(p[0], self.gene(p[2])) and \
                self.hasLink(p[1], self.gene(p[2]));
            if hasM: comSet.add(d)
        return comSet

    def createFFL(self):
        testSet = set(self.proteinList)
        fflSet = set()
        fblSet = set()
        for e in self.links:
            if self.isGene(e[1]):
                v = e[0]
                w = self.protein(e[1])
                if v != w and v in self.pCand and w in self.gCand: 
                    c = self.pCand[v].intersection(self.gCand[w]);
                    for u in c:
                        if v != u and w != u and (v in testSet or w in testSet or u in testSet): 
                            fflId = v + self.delim + u + self.delim + w
                            if fflId in fflSet:
                                print (e, u, fflId)
                                exit()
                            fflSet.add(fflId)
        self.motifSets["FFL"] = fflSet
    def createFBLD(self):
        fblSet = set()
        fflSet = self.motifSets["FFL"]
        for m in fflSet:
            p = m.split(self.delim)
            n = p[0] + self.delim + p[2] + self.delim + p[1]
            if n in fflSet:
                if p[1] > p[2]: m = n
                fblSet.add(m)
        self.motifSets["FBLD"] = fblSet
    def createFBLU(self):
        fblSet = set()
        fflSet = self.motifSets["FFL"]
        for m in fflSet:
            p = m.split(self.delim)
            n = p[1] + self.delim + p[0] + self.delim + p[2]
            if n in fflSet:
                if p[0] > p[1]: m = n
                fblSet.add(m)
        self.motifSets["FBLU"] = fblSet
    def createCycles(self):
        for t in ["C2", "C3", "C4"]:
            self.motifSets[t] = set()
            for po in self.proteinList:
                c = set()
                if t == "C2": c = self.getC2I(po)
                elif t == "C3": c = self.getC3I(po)
                elif t == "C4": c = self.getC4I(po)
                self.motifSets[t] |= self.addCX(po, c)

    def initRandomize(self):
        self.createCycles()
        self.createFFL()
        self.createFBLD()
        self.createFBLU()
        self.randCount = 0
        self.motifCount = {}
        for t in self.motifSets.keys():
            self.motifCount[t] = {}
            for c in self.motifSets[t]:
                self.motifCount[t][c] = 0
        # print ("FBLU diff FFL", len(set(self.motifCount["FBLU"].keys()).difference(set(self.motifCount["FFL"].keys()))))
        # print ("FFL diff FBLU", len(set(self.motifCount["FFL"].keys()).difference(set(self.motifCount["FBLU"].keys()))))
        # print ("FBLU inter FFL", len(set(self.motifCount["FBLU"].keys()).intersection(set(self.motifCount["FFL"].keys()))))
        # print ("FFL", len(self.motifCount["FFL"].keys()))
        # print ("FBLU", len(self.motifCount["FBLU"].keys()))
        # exit()
        self.motifSets = None

    def countMotifs(self, cs, t):
        for c in cs: self.motifCount[t][c] += 1
    def testMotifs(self):
        self.countMotifs(self.testCX(self.motifCount["C2"].keys()), "C2")
        self.countMotifs(self.testCX(self.motifCount["C3"].keys()), "C3")
        self.countMotifs(self.testCX(self.motifCount["C4"].keys()), "C4")
        self.countMotifs(self.testFFL(self.motifCount["FFL"].keys()), "FFL")
        self.countMotifs(self.testFBLD(self.motifCount["FBLD"].keys()), "FBLD")
        self.countMotifs(self.testFBLU(self.motifCount["FBLU"].keys()), "FBLU")
    def runRadomize(self, count):
        cSet = set()
        for k in range(count):
            random.seed(time.time())
            self.randomizeGraph(0.45)
            self.testMotifs()
            self.randCount += 1
    def randomizeGraph(self, kk):
        nlinkSet = set()
        swCount = 0
        eN = self.pgLinkCount
        iters = int(eN * kk)
        for k in range(iters):
            i1 = random.randint(0, eN - 1)
            src = True
            while src:
                i2 = random.randint(0, eN - 2)
                if i2 == i1: i2 += 1
                e1 = self.links[i1]
                e2 = self.links[i2]
                eid12 = self.edgeDId(e1[0], e2[1])
                eid21 = self.edgeDId(e2[0], e1[1])
                if e1[0] != e2[0] and e1[1] != e2[1] and eid12 not in self.linkSet and eid21 not in self.linkSet and eid12 not in nlinkSet and eid21 not in nlinkSet:
                    self.linkSet.remove(self.edgeDId(e1[0], e1[1]))
                    # if self.edgeDId(e2[0], e2[1]) not in self.linkSet:
                    #     print ("not in linkSet", i1, i2, eN, e2, len(self.linkSet), swCount)
                    self.linkSet.remove(self.edgeDId(e2[0], e2[1]))
                    v = e1[0]
                    e1[0] = e2[0]
                    e2[0] = v
                    self.linkSet.add(eid12)
                    self.linkSet.add(eid21)
                    nlinkSet.add(eid12)
                    nlinkSet.add(eid12)
                    swCount += 1
                    src = False
        self.linkSet = set()
        for e in self.links: self.linkSet.add(self.edgeDId(e[0], e[1]))
    def saveRandomize(self):
        self.randomizedData = {}
        self.randomizedData["randCount"] = self.randCount
        self.randomizedData["motifCount"] = self.motifCount
        self.randomizedData["proteinList"] = self.proteinList
        Utils.utilsReader.saveJson(self.template["fileOutput"] , self.randomizedData)
    def readRandomize(self):
        randomizedData = Utils.utilsReader.readJson(self.template["fileOutput"])
        self.randCount = randomizedData["randCount"]
        self.motifCount = randomizedData["motifCount"]
        self.links = Utils.utilsReader.readJson(self.template["fileInput"])
        self.proteinList = Utils.utilsReader.readJson(self.template["fileInputProteins"])
    def executeRandomize(self):
        if os.path.isfile(self.template["fileOutput"]):
            self.readRandomize()
            if self.randCount >= self.template["randomizeCount"]:
                print ("Randomization has already been completed", self.randCount)
                exit()
            else:
                print ("Continue randomize from", self.randCount)
                self.updateCand()
        else:
            self.initData()
            self.initRandomize()
        while self.randCount < self.template["randomizeCount"]:
            self.runRadomize(self.template["saveAfterCount"])
            self.saveRandomize()
            print ("After save rand count", self.randCount)
        self.saveRandomize()
        print ("Randomization is completed", self.randCount)