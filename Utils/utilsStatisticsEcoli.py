import os.path, math, random, time
import Utils.utilsReader
# import Utils.DGDVS
class Dgdvs:
    def __init__(self, orbNum):
        self.no = orbNum
        # self.no = 13
        q = [[0,1], [2, 4, 5, 6, 7, 8, 10, 11], [3, 9, 12, 13, 16, 17, 20, 21, 24, 25, 28, 29, 30, 31, 32, 33, 36, 60, 64, 68, 72, 118, 120, 125, 128], [19, 22, 23, 27, 35, 38, 44, 45, 47, 48, 52, 56, 57, 61, 63, 66, 67, 70, 73, 74, 76, 77, 78, 80, 85, 87, 91, 92, 95, 96, 100, 105, 113, 115], [14, 15, 18, 26, 34, 37, 42, 43, 49, 50, 53, 54, 58, 59, 62, 65, 69, 71, 88, 90, 93, 94, 99, 103, 107, 111, 121, 122, 126, 127], [39, 40, 41, 46, 51, 55, 75, 79, 81, 82, 83, 84, 86, 89, 97, 98, 101, 102, 104, 106, 108, 109, 110, 112, 114, 116,117, 119, 123, 124]]
        oq = {}
        for i in range(6):
            for o in q[i]: oq[o] = i 
        self.w = []
        self.ws = 0 
        for i in range(self.no):
            self.w.append(1 - math.log(oq[i] + 1) / math.log(self.no))
            self.ws += self.w[i]
    def getValue(self, u, v):
        sd = 0.0
        sw = 0.0
        for i in range(self.no):
            r = self.w[i] * abs(math.log(u[i] + 1) - math.log(v[i] + 1)) / math.log(max(u[i], v[i]) + 2)
            sd += r
        return 1 - sd / self.ws
class utilsGraphlets:
    def __init__(self, templStats):
        self.delim = "#"
        self.orbitNumb = 13
        self.templStats = templStats
    def eidD(self, v, w):
        return w + (v << self.dim)
    def eidDDir(self, v, w, dir):
        if dir == 0: return self.eidD(v, w)
        if dir == 1: return self.eidD(w, v)
    def opositdir(self, dir):
        return 3 - dir
    def eidU(self, v, w):
        if v > w: return self.eidD(w, v)
        return self.eidD(v, w)
    def eidFrom(self, eid):
        return eid >> self.dim
    def eidTo(self, eid):
        return eid & self.maskdim
    def eidEdge(self, eid):
        return [self.eidFrom(eid), self.eidTo(eid)]
    def oposite(self, eid):
        e = self.eidEdge(eid)
        return self.eidD(e[1], e[0])
    def undirect(self, eid):
        e = self.eidEdge(eid)
        return self.eidU(e[0], e[1])
    def utSet(self, v, w): 
        eid = self.eidU(v, w)
        if eid not in self.utadj: return set()
        return self.utadj[eid]
    def dtSet(self, v, w, i):
        eid = self.eidU(v, w)
        if eid not in self.dtadj: return set()
        return self.dtadj[eid][i]
    def setEdges(self, edges): 
        nodes = set()
        for e in edges: 
            nodes.add(e[0])
            nodes.add(e[1])
        self.nodes = sorted(list(nodes))
        self.N = len(self.nodes)
        self.M = len(edges)
        self.orbit = [[0 for j in range(self.orbitNumb)] for i in range(self.N)]
        self.graphlets = {}

        self.nodeInds = {}
        n = len(self.nodes)
        for i in range(n): self.nodeInds[self.nodes[i]] = i
        self.dim = 1
        while 1 << self.dim < n + 1: self.dim += 1
        self.maskdim = (1 << self.dim) - 1
        self.adj = {}
        self.inc = {}
        self.uadj = {}
        self.uinc = {}
        self.dtinc = {}
        self.utinc = {}
        for v in range(n): 
            self.adj[v] = [set(), set()]
            self.inc[v] = [set(), set()]
            self.uadj[v] = set()
            self.uinc[v] = set()
            self.dtinc[v] = [set(), set(), set(), set()]
            self.utinc[v] = set()
        self.dedgeSet = set()
        self.uedgeSet = set()
        for es in edges:
            e = [self.nodeInds[es[0]], self.nodeInds[es[1]]]
            eid = self.eidD(e[0], e[1])
            self.dedgeSet.add(eid)
            self.adj[e[0]][0].add(e[1])
            self.adj[e[1]][1].add(e[0])
            self.inc[e[0]][0].add(eid)
            self.inc[e[1]][1].add(eid)
            ueid = self.eidU(e[0], e[1])
            self.uedgeSet.add(ueid)
            self.uadj[e[0]].add(e[1])
            self.uadj[e[1]].add(e[0])
            self.uinc[e[0]].add(ueid)
            self.uinc[e[1]].add(ueid)
    def countPairOrbit_G0(self, p):
        k = self.pairOrbitLim[0]
        for i in [0, 1]:
            for j in [0, 1]: 
                self.pairOrbit[i + j * 2] = len(self.adj[p[i]][j] - {p[1 - i]})
            if (p[1 - i] in self.adj[p[i]][0]): 
                self.pairOrbit[4 + i] += 1
            print(self.adj[p[0]][i], "&", self.adj[p[1]][i], "-", {p[0], p[1]}, (self.adj[p[0]][i] & self.adj[p[1]][i] - {p[0], p[1]}))
            self.pairOrbit[6 + i] = len(self.adj[p[0]][i] & self.adj[p[1]][i] - {p[0], p[1]})
    def countPairOrbit_G1(self, p):
        k = self.pairOrbitLim[1]
        for i in [0, 1]:
            for v in self.adj[p[i]][0] - {p[1 - i]}:
                for w in self.adj[v][0] - {p[0], p[1]}:
                    if self.eidU(p[i], w) not in self.uedgeSet:
                        self.pairOrbit[k + i] += 1
            for v in self.adj[p[i]][1] - {p[1 - i]}:
                for w in self.adj[v][1] - {p[0], p[1]}:
                    if self.eidU(p[i], w) not in self.uedgeSet:
                        self.pairOrbit[k + 2 + i] += 1
                for w in self.adj[p[i]][0] - {p[1 - i]}:
                    if self.eidU(v, w) not in self.uedgeSet:
                        self.pairOrbit[k + 4 + i] += 1
                if p[1 - i] in self.adj[p[i]][0]:
                    if self.eidU(v, p[1 - i]) not in self.uedgeSet:
                        self.pairOrbit[k + 8 + i] += 1
            if p[1 - i] in self.adj[p[i]][0]:
                for w in self.adj[p[1 - i]][0] - {p[i]}:
                    if self.eidU(p[i], w) not in self.uedgeSet:
                        self.pairOrbit[k + 6 + i] += 1
            for v in self.adj[p[0]][i] & self.adj[p[1]][i]:
                for w in self.adj[v][i] - {p[0], p[1]}:
                    if self.eidU(p[0], w) not in self.uedgeSet and self.eidU(p[1], w) not in self.uedgeSet:
                        self.pairOrbit[k + 10 + i] += 1

        for v in self.adj[p[0]][0] & self.adj[p[1]][0]:
            for w in self.adj[v][0] - {p[0], p[1]}:
                if self.eidU(p[0], w) not in self.uedgeSet and self.eidU(p[1], w) not in self.uedgeSet:
                    self.pairOrbit[k + 12] += 1
        for v in self.adj[p[0]][1] & self.adj[p[1]][1]:
            for w in self.adj[v][1] - {p[0], p[1]}:
                if self.eidU(w, p[0]) not in self.uedgeSet and self.eidU(w, p[1]) not in self.uedgeSet:
                    self.pairOrbit[k + 13] += 1
        for v in self.adj[p[0]][1] & self.adj[p[1]][1]:
            for w in self.adj[p[0]][0] & self.adj[p[1]][0]:
                if self.eidU(v, w) not in self.uedgeSet:
                    self.pairOrbit[k + 14] += 1
    def countPairOrbit_G2_3(self, p, m):
        k = self.pairOrbitLim[2 + m]
        for i in [0, 1]:
            for v in self.adj[p[i]][1 - m] - {p[1 - i]}:
                for w in self.adj[v][m] - {p[0], p[1]}:
                    if self.eidU(p[i], w) not in self.uedgeSet:
                        self.pairOrbit[k + i] += 1
            for v in self.adj[p[i]][m] - {p[1 - i]}:
                for w in self.adj[p[i]][m] - {p[1 - i], v}:
                    if self.eidU(v, w) not in self.uedgeSet:
                        self.pairOrbit[k + 2 + i] += 1
                if p[1 - i] in self.adj[p[i]][m]:
                    if self.eidU(v, p[1 - i]) not in self.uedgeSet:
                        self.pairOrbit[k + 4 + i] += 1
        for v in self.adj[p[0]][1 - m] & self.adj[p[1]][1 - m]:
            for w in self.adj[v][m] - {p[0], p[1]}:
                if self.eidU(p[0], w) not in self.uedgeSet and self.eidU(p[1], w) not in self.uedgeSet:
                    self.pairOrbit[k + 6] += 1
        a = list(self.adj[p[0]][m] & self.adj[p[1]][m])
        for i in range(len(a) - 1):
            for j in range(i + 1, len(a)):
                if self.eidU(a[i], a[j]) not in self.uedgeSet:
                    self.pairOrbit[k + 7] += 1
    def countPairOrbit_G4(self, p):
        k = self.pairOrbitLim[4]
        for i in [0, 1]:
            for v in self.adj[p[i]][0] - {p[1 - i]}:
                self.pairOrbit[k + i] = len(self.adj[v][0] & self.adj[p[i]][1] - {p[1 - i]})
            if p[1 - i] in self.adj[p[i]][0]:
                self.pairOrbit[k + 2 + i] = len(self.adj[p[1 - i]][0] & self.adj[p[i]][1])
        for v in self.adj[p[0]][1] & self.adj[p[1]][1]:
            for w in self.adj[p[0]][0] & self.adj[p[1]][0]:
                self.pairOrbit[k + 4] = len(self.adj[p[0]][0] & self.adj[p[1]][0] & self.adj[v][1])
    def countPairOrbit_G5(self, p):
        k = self.pairOrbitLim[5]
        for i in [0, 1]:
            for v in self.adj[p[i]][1] - {p[1 - i]}:
                self.pairOrbit[k + i] = len(self.adj[v][1] & self.adj[p[i]][1] - {p[1 - i]})
            for v in self.adj[p[i]][0] - {p[1 - i]}:
                self.pairOrbit[k + 2 + i] = len(self.adj[v][1] & self.adj[p[i]][1] - {p[1 - i]})
            for v in self.adj[p[i]][0] - {p[1 - i]}:
                self.pairOrbit[k + 4 + i] = len(self.adj[v][0] & self.adj[p[i]][0] - {p[1 - i]})
            if p[1 - i] in self.adj[p[i]][1]:
                self.pairOrbit[k + 6 + i] = len(self.adj[p[1 - i]][1] & self.adj[p[i]][1])
            if p[1 - i] in self.adj[p[i]][1]:
                self.pairOrbit[k + 8 + i] = len(self.adj[p[1 - i]][0] & self.adj[p[i]][1])
            if p[1 - i] in self.adj[p[i]][1]:
                self.pairOrbit[k + 10 + i] = len(self.adj[p[1 - i]][0] & self.adj[p[i]][0])
        for v in self.adj[p[i]][1] & self.adj[p[1 - i]][1]:
            self.pairOrbit[k + 12] = len(self.adj[v][1] & self.adj[p[i]][1] & self.adj[p[1 - i]][1])
        for v in self.adj[p[i]][0] & self.adj[p[1 - i]][0]:
            self.pairOrbit[k + 13] = len(self.adj[v][0] & self.adj[p[i]][0] & self.adj[p[1 - i]][0])
        for v in self.adj[p[i]][0] & self.adj[p[1 - i]][0]:
            self.pairOrbit[k + 14] = len(self.adj[v][1] & self.adj[p[i]][1] & self.adj[p[1 - i]][1])
    # ############################################################################################################
    def countGenPairOrbits(self):
        for p in self.ohn:
            if p[0] not in self.nodeInds or p[1] not in self.nodeInds:
                print ("\n***** pair is not in gene list:", p)
                exit()
        self.pairOrbitLim = [0, 8, 23, 31, 39, 44, 59]
        self.pairNumb = len(self.ohn)
        self.pairOrbits = []
        self.pairGraphlets = {}
        for pind in range(len(self.ohn)):
            p = self.ohn[pind]
            self.pairOrbit = [0 for j in range(self.pairOrbitLim[6])]
            v = [self.nodeInds[p[0]], self.nodeInds[p[1]]]
            self.countPairOrbit_G0(v)
            self.countPairOrbit_G1(v)
            self.countPairOrbit_G2_3(v, 0)
            self.countPairOrbit_G2_3(v, 1)
            self.countPairOrbit_G4(v)
            self.countPairOrbit_G5(v)
            self.pairOrbits.append(self.pairOrbit)
    # ############################################################################################################
    def getResultsWithDGDVS(self, oNumb):
        dg = Dgdvs(oNumb)
        header = ["ohnolog", "other", "value", "type", "key", "DGDVS"]
        for i in range(len(self.pairOrbitLim) - 1):
            for j in range(self.pairOrbitLim[i + 1] - self.pairOrbitLim[i]):
                header.append("PG" + str(i) + "-" + str(j))
        rez = [header]
        orbInds = [0, 2, 8, 12, 10, 23, 25, 31, 33, 39, 44, 48, 46]
        orbits = [[0 for i in range(oNumb)], [0 for i in range(oNumb)]]
        for i in range(len(self.pairOrbits)):
            for k in range(oNumb):
                orbits[0][k] = self.pairOrbits[i][orbInds[k]]
                orbits[1][k] = self.pairOrbits[i][orbInds[k] + 1]
            dgdvs = dg.getValue(orbits[0], orbits[1])
            rez.append(self.ohn[i] + [dgdvs] + self.pairOrbits[i])
        return rez
    def finalStats(self):
        self.ohn = Utils.utilsReader.readJson(self.templStats["fileInputGenePairs"])
        for d in self.ohn: d += ["", "", ""]
        self.setEdges(Utils.utilsReader.readJson(self.templStats["fileInputLinks"]))
        self.countGenPairOrbits()
        Utils.utilsReader.saveCsv(self.templStats["fileOutputStats"], self.getResultsWithDGDVS(13))
