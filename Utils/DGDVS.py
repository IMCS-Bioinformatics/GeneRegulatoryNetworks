import math

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