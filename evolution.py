import numpy as np
import math
import matplotlib.pyplot as plt
from random import randint
from operator import itemgetter
from itertools import islice
import heapq

class Evolution:
    '''
    遗传算法
    默认用于解决迷宫问题，同时可配置于求解其他问题
    usage:
        e = Evolution()
        e.efficacy_evaluate()
    '''
    def __init__(self, evaluate=None, geneNum=25, genNum=100, genCount=1000):
        '''
        evaluate 评估函数
        geneNum 基因数目
        genNum 初始个体数
        genCount 迭代次数
        '''
        if evaluate is None:
            mm = Evolution.Map()
            self.evaluate = mm.evaluate
        else:
            self.evaluate = evaluate
        self.geneNum = geneNum
        self.genNum = genNum
        self.genCount = genCount
        pass
    def getRandom(self):
        return list(map(lambda x:randint(0,3),[0]*self.geneNum))
    def getRandomGen(self):
        return list(map(lambda x:[(self.evaluate(x),x) for x in [self.getRandom()]][0], [0]*self.genNum))
    def cross(self, gene1, gene2):
        index = randint(1, self.geneNum-2)
        gene1, gene2 = gene1[:index]+gene2[index:], gene1[index:]+gene2[:index]
        return gene1, gene2
    def evolve_once(self, gen):
        gen.sort(key=itemgetter(0))
        newgen = []
        for x,y in zip(islice(gen,None,None,2), islice(gen,1,None,2)):
            newgene1, newgene2 = self.cross(x[1],y[1])
            newgen.append((self.evaluate(newgene1),newgene1))
            newgen.append((self.evaluate(newgene2),newgene2))
        gen.extend(newgen)
        gen.sort(key=itemgetter(0))
        return gen[:self.genNum]
    def evolve(self):
        gen = e.getRandomGen()
        for x in range(self.genCount):
            gen = self.evolve_once(gen)
            if gen[0][0]==0:
                break
        return (x, gen)
    def efficacy_evaluate(self, geneNum=25, genNum=100, evalNum=100, maxGenCount=1000):
        '''
        对遗传算法进行效率评估
        geneNum: 基因数目
        genNum: 初始个体数
        evalNum: 运行算法的次数
        maxGenCount: 每次运行算法的最大迭代数
        '''
        self.genCount = maxGenCount
        self.geneNum = geneNum
        self.genNum = genNum

        cnt = 0
        genCnt = 0
        failedTimes = 0
        while cnt<evalNum:
            cnt += 1
            print('\r第{0}次运行'.format(cnt), end='')
            times,_ = self.evolve()
            if times>=maxGenCount-1:
                failedTimes += 1
            else:
                genCnt += times
        print('\r评估完成。共运行{0}次，成功{3}次，失败{2}次，找到可行解平均需要{1:.1f}次迭代。'
              .format(cnt, genCnt/(cnt-failedTimes), failedTimes, cnt-failedTimes))
    
    class Map:
        def __init__(self):
            self.map = [1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1,
                1,0,1,0,0, 0,0,0,1,1, 1,0,0,0,1,
                8,0,0,0,0, 0,0,0,1,1, 1,0,0,0,1,
                1,0,0,0,1, 1,1,0,0,1, 0,0,0,0,1,
                1,0,0,0,1, 1,1,0,0,0, 0,0,1,0,1,

                1,1,0,0,1, 1,1,0,0,0, 0,0,1,0,1,
                1,0,0,0,0, 1,0,0,0,0, 1,1,1,0,1,
                1,0,1,1,0, 0,0,1,0,0, 0,0,0,0,5,
                1,0,1,1,0, 0,0,1,0,0, 0,0,0,0,1,
                1,1,1,1,1, 1,1,1,1,1, 1,1,1,1,1]
            self.width = 15
            self.height = 10
            self.start = [7,14]
            self.end = [2,0]
        def get(self,x,y):
            if x<0 or x>=self.height or y<0 or y>=self.width:
                #print(x,y,'return None')
                return None
            return self.map[x*self.width+y]
        def evaluate(self, gene):
            return self.distance(self.end, self.getDest(gene))
        def distance(self, posa, posb):
            dx = posa[0]-posb[0]
            dy = posa[1]-posb[1]
            return math.sqrt(dx**2+dy**2)
        def getDest(self, gene):
            dx = [0,0,-1,1]
            dy = [1,-1,0,0]
            pos = self.start[:]
            for x in gene:
                nx = pos[0]+dx[x]
                ny = pos[1]+dy[x]
                val = self.get(nx, ny)
                if val is None or val==1:
                    #print(x,pos,'jump')
                    continue
                elif val==8:
                    return [nx,ny]
                else:
                    #print(x,pos,'continue')
                    pos[0] = nx
                    pos[1] = ny
            return pos


if __name__ == '__main__':
    e = Evolution()
    e.efficacy_evaluate()