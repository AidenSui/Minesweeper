# ==============================CS-199==================================
# FILE:         MyAI.py
#
# AUTHOR:       Justin Chung
#
# DESCRIPTION:  This file contains the MyAI class. You will implement your
#                agent in this file. You will write the 'getAction' function,
#                the constructor, and any additional helper functions.
#
# NOTES:        MyAI inherits from the abstract AI class in AI.py.
#
#               - DO NOT MAKE CHANGES TO THIS FILE.
# ==============================CS-199==================================

from AI import AI
from Action import Action
from collections import defaultdict
from itertools import combinations
import time


class MyAI( AI ):
    def __init__(self, rowDimension, colDimension, totalMines, startX, startY):
        self.__rowDimension = rowDimension
        self.__colDimension = colDimension
        self.__board = [["*" for j in range(colDimension)] for i in range(rowDimension)]
        # We use a list of length three to store three informations for each tile - 1. face number 2. effective label 3. number of adjacent uncovered tiles
        # 1. there are three candidates for a tile's face number: * -> uncovered; M -> Marked; int -> covered
        # 2. a tile's effective label = face number - number of adjacent marked tiles (-1 for uncovered and marked tiles)
        # 3. the number of adjacent uncovered tiles stores how many adjacent tiles have face number "*" 
        # Initialize the board 
        for i in range(rowDimension):
            for j in range(colDimension):
                if i == 0 or i == rowDimension - 1 or j == 0 or j == colDimension - 1:
                    if (i == 0 and j == 0) or (i == 0 and j == colDimension - 1) or (i == rowDimension - 1 and j == 0) or (i == rowDimension - 1 and j == colDimension - 1):
                        self.__board[i][j] = ["*", -1, 3]  # A corner tile has 3 neighbors
                    else:
                        self.__board[i][j] = ["*", -1, 5]  # An edge tile has 5 neighbors
                else:
                    self.__board[i][j] = ["*", -1, 8] # A core tile has 8 neighbors
        self.__toUncover = [] # This variable stores all the coordinates we choose to uncover
        self.__numOfMinesLeft = totalMines # This variable keeps track of how many mines are left, which is equal to totalMines initially
        self.__X = startY  # These two variables keeps track of the coordinate of current tile we choose to uncover 
        self.__Y = startX
        self.__time = time.time()

    # helper functions
    def getValidNeighbors(self, x, y) -> "list[tuple]":
        def isValid(x, y) -> bool:
            if x <= self.__rowDimension - 1 and x >= 0 and y <= self.__colDimension - 1 and y >= 0:
                return True
            return False
        allNeighbors = [(x - 1, y - 1), 
                     (x    , y - 1), 
                     (x + 1, y - 1), 
                     (x - 1, y), 
                     (x + 1, y), 
                     (x - 1, y + 1), 
                     (x    , y + 1), 
                     (x + 1, y + 1) ]
        return {n for n in allNeighbors if isValid(n[0], n[1]) }

    def updateBoard(self, x, y, label) -> None:
        # This function updates the variable self.__board based on one given information - the label of the tile at (x, y)
        #print("Now updating board with info given on tile", x, y)
        self.__board[x][y][0] = label
        neighbors = self.getValidNeighbors(x, y)
        # print("The neighbors of the current tile is", neighbors)
        mineCount = 0
        for coord in neighbors:           
            coord_x = coord[0]
            coord_y = coord[1]
            if label != "M":
            # if the current tile is not a mine
                if self.__board[coord_x][coord_y][0] == "M":
                    mineCount += 1 # increment the count of mines in neighborhood by 1
            else:
            # if the current tile is a mine
                if self.__board[coord_x][coord_y][0] != "M" and self.__board[coord_x][coord_y][0] != "*":
                    self.__board[coord_x][coord_y][1] -= 1 # decrement the effective label of a neighbor by 1
                
            # print("Decrementing the amount of unknown neighbors of", coord, "by 1")           
            self.__board[coord_x][coord_y][2] -= 1 # decrement the num of uncovered neighbor of a neighbor by 1
        if type(label) == int:
            self.__board[x][y][1] = label - mineCount # decrement the effective label by the count of flagged mines in the neighborhood
        else:
            self.__numOfMinesLeft -= 1


    def getFrontierU(self) -> "list[tuple]":
        # This function gives a list of all the current "Uncovered" frontiers based on the current board condition
        # Return a list of frontiers, each frontier is a list of tile coordinates. e.g. frontier = [(1, 3), (0, 4), (0, 5), (1, 6)]
              
        allFrontiers = []
        for i in range(self.__rowDimension):
            for j in range(self.__colDimension):
                if type(self.__board[i][j][0]) == int and self.__board[i][j][2] != 0:
                    allFrontiers.append((i, j))

        return allFrontiers

    def getFrontierC(self, FU) -> "list[tuple]":
        # This function computes an uncovered frontier based on the current board condition and a given covered frontier
        # A given covered frontier is needed as multiple frontiers can exist at the same time
        frontier = []
        for f in FU:
            for coord in self.getValidNeighbors(f[0], f[1]):
                if self.__board[coord[0]][coord[1]][0] == "*" and coord not in frontier:
                    frontier.append(coord)
        return frontier
    
    def splitFrontier(self, FC, FU) -> "splittedFC, splittedFU":
        # returns two splitted frontiers
        # A covered frontier is splitted once its length goes over 15 or has no overlapping tiles with the others
        # A covered frontier will be splitted by looping through the following steps:
        # 1. A new sub-frontier a initiated by finding a head and append it to a empty list
        # 2. We look for its neighbors that are also in the covered frontier and add them to the sub-frontier
        # 3. If no more connected tiles can be found or the length exceed 15, we append the sub-frontier to splittedFC and initiate a new sub-frontier
        #    the corresponding uncovered frontier is formed by joining all the neighbors for tiles in our covered sub-frontier, and is append to the splittedFU
        # 4. If no more element in FC, we stop the loop and return splittedFC and splitted FU
        
        def getUncoveredNeighborDict(FC, FU) -> "{tuple: set(tuple)}":
            # return a dictionary that stores all its neighbors in FU for a tile in FC
            # key: tile in FC
            # value: the set of all the valid neighbors of the key that are in FU
            unbd = dict()
            for t in FC:
                unbd[t] = self.getValidNeighbors(t[0], t[1]).intersection(set(FU))
            return unbd

        def getCoveredNeighborDict(uncoveredNeighborDict) -> "{tuple: set(tuple)}":
            # return a dictionary that stores all its neighbors in FC that have mutual neighbors in FU for a tile in FC
            # key: tile in FC
            # value: the set of all the neighbors of the key that are in FC and have mutual neighbors with the key tile
            cnbd = dict()
            for tile, unb in uncoveredNeighborDict.items():
                cnb = set()
                for k, v in uncoveredNeighborDict.items():
                    if tile != k and unb.intersection(v) != set():
                        cnb.add(k)
                cnbd[tile] = cnb
            return cnbd

        def getHead(coveredNeighborDict) -> tuple:
            # return the head of a frontier
            # the head is the tile with the fewest neighbors in the frontier
            maxLen = 8
            head = (-1, -1)
            for tile, cnb in coveredNeighborDict.items():
                if len(cnb) < maxLen:
                    head = tile
            return head

        def getSplittedFU(frontier: "list[tuple]") -> "list[tuple]":
            # this function get the corresponding covered sub-frontier from a given uncovered sub-frontier
            # note that there might be duplicated tiles for different covered sub-frontier and this wouldn't be problematic 
            # because the uncovered frontier is only used to check for condition satisfaction
            out = set()
            for f in frontier:
                out = out.union(unbd[f])
            return list(out)

        def addFriends(someFriends, frontier):
            for friend in someFriends:
                if friend not in frontier:
                    #print("removing friend:", friend)
                    FC.remove(friend)
                    frontier.append(friend)
                    for tile in cnbd[friend]:
                        #print("removing", friend, "from dict[" + str(tile) + "]" )
                        cnbd[tile].discard(friend)
            #print(cnbd)

        splittedFC = []
        splittedFU = []
        unbd = getUncoveredNeighborDict(FC, FU)
        #print("unbd:", unbd)
        cnbd = getCoveredNeighborDict(unbd)
        #print("cnbd:", cnbd)
        head = getHead(cnbd)
        frontierStack = [head]
        FC.remove(head)
        idx = 0
        while True:
            #print("head:", head)
            #print("idx:", idx)
            friends = cnbd[head]
            closestFriends = {(head[0], head[1] + 1), (head[0], head[1] - 1), (head[0] + 1, head[1]), (head[0] - 1, head[1]) }.intersection(friends)
            #print("closestFriends:",closestFriends)
            secondClosestFriends = {(head[0] + 1, head[1] + 1), (head[0] - 1, head[1] - 1), (head[0] + 1, head[1] - 1), (head[0] - 1, head[1] + 1) }.intersection(friends)
            #print("secondClosestFriends:",secondClosestFriends)
            leastClosestFriends = friends - closestFriends - secondClosestFriends
            #print("leastClosestFriends:",leastClosestFriends)
            addFriends(closestFriends, frontierStack)
            addFriends(secondClosestFriends, frontierStack)
            addFriends(leastClosestFriends, frontierStack)
            #print("frontierStack:", frontierStack)
            #print("FC:", FC)
            if len(FC) == 0:
                splittedFC.append(frontierStack)
                splittedFU.append(getSplittedFU(frontierStack))
                break
            if len(frontierStack) >= 15 or idx == len(frontierStack) - 1:
                #print("splitting the frontier ...")
                splittedFC.append(frontierStack)
                splittedFU.append(getSplittedFU(frontierStack))
                for tile in frontierStack:
                    del unbd[tile]
                    del cnbd[tile]
                head = getHead(cnbd)
                #print("new head:", head)
                frontierStack = [head]
                FC.remove(head)
                idx = 0
                continue
            idx += 1
            head = frontierStack[idx]
        
        #print("Splitted froontierC:", splittedFC)
        #print("Splitted froontierU:", splittedFU)
        return splittedFC, splittedFU
               


    def analyzeFrontier(self, FU, FC) -> "list[float]":
        # This function analyzes a given pair of "covered" and "uncovered" frontier and the board situation to give a list of floats for
        # each coordinate in the "uncovered" frontier. Each float represents the probability that there is a mine at the coordinate
        # e.g. FU = [(2, 2), (2, 3), (3, 2), (4, 2)], res = [0.0, 1.0, 0.6, 0.25]
        #print("FC", FC)
        #print("FU", FU)
        def isValidComb(comb, fc) -> bool:
            #print(relationDict)
            #print("checking", comb)
            if len(fc) == unknownCount:
                if comb.count("1") != self.__numOfMinesLeft:
                    return False
            for k in relationDict.keys():
            # k -> a uncovered tile in the uncovered frontier
            # v -> a set containing all the covered tiles in the covered frontier that are next to k
                effectiveLabel = self.__board[k[0]][k[1]][1]
                
                #print("effectiveLabel for", k, "is", effectiveLabel)
                #print("comb:", comb, "k:", k, "v", v, "effectiveLabel:", effectiveLabel)
                #print("FC", FC, "comb", comb)
                #print("neighbors:", relationDict[k])
                for ct in relationDict[k]:
                    #print("ct", ct)
                    
                    i = fc.index(ct)
                    # if len(comb) >= 9 and comb[8] == "1" and comb[9] == "1":
                    #     print("FC:", FC)
                    #     print("i", i)
                    #     print("comb:", comb)
                    #     print("k", k)
                    #     print("ct", ct, "   effectivelabel:", effectiveLabel)
                    #     print()
                    #print("i:", i)
                    
                    if comb[i] == "1":
                        effectiveLabel -= 1
                if effectiveLabel != 0:
                    return False
            
            return True       

        def generateAllCases(fc, fu) -> "list[str]":
            print("fc", fc)
            out = []

            
                # if n >= remainingMines:
                #     idxOfOnes = combinations(range(n), remainingMines)
                #     for idxs in idxOfOnes:
                #         combination = ""
                #         for i in range(n):
                #             if i in idxs:
                #                 combination += "1"
                #             else:
                #                 combination += "0"
                #         #print("combination", combination, "is", isValidComb(combination), "valid")
                #         if isValidComb(combination, relationDict):
                #             out.append(combination)
                # else:
            n = len(fc)
            for i in range(2**n):
                combination = bin(i).replace("0b", "")
                combination = (n - len(combination)) * "0" + combination
                #print("combination", combination, "is", isValidComb(combination), "valid")

                if isValidComb(combination, fc):
                    out.append(combination)

            #print("all cases:", out)
            return out
            
 
         

        def generateReport(inputCombs: "list[str]", fc) -> "list[float]":
            rep = [0 for _ in range(len(fc))]
            for comb in inputCombs:
                for i in range(len(comb)):
                    if comb[i] == "1":
                        rep[i] += 1
            for j in range(len(rep)):
                rep[j] /= len(inputCombs)
            return rep

        

        easyReport = 0
        zeros = []
        ones = []
        for fu in FU:
            for tile in fu:
                if self.__board[tile[0]][tile[1]][1] == 0:
                    easyReport = 1
                    for nb in self.getValidNeighbors(tile[0], tile[1]):
                        if nb not in zeros and self.__board[nb[0]][nb[1]][0] != "M":
                            zeros.append(nb)
                elif self.__board[tile[0]][tile[1]][1] == self.__board[tile[0]][tile[1]][2]:
                    easyReport = 1
                    for nb in self.getValidNeighbors(tile[0], tile[1]):
                        if nb not in ones and self.__board[nb[0]][nb[1]][0] != "M":
                            ones.append(nb)
        #print("easyReport:", easyReport)
        #print("zeros", zeros)
        #print("ones", ones)
        if easyReport == 1:
            for fc in FC:
                for tile in fc:
                    if tile in zeros:
                        self.__toUncover.append(tile)
                    elif tile in ones:
                        self.updateBoard(tile[0], tile[1], "M")
            return "break"
        else:
            
            unknownCount = 0
            for row in self.__board:
                for tile in row:
                    if tile[0] == "*":
                        unknownCount += 1


            finalReport = []
            for i in range(len(FU)):
                relationDict = {}
                for tile in FU[i]:
                    relationDict[tile] = self.getValidNeighbors(tile[0], tile[1]).intersection(set(FC[i]))

                #print("relation:", relationDict)
                allCombs = generateAllCases(FC[i], FU[i])
                print("validCombs:", allCombs)
                finalReport.append(generateReport(allCombs, FC[i])) 
            return finalReport

    def handleReport(self, report, FC) -> None:
        # This function reads the report and does three jobs:
        # 1. Check if there is any tile with 0.0 probability to be a mine. If so, push (append) all of them into self.__toUncover (could try: if
        # self.__numOfMinesLeft == 0, push (append) all uncovered tiles into self.__toUncover).
        # 2. Check if there is any tile with 1.0 probability to be a mine. If so, change all of their labels in self.__board accordingly
        # and subtract self.__numOfMinesLeft by the number of mines detected
        # 3. If there is neither 0.0 nor 1.0 in the report, choose the tile with the lowest probability to explode to uncover
        mineIdx = []
        safeIdx = []
        minProb = 1
        minProbIdx = None
        for i in range(len(report)):
            for j in range(len(report[i])):
                prob = report[i][j]
                if prob == 1.0:
                    mineIdx.append((i, j))
                elif prob == 0.0:
                    safeIdx.append((i, j))
                else:
                    if prob < minProb:
                        minProb = prob
                        minProbIdx = (i, j)
        if len(mineIdx) == 0 and len(safeIdx) == 0:
            safeIdx.append(minProbIdx)
        #print("mineIdx", mineIdx)
        #print("safeIdx", safeIdx)
        for idx in mineIdx:
            self.updateBoard(FC[idx[0]][idx[1]][0], FC[idx[0]][idx[1]][1], "M")
        for idx in safeIdx:
            self.__toUncover.append(FC[idx[0]][idx[1]])

    def printBoardInfo(self) -> None:
        print(self.__numOfMinesLeft, "mines left")
        r = len(self.__board)
        c = len(self.__board[0])
        for row in self.__board:
            for col in row:
                print(col[0], end=" |")
            print()
            print("-"*3*self.__colDimension)
        #print("there are still", self.__numOfMinesLeft, "mines left")
        
        #print("We have", r, "rows and", c, "columns in the board")
        #print("There are", self.__numOfMinesLeft, "mines left in the board")

    def getAction(self, number: int) -> "Action Object":
        #print("stack:", self.__toUncover)
        def analyzeBoard():
            frontierU = self.getFrontierU()    
            frontierC = self.getFrontierC(frontierU)
            #print("frontierU:", frontierU)            
            splittedfrontierC, splittedfrontierU = self.splitFrontier(frontierC, frontierU)
            print("splittedfrontierC:", splittedfrontierC)
            report = self.analyzeFrontier(splittedfrontierU, splittedfrontierC)
            if report == "break":
                return
            print("report:", report)
            self.handleReport(report, splittedfrontierC)

        self.updateBoard(self.__X, self.__Y, number) # After getting the face number back, we update the board (knowledge base) with this new knowledge
        while len(self.__toUncover) == 0:
            self.printBoardInfo()
            print()
            if self.__numOfMinesLeft == 0:
                timeUsed = time.time() - self.__time
                print(timeUsed, "seconds used.")
                return Action(AI.Action.LEAVE)
            analyzeBoard()
            #print("uncover", self.__toUncover)
            if self.__numOfMinesLeft == 0:
                for i in range(self.__rowDimension):
                    for j in range(self.__colDimension):
                        if self.__board[i][j][0] == "*" and (i, j) not in self.__toUncover:
                            self.__toUncover.append((i, j))
                if len(self.__toUncover) == 0:
                    timeUsed = time.time() - self.__time
                    print( timeUsed, "seconds used.")
                    return Action(AI.Action.LEAVE)
        nextTile = self.__toUncover.pop(0)
        #self.printBoardInfo()
        #print()
        #print("next move:", nextTile)
        self.__X, self.__Y = nextTile[0], nextTile[1]
        #print("It has taken", time.time() - self.__time, "seconds until now.")
        return Action(AI.Action.UNCOVER, nextTile[1], nextTile[0])


        


        
        
            