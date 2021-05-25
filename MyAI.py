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
        self.__X = startX  # These two variables keeps track of the coordinate of current tile we choose to uncover 
        self.__Y = startY
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
    
    def splitFrontierC(self, FC, FU) -> "list[[tuple]]":
        allFrontiers = []
        allNeighbors = []
        #print("FC:", FC)
        for f in FC:
            #print("f:", f)
            newNeighbors = {nb for nb in self.getValidNeighbors(f[0], f[1]) if type(self.__board[nb[0]][nb[1]][0]) == int}
            newFrontiers = {f}
            for i in range(len(allFrontiers)):
                #print("The intersection between tile", allFrontiers[i], "and newFrontiers", newFrontiers, "is:", allNeighbors[i].intersection(newNeighbors))
                if allNeighbors[i].intersection(newNeighbors) != set():
                    newFrontiers = newFrontiers.union(allFrontiers[i])
                    allFrontiers[i] = set()
                    newNeighbors = newNeighbors.union(allNeighbors[i])
                    allNeighbors[i] = set()
            while set() in allFrontiers:
                allFrontiers.remove(set())
                allNeighbors.remove(set())
            #print("new:", newFrontiers)
            allFrontiers.append(newFrontiers)
            allNeighbors.append(newNeighbors)
        #print(allFrontiers)
        for i in range(len(allFrontiers)):
            allFrontiers[i] = list(allFrontiers[i])
        
        return allFrontiers                

    def splitFrontierU(self, FC, FU) -> "list[[tuple]]":
        out = [[] for i in range(len(FC))]
        for f in FU:
            for i in range(len(FC)):
                if self.getValidNeighbors(f[0], f[1]).intersection(set(FC[i])) != set():
                    out[i].append(f)
        return out


    def analyzeFrontier(self, FU, FC) -> "list[float]":
        # This function analyzes a given pair of "covered" and "uncovered" frontier and the board situation to give a list of floats for
        # each coordinate in the "uncovered" frontier. Each float represents the probability that there is a mine at the coordinate
        # e.g. FU = [(2, 2), (2, 3), (3, 2), (4, 2)], res = [0.0, 1.0, 0.6, 0.25]
        def generateAllCases(n) -> "list[str]":
            if n > 20:
                print("large n:", n)
            remainingMines = self.__numOfMinesLeft
            out = []
            while remainingMines >= 1:
                idxOfOnes = combinations(range(n), remainingMines)
                for idxs in idxOfOnes:
                    combination = ""
                    for i in range(n):
                        if i in idxs:
                            combination += "1"
                        else:
                            combination += "0"
                    out.append(combination)
                remainingMines -= 1
                #print("all cases:", out)
            return out
            
        relationDict = {}
        for tile in FU:
            inter = self.getValidNeighbors(tile[0], tile[1]).intersection(set(FC))
            if inter != set():
                relationDict[tile] = set(inter)
        #print("relation:", relationDict)
        def isValidComb(comb: str) -> bool:
            for k, v in relationDict.items():
            # k -> a uncovered tile in the uncovered frontier
            # v -> a set containing all the covered tiles in the covered frontier that are next to k
                effectiveLabel = self.__board[k[0]][k[1]][1]
                #print("comb:", comb, "k:", k, "v", v, "effectiveLabel:", effectiveLabel)
                #print("FC", FC, "comb", comb)
                for ct in v:
                    i = FC.index(ct)
                    if comb[i] == "1":
                        effectiveLabel -= 1
                        if effectiveLabel < 0:
                            return False
                if effectiveLabel != 0:
                    return False
            return True        

        def generateReport(inputCombs: "list[str]") -> "list[float]":
            rep = [0 for _ in range(len(FC))]
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
        for tile in FU:
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
        if easyReport == 1:
            for tile in FC:
                if tile in zeros:
                    self.__toUncover.append(tile)
                elif tile in ones:
                    self.updateBoard(tile[0], tile[1], "M")
            return []
        else:
            allCombs = generateAllCases(len(FC))
            validCombs = [comb for comb in allCombs if isValidComb(comb)]
            #print("validCombs:", validCombs)
            finalReport = generateReport(validCombs)
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
            #print("111")
            frontierU = self.getFrontierU()
            #print("111")
            frontierC = self.getFrontierC(frontierU)
            #print("111")
            frontierC = self.splitFrontierC(frontierC, frontierU)
            #print("111")
            frontierU = self.splitFrontierU(frontierC, frontierU)
            #print("frontierU:", frontierU)
            #print("frontierC:", frontierC)
            report = [self.analyzeFrontier(frontierU[i], frontierC[i]) for i in range(len(frontierC))]
            #print("report:", report)
            if [] in report:
                return
            self.handleReport(report, frontierC)

        self.updateBoard(self.__X, self.__Y, number) # After getting the face number back, we update the board (knowledge base) with this new knowledge
        while len(self.__toUncover) == 0:
            self.printBoardInfo()
            #print()
            if self.__numOfMinesLeft == 0:
                timeUsed = time.time() - self.__time
                print("Result 1: It taked", timeUsed, "seconds to solve this board.")
                return Action(AI.Action.LEAVE)
            analyzeBoard()
            if self.__numOfMinesLeft == 0:
                for i in range(self.__rowDimension):
                    for j in range(self.__colDimension):
                        if self.__board[i][j][0] == "*" and (i, j) not in self.__toUncover:
                            self.__toUncover.append((i, j))
                if len(self.__toUncover) == 0:
                    timeUsed = time.time() - self.__time
                    print("Result 2: It taked", timeUsed, "seconds to solve this board.")
                    return Action(AI.Action.LEAVE)
        nextTile = self.__toUncover.pop(0)
        #self.printBoardInfo()
        #print()
        print("next move:", nextTile)
        self.__X, self.__Y = nextTile[0], nextTile[1]
        #print("It has taken", time.time() - self.__time, "seconds until now.")
        return Action(AI.Action.UNCOVER, nextTile[0], nextTile[1])


        


        
        
            