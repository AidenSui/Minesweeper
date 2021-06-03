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
import time


class MyAI( AI ):
    
    class AllCombs:
        def __init__(self):
            self._frontier = []
            self._comb_list = []
            self._look_up = dict()


        def addNewTile(self, x, y):
            def addAtTheEnd(comb_list):
                if len(comb_list) == 0:
                    comb_list.append("0")
                    comb_list.append("1")
                    return

                new_elements = []
                for comb in comb_list:
                    new_elements.append(comb + "0")
                    new_elements.append(comb + "1")
                for new_comb in new_elements:
                    comb_list.append(new_comb)

            self._frontier.append((x, y))
            addAtTheEnd(self._comb_list)
        
        def removeTile(self, x, y):
            def remove(target, x, y):
                idx = 0
                while target[idx] != (x, y):
                    idx += 1
                target.pop(idx)
                return idx
        
            index = remove(self._frontier, x, y)
            for i in range(len(self._comb_list)):
                comb = self._comb_list[i]
                self._comb_list[i] = comb[:index] + comb[index + 1:]


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
        self.__allCombs = self.AllCombs()
        self.__uncheckedLocation = []

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
            self.__uncheckedLocation.append((x, y))
            self.__board[x][y][1] = label - mineCount # decrement the effective label by the count of flagged mines in the neighborhood
        else:
            self.__numOfMinesLeft -= 1

    def getFrontierU(self) -> "list[tuple]":
        # This function gives a list of all the current "Uncovered" frontiers based on the current board condition
        # Return a list of frontiers, each frontier is a list of tile coordinates. e.g. frontier = [(1, 3), (0, 4), (0, 5), (1, 6)]
              
        frontier = []
        for i in range(self.__rowDimension):
            for j in range(self.__colDimension):
                if type(self.__board[i][j][0]) == int and self.__board[i][j][2] != 0:
                    frontier.append((i, j))

        return frontier

    def getFrontierC(self, FU) -> "list[tuple]":
        # This function computes an uncovered frontier based on the current board condition and a given covered frontier
        # A given covered frontier is needed as multiple frontiers can exist at the same time
        frontier = []
        for f in FU:
            for coord in self.getValidNeighbors(f[0], f[1]):
                if self.__board[coord[0]][coord[1]][0] == "*" and coord not in frontier:
                    frontier.append(coord)
        return frontier
    
    def splitFrontier(self, FC, FU):
        # need to complete
        return 1
               

    def analyzeFrontier(self, FU, FC) -> "list[float]":
        # need to complete
        return 1

    def handleReport(self, report, FC) -> None:
        # need to complete
        return

    def printBoardInfo(self) -> None:
        print(self.__numOfMinesLeft, "mines left")
        r = len(self.__board)
        c = len(self.__board[0])
        for row in self.__board:
            for col in row:
                print(col[0], end=" |")
            print()
            print("-"*3*self.__colDimension)

    def getAction(self, number: int) -> "Action Object":
        
        def analyzeBoard():
            frontierU = self.getFrontierU()    
            frontierC = self.getFrontierC(frontierU)
        
        # need to complete
        
        def uncoverAll():
            for i in range(self.__rowDimension):
                for j in range(self.__colDimension):
                    if self.__board[i][j][0] == "*" and (i, j) not in self.__toUncover:
                        self.__toUncover.append((i, j))

        self.updateBoard(self.__X, self.__Y, number) # After getting the face number back, we update the board (knowledge base) with this new knowledge
        while len(self.__toUncover) == 0:
            #self.printBoardInfo()
            #print()
            if self.__numOfMinesLeft == 0:
                timeUsed = time.time() - self.__time
                if timeUsed > 10:
                    print("%.2f seconds used" % timeUsed)
                return Action(AI.Action.LEAVE)
            
            action = analyzeBoard()
            if action == "leave":
                return Action(AI.Action.LEAVE)
            
            #print("uncover", self.__toUncover)
            if self.__numOfMinesLeft == 0:
                uncoverAll()
                if len(self.__toUncover) == 0:
                    timeUsed = time.time() - self.__time
                    if timeUsed > 10:
                        print("%.2f seconds used" % timeUsed)
                    return Action(AI.Action.LEAVE)
        nextTile = self.__toUncover.pop(0)
        #self.printBoardInfo()
        #print()
        #print("next move:", nextTile)
        self.__X, self.__Y = nextTile[0], nextTile[1]
        #print("It has taken", time.time() - self.__time, "seconds until now.")
        return Action(AI.Action.UNCOVER, nextTile[1], nextTile[0])


        


        
        
            