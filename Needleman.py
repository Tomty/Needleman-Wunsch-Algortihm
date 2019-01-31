from __future__ import print_function
import re

class Sequence:

	def __init__(self, name, vec):
		self.name = name
		self.vec = vec

class NeedlemanWunsh:

	def __init__(self, seqFile, subMatrixFile, gapPenalty):
		self.subMatrix = []
		self.subMatrixLab = []
		self.scoreMatrix = []
		self.pathMatrix = []
		self.seq1 = Sequence("", [])
		self.seq2 = Sequence("", [])
		self.gapPenalty = gapPenalty
		self.endGapPenalty = 0
		self.align1 = []
		self.align2 = []
		self.score = 0

		self.readFileSeq(seqFile)
		self.readFileSubMatrix(subMatrixFile)

		self.compute()
		#self.print_scoreMatrix(self.scoreMatrix, self.seq1.vec, self.seq2.vec)
		self.traceback()
		self.print_align()

	def compute(self):

		#Init score matrix and path matrix
		self.scoreMatrix = [[0 for x in range(len(self.seq2.vec) + 1)] for y in range(len(self.seq1.vec) + 1)]
		self.pathMat = self.init_pathMat()

		#Compute de score matrix
		for i in range(1,len(self.seq1.vec)+1):
			lVal = self.seq1.vec[i-1]
			for j in range(1,len(self.seq2.vec)+1):
				cVal = self.seq2.vec[j - 1]
				tempScore = []

				#Gap up
				tempScore.append(self.scoreMatrix[i-1][j] + self.gapPenalty)
				#Gap left / no end gap penalty
				if(j==len(self.seq2.vec)):
					tempScore.append(self.scoreMatrix[i][j - 1])
				else:
					tempScore.append(self.scoreMatrix[i][j - 1] + self.gapPenalty)
				#Match/Mismatch
				tempScore.append((self.scoreMatrix[i-1][j-1] + self.getSubScore(lVal,cVal)))
				#Set max Score
				self.scoreMatrix[i][j] = max(tempScore)

				x = tempScore.index(max(tempScore))

				#Traceback Matrix
				if (x == 0):
					self.pathMat[i][j] = 'h'
				elif (x == 1):
					self.pathMat[i][j] = 'g'
				elif (x==2):
					self.pathMat[i][j] = 'd'

		self.score = self.scoreMatrix[i][j]

	def traceback(self):

		l = len(self.seq1.vec)
		c = len(self.seq2.vec)
		elem = self.pathMat[l][c]

		while (l > 0 and c > 0):
			if (elem == 'g'):
				c -= 1
				elem = self.pathMat[l][c]
				self.align1.insert(0, "-")
				self.align2.insert(0, self.seq2.vec[c])
			elif (elem == 'h'):
				l -= 1
				elem = self.pathMat[l][c]
				self.align2.insert(0, "-")
				self.align1.insert(0, self.seq1.vec[l])
			elif (elem == 'd'):
				c -= 1
				l -= 1
				elem = self.pathMat[l][c]
				self.align1.insert(0, self.seq1.vec[l])
				self.align2.insert(0, self.seq2.vec[c])

		while (c > 0):
			c -= 1
			self.align1.insert(0, "-")
			self.align2.insert(0, self.seq2.vec[c])

		while (l > 0):
			l -= 1
			self.align2.insert(0, "-")
			self.align1.insert(0, self.seq1.vec[l])

	def print_align(self):
		print("\n")
		for e in self.align1:
			print(e, end=" ")
		print("")
		for e in self.align2:
			print(e, end=" ")
		print("\n")
		print("Score: "+str(self.score))

	def init_pathMat(self):
		pathMat = [[0 for x in range(len(self.seq2.vec) + 1)] for y in range(len(self.seq1.vec) + 1)]

		for i in range(1, len(self.seq2.vec) + 1):
			pathMat[0][i] = 'g'
		for i in range(1, len(self.seq1.vec) + 1):
			pathMat[i][0] = 'h'

		pathMat[0][0] = 's'

		return pathMat

	def print_scoreMatrix(self,mat,seq1,seq2):

		z = -1
		print("  \t\t", end="")
		for elem in seq2:
			x = str(elem) + "  \t"
			print(x, end="")
		for line in mat:
			print("\n")
			if (z >= 0):
				print(seq1[z], end="\t")
			else:
				print("\t", end="")
			z += 1
			for elem in line:
				x = str(elem) + "\t"
				print(x, end="")

	def getSubScore(self, s1, s2):
		return self.subMatrix[self.subMatrixLab.index(s1)][self.subMatrixLab.index(s2)]

	def readFileSeq(self, file):
		vec1 = []
		vec2 = []
		with open(file, 'r') as f:
			name = f.readline().split('>')
			line = f.readline()
			for ch in line:
				vec1.append(ch)
				self.seq1 = Sequence(name[1],vec1)

			self.removeThis(self.seq1.vec, '\n')

			name = f.readline().split('>')
			line = f.readline()
			for ch in line:
				vec2.append(ch)
				self.seq2 = Sequence(name[1], vec2)

			self.removeThis(self.seq2.vec, '\n')

	def removeThis(self, list, this):
		for item in range(list.count(this)):
			list.remove(this)

	def readFileSubMatrix(self, file):
		with open(file, 'r') as f:
			vec1 = re.split(r'[\t\n\s]',f.readline())
			self.removeThis(vec1, '')
			self.subMatrixLab = vec1
			lines = f.readlines()

			for i in range(0,len(lines)):
				line = []
				vec = re.split(r'[^-?\d*\.{0,1}\d+$]',lines[i])
				self.removeThis(vec, '')

				for item in vec:
					line.append(int(item))

				self.subMatrix.append(line)



NeedlemanWunsh("seq.txt", "pam250.tab", -6)

