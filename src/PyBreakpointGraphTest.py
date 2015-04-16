import PyBreakpointGraph
import unittest
import random
import math


def cycle_search(bg, id, maxvisit=1000):
	(start,end) = bg.SetCycleSearch(id)
	path = []
	pathScore = PyBreakpointGraph.ShortestPath(bg, start, end, maxvisit, path)
	clusterIDs = []
	for vertex in path:
		clusterID = bg.GetClusterID(vertex)
		if len(clusterIDs) == 0 or clusterID != clusterIDs[-1]:
			clusterIDs.append(bg.GetClusterID(vertex))
	return (pathScore,clusterIDs)

def path_search(bg, ref1, str1, pos1, ref2, str2, pos2, maxvisit=1000):
	(start,end) = bg.SetPathSearch(ref1, str1, pos1, ref2, str2, pos2)
	path = []
	pathScore = PyBreakpointGraph.ShortestPath(bg, start, end, maxvisit, path)
	clusterIDs = []
	for vertex in path:
		if vertex != start and vertex != end:
			clusterID = bg.GetClusterID(vertex)
			if len(clusterIDs) == 0 or clusterID != clusterIDs[-1]:
				clusterIDs.append(bg.GetClusterID(vertex))
	return (pathScore,clusterIDs)

class PyBreakpointGraphTest(unittest.TestCase):
	
	def CreateRandomGraph(self, seed, numBr, delLam, insLam):
		bg = PyBreakpointGraph.BreakpointGraph(delLam,insLam)
		
		random.seed(seed)
		for br in range(numBr):
			ref1 = random.randint(0,9)
			ref2 = random.randint(0,9)
			str1 = random.randint(0,1)
			str2 = random.randint(0,1)
			pos1 = random.randint(1,1000000)
			pos2 = random.randint(1,1000000)
			score = math.log(1.0 - random.random())
			bg.AddBreakpoint(br, ref1, str1, pos1, ref2, str2, pos2, score)
		bg.ConstructGraph()
		
		return bg;
		
	def testSimplePath1(self):
		bg = PyBreakpointGraph.BreakpointGraph(10.0,1.0)
		
		bg.AddBreakpoint(1, 1, 0, 200, 2, 1, 100, 1)
		bg.AddBreakpoint(2, 2, 0, 200, 3, 1, 100, 1)
		bg.ConstructGraph()
		
		result = path_search(bg, 1, 0, 100, 3, 1, 200)
		
		self.assertEqual(result, (32.0, [1, 2]), 'incorrect path')
		
	def testSimplePath2(self):
		bg = PyBreakpointGraph.BreakpointGraph(2000.0,200.0)
		
		bg.AddBreakpoint(12511188, 11, 0, 70596341, 11, 0, 72732885, 4.97534e-07)
		bg.ConstructGraph()
		
		result = path_search(bg, 11, 0, 70596485, 11, 0, 72732835)
		
		self.assertEqual(result, (0.74500048160552979, [12511188]), 'incorrect path')
		
	def testSimpleCycle1(self):
		bg = PyBreakpointGraph.BreakpointGraph(10.0,1.0)
		
		bg.AddBreakpoint(0, 1, 0, 200, 2, 1, 100, 1)
		bg.AddBreakpoint(1, 2, 0, 200, 3, 1, 100, 1)
		bg.AddBreakpoint(2, 3, 0, 200, 4, 1, 100, 1)
		bg.AddBreakpoint(3, 4, 0, 200, 1, 1, 100, 1)
		bg.ConstructGraph()
		
		result = cycle_search(bg, 0)
		
		self.assertEqual(result, (43.0, [0, 3, 2, 1, 0]), 'incorrect cycle')
		
	def testSimpleCycle2(self):
		bg = PyBreakpointGraph.BreakpointGraph(10.0,10.0)
		
		bg.AddBreakpoint(0, 6, 0, 17782176, 4, 1, 59212722, 0.0000000001)
		bg.AddBreakpoint(1, 6, 0, 17782176, 4, 1, 59212722, 0.00000000001)
		bg.AddBreakpoint(2, 6, 0, 17782176, 4, 1, 59212722, 0.0000000000000001)
		bg.AddBreakpoint(3, 6, 0, 17782176, 3, 1,  59212722, 1)
		bg.AddBreakpoint(4, 2, 0, 11861249, 9, 0, 134085264, 2)
		bg.AddBreakpoint(5, 2, 1, 11861626, 3, 0,  59212565, 3)
		bg.AddBreakpoint(6, 6, 1, 17782626, 9, 1, 134085420, 2)
		bg.ConstructGraph()
		
		result = cycle_search(bg, 5)
		
		self.assertEqual(result, (119.0, [5, 4, 6, 3, 5]), 'incorrect cycle')
		
	def testRandomCycle(self):
		bg = self.CreateRandomGraph(0, 10000, 10.0,10.0)
		
		result = cycle_search(bg, 0)
		
	def testRandomPath(self):
		bg = self.CreateRandomGraph(0, 10000, 10.0,10.0)
		
		result = path_search(bg, 1, 0, 500000, 2, 0, 500000)
		

		

if __name__ == '__main__':
	unittest.main()

