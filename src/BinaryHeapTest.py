import PyBinaryHeap
import unittest
import random
import math
import operator

class PyBinaryMinHeap:
	def __init__(self):
		self.data = dict()
	def Push(self,id,priority):
		self.data[id] = priority
	def Pop(self):
		top = sorted(self.data.iteritems(), key=operator.itemgetter(1))[0]
		del self.data[top[0]]
		return top[0]
	def Min(self):
		top = sorted(self.data.iteritems(), key=operator.itemgetter(1))[0]
		return top[1]
	def Remove(self,id):
		del self.data[id]
	def ReplaceKey(self,id,priority):
		self.data[id] = priority
	def DecreaseKey(self,id,priority):
		self.data[id] = priority
	def IncreaseKey(self,id,priority):
		self.data[id] = priority
	def Empty(self):
		return len(self.data) == 0

class BinaryHeapTest(unittest.TestCase):
	
	def testRandom(self):
		for seed in range(1000):
			
			base = PyBinaryMinHeap()
			test = PyBinaryHeap.BinaryMinHeap()
			
			random.seed(seed)
			
			operations = 0
			contains = set()
			fail = 0
			for i in range(1000):
				self.assertEqual(base.Empty(), test.Empty(), 'bad empty flag')
				if len(contains) > 0:
					self.assertEqual(base.Min(), test.Min(), 'incorrect min')
				operation = "push"
				if len(contains) > 0:
					operation = random.choice(("pop","push","replace","decrease","increase","remove"))
				if operation == "pop":
					a = base.Pop()
					b = test.Pop()
					contains.remove(a)
					operations += 1
					self.assertEqual(a, b, 'incorrect pop')
				elif operation == "replace":
					priority = random.random()
					id = list(contains)[random.randint(0, len(contains)-1)]
					a = base.ReplaceKey(id,priority)
					b = test.ReplaceKey(id,priority)
					operations += 1
				elif operation == "decrease":
					id = list(contains)[random.randint(0, len(contains)-1)]
					priority = random.uniform(0.0, base.data[id])
					a = base.DecreaseKey(id,priority)
					b = test.DecreaseKey(id,priority)
					operations += 1
				elif operation == "increase":
					id = list(contains)[random.randint(0, len(contains)-1)]
					priority = random.uniform(base.data[id], 1.0)
					a = base.IncreaseKey(id,priority)
					b = test.IncreaseKey(id,priority)
					operations += 1
				elif operation == "remove":
					id = list(contains)[random.randint(0, len(contains)-1)]
					base.Remove(id)
					test.Remove(id)
					contains.remove(id)
					operations += 1
				elif operation == "push":
					priority = random.random()
					id = random.randint(0, 100)
					base.Push(id,priority)
					test.Push(id,priority)
					contains.add(id)
					operations += 1
			self.assertEqual(a,a)

if __name__ == '__main__':
	unittest.main()

