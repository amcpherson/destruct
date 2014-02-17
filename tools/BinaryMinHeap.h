/*
 *  BinaryMinHeap.h
 */

#ifndef BINARYMINHEAP_H_
#define BINARYMINHEAP_H_

#include "Common.h"
#include "DebugCheck.h"

#include <iostream>
#include <string>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


class BinaryMinHeap
{
public:
	void Push(int id, double priority)
	{
		if (mPositions.find(id) != mPositions.end())
		{
			ReplaceKey(id, priority);
			return;
		}
		
		int index = mPriorities.size();
		
		mPriorities.push_back(priority);
		mPositions[id] = index;
		mIDs.push_back(id);
		
		BubbleUp(index);
	}
	
	int Pop()
	{
		swap(mPriorities.front(), mPriorities.back());
		swap(mPositions[mIDs.front()], mPositions[mIDs.back()]);
		swap(mIDs.front(), mIDs.back());
		
		int minID = mIDs.back();
		mIDs.pop_back();
		mPriorities.pop_back();
		mPositions.erase(minID);
		
		BubbleDown(0);
		
		return minID;
	}
	
	double MinPriority() const
	{
		if (mPriorities.empty())
		{
			cerr << "Error: Attempt to find minimum of empty heap" << endl;
			exit(1);
		}
		
		return mPriorities.front();
	}
	
	int MinID() const
	{
		if (mPriorities.empty())
		{
			cerr << "Error: Attempt to find minimum of empty heap" << endl;
			exit(1);
		}
		
		return mIDs.front();
	}
	
	void Remove(int id)
	{
		if (mPositions.find(id) == mPositions.end())
		{
			cerr << "Error: Attempt to remove non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		swap(mPriorities[index], mPriorities.back());
		swap(mPositions[id], mPositions[mIDs.back()]);
		swap(mIDs[index], mIDs.back());
		
		mIDs.pop_back();
		mPriorities.pop_back();
		mPositions.erase(id);
		
		BubbleUp(index);
		BubbleDown(index);
	}
	
	void ReplaceKey(int id, double priority)
	{
		if (mPositions.find(id) == mPositions.end())
		{
			cerr << "Error: Attempt to replace non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		mPriorities[index] = priority;
		
		BubbleUp(index);
		BubbleDown(index);
	}
	
	bool DecreaseKey(int id, double priority)
	{
		if (mPositions.find(id) == mPositions.end())
		{
			cerr << "Error: Attempt to decrease key for non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		if (priority > mPriorities[index])
		{
			cerr << "Error: Attempt to decrease key to larger value" << endl;
			exit(1);
		}
		
		mPriorities[index] = priority;
		BubbleUp(index);
		
		return true;
	}
	
	bool IncreaseKey(int id, double priority)
	{
		if (mPositions.find(id) == mPositions.end())
		{
			cerr << "Error: Attempt to increase key non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		if (priority < mPriorities[index])
		{
			cerr << "Error: Attempt to increase key to smaller value" << endl;
			exit(1);
		}
		
		mPriorities[index] = priority;
		BubbleDown(index);
		
		return true;
	}
	
	bool Empty() const
	{
		return mPriorities.empty();
	}
	
	int Size() const
	{
		return mPriorities.size();
	}
	
	void Print() const
	{
		for (int index = 0; index < (int)mPositions.size(); index++)
		{
			cerr << mIDs[index] << ":" << mPriorities[index] << "\t";
		}
		cerr << endl;
	}
	
private:
	void BubbleUp(int index)
	{
		while (index > 0)
		{
			int parent = (index-1) / 2;
			
			if (mPriorities[parent] > mPriorities[index])
			{
				swap(mPriorities[index], mPriorities[parent]);
				swap(mPositions[mIDs[index]], mPositions[mIDs[parent]]);
				swap(mIDs[index], mIDs[parent]);
				index = parent;
			}
			else
			{
				break;
			}
		}
	}
	
	void BubbleDown(int index)
	{
		while (index < (int)mPriorities.size())
		{
			int child1 = 2 * index + 1;
			int child2 = 2 * index + 2;
			
			int smallest = index;
			if (child1 < (int)mPriorities.size() && mPriorities[child1] < mPriorities[index])
			{
				smallest = child1;
			}
			if (child2 < (int)mPriorities.size() && mPriorities[child2] < mPriorities[smallest])
			{
				smallest = child2;
			}
			
			if (smallest != index)
			{
				swap(mPriorities[index], mPriorities[smallest]);
				swap(mPositions[mIDs[index]], mPositions[mIDs[smallest]]);
				swap(mIDs[index], mIDs[smallest]);
				index = smallest;
			}
			else
			{
				break;
			}
		}
	}
	
	vector<double> mPriorities;
	unordered_map<int,int> mPositions;
	vector<int> mIDs;
};

class BinaryMinHeapSharedIndex
{
public:
	explicit BinaryMinHeapSharedIndex(vector<int>& sharedIndex) : mPositions(sharedIndex) {}
	
	~BinaryMinHeapSharedIndex()
	{
		for (vector<int>::const_iterator idIter = mIDs.begin(); idIter != mIDs.end(); idIter++)
		{
			mPositions[*idIter] = -1;
		}
	}
	
	void Push(int id, double priority)
	{
		if (mPositions[id] >= 0)
		{
			ReplaceKey(id, priority);
			return;
		}
		
		int index = mPriorities.size();
		
		mPriorities.push_back(priority);
		mPositions[id] = index;
		mIDs.push_back(id);
		
		BubbleUp(index);
	}
	
	int Pop()
	{
		swap(mPriorities.front(), mPriorities.back());
		swap(mPositions[mIDs.front()], mPositions[mIDs.back()]);
		swap(mIDs.front(), mIDs.back());
		
		int minID = mIDs.back();
		mIDs.pop_back();
		mPriorities.pop_back();
		mPositions[minID] = -1;
		
		BubbleDown(0);
		
		return minID;
	}
	
	double MinPriority() const
	{
		if (mPriorities.empty())
		{
			cerr << "Error: Attempt to find minimum of empty heap" << endl;
			exit(1);
		}
		
		return mPriorities.front();
	}
	
	int MinID() const
	{
		if (mPriorities.empty())
		{
			cerr << "Error: Attempt to find minimum of empty heap" << endl;
			exit(1);
		}
		
		return mIDs.front();
	}
	
	void Remove(int id)
	{
		if (mPositions[id] < 0)
		{
			cerr << "Error: Attempt to remove non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		swap(mPriorities[index], mPriorities.back());
		swap(mPositions[id], mPositions[mIDs.back()]);
		swap(mIDs[index], mIDs.back());
		
		mIDs.pop_back();
		mPriorities.pop_back();
		mPositions[id] = -1;
		
		BubbleUp(index);
		BubbleDown(index);
	}
	
	void ReplaceKey(int id, double priority)
	{
		if (mPositions[id] < 0)
		{
			cerr << "Error: Attempt to replace non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		mPriorities[index] = priority;
		
		BubbleUp(index);
		BubbleDown(index);
	}
	
	bool DecreaseKey(int id, double priority)
	{
		if (mPositions[id] < 0)
		{
			cerr << "Error: Attempt to decrease key for non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		if (priority > mPriorities[index])
		{
			cerr << "Error: Attempt to decrease key to larger value" << endl;
			exit(1);
		}
		
		mPriorities[index] = priority;
		BubbleUp(index);
		
		return true;
	}
	
	bool IncreaseKey(int id, double priority)
	{
		if (mPositions[id] < 0)
		{
			cerr << "Error: Attempt to increase key non-existent element" << endl;
			exit(1);
		}
		
		int index = mPositions[id];
		
		if (priority < mPriorities[index])
		{
			cerr << "Error: Attempt to increase key to smaller value" << endl;
			exit(1);
		}
		
		mPriorities[index] = priority;
		BubbleDown(index);
		
		return true;
	}
	
	bool Empty() const
	{
		return mPriorities.empty();
	}
	
	int Size() const
	{
		return mPriorities.size();
	}
	
	void Print() const
	{
		for (int index = 0; index < (int)mPositions.size(); index++)
		{
			cerr << mIDs[index] << ":" << mPriorities[index] << "\t";
		}
		cerr << endl;
	}
	
private:
	void BubbleUp(int index)
	{
		while (index > 0)
		{
			int parent = (index-1) / 2;
			
			if (mPriorities[parent] > mPriorities[index])
			{
				swap(mPriorities[index], mPriorities[parent]);
				swap(mPositions[mIDs[index]], mPositions[mIDs[parent]]);
				swap(mIDs[index], mIDs[parent]);
				index = parent;
			}
			else
			{
				break;
			}
		}
	}
	
	void BubbleDown(int index)
	{
		while (index < (int)mPriorities.size())
		{
			int child1 = 2 * index + 1;
			int child2 = 2 * index + 2;
			
			int smallest = index;
			if (child1 < (int)mPriorities.size() && mPriorities[child1] < mPriorities[index])
			{
				smallest = child1;
			}
			if (child2 < (int)mPriorities.size() && mPriorities[child2] < mPriorities[smallest])
			{
				smallest = child2;
			}
			
			if (smallest != index)
			{
				swap(mPriorities[index], mPriorities[smallest]);
				swap(mPositions[mIDs[index]], mPositions[mIDs[smallest]]);
				swap(mIDs[index], mIDs[smallest]);
				index = smallest;
			}
			else
			{
				break;
			}
		}
	}
	
	vector<double> mPriorities;
	vector<int> mIDs;
	vector<int>& mPositions;
};

class BinaryMaxHeap
{
public:
	void Push(int id, double priority)
	{
		mHeap.Push(id, -priority);
	}
	
	int Pop()
	{
		return mHeap.Pop();
	}
	
	double MaxPriority() const
	{
		return -mHeap.MinPriority();
	}
	
	int MaxID() const
	{
		return mHeap.MinID();
	}
	
	void Remove(int id)
	{
		mHeap.Remove(id);
	}
	
	void ReplaceKey(int id, double priority)
	{
		mHeap.ReplaceKey(id, -priority);
	}
	
	bool DecreaseKey(int id, double priority)
	{
		return mHeap.IncreaseKey(id, -priority);
	}
	
	bool IncreaseKey(int id, double priority)
	{
		return mHeap.DecreaseKey(id, -priority);
	}
	
	bool Empty() const
	{
		return mHeap.Empty();
	}
	
	int Size() const
	{
		return mHeap.Size();
	}
	
	void Print() const
	{
		mHeap.Print();
	}
	
private:
	
	BinaryMinHeap mHeap;
};

class BinaryMaxHeapSharedIndex
{
public:
	explicit BinaryMaxHeapSharedIndex(vector<int>& sharedIndex) : mHeap(sharedIndex) {}
	
	void Push(int id, double priority)
	{
		mHeap.Push(id, -priority);
	}
	
	int Pop()
	{
		return mHeap.Pop();
	}
	
	double MaxPriority() const
	{
		return -mHeap.MinPriority();
	}
	
	int MaxID() const
	{
		return mHeap.MinID();
	}
	
	void Remove(int id)
	{
		mHeap.Remove(id);
	}
	
	void ReplaceKey(int id, double priority)
	{
		mHeap.ReplaceKey(id, -priority);
	}
	
	bool DecreaseKey(int id, double priority)
	{
		return mHeap.IncreaseKey(id, -priority);
	}
	
	bool IncreaseKey(int id, double priority)
	{
		return mHeap.DecreaseKey(id, -priority);
	}
	
	bool Empty() const
	{
		return mHeap.Empty();
	}
	
	int Size() const
	{
		return mHeap.Size();
	}
	
	void Print() const
	{
		mHeap.Print();
	}
	
private:
	
	BinaryMinHeapSharedIndex mHeap;
};

#endif
