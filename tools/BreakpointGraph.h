/*
 *  BreakpointGraph.h
 */

#ifndef BREAKPOINTGRAPH_H_
#define BREAKPOINTGRAPH_H_

#include "Common.h"
#include "DebugCheck.h"

#include <iostream>
#include <string>
#include <boost/unordered_map.hpp>

using namespace boost;
using namespace std;


struct PackedVertexData
{
public:
	enum EVertexType
	{
		eOutgoing = 0,
		eIncoming = 1,
	};
	
	PackedVertexData(unsigned int clusterID, unsigned int clusterEnd, unsigned int vertexType) : mClusterID(clusterID + eFirstClusterID), mClusterEnd(clusterEnd), mVertexType(vertexType)
	{}
	
	PackedVertexData(unsigned int packed) : mPacked(packed)
	{}
	
	inline unsigned int ClusterID() { return mClusterID - eFirstClusterID; }
	inline unsigned int ClusterEnd() { return mClusterEnd; }
	inline unsigned int VertexType() { return mVertexType; }
	inline unsigned int Packed() { return mPacked; }
	
	inline static PackedVertexData StartVertex()
	{
		static PackedVertexData sStartVertex(PackedVertexData::eStartClusterID - eFirstClusterID,0,0);
		return sStartVertex;
	}
	
	inline static PackedVertexData EndVertex()
	{
		static PackedVertexData sEndVertex(PackedVertexData::eEndClusterID - eFirstClusterID,0,0);
		return sEndVertex;
	}
	
	inline bool IsStartVertex() { return mClusterID == eStartClusterID; }
	inline bool IsEndVertex() { return mClusterID == eEndClusterID; }
	inline bool IsSpecialVertex() { return IsStartVertex() || IsEndVertex(); }
	
private:
	enum ESpecialClusterID
	{
		eStartClusterID = 0,
		eEndClusterID = 1,
		eFirstClusterID = 2,
	};
	
	union
	{
		struct
		{
			unsigned mClusterID : 30;
			unsigned mClusterEnd : 1;
			unsigned mVertexType : 1;
		};
		
		unsigned int mPacked;
	};
};

class BreakpointChildIterator : public IChildIterator
{
public:
	BreakpointChildIterator(const vector<pair<int,unsigned int> >& breakpoints, double toScore, int position, int leftIndex, int rightIndex, float leftLambda, float rightLambda)
	: mBreakpoints(breakpoints), mToScore(toScore), mPosition(position), mLeftIndex(leftIndex), mRightIndex(rightIndex), mLeftLambda(leftLambda), mRightLambda(rightLambda)
	{
		Next();
	}
	
	virtual bool Valid()
	{
		return mCurrentIndex >= 0;
	}
	
	virtual void Next()
	{
		mCurrentIndex = NextIndex();
		
		if (mCurrentIndex >= 0 && mCurrentIndex < (int)mBreakpoints.size())
		{
			mChildVertex = mBreakpoints[mCurrentIndex].second;
			mChildScore = ToBreakpointScore(mBreakpoints[mCurrentIndex].first);
		}
	}
	
	virtual unsigned int ChildID()
	{
		return mChildVertex;
	}
	
	virtual double ToChildScore()
	{
		return mChildScore + mToScore;
	}
	
protected:
	BreakpointChildIterator(const vector<pair<int,unsigned int> >& breakpoints, double toScore, int position, int leftIndex, int rightIndex, float leftLambda, float rightLambda, int noNext)
	: mBreakpoints(breakpoints), mToScore(toScore), mPosition(position), mLeftIndex(leftIndex), mRightIndex(rightIndex), mLeftLambda(leftLambda), mRightLambda(rightLambda)
	{
	}
	
	double ToBreakpointScore(int position)
	{
		if (position < mPosition)
		{
			return ((double)mPosition - (double)position) / (double)mLeftLambda;
		}
		else
		{
			return ((double)position - (double)mPosition) / (double)mRightLambda;
		}
	}
	
	int NextIndex()
	{
		int index;
		
		if (mLeftIndex < 0 && mRightIndex >= (int)mBreakpoints.size())
		{
			index = -1;
		}
		else if (mLeftIndex < 0)
		{
			index = mRightIndex;
			mRightIndex++;
		}
		else if (mRightIndex >= (int)mBreakpoints.size())
		{
			index = mLeftIndex;
			mLeftIndex--;
		}
		else if (ToBreakpointScore(mBreakpoints[mLeftIndex].first) < ToBreakpointScore(mBreakpoints[mRightIndex].first))
		{
			index = mLeftIndex;
			mLeftIndex--;
		}
		else
		{
			index = mRightIndex;
			mRightIndex++;
		}
		
		return index;
	}
	
	const vector<pair<int,unsigned int> >& mBreakpoints;
	
	float mToScore;
	int mPosition;
	
	int mLeftIndex;
	int mRightIndex;
	
	float mLeftLambda;
	float mRightLambda;
	
	int mCurrentIndex;
	unsigned int mChildVertex;
	float mChildScore;
};



class PathEndBreakpointChildIterator : public BreakpointChildIterator
{
public:
	PathEndBreakpointChildIterator(const vector<pair<int,unsigned int> >& breakpoints, double toScore, int position, int leftIndex, int rightIndex, float leftLambda, float rightLambda, int endPosition)
	: BreakpointChildIterator(breakpoints,toScore,position,leftIndex,rightIndex,leftLambda,rightLambda,0), mEndPosition(endPosition)
	{
		Next();
	}
	
	virtual void Next()
	{
		mCurrentIndex = NextIndex();
		
		if (mCurrentIndex < 0 || mCurrentIndex >= (int)mBreakpoints.size() || ToBreakpointScore(mEndPosition) < ToBreakpointScore(mBreakpoints[mCurrentIndex].first))
		{
			mChildVertex = PackedVertexData::EndVertex().Packed();
			mChildScore = ToBreakpointScore(mEndPosition);
			
			mCurrentIndex = 0;
			mLeftIndex = -1;
			mRightIndex = (int)mBreakpoints.size();
		}
		else if (mCurrentIndex >= 0 && mCurrentIndex < (int)mBreakpoints.size())
		{
			mChildVertex = mBreakpoints[mCurrentIndex].second;
			mChildScore = ToBreakpointScore(mBreakpoints[mCurrentIndex].first);
		}
	}
	
private:
	int mEndPosition;
};

class MatePairChildIterator : public IChildIterator
{
public:
	MatePairChildIterator(unsigned int mateVertex, double toMateScore) : mValid(true), mMateVertex(mateVertex), mToMateScore(toMateScore)
	{}
	
	virtual bool Valid()
	{
		return mValid;
	}
	
	virtual void Next()
	{
		mValid = false;
	}
	
	virtual unsigned int ChildID()
	{
		return mMateVertex;
	}
	
	virtual double ToChildScore()
	{
		return mToMateScore;
	}
	
private:
	bool mValid;
	unsigned int mMateVertex;
	double mToMateScore;
};

class NullChildIterator : public IChildIterator
{
public:
	NullChildIterator()
	{}
	
	virtual bool Valid()
	{
		return false;
	}
	
	virtual void Next()
	{
	}
	
	virtual unsigned int ChildID()
	{
		return -1;
	}
	
	virtual double ToChildScore()
	{
		return 0.0;
	}
};

class BreakpointGraph : public IBreakpointGraph
{
public:
	BreakpointGraph(double insertionLambda, double deletionLambda) : mInsertionLambda(insertionLambda), mDeletionLambda(deletionLambda), mConstructed(false), mCycleSearch(false), mPathSearch(false) {}
	
	IChildIterator* BeginChildIteration(unsigned int vertex, double toScore) const
	{
		DebugCheck(mConstructed);
		
		PackedVertexData vertexData(vertex);
		
		DebugCheck(vertexData.IsSpecialVertex() || mBreakpointScores.find(vertexData.ClusterID()) != mBreakpointScores.end());
		
		IChildIterator* childIter = 0;
		if (mBlockedVertices.find(vertex) != mBlockedVertices.end())
		{
			childIter = new NullChildIterator();
		}
		else if (vertexData.VertexType() == PackedVertexData::eIncoming)
		{
			PackedVertexData mateVertexData(vertexData.ClusterID(), OtherClusterEnd(vertexData.ClusterEnd()), PackedVertexData::eOutgoing);
			
			childIter = new MatePairChildIterator(mateVertexData.Packed(), toScore + mBreakpointScores.find(mateVertexData.ClusterID())->second);
		}
		else if (vertexData.VertexType() == PackedVertexData::eOutgoing)
		{
			if (mPathSearch && vertexData.IsStartVertex())
			{
				const vector<pair<int,unsigned int> >& breakpoints = mBreakpoints.find(mPathStartRefStrand.id)->second;
				float leftLambda = (mPathStartRefStrand.strand == PlusStrand) ? mDeletionLambda : mInsertionLambda;
				float rightLambda = (mPathStartRefStrand.strand == PlusStrand) ? mInsertionLambda : mDeletionLambda;
				
				childIter = new BreakpointChildIterator(breakpoints, toScore, mPathStartPosition, mPathStartLeftIndex, mPathStartLeftIndex + 1, leftLambda, rightLambda);				
			}
			else
			{
				RefStrand adjacentRefStrand = mAdjacentRefStrands.find(vertexData.Packed())->second;
				int position = mPositions.find(vertexData.Packed())->second;
				int leftIndex = mLeftIndices.find(vertexData.Packed())->second;
				float leftLambda = (adjacentRefStrand.strand == PlusStrand) ? mDeletionLambda : mInsertionLambda;
				float rightLambda = (adjacentRefStrand.strand == PlusStrand) ? mInsertionLambda : mDeletionLambda;
				
				const vector<pair<int,unsigned int> >& breakpoints = mBreakpoints.find(adjacentRefStrand.id)->second;
				
				if (mPathSearch && adjacentRefStrand.id == mPathEndRefStrand.id)
				{
					childIter = new PathEndBreakpointChildIterator(breakpoints, toScore, position, leftIndex, leftIndex + 1, leftLambda, rightLambda, mPathEndPosition);				
				}
				else
				{
					childIter = new BreakpointChildIterator(breakpoints, toScore, position, leftIndex, leftIndex + 1, leftLambda, rightLambda);				
				}
			}			
		}
		
		DebugCheck(mChildIterators.find(vertex) == mChildIterators.end());
		mChildIterators[vertex] = childIter;
		
		return childIter;			
	}
	
	void EndChildIteration(unsigned int vertex) const
	{
		DebugCheck(mConstructed);
		DebugCheck(mChildIterators.find(vertex) != mChildIterators.end());
		
		delete mChildIterators[vertex];
		mChildIterators.erase(vertex);
	}
	
	void AddBreakpoint(unsigned int clusterID, int refID1, int strand1, int position1, int refID2, int strand2, int position2, double score)
	{
		mConstructed = false;
		
		PackedVertexData vertex1(clusterID, 0, PackedVertexData::eOutgoing);
		PackedVertexData vertex2(clusterID, 1, PackedVertexData::eOutgoing);
		
		unsigned int vertices[2];
		int refIDs[2];
		int strands[2];
		int positions[2];
		
		vertices[0] = vertex1.Packed();
		refIDs[0] = refID1;
		strands[0] = strand1;
		positions[0] = position1;
		
		vertices[1] = vertex2.Packed();
		refIDs[1] = refID2;
		strands[1] = strand2;
		positions[1] = position2;
		
		for (int clusterEnd = 0; clusterEnd <= 1; clusterEnd++)
		{
			RefStrand refStrand;
			refStrand.referenceIndex = refIDs[clusterEnd];
			refStrand.strand = strands[clusterEnd];
			
			RefStrand refOtherStrand;
			refOtherStrand.referenceIndex = refIDs[clusterEnd];
			refOtherStrand.strand = OtherStrand(strands[clusterEnd]);
			
			mRefIDs.insert(refIDs[clusterEnd]);
			mAdjacentRefStrands[vertices[clusterEnd]] = refOtherStrand;
			mBreakpoints[refStrand.id].push_back(make_pair(positions[clusterEnd], vertices[clusterEnd]));
			mPositions[vertices[clusterEnd]] = positions[clusterEnd];
		}
		
		mBreakpointScores[clusterID] = score;
	}
	
	static void PrintVertex(string pref, int vertex)
	{
		PackedVertexData vertexData(vertex);
		cerr << pref << vertexData.ClusterID() << "\t" << vertexData.ClusterEnd() << "\t" << vertexData.VertexType() << "\t" << vertexData.Packed() << endl;
	}
	
	void ConstructGraph()
	{
		for (unordered_map<int,vector<pair<int,unsigned int> > >::iterator breakpointsIter = mBreakpoints.begin(); breakpointsIter != mBreakpoints.end(); breakpointsIter++)
		{
			sort(breakpointsIter->second.begin(), breakpointsIter->second.end());
		}
		
		for (unordered_set<int>::const_iterator refIDIter = mRefIDs.begin(); refIDIter != mRefIDs.end(); refIDIter++)
		{
			RefStrand refPlus;
			RefStrand refMinus;
			
			refPlus.referenceIndex = *refIDIter;
			refPlus.strand = PlusStrand;
			
			refMinus.referenceIndex = *refIDIter;
			refMinus.strand = MinusStrand;
			
			const vector<pair<int,unsigned int> >& breapointsPlus = mBreakpoints[refPlus.id];
			const vector<pair<int,unsigned int> >& breapointsMinus = mBreakpoints[refMinus.id];
			
			int indexPlus = 0;
			int indexMinus = 0;
			int leftIndexPlus = -1;
			int leftIndexMinus = -1;
			
			while (indexPlus < (int)breapointsPlus.size() || indexMinus < (int)breapointsMinus.size())
			{
				if (indexPlus < (int)breapointsPlus.size() && (indexMinus >= (int)breapointsMinus.size() || breapointsPlus[indexPlus].first <= breapointsMinus[indexMinus].first))
				{
					leftIndexPlus = indexPlus;
					mLeftIndices[breapointsPlus[indexPlus].second] = leftIndexMinus;
					indexPlus++;
				}
				else if (indexMinus < (int)breapointsMinus.size() && (indexPlus >= (int)breapointsPlus.size() || breapointsPlus[indexPlus].first > breapointsMinus[indexMinus].first))
				{
					leftIndexMinus = indexMinus;
					mLeftIndices[breapointsMinus[indexMinus].second] = leftIndexPlus;
					indexMinus++;
				}
			}
		}
		
		for (unordered_map<int,vector<pair<int,unsigned int> > >::iterator breakpointsIter = mBreakpoints.begin(); breakpointsIter != mBreakpoints.end(); breakpointsIter++)
		{
			for (vector<pair<int,unsigned int> >::iterator breakpointIter = breakpointsIter->second.begin(); breakpointIter != breakpointsIter->second.end(); breakpointIter++)
			{
				PackedVertexData vertexData(breakpointIter->second);
				vertexData = PackedVertexData(vertexData.ClusterID(), vertexData.ClusterEnd(), PackedVertexData::eIncoming);
				breakpointIter->second = vertexData.Packed();
			}
		}
		
		mConstructed = true;
	}
	
	bool SetCycleSearch(unsigned int clusterID, unsigned int& startVertex, unsigned int& endVertex)
	{
		if (mBreakpointScores.find(clusterID) == mBreakpointScores.end())
		{
			return false;
		}
		
		mCycleSearch = true;
		
		PackedVertexData startVertexData(clusterID, 0, PackedVertexData::eOutgoing);
		PackedVertexData endVertexData(clusterID, 1, PackedVertexData::eIncoming);	
		startVertex = startVertexData.Packed();
		endVertex = endVertexData.Packed();
		
		BlockVertex(clusterID, 0);
		
		return true;
	}
	
	bool SetPathSearch(int refID1, int strand1, int position1, int refID2, int strand2, int position2, unsigned int& startVertex, unsigned int& endVertex)
	{
		DebugCheck(mConstructed);
		
		mPathStartRefStrand.referenceIndex = refID1;
		mPathStartRefStrand.strand = strand1;
		mPathStartPosition = position1;
		if (!PathFindLeftIndex(refID1, strand1, position1, mPathStartLeftIndex))
		{
			return false;
		}
		
		mPathEndRefStrand.referenceIndex = refID2;
		mPathEndRefStrand.strand = OtherStrand(strand2);
		mPathEndPosition = position2;
		
		startVertex = PackedVertexData::StartVertex().Packed();
		endVertex = PackedVertexData::EndVertex().Packed();
		
		mPathSearch = true;
		
		return true;
	}
	
	void BlockVertex(unsigned int clusterID, unsigned int clusterEnd)
	{
		PackedVertexData blockedVertexData(clusterID, clusterEnd, PackedVertexData::eIncoming);
		mBlockedVertices.insert(blockedVertexData.Packed());
	}
	
	void Reset()
	{
		for (unordered_map<unsigned int,IChildIterator*>::iterator iter = mChildIterators.begin(); iter != mChildIterators.end(); iter++)
		{
			delete iter->second;
		}
		mChildIterators.clear();
		
		mBlockedVertices.clear();
		
		mCycleSearch = false;
		mPathSearch = false;
	}
	
	unsigned int GetClusterID(unsigned int vertex) const
	{
		PackedVertexData vertexData(vertex);
		return vertexData.ClusterID();
	}
	
	int GetClusterEnd(unsigned int vertex) const
	{
		PackedVertexData vertexData(vertex);
		return vertexData.ClusterEnd();
	}
	
private:
	bool PathFindLeftIndex(int refID, int strand, int position, int& leftIndex)
	{
		RefStrand refStrand;
		refStrand.referenceIndex = refID;
		refStrand.strand = strand;
		
		const vector<pair<int,unsigned int> >& breakpoints = mBreakpoints.find(refStrand.id)->second;
		
		if (breakpoints.empty())
		{
			return false;
		}
		
		vector<pair<int,unsigned int> >::const_iterator breakpointIter = upper_bound(breakpoints.begin(), breakpoints.end(), pair<int,unsigned int>(position,numeric_limits<unsigned int>::max()));
		
		if (breakpointIter == breakpoints.end())
		{
			leftIndex = breakpoints.size() - 1;
		}
		else
		{
			leftIndex = breakpointIter - breakpoints.begin() - 1;
		}
		
		DebugCheck(leftIndex < 0 || leftIndex >= breakpoints.size() || breakpoints[leftIndex].first <= position);
		DebugCheck(leftIndex + 1 < 0 || leftIndex + 1 >= breakpoints.size() || breakpoints[leftIndex + 1].first >= position);
		
		return true;
	}
	
	const float mInsertionLambda;
	const float mDeletionLambda;
	bool mConstructed;
	
	bool mCycleSearch;
	bool mPathSearch;
	RefStrand mPathStartRefStrand;
	int mPathStartPosition;
	int mPathStartLeftIndex;
	RefStrand mPathEndRefStrand;
	int mPathEndPosition;
		
	unordered_set<int> mRefIDs;
	unordered_map<unsigned int,RefStrand> mAdjacentRefStrands;
	unordered_map<int,vector<pair<int,unsigned int> > > mBreakpoints;
	unordered_map<unsigned int,double> mBreakpointScores;
	unordered_map<unsigned int,int> mPositions;
	
	unordered_map<unsigned int,int> mLeftIndices;
	
	mutable unordered_map<unsigned int,IChildIterator*> mChildIterators;
	
	unordered_set<unsigned int> mBlockedVertices;
};


#endif
