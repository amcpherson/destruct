/*
 *  DiskPriorityQueue.h
 *
 *  Created by Andrew McPherson on 06/12/13.
 *
 */

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <boost/shared_ptr.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>


using namespace boost;
using namespace std;


template <typename TType, class Compare>
struct FileStack
{
	FileStack(const string& filename) : mLen(0), mIdx(0)
	{
		mFile = boost::shared_ptr<ifstream>(new ifstream(filename.c_str(), ios::in|ios::binary));
		if (mFile->good())
		{
			// Stream via a gzip decompressor
			mStream = boost::shared_ptr<iostreams::filtering_streambuf<iostreams::input> >(new iostreams::filtering_streambuf<iostreams::input>());
			mStream->push(iostreams::gzip_decompressor());
			mStream->push(*mFile);

			// Input archive in binary from the stream
			mArchive = boost::shared_ptr<archive::binary_iarchive>(new archive::binary_iarchive(*mStream));

			// Number of elements at the top of the archive
			(*mArchive) >> mLen;
		}

		if (mFile->good())
		{
			(*mArchive) >> Top;

			mIdx++;
		}
	}

	void Pop()
	{
		if (Good())
		{
			(*mArchive) >> Top;

			mIdx++;
		}
	}

	bool Good() const
	{
		return mFile->good() && mIdx < mLen;
	}

	bool operator<(const FileStack<TType,Compare>& other) const
	{
		return !Compare()(Top, other.Top);
	}

	boost::shared_ptr<ifstream> mFile;
	boost::shared_ptr<iostreams::filtering_streambuf<iostreams::input> > mStream;
	boost::shared_ptr<archive::binary_iarchive> mArchive;
	size_t mLen;
	size_t mIdx;

	TType Top;
};

template <typename TType, class Compare=less<TType> >
struct DiskPriorityQueue
{
	DiskPriorityQueue(const string& prefix, int maxSize) : mPrefix(prefix), mMaxSize(maxSize), mSuffix(0)
	{
	}

	~DiskPriorityQueue()
	{
		// Delete alll temporary files
		for (vector<string>::const_iterator iter = mFilenames.begin(); iter != mFilenames.end(); iter++)
		{
			if(remove(iter->c_str()) != 0)
			{
				string errorMsg = "Error deleting file " + *iter;
				perror(errorMsg.c_str());
			}
		}
	}

	void Push(const TType& a)
	{
		mBuffer.push_back(a);

		if (mBuffer.size() >= mMaxSize)
		{
			Flush();
		}
	}

	void Finalize()
	{
		Flush();

		for (vector<string>::const_iterator iter = mFilenames.begin(); iter != mFilenames.end(); iter++)
		{
			mQueue.push(FileStack<TType,Compare>(*iter));
		}
	}
	
	void Flush()
	{
		sort(mBuffer.begin(), mBuffer.end(), Compare());

		// Create filename with numeric suffix
		ostringstream suffixConvert;
		suffixConvert << mSuffix;
		string outFilename = mPrefix + suffixConvert.str();
		mSuffix++;

		// Maintain a list of temp filenames, deleted by destructor
		mFilenames.push_back(outFilename);

		ofstream outFile(outFilename.c_str(), ios::out|ios::binary|ios::trunc);

		// Stream via a gzip compressor
		iostreams::filtering_streambuf<iostreams::output> outStream;
		outStream.push(iostreams::gzip_compressor());
		outStream.push(outFile);

		// Archive in binary to the stream
		archive::binary_oarchive outArchive(outStream);
		
		// Number of elements at the top of the archive
		size_t len = mBuffer.size();
		outArchive << len;

		// Archive sorted list
		for (typename vector<TType>::const_iterator iter = mBuffer.begin(); iter != mBuffer.end(); iter++)
		{
			outArchive << *iter;
		}

		mBuffer.clear();
	}

	const TType& Top()
	{
		return mQueue.top().Top;
	}

	void Pop()
	{
		if (mQueue.top().Good())
		{
			// Pop from the priority queue but store as temp
			FileStack<TType,Compare> tempStack = mQueue.top();
			mQueue.pop();

			// Pop from the file
			tempStack.Pop();

			// Add back to the priority queue
			mQueue.push(tempStack);
		}
		else
		{
			// End of file, just pop from the priority queue
			mQueue.pop();
		}
	}

	bool Empty()
	{
		return mQueue.empty();
	}

	const string mPrefix;
	int mMaxSize;
	int mSuffix;
	vector<TType> mBuffer;
	vector<string> mFilenames;
	priority_queue<FileStack<TType,Compare> > mQueue;
};


// int main(int argc, char* argv[])
// {
// 	DiskPriorityQueue<int> dpq("dpq", 3);

// 	for (int i = 100; i >= 0; i--)
// 	{
// 		dpq.Push(i);
// 	}

// 	dpq.Finalize();
	
// 	cout << "pop" << endl;

// 	while (!dpq.Empty())
// 	{
// 		cout << dpq.Top() << endl; dpq.Pop();
// 	}

// 	return 0;
// }
	
