// Copyright Ralf W. Grosse-Kunstleve 2002-2004. Distributed under the Boost
// Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/manage_new_object.hpp>
#include <boost/python/tuple.hpp>

#include "ShortestPath.h"
#include "BreakpointGraph.h"

using namespace boost;
using namespace std;


class BreakpointGraphWrapper : public BreakpointGraph
{
public:
	BreakpointGraphWrapper(double insertionLambda, double deletionLambda) : BreakpointGraph(insertionLambda, deletionLambda) {}
	
	python::object SetCycleSearchWrapper(unsigned int clusterID)
	{
		unsigned int startVertex;
		unsigned int endVertex;
		
		if (!BreakpointGraph::SetCycleSearch(clusterID, startVertex, endVertex))
		{
			return python::object();
		}
		
		return python::make_tuple(startVertex,endVertex);
	}
	
	python::object SetPathSearchWrapper(int refID1, int strand1, int position1, int refID2, int strand2, int position2)
	{
		unsigned int startVertex;
		unsigned int endVertex;
		
		if (!BreakpointGraph::SetPathSearch(refID1, strand1, position1, refID2, strand2, position2, startVertex, endVertex))
		{
			return python::object();
		}
		
		return python::make_tuple(startVertex,endVertex);
	}
};

python::object ShortestPathWrapper(const BreakpointGraphWrapper* breakpointGraph, unsigned int startVertex, unsigned int endVertex, int visitLimit, python::list pyPath)
{
	vector<unsigned int> path;
	double pathScore = 0.0;
	if (!ShortestPath(breakpointGraph, startVertex, endVertex, visitLimit, path, pathScore))
	{
		return python::object();
	}
	
	for (vector<unsigned int>::const_iterator pathIter = path.begin(); pathIter != path.end(); pathIter++)
	{
		pyPath.append(*pathIter);
	}
	
	return python::object(pathScore);
}

BOOST_PYTHON_MODULE(PyBreakpointGraph)
{
	using namespace python;
	
	def("ShortestPath", ShortestPathWrapper);
	
	class_<BreakpointGraphWrapper>("BreakpointGraph", init<double,double>())
		.def("AddBreakpoint", &BreakpointGraphWrapper::AddBreakpoint)
		.def("ConstructGraph", &BreakpointGraphWrapper::ConstructGraph)
		.def("SetCycleSearch", &BreakpointGraphWrapper::SetCycleSearchWrapper)
		.def("SetPathSearch", &BreakpointGraphWrapper::SetPathSearchWrapper)
		.def("GetClusterID", &BreakpointGraphWrapper::GetClusterID)
		;
}

