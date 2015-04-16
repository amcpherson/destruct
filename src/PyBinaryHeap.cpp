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

#include "BinaryMinHeap.h"

using namespace boost;
using namespace std;


BOOST_PYTHON_MODULE(PyBinaryHeap)
{
	using namespace python;
	
	class_<BinaryMinHeap>("BinaryMinHeap")
		.def("Push", &BinaryMinHeap::Push)
		.def("Pop", &BinaryMinHeap::Pop)
		.def("Min", &BinaryMinHeap::Min)
		.def("Remove", &BinaryMinHeap::Remove)
		.def("ReplaceKey", &BinaryMinHeap::ReplaceKey)
		.def("DecreaseKey", &BinaryMinHeap::DecreaseKey)
		.def("IncreaseKey", &BinaryMinHeap::IncreaseKey)
		.def("Empty", &BinaryMinHeap::Empty)
		.def("Print", &BinaryMinHeap::Print)
		;
}

