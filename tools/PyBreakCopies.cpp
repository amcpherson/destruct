//
//  PyBreakCopies.cpp
//

#include <boost/python/class.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/overloads.hpp>

#include <string>
#include <fstream>
#include <boost/unordered_map.hpp>
#include <boost/lexical_cast.hpp>

#include <coin/IpIpoptApplication.hpp>

#include <Common.h>

using namespace std;
using namespace boost;
using namespace Ipopt;

class BreakCopies : public TNLP
{
public:
	// constructor
	BreakCopies()
	{}

	// destructor
	~BreakCopies()
	{}

	struct EdgeInfo
	{
		double read_count;
		double region_length;
		double obj_coeff;
		double upper_bound;
	};

	struct NodeInfo
	{
		int edge_idx;
		double edge_coeff;
	};

	vector<EdgeInfo> mEdges;
	unordered_map<int,vector<NodeInfo> > mNodes;
	int mIncidences;

	vector<double> mCopies;

	// read the problem from files
	void ReadProblem(const string& edge_filename, const string& node_filename)
	{
		{
			ifstream edge_file(edge_filename.c_str());
			
			if (!edge_file.good())
			{
				throw runtime_error("Error: Unable to open file " + edge_filename);
			}

			mEdges.clear();
			
			StringVec fields;
			while (ReadTSV(edge_file, fields))
			{
				EdgeInfo edge_info;
				edge_info.read_count = SAFEPARSE(double, fields[0]);
				edge_info.region_length = SAFEPARSE(double, fields[1]);
				edge_info.obj_coeff = SAFEPARSE(double, fields[2]);
				edge_info.upper_bound = SAFEPARSE(double, fields[3]);

				mEdges.push_back(edge_info);
			}
		}

		{
			ifstream node_file(node_filename.c_str());
			
			if (!node_file.good())
			{
				throw runtime_error("Error: Unable to open file " + node_filename);
			}

			mNodes.clear();
			mIncidences = 0;
			
			StringVec fields;
			while (ReadTSV(node_file, fields))
			{
				int node_idx = SAFEPARSE(int, fields[0]);

				NodeInfo node_info;
				node_info.edge_idx = SAFEPARSE(int, fields[1]);
				node_info.edge_coeff = SAFEPARSE(double, fields[2]);

				mNodes[node_idx].push_back(node_info);
				mIncidences++;
			}
		}
	}

	// initialize the problem from python lists
	void InitProblem(python::list edges, python::list nodes)
	{
		{
			mEdges.clear();
			
			python::ssize_t num_edges = python::len(edges);

			for (python::ssize_t i = 0; i < num_edges; i++)
			{
				python::tuple edge_tuple = python::extract<python::tuple>(edges[i]);

				EdgeInfo edge_info;
				edge_info.read_count = python::extract<double>(edge_tuple[0]);
				edge_info.region_length = python::extract<double>(edge_tuple[1]);
				edge_info.obj_coeff = python::extract<double>(edge_tuple[2]);
				edge_info.upper_bound = python::extract<double>(edge_tuple[3]);

				mEdges.push_back(edge_info);
			}
		}

		{
			mNodes.clear();
			mIncidences = 0;
			
			python::ssize_t num_nodes = python::len(nodes);

			for (python::ssize_t i = 0; i < num_nodes; i++)
			{
				python::tuple node_tuple = python::extract<python::tuple>(nodes[i]);

				int node_idx = python::extract<int>(node_tuple[0]);

				NodeInfo node_info;
				node_info.edge_idx = python::extract<int>(node_tuple[1]);
				node_info.edge_coeff = python::extract<double>(node_tuple[2]);

				mNodes[node_idx].push_back(node_info);
				mIncidences++;
			}
		}
	}

	// returns the size of the problem
	bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
	                  Index& nnz_h_lag, IndexStyleEnum& index_style)
	{
		// One variable per edge 
		n = mEdges.size();

		// One equality constraint per node
		m = mNodes.size();

		// One nonzero entry per node edge incidence
		nnz_jac_g = mIncidences;

		// Only the objective contributes to the hessian of the lagrangian
		nnz_h_lag = mEdges.size();

		// use the C style indexing (0-based)
		index_style = TNLP::C_STYLE;

		return true;
	}

	// returns the variable bounds
	bool get_bounds_info(Index n, Number* x_l, Number* x_u,
	                     Index m, Number* g_l, Number* g_u)
	{
		// Edge variables have lower bound 0, upper bound specified as part of the problem
		for (Index i = 0; i < n; i++)
		{
			x_l[i] = 0.0;
			x_u[i] = mEdges[i].upper_bound;
		}

		// Node constraints equal 0
		for (Index i = 0; i < m; i++)
		{
			g_l[i] = 0.0;
			g_u[i] = 0.0;
		}

		return true;
	}

	// returns the initial point for the problem
	bool get_starting_point(Index n, bool init_x, Number* x,
                            bool init_z, Number* z_L, Number* z_U,
                            Index m, bool init_lambda,
                            Number* lambda)
	{
		// Here, we assume we only have starting values for x, if you code
		// your own NLP, you can provide starting values for the dual variables
		// if you wish
		assert(init_x == true);
		assert(init_z == false);
		assert(init_lambda == false);

		// initialize to the given starting point at 0
		for (Index i = 0; i < n; i++)
		{
			x[i] = 0.0;
		}

		return true;
	}

	// returns the value of the objective function
	bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
	{
		obj_value = 0.0;
		for (Index i = 0; i < n; i++)
		{
			// Add likelihood if we have at least one read
			if (mEdges[i].read_count > 0)
			{
				obj_value += mEdges[i].region_length * x[i] - mEdges[i].read_count * log(x[i] + 1e-20);
			}

			// Add regularization penalty
			obj_value += mEdges[i].obj_coeff * x[i];
		}

		return true;
	}

	// return the gradient of the objective function grad_{x} f(x)
	bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
	{
		for (Index i = 0; i < n; i++)
		{
			grad_f[i] = 0.0;

			// Likelihood gradiant if we have at least one read
			if (mEdges[i].read_count > 0)
			{
				grad_f[i] += mEdges[i].region_length - mEdges[i].read_count / (x[i] + 1e-20);
			}

			// Add regularization penalty
			grad_f[i] += mEdges[i].obj_coeff;
		}

		return true;
	}

	// return the value of the constraints: g(x)
	bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
	{
		for (Index i = 0; i < m; i++)
		{
			g[i] = 0.0;

			for (vector<NodeInfo>::const_iterator iter = mNodes[i].begin(); iter != mNodes[i].end(); iter++)
			{
				g[i] += x[iter->edge_idx] * iter->edge_coeff;
			}
		}

		return true;
	}

	// return the structure or values of the jacobian
	bool eval_jac_g(Index n, const Number* x, bool new_x,
                    Index m, Index nele_jac, Index* iRow, Index *jCol,
                    Number* values)
	{
		if (values == NULL)
		{
			int incidence = 0;

			// return the structure of the jacobian
			for (Index i = 0; i < m; i++)
			{
				for (vector<NodeInfo>::const_iterator iter = mNodes[i].begin(); iter != mNodes[i].end(); iter++)
				{
					iRow[incidence] = i;
					jCol[incidence] = iter->edge_idx;
					incidence++;
				}
			}

			assert(incidence = mIncidences);
		}
		else
		{
			int incidence = 0;

			// return the values of the jacobian of the constraints
			for (Index i = 0; i < m; i++)
			{
				for (vector<NodeInfo>::const_iterator iter = mNodes[i].begin(); iter != mNodes[i].end(); iter++)
				{
					values[incidence] = iter->edge_coeff;
					incidence++;
				}
			}

			assert(incidence = mIncidences);
		}

		return true;
	}

	//return the structure or values of the hessian
	bool eval_h(Index n, const Number* x, bool new_x,
                Number obj_factor, Index m, const Number* lambda,
                bool new_lambda, Index nele_hess, Index* iRow,
                Index* jCol, Number* values)
	{
		if (values == NULL)
		{
			// the hessian is the identity matrix of the objective
			for (Index i = 0; i < n; i++)
			{
				iRow[i] = i;
				jCol[i] = i;
			}
		}
		else
		{
			// return the values of the hessian of the lagrangian
			for (Index i = 0; i < n; i++)
			{
				values[i] = obj_factor * mEdges[i].read_count / ((x[i] + 1e-20) * (x[i] + 1e-20));
			}
		}

		return true;
	}

	void finalize_solution(SolverReturn status,
						   Index n, const Number* x, const Number* z_L, const Number* z_U,
						   Index m, const Number* g, const Number* lambda,
						   Number obj_value,
						   const IpoptData* ip_data,
						   IpoptCalculatedQuantities* ip_cq)
	{
		mCopies.resize(n);
		for (Index i = 0; i < n; i++)
		{
			mCopies[i] = x[i];
		}

		/*
	  // here is where we would store the solution to variables, or write to a file, etc
	  // so we could use the solution.

	  // For this example, we write the solution to the console
	  std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
	  for (Index i=0; i<n; i++) {
	     std::cout << "x[" << i << "] = " << x[i] << std::endl;
	  }

	  std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
	  for (Index i=0; i<n; i++) {
	    std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
	  }
	  for (Index i=0; i<n; i++) {
	    std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
	  }

	  std::cout << std::endl << std::endl << "Objective value" << std::endl;
	  std::cout << "f(x*) = " << obj_value << std::endl;

	  std::cout << std::endl << "Final value of the constraints:" << std::endl;
	  for (Index i=0; i<m ;i++) {
	    std::cout << "g(" << i << ") = " << g[i] << std::endl;
	  }
	  */
	}
};


/*
int main(int argv, char* argc[])
{
	// Create a new instance of your nlp 
	//  (use a SmartPtr, not raw)
	SmartPtr<BreakCopies> mynlp = new BreakCopies();

	mynlp->ReadProblem("edges.txt", "nodes.txt");

	// Create a new instance of IpoptApplication
	//  (use a SmartPtr, not raw)
	// We are using the factory, since this allows us to compile this
	// example with an Ipopt Windows DLL
	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();

	// Change some options
	// Note: The following choices are only examples, they might not be
	//       suitable for your optimization problem.
	app->Options()->SetNumericValue("tol", 1e-9);
	app->Options()->SetIntegerValue("print_level", 0);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetStringValue("output_file", "ipopt.out");
	app->Options()->SetStringValue("derivative_test", "second-order");
	app->Options()->SetStringValue("derivative_test_print_all", "yes");
	app->Options()->SetStringValue("check_derivatives_for_naninf", "yes");

	// Intialize the IpoptApplication and process the options
	ApplicationReturnStatus status;
	status = app->Initialize();
	if (status != Solve_Succeeded)
	{
		printf("\n\n*** Error during initialization!\n");
		return (int) status;
	}

	// Ask Ipopt to solve the problem
	status = app->OptimizeTNLP(mynlp);

	if (status == Solve_Succeeded)
	{
		printf("\n\n*** The problem solved!\n");
	}
	else
	{
		printf("\n\n*** The problem FAILED!\n");
	}

	// As the SmartPtrs go out of scope, the reference count
	// will be decremented and the objects will automatically 
	// be deleted.

	return (int) status;
}
*/

python::list solve(python::list edges, python::list nodes, int verbose=0)
{
	SmartPtr<BreakCopies> mynlp = new BreakCopies();

	mynlp->InitProblem(edges, nodes);

	SmartPtr<IpoptApplication> app = IpoptApplicationFactory();
	app->Options()->SetNumericValue("tol", 1e-9);
	app->Options()->SetStringValue("mu_strategy", "adaptive");
	app->Options()->SetIntegerValue("print_level", verbose);

	ApplicationReturnStatus status;

	status = app->Initialize();

	if (status != Solve_Succeeded)
	{
		throw runtime_error("Error: Initialization failed with error code " + lexical_cast<string>(status));
	}

	status = app->OptimizeTNLP(mynlp);

	if (status != Solve_Succeeded)
	{
		throw runtime_error("Error: Solve failed with error code " + lexical_cast<string>(status));
	}

	python::list copies;
	for (size_t i = 0; i < mynlp->mCopies.size(); i++)
	{
		copies.append(mynlp->mCopies[i]);
	}

	return copies;
}

BOOST_PYTHON_FUNCTION_OVERLOADS(solve_overloads, solve, 2, 3);

BOOST_PYTHON_MODULE(pybreakcopies)
{
	using namespace python;
	
	def("solve", solve, solve_overloads());
}


