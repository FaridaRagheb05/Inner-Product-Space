/* verifies whether a function defines an inner product on R^n by testing the three axioms:
    1. Positive Definiteness: <v,v> >= 0 for all v, with equality iff v = 0
    2. Symmetry: <u,v> = <v,u> for all u,v
    3. Linearity: <au+bv,w> = a<u,w> + b<v,w> for all u,v,w and scalars a,b
NOTE: Numerical testing cannot truly prove that a function is an inner product since the axioms must hold for ALL vectors
This program just performs randomized validation: 
it can definitively, so it can disprove an axiom (if any test fails), but passing all tests is only strong empirical evidence, but it is not a formal proof
*/

#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <string>
#include <functional>
#include <ctime>
#include <cstdlib>
#include <stdexcept>

using namespace std;


const double DEFAULT_EPSILON = 1e-9;

// Global state for current test
struct TestState {
    int positivity_tests = 0;
    int positivity_passed = 0;
    
    int symmetry_tests = 0;
    int symmetry_passed = 0;
    
    int linearity_tests = 0;
    int linearity_passed = 0;
    
    bool all_axioms_pass() const 
    {
        return (positivity_passed == positivity_tests) && (symmetry_passed == symmetry_tests) && (linearity_passed == linearity_tests) &&
               (positivity_tests > 0 && symmetry_tests > 0 && linearity_tests > 0);
    }
};

TestState test_state;

// helper functions

// Checks if two doubles are approximately equal within a specified epsilon tolerance
bool approx_equal(double a, double b, double epsilon) 
{
    return fabs(a - b) <= epsilon;
}


// Computes the squared norm of a vector (useful for checking posiive )
double vector_norm_squared(const vector<double>& v) 
{
    double sum = 0.0;
    for (double x : v) 
    {
        sum += x * x;
    }
    return sum;
}

// Inner Product Functions

// Example 1: Standard Dot Product
double standard_dot_product(const vector<double>& u, const vector<double>& v) 
{
    if (u.size() != v.size()) 
    {
        throw invalid_argument("Vector dimensions must match");
    }
    double sum = 0.0;
    for (size_t i = 0; i < u.size(); ++i) 
    {
        sum += u[i] * v[i];
    }
    return sum;
}

// Example 2: Weighted Dot Product
// <u,v>_w = sum(w_i * u_i * v_i) where all w_i > 0
vector<double> global_weights;  // populated by argument parser

double weighted_dot_product(const vector<double>& u, const vector<double>& v) {
    if (u.size() != v.size()) {
        throw invalid_argument("Vector dimensions must match");
    }
    if (u.size() != global_weights.size()) {
        throw invalid_argument("Weights dimension must match vector dimension");
    }
    double sum = 0.0;
    for (size_t i = 0; i < u.size(); ++i) {
        sum += global_weights[i] * u[i] * v[i];
    }
    return sum;
}
// Example 3: Missing Symmetry 
// <u,v> = u_1*v_2 - u_2*v_1 (only for dimension >= 2)
double antisymmetric_invalid(const vector<double>& u, const vector<double>& v) 
{
    if (u.size() < 2) 
    {
        throw invalid_argument("This example requires dimension >= 2");
    }
    return u[0] * v[1] - u[1] * v[0];
}

// Example 4: Negative Weights
// <u,v> = -u_1*v_1 + u_2*v_2
// Fails positive definiteness because <e_1, e_1> = -1
double negative_weights_invalid(const vector<double>& u, const vector<double>& v) 
{
    if (u.size() != v.size()) 
    {
        throw invalid_argument("Vector dimensions must match");
    }
    double sum = 0.0;
    for (size_t i = 0; i < u.size(); ++i) 
    {
        double weight = (i == 0) ? -1.0 : 1.0;
        sum += weight * u[i] * v[i];
    }
    return sum;
}

// Example 5: Nonlinear 
// <u,v> = (u_1*v_1)^2 + u_2*v_2
// fails linearity
double nonlinear_invalid(const vector<double>& u, const vector<double>& v) 
{
    if (u.size() != v.size()) 
    {
        throw invalid_argument("Vector dimensions must match");
    }
    if (u.size() < 2) 
    {
        throw invalid_argument("Nonlinear example requires dimension >= 2");
    }
    double first_term = (u.size() > 0) ? pow(u[0] * v[0], 2.0) : 0.0;
    double rest = 0.0;
    for (size_t i = 1; i < u.size(); ++i) 
    {
        rest += u[i] * v[i];
    }
    return first_term + rest;
}



// Test implementations

void print_header(const string& function_name, int dim, double epsilon, const string& weights_str = "") 
{
    cout << "\n";
    cout << "+" << string(78, '=') << "+\n";
    cout << "| " << left << setw(76) << function_name << " |\n";
    cout << "+---" << string(75, '-') << "+\n";
    cout << "| Dimension: " << dim << setw(65) << " " << "|\n";
    cout << "| Tolerance: " << scientific << setprecision(2) << epsilon
         << setw(63) << " " << "|\n";
    if (!weights_str.empty()) 
    {
        cout << "| Weights: " << left << setw(67) << weights_str << "|\n";
    }
    cout << "+" << string(78, '=') << "+\n\n";
}

// AXIOM 1: Positive Definiteness Test

void test_positivity(function<double(const vector<double>&, const vector<double>&)> inner_prod, int dim, double epsilon) 
{
    cout << "\nAXIOM 1: Positive Definiteness: <v,v> >= 0 for all v, and <v,v> = 0 iff v = 0\n";
    cout << string(80, '-') << "\n";

    int passed_count = 0;
    int prev_positivity_tests = test_state.positivity_tests;  // Baseline for local display

    // Test: Zero vector (required: <0,0> = 0)
    cout << "\nTest 1.1 (Zero vector): v = (0, 0, ..., 0)\n";
    vector<double> zero(dim, 0.0);
    double zero_result = inner_prod(zero, zero);
    bool zero_test_pass = approx_equal(zero_result, 0.0, epsilon);
    cout << "  <0, 0> = " << fixed << setprecision(10) << zero_result << "\n";
    cout << "  Expected: 0\n";
    cout << "  Difference: " << fabs(zero_result - 0.0) << "\n";
    cout << "  Result: " << (zero_test_pass ? "PASS" : "FAIL") << "\n";
    if (zero_test_pass) passed_count++;
    test_state.positivity_tests++;
    test_state.positivity_passed += (zero_test_pass ? 1 : 0);

    // Tests: Standard basis vectors
    // e_i isolates each component of the function
    // A negative weight on any component is guaranteed to fail here
    for (int i = 0; i < min(dim, 2); ++i) 
    {
        cout << "\nTest 1." << (i + 2) << " (Standard basis vector e" << (i + 1) << "):\n";
        vector<double> e(dim, 0.0);
        e[i] = 1.0;
        cout << "  v = e" << (i + 1) << "  (1 in position " << (i + 1) << ", 0 elsewhere)\n";
        double result = inner_prod(e, e);
        bool pass = (result >= -epsilon);
        cout << "  <e" << (i + 1) << ", e" << (i + 1) << "> = " << fixed << setprecision(10) << result << "\n";
        cout << "  Requirement: <v, v> > 0 for nonzero v\n";
        cout << "  Result: " << (pass ? "PASS (positive)" : "FAIL (non-positive)") << "\n";
        if (pass) passed_count++;
        test_state.positivity_tests++;
        test_state.positivity_passed += (pass ? 1 : 0);
    }

    // Tests: Random nonzero vectors (seed set in main)
    for (int i = 0; i < 3; ++i) 
    {
        cout << "\nTest 1." << (min(dim, 2) + 2 + i) << " (Random nonzero vector):\n";
        vector<double> v(dim);
        double norm_sq = 0.0;
        cout << "  v = (";
        for (int j = 0; j < dim; ++j) 
        {
            v[j] = (2.0 * rand() / RAND_MAX) - 1.0;
            norm_sq += v[j] * v[j];
            if (j < 3)
            {
                cout << fixed << setprecision(3) << v[j] << (j < min(dim-1, 2) ? ", " : "");
            }
        }
        if (dim > 3) cout << ", ...";
        cout << ")\n";

        if (norm_sq < 1e-14) 
        {
            cout << "  [Skipped: numerically zero]\n";
            continue;
        }

        double result = inner_prod(v, v);

        bool pass;

        if (norm_sq <= epsilon)
        {
            pass = approx_equal(result, 0.0, epsilon);
        }
        else
        {
            pass = (result > epsilon);
        }
        cout << "  <v, v> = " << fixed << setprecision(10) << result << "\n";
        cout << "  Requirement: <v, v> > 0 for nonzero v\n";
        cout << "  Result: " << (pass ? "PASS" : "FAIL");
        if (!pass) cout << " -> violates positive definiteness";
        cout << "\n";
        if (pass) passed_count++;
        test_state.positivity_tests++;
        test_state.positivity_passed += (pass ? 1 : 0);
    }

    int local_test_count = test_state.positivity_tests - prev_positivity_tests;
    cout << "\n" << string(80, '-') << "\n";
    cout << "Positivity: " << passed_count << "/" << local_test_count << " tests passed";
    cout << (passed_count == local_test_count ? " [PASS]" : " [FAIL]") << "\n";
}

// AXIOM 2: Symmetry Test

void test_symmetry(function<double(const vector<double>&, const vector<double>&)> inner_prod,int dim, double epsilon
) {
    cout << "\nAXIOM 2: Symmetry: <u,v> = <v,u> for all u, v\n";
    cout << string(80, '-') << "\n";

    int passed_count = 0;
    int prev_symmetry_tests = test_state.symmetry_tests;

    for (int test = 0; test < 7; ++test) 
    {
        cout << "\nTest 2." << (test + 1) << ":\n";
        
        vector<double> u(dim), v(dim);
        cout << "  u = (";
        for (int i = 0; i < dim; ++i) 
        {
            u[i] = (2.0 * rand() / RAND_MAX) - 1.0;
            if (i < 3) cout << fixed << setprecision(3) << u[i] << ", ";
        }
        if (dim > 3) cout << "...";
        cout << ")\n";
        
        cout << "  v = (";
        for (int i = 0; i < dim; ++i) 
        {
            v[i] = (2.0 * rand() / RAND_MAX) - 1.0;
            if (i < 3) cout << fixed << setprecision(3) << v[i] << ", ";
        }
        if (dim > 3) cout << "...";
        cout << ")\n";

        double uv = inner_prod(u, v);
        double vu = inner_prod(v, u);
        double diff = fabs(uv - vu);
        bool pass = approx_equal(uv, vu, epsilon);

        cout << "  <u, v> = " << fixed << setprecision(10) << uv << "\n";
        cout << "  <v, u> = " << fixed << setprecision(10) << vu << "\n";
        cout << "  Difference: " << scientific << setprecision(2) << diff << "\n";
        cout << "  Tolerance: " << scientific << setprecision(2) << epsilon << "\n";
        cout << "  Result: " << (pass ? "PASS" : "FAIL");
        if (!pass) cout << " -> violates symmetry axiom";
        cout << "\n";
        
        if (pass) passed_count++;
        test_state.symmetry_tests++;
        test_state.symmetry_passed += (pass ? 1 : 0);
    }

    int local_sym_count = test_state.symmetry_tests - prev_symmetry_tests;
    cout << "\n" << string(80, '-') << "\n";
    cout << "Symmetry: " << passed_count << "/" << local_sym_count << " tests passed";
    cout << (passed_count == local_sym_count ? " [PASS]" : " [FAIL]") << "\n";
}

// AXIOM 3:  Linearity Test

void test_linearity(function<double(const vector<double>&, const vector<double>&)> inner_prod, int dim, double epsilon) 
{
    cout << "\nAXIOM 3: Linearity: <au+bv,w> = a<u,w> + b<v,w> for all u,v,w,a,b\n";
    cout << string(80, '-') << "\n";

    int passed_count = 0;
    int prev_linearity_tests = test_state.linearity_tests;

    // Test with specific scalars: positive, negative, and ZERO
    vector<pair<double, double>> scalar_pairs = 
    {

        {2.0, 3.0},
        {-1.0, 2.5},
        {0.0, 1.0}, // Zero scalar test (a=0)
        {1.5, 0.0}, // Zero scalar test (b=0)
        {-0.5, -2.0},
        {3.3, -1.7},
        {-2.2, 0.0}
    };

    for (size_t test = 0; test < scalar_pairs.size(); ++test) 
    {
        double a = scalar_pairs[test].first;
        double b = scalar_pairs[test].second;

        cout << "\nTest 3." << (test + 1) << " (a=" << fixed << setprecision(1) << a << ", b=" << b << "):\n";

        vector<double> u(dim), v(dim), w(dim);
        for (int i = 0; i < dim; ++i) 
        {
            u[i] = (2.0 * rand() / RAND_MAX) - 1.0;
            v[i] = (2.0 * rand() / RAND_MAX) - 1.0;
            w[i] = (2.0 * rand() / RAND_MAX) - 1.0;
        }

        // Compute au + bv
        vector<double> linear_comb(dim);
        for (int i = 0; i < dim; ++i) 
        {
            linear_comb[i] = a * u[i] + b * v[i];
        }

        // LHS: <au+bv, w>
        double lhs = inner_prod(linear_comb, w);

        // RHS: a<u,w> + b<v,w>
        double uw = inner_prod(u, w);
        double vw = inner_prod(v, w);
        double rhs = a * uw + b * vw;

        double diff = fabs(lhs - rhs);
        bool pass = approx_equal(lhs, rhs, epsilon);

        cout << "  <au+bv, w> = " << fixed << setprecision(10) << lhs << "\n";
        cout << "  a<u,w> + b<v,w> = " << fixed << setprecision(10) << rhs << "\n";
        cout << "  Difference: " << scientific << setprecision(2) << diff << "\n";
        cout << "  Tolerance: " << scientific << setprecision(2) << epsilon << "\n";
        cout << "  Result: " << (pass ? "PASS" : "FAIL");
        if (!pass) cout << " -> violates linearity axiom";
        cout << "\n";

        if (pass) passed_count++;
        test_state.linearity_tests++;
        test_state.linearity_passed += (pass ? 1 : 0);
    }

    int local_lin_count = test_state.linearity_tests - prev_linearity_tests;
    cout << "\n" << string(80, '-') << "\n";
    cout << "Linearity: " << passed_count << "/" << local_lin_count << " tests passed";
    cout << (passed_count == local_lin_count ? " [PASS]" : " [FAIL]") << "\n";
}

// summary

void print_summary() 
{
    cout << "\n\n";
    cout << "+" << string(78, '=') << "+\n";
    cout << "| " << left << setw(76) << "SUMMARY" << " |\n";
    cout << "+---" << string(75, '-') << "+\n";

    bool positivity_pass = (test_state.positivity_passed == test_state.positivity_tests);
    bool symmetry_pass = (test_state.symmetry_passed == test_state.symmetry_tests);
    bool linearity_pass = (test_state.linearity_passed == test_state.linearity_tests);

    cout << "| Positivity (Axiom 1): " << test_state.positivity_passed << "/" << test_state.positivity_tests << " tests passed";
    cout << setw(46 - to_string(test_state.positivity_passed).length() - to_string(test_state.positivity_tests).length()) << " " << (positivity_pass ? "[PASS]" : "[FAIL]") << " |\n";

    cout << "| Symmetry (Axiom 2):   " << test_state.symmetry_passed << "/" << test_state.symmetry_tests << " tests passed";
    cout << setw(46 - to_string(test_state.symmetry_passed).length() - to_string(test_state.symmetry_tests).length()) << " " << (symmetry_pass ? "[PASS]" : "[FAIL]") << " |\n";

    cout << "| Linearity (Axiom 3):  " << test_state.linearity_passed << "/" << test_state.linearity_tests << " tests passed";
    cout << setw(46 - to_string(test_state.linearity_passed).length() -  to_string(test_state.linearity_tests).length()) << " " << (linearity_pass ? "[PASS]" : "[FAIL]") << " |\n";

    cout << "+---" << string(75, '-') << "+\n";

    if (test_state.all_axioms_pass()) 
    {
        cout << "| " << left << setw(76) << "VERDICT: VALID INNER PRODUCT" << " |\n";
    } 
    else 
    {
        cout << "| " << left << setw(76) << "VERDICT: NOT A VALID INNER PRODUCT" << " |\n";
        cout << "| ";
        if (!positivity_pass) cout << "Failed: Positivity. ";
        if (!symmetry_pass) cout << "Failed: Symmetry. ";
        if (!linearity_pass) cout << "Failed: Linearity. ";
        cout << setw(76 - ((!positivity_pass ? 19 : 0) + (!symmetry_pass ? 17 : 0) + (!linearity_pass ? 17 : 0))) << " " << "|\n";
    }

    cout << "+---" << string(75, '-') << "+\n";
    int total_tests = test_state.positivity_tests + test_state.symmetry_tests + test_state.linearity_tests;
    cout << "| Total tests performed: " << total_tests;    
    cout << setw(53) << " " << "|\n";
    cout << "+" << string(78, '=') << "+\n\n";
}

// CLI Parsing

void print_usage(const string& program_name) 
{
    cout << "\nUsage: " << program_name << " [options]\n\n";
    cout << "Options:\n";
    cout << "  -d, --dimension DIM      Dimension of vector space (default: 3)\n";
    cout << "  -f, --function TYPE      Function type: dot, weighted, negative, nonlinear, antisymmetric\n";
    cout << "                           (default: dot)\n";
    cout << "  -w, --weights W1 W2 ... Weights for weighted inner product\n";
    cout << "  -t, --tolerance EPS      Floating-point tolerance (default: 1e-9)\n";
    cout << "  -h, --help               Show this help message\n\n";
    cout << "Examples:\n";
    cout << "  " << program_name << " -d 3 -f dot\n";
    cout << "  " << program_name << " -d 4 -f weighted -w 1 2 3 4\n";
    cout << "  " << program_name << " -d 2 -f nonlinear\n";
}

struct Arguments 
{
    int dimension = 3;
    string function_type = "dot";
    vector<double> weights;
    double tolerance = DEFAULT_EPSILON;
    bool help_requested = false;
};

Arguments parse_arguments(int argc, char* argv[]) 
{
    Arguments args;

    for (int i = 1; i < argc; ++i) 
    {
        string arg = argv[i];

        if (arg == "-h" || arg == "--help") 
        {
            args.help_requested = true;
            return args;
        } 
        else if (arg == "-d" || arg == "--dimension") 
        {
            if (i + 1 < argc) 
            {
                args.dimension = stoi(argv[++i]);
                if (args.dimension < 1 || args.dimension > 100) 
                {
                    throw invalid_argument("Dimension must be between 1 and 100");
                }
            }
        } 
        else if (arg == "-f" || arg == "--function") 
        {
            if (i + 1 < argc) 
            {
                args.function_type = argv[++i];
            }
        } 
        else if (arg == "-w" || arg == "--weights") 
        {
            while (i + 1 < argc && argv[i + 1][0] != '-') 
            {
                args.weights.push_back(stod(argv[++i]));
            }
        } else if (arg == "-t" || arg == "--tolerance") 
        {
            if (i + 1 < argc) 
            {
                args.tolerance = stod(argv[++i]);
                if (args.tolerance <= 0) 
                {
                    throw invalid_argument("Tolerance must be positive");
                }
            }
        }
    }

    return args;
}


int main(int argc, char* argv[]) 
{
    try 
    {
        Arguments args = parse_arguments(argc, argv);

        if (args.help_requested) {
            print_usage(argv[0]);
            return 0;
        }

        // Seed random number generator 
        srand(static_cast<unsigned>(time(nullptr)));

        // Validate dimension for examples that require dim >= 2
        if ((args.function_type == "antisymmetric" || args.function_type == "nonlinear") && args.dimension < 2) 
        {
            throw invalid_argument("Function type '" + args.function_type + "' requires dimension >= 2");
        }

        // Select function
        function<double(const vector<double>&, const vector<double>&)> inner_prod;
        string function_name;
        string weights_str = "";

        if (args.function_type == "dot") 
        {
            inner_prod = standard_dot_product;
            function_name = "Standard Dot Product: <u,v> = u₁v₁ + u₂v₂ + ... + uₙvₙ";
        } 
        else if (args.function_type == "weighted") 
        {
            if (args.weights.empty()) 
            {
                // Default weights
                for (int i = 1; i <= args.dimension; ++i) 
                {
                    args.weights.push_back(static_cast<double>(i));
                }
            }
            if (args.weights.size() != static_cast<size_t>(args.dimension)) 
            {
                throw invalid_argument("Number of weights must match dimension");
            }
            // Weights must all be positive for a valid inner product
            bool has_nonpositive = false;

            for (double w : args.weights) 
            {
                if (w <= 0.0) has_nonpositive = true;
            }

            if (has_nonpositive) 
            {
                cout << "\nWARNING: Non-positive weights detected -> expected to FAIL positivity axiom\n";
            }
            global_weights = args.weights;
            inner_prod = weighted_dot_product;
            
            function_name = "Weighted Dot Product: <u,v> = Σ wᵢuᵢvᵢ (wᵢ > 0)";
            weights_str = "[";
            for (size_t i = 0; i < args.weights.size(); ++i) 
            {
                if (i > 0) weights_str += ", ";
                weights_str += to_string(args.weights[i]);
            }
            weights_str += "]";
        } 
        else if (args.function_type == "negative") 
        {
            inner_prod = negative_weights_invalid;
            function_name = "INVALID: Negative Weights: <u,v> = -u₁v₁ + u₂v₂ + ...";
        } 
        else if (args.function_type == "nonlinear") 
        {
            inner_prod = nonlinear_invalid;
            function_name = "INVALID: Nonlinear: <u,v> = (u₁v₁)² + u₂v₂ + ...";
        } 
        else if (args.function_type == "antisymmetric") 
        {
            inner_prod = antisymmetric_invalid;
            function_name = "INVALID: Antisymmetric: <u,v> = u₁v₂ - u₂v₁";
        } 
        else 
        {
            throw invalid_argument("Unknown function type: " + args.function_type);
        }

        // Print header
        print_header(function_name, args.dimension, args.tolerance, weights_str);

        // run tests
        test_positivity(inner_prod, args.dimension, args.tolerance);
        test_symmetry(inner_prod, args.dimension, args.tolerance);
        test_linearity(inner_prod, args.dimension, args.tolerance);

        // print summary
        print_summary();

        return 0;

    } 
    catch (const exception& e) 
    {
        cerr << "\nERROR: " << e.what() << "\n";
        print_usage(argv[0]);
        return 1;
    }
}