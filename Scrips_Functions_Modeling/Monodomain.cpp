/*
 ============================================================================

 .______    _______     ___   .___________.    __  .___________.
 |   _  \  |   ____|   /   \  |           |   |  | |           |
 |  |_)  | |  |__     /  ^  \ `---|  |----`   |  | `---|  |----`
 |   _  <  |   __|   /  /_\  \    |  |        |  |     |  |
 |  |_)  | |  |____ /  _____  \   |  |        |  |     |  |
 |______/  |_______/__/     \__\  |__|        |__|     |__|

 BeatIt - code for cardiovascular simulations
 Copyright (C) 2016 Simone Rossi

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ============================================================================
 */
#include <chrono>

#include "libmesh/exodusII_io.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_modification.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/mesh_tools.h"
#include "libmesh/parallel_mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/elem.h"
#include "libmesh/analytic_function.h"

#include "libmesh/equation_systems.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/dof_map.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/quadrature_gauss.h"

#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/explicit_system.h"
#include "libmesh/dtk_solution_transfer.h"

#include <sys/stat.h>
#include "libmesh/getpot.h"

#include "libmesh/vtk_io.h"
#include <iomanip>      // std::setw

#include "libmesh/mesh_serializer.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <iterator>
#include <numeric>

#include <algorithm>

#include "petscksp.h"
PetscErrorCode  KSPSetNullSpace(KSP ksp,MatNullSpace nullsp);




using namespace libMesh;

double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::floor(value * multiplier) / multiplier;
}

std::string to_string_with_precision(const double a_value, const int n = 4)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}

// Domain numbering
enum class Subdomain : unsigned int
{
    TISSUE = 1, FIBROSIS = 0
};
// Time integrators
enum class TimeIntegrator : unsigned int
{
    SBDF1 = 1, SBDF2 = 2, SBDF3 = 3, EXPLICIT_EXTRACELLULAR = 4, EXPLICIT_INTRACELLULAR = 5, SEMI_IMPLICIT = 6, SEMI_IMPLICIT_HALF_STEP = 7
};
enum class EquationType : unsigned int
{
    PARABOLIC = 1
};

// Store here details about time stepping
struct TimeData
{
    // constructor
    TimeData(GetPot &data);
    // show parameters
    void print();
    // parameters
    double time;      // Time
    double end_time;      // Final Time
    double dt;            // Timestep
    int timestep;         // timestep coutner
    int export_timesteps; // interval at which export the solution
};


double initial_condition_V(const libMesh::Point& p, const double time);
double exact_solution_V(const libMesh::Point& p, const double time);



// Save vector as matrix type in Octave (.mat) stm.
/*
template <typename T>
std::ostream& save_vector_as_matrix( const std::string& name, const std::vector<T>& matrix, std::ofstream& stm )
{
  stm << name_tag << name << '\n' << type_tag << '\n' << rows_tag << '\n'
      << "# columns: " << matrix.size() << '\n' ;

    std::copy( matrix.begin(), matrix.end(), std::ostream_iterator<T>( stm, " " ) ) ;
  return stm << "\n\n\n" ;
}
*/


// Read mesh as specified in the input file
void read_mesh(const GetPot &data, libMesh::ParallelMesh &mesh, libMesh::SerialMesh &mesh2);
// read BC sidesets from string: e.g. bc = "1 2 3 5", or bc = "1, 55, 2, 33"
void read_bc_list(std::string &bc, std::set<int> &bc_sidesets);
// Assemble matrices
void assemble_matrices(libMesh::EquationSystems &es, const TimeData &timedata, TimeIntegrator time_integrator, libMesh::Order p_order, const GetPot &data, int rankLoc);

//double calculate_PhiE(libMesh::EquationSystems &es, const TimeData &timedata, TimeIntegrator time_integrator, libMesh::Order p_order, const GetPot &data);
// Solve ionic model / evaluate ionic currents

// void assemble RHS
void assemble_rhs(  libMesh::EquationSystems &es,
                    const TimeData &timedata,
                    TimeIntegrator time_integrator,
                    libMesh::Order p_order,
                    const GetPot &data,
                    EquationType type = EquationType::PARABOLIC );


void ForcingTermConvergence(libMesh::EquationSystems &es, const double dt, const double Beta, const double Cm, const double SigmaSI, const double SigmaSE, const double SigmaBath, const double SigmaTorso, const double CurTime, const double xDim, const double v0, const double v1, const double v2, const double kcubic, const std::string integrator );

void init_cd_exact (EquationSystems & es, const double xDim, std::string integrator,  const GetPot &data);

void SolverRecovery (EquationSystems & es, const GetPot &data, TimeData& datatime);
//void SolverRecoveryMono (EquationSystems & es, const GetPot &data, TimeData& datatime);

double exact_solutionV_all (const double x,
              const double y,
              const int typePlace,
              const double time,
              const double z,
              const double xDim){

    double value = 0.;
    //double pi = acos(-1);
    double y0 = -.5;
    double c = .125;
    double sigma = .125;
    double delta = .1;
    double alpha = 50.;
    double a1 = 8.0/(xDim*xDim/4.0);
    double b1 = 16.0/(xDim/2.0);
    double c1 = 8.0;
    double a2 = -8.0/(xDim*xDim/4.0);
    double b2 = -8.0/(xDim/2.0);
    double c2 = -1.0;


    double dummy = y0 - y + c*time;

    double V = tanh(alpha*(dummy))*.5 + .5;
    double Vx = 0.;
    double Vxx = 0.;
    double Vy = .5*alpha*(pow(tanh(alpha*(dummy)),2) - 1.0);
    double Vyy = alpha*alpha*tanh(alpha*(dummy))*(pow(tanh(alpha*(dummy)),2) - 1.0);
    double Vxy = 0.0;
    double Vt = -(.5*alpha*c*(pow(tanh(alpha*dummy),2) - 1.0));

    double sigma6 = pow(sigma,6);
    double innerPow1 = delta - dummy;
    double innerPow1y = 1.0;
    double innerPow1t = -c;
    double innerPow2 = delta + dummy;
    double innerPow2y = -1.0;
    double innerPow2t = c;

    double exp1 = exp(-pow((innerPow1),6)/sigma6);
    double exp1y = -6.0*exp1*pow(innerPow1,5)*innerPow1y/sigma6;
    double exp1yy =  -30.0*exp1*pow(innerPow1,4)*innerPow1y*innerPow1y/sigma6 - 6.0*exp1y*pow(innerPow1,5)*innerPow1y/sigma6;
    double exp2 = exp(-pow((innerPow2),6)/sigma6);
    double exp2y = -6.0*exp2*pow(innerPow2,5)*innerPow2y/sigma6;
    double exp2yy =  -30.0*exp2*pow(innerPow2,4)*innerPow2y*innerPow2y/sigma6 - 6.0*exp2y*pow(innerPow2,5)*innerPow2y/sigma6;

    double f = exp1 - exp2;
    double fy = exp1y - exp2y;
    double fyy = exp1yy - exp2yy;
    if(x > 0){
      b1 = -b1;
      b2 = -b2;
    }

    double g1 = a1*x*x + b1*x + c1;
    double g1x = 2.0*a1*x + b1;
    double g1xx = 2.0*a1;
    double g2 = a2*x*x + b2*x + c2;
    double g2x = 2.0*a2*x + b2;
    double g2xx = 2.0*a2;

    double Ve = f;
    double Vex = 0.0;
    double Vexx = 0.0;
    double Vey = fy;
    double Veyy = fyy;
    double g, gx, gxx;

    if(x <= (-3.0/4.0)*xDim/2.0){
      g = g1;
      gx = g1x;
      gxx = g1xx;
    }
    else if(x > (-3.0/4.0)*xDim/2.0 && x < -xDim/4.0){
      g = g2;
      gx = g2x;
      gxx = g2xx;
    }
    else if(x >= (3.0/4.0)*xDim/2.0){
      g = g1;
      gx = g1x;
      gxx = g1xx;
    }
    else{
      g = g2;
      gx = g2x;
      gxx = g2xx;
    }


    double Vb = f*g;
    double Vbx = f*gx;
    double Vby = fy*g;
    double Vbxx = f*gxx;
    double Vbyy = fyy*g;
    double Vbxy = fy*gx;



    // Rest of tissue part: V
    if(typePlace == 2){
      //libMesh::out << .5 -.5*tanh(3.*pi*std::abs(.7)/(yDim/2.0) - pi - (1.0/5.0)*time)<< " " << std::abs(.7) << std::endl;
      //value = .5 -.5*tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time);
      value = V;
      //libMesh::out << pow(.5,2) << std::endl;
    }
    // Rest of tissue part: V time derivative
  else if(typePlace == 22){
    //value = (-1./10.0)*(1 - pow(tanh(-3.*pi*(std::abs(y))/(yDim/2.0) + pi + (1.0/5.0)*time),2));
    value = Vt;
  }
    // Rest of tissue part: V first derivative
  else if(typePlace == 20){
    double dummy2 = tanh(dummy);
    value = Vy;
  }
    // Rest of tissue part: V second derivative
  else if(typePlace == 200){
    //value = 9.*(pow(pi,2))*(1. - pow(tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time),2))*(tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time));
    value = Vyy;
  }



    // Rest of tissue part: Ve and Vb
    else if(typePlace == 3){
      //value = -.5 +.5*tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time);
      //double expo1 = pow(((-dummy+delta)/sigma),6);
      //double expo2 = pow(((dummy+delta)/sigma),6);
      //value = exp(expo1) - exp(expo2);
      value = Ve;

    }

    else if(typePlace == 4){
     value = Vb;
    }


    // Rest of tissue part: Ve and Vb first derivative
    else if(typePlace == 30){
      //value = (-.25*exp(-pow(((x - time*.15 + 2.0)/(2.0/8.0)),2)) + .25*exp(-pow(((x - time*.15 + 1.5)/(2.0/8.0)),2)))*(1.0 - abs(y)*2.0);
      //value = (3.*pi*y/(2.*(std::abs(y))))*(1. - pow(tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time),2));
      //double expo1 = pow(((dummy+delta)/sigma),5);
      //double expo2 = pow(((dummy-delta)/sigma),5);
      value = Vey;
    }
    else if(typePlace == 40){

    value = Vby;
  }

    // Rest of tissue part: Ve and Vb second derivative
    else if(typePlace == 300){
      //value = (-.25*exp(-pow(((x - time*.15 + 2.0)/(2.0/8.0)),2)) + .25*exp(-pow(((x - time*.15 + 1.5)/(2.0/8.0)),2)))*(1.0 - abs(y)*2.0);
      //value = -9.*(pow(pi,2))*(1. - pow(tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time),2))*(tanh(3.*pi*(std::abs(y))/(yDim/2.0) - pi - (1.0/5.0)*time));
      value = Veyy;
      //libMesh::out << value << "     place: " << typePlace  << " y: " << y << std::endl;
    }
    //summation of both the derivative in respect to y and x of Vb
    else if(typePlace == 4004){
      value = Vbyy + Vbxx;
    }




    //libMesh::out << value << "     place: " << typePlace<< std::endl;
    return value;
}


double CalculateF(const double x,
          const double y,
          const int typeEqn,
          const double time,
          const double z,
          const double dt1,
          const double Beta,
          const double Cm,
          const double SigmaSI,
          const double SigmaSE,
          const double SigmaBath,
          const double SigmaTorso,
          const double xDim,
          const double v0,
          const double v1,
          const double v2,
          const double kcubic){

  double value;

  //if(alphaTime != 1){

    //dt1 = alphaTime*(yDim/numEl)/(.125);
  //}

  // 0 is for V, 1 for Ve, 2 is for Vbath in blood, 3 is for Vbath in torso
  //USUALLY BLOOD BATH ON ONE SIDE OF THE TISSUE AND TORSO ON THE OTHER SIDE

  //V equation, parabolic, so Fv
  if(typeEqn == 0){
    double prevTime = time - dt1;
    double prevprevTime = time - dt1*2.0;
    double u_current = exact_solutionV_all(x,y, 2, time, 0.0, xDim);
    double u_prev = exact_solutionV_all(x,y, 2, prevTime, 0.0, xDim);
    double u_prev_prev = exact_solutionV_all(x,y, 2, prevprevTime, 0.0, xDim);

    value = (1*Cm)*exact_solutionV_all(x,y, 22, time, 0.0, xDim) - 1.*( -kcubic*((u_current - v0)*(u_current - v1)*(u_current - v2)) ) - SigmaSI*exact_solutionV_all(x,y, 200, time, 0.0, xDim)/Beta - SigmaSI*exact_solutionV_all(x,y, 300, time, 0.0, xDim)/Beta;

  }
  //Ve equation, elliptic, so Fve
  else if(typeEqn == 1){
    value = -SigmaSI*exact_solutionV_all(x,y, 200, time, 0.0, xDim)/Beta - (SigmaSI + SigmaSE)*exact_solutionV_all(x,y, 300, time, 0.0, xDim)/Beta;
  }
  //Bath equation, Vb, so Fvb
  else if(typeEqn == 2){
    value = -SigmaBath*exact_solutionV_all(x,y, 4004, time, 0.0, xDim)/Beta;
  }
  else if(typeEqn == 3){
      value = -SigmaTorso*exact_solutionV_all(x,y, 4004, time, 0.0, xDim)/Beta;
    }

  return value;
}




double HofX(double uvalue, double otherValue){
  double Hvalue;

  if(uvalue - otherValue > 0){
    Hvalue = 1.;
  }
  else{
    Hvalue = 0.;
  }

  return Hvalue;
}




using namespace std;

int main(int argc, char **argv)
{

    // Read input file
    GetPot cl(argc, argv);
    std::string datafile_name = cl.follow("data.beat", 2, "-i", "--input");
    GetPot data(datafile_name);
    // Use libMesh stuff without writing libMesh:: everytime
    using namespace libMesh;
    // Initialize libMesh
    LibMeshInit init(argc, argv, MPI_COMM_WORLD);

    chrono::steady_clock sc;   // create an object of `steady_clock` class
    auto start = sc.now();     // start timer

    //#if !defined(LIBMESH_TRILINOS_HAVE_DTK)
      // Skip this example (and use a different return code) if libMesh
      // was compiled without Trilinos+DTK support.
      //libmesh_example_requires(true, "--enable-trilinos");

    //#else

    // Create folde in which we save the output of the solution
    std::string output_folder = data("output_folder", "Output");
    std::string stimulus_type = data("stimulus_type", "Transmembrane");
    double IstimD = data("stimulus_duration",2.);
    bool exportVe_time = data("exportVe_time",false);
    bool exportVe_mono = data("exportVe_mono", false);
    bool exportVMAT_time = data("exportVMAT_time", false);
    double exportVe_timeX = data("exportVe_timeX", 0.);
    double exportVe_timeY = data("exportVe_timeY", 0.);
    double exportVe_timeZ = data("exportVe_timeZ", 0.);

    bool checkSoFar = data("checkSoFar", true);

    bool FK_in_mV = data("FK_in_mV", false);

    bool convergence_test_save = data("convergence_test_save", false);
    bool monodomainType = data("monodomainType", false);
    std::string cellModel = data("cellModel", "Fenton-Karma");

    struct stat out_dir;
    if (stat(&output_folder[0], &out_dir) != 0)
    {
        if (init.comm().rank() == 0)
        {
            mkdir(output_folder.c_str(), 0777);
        }
    }

    // Create empty mesh
    ParallelMesh mesh(init.comm());
    SerialMesh mesh2(init.comm());

    int nelxS = data("nelxS", 10);
    int nelyS = data("nelyS", 10);
    int nelzS = data("nelzS", 10);

    bool extraMesh = data("extraMesh", false);
    double dXmesh = data("dXmesh", 0.015);
    double minxS = exportVe_timeX - dXmesh*(nelxS/2.);
    double minyS = exportVe_timeY - dXmesh*(nelyS/2.);
    double minzS = exportVe_timeZ - dXmesh*(nelzS/2.);
    double maxxS = exportVe_timeX + dXmesh*(nelxS/2.);
    double maxyS = exportVe_timeY + dXmesh*(nelyS/2.);
    double maxzS = exportVe_timeZ + dXmesh*(nelzS/2.);


    bool prescribedFibers = data("prescribedFibers",false);

    read_mesh(data, mesh, mesh2);
    // output the details about the mesh
    mesh.print_info();

    // The dimension that we are running
    const unsigned int dim = mesh.mesh_dimension();

    // Time stuff
    TimeData datatime(data);
    datatime.print();

    // Pick time integrator
    // This will define if we use 1 or 2 systems
    std::map<std::string, TimeIntegrator> time_integrator_map =
    {
    { "SBDF1", TimeIntegrator::SBDF1 },
    { "SBDF2", TimeIntegrator::SBDF2 },
    { "SBDF3", TimeIntegrator::SBDF3 },
    { "EXPLICIT_EXTRACELLULAR", TimeIntegrator::EXPLICIT_EXTRACELLULAR },
    { "EXPLICIT_INTRACELLULAR", TimeIntegrator::EXPLICIT_INTRACELLULAR },
    { "SEMI_IMPLICIT", TimeIntegrator::SEMI_IMPLICIT },
    { "SEMI_IMPLICIT_HALF_STEP", TimeIntegrator::SEMI_IMPLICIT_HALF_STEP } };

    std::string integrator = data("integrator", "SBDF1");
    auto it = time_integrator_map.find(integrator);
    TimeIntegrator time_integrator = it->second;
    bool using_implicit_time_integrator = false;
    if (time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3)
    {
        using_implicit_time_integrator = true;
    }
    bool convergence_test = data("convergence_test", false);
    bool ICstim = data("ICstim", false);
    bool ICconditions = data("ICconditions", false);

    if (time_integrator == TimeIntegrator::SEMI_IMPLICIT || time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP)
    {
        libMesh::out << "SBDF1 and SEMI-IMPLICIT are the same thing. SEMI-IMPLICIT HALF-STEP no longer valid since Monodomain doesn't need equation decoupling." << std::endl;
        exit(1);
    }

    // Create libMesh Equations systems
    // This will hold the mesh and create the corresponding
    // finite element spaces
    libMesh::EquationSystems es(mesh);
    //DTKSolutionTransfer dtk_transfer(init.comm());

    // Define tissue active subdomain for subdomain restricted variables
    std::set < libMesh::subdomain_id_type > tissue_subdomains, fibrosis_subdomains;
    tissue_subdomains.insert(static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE));
    fibrosis_subdomains.insert(static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS));
    // Define finite element ORDER
    int p_order = data("p_order", 1);
    bool fibrosisBool = data("fibrosisBool", false);
    std::string fibrosisMethod = data("fibrosisMethod","Conductivities");
    std::cout << "P order: " << p_order << std::endl;
    Order order = FIRST;

    if (p_order == 2)
        order = SECOND;

    std::set<std::string> exported_variables;
    //equation_systems.parameters.set<DenseMatrix <Number> >("rlocal");
    
    // Bidomain System
    //switch (time_integrator)
    //{
        //case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
        //case TimeIntegrator::EXPLICIT_INTRACELLULAR:
        //case TimeIntegrator::SEMI_IMPLICIT:
        //case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
        //{
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.add_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            TransientExplicitSystem & recovery = es.add_system <TransientExplicitSystem> ("Recovery");

            //TransientExplicitSystem & fibrosis = es.add_system <TransientExplicitSystem> ("Fibrosis");

            std::cout << "Using element of order p = " << order << std::endl;
            //libMesh::out << "ALL GOOD TILL HERE 0" << std::endl;
            if(convergence_test_save){
                LinearImplicitSystem & solutionSystem = es.add_system <LinearImplicitSystem> ("Solution");
                solutionSystem.add_variable("V_exact",order, LAGRANGE, &tissue_subdomains);
            }

            if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
              LinearImplicitSystem & fibrosis_system = es.add_system <LinearImplicitSystem> ("Fibrosis");
              fibrosis_system.add_variable("Fibrosis_region",order, LAGRANGE, &tissue_subdomains);
              parabolic_system.add_vector("Fibrosis_Check");
            }

            if(exportVMAT_time){
                LinearImplicitSystem & LATSystem = es.add_system <LinearImplicitSystem> ("LATSystem");
                LATSystem.add_variable("LAT",order, LAGRANGE, &tissue_subdomains);
            }

           

            parabolic_system.add_variable("V", order, LAGRANGE, &tissue_subdomains);
            //parabolic_system.add_variable("Vfib", order, LAGRANGE, &fibrosis_subdomains);
            parabolic_system.add_matrix("Ki");
            parabolic_system.add_matrix("K");
            parabolic_system.add_matrix("Kmono");
            parabolic_system.add_matrix("Ke");
            parabolic_system.add_matrix("Kg");
            parabolic_system.add_matrix("M");
            parabolic_system.add_vector("ML");
            parabolic_system.add_vector("R_dist");
            parabolic_system.add_vector("In");
            parabolic_system.add_vector("FV"); // parabolic forcing term for exact solution
            parabolic_system.add_vector("F"); // parabolic forcing term for exact solution
            parabolic_system.add_vector("aux1"); // auxilliary vector for assembling the RHS
            parabolic_system.add_vector("ForcingV");
            parabolic_system.add_vector("ForcingConv");
            parabolic_system.add_vector("Inm1");
            parabolic_system.add_vector("Inm2");
            parabolic_system.add_vector("Vnm1");
            parabolic_system.add_vector("Vnm2");
            parabolic_system.add_vector("aux2"); // auxilliary vector for assembling the RHS
            parabolic_system.add_vector("aux3"); // auxilliary vector for assembling the RHS
            parabolic_system.add_vector("Vcurr"); // auxilliary vector for assembling the RHS
            parabolic_system.add_vector("Vprev"); // auxilliary vector for assembling the RHS
            parabolic_system.add_vector("dVmdt_max");
            parabolic_system.add_vector("dVmdt_max_prev");
            parabolic_system.add_vector("dVmdt_max_check");
            parabolic_system.add_vector("LATvec");


            if(exportVe_mono){
              parabolic_system.add_vector("Beta_Im");
              
            }

            if(cellModel.compare(0,12,"Fenton-Karma") == 0){
              recovery.add_variable("v", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("w", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("s", order, LAGRANGE, &tissue_subdomains);
              recovery.add_vector("v_prev");
              recovery.add_vector("w_prev");
              recovery.add_vector("s_prev");
              recovery.add_vector("v_prev_prev");
              recovery.add_vector("w_prev_prev");
              recovery.add_vector("s_prev_prev");   
            }
            else if(cellModel.compare(0,12,"Courtemanche") == 0){
              recovery.add_variable("v", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("w", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("xs", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("u", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("m", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("h", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("d", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("xr", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Nai", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Ki", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Carel", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("oa", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("ui", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("oi", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("f", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("j", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Cai", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Caup", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("ua", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("fCa", order, LAGRANGE, &tissue_subdomains); 

              //BacNav variables
              recovery.add_variable("mBN", order, LAGRANGE, &tissue_subdomains); 
              recovery.add_variable("hBN", order, LAGRANGE, &tissue_subdomains); 
            }

            if(prescribedFibers){
              
              LinearImplicitSystem & fiberSoln = es.add_system <LinearImplicitSystem> ("fiberSolution");
              fiberSoln.add_variable("fibx",order, LAGRANGE);
              fiberSoln.add_variable("fiby",order, LAGRANGE);
              fiberSoln.add_variable("fibz",order, LAGRANGE);

              fiberSoln.add_variable("sheetx",order, LAGRANGE);
              fiberSoln.add_variable("sheety",order, LAGRANGE);
              fiberSoln.add_variable("sheetz",order, LAGRANGE);

              fiberSoln.add_variable("crossfx",order, LAGRANGE);
              fiberSoln.add_variable("crossfy",order, LAGRANGE);
              fiberSoln.add_variable("crossfz",order, LAGRANGE);
              
            }
            

            es.init();
            //libMesh::out << "ALL GOOD TILL HERE 1" << std::endl;

            //break;
        //}
        /*
        case TimeIntegrator::SBDF1:
        case TimeIntegrator::SBDF2:
        case TimeIntegrator::SBDF3:
        default:
        {
            libMesh::TransientLinearImplicitSystem &system = es.add_system < libMesh::TransientLinearImplicitSystem > ("bidomain");
            std::cout << "Using element of order p = " << order << std::endl;

            TransientExplicitSystem & recovery = es.add_system <TransientExplicitSystem> ("Recovery");

            if(convergence_test_save){
                LinearImplicitSystem & solutionSystem = es.add_system <LinearImplicitSystem> ("Solution");
                solutionSystem.add_variable("V_exact",order, LAGRANGE, &tissue_subdomains);
            } 

            if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
              LinearImplicitSystem & fibrosis_system = es.add_system <LinearImplicitSystem> ("Fibrosis");
              fibrosis_system.add_variable("Fibrosis_region",order, LAGRANGE, &tissue_subdomains);
              fibrosis_system.add_vector("Fibrosis_Check");
            }

            system.add_variable("V", order, LAGRANGE, &tissue_subdomains);
            system.add_matrix("K");
            system.add_matrix("Ki");
            system.add_matrix("Kmono");
            system.add_matrix("M");
            system.add_vector("ML");
            system.add_vector("aux1");
            system.add_vector("aux2");
            system.add_vector("R_dist");
            system.add_vector("In");
            system.add_vector("Inm1");
            system.add_vector("Inm2");
            system.add_vector("Vnm1");
            system.add_vector("Vnprev");
            system.add_vector("Vnm2");
            system.add_vector("F"); // forcing term for exact solution
            system.add_vector("aux1"); // auxilliary vector for assembling the RHS
            system.add_vector("ForcingConv");
            system.add_vector("I_extra_stim");
            system.add_vector("DiffVecV");
            system.add_vector("DiffVecVe");
            system.add_vector("DiffVecVb");
            system.add_vector("Vcurr"); // auxilliary vector for assembling the RHS
            system.add_vector("Vprev"); // auxilliary vector for assembling the RHS

            if(exportVe_mono){
              system.add_vector("Beta_Im");
            }


           

            if(cellModel.compare(0,12,"Fenton-Karma") == 0){
              recovery.add_variable("v", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("w", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("s", order, LAGRANGE, &tissue_subdomains);
              recovery.add_vector("v_prev");
              recovery.add_vector("w_prev");
              recovery.add_vector("s_prev");
              recovery.add_vector("v_prev_prev");
              recovery.add_vector("w_prev_prev");
              recovery.add_vector("s_prev_prev");   
            }
            else if(cellModel.compare(0,12,"Courtemanche") == 0){
              recovery.add_variable("v", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("w", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("xs", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("u", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("m", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("h", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("d", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("xr", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Nai", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Ki", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Carel", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("oa", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("ui", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("oi", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("f", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("j", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Cai", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("Caup", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("ua", order, LAGRANGE, &tissue_subdomains);
              recovery.add_variable("fCa", order, LAGRANGE, &tissue_subdomains); 

              
              recovery.add_vector("v_prev");
              recovery.add_vector("w_prev");
              recovery.add_vector("xs_prev");
              recovery.add_vector("u_prev");
              recovery.add_vector("m_prev");
              recovery.add_vector("h_prev");
              recovery.add_vector("d_prev");
              recovery.add_vector("xr_prev");
              recovery.add_vector("Nai_prev");
              recovery.add_vector("Ki_prev");
              recovery.add_vector("Carel_prev");
              recovery.add_vector("oa_prev");
              recovery.add_vector("ui_prev");
              recovery.add_vector("oi_prev");
              recovery.add_vector("f_prev");
              recovery.add_vector("j_prev");
              recovery.add_vector("Cai_prev");
              recovery.add_vector("Caup_prev");
              recovery.add_vector("ua_prev");
              recovery.add_vector("fCa_prev");

              recovery.add_vector("v_prev_prev");
              recovery.add_vector("w_prev_prev");
              recovery.add_vector("xs_prev_prev");
              recovery.add_vector("u_prev_prev");
              recovery.add_vector("m_prev_prev");
              recovery.add_vector("h_prev_prev");
              recovery.add_vector("d_prev_prev");
              recovery.add_vector("xr_prev_prev");
              recovery.add_vector("Nai_prev_prev");
              recovery.add_vector("Ki_prev_prev");
              recovery.add_vector("Carel_prev_prev");
              recovery.add_vector("oa_prev_prev");
              recovery.add_vector("ui_prev_prev");
              recovery.add_vector("oi_prev_prev");
              recovery.add_vector("f_prev_prev");
              recovery.add_vector("j_prev_prev");
              recovery.add_vector("Cai_prev_prev");
              recovery.add_vector("Caup_prev_prev");
              recovery.add_vector("ua_prev_prev");
              recovery.add_vector("fCa_prev_prev");
              
            }

            es.init();



            break;
        }
        */

    //}

    es.print_info();
    //libMesh::out << "ALL GOOD TILL HERE 2" << std::endl;

    // setup pacing protocol
    //Pacing pacing(data);
    //pacing.print();
    //IonicModel ionic_model;
    //onic_model.setup(data);
    //ionic_model.print();

    //bool ICconditions = data("ICconditions", false);
    int ICsave_iter_start = data("ICsave_iter_start", 0);
    double ICtime = data("ICtime", 0.);
    double xDim = data("maxx",1.) - data("minx", -1.);

    //libMesh::out << "ALL GOOD TILL HERE 3" << std::endl;
    
    exported_variables.insert("parabolic");
          
    if(prescribedFibers){
        exported_variables.insert("fiberSolution");
     }
     if(exportVMAT_time){
        exported_variables.insert("LATSystem");
     }
        
    //if(cellModel.compare(0,12,"Courtemanche") == 0 || cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
    if(cellModel.compare(0,12,"Courtemanche") == 0){
       init_cd_exact(es, xDim, integrator, data);
       libMesh::out << "Setting to resting membrane voltage..." << std::endl;
    }
    if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0 && FK_in_mV){
       init_cd_exact(es, xDim, integrator, data);
       libMesh::out << "Setting to resting membrane voltage..." << std::endl;
    }

    int save_iter;
    int save_iter2;
    //libMesh::out << "ALL GOOD TILL HERE 4" << std::endl;
      if(ICconditions){

        save_iter = ICsave_iter_start;
        datatime.time = ICtime;

        std::cout << "Time: " << datatime.time << ", " << std::flush;
        std::cout << "Not exporting because of initial conditions  ... " << std::flush;
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << save_iter;
        std::string step_str = ss.str();
        

        libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es, &exported_variables);


        std::cout << "done " << std::endl;
      }
      else{
        
        save_iter = 0;
        save_iter2 = 0;

        std::cout << "Time: " << datatime.time << ", " << std::flush;
        std::cout << "Exporting  ... " << std::flush;
        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << save_iter;
        std::string step_str = ss.str();
        

        libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es, &exported_variables);


        std::cout << "done " << std::endl;
        
      }
        
    //libMesh::out << "ALL GOOD TILL HERE 5" << std::endl;

    // Start loop in time

    libMesh::out << "Outside the time is: " << datatime.time << "   with save iter being: " <<  save_iter<< std::endl;

    double totalErrorInSpaceTimeV = 0.;
    double totalErrorInSpaceTimeVe = 0.;
    double totalErrorInSpaceTimeVb = 0.;
    double L_inf_Total_space_time = 0.;
    double L_inf_V_space_time = 0.;
    double L_inf_Ve_space_time = 0.;
    double L_inf_Vb_space_time = 0.;
    

    std::vector<double> error(3);
    int numValues = data("numValues",340);
    int numElectrodeValues = data("numElectrodeValues",1);
    int numVMATValues = data("numVMATValues",1);
    //double * Ve_Cath_vector = malloc(numValues*sizeof(double));
    //double Ve_Cath_vector [numValues] = {};


    //std::vector < double > Ve_Cath_vector(numValues);

    std::vector < double > processorContribution ( mesh.n_processors() );

    std::vector < double > Ve_Cath_vector(numValues);

    std::vector < double > Ve_Cath_Mono_vector(numValues);
    std::vector < double > V_Cath_vector(numValues);

    std::vector < double > Ve_Cath_Mono_vector_curr(numElectrodeValues);
    std::vector < double > Ve_Cath_vector_curr(numElectrodeValues);

    std::vector < std::string > Ve_Cath_Mono_vector_string(numValues);
    std::vector < std::string > Ve_Cath_vector_string(numValues);

    vector<std::string> electrode_coords;
    vector<double> electrodeXvec;
    vector<double> electrodeYvec;
    vector<double> electrodeZvec;

    vector<std::string> VMAT_coords;
    vector<double> VMATXvec;
    vector<double> VMATYvec;
    vector<double> VMATZvec;
    std::vector < double > VMAT_vector_curr(numVMATValues);
    std::vector < std::string > VMAT_vector_string(numValues);

    if( numElectrodeValues > 1 ){
      std::string electrode_coords_matrix = data("electrode_coords_matrix", "fibrosisCoords.csv");
      std::string line, word;
      fstream file (electrode_coords_matrix, ios::in);

      //PULLING ELECTRODE COORDINATES INTO VECTOR
      if(file.is_open())
      {
        while(getline(file, line))
        {
          stringstream str(line);
          electrode_coords.push_back(line);

          std::vector<std::string> CoordsInfoVec;
          stringstream ss(line);
               
          while (ss.good()){
              string substr;
              getline(ss, substr, ',');
              CoordsInfoVec.push_back(substr);
          }
          electrodeXvec.push_back( std::stod(CoordsInfoVec[0]) );
          electrodeYvec.push_back( std::stod(CoordsInfoVec[1]) );
          electrodeZvec.push_back( std::stod(CoordsInfoVec[2]) );
        }
        libMesh::out << "Loading electrode points to calculate signals at: File: " << electrode_coords_matrix << " found successfully... Creating grid of electrodes..." << std::endl;
      }
      else{
        libMesh::out << "Could not open the electrodes coodinates file..." << std::endl;
      }
    }


    //std::vector < double > t_vector;
    int counterVe = 0;
    int counterVeMono = 0;
    int counterVMAT = 0;

    //es.parameters.set< std::vector < double > >("systemVectorContribution") = processorContribution;

    libMesh::out << "Cell model to use: " << cellModel << std::endl;
    int extraMeshInt = data("extraMeshInt", 1);
    int counterExport = 0;

    //libMesh::out << "ALL GOOD TILL HERE 6" << std::endl;

    double whereX = 0., whereY = 0., whereZ = 0.;

     std::string nameOfVFileExt = data("nameOfVFileExt", "apnd.txt");
     std::string nameOfImFileExt = data("nameOfImFileExt", "apnd.txt");
     std::string nameOfxFileExt = data("nameOfxFileExt", "apnd.txt");
     std::string nameOfyFileExt = data("nameOfyFileExt", "apnd.txt");
     std::string nameOfzFileExt = data("nameOfzFileExt", "apnd.txt");

     

     bool dummtrackvar = true;

     bool processorOfIntBi = false;
     std::vector<dof_id_type> dof_indicesSet;
     int varInterestSet = 0;

     //libMesh::out << "ALL GOOD TILL HERE 7" << std::endl;

     if(exportVMAT_time){
        libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
        system.get_vector("dVmdt_max_prev").zero();
        system.get_vector("dVmdt_max_check").zero();
      }

    for (; datatime.time < datatime.end_time; )
    {

        // advance time
        datatime.time += datatime.dt;
        datatime.timestep++;
        es.parameters.set< double >("phiEvalue") = 0.;
        //std::cout << "Time: " << datatime.time << std::endl;
        // advance vectors and solve
        // Bidomain System
        switch (time_integrator)
        {
            case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
            case TimeIntegrator::EXPLICIT_INTRACELLULAR:
            case TimeIntegrator::SEMI_IMPLICIT:
            case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
            {
                if(datatime.timestep == 1)
                {
                    int locRankk = init.comm().rank();
                    assemble_matrices(es, datatime, time_integrator, order, data, locRankk);
                    if(convergence_test){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICstim){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICconditions){
                      init_cd_exact(es, xDim, integrator, data);
                    }

                    es.get_system("parabolic").get_vector("R_dist").zero();
                    
                }

                libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
                //libMesh::TransientExplicitSystem & fibrosis_system = es.get_system < libMesh::TransientExplicitSystem > ("Fibrosis");

                *parabolic_system.old_local_solution = *parabolic_system.solution;
                parabolic_system.get_vector("Vcurr").zero();
                parabolic_system.get_vector("Vprev").zero();

                parabolic_system.get_vector("Vprev") = *parabolic_system.solution;

                //solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                SolverRecovery (es, data, datatime);
                // Solve Parabolic
                assemble_rhs(es, datatime, time_integrator, order, data, EquationType::PARABOLIC);
                if(time_integrator == TimeIntegrator::SEMI_IMPLICIT ||
                   time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP)
                {
                    parabolic_system.solve();
                }
                else
                {
                    *parabolic_system.solution += *parabolic_system.rhs;
                }
                parabolic_system.update();
                // Solve Elliptic
                //libMesh::out << "TRACKING BUG 0" << std::endl;
                
                //libMesh::out << "HERE" << std::endl;
                if(fibrosisBool && datatime.timestep < 5 && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                  libMesh::LinearImplicitSystem &fibrosis_system = es.get_system < libMesh::LinearImplicitSystem > ("Fibrosis");
                  *fibrosis_system.solution += parabolic_system.get_vector("Fibrosis_Check");
                }

                parabolic_system.get_vector("Vcurr") = *parabolic_system.solution;

                break;
            }
            case TimeIntegrator::SBDF3:
            {
                if(datatime.timestep == 1)
                {
                    int locRankk = init.comm().rank();
                    assemble_matrices(es, datatime, TimeIntegrator::SBDF1, order, data, locRankk);
                    if(convergence_test){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICstim){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICconditions){
                      init_cd_exact(es, xDim, integrator, data);
                    }

                    es.get_system("parabolic").get_vector("R_dist").zero();
                    
                }
                if(datatime.timestep == 2)
                {
                    int locRankk = init.comm().rank();
                    assemble_matrices(es, datatime, TimeIntegrator::SBDF2, order, data, locRankk);
                }
                if(datatime.timestep == 3)
                {
                    int locRankk = init.comm().rank();
                    assemble_matrices(es, datatime, time_integrator, order, data, locRankk);
                }

                libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");

                system.get_vector("Inm2") = system.get_vector("Inm1");
                system.get_vector("Inm1") = system.get_vector("In");
                system.get_vector("Vnm2") = *system.older_local_solution;
                system.get_vector("Vnm1") = *system.old_local_solution;
                *system.older_local_solution = *system.old_local_solution;
                *system.old_local_solution = *system.solution;
                //solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                SolverRecovery (es, data, datatime);
                if(datatime.timestep == 1) assemble_rhs(es, datatime, TimeIntegrator::SBDF1, order, data);
                else if(datatime.timestep == 2) assemble_rhs(es, datatime, TimeIntegrator::SBDF2, order, data);
                else assemble_rhs(es, datatime, time_integrator, order, data);

                system.get_vector("Vcurr").zero();
                system.get_vector("Vprev").zero();

                system.get_vector("Vprev") = *system.solution;

              
                system.solve();
                system.update();

                if(fibrosisBool && datatime.timestep < 5 && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                   libMesh::LinearImplicitSystem &fibrosis_system = es.get_system < libMesh::LinearImplicitSystem > ("Fibrosis");
                   *fibrosis_system.solution += fibrosis_system.get_vector("Fibrosis_Check");
                }
                system.get_vector("Vcurr") = *system.solution;
                

                break;
            }
            case TimeIntegrator::SBDF2:
            {
                if(datatime.timestep == 1)
                {
                    int locRankk = init.comm().rank();
                    assemble_matrices(es, datatime, TimeIntegrator::SBDF1, order, data, locRankk);
                    if(convergence_test){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICstim){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICconditions){
                      init_cd_exact(es, xDim, integrator, data);
                    }
                    es.get_system("parabolic").get_vector("R_dist").zero();
                    
                }
                if(datatime.timestep == 2)
                {
                    int locRankk = init.comm().rank();
                    assemble_matrices(es, datatime, time_integrator, order, data, locRankk);
                }

                libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");

                system.get_vector("Inm2") = system.get_vector("Inm1");
                system.get_vector("Inm1") = system.get_vector("In");
                system.get_vector("Vnm2") = *system.older_local_solution;
                system.get_vector("Vnm1") = *system.old_local_solution;
                *system.older_local_solution = *system.old_local_solution;
                *system.old_local_solution = *system.solution;
                //solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                SolverRecovery (es, data, datatime);

                

                if(datatime.timestep == 1) assemble_rhs(es, datatime, TimeIntegrator::SBDF1, order, data);
                else assemble_rhs(es, datatime, time_integrator, order, data);
                //system.rhs->print();
                //system.get_vector("ML").print();
                //system.get_matrix("M").print();
                //system.matrix->print();

                system.get_vector("Vcurr").zero();
                system.get_vector("Vprev").zero();

                system.get_vector("Vprev") = *system.solution;

                system.solve();
                //system.solution->print();
                system.update();

                if(fibrosisBool && datatime.timestep < 5 && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                  libMesh::LinearImplicitSystem &fibrosis_system = es.get_system < libMesh::LinearImplicitSystem > ("Fibrosis");
                  *fibrosis_system.solution += fibrosis_system.get_vector("Fibrosis_Check");
                }
                system.get_vector("Vcurr") = *system.solution;
                

                break;
            }
            case TimeIntegrator::SBDF1:
            default:
            {
                //libMesh::out << "ALL GOOD TILL HERE 8" << std::endl;
                if(datatime.timestep == 1){
                    int locRankk = init.comm().rank();
                    //libMesh::out << "ALL GOOD TILL HERE 9" << std::endl;
                    assemble_matrices(es, datatime, time_integrator, order, data, locRankk);
                    if(convergence_test){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICstim){
                       init_cd_exact(es, xDim, integrator, data);
                    }
                    if(ICconditions){
                      init_cd_exact(es, xDim, integrator, data);
                    }

                    es.get_system("parabolic").get_vector("R_dist").zero();
                    
                }

                libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
                system.get_vector("Inm2") = system.get_vector("Inm1");
                system.get_vector("Inm1") = system.get_vector("In");
                system.get_vector("Vnm2") = *system.older_local_solution;
                system.get_vector("Vnm1") = *system.old_local_solution;
                *system.older_local_solution = *system.old_local_solution;
                *system.old_local_solution = *system.solution;
                //solve_ionic_model_evaluate_ionic_currents(es, ionic_model, pacing, datatime, time_integrator);
                SolverRecovery (es, data, datatime);
                assemble_rhs(es, datatime, time_integrator, order, data);

                system.get_vector("Vcurr").zero();
                system.get_vector("Vprev").zero();

                system.get_vector("Vprev") = *system.solution;

                system.solve();
                system.update();

                if(fibrosisBool && datatime.timestep < 5 && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                  libMesh::LinearImplicitSystem &fibrosis_system = es.get_system < libMesh::LinearImplicitSystem > ("Fibrosis");
                  *fibrosis_system.solution += fibrosis_system.get_vector("Fibrosis_Check");
                }
                system.get_vector("Vcurr") = *system.solution;
                

                break;
            }
        }


      
        if(checkSoFar && 0 == datatime.timestep % extraMeshInt){

          if( time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3 ){
             libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
             system.get_vector("aux3").zero();
             system.get_vector("aux3") = system.get_vector("Vcurr");
             system.get_vector("aux3") -= system.get_vector("Vprev");

             double linfinityCurr = system.get_vector("aux3").linfty_norm();
             libMesh::out << "L-infinity norm for dV is " << linfinityCurr << std::endl;
          }
          else{
             libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
             system.get_vector("aux3").zero();
             system.get_vector("aux3") = system.get_vector("Vcurr");
             system.get_vector("aux3") -= system.get_vector("Vprev");

             double linfinityCurr = system.get_vector("aux3").linfty_norm();
             libMesh::out << "L-infinity norm for dV is " << linfinityCurr << std::endl;
          }
           
        }

        double exportInterval = data("exportInterval", 0.);
        //beginning of Monodomain approach for Ve calculation

        //I_m = Cm*dVm/dt + I_ion
        //phi_E(t) = (1/(4*pi*sigma_B))*integral(Beta*I_m/r)
        //EXPORTING THE SPECIFIC POINT DESIRED
        //Basically similar to the other exportVe method. We create an aux vector of Cm*dVm/dt + I_ion. Then, per processor, we
        //iterate through all the nodes and calculate that local Ve from those local nodes using the forward calculation approach
        //and then in processor 0 we add all of the individual local Ve's coming from all the processors to create and export a final Ve.
        
        if(exportVMAT_time && 0 == datatime.timestep % extraMeshInt){
        //if(exportVe_mono && 0. == ( (lround((datatime.timestep*datatime.dt*100))) % (lround(exportInterval*100)) ) ){
            counterVMAT++;
            libMesh::out << "Checking for VMAT at t = " << datatime.time << std::endl;
            double Ve, x, y, z, diffX, diffY, diffZ;
            double threshDist = data("threshDist", 0.1);
            std::vector<double> locVe, locDist;
            int nelz = data("nelz",0);


            std::vector<dof_id_type> dof_indices;
            std::vector<dof_id_type> dof_indicesP;
            std::vector<dof_id_type> dof_indicesP2;
            std::vector<dof_id_type> dof_indicesR;
            std::vector<dof_id_type> dof_indicesESmall;
            std::vector<dof_id_type> dof_indicesPSmall;

            double phi_E = 0.;

            //libMesh::out << "HERE exporting Ve" <<std::endl;
            if( time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3 ){
              //const unsigned int dim = mesh.mesh_dimension();
              const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();
              auto femSoluV = es.get_system("parabolic").variable_number("V");
              libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
              system.get_vector("dVmdt_max").zero();
              system.get_vector("aux1").zero();
              system.get_vector("aux1") = system.get_vector("Vcurr");
              system.get_vector("aux1") -= system.get_vector("Vprev");    
              system.get_vector("dVmdt_max").pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));

                    for (auto node : mesh.local_node_ptr_range())
                    {
                          dof_mapP.dof_indices(node, dof_indicesP);
                          int parabolic_ndofs = dof_indicesP.size();
                          //libMesh::out << "Here outside of if within node with elliptic size " << elliptic_ndofs << " and parabolic size of " << parabolic_ndofs << std::endl;

                          if (parabolic_ndofs > 0)
                          {
                                  double dVmdt_curr = (system.get_vector("dVmdt_max"))  (dof_indicesP[femSoluV]);
                                  double dVmdt_prev = (system.get_vector("dVmdt_max_prev"))  (dof_indicesP[femSoluV]);
                                  double dVmdt_check = (system.get_vector("dVmdt_max_check"))  (dof_indicesP[femSoluV]);
                                  if( dVmdt_curr > dVmdt_check ){
                                    system.get_vector("dVmdt_max_check").set( dof_indicesP[femSoluV], dVmdt_curr );
                                    system.get_vector("LATvec").set( dof_indicesP[femSoluV], datatime.time );
                                  }
                          }
                    }
              system.get_vector("dVmdt_max_prev") = system.get_vector("dVmdt_max");
           }
            else{
              //const unsigned int dim = mesh.mesh_dimension();
              const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();
              auto femSoluV = es.get_system("parabolic").variable_number("V");
              libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
              system.get_vector("dVmdt_max").zero();
              system.get_vector("aux1").zero();
              system.get_vector("aux1") = system.get_vector("Vcurr");
              system.get_vector("aux1") -= system.get_vector("Vprev");    
              system.get_vector("dVmdt_max").pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));

                    for (auto node : mesh.local_node_ptr_range())
                    {
                          dof_mapP.dof_indices(node, dof_indicesP);
                          int parabolic_ndofs = dof_indicesP.size();
                          if (parabolic_ndofs > 0)
                          {
                                  double dVmdt_curr = (system.get_vector("dVmdt_max"))  (dof_indicesP[femSoluV]);
                                  double dVmdt_prev = (system.get_vector("dVmdt_max_prev"))  (dof_indicesP[femSoluV]);
                                  double dVmdt_check = (system.get_vector("dVmdt_max_check"))  (dof_indicesP[femSoluV]);
                                  if( dVmdt_curr > dVmdt_check ){
                                    system.get_vector("dVmdt_max_check").set( dof_indicesP[femSoluV], dVmdt_curr );
                                    system.get_vector("LATvec").set( dof_indicesP[femSoluV], datatime.time );
                                  }
                          }
                    }
              system.get_vector("dVmdt_max_prev") = system.get_vector("dVmdt_max");              
                
            }
        }




        if( exportVMAT_time && datatime.time + datatime.dt  >= datatime.end_time ){
            libMesh::out << "Exporting Local Activation Times at t = " << datatime.time << std::endl;
            
            std::vector<dof_id_type> dof_indices;
            std::vector<dof_id_type> dof_indicesP;

            //libMesh::out << "HERE exporting Ve" <<std::endl;
            if( time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3 )
            {
              const DofMap & dof_map2 = es.get_system("LATSystem").get_dof_map();
              const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();
              //std::vector<dof_id_type> dof_indices;
              auto femSolu2 = es.get_system("LATSystem").variable_number("LAT");
              auto femSoluV = es.get_system("parabolic").variable_number("V");

              libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
              libMesh::LinearImplicitSystem &LATSystem = es.get_system < libMesh::LinearImplicitSystem > ("LATSystem");
                
                  for (auto node : mesh.local_node_ptr_range()){
                        dof_mapP.dof_indices(node, dof_indicesP);
                        dof_map2.dof_indices(node, dof_indices);
                        if( dof_indices.size() > 0 && dof_indicesP.size() > 0  ){
                                double LATvalue = (system.get_vector("LATvec"))  (dof_indicesP[femSoluV]);
                                LATSystem.solution -> set( dof_indices[femSolu2], LATvalue );
                        }
                  }  
           }
           else{
              const DofMap & dof_map2 = es.get_system("LATSystem").get_dof_map();
              const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();
              //std::vector<dof_id_type> dof_indices;
              auto femSolu2 = es.get_system("LATSystem").variable_number("LAT");
              auto femSoluV = es.get_system("parabolic").variable_number("V");

              libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
              libMesh::LinearImplicitSystem &LATSystem = es.get_system < libMesh::LinearImplicitSystem > ("LATSystem");
                
                  for (auto node : mesh.local_node_ptr_range()){
                        dof_mapP.dof_indices(node, dof_indicesP);
                        dof_map2.dof_indices(node, dof_indices);
                        if( dof_indices.size() > 0 && dof_indicesP.size() > 0 ){
                                double LATvalue = (system.get_vector("LATvec"))  (dof_indicesP[femSoluV]);
                                LATSystem.solution -> set( dof_indices[femSolu2], LATvalue );
                        }
                  }   
           }

           /*
            std::cout << "Final Time: " << datatime.time << ", " << std::flush;
            std::cout << "Exporting at final point  ... " << std::flush;
            // update output file counter
            //save_iter++;//
            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << save_iter;
            std::string step_str = ss.str();
            libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/Finalbath_" + step_str + ".pvtu", es, &exported_variables);
            std::cout << "done " << std::endl;
            */
        }



        //Export the solution if at the right timestep
        if (0 == datatime.timestep % datatime.export_timesteps)
        {
            std::cout << "Time: " << datatime.time << ", " << std::flush;
            std::cout << "Exporting  ... " << std::flush;
            // update output file counter
            save_iter++;//
            std::ostringstream ss;
            ss << std::setw(4) << std::setfill('0') << save_iter;
            std::string step_str = ss.str();
            libMesh::VTKIO(mesh).write_equation_systems(output_folder+"/bath_" + step_str + ".pvtu", es, &exported_variables);
            std::cout << "done " << std::endl;
        }

        if (0 == datatime.timestep % 3)
        {
            std::cout << "At time: " << datatime.time << " ms: " << std::flush;
            std::cout << "done " << std::endl;
        }



        if(exportVe_mono && 0 == datatime.timestep % extraMeshInt){
        //if(exportVe_mono && 0. == ( (lround((datatime.timestep*datatime.dt*100))) % (lround(exportInterval*100)) ) ){
            counterExport++;
            libMesh::out << "Exporting Ve via Monodomain approach at t = " << datatime.time << std::endl;
            double Ve, x, y, z, diffX, diffY, diffZ;
            double threshDist = data("threshDist", 0.1);
            std::vector<double> locVe, locDist;
            int nelz = data("nelz",0);

            double xE1 = data("exportVe_timeX", 0.);
            //xE1 = round_up(xE1, 8);
            double yE1 = data("exportVe_timeY", 0.);
            //yE1 = round_up(yE1, 8);
            double zE1 = data("exportVe_timeZ", 0.);
            //zE1 = round_up(zE1, 8);
            double lowBoxX = xE1 - threshDist;
            double lowBoxY = yE1 - threshDist;
            double highBoxX = xE1 + threshDist;
            double highBoxY = yE1 + threshDist;

            double lowBoxZ, highBoxZ;

            if(nelz == 0){
              lowBoxZ = 0.;
              highBoxZ = 0.;
            }
            else{
              lowBoxZ = zE1 - threshDist;
              highBoxZ = zE1 + threshDist;
            }


            //zE3 = round_up(zE3, 8);
            
            std::vector<dof_id_type> dof_indices;
            std::vector<dof_id_type> dof_indicesP;
            std::vector<dof_id_type> dof_indicesP2;
            std::vector<dof_id_type> dof_indicesR;
            std::vector<dof_id_type> dof_indicesESmall;
            std::vector<dof_id_type> dof_indicesPSmall;

            /*
            if(es.parameters.get<double>("specX") != 1000.){
              xE1 = es.parameters.get<double>("specX");
            }
            if(es.parameters.get<double>("specY") != 1000.){
              yE1 = es.parameters.get<double>("specY");
            }
            if(es.parameters.get<double>("specZ") != 1000.){
              zE1 = es.parameters.get<double>("specZ");
            }
            */

            double phi_E = 0.;
            double Vm_local = 0.;

            //libMesh::out << "HERE exporting Ve" <<std::endl;
            if( time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3 ){
              
              //const unsigned int dim = mesh.mesh_dimension();
              const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();

              auto femSoluV = es.get_system("parabolic").variable_number("V");

              libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
              
              double Cm = data("Cm", 1.);
              double sigma_b = data("sigma_b", 6.0);
              double chi = data("chi", 1e3);
              double pi = 3.14159265;
              double sigma_f_e = data("sigma_f_e", 1.5448);
              double sigma_s_e = data("sigma_s_e", 1.0438);  // sigma_t in the paper
              double sigma_n_e = data("sigma_n_e", 0.37222);

              double sigma_eTot = (sigma_f_e + sigma_s_e + sigma_n_e)/3.;

              double tissue_minz = data("tissue_minz", 0.37222);
              double tissue_maxz = data("tissue_maxz", 0.37222);

              system.get_vector("Beta_Im").zero();
              system.get_vector("aux1").zero();
              system.get_vector("aux2").zero();


              //system.get_vector("aux3").zero();

              //Im based on LHS
              //system.get_vector("aux1").add_vector(*system.solution, system.get_matrix("Kmono"));
              //system.get_vector("Beta_Im") = system.get_vector("aux1");
              //system.get_vector("Beta_Im").scale( -1.*chi / (4. * pi * sigma_b) );
              //system.get_vector("aux1").pointwise_mult(system.get_vector("Beta_Im"), system.get_vector("R_dist"));
              //system.get_vector("Beta_Im") /= system.get_vector("ML2");
              //system.get_vector("Beta_Im").pointwise_mult(system.get_vector("aux1"), system.get_vector("ML2"));
             
              //Im based on RHS
              system.get_vector("Vcurr") -= system.get_vector("Vprev");    
              system.get_vector("aux1").pointwise_mult(system.get_vector("Vcurr"), system.get_vector("ML"));
              //system.get_vector("Vcurr").scale(Cm/datatime.dt);
              system.get_vector("aux2").add_vector(system.get_vector("In"), system.get_matrix("M"));
              system.get_vector("aux1") -= system.get_vector("aux2");
              
              if( numElectrodeValues == 1 ){
                if(datatime.timestep < 500){
                    for (auto node : mesh.local_node_ptr_range())
                    {
                        dof_mapP.dof_indices(node, dof_indicesP);
                        int parabolic_ndofs = dof_indicesP.size();

                        if(dof_indicesP.size() > 0){

                                double x = (*node)(0);
                                double y = (*node)(1);
                                double z = (*node)(2);

                                double r_distLoc = ( sqrt( (xE1 - x)*(xE1 - x) + (yE1 - y)*(yE1 - y) + (zE1 - z)*(zE1 - z) ) );
                                double BImloc = (system.get_vector("aux1"))  (dof_indicesP[femSoluV]);//( 1./( 4.*pi*sigma_b*chi ) )*( Cm*(Vcurloc - Vprevloc)/datatime.dt - Inloc );//(system.get_vector("aux1"))  (dof_indicesP[femSoluV]);

                                //phi_E += ( chi/( 4.*pi*sigma_b ) ) * (BImloc / r_distLoc);
                                system.get_vector("R_dist").set( dof_indicesP[femSoluV], ((chi/(4.*pi*sigma_b)) / r_distLoc) );
                                //system.get_vector("R_dist").set( dof_indicesP[femSoluV], ( chi/( 4.*pi*sigma_b ) ) * (1. / r_distLoc) );
                        }    
                    }
                }
                phi_E = system.get_vector("aux1").dot(system.get_vector("R_dist"));

                Point point_D(xE1, yE1, zE1);
                Vm_local = system.point_value(femSoluV, point_D);

              }
              else{
                for(int i = 0; i < numElectrodeValues; i++){
                    xE1 = electrodeXvec[i];
                    yE1 = electrodeYvec[i];
                    zE1 = electrodeZvec[i];

                    for (auto node : mesh.local_node_ptr_range())
                    {
                          dof_mapP.dof_indices(node, dof_indicesP);
                          int parabolic_ndofs = dof_indicesP.size();
                          //libMesh::out << "Here outside of if within node with elliptic size " << elliptic_ndofs << " and parabolic size of " << parabolic_ndofs << std::endl;

                          if (parabolic_ndofs > 0)
                          {
                                  double x = (*node)(0);
                                  double y = (*node)(1);
                                  double z = (*node)(2);

                                  double r_distLoc = ( sqrt( (xE1 - x)*(xE1 - x) + (yE1 - y)*(yE1 - y) + (zE1 - z)*(zE1 - z) ) );
                                  double BImloc = (system.get_vector("aux1"))  (dof_indicesP[femSoluV]);//( 1./( 4.*pi*sigma_b*chi ) )*( Cm*(Vcurloc - Vprevloc)/datatime.dt - Inloc );//(system.get_vector("aux1"))  (dof_indicesP[femSoluV]);
                                  phi_E += (chi/(4.*pi*sigma_b)) * (BImloc / r_distLoc);  
                          }
                    }

                    // Gather all partial averages down to the root process
                    double *vecPhiEval = NULL;

                    if (mesh.processor_id() == 0) {
                      int numprocessors = mesh.n_processors();
                      vecPhiEval = (double *) malloc(numprocessors*sizeof(double)) ;
                    }

                    MPI_Gather(&phi_E, 1, MPI_DOUBLE, vecPhiEval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                    if(mesh.processor_id() == 0){
                        double localvaluePhiE = 0.;
                        for(int i = 0; i < mesh.n_processors() ; i++){
                          localvaluePhiE += vecPhiEval[i];
                        }
                      Ve_Cath_Mono_vector_curr[i] = localvaluePhiE;// / mesh.n_processors();//std::accumulate(processorContribution.begin(), processorContribution.end(), 0.0);
                      //numElectrodeLocal++;
                    }
                    phi_E = 0.;
                }
              }




           }
            else{
              //const unsigned int dim = mesh.mesh_dimension();
              const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();

              auto femSoluV = es.get_system("parabolic").variable_number("V");

              libMesh::TransientLinearImplicitSystem &system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
              
              double Cm = data("Cm", 1.);
              double sigma_b = data("sigma_b", 6.0);
              double chi = data("chi", 1e3);
              double pi = 3.14159265;
              double sigma_f_e = data("sigma_f_e", 1.5448);
              double sigma_s_e = data("sigma_s_e", 1.0438);  // sigma_t in the paper
              double sigma_n_e = data("sigma_n_e", 0.37222);

              double sigma_eTot = (sigma_f_e + sigma_s_e + sigma_n_e)/3.;

              double tissue_minz = data("tissue_minz", 0.37222);
              double tissue_maxz = data("tissue_maxz", 0.37222);

              system.get_vector("Beta_Im").zero();
              system.get_vector("aux1").zero();
              system.get_vector("aux2").zero();
              system.get_vector("aux3").zero();

              //Im based on LHS
              //system.get_vector("aux1").add_vector(*system.solution, system.get_matrix("Kmono"));
              //system.get_vector("Beta_Im") = system.get_vector("aux1");
              //system.get_vector("Beta_Im").scale( -1.*chi / (4. * pi * sigma_b) );
              //system.get_vector("aux1").pointwise_mult(system.get_vector("Beta_Im"), system.get_vector("R_dist"));
              //system.get_vector("Beta_Im") /= system.get_vector("ML2");
              //system.get_vector("Beta_Im").pointwise_mult(system.get_vector("aux1"), system.get_vector("ML2"));
             
              //Im based on RHS
              system.get_vector("Vcurr") -= system.get_vector("Vprev");    
              system.get_vector("aux1").pointwise_mult(system.get_vector("Vcurr"), system.get_vector("ML"));
              //system.get_vector("Vcurr").scale(Cm/datatime.dt);
              system.get_vector("aux2").add_vector(system.get_vector("In"), system.get_matrix("M"));
              system.get_vector("aux1") -= system.get_vector("aux2");
              
              if( numElectrodeValues == 1 ){
                if(datatime.timestep < 500){
                    for (auto node : mesh.local_node_ptr_range())
                    {
                        dof_mapP.dof_indices(node, dof_indicesP);
                        int parabolic_ndofs = dof_indicesP.size();

                        if(dof_indicesP.size() > 0){

                                double x = (*node)(0);
                                double y = (*node)(1);
                                double z = (*node)(2);

                                double r_distLoc = ( sqrt( (xE1 - x)*(xE1 - x) + (yE1 - y)*(yE1 - y) + (zE1 - z)*(zE1 - z) ) );
                                double BImloc = (system.get_vector("aux1"))  (dof_indicesP[femSoluV]);//( 1./( 4.*pi*sigma_b*chi ) )*( Cm*(Vcurloc - Vprevloc)/datatime.dt - Inloc );//(system.get_vector("aux1"))  (dof_indicesP[femSoluV]);

                                //phi_E += ( chi/( 4.*pi*sigma_b ) ) * (BImloc / r_distLoc);
                                //phi_E += (BImloc / r_distLoc);
                                system.get_vector("R_dist").set( dof_indicesP[femSoluV], ((chi/(4.*pi*sigma_b)) / r_distLoc) );
                        }
                        
                    }
                }
                phi_E = system.get_vector("aux1").dot(system.get_vector("R_dist"));

                Point point_D(xE1, yE1, zE1);
                Vm_local = system.point_value(femSoluV, point_D);

              }
              else{
                for(int i = 0; i < numElectrodeValues; i++){
                    xE1 = electrodeXvec[i];
                    yE1 = electrodeYvec[i];
                    zE1 = electrodeZvec[i];

                    for (auto node : mesh.local_node_ptr_range())
                    {
                          dof_mapP.dof_indices(node, dof_indicesP);
                          int parabolic_ndofs = dof_indicesP.size();
                          //libMesh::out << "Here outside of if within node with elliptic size " << elliptic_ndofs << " and parabolic size of " << parabolic_ndofs << std::endl;

                          if (parabolic_ndofs > 0)
                          {
                                  double x = (*node)(0);
                                  double y = (*node)(1);
                                  double z = (*node)(2);

                                  double r_distLoc = ( sqrt( (xE1 - x)*(xE1 - x) + (yE1 - y)*(yE1 - y) + (zE1 - z)*(zE1 - z) ) );
                                  double BImloc = (system.get_vector("aux1"))  (dof_indicesP[femSoluV]);//( 1./( 4.*pi*sigma_b*chi ) )*( Cm*(Vcurloc - Vprevloc)/datatime.dt - Inloc );//(system.get_vector("aux1"))  (dof_indicesP[femSoluV]);
                                  phi_E += (chi/(4.*pi*sigma_b)) * (BImloc / r_distLoc);  
                          }
                    }

                    // Gather all partial averages down to the root process
                    double *vecPhiEval = NULL;

                    if (mesh.processor_id() == 0) {
                      int numprocessors = mesh.n_processors();
                      vecPhiEval = (double *) malloc(numprocessors*sizeof(double)) ;
                    }

                    MPI_Gather(&phi_E, 1, MPI_DOUBLE, vecPhiEval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

                    if(mesh.processor_id() == 0){
                        double localvaluePhiE = 0.;
                        for(int i = 0; i < mesh.n_processors() ; i++){
                          localvaluePhiE += vecPhiEval[i];
                        }
                      Ve_Cath_Mono_vector_curr[i] = localvaluePhiE;// / mesh.n_processors();//std::accumulate(processorContribution.begin(), processorContribution.end(), 0.0);
                      //numElectrodeLocal++;
                    }
                    phi_E = 0.;
                }
              }
                
            }

            // Gather all partial averages down to the root process
            if( numElectrodeValues == 1 ){
              double *vecPhiEval = NULL;
              if (mesh.processor_id() == 0) {
                int numprocessors = mesh.n_processors();
                vecPhiEval = (double *) malloc(numprocessors*sizeof(double)) ;
              }

              MPI_Gather(&phi_E, 1, MPI_DOUBLE, vecPhiEval, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

              if(mesh.processor_id() == 0){
                int numNonzero = 0;
                //if(mesh.processor_id() == 0){
                  double localvaluePhiE = 0.;
                  for(int i = 0; i < mesh.n_processors() ; i++){
                    localvaluePhiE += vecPhiEval[i];
                  }
                Ve_Cath_Mono_vector[counterVeMono] = phi_E;//std::accumulate(processorContribution.begin(), processorContribution.end(), 0.0);
                V_Cath_vector[counterVeMono] = Vm_local;
                counterVeMono++;
              }
            }
            else{

              if(mesh.processor_id() == 0){
                std::string currentTotalPhiE;
                
                for(int i = 0; i < Ve_Cath_Mono_vector_curr.size(); i++){
                    if( i != Ve_Cath_Mono_vector_curr.size()-1 ){
                        currentTotalPhiE += to_string( Ve_Cath_Mono_vector_curr[i] ) + ",";
                    }
                    else{
                        currentTotalPhiE += to_string( Ve_Cath_Mono_vector_curr[i] );
                    }
                }

                Ve_Cath_Mono_vector_string[counterVeMono] = currentTotalPhiE;//std::accumulate(processorContribution.begin(), processorContribution.end(), 0.0);
                counterVeMono++;
              }

            }

            
            //libMesh::out << "The forward-calculated extracellular signal value is " << phi_E << " mV" << std::endl;

        }

        //end of Monodomain approach for Ve calculation


        //if(counterExport >= 1){
          //dummtrackvar = false;
        //}
        
    /*
     //START OF ERROR CALCULATION
      if(convergence_test){

          double totalErrorInSpaceV = 0.;
          double totalErrorInSpaceVe = 0.;
          double totalErrorInSpaceVb = 0.;
          double L_inf_V_space = 0.;
          double L_inf_Ve_space = 0.;
          double L_inf_Vb_space = 0.;
          double L_inf_Total_space = 0.;
          double dummyVar = 0.;
          double u_h, u_exact, ue_h, ue_exact, ub_h, ub_exact;
          const unsigned int dim = mesh.mesh_dimension();
          const DofMap & dof_map2 = es.get_system("Recovery").get_dof_map();
          std::vector<dof_id_type> dof_indices;
          std::vector<dof_id_type> dof_indices2;
          std::vector<dof_id_type> dof_indices3;



          if( time_integrator == TimeIntegrator::SBDF1 || time_integrator == TimeIntegrator::SBDF2 || time_integrator == TimeIntegrator::SBDF3 )
          {
              auto femSolu = es.get_system("bidomain").variable_number("V");
              FEType fe_type = dof_map2.variable_type(femSolu);
              std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
              QGauss qrule (dim, FIFTH);
              // Tell the finite element object to use our quadrature rule.
              fe->attach_quadrature_rule (&qrule);
              const std::vector<Point> & q_point = fe->get_xyz();
              const std::vector<Real> & JxW = fe->get_JxW();
              const std::vector<std::vector<Real>> & phi = fe->get_phi();

              for(const auto & elem : mesh.active_local_element_ptr_range()){

              fe->reinit (elem);

              for (unsigned int qp=0; qp<qrule.n_points(); qp++){

                    if(elem -> subdomain_id() == 1){
                      u_h += es.get_system("bidomain").point_value(femSolu, q_point[qp], elem);
                      u_exact += exact_solutionV_all(q_point[qp](0), q_point[qp](1), 2, datatime.time, 0.0, xDim);
                      ue_h += es.get_system("bidomain").point_value(femSolu2, q_point[qp], elem);
                      ue_exact += exact_solutionV_all(q_point[qp](0),q_point[qp](1), 3, datatime.time, 0.0, xDim);
                      ub_h += 0.;
                      ub_exact += 0.;

                      //equation_systems.get_system("Bidomain").point_value(femSoluV_exact, q_point[qp], elem) = u_exact;
                      //equation_systems.get_system("Bidomain").point_value(femSoluVe_exact, q_point[qp], elem) = ue_exact;

                    }
                    else{
                      u_h += 0.;
                      u_exact += 0.;
                      ue_h += 0.;
                      ue_exact += 0.;
                      ub_h += es.get_system("bidomain").point_value(femSolu2, q_point[qp], elem);
                      ub_exact += exact_solutionV_all(q_point[qp](0),q_point[qp](1), 4, datatime.time, 0.0, xDim);

                      //equation_systems.get_system("Bidomain").point_value(femSoluVe_exact, q_point[qp], elem) = ue_exact;
                    }


                    totalErrorInSpaceV += JxW[qp]*((u_h - u_exact)*(u_h - u_exact));
                    totalErrorInSpaceVe += JxW[qp]*((ue_h - ue_exact)*(ue_h - ue_exact));
                    totalErrorInSpaceVb += JxW[qp]*((ub_h - ub_exact)*(ub_h - ub_exact));


                    dummyVar = 1.*std::abs(u_h - u_exact);
                    //libMesh::out << "Difference: " << (u_h - u_exact)<< std::endl;

                    if(dummyVar > L_inf_V_space){
                     L_inf_V_space = dummyVar;
                    }

                    dummyVar = 1.*std::abs(ue_h - ue_exact);

                    if(dummyVar > L_inf_Ve_space){
                      L_inf_Ve_space = dummyVar;
                    }

                    dummyVar = 1.*std::abs(ub_h - ub_exact);

                    if(dummyVar > L_inf_Vb_space){
                      L_inf_Vb_space = dummyVar;
                    }
                    //rowsn = double(elem);
                    //colsn = qp;
                    //cout << qp << '\n';

                    //if(qp == 0){
                      //rowsn = rowsn + 1;
                      //cout << 'this is rows: ' << rowsn << '\n';
                    //}

                    //cout << 'element: ' << elem << ' with qp: ' << q_point[qp] << endl;
                    //cout << typeof(elem).name() << endl;
                  }

              }
           }
          else{
              auto femSolu = es.get_system("parabolic").variable_number("V");
              auto femSolu2 = es.get_system("elliptic").variable_number("Ve");
              FEType fe_type = dof_map2.variable_type(femSolu);
              std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
              QGauss qrule (dim, FIFTH);
              // Tell the finite element object to use our quadrature rule.
              fe->attach_quadrature_rule (&qrule);
              const std::vector<Point> & q_point = fe->get_xyz();
              const std::vector<Real> & JxW = fe->get_JxW();
              const std::vector<std::vector<Real>> & phi = fe->get_phi();

              for(const auto & elem : mesh.active_local_element_ptr_range()){

              fe->reinit (elem);

              for (unsigned int qp=0; qp<qrule.n_points(); qp++){

                    if(elem -> subdomain_id() == 1){
                      u_h += es.get_system("parabolic").point_value(femSolu, q_point[qp], elem);
                      u_exact += exact_solutionV_all(q_point[qp](0), q_point[qp](1), 2, datatime.time, 0.0, xDim);
                      ue_h += es.get_system("elliptic").point_value(femSolu2, q_point[qp], elem);
                      ue_exact += exact_solutionV_all(q_point[qp](0),q_point[qp](1), 3, datatime.time, 0.0, xDim);
                      ub_h += 0.;
                      ub_exact += 0.;

                      //equation_systems.get_system("Bidomain").point_value(femSoluV_exact, q_point[qp], elem) = u_exact;
                      //equation_systems.get_system("Bidomain").point_value(femSoluVe_exact, q_point[qp], elem) = ue_exact;

                    }
                    else{
                      u_h += 0.;
                      u_exact += 0.;
                      ue_h += 0.;
                      ue_exact += 0.;
                      ub_h += es.get_system("elliptic").point_value(femSolu2, q_point[qp], elem);
                      ub_exact += exact_solutionV_all(q_point[qp](0),q_point[qp](1), 4, datatime.time, 0.0, xDim);

                      //equation_systems.get_system("Bidomain").point_value(femSoluVe_exact, q_point[qp], elem) = ue_exact;
                    }

                    totalErrorInSpaceV += JxW[qp]*((u_h - u_exact)*(u_h - u_exact));
                    totalErrorInSpaceVe += JxW[qp]*((ue_h - ue_exact)*(ue_h - ue_exact));
                    totalErrorInSpaceVb += JxW[qp]*((ub_h - ub_exact)*(ub_h - ub_exact));



                    dummyVar = 1.*std::abs(u_h - u_exact);
                    //libMesh::out << "Difference: " << (u_h - u_exact)<< std::endl;

                    if(dummyVar > L_inf_V_space){
                     L_inf_V_space = dummyVar;
                    }

                    dummyVar = 1.*std::abs(ue_h - ue_exact);

                    if(dummyVar > L_inf_Ve_space){
                      L_inf_Ve_space = dummyVar;
                    }

                    dummyVar = 1.*std::abs(ub_h - ub_exact);

                    if(dummyVar > L_inf_Vb_space){
                      L_inf_Vb_space = dummyVar;
                    }
                    //rowsn = double(elem);
                    //colsn = qp;
                    //cout << qp << '\n';

                    //if(qp == 0){
                      //rowsn = rowsn + 1;
                      //cout << 'this is rows: ' << rowsn << '\n';
                    //}

                    //cout << 'element: ' << elem << ' with qp: ' << q_point[qp] << endl;
                    //cout << typeof(elem).name() << endl;
                  }

              }
           }



      totalErrorInSpaceTimeV += datatime.dt*std::sqrt(totalErrorInSpaceV);
      totalErrorInSpaceTimeVe += datatime.dt*std::sqrt(totalErrorInSpaceVe);
      totalErrorInSpaceTimeVb += datatime.dt*std::sqrt(totalErrorInSpaceVb);

      L_inf_V_space_time += datatime.dt*L_inf_V_space;
      L_inf_Ve_space_time += datatime.dt*L_inf_Ve_space;
      L_inf_Vb_space_time += datatime.dt*L_inf_Vb_space;

      libMesh::out << "The L2 norm for V in space is: " << std::sqrt(totalErrorInSpaceV)<< std::endl;
      libMesh::out << "The L2 norm for Ve in space is: " << std::sqrt(totalErrorInSpaceVe)<< std::endl;
      libMesh::out << "The L2 norm for Vb in space is: " << std::sqrt(totalErrorInSpaceVb)<< std::endl;

      libMesh::out << "The L_infinity norm for the V in space is: " << L_inf_V_space << std::endl;
      libMesh::out << "The L_infinity norm for the Ve in space is: " << L_inf_Ve_space << std::endl;
      libMesh::out << "The L_infinity norm for the Vb in space is: " << L_inf_Vb_space << std::endl;

      //libMesh::out << "The CUMULATIVE L2 norm for V in space and time is: " << (totalErrorInSpaceTimeV)<< std::endl;

      //START OF SAVING SOLUTION SYSTEM
        if(convergence_test_save){
          auto femSoluV_exact = es.get_system("Solution").variable_number("V_exact");
          auto femSoluVe_exact = es.get_system("Solution").variable_number("Ve_exact");

          const DofMap & dof_map = es.get_system("Solution").get_dof_map();

          for(const auto & node : mesh.local_node_ptr_range()){

            dof_map.dof_indices (node, dof_indices);
            dof_map2.dof_indices (node, dof_indices2);

            const Real x = (*node)(0);
            const Real y = (*node)(1);

            if(dof_indices2.size() > 0){

            u_exact = exact_solutionV_all(x, y, 2, datatime.time, 0.0, xDim);
            es.get_system("Solution").solution -> set(dof_indices[femSoluV_exact],u_exact);
            ue_exact = exact_solutionV_all(x, y, 3, datatime.time, 0.0, xDim);
            es.get_system("Solution").solution -> set(dof_indices[femSoluVe_exact],ue_exact);

            }
            else{
            ue_exact = exact_solutionV_all(x, y, 4, datatime.time, 0.0, xDim);
            es.get_system("Solution").solution -> set(dof_indices[femSoluVe_exact],ue_exact);
            }

          }
        }
        //END OF SAVING SOLUTION SYSTEM

    }
    */
    //END OF ERROR CALCULATION

    }
    //END OF TIME LOOP

    auto end = sc.now();       // end timer (starting & ending is done by measuring the time at the moment the process started & ended respectively)
    auto time_span = static_cast<chrono::duration<double>>(end - start);   // measure time span between start & end
    cout<<"Operation took: "<<time_span.count()<<" seconds"<<std::endl;



    //bool ICsave_state = data("ICsave_state", false);


    /*
    if(ICsave_state){


        //const unsigned int dim = mesh.mesh_dimension();
        //libMesh::out << "DEBUGGING 0" << std::endl;
        const DofMap & dof_mapR = es.get_system("Recovery").get_dof_map();
        const DofMap & dof_mapE = es.get_system("elliptic").get_dof_map();
        const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();
        std::vector<dof_id_type> dof_indicesP;
        std::vector<dof_id_type> dof_indicesE;
        std::vector<dof_id_type> dof_indicesR;

        //libMesh::out << "DEBUGGING 1" << std::endl;
        auto femSoluV = es.get_system("parabolic").variable_number("V");
        auto femSoluVe = es.get_system("elliptic").variable_number("Ve");
        auto femSoluv = es.get_system("Recovery").variable_number("v");
        auto femSoluw = es.get_system("Recovery").variable_number("w");
        auto femSolus = es.get_system("Recovery").variable_number("s");

        //libMesh::out << "DEBUGGING 2" << std::endl;
        std::vector<int> IC_Pindices_vec;
        std::vector<int> IC_Eindices_vec;
        std::vector<double> IC_V_vec;
        std::vector<double> IC_Ve_vec;
        std::vector<double> IC_v_vec;
        std::vector<double> IC_w_vec;
        std::vector<double> IC_s_vec;

        libMesh::LinearImplicitSystem &elliptic_system = es.get_system < libMesh::LinearImplicitSystem > ("elliptic");
        libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
        libMesh::TransientExplicitSystem & recovery_system = es.get_system < libMesh::TransientExplicitSystem > ("Recovery");

        
        //libMesh::out << "DEBUGGING 3" << std::endl;
        for(const auto & node : mesh.local_node_ptr_range()){

          dof_mapP.dof_indices (node, dof_indicesP);
          dof_mapR.dof_indices (node, dof_indicesR);
          dof_mapE.dof_indices (node, dof_indicesE);

          //libMesh::out << "DEBUGGING 4" << std::endl;

          const int nodeid = node->id();

          //libMesh::out << "DEBUGGING 5" << std::endl;
          if(dof_indicesP.size() == dof_indicesE.size()){

            //libMesh::out << "DEBUGGING 6 P" << std::endl;
            IC_V_vec.push_back( (*parabolic_system.solution)(dof_indicesP[femSoluV]) );
            IC_Ve_vec.push_back( (*elliptic_system.solution)(dof_indicesE[femSoluVe]) );
            IC_Pindices_vec.push_back( nodeid );
            IC_Eindices_vec.push_back( nodeid );
            IC_v_vec.push_back( (*recovery_system.solution)(dof_indicesR[femSoluv]) );
            IC_w_vec.push_back( (*recovery_system.solution)(dof_indicesR[femSoluw]) );
            IC_s_vec.push_back( (*recovery_system.solution)(dof_indicesR[femSolus]) );
            //libMesh::out << "DEBUGGING 7 P" << std::endl;

          }
          else{
            //libMesh::out << "DEBUGGING 6 E" << std::endl;
            IC_Ve_vec.push_back( (*elliptic_system.solution)(dof_indicesE[femSoluVe]) );
            IC_Eindices_vec.push_back( nodeid );
            //libMesh::out << "DEBUGGING 7 E" << std::endl;

          }

        }

        //Export here
        //libMesh::out << "DEBUGGING 8" << std::endl;
        std::string nameOfPindicesFile, nameOfEindicesFile, nameOfVFile, nameOfVeFile, nameOfvFile, nameOfwFile, nameOfsFile;
        nameOfPindicesFile = "IC_Pindices_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        nameOfEindicesFile = "IC_Eindices_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        nameOfVFile = "IC_V_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        nameOfVeFile = "IC_Ve_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        nameOfvFile = "IC_v_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        nameOfwFile = "IC_w_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        nameOfsFile = "IC_s_"+integrator+"_"+std::to_string(datatime.time - datatime.dt)+"ms.txt";
        //nameOfTimeFile = "t_vector_"+integrator+"_"+exportVe_Name+"_"+std::to_string(counterVe)+"points.txt";
        //outputFile2.open (nameOfOutputFile2);
        //outputFile2 << (Ve_Cath_vector)<< std::endl;
        //outputFile2.close();
        //libMesh::out << "DEBUGGING 9" << std::endl;
        ofstream output_filePI(nameOfPindicesFile);
        std::ostream_iterator<double> output_iteratorPI(output_filePI, "\n");
        std::copy(IC_Pindices_vec.begin(), IC_Pindices_vec.end(), output_iteratorPI);

        ofstream output_fileEI(nameOfEindicesFile);
        std::ostream_iterator<double> output_iteratorEI(output_fileEI, "\n");
        std::copy(IC_Eindices_vec.begin(), IC_Eindices_vec.end(), output_iteratorEI);

        ofstream output_fileV(nameOfVFile);
        std::ostream_iterator<double> output_iteratorV(output_fileV, "\n");
        std::copy(IC_V_vec.begin(), IC_V_vec.end(), output_iteratorV);

        ofstream output_fileVe(nameOfVeFile);
        std::ostream_iterator<double> output_iteratorVe(output_fileVe, "\n");
        std::copy(IC_Ve_vec.begin(), IC_Ve_vec.end(), output_iteratorVe);

        ofstream output_filev(nameOfvFile);
        std::ostream_iterator<double> output_iteratorv(output_filev, "\n");
        std::copy(IC_v_vec.begin(), IC_v_vec.end(), output_iteratorv);

        ofstream output_filew(nameOfwFile);
        std::ostream_iterator<double> output_iteratorw(output_filew, "\n");
        std::copy(IC_w_vec.begin(), IC_w_vec.end(), output_iteratorw);

        ofstream output_files(nameOfsFile);
        std::ostream_iterator<double> output_iterators(output_files, "\n");
        std::copy(IC_s_vec.begin(), IC_s_vec.end(), output_iterators);


    }

    */

    /*
    finputV.close();
    foutputV.close(); 

    finputIm.close();
    foutputIm.close(); 
    */
    

    if(exportVe_mono){

      if(mesh.processor_id() == 0){

        if( numElectrodeValues == 1 ){
          //ofstream outputFile2;
          std::string nel2 = data("nelx", "50");
          std::string exportVe_Name = data("exportVe_Name","ExtracellularNegative_NoFibrosis");
          std::string nameOfOutputFile2, nameOfOutputFileCoords, nameOfOutputFileV1;
          std::string nameOfTimeFile;
          nameOfOutputFile2 = "Ve_Catheter_MonoAp_"+integrator+"_"+exportVe_Name+"_"+std::to_string(counterVeMono)+"points.txt";
          nameOfOutputFileCoords = "Ve_Catheter_Coords_MonoAp_"+integrator+"_"+exportVe_Name+"_"+std::to_string(counterVeMono)+"points.txt";

          ofstream output_file2(nameOfOutputFile2);
          std::ostream_iterator<double> output_iterator(output_file2, "\n");
          std::copy(Ve_Cath_Mono_vector.begin(), Ve_Cath_Mono_vector.end(), output_iterator);

          nameOfOutputFileV1 = "V_Catheter_"+integrator+"_"+exportVe_Name+"_"+std::to_string(counterVe)+"points.csv";
          ofstream output_fileV1(nameOfOutputFileV1);
          std::ostream_iterator<double> output_iteratorV1(output_fileV1, "\n");
          std::copy(V_Cath_vector.begin(), V_Cath_vector.end(), output_iteratorV1);

          //ofstream output_file3(nameOfOutputFileCoords);
          //output_file3 << es.parameters.get<double>("specX") << " , " << es.parameters.get<double>("specY") << " , " << es.parameters.get<double>("specZ") << std::endl;
          //output_file3.close();
        }
        else{
          std::string nel2 = data("nelx", "50");
          std::string exportVe_Name = data("exportVe_Name","ExtracellularNegative_NoFibrosis");
          std::string nameOfOutputFile2;
          std::string nameOfTimeFile;
          nameOfOutputFile2 = "Ve_Catheter_MonoAp_"+integrator+"_"+exportVe_Name+"_"+std::to_string(counterVeMono)+"points.csv";
          
          ofstream output_file2(nameOfOutputFile2);
          std::ostream_iterator< std::string > output_iterator(output_file2, "\n");
          std::copy(Ve_Cath_Mono_vector_string.begin(), Ve_Cath_Mono_vector_string.end(), output_iterator);
        }

      }

    }
    /*
    if(convergence_test)
    {
        libMesh::out << "The L2 norm for V in space AND time is: " << (totalErrorInSpaceTimeV)<< std::endl;
        libMesh::out << "The L2 norm for Ve in space AND time is: " << (totalErrorInSpaceTimeVe)<< std::endl;
        libMesh::out << "The L2 norm for Vb in space AND time is: " << (totalErrorInSpaceTimeVb)<< std::endl;

        libMesh::out << "The L_infinity norm for V in space AND time is: " << (L_inf_V_space_time)<< std::endl;
        libMesh::out << "The L_infinity norm for Ve in space AND time is: " << (L_inf_Ve_space_time)<< std::endl;
        libMesh::out << "The L_infinity norm for Vb in space AND time is: " << (L_inf_Vb_space_time)<< std::endl;

        ofstream outputFile;
        std::string nel2 = data("nelx", "50");
        std::string eltype2 = data("eltype", "QUAD4");

        std::string eltype = data("eltype", "simplex");
        std::string nameOfOutputFile = "outPut_for_"+integrator+"_"+nel2+"elems_"+eltype2+"type.txt";
        outputFile.open (nameOfOutputFile);
        outputFile << "The L2 norm for V in space AND time is: " << (totalErrorInSpaceTimeV)<< std::endl;
        outputFile << "The L2 norm for Ve in space AND time is: " << (totalErrorInSpaceTimeVe)<< std::endl;
        outputFile << "The L2 norm for Vb in space AND time is: " << (totalErrorInSpaceTimeVb)<< std::endl;
        outputFile << "The L_infinity norm for V in space AND time is: " << (L_inf_V_space_time)<< std::endl;
        outputFile << "The L_infinity norm for Ve in space AND time is: " << (L_inf_Ve_space_time)<< std::endl;
        outputFile << "The L_infinity norm for Vb in space AND time is: " << (L_inf_Vb_space_time)<< std::endl;
        outputFile << "Operation took: "<<time_span.count()<<" seconds"<<std::endl;
        outputFile.close();
    }
    */

    //#endif

    return 0;
}



void read_mesh(const GetPot &data, libMesh::ParallelMesh &mesh, libMesh::SerialMesh &meshRef)
{
    using namespace libMesh;
    // read mesh from file?
    std::string mesh_filename = data("mesh_filename", "NONE");
    if (mesh_filename.compare("NONE") != 0)
    {
        // READ MESH
        std::cout << "I will read the mesh here: " << mesh_filename << std::endl;
        mesh.read(mesh_filename);
        //meshRef.read("block11.e");
        //libMesh::ExodusII_IO_Helper ehh(mesh);
        std::cout << "I read the mesh " << mesh_filename << " and I am about to prepare it for use." << std::endl;

        //MESH BLOCKS
       // mesh.get_subdomain_name_map();
        //libMesh::out << "The number of blocks is " <<  << std::endl;

        //MESH REFINEMENT
        int mesh_refinement_factor_tissue = data("mesh_refinement_factor_tissue", 0);
        int mesh_refinement_factor_bath = data("mesh_refinement_factor_bath", 0);
        double mesh_refinement_desired_h = data("mesh_refinement_desired_h", 0.03); //cm
        bool uniformRefinement = data("uniformRefinement", false);

        if(mesh_refinement_factor_tissue > 0){
          MeshRefinement meshR(mesh);
          //MeshRefinement meshR2(meshRef);
          libMesh::out << "HERE before mesh refinement" << std::endl;
          
          for(int iMesh = 0; iMesh < mesh_refinement_factor_tissue; iMesh++){

            if(uniformRefinement){
              meshR.uniformly_refine();
              //meshR2.uniformly_refine();
              mesh.prepare_for_use(false);
              //meshRef.prepare_for_use(false);
            }
            else{
              for(auto el: mesh.active_local_element_ptr_range()){

                double distanceElem = ( (el -> hmax()) + (el -> hmin()) ) /2.;
                if( distanceElem >= mesh_refinement_desired_h*2. && el -> hmin() > mesh_refinement_desired_h/3. ){
                  el -> set_refinement_flag(libMesh::Elem::REFINE);
                }
                else{
                  el -> set_refinement_flag(libMesh::Elem::DO_NOTHING);
                }
              }
              meshR.refine_elements();
              mesh.prepare_for_use(false);
            }
          }
          libMesh::out << "HERE after mesh refinement" << std::endl;
        }


        //IMPORTING FIBROSIS INFO
        bool fibrosisBool = data("fibrosisBool", false);
        std::string fibrosisMethod = data("fibrosisMethod","Conductivities");
        double thresholdPerElement = data("thresholdPerElement", .5);

        vector<std::string> fibrosis_coords2;
        double threshDist = data("threshDist", .1);
        vector<double> fibrosisNodeXvec2;
        vector<double> fibrosisNodeYvec2;
        vector<double> fibrosisNodeZvec2;

        if(fibrosisBool && fibrosisMethod.compare(0,11,"Percolation") == 0){
          std::string fibrosis_coords_matrix = data("fibrosis_coords_matrix", "fibrosisCoords.csv");
          std::string line, word;
          fstream file (fibrosis_coords_matrix, ios::in);

          //PULLING FIBROSIS COORDINATES INTO VECTOR
          if(file.is_open())
          {
            while(getline(file, line))
            {
              stringstream str(line);
              fibrosis_coords2.push_back(line);

              std::vector<std::string> CoordsInfoVec;
              stringstream ss(line);
                   
              while (ss.good()){
                  string substr;
                  getline(ss, substr, ',');
                  CoordsInfoVec.push_back(substr);
              }
              fibrosisNodeXvec2.push_back( std::stod(CoordsInfoVec[0]) );
              fibrosisNodeYvec2.push_back( std::stod(CoordsInfoVec[1]) );
              fibrosisNodeZvec2.push_back( std::stod(CoordsInfoVec[2]) );
            }
            libMesh::out << "Fibrosis case: File: " << fibrosis_coords_matrix << " found successfully... Eliminating elements from tissue domain..." << std::endl;
            libMesh::out << "Fibrosis method used: " << fibrosisMethod << std::endl;
          }
          else{
            libMesh::out << "Could not open the file" << std::endl;
          }
        }


        //SUBDOMAIN DEFINITION
        //for(auto elRef: meshRef.element_ptr_range()){
            //auto c = elRef->centroid();
        int numFibrosisPoints = 0, numTissuePoints = 0;
            for(auto el: mesh.element_ptr_range()){
                if(fibrosisBool && fibrosisMethod.compare(0,11,"Percolation") == 0){
                  //check for fibrosis
                  auto localNodes = el -> get_nodes();
                  auto numNodes = el -> n_nodes();
                  std::vector<int> fibBool = {2, 2, 2, 2};
                  int max = 0, indToUse = 0, mostvalue = 2;
                
                  for(int locnode = 0; locnode < numNodes; locnode++){

                    double x = (*localNodes[locnode])(0);
                    double y = (*localNodes[locnode])(1);
                    double z = (*localNodes[locnode])(2);

                    //libMesh::out << std::fixed;
                    //libMesh::out << std::setprecision(4);

                    std::string check_Presence;

                    check_Presence = to_string(round_up(x,4))+","+to_string(round_up(y,4))+","+to_string(round_up(z,4));

                    auto icheck = std::find(fibrosis_coords2.begin(), fibrosis_coords2.end(), check_Presence);
                    //libMesh::out << check_Presence << std::endl;
                    //libMesh::out << x << "," << y << "," << z << std::endl;
                    if(fibrosis_coords2.end() != icheck){
                      //libMesh::out << "Found  first   " << check_Presence << std::endl;
                      //libMesh::out << x << "," << y << "," << z << std::endl;
                      //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);
                      fibBool[locnode] = 1;
                    }
                    else{
                      //fibBool[locnode] = 0;
                      
                      for(int checkCoords = 0; checkCoords < fibrosisNodeXvec2.size(); checkCoords++){
                        double distX = std::abs( round_up(x,8) - fibrosisNodeXvec2[checkCoords] );
                        double distY = std::abs( round_up(y,8) - fibrosisNodeYvec2[checkCoords] );
                        double distZ = std::abs( round_up(z,8) - fibrosisNodeZvec2[checkCoords] );
                        if( distX < .005 && distY < .005 && distZ < .005 ){
                          //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);
                          fibBool[locnode] = 1;
                          //break;
                          //libMesh::out << "Found  second   " << check_Presence << std::endl;
                          //libMesh::out << x << "," << y << "," << z << std::endl;
                        }
                        else{
                          if( fibBool[locnode] == 1 ){
                            //NOTHING TO DO BECAUSE IT'S ALREADY FIBROSIS
                          }
                          else if( fibBool[locnode] == 2 ){
                            //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE); 
                            fibBool[locnode] = 0;
                          }                         
                        }
                      }
                      
                    }
                  }

                  double avgFib=0., numNonTwo=0.;
                  for (int i = 0;i < fibBool.size(); i++){
                    if(fibBool[i] != 2){
                      avgFib += (double) fibBool[i];
                      numNonTwo = numNonTwo + 1.;
                    }
                  }
                  if( avgFib/numNonTwo >= thresholdPerElement ){
                    el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);
                    numFibrosisPoints += 1;
                  }
                  else{
                    el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
                    numTissuePoints += 1;
                  }

                }
                else{
                  el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
                }
                
            }
        //}


        libMesh::out << "Fibrosis density in percentage is: " << " with " << numFibrosisPoints << " fibrosis points and " << numTissuePoints << " tissue points..." << std::endl;
        if(mesh_refinement_factor_tissue > 0){
          mesh.prepare_for_use(true);
          std::cout << "Mesh was refined and subdomains specified..." << std::endl;
        }
        std::cout << "Mesh is prepared..." << std::endl;
    }
    else
    {
        // Create Mesh
        // number of elements in the x,y and z direction
        int nelx = data("nelx", 10);
        int nely = data("nely", 10);
        // if nelz = 0 we run in 2D
        int nelz = data("nelz", 0);

        // the cube dimensions are defined as [minx, maxx] x [miny, maxy] x [minz, maxz]
        double maxx = data("maxx", 1.0);
        double maxy = data("maxy", 1.0);
        double maxz = data("maxz", 1.0);
        double minx = data("minx", 0.0);
        double miny = data("miny", 0.0);
        double minz = data("minz", 0.0);

        //IMPORTING FIBROSIS INFO
        bool fibrosisBool = data("fibrosisBool", false);
        std::string fibrosisMethod = data("fibrosisMethod","Conductivities");
        double thresholdPerElement = data("thresholdPerElement", .5);

        vector<std::string> fibrosis_coords2;
        double threshDist = data("threshDist", .1);
        vector<double> fibrosisNodeXvec2;
        vector<double> fibrosisNodeYvec2;
        vector<double> fibrosisNodeZvec2;

        int numPointsToFibrosis = 0;
        if(fibrosisBool && fibrosisMethod.compare(0,11,"Percolation") == 0){
          std::string fibrosis_coords_matrix = data("fibrosis_coords_matrix", "fibrosisCoords.csv");
          std::string line, word;
          fstream file (fibrosis_coords_matrix, ios::in);

          //PULLING FIBROSIS COORDINATES INTO VECTOR
          if(file.is_open())
          {
            while(getline(file, line))
            {
              stringstream str(line);
              fibrosis_coords2.push_back(line);

              std::vector<std::string> CoordsInfoVec;
              stringstream ss(line);
                   
              while (ss.good()){
                  string substr;
                  getline(ss, substr, ',');
                  CoordsInfoVec.push_back(substr);
              }
              fibrosisNodeXvec2.push_back( std::stod(CoordsInfoVec[0]) );
              fibrosisNodeYvec2.push_back( std::stod(CoordsInfoVec[1]) );
              fibrosisNodeZvec2.push_back( std::stod(CoordsInfoVec[2]) );

              numPointsToFibrosis++;
            }
            libMesh::out << "Fibrosis case: File: " << fibrosis_coords_matrix << " found successfully... Eliminating elements from tissue domain..." << std::endl;
            libMesh::out << "Fibrosis method used: " << fibrosisMethod << std::endl;
            libMesh::out << "Fibrosis coordinate points to check: " << numPointsToFibrosis << std::endl;
          }
          else{
            libMesh::out << "Could not open the file" << std::endl;
          }
        }




        bool convergence_test = data("convergence_test", false);

        // Get mesh parameters
        std::string eltype = data("eltype", "simplex");
        std::map < std::string, ElemType > elem_type_map =
        {
        { "TRI3", TRI3 },
        { "TRI6", TRI6 },
        { "QUAD4", QUAD4 },
        { "QUAD9", QUAD9 },
        { "TET4", TET4 },
        { "HEX8", HEX8 },
        {"HEX20",HEX20},
        {"HEX27",HEX27},
        {"EDGE3",EDGE3},
        {"EDGE4",EDGE4} };
        auto elType = TRI3;
        auto elem_type_it = elem_type_map.find(eltype);
        if (elem_type_it != elem_type_map.end())
            elType = elem_type_it->second;

        if(nelz == 0){
          std::cout << "Creating the box [" << minx << ", " << maxx << "] x [" << miny << ", " << maxy << "] " << std::endl;
          std::cout << "Using " << nelx << " x " << nely << " elements." << std::endl;
        }
        else{
          std::cout << "Creating the cube [" << minx << ", " << maxx << "] x [" << miny << ", " << maxy << "] x [" << minz << ", " << maxz << "] " << std::endl;
          std::cout << "Using " << nelx << " x " << nely << " x " << nelz << " elements." << std::endl;
        }
        std::cout << "Element type " << elem_type_it->first << std::endl;

        // Create mesh
        MeshTools::Generation::build_cube(mesh, nelx, nely, nelz, minx, maxx, miny, maxy, minz, maxz, elType);

        // setup subdomains
        // tissue domain box
        double tissue_maxx = data("tissue_maxx", 1.0);
        double tissue_maxy = data("tissue_maxy", 1.0);
        double tissue_maxz = data("tissue_maxz", 1.0);
        double tissue_minx = data("tissue_minx", 0.0);
        double tissue_miny = data("tissue_miny", 0.0);
        double tissue_minz = data("tissue_minz", 0.0);
        libMesh::BoundingBox box(libMesh::Point(tissue_minx, tissue_miny, tissue_minz), libMesh::Point(tissue_maxx, tissue_maxy, tissue_maxz));


        
        //MESH REFINEMENT
        int mesh_refinement_factor_tissue = data("mesh_refinement_factor_tissue", 0);
        int mesh_refinement_factor_bath = data("mesh_refinement_factor_bath", 0);
        MeshRefinement meshR(mesh);

        if(mesh_refinement_factor_tissue > 0){
          libMesh::out << "HERE before meshR" << std::endl;
          for(int iMesh = 0; iMesh < mesh_refinement_factor_tissue; iMesh++){
            for(auto el: mesh.element_ptr_range()){
              auto c = el->centroid();

              if(nelz == 0){
                if (c(1) >= tissue_maxy + ((maxy - miny)/nely)*1.5 && c(1) <= tissue_maxy + ((maxy - miny)/nely)*3.){
                  el -> set_refinement_flag(libMesh::Elem::REFINE);
                }
                else if(c(1) <= tissue_miny - ((maxy - miny)/nely)*1.5 && c(1) >= tissue_miny - ((maxy - miny)/nely)*3.){
                  el -> set_refinement_flag(libMesh::Elem::REFINE);
                }
                else{
                  el -> set_refinement_flag(libMesh::Elem::DO_NOTHING);
                }
              }
              else{
                if (c(2) <= tissue_maxz + ((maxz - minz)/nelz) && c(2) >= tissue_minz - ((maxz - minz)/nelz)){
                  el -> set_refinement_flag(libMesh::Elem::REFINE);
                }
                //else if(c(2) <= tissue_minz - ((maxz - minz)/nelz)*1.5 && c(2) >= tissue_minz - ((maxz - minz)/nelz)*3.){
                  //el -> set_refinement_flag(libMesh::Elem::REFINE);
                //}
                else{
                  el -> set_refinement_flag(libMesh::Elem::DO_NOTHING);
                }
              }
            }
            meshR.refine_elements();
            mesh.prepare_for_use(false);
          }
        }

      

        int numFibrosisPoints = 0, numTissuePoints = 0;
        libMesh::out << "HERE before mesh elemnt ptr range" << std::endl;
        // right_interface
        for (auto el : mesh.element_ptr_range())
        {
            auto c = el->centroid();
            // Are we in the tissue region?
            if (box.contains_point(c))
            {
              if(fibrosisBool && fibrosisMethod.compare(0,11,"Percolation") == 0){
                //check for fibrosis
                auto localNodes = el -> get_nodes();
                auto numNodes = el -> n_nodes();
                std::vector<int> fibBool;// = {2, 2, 2, 2};

                if( elType == HEX8 ){
                  fibBool = {2, 2, 2, 2, 2, 2, 2, 2};
                }
                else if( elType == HEX20 ){
                  fibBool = {2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2};
                }
                else if( elType == TET4 ){
                  fibBool = {2, 2, 2, 2};
                }
                else if( elType == QUAD4 ){
                  fibBool = {2, 2, 2, 2};
                }
                else if( elType == QUAD9 ){
                  fibBool = {2, 2, 2, 2, 2, 2, 2, 2, 2};
                }


                int max = 0, indToUse = 0, mostvalue = 2;
              
                for(int locnode = 0; locnode < numNodes; locnode++){

                  double x = (*localNodes[locnode])(0);
                  double y = (*localNodes[locnode])(1);
                  double z = (*localNodes[locnode])(2);

                  //libMesh::out << std::fixed;
                  //libMesh::out << std::setprecision(4);

                  std::string check_Presence;

                  check_Presence = to_string(round_up(x,4))+","+to_string(round_up(y,4))+","+to_string(round_up(z,4));

                  auto icheck = std::find(fibrosis_coords2.begin(), fibrosis_coords2.end(), check_Presence);
                  //libMesh::out << check_Presence << std::endl;
                  //libMesh::out << x << "," << y << "," << z << std::endl;
                  if(fibrosis_coords2.end() != icheck){
                    //libMesh::out << "Found  first   " << check_Presence << std::endl;
                    //libMesh::out << x << "," << y << "," << z << std::endl;
                    //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);
                    fibBool[locnode] = 1;
                  }
                  else{
                    //fibBool[locnode] = 0;
                    

                    for(int checkCoords = 0; checkCoords < fibrosisNodeXvec2.size(); checkCoords++){
                      double distX = std::abs( round_up(x,8) - fibrosisNodeXvec2[checkCoords] );
                      double distY = std::abs( round_up(y,8) - fibrosisNodeYvec2[checkCoords] );
                      double distZ = std::abs( round_up(z,8) - fibrosisNodeZvec2[checkCoords] );
                      if( distX < .005 && distY < .005 && distZ < .005 ){
                        //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);
                        fibBool[locnode] = 1;
                        //break;
                        //libMesh::out << "Found  second   " << check_Presence << std::endl;
                        //libMesh::out << x << "," << y << "," << z << std::endl;
                      }
                      else{
                        if( fibBool[locnode] == 1 ){
                          //NOTHING TO DO BECAUSE IT'S ALREADY FIBROSIS
                        }
                        else if( fibBool[locnode] == 2 ){
                          //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE); 
                          fibBool[locnode] = 0;
                        }                         
                      }
                    }
                    

                    
                  }
                }

                double avgFib=0., numNonTwo=0.;
                for (int i = 0;i < fibBool.size(); i++){
                  if(fibBool[i] != 2){
                    avgFib += (double) fibBool[i];
                    numNonTwo = numNonTwo + 1.;
                  }
                }
                if( avgFib/numNonTwo >= thresholdPerElement ){
                  el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);
                  numFibrosisPoints += 1;
                }
                else{
                  el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
                  numTissuePoints += 1;
                }

              }
              else{
                  el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
                  //libMesh::out << "x: " << c(0) << "     y: " << c(1) << "     z: " << c(2) << std::endl;
              }
            }
            
        }
        libMesh::out << "Fibrosis density in percentage is: " << " with " << numFibrosisPoints << " fibrosis points and " << numTissuePoints << " tissue points..." << std::endl;
    }

    //mesh.prepare_for_use(true);
    libMesh::out << "HERE after element ptr range" << std::endl;

}


// Initialize TimeData from input
TimeData::TimeData(GetPot &data) :
        time(0.0),
        end_time(data("end_time", 1.0)),
        dt(data("dt", 0.125)),
        timestep(0),
        export_timesteps(data("export_timesteps", 1))
{
}
// Initialize TimeData from input
void TimeData::print()
{
    std::cout << "**************************************** " << std::endl;
    std::cout << "TimeData Parameters: " << std::endl;
    std::cout << "End time: " << end_time << std::endl;
    std::cout << "Dt: " << dt << std::endl;
    std::cout << "Current Timestep: " << timestep << std::endl;
    std::cout << "Export Timesteps: " << export_timesteps << std::endl;
    std::cout << "**************************************** " << std::endl;
}

void assemble_matrices(libMesh::EquationSystems &es, const TimeData &timedata, TimeIntegrator time_integrator, libMesh::Order p_order, const GetPot &data, int rankLoc)
{
    using namespace libMesh;
    // Create vector of BC sidesets ids
    std::set<int> bc_sidesets;
    std::string bcs = data("bcs", "");
    read_bc_list(bcs, bc_sidesets);

    bool fibrosisBool = data("fibrosisBool", false);
    std::string fibrosisMethod = data("fibrosisMethod","Conductivities");

    //vector<vector<double>> fibrosis_coords;
    //vector<vector<double>> tissue_coords;

    vector<std::string> fibrosis_coords2;
    //vector<std::string> tissue_coords2;

    if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
      std::string fibrosis_coords_matrix = data("fibrosis_coords_matrix", "fibrosisCoords.csv");

      std::string line, word;
      fstream file (fibrosis_coords_matrix, ios::in);

      //PULLING FIBROSIS COORDINATES INTO VECTOR
      if(file.is_open())
      {
      while(getline(file, line))
      {
        stringstream str(line);

        fibrosis_coords2.push_back(line);
      }

      libMesh::out << "Fibrosis case: File: " << fibrosis_coords_matrix << " found successfully... Changing conductivities..." << std::endl;
      libMesh::out << "Fibrosis method used: " << fibrosisMethod << std::endl;
      }
      else{
      libMesh::out << "Could not open the file" << std::endl;
      }
      /*
      //PULLING TISSUE COORDINATES INTO VECTOR
      if(file2.is_open())
      {
      while(getline(file2, line))
      {
        local_row.clear();

        stringstream str(line);

        tissue_coords2.push_back(line);

        while(getline(str, word, ',')){
        local_row.push_back(std::stod(word));
        }
        tissue_coords.push_back(local_row);
      }
      }
      else{
      libMesh::out << "Could not open the file" << std::endl;
      }

      libMesh::out << "FIBROTIC COORDINATES" << std::endl;
      for(int i=0;i<fibrosis_coords.size();i++)
      {
      for(int j=0;j<fibrosis_coords[i].size();j++)
      {
        libMesh::out << fibrosis_coords[i][j] << " ";
      }
      libMesh::out << std::endl;
      }

      libMesh::out << "TISSUE COORDINATES" << std::endl;
      for(int i=0;i<tissue_coords.size();i++)
      {
      for(int j=0;j<tissue_coords[i].size();j++)
      {
        libMesh::out << tissue_coords[i][j] << " ";
      }
      libMesh::out << std::endl;
      }

      libMesh::out << "CHECK CHECK CHECK" << std::endl;
      for(int i=0;i<fibrosis_coords2.size();i++)
      {
      libMesh::out << fibrosis_coords2[i] << std::endl;
      }
      */
    }



    bool prescribedFibers = data("prescribedFibers",false);
    
    vector<std::string> CoordsPoints;
    vector<double> CoordsPointsX;
    vector<double> CoordsPointsY;
    vector<double> CoordsPointsZ;

    if(prescribedFibers){
      std::string CoordsPointsCSV = data("CoordsPointsCSV", "fibrosisCoords.csv");

      std::string line, word;
      fstream file (CoordsPointsCSV, ios::in);

      //PULLING POINT COORDINATES INTO VECTOR
      if(file.is_open())
      {
      while(getline(file, line))
      {
        stringstream str(line);
        CoordsPoints.push_back(line);

        std::vector<std::string> CoordsInfoVec;
        stringstream ss(line);
                 
        while (ss.good()){
            string substr;
            getline(ss, substr, ',');
            CoordsInfoVec.push_back(substr);
        }
        CoordsPointsX.push_back( std::stod(CoordsInfoVec[0]) );
        CoordsPointsY.push_back( std::stod(CoordsInfoVec[1]) );
        CoordsPointsZ.push_back( std::stod(CoordsInfoVec[2]) );


      }
      libMesh::out << "Importing the coordinate points to map the fibers: File: " << CoordsPointsCSV << " found successfully..." << std::endl;
      }
      else{
      libMesh::out << "Could not open the file..." << std::endl;
      }
    }



    const MeshBase &mesh = es.get_mesh();
    const unsigned int dim = mesh.mesh_dimension();

    // volume element
    FEType fe_type(p_order);
    std::unique_ptr < FEBase > fe(FEBase::build(dim, fe_type));
    Order qp_order = THIRD;
    if (p_order > 1)
        qp_order = FIFTH;
    libMesh::QGauss qrule(dim, qp_order);
    fe->attach_quadrature_rule(&qrule);
    // surface element
    std::unique_ptr < FEBase > fe_face(FEBase::build(dim, fe_type));
    QGauss qface(dim - 1, qp_order);
    fe_face->attach_quadrature_rule(&qface);

    // quantities for volume integration
    const std::vector<Real> &JxW = fe->get_JxW();
    const std::vector<Point> &q_point = fe->get_xyz();
    const std::vector<std::vector<Real>> &phi = fe->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

    // quantities for surface integration
    const std::vector<Real> &JxW_face = fe_face->get_JxW();
    const std::vector<Point> &q_point_face = fe_face->get_xyz();
    const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
    const std::vector<std::vector<RealGradient>> &dphi_face = fe_face->get_dphi();
    const std::vector<Point> &normal = fe_face->get_normals();

    // define fiber field
    double fx = data("fx", 1.0), fy = data("fy", 0.0), fz = data("fz", 0.0);
    double sx = data("sx", 0.0), sy = data("sy", 1.0), sz = data("sz", 0.0);
    double nx = data("nx", 0.0), ny = data("ny", 0.0), nz = data("nz", 1.0);
    VectorValue<Number>  f0(fx, fy, fz);
    VectorValue<Number>  s0(sx, sy, sz);
    VectorValue<Number>  n0(nx, ny, nz);
    // setup conductivities:
    // Default parameters from
    // Cardiac excitation mechanisms, wavefront dynamics and strength�interval
    // curves predicted by 3D orthotropic bidomain simulations

    // read parameters
    double exportVe_timeX = data("exportVe_timeX", 2.3172);
    double exportVe_timeY = data("exportVe_timeY", 2.3172);
    double exportVe_timeZ = data("exportVe_timeZ", 2.3172);
    double dXmesh = data("dXmesh", .0125);
    double sigma_f_i = data("sigma_f_i", 2.3172);
    double sigma_s_i = data("sigma_s_i", 0.2435); // sigma_t in the paper
    double sigma_n_i = data("sigma_n_i", 0.0569);
    double sigma_f_e = data("sigma_f_e", 1.5448);
    double sigma_s_e = data("sigma_s_e", 1.0438);  // sigma_t in the paper
    double sigma_n_e = data("sigma_n_e", 0.37222);
    double sigma_b_ie = data("sigma_b", 6.0);
    double sigma_torso_ie = data("sigma_torso", 6.0);
    double chi = data("chi", 1e3);
    double Cm = data("Cm", 1.5);
    double penalty = data("penalty", 1e8);
    double interface_penalty = data("interface_penalty", 1e4);
    double zDim = data("nelz",0);
    double nelz = data("nelz",0);
    double stimulus_maxx = data("stimulus_maxx", .5);
    double stimulus_minx = data("stimulus_minx", -.5);
    double stimulus_maxy = data("stimulus_maxy", 1.5);
    double stimulus_miny = data("stimulus_miny", 1.3);
    double stimulus_maxz = data("stimulus_maxz", .05);
    double stimulus_minz = data("stimulus_minz", -.05);
    double tissue_maxx = data("tissue_maxx", .5);
    double tissue_maxy = data("tissue_maxy", .5);
    double tissue_maxz = data("tissue_maxz", .5);
    double tissue_minx = data("tissue_minx", .5);
    double tissue_miny = data("tissue_miny", .5);
    double tissue_minz = data("tissue_minz", .5);
    double maxx = data("maxx", .5);
    double maxy = data("maxy", .5);
    double maxz = data("maxz", .5);
    double minx = data("minx", .5);
    double miny = data("miny", .5);
    double minz = data("minz", .5);
    //double dXmesh = data("dXmesh", .0125);
    bool exportVe_time = data("exportVe_time", false);
    double pi = 3.14159265;

    double IstimV = data("stimulus_amplitude", -1.);
    bool differentFibers = data("differentFibers", false);
    double differentFibersAngle = data("differentFibersAngle", 8.);

    bool monodomainType = data("monodomainType", false);

    double fx2, fy2, sx2, sy2;
    if(differentFibers){
      fx2 = fx*cos(differentFibersAngle*acos(-1.)/180.) + fy*sin(differentFibersAngle*acos(-1.)/180.);
      fy2 = -fx*sin(differentFibersAngle*acos(-1.)/180.) + fy*cos(differentFibersAngle*acos(-1.)/180.);

      sx2 = sx*cos(differentFibersAngle*acos(-1.)/180.) + sy*sin(differentFibersAngle*acos(-1.)/180.);
      sy2 = -sx*sin(differentFibersAngle*acos(-1.)/180.) + sy*cos(differentFibersAngle*acos(-1.)/180.);
    }

    double sigma_f_M = (sigma_f_i*sigma_f_e)/(sigma_f_i + sigma_f_e);
    double sigma_s_M = (sigma_s_i*sigma_s_e)/(sigma_s_i + sigma_s_e);
    double sigma_n_M = (sigma_n_i*sigma_n_e)/(sigma_n_i + sigma_n_e); 

    //THIS DOES NOT MATTER BECAUSE IF MONODOMAIN, THE CODE WILL JUST TAKE THE CONDUCTIVITIES WITH AN i, for intracellular.

    // setup tensors parameters
    // f0 \otimes f0
    TensorValue<Number> fof(fx * fx, fx * fy, fx * fz, fy * fx, fy * fy, fy * fz, fz * fx, fz * fy, fz * fz);
    // s0 \otimes s0
    TensorValue<Number> sos(sx * sx, sx * sy, sx * sz, sy * sx, sy * sy, sy * sz, sz * sx, sz * sy, sz * sz);

    TensorValue<Number> fof2(fx2 * fx2, fx2 * fy2, fx2 * fz, fy2 * fx2, fy2 * fy2, fy2 * fz, fz * fx2, fz * fy2, fz * fz);
    // s0 \otimes s0
    TensorValue<Number> sos2(sx2 * sx2, sx2 * sy2, sx2 * sz, sy2 * sx2, sy2 * sy2, sy2 * sz, sz * sx2, sz * sy2, sz * sz);
    // n0 \otimes n0
    TensorValue<Number> non(nx * nx, nx * ny, nx * nz, ny * nx, ny * ny, ny * nz, nz * nx, nz * ny, nz * nz);

    TensorValue<Number> sigma_i = ( sigma_f_i * fof + sigma_s_i * sos + sigma_n_i * non ) /chi;
    TensorValue<Number> sigma_e = ( sigma_f_e * fof + sigma_s_e * sos + sigma_n_e * non) / chi;
    TensorValue<Number> sigma_i2 = ( sigma_f_i * fof2 + sigma_s_i * sos2 + sigma_n_i * non ) /chi;
    TensorValue<Number> sigma_e2 = ( sigma_f_e * fof2 + sigma_s_e * sos2 + sigma_n_e * non) / chi;
    TensorValue<Number> sigma_b(sigma_b_ie/chi, 0.0, 0.0, 0.0, sigma_b_ie/chi, 0.0, 0.0, 0.0, sigma_b_ie/chi);
    TensorValue<Number> sigma_torso(sigma_torso_ie/chi, 0.0, 0.0, 0.0, sigma_torso_ie/chi, 0.0, 0.0, 0.0, sigma_torso_ie/chi);
    TensorValue<Number> sigma_M = ( sigma_f_M * fof + sigma_s_M * sos + sigma_n_M * non ) /chi;
    TensorValue<Number> sigma_IMono = ( sigma_f_i * fof + sigma_s_i * sos + sigma_n_i * non );
    TensorValue<Number> sigma_IMono2 = ( sigma_f_i * fof2 + sigma_s_i * sos2 + sigma_n_i * non );
    TensorValue<Number> sigma_M2 = ( sigma_f_M * fof2 + sigma_s_M * sos2 + sigma_n_M * non ) /chi;

    //IT'S ALWATS GOING TO BE FOUR VALUES JUST BECAUSE THE MESHES TO USE WILL ALWAYS BE TETRAHEDRONS (4 POINTS)
    std::vector< TensorValue<Number> > sigma_I_fibersV = {sigma_i, sigma_i, sigma_i, sigma_i};
    std::vector< TensorValue<Number> > sigma_E_fibersV = {sigma_e, sigma_e, sigma_e, sigma_e};
    std::vector< TensorValue<Number> > sigma_M_fibersV = {sigma_M, sigma_M, sigma_M, sigma_M};

    DenseMatrix < Number > K;
    DenseMatrix < Number > M;
    DenseVector < Number > M_lumped;
    DenseVector < Number > M_lumped2;
    //DenseVector < Number > R_dist;
    DenseVector < Number > I_extra_stim;
    DenseVector < Number > Fibrosis_CheckV;

    std::vector < dof_id_type > dof_indices;
    std::vector < dof_id_type > parabolic_dof_indices;
    std::vector < dof_id_type > node_dof_indices;

    double scfa = 1.;

    std::vector<std::string> FiberInfo;
    std::vector<double> decidingWhatToUse = {1., 1., 1., 1.};

    if(prescribedFibers){
      
      std::string FiberInfoCSV = data("FiberInfoCSV", "fibrosisCoords.csv");

      std::string line, word;
      fstream file (FiberInfoCSV, ios::in);

      if(file.is_open())
      {
        while(getline(file, line))
        {
          stringstream str(line);
          FiberInfo.push_back(line);
        }
        libMesh::out << "Importing the fiber orientations for each coordinate point: File: " << FiberInfoCSV << " found successfully..." << std::endl;
      }
      else{
        libMesh::out << "Could not open the file..." << std::endl;
      }
    }

    //libMesh::out << "ALL GOOD TILL HERE 10" << std::endl;
    // Assemble matrices differently based on time integrator
    switch (time_integrator)
    {
        case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
        case TimeIntegrator::EXPLICIT_INTRACELLULAR:
        case TimeIntegrator::SEMI_IMPLICIT:
        case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
        {


            std::cout << "Assembling matrices for EXPLICIT EXTRACELLULAR: " << static_cast<int>(time_integrator) << std::endl;
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            parabolic_system.zero_out_matrix_and_rhs = false;
            parabolic_system.assemble_before_solve = false;
            parabolic_system.matrix->zero();
            parabolic_system.get_matrix("Ki").zero();
            parabolic_system.get_matrix("Kmono").zero();
            parabolic_system.get_matrix("Ke").zero();
            parabolic_system.get_matrix("Kg").zero();
            parabolic_system.get_matrix("M").zero();
            parabolic_system.get_vector("ML").zero();
            //parabolic_system.get_vector("ML2").zero();
            //parabolic_system.get_vector("R_dist").zero();
            DenseMatrix < Number > Ki;
            DenseMatrix < Number > Kmono;
            DenseMatrix < Number > Ke;
            DenseMatrix < Number > Kg;
            DenseMatrix < Number > Ke_interface;
            DenseMatrix < Number > S;

            bool checkMatrixV = false;
            const DofMap & parabolic_dof_map = parabolic_system.get_dof_map();

            for (auto elem : mesh.active_local_element_ptr_range())
            {
                auto c = elem -> centroid();
                parabolic_dof_map.dof_indices(elem, parabolic_dof_indices);
                fe->reinit(elem);
                int parabolic_ndofs = parabolic_dof_indices.size();

                std::vector<double> localfx = {2.,2.,2.,2.}, localfy = {2.,2.,2.,2.}, localfz = {2.,2.,2.,2.};
                int max = 0, indToUse = 0;

                if(prescribedFibers){
                  std::string check_Presence;

                  auto localNodes = elem -> get_nodes();
                  auto numNodes = elem -> n_nodes();
                  //libMesh::out << "size is: " << sizeof(localNodes) << " and type is: " << typeid(localNodes).name() << std::endl;
                  //libMesh::out << "Checking how many per element " << numNodes << std::endl;
                  for(int locnode = 0; locnode < numNodes; locnode++){

                    double x = (*localNodes[locnode])(0);
                    double y = (*localNodes[locnode])(1);
                    double z = (*localNodes[locnode])(2);

                    check_Presence = to_string( round_up(x,4) )+","+to_string( round_up(y,4) )+","+to_string( round_up(z,4));
                    //libMesh::out << "Local point:" << std::endl;
                    //libMesh::out << check_Presence << std::endl;
                    //libMesh::out << "Point in CSV:" << std::endl;
                    //libMesh::out << CoordsPoints[0] << std::endl;

                    int elementToUse;

                    auto icheck = std::find(CoordsPoints.begin(), CoordsPoints.end(), check_Presence);
                    if(CoordsPoints.end() != icheck){
                      elementToUse = ( icheck - CoordsPoints.begin() );
                      //libMesh::out << "CSV points: " << to_string(CoordsPointsX[elementToUse]) << "," << to_string(CoordsPointsY[elementToUse]) << "," << to_string(CoordsPointsZ[elementToUse])  << std::endl;
                      //libMesh::out << "Local point: " << (check_Presence)  << std::endl;
                    }
                    else{
                      for(int checkCoords = 0; checkCoords < CoordsPointsX.size(); checkCoords++){
                        double distX = std::abs( round_up(x,8) - CoordsPointsX[checkCoords] );
                        double distY = std::abs( round_up(y,8) - CoordsPointsY[checkCoords] );
                        double distZ = std::abs( round_up(z,8) - CoordsPointsZ[checkCoords] );
                        if( distX < .005 && distY < .005 && distZ < .005 ){
                          elementToUse = checkCoords;
                          //libMesh::out << "Distance: " << to_string(distX) << "," << to_string(distY) << "," << to_string(distZ) << " with node being: " << to_string(elementToUse)  << std::endl;
                          //libMesh::out << "CSV points: " << to_string(CoordsPointsX[checkCoords]) << "," << to_string(CoordsPointsY[checkCoords]) << "," << to_string(CoordsPointsZ[checkCoords])  << std::endl;
                          //libMesh::out << "Local point: " << (check_Presence)  << std::endl;
                        }
                        else{
                          //libMesh::out << "Points not found at all... something is wrong..." << std::endl;
                          
                        }

                      }
                      //libMesh::out << "Point not found at all... something is wrong..." << std::endl;
                      //libMesh::out << check_Presence << std::endl;
                      //libMesh::out << CoordsPointsX.size() << std::endl;
                      //libMesh::out << CoordsPointsY.size() << std::endl;
                      //libMesh::out << CoordsPointsZ.size() << std::endl;
                    }

                    std::string locFiberInfo;
                    //if( elementToUse >= FiberInfo.size() ){
                    //}
                    //else{
                      locFiberInfo = FiberInfo[elementToUse];
                    //}

                    std::vector<std::string> FiberInfoVec;
                    stringstream ss(locFiberInfo);
                 
                    while (ss.good()) {
                        string substr;
                        getline(ss, substr, ',');
                        FiberInfoVec.push_back(substr);
                    }

                    //if(FiberInfoVec.size() > 0){
                      fx = std::stod(FiberInfoVec[0]);
                      fy = std::stod(FiberInfoVec[1]);
                      fz = std::stod(FiberInfoVec[2]);

                      nx = std::stod(FiberInfoVec[3]);
                      ny = std::stod(FiberInfoVec[4]);
                      nz = std::stod(FiberInfoVec[5]);

                      sx = std::stod(FiberInfoVec[6]);
                      sy = std::stod(FiberInfoVec[7]);
                      sz = std::stod(FiberInfoVec[8]);
                    //}

                    decidingWhatToUse[locnode] = fx;
                    localfx[locnode] = fx; localfy[locnode] = fy; localfz[locnode] = fz;


                    TensorValue<Number> fofp(fx * fx, fx * fy, fx * fz, fy * fx, fy * fy, fy * fz, fz * fx, fz * fy, fz * fz);
                    TensorValue<Number> sosp(sx * sx, sx * sy, sx * sz, sy * sx, sy * sy, sy * sz, sz * sx, sz * sy, sz * sz);
                    TensorValue<Number> nonp(nx * nx, nx * ny, nx * nz, ny * nx, ny * ny, ny * nz, nz * nx, nz * ny, nz * nz);

                    sigma_I_fibersV[locnode] = ( sigma_f_i * fofp + sigma_s_i * sosp + sigma_n_i * nonp ) /chi;
                    sigma_E_fibersV[locnode] = ( sigma_f_e * fofp + sigma_s_e * sosp + sigma_n_e * nonp ) / chi;
                    sigma_M_fibersV[locnode] = ( sigma_f_M * fofp + sigma_s_M * sosp + sigma_n_M * nonp ) /chi;


                    //CHECK WHAT IS THE MOST COMMON VALUE? MAKE THAT THE ELEMENT CONDUCTIVITY TENSOR, NOT A DIFFERENT ONE BASED ON VECTOR LOCALITY.

                    
                    libMesh::LinearImplicitSystem & fiberSoln = es.add_system <libMesh::LinearImplicitSystem> ("fiberSolution");
                    const DofMap & recovery_dof_map = fiberSoln.get_dof_map();

                    recovery_dof_map.dof_indices(localNodes[locnode], node_dof_indices);

                    const unsigned int fibx_var = fiberSoln.variable_number ("fibx");
                    const unsigned int fiby_var = fiberSoln.variable_number ("fiby");
                    const unsigned int fibz_var = fiberSoln.variable_number ("fibz");
                    const unsigned int sheetx_var = fiberSoln.variable_number ("sheetx");
                    const unsigned int sheety_var = fiberSoln.variable_number ("sheety");
                    const unsigned int sheetz_var = fiberSoln.variable_number ("sheetz");
                    const unsigned int crossfx_var = fiberSoln.variable_number ("crossfx");
                    const unsigned int crossfy_var = fiberSoln.variable_number ("crossfy");
                    const unsigned int crossfz_var = fiberSoln.variable_number ("crossfz");
                    fiberSoln.solution -> set(node_dof_indices[fibx_var],fx);
                    fiberSoln.solution -> set(node_dof_indices[fiby_var],fy);
                    fiberSoln.solution -> set(node_dof_indices[fibz_var],fz);

                    fiberSoln.solution -> set(node_dof_indices[sheetx_var],nx);
                    fiberSoln.solution -> set(node_dof_indices[sheety_var],ny);
                    fiberSoln.solution -> set(node_dof_indices[sheetz_var],nz);

                    fiberSoln.solution -> set(node_dof_indices[crossfx_var],sx);
                    fiberSoln.solution -> set(node_dof_indices[crossfy_var],sy);
                    fiberSoln.solution -> set(node_dof_indices[crossfz_var],sz);
                    

                  }

                  for (int i = 0;i < decidingWhatToUse.size(); i++)
                  {
                      int co = (int)count(decidingWhatToUse.begin(), decidingWhatToUse.end(), decidingWhatToUse[i]);
                      if (co > max)
                      {
                          indToUse = i;
                          max = co;
                          //mostvalue = decidingWhatToUse[i];
                          //cout << "Most value: " + to_string(mostvalue) << endl;
                      }
                  }

                }

                
                double xcen = c(0), ycen = c(1), zcen = c(2);
                //libMesh::out << "At centroid point: (" << xcen << "," << ycen << "," << zcen << ") the fiber direction is: <" << (localfx[indToUse]) << "," << (localfy[indToUse]) << "," << (localfz[indToUse]) << ">  with conductivity tensor being: " << (sigma_I_fibersV[indToUse]) << std::endl;
                //libMesh::out << "At centroid point: (" << xcen << "," << ycen << "," << zcen << ") the fiber direction is: <" << (localfx[indToUse]) << "," << (localfy[indToUse]) << "," << (localfz[indToUse]) << ">  with conductivity tensor being: " << (sigma_I_fibersV[indToUse]) << std::endl;
                //libMesh::out << "Previous tensor: " << (sigma_i) << std::endl;

                // resize local elemental matrices
                Ki.resize(parabolic_ndofs, parabolic_ndofs);
                Kmono.resize(parabolic_ndofs, parabolic_ndofs);
                Ke.resize(parabolic_ndofs, parabolic_ndofs);
                Kg.resize(parabolic_ndofs, parabolic_ndofs);
                M.resize(parabolic_ndofs, parabolic_ndofs);
                K.resize(parabolic_ndofs, parabolic_ndofs);
                M_lumped.resize(parabolic_ndofs);
                M_lumped2.resize(parabolic_ndofs);
                //R_dist.resize(parabolic_ndofs);
                if(fibrosisMethod.compare(0,14,"Conductivities") == 0 && fibrosisBool){
                  Fibrosis_CheckV.resize(parabolic_ndofs);
                }



                auto subdomain_id = elem->subdomain_id();

                //std::cout << "Loop over volume: " << std::flush;

                if(subdomain_id == static_cast<short unsigned int>(Subdomain::TISSUE))
                {
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i != parabolic_ndofs; i++)
                        {
                            //libMesh::out << parabolic_ndofs << std::endl;
                            double x = c(0);
                            double y = c(1);
                            double z = c(2);

                            const double xT = q_point[qp](0);
                            const double yT = q_point[qp](1);
                            const double zT = q_point[qp](2);

                            /*
                            if( distanceValue < threshDist*2. ){
                              scfa = 1.;
                            }
                            else if( (distanceValue >= threshDist*2.) && (distanceValue < threshDist*4.) ){
                              scfa = 2.;
                            }
                            else if( (distanceValue >= threshDist*4.) && (distanceValue < threshDist*7.) ){
                              scfa = 4.;
                            }
                            else if( (distanceValue >= threshDist*7.) ){
                              scfa = 8.;
                            }
                            */
                            
                            double locdistance = ( sqrt( (exportVe_timeX - xT)*(exportVe_timeX - xT) + (exportVe_timeY - yT)*(exportVe_timeY - yT) + (exportVe_timeZ - zT)*(exportVe_timeZ - zT) ) );

                            if(locdistance < .075){
                              checkMatrixV = true;
                            }

                            //R_dist(i) += ( 1. / ( sqrt( (exportVe_timeX - xT)*(exportVe_timeX - xT) + (exportVe_timeY - yT)*(exportVe_timeY - yT) + (exportVe_timeZ - zT)*(exportVe_timeZ - zT) ) ) );
                            
                            for (unsigned int j = 0; j != parabolic_ndofs; j++)
                            {
                              double factorDiv = 1.;
                              auto locSigma_i = sigma_i;
                              auto locSigma_e = sigma_e;

                              if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                                std::string check_Presence;

                                if(nelz == 0){
                                  check_Presence = to_string(round_up(x,8))+","+to_string(round_up(y,8));
                                }
                                else{
                                  check_Presence = to_string(round_up(x,8))+","+to_string(round_up(y,8))+","+to_string(round_up(z,8));
                                }

                                auto icheck = std::find(fibrosis_coords2.begin(), fibrosis_coords2.end(), check_Presence);
                                //libMesh::out << check_Presence << std::endl;
                                if(fibrosis_coords2.end() != icheck){
                                //libMesh::out << "Found     " << check_Presence << "    at row " << " col " << i-fibrosis_coords2.begin() << '\n';
                                //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);

                                  factorDiv = 100000.;
                                  Fibrosis_CheckV(i) += JxW[qp]*(     phi[i][qp]*phi[j][qp]*((10./1.)*IstimV)     );

                                  //locSigma_e = ( 0. * fof + 0. * sos + 0. * non ) /chi;
                                  //locSigma_i = ( 0. * fof + 0. * sos + 0. * non ) /chi;

                                }
                                else{
                                //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
                                  //Fibrosis_CheckV(i) += JxW[qp]*(     0.    );
                                }
                              }

                                // Elliptic equation matrix
                              if(differentFibers){
                                if(z > 0){
                                  K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_i/factorDiv * dphi[j][qp]);
                                  Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i/factorDiv) * dphi[j][qp]);
                                  Kmono(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_M/factorDiv) * dphi[j][qp]);
                                  Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_e/factorDiv) * dphi[j][qp]);
                                  
                                }
                                else{
                                  K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_i2/factorDiv * dphi[j][qp]);
                                  Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i2/factorDiv) * dphi[j][qp]);
                                  Kmono(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_M2/factorDiv) * dphi[j][qp]);
                                  Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_e2/factorDiv) * dphi[j][qp]);
                                }
                              }
                              else{
                                K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_I_fibersV[indToUse]/factorDiv * dphi[j][qp]);
                                Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_I_fibersV[indToUse]/factorDiv) * dphi[j][qp]);
                                //Kmono(i, j) += JxW[qp] * (( 1. / ( phi[i][qp] * sqrt( (pow((exportVe_timeX - xT),2)) + (pow((exportVe_timeY - yT),2)) + (pow((exportVe_timeZ - zT),2)) ) ) )) * dphi[i][qp] * ((sigma_IMono/factorDiv) * dphi[j][qp]);
                                Kmono(i, j) += JxW[qp]  * dphi[i][qp] * ((sigma_M_fibersV[indToUse]) * dphi[j][qp]);
                                Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_E_fibersV[indToUse]/factorDiv) * dphi[j][qp]);
                              }

                                // Parabolic matrices
                                // add mass matrix
                                M_lumped(i) += JxW[qp] * Cm / timedata.dt * phi[i][qp] * phi[j][qp];
                                M_lumped2(i) += JxW[qp] * phi[i][qp] * phi[j][qp];


                                // stiffness matrix
                                M(i,  j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                                K(i, i)  += JxW[qp] * Cm / timedata.dt * phi[i][qp] * phi[j][qp];
                            }
                        }
                    }
                    //IF FOR SOME REASON THERE NEEDS TO BE A GROUNDED NODE (VM = 0)
                    if( maxz - tissue_maxz == 0. || minz - tissue_minz == 0. ){
                      for (auto side : elem->side_index_range())
                      {
                          if (elem->neighbor_ptr(side) == nullptr)
                          {
                              auto boundary_id = mesh.get_boundary_info().boundary_id(elem, side);
                              auto it = bc_sidesets.find(static_cast<int>(boundary_id));
                              if (it != bc_sidesets.end())
                              {
                                  //std::cout << "found_boundary " << std::endl;
                                  fe_face->reinit(elem, side);
                                  for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                                  {
                                      for (unsigned int i = 0; i != parabolic_ndofs; i++)
                                      {
                                          for (unsigned int j = 0; j != parabolic_ndofs; j++)
                                          {
                                              K(i, j) += JxW_face[qp] * penalty * phi_face[j][qp] * phi_face[i][qp];
                                          }
                                      }
                                  }
                              }
                          }
                      }
                    }

                }


                /*
                else if(subdomain_id == static_cast<short unsigned int>(Subdomain::FIBROSIS))
                {
                  for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                  {
                    for (unsigned int i = 0; i != elliptic_ndofs; i++)
                    {
                      for (unsigned int j = 0; j != elliptic_ndofs; j++)
                      {
                        
                        
                      }
                    }
                  }
                }
                */
                //std::cout << ", Loop over interface: " << std::flush;

                // interface
                //if( time_integrator == TimeIntegrator::EXPLICIT_EXTRACELLULAR )
                //{
                //}
                // boundaries
                //std::cout << ", Loop over boundaries: " << std::flush;
                
                //std::cout << "Add matrix " << std::endl;
                parabolic_system.matrix->add_matrix(K, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("M").add_matrix(M, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_vector("ML").add_vector(M_lumped, parabolic_dof_indices);
                //parabolic_system.get_vector("ML2").add_vector(M_lumped2, parabolic_dof_indices);
                //parabolic_system.get_vector("R_dist").add_vector(R_dist, parabolic_dof_indices);
                if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                  parabolic_system.get_vector("Fibrosis_Check").add_vector(Fibrosis_CheckV, parabolic_dof_indices);
                }
                parabolic_system.get_matrix("Ki").add_matrix(Ki, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Kmono").add_matrix(Kmono, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Ke").add_matrix(Ke, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Kg").add_matrix(Kg, parabolic_dof_indices, parabolic_dof_indices);

                if(checkMatrixV){
                  //libMesh::out << "For element " << *elem << "with size " << parabolic_ndofs << "x" << parabolic_ndofs << " at point within the tissue the matrix is: " << std::endl;
                  //libMesh::out << Kmono << std::endl;
                }

            }
            std::cout << "closing" << std::endl;
            if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
              parabolic_system.get_vector("Fibrosis_Check").close();
            }
            parabolic_system.matrix->close();
            parabolic_system.get_matrix("M").close();
            parabolic_system.get_vector("ML").close();
            //parabolic_system.get_vector("ML2").close();
            //parabolic_system.get_vector("R_dist").close();
            parabolic_system.get_matrix("Ki").close();
            parabolic_system.get_matrix("Kmono").close();
            parabolic_system.get_matrix("Ke").close();
            parabolic_system.get_matrix("Kg").close();
            std::cout << "done" << std::endl;

            if(prescribedFibers){
              libMesh::LinearImplicitSystem & fiberSoln = es.get_system < libMesh::LinearImplicitSystem > ("fiberSolution");
              fiberSoln.solution -> close();
              fiberSoln.update();
            }



            break;
        }
        case TimeIntegrator::SBDF1:
            // Cm * M * ( V^n+1 - Vn ) / dt + Ki * V^n+1 + Ki * Ve^n+1 = -I^n
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //         -         -
            //        | Cm /dt * M + Ki,   Ki    |
            //   K  = | Ki,                Ki+Ke |;
            //         -         -
        case TimeIntegrator::SBDF2:
            // Cm * M * ( 3 * V^n+1 - 4 * Vn + V^n-1 ) / (2*dt) + Ki * V^n+1 + Ki * Ve^n+1 = -2*I^n + I^n-1
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //         -         -
            //        | 3/2 * Cm /dt * M + Ki,   Ki    |
            //   K  = | Ki,                      Ki+Ke |;
            //         -         -
        case TimeIntegrator::SBDF3:
            // Cm * M * ( 11/6 * V^n+1 - 3 * Vn + 3/2 V^n-1 -1/3 V^n-2) / (dt) + Ki * V^n+1 + Ki * Ve^n+1 = -3*I^n + 3*I^n-1 - I^n-2
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //         -         -
            //        | 11/6 * Cm /dt * M + Ki,   Ki    |
            //   K  = | Ki,                      Ki+Ke |;
            //         -         -
        default:
        {

            std::cout << "Assembling matrices for SBDF" << static_cast<int>(time_integrator) << std::endl;
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            parabolic_system.zero_out_matrix_and_rhs = false;
            parabolic_system.assemble_before_solve = false;
            parabolic_system.matrix->zero();
            parabolic_system.get_matrix("Ki").zero();
            parabolic_system.get_matrix("Kmono").zero();
            parabolic_system.get_matrix("Ke").zero();
            parabolic_system.get_matrix("Kg").zero();
            parabolic_system.get_matrix("M").zero();
            parabolic_system.get_vector("ML").zero();
            //parabolic_system.get_vector("ML2").zero();
            //parabolic_system.get_vector("R_dist").zero();
            DenseMatrix < Number > Ki;
            DenseMatrix < Number > Kmono;
            DenseMatrix < Number > Ke;
            DenseMatrix < Number > Kg;
            DenseMatrix < Number > Ke_interface;
            DenseMatrix < Number > S;

            double coefficient = Cm / timedata.dt;

            if( time_integrator == TimeIntegrator::SBDF2 ){
              coefficient *= 1.5;
            }
            else if( time_integrator == TimeIntegrator::SBDF3 ){
              coefficient *= 11. / 6.;
            }

            bool checkMatrixV = false;
            const DofMap & parabolic_dof_map = parabolic_system.get_dof_map();

            for (auto elem : mesh.active_local_element_ptr_range())
            {
                auto c = elem -> centroid();
                parabolic_dof_map.dof_indices(elem, parabolic_dof_indices);
                fe->reinit(elem);
                int parabolic_ndofs = parabolic_dof_indices.size();


                std::vector<double> localfx = {2.,2.,2.,2.}, localfy = {2.,2.,2.,2.}, localfz = {2.,2.,2.,2.};
                int max = 0, indToUse = 0;

                if(prescribedFibers){
                  std::string check_Presence;

                  auto localNodes = elem -> get_nodes();
                  auto numNodes = elem -> n_nodes();
                  //libMesh::out << "size is: " << sizeof(localNodes) << " and type is: " << typeid(localNodes).name() << std::endl;
                  //libMesh::out << "Checking how many per element " << numNodes << std::endl;
                  for(int locnode = 0; locnode < numNodes; locnode++){

                    double x = (*localNodes[locnode])(0);
                    double y = (*localNodes[locnode])(1);
                    double z = (*localNodes[locnode])(2);

                    check_Presence = to_string( round_up(x,4) )+","+to_string( round_up(y,4) )+","+to_string( round_up(z,4) );
                    //libMesh::out << "Local point:" << std::endl;
                    //libMesh::out << check_Presence << std::endl;
                    //libMesh::out << "Point in CSV:" << std::endl;
                    //libMesh::out << CoordsPoints[0] << std::endl;

                    int elementToUse;

                    auto icheck = std::find(CoordsPoints.begin(), CoordsPoints.end(), check_Presence);
                    if(CoordsPoints.end() != icheck){
                      elementToUse = ( icheck - CoordsPoints.begin() );
                      //libMesh::out << "CSV points: " << to_string(CoordsPointsX[elementToUse]) << "," << to_string(CoordsPointsY[elementToUse]) << "," << to_string(CoordsPointsZ[elementToUse])  << std::endl;
                      //libMesh::out << "Local point: " << (check_Presence)  << std::endl;
                    }
                    else{
                      for(int checkCoords = 0; checkCoords < CoordsPointsX.size(); checkCoords++){
                        double distX = std::abs( round_up(x,8) - CoordsPointsX[checkCoords] );
                        double distY = std::abs( round_up(y,8) - CoordsPointsY[checkCoords] );
                        double distZ = std::abs( round_up(z,8) - CoordsPointsZ[checkCoords] );
                        if( distX < .005 && distY < .005 && distZ < .005 ){
                          elementToUse = checkCoords;
                          //libMesh::out << "Distance: " << to_string(distX) << "," << to_string(distY) << "," << to_string(distZ) << " with node being: " << to_string(elementToUse)  << std::endl;
                          //libMesh::out << "CSV points: " << to_string(CoordsPointsX[checkCoords]) << "," << to_string(CoordsPointsY[checkCoords]) << "," << to_string(CoordsPointsZ[checkCoords])  << std::endl;
                          //libMesh::out << "Local point: " << (check_Presence)  << std::endl;
                        }
                        else{
                          //libMesh::out << "Points not found at all... something is wrong..." << std::endl;
                          
                        }

                      }
                      //libMesh::out << "Point not found at all... something is wrong..." << std::endl;
                      //libMesh::out << check_Presence << std::endl;
                      //libMesh::out << CoordsPointsX.size() << std::endl;
                      //libMesh::out << CoordsPointsY.size() << std::endl;
                      //libMesh::out << CoordsPointsZ.size() << std::endl;
                    }

                    std::string locFiberInfo;
                    //if( elementToUse >= FiberInfo.size() ){
                    //}
                    //else{
                      locFiberInfo = FiberInfo[elementToUse];
                    //}

                    std::vector<std::string> FiberInfoVec;
                    stringstream ss(locFiberInfo);
                 
                    while (ss.good()) {
                        string substr;
                        getline(ss, substr, ',');
                        FiberInfoVec.push_back(substr);
                    }

                    //if(FiberInfoVec.size() > 0){
                      fx = std::stod(FiberInfoVec[0]);
                      fy = std::stod(FiberInfoVec[1]);
                      fz = std::stod(FiberInfoVec[2]);

                      nx = std::stod(FiberInfoVec[3]);
                      ny = std::stod(FiberInfoVec[4]);
                      nz = std::stod(FiberInfoVec[5]);

                      sx = std::stod(FiberInfoVec[6]);
                      sy = std::stod(FiberInfoVec[7]);
                      sz = std::stod(FiberInfoVec[8]);
                    //}

                    decidingWhatToUse[locnode] = fx;
                    localfx[locnode] = fx; localfy[locnode] = fy; localfz[locnode] = fz;


                    TensorValue<Number> fofp(fx * fx, fx * fy, fx * fz, fy * fx, fy * fy, fy * fz, fz * fx, fz * fy, fz * fz);
                    TensorValue<Number> sosp(sx * sx, sx * sy, sx * sz, sy * sx, sy * sy, sy * sz, sz * sx, sz * sy, sz * sz);
                    TensorValue<Number> nonp(nx * nx, nx * ny, nx * nz, ny * nx, ny * ny, ny * nz, nz * nx, nz * ny, nz * nz);

                    sigma_I_fibersV[locnode] = ( sigma_f_i * fofp + sigma_s_i * sosp + sigma_n_i * nonp ) /chi;
                    sigma_E_fibersV[locnode] = ( sigma_f_e * fofp + sigma_s_e * sosp + sigma_n_e * nonp ) / chi;
                    sigma_M_fibersV[locnode] = ( sigma_f_M * fofp + sigma_s_M * sosp + sigma_n_M * nonp ) /chi;


                    //CHECK WHAT IS THE MOST COMMON VALUE? MAKE THAT THE ELEMENT CONDUCTIVITY TENSOR, NOT A DIFFERENT ONE BASED ON VECTOR LOCALITY.

                    
                    libMesh::LinearImplicitSystem & fiberSoln = es.add_system <libMesh::LinearImplicitSystem> ("fiberSolution");
                    const DofMap & recovery_dof_map = fiberSoln.get_dof_map();

                    recovery_dof_map.dof_indices(localNodes[locnode], node_dof_indices);

                    const unsigned int fibx_var = fiberSoln.variable_number ("fibx");
                    const unsigned int fiby_var = fiberSoln.variable_number ("fiby");
                    const unsigned int fibz_var = fiberSoln.variable_number ("fibz");
                    const unsigned int sheetx_var = fiberSoln.variable_number ("sheetx");
                    const unsigned int sheety_var = fiberSoln.variable_number ("sheety");
                    const unsigned int sheetz_var = fiberSoln.variable_number ("sheetz");
                    const unsigned int crossfx_var = fiberSoln.variable_number ("crossfx");
                    const unsigned int crossfy_var = fiberSoln.variable_number ("crossfy");
                    const unsigned int crossfz_var = fiberSoln.variable_number ("crossfz");
                    fiberSoln.solution -> set(node_dof_indices[fibx_var],fx);
                    fiberSoln.solution -> set(node_dof_indices[fiby_var],fy);
                    fiberSoln.solution -> set(node_dof_indices[fibz_var],fz);

                    fiberSoln.solution -> set(node_dof_indices[sheetx_var],nx);
                    fiberSoln.solution -> set(node_dof_indices[sheety_var],ny);
                    fiberSoln.solution -> set(node_dof_indices[sheetz_var],nz);

                    fiberSoln.solution -> set(node_dof_indices[crossfx_var],sx);
                    fiberSoln.solution -> set(node_dof_indices[crossfy_var],sy);
                    fiberSoln.solution -> set(node_dof_indices[crossfz_var],sz);
                    

                  }

                  for (int i = 0;i < decidingWhatToUse.size(); i++)
                  {
                      int co = (int)count(decidingWhatToUse.begin(), decidingWhatToUse.end(), decidingWhatToUse[i]);
                      if (co > max)
                      {
                          indToUse = i;
                          max = co;
                          //mostvalue = decidingWhatToUse[i];
                          //cout << "Most value: " + to_string(mostvalue) << endl;
                      }
                  }

                }

                
                double xcen = c(0), ycen = c(1), zcen = c(2);
                //libMesh::out << "At centroid point: (" << xcen << "," << ycen << "," << zcen << ") the fiber direction is: <" << (localfx[indToUse]) << "," << (localfy[indToUse]) << "," << (localfz[indToUse]) << ">  with conductivity tensor being: " << (sigma_I_fibersV[indToUse]) << std::endl;
                //libMesh::out << "At centroid point: (" << xcen << "," << ycen << "," << zcen << ") the fiber direction is: <" << (localfx[indToUse]) << "," << (localfy[indToUse]) << "," << (localfz[indToUse]) << ">  with conductivity tensor being: " << (sigma_I_fibersV[indToUse]) << std::endl;
                //libMesh::out << "Previous tensor: " << (sigma_i) << std::endl;

                // resize local elemental matrices
                Ki.resize(parabolic_ndofs, parabolic_ndofs);
                Kmono.resize(parabolic_ndofs, parabolic_ndofs);
                Ke.resize(parabolic_ndofs, parabolic_ndofs);
                Kg.resize(parabolic_ndofs, parabolic_ndofs);
                M.resize(parabolic_ndofs, parabolic_ndofs);
                K.resize(parabolic_ndofs, parabolic_ndofs);
                M_lumped.resize(parabolic_ndofs);
                M_lumped2.resize(parabolic_ndofs);
                //R_dist.resize(parabolic_ndofs);
                if(fibrosisMethod.compare(0,14,"Conductivities") == 0 && fibrosisBool){
                  Fibrosis_CheckV.resize(parabolic_ndofs);
                }




                auto subdomain_id = elem->subdomain_id();

                //std::cout << "Loop over volume: " << std::flush;

                if(subdomain_id == static_cast<short unsigned int>(Subdomain::TISSUE))
                {
                    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                    {
                        for (unsigned int i = 0; i != parabolic_ndofs; i++)
                        {

                            double x = c(0);
                            double y = c(1);
                            double z = c(2);

                            const double xT = q_point[qp](0);
                            const double yT = q_point[qp](1);
                            const double zT = q_point[qp](2);

                            /*
                            if( distanceValue < threshDist*2. ){
                              scfa = 1.;
                            }
                            else if( (distanceValue >= threshDist*2.) && (distanceValue < threshDist*4.) ){
                              scfa = 2.;
                            }
                            else if( (distanceValue >= threshDist*4.) && (distanceValue < threshDist*7.) ){
                              scfa = 4.;
                            }
                            else if( (distanceValue >= threshDist*7.) ){
                              scfa = 8.;
                            }
                            */
                            
                            double locdistance = ( sqrt( (exportVe_timeX - xT)*(exportVe_timeX - xT) + (exportVe_timeY - yT)*(exportVe_timeY - yT) + (exportVe_timeZ - zT)*(exportVe_timeZ - zT) ) );

                            if(locdistance < .075){
                              checkMatrixV = true;
                            }

                            //R_dist(i) += ( 1. / ( sqrt( (exportVe_timeX - xT)*(exportVe_timeX - xT) + (exportVe_timeY - yT)*(exportVe_timeY - yT) + (exportVe_timeZ - zT)*(exportVe_timeZ - zT) ) ) );
                            
                            for (unsigned int j = 0; j != parabolic_ndofs; j++)
                            {
                              double factorDiv = 1.;
                              auto locSigma_i = sigma_i;
                              auto locSigma_e = sigma_e;

                              if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                                std::string check_Presence;

                                if(nelz == 0){
                                  check_Presence = to_string(round_up(x,8))+","+to_string(round_up(y,8));
                                }
                                else{
                                  check_Presence = to_string(round_up(x,8))+","+to_string(round_up(y,8))+","+to_string(round_up(z,8));
                                }

                                auto icheck = std::find(fibrosis_coords2.begin(), fibrosis_coords2.end(), check_Presence);
                                //libMesh::out << check_Presence << std::endl;
                                if(fibrosis_coords2.end() != icheck){
                                //libMesh::out << "Found     " << check_Presence << "    at row " << " col " << i-fibrosis_coords2.begin() << '\n';
                                //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::FIBROSIS);

                                  factorDiv = 10000.;
                                  Fibrosis_CheckV(i) += JxW[qp]*(     phi[i][qp]*phi[j][qp]*((10./1.)*IstimV)     );

                                  //locSigma_e = ( 0. * fof + 0. * sos + 0. * non ) /chi;
                                  //locSigma_i = ( 0. * fof + 0. * sos + 0. * non ) /chi;

                                }
                                else{
                                //el->subdomain_id() = static_cast<libMesh::subdomain_id_type>(Subdomain::TISSUE);
                                  //Fibrosis_CheckV(i) += JxW[qp]*(     0.    );
                                }
                              }

                                // Elliptic equation matrix
                              if(differentFibers){
                                if(z > 0){
                                  K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_i/factorDiv * dphi[j][qp]);
                                  Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i/factorDiv) * dphi[j][qp]);
                                  Kmono(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_M/factorDiv) * dphi[j][qp]);
                                  Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_e/factorDiv) * dphi[j][qp]);
                                  
                                }
                                else{
                                  K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_i2/factorDiv * dphi[j][qp]);
                                  Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_i2/factorDiv) * dphi[j][qp]);
                                  Kmono(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_M2/factorDiv) * dphi[j][qp]);
                                  Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_e2/factorDiv) * dphi[j][qp]);
                                }
                              }
                              else{
                                K(i, j) += JxW[qp] * dphi[i][qp] * (sigma_I_fibersV[indToUse]/factorDiv * dphi[j][qp]);
                                Ki(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_I_fibersV[indToUse]/factorDiv) * dphi[j][qp]);
                                //Kmono(i, j) += JxW[qp] * (( 1. / ( phi[i][qp] * sqrt( (pow((exportVe_timeX - xT),2)) + (pow((exportVe_timeY - yT),2)) + (pow((exportVe_timeZ - zT),2)) ) ) )) * dphi[i][qp] * ((sigma_IMono/factorDiv) * dphi[j][qp]);
                                Kmono(i, j) += JxW[qp]  * dphi[i][qp] * ((sigma_M_fibersV[indToUse]) * dphi[j][qp]);
                                Ke(i, j) += JxW[qp] * dphi[i][qp] * ((sigma_E_fibersV[indToUse]/factorDiv) * dphi[j][qp]);
                              }

                                // Parabolic matrices
                                // add mass matrix
                                M_lumped(i) += JxW[qp] * Cm / timedata.dt * phi[i][qp] * phi[j][qp];
                                M_lumped2(i) += JxW[qp] * phi[i][qp] * phi[j][qp];


                                // stiffness matrix
                                M(i,  j) += JxW[qp] * phi[i][qp] * phi[j][qp];
                                K(i, i)  += JxW[qp] * coefficient * phi[i][qp] * phi[j][qp];
                            }
                        }
                    }
                    //IF FOR SOME REASON THERE NEEDS TO BE A GROUNDED NODE (VM = 0)
                    if( maxz - tissue_maxz == 0. || minz - tissue_minz == 0. ){
                      for (auto side : elem->side_index_range())
                      {
                          if (elem->neighbor_ptr(side) == nullptr)
                          {
                              auto boundary_id = mesh.get_boundary_info().boundary_id(elem, side);
                              auto it = bc_sidesets.find(static_cast<int>(boundary_id));
                              if (it != bc_sidesets.end())
                              {
                                  //std::cout << "found_boundary " << std::endl;
                                  fe_face->reinit(elem, side);
                                  for (unsigned int qp = 0; qp < qface.n_points(); qp++)
                                  {
                                      for (unsigned int i = 0; i != parabolic_ndofs; i++)
                                      {
                                          for (unsigned int j = 0; j != parabolic_ndofs; j++)
                                          {
                                              K(i, j) += JxW_face[qp] * penalty * phi_face[j][qp] * phi_face[i][qp];
                                          }
                                      }
                                  }
                              }
                          }
                      }
                    }

                }


                /*
                else if(subdomain_id == static_cast<short unsigned int>(Subdomain::FIBROSIS))
                {
                  for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
                  {
                    for (unsigned int i = 0; i != elliptic_ndofs; i++)
                    {
                      for (unsigned int j = 0; j != elliptic_ndofs; j++)
                      {
                        
                        
                      }
                    }
                  }
                }
                */
                //std::cout << ", Loop over interface: " << std::flush;

                // interface
                //if( time_integrator == TimeIntegrator::EXPLICIT_EXTRACELLULAR )
                //{
                //}
                // boundaries
                //std::cout << ", Loop over boundaries: " << std::flush;
                
                //std::cout << "Add matrix " << std::endl;
                parabolic_system.matrix->add_matrix(K, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("M").add_matrix(M, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_vector("ML").add_vector(M_lumped, parabolic_dof_indices);
                //parabolic_system.get_vector("ML2").add_vector(M_lumped2, parabolic_dof_indices);
                //parabolic_system.get_vector("R_dist").add_vector(R_dist, parabolic_dof_indices);
                if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
                  parabolic_system.get_vector("Fibrosis_Check").add_vector(Fibrosis_CheckV, parabolic_dof_indices);
                }
                parabolic_system.get_matrix("Ki").add_matrix(Ki, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Kmono").add_matrix(Kmono, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Ke").add_matrix(Ke, parabolic_dof_indices, parabolic_dof_indices);
                parabolic_system.get_matrix("Kg").add_matrix(Kg, parabolic_dof_indices, parabolic_dof_indices);

                if(checkMatrixV){
                  //libMesh::out << "For element " << *elem << "with size " << parabolic_ndofs << "x" << parabolic_ndofs << " at point within the tissue the matrix is: " << std::endl;
                  //libMesh::out << Kmono << std::endl;
                }

            }
            std::cout << "closing" << std::endl;
            if(fibrosisBool && fibrosisMethod.compare(0,14,"Conductivities") == 0){
              parabolic_system.get_vector("Fibrosis_Check").close();
            }
            parabolic_system.matrix->close();
            parabolic_system.get_matrix("M").close();
            parabolic_system.get_vector("ML").close();
            //parabolic_system.get_vector("ML2").close();
            //parabolic_system.get_vector("R_dist").close();
            parabolic_system.get_matrix("Ki").close();
            parabolic_system.get_matrix("Kmono").close();
            parabolic_system.get_matrix("Ke").close();
            parabolic_system.get_matrix("Kg").close();
            std::cout << "done" << std::endl;

            if(prescribedFibers){
              libMesh::LinearImplicitSystem & fiberSoln = es.get_system < libMesh::LinearImplicitSystem > ("fiberSolution");
              fiberSoln.solution -> close();
              fiberSoln.update();
            }

            break;
        }

    }
    //libMesh::out << "HERE" << std::endl;
}



void assemble_rhs(  libMesh::EquationSystems &es,
                    const TimeData &timedata,
                    TimeIntegrator time_integrator,
                    libMesh::Order p_order,
                    const GetPot &data,
                    EquationType type)
{
    using namespace libMesh;
    double Cm = data("Cm", 1.0);
    double chi = data("chi", 1000.0);
    double sigma_s_i = data("sigma_s_i", 0.2435); // sigma_t in the paper
    double sigma_s_e = data("sigma_s_e", 1.0438);  // sigma_t in the paper
    double sigma_b_ie = data("sigma_b", 6.0);
    double sigma_torso_ie = data("sigma_torso", 6.0);
    bool convergence_test = data("convergence_test", false);
    double xDim = data("maxx", 1.) - data("minx", -1.);
    double v0 = data("v0", 0.);
    double v1 = data("v1", 0.05);
    double v2 = data("v2", 1.);
    double kcubic = data("k", 8.);
    std::string integrator = data("integrator","SBDF1");

    bool monodomainType = data("monodomainType", false);
    bool exportVe_mono = data("exportVe_mono", false);

    double xE1 = data("exportVe_timeX", 2.3172);
    double yE1 = data("exportVe_timeY", 2.3172);
    double zE1 = data("exportVe_timeZ", 2.3172);

    //libMesh::out << "Assembling RHS" << std::endl;
    //std::cout << " assemble_rhs " << std::endl;
    switch (time_integrator)
    {
        case TimeIntegrator::EXPLICIT_EXTRACELLULAR:
        case TimeIntegrator::EXPLICIT_INTRACELLULAR:
        case TimeIntegrator::SEMI_IMPLICIT:
        case TimeIntegrator::SEMI_IMPLICIT_HALF_STEP:
         {
            const MeshBase& mesh = es.get_mesh();
            libMesh::TransientLinearImplicitSystem & parabolic_system = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            const DofMap & parabolic_dof_map = parabolic_system.get_dof_map();
            std::vector < dof_id_type > parabolic_dof_indices;

            if(convergence_test){
                parabolic_system.get_vector("ForcingV").zero();
                ForcingTermConvergence ( es, timedata.dt, chi, Cm, sigma_s_i, sigma_s_e, sigma_b_ie, sigma_torso_ie, timedata.time, xDim, v0, v1, v2, kcubic, integrator);
            }

            if(type == EquationType::PARABOLIC)
            {
                //parabolic_system.get_vector("aux1").zero();
                //parabolic_system.get_vector("R_dist").zero();
                parabolic_system.rhs->zero();
                //parabolic_system.get_vector("R_dist").close();
                if(time_integrator == TimeIntegrator::EXPLICIT_INTRACELLULAR)
                {
                    // Assemble RHS:
                    parabolic_system.rhs->add_vector(*parabolic_system.solution, parabolic_system.get_matrix("Ki"));
                    parabolic_system.rhs->scale(-1.0);
                }
                parabolic_system.get_vector("aux1").zero();
                parabolic_system.get_vector("aux1").add(1.0, parabolic_system.get_vector("In"));
                // add  M * In
                parabolic_system.rhs->add_vector(parabolic_system.get_vector("aux1"), parabolic_system.get_matrix("M"));

                if(convergence_test){
                    parabolic_system.rhs->add(parabolic_system.get_vector("ForcingV"));
                }

                if(time_integrator == TimeIntegrator::SEMI_IMPLICIT ||
                   time_integrator == TimeIntegrator::SEMI_IMPLICIT_HALF_STEP)
                {
                    parabolic_system.get_vector("aux1").zero();
                    parabolic_system.get_vector("aux1").pointwise_mult(*parabolic_system.old_local_solution, parabolic_system.get_vector("ML"));
                    *parabolic_system.rhs += parabolic_system.get_vector("aux1");
                }
                else
                {
                    (*parabolic_system.rhs) /= parabolic_system.get_vector("ML");
                }
                // add forcing



            }
            
            break;
        }
        case TimeIntegrator::SBDF3:
            // Cm * M * ( 11/6 * V^n+1 - 3 * Vn + 3/2 V^n-1 -1/3 V^n-2) / (dt) + Ki * V^n+1 + Ki * Ve^n+1 = -3*I^n + 3*I^n-1 - I^n-2
            // Ki * V^n+1 + Kie * Ve^n+1 = 0            //
            //   RHS = 3 Cm /dt * M * Vn - 3/2 Cm /dt * M * Vnm1 +Cm/dt/3 M Vnm2 - 3 * M In + 3 Inm1 - Inm2

        {
            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("parabolic");
            //system.get_vector("In").print();
            system.rhs->zero();
            system.get_vector("aux1").zero();
            //system.get_vector("aux1").print();

            if(p_order == SECOND)
            {
                // eval: 3 Cm /dt * M * Vn
                system.get_vector("aux1").add(3.0*Cm/timedata.dt, *system.old_local_solution);
                // add: -1.5 Cm /dt * M * Vnm1
                system.get_vector("aux1").add(-1.5*Cm/timedata.dt, *system.older_local_solution);
                // add: +1/3 Cm /dt * M * Vnm2
                system.get_vector("aux1").add(Cm/(3.*timedata.dt), system.get_vector("Vnm2"));
                system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            }
            else
            {
                // eval: 3 Cm /dt * M * Vn
                system.get_vector("aux1").add(3.0, *system.old_local_solution);
                // add: -1.5 Cm /dt * M * Vnm1
                system.get_vector("aux1").add(-1.5, *system.older_local_solution);
                // add: +1/3 Cm /dt * M * Vnm2
                system.get_vector("aux1").add(1./(3.), system.get_vector("Vnm2"));
                system.rhs->pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            }
            // add  M * (2*In-Inm1)
            system.get_vector("aux1").zero();
            system.get_vector("aux1").add(3.0, system.get_vector("In"));
            system.get_vector("aux1").add(-3.0, system.get_vector("Inm1"));
            system.get_vector("aux1").add(1.0, system.get_vector("Inm2"));
            //system.get_vector("aux1").print();

            // add  M * In
            system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            // add forcing
            if(convergence_test){
                system.get_vector("ForcingConv").zero();
                ForcingTermConvergence ( es, timedata.dt, chi, Cm, sigma_s_i, sigma_s_e, sigma_b_ie, sigma_torso_ie, timedata.time, xDim, v0, v1, v2, kcubic, integrator);
                system.rhs->add(system.get_vector("ForcingConv"));
            }
            //system.rhs->add(system.get_vector("F"));



            break;
        }

        case TimeIntegrator::SBDF2:
            // Cm * M * ( 3 * V^n+1 - 4 * Vn + V^n-1 ) / (2*dt) + Ki * V^n+1 + Ki * Ve^n+1 = -2*I^n + I^n-1
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //
            //   RHS = 2 Cm /dt * M * Vn - 1/2 Cm /dt * M * Vnm1 + 2 * M In - Inm1
            // Cm /dt M ( 2 * Vn - 1/2  Vnm1)
        {
            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("parabolic");
            system.rhs->zero();
            system.get_vector("aux1").zero();
            //system.get_vector("In").print();
            if(p_order == SECOND)
            {
                // eval: 3 Cm /dt * M * Vn
                system.get_vector("aux1").add(2.0*Cm/timedata.dt, *system.old_local_solution);
                // add: -1.5 Cm /dt * M * Vnm1
                system.get_vector("aux1").add(-.5*Cm/timedata.dt, *system.older_local_solution);
                system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            }
            else
            {
                // eval: 3 Cm /dt * M * Vn
                system.get_vector("aux1").add(2.0, *system.old_local_solution);
                // add: -1.5 Cm /dt * M * Vnm1
                system.get_vector("aux1").add(-.5, *system.older_local_solution);
                system.rhs->pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            }
            // add  M * (2*In-Inm1)
            system.get_vector("aux1").zero();
            system.get_vector("aux1").add(2.0, system.get_vector("In"));
            system.get_vector("aux1").add(-1.0, system.get_vector("Inm1"));
            //system.get_vector("aux1").print();
            // add  M * In
            system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            //system.get_vector("aux1").pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            //*system.rhs += system.get_vector("aux1");
            // add forcing
            if(convergence_test){
                system.get_vector("ForcingConv").zero();
                ForcingTermConvergence ( es, timedata.dt, chi, Cm, sigma_s_i, sigma_s_e, sigma_b_ie, sigma_torso_ie, timedata.time, xDim, v0, v1, v2, kcubic, integrator);
                system.rhs->add(system.get_vector("ForcingConv"));
            }
            //system.rhs->add(system.get_vector("F"));

           

            break;
        }
        case TimeIntegrator::SBDF1:
            // Cm * M * ( V^n+1 - Vn ) / dt + Ki * V^n+1 + Ki * Ve^n+1 = -I^n
            // Ki * V^n+1 + Kie * Ve^n+1 = 0
            //
            //   RHS = Cm /dt * M * Vn + M In
        default:
        {

            TransientLinearImplicitSystem &system = es.get_system < TransientLinearImplicitSystem > ("parabolic");
            system.rhs->zero();
            system.get_vector("aux1").zero();

            if(p_order == SECOND)
            {
                // eval: Cm /dt * M * Vn
                system.get_vector("aux1").add(Cm/timedata.dt, *system.old_local_solution);
                system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            }
            else
            {
                // eval: Cm /dt * M * Vn
                system.get_vector("aux1").add(1., *system.old_local_solution);
                system.rhs->pointwise_mult(system.get_vector("aux1"), system.get_vector("ML"));
            }
            // add  M * In
            system.get_vector("aux1").zero();
            system.get_vector("aux1").add(1.0, system.get_vector("In"));
            // add  M * In
            system.rhs->add_vector(system.get_vector("aux1"), system.get_matrix("M"));
            // add forcing
            if(convergence_test){
                system.get_vector("ForcingConv").zero();
                ForcingTermConvergence ( es, timedata.dt, chi, Cm, sigma_s_i, sigma_s_e, sigma_b_ie, sigma_torso_ie, timedata.time, xDim, v0, v1, v2, kcubic, integrator);
                system.rhs->add(system.get_vector("ForcingConv"));
            }
            //system.rhs->add(system.get_vector("F"));



            break;
        }
    }

    //libMesh::out << "Done Assembling RHS" << std::endl;
}


// read BC sidesets from string: e.g. bc = "1 2 3 5", or bc = "1, 55, 2, 33"
void read_bc_list(std::string &bc, std::set<int> &bc_sidesets)
{
    std::string number;
    for (auto it = bc.begin(); it != bc.end(); it++)
    {
        auto character = *it;
        if (std::isdigit(character))
            number += character;
        else
        {
            if (number.size() > 0)
            {
                bc_sidesets.insert(std::stoi(number));
                number.clear();
            }
        }
        if (it + 1 == bc.end())
        {
            bc_sidesets.insert(std::stoi(number));
        }
    }
    std::cout << "BC sideset list: " << std::flush;
    for (auto &&sideset : bc_sidesets)
        std::cout << sideset << ", " << std::flush;
    std::cout << std::endl;
}



void SolverRecovery (EquationSystems & es, const GetPot &data, TimeData& datatime)
{
  auto time = datatime.time;
  const Real dt = datatime.dt;
  std::string integrator = data("integrator", "SBDF1");
  double zDim = data("maxz", 0.) - data("minz", 0.);
  double xDim = data("maxx", 0.) - data("minx", 0.);
  int nelx = data("nelx", 40);
  double IstimD = data("stimulus_duration", 2.);
  double IstimV = data("stimulus_amplitude", -1.);
  double tissue_maxx = data("tissue_maxx", .5);
  double tissue_minx = data("tissue_minx", .5);
  double tissue_maxy = data("tissue_maxy", .5);
  double tissue_miny = data("tissue_miny", .5);
  double tissue_maxz = data("tissue_maxz", .5);
  double tissue_minz = data("tissue_minz", .5);
  double stimulus_maxx = data("stimulus_maxx", .5);
  double stimulus_maxy = data("stimulus_maxy", .5);
  double stimulus_minx = data("stimulus_minx", -.5);
  double stimulus_miny = data("stimulus_miny", .85);
  double stimulus_maxz = data("stimulus_maxz", 0.0);
  double stimulus_minz = data("stimulus_minz", 0.0);
  double stimulus_start_time = data("stimulus_start_time", 0.);
  bool convergence_test = data("convergence_test", false);
  bool ICstim = data("ICstim", false);
  bool ICconditions = data("ICconditions", false);
  bool FK_in_mV = data("FK_in_mV", false);
  double u0 = data("v0", 0.);//-85;
  double u1 = data("v1",.05);//-57.6;
  double u2 = data("v2",1.0);//30;
  double kcubic = data("k", 8.);
  double Cm = data("Cm", 1.);
  std::string StimPlace = data("stimulus_type", "Transmembrane");

  std::string StimulusLocation = data("StimulusLocation", "Cube");// Cube or Points
  std::string StimulusLocation_File = data("StimulusLocation_File", "Cube");//If Points, then load file with points to stimulate
  double StimulusLocation_Radius = data("StimulusLocation_Radius", .2);

  
  int SpiralBool = data("SpiralBool", 0);
  double SpiralS2duration = data("SpiralS2duration", 4.);
  double SpiralS2time = data("SpiralS2time", 240.);
  double SpiralS2strength = data("SpiralS2strength", 240.);

  std::string cellModel = data("cellModel", "Fenton-Karma");

  bool monodomainType = data("monodomainType", false);

  //BacNav stuff
  bool boolBacNav = data("boolBacNav", false);
  
  TransientExplicitSystem & Rsystem = es.get_system<TransientExplicitSystem> ("Recovery");
  const DofMap & dof_map2 = Rsystem.get_dof_map();

  MeshBase & mesh = es.get_mesh();

  const unsigned int dim = mesh.mesh_dimension();

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices2;


  vector<double> CoordsPointsX;
  vector<double> CoordsPointsY;
  vector<double> CoordsPointsZ;


  if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
    if(StimulusLocation.compare(0,6,"Points") == 0){
      std::string line, word;
      fstream file (StimulusLocation_File, ios::in);
      //PULLING POINT COORDINATES INTO VECTOR
      if(file.is_open())
      {
      while(getline(file, line))
      {
        std::vector<std::string> CoordsInfoVec;
        stringstream ss(line);     
        while (ss.good()){
            string substr;
            getline(ss, substr, ',');
            CoordsInfoVec.push_back(substr);
        }
        CoordsPointsX.push_back( std::stod(CoordsInfoVec[0]) );
        CoordsPointsY.push_back( std::stod(CoordsInfoVec[1]) );
        CoordsPointsZ.push_back( std::stod(CoordsInfoVec[2]) );
      }
      libMesh::out << "Stimulating at time: " << time << " - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
      }
      else{
      libMesh::out << "Could not open the file..." << std::endl;
      }
    }
    else if(StimulusLocation.compare(0,4,"Cube") == 0){
      libMesh::out << "Stimulating at time: " << time << " - Using the cube with limits: [x: {" << stimulus_minx << "," << stimulus_maxx << "}; y: {" << stimulus_miny << "," << stimulus_maxy << "}; z: {"<< stimulus_minz << "," << stimulus_maxz << "} ]" << std::endl;
    }
  }




  //Fenton-Karma Modified starts here
  if(cellModel.compare(0,12,"Fenton-Karma") == 0){

    const unsigned int v_var = Rsystem.variable_number ("v");
    const unsigned int w_var = Rsystem.variable_number ("w");
    const unsigned int s_var = Rsystem.variable_number ("s");

    double v_new, w_new, s_new, u_old, v_old, w_old, s_old, u_old_old, v_old_old, w_old_old, s_old_old, u_old_old_old, v_old_old_old, w_old_old_old, s_old_old_old;

    double r_new;

    //Bsystem.get_vector("I_ion").zero();
    //Rsystem.solution->zero();

    //Bsystem.get_vector("I_ion").close();

    //libMesh::out << "I AM HERE too INSIDE SOLVER RECOVERY" << std::endl;

    double tau_v_plus = 3.33;
    double tau_v1_minus = 19.2;
    double tau_v2_minus = 10.0;
    double tau_w_plus = 160.0;
    double tau_w1_minus = 75.0;
    double tau_w2_minus = 75.0;
    double tau_d = .065;
    double tau_si = 31.8364;
    double tau_o = 39.0;
    double tau_a = .009;
    double u_c = .23;
    double u_v = .055;
    double u_w = .146;
    double u_o = 0.0;
    double u_m = 1.0;
    double u_csi = .8;
    double u_so = .3;
    double r_s_plus = .02;
    double r_s_minus = 1.2;
    double k = 3.0;
    double a_so = .115;
    double b_so = .84;
    double c_so = .02;

    //Specifically for Original Fenton-Karma
    double g_fi = 4.;
    double tau_r = 50.;
    double tau_0 = 8.3;
    double tau_w_minus = 11.;
    double V_fi = 15.;
    double V_rest = -85.;

    if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
      tau_si = 44.84;
      tau_v_plus = 3.33;
      tau_v1_minus = 1000.;
      tau_v2_minus = 19.2;
      tau_w_plus = 667.;
      u_c = .13;
      u_v = .055;
      u_csi = .85;
      tau_d = Cm / g_fi;
      k = 10.;
    }


    double phi [1000] = {};
    double x1 [1000] = {};
    double y1 [1000] = {};
    double threshSpiral = .25;//(xDim*10./ (double)nelx); //This is just to make sure that the spiral is thick enough

    if(SpiralBool == 1 || SpiralBool == 2){
        //spiral wave change for variables: parameter set 1 Fenton et al 2002
        tau_v_plus = 3.33;
        tau_v1_minus = 19.6;
        tau_v2_minus = 1000.0;
        tau_w_plus = 667.0;
        tau_w1_minus = 11.0;
        tau_w2_minus = 11.0;
        tau_d = .25;
        tau_si = 45.;
        tau_o = 8.3;
        tau_a = .02; //tau_a = 1/tau_r
        k = 10.;
        u_csi = .85;
        u_c = .13;
        u_v = .055;

      }


    
    if(integrator.compare(0,4,"SBDF") == 0){

      TransientLinearImplicitSystem &Bsystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
      Bsystem.update();

    

      const unsigned int u_var = Bsystem.variable_number ("V");
      const DofMap & dof_map = Bsystem.get_dof_map();


      for(const auto & node : mesh.local_node_ptr_range()){

          dof_map.dof_indices (node, dof_indices);
          dof_map2.dof_indices (node, dof_indices2);

          if(dof_indices.size() > 0){

          const Real x = (*node)(0);
          const Real y = (*node)(1);
          const Real z = (*node)(2);

              if(datatime.timestep == 1){
                  v_old = 1.;
                  w_old = 1.;
                  s_old = 0.;
                  u_old = (*Bsystem.current_local_solution)  (dof_indices[u_var]);

                  if(SpiralBool == 1){
                      if( y < 0){
                          v_old = 0.;
                          w_old = 0.;
                          s_old = 1.;
                          if(FK_in_mV){
                            u_old = V_fi;
                          }
                          else{
                            u_old = 1.;
                          }
                      }
                  }

                  else if(SpiralBool == 2){

                      double a = .5;
                      double k = .09;
                      int thickness = 15;
                      int finalDim = thickness*1000; //15,000
                      double bathRegion = xDim/2.0 - tissue_maxx;
                      double factorAxis = 3.0*(xDim/2.-bathRegion);
                      phi[0] = 0.;
                      x1[0] = -a*exp(k*phi[0])*cos(phi[0]);
                      y1[0] = -a*exp(k*phi[0])*sin(phi[0]);
                      for(int i = 1; i < 1000; i++){
                          phi[i] = phi[i-1] + ((factorAxis*acos(-1) - phi[0])/1000.);
                          x1[i] = -a*exp(k*phi[i])*cos(phi[i]);
                          y1[i] = -a*exp(k*phi[i])*sin(phi[i]);
                      }
                      for(int i = 0; i < 1000; i++){
                          if( std::abs(x - x1[i]) <= threshSpiral && std::abs(y - y1[i]) <= threshSpiral){
                              //conditions in notes seen from paraview
                              /*
                               * initial conditions for multiple spirals to happen:
                                  u - spiral of 1. and rest of 0
                                  v - inverse of spiral of 1. and rest of 0
                                  w - spiral of .8 and rest of 1.
                                  s - same spiral of .1 and rest of 0
                              */
                              v_old = 0.;
                              w_old = 0.;
                              s_old = 1.;
                              if(FK_in_mV){
                                u_old = V_fi;
                              }
                              else{
                                u_old = 1.;
                              }
                              //libMesh::out << "HERE --> x = " << x << ";     y = " << y << ";    with threshold of the spiral: " << threshSpiral << std::endl;
                          }
                      }
                  }


              }
              else{
                  v_old = (*Rsystem.current_local_solution) (dof_indices2[v_var]);
                  w_old = (*Rsystem.current_local_solution) (dof_indices2[w_var]);        
                  s_old = (*Rsystem.current_local_solution) (dof_indices2[s_var]);        
                  u_old = (*Bsystem.current_local_solution)  (dof_indices[u_var]);
                  u_old_old = (Bsystem.get_vector("Vnm1"))  (dof_indices[0]);
                  u_old_old_old = (Bsystem.get_vector("Vnm2"))  (dof_indices[0]);
                  v_old_old = (Rsystem.get_vector("v_prev")) (dof_indices2[0]);
                  v_old_old_old = (Rsystem.get_vector("v_prev_prev")) (dof_indices2[0]);
                  w_old_old = (Rsystem.get_vector("w_prev")) (dof_indices2[0]);
                  w_old_old_old = (Rsystem.get_vector("w_prev_prev")) (dof_indices2[0]);
                  s_old_old = (Rsystem.get_vector("s_prev")) (dof_indices2[0]);
                  s_old_old_old = (Rsystem.get_vector("s_prev_prev")) (dof_indices2[0]);
              }


              if(integrator.compare(0,5,"SBDF3") == 0){

                  if(datatime.timestep == 1){
                    if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                        double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                        double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                        v_new = v_old + dt*RhsV;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        w_new = w_old + dt*RhsW;
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                        s_new = s_old + dt*RhsS;
                        Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    }
                    else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                        if(FK_in_mV){
                          u_old = (u_old - V_rest)/(V_fi - V_rest);
                        }

                        double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                        v_new = v_old + dt*RhsV;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        w_new = w_old + dt*RhsW;
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    }

                  }
                  else if(datatime.timestep == 2){
                    if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                        double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                        double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                        double RhsV_old = ( ((1.0 - HofX(u_old_old,u_c))*(1.0 - v_old_old))/(tau_v2_minus*HofX(u_old_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old_old,u_v)))) - ((HofX(u_old_old,u_c)*v_old_old)/(tau_v_plus));
                        double RhsW_old = ( ((1.0 - HofX(u_old_old,u_c))*(1.0 - w_old_old))/(tau_w2_minus*HofX(u_old_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old_old,u_w)))) - ((HofX(u_old_old,u_c)*w_old_old)/(tau_w_plus));
                        double RhsS_old = (r_s_plus*(HofX(u_old_old,u_c)) + r_s_minus*(1.0 - HofX(u_old_old,u_c))) * (.5*(1.0 + tanh(k*(u_old_old - u_csi))) - s_old_old);

                        //v_new = (4.0/3.)*v_old - (1.0/3.)*v_old_old + (4.0/3.)*dt*RhsV -(2.0/3.)*dt*RhsV_old;
                        
                        v_new = v_old + (3./2.)*dt*RhsV - (1./2.)*dt*RhsV_old;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        Rsystem.get_vector("v_prev").set(dof_indices2[0],v_old);
                        w_new = w_old + (3./2.)*dt*RhsW - (1./2.)*dt*RhsW_old;
                        Rsystem.get_vector("w_prev").set(dof_indices2[0],w_old);
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                        s_new = s_old + (3./2.)*dt*RhsS - (1./2.)*dt*RhsS_old;
                        Rsystem.get_vector("s_prev").set(dof_indices2[0],s_old);
                        Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    }
                    else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                        if(FK_in_mV){
                          u_old = (u_old - V_rest)/(V_fi - V_rest);
                          u_old_old = (u_old_old - V_rest)/(V_fi - V_rest);
                        }

                        double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                        double RhsV_old = ( (HofX(u_c,u_old_old)*(1.0 - v_old_old))/( tau_v2_minus*HofX(u_v,u_old_old) + tau_v1_minus*HofX(u_old_old,u_v) )) - ((HofX(u_old_old,u_c)*v_old_old)/(tau_v_plus));
                        double RhsW_old = ( (HofX(u_c,u_old_old)*(1.0 - w_old_old))/( tau_w_minus )) - ((HofX(u_old_old,u_c)*w_old_old)/(tau_w_plus));

                        v_new = v_old + (3./2.)*dt*RhsV - (1./2.)*dt*RhsV_old;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        Rsystem.get_vector("v_prev").set(dof_indices2[0],v_old);
                        w_new = w_old + (3./2.)*dt*RhsW - (1./2.)*dt*RhsW_old;
                        Rsystem.get_vector("w_prev").set(dof_indices2[0],w_old);
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    }
                    /*
                    v_new = v_old + dt*RhsV;
                    Rsystem.solution -> set(dof_indices2[v_var],v_new);
                    w_new = w_old + dt*RhsW;
                    Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    s_new = s_old + dt*RhsS;
                    Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    */
                  }
                  else{
                    if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                        double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                        double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                        double RhsV_old = ( ((1.0 - HofX(u_old_old,u_c))*(1.0 - v_old_old))/(tau_v2_minus*HofX(u_old_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old_old,u_v)))) - ((HofX(u_old_old,u_c)*v_old_old)/(tau_v_plus));
                        double RhsW_old = ( ((1.0 - HofX(u_old_old,u_c))*(1.0 - w_old_old))/(tau_w2_minus*HofX(u_old_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old_old,u_w)))) - ((HofX(u_old_old,u_c)*w_old_old)/(tau_w_plus));
                        double RhsS_old = (r_s_plus*(HofX(u_old_old,u_c)) + r_s_minus*(1.0 - HofX(u_old_old,u_c))) * (.5*(1.0 + tanh(k*(u_old_old - u_csi))) - s_old_old);

                        double RhsV_old_old = ( ((1.0 - HofX(u_old_old_old,u_c))*(1.0 - v_old_old_old))/(tau_v2_minus*HofX(u_old_old_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old_old_old,u_v)))) - ((HofX(u_old_old_old,u_c)*v_old_old_old)/(tau_v_plus));
                        double RhsW_old_old = ( ((1.0 - HofX(u_old_old_old,u_c))*(1.0 - w_old_old_old))/(tau_w2_minus*HofX(u_old_old_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old_old_old,u_w)))) - ((HofX(u_old_old_old,u_c)*w_old_old_old)/(tau_w_plus));
                        double RhsS_old_old = (r_s_plus*(HofX(u_old_old_old,u_c)) + r_s_minus*(1.0 - HofX(u_old_old_old,u_c))) * (.5*(1.0 + tanh(k*(u_old_old_old - u_csi))) - s_old_old_old);

                        
                        v_new = v_old + (23./12.)*dt*RhsV - (16./12.)*dt*RhsV_old + (5./12.)*dt*RhsV_old_old;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        Rsystem.get_vector("v_prev").set(dof_indices2[0],v_old);
                        Rsystem.get_vector("v_prev_prev").set(dof_indices2[0],v_old_old);
                        w_new = w_old + (23./12.)*dt*RhsW - (16./12.)*dt*RhsW_old + (5./12.)*dt*RhsW_old_old;
                        Rsystem.get_vector("w_prev").set(dof_indices2[0],w_old);
                        Rsystem.get_vector("w_prev_prev").set(dof_indices2[0],w_old_old);
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                        s_new = s_old + (23./12.)*dt*RhsS - (16./12.)*dt*RhsS_old + (5./12.)*dt*RhsS_old_old;
                        Rsystem.get_vector("s_prev").set(dof_indices2[0],s_old);
                        Rsystem.get_vector("s_prev_prev").set(dof_indices2[0],s_old_old);
                        Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    }
                    else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                        if(FK_in_mV){
                          u_old = (u_old - V_rest)/(V_fi - V_rest);
                          u_old_old = (u_old_old - V_rest)/(V_fi - V_rest);
                          u_old_old_old = (u_old_old_old - V_rest)/(V_fi - V_rest);
                        }

                        double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                        double RhsV_old = ( (HofX(u_c,u_old_old)*(1.0 - v_old_old))/( tau_v2_minus*HofX(u_v,u_old_old) + tau_v1_minus*HofX(u_old_old,u_v) )) - ((HofX(u_old_old,u_c)*v_old_old)/(tau_v_plus));
                        double RhsW_old = ( (HofX(u_c,u_old_old)*(1.0 - w_old_old))/( tau_w_minus )) - ((HofX(u_old_old,u_c)*w_old_old)/(tau_w_plus));

                        double RhsV_old_old = ( (HofX(u_c,u_old_old_old)*(1.0 - v_old_old_old))/( tau_v2_minus*HofX(u_v,u_old_old_old) + tau_v1_minus*HofX(u_old_old_old,u_v) )) - ((HofX(u_old_old_old,u_c)*v_old_old_old)/(tau_v_plus));
                        double RhsW_old_old = ( (HofX(u_c,u_old_old_old)*(1.0 - w_old_old_old))/( tau_w_minus )) - ((HofX(u_old_old_old,u_c)*w_old_old_old)/(tau_w_plus));

                        v_new = v_old + (23./12.)*dt*RhsV - (16./12.)*dt*RhsV_old + (5./12.)*dt*RhsV_old_old;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        Rsystem.get_vector("v_prev").set(dof_indices2[0],v_old);
                        Rsystem.get_vector("v_prev_prev").set(dof_indices2[0],v_old_old);
                        w_new = w_old + (23./12.)*dt*RhsW - (16./12.)*dt*RhsW_old + (5./12.)*dt*RhsW_old_old;
                        Rsystem.get_vector("w_prev").set(dof_indices2[0],w_old);
                        Rsystem.get_vector("w_prev_prev").set(dof_indices2[0],w_old_old);
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    }
                    

                    /*
                    v_new = v_old + dt*RhsV;
                    Rsystem.solution -> set(dof_indices2[v_var],v_new);
                    w_new = w_old + dt*RhsW;
                    Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    s_new = s_old + dt*RhsS;
                    Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    */
                  }

              }
              else if(integrator.compare(0,5,"SBDF2") == 0){

                  if(datatime.timestep == 1){
                    if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                        double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                        double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                        v_new = v_old + dt*RhsV;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        w_new = w_old + dt*RhsW;
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                        s_new = s_old + dt*RhsS;
                        Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    }
                    else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                        if(FK_in_mV){
                          u_old = (u_old - V_rest)/(V_fi - V_rest);
                        }

                        double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                        v_new = v_old + dt*RhsV;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        w_new = w_old + dt*RhsW;
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    }
                  }
                  else{
                    if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                        double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                        double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                        double RhsV_old = ( ((1.0 - HofX(u_old_old,u_c))*(1.0 - v_old_old))/(tau_v2_minus*HofX(u_old_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old_old,u_v)))) - ((HofX(u_old_old,u_c)*v_old_old)/(tau_v_plus));
                        double RhsW_old = ( ((1.0 - HofX(u_old_old,u_c))*(1.0 - w_old_old))/(tau_w2_minus*HofX(u_old_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old_old,u_w)))) - ((HofX(u_old_old,u_c)*w_old_old)/(tau_w_plus));
                        double RhsS_old = (r_s_plus*(HofX(u_old_old,u_c)) + r_s_minus*(1.0 - HofX(u_old_old,u_c))) * (.5*(1.0 + tanh(k*(u_old_old - u_csi))) - s_old_old);

                        //v_new = (4.0/3.)*v_old - (1.0/3.)*v_old_old + (4.0/3.)*dt*RhsV -(2.0/3.)*dt*RhsV_old;
                        
                        v_new = v_old + (3./2.)*dt*RhsV - (1./2.)*dt*RhsV_old;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        Rsystem.get_vector("v_prev").set(dof_indices2[0],v_old);
                        w_new = w_old + (3./2.)*dt*RhsW - (1./2.)*dt*RhsW_old;
                        Rsystem.get_vector("w_prev").set(dof_indices2[0],w_old);
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                        s_new = s_old + (3./2.)*dt*RhsS - (1./2.)*dt*RhsS_old;
                        Rsystem.get_vector("s_prev").set(dof_indices2[0],s_old);
                        Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    }
                    else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                        if(FK_in_mV){
                          u_old = (u_old - V_rest)/(V_fi - V_rest);
                          u_old_old = (u_old_old - V_rest)/(V_fi - V_rest);
                        }

                        double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                        double RhsV_old = ( (HofX(u_c,u_old_old)*(1.0 - v_old_old))/( tau_v2_minus*HofX(u_v,u_old_old) + tau_v1_minus*HofX(u_old_old,u_v) )) - ((HofX(u_old_old,u_c)*v_old_old)/(tau_v_plus));
                        double RhsW_old = ( (HofX(u_c,u_old_old)*(1.0 - w_old_old))/( tau_w_minus )) - ((HofX(u_old_old,u_c)*w_old_old)/(tau_w_plus));

                        v_new = v_old + (3./2.)*dt*RhsV - (1./2.)*dt*RhsV_old;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        Rsystem.get_vector("v_prev").set(dof_indices2[0],v_old);
                        w_new = w_old + (3./2.)*dt*RhsW - (1./2.)*dt*RhsW_old;
                        Rsystem.get_vector("w_prev").set(dof_indices2[0],w_old);
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    }
                    

                    /*
                    v_new = v_old + dt*RhsV;
                    Rsystem.solution -> set(dof_indices2[v_var],v_new);
                    w_new = w_old + dt*RhsW;
                    Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    s_new = s_old + dt*RhsS;
                    Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    */

                  }

              }
              else{
                    if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                        double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                        double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                        v_new = v_old + dt*RhsV;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        w_new = w_old + dt*RhsW;
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                        s_new = s_old + dt*RhsS;
                        Rsystem.solution -> set(dof_indices2[s_var],s_new);
                    }
                    else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                        if(FK_in_mV){
                          u_old = (u_old - V_rest)/(V_fi - V_rest);
                        }

                        double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                        double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                        v_new = v_old + dt*RhsV;
                        Rsystem.solution -> set(dof_indices2[v_var],v_new);
                        w_new = w_old + dt*RhsW;
                        Rsystem.solution -> set(dof_indices2[w_var],w_new);
                    }
              }

              double freact = 0.0;

              double Ifival = (-v_old*HofX(u_old,u_c)*(u_old-u_c)*(u_m-u_old)) / (tau_d);
              double Isival, Isoval;
              if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                Isival = (-w_old*s_old) / (tau_si);
                Isoval = (((u_old-u_o)*(1.0 - HofX(u_old,u_so))) / tau_o) + HofX(u_old,u_so)*tau_a + .5*(a_so-tau_a)*(1.0 + tanh((u_old - b_so)/(c_so)));
              }
              else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                Isival = (-w_old / (2.*tau_si))*( 1. + tanh(k*(u_old - u_csi)) );
                Isoval = (u_old/tau_0)*(HofX(u_c, u_old)) + (1./tau_r)*(HofX(u_old,u_c));
              }
              double Istim = 0.0;
              double Istim2 = 0.0;

              double IonCurrents;

              if(FK_in_mV){
                IonCurrents = (Ifival + Isival + Isoval)*( Cm*( V_fi - V_rest ) );
              }
              else{
                IonCurrents = Ifival + Isival + Isoval;
              }

              if(SpiralBool == 0){
                  if( zDim == 0. ){
                    if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){

                      if(StimulusLocation.compare(0,4,"Cube") == 0){
                        if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx ){
                          Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                          //libMesh::out << x << std::endl;
                        }
                      }
                      else if(StimulusLocation.compare(0,6,"Points") == 0){
                        //READ POINTS
                        //LOOP OVER POINTS
                        //SET STIMULUS

                        for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                         if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius ){
                            Istim = IstimV;
                            //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                            break;
                         }
                         else{}
                        }

                      }

                    }
                  }
                  else{
                    if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){


                      if(StimulusLocation.compare(0,4,"Cube") == 0){
                        if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx && z > stimulus_minz && z < stimulus_maxz ){
                          Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                          //libMesh::out << x << std::endl;
                        }
                      }
                      else if(StimulusLocation.compare(0,6,"Points") == 0){
                        //READ POINTS
                        //LOOP OVER POINTS
                        //SET STIMULUS
                        
                        for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                         if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius && std::abs(z - CoordsPointsZ[iter]) < StimulusLocation_Radius ){
                            Istim = IstimV;
                            //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                            break;
                         }
                         else{}
                        }

                      }


                    }
                  }
                  if(convergence_test){
                    freact = -kcubic*((u_old - u0)*(u_old - u1)*(u_old - u2));// - 1*r_new*(u_old-u0)) - Istim  - Istim2;
                  }
                  else{
                      if(StimPlace.compare(0,13,"Transmembrane") == 0){
                          if(ICstim){
                            freact = (-(IonCurrents));
                          }
                          else if(ICconditions){
                            freact = (-(IonCurrents));
                          }
                          else{
                            freact = (-(IonCurrents)) - Istim  - Istim2;
                          }
                          
                      }
                      else{
                          freact = (-(IonCurrents));
                      }
                  }
              }

              else if(SpiralBool == 1){
                  if(time < IstimD && x < 0. && y > 0.){
                      Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                  //libMesh::out << Istim << std::endl;
                  }
                  freact = (-(IonCurrents)) - Istim  - Istim2;
              }

              else if(SpiralBool == 2){
                  if(time < IstimD &&  y > tissue_maxy - 0.25 ){
                  //if(time < IstimD &&  x < tissue_minx + 0.25 && y > 0 ){
                      Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  }
                  if(time < IstimD &&  x > tissue_maxx - 0.15 ){
                      Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  }
                  //SECOND STIMULUS
                  //if(time > 150. && time < 152. &&  x > tissue_maxx - 0.15 && y < 0. ){
                      //Istim = IstimV*.5;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  //}

                  freact = (-(IonCurrents)) - Istim  - Istim2;
              }



              
              //libMesh::out << freact << "    " << x << "    " << dof_indices[0] << std::endl;
              //libMesh::out << x << std::endl;


              Bsystem.get_vector("In").set(dof_indices[u_var], freact);


           
          }

          

      }
      Bsystem.get_vector("In").close();
      Bsystem.get_vector("Inm1").close();

      Rsystem.solution -> close();
      Rsystem.update(); 

    //Rsystem.get_vector("Rvalues") = Rnew;
    

    }




    if(integrator.compare(0,4,"SBDF") != 0){

      TransientLinearImplicitSystem &Bsystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
      Bsystem.update();
   
      const unsigned int u_var = Bsystem.variable_number ("V");
      const DofMap & dof_map = Bsystem.get_dof_map();


      for(const auto & node : mesh.local_node_ptr_range()){

          dof_map.dof_indices (node, dof_indices);
          dof_map2.dof_indices (node, dof_indices2);

          if(dof_indices.size() > 0){

          const Real x = (*node)(0);
          const Real y = (*node)(1);
          const Real z = (*node)(2);

            if(datatime.timestep == 1){
              v_old = 1.;
              w_old = 1.;
              s_old = 0.;
              u_old = (*Bsystem.current_local_solution)  (dof_indices[u_var]);

              if(SpiralBool == 1){
                  if( y < 0){
                      v_old = 0.;
                      w_old = 0.;
                      s_old = 1.;
                      if(FK_in_mV){
                        u_old = V_fi;
                      }
                      else{
                        u_old = 1.;
                      }
                  }
              }

              else if(SpiralBool == 2){

                  double a = .5;
                  double k = .09;
                  int thickness = 15;
                  int finalDim = thickness*1000; //15,000
                  double bathRegion = xDim/2.0 - tissue_maxx;
                  double factorAxis = 3.0*(xDim/2.-bathRegion);
                  phi[0] = 0.;
                  x1[0] = -a*exp(k*phi[0])*cos(phi[0]);
                  y1[0] = -a*exp(k*phi[0])*sin(phi[0]);
                  for(int i = 1; i < 1000; i++){
                      phi[i] = phi[i-1] + ((factorAxis*acos(-1) - phi[0])/1000.);
                      x1[i] = -a*exp(k*phi[i])*cos(phi[i]);
                      y1[i] = -a*exp(k*phi[i])*sin(phi[i]);
                  }
                  for(int i = 0; i < 1000; i++){
                      if( std::abs(x - x1[i]) <= threshSpiral && std::abs(y - y1[i]) <= threshSpiral){
                          //conditions in notes seen from paraview
                          /*
                           * initial conditions for multiple spirals to happen:
                              u - spiral of 1. and rest of 0
                              v - inverse of spiral of 1. and rest of 0
                              w - spiral of .8 and rest of 1.
                              s - same spiral of .1 and rest of 0
                          */
                          v_old = 0.;
                          w_old = 0.;
                          s_old = 1.;
                          if(FK_in_mV){
                            u_old = V_fi;
                          }
                          else{
                            u_old = 1.;
                          }
                          //libMesh::out << "HERE --> x = " << x << ";     y = " << y << ";    with threshold of the spiral: " << threshSpiral << std::endl;
                      }
                  }
              }

            }
            else{
              v_old = (*Rsystem.current_local_solution) (dof_indices2[v_var]);
              w_old = (*Rsystem.current_local_solution) (dof_indices2[w_var]);        
              s_old = (*Rsystem.current_local_solution) (dof_indices2[s_var]);        
              u_old = (*Bsystem.current_local_solution)  (dof_indices[u_var]);
            }

            if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
                double RhsV = ( ((1.0 - HofX(u_old,u_c))*(1.0 - v_old))/(tau_v2_minus*HofX(u_old,u_v) + tau_v1_minus*(1.0 - HofX(u_old,u_v)))) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                double RhsW = ( ((1.0 - HofX(u_old,u_c))*(1.0 - w_old))/(tau_w2_minus*HofX(u_old,u_w) + tau_w1_minus*(1.0 - HofX(u_old,u_w)))) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));
                double RhsS = (r_s_plus*(HofX(u_old,u_c)) + r_s_minus*(1.0 - HofX(u_old,u_c))) * (.5*(1.0 + tanh(k*(u_old - u_csi))) - s_old);

                v_new = v_old + dt*RhsV;
                Rsystem.solution -> set(dof_indices2[v_var],v_new);
                w_new = w_old + dt*RhsW;
                Rsystem.solution -> set(dof_indices2[w_var],w_new);
                s_new = s_old + dt*RhsS;
                Rsystem.solution -> set(dof_indices2[s_var],s_new);
            }
            else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                if(FK_in_mV){
                  u_old = (u_old - V_rest)/(V_fi - V_rest);
                }

                double RhsV = ( (HofX(u_c,u_old)*(1.0 - v_old))/( tau_v2_minus*HofX(u_v,u_old) + tau_v1_minus*HofX(u_old,u_v) )) - ((HofX(u_old,u_c)*v_old)/(tau_v_plus));
                double RhsW = ( (HofX(u_c,u_old)*(1.0 - w_old))/( tau_w_minus )) - ((HofX(u_old,u_c)*w_old)/(tau_w_plus));

                v_new = v_old + dt*RhsV;
                Rsystem.solution -> set(dof_indices2[v_var],v_new);
                w_new = w_old + dt*RhsW;
                Rsystem.solution -> set(dof_indices2[w_var],w_new);
            }
      

            double freact = 0.0;
            double Ifival = (-v_old*HofX(u_old,u_c)*(u_old-u_c)*(u_m-u_old)) / (tau_d);
            double Isival, Isoval;
            if(cellModel.compare(0,15,"Fenton-KarmaMod") == 0){
              Isival = (-w_old*s_old) / (tau_si);
              Isoval = (((u_old-u_o)*(1.0 - HofX(u_old,u_so))) / tau_o) + HofX(u_old,u_so)*tau_a + .5*(a_so-tau_a)*(1.0 + tanh((u_old - b_so)/(c_so)));
            }
            else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
              Isival = (-w_old / (2.*tau_si))*( 1. + tanh(k*(u_old - u_csi)) );
              Isoval = (u_old/tau_0)*(HofX(u_c, u_old)) + (1./tau_r)*(HofX(u_old,u_c));
            }
            double Istim = 0.0;
            double Istim2 = 0.0;

            double IonCurrents;

            if(FK_in_mV){
              IonCurrents = (Ifival + Isival + Isoval)*( Cm*( V_fi - V_rest ) );
            }
            else{
              IonCurrents = Ifival + Isival + Isoval;
            }

            if(SpiralBool == 0){
                  if( zDim == 0. ){
                    if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
                      

                      if(StimulusLocation.compare(0,4,"Cube") == 0){
                        if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx ){
                          Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                          //libMesh::out << x << std::endl;
                        }
                      }
                      else if(StimulusLocation.compare(0,6,"Points") == 0){
                        //READ POINTS
                        //LOOP OVER POINTS
                        //SET STIMULUS
                        
                        for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                         if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius ){
                            Istim = IstimV;
                            //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                            break;
                         }
                         else{}
                        }

                      }


                    }
                  }
                  else{
                    if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
                      


                      if(StimulusLocation.compare(0,4,"Cube") == 0){
                        if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx && z > stimulus_minz && z < stimulus_maxz ){
                          Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                          //libMesh::out << x << std::endl;
                        }
                      }
                      else if(StimulusLocation.compare(0,6,"Points") == 0){
                        //READ POINTS
                        //LOOP OVER POINTS
                        //SET STIMULUS
                        
                        for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                         if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius && std::abs(z - CoordsPointsZ[iter]) < StimulusLocation_Radius ){
                            Istim = IstimV;
                            //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                            break;
                         }
                         else{}
                        }

                      }



                    }
                  }
                  if(convergence_test){
                    freact = -kcubic*((u_old - u0)*(u_old - u1)*(u_old - u2));// - 1*r_new*(u_old-u0)) - Istim  - Istim2;
                  }
                  else{
                      if(StimPlace.compare(0,13,"Transmembrane") == 0){
                          if(ICstim){
                            freact = (-(IonCurrents));
                          }
                          else if(ICconditions){
                            freact = (-(IonCurrents));
                          }
                          else{
                            freact = (-(IonCurrents)) - Istim  - Istim2;
                          }
                      }
                      else{
                          freact = (-(IonCurrents));
                      }
                  }
              }


              else if(SpiralBool == 1){
                  if(time < IstimD && x < 0. && y > 0.){
                      Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                  //libMesh::out << Istim << std::endl;
                  }
                  freact = (-(IonCurrents)) - Istim  - Istim2;
              }

              else if(SpiralBool == 2){
                  //if(time < IstimD &&  x < tissue_minx + 0.5 ){
                      //Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  //}
                  if(time < IstimD &&  y > tissue_maxy - 0.25 ){
                  //if(time < IstimD &&  x < tissue_minx + 0.25 && y > 0 ){
                      Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  }
                  if(time < IstimD &&  x > tissue_maxx - 0.15 ){
                      Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  }
                  //SECOND STIMULUS
                  //if(time > 150. && time < 152. &&  x > tissue_maxx - 0.15 && y < 0. ){
                      //Istim = IstimV*.5;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                      //libMesh::out << Istim << std::endl;
                  //}

                  freact = (-(IonCurrents)) - Istim  - Istim2;
              }

            //libMesh::out << freact << "    " << x << "    " << dof_indices[0] << std::endl;
            //libMesh::out << x << std::endl;


            Bsystem.get_vector("In").set(dof_indices[u_var], freact);

          
        }


      }
      Bsystem.get_vector("In").close();
      //Rsystem.get_vector("Rvalues") = Rnew;

      Rsystem.solution -> close();
      Rsystem.update();


    
      
    }
  }
  //end of Fenton-Karma Modified

  //start of Courtemanche

  else if(cellModel.compare(0,12,"Courtemanche") == 0){
        
      const unsigned int v_var = Rsystem.variable_number ("v");
      const unsigned int w_var = Rsystem.variable_number ("w");
      const unsigned int xs_var = Rsystem.variable_number ("xs");
      const unsigned int u_var = Rsystem.variable_number ("u");
      const unsigned int m_var = Rsystem.variable_number ("m");
      const unsigned int h_var = Rsystem.variable_number ("h");
      const unsigned int d_var = Rsystem.variable_number ("d");
      const unsigned int xr_var = Rsystem.variable_number ("xr");
      const unsigned int Nai_var = Rsystem.variable_number ("Nai");
      const unsigned int Ki_var = Rsystem.variable_number ("Ki");
      const unsigned int Carel_var = Rsystem.variable_number ("Carel");
      const unsigned int oa_var = Rsystem.variable_number ("oa");
      const unsigned int ui_var = Rsystem.variable_number ("ui");
      const unsigned int oi_var = Rsystem.variable_number ("oi");
      const unsigned int f_var = Rsystem.variable_number ("f");
      const unsigned int j_var = Rsystem.variable_number ("j");
      const unsigned int Cai_var = Rsystem.variable_number ("Cai");
      const unsigned int Caup_var = Rsystem.variable_number ("Caup");
      const unsigned int ua_var = Rsystem.variable_number ("ua");
      const unsigned int fCa_var = Rsystem.variable_number ("fCa");

      const unsigned int mBN_var = Rsystem.variable_number ("mBN"); //BacNav m
      const unsigned int hBN_var = Rsystem.variable_number ("hBN"); //BacNav m

      double v_new, u_new, w_new, xs_new, vm_old, v_old, w_old, xs_old, vm_old_old, v_old_old, w_old_old, xs_old_old, vm_old_old_old, v_old_old_old, w_old_old_old, xs_old_old_old, u_old, u_old_old, u_old_old_old;
      double h_new, d_new, j_new, f_new, h_old, d_old, j_old, f_old, h_old_old, d_old_old, j_old_old, f_old_old, h_old_old_old, d_old_old_old, j_old_old_old, f_old_old_old;
      double xr_new, Nai_new, Ki_new, Carel_new, xr_old, Nai_old, Ki_old, Carel_old, xr_old_old, Nai_old_old, Ki_old_old, Carel_old_old, xr_old_old_old, Nai_old_old_old, Ki_old_old_old, Carel_old_old_old;
      double oi_new, ui_new, m_new, Cai_new, oi_old, ui_old, m_old, Cai_old, oi_old_old, ui_old_old, m_old_old, Cai_old_old, oi_old_old_old, ui_old_old_old, m_old_old_old, Cai_old_old_old;
      double Caup_new, oa_new, ua_new, fCa_new, Caup_old, oa_old, ua_old, fCa_old, Caup_old_old, oa_old_old, ua_old_old, fCa_old_old, Caup_old_old_old, oa_old_old_old, ua_old_old_old, fCa_old_old_old;
      double hBN_new, hBN_old, hBN_old_old, hBN_old_old_old, mBN_new, mBN_old, mBN_old_old, mBN_old_old_old; 
      //Bsystem.get_vector("I_ion").zero();
      //Rsystem.solution->zero();
      //Bsystem.get_vector("I_ion").close();
      //libMesh::out << "I AM HERE too INSIDE SOLVER RECOVERY" << std::endl;
      double R = 8.3143; //J/(K*mol) gas constant
      double T = 310.; //K normal body temperature
      double F = 96.4867; //C/mmol Faraday's constant
      //double Cm = 1.; //pF membrane capacitance
      double Vol_cell = 20100.; //um^3  cell volume
      double Vol_i = 13668.; //um^3 intracellular volume
      double Vol_up = 1109.52; //um^3 Sarcoplasmic Reticulum (SR) uptake compartment volume
      double Vol_rel = 96.48; //um^3 SR release compartment volume
      double K_o = 5.4; //mM extracellular potassium concentration
      double Na_o = 140.; //mM extracellular sodium concentration
      double Ca_o = 1.8; //mM extracellular calcium concentration
      double g_Na = 7.8; //nS/pF max I_Na conductance
      double g_K1 = .09; //nS/pF max I_K1 conductance
      double g_to = .1652; //nS/pF max I_to conductance
      double g_Kr = .0294; //nS/pF max I_Kr conductance
      double g_Ks = .129; //nS/pF max I_Ks conductance
      double g_CaL = .1238; //nS/pF max I_CaL conductance
      double g_bCa = .00113; //nS/pF max I_bCa conductance
      double g_bNa = .000674; //nS/pF max I_bNa conductance
      double I_NaK_max = .60; //pA/pF max I_NaK
      double I_NaCa_max = 1600.; //pA/pF max I_NaCa
      double I_pCa_max = .275; //pA/pF max I_pCa calcium pump current
      double I_up_max = .005; //mM/ms max I_up
      double K_Q10 = 3.; //Temperature scaling factor for I_Kur and I_to kinetics
      double gamma = .35; //Voltage-dependence parameter for I_NaCa
      double K_m_Na_i = 10.; //mM  Na_intra concentration half-saturation constant for I_NaK
      double K_m_K_o = 1.5; //mM  K_extra concentration half-saturation constant for I_NaK
      double K_m_Na = 87.5; //mM   Na_extra concentration half-saturation constant for I_NaCa
      double K_m_Ca = 1.38; //mM   Ca_extra concentration half-saturation constant for I_NaCa
      double k_sat = .1; //saturation factor for I_NaCa
      double k_rel = 30.; //ms^-1    maximal release rate for I_rel
      double K_up = .00092; //mM   Ca_intra concentration half-saturation constant for I_up
      double Ca_up_max = 15.; //mM    max Ca concentration in uptake compartment
      double Cmdn_max = .05; //mM   total Calmodulin concentration in myoplasm
      double Trpn_max = .07; //mM   total Troponin concentration in myoplasm
      double Csqn_max = 10.; //mM   total Calsequestrin concentration in SR release compartment
      double Km_Cmdn = .00238; //mM  Ca_intra concentration half-saturation constant for Calmodulin
      double Km_Trpn = .0005; //mM  Ca_intra concentration half-saturation constant for Troponin
      double Km_Csqn = .8; //mM  Ca_rel concentration half-saturation constant for I_up
      double cell_length = 100.; //um
      double cell_diameter = 16.; //um
      double tau_tr = 180.;
      double tau_fCa = 2.;
      double sigma = (1./7.)*(exp(Na_o/67.3) - 1.);
      double I_bK = 0.;
      double tau_u = 8.;

      double decreaseFacgNa = data("decreaseFacgNa", 1.); //factor to decrease conductance of Na, 1 if normal conductance
      double increaseFacgto = data("increaseFacgto", 1.); //factor to increase conductance of transient outwards, 1 if normal conductance
      double rateOfInact = data("rateOfInact", 1.);//factor to control speed of fast inactivation for INa, 1 if normal speed
      double percBacNav = data("percBacNav", 0.);//factor to control how much of BacNav current will be added, 0 if no BacNav at all
      double g_Na_BacNav = data("g_Na_BacNav", 10.); //nS/pF max I_Na_BacNav conductance   (21.6399 is original)  instead of 10 to affect peak of notch in AP (the higher it is, the higher the peak is)
      double tauhParam = data("tauhParam", 90.);//parameter within calculation of tau_h_BN to affect AP shape, the higher it is the shorter the APD is... original is 84.6609 instead of 90.
      double hBN_IC = 0.8231; //IC
      double mBN_IC = 0.0000094; //IC

      double E_Ca, E_Na, E_K, I_Na, I_K1, I_to, I_Kur, I_Kr, I_Ks, I_CaL, I_pCa, I_NaK, I_NaCa, I_bNa, I_bCa, I_rel, I_tr, I_up, I_up_leak;
      double alpha_m, alpha_h, alpha_j, alpha_oa, alpha_oi, alpha_ua, alpha_ui, alpha_xr, alpha_xs;
      double beta_m, beta_h, beta_j, beta_oa, beta_oi, beta_ua, beta_ui, beta_xr, beta_xs;
      double tau_d, tau_f, tau_v, tau_w, tau_m, tau_h, tau_j, tau_oa, tau_oi, tau_ua, tau_ui, tau_xr, tau_xs;
      double u_inf, fCa_inf, d_inf, f_inf, v_inf, w_inf, m_inf, h_inf, j_inf, oa_inf, oi_inf, ua_inf, ui_inf, xr_inf, xs_inf;
      double g_Kur, f_NaK, F_n, B1, B2;

      double E_Ca_old, E_Na_old, E_K_old, I_Na_old, I_K1_old, I_to_old, I_Kur_old, I_Kr_old, I_Ks_old, I_CaL_old, I_pCa_old, I_NaK_old, I_NaCa_old, I_bNa_old, I_bCa_old, I_rel_old, I_tr_old, I_up_old, I_up_leak_old;
      double alpha_m_old, alpha_h_old, alpha_j_old, alpha_oa_old, alpha_oi_old, alpha_ua_old, alpha_ui_old, alpha_xr_old, alpha_xs_old;
      double beta_m_old, beta_h_old, beta_j_old, beta_oa_old, beta_oi_old, beta_ua_old, beta_ui_old, beta_xr_old, beta_xs_old;
      double tau_d_old, tau_f_old, tau_v_old, tau_w_old, tau_m_old, tau_h_old, tau_j_old, tau_oa_old, tau_oi_old, tau_ua_old, tau_ui_old, tau_xr_old, tau_xs_old;
      double u_inf_old, fCa_inf_old, d_inf_old, f_inf_old, v_inf_old, w_inf_old, m_inf_old, h_inf_old, j_inf_old, oa_inf_old, oi_inf_old, ua_inf_old, ui_inf_old, xr_inf_old, xs_inf_old;
      double g_Kur_old, f_NaK_old, F_n_old, B1_old, B2_old;

      double E_Ca_old_old, E_Na_old_old, E_K_old_old, I_Na_old_old, I_K1_old_old, I_to_old_old, I_Kur_old_old, I_Kr_old_old, I_Ks_old_old, I_CaL_old_old, I_pCa_old_old, I_NaK_old_old, I_NaCa_old_old, I_bNa_old_old, I_bCa_old_old, I_rel_old_old, I_tr_old_old, I_up_old_old, I_up_leak_old_old;
      double alpha_m_old_old, alpha_h_old_old, alpha_j_old_old, alpha_oa_old_old, alpha_oi_old_old, alpha_ua_old_old, alpha_ui_old_old, alpha_xr_old_old, alpha_xs_old_old;
      double beta_m_old_old, beta_h_old_old, beta_j_old_old, beta_oa_old_old, beta_oi_old_old, beta_ua_old_old, beta_ui_old_old, beta_xr_old_old, beta_xs_old_old;
      double tau_d_old_old, tau_f_old_old, tau_v_old_old, tau_w_old_old, tau_m_old_old, tau_h_old_old, tau_j_old_old, tau_oa_old_old, tau_oi_old_old, tau_ua_old_old, tau_ui_old_old, tau_xr_old_old, tau_xs_old_old;
      double u_inf_old_old, fCa_inf_old_old, d_inf_old_old, f_inf_old_old, v_inf_old_old, w_inf_old_old, m_inf_old_old, h_inf_old_old, j_inf_old_old, oa_inf_old_old, oi_inf_old_old, ua_inf_old_old, ui_inf_old_old, xr_inf_old_old, xs_inf_old_old;
      double g_Kur_old_old, f_NaK_old_old, F_n_old_old, B1_old_old, B2_old_old;

      double tau_m_BN, tau_h_BN, m_inf_BN, h_inf_BN, I_BacNav;
      double tau_m_BN_old, tau_h_BN_old, m_inf_BN_old, h_inf_BN_old, I_BacNav_old;
      double tau_m_BN_old_old, tau_h_BN_old_old, m_inf_BN_old_old, h_inf_BN_old_old, I_BacNav_old_old;

      double phi [1000] = {};
      double x1 [1000] = {};
      double y1 [1000] = {};
      double threshSpiral = .25;//(xDim*10./ (double)nelx); //This is just to make sure that the spiral is thick enough
      
      if(integrator.compare(0,4,"SBDF") == 0){

        TransientLinearImplicitSystem &Bsystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
        Bsystem.update();

        const unsigned int vm_var = Bsystem.variable_number ("V");
        const DofMap & dof_map = Bsystem.get_dof_map();

        for(const auto & node : mesh.local_node_ptr_range()){

            dof_map.dof_indices (node, dof_indices);
            dof_map2.dof_indices (node, dof_indices2);

          if(dof_indices.size() > 0){

            const Real x = (*node)(0);
            const Real y = (*node)(1);
            const Real z = (*node)(2);

                if(datatime.timestep == 1){
                    v_old = 1.;
                    w_old = .999;
                    xs_old = 0.0187;
                    h_old = .965;
                    d_old = .000137;
                    xr_old = .0000329;
                    Nai_old = 11.2;
                    Ki_old = 139.;
                    Carel_old = 1.49;
                    oi_old = .999;
                    ui_old = .999;
                    m_old = .00291;
                    j_old = .978;
                    f_old = .999;
                    Cai_old = .000102;
                    Caup_old = 1.49;
                    oa_old = .0304;
                    ua_old = .00496;
                    fCa_old = .775;
                    u_old = 0.0;
                    vm_old = -81.2;//(*Bsystem.current_local_solution)  (dof_indices[vm_var]);

                    //BacNav
                    hBN_old = hBN_IC;
                    mBN_old = mBN_IC;

                    if(SpiralBool == 1){
                        if( x < (tissue_maxx + tissue_minx)/2. ){
                          //values at t = 285 after a stimulus, so in refractory period still.
                          //Basically, like if propagation wave is halfway through
                            v_old = 1.;
                            w_old = .9979;
                            xs_old = 0.1802;
                            h_old = .0895;
                            d_old = .001;
                            xr_old = .6142;
                            Nai_old = 11.1997;
                            Ki_old = 138.9995;
                            Carel_old = 1.4905;
                            oi_old = .3322;
                            ui_old = .9774;
                            m_old = .0357;
                            j_old = .0423;
                            f_old = .4576;
                            Cai_old = .00010355;
                            Caup_old = 1.4921;
                            oa_old = .0828;
                            ua_old = .0344;
                            fCa_old = .7717;
                            u_old = 0.0;
                            vm_old = -65.3708;

                            //BacNav
                            hBN_old = .0895;
                            mBN_old = .0357;

                            
                        }
                    }
                    else if(SpiralBool == 2){

                        double a = .5;
                        double k = .09;
                        int thickness = 15;
                        int finalDim = thickness*1000; //15,000
                        double bathRegion = xDim/2.0 - tissue_maxx;
                        double factorAxis = 3.0*(xDim/2.-bathRegion);
                        phi[0] = 0.;
                        x1[0] = -a*exp(k*phi[0])*cos(phi[0]);
                        y1[0] = -a*exp(k*phi[0])*sin(phi[0]);
                        for(int i = 1; i < 1000; i++){
                            phi[i] = phi[i-1] + ((factorAxis*acos(-1) - phi[0])/1000.);
                            x1[i] = -a*exp(k*phi[i])*cos(phi[i]);
                            y1[i] = -a*exp(k*phi[i])*sin(phi[i]);
                        }
                        for(int i = 0; i < 1000; i++){
                            if( std::abs(x - x1[i]) <= threshSpiral && std::abs(y - y1[i]) <= threshSpiral){
                                //conditions in notes seen from paraview
                                /*
                                 * initial conditions for multiple spirals to happen:
                                    u - spiral of 1. and rest of 0
                                    v - inverse of spiral of 1. and rest of 0
                                    w - spiral of .8 and rest of 1.
                                    s - same spiral of .1 and rest of 0
                                */
                                v_old = 1.;
                                w_old = .9979;
                                xs_old = 0.1802;
                                h_old = .0895;
                                d_old = .001;
                                xr_old = .6142;
                                Nai_old = 11.1997;
                                Ki_old = 138.9995;
                                Carel_old = 1.4905;
                                oi_old = .3322;
                                ui_old = .9774;
                                m_old = .0357;
                                j_old = .0423;
                                f_old = .4576;
                                Cai_old = .00010355;
                                Caup_old = 1.4921;
                                oa_old = .0828;
                                ua_old = .0344;
                                fCa_old = .7717;
                                u_old = 0.0;
                                vm_old = -65.3708;

                                //BacNav
                                hBN_old = .0895;
                                mBN_old = .0357;
                                //libMesh::out << "HERE --> x = " << x << ";     y = " << y << ";    with threshold of the spiral: " << threshSpiral << std::endl;
                            }
                        }
                    }
                }
                else{
                    v_old = (*Rsystem.current_local_solution) (dof_indices2[v_var]);
                    w_old = (*Rsystem.current_local_solution) (dof_indices2[w_var]);        
                    xs_old = (*Rsystem.current_local_solution) (dof_indices2[xs_var]);
                    h_old = (*Rsystem.current_local_solution) (dof_indices2[h_var]);
                    xr_old = (*Rsystem.current_local_solution) (dof_indices2[xr_var]);        
                    d_old = (*Rsystem.current_local_solution) (dof_indices2[d_var]);
                    Nai_old = (*Rsystem.current_local_solution) (dof_indices2[Nai_var]);
                    Ki_old = (*Rsystem.current_local_solution) (dof_indices2[Ki_var]);        
                    Carel_old = (*Rsystem.current_local_solution) (dof_indices2[Carel_var]);     
                    oi_old = (*Rsystem.current_local_solution) (dof_indices2[oi_var]);
                    ui_old = (*Rsystem.current_local_solution) (dof_indices2[ui_var]);        
                    m_old = (*Rsystem.current_local_solution) (dof_indices2[m_var]);    
                    j_old = (*Rsystem.current_local_solution) (dof_indices2[j_var]);
                    f_old = (*Rsystem.current_local_solution) (dof_indices2[f_var]);        
                    Cai_old = (*Rsystem.current_local_solution) (dof_indices2[Cai_var]);     
                    Caup_old = (*Rsystem.current_local_solution) (dof_indices2[Caup_var]);
                    oa_old = (*Rsystem.current_local_solution) (dof_indices2[oa_var]);        
                    ua_old = (*Rsystem.current_local_solution) (dof_indices2[ua_var]);     
                    fCa_old = (*Rsystem.current_local_solution) (dof_indices2[fCa_var]);
                    u_old = (*Rsystem.current_local_solution) (dof_indices2[u_var]);           
                    vm_old = (*Bsystem.current_local_solution)  (dof_indices[vm_var]);

                    hBN_old = (*Rsystem.current_local_solution) (dof_indices2[hBN_var]);
                    mBN_old = (*Rsystem.current_local_solution) (dof_indices2[mBN_var]);

                    /*
                    v_old_old = (Rsystem.get_vector("v_prev")) (dof_indices2[0]);
                    w_old_old = (Rsystem.get_vector("w_prev")) (dof_indices2[0]);        
                    xs_old_old = (Rsystem.get_vector("xs_prev")) (dof_indices2[0]);
                    h_old_old = (Rsystem.get_vector("h_prev")) (dof_indices2[0]);
                    xr_old_old = (Rsystem.get_vector("xr_prev")) (dof_indices2[0]);        
                    d_old_old = (Rsystem.get_vector("d_prev")) (dof_indices2[0]);
                    Nai_old_old = (Rsystem.get_vector("Nai_prev")) (dof_indices2[0]);
                    Ki_old_old = (Rsystem.get_vector("Ki_prev")) (dof_indices2[0]);       
                    Carel_old_old = (Rsystem.get_vector("Carel_prev")) (dof_indices2[0]);     
                    oi_old_old = (Rsystem.get_vector("oi_prev")) (dof_indices2[0]);
                    ui_old_old = (Rsystem.get_vector("ui_prev")) (dof_indices2[0]);       
                    m_old_old = (Rsystem.get_vector("m_prev")) (dof_indices2[0]);   
                    j_old_old = (Rsystem.get_vector("j_prev")) (dof_indices2[0]);
                    f_old_old = (Rsystem.get_vector("f_prev")) (dof_indices2[0]);        
                    Cai_old_old = (Rsystem.get_vector("Cai_prev")) (dof_indices2[0]);     
                    Caup_old_old = (Rsystem.get_vector("Caup_prev")) (dof_indices2[0]);
                    oa_old_old = (Rsystem.get_vector("oa_prev")) (dof_indices2[0]);     
                    ua_old_old = (Rsystem.get_vector("ua_prev")) (dof_indices2[0]);     
                    fCa_old_old = (Rsystem.get_vector("fCa_prev")) (dof_indices2[0]);
                    u_old_old = (Rsystem.get_vector("u_prev")) (dof_indices2[0]);          
                    vm_old_old = (Bsystem.get_vector("Vnm1"))  (dof_indices[0]);

                    v_old_old_old = (Rsystem.get_vector("v_prev_prev")) (dof_indices2[0]);
                    w_old_old_old = (Rsystem.get_vector("w_prev_prev")) (dof_indices2[0]);        
                    xs_old_old_old = (Rsystem.get_vector("xs_prev_prev")) (dof_indices2[0]);
                    h_old_old_old = (Rsystem.get_vector("h_prev_prev")) (dof_indices2[0]);
                    xr_old_old_old = (Rsystem.get_vector("xr_prev_prev")) (dof_indices2[0]);        
                    d_old_old_old = (Rsystem.get_vector("d_prev_prev")) (dof_indices2[0]);
                    Nai_old_old_old = (Rsystem.get_vector("Nai_prev_prev")) (dof_indices2[0]);
                    Ki_old_old_old = (Rsystem.get_vector("Ki_prev_prev")) (dof_indices2[0]);       
                    Carel_old_old_old = (Rsystem.get_vector("Carel_prev_prev")) (dof_indices2[0]);     
                    oi_old_old_old = (Rsystem.get_vector("oi_prev_prev")) (dof_indices2[0]);
                    ui_old_old_old = (Rsystem.get_vector("ui_prev_prev")) (dof_indices2[0]);       
                    m_old_old_old = (Rsystem.get_vector("m_prev_prev")) (dof_indices2[0]);   
                    j_old_old_old = (Rsystem.get_vector("j_prev_prev")) (dof_indices2[0]);
                    f_old_old_old = (Rsystem.get_vector("f_prev_prev")) (dof_indices2[0]);        
                    Cai_old_old_old = (Rsystem.get_vector("Cai_prev_prev")) (dof_indices2[0]);     
                    Caup_old_old_old = (Rsystem.get_vector("Caup_prev_prev")) (dof_indices2[0]);
                    oa_old_old_old = (Rsystem.get_vector("oa_prev_prev")) (dof_indices2[0]);     
                    ua_old_old_old = (Rsystem.get_vector("ua_prev_prev")) (dof_indices2[0]);     
                    fCa_old_old_old = (Rsystem.get_vector("fCa_prev_prev")) (dof_indices2[0]);
                    u_old_old_old = (Rsystem.get_vector("u_prev_prev")) (dof_indices2[0]);          
                    vm_old_old_old = (Bsystem.get_vector("Vnm2"))  (dof_indices[0]);
                    */
                }
                
                /*
                if(integrator.compare(0,5,"SBDF3") == 0){
                    if(datatime.timestep == 1){
                      //Nernst Potential Calculation
                      E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
                      E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
                      E_K = (R*T/(1.*F))*log(K_o/Ki_old);
                      //Fast Sodium Current
                      I_Na = g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
                      if(vm_old == -47.13){
                          alpha_m = 3.2;
                      }
                      else{
                          alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
                      }
                      beta_m = .08*exp(-vm_old/11.);
                      if(vm_old >= -40.){
                          alpha_h = 0.;
                      }
                      else{
                          alpha_h = .135*exp((vm_old + 80.)/(-6.8));
                      }
                      if(vm_old >= -40.){
                          beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
                      }
                      else{
                          beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
                      }
                      if(vm_old >= -40.){
                          alpha_j = 0.;
                      }
                      else{
                          alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
                      }
                      if(vm_old >= -40.){
                          beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
                      }
                      else{
                          beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
                      }
                      tau_m = 1./(alpha_m + beta_m);
                      tau_h = 1./(alpha_h + beta_h);
                      tau_j = 1./(alpha_j + beta_j);
                      m_inf = alpha_m/(alpha_m + beta_m);
                      h_inf = alpha_h/(alpha_h + beta_h);
                      j_inf = alpha_j/(alpha_j + beta_j);
                      //Time-independent K current
                      I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
                      //Transient outward K current
                      I_to = g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
                      oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
                      alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
                      beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
                      tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
                      oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
                      //Ultrarapid rapid delayed rectifier K current
                      g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
                      I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
                      alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
                      ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
                      alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
                      beta_ui = exp((vm_old - 158.)/16.);
                      tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
                      ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
                      //Rapid delayed outward rectifier K current
                      I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
                      alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
                      beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
                      tau_xr = 1./(alpha_xr + beta_xr);
                      xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
                      //Slow delayed outward rectifier K current
                      I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
                      alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
                      beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
                      tau_xs = 1./(2.*(alpha_xs + beta_xs));
                      xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
                      //L-type Ca current
                      I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
                      d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
                      tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
                      f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
                      fCa_inf = 1./(1. + Cai_old/.00035);
                      //NaK pump current
                      f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
                      I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
                      //NaCa exchanger current
                      I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
                      //Background currents
                      I_bCa = g_bCa*(vm_old - E_Ca);
                      I_bNa = g_bNa*(vm_old - E_Na);
                      //Ca pump current
                      I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
                      //Ca release current from JSR
                      I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
                      F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
                      u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
                      tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
                      w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
                      //Transfer current from NSR to JSR
                      I_tr = (Caup_old - Carel_old)/(tau_tr); 
                      //Ca uptake current by NSR
                      I_up = I_up_max/(1. + (K_up/Cai_old));
                      //Ca leak current by NSR
                      I_up_leak = I_up_max*Caup_old/(Ca_up_max);

                      //RHS of Concentrations
                      double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
                      double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
                      B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
                      B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
                      double Rhs_Cai = B1/B2;
                      double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
                      double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v = (v_inf - v_old)/tau_v;
                      double Rhs_w = (w_inf - w_old)/tau_w;
                      double Rhs_xs = (xs_inf - xs_old)/tau_xs;
                      double Rhs_h = (h_inf - h_old)/tau_h;
                      double Rhs_xr = (xr_inf - xr_old)/tau_xr;
                      double Rhs_d = (d_inf - d_old)/tau_d;
                      double Rhs_oi = (oi_inf - oi_old)/tau_oi;
                      double Rhs_ui = (ui_inf - ui_old)/tau_ui;
                      double Rhs_m = (m_inf - m_old)/tau_m;
                      double Rhs_j = (j_inf - j_old)/tau_j;
                      double Rhs_f = (f_inf - f_old)/tau_f;
                      double Rhs_u = (u_inf - u_old)/tau_u;
                      double Rhs_oa = (oa_inf - oa_old)/tau_oa;
                      double Rhs_ua = (ua_inf - ua_old)/tau_ua;
                      double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

                      //Integrating over time
                      Nai_new = Nai_old + dt*Rhs_Nai;
                      Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
                      Ki_new = Ki_old + dt*Rhs_Ki;
                      Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
                      Cai_new = Cai_old + dt*Rhs_Cai;
                      Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
                      Caup_new = Caup_old + dt*Rhs_Caup;
                      Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
                      Carel_new = Carel_old + dt*Rhs_Carel;
                      Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
                      v_new = v_old + dt*Rhs_v;
                      Rsystem.solution -> set(dof_indices2[v_var],v_new);
                      w_new = w_old + dt*Rhs_w;
                      Rsystem.solution -> set(dof_indices2[w_var],w_new);
                      xs_new = xs_old + dt*Rhs_xs;
                      Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
                      h_new = h_old + dt*Rhs_h;
                      Rsystem.solution -> set(dof_indices2[h_var],h_new);
                      xr_new = xr_old + dt*Rhs_xr;
                      Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
                      d_new = d_old + dt*Rhs_d;
                      Rsystem.solution -> set(dof_indices2[d_var],d_new);
                      oi_new = oi_old + dt*Rhs_oi;
                      Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
                      ui_new = ui_old + dt*Rhs_ui;
                      Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
                      m_new = m_old + dt*Rhs_m;
                      Rsystem.solution -> set(dof_indices2[m_var],m_new);
                      j_new = j_old + dt*Rhs_j;
                      Rsystem.solution -> set(dof_indices2[j_var],j_new);
                      f_new = f_old + dt*Rhs_f;
                      Rsystem.solution -> set(dof_indices2[f_var],f_new);
                      u_new = u_old + dt*Rhs_u;
                      Rsystem.solution -> set(dof_indices2[u_var],u_new);
                      oa_new = oa_old + dt*Rhs_oa;
                      Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
                      ua_new = ua_old + dt*Rhs_ua;
                      Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
                      fCa_new = fCa_old + dt*Rhs_fCa;
                      Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);
                    }
                    else if(datatime.timestep == 2){
                      //Nernst Potential Calculation
                      E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
                      E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
                      E_K = (R*T/(1.*F))*log(K_o/Ki_old);
                      //Fast Sodium Current
                      I_Na = g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
                      if(vm_old == -47.13){
                          alpha_m = 3.2;
                      }
                      else{
                          alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
                      }
                      beta_m = .08*exp(-vm_old/11.);
                      if(vm_old >= -40.){
                          alpha_h = 0.;
                      }
                      else{
                          alpha_h = .135*exp((vm_old + 80.)/(-6.8));
                      }
                      if(vm_old >= -40.){
                          beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
                      }
                      else{
                          beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
                      }
                      if(vm_old >= -40.){
                          alpha_j = 0.;
                      }
                      else{
                          alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
                      }
                      if(vm_old >= -40.){
                          beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
                      }
                      else{
                          beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
                      }
                      tau_m = 1./(alpha_m + beta_m);
                      tau_h = 1./(alpha_h + beta_h);
                      tau_j = 1./(alpha_j + beta_j);
                      m_inf = alpha_m/(alpha_m + beta_m);
                      h_inf = alpha_h/(alpha_h + beta_h);
                      j_inf = alpha_j/(alpha_j + beta_j);
                      //Time-independent K current
                      I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
                      //Transient outward K current
                      I_to = g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
                      oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
                      alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
                      beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
                      tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
                      oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
                      //Ultrarapid rapid delayed rectifier K current
                      g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
                      I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
                      alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
                      ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
                      alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
                      beta_ui = exp((vm_old - 158.)/16.);
                      tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
                      ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
                      //Rapid delayed outward rectifier K current
                      I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
                      alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
                      beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
                      tau_xr = 1./(alpha_xr + beta_xr);
                      xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
                      //Slow delayed outward rectifier K current
                      I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
                      alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
                      beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
                      tau_xs = 1./(2.*(alpha_xs + beta_xs));
                      xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
                      //L-type Ca current
                      I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
                      d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
                      tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
                      f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
                      fCa_inf = 1./(1. + Cai_old/.00035);
                      //NaK pump current
                      f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
                      I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
                      //NaCa exchanger current
                      I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
                      //Background currents
                      I_bCa = g_bCa*(vm_old - E_Ca);
                      I_bNa = g_bNa*(vm_old - E_Na);
                      //Ca pump current
                      I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
                      //Ca release current from JSR
                      I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
                      F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
                      u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
                      tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
                      w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
                      //Transfer current from NSR to JSR
                      I_tr = (Caup_old - Carel_old)/(tau_tr); 
                      //Ca uptake current by NSR
                      I_up = I_up_max/(1. + (K_up/Cai_old));
                      //Ca leak current by NSR
                      I_up_leak = I_up_max*Caup_old/(Ca_up_max);
                      //RHS of Concentrations
                      double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
                      double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
                      B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
                      B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
                      double Rhs_Cai = B1/B2;
                      double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
                      double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v = (v_inf - v_old)/tau_v;
                      double Rhs_w = (w_inf - w_old)/tau_w;
                      double Rhs_xs = (xs_inf - xs_old)/tau_xs;
                      double Rhs_h = (h_inf - h_old)/tau_h;
                      double Rhs_xr = (xr_inf - xr_old)/tau_xr;
                      double Rhs_d = (d_inf - d_old)/tau_d;
                      double Rhs_oi = (oi_inf - oi_old)/tau_oi;
                      double Rhs_ui = (ui_inf - ui_old)/tau_ui;
                      double Rhs_m = (m_inf - m_old)/tau_m;
                      double Rhs_j = (j_inf - j_old)/tau_j;
                      double Rhs_f = (f_inf - f_old)/tau_f;
                      double Rhs_u = (u_inf - u_old)/tau_u;
                      double Rhs_oa = (oa_inf - oa_old)/tau_oa;
                      double Rhs_ua = (ua_inf - ua_old)/tau_ua;
                      double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

                      //RHS of Concentrations
                      double Rhs_Nai_old = Nai_old_old;// (-3.*I_NaK_old - 3.*I_NaCa_old - I_bNa_old - I_Na_old)/(F*Vol_i);
                      double Rhs_Ki_old = Ki_old_old;// (2.*I_NaK_old - I_K1_old - I_to_old - I_Kur_old - I_Kr_old - I_Ks_old - I_bK)/(F*Vol_i);
                      //B1_old = ((2.*I_NaCa_old - I_pCa_old - I_CaL_old - I_bCa_old)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak_old - I_up_old) + I_rel_old*Vol_rel)/Vol_i);
                      //B2_old = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old_old + Km_Cmdn),2)));
                      double Rhs_Cai_old = Cai_old_old;// B1_old/B2_old;
                      double Rhs_Caup_old = Caup_old_old;// I_up_old - I_up_leak_old -I_tr_old*Vol_rel/Vol_up;
                      double Rhs_Carel_old = Carel_old_old;// (I_tr_old - I_rel_old)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v_old = v_old_old;// (v_inf_old - v_old_old)/tau_v_old;
                      double Rhs_w_old = w_old_old;// (w_inf_old - w_old_old)/tau_w_old;
                      double Rhs_xs_old = xs_old_old;// (xs_inf_old - xs_old_old)/tau_xs_old;
                      double Rhs_h_old = h_old_old;// (h_inf_old - h_old_old)/tau_h_old;
                      double Rhs_xr_old = xr_old_old;// (xr_inf_old - xr_old_old)/tau_xr_old;
                      double Rhs_d_old = d_old_old;// (d_inf_old - d_old_old)/tau_d_old;
                      double Rhs_oi_old = oi_old_old;// (oi_inf_old - oi_old_old)/tau_oi_old;
                      double Rhs_ui_old = ui_old_old;// (ui_inf_old - ui_old_old)/tau_ui_old;
                      double Rhs_m_old = m_old_old;// (m_inf_old - m_old_old)/tau_m_old;
                      double Rhs_j_old = j_old_old;// (j_inf_old - j_old_old)/tau_j_old;
                      double Rhs_f_old = f_old_old;// (f_inf_old - f_old_old)/tau_f_old;
                      double Rhs_u_old = u_old_old;// (u_inf_old - u_old_old)/tau_u;
                      double Rhs_oa_old = oa_old_old;// (oa_inf_old - oa_old_old)/tau_oa_old;
                      double Rhs_ua_old = ua_old_old;// (ua_inf_old - ua_old_old)/tau_ua_old;
                      double Rhs_fCa_old = fCa_old_old;// (fCa_inf_old - fCa_old_old)/tau_fCa;

                      //Integrating over time
                      v_new = v_old + (3./2.)*dt*Rhs_v - (1./2.)*dt*Rhs_v_old;
                      Rsystem.get_vector("v_prev").set(dof_indices2[0],Rhs_v);
                      Rsystem.solution -> set(dof_indices2[v_var],v_new);
                      w_new = w_old + (3./2.)*dt*Rhs_w - (1./2.)*dt*Rhs_w_old;
                      Rsystem.get_vector("w_prev").set(dof_indices2[0],Rhs_w);
                      Rsystem.solution -> set(dof_indices2[w_var],w_new);
                      xs_new = xs_old + (3./2.)*dt*Rhs_xs - (1./2.)*dt*Rhs_xs_old;
                      Rsystem.get_vector("xs_prev").set(dof_indices2[0],Rhs_xs);
                      Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
                      h_new = h_old + (3./2.)*dt*Rhs_h - (1./2.)*dt*Rhs_h_old;
                      Rsystem.get_vector("h_prev").set(dof_indices2[0],Rhs_h);
                      Rsystem.solution -> set(dof_indices2[h_var],h_new);
                      xr_new = xr_old + (3./2.)*dt*Rhs_xr - (1./2.)*dt*Rhs_xr_old;
                      Rsystem.get_vector("xr_prev").set(dof_indices2[0],Rhs_xr);
                      Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
                      d_new = d_old + (3./2.)*dt*Rhs_d - (1./2.)*dt*Rhs_d_old;
                      Rsystem.get_vector("d_prev").set(dof_indices2[0],Rhs_d);
                      Rsystem.solution -> set(dof_indices2[d_var],d_new);
                      oi_new = oi_old + (3./2.)*dt*Rhs_oi - (1./2.)*dt*Rhs_oi_old;
                      Rsystem.get_vector("oi_prev").set(dof_indices2[0],Rhs_oi);
                      Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
                      ui_new = ui_old + (3./2.)*dt*Rhs_ui - (1./2.)*dt*Rhs_ui_old;
                      Rsystem.get_vector("ui_prev").set(dof_indices2[0],Rhs_ui);
                      Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
                      m_new = m_old + (3./2.)*dt*Rhs_m - (1./2.)*dt*Rhs_m_old;
                      Rsystem.get_vector("m_prev").set(dof_indices2[0],Rhs_m);
                      Rsystem.solution -> set(dof_indices2[m_var],m_new);
                      j_new = j_old + (3./2.)*dt*Rhs_j - (1./2.)*dt*Rhs_j_old;
                      Rsystem.get_vector("j_prev").set(dof_indices2[0],Rhs_j);
                      Rsystem.solution -> set(dof_indices2[j_var],j_new);
                      f_new = f_old + (3./2.)*dt*Rhs_f - (1./2.)*dt*Rhs_f_old;
                      Rsystem.get_vector("f_prev").set(dof_indices2[0],Rhs_f);
                      Rsystem.solution -> set(dof_indices2[f_var],f_new);
                      u_new = u_old + (3./2.)*dt*Rhs_u - (1./2.)*dt*Rhs_u_old;
                      Rsystem.get_vector("u_prev").set(dof_indices2[0],Rhs_u);
                      Rsystem.solution -> set(dof_indices2[u_var],u_new);
                      oa_new = oa_old + (3./2.)*dt*Rhs_oa - (1./2.)*dt*Rhs_oa_old;
                      Rsystem.get_vector("oa_prev").set(dof_indices2[0],Rhs_oa);
                      Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
                      ua_new = ua_old + (3./2.)*dt*Rhs_ua - (1./2.)*dt*Rhs_ua_old;
                      Rsystem.get_vector("ua_prev").set(dof_indices2[0],Rhs_ua);
                      Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
                      fCa_new = fCa_old + (3./2.)*dt*Rhs_fCa - (1./2.)*dt*Rhs_fCa_old;
                      Rsystem.get_vector("fCa_prev").set(dof_indices2[0],Rhs_fCa);
                      Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);
                      Nai_new = Nai_old + (3./2.)*dt*Rhs_Nai - (1./2.)*dt*Rhs_Nai_old;
                      Rsystem.get_vector("Nai_prev").set(dof_indices2[0],Rhs_Nai);
                      Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
                      Ki_new = Ki_old + (3./2.)*dt*Rhs_Ki - (1./2.)*dt*Rhs_Ki_old;
                      Rsystem.get_vector("Ki_prev").set(dof_indices2[0],Rhs_Ki);
                      Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
                      Cai_new = Cai_old + (3./2.)*dt*Rhs_Cai - (1./2.)*dt*Rhs_Cai_old;
                      Rsystem.get_vector("Cai_prev").set(dof_indices2[0],Rhs_Cai);
                      Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
                      Caup_new = Caup_old + (3./2.)*dt*Rhs_Caup - (1./2.)*dt*Rhs_Caup_old;
                      Rsystem.get_vector("Caup_prev").set(dof_indices2[0],Rhs_Caup);
                      Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
                      Carel_new = Carel_old + (3./2.)*dt*Rhs_Carel - (1./2.)*dt*Rhs_Carel_old;
                      Rsystem.get_vector("Carel_prev").set(dof_indices2[0],Rhs_Carel);
                      Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
                    }
                    else{
                      //Nernst Potential Calculation
                      E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
                      E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
                      E_K = (R*T/(1.*F))*log(K_o/Ki_old);
                      //Fast Sodium Current
                      I_Na = g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
                      if(vm_old == -47.13){
                          alpha_m = 3.2;
                      }
                      else{
                          alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
                      }
                      beta_m = .08*exp(-vm_old/11.);
                      if(vm_old >= -40.){
                          alpha_h = 0.;
                      }
                      else{
                          alpha_h = .135*exp((vm_old + 80.)/(-6.8));
                      }
                      if(vm_old >= -40.){
                          beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
                      }
                      else{
                          beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
                      }
                      if(vm_old >= -40.){
                          alpha_j = 0.;
                      }
                      else{
                          alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
                      }
                      if(vm_old >= -40.){
                          beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
                      }
                      else{
                          beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
                      }
                      tau_m = 1./(alpha_m + beta_m);
                      tau_h = 1./(alpha_h + beta_h);
                      tau_j = 1./(alpha_j + beta_j);
                      m_inf = alpha_m/(alpha_m + beta_m);
                      h_inf = alpha_h/(alpha_h + beta_h);
                      j_inf = alpha_j/(alpha_j + beta_j);
                      //Time-independent K current
                      I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
                      //Transient outward K current
                      I_to = g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
                      oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
                      alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
                      beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
                      tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
                      oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
                      //Ultrarapid rapid delayed rectifier K current
                      g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
                      I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
                      alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
                      ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
                      alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
                      beta_ui = exp((vm_old - 158.)/16.);
                      tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
                      ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
                      //Rapid delayed outward rectifier K current
                      I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
                      alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
                      beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
                      tau_xr = 1./(alpha_xr + beta_xr);
                      xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
                      //Slow delayed outward rectifier K current
                      I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
                      alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
                      beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
                      tau_xs = 1./(2.*(alpha_xs + beta_xs));
                      xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
                      //L-type Ca current
                      I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
                      d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
                      tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
                      f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
                      fCa_inf = 1./(1. + Cai_old/.00035);
                      //NaK pump current
                      f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
                      I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
                      //NaCa exchanger current
                      I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
                      //Background currents
                      I_bCa = g_bCa*(vm_old - E_Ca);
                      I_bNa = g_bNa*(vm_old - E_Na);
                      //Ca pump current
                      I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
                      //Ca release current from JSR
                      I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
                      F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
                      u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
                      tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
                      w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
                      //Transfer current from NSR to JSR
                      I_tr = (Caup_old - Carel_old)/(tau_tr); 
                      //Ca uptake current by NSR
                      I_up = I_up_max/(1. + (K_up/Cai_old));
                      //Ca leak current by NSR
                      I_up_leak = I_up_max*Caup_old/(Ca_up_max);
                      //RHS of Concentrations
                      double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
                      double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
                      B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
                      B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
                      double Rhs_Cai = B1/B2;
                      double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
                      double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v = (v_inf - v_old)/tau_v;
                      double Rhs_w = (w_inf - w_old)/tau_w;
                      double Rhs_xs = (xs_inf - xs_old)/tau_xs;
                      double Rhs_h = (h_inf - h_old)/tau_h;
                      double Rhs_xr = (xr_inf - xr_old)/tau_xr;
                      double Rhs_d = (d_inf - d_old)/tau_d;
                      double Rhs_oi = (oi_inf - oi_old)/tau_oi;
                      double Rhs_ui = (ui_inf - ui_old)/tau_ui;
                      double Rhs_m = (m_inf - m_old)/tau_m;
                      double Rhs_j = (j_inf - j_old)/tau_j;
                      double Rhs_f = (f_inf - f_old)/tau_f;
                      double Rhs_u = (u_inf - u_old)/tau_u;
                      double Rhs_oa = (oa_inf - oa_old)/tau_oa;
                      double Rhs_ua = (ua_inf - ua_old)/tau_ua;
                      double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

                      //RHS of Concentrations
                      double Rhs_Nai_old = Nai_old_old;// (-3.*I_NaK_old - 3.*I_NaCa_old - I_bNa_old - I_Na_old)/(F*Vol_i);
                      double Rhs_Ki_old = Ki_old_old;// (2.*I_NaK_old - I_K1_old - I_to_old - I_Kur_old - I_Kr_old - I_Ks_old - I_bK)/(F*Vol_i);
                      //B1_old = ((2.*I_NaCa_old - I_pCa_old - I_CaL_old - I_bCa_old)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak_old - I_up_old) + I_rel_old*Vol_rel)/Vol_i);
                      //B2_old = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old_old + Km_Cmdn),2)));
                      double Rhs_Cai_old = Cai_old_old;// B1_old/B2_old;
                      double Rhs_Caup_old = Caup_old_old;// I_up_old - I_up_leak_old -I_tr_old*Vol_rel/Vol_up;
                      double Rhs_Carel_old = Carel_old_old;// (I_tr_old - I_rel_old)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v_old = v_old_old;// (v_inf_old - v_old_old)/tau_v_old;
                      double Rhs_w_old = w_old_old;// (w_inf_old - w_old_old)/tau_w_old;
                      double Rhs_xs_old = xs_old_old;// (xs_inf_old - xs_old_old)/tau_xs_old;
                      double Rhs_h_old = h_old_old;// (h_inf_old - h_old_old)/tau_h_old;
                      double Rhs_xr_old = xr_old_old;// (xr_inf_old - xr_old_old)/tau_xr_old;
                      double Rhs_d_old = d_old_old;// (d_inf_old - d_old_old)/tau_d_old;
                      double Rhs_oi_old = oi_old_old;// (oi_inf_old - oi_old_old)/tau_oi_old;
                      double Rhs_ui_old = ui_old_old;// (ui_inf_old - ui_old_old)/tau_ui_old;
                      double Rhs_m_old = m_old_old;// (m_inf_old - m_old_old)/tau_m_old;
                      double Rhs_j_old = j_old_old;// (j_inf_old - j_old_old)/tau_j_old;
                      double Rhs_f_old = f_old_old;// (f_inf_old - f_old_old)/tau_f_old;
                      double Rhs_u_old = u_old_old;// (u_inf_old - u_old_old)/tau_u;
                      double Rhs_oa_old = oa_old_old;// (oa_inf_old - oa_old_old)/tau_oa_old;
                      double Rhs_ua_old = ua_old_old;// (ua_inf_old - ua_old_old)/tau_ua_old;
                      double Rhs_fCa_old = fCa_old_old;// (fCa_inf_old - fCa_old_old)/tau_fCa;

                      //RHS of Concentrations
                      double Rhs_Nai_old_old = Nai_old_old_old;// (-3.*I_NaK_old - 3.*I_NaCa_old - I_bNa_old - I_Na_old)/(F*Vol_i);
                      double Rhs_Ki_old_old = Ki_old_old_old;// (2.*I_NaK_old - I_K1_old - I_to_old - I_Kur_old - I_Kr_old - I_Ks_old - I_bK)/(F*Vol_i);
                      //B1_old = ((2.*I_NaCa_old - I_pCa_old - I_CaL_old - I_bCa_old)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak_old - I_up_old) + I_rel_old*Vol_rel)/Vol_i);
                      //B2_old = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old_old + Km_Cmdn),2)));
                      double Rhs_Cai_old_old = Cai_old_old_old;// B1_old/B2_old;
                      double Rhs_Caup_old_old = Caup_old_old_old;// I_up_old - I_up_leak_old -I_tr_old*Vol_rel/Vol_up;
                      double Rhs_Carel_old_old = Carel_old_old_old;// (I_tr_old - I_rel_old)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v_old_old = v_old_old_old;// (v_inf_old - v_old_old)/tau_v_old;
                      double Rhs_w_old_old = w_old_old_old;// (w_inf_old - w_old_old)/tau_w_old;
                      double Rhs_xs_old_old = xs_old_old_old;// (xs_inf_old - xs_old_old)/tau_xs_old;
                      double Rhs_h_old_old = h_old_old_old;// (h_inf_old - h_old_old)/tau_h_old;
                      double Rhs_xr_old_old = xr_old_old_old;// (xr_inf_old - xr_old_old)/tau_xr_old;
                      double Rhs_d_old_old = d_old_old_old;// (d_inf_old - d_old_old)/tau_d_old;
                      double Rhs_oi_old_old = oi_old_old_old;// (oi_inf_old - oi_old_old)/tau_oi_old;
                      double Rhs_ui_old_old = ui_old_old_old;// (ui_inf_old - ui_old_old)/tau_ui_old;
                      double Rhs_m_old_old = m_old_old_old;// (m_inf_old - m_old_old)/tau_m_old;
                      double Rhs_j_old_old = j_old_old_old;// (j_inf_old - j_old_old)/tau_j_old;
                      double Rhs_f_old_old = f_old_old_old;// (f_inf_old - f_old_old)/tau_f_old;
                      double Rhs_u_old_old = u_old_old_old;// (u_inf_old - u_old_old)/tau_u;
                      double Rhs_oa_old_old = oa_old_old_old;// (oa_inf_old - oa_old_old)/tau_oa_old;
                      double Rhs_ua_old_old = ua_old_old_old;// (ua_inf_old - ua_old_old)/tau_ua_old;
                      double Rhs_fCa_old_old = fCa_old_old_old;// (fCa_inf_old - fCa_old_old)/tau_fCa;
                      
                      v_new = v_old + (23./12.)*dt*Rhs_v - (16./12.)*dt*Rhs_v_old + (5./12.)*dt*Rhs_v_old_old;
                      Rsystem.get_vector("v_prev").set(dof_indices2[0],Rhs_v);
                      Rsystem.get_vector("v_prev_prev").set(dof_indices2[0],Rhs_v_old);
                      Rsystem.solution -> set(dof_indices2[v_var],v_new);
                      w_new = w_old + (23./12.)*dt*Rhs_w - (16./12.)*dt*Rhs_w_old + (5./12.)*dt*Rhs_w_old_old;
                      Rsystem.get_vector("w_prev").set(dof_indices2[0],Rhs_w);
                      Rsystem.get_vector("w_prev_prev").set(dof_indices2[0],Rhs_w_old);
                      Rsystem.solution -> set(dof_indices2[w_var],w_new);
                      xs_new = xs_old + (23./12.)*dt*Rhs_xs - (16./12.)*dt*Rhs_xs_old + (5./12.)*dt*Rhs_xs_old_old;
                      Rsystem.get_vector("xs_prev").set(dof_indices2[0],Rhs_xs);
                      Rsystem.get_vector("xs_prev_prev").set(dof_indices2[0],Rhs_xs_old);
                      Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
                      h_new = h_old + (23./12.)*dt*Rhs_h - (16./12.)*dt*Rhs_h_old + (5./12.)*dt*Rhs_h_old_old;
                      Rsystem.get_vector("h_prev").set(dof_indices2[0],Rhs_h);
                      Rsystem.get_vector("h_prev_prev").set(dof_indices2[0],Rhs_h_old);
                      Rsystem.solution -> set(dof_indices2[h_var],h_new);
                      d_new = d_old + (23./12.)*dt*Rhs_d - (16./12.)*dt*Rhs_d_old + (5./12.)*dt*Rhs_d_old_old;
                      Rsystem.get_vector("d_prev").set(dof_indices2[0],Rhs_d);
                      Rsystem.get_vector("d_prev_prev").set(dof_indices2[0],Rhs_d_old);
                      Rsystem.solution -> set(dof_indices2[d_var],d_new);
                      xr_new = xr_old + (23./12.)*dt*Rhs_xr - (16./12.)*dt*Rhs_xr_old + (5./12.)*dt*Rhs_xr_old_old;
                      Rsystem.get_vector("xr_prev").set(dof_indices2[0],Rhs_xr);
                      Rsystem.get_vector("xr_prev_prev").set(dof_indices2[0],Rhs_xr_old);
                      Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
                      oi_new = oi_old + (23./12.)*dt*Rhs_oi - (16./12.)*dt*Rhs_oi_old + (5./12.)*dt*Rhs_oi_old_old;
                      Rsystem.get_vector("oi_prev").set(dof_indices2[0],Rhs_oi);
                      Rsystem.get_vector("oi_prev_prev").set(dof_indices2[0],Rhs_oi_old);
                      Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
                      ui_new = ui_old + (23./12.)*dt*Rhs_ui - (16./12.)*dt*Rhs_ui_old + (5./12.)*dt*Rhs_ui_old_old;
                      Rsystem.get_vector("ui_prev").set(dof_indices2[0],Rhs_ui);
                      Rsystem.get_vector("ui_prev_prev").set(dof_indices2[0],Rhs_ui_old);
                      Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
                      m_new = m_old + (23./12.)*dt*Rhs_m - (16./12.)*dt*Rhs_m_old + (5./12.)*dt*Rhs_m_old_old;
                      Rsystem.get_vector("m_prev").set(dof_indices2[0],Rhs_m);
                      Rsystem.get_vector("m_prev_prev").set(dof_indices2[0],Rhs_m_old);
                      Rsystem.solution -> set(dof_indices2[m_var],m_new);
                      j_new = j_old + (23./12.)*dt*Rhs_j - (16./12.)*dt*Rhs_j_old + (5./12.)*dt*Rhs_j_old_old;
                      Rsystem.get_vector("j_prev").set(dof_indices2[0],Rhs_j);
                      Rsystem.get_vector("j_prev_prev").set(dof_indices2[0],Rhs_j_old);
                      Rsystem.solution -> set(dof_indices2[j_var],j_new);
                      f_new = f_old + (23./12.)*dt*Rhs_f - (16./12.)*dt*Rhs_f_old + (5./12.)*dt*Rhs_f_old_old;
                      Rsystem.get_vector("f_prev").set(dof_indices2[0],Rhs_f);
                      Rsystem.get_vector("f_prev_prev").set(dof_indices2[0],Rhs_f_old);
                      Rsystem.solution -> set(dof_indices2[f_var],f_new);
                      u_new = u_old + (23./12.)*dt*Rhs_u - (16./12.)*dt*Rhs_u_old + (5./12.)*dt*Rhs_u_old_old;
                      Rsystem.get_vector("u_prev").set(dof_indices2[0],Rhs_u);
                      Rsystem.get_vector("u_prev_prev").set(dof_indices2[0],Rhs_u_old);
                      Rsystem.solution -> set(dof_indices2[u_var],u_new);
                      oa_new = oa_old + (23./12.)*dt*Rhs_oa - (16./12.)*dt*Rhs_oa_old + (5./12.)*dt*Rhs_oa_old_old;
                      Rsystem.get_vector("oa_prev").set(dof_indices2[0],Rhs_oa);
                      Rsystem.get_vector("oa_prev_prev").set(dof_indices2[0],Rhs_oa_old);
                      Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
                      ua_new = ua_old + (23./12.)*dt*Rhs_ua - (16./12.)*dt*Rhs_ua_old + (5./12.)*dt*Rhs_ua_old_old;
                      Rsystem.get_vector("ua_prev").set(dof_indices2[0],Rhs_ua);
                      Rsystem.get_vector("ua_prev_prev").set(dof_indices2[0],Rhs_ua_old);
                      Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
                      fCa_new = fCa_old + (23./12.)*dt*Rhs_fCa - (16./12.)*dt*Rhs_fCa_old + (5./12.)*dt*Rhs_fCa_old_old;
                      Rsystem.get_vector("fCa_prev").set(dof_indices2[0],Rhs_fCa);
                      Rsystem.get_vector("fCa_prev_prev").set(dof_indices2[0],Rhs_fCa_old);
                      Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);
                      Cai_new = Cai_old + (23./12.)*dt*Rhs_Cai - (16./12.)*dt*Rhs_Cai_old + (5./12.)*dt*Rhs_Cai_old_old;
                      Rsystem.get_vector("Cai_prev").set(dof_indices2[0],Rhs_Cai);
                      Rsystem.get_vector("Cai_prev_prev").set(dof_indices2[0],Rhs_Cai_old);
                      Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
                      Caup_new = Caup_old + (23./12.)*dt*Rhs_Caup - (16./12.)*dt*Rhs_Caup_old + (5./12.)*dt*Rhs_Caup_old_old;
                      Rsystem.get_vector("Caup_prev").set(dof_indices2[0],Rhs_Caup);
                      Rsystem.get_vector("Caup_prev_prev").set(dof_indices2[0],Rhs_Caup_old);
                      Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
                      Carel_new = Carel_old + (23./12.)*dt*Rhs_Carel - (16./12.)*dt*Rhs_Carel_old + (5./12.)*dt*Rhs_Carel_old_old;
                      Rsystem.get_vector("Carel_prev").set(dof_indices2[0],Rhs_Carel);
                      Rsystem.get_vector("Carel_prev_prev").set(dof_indices2[0],Rhs_Carel_old);
                      Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
                      Nai_new = Nai_old + (23./12.)*dt*Rhs_Nai - (16./12.)*dt*Rhs_Nai_old + (5./12.)*dt*Rhs_Nai_old_old;
                      Rsystem.get_vector("Nai_prev").set(dof_indices2[0],Rhs_Nai);
                      Rsystem.get_vector("Nai_prev_prev").set(dof_indices2[0],Rhs_Nai_old);
                      Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
                      Ki_new = Ki_old + (23./12.)*dt*Rhs_Ki - (16./12.)*dt*Rhs_Ki_old + (5./12.)*dt*Rhs_Ki_old_old;
                      Rsystem.get_vector("Ki_prev").set(dof_indices2[0],Rhs_Ki);
                      Rsystem.get_vector("Ki_prev_prev").set(dof_indices2[0],Rhs_Ki_old);
                      Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
                      
                    }

                }
                else if(integrator.compare(0,5,"SBDF2") == 0){

                    if(datatime.timestep == 1){
                      //Nernst Potential Calculation
                      E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
                      E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
                      E_K = (R*T/(1.*F))*log(K_o/Ki_old);
                      //Fast Sodium Current
                      I_Na = g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
                      if(vm_old == -47.13){
                          alpha_m = 3.2;
                      }
                      else{
                          alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
                      }
                      beta_m = .08*exp(-vm_old/11.);
                      if(vm_old >= -40.){
                          alpha_h = 0.;
                      }
                      else{
                          alpha_h = .135*exp((vm_old + 80.)/(-6.8));
                      }
                      if(vm_old >= -40.){
                          beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
                      }
                      else{
                          beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
                      }
                      if(vm_old >= -40.){
                          alpha_j = 0.;
                      }
                      else{
                          alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
                      }
                      if(vm_old >= -40.){
                          beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
                      }
                      else{
                          beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
                      }
                      tau_m = 1./(alpha_m + beta_m);
                      tau_h = 1./(alpha_h + beta_h);
                      tau_j = 1./(alpha_j + beta_j);
                      m_inf = alpha_m/(alpha_m + beta_m);
                      h_inf = alpha_h/(alpha_h + beta_h);
                      j_inf = alpha_j/(alpha_j + beta_j);
                      //Time-independent K current
                      I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
                      //Transient outward K current
                      I_to = g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
                      oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
                      alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
                      beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
                      tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
                      oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
                      //Ultrarapid rapid delayed rectifier K current
                      g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
                      I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
                      alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
                      ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
                      alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
                      beta_ui = exp((vm_old - 158.)/16.);
                      tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
                      ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
                      //Rapid delayed outward rectifier K current
                      I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
                      alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
                      beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
                      tau_xr = 1./(alpha_xr + beta_xr);
                      xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
                      //Slow delayed outward rectifier K current
                      I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
                      alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
                      beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
                      tau_xs = 1./(2.*(alpha_xs + beta_xs));
                      xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
                      //L-type Ca current
                      I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
                      d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
                      tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
                      f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
                      fCa_inf = 1./(1. + Cai_old/.00035);
                      //NaK pump current
                      f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
                      I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
                      //NaCa exchanger current
                      I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
                      //Background currents
                      I_bCa = g_bCa*(vm_old - E_Ca);
                      I_bNa = g_bNa*(vm_old - E_Na);
                      //Ca pump current
                      I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
                      //Ca release current from JSR
                      I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
                      F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
                      u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
                      tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
                      w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
                      //Transfer current from NSR to JSR
                      I_tr = (Caup_old - Carel_old)/(tau_tr); 
                      //Ca uptake current by NSR
                      I_up = I_up_max/(1. + (K_up/Cai_old));
                      //Ca leak current by NSR
                      I_up_leak = I_up_max*Caup_old/(Ca_up_max);

                      //RHS of Concentrations
                      double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
                      double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
                      B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
                      B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
                      double Rhs_Cai = B1/B2;
                      double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
                      double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v = (v_inf - v_old)/tau_v;
                      double Rhs_w = (w_inf - w_old)/tau_w;
                      double Rhs_xs = (xs_inf - xs_old)/tau_xs;
                      double Rhs_h = (h_inf - h_old)/tau_h;
                      double Rhs_xr = (xr_inf - xr_old)/tau_xr;
                      double Rhs_d = (d_inf - d_old)/tau_d;
                      double Rhs_oi = (oi_inf - oi_old)/tau_oi;
                      double Rhs_ui = (ui_inf - ui_old)/tau_ui;
                      double Rhs_m = (m_inf - m_old)/tau_m;
                      double Rhs_j = (j_inf - j_old)/tau_j;
                      double Rhs_f = (f_inf - f_old)/tau_f;
                      double Rhs_u = (u_inf - u_old)/tau_u;
                      double Rhs_oa = (oa_inf - oa_old)/tau_oa;
                      double Rhs_ua = (ua_inf - ua_old)/tau_ua;
                      double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

                      //Integrating over time
                      Nai_new = Nai_old + dt*Rhs_Nai;
                      Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
                      Ki_new = Ki_old + dt*Rhs_Ki;
                      Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
                      Cai_new = Cai_old + dt*Rhs_Cai;
                      Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
                      Caup_new = Caup_old + dt*Rhs_Caup;
                      Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
                      Carel_new = Carel_old + dt*Rhs_Carel;
                      Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
                      v_new = v_old + dt*Rhs_v;
                      Rsystem.solution -> set(dof_indices2[v_var],v_new);
                      w_new = w_old + dt*Rhs_w;
                      Rsystem.solution -> set(dof_indices2[w_var],w_new);
                      xs_new = xs_old + dt*Rhs_xs;
                      Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
                      h_new = h_old + dt*Rhs_h;
                      Rsystem.solution -> set(dof_indices2[h_var],h_new);
                      xr_new = xr_old + dt*Rhs_xr;
                      Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
                      d_new = d_old + dt*Rhs_d;
                      Rsystem.solution -> set(dof_indices2[d_var],d_new);
                      oi_new = oi_old + dt*Rhs_oi;
                      Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
                      ui_new = ui_old + dt*Rhs_ui;
                      Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
                      m_new = m_old + dt*Rhs_m;
                      Rsystem.solution -> set(dof_indices2[m_var],m_new);
                      j_new = j_old + dt*Rhs_j;
                      Rsystem.solution -> set(dof_indices2[j_var],j_new);
                      f_new = f_old + dt*Rhs_f;
                      Rsystem.solution -> set(dof_indices2[f_var],f_new);
                      u_new = u_old + dt*Rhs_u;
                      Rsystem.solution -> set(dof_indices2[u_var],u_new);
                      oa_new = oa_old + dt*Rhs_oa;
                      Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
                      ua_new = ua_old + dt*Rhs_ua;
                      Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
                      fCa_new = fCa_old + dt*Rhs_fCa;
                      Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);
                    }
                    else{
                      //Nernst Potential Calculation
                      E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
                      E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
                      E_K = (R*T/(1.*F))*log(K_o/Ki_old);
                      //Fast Sodium Current
                      I_Na = g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
                      if(vm_old == -47.13){
                          alpha_m = 3.2;
                      }
                      else{
                          alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
                      }
                      beta_m = .08*exp(-vm_old/11.);
                      if(vm_old >= -40.){
                          alpha_h = 0.;
                      }
                      else{
                          alpha_h = .135*exp((vm_old + 80.)/(-6.8));
                      }
                      if(vm_old >= -40.){
                          beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
                      }
                      else{
                          beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
                      }
                      if(vm_old >= -40.){
                          alpha_j = 0.;
                      }
                      else{
                          alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
                      }
                      if(vm_old >= -40.){
                          beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
                      }
                      else{
                          beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
                      }
                      tau_m = 1./(alpha_m + beta_m);
                      tau_h = 1./(alpha_h + beta_h);
                      tau_j = 1./(alpha_j + beta_j);
                      m_inf = alpha_m/(alpha_m + beta_m);
                      h_inf = alpha_h/(alpha_h + beta_h);
                      j_inf = alpha_j/(alpha_j + beta_j);
                      //Time-independent K current
                      I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
                      //Transient outward K current
                      I_to = g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
                      oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
                      alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
                      beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
                      tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
                      oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
                      //Ultrarapid rapid delayed rectifier K current
                      g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
                      I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
                      alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
                      ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
                      alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
                      beta_ui = exp((vm_old - 158.)/16.);
                      tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
                      ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
                      //Rapid delayed outward rectifier K current
                      I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
                      alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
                      beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
                      tau_xr = 1./(alpha_xr + beta_xr);
                      xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
                      //Slow delayed outward rectifier K current
                      I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
                      alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
                      beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
                      tau_xs = 1./(2.*(alpha_xs + beta_xs));
                      xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
                      //L-type Ca current
                      I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
                      d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
                      tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
                      f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
                      fCa_inf = 1./(1. + Cai_old/.00035);
                      //NaK pump current
                      f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
                      I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
                      //NaCa exchanger current
                      I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
                      //Background currents
                      I_bCa = g_bCa*(vm_old - E_Ca);
                      I_bNa = g_bNa*(vm_old - E_Na);
                      //Ca pump current
                      I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
                      //Ca release current from JSR
                      I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
                      F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
                      u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
                      tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
                      w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
                      //Transfer current from NSR to JSR
                      I_tr = (Caup_old - Carel_old)/(tau_tr); 
                      //Ca uptake current by NSR
                      I_up = I_up_max/(1. + (K_up/Cai_old));
                      //Ca leak current by NSR
                      I_up_leak = I_up_max*Caup_old/(Ca_up_max);
                      //RHS of Concentrations
                      double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
                      double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
                      B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
                      B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
                      double Rhs_Cai = B1/B2;
                      double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
                      double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v = (v_inf - v_old)/tau_v;
                      double Rhs_w = (w_inf - w_old)/tau_w;
                      double Rhs_xs = (xs_inf - xs_old)/tau_xs;
                      double Rhs_h = (h_inf - h_old)/tau_h;
                      double Rhs_xr = (xr_inf - xr_old)/tau_xr;
                      double Rhs_d = (d_inf - d_old)/tau_d;
                      double Rhs_oi = (oi_inf - oi_old)/tau_oi;
                      double Rhs_ui = (ui_inf - ui_old)/tau_ui;
                      double Rhs_m = (m_inf - m_old)/tau_m;
                      double Rhs_j = (j_inf - j_old)/tau_j;
                      double Rhs_f = (f_inf - f_old)/tau_f;
                      double Rhs_u = (u_inf - u_old)/tau_u;
                      double Rhs_oa = (oa_inf - oa_old)/tau_oa;
                      double Rhs_ua = (ua_inf - ua_old)/tau_ua;
                      double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

                      //RHS of Concentrations
                      double Rhs_Nai_old = Nai_old_old;// (-3.*I_NaK_old - 3.*I_NaCa_old - I_bNa_old - I_Na_old)/(F*Vol_i);
                      double Rhs_Ki_old = Ki_old_old;// (2.*I_NaK_old - I_K1_old - I_to_old - I_Kur_old - I_Kr_old - I_Ks_old - I_bK)/(F*Vol_i);
                      //B1_old = ((2.*I_NaCa_old - I_pCa_old - I_CaL_old - I_bCa_old)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak_old - I_up_old) + I_rel_old*Vol_rel)/Vol_i);
                      //B2_old = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old_old + Km_Cmdn),2)));
                      double Rhs_Cai_old = Cai_old_old;// B1_old/B2_old;
                      double Rhs_Caup_old = Caup_old_old;// I_up_old - I_up_leak_old -I_tr_old*Vol_rel/Vol_up;
                      double Rhs_Carel_old = Carel_old_old;// (I_tr_old - I_rel_old)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v_old = v_old_old;// (v_inf_old - v_old_old)/tau_v_old;
                      double Rhs_w_old = w_old_old;// (w_inf_old - w_old_old)/tau_w_old;
                      double Rhs_xs_old = xs_old_old;// (xs_inf_old - xs_old_old)/tau_xs_old;
                      double Rhs_h_old = h_old_old;// (h_inf_old - h_old_old)/tau_h_old;
                      double Rhs_xr_old = xr_old_old;// (xr_inf_old - xr_old_old)/tau_xr_old;
                      double Rhs_d_old = d_old_old;// (d_inf_old - d_old_old)/tau_d_old;
                      double Rhs_oi_old = oi_old_old;// (oi_inf_old - oi_old_old)/tau_oi_old;
                      double Rhs_ui_old = ui_old_old;// (ui_inf_old - ui_old_old)/tau_ui_old;
                      double Rhs_m_old = m_old_old;// (m_inf_old - m_old_old)/tau_m_old;
                      double Rhs_j_old = j_old_old;// (j_inf_old - j_old_old)/tau_j_old;
                      double Rhs_f_old = f_old_old;// (f_inf_old - f_old_old)/tau_f_old;
                      double Rhs_u_old = u_old_old;// (u_inf_old - u_old_old)/tau_u;
                      double Rhs_oa_old = oa_old_old;// (oa_inf_old - oa_old_old)/tau_oa_old;
                      double Rhs_ua_old = ua_old_old;// (ua_inf_old - ua_old_old)/tau_ua_old;
                      double Rhs_fCa_old = fCa_old_old;// (fCa_inf_old - fCa_old_old)/tau_fCa;

                      //Integrating over time
                      v_new = v_old + (3./2.)*dt*Rhs_v - (1./2.)*dt*Rhs_v_old;
                      Rsystem.get_vector("v_prev").set(dof_indices2[0],Rhs_v);
                      Rsystem.solution -> set(dof_indices2[v_var],v_new);
                      w_new = w_old + (3./2.)*dt*Rhs_w - (1./2.)*dt*Rhs_w_old;
                      Rsystem.get_vector("w_prev").set(dof_indices2[0],Rhs_w);
                      Rsystem.solution -> set(dof_indices2[w_var],w_new);
                      xs_new = xs_old + (3./2.)*dt*Rhs_xs - (1./2.)*dt*Rhs_xs_old;
                      Rsystem.get_vector("xs_prev").set(dof_indices2[0],Rhs_xs);
                      Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
                      h_new = h_old + (3./2.)*dt*Rhs_h - (1./2.)*dt*Rhs_h_old;
                      Rsystem.get_vector("h_prev").set(dof_indices2[0],Rhs_h);
                      Rsystem.solution -> set(dof_indices2[h_var],h_new);
                      xr_new = xr_old + (3./2.)*dt*Rhs_xr - (1./2.)*dt*Rhs_xr_old;
                      Rsystem.get_vector("xr_prev").set(dof_indices2[0],Rhs_xr);
                      Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
                      d_new = d_old + (3./2.)*dt*Rhs_d - (1./2.)*dt*Rhs_d_old;
                      Rsystem.get_vector("d_prev").set(dof_indices2[0],Rhs_d);
                      Rsystem.solution -> set(dof_indices2[d_var],d_new);
                      oi_new = oi_old + (3./2.)*dt*Rhs_oi - (1./2.)*dt*Rhs_oi_old;
                      Rsystem.get_vector("oi_prev").set(dof_indices2[0],Rhs_oi);
                      Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
                      ui_new = ui_old + (3./2.)*dt*Rhs_ui - (1./2.)*dt*Rhs_ui_old;
                      Rsystem.get_vector("ui_prev").set(dof_indices2[0],Rhs_ui);
                      Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
                      m_new = m_old + (3./2.)*dt*Rhs_m - (1./2.)*dt*Rhs_m_old;
                      Rsystem.get_vector("m_prev").set(dof_indices2[0],Rhs_m);
                      Rsystem.solution -> set(dof_indices2[m_var],m_new);
                      j_new = j_old + (3./2.)*dt*Rhs_j - (1./2.)*dt*Rhs_j_old;
                      Rsystem.get_vector("j_prev").set(dof_indices2[0],Rhs_j);
                      Rsystem.solution -> set(dof_indices2[j_var],j_new);
                      f_new = f_old + (3./2.)*dt*Rhs_f - (1./2.)*dt*Rhs_f_old;
                      Rsystem.get_vector("f_prev").set(dof_indices2[0],Rhs_f);
                      Rsystem.solution -> set(dof_indices2[f_var],f_new);
                      u_new = u_old + (3./2.)*dt*Rhs_u - (1./2.)*dt*Rhs_u_old;
                      Rsystem.get_vector("u_prev").set(dof_indices2[0],Rhs_u);
                      Rsystem.solution -> set(dof_indices2[u_var],u_new);
                      oa_new = oa_old + (3./2.)*dt*Rhs_oa - (1./2.)*dt*Rhs_oa_old;
                      Rsystem.get_vector("oa_prev").set(dof_indices2[0],Rhs_oa);
                      Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
                      ua_new = ua_old + (3./2.)*dt*Rhs_ua - (1./2.)*dt*Rhs_ua_old;
                      Rsystem.get_vector("ua_prev").set(dof_indices2[0],Rhs_ua);
                      Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
                      fCa_new = fCa_old + (3./2.)*dt*Rhs_fCa - (1./2.)*dt*Rhs_fCa_old;
                      Rsystem.get_vector("fCa_prev").set(dof_indices2[0],Rhs_fCa);
                      Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);
                      Nai_new = Nai_old + (3./2.)*dt*Rhs_Nai - (1./2.)*dt*Rhs_Nai_old;
                      Rsystem.get_vector("Nai_prev").set(dof_indices2[0],Rhs_Nai);
                      Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
                      Ki_new = Ki_old + (3./2.)*dt*Rhs_Ki - (1./2.)*dt*Rhs_Ki_old;
                      Rsystem.get_vector("Ki_prev").set(dof_indices2[0],Rhs_Ki);
                      Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
                      Cai_new = Cai_old + (3./2.)*dt*Rhs_Cai - (1./2.)*dt*Rhs_Cai_old;
                      Rsystem.get_vector("Cai_prev").set(dof_indices2[0],Rhs_Cai);
                      Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
                      Caup_new = Caup_old + (3./2.)*dt*Rhs_Caup - (1./2.)*dt*Rhs_Caup_old;
                      Rsystem.get_vector("Caup_prev").set(dof_indices2[0],Rhs_Caup);
                      Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
                      Carel_new = Carel_old + (3./2.)*dt*Rhs_Carel - (1./2.)*dt*Rhs_Carel_old;
                      Rsystem.get_vector("Carel_prev").set(dof_indices2[0],Rhs_Carel);
                      Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
                    }

                }
                */
                //else{
                      //Nernst Potential Calculation
                      E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
                      E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
                      E_K = (R*T/(1.*F))*log(K_o/Ki_old);

                      //BacNav Na current
                      I_BacNav = percBacNav*g_Na_BacNav*(mBN_old*mBN_old*mBN_old)*hBN_old*( vm_old - E_Na );

                      tau_m_BN = ((34.65)/(exp((vm_old + 43.47)/14.36) + exp((vm_old + 15.75)/(-.2351)))) + 1.66;//(4.2451./(exp((Vm(i) - -38.3561)/11.43877) + exp(-(Vm(i) - -34.4288)/1.)) + 0.14824);
                      tau_h_BN = tauhParam+(0.01 - 84.6609)*(1.0/(1.0+exp((-18.9945-vm_old)/2.4304))) + 0.01 + (12.47004 - 0.01)*(1.0/(1.0+exp((-40. - vm_old )/1.)));//((107.8)./(exp((Vm(i) + 27.15)./.1281) + exp((Vm(i) + 25.63)./(-25.19)))) + 9.593;
                      //84.6609 instead of 90. above to go back to original (only change is in tauh in bacnav)
                      m_inf_BN = 1.0/(1.0 + exp((-22.1573 - vm_old)/8.1769));//1./(1 + exp((-22.5-Vm(i))./2.704));
                      h_inf_BN = 1.-1.0/(1.0+exp((-76.7507-vm_old )/10.4215));//1./(1 + exp((77.05+Vm(i))./10.64));

                      //Fast Sodium Current
                      I_Na = decreaseFacgNa*g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
                      if(vm_old == -47.13){
                          alpha_m = 3.2;
                      }
                      else{
                          alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
                      }
                      beta_m = .08*exp(-vm_old/11.);
                      if(vm_old >= -40.){
                          alpha_h = 0.;
                      }
                      else{
                          alpha_h = .135*exp((vm_old + 80.)/(-6.8));
                      }
                      if(vm_old >= -40.){
                          beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
                      }
                      else{
                          beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
                      }
                      if(vm_old >= -40.){
                          alpha_j = 0.;
                      }
                      else{
                          alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
                      }
                      if(vm_old >= -40.){
                          beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
                      }
                      else{
                          beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
                      }
                      tau_m = 1./(alpha_m + beta_m);
                      tau_h = 1./(alpha_h + beta_h);
                      tau_j = 1./(alpha_j + beta_j);
                      m_inf = alpha_m/(alpha_m + beta_m);
                      h_inf = alpha_h/(alpha_h + beta_h);
                      j_inf = alpha_j/(alpha_j + beta_j);
                      //Time-independent K current
                      I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));

                      //Seemann et. al. 2010
                      if(SpiralBool > 0){
                        I_K1 = (2.1*g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
                      }
                      //Transient outward K current
                      I_to = increaseFacgto*g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      //Seemann et. al. 2010
                      if(SpiralBool > 0){
                        I_to = .35*increaseFacgto*g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
                      }

                      alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
                      oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
                      alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
                      beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
                      tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
                      oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
                      //Ultrarapid rapid delayed rectifier K current
                      g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
                      I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
                      alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
                      beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
                      tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
                      ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
                      alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
                      beta_ui = exp((vm_old - 158.)/16.);
                      tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
                      ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
                      //Rapid delayed outward rectifier K current
                      I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
                      alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
                      beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
                      tau_xr = 1./(alpha_xr + beta_xr);
                      xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
                      //Slow delayed outward rectifier K current
                      I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
                      alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
                      beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
                      tau_xs = 1./(2.*(alpha_xs + beta_xs));
                      xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
                      //L-type Ca current
                      I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      //Seemann et. al. 2010
                      if(SpiralBool > 0){
                        I_CaL = .35*g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
                      }

                      tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
                      d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
                      tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
                      f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
                      fCa_inf = 1./(1. + Cai_old/.00035);
                      //NaK pump current
                      f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
                      I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
                      //NaCa exchanger current
                      I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
                      //Background currents
                      I_bCa = g_bCa*(vm_old - E_Ca);
                      I_bNa = g_bNa*(vm_old - E_Na);
                      //Ca pump current
                      I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
                      //Ca release current from JSR
                      I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
                      F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
                      u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
                      v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
                      tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
                      w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
                      //Transfer current from NSR to JSR
                      I_tr = (Caup_old - Carel_old)/(tau_tr); 
                      //Ca uptake current by NSR
                      I_up = I_up_max/(1. + (K_up/Cai_old));
                      //Ca leak current by NSR
                      I_up_leak = I_up_max*Caup_old/(Ca_up_max);

                      //RHS of Concentrations
                      double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
                      double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
                      B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
                      B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
                      double Rhs_Cai = B1/B2;
                      double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
                      double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
                      //RHS of gating variables
                      double Rhs_v = (v_inf - v_old)/tau_v;
                      double Rhs_w = (w_inf - w_old)/tau_w;
                      double Rhs_xs = (xs_inf - xs_old)/tau_xs;
                      double Rhs_h = (h_inf - h_old)/tau_h;
                      double Rhs_xr = (xr_inf - xr_old)/tau_xr;
                      double Rhs_d = (d_inf - d_old)/tau_d;
                      double Rhs_oi = (oi_inf - oi_old)/tau_oi;
                      double Rhs_ui = (ui_inf - ui_old)/tau_ui;
                      double Rhs_m = (m_inf - m_old)/tau_m;
                      double Rhs_j = (j_inf - j_old)/tau_j;
                      double Rhs_f = (f_inf - f_old)/tau_f;
                      double Rhs_u = (u_inf - u_old)/tau_u;
                      double Rhs_oa = (oa_inf - oa_old)/tau_oa;
                      double Rhs_ua = (ua_inf - ua_old)/tau_ua;
                      double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

                      //BacNav
                      double Rhs_mBN = (m_inf_BN - mBN_old)/tau_m_BN;
                      double Rhs_hBN = (h_inf_BN - hBN_old)/tau_h_BN;

                      //Integrating over time
                      Nai_new = Nai_old + dt*Rhs_Nai;
                      Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
                      Ki_new = Ki_old + dt*Rhs_Ki;
                      Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
                      Cai_new = Cai_old + dt*Rhs_Cai;
                      Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
                      Caup_new = Caup_old + dt*Rhs_Caup;
                      Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
                      Carel_new = Carel_old + dt*Rhs_Carel;
                      Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
                      v_new = v_old + dt*Rhs_v;
                      Rsystem.solution -> set(dof_indices2[v_var],v_new);
                      w_new = w_old + dt*Rhs_w;
                      Rsystem.solution -> set(dof_indices2[w_var],w_new);
                      xs_new = xs_old + dt*Rhs_xs;
                      Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
                      h_new = h_old + dt*rateOfInact*Rhs_h;
                      Rsystem.solution -> set(dof_indices2[h_var],h_new);
                      xr_new = xr_old + dt*Rhs_xr;
                      Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
                      d_new = d_old + dt*Rhs_d;
                      Rsystem.solution -> set(dof_indices2[d_var],d_new);
                      oi_new = oi_old + dt*Rhs_oi;
                      Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
                      ui_new = ui_old + dt*Rhs_ui;
                      Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
                      m_new = m_old + dt*Rhs_m;
                      Rsystem.solution -> set(dof_indices2[m_var],m_new);
                      j_new = j_old + dt*Rhs_j;
                      Rsystem.solution -> set(dof_indices2[j_var],j_new);
                      f_new = f_old + dt*Rhs_f;
                      Rsystem.solution -> set(dof_indices2[f_var],f_new);
                      u_new = u_old + dt*Rhs_u;
                      Rsystem.solution -> set(dof_indices2[u_var],u_new);
                      oa_new = oa_old + dt*Rhs_oa;
                      Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
                      ua_new = ua_old + dt*Rhs_ua;
                      Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
                      fCa_new = fCa_old + dt*Rhs_fCa;
                      Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);

                      mBN_new = mBN_old + dt*Rhs_mBN;
                      Rsystem.solution -> set(dof_indices2[mBN_var],mBN_new);
                      hBN_new = hBN_old + dt*Rhs_hBN;
                      Rsystem.solution -> set(dof_indices2[hBN_var],hBN_new);
                //}

                double freact = 0.0;
                double I_ion = I_BacNav + I_Na + I_K1 + I_to + I_Kur + I_Kr + I_Ks + I_CaL + I_pCa + I_NaK + I_NaCa + I_bNa + I_bCa;
                double Istim = 0.0;
                double Istim2 = 0.0;

                if(SpiralBool == 0){
                    if( zDim == 0. ){
                      if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
                        

                        if(StimulusLocation.compare(0,4,"Cube") == 0){
                          if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx ){
                            Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                            //libMesh::out << x << std::endl;
                          }
                        }
                        else if(StimulusLocation.compare(0,6,"Points") == 0){
                          //READ POINTS
                          //LOOP OVER POINTS
                          //SET STIMULUS
                          
                          for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                           if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius ){
                              Istim = IstimV;
                              //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                              break;
                           }
                           else{}
                          }

                        }


                      }
                    }
                    else{
                      if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
                        


                        if(StimulusLocation.compare(0,4,"Cube") == 0){
                          if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx && z > stimulus_minz && z < stimulus_maxz ){
                            Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                            //libMesh::out << x << std::endl;
                          }
                        }
                        else if(StimulusLocation.compare(0,6,"Points") == 0){
                          //READ POINTS
                          //LOOP OVER POINTS
                          //SET STIMULUS
                          
                          for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                           if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius && std::abs(z - CoordsPointsZ[iter]) < StimulusLocation_Radius ){
                              Istim = IstimV;
                              //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                              break;
                           }
                           else{}
                          }

                        }



                      }
                    }
                    if(convergence_test){
                      freact = -kcubic*((u_old - u0)*(u_old - u1)*(u_old - u2));// - 1*r_new*(u_old-u0)) - Istim  - Istim2;
                    }
                    else{
                        if(StimPlace.compare(0,13,"Transmembrane") == 0){
                            if(ICstim){
                              freact = (-(I_ion));
                            }
                            else if(ICconditions){
                              freact = (-(I_ion));
                            }
                            else{
                              freact = (-(I_ion)) - Istim  - Istim2;
                            }
                            
                        }
                        else{
                            freact = (-(I_ion));
                        }
                    }
                }

                else if(SpiralBool == 1){
                    if( time < SpiralS2time + SpiralS2duration && time > SpiralS2time && x < (tissue_maxx + tissue_minx)/2. && y > (tissue_maxy + tissue_miny)/2. && z > (tissue_maxz + tissue_minz)/2. ){
                        Istim = SpiralS2strength;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                    //libMesh::out << Istim << std::endl;
                    }
                    freact = (-(I_ion)) - Istim  - Istim2;
                }

                else if(SpiralBool == 2){
                    if(time < IstimD &&  y > tissue_maxy - 0.25 ){
                    //if(time < IstimD &&  x < tissue_minx + 0.25 && y > 0 ){
                        Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    }
                    if(time < IstimD &&  x > tissue_maxx - 0.15 ){
                        Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    }
                    //SECOND STIMULUS
                    //if(time > 150. && time < 152. &&  x > tissue_maxx - 0.15 && y < 0. ){
                        //Istim = IstimV*.5;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    //}

                    freact = (-(I_ion)) - Istim  - Istim2;
                }



                
                //libMesh::out << freact << "    " << x << "    " << dof_indices[0] << std::endl;
                //libMesh::out << x << std::endl;


                Bsystem.get_vector("In").set(dof_indices[vm_var], freact);


             

          }
            

        }
        Bsystem.get_vector("In").close();
        Bsystem.get_vector("Inm1").close();

        Rsystem.solution -> close();
        Rsystem.update(); 

      //Rsystem.get_vector("Rvalues") = Rnew;
      

      }




      if(integrator.compare(0,4,"SBDF") != 0){

        TransientLinearImplicitSystem &Bsystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
        Bsystem.update();
     
        const unsigned int vm_var = Bsystem.variable_number ("V");
        const DofMap & dof_map = Bsystem.get_dof_map();


        for(const auto & node : mesh.local_node_ptr_range()){

            dof_map.dof_indices (node, dof_indices);
            dof_map2.dof_indices (node, dof_indices2);

          if(dof_indices.size() > 0){

            const Real x = (*node)(0);
            const Real y = (*node)(1);
            const Real z = (*node)(2);

              if(datatime.timestep == 1){
                v_old = 1.;
                w_old = .999;
                xs_old = 0.0187;
                h_old = .965;
                d_old = .000137;
                xr_old = .0000329;
                Nai_old = 11.2;
                Ki_old = 139.;
                Carel_old = 1.49;
                oi_old = .999;
                ui_old = .999;
                m_old = .00291;
                j_old = .978;
                f_old = .999;
                Cai_old = .000102;
                Caup_old = 1.49;
                oa_old = .0304;
                ua_old = .00496;
                fCa_old = .775;
                u_old = 0.0;
                vm_old = -81.2;//(*Bsystem.current_local_solution)  (dof_indices[vm_var]);

                //BacNav
                hBN_old = hBN_IC;
                mBN_old = mBN_IC;

                if(SpiralBool == 1){
                    if( x < (tissue_maxx + tissue_minx)/2. ){
                        v_old = 1.;
                        w_old = .9979;
                        xs_old = 0.1802;
                        h_old = .0895;
                        d_old = .001;
                        xr_old = .6142;
                        Nai_old = 11.1997;
                        Ki_old = 138.9995;
                        Carel_old = 1.4905;
                        oi_old = .3322;
                        ui_old = .9774;
                        m_old = .0357;
                        j_old = .0423;
                        f_old = .4576;
                        Cai_old = .00010355;
                        Caup_old = 1.4921;
                        oa_old = .0828;
                        ua_old = .0344;
                        fCa_old = .7717;
                        u_old = 0.0;
                        vm_old = -65.3708;

                        //BacNav
                        hBN_old = .0895;
                        mBN_old = .0357;
                    }
                }

                else if(SpiralBool == 2){

                    double a = .5;
                    double k = .09;
                    int thickness = 15;
                    int finalDim = thickness*1000; //15,000
                    double bathRegion = xDim/2.0 - tissue_maxx;
                    double factorAxis = 3.0*(xDim/2.-bathRegion);
                    phi[0] = 0.;
                    x1[0] = -a*exp(k*phi[0])*cos(phi[0]);
                    y1[0] = -a*exp(k*phi[0])*sin(phi[0]);
                    for(int i = 1; i < 1000; i++){
                        phi[i] = phi[i-1] + ((factorAxis*acos(-1) - phi[0])/1000.);
                        x1[i] = -a*exp(k*phi[i])*cos(phi[i]);
                        y1[i] = -a*exp(k*phi[i])*sin(phi[i]);
                    }
                    for(int i = 0; i < 1000; i++){
                        if( std::abs(x - x1[i]) <= threshSpiral && std::abs(y - y1[i]) <= threshSpiral){
                            //conditions in notes seen from paraview
                            /*
                             * initial conditions for multiple spirals to happen:
                                u - spiral of 1. and rest of 0
                                v - inverse of spiral of 1. and rest of 0
                                w - spiral of .8 and rest of 1.
                                s - same spiral of .1 and rest of 0
                            */
                            v_old = 1.;
                            w_old = .9979;
                            xs_old = 0.1802;
                            h_old = .0895;
                            d_old = .001;
                            xr_old = .6142;
                            Nai_old = 11.1997;
                            Ki_old = 138.9995;
                            Carel_old = 1.4905;
                            oi_old = .3322;
                            ui_old = .9774;
                            m_old = .0357;
                            j_old = .0423;
                            f_old = .4576;
                            Cai_old = .00010355;
                            Caup_old = 1.4921;
                            oa_old = .0828;
                            ua_old = .0344;
                            fCa_old = .7717;
                            u_old = 0.0;
                            vm_old = -65.3708;

                            //BacNav
                            hBN_old = .0895;
                            mBN_old = .0357;
                            //libMesh::out << "HERE --> x = " << x << ";     y = " << y << ";    with threshold of the spiral: " << threshSpiral << std::endl;
                        }
                    }
                }

              }
              else{
                  v_old = (*Rsystem.current_local_solution) (dof_indices2[v_var]);
                  w_old = (*Rsystem.current_local_solution) (dof_indices2[w_var]);        
                  xs_old = (*Rsystem.current_local_solution) (dof_indices2[xs_var]);
                  h_old = (*Rsystem.current_local_solution) (dof_indices2[h_var]);
                  xr_old = (*Rsystem.current_local_solution) (dof_indices2[xr_var]);        
                  d_old = (*Rsystem.current_local_solution) (dof_indices2[d_var]);
                  Nai_old = (*Rsystem.current_local_solution) (dof_indices2[Nai_var]);
                  Ki_old = (*Rsystem.current_local_solution) (dof_indices2[Ki_var]);        
                  Carel_old = (*Rsystem.current_local_solution) (dof_indices2[Carel_var]);     
                  oi_old = (*Rsystem.current_local_solution) (dof_indices2[oi_var]);
                  ui_old = (*Rsystem.current_local_solution) (dof_indices2[ui_var]);        
                  m_old = (*Rsystem.current_local_solution) (dof_indices2[m_var]);    
                  j_old = (*Rsystem.current_local_solution) (dof_indices2[j_var]);
                  f_old = (*Rsystem.current_local_solution) (dof_indices2[f_var]);        
                  Cai_old = (*Rsystem.current_local_solution) (dof_indices2[Cai_var]);     
                  Caup_old = (*Rsystem.current_local_solution) (dof_indices2[Caup_var]);
                  oa_old = (*Rsystem.current_local_solution) (dof_indices2[oa_var]);        
                  ua_old = (*Rsystem.current_local_solution) (dof_indices2[ua_var]);     
                  fCa_old = (*Rsystem.current_local_solution) (dof_indices2[fCa_var]);
                  u_old = (*Rsystem.current_local_solution) (dof_indices2[u_var]);           
                  vm_old = (*Bsystem.current_local_solution)  (dof_indices[vm_var]);

                  hBN_old = (*Rsystem.current_local_solution) (dof_indices2[hBN_var]);
                  mBN_old = (*Rsystem.current_local_solution) (dof_indices2[mBN_var]);
              }

              //Nernst Potential Calculation
              E_Ca = (R*T/(2.*F))*log(Ca_o/Cai_old);
              E_Na = (R*T/(1.*F))*log(Na_o/Nai_old);
              E_K = (R*T/(1.*F))*log(K_o/Ki_old);

              //BacNav Na current
              I_BacNav = percBacNav*g_Na_BacNav*(mBN_old*mBN_old*mBN_old)*hBN_old*( vm_old - E_Na );

              tau_m_BN = ((34.65)/(exp((vm_old + 43.47)/14.36) + exp((vm_old + 15.75)/(-.2351)))) + 1.66;//(4.2451./(exp((Vm(i) - -38.3561)/11.43877) + exp(-(Vm(i) - -34.4288)/1.)) + 0.14824);
              tau_h_BN = tauhParam+(0.01 - 84.6609)*(1.0/(1.0+exp((-18.9945-vm_old)/2.4304))) + 0.01 + (12.47004 - 0.01)*(1.0/(1.0+exp((-40. - vm_old )/1.)));//((107.8)./(exp((Vm(i) + 27.15)./.1281) + exp((Vm(i) + 25.63)./(-25.19)))) + 9.593;
              //84.6609 instead of 90. above to go back to original (only change is in tauh in bacnav)
              m_inf_BN = 1.0/(1.0 + exp((-22.1573 - vm_old)/8.1769));//1./(1 + exp((-22.5-Vm(i))./2.704));
              h_inf_BN = 1.-1.0/(1.0+exp((-76.7507-vm_old )/10.4215));//1./(1 + exp((77.05+Vm(i))./10.64));

              //Fast Sodium Current
              I_Na = decreaseFacgNa*g_Na*m_old*m_old*m_old*h_old*j_old*( vm_old - E_Na );
              if(vm_old == -47.13){
                  alpha_m = 3.2;
              }
              else{
                  alpha_m = .32*(vm_old + 47.13)/(1. - exp(-.1*(vm_old + 47.13)));
              }
              beta_m = .08*exp(-vm_old/11.);
              if(vm_old >= -40.){
                  alpha_h = 0.;
              }
              else{
                  alpha_h = .135*exp((vm_old + 80.)/(-6.8));
              }
              if(vm_old >= -40.){
                  beta_h = 1./(.13*(1. + exp((vm_old + 10.66)/-11.1)));
              }
              else{
                  beta_h = 3.56*exp(.079*vm_old) + 310000.*exp(.35*vm_old);
              }
              if(vm_old >= -40.){
                  alpha_j = 0.;
              }
              else{
                  alpha_j = (-127140.*exp(.2444*vm_old) - .00003474*exp(-.04391*vm_old))*((vm_old + 37.78)/(1. + exp(.311*(vm_old + 79.23))));
              }
              if(vm_old >= -40.){
                  beta_j = (.3*exp(-.0000002535*vm_old))/(1. + exp(-.1*(vm_old + 32.)));
              }
              else{
                  beta_j = (.1212*exp(-.01052*vm_old))/(1. + exp(-.1378*(vm_old + 40.14)));
              }
              tau_m = 1./(alpha_m + beta_m);
              tau_h = 1./(alpha_h + beta_h);
              tau_j = 1./(alpha_j + beta_j);
              m_inf = alpha_m/(alpha_m + beta_m);
              h_inf = alpha_h/(alpha_h + beta_h);
              j_inf = alpha_j/(alpha_j + beta_j);
              //Time-independent K current
              I_K1 = (g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));

              //Seemann et. al. 2010
              if(SpiralBool > 0){
                I_K1 = (2.1*g_K1*(vm_old - E_K))/(1. + exp(.07*(vm_old + 80.)));
              }
              //Transient outward K current
              I_to = increaseFacgto*g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
              //Seemann et. al. 2010
              if(SpiralBool > 0){
                I_to = .35*increaseFacgto*g_to*(oa_old*oa_old*oa_old)*(oi_old)*(vm_old - E_K);
              }

              alpha_oa = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
              beta_oa = .65/(2.5 + exp((vm_old + 82.)/17.));
              tau_oa = 1./(K_Q10*(alpha_oa + beta_oa));
              oa_inf = 1./(1. + exp((vm_old + 20.47)/-17.54));
              alpha_oi = 1./(18.53 + exp((vm_old + 113.7)/10.95));
              beta_oi = 1./(35.56 + exp((vm_old + 1.26)/-7.44));
              tau_oi = 1./(K_Q10*(alpha_oi + beta_oi));
              oi_inf = 1./(1. + exp((vm_old + 43.1)/5.3));
              //Ultrarapid rapid delayed rectifier K current
              g_Kur = .005 + (.05)/(1. + exp((vm_old - 15.)/-13.));
              I_Kur = g_Kur*(ua_old*ua_old*ua_old)*(ui_old)*(vm_old - E_K);
              alpha_ua = .65/(exp((vm_old + 10.)/-8.5) + exp((vm_old - 30.)/-59.));
              beta_ua = .65/(2.5 + exp((vm_old + 82.)/17.));
              tau_ua = 1./(K_Q10*(alpha_ua + beta_ua));
              ua_inf = 1./(1. + exp((vm_old + 30.3)/-9.6));
              alpha_ui = 1./(21. + exp((vm_old - 185.)/-28.));
              beta_ui = exp((vm_old - 158.)/16.);
              tau_ui = 1./(K_Q10*(alpha_ui + beta_ui));
              ui_inf = 1./(1. + exp((vm_old - 99.45)/27.48));
              //Rapid delayed outward rectifier K current
              I_Kr = (g_Kr*xr_old*(vm_old - E_K))/(1. + exp((vm_old + 15.)/22.4));
              alpha_xr = (.0003*(vm_old + 14.1))/(1. - exp((vm_old + 14.1)/-5.));
              beta_xr = (.000073898*(vm_old - 3.3328))/(-1. + exp((vm_old - 3.3328)/5.1237));
              tau_xr = 1./(alpha_xr + beta_xr);
              xr_inf = 1./(1. + exp((vm_old + 14.1)/-6.5));
              //Slow delayed outward rectifier K current
              I_Ks = g_Ks*(xs_old*xs_old)*(vm_old - E_K);
              alpha_xs = (.00004*(vm_old - 19.9))/(1. - exp((vm_old - 19.9)/-17.));
              beta_xs = (.000035*(vm_old - 19.9))/(-1. + exp((vm_old - 19.9)/9.));
              tau_xs = 1./(2.*(alpha_xs + beta_xs));
              xs_inf = 1./(sqrt(1. + exp((vm_old - 19.9)/-12.7)));
              //L-type Ca current
              I_CaL = g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
              //Seemann et. al. 2010
              if(SpiralBool > 0){
                I_CaL = .35*g_CaL*d_old*f_old*fCa_old*(vm_old - 65.);
              }

              tau_d = (1. - exp((vm_old + 10.)/-6.24))/(.035*(vm_old + 10.)*(1. + exp((vm_old+ 10.)/-6.24)));
              d_inf = 1./(1. + exp((vm_old + 10.)/-8.));
              tau_f = 9./(.02 + .0197*exp(-(.0337*.0337)*pow((vm_old + 10.),2)));
              f_inf = 1./(1. + exp((vm_old + 28.)/6.9));
              fCa_inf = 1./(1. + Cai_old/.00035);
              //NaK pump current
              f_NaK = 1./(1. + (.1245*exp(-.1*F*vm_old/(R*T))) + (.0365*sigma*exp(-F*vm_old/(R*T))));
              I_NaK = I_NaK_max*f_NaK*(1./(1. + (pow((K_m_Na_i/Nai_old),1.5))))*(K_o/(K_o + K_m_K_o));
              //NaCa exchanger current
              I_NaCa = ( I_NaCa_max*(exp(gamma*F*vm_old/(R*T))*(Nai_old*Nai_old*Nai_old)*(Ca_o) - exp((gamma-1.)*F*vm_old/(R*T))*(Na_o*Na_o*Na_o)*(Cai_old)) )/( (pow(K_m_Na,3.) + pow(Na_o,3.))*(K_m_Ca + Ca_o)*(1. + k_sat*exp((gamma - 1.)*F*vm_old/(R*T))) );
              //Background currents
              I_bCa = g_bCa*(vm_old - E_Ca);
              I_bNa = g_bNa*(vm_old - E_Na);
              //Ca pump current
              I_pCa = (I_pCa_max*Cai_old)/(.0005 + Cai_old);
              //Ca release current from JSR
              I_rel = k_rel*(u_old*u_old)*v_old*w_old*(Carel_old - Cai_old);
              F_n = (pow(10.,(-12)))*Vol_rel*I_rel - ((5.0e-13)/F)*(.5*I_CaL - .2*I_NaCa);
              u_inf = 1./(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
              tau_v = 1.91 + 2.09/(1. + exp((F_n - (3.4175e-13))/(-(13.67e-16))));
              v_inf = 1. - 1./(1. + exp((F_n - (6.835e-14))/(-(13.67e-16))));
              tau_w = 6.*(1. - exp((vm_old - 7.9)/-5.))/((1. + .3*exp((vm_old - 7.9)/-5.))*(vm_old - 7.9));
              w_inf = 1. - 1./(1. + exp((vm_old - 40.)/(-(17.))));
              //Transfer current from NSR to JSR
              I_tr = (Caup_old - Carel_old)/(tau_tr); 
              //Ca uptake current by NSR
              I_up = I_up_max/(1. + (K_up/Cai_old));
              //Ca leak current by NSR
              I_up_leak = I_up_max*Caup_old/(Ca_up_max);

              //RHS of Concentrations
              double Rhs_Nai = (-3.*I_NaK - 3.*I_NaCa - I_bNa - I_Na)/(F*Vol_i);
              double Rhs_Ki = (2.*I_NaK - I_K1 - I_to - I_Kur - I_Kr - I_Ks - I_bK)/(F*Vol_i);
              B1 = ((2.*I_NaCa - I_pCa - I_CaL - I_bCa)/(2.*F*Vol_i)) + ((Vol_up*(I_up_leak - I_up) + I_rel*Vol_rel)/Vol_i);
              B2 = 1. + ((Trpn_max*Km_Trpn)/(pow((Cai_old + Km_Trpn),2))) + ((Cmdn_max*Km_Cmdn)/(pow((Cai_old + Km_Cmdn),2)));
              double Rhs_Cai = B1/B2;
              double Rhs_Caup = I_up - I_up_leak -I_tr*Vol_rel/Vol_up;
              double Rhs_Carel = (I_tr - I_rel)/(1. + (((Csqn_max*Km_Csqn)/(pow((Carel_old + Km_Csqn),2)))));
              //RHS of gating variables
              double Rhs_v = (v_inf - v_old)/tau_v;
              double Rhs_w = (w_inf - w_old)/tau_w;
              double Rhs_xs = (xs_inf - xs_old)/tau_xs;
              double Rhs_h = (h_inf - h_old)/tau_h;
              double Rhs_xr = (xr_inf - xr_old)/tau_xr;
              double Rhs_d = (d_inf - d_old)/tau_d;
              double Rhs_oi = (oi_inf - oi_old)/tau_oi;
              double Rhs_ui = (ui_inf - ui_old)/tau_ui;
              double Rhs_m = (m_inf - m_old)/tau_m;
              double Rhs_j = (j_inf - j_old)/tau_j;
              double Rhs_f = (f_inf - f_old)/tau_f;
              double Rhs_u = (u_inf - u_old)/tau_u;
              double Rhs_oa = (oa_inf - oa_old)/tau_oa;
              double Rhs_ua = (ua_inf - ua_old)/tau_ua;
              double Rhs_fCa = (fCa_inf - fCa_old)/tau_fCa;

              //BacNav
              double Rhs_mBN = (m_inf_BN - mBN_old)/tau_m_BN;
              double Rhs_hBN = (h_inf_BN - hBN_old)/tau_h_BN;

              //Integrating over time
              Nai_new = Nai_old + dt*Rhs_Nai;
              Rsystem.solution -> set(dof_indices2[Nai_var],Nai_new);
              Ki_new = Ki_old + dt*Rhs_Ki;
              Rsystem.solution -> set(dof_indices2[Ki_var],Ki_new);
              Cai_new = Cai_old + dt*Rhs_Cai;
              Rsystem.solution -> set(dof_indices2[Cai_var],Cai_new);
              Caup_new = Caup_old + dt*Rhs_Caup;
              Rsystem.solution -> set(dof_indices2[Caup_var],Caup_new);
              Carel_new = Carel_old + dt*Rhs_Carel;
              Rsystem.solution -> set(dof_indices2[Carel_var],Carel_new);
              v_new = v_old + dt*Rhs_v;
              Rsystem.solution -> set(dof_indices2[v_var],v_new);
              w_new = w_old + dt*Rhs_w;
              Rsystem.solution -> set(dof_indices2[w_var],w_new);
              xs_new = xs_old + dt*Rhs_xs;
              Rsystem.solution -> set(dof_indices2[xs_var],xs_new);
              h_new = h_old + dt*rateOfInact*Rhs_h;
              Rsystem.solution -> set(dof_indices2[h_var],h_new);
              xr_new = xr_old + dt*Rhs_xr;
              Rsystem.solution -> set(dof_indices2[xr_var],xr_new);
              d_new = d_old + dt*Rhs_d;
              Rsystem.solution -> set(dof_indices2[d_var],d_new);
              oi_new = oi_old + dt*Rhs_oi;
              Rsystem.solution -> set(dof_indices2[oi_var],oi_new);
              ui_new = ui_old + dt*Rhs_ui;
              Rsystem.solution -> set(dof_indices2[ui_var],ui_new);
              m_new = m_old + dt*Rhs_m;
              Rsystem.solution -> set(dof_indices2[m_var],m_new);
              j_new = j_old + dt*Rhs_j;
              Rsystem.solution -> set(dof_indices2[j_var],j_new);
              f_new = f_old + dt*Rhs_f;
              Rsystem.solution -> set(dof_indices2[f_var],f_new);
              u_new = u_old + dt*Rhs_u;
              Rsystem.solution -> set(dof_indices2[u_var],u_new);
              oa_new = oa_old + dt*Rhs_oa;
              Rsystem.solution -> set(dof_indices2[oa_var],oa_new);
              ua_new = ua_old + dt*Rhs_ua;
              Rsystem.solution -> set(dof_indices2[ua_var],ua_new);
              fCa_new = fCa_old + dt*Rhs_fCa;
              Rsystem.solution -> set(dof_indices2[fCa_var],fCa_new);

              mBN_new = mBN_old + dt*Rhs_mBN;
              Rsystem.solution -> set(dof_indices2[mBN_var],mBN_new);
              hBN_new = hBN_old + dt*Rhs_hBN;
              Rsystem.solution -> set(dof_indices2[hBN_var],hBN_new);
        

              double freact = 0.0;
              double I_ion = I_BacNav + I_Na + I_K1 + I_to + I_Kur + I_Kr + I_Ks + I_CaL + I_pCa + I_NaK + I_NaCa + I_bNa + I_bCa;
              double Istim = 0.0;
              double Istim2 = 0.0;

              if(SpiralBool == 0){
                    if( zDim == 0. ){
                      if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
                        

                        if(StimulusLocation.compare(0,4,"Cube") == 0){
                          if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx ){
                            Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                            //libMesh::out << x << std::endl;
                          }
                        }
                        else if(StimulusLocation.compare(0,6,"Points") == 0){
                          //READ POINTS
                          //LOOP OVER POINTS
                          //SET STIMULUS
                          
                          for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                           if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius ){
                              Istim = IstimV;
                              //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                              break;
                           }
                           else{}
                          }

                        }


                      }
                    }
                    else{
                      if( time < (stimulus_start_time + IstimD) && time > stimulus_start_time ){
                        


                        if(StimulusLocation.compare(0,4,"Cube") == 0){
                          if(  y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx && z > stimulus_minz && z < stimulus_maxz ){
                            Istim = IstimV;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                            //libMesh::out << x << std::endl;
                          }
                        }
                        else if(StimulusLocation.compare(0,6,"Points") == 0){
                          //READ POINTS
                          //LOOP OVER POINTS
                          //SET STIMULUS
                          
                          for(int iter = 0; iter < CoordsPointsX.size(); iter++){
                           if( std::abs(x - CoordsPointsX[iter]) < StimulusLocation_Radius && std::abs(y - CoordsPointsY[iter]) < StimulusLocation_Radius && std::abs(z - CoordsPointsZ[iter]) < StimulusLocation_Radius ){
                              Istim = IstimV;
                              //libMesh::out << "Stimulating at time: " << time << " ms - Using the points to stimulate: File: " << StimulusLocation_File << " found successfully..." << std::endl;
                              break;
                           }
                           else{}
                          }

                        }



                      }
                    }
                    if(convergence_test){
                      freact = -kcubic*((u_old - u0)*(u_old - u1)*(u_old - u2));// - 1*r_new*(u_old-u0)) - Istim  - Istim2;
                    }
                    else{
                        if(StimPlace.compare(0,13,"Transmembrane") == 0){
                            if(ICstim){
                              freact = (-(I_ion));
                            }
                            else if(ICconditions){
                              freact = (-(I_ion));
                            }
                            else{
                              freact = (-(I_ion)) - Istim  - Istim2;
                            }
                        }
                        else{
                            freact = (-(I_ion));
                        }
                    }
                }


                else if(SpiralBool == 1){
                    if( time < SpiralS2time + SpiralS2duration && time > SpiralS2time && x < (tissue_maxx + tissue_minx)/2. && y > (tissue_maxy + tissue_miny)/2. && z > (tissue_maxz + tissue_minz)/2. ){
                        Istim = SpiralS2strength;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                    //libMesh::out << Istim << std::endl;
                    }
                    freact = (-(I_ion)) - Istim  - Istim2;
                }

                else if(SpiralBool == 2){
                    //if(time < IstimD &&  x < tissue_minx + 0.5 ){
                        //Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    //}
                    if(time < IstimD &&  y > tissue_maxy - 0.25 ){
                    //if(time < IstimD &&  x < tissue_minx + 0.25 && y > 0 ){
                        Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    }
                    if(time < IstimD &&  x > tissue_maxx - 0.15 ){
                        Istim = IstimV*1.;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    }
                    //SECOND STIMULUS
                    //if(time > 150. && time < 152. &&  x > tissue_maxx - 0.15 && y < 0. ){
                        //Istim = IstimV*.5;//with TIME < 1.2, -.08 is the minimum current stimulus to initiate AP propagation
                        //libMesh::out << Istim << std::endl;
                    //}

                    freact = (-(I_ion)) - Istim  - Istim2;
                }

              //libMesh::out << freact << "    " << x << "    " << dof_indices[0] << std::endl;
              //libMesh::out << x << std::endl;


              Bsystem.get_vector("In").set(dof_indices[vm_var], freact);

            
          }


        }
        Bsystem.get_vector("In").close();
        //Rsystem.get_vector("Rvalues") = Rnew;

        Rsystem.solution -> close();
        Rsystem.update();

      }


  }

  //end of Courtemanche

//end of function
}








void init_cd_exact (EquationSystems & es, const double xDim, std::string integrator,  const GetPot &data)
{
      libMesh::out << "Initial conditions..." << std::endl;
      MeshBase & mesh = es.get_mesh();
      const DofMap & dof_map2 = es.get_system("Recovery").get_dof_map();
      const unsigned int dim = mesh.mesh_dimension();
      QGauss qrule (dim, FIFTH);
      double u_h, u_exact, ue_h, ue_exact;
      double v, w, s;
      std::vector<dof_id_type> dof_indices;
      std::vector<dof_id_type> dof_indices2;
      std::vector<dof_id_type> dof_indicesP;

      double stimulus_maxx = data("stimulus_maxx", .5);
      double stimulus_maxy = data("stimulus_maxy", .5);
      double stimulus_minx = data("stimulus_minx", -.5);
      double stimulus_miny = data("stimulus_miny", .85);
      double stimulus_maxz = data("stimulus_maxz", 0.0);
      double stimulus_minz = data("stimulus_minz", 0.0);
      double zDim = data("maxz", 0.) - data("minz", 0.);
      //double stimulus_start_time = data("stimulus_start_time", 0.);
      bool convergence_test = data("convergence_test", false);

      bool monodomainType = data("monodomainType", false);
      bool ICstim = data("ICstim", false);

      bool ICconditions = data("ICconditions", false);
      int ICsave_iter_start = data("ICsave_iter_start", 0);
      double ICtime = data("ICtime", 0.);


        std::string ICV_coordinates = data("ICV_coordinates", "fibrosisCoords.csv");
        std::string ICVe_coordinates = data("ICVe_coordinates", "IC_Pindices_SBDF1.txt");
        std::string ICV_values = data("ICV_values", "IC_Pindices_SBDF1.txt");
        std::string ICVe_values = data("ICVe_values", "IC_Pindices_SBDF1.txt");
        std::string ICv_values = data("ICv_values", "IC_Pindices_SBDF1.txt");
        std::string ICw_values = data("ICw_values", "IC_Pindices_SBDF1.txt");
        std::string ICs_values = data("ICs_values", "IC_Pindices_SBDF1.txt");

        std::string cellModel = data("cellModel", "Fenton-Karma");


        std::vector<int> ICV_coordinates_vec, ICVe_coordinates_vec;
        std::vector<double> ICV_values_vec, ICVe_values_vec, ICv_values_vec, ICw_values_vec, ICs_values_vec;
        ifstream inputFile(ICV_coordinates);
        ifstream inputFile2(ICVe_coordinates);
        ifstream inputFile3(ICV_values);
        ifstream inputFile4(ICVe_values);
        ifstream inputFile5(ICv_values);
        ifstream inputFile6(ICw_values);
        ifstream inputFile7(ICs_values);



        if(ICconditions){
            // V coordinates  
            if (inputFile) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile >> value ) {
                  ICV_coordinates_vec.push_back( (int) value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "V coordinates file for initial conditions: " << ICV_coordinates << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
            // Ve coordinates  
            if (inputFile2) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile2 >> value ) {
                  ICVe_coordinates_vec.push_back( (int) value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "Ve coordinates file for initial conditions: " << ICVe_coordinates << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
            // V values  
            if (inputFile3) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile3 >> value ) {
                  ICV_values_vec.push_back(value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "V values file for initial conditions: " << ICV_values << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
            // Ve values  
            if (inputFile4) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile4 >> value ) {
                  ICVe_values_vec.push_back(value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "Ve values file for initial conditions: " << ICVe_values << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
            // v values  
            if (inputFile5) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile5 >> value ) {
                  ICv_values_vec.push_back(value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "v values file for initial conditions: " << ICv_values << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
            // w values  
            if (inputFile6) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile6 >> value ) {
                  ICw_values_vec.push_back(value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "w values file for initial conditions: " << ICw_values << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
            // s values  
            if (inputFile7) {        
              double value;
              // read the elements in the file into a vector  
              while ( inputFile7 >> value ) {
                  ICs_values_vec.push_back(value);
                  //libMesh::out << value << std::endl;
              }
              libMesh::out << "s values file for initial conditions: " << ICs_values << " found successfully..." << std::endl;
            }
            else{
              libMesh::out << "Could not open the file" << std::endl;
            }
        }



          //libMesh::out << ICV_coordinates_vec[20] <<std::endl;




      if(integrator.compare(0,4,"SBDF") == 0){
            TransientLinearImplicitSystem & Bsystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
            Bsystem.update();

           

            auto femSolu = es.get_system("parabolic").variable_number("V");
            const DofMap & dof_map = es.get_system("parabolic").get_dof_map();
            FEType fe_type = dof_map.variable_type(femSolu);
            std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
            fe->attach_quadrature_rule (&qrule);
            const std::vector<Point> & q_point = fe->get_xyz();
            const std::vector<Real> & JxW = fe->get_JxW();
            const std::vector<std::vector<Real>> & phi = fe->get_phi();

            for(const auto & node : mesh.local_node_ptr_range()){

              dof_map.dof_indices (node, dof_indices);

              const Real x = (*node)(0);
              const Real y = (*node)(1);
              const Real z = (*node)(2);

              if(cellModel.compare(0,12,"Courtemanche") == 0){

                if(dof_indices.size() > 0){

                  u_exact = -81.2;
                  Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                  //ue_exact = exact_solutionV_all(x, y, 3, 0.0, 0.0, xDim);
                  //Bsystem.solution -> set(dof_indices[femSolu2],ue_exact);

                }
                else{
                  //ue_exact = exact_solutionV_all(x, y, 4, 0.0, 0.0, xDim);
                  //Bsystem.solution -> set(dof_indices[femSolu2],ue_exact);
                }

              }
              else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){

                if(dof_indices.size() > 0){

                  u_exact = -85.;
                  Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                  //ue_exact = exact_solutionV_all(x, y, 3, 0.0, 0.0, xDim);
                  //Bsystem.solution -> set(dof_indices[femSolu2],ue_exact);

                }
                else{
                  //ue_exact = exact_solutionV_all(x, y, 4, 0.0, 0.0, xDim);
                  //Bsystem.solution -> set(dof_indices[femSolu2],ue_exact);
                }

              }


              if(dof_indices.size() > 0){
                if(convergence_test){
                    u_exact = exact_solutionV_all(x, y, 2, 0.0, 0.0, xDim);
                    Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                }
                if(ICstim){

                    if( zDim == 0. ){
                      if(y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx){
                        u_exact = 1.;
                        Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                      }
                      else{
                        u_exact = 0.;
                        Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                      }

                      

                    }
                    else{
                        if(y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx && z > stimulus_minz && z < stimulus_maxz){
                          u_exact = 1.;
                          Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                        }
                        else{
                          u_exact = 0.;
                          Bsystem.solution -> set(dof_indices[femSolu],u_exact);
                        }
                        

                    }

                }
              }
              
            }

          Bsystem.solution -> close();
      }
      else{
            libMesh::TransientLinearImplicitSystem & Psystem = es.get_system < libMesh::TransientLinearImplicitSystem > ("parabolic");
            libMesh::TransientExplicitSystem & Rsystem = es.get_system < libMesh::TransientExplicitSystem > ("Recovery");
            Psystem.update();
            Rsystem.update();
            auto femSolu = es.get_system("parabolic").variable_number("V");
            //auto femSoluv = es.get_system("Recovery").variable_number("v");
            //auto femSoluw = es.get_system("Recovery").variable_number("w");
            //auto femSolus = es.get_system("Recovery").variable_number("s");

            const DofMap & dof_mapP = es.get_system("parabolic").get_dof_map();
            FEType fe_type = dof_mapP.variable_type(femSolu);
            std::unique_ptr<FEBase> fe (FEBase::build(dim, fe_type));
            fe->attach_quadrature_rule (&qrule);
            const std::vector<Point> & q_point = fe->get_xyz();
            const std::vector<Real> & JxW = fe->get_JxW();
            const std::vector<std::vector<Real>> & phi = fe->get_phi();
            libMesh::out << "Inside Initial Conditions component..." << std::endl;
            for(const auto & node : mesh.local_node_ptr_range()){

                dof_mapP.dof_indices (node, dof_indicesP);

                const Real x = (*node)(0);
                const Real y = (*node)(1);
                const Real z = (*node)(2);

                auto nodeid = node -> id();
                //libMesh::out << "the node id is: " << nodeid << std::endl;

                //libMesh::out << "CHECKING FOR ERROR" << std::endl;
                if(dof_indicesP.size() > 0){
                  if(cellModel.compare(0,12,"Courtemanche") == 0){
                      u_exact = -81.2;
                      Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                      //ue_exact = exact_solutionV_all(x, y, 3, 0.0, 0.0, xDim);
                      //Bsystem.solution -> set(dof_indices[femSolu2],ue_exact);
                  }
                  else if(cellModel.compare(0,14,"Fenton-KarmaOG") == 0){
                      u_exact = -85.;
                      Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                      //ue_exact = exact_solutionV_all(x, y, 3, 0.0, 0.0, xDim);
                      //Bsystem.solution -> set(dof_indices[femSolu2],ue_exact);
                  }

                  if(convergence_test){
                      //libMesh::out << "INSIDE" << std::endl;
                      u_exact = exact_solutionV_all(x, y, 2, 0.0, 0.0, xDim);
                      Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                  }
                  if(ICstim){
                        if( zDim == 0. ){
                          if(y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx){
                            u_exact = 1.;
                            Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                          }
                          else{
                            u_exact = 0.;
                            Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                          }
                        }
                        else{
                            if(y > stimulus_miny && y < stimulus_maxy && x > stimulus_minx && x < stimulus_maxx && z > stimulus_minz && z < stimulus_maxz){
                              u_exact = 1.;
                              Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                            }
                            else{
                              u_exact = 0.;
                              Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                            }
                        }
                  }
                }

                /*
                if(ICconditions){

                  if(dof_indicesP.size() == dof_indices.size()){
                      //libMesh::out << "INSIDE" << std::endl;

                    auto Vcoordscheck = std::find(ICV_coordinates_vec.begin(), ICV_coordinates_vec.end(), nodeid);
                    if(ICV_coordinates_vec.end() != Vcoordscheck){
                      int indexCoords = Vcoordscheck - ICV_coordinates_vec.begin();


                      u_exact = ICV_values_vec[indexCoords];
                      Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                      v = ICv_values_vec[indexCoords];
                      Rsystem.solution -> set(dof_indices2[femSoluv], v);
                      w = ICw_values_vec[indexCoords];
                      Rsystem.solution -> set(dof_indices2[femSoluw], w);
                      s = ICs_values_vec[indexCoords];
                      Rsystem.solution -> set(dof_indices2[femSolus], s);

                    }

                    for(int coordCheck=0; coordCheck < ICV_coordinates_vec.length(): coordCheck++){
                        if(ICV_coordinates_vec[coordCheck] == nodeid){
                          //something happens
                          u_exact = ICV_values_vec[nodeid];
                          Psystem.solution -> set(dof_indicesP[femSolu],u_exact);
                          v = ICv_values_vec[nodeid];
                          Rsystem.solution -> set(dof_indices2[femSoluv], v);
                          w = ICw_values_vec[nodeid];
                          Rsystem.solution -> set(dof_indices2[femSoluw], w);
                          s = ICs_values_vec[nodeid];
                          Rsystem.solution -> set(dof_indices2[femSolus], s);
                        }
                      }



                    auto Vecoordscheck = std::find(ICVe_coordinates_vec.begin(), ICVe_coordinates_vec.end(), nodeid);
                    if(ICVe_coordinates_vec.end() != Vecoordscheck){
                      int indexCoords = Vecoordscheck - ICVe_coordinates_vec.begin();

                      ue_exact = ICVe_values_vec[indexCoords];
                      Esystem.solution -> set(dof_indices[femSolu2],ue_exact);
                    }

                  }
                  else{
                    
                    auto Vecoordscheck = std::find(ICVe_coordinates_vec.begin(), ICVe_coordinates_vec.end(), nodeid);
                    if(ICVe_coordinates_vec.end() != Vecoordscheck){
                      int indexCoords = Vecoordscheck - ICVe_coordinates_vec.begin();

                      ue_exact = ICVe_values_vec[indexCoords];
                      Esystem.solution -> set(dof_indices[femSolu2],ue_exact);
                    }

                  }
                    
                    
                }
                */



              
            }

          Psystem.solution -> close();
      }
      libMesh::out << "End of Initial conditions..." << std::endl;

}


void ForcingTermConvergence (libMesh::EquationSystems &es, const double dt, const double Beta, const double Cm, const double SigmaSI, const double SigmaSE, const double SigmaBath, const double SigmaTorso, const double CurTime, const double xDim, const double v0, const double v1, const double v2, const double kcubic, const std::string integrator )
{

  const MeshBase & mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  QGauss qrule (dim,  FIFTH);
  QGauss qface (dim-1, SECOND);
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_V;

  if(integrator.compare(0,4,"SBDF") == 0){

      TransientLinearImplicitSystem & Bsystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
      FEType fe_type_V = Bsystem.variable_type(0);
      std::unique_ptr<FEBase> fe      (FEBase::build(dim, fe_type_V));
      std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type_V));
      fe->attach_quadrature_rule      (&qrule);
      fe_face->attach_quadrature_rule (&qface);
      const std::vector<Real> & JxW      = fe->get_JxW();
      const std::vector<Real> & JxW_face = fe_face->get_JxW();
      const std::vector<std::vector<Real>> & phi = fe->get_phi();
      const std::vector<std::vector<Real>> & psi = fe_face->get_phi();
      const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
      const std::vector<Point> & q_point = fe->get_xyz();
      const std::vector<Point> & qface_points = fe_face->get_xyz();
      const DofMap & dof_map = Bsystem.get_dof_map();

      DenseVector<Number> FV;
      
      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          dof_map.dof_indices (elem, dof_indices_V);

          const unsigned int n_V_dofsP = dof_indices_V.size();

          fe->reinit (elem);

          FV.resize (n_V_dofsP);

            for (unsigned int qp=0; qp<qrule.n_points(); qp++)
            {
              for (unsigned int i=0; i<n_V_dofsP; i++){
                FV(i) += JxW[qp]*((    CalculateF( q_point[qp](0), q_point[qp](1), 0, CurTime, 0.0, dt, Beta, Cm, SigmaSI, SigmaSE, SigmaBath, SigmaTorso, xDim, v0, v1, v2, kcubic)      )*(phi[i][qp]));
              }
            }
            {
                    //BOUNDARIES
                              for (auto s : elem->side_index_range()){

                                //libMesh::out << elem << std::endl;
                                            if (elem->neighbor_ptr(s) == nullptr)
                                              {
                                              //libMesh::out << "SOMETHING        " << s << std::endl;

                                                std::unique_ptr<const Elem> side (elem->build_side_ptr(s,false));
                                                //libMesh::out << side.get() << std::endl;
                                                // Loop over the nodes on the side.
                                                for (auto ns : side->node_index_range())
                                                  {
                                                    // The location on the boundary of the current
                                                    // node.
                                                  //libMesh::out << ns << std::endl;

                                                    const Real xf = side->point(ns)(0);
                                                    const Real yf = side->point(ns)(1);

                                                    // The penalty value.  \f$ \frac{1}{\epsilon \f$
                                                    const Real penalty = 1.e10;

                                                    // The boundary values.

                                                    // Set v = 0 everywhere
                                                    // Find the node on the element matching this node on
                                                    // the side.  That defined where in the element matrix
                                                    // the boundary condition will be applied.
                                                    //libMesh::out << "  BEFORE CONDITIONS " << std::endl;

                                                    //libMesh::out << "  AFTER CONDITIONS " << std::endl;

                                                  } // end face node loop
                                              }
                              }// end if (elem->neighbor(side) == nullptr)

             }

            Bsystem.get_vector("ForcingV").add_vector    (FV, dof_indices_V);
        }

          Bsystem.get_vector("ForcingV").close();


  }

  else{
      TransientLinearImplicitSystem & Psystem = es.get_system <TransientLinearImplicitSystem> ("parabolic");
      FEType fe_type_V = Psystem.variable_type(0);
      std::unique_ptr<FEBase> fe      (FEBase::build(dim, fe_type_V));
      std::unique_ptr<FEBase> fe_face (FEBase::build(dim, fe_type_V));
      fe->attach_quadrature_rule      (&qrule);
      fe_face->attach_quadrature_rule (&qface);
      const std::vector<Real> & JxW      = fe->get_JxW();
      const std::vector<Real> & JxW_face = fe_face->get_JxW();
      const std::vector<std::vector<Real>> & phi = fe->get_phi();
      const std::vector<std::vector<Real>> & psi = fe_face->get_phi();
      const std::vector<std::vector<RealGradient>> & dphi = fe->get_dphi();
      const std::vector<Point> & q_point = fe->get_xyz();
      const std::vector<Point> & qface_points = fe_face->get_xyz();
      const DofMap & dof_mapP = Psystem.get_dof_map();

      DenseVector<Number> FV;

      for (const auto & elem : mesh.active_local_element_ptr_range())
        {
          dof_mapP.dof_indices (elem, dof_indices_V);

          const unsigned int n_V_dofsP = dof_indices_V.size();

          fe->reinit (elem);

          FV.resize (n_V_dofsP);

            for (unsigned int qp=0; qp<qrule.n_points(); qp++)
            {
              for (unsigned int i=0; i<n_V_dofsP; i++){
                FV(i) += JxW[qp]*((    CalculateF( q_point[qp](0), q_point[qp](1), 0, CurTime, 0.0, dt, Beta, Cm, SigmaSI, SigmaSE, SigmaBath, SigmaTorso, xDim, v0, v1, v2, kcubic)      )*(phi[i][qp]));
              }
            }
            {
                    //BOUNDARIES
                              for (auto s : elem->side_index_range()){

                                //libMesh::out << elem << std::endl;
                                            if (elem->neighbor_ptr(s) == nullptr)
                                              {
                                              //libMesh::out << "SOMETHING        " << s << std::endl;

                                                std::unique_ptr<const Elem> side (elem->build_side_ptr(s,false));
                                                //libMesh::out << side.get() << std::endl;
                                                // Loop over the nodes on the side.
                                                for (auto ns : side->node_index_range())
                                                  {
                                                    // The location on the boundary of the current
                                                    // node.
                                                  //libMesh::out << ns << std::endl;

                                                    const Real xf = side->point(ns)(0);
                                                    const Real yf = side->point(ns)(1);

                                                    // The penalty value.  \f$ \frac{1}{\epsilon \f$
                                                    const Real penalty = 1.e10;

                                                    // The boundary values.

                                                    // Set v = 0 everywhere
                                                    // Find the node on the element matching this node on
                                                    // the side.  That defined where in the element matrix
                                                    // the boundary condition will be applied.
                                                    //libMesh::out << "  BEFORE CONDITIONS " << std::endl;

                                                    //libMesh::out << "  AFTER CONDITIONS " << std::endl;

                                                  } // end face node loop
                                              }
                              }// end if (elem->neighbor(side) == nullptr)

             }

            Psystem.get_vector("ForcingV").add_vector    (FV, dof_indices_V);
        }

          Psystem.get_vector("ForcingV").close();


    }


}




