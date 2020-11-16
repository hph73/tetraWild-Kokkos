#include <Eigen/Dense>

int runfTetWild(std::string input_path, double ideal_edge_length_rel, double eps_rel, double stop_energy, int max_its, 
  bool skip_simplify, bool is_smooth_open_boundary, bool use_floodfill, bool use_general_wn,
  Eigen::MatrixXd &verticesInput, Eigen::MatrixXi &facetsInput, std::string outputPath,
  Eigen::MatrixXd &points, Eigen::MatrixXi &elements, Eigen::VectorXi &material, Eigen::MatrixXd &V_sf, Eigen::MatrixXi &F_sf);

#ifdef ENABLE_MPI
#include <mpi.h>
#endif
int main(int argc, char *argv[]){
  // generate tetra elements in the bounding box
  double ideal_edge_length_rel = 1.0/20.0;
  double eps_rel = 1e-3;
  double stop_energy = 10.0;
  int max_its = 80;
  bool skip_simplify = false;
  bool use_floodfill = false;
  bool use_general_wn = false;
  bool smooth_open_boundary = false;
  Eigen::MatrixXd vertexInput; 
  Eigen::MatrixXi facetInput;
  Eigen::MatrixXd pointsTemp;
  Eigen::MatrixXi elementsTemp;
  Eigen::VectorXi materialTemp;
  Eigen::MatrixXd vertexTemp;
  Eigen::MatrixXi facetTemp;
  runfTetWild(
    "./cube10.stl", ideal_edge_length_rel, eps_rel, stop_energy, max_its, skip_simplify, smooth_open_boundary, use_floodfill, use_general_wn,
    vertexInput, facetInput, "./output.msh", 
    pointsTemp, elementsTemp, materialTemp, vertexTemp, facetTemp);
  return 0;
}
