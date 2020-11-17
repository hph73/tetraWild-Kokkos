#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <floattetwild/Mesh.hpp>
#include <floattetwild/MeshIO.hpp>
#include <floattetwild/FloatTetDelaunay.h>
#include <floattetwild/LocalOperations.h>
#include <floattetwild/MeshImprovement.h>
#include <floattetwild/Simplification.h>
#include <floattetwild/AABBWrapper.h>
#include <floattetwild/Statistics.h>
#include <floattetwild/TriangleInsertion.h>
#include <floattetwild/CSGTreeParser.hpp>

#include <floattetwild/Logger.hpp>
#include <Eigen/Dense>

#include <igl/Timer.h>
#include <igl/write_triangle_mesh.h>

#include <geogram/basic/logger.h>
#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <geogram/mesh/mesh.h>

#include <bitset>

using namespace floatTetWild;
using namespace Eigen;

#include <geogram/basic/common.h>
#include <geogram/basic/numeric.h>
#include <geogram/basic/geometry.h>
#include <floattetwild/Predicates.hpp>

#include <geogram/mesh/mesh_AABB.h>
#include <floattetwild/MshLoader.h>

// params.input_path "Input surface mesh INPUT in .off/.obj/.stl/.ply format. (string, required)"
// params.output_path "Output tetmesh OUTPUT in .msh format. (string, optional, default: input_file+postfix+'.msh')"
// csg_file "json file containg a csg tree"
// boolean_op "Boolean operation: 0: union, 1: intersection, 2: difference."
// params.ideal_edge_length_rel "ideal_edge_length = diag_of_bbox * L. (double, optional, default: 0.05)"
// params.eps_rel "epsilon = diag_of_bbox * EPS. (double, optional, default: 1e-3)"
// params.stop_energy "Stop optimization when max energy is lower than this."
// params.log_path "Log info to given file."
// params.log_level "Log level (0 = most verbose, 6 = off)."
// params.is_quiet "Mute console output. (optional)"
// skip_simplify "skip preprocessing."
// nobinary "export meshes as ascii"
// nocolor "don't export color"
// params.smooth_open_boundary, "Smooth the open boundary."
// params.manifold_surface, "Force the output to be manifold."
// params.coarsen, "Coarsen the output as much as possible."
// params.disable_filtering, "Disable filtering out outside elements."
// params.use_floodfill, "Use flood-fill to extract interior volume."
// params.use_general_wn, "Use general winding number."
// params.use_input_for_wn, "Use input surface for winding number."
int runfTetWild(std::string input_path, double ideal_edge_length_rel, double eps_rel, double stop_energy, int max_its, 
  bool skip_simplify, bool is_smooth_open_boundary, bool use_floodfill, bool use_general_wn,
  Eigen::MatrixXd &verticesInput, Eigen::MatrixXi &facetsInput, std::string outputPath,
  Eigen::MatrixXd &points, Eigen::MatrixXi &elements, Eigen::VectorXi &material, Eigen::MatrixXd &V_sf, Eigen::MatrixXi &F_sf) {

  GEO::initialize();

  // Import standard command line arguments, and custom ones
  GEO::CmdLine::import_arg_group("standard");
  GEO::CmdLine::import_arg_group("pre");
  GEO::CmdLine::import_arg_group("algo");

  unsigned int max_threads = std::numeric_limits<unsigned int>::max();

  Mesh mesh;
  Parameters &params = mesh.params;

  bool nobinary = false;
  bool nocolor = false;
  params.ideal_edge_length_rel = ideal_edge_length_rel;
  params.eps_rel = eps_rel;
  params.stop_energy = stop_energy;
  params.max_its = max_its;
  params.input_path = input_path;
  params.output_path = outputPath;
  params.disable_filtering = false;
  params.use_floodfill = use_floodfill;
  params.use_general_wn = use_general_wn;
  params.smooth_open_boundary = is_smooth_open_boundary;
  int boolean_op = -1;
  std::string csg_file="";
  std::string background_mesh = "";

  #ifdef USE_TBB
    const size_t MB = 1024 * 1024;
    const size_t stack_size = 64 * MB;
    unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    num_threads = std::min(max_threads, num_threads);
    params.num_threads = num_threads;
    std::cout << "TBB threads " << num_threads << std::endl;
    tbb::task_scheduler_init scheduler(num_threads, stack_size);
  #endif

  Logger::init(!params.is_quiet, params.log_path);
  params.log_level = std::max(0, std::min(6, params.log_level));
  spdlog::set_level(static_cast<spdlog::level::level_enum>(params.log_level));
  spdlog::flush_every(std::chrono::seconds(3));

  if (params.output_path.empty())
    params.output_path = params.input_path;
  if (params.log_path.empty())
    params.log_path = params.output_path;

  std::string output_mesh_name = params.output_path;
  if (params.output_path.size() > 3
    && params.output_path.substr(params.output_path.size() - 3, params.output_path.size()) == "msh")
    output_mesh_name = params.output_path;
  else if (params.output_path.size() > 4
    && params.output_path.substr(params.output_path.size() - 4, params.output_path.size()) == "mesh")
    output_mesh_name = params.output_path;
  else
    output_mesh_name = params.output_path + "_" + params.postfix + ".msh";

  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ///set input tage
  std::vector<Vector3> input_vertices;
  std::vector<Vector3i> input_faces;
  std::vector<int> input_tags;

  input_vertices.resize(verticesInput.rows());
  for (size_t i = 0; i < verticesInput.rows(); i++)
    input_vertices[i] << verticesInput(i, 0), verticesInput(i, 1), verticesInput(i, 2);

  input_faces.resize(facetsInput.rows());
  for (size_t i = 0; i < facetsInput.rows(); i++)
    input_faces[i] << facetsInput(i, 0), facetsInput(i, 1), facetsInput(i, 2);

  ///set envelope
  igl::Timer timer;
  GEO::Mesh sf_mesh;
  json tree_with_ids;
  std::vector<std::string> meshes;

  std::cout << "================== loading mesh =====================" << std::endl;
  if (params.input_path.length() > 0) {
    if (!MeshIO::load_mesh(params.input_path, input_vertices, input_faces, sf_mesh, input_tags)) {
        logger().error("Unable to load mesh at {}", params.input_path);
        MeshIO::write_mesh(output_mesh_name, mesh, false);
        return EXIT_FAILURE;
    } else if (input_vertices.empty() || input_faces.empty()) {
        MeshIO::write_mesh(output_mesh_name, mesh, false);
        return EXIT_FAILURE;
    }
  } else {
    MeshIO::load_mesh(input_vertices, input_faces, sf_mesh, input_tags);
  }

  if (input_tags.size() != input_faces.size()) {
    input_tags.resize(input_faces.size());
    std::fill(input_tags.begin(), input_tags.end(), 0);
  }
  
  AABBWrapper tree(sf_mesh);
  if (!params.init(tree.get_sf_diag())) {
    return EXIT_FAILURE;
  }

  stats().record(StateInfo::init_id, 0, input_vertices.size(), input_faces.size(), -1, -1);

  timer.start();
  std::cout << "================== simplify =====================" << std::endl;
  simplify(input_vertices, input_faces, input_tags, tree, params, skip_simplify);
  tree.init_b_mesh_and_tree(input_vertices, input_faces, mesh);
  logger().info("preprocessing {}s", timer.getElapsedTimeInSec());
  logger().info("");
  stats().record(StateInfo::preprocessing_id, timer.getElapsedTimeInSec(), input_vertices.size(),
    input_faces.size(), -1, -1);
  if (params.log_level <= 1)
    output_component(input_vertices, input_faces, input_tags);

  timer.start();
  std::cout << "================== tetrahedralize =====================" << std::endl;
  std::vector<bool> is_face_inserted(input_faces.size(), false);
  FloatTetDelaunay::tetrahedralize(input_vertices, input_faces, tree, mesh, is_face_inserted);
  logger().info("#v = {}", mesh.get_v_num());
  logger().info("#t = {}", mesh.get_t_num());
  logger().info("tetrahedralizing {}s", timer.getElapsedTimeInSec());
  logger().info("");
  stats().record(StateInfo::tetrahedralization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(), -1, -1);

  timer.start();
  std::cout << "================== insert_triangles =====================" << std::endl;
  insert_triangles(input_vertices, input_faces, input_tags, mesh, is_face_inserted, tree, false);
  logger().info("cutting {}s", timer.getElapsedTimeInSec());
  logger().info("");
  stats().record(StateInfo::cutting_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
    mesh.get_max_energy(), mesh.get_avg_energy(),
    std::count(is_face_inserted.begin(), is_face_inserted.end(), false));

  timer.start();
  std::cout << "================== optimization =====================" << std::endl;
  optimization(input_vertices, input_faces, input_tags, is_face_inserted, mesh, tree, {{1, 1, 1, 1}});
  logger().info("mesh optimization {}s", timer.getElapsedTimeInSec());
  logger().info("");
  stats().record(StateInfo::optimization_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
    mesh.get_max_energy(), mesh.get_avg_energy());

  timer.start();
  std::cout << "================== correct_tracked_surface_orientation =====================" << std::endl;
  correct_tracked_surface_orientation(mesh, tree);
  logger().info("correct_tracked_surface_orientation done");

  if (params.smooth_open_boundary) {
    smooth_open_boundary(mesh, tree);
    for (auto &t: mesh.tets) {
      if (t.is_outside)
        t.is_removed = true;
    }
  } else {
    std::cout << "================== filter_outside =====================" << std::endl;
    elements.resize(mesh.tets.size(), 4);
    material.resize(mesh.tets.size()); material.fill(0);
    for (size_t i = 0; i < mesh.tets.size(); i++) {
      elements(i, 0) = mesh.tets[i][0];
      elements(i, 1) = mesh.tets[i][1];
      elements(i, 2) = mesh.tets[i][2];
      elements(i, 3) = mesh.tets[i][3];
      if (mesh.tets[i].is_removed) material(i) = -1;
      // std::cout << "elements.row(i)=" << elements.row(i) << std::endl;
      // std::cout << "material(i)=" << material(i) << std::endl;
    }

    if(!params.disable_filtering) {
      if(params.use_floodfill) {
        filter_outside_floodfill(mesh);
      } else if(params.use_input_for_wn){
        filter_outside(mesh, input_vertices, input_faces);
      } else
        filter_outside(mesh);
    }

    for (size_t i = 0; i < mesh.tets.size(); i++) {
      if (!mesh.tets[i].is_removed) material(i) = 1;
    }
    points.resize( mesh.tet_vertices.size(), 3);
    for (size_t i = 0; i < mesh.tet_vertices.size(); i++) {
      points(i, 0) = mesh.tet_vertices[i][0];
      points(i, 1) = mesh.tet_vertices[i][1];
      points(i, 2) = mesh.tet_vertices[i][2];
    }
  }
  
  if(params.manifold_surface){
    manifold_surface(mesh, V_sf, F_sf);
  } else {
    get_surface(mesh, V_sf, F_sf);
  }
  stats().record(StateInfo::wn_id, timer.getElapsedTimeInSec(), mesh.get_v_num(), mesh.get_t_num(),
    mesh.get_max_energy(), mesh.get_avg_energy());
  logger().info("after winding number");
  logger().info("#v = {}", mesh.get_v_num());
  logger().info("#t = {}", mesh.get_t_num());
  logger().info("winding number {}s", timer.getElapsedTimeInSec());
  logger().info("");

  //fortest
  std::vector<Scalar> colors;
  if (!nocolor) {
    colors.resize(mesh.tets.size(), -1);
    for (int i = 0; i < mesh.tets.size(); i++) {
      if (mesh.tets[i].is_removed)
        continue;
      colors[i] = mesh.tets[i].quality;
    }
  }
  //fortest
  MeshIO::write_mesh(output_mesh_name, mesh, false, colors, !nobinary, !csg_file.empty());
  igl::write_triangle_mesh(params.output_path + "_" + params.postfix + "_sf.stl", V_sf, F_sf);

  return EXIT_SUCCESS;
}