
#include <memory>
#include <acado_optimal_control.hpp>
#include <acado_code_generation.hpp>
#include <acado_gnuplot.hpp>
#include <acado_optimal_control.hpp>
#include <cmath>
#include <iostream>
#include <reference_trajectory/simulated_reference_trajectory.hpp>
#include <iostream>
#include <acado_toolkit.hpp>
//#include <acado_gnuplot.hpp>
//#include <acado/utils/acado_utils.hpp>
//#include <acado/matrix_vector/matrix_vector.hpp>
//#include <bits/stdc++.h>
//#include <acado_code_generation.hpp>
////#include <boost/container/vector.hpp>
//#include <boost/numeric/ublas/matrix.hpp>

int main()
{
  // Use Acado
  USING_NAMESPACE_ACADO

  //------------------------------------------Defining System States------------------------------------------------//
  DifferentialState X;       //x position of the robot
  DifferentialState Y;       // y position of the robot
  DifferentialState theta_c; // heading angle

  Control omega_1; // Wheel 1 rotational speed
  Control omega_2; // Wheel 2 rotational speed
  Control omega_3; // Wheel 3 rotational speed

  OnlineData N_phi_x; // North point cloud sum
  OnlineData S_phi_x; // South point cloud sum
  OnlineData W_phi_x; // West point cloud
  OnlineData E_phi_x; // East point cloud

  OnlineData N_phi_y; // North point cloud sum
  OnlineData S_phi_y; // South point cloud sum
  OnlineData W_phi_y; // West point cloud
  OnlineData E_phi_y; // East point cloud

  OnlineData NE_phi_x; // Northeast
  OnlineData NW_phi_x; // Northwest
  OnlineData SE_phi_x; // Southeast
  OnlineData SW_phi_x; // Southwest

  OnlineData NE_phi_y; // Northeast
  OnlineData NW_phi_y; // Northwest
  OnlineData SE_phi_y; // Southeast
  OnlineData SW_phi_y;

  OnlineData x_ref;
  OnlineData y_ref;
  OnlineData theta_ref;

  DifferentialEquation f; // the right-hand side of the model
  Function h, hN;         // costs for the control window and at the of  the control  window

  IntermediateState attractive_pf, repulsive_pf;

  //-------------------------------------------Constant robot geometry------------------------------------------------//

  double PI = 2 * acos(0.0);
  const float radi = 0.07; // radius of the mecanuum wheels (m)
  const float L = 0.25;    // front arm (m)
  const float theta = 60 * PI / 180;

  //----------------------------------------------Kinematics----------------------------------------------------------//

  const double j_inv[3][3]{// inverse jacobian matrix
                           {0, 1 / radi, -L / radi},
                           {sin(theta) / radi, -cos(theta) / radi, -L / radi},
                           {-sin(theta) / radi, -cos(theta) / radi, -L / radi}};

  const double j[3][3]{// forward jacobian matrix
                       {0, -radi / (sqrt(3)), +radi / (sqrt(3))},
                       {+2 * radi / 3, -1 * radi / 3, -radi / 3},
                       {1 * radi / 3 / L, radi / 3 / L, radi / 3 / L}};
  //------------------------------------------Solver Relevant Parameters----------------------------------------------//

  const double t_start = 0.0;      // Initial time [s]
  const double t_end = 15.0;       // Time horizon [s]
  const double dt = 0.1;           // Discretization time [s]
  const int N = round(t_end / dt); // Number of nodes
  const double w_max = 200;        // Maximal yaw rate [rad/s]
  const double max_velocity = 2;
  const double epsilon = 0.000001; // Bias to prevent division by zero.
  const double alpha = 1;          //attractive potential weight
  const double beta = 1;           // repulsive potential weight

  //------------------------------------------System Dynamics and Equations-------------------------------------------//

  f << dot(X) == j[0][0] * omega_1 + j[0][1] * omega_2 + j[0][2] * omega_3;
  f << dot(Y) == j[1][0] * omega_1 + j[1][1] * omega_2 + j[1][2] * omega_3;
  f << dot(theta_c) == j[2][0] * omega_1 + j[2][1] * omega_2 + j[2][2] * omega_3;

  //    IntermediateState V_X = j[0][0] * omega_1 + j[0][1] * omega_2 + j[0][2] * omega_3;
  //    IntermediateState V_Y = j[1][0] * omega_1 + j[1][1] * omega_2 + j[1][2] * omega_3;
  //    IntermediateState omega_c = j[2][0] * omega_1 + j[2][1] * omega_2 + j[2][2] * omega_3;

  attractive_pf = alpha * ((X - x_ref) * (X - x_ref) + (Y - y_ref) * (Y - y_ref));
  repulsive_pf = beta * (1 / ((N_phi_x - X) * (N_phi_x - X) + (N_phi_y - Y) * (N_phi_y - Y) + epsilon) + 1 / ((S_phi_x - X) * (S_phi_x - X) + (S_phi_y - Y) * (S_phi_y - Y) + epsilon) + 1 / ((W_phi_x - X) * (W_phi_x - X) + (W_phi_y - Y) * (W_phi_y - Y) + epsilon) + 1 / ((E_phi_x - X) * (E_phi_x - X) + (E_phi_y - Y) * (E_phi_y - Y) + epsilon) + 1 / ((NW_phi_x - X) * (NW_phi_x - X) + (NW_phi_y - Y) * (NW_phi_y - Y) + epsilon) + 1 / ((NE_phi_x - X) * (NE_phi_x - X) + (NE_phi_y - Y) * (NE_phi_y - Y) + epsilon) + 1 / ((SW_phi_x - X) * (SW_phi_x - X) + (SW_phi_y - Y) * (SW_phi_y - Y) + epsilon) + 1 / ((SE_phi_x - X) * (SE_phi_x - X) + (SE_phi_y - Y) * (SE_phi_y - Y) + epsilon));

  attractive_pf = alpha * ((goal_point[0] - X) * (goal_point[0] - X) + (goal_point[1] - Y) * (goal_point[1] - Y));

  //-------------------------------------------------Cost Functions----------------------------------------------------//

  // End cost vector consists of all states (no inputs at last state).
  h << X;
  h << Y;
  h << theta_c;

  h << omega_1;
  h << omega_2;
  h << omega_3;
  h << repulsive_pf;
  //    h << attractive_pf;

  hN << X;
  hN << Y;
  hN << theta_c;
  hN << repulsive_pf;
  //    hN << attractive_pf;

  // Running cost weight matrix
  BMatrix Q = eye<bool>(h.getDim());

  // End cost weight matrix
  BMatrix QN = eye<bool>(hN.getDim());

  //-----------------------------------Reference Definition (for offline optimization)---------------------------------//

  //  // Set a reference for the analysis (if CODE_GEN is false).
  DVector r(h.getDim()); // Running cost reference
  r.setZero();
  r(0) = 12;
  r(1) = 0;
  r(2) = 0.0;

  DVector rN(hN.getDim()); // End cost reference
  rN.setZero();
  rN(0) = r(0);
  rN(1) = r(1);
  rN(2) = r(2);

  //-------------------------------------------------------------------------------------------------------------------//
  // DEFINE AN OPTIMAL CONTROL PROBLEM:
  // ----------------------------------
  OCP ocp(t_start, t_end, N);
  ocp.setNOD(19);
  ocp.subjectTo(f);

  ocp.minimizeLSQ(Q, h);
  ocp.minimizeLSQEndTerm(QN, hN);
  //    ocp.subjectTo(sqrt((N_phi_x - X)*(N_phi_x - X) + (N_phi_y - Y)*(N_phi_y - Y) ) - 1.2 >= 0);
  //	ocp.subjectTo(sqrt((S_phi_x - X)*(S_phi_x - X) + (S_phi_y - Y)*(S_phi_y - Y) ) - 1.2 >= 0);
  //	ocp.subjectTo(sqrt((W_phi_x - X)*(W_phi_x - X) + (W_phi_y - Y)*(W_phi_y - Y) ) - 1.2 >= 0);
  //	ocp.subjectTo(sqrt((E_phi_x - X)*(E_phi_x - X) + (E_phi_y - Y)*(E_phi_y - Y) ) - 1.2 >= 0);
  //	ocp.subjectTo(sqrt((NW_phi_x - X)*(NW_phi_x - X) + (NW_phi_y - Y)*(NW_phi_y - Y) )- 1.2 >= 0);
  //	ocp.subjectTo(sqrt((NE_phi_x - X)*(NE_phi_x - X) + (NE_phi_y - Y)*(NE_phi_y - Y) ) - 1.2 >= 0);
  //	ocp.subjectTo(sqrt((SW_phi_x - X)*(SW_phi_x - X) + (SW_phi_y - Y)*(SW_phi_y - Y) ) - 1.2 >= 0);
  //	ocp.subjectTo(sqrt((SE_phi_x - X)*(SE_phi_x - X) + (SE_phi_y - Y)*(SE_phi_y - Y) ) - 1.2 >= 0);

  ocp.subjectTo(-30.0 <= omega_1 <= 30.0);
  ocp.subjectTo(-30.0 <= omega_2 <= 30.0);
  ocp.subjectTo(-30.0 <= omega_3 <= 30.0);
  //	ocp.subjectTo( ((X - world_1[0][0])*(X - world_1[0][0]) + (Y - world_1[0][1])*(Y - world_1[0][1]) >= 1 );

  //    repulsive_pf = beta * ( 1 / ((N_phi_x - X)*(N_phi_x - X) + (N_phi_y - Y)*(N_phi_y - Y) + epsilon )
  //                          + 1 / ((S_phi_x - X)*(S_phi_x - X) + (S_phi_y - Y)*(S_phi_y - Y) + epsilon )
  //                          + 1 / ((W_phi_x - X)*(W_phi_x - X) + (W_phi_y - Y)*(W_phi_y - Y) + epsilon )
  //                          + 1 / ((E_phi_x - X)*(E_phi_x - X) + (E_phi_y - Y)*(E_phi_y - Y) + epsilon )
  //                          + 1 / ((NW_phi_x - X)*(NW_phi_x - X) + (NW_phi_y - Y)*(NW_phi_y - Y) + epsilon )
  //                          + 1 / ((NE_phi_x - X)*(NE_phi_x - X) + (NE_phi_y - Y)*(NE_phi_y - Y) + epsilon )
  //                          + 1 / ((SW_phi_x - X)*(SW_phi_x - X) + (SW_phi_y - Y)*(SW_phi_y - Y) + epsilon )
  //                          + 1 / ((SE_phi_x - X)*(SE_phi_x - X) + (SE_phi_y - Y)*(SE_phi_y - Y) + epsilon ));
  //    StaticReferenceTrajectory reference;
  //    Controller controller;

  //	// Export the code:
  OCPexport mpc(ocp);
  mpc.set(HESSIAN_APPROXIMATION, GAUSS_NEWTON);
  mpc.set(DISCRETIZATION_TYPE, SINGLE_SHOOTING);
  mpc.set(INTEGRATOR_TYPE, INT_RK4);
  mpc.set(NUM_INTEGRATOR_STEPS, 30);
  mpc.set(CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);

  mpc.set(QP_SOLVER, QP_QPOASES);
  mpc.set(FIX_INITIAL_STATE, YES);
  mpc.set(GENERATE_TEST_FILE, YES);
  mpc.set(GENERATE_MAKE_FILE, YES);
  mpc.set(GENERATE_MATLAB_INTERFACE, YES);
  mpc.set(GENERATE_SIMULINK_INTERFACE, YES);
  //  mpc.set(CG_USE_VARIABLE_WEIGHTING_MATRIX, YES);

  if (mpc.exportCode("dataset_6") != SUCCESSFUL_RETURN)
    exit(EXIT_FAILURE);

  mpc.printDimensionsQP();
  //    alg.solve();
  //    window.plot();

  return EXIT_SUCCESS;
}