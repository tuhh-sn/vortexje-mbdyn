//
// Vortexje -- Simple VAWT example.
// Kinematics driven by MBDyn, http://www.aero.polimi.it/mbdyn/
//
// Copyright (C) 2012 - 2014 Baayen & Heinz GmbH.
//
// Authors: Jorn Baayen <jorn.baayen@baayen-heinz.com>
// Editors: Sönke Neumann <soenke.neumann@tuhh.de>
//

#include <vortexje/solver.hpp>
#include <vortexje/lifting-surface-builder.hpp>
#include <vortexje/shape-generators/airfoils/naca4-airfoil-generator.hpp>
#include <vortexje/shape-generators/ellipse-generator.hpp>
#include <vortexje/surface-writers/vtk-surface-writer.hpp>
#include <vortexje/boundary-layers/dummy-boundary-layer.hpp>
#include <vortexje/empirical-wakes/ramasamy-leishman-wake.hpp>

#include <iostream>
#include <fstream>

#include <mbdiface.h>

using namespace std;
using namespace Eigen;
using namespace Vortexje;

typedef Map<Matrix<double,Dynamic,Dynamic> > MapMatDynamic;
typedef Map<Matrix<double,Dynamic,1> > MapVecDynamic;


#define N_BLADES        2
#define MILL_RADIUS     2.5
#define TIP_SPEED_RATIO 5
#define WIND_VELOCITY   6
//#define INCLUDE_TOWER
#define END_TIME 0.5

class Blade : public LiftingSurface
{
public:
    // Constructor:
    Blade()
    {
        // Create blade:
        LiftingSurfaceBuilder surface_builder(*this);
        
        const int n_points_per_airfoil = 32;
        const int n_airfoils = 21;
        
        const double chord = 0.75;
        const double span = 4.5;
        
        int trailing_edge_point_id;
        vector<int> prev_airfoil_nodes;
        
        vector<vector<int> > node_strips;
        vector<vector<int> > panel_strips;
        
        for (int i = 0; i < n_airfoils; i++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > airfoil_points =
                NACA4AirfoilGenerator::generate(0, 0, 0.12, true, chord, n_points_per_airfoil, trailing_edge_point_id);
            for (int j = 0; j < (int) airfoil_points.size(); j++)
                airfoil_points[j](2) += i * span / (double) (n_airfoils - 1);
                
            vector<int> airfoil_nodes = surface_builder.create_nodes_for_points(airfoil_points);
            node_strips.push_back(airfoil_nodes);
            
            if (i > 0) {
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(airfoil_nodes, prev_airfoil_nodes, trailing_edge_point_id);
                panel_strips.push_back(airfoil_panels);
            }
                
            prev_airfoil_nodes = airfoil_nodes;
        }

        surface_builder.finish(node_strips, panel_strips, trailing_edge_point_id);
        
        // Translate and rotate into the canonical coordinate system:
        Vector3d translation(-chord / 3.0, 0.0, -span / 2.0);
        translate(translation);
        
        rotate(Vector3d::UnitZ(), -M_PI / 2.0);
    }
};

class Tower : public Surface
{
public:
    // Constructor:
    Tower()
    {
        // Create cylinder:      
        SurfaceBuilder surface_builder(*this);
        
        const double r = 0.1;
        const double h = 4.5;
        
        const int n_points = 32;
        const int n_layers = 21;
        
        vector<int> prev_nodes;
        
        for (int i = 0; i < n_layers; i++) {
            vector<Vector3d, Eigen::aligned_allocator<Vector3d> > points =
                EllipseGenerator::generate(r, r, n_points);
            for (int j = 0; j < (int) points.size(); j++)
                points[j](2) += i * h / (double) (n_layers - 1);
                 
            vector<int> nodes = surface_builder.create_nodes_for_points(points);
            
            if (i > 0)
                vector<int> airfoil_panels = surface_builder.create_panels_between_shapes(nodes, prev_nodes);
                
            prev_nodes = nodes;
        }

        surface_builder.finish();
        
        // Translate into the canonical coordinate system:
        Vector3d translation(0.0, 0.0, -h / 2.0);
        translate(translation);
    }
};

class VAWT : public Body
{
public:
    double rotor_radius;
    
    // Constructor:
    VAWT(string   id,
         double   rotor_radius,
         int      n_blades,
         Vector3d position,
         Vector3d velocity,
         double theta_0,
         Vector3d dthetadt) :
         Body(id), rotor_radius(rotor_radius)
    {
        // Initialize kinematics:
        this->position = position;
        this->velocity = velocity;
        this->attitude = AngleAxis<double>(theta_0, Vector3d::UnitZ());
        this->rotational_velocity = dthetadt;
        
#ifdef INCLUDE_TOWER
        // Initialize tower:
        Surface *tower = new Tower();
        add_non_lifting_surface(*tower);
        
        allocated_surfaces.push_back(tower);
#endif
        
        // Initialize blades:
        for (int i = 0; i < n_blades; i++) {
            Blade *blade = new Blade();
            
            Vector3d translation(rotor_radius, 0, 0);
            blade->translate(translation);
            
            double theta_blade = theta_0 + 2 * M_PI / n_blades * i;
            blade->rotate(Vector3d::UnitZ(), theta_blade);
            
            blade->translate(position);
            
            BoundaryLayer *boundary_layer = new DummyBoundaryLayer();

            Wake *wake = new RamasamyLeishmanWake(*blade);

            add_lifting_surface(*blade, *boundary_layer, *wake);

            allocated_boundary_layers.push_back(boundary_layer);
            allocated_surfaces.push_back(blade);
            allocated_surfaces.push_back(wake);
        }
    }
    
};

int
main (int argc, char **argv)
{    
    // Set simulation parameters:
    Parameters::convect_wake = true;

    // Init force and moment vector
    Vector3d F = Vector3d::Zero();
    Vector3d M = Vector3d::Zero();
    Vector3d F0, M0, position0;

    // Residuum vector for forces and position
    VectorXd r_k(9);
    double conv_tolerance = 1e-3;   // convergence tolerance
    double relax = 1.0;             // relaxation factor w, 0 < w <= 1

    /*** Connect to MBDyn: ***/
    cout << "Connecting to MBDyn..." << endl;
    // Init
    mbdInterface mbd1;
    mbd1.nodes = 1;
    int socketstate = 0;
    mbc_nodal_t* mbc = mbd1.Init(); // mbc node handle
    cout << "connected." << endl;

    // Send initial forces and moments
    mbd1.setForce(F,M,0);
    
    // Set up VAWT:
    // First, get position from MBDyn:
    mbd1.getMotion();                           // update MBDyn state vectors
    MapVecDynamic position(MBC_N_X(mbc),3);     // 3D vector mapping for position of the 1st node
    MapVecDynamic velocity(MBC_N_XP(mbc),3);    // Velocity
    Map<Matrix<double,3,3> > R1(MBC_N_R(mbc));  // don't use MapType because we need fixed size 3d matrix
    MapVecDynamic Omega(MBC_N_OMEGA(mbc),3);    // Rotational velocity
    AngleAxisd Theta;                           // Create theta in angle-axis definitions
    Theta = R1;                                 // explicit conversion from matrix-type rotation
//    cout << Theta.axis() << endl;             // z-axis in this case

    
    VAWT vawt(string("vawt"),
              MILL_RADIUS,
              N_BLADES,
              position,
              velocity,
              Theta.angle(),    //  the angle of the initial z-rotation
              Omega);
    
    // Set up solver:
    Solver solver("vawt-mbdyn-log");
    solver.add_body(vawt);
    
    Vector3d freestream_velocity(WIND_VELOCITY, 0, 0);
    solver.set_freestream_velocity(freestream_velocity);
    
    double fluid_density = 1.2;
    solver.set_fluid_density(fluid_density);
    
    // Set up file format for logging:
    VTKSurfaceWriter surface_writer;
    
    // Log shaft moments:
    ofstream f;
    f.open("vawt-mbdyn-log/shaft_moment.txt");
    
    // Run simulation:
    double t = 0.0;
    double dt = 0.0033;
    int step_number = 0;
    bool conv;              // flag for convergence
    int qt;                 // convergence step counter
    
    solver.initialize_wakes(dt);
    while (t < END_TIME) {

        /*** CONVERGENCE LOOP ***/
        // Resets
        conv = 0;
        qt = 0;
        while( !conv ){
            cout << "convergence step " << qt << "..." << endl;

            // Update position and rotate blades:
            position0 = position;   // cycle old states
            mbd1.getMotion();
            vawt.set_attitude((Quaterniond) R1);
            vawt.set_position(position);
            // update velocities
            vawt.set_rotational_velocity(Omega);
            vawt.set_velocity(velocity);

            // Log shaft moment:
            F0 = F; // cycle old states
            M0 = M;
            M = solver.moment(vawt, position);
//            f << qt << "," << M(2) << endl;
            F = solver.force(vawt);
//            cout << "F: " << F << endl;
//            cout << "M: " << M << endl;

            // Convergence check
            r_k.head(3) = position - position0;
            r_k.segment(3,3) = F - F0;
            r_k.tail(3) = M - M0;
            cout << "residuum: " << r_k.norm()/sqrt(9.) << endl;
            conv = r_k.norm()/sqrt(9.) < conv_tolerance;

            // underrelaxation
            F = F0 + relax*r_k.segment(3,3);
            M = M0 + relax*r_k.tail(3);

            socketstate = mbd1.setForce(F,M,conv);        // 1, convergence achieved for this time step
            if (socketstate < 0) {
                cout << "Error while sending forces. MBDyn not connected anymore?" << endl;
                exit(EXIT_FAILURE);
            }

            // Solve:
//            solver.solve(dt);
            // For strong coupling convergence is checked here and
            // the solver called with or without propagation:
             solver.solve(dt,conv); // ---> GOTO next convergence step if 0

             qt++;
        }

        // Write forces and rotational speed
        f << M(2) << "," << Omega(2) << endl;

        // Update wakes:
        solver.update_wakes(dt);

        // Log coefficients:
        solver.log(step_number, surface_writer);
        
        // Step time:
        t += dt;
        step_number++;
    }
    
    // Close shaft log file:
    f.close();
    
    // Done:
    return 0;
}
