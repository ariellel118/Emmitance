#include "epot_bicgstabsolver.hpp"
#include "geometry.hpp"
#include "func_solid.hpp"
#include "epot_efield.hpp"
#include "meshvectorfield.hpp"
#include "particledatabase.hpp"
#include "geomplotter.hpp"
#include "ibsimu.hpp"
#include "error.hpp"
#include "particlediagplotter.hpp"

bool solid1( double x, double y, double z )
{
    return( x >= 0.003556 && x<=0.004064 && y>=0.0015875 );
}


bool solid2( double x, double y, double z )
{
    return( x >= 0.00762 && x <= 0.008128 && y >= 0.0015875 );
}


bool solid3( double x, double y, double z )
{
    return( x >= 0.011684 && x<= 0.012192 && y >= 0.0015875);
}


void simu( void )
{
    Geometry geom( MODE_CYL, Int3D(341,241,1), Vec3D(0,0,0), 0.00005 );

    Solid *s1 = new FuncSolid( solid1 );
    geom.set_solid( 7, s1 );
    Solid *s2 = new FuncSolid( solid2 );
    geom.set_solid( 8, s2 );
    Solid *s3 = new FuncSolid( solid3 );
    geom.set_solid( 9, s3 );
    
    geom.set_boundary( 1, Bound(BOUND_DIRICHLET,  0.0) );
    geom.set_boundary( 2, Bound(BOUND_DIRICHLET,  400.0) );
    geom.set_boundary( 3, Bound(BOUND_NEUMANN,     0.0  ) );
    geom.set_boundary( 4, Bound(BOUND_NEUMANN,     0.0  ) );
    geom.set_boundary( 7, Bound(BOUND_DIRICHLET,  400) );
    geom.set_boundary( 8, Bound(BOUND_DIRICHLET, -341) );
    geom.set_boundary( 9, Bound(BOUND_DIRICHLET, 400) );
    geom.build_mesh();

    EpotField epot( geom );
    MeshScalarField scharge( geom );
    MeshVectorField bfield;
    EpotEfield efield( epot );
    field_extrpl_e efldextrpl[6] = { FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE, 
	     FIELD_SYMMETRIC_POTENTIAL, FIELD_EXTRAPOLATE,
	          FIELD_EXTRAPOLATE, FIELD_EXTRAPOLATE };
    efield.set_extrapolation( efldextrpl );

    EpotBiCGSTABSolver solver( geom );

    ParticleDataBaseCyl pdb( geom );
    bool pmirror[6] = { false, false, true, false, false, false };
    pdb.set_mirror( pmirror );
    pdb.set_polyint(true);

    for( size_t i = 0; i < 5; i++ ) {
	solver.solve( epot, scharge );
	efield.recalculate();
	pdb.clear();
	pdb.add_2d_beam_with_energy( 1500, 5e-9, -1.0, 5.485799e-4, 
				     500, 0.0, 0.0, 
				     0.0, 0.0, 
				     0.0, 0.000012 );
	pdb.iterate_trajectories( scharge, efield, bfield );
    }

    GeomPlotter geomplotter( geom );
    geomplotter.set_size( 750, 750 );
    geomplotter.set_epot( &epot );
    geomplotter.set_particle_database( &pdb );
    geomplotter.plot_png( "plot3.png" );

    ParticleDiagPlotter partplotter(geom, pdb, AXIS_X,0.015, PARTICLE_DIAG_PLOT_SCATTER, DIAG_R, DIAG_EK);
    partplotter.set_size(750,750);
    partplotter.plot_png("plot4.png");

    /* GTKPlotter plotter( argc, argv );
    plotter.set_geometry( &geom );
    plotter.set_epot( &epot );
    plotter.set_scharge( &scharge );
    plotter.set_particledatabase( &pdb );
    plotter.new_geometry_plot_window();
    plotter.run();*/
}


int main( int argc, char **argv )
{
    try {
	ibsimu.set_message_threshold( MSG_VERBOSE, 1 );
	ibsimu.set_thread_count( 4 );
        simu();
    } catch ( Error e ) {
        e.print_error_message( std::cout );
        exit( 1 );
    }

    return( 0 );
}
