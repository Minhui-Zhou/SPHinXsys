#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"
using namespace SPH;

Real DH = 0.0058;
Real DL1 = 4 * DH;
Real DL2 = 2 * DH;
Real DL3 = 2 * DH;
Real DL = DL1 + DL2 + DL3;
Real BW = 0.55e-3;
Real resolution_ref = BW / 4.0;
Real max_narrowing = 0.3; // 0.3 0.5 0.7
Real interpolationNum = 200;

//Real DL = 0.4;
//Real DH = 0.1;
//Real resolution_ref = DH / 20.0;

Real rho0_s = 1060.0;
Real Youngs_modulus = 0.4e6;
Real poisson = 0.45;
Real gravity_g = 1.0e-1;

//std::vector<Vecd> createBodyShape()
//{
//    // geometry
//    std::vector<Vecd> body_shape;
//    body_shape.push_back(Vecd(0.0, 0.0));
//    body_shape.push_back(Vecd(0.0, DH));
//    body_shape.push_back(Vecd(DL, DH));
//    body_shape.push_back(Vecd(DL, 0.0));
//    body_shape.push_back(Vecd(0.0, 0.0));
//    return body_shape;
//}
//
//std::vector<Vecd> createBodyShape1()
//{
//    // geometry
//    std::vector<Vecd> body_shape1;
//    body_shape1.push_back(Vecd(0.0, 2.0 * DH));
//    body_shape1.push_back(Vecd(0.0, 3.0 * DH));
//    body_shape1.push_back(Vecd(DL, 3.0 * DH));
//    body_shape1.push_back(Vecd(DL, 2.0 * DH));
//    body_shape1.push_back(Vecd(0.0, 2.0 * DH));
//    return body_shape1;
//}
//
//class TestBody : public ComplexShape
//{
//  public:
//    explicit TestBody(const std::string &shape_name) : ComplexShape(shape_name)
//    {
//        MultiPolygon test_body(createBodyShape());
//        MultiPolygon test_body_1(createBodyShape1());
//        add<MultiPolygonShape>(test_body, "TestBody");
//        add<MultiPolygonShape>(test_body_1, "TestBody1");
//    }
//};

Real outline(Real x_rel, Real a, Real X0)
{
    if (std::abs(x_rel) <= X0)
    {
        return 1.0 - 0.5 * a * (1.0 + std::cos(Pi * x_rel / X0));
    }
    else
    {
        return 1.0;
    }
}
std::vector<Vecd> createStenosisUpper(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    Real dx = (2.0 * X0) / N; // From -X0 to +X0 split into N sections
    std::vector<Vecd> stenosis_upper;
    stenosis_upper.reserve(N + 1);

    for (int i = 0; i <= N; ++i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = (D * 0.5) * y_rel;
        stenosis_upper.push_back(Vecd(x_rel, +y));
    }
    stenosis_upper.push_back(stenosis_upper.front());
    return stenosis_upper;
}

std::vector<Vecd> createStenosisLower(Real a,
                                      Real X0,
                                      Real D, // diam of normal section
                                      int N)  // numbers of iteration
{
    Real dx = (2.0 * X0) / N;
    std::vector<Vecd> stenosis_lower;
    stenosis_lower.reserve(N + 1);
    for (int i = 0; i <= N; ++i)
    {
        Real x_rel = -X0 + i * dx;
        Real y_rel = outline(x_rel, a, X0);
        Real y = -(D * 0.5) * y_rel;
        stenosis_lower.push_back(Vecd(x_rel, y));
    }
    stenosis_lower.push_back(stenosis_lower.front());
    return stenosis_lower;
}
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH + BW));
    outer_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH + BW));
    outer_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH - BW));
    outer_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    return outer_wall_shape;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    inner_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH));
    inner_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH));
    inner_wall_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH));
    inner_wall_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    return inner_wall_shape;
}
class TestBody : public ComplexShape
{
  public:
    explicit TestBody(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon boundary_outer_wall_shape(createOuterWallShape());
        add<MultiPolygonShape>(boundary_outer_wall_shape, "OuterWallBoundary");
        MultiPolygon boundary_inner_wall_shape(createInnerWallShape());
        subtract<MultiPolygonShape>(boundary_inner_wall_shape, "OuterWallBoundary");

        MultiPolygon stenosisUpper(createStenosisUpper(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        add<MultiPolygonShape>(stenosisUpper);
        MultiPolygon stenosisLower(createStenosisLower(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        add<MultiPolygonShape>(stenosisLower);
    }
};

// the width of constraint: 4 * resolution
//std::vector<Vecd> createLeftConstraint()
//{
//    std::vector<Vecd> left_constraint_shape;
//    left_constraint_shape.push_back(Vecd(0.0, 0.0));
//    left_constraint_shape.push_back(Vecd(0.0, 3.0 * DH));
//    left_constraint_shape.push_back(Vecd(resolution_ref * 4.0, 3.0 * DH));
//    left_constraint_shape.push_back(Vecd(resolution_ref * 4.0, 0.0));
//    left_constraint_shape.push_back(Vecd(0.0, 0.0));
//
//    return left_constraint_shape;
//}
//
//std::vector<Vecd> createRightConstraint()
//{
//    std::vector<Vecd> right_constraint_shape;
//    right_constraint_shape.push_back(Vecd(DL - resolution_ref * 4.0, 0.0));
//    right_constraint_shape.push_back(Vecd(DL - resolution_ref * 4.0, 3.0 * DH));
//    right_constraint_shape.push_back(Vecd(DL, 3.0 * DH));
//    right_constraint_shape.push_back(Vecd(DL, 0.0));
//    right_constraint_shape.push_back(Vecd(DL - resolution_ref * 4.0, 0.0));
//
//    return right_constraint_shape;
//}
//
//MultiPolygon conbimedConstraint()
//{
//    MultiPolygon multi_polygon;
//    multi_polygon.addAPolygon(createLeftConstraint(), ShapeBooleanOps::add);
//    multi_polygon.addAPolygon(createRightConstraint(), ShapeBooleanOps::add);
//
//    return multi_polygon;
//}

std::vector<Vecd> createLeftConstraint()
{
    std::vector<Vecd> left_constraint_shape;
    left_constraint_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));
    left_constraint_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH + BW));
    left_constraint_shape.push_back(Vecd(-DL1 - 0.5 * DL2 + BW, 0.5 * DH + BW));
    left_constraint_shape.push_back(Vecd(-DL1 - 0.5 * DL2 + BW, -0.5 * DH - BW));
    left_constraint_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH - BW));

    return left_constraint_shape;
}

std::vector<Vecd> createRightConstraint()
{
    std::vector<Vecd> right_constraint_shape;
    right_constraint_shape.push_back(Vecd(0.5 * DL2 + DL3 - BW, -0.5 * DH - BW));
    right_constraint_shape.push_back(Vecd(0.5 * DL2 + DL3 - BW, 0.5 * DH + BW));
    right_constraint_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH + BW));
    right_constraint_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH - BW));
    right_constraint_shape.push_back(Vecd(0.5 * DL2 + DL3 - BW, -0.5 * DH - BW));

    return right_constraint_shape;
}

MultiPolygon conbimedConstraint()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(createLeftConstraint(), ShapeBooleanOps::add);
    multi_polygon.addAPolygon(createRightConstraint(), ShapeBooleanOps::add);

    return multi_polygon;
}

// The following are classes with different constitutions: Linear, Saint Venat-Kirchhoff, Neo-Hookean constitutive model
class LinearElasticSolidWithPrestretch : public LinearElasticSolid
{
  public:
    explicit LinearElasticSolidWithPrestretch(Real rho0, Real youngs_modulus, Real poisson_ratio)
        : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio)
    {
        material_type_name_ = "LinearElasticSolidWithPrestretch";
        F_pre_(0, 0) = 1.5;
        F_pre_(0, 1) = 0.0;
        F_pre_(1, 0) = 0.0;
        F_pre_(1, 1) = 1.0 / 1.5;
    };
    virtual ~LinearElasticSolidWithPrestretch() {};

    virtual Matd StressPK2(Matd &F, size_t index_i) override
    {
        Matd strain = 0.5 * (F.transpose() + F) - Matd::Identity() +
                      0.5 * (F_pre_.transpose() + F_pre_) - Matd::Identity();
        return lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain;
    }

  private:
    Matd F_pre_;
};
//----------------------------------------------------------------------
//	Solid material definition.
//----------------------------------------------------------------------
 class WallBoundaryComposite : public CompositeSolid
{
  public:
    WallBoundaryComposite() : CompositeSolid(rho0_s)
    {
        add<LinearElasticSolid>(rho0_s, 0.4e6, poisson); // id = 0
        //add<LinearElasticSolid>(rho0_s, 0.4e6, poisson); // id = 1
    };
};

 class WallMaterialInitialization : public MaterialIdInitialization
{
  public:
    explicit WallMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body) {};
    void update(size_t index_i, Real dt = 0.0)
    {
        //Real x = pos_[index_i][0];
        //Real y = pos_[index_i][1];

        //if ((std::abs(x) <= (0.5 * DL2)) && (y > (outline(x, max_narrowing, 0.5 * DL2) * (DH * 0.5))) && (y < (0.5 * DH)))
        //{
        //    material_id_[index_i] = 1;
        //}
        //else if ((std::abs(x) <= (0.5 * DL2)) && (y < (-outline(x, max_narrowing, 0.5 * DL2) * (DH * 0.5))) && (y > (-0.5 * DH)))
        //{
        //    material_id_[index_i] = 1;
        //}
        //else
        //{
        //    material_id_[index_i] = 0;
        //}
        material_id_[index_i] = 0;
    };
};

// class PrestressedSolid : public SaintVenantKirchhoffSolid
//{
//   public:
//     explicit PrestressedSolid(Real rho0, Real youngs_modulus, Real poisson_ratio)
//         : SaintVenantKirchhoffSolid(rho0, youngs_modulus, poisson_ratio)
//     {
//         material_type_name_ = "PrestressedSolid";
//         F_p_inv(0, 0) = 1.0 / 1.2;
//         F_p_inv(0, 1) = 0.0;
//         F_p_inv(1, 0) = 0.0;
//         F_p_inv(1, 1) = 1.2;
//     };
//     virtual ~PrestressedSolid() {};
//
//     /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
//     virtual Matd StressPK2(Matd &F, size_t index_i) override
//     {
//         // Fe = F * Fp'
//         //Matd F_e = F * F_p_inv;
//         //Matd strain = 0.5 * (F_e.transpose() * F_e - Matd::Identity());
//         //return lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain;
//
//         Matd strain = 0.5 * (F.transpose() * F - Matd::Identity());
//         Matd pre_strain = 0.5 * (F_p_inv.transpose() * F_p_inv - Matd::Identity());
//         pre_strain -= pre_strain.trace() / Dimensions * Matd::Identity();
//         return lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain + 2.0 * G0_ * pre_strain;
//     };
//
//   private:
//     Matd F_p_inv;
// };

// class PrestressedSolidBasedOnNeoHookean : public NeoHookeanSolid
//{
//   public:
//     explicit PrestressedSolidBasedOnNeoHookean(Real rho0, Real youngs_modulus, Real poisson_ratio)
//         : NeoHookeanSolid(rho0, youngs_modulus, poisson_ratio)
//     {
//         material_type_name_ = "PrestressedSolidBasedOnNeoHookean";
//         F_p_inv(0, 0) = 1.0 / 1.2;
//         F_p_inv(0, 1) = 0.0;
//         F_p_inv(1, 0) = 0.0;
//         F_p_inv(1, 1) = 1.2;
//     };
//     virtual ~PrestressedSolidBasedOnNeoHookean() {};
//
//     /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
//     virtual Matd StressPK2(Matd &F, size_t index_i) override
//     {
//         Matd Fe = F * F_p_inv;
//         Matd right_cauchy = Fe.transpose() * Fe;
//         Real J = F.determinant();
//         return G0_ * Matd::Identity() + (lambda0_ * (J - 1.0) - G0_) * J * right_cauchy.inverse();
//     };
//     virtual Matd StressCauchy(Matd &almansi_strain, size_t particle_index_i) override
//     {
//         Matd B = (-2.0 * almansi_strain + Matd::Identity()).inverse();
//         Real J = sqrt(B.determinant());
//         Matd cauchy_stress = 0.5 * K0_ * (J - 1.0 / J) * Matd::Identity() +
//                              G0_ * pow(J, -2.0 * OneOverDimensions - 1.0) *
//                                  (B - OneOverDimensions * B.trace() * Matd::Identity());
//         return cauchy_stress;
//     };
//     /** Volumetric Kirchhoff stress from determinate */
//     virtual Real VolumetricKirchhoff(Real J) override
//     {
//         return 0.5 * K0_ * (J * J - 1);
//     };
//     /** Define the calculation of the stress matrix for postprocessing */
//     virtual std::string getRelevantStressMeasureName() override { return "Cauchy"; };
//
//   private:
//     Matd F_p_inv;
// };

//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL, -DH), Vec2d(DL * 2.0, 4.0 * DH)); // Expand the domain to prepare for stretching

    //BoundingBox system_domain_bounds(Vec2d(-DL1 - 0.5 * DL2 - BW, -0.5 * DH - BW), Vec2d(0.5 * DL2 + DL3, 0.5 * DH + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);

    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody test_body(sph_system, makeShared<TestBody>("TestBody"));
    //test_body.defineMaterial<LinearElasticSolidWithPrestretch>(rho0_s, Youngs_modulus, poisson);
    test_body.defineMaterial<WallBoundaryComposite>();
    test_body.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation test_body_inner(test_body);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(test_body, gravity); // Apply gravitational force to solids

    SimpleDynamics<WallMaterialInitialization> composite_material_id(test_body);

    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> test_body_inner_corrected_configuration(test_body_inner);

    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> stress_relaxation_first_half(test_body_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(test_body_inner);

    /** Constrain the Boundary. */
    BodyRegionByParticle fixed_test_body(test_body, makeShared<MultiPolygonShape>(conbimedConstraint()));
    SimpleDynamics<FixBodyPartConstraint> constraint_base(fixed_test_body);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(test_body);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    SimpleDynamics<VonMisesStress> wall_stress(test_body);
    write_real_body_states.addToWrite<Real>(test_body, "VonMisesStress");
    write_real_body_states.addToWrite<Matd>(test_body, "DeformationGradient");
    write_real_body_states.addToWrite<int>(test_body, "MaterialID");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    test_body_inner_corrected_configuration.exec();
    constant_gravity.exec();
    composite_material_id.exec();

    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 2.0;
    Real output_interval = end_time / 100.0;
    Real Dt = 0.1 * output_interval;
    Real dt = 0.0;

    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integral_time = 0.0;
        // integrate time (loop) until the next output time
        while (integral_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                stress_relaxation_first_half.exec(dt);
                constraint_base.exec();
                stress_relaxation_second_half.exec(dt);

                ite++;
                dt = computing_time_step_size.exec();
                relaxation_time += dt;
                integral_time += dt;
                physical_time += dt;

                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }
            }
        }

        TickCount t2 = TickCount::now();
        wall_stress.exec();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
}