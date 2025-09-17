#include "bidirectional_buffer.h"
#include "density_correciton.h"
#include "density_correciton.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "sphinxsys.h"

using namespace SPH;

std::string womersley_velocity_profile_csv = "./input/womersley_velocity_profile.csv";
std::string outlet_pressure_csv = "./input/outlet_pressure.csv";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 0.0058;
Real DL1 = 4 * DH;
Real DL2 = 2 * DH;
Real DL3 = 2 * DH;
Real DL = DL1 + DL2 + DL3;   
Real BW = 0.55e-3;              
Real resolution_ref = BW / 4.0; 
Real max_narrowing = 0.3;       // 0.3 0.5 0.7
Real interpolationNum = 200;
BoundingBox system_domain_bounds(Vec2d(-DL1 - 0.5 * DL2 - BW, -0.5 * DH - BW), Vec2d(0.5 * DL2 + DL3, 0.5 * DH + BW));
//----------------------------------------------------------------------
//	Material parameters of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1040.0;
Real Re = 50.0;
Real mu_f = 0.004;
Real U_f = 1.2;
Real c_f = 10.0 * U_f;
//----------------------------------------------------------------------
//	Material parameters of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1060.0;
Real Youngs_modulus1 = 0.4e6;
Real Youngs_modulus2 = 0.4e6;
Real poisson = 0.49;
//----------------------------------------------------------------------
//	Buffer parameters
//----------------------------------------------------------------------
Vec2d bidirectional_buffer_halfsize = Vec2d(2 * resolution_ref, 0.5 * DH);
Vec2d left_bidirectional_translation = Vec2d(-DL1 - 0.5 * DL2 + 2 * resolution_ref, 0);
Vec2d right_bidirectional_translation = Vec2d(0.5 * DL2 + DL3 - 2 * resolution_ref, 0);
Vec2d normal = Vec2d(1.0, 0.0);
//----------------------------------------------------------------------
//	Pressure boundary definition.
//----------------------------------------------------------------------
struct LeftInflowPressure
{
    template <class BoundaryConditionType>
    LeftInflowPressure(BoundaryConditionType &boundary_condition) {}

    Real operator()(Real p, Real physical_time)
    {
        return p;
    }
};
//struct RightInflowPressure
//{
//    template <class BoundaryConditionType>
//    RightInflowPressure(BoundaryConditionType &boundary_condition) {}
//
//    Real operator()(Real p, Real physical_time)
//    {
//        return 0.1;
//    }
//};
 class OutletPressureCSV
{
  public:
    std::vector<Real> times;
    std::vector<Real> pressures;
    Real period{0.0};

    explicit OutletPressureCSV(const std::string &csv_file)
    {
        std::ifstream file(csv_file);
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Cannot open file: " << csv_file << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::string line, tok_t, tok_p;

        // 逐行读取：两列格式为 t, p
        while (std::getline(file, line))
        {
            if (line.empty())
                continue;
            std::stringstream ss(line);

            if (!std::getline(ss, tok_t, ','))
                continue;
            if (!std::getline(ss, tok_p, ','))
                continue;

            // 跳过表头或无法转成数值的行
            try
            {
                Real t = static_cast<Real>(std::stod(tok_t));
                Real p = static_cast<Real>(std::stod(tok_p));
                times.push_back(t);
                pressures.push_back(p);
            }
            catch (...)
            {
                continue;
            }
        }

        if (!times.empty())
            period = times.back() - times.front();

        std::cout << "[INFO] Loaded " << times.size()
                  << " pressure samples." << std::endl;
    }

    static size_t nearestIndex(const std::vector<Real> &v, Real x)
    {
        auto it = std::lower_bound(v.begin(), v.end(), x);
        if (it == v.begin())
            return 0;
        if (it == v.end())
            return v.size() - 1;
        size_t hi = static_cast<size_t>(std::distance(v.begin(), it));
        size_t lo = hi - 1;
        return (std::fabs(x - v[lo]) <= std::fabs(x - v[hi])) ? lo : hi;
    }

    Real value(Real t) const
    {
        if (times.empty())
            return Real(0);

        // wrap t into one period [t0, t0+T)
        if (period > Real(0))
        {
            Real t0 = times.front();
            Real dt = t - t0;
            Real m = std::fmod(dt, period);
            if (m < Real(0))
                m += period;
            t = t0 + m;
        }

        size_t it = nearestIndex(times, t);
        return pressures[it];
    }
};

 struct RightInflowPressure
{
    OutletPressureCSV pressure;

    template <class BoundaryConditionType>
    RightInflowPressure(BoundaryConditionType &boundary_condition,
                        const std::string &csv_path = outlet_pressure_csv)
        : pressure(csv_path) {}

    Real operator()(Real p, Real physical_time)
        {
            return 1000.0 * pressure.value(physical_time);
        }
};

//----------------------------------------------------------------------
//	inflow velocity definition.
//----------------------------------------------------------------------
//struct InflowVelocity
//{
//    template <class BoundaryConditionType>
//    InflowVelocity(BoundaryConditionType &boundary_condition) {}
//
//    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
//    {
//        Vecd target_velocity = Vecd::Zero();
//
//        target_velocity[0] = 0.5;
//        target_velocity[1] = 0.0;
//
//        return target_velocity;
//    }
//};
class WomersleyProfileCSV
{
  public:
    std::vector<Real> times;                   // Nt
    std::vector<Real> radii;                   // Nr
    std::vector<std::vector<Real>> velocities; // Nt x Nr
    Real period{0.0};

    explicit WomersleyProfileCSV(const std::string &csv_file)
    {
        std::ifstream file(csv_file); // input file stream， read data from file
        if (!file.is_open())
        {
            std::cerr << "[ERROR] Cannot open file: " << csv_file << std::endl;
            std::exit(EXIT_FAILURE);
        }

        std::string line, token; // token - 最小词元

        // 第一行：NaN, r1, r2, ... 因为有NaN所以单独处理 第一行存半径
        std::getline(file, line); // read the entire line from file and save it to line
        {
            std::stringstream ss(line);   // enclose this line within a stringstream enabling each field to be retrieved sequentially as if from a stream
            std::getline(ss, token, ','); // skip first "NaN", extract character until the given character "," is found
            while (std::getline(ss, token, ','))
                radii.push_back(static_cast<Real>(std::stod(token))); // stod: string to double，save in vector radii
        }

        // 其余行：t, w(t,r1)...w(t,rNr) 其余行存的 时间 + 各半径的速度值
        while (std::getline(file, line))
        {
            if (line.empty())
                continue;               // 跳过空行
            std::stringstream ls(line); // ls - line stream 当前行字符串流

            std::getline(ls, token, ','); // time 第一个字段 = 时间
            if (token == "NaN" || token.empty())
                continue;
            times.push_back(static_cast<Real>(std::stod(token))); // 第一个字段存入times向量

            std::vector<Real> row;
            while (std::getline(ls, token, ',')) // 读后续字段
                row.push_back(static_cast<Real>(std::stod(token)));
            velocities.push_back(std::move(row)); // 存到velocities
        }

        if (!times.empty())
            period = times.back();

        std::cout << "[INFO] Loaded " << times.size()
                  << " time steps and " << radii.size() << " radii.\n";
    }

    static size_t nearestIndex(const std::vector<Real> &v, Real x)
    {
        auto it = std::lower_bound(v.begin(), v.end(), x); // return the POSITION of first element that is not less than x
        if (it == v.begin())                               // 第一个就比x大
            return 0;
        if (it == v.end()) // 最后一个都比x小
            return v.size() - 1;
        size_t hi = std::distance(v.begin(), it);
        size_t lo = hi - 1;
        return (std::fabs(x - v[lo]) <= std::fabs(x - v[hi])) ? lo : hi;
    }

    Real value(Real r, Real t) const
    {
        // clamp r
        if (r <= radii.front())
            r = radii.front();
        else if (r >= radii.back())
            r = radii.back();

        // wrap t into one period [0, T)
        if (period > 0.0)
        {
            t = std::fmod(t, period);
            if (t < 0)
                t += period;
        }

        size_t it = nearestIndex(times, t);
        size_t ir = nearestIndex(radii, r);
        return velocities[it][ir];
    }
};

struct InflowVelocity
{
    WomersleyProfileCSV profile;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition,
                   const std::string &csv_path = womersley_velocity_profile_csv)
        : profile(csv_path) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = Vecd::Zero();

        Real r = std::fabs(position[1]); // 2D：管道中心在原点
        Real speed = profile.value(r, current_time);
        target_velocity[0] = speed; // x 方向为入口轴向
        return target_velocity;
    }
};

//----------------------------------------------------------------------
//	stenosis definition
//----------------------------------------------------------------------
// @param a: maximaum narrowing, X0 - halflength of stenosis
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
//----------------------------------------------------------------------
//	Fluid body definition.
//----------------------------------------------------------------------
std::vector<Vecd> createBloodShape()
{
    // geometry
    std::vector<Vecd> blood_shape;
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, 0.5 * DH));
    blood_shape.push_back(Vecd(0.5 * DL2 + DL3, 0.5 * DH));
    blood_shape.push_back(Vecd(0.5 * DL2 + DL3, -0.5 * DH));
    blood_shape.push_back(Vecd(-DL1 - 0.5 * DL2, -0.5 * DH));
    return blood_shape;
}
class Blood : public ComplexShape
{
  public:
    explicit Blood(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon blood_outer_boundary(createBloodShape());
        add<MultiPolygonShape>(blood_outer_boundary, "BloodOuterBoundary");

        MultiPolygon stenosisUpper(createStenosisUpper(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(stenosisUpper);
        MultiPolygon stenosisLower(createStenosisLower(max_narrowing, 0.5 * DL2, DH, interpolationNum));
        subtract<MultiPolygonShape>(stenosisLower);
    }
};

//----------------------------------------------------------------------
//	Wall boundary body definition.
//----------------------------------------------------------------------
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
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
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
//----------------------------------------------------------------------
//	Prestretch
//----------------------------------------------------------------------
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
        Matd strain = 0.5 * (F.transpose() + F) - Matd::Identity();
        Matd strain_pre = 0.5 * (F_pre_.transpose() + F_pre_) - Matd::Identity();
        return lambda0_ * strain.trace() * Matd::Identity() + 2.0 * G0_ * strain +
               lambda0_ * strain_pre.trace() * Matd::Identity() + 2.0 * G0_ * strain_pre;
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
        add<LinearElasticSolidWithPrestretch>(rho0_s, Youngs_modulus1, poisson); // id = 0
        add<LinearElasticSolidWithPrestretch>(rho0_s, Youngs_modulus2, poisson); // id = 1
    };
};

class WallMaterialInitialization : public MaterialIdInitialization
{
  public:
    explicit WallMaterialInitialization(SolidBody &solid_body)
        : MaterialIdInitialization(solid_body) {};
    void update(size_t index_i, Real dt = 0.0)
    {
        Real x = pos_[index_i][0];
        Real y = pos_[index_i][1];

        if ((std::abs(x) <= (0.5 * DL2)) && (y > (outline(x, max_narrowing, 0.5 * DL2) * (DH * 0.5))) && (y < (0.5 * DH)))
        {
            material_id_[index_i] = 1;
        }
        else if ((std::abs(x) <= (0.5 * DL2)) && (y < (-outline(x, max_narrowing, 0.5 * DL2) * (DH * 0.5))) && (y > (-0.5 * DH)))
        {
            material_id_[index_i] = 1;
        }
        else
        {
            material_id_[index_i] = 0;
        }
    };
};

//----------------------------------------------------------------------
//	Constraint definition.
//----------------------------------------------------------------------
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
// For boundary_geometry
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry() {};

  private:
    bool tagManually(size_t index_i)
    {
        return 1;
    };
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem and IO environment.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setGenerateRegressionData(false);

    //sph_system.setRunParticleRelaxation(true);
    //sph_system.setReloadParticles(false);

     sph_system.setRunParticleRelaxation(false);
     sph_system.setReloadParticles(true);

    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody blood(sph_system, makeShared<Blood>("Blood"));
    blood.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(6.0);
    blood.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineAdaptationRatios(1.3, 1.0);
    wall_boundary.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(sph_system);
    //wall_boundary.defineMaterial<LinearElasticSolidWithPrestretch>(rho0_s, Youngs_modulus1, poisson);
    wall_boundary.defineMaterial<WallBoundaryComposite>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<BaseParticles, Reload>(wall_boundary.getName())
        : wall_boundary.generateParticles<BaseParticles, Lattice>();
    //wall_boundary.defineMaterial<LinearElasticSolidWithPrestretch>(rho0_s, Youngs_modulus1, poisson);
    //wall_boundary.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation blood_inner(blood);
    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelation blood_contact(blood, {&wall_boundary});
    ContactRelation wall_boundary_contact(wall_boundary, {&blood});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation blood_complex(blood_inner, blood_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
        BodyStatesRecordingToVtp write_wall_boundary(wall_boundary);
        ReloadParticleIO write_particle_reload_files(wall_boundary);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(wall_boundary_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_wall_boundary_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_wall_boundary.writeToFile();

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_wall_boundary.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        write_particle_reload_files.writeToFile();
        return 0;
    }
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<NablaWVComplex> kernel_summation(blood_inner, blood_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> wall_boundary_corrected_configuration(wall_boundary_inner);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> boundary_indicator(blood_inner, blood_contact);
    //InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(DynamicsArgs(blood_inner, 1.0), blood_contact); // for Integration1stHalfCorrectionWithWallRiemann

    SimpleDynamics<WallMaterialInitialization> composite_material_id(wall_boundary);

    Dynamics1Level<solid_dynamics::DecomposedIntegration1stHalf> wall_boundary_stress_relaxation_first_half(wall_boundary_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> wall_boundary_stress_relaxation_second_half(wall_boundary_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> wall_boundary_computing_time_step_size(wall_boundary);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(blood_inner, blood_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(blood_inner, blood_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(blood_inner, blood_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(blood_inner, blood_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(blood, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(blood);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(wall_boundary);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> wall_boundary_update_normal(wall_boundary);
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(wall_boundary_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(wall_boundary_contact);
    //----------------------------------------------------------------------
    //  Buffer
    //----------------------------------------------------------------------
    AlignedBox left_emitter_shape(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize);
    AlignedBoxByCell left_emitter(blood, left_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);

    AlignedBox right_emitter_shape(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize);
    AlignedBoxByCell right_emitter(blood, right_emitter_shape);
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::DensitySummationPressureComplex> update_fluid_density(blood_inner, blood_contact);

    SimpleDynamics<fluid_dynamics::PressureCondition<LeftInflowPressure>> left_inflow_pressure_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::PressureCondition<RightInflowPressure>> right_inflow_pressure_condition(right_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    //----------------------------------------------------------------------
    //	Constrain the Boundary.
    //----------------------------------------------------------------------
    BodyRegionByParticle wall_base(wall_boundary, makeShared<MultiPolygonShape>(conbimedConstraint()));
    SimpleDynamics<FixBodyPartConstraint> constraint_base(wall_base);
    BoundaryGeometry boundary_geometry(wall_boundary);
    SimpleDynamics<FixedInAxisDirection> constrain_axial(boundary_geometry, Vecd(0.0, 1.0)); // 1 - Limit the direction of movement
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(blood);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(blood, "Pressure");
    body_states_recording.addToWrite<int>(blood, "Indicator");
    body_states_recording.addToWrite<Real>(blood, "Density");
    body_states_recording.addToWrite<int>(blood, "BufferIndicator");
    SimpleDynamics<VonMisesStress> wall_stress(wall_boundary);
    body_states_recording.addToWrite<Real>(wall_boundary, "VonMisesStress");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    wall_boundary_normal_direction.exec();
    wall_boundary_corrected_configuration.exec();
    composite_material_id.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 3.0;     /**< End time. */
    Real Output_Time = 0.01; /**< Time stamps for output of body states. */
    Real dt = 0.0;           /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < Output_Time)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_from_fluid.exec();
            /** Update normal direction on elastic body.*/
            wall_boundary_update_normal.exec();
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                pressure_force_from_fluid.exec();

                kernel_summation.exec();
                left_inflow_pressure_condition.exec(dt);
                right_inflow_pressure_condition.exec(dt);
                inflow_velocity_condition.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = SMIN(wall_boundary_computing_time_step_size.exec(), dt - dt_s_sum);
                    wall_boundary_stress_relaxation_first_half.exec(dt_s);
                    constraint_base.exec();
                    constrain_axial.exec();
                    wall_boundary_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                inner_ite_dt++;
            }
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                            << physical_time
                            << "	Dt = " << Dt << "\n";
            }
            number_of_iterations++;

            // first do injection for all buffers
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            blood.updateCellLinkedList();
            blood_complex.updateConfiguration();
            wall_boundary.updateCellLinkedList();
            wall_boundary_contact.updateConfiguration();
            boundary_indicator.exec();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
            wall_stress.exec();
            body_states_recording.writeToFile();
        }
        TickCount t2 = TickCount::now();
        //wall_stress.exec();
        //body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    return 0;
}
