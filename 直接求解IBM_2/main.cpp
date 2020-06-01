#include "IBMesh.h"
#include "IBInterpolation.h"

#include <dolfin.h>
#include "TentativeVelocity.h"
#include "PressureUpdate.h"
#include "VelocityUpdate.h"
#include "Poisson.h"

using namespace dolfin;

// Define noslip domain
class NoslipDomain : public SubDomain
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        return near(x[0], 0) || near(x[0], 1.0) || near(x[1], 0);
    }
};

// Define inflow domain
class InflowDomain : public SubDomain
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        return near(x[1], 1);
    }
};

// Define inflow domain
class Pinpoint : public SubDomain
{
    bool inside(const Array<double> &x, bool on_boundary) const
    {
        return near(x[1], 0) && near(x[0], 0);
    }
};

// Define pressure boundary value at inflow
class InflowVelocity : public Expression
{
public:
    // Constructor
    InflowVelocity() : Expression(2) {}

    // Evaluate pressure at inflow
    void eval(Array<double> &values, const Array<double> &x) const
    {
        values[0] = 1.0;
        values[1] = 0.0;
    }
};
class Position : public Expression
{
public:
    // Constructor
    Position() : Expression(2) {}

    // Evaluate pressure at inflow
    void eval(Array<double> &values, const Array<double> &x) const
    {
        values[0] = x[0];
        values[1] = x[1];
    }
};

void my_move(Mesh &mesh, Function &displacement)
{
    std::vector<double> vertex_values;
    const std::size_t N = mesh.num_vertices();
    displacement.compute_vertex_values(vertex_values);

    // Move vertex coordinates
    MeshGeometry &geometry = mesh.geometry();
    std::vector<double> x(2);
    for (std::size_t i = 0; i < N; i++)
    {
        for (std::size_t j = 0; j < 2; j++)
            x[j] = vertex_values[j * N + i];
        geometry.set(i, x.data());
    }
}

int main()
{
    // Create chanel mesh
    size_t nnn = 128;
    Point point0(0, 0, 0);
    Point point1(1.0, 1.0, 0);
    IBMesh ba({point0, point1}, {nnn, nnn});

    // Create circle mesh
    auto circle = std::make_shared<Mesh>("./circle.xml.gz");
    auto U = std::make_shared<Poisson::FunctionSpace>(circle);
    auto body_velocity = std::make_shared<Function>(U);
    auto body_position = std::make_shared<Function>(U);
    auto body_disp = std::make_shared<Function>(U);

    auto position = std::make_shared<Position>();
    body_position->interpolate(*position);
    body_disp->interpolate(*position);

    // Create function spaces
    auto V = std::make_shared<VelocityUpdate::FunctionSpace>(ba.mesh());
    auto Q = std::make_shared<PressureUpdate::FunctionSpace>(ba.mesh());

    // Set parameter values
    double dt = 0.125 / static_cast<double>(nnn);
    double T = 100;

    // Define values for boundary conditions
    auto v_in = std::make_shared<InflowVelocity>();
    auto zero = std::make_shared<Constant>(0.0);
    auto zero_vector = std::make_shared<Constant>(0.0, 0.0);

    // Define subdomains for boundary conditions
    auto noslip_domain = std::make_shared<NoslipDomain>();
    auto inflow_domain = std::make_shared<InflowDomain>();
    auto pinpoint_domain = std::make_shared<Pinpoint>();

    // Define boundary conditions
    DirichletBC noslip(V, zero_vector, noslip_domain);
    DirichletBC inflow(V, v_in, inflow_domain);
    DirichletBC pinpoint(Q, zero, pinpoint_domain, "pointwise");
    std::vector<DirichletBC*> bcu = {{&inflow, &noslip}};
    std::vector<DirichletBC*> bcp = {&pinpoint};

    // Create functions
    auto u0 = std::make_shared<Function>(V);
    auto u1 = std::make_shared<Function>(V);
    auto p0 = std::make_shared<Function>(Q);
    auto p1 = std::make_shared<Function>(Q);
    auto domain_force = std::make_shared<Function>(V);
    
    // Create coefficients
    auto k = std::make_shared<Constant>(dt);
    auto f = std::make_shared<Constant>(0, 0);

    // Create forms
    TentativeVelocity::BilinearForm a1(V, V);
    TentativeVelocity::LinearForm L1(V);
    PressureUpdate::BilinearForm a2(Q, Q);
    PressureUpdate::LinearForm L2(Q);
    VelocityUpdate::BilinearForm a3(V, V);
    VelocityUpdate::LinearForm L3(V);

    // Set coefficients
    a1.k = k;
    L1.k = k;
    L1.u0 = u0;
    L1.f = f;

    L2.k = k;
    L2.u1 = u1;

    L3.k = k;
    L3.u1 = u1;
    L3.p1 = p1;

    // Assemble matrices
    Matrix A1, A2, A3;
    assemble(A1, a1);
    assemble(A2, a2);
    assemble(A3, a3);

    // Create vectors
    Vector b1, b2, b3;

  // Use amg preconditioner if available
  const std::string prec(has_krylov_solver_preconditioner("amg") ? "amg" : "default");

    // Create files for storing solution
    File ufile("results/velocity.pvd");
    File pfile("results/pressure.pvd");
    File bfile("results/body.pvd");
    File ffile("results/domain_force.pvd");


    // Time-stepping
    double t = dt;
    std::vector<double> force(V->dim());

    while (t < T + DOLFIN_EPS)
    {
        // Compute tentative velocity step
        begin("Computing tentative velocity");
        assemble(b1, L1);
        /// apply the source term.
        
        for (size_t i = 0; i < b1.size(); i++)
            b1.setitem(i, b1.getitem(i) - force[i]);
        
        for (std::size_t i = 0; i < bcu.size(); i++)
            bcu[i]->apply(A1, b1);
        
        solve(A1, *u1->vector(), b1, "gmres", "default");
        end();

        // Pressure correction
        begin("Computing pressure correction");
        assemble(b2, L2);
        for (std::size_t i = 0; i < bcp.size(); i++)
        {
            bcp[i]->apply(A2, b2);
            bcp[i]->apply(*p1->vector());
        }
        solve(A2, *p1->vector(), b2, "bicgstab", prec);
        end();

        // Velocity correction
        begin("Computing velocity correction");
        assemble(b3, L3);
        for (std::size_t i = 0; i < bcu.size(); i++)
            bcu[i]->apply(A3, b3);
        solve(A3, *u1->vector(), b3, "gmres", "default");
        end();

        // force calculation
        begin("Computing elastic force");
        auto temp_disp = std::make_shared<Function>(U);
        
        /// 1. 进行速度插值
        body_velocity->interpolate(*u1);

        /// 2. 更新坐标
        *temp_disp = FunctionAXPY(body_velocity, dt) + body_disp;
        *body_disp = *temp_disp;

        /// 3. 移动网格
        my_move(*circle, *body_position);

        /// 4. 计算出实时坐标下的高斯点、高斯权重、高斯点上的 P = mu*F*F^{T}
        std::vector<double> weights; /// shape = (1,)
        std::vector<std::vector<double>> points; /// shape = (2,)
        std::vector<std::vector<double>> values; /// shape = (2,2)
        calculate_values_at_gauss_points(*u1, weights, points, values);
        
        std::cout << weights.size() << std::endl;
        std::cout << points.size() << std::endl;
        std::cout << values.size() << std::endl;

        /*
        for (size_t i = 0; i < weights.size(); i++)
        {
            std::cout<< points[i].size()*points[i].size() << ", "<< values[i].size() <<  std::endl;
        }
        */

        /*
        std::cout<< values.size() <<std::endl;
        for (size_t i = 0; i < values.size(); i++)
        {
            std::cout<< values[i].size() <<std::endl;
            for (size_t j = 0; j < values[i].size(); j++)
            {
                std::cout<< values[i][j] <<std::endl;
            }
            
        }*/
        

        /// 5. 组装 \int mu*FF^{T}:grad(v) dx
        force = source_assemble(weights, points, values, *(u1->function_space()), ba);
        /*
        double force_max = 0.0;
        for (size_t i = 0; i < force.size(); i++)
        {
            force_max = force_max > force[i] ? force_max : force[i];
        }
        std::cout<<"force_max: "<< force_max << std::endl;
        
        std::cout << force.size() <<std::endl;
        end();
        */
       
        // Save to file
        domain_force->vector()->set_local(force);
        ffile << *domain_force;
        ufile << *u1;
        pfile << *p1;
        bfile << *body_disp;

        // Move to next time step
        *u0 = *u1;
        *p0 = *p1;
        t += dt;
        cout << "t = " << t << endl;
    }

    return 0;
}
