#include "scattering_2d.h"
#include <algorithm>
#include <basix/finite-element.h>
#include <cmath>
#include <dolfinx.h>
#include <dolfinx/fem/Constant.h>
#include <dolfinx/fem/CoordinateElement.h>
#include <dolfinx/fem/petsc.h>
#include <dolfinx/io/VTKFile.h>
#include <dolfinx/io/XDMFFile.h>
#include <dolfinx/la/petsc.h>
#include <dolfinx/mesh/MeshTags.h>
#include <petscmat.h>
#include <petscsys.h>
#include <memory>
#include <string>
#include <utility>
#include <vector>

using namespace dolfinx;
using T = PetscScalar;
using U = typename dolfinx::scalar_value_t<T>;

namespace
{
struct PhysicalConfig
{
  double k0 = 4.0;
  double n_particle = 1.5;
  double n_background = 1.0;
  double beam_waist = 0.35;
  double beam_center_y = 0.0;
  double absorbing_alpha_factor = 1.0;

  double k_background() const { return k0 * n_background; }
  double absorbing_alpha() const
  {
    return absorbing_alpha_factor * k_background();
  }
};

enum class IncidentKind
{
  plane_wave,
  gaussian_beam
};

void assign_cell_value(fem::Function<T>& fn, const mesh::MeshTags<std::int32_t>& tags,
                       std::int32_t marker, T value);

std::shared_ptr<fem::Function<T>>
build_refractive_index_squared(const std::shared_ptr<fem::FunctionSpace<U>>& Q,
                               const mesh::MeshTags<std::int32_t>& cell_tags,
                               const PhysicalConfig& cfg)
{
  auto n2 = std::make_shared<fem::Function<T>>(Q);
  std::ranges::fill(n2->x()->array(), T(cfg.n_background * cfg.n_background));
  assign_cell_value(*n2, cell_tags, 1, T(cfg.n_particle * cfg.n_particle));
  assign_cell_value(*n2, cell_tags, 2, T(cfg.n_background * cfg.n_background));
  n2->x()->scatter_fwd();
  return n2;
}

std::shared_ptr<fem::Function<T>>
build_incident_field(const std::shared_ptr<fem::FunctionSpace<U>>& V,
                     const PhysicalConfig& cfg, IncidentKind kind)
{
  auto u_inc = std::make_shared<fem::Function<T>>(V);
  u_inc->interpolate(
      [cfg, kind](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>>
      {
        std::vector<T> values;
        values.reserve(x.extent(1));
        const double kb = cfg.k_background();
        for (std::size_t p = 0; p < x.extent(1); ++p)
        {
          const double phase = kb * x(0, p);
          if (kind == IncidentKind::plane_wave)
          {
            values.push_back(std::cos(phase));
            continue;
          }

          const double y_shift = x(1, p) - cfg.beam_center_y;
          const double envelope
              = std::exp(-(y_shift * y_shift)
                         / (cfg.beam_waist * cfg.beam_waist));
          values.push_back(envelope * std::cos(phase));
        }
        return {values, {values.size()}};
      });
  u_inc->x()->scatter_fwd();
  return u_inc;
}

std::shared_ptr<fem::Function<T>>
build_background_residual(const std::shared_ptr<fem::FunctionSpace<U>>& V,
                          const PhysicalConfig& cfg, IncidentKind kind)
{
  auto source_bg = std::make_shared<fem::Function<T>>(V);
  source_bg->interpolate(
      [cfg, kind](auto x) -> std::pair<std::vector<T>, std::vector<std::size_t>>
      {
        std::vector<T> values;
        values.reserve(x.extent(1));
        if (kind == IncidentKind::plane_wave)
        {
          values.assign(x.extent(1), T(0.0));
          return {values, {values.size()}};
        }

        const double kb = cfg.k_background();
        const double inv_w2 = 1.0 / (cfg.beam_waist * cfg.beam_waist);
        const double inv_w4 = inv_w2 * inv_w2;
        for (std::size_t p = 0; p < x.extent(1); ++p)
        {
          const double y_shift = x(1, p) - cfg.beam_center_y;
          const double envelope = std::exp(-(y_shift * y_shift) * inv_w2);
          const double envelope_dd
              = ((4.0 * y_shift * y_shift) * inv_w4 - 2.0 * inv_w2)
                * envelope;
          const double phase = kb * x(0, p);
          values.push_back(envelope_dd * std::cos(phase));
        }
        return {values, {values.size()}};
      });
  source_bg->x()->scatter_fwd();
  return source_bg;
}

std::shared_ptr<fem::Function<T>>
build_sum_field(const std::shared_ptr<fem::Function<T>>& lhs,
                const std::shared_ptr<fem::Function<T>>& rhs)
{
  auto sum = std::make_shared<fem::Function<T>>(lhs->function_space());
  auto& result = sum->x()->array();
  const auto& lhs_values = lhs->x()->array();
  const auto& rhs_values = rhs->x()->array();
  for (std::size_t i = 0; i < result.size(); ++i)
    result[i] = lhs_values[i] + rhs_values[i];
  sum->x()->scatter_fwd();
  return sum;
}

void write_field(const std::shared_ptr<fem::Function<T>>& u,
                 const std::string& filename)
{
  io::VTKFile file(MPI_COMM_WORLD, filename, "w");
  file.write<T>({*u}, 0.0);
}

void assign_cell_value(fem::Function<T>& fn, const mesh::MeshTags<std::int32_t>& tags,
                       std::int32_t marker, T value)
{
  const auto cells = tags.find(marker);
  auto dofmap = fn.function_space()->dofmap();
  for (const std::int32_t cell : cells)
  {
    for (const std::int32_t dof : dofmap->cell_dofs(cell))
      fn.x()->array()[dof] = value;
  }
}

void solve_case(const std::shared_ptr<fem::FunctionSpace<U>>& V,
                const std::shared_ptr<fem::Function<T>>& n2,
                const std::string& case_prefix, const PhysicalConfig& cfg,
                IncidentKind kind)
{
  auto k0_sq = std::make_shared<fem::Constant<T>>(cfg.k0 * cfg.k0);
  auto n_background_sq
      = std::make_shared<fem::Constant<T>>(cfg.n_background * cfg.n_background);
  // Real-valued Robin damping on the truncated outer boundary.
  auto alpha = std::make_shared<fem::Constant<T>>(cfg.absorbing_alpha());

  fem::Form<T> a = fem::create_form<T>(*form_scattering_2d_a, {V, V},
                                       {{"n2", n2}},
                                       {{"k0_sq", k0_sq}, {"alpha", alpha}},
                                       {}, {});

  auto u_inc = build_incident_field(V, cfg, kind);
  // Background residual is zero for the plane wave and nonzero for the
  // Gaussian beam because the latter is only an approximate background mode.
  auto source_bg = build_background_residual(V, cfg, kind);
  fem::Form<T> L = fem::create_form<T>(*form_scattering_2d_L, {V},
                                       {{"n2", n2},
                                        {"u_inc", u_inc},
                                        {"source_bg", source_bg}},
                                       {{"k0_sq", k0_sq},
                                        {"n_background_sq", n_background_sq}},
                                       {}, {});

  auto u_scattered = std::make_shared<fem::Function<T>>(V);
  la::petsc::Matrix A(fem::petsc::create_matrix(a), false);
  la::Vector<T> b(L.function_spaces()[0]->dofmap()->index_map,
                  L.function_spaces()[0]->dofmap()->index_map_bs());

  MatZeroEntries(A.mat());
  fem::assemble_matrix(la::petsc::Matrix::set_block_fn(A.mat(), ADD_VALUES), a,
                       {});
  MatAssemblyBegin(A.mat(), MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A.mat(), MAT_FINAL_ASSEMBLY);

  std::ranges::fill(b.array(), 0);
  fem::assemble_vector(b.array(), L);
  b.scatter_rev(std::plus<T>());

  la::petsc::KrylovSolver lu(MPI_COMM_WORLD);
  la::petsc::options::set("ksp_type", "preonly");
  la::petsc::options::set("pc_type", "lu");
  lu.set_from_options();
  lu.set_operator(A.mat());

  la::petsc::Vector wrapped_u(
      la::petsc::create_vector_wrap(*u_scattered->x()), false);
  la::petsc::Vector wrapped_b(la::petsc::create_vector_wrap(b), false);
  lu.solve(wrapped_u.vec(), wrapped_b.vec());
  u_scattered->x()->scatter_fwd();

  auto u_total = build_sum_field(u_scattered, u_inc);

  write_field(u_total, case_prefix + "_total.pvd");
  write_field(u_inc, case_prefix + "_incident.pvd");
  write_field(u_scattered, case_prefix + "_scattered.pvd");
  write_field(source_bg, case_prefix + "_background_source.pvd");
}
} // namespace

int main(int argc, char* argv[])
{
  dolfinx::init_logging(argc, argv);
  PetscInitialize(&argc, &argv, nullptr, nullptr);

  {
    io::XDMFFile mesh_file(MPI_COMM_WORLD, MIE_MESH_XDMF, "r");
    fem::CoordinateElement<double> cmap(mesh::CellType::triangle, 1);

    auto msh = std::make_shared<mesh::Mesh<U>>(mesh_file.read_mesh(
        cmap, mesh::GhostMode::shared_facet, "circle_2d"));

    auto cell_tags = std::make_shared<mesh::MeshTags<std::int32_t>>(
        mesh_file.read_meshtags(*msh, "circle_2d_cells", std::nullopt));

    auto element_u = basix::create_element<U>(
        basix::element::family::P, basix::cell::type::triangle, 1,
        basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, false);
    auto element_q = basix::create_element<U>(
        basix::element::family::P, basix::cell::type::triangle, 0,
        basix::element::lagrange_variant::unset,
        basix::element::dpc_variant::unset, true);

    auto V = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace<U>(
        msh, std::make_shared<fem::FiniteElement<U>>(element_u)));
    auto Q = std::make_shared<fem::FunctionSpace<U>>(fem::create_functionspace<U>(
        msh, std::make_shared<fem::FiniteElement<U>>(element_q)));

    PhysicalConfig cfg;
    auto n2 = build_refractive_index_squared(Q, *cell_tags, cfg);

    write_field(n2, std::string(MIE_RESULTS_DIR) + "/material_n2_2d.pvd");

    solve_case(V, n2, std::string(MIE_RESULTS_DIR) + "/plane_wave_2d", cfg,
               IncidentKind::plane_wave);
    solve_case(V, n2, std::string(MIE_RESULTS_DIR) + "/gaussian_beam_2d", cfg,
               IncidentKind::gaussian_beam);
  }

  PetscFinalize();
  return 0;
}
