#ifndef _IBMesh_H_
#define _IBMesh_H_
#include <iostream>
#include <dolfin.h>
using namespace dolfin;

class IBMesh
{
public:
	IBMesh(std::array<dolfin::Point, 2> points, std::vector<size_t> dims, CellType::Type cell_type = CellType::Type::triangle)
	{

		// toplogy dimesion
		top_dim = dims.size();

		// generate mesh
		auto _mesh = RectangleMesh::create(points, {dims[0], dims[1]}, CellType::Type::triangle, "left/right");
		mesh_ptr = std::make_shared<Mesh>(_mesh);

		// check again.
		if (top_dim != 2 && top_dim != 3)
			dolfin_error("the size of dims must be 2 & 3.", ".", ".");

		nx = dims[0];
		ny = dims[1];
		nz = 0;

		x0 = points[0].x();
		x1 = points[1].x();
		y0 = points[0].y();
		y1 = points[1].y();
		z0 = 0.0;
		z1 = 0.0;

		if (top_dim == 3)
		{
			nz = dims[2];
			z0 = points[0].z();
			z1 = points[1].z();
		}
		mpi_rank = dolfin::MPI::rank(mesh_ptr->mpi_comm());
		index_mesh();
	}

	// global index and cell center
	std::array<size_t, 2> map(size_t i)
	{
		return global_map[i];
	}
	std::shared_ptr<Mesh> mesh()
	{
		return mesh_ptr;
	}

	std::vector<double> side_length()
	{
		if (top_dim != 2)
			dolfin_error("the size of dims must be 2 & 3.", "side_length()", ".");
		std::vector<double> side_lengths;
		side_lengths.push_back((x1 - x0) / nx);
		side_lengths.push_back((y1 - y0) / ny);
		return side_lengths;
	}

	void index_mesh()
	{
		// The local map local to global
		std::vector<size_t> local_map;
		for (dolfin::CellIterator e(*mesh_ptr); !e.end(); ++e)
		{
			auto center = e->midpoint();
			local_map.push_back(e->global_index());
			local_map.push_back(mpi_rank);
			local_map.push_back(e->index());
		}
		// send local map to every peocess.
		std::vector<std::vector<size_t>> mpi_collect(dolfin::MPI::size(mesh_ptr->mpi_comm()));
		dolfin::MPI::all_gather(mesh_ptr->mpi_comm(), local_map, mpi_collect);
		// alloc memory for global map.
		auto num_cell_global = mesh_ptr->num_entities_global(2);
		global_map.resize(num_cell_global);
		for (auto iter = mpi_collect.cbegin(); iter != mpi_collect.cend(); iter++)
		{
			for (auto jter = iter->begin(); jter != iter->cend();)
			{
				size_t cell_index = *jter;
				jter++;
				global_map[cell_index][0] = *jter; /// mpi_rank;size_t
				jter++;
				global_map[cell_index][1] = *jter; /// cell local index
				jter++;
			}
		}
	}

	size_t hash(const dolfin::Point &point) const
	{
		if (top_dim != 2 && top_dim != 3)
			dolfin_error("the size of dims must be 2 and 3.", ".", ".");

		double x = point.x();
		double y = point.y();
		double z = point.z();

		double dx = (x1 - x0) / static_cast<double>(nx);
		double dy = (y1 - y0) / static_cast<double>(ny);
		double dz = 0.0;

		size_t i = static_cast<size_t>((x - x0) / dx);
		size_t j = static_cast<size_t>((y - y0) / dy);
		size_t k = 0;

		if (top_dim == 3)
		{
			dz = (z1 - z0) / static_cast<double>(nz);
			k = static_cast<size_t>((z - z0) / dz);
		}

		return k * nx * ny + j * nx + i;
	}

	double x0, x1, y0, y1, z0, z1;
	size_t nx, ny, nz;
	size_t top_dim;
	size_t mpi_rank;
	// The map of global index to hash index for cells.
	std::vector<std::array<size_t, 2>> global_map;
	std::shared_ptr<Mesh> mesh_ptr;
};
#endif