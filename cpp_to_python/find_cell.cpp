#include <dolfin.h>
using namespace dolfin;

	size_t valid_global_index(size_t global_index){
		if (global_map.size() > global_index){
			if (global_map[global_index][0] == mpi_rank){
				return global_map[global_index][1];
			}
		} else {
			return static_cast<size_t>(-1);
		}
	}

	bool find_cell(const Point &point, Cell &cell)
	{
    	auto index_1 = valid_global_index(2 * mesh_ptr->hash(point));
		auto index_2 = valid_global_index(2 * mesh_ptr->hash(point)+1);

		/// if index is not inside local mesh, assign it with 0.
		index_1 == mesh_ptr->entity(top_dim) > index_1 ? index_1 : 0;
		index_2 == mesh_ptr->entity(top_dim) > index_2 ? index_2 : 0;

		Cell cell_1(*mesh_ptr, index_1);
    	Cell cell_2(*mesh_ptr, index_2);

		cell = cell_1.contains(point) ? cell_1 : cell_2;
		
		if (cell.contains(point)){
			return true;
		} else {
			return false;
		}
	}