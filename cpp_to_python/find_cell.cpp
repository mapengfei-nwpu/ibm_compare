#include <dolfin.h>
using namespace dolfin;

Cell find_cell(const Point &point, const IBMesh &um)
{
    auto mesh = um.mesh_ptr;

    auto index_1 = 2 * um.hash(point);
    auto index_2 = 2 * um.hash(point) + 1;
    Cell cell_1(*mesh, index_1);
    Cell cell_2(*mesh, index_2);

    return cell_1.contains(point) ? cell_1 : cell_2;
}