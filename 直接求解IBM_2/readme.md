ffc -l dolfin ElasticStructure.ufl 
ffc -l dolfin PressureUpdate.ufl 
ffc -l dolfin VelocityUpdate.ufl 
ffc -l dolfin TentativeVelocity.ufl

mkdir build && cd build

cmake ..

make

python3 ../generate.py

./demo_navier-stokes

