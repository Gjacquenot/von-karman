main:
	g++ -std=c++11 main.cpp vtkhdf5.cpp -o main -isystem /usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5

demo: demo_vtkhdf.cpp
	g++ -std=c++11 demo_vtkhdf.cpp -o demo -isystem /usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_cpp -lhdf5

vtk_output_temporal_cpp.hdf: demo
	./demo
	@echo "vtk_output_temporal_cpp.hdf created"

clean:
	@rm -f main
	@rm -f demo
	@rm -f *.hdf