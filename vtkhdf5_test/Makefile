
docker_image:
	docker build -t vtkhdf_test .

docker_compile: vtk_output_temporal_cpp.hdf
	docker run --rm -v $(PWD):/src -w /src vtkhdf_test \
		bash -c "make compile"

compile: vtk_output_temporal_cpp.hdf
	mkdir -p build && \
	cd build && \
	cmake .. && \
	make && \
	./ProbeHDFVelocity ../vtk_output_temporal_cpp.hdf && \
	./ProbeHDFVelocityDynamic ../vtk_output_temporal_cpp.hdf

vtk_output_temporal_cpp.hdf:
	cp ../vtkhdf5/vtk_output_temporal_cpp.hdf .

clean:
	@rm -rf build
	@rm -rf vtk_output_temporal_cpp.hdf
	@rm -rf vtk_output_temporal_cpp_3d.hdf
	@rm -rf vtk_output_temporal_cpp_3d.hdf5
