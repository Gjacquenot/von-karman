
docker_image:
	docker build -t vtkhdf_test .

docker_compile:
	docker run --rm -v $(PWD):/src -w /src vtkhdf_test \
		bash -c "make compile"

compile:
	mkdir -p build && \
	cd build && \
	cmake .. && \
	make && \
	./ProbeHDFVelocity ../vtk_output_temporal_cpp_3d.hdf

clean:
	@rm -rf build
	@rm -rf vtk_output_temporal_cpp_3d.hdf
	@rm -rf vtk_output_temporal_cpp_3d.hdf5
