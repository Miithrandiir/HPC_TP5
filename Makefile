create:
	cd build && cmake .. && make -j16 && mpirun -n 6 ./TP5.3