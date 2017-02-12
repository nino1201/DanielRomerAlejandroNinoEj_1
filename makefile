
placas: placas.pdf  
	evince placas.pdf 
placas.pdf: grafica.py data.dat	
	python grafica.py
data.dat: placas.c
	mpicc -o placas placas.c 
	mpiexec -n 4 ./placas > data.dat
ejec: 
	mpicc -o placas placas.c
	qsub submit_job.sh
clean:
	rm -f data.dat placas placas.pdf