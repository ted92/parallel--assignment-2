all: tsp-seq tsp-par

tsp-seq: tsp-seq.c
	gcc tsp-seq.c -o tsp-seq -lm && ls

tsp-seq-prof: tsp-seq.c
	gcc tsp-seq.c -pg -o tsp-seq-prof -lm && ./tsp-seq-prof && gprof tsp-seq-prof gmon.out >> seq-prof.txt

tsp-par-omp: tsp-par.c
	kinst-ompp gcc -fopenmp tsp-par.c -lm -o tsp-par-omp

tsp-par: tsp-par.c
	gcc -fopenmp tsp-par.c -lm -o tsp-par

tsp-par-prof: tsp-par.c
	gcc -fopenmp tsp-par.c -pg -o tsp-par-prof -lm && ./tsp-par-prof && gprof tsp-par-prof gmon.out >> par-prof.txt

clean:
	rm -f tsp-seq tsp-seq-prof gmon.out seq-prof.txt tsp-par tsp-par-omp tsp-par-prof par-prof.txt && ls
