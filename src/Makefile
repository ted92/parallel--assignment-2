all: tsp-seq

tsp-seq: tsp-seq.c
	gcc tsp-seq.c -o tsp-seq -lm && ls

tsp-seq-prof: tsp-seq.c
	gcc tsp-seq.c -pg -o tsp-seq-prof -lm && ./tsp-seq-prof && gprof tsp-seq-prof gmon.out >> gmon.txt

tsp-par: tsp-par.c
	gcc -fopenmp tsp-par.c -lm -o tsp-par

clean:
	rm -f tsp-seq tsp-seq-prof gmon.out gmon.txt tsp-par && ls
