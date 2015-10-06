all: tsp-seq

tsp-seq: tsp-seq.c
	gcc tsp-seq.c -o tsp-seq -lm && ls

tsp-seq-prof: tsp-seq.c
	gcc tsp-seq.c -pg -o tsp-seq-prof -lm && gprof tsp-seq-prof gmon.out >> gmon.txt

clean:
	rm -f tsp-seq tsp-seq-prof gmon.out gmon.txt && ls
