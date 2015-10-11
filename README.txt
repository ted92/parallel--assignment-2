Enrico Tedeschi, ete011
Assignment 2 - 3201 Parallel Programming

doc:
	-pictures used in the report
	-report -> report_assignment_ete011.pdf

src:
	Makefile: all, clean, tsp-seq, tsp-par
	route.dat: fixed coordinates for cities
	tsp-par.c: parallel version of tsp-seq
	tsp-seq.c: changed sequential version
	tsp-seq-initial-code.c: given code, intial sequential version
	
---run---
	make all
	./tsp-seq [num_cities]
	./tsp-par [num_cities][num_processor]
	
	--no arg given--
		num_cities = 14
		num_processors = 4