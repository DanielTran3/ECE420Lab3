all: gje_seq.c gje_parallel.c gje_parallel2.c gje_taskfunction.c
	gcc gje_seq.c -o gje_seq Lab3IO.c
	gcc gje_parallel.c -o gje_parallel -g -Wall -fopenmp Lab3IO.c
	gcc gje_parallel2.c -o gje_parallel2 -g -Wall -fopenmp Lab3IO.c
	gcc gje_taskfunction.c -o gje_taskfunction -g -Wall -fopenmp Lab3IO.c
	gcc gje_sections.c -o gje_tasksections -g -Wall -fopenmp Lab3IO.c

datagen: datagen.c
	gcc datagen.c -o datagen Lab3IO.c
gje_seq: gje_seq.c
	gcc gje_seq.c -o gje_seq Lab3IO.c
gje_parallel: gje_parallel.c
	gcc gje_parallel.c -o gje_parallel -g -Wall -fopenmp Lab3IO.c
gje_parallel2: gje_parallel2.c
	gcc gje_parallel2.c -o gje_parallel2 -g -Wall -fopenmp Lab3IO.c
gje_taskfunction: gje_taskfunction.c
	gcc gje_taskfunction.c -o gje_taskfunction -g -Wall -fopenmp Lab3IO.c
gje_sections: gje_sections.c
	gcc gje_sections.c -o gje_tasksections -g -Wall -fopenmp Lab3IO.c
