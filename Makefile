CC= gcc
LFLAGS=  -Wall -lm -L/opt/local/lib -lgsl -lgslcblas 
CFLAGS=  -Wall   -I/opt/local/include

.o: %.c
	$(CC) $(CFLAGS) -c $< 

ffoptim:	ffoptim.o misc.o top_tools.o io.o relaxation.o
	$(CC) ffoptim.o misc.o top_tools.o io.o relaxation.o -o ffoptim.x $(LFLAGS)
clean:
	rm *.o *.x

