CC = g++ -O3

EFILE = tis.x
OBJS = tis.o energy.o param.o global.o io.o random_generator.o misc.o

tis.x: $(OBJS)	
	@echo "linking ..."	
	$(CC) -o $(EFILE) $(OBJS)

tis.o: tis.h random_generator.h energy.h param.h global.h io.h random_generator.h misc.h
	$(CC) -c tis.cpp -o tis.o

energy.o: tis.h energy.h global.h random_generator.h misc.h
	$(CC) -c energy.cpp -o energy.o

param.o: param.h global.h
	$(CC) -c param.cpp -o param.o

global.o: global.h
	$(CC) -c global.cpp -o global.o

io.o: io.h global.h
	$(CC) -c io.cpp -o io.o

random_generator.o: random_generator.h
	$(CC) -c random_generator.cpp -o random_generator.o

misc.o: misc.h
	$(CC) -c misc.cpp -o misc.o

clean:
	rm -f $(OBJS) $(EFILE)
