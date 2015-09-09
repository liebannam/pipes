include Sysdep

headers = $(shell ls ../*.h)

testjunction1s: testjunction1s.o setupandrun.o mp_mat.o levmar.o mp_mat_double.o str_double.o channel.o network.o file_output.o libcla.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

optimizeit:  optimizeit.o  setupandrun.o mp_mat.o levmar.o mp_mat_double.o str_double.o channel.o network.o file_output.o libcla.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

justrunit:  justrunit.o setupandrun.o  mp_mat.o levmar.o mp_mat_double.o str_double.o channel.o network.o file_output.o	libcla.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_levmar : test_levmar.o levmar.o mp_mat.o mp_mat_double.o \
	str_double.o libcla.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

test_dgesvd : test_dgesvd.o libcla.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)


%: %.o
	$(CC) -o $@ $(CFLAGS) $^ $(LIBS)

%.o : ../%.f
	$(f77) -c $(CFAST) $<

%.o : ../%.c $(headers)
	$(cc) -c $(INC) $(cFLAGS) $<

%.o: ../%.cpp $(headers)
	$(CC) -c $(INC) $(CFLAGS) $<

%.o: ../%.cxx $(headers)
	$(CC) -c $(INC) $(CFLAGS) $<

%.o : ../%.cc $(headers)
	$(MPICC) -c $(INC) $(CFLAGS) $<

%.o : ../%.cu $(headers)
	$(NVCC) -c $(CUINC) $(CUFLAGS) $<

%.do: ../%.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGSD) $<

%.ddo: ../%.cpp $(headers)
	$(CC) -c -o $@ $(INC) $(CFLAGS) -DDDDD=1 $<

%.ddo: ../%.cc $(headers)
	$(MPICC) -c -o $@ $(INC) $(CFLAGS) -DDDDD=1 $<

%.ddo: ../%.cu $(headers)
	$(NVCC) -c -o $@ $(CUINC) $(CUFLAGS) -DDDDD=1 $<

.PHONY : always
%.ps : always
	latex ../template; dvips -o $*.ps template

.PHONY : clean
clean:
	rm -f *.o *.ddo *.do a.out
