subsystem:
	 $(MAKE) -C src
	 mv src/runsofun .

.PHONY: clean
clean:
	rm $(EXE) src/*.o src/*.mod

