default: tilings

tilings: compute_tilings.c
	gcc compute_tilings.c -std=c11 -fopenmp -o tilings

clean:
	-rm -f tilings	
