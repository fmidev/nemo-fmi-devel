#!/bin/awk

BEGIN{}

func relative_path() { sub(/.*NEMOGCM\/[[:upper:]]*\//,//); return substr($0,2) }

func array_loop(action,nb,array) {
    if ( nb != 0 ) {
	printf("\n%-10s: %3d\n\t",action,nb)
	j = 0
	for (i in array) {
	    j += 1
	    printf("%-45s\t",array[i])
	    if ( j % 3 == 0 && j != nb ) { printf"\n\t" }
	}
    }
}

/^A[[:upper:] ]/ { acount += 1; aroutines[acount] = relative_path() }
/^C[[:upper:] ]/ { ccount += 1; croutines[ccount] = relative_path() }
/^D[[:upper:] ]/ { dcount += 1; droutines[dcount] = relative_path() }
/^E[[:upper:] ]/ { ecount += 1; eroutines[ecount] = relative_path() }
/^G[[:upper:] ]/ { gcount += 1; groutines[gcount] = relative_path() }
/^U[[:upper:] ]/ { ucount += 1; uroutines[ucount] = relative_path() }

END{
    array_loop("added",     acount,aroutines)
    array_loop("conflicted",ccount,croutines)
    array_loop("deleted",   dcount,droutines)
    array_loop("edited",    ecount,eroutines)
    array_loop("merged",    gcount,groutines)
    array_loop("updated",   ucount,uroutines)
    printf "\n"
}
