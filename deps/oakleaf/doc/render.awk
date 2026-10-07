{
    if(/<!-- #include .* -->/) {
	split($0,A,"\"");
	while( (getline line < A[2]) > 0) {
	    print line;
	}
    }
    else {
	print;
    }
}
