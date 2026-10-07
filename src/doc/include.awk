{
if( /<!-- #include virtual=.* -->/ ) {
  split($0,A,"\"");
  filename = "." A[2];
  while((getline line < filename) > 0) {
      print line;
  }
}
else {
    print;
}

}
