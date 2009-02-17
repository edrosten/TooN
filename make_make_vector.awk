BEGIN{
	N=30
	print "//Mechanically generate code. Do not edit."
	print "//modify make_make_vector.awk and run make internal/make_vector.hh"
	for(i=1; i <= N; i++)
	{
		printf "inline Vector<" i "> makeVector(double x1"
		
		for(j=2; j<= i; j++)
			printf ", double x"j

		print "){"
		print  "  Vector<"i"> vals;"
		for(j=1; j<= i; j++)
			print "  vals["j-1"] = x"j";"

		print "  return vals;"
		print "}"
	}

}
