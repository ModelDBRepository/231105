Doub gammln(const Doub xx);
/*
Doub factrl(const Int n) {
	static VecDoub a(171);
	static Bool init=true;
	if (init) {
		init = false;
		a[0] = 1.;
		for (Int i=1;i<171;i++) a[i] = i*a[i-1];
	}
	if (n < 0 || n > 170) throw("factrl out of range");
	return a[n];
}
Doub factln(const Int n) {
	static const Int NTOP=2000;
	static VecDoub a(NTOP);
	static Bool init=true;
	if (init) {
		init = false;
		for (Int i=0;i<NTOP;i++) a[i] = gammln(i+1.);
	}
	if (n < 0) throw("negative arg in factln");
	if (n < NTOP) return a[n];
	return gammln(n+1.);
}
Doub bico(const Int n, const Int k) {
	if (n<0 || k<0 || k>n) throw("bad args in bico");
	if (n<171) return floor(0.5+factrl(n)/(factrl(k)*factrl(n-k)));
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}
Doub beta(const Doub z, const Doub w) {
	return exp(gammln(z)+gammln(w)-gammln(z+w));
}*/
