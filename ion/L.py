def L(obj, I=0):
	# L products of acidity constants, for use in the pH calculating routine.
	# It can use ionic strength correction if an ionic strength is specified. Otherwise. It uses
	# uncorrected acidity coefficients.

	L=obj.z0;
	Ka=obj.Ka_eff(I);

	index_0=find(L==0);
	L(index_0)=1;

	if index_0~=1
		for i=(index_0-1):-1:1
			L(i)= L(i+1)*Ka(i);


	if index_0~=length(L)
		for i=(index_0+1):length(L)
			L(i)= L(i-1)/Ka(i-1);
	return L
