import numpy


def calc_pH(obj, I=0):
    """Return the pH of the object.
    If an ionic strength is specified, uses the corrected acidity constants.
    This function should be used only when finding the equilibrium state.
    After that, the value should be pulled from obj.pH.

    If ionic strength does not exist, assume it is zero.
    This function is used to find the equilibrium state,
    so it cannot pull the ionic strength from the object.
    """

    # Find the order of the polynomial. This is the maximum
    # size of the list of charge states in an ion.
    MaxCol = -inf
    for i in len(range(obj.concentrations)):
        MaxCol = max(MaxCol, max(obj.ions{i}.z)-min(obj.ions{i}.z)+1);
    end

    # Set up the matrix of Ls, the multiplication
    # of acidity coefficients for each ion.
    LMat = zeros(length(obj.ions), MaxCol)

    for i = range(len(obj.ions)):
        LMat(i,1:length(obj.ions{i}.z)+1)=obj.ions{i}.L(I)

    # Construct Q vector.
    Q=1
    # Convolve each line of the L matrix.
    for j in  1:size(LMat,1):
        Q=conv(Q,LMat(j,:))

    #Convolve with water dissociation.
    Q=conv(Q, [-obj.Kw_eff(I) 0 1]);

    #Construct P matrix
    for i=1:length(obj.concentrations)
        z_list=obj.ions{i}.z0

        tmp=zeros(1,size(LMat,2))
        tmp(1:length(z_list))=z_list
        Mmod=LMat;     Mmod(i,:)=Mmod(i,:).*tmp

        Pi = 1
        for kl=1:size(Mmod,1)
            Pi=conv(Pi,Mmod(kl,:))
        end

        Pi = conv([0 1],Pi)  # Convolve with P2
        PMat(i,:) = Pi

    # Multiply P matrix by concentrations, and sum.
    C = repmat(obj.concentrations', 1, size(PMat, 2));
    P = sum(PMat.*C, 1)

    # Pad whichever is smaller, P or Q
    SizeDiff = size(Q,2)-size(PMat,2)
    if SizeDiff>0:
        P = [P,repmat(0,1,SizeDiff)];
    elif SizeDiff<0:
        Q = [Q,repmat(0,1,SizeDiff)];

    # Construct polynomial.
    Poly = zeros(1,max(length(P), length(Q)))
    Poly(1:length(P)) = Poly(1:length(P))+P
    Poly(1:length(Q)) = Poly(1:length(Q))+Q # from QMat

    Poly = fliplr(Poly)

    # Solve Polynomial for concentration
    roo = numpy.roots(Poly)
    cH = roo(imag(roo)==0 & roo>0)

    # Convert to pH. Use the activity to calculate properly.
    pH = -log10(cH*obj.H.activity_coefficient(I, 1))
    return pH
