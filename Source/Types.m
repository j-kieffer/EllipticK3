
/* Elliptic K3 surface with parameters */

declare type EllK3;

declare attributes EllK3: BaseField,
	EllParam,
	Coeffs, //a2, a4, a6
	RedFib,
	CmpGroup, //Abstract finite group
	NSLat, //Even positive-definite lattice
	RootSublattice,
	RootQuotient, //Quotient L^*/L where L is the root sublattice
	MW;

/* Reducible fiber */

declare type EllK3RedFib;

declare attributes EllK3RedFib: Pl,
	Kodaira,
	RootType,
	RatComps,
	Group;

/* Mordell--Weil group */

declare type EllK3MW[EllK3MWElt]: Grp;

declare attributes EllK3MW: K3,
	Grp,
	Basis;

declare attributes EllK3MWElt: Parent,
	x,
	y;



