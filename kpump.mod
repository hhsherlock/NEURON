TITLE kpump

NEURON {
	SUFFIX kpump
	USEION k READ  ki, ko WRITE ik
    USEION na READ  nai	WRITE nai
	 RANGE Kp, Krest, imax
    :RANGE Krest, imax	
}
UNITS {
    (molar) = (1/liter)
    (mV) =	(millivolt)    
    (mM) =	(millimolar)
	(mA) = (milliamp)
	
}

PARAMETER {
	Kp = 0.1 (mA/cm2)
	Krest = 110 (mM)
	imax = 0.3 (mA/cm2)
	
	
}

ASSIGNED {
	
	ik (mA/cm2)
    ki       (mM)	
	ko       (mM)
	nai      (mM) 
}

BREAKPOINT {
    : nai=10
	: ik = Kp*(ki/Krest - 1)
	 ik= imax*(1/(1 + 7.3/ko)^2)*(1/(1 + 10/nai)^3) *(ki/Krest - 1)

}

