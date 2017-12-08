package ru.norgorn.fiz1;

public class FzDerivatives extends FzIterative {

	protected double firstDerX(double[][] c, int j, int m) {
		if(j == 0)
			return (c[1][m]-c[0][m])/dx;
		if (j == lastXInd)
			return (c[lastXInd][m]-c[lastXInd-1][m])/dx;
		return (c[j+1][m]-c[j-1][m])/(2*dx); 
	}
	
	protected double secondDerX(double[][] c, int j, int m) {
		if(j == 0)
			return (2*c[0][m]-5*c[1][m]+4*c[2][m]-c[3][m])/dx/dx;
		if (j == lastXInd)
			return (2*c[lastXInd][m]-5*c[lastXInd-1][m]+4*c[lastXInd-2][m]-c[lastXInd-3][m])/dx/dx;
		double r = (c[j+1][m]-2*c[j][m]+c[j-1][m]) / dx/dx;
		return r;
	}
	
	protected double firstDerZ(double[][] c, int j, int m) {
		if(m == 0)
			return (c[j][1]-c[j][0])/dz;
		if (m == lastZInd)
			return (c[j][lastZInd]-c[j][lastZInd-1])/dz;
		return (c[j][m+1]-c[j][m-1])/(2*dz); 
	}
	
	protected double secondDerZ(double[][] c, int j, int m) {
		if( m ==0)
			return (2*c[j][0]-5*c[j][1]+4*c[j][2]-c[j][3])/dz/dz;
		if (m == lastZInd)
			return (2*c[j][lastZInd]-5*c[j][lastZInd-1]+4*c[j][lastZInd-2]-c[j][lastZInd-3])/dz/dz;
		double r = (c[j][m+1]-2*c[j][m]+c[j][m-1]) / dz/dz;
		return r;
	}
}
