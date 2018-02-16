package ru.norgorn.fiz1;

import java.util.Arrays;

public class FzPoisson2DExplicitPsi extends FzPoisson2DExplicit {

	public FzPoisson2DExplicitPsi(Fz2Values previousValues, Fz2Values currentValues, double Pe, double Rp, double C0,
			double beta) {
		super(previousValues, currentValues, Pe, Rp, C0, beta);
	}
	
	@Override
	protected double solve(Integer j, Integer m) {
		int jj = j;
		double[][] p = previousValues.p;
		double[][] c = previousValues.c;
		double[][] k = currentValues.k;
		double kk = k[j][m];
		double cc = c[j][m];
		double dkx = firstDerX(k, j, m);
		double dkz = firstDerZ(k, j, m);
		double dcx = firstDerX(c, j, m);
		double dcz = firstDerZ(c, j, m);
		
		
		double a1 = 0;
		double a2 =  - Rp/Pe*(1/kk*dkx*cc + dcx + 1/kk*dkz*cc + dcz);
		double a3 = + (p[j+1][m]+p[j-1][m])/dx/dx;
		double a4 = + (p[j][m+1]+p[j][m-1])/dz/dz;
		double val = a1 + a2 + a3 + a4;
		
		double cur = dz*dz*dx*dx/2/(dz*dz+dx*dx) * val;
		return cur;
	}

	protected void leftBoundaryCondition() {
		currentValues.p[0] = new double[stepsZ];
		Arrays.fill(currentValues.p[0], 1);
	}
	
	protected void rightBoundaryCondition() {
		currentValues.p[lastXInd] = new double[stepsZ];
	}
	
	protected void bottomBoundaryCondition(int j) {
		double[][] c = previousValues.c;
		currentValues.p[j][0] = currentValues.p[j][1] - dz*Rp/Pe*(c[j][0]);
		if(!simpleCase){
			currentValues.p[j][0] += -dz*Rp/Pe*0*firstDerZ(c, j, lastZInd);
		}
	}
	
	protected void topBoundaryCondition(int j) {
		double[][] c = previousValues.c;
		currentValues.p[j][lastZInd] = dz*Rp/Pe*(c[j][lastZInd] + 1) - currentValues.p[j][lastZInd-1];
		if(!simpleCase){
			currentValues.p[j][lastZInd] += dz*Rp/Pe*maxZ*firstDerZ(c, j, lastZInd);
		}
	}
}
